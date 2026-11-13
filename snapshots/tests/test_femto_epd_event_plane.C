// ROOT5/STAR regression for the historical first-order EPD mixing-bin convention.
int test_femto_epd_event_plane(const char* mainconf);
#ifndef __CINT__
#include "FemtoEpdEventPlane.h"
#include "Rtypes.h"
#include "TMath.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace EpdBinTest {
unsigned checks=0;
void Need(bool pass,const char* label) { ++checks; if(!pass)throw std::runtime_error(label); }
// Freeze the original two-stage Float_t interface: the stored double first
// becomes a float, the caller halves it, then getPsibin receives another float.
int Original(Float_t storedPsi) {
  const Float_t argument=storedPsi/2;
  int bin=-1;
  Double_t limits[3]={-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0};
  const Float_t width=(limits[1]-limits[0])/Float_t(8);
  limits[2]=width;
  if(argument>limits[0] && argument<limits[1])
    bin=Int_t((argument-limits[0])/limits[2]);
  return bin;
}
void Check(Double_t angle) {
  Need(femto_epd::XialeiBin(angle)==Original(Float_t(angle)),
       "first-order bin differs from original stored-Float_t convention");
}
void Run() {
  const double pi=TMath::Pi();
  for(int i=0;i<8;++i) {
    const double middle=-pi+(i+.5)*(2*pi/8.);
    Need(femto_epd::XialeiBin(middle)==i,"all eight first-order angular sectors are distinct");
  }
  // Original float bin width places exactly zero in bin 3, not ideal-real bin 4.
  Need(femto_epd::XialeiBin(0.)==3,"zero preserves historical Float_t boundary arithmetic");
  const double invalid[]={-pi,pi,-2*pi,2*pi,-100.,100.,
    std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::quiet_NaN()};
  for(unsigned i=0;i<sizeof(invalid)/sizeof(invalid[0]);++i) {
    Need(femto_epd::XialeiBin(invalid[i])==-1,"invalid or open endpoints are rejected");
    Check(invalid[i]);
  }
  for(int edge=0;edge<=8;++edge) {
    const double a=-pi+edge*(2*pi/8.);
    const Float_t f=Float_t(a);
    const double neighbours[]={a,std::nextafter(a,-std::numeric_limits<double>::infinity()),
      std::nextafter(a,std::numeric_limits<double>::infinity()),double(f),
      double(std::nextafter(f,-std::numeric_limits<float>::infinity())),
      double(std::nextafter(f,std::numeric_limits<float>::infinity()))};
    for(unsigned i=0;i<sizeof(neighbours)/sizeof(neighbours[0]);++i)Check(neighbours[i]);
  }
  // A deterministic arithmetic grid exercises conversion roundoff without
  // consuming ROOT's global random stream used by event reconstruction.
  unsigned state=314159u;
  for(unsigned i=0;i<20000;++i) {
    state=1664525u*state+1013904223u;
    const double fraction=double(state)/double(std::numeric_limits<unsigned>::max());
    Check(-1.1*pi+2.2*pi*fraction);
  }
}
struct PicoFixture {
  TClonesArray events,tracks,epdHits;TClonesArray* arrays[StPicoArrays::NAllPicoArrays];StPicoDst dst;
  PicoFixture():events("StPicoEvent"),tracks("StPicoTrack"),epdHits("StPicoEpdHit") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;
    arrays[StPicoArrays::EpdHit]=&epdHits;
    new(events[0])StPicoEvent();StPicoDst::set(arrays);
  }
  ~PicoFixture(){StPicoDst::unset();}
};
struct RandomFixture {
  TRandom* previous;TRandom3 random;
  RandomFixture():previous(gRandom),random(13579){gRandom=&random;}
  ~RandomFixture(){gRandom=previous;}
};
bool SameResult(const femto_epd::Result& a,const femto_epd::Result& b) {
  return a.valid==b.valid&&a.status==b.status&&a.bin==b.bin&&a.runIndex==b.runIndex&&
    a.calibrationCent==b.calibrationCent&&a.hits==b.hits&&a.psi1==b.psi1&&
    a.rawQx==b.rawQx&&a.rawQy==b.rawQy&&a.recenteredQx==b.recenteredQx&&
    a.recenteredQy==b.recenteredQy&&a.weight==b.weight;
}
void ProviderChecks(const char* mainconf) {
  femto_epd::Settings settings;settings.Load(mainconf);
  Need(settings.enabled,"explicit calibrated EPD configuration required");
  for(int file=0;file<3;++file) {
    femto_epd::Settings bad=settings;
    if(file==0)bad.recenterSha256=std::string(64,'0');
    else if(file==1)bad.shiftSha256=std::string(64,'0');
    else bad.runListSha256=std::string(64,'0');
    bool rejected=false;try {femto_epd::Provider p;p.Init(bad,true);}
    catch(const std::exception& e){rejected=std::string(e.what()).find("SHA256 mismatch")!=std::string::npos;}
    Need(rejected,"changed recenter, shift and run-list inputs are rejected");
  }
  const std::string emptyPath=std::string(gSystem->WorkingDirectory())+Form("/tmp/epd-empty-calibration-%d.root",gSystem->GetPid());
  {TFile file(emptyPath.c_str(),"CREATE");Need(!file.IsZombie(),"create tiny missing-key calibration fixture");file.Close();}
  for(int file=0;file<2;++file) {
    femto_epd::Settings bad=settings;const std::string hash=femto_epd::Sha256File(emptyPath);
    if(file==0){bad.recenterFile=emptyPath;bad.recenterSha256=hash;}
    else {bad.shiftFile=emptyPath;bad.shiftSha256=hash;}
    bool rejected=false;try {femto_epd::Provider p;p.Init(bad,true);}
    catch(const std::exception& e){rejected=std::string(e.what()).find("missing/wrong calibration axes")!=std::string::npos;}
    Need(rejected,"valid ROOT file and SHA cannot mask missing recenter/shift keys");
  }
  femto_epd::Provider provider;provider.Init(settings,true);
  Need(!provider.Calculate(0).valid,"missing PicoDst is an invalid plane");
  PicoFixture f;StPicoEvent* event=f.dst.event();
  Need(!provider.Calculate(&f.dst).valid&&provider.Calculate(&f.dst).status==1,"missing calibration trigger invalidates plane");
  event->setTriggerId(settings.triggerIds[0]);event->setPrimaryVertexPosition(0.f,0.f,0.f);
  Need(provider.Calculate(&f.dst).status==2,"calibration event vertex selection is enforced");
  event->setPrimaryVertexPosition(settings.vertexCenterX,settings.vertexCenterY,.5*(settings.minVz+settings.maxVz));
  event->setRunId(1);bool unknown=false;
  try {provider.Calculate(&f.dst);}catch(const std::exception& e){unknown=std::string(e.what()).find("run not covered")!=std::string::npos;}
  Need(unknown,"unknown run cannot silently use another calibration bin");
  std::ifstream runs(settings.runListFile.c_str());int run=0;runs>>run;Need(run>0,"fixture selects first configured run");event->setRunId(run);
  const femto_epd::Result empty=provider.Calculate(&f.dst);
  Need(!empty.valid&&empty.status==4&&empty.runIndex==0,"known run index and too-low primary multiplicity");
  // Five primary candidates give calibration cent0 while zero fitted hits make
  // both TPC subevents empty. This path must not invent an EPD angle.
  for(int i=0;i<settings.centralityThresholds[0];++i) {
    StPicoTrack* t=new(f.tracks[i])StPicoTrack();t->setPrimaryMomentum(1.,0.,1.);
  }
  const femto_epd::Result noTpc=provider.Calculate(&f.dst);
  Need(!noTpc.valid&&noTpc.status==5&&noTpc.calibrationCent==0,"original all-plane TPC validity rejects an incomplete event");
  // Select an actually populated cent0 run from the frozen reference calibration,
  // so a missing EPD subevent can be tested after valid TPC recentering.
  std::vector<int> runIds;runIds.push_back(run);while(runs>>run)runIds.push_back(run);
  int selected=-1;TFile reference(settings.recenterFile.c_str(),"READ"),shifts(settings.shiftFile.c_str(),"READ");
  for(unsigned r=0;r<runIds.size()&&selected<0;++r) {
    bool populated=true;
    for(int harmonic=0;harmonic<2;++harmonic)for(int group=0;group<2;++group)for(int coordinate=0;coordinate<2;++coordinate) {
      TProfile2D* p=dynamic_cast<TProfile2D*>(reference.Get(Form("TPCQ%c_%c_recen%s",coordinate?'y':'x',group?'B':'A',harmonic?"_psi2":"")));
      populated=populated&&p&&p->GetBinEntries(p->GetBin(1,r+1))>0;
    }
    for(int harmonic=0;harmonic<2;++harmonic)for(int group=0;group<7;++group)for(int coordinate=0;coordinate<2;++coordinate) {
      TProfile2D* p=dynamic_cast<TProfile2D*>(reference.Get(Form("EPDQ%c_ring%d_recen%s",coordinate?'y':'x',group,harmonic?"_psi2":"")));
      populated=populated&&p&&p->GetBinEntries(p->GetBin(1,r+1))>0;
    }
    for(int coordinate=0;coordinate<2;++coordinate)for(int term=0;term<settings.shiftTerms;++term) {
      TProfile3D* p=dynamic_cast<TProfile3D*>(shifts.Get(coordinate?"EPDshiftpar_cos_ring4":"EPDshiftpar_sin_ring4"));
      populated=populated&&p&&p->GetBinEntries(p->GetBin(1,term+1,r+1))>0;
    }
    if(populated)selected=r;
  }
  Need(selected>=0,"reference calibration contains a fully populated cent0 TPC run");
  event->setRunId(runIds[selected]);
  const double pt=.5*(settings.tpcMinPt+settings.tpcMaxPt);
  for(int i=0;i<settings.centralityThresholds[0];++i) {
    StPicoTrack* t=f.dst.track(i);const int group=i<settings.centralityThresholds[0]/2?0:1;
    const double eta=.5*(settings.tpcEtaRanges[group][0]+settings.tpcEtaRanges[group][1]),phi=.21+.17*i;
    const double px=pt*std::cos(phi),py=pt*std::sin(phi),pz=pt*std::sinh(eta);
    t->setPrimaryMomentum(px,py,pz);t->setGlobalMomentum(px,py,pz);
    t->setOrigin(event->primaryVertex());t->setNHitsFit(25);t->setNHitsMax(40);
  }
  const femto_epd::Result noEpd=provider.Calculate(&f.dst);
  Need(!noEpd.valid&&noEpd.status==6&&noEpd.runIndex==selected,"missing EPD group invalidates an otherwise valid TPC event");
  for(int group=0;group<4;++group)for(int hit=0;hit<2;++hit) {
    const int row=settings.epdRowGroups[group][0]+1,tile=2*(row-1)+hit;
    new(f.epdHits[2*group+hit])StPicoEpdHit(1+2*group+hit,tile,-1,100,0,0,kFALSE,Float_t(hit?1.3:.8),kTRUE);
  }
  RandomFixture random;TRandom3 untouched(13579);
  const femto_epd::Result first=provider.Calculate(&f.dst),second=provider.Calculate(&f.dst);
  Need(first.valid&&first.hits==4&&first.weight>0,"nonempty EPD geometry exercises successful corrected plane");
  Need(SameResult(first,second),"repeated event has bitwise-identical sampled geometry and plane");
  Need(gRandom==&random.random,"provider restores caller global random pointer");
  for(int draw=0;draw<4;++draw)Need(random.random.Rndm()==untouched.Rndm(),"provider preserves caller random stream state");
  for(int draw=0;draw<37;++draw)random.random.Rndm();
  Need(SameResult(first,provider.Calculate(&f.dst)),"EP sampling is independent of unrelated global random draws");
}
}
#endif
int test_femto_epd_event_plane(const char* mainconf) {
#ifndef __CINT__
  try { EpdBinTest::checks=0;EpdBinTest::Run();
    if(mainconf&&*mainconf)EpdBinTest::ProviderChecks(mainconf);
    std::cout<<"FEMTO_EPD_EVENT_PLANE_TEST_PASS checks="<<EpdBinTest::checks<<std::endl;return 0;
  } catch(const std::exception& e) {
    std::cerr<<"FEMTO_EPD_EVENT_PLANE_TEST_FAIL after "<<EpdBinTest::checks<<" checks: "<<e.what()<<std::endl;return 1;
  }
#else
  return 1;
#endif
}
