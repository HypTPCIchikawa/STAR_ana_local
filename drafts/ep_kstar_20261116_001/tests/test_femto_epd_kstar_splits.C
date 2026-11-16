// Synthetic ROOT/STAR regression; momenta and event counts are test fixtures.
int test_femto_epd_kstar_splits(const char* ep8Mainconf,const char* ep1Mainconf);
#ifndef __CINT__
#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "StMaker.h"
#include "FemtoCandidate.h"
#include "FemtoEpdEventPlane.h"
#include "ConfigManager.h"
#include "cuts/MixingConfig.h"
// Match the existing pair-diagnostics test access without changing production API.
#define private public
#include "HistManager.h"
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#undef private
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TMemFile.h"
#include "TString.h"

namespace EpdKstarTest {
unsigned checks=0;
void Need(bool pass,const char* label) { ++checks;if(!pass)throw std::runtime_error(label); }
bool Close(double a,double b) { return std::fabs(a-b)<=1.e-10*std::max(1.,std::max(std::fabs(a),std::fabs(b))); }
const char* Suffix(int region) {
  const char* suffix[]={"","_signal","_leftSB","_rightSB"};return suffix[region];
}
std::string Name(bool mixed,int region,bool ep) {
  return std::string("hKstar")+(mixed?"ME":"SE")+(ep?"VsEpBinVsCent":"VsCent")+
    "_lambda_he4"+Suffix(region);
}
TH3D* Split(FemtoLambdaLegacy& helper,bool mixed,int region) {
  TH3D* h=dynamic_cast<TH3D*>(helper.RootHistograms()->Get(Name(mixed,region,true).c_str()));
  Need(h!=0,"EP diagnostic exists and is TH3D");return h;
}
TH2* Inclusive(FemtoLambdaLegacy& helper,bool mixed,int region) {
  TH2* h=dynamic_cast<TH2*>(helper.RootHistograms()->Get(Name(mixed,region,false).c_str()));
  Need(h!=0,"inclusive centrality reference exists");return h;
}
void SameAxis(const TAxis& a,const TAxis& b) {
  Need(a.GetNbins()==b.GetNbins(),"projection axis bin counts agree");
  for(int i=1;i<=a.GetNbins()+1;++i)Need(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),"projection axis edges agree");
}
void Closure(FemtoLambdaLegacy& helper,int classes) {
  for(int mode=0;mode<2;++mode)for(int region=0;region<4;++region) {
    TH3D* h=Split(helper,mode,region);TH2* ref=Inclusive(helper,mode,region);
    SameAxis(*h->GetXaxis(),*ref->GetXaxis());SameAxis(*h->GetZaxis(),*ref->GetYaxis());
    Need(h->GetNbinsY()==classes&&h->GetYaxis()->GetXmin()==-.5&&
         h->GetYaxis()->GetXmax()==classes-.5,"Y axis describes actual mixing classes");
    Need(h->GetEntries()==ref->GetEntries(),"all EP classes preserve inclusive entry count");
    for(int x=0;x<=h->GetNbinsX()+1;++x)for(int z=0;z<=h->GetNbinsZ()+1;++z) {
      double value=0.,variance=0.;
      for(int y=1;y<=classes;++y) {
        const int cell=h->GetBin(x,y,z);value+=h->GetBinContent(cell);
        variance+=std::pow(h->GetBinError(cell),2);
      }
      const int cell=ref->GetBin(x,z);
      Need(value==ref->GetBinContent(cell),"EP sum reproduces inclusive content including kstar/cent flows");
      Need(Close(variance,std::pow(ref->GetBinError(cell),2)),"EP sum reproduces inclusive variance");
      Need(h->GetBinContent(h->GetBin(x,0,z))==0.&&
           h->GetBinContent(h->GetBin(x,classes+1,z))==0.,"no rejected class leaks to EP flow bins");
    }
  }
}
typedef std::map<std::string,double> Entries;
Entries Snapshot(const HistManager* manager) {
  Entries result;
  for(std::map<std::string,TH1*>::const_iterator i=manager->m_histograms.begin();i!=manager->m_histograms.end();++i)
    result[i->first]=i->second->GetEntries();
  return result;
}
FemtoCandidate Lambda(int event,int variation,double mass) {
  FemtoCandidate c;c.source=kFemtoCandResonance;c.speciesKey="lambda";c.eventIndex=event;
  c.reso.dau1EventIndex=event;c.reso.dau2EventIndex=event;c.reso.dau1Index=0;c.reso.dau2Index=1;
  TLorentzVector p;p.SetXYZM(.31+.13*variation,.12-.04*variation,-.43+.17*variation,1.115683);
  c.SetP4(p);c.reso.invMass=mass;return c;
}
FemtoCandidate Nucleus(FemtoLambdaLegacy& helper,int event,int variation) {
  FemtoCandidate c;c.source=kFemtoCandTrack;c.speciesKey="he4";c.eventIndex=event;
  c.trk.trackIndex=9;c.charge=2;
  c.SetP4(helper.NuclearP4(TVector3(.22+.08*variation,-.19+.03*variation,-.51-.12*variation),3));return c;
}
double Kstar(const FemtoCandidate& a,const FemtoCandidate& b) {
  TLorentzVector p=a.P4();const TLorentzVector pair=p+b.P4();p.Boost(-pair.BoostVector());return p.P();
}
struct Oracle {
  TH3D* hist[2][4];
  explicit Oracle(FemtoLambdaLegacy& helper) {
    for(int mode=0;mode<2;++mode)for(int r=0;r<4;++r) {
      hist[mode][r]=static_cast<TH3D*>(Split(helper,mode,r)->Clone(Form("oracle_%d_%d",mode,r)));
      hist[mode][r]->SetDirectory(0);hist[mode][r]->Reset();
    }
  }
  ~Oracle(){for(int mode=0;mode<2;++mode)for(int r=0;r<4;++r)delete hist[mode][r];}
  void Fill(bool mixed,int region,double k,int ep,int cent) {
    hist[mixed][0]->Fill(k,ep,cent);if(region)hist[mixed][region]->Fill(k,ep,cent);
  }
  void Check(FemtoLambdaLegacy& helper) {
    for(int mode=0;mode<2;++mode)for(int r=0;r<4;++r) {
      const TH3D* actual=Split(helper,mode,r);const TH3D* expected=hist[mode][r];
      Need(actual->GetEntries()==expected->GetEntries(),"oracle exact routed entries");
      const int cells=(actual->GetNbinsX()+2)*(actual->GetNbinsY()+2)*(actual->GetNbinsZ()+2);
      for(int i=0;i<cells;++i) {
        Need(actual->GetBinContent(i)==expected->GetBinContent(i),"oracle exact per-class/mass/flow routing");
        Need(Close(std::pow(actual->GetBinError(i),2),std::pow(expected->GetBinError(i),2)),"oracle saved variance");
      }
    }
  }
};
void Load(const char* mainconf,int classes) {
  Need(mainconf&&*mainconf,"explicit mainconf required");
  Need(ConfigManager::GetInstance().LoadConfig(mainconf),"mainconf loads");
  Need(ConfigManager::GetInstance().GetMixingConfig().nEventPlaneBins==classes,"fixture mixing class count");
}
void DirectPairs(const char* mainconf,int classes) {
  Load(mainconf,classes);FemtoLambdaLegacy helper;
  Need(helper.Init(mainconf,std::cerr)&&helper.SpeciesKey()=="he4","He4 legacy initializes");
  Oracle oracle(helper);
  const double side=.5*(1.+helper.mOuterFactor)*helper.mWindow;
  const double mass[]={helper.mMean,helper.mMean-side,helper.mMean+side,
                       helper.mMean+(helper.mOuterFactor+1.)*helper.mWindow};
  const int expectedRegion[]={1,2,3,0};
  const double low=Split(helper,false,0)->GetXaxis()->GetXmin();
  const double high=Split(helper,false,0)->GetXaxis()->GetXmax();
  const double k[]={low-.1,low,low+.217*(high-low),high-1.e-7,high,high+.1};
  const int centralities[]={0,4,8};
  for(int ep=0;ep<classes;++ep)for(int mode=0;mode<3;++mode)for(int region=0;region<4;++region)
    for(unsigned c=0;c<sizeof(centralities)/sizeof(centralities[0]);++c)
      for(unsigned x=0;x<sizeof(k)/sizeof(k[0]);++x) {
        FemtoCandidate l=Lambda(mode==2?101:100,region,mass[region]);
        FemtoCandidate n=Nucleus(helper,mode==1?101:100,region);
        if(mode)n.trk.trackIndex=0; // Reused local daughter index is allowed across events.
        Need(helper.MassRegion(l.reso.invMass)==expectedRegion[region],"fixture covers signal/left/right/outside");
        Need(helper.FillPair(l,n,k[x],.2,centralities[c],mode!=0,0,ep),"SE and both ME directions are retained");
        oracle.Fill(mode!=0,expectedRegion[region],k[x],ep,centralities[c]);
      }
  oracle.Check(helper);Closure(helper,classes);
  FemtoCandidate l=Lambda(200,0,helper.mMean),n=Nucleus(helper,200,0);
  const int invalid[]={-1,classes};
  for(int mode=0;mode<2;++mode)for(unsigned i=0;i<sizeof(invalid)/sizeof(invalid[0]);++i) {
    n.eventIndex=mode?201:200;
    const Entries rootBefore=Snapshot(helper.RootHistograms()),nuclearBefore=Snapshot(helper.NuclearHistograms());
    bool rejected=false;try{helper.FillPair(l,n,.2,.2,4,mode!=0,0,invalid[i]);}catch(const std::exception&){rejected=true;}
    Need(rejected,"invalid diagnostic mixing class fails before physical histogram fills");
    Need(Snapshot(helper.RootHistograms())==rootBefore&&Snapshot(helper.NuclearHistograms())==nuclearBefore,
         "invalid class leaves every histogram unchanged");
  }
  TMemFile output(Form("ep-kstar-write-%d",classes),"RECREATE");
  Need(helper.Write(&output),"diagnostics write succeeds");
  for(int mode=0;mode<2;++mode)for(int r=0;r<4;++r)
    Need(dynamic_cast<TH3D*>(output.Get(Name(mode,r,true).c_str()))!=0,"all eight TH3D keys persist");
}
size_t PoolEvents(const StFemtoMaker& maker) {
  size_t result=0;
  for(std::map<Int_t,std::deque<StFemtoMaker::FemtoMixingEvent> >::const_iterator i=maker.m_mixingPool.begin();i!=maker.m_mixingPool.end();++i)
    result+=i->second.size();
  return result;
}
void MakerPairs(const char* mainconf,int classes) {
  Load(mainconf,classes);StFemtoMaker maker("ep-kstar-synthetic",0,"");
  maker.mLambdaLegacy=new FemtoLambdaLegacy();
  Need(maker.mLambdaLegacy->Init(mainconf,std::cerr),"Maker legacy fixture initializes");
  maker.m_histManager=maker.mLambdaLegacy->RootHistograms();
  maker.mLambdaEpd=new femto_epd::Provider();maker.mLambdaEpdResult=new femto_epd::Result();
  MixingConfig& mix=ConfigManager::GetInstance().GetMixingConfig();
  Need(mix.mixBothDirections&&mix.bufferSize>=3,"fixture retains both directions and at least three events");
  maker.mLambdaMaxMixEvents=mix.bufferSize;maker.mLambdaMixMinVz=0.;maker.mLambdaMixMaxVz=1.;maker.m_cent9=6;
  Oracle oracle(*maker.mLambdaLegacy);
  for(int ep=0;ep<classes;++ep) {
    std::vector<FemtoCandidate> pastLambda,pastNucleus;
    for(int step=0;step<3;++step) {
      const int event=1000+ep*3+step;
      const int definitionBin=classes==1?(step+3)%8:ep;
      maker.mLambdaEpdResult->valid=true;maker.mLambdaEpdResult->status=0;maker.mLambdaEpdResult->bin=definitionBin;
      maker.mLambdaEpdResult->psi1=-TMath::Pi()+(definitionBin+.5)*(2.*TMath::Pi()/8.);
      maker.mEventCounter=event;maker.m_eventCandidates.clear();
      const FemtoCandidate l=Lambda(event,step,maker.mLambdaLegacy->mMean),n=Nucleus(*maker.mLambdaLegacy,event,step);
      maker.m_eventCandidates["lambda"].push_back(l);maker.m_eventCandidates["he4"].push_back(n);
      const int key=maker.GetMixingBin(.5,6,-1.);
      Need(key>=0&&key/(mix.nVzBins*mix.nCentralityBins)==ep,"pool key uses actual class; EP1 folds physical bins to zero");
      const double oldME=Inclusive(*maker.mLambdaLegacy,true,0)->GetEntries();
      maker.ProcessLambdaPairs(.5,6,-1.);
      Need(Inclusive(*maker.mLambdaLegacy,true,0)->GetEntries()==oldME+2*step,
           "both ME directions mix only prior events of this class");
      oracle.Fill(false,1,Kstar(l,n),ep,6);
      for(size_t p=0;p<pastLambda.size();++p) {
        oracle.Fill(true,1,Kstar(l,pastNucleus[p]),ep,6);
        oracle.Fill(true,1,Kstar(pastLambda[p],n),ep,6);
      }
      pastLambda.push_back(l);pastNucleus.push_back(n);
    }
  }
  Need(PoolEvents(maker)==size_t(3*classes),"separate pools contain each valid event once");
  Need(maker.m_mixingPool.size()==size_t(classes),"EP8 classes stay separate; EP1 has one pool");
  oracle.Check(*maker.mLambdaLegacy);Closure(*maker.mLambdaLegacy,classes);
  const size_t oldPool=PoolEvents(maker);
  const double oldSE=Inclusive(*maker.mLambdaLegacy,false,0)->GetEntries(),oldME=Inclusive(*maker.mLambdaLegacy,true,0)->GetEntries();
  maker.mLambdaEpdResult->valid=false;maker.mLambdaEpdResult->bin=-1;maker.mLambdaEpdResult->status=6;
  const Long64_t rejected=maker.mLambdaEpdRejected;
  maker.ProcessLambdaPairs(.5,6,-1.);
  Need(maker.mLambdaEpdRejected==rejected+1,"invalid plane is explicitly rejected in EP1 and EP8");
  Need(oldSE==Inclusive(*maker.mLambdaLegacy,false,0)->GetEntries()&&
       oldME==Inclusive(*maker.mLambdaLegacy,true,0)->GetEntries()&&PoolEvents(maker)==oldPool,
       "invalid plane contributes neither pairs nor pool events");
  oracle.Check(*maker.mLambdaLegacy);Closure(*maker.mLambdaLegacy,classes);
}
}
#endif
int test_femto_epd_kstar_splits(const char* ep8Mainconf,const char* ep1Mainconf) {
#ifndef __CINT__
  try {
    EpdKstarTest::checks=0;
    EpdKstarTest::DirectPairs(ep8Mainconf,8);EpdKstarTest::MakerPairs(ep8Mainconf,8);
    EpdKstarTest::DirectPairs(ep1Mainconf,1);EpdKstarTest::MakerPairs(ep1Mainconf,1);
    std::cout<<"FEMTO_EPD_KSTAR_SPLITS_TEST_PASS checks="<<EpdKstarTest::checks<<std::endl;return 0;
  }catch(const std::exception& e) {
    std::cerr<<"FEMTO_EPD_KSTAR_SPLITS_TEST_FAIL after "<<EpdKstarTest::checks<<" checks: "<<e.what()<<std::endl;return 1;
  }
#else
  return 1;
#endif
}
