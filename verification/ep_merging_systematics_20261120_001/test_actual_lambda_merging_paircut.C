// Small deterministic integration test of the production FemtoLambdaLegacy.
// No input event loop or job submission. Caller supplies a cut-on study mainconf
// (He4/proton, strict 0.05 rectangle, writeQA=true, He3 control enabled).
int test_actual_lambda_merging_paircut(const char* mainconf);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "HistManager.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ActualLambdaMergingTest {
using femto_lambda_pair_geometry::EventGeometry;
using femto_lambda_pair_geometry::TrackPoint;
unsigned checks=0,calls=0;
void Need(bool good,const char* reason) {++checks;if(!good)throw std::runtime_error(reason);}
FemtoCandidate Lambda(int event) {
  FemtoCandidate c;c.speciesKey="lambda";c.source=kFemtoCandResonance;
  c.eventIndex=event;c.reso.dau1EventIndex=c.reso.dau2EventIndex=event;
  c.reso.dau1Index=0;c.reso.dau2Index=1;c.reso.invMass=1.11596f;
  TLorentzVector p;p.SetXYZM(.32,.14,.8,1.115683);c.SetP4(p);return c;
}
FemtoCandidate Nucleus(FemtoLambdaLegacy& helper,int event,int type,int index=2) {
  FemtoCandidate c;c.speciesKey=type==3?"he4":"he3";c.source=kFemtoCandTrack;
  c.eventIndex=event;c.trk.trackIndex=index;c.charge=2;
  c.SetP4(helper.NuclearP4(TVector3(.55,.08,.18),type));return c;
}
TrackPoint Point(double phi,double eta,bool valid=true) {
  TrackPoint p;p.phi=phi;p.eta=eta;p.valid=valid;return p;
}
EventGeometry Geometry(int event,double pphi=.01,double peta=.01,
                       double piphi=.15,double pieta=.15) {
  EventGeometry e;e.eventIndex=event;e.points[0]=Point(pphi,peta);
  e.points[1]=Point(piphi,pieta);e.points[2]=Point(0.,0.);return e;
}
TH2* Qa(FemtoLambdaLegacy& h) {
  TH2* q=dynamic_cast<TH2*>(h.NuclearHistograms()->Get("hDphiDetaPairCutQA"));
  Need(q!=0,"dedicated pair-cut counter exists");return q;
}
double Counter(FemtoLambdaLegacy& h,int x,bool mixed) {
  TH2* q=Qa(h);return q->GetBinContent(q->FindBin(x,mixed?1.:0.));
}
double MapEntries(FemtoLambdaLegacy& h,const char* daughter,bool mixed,bool after) {
  const std::string key=std::string("hDphiDetaPair_")+daughter+"_"+
      (mixed?"Mixed_":"")+"4He_"+(after?"after":"before");
  TH1* map=h.NuclearHistograms()->Get(key.c_str());Need(map!=0,"before/after daughter map exists");
  return map->GetEntries();
}
double CfEntries(FemtoLambdaLegacy& h,const char* species,bool mixed) {
  TH1* cf=h.RootHistograms()->Get((std::string("hKstar")+(mixed?"ME":"SE")+
      "_lambda_"+species+"_signal").c_str());Need(cf!=0,"signal kstar histogram exists");return cf->GetEntries();
}
void Call(FemtoLambdaLegacy& h,const FemtoCandidate& l,const FemtoCandidate& n,bool mixed,
          const EventGeometry& lg,const EventGeometry& ng,bool reject,bool protonValid=true,bool pionValid=true) {
  const double cf=CfEntries(h,"he4",mixed),considered=Counter(h,0,mixed),
    passed=Counter(h,1,mixed),rejected=Counter(h,2,mixed),
    beforeP=MapEntries(h,"proton",mixed,false),afterP=MapEntries(h,"proton",mixed,true),
    beforePi=MapEntries(h,"pion",mixed,false),afterPi=MapEntries(h,"pion",mixed,true),
    validP=Counter(h,3,mixed),invalidP=Counter(h,4,mixed),
    validPi=Counter(h,5,mixed),invalidPi=Counter(h,6,mixed),
    signalAll=Counter(h,8,mixed),signalPass=Counter(h,9,mixed),signalReject=Counter(h,10,mixed),
    signalInvalid=Counter(h,11,mixed);
  bool reason=!reject;
  const bool accepted=h.FillPair(l,n,.217,.31,4,mixed,0,&lg,&ng,&reason);++calls;
  Need(accepted==!reject&&reason==reject,"actual FillPair result and merging rejection reason");
  Need(CfEntries(h,"he4",mixed)==cf+(reject?0:1),"physical signal kstar changes only for accepted pairs");
  Need(Counter(h,0,mixed)==considered+1&&Counter(h,1,mixed)==passed+(reject?0:1)&&
       Counter(h,2,mixed)==rejected+(reject?1:0),"all-mass cut counter accepted/rejected closure");
  Need(Counter(h,8,mixed)==signalAll+1&&Counter(h,9,mixed)==signalPass+(reject?0:1)&&
       Counter(h,10,mixed)==signalReject+(reject?1:0),"signal counter accepted/rejected closure");
  Need(Counter(h,3,mixed)==validP+(protonValid?1:0)&&Counter(h,4,mixed)==invalidP+(protonValid?0:1)&&
       Counter(h,5,mixed)==validPi+(pionValid?1:0)&&Counter(h,6,mixed)==invalidPi+(pionValid?0:1)&&
       Counter(h,11,mixed)==signalInvalid+(protonValid?0:1),"valid/unreachable counters follow each daughter");
  Need(MapEntries(h,"proton",mixed,false)==beforeP+(protonValid?1:0)&&
       MapEntries(h,"proton",mixed,true)==afterP+(!reject&&protonValid?1:0)&&
       MapEntries(h,"pion",mixed,false)==beforePi+(pionValid?1:0)&&
       MapEntries(h,"pion",mixed,true)==afterPi+(!reject&&pionValid?1:0),
       "signal QA maps retain before-cut pairs and gate after-cut pairs");
}
void EqualAliases(FemtoLambdaLegacy& h) {
  const char* daughters[]={"proton","pion"};
  for(int d=0;d<2;++d) {
    const std::string base=std::string("hDphiDeta_")+daughters[d]+"_4He";
    const TH1* a=h.NuclearHistograms()->Get(base.c_str());
    const TH1* b=h.NuclearHistograms()->Get((base+"_EventField").c_str());
    Need(a&&b&&a->GetEntries()==b->GetEntries(),"legacy aliases have identical entries");
    const int cells=(a->GetNbinsX()+2)*(a->GetNbinsY()+2);
    for(int bin=0;bin<cells;++bin)Need(a->GetBinContent(bin)==b->GetBinContent(bin)&&
      a->GetBinError(bin)==b->GetBinError(bin),"legacy aliases preserve exact bins/errors including flow");
  }
}
std::vector<FemtoLambdaLegacy::TrackGeometry> TrackBank() {
  std::vector<FemtoLambdaLegacy::TrackGeometry> tracks(3);
  for(int i=0;i<3;++i) {
    tracks[i].trackIndex=i;tracks[i].charge=i==1?-1:1;
    tracks[i].origin=TVector3(0.,0.,0.);
    tracks[i].globalMomentum=i==1?TVector3(.35,-.18,.09):TVector3(.55,.08,.18);
  }
  return tracks;
}
void OwnEventFieldAndOwnership(FemtoLambdaLegacy& h) {
  const std::vector<FemtoLambdaLegacy::TrackGeometry> tracks=TrackBank();
  EventGeometry cached[2];const double fields[]={5.,-5.};
  for(int e=0;e<2;++e) {
    FemtoCandidateStore store;store["lambda"].push_back(Lambda(10+10*e));
    store["he4"].push_back(Nucleus(h,10+10*e,3));
    h.SetReplayTrackGeometry(tracks,fields[e]);h.BuildPairGeometry(10+10*e,store,0,cached[e]);
    Need(cached[e].eventIndex==10+10*e&&cached[e].points.size()==3,"cache has exact event identity and required tracks only");
    for(int i=0;i<3;++i) {
      TVector3 pos;const bool reached=FemtoLambdaLegacy::MergingTrackPositionAtRadius(tracks[i],fields[e],1.4,210.,1.e-3,pos);
      Need(reached&&cached[e].points[i].valid,"synthetic saved helix reaches common TPC cylinder");
      Need(std::fabs(TVector2::Phi_mpi_pi(cached[e].points[i].phi-pos.Phi()))<1.e-12&&
           std::fabs(cached[e].points[i].eta-tracks[i].globalMomentum.Eta())<1.e-12,
           "cached point uses own event field, common radius and global momentum eta");
    }
  }
  Need(std::fabs(TVector2::Phi_mpi_pi(cached[0].points[0].phi-cached[1].points[0].phi))>.01,
       "opposite-field fixture discriminates event field");
  EventGeometry old=cached[0];cached[0].Clear();h.ClearReplayTrackGeometry();
  Need(old.eventIndex==10&&old.points[0].valid&&old.points.size()==3,"pool geometry owns coordinates after current replay bank is cleared");
  // Restore a valid native-equivalent bank for the historical SE QA path.
  h.SetReplayTrackGeometry(tracks,5.);
  FemtoCandidate l=Lambda(10),n=Nucleus(h,10,3);
  const double before=h.NuclearHistograms()->Get("hDphiDeta_proton_4He")->GetEntries();
  Call(h,l,n,false,old,old,true);
  Need(h.NuclearHistograms()->Get("hDphiDeta_proton_4He")->GetEntries()==before+1,
       "historical SE raw map includes the rejected close pair before the gate");
  EqualAliases(h);
}
void MissingCacheRefuses(FemtoLambdaLegacy& h) {
  FemtoCandidate l=Lambda(10),n=Nucleus(h,20,3);EventGeometry a=Geometry(10),b=Geometry(20);
  const double cf=CfEntries(h,"he4",true),q=Qa(h)->GetEntries();
  bool rejected=false,threw=false;
  try {h.FillPair(l,n,.217,.31,4,true,0,&b,&b,&rejected);}catch(const std::runtime_error&){threw=true;}
  Need(threw&&!rejected,"wrong Lambda event cache fails before a cut decision");
  threw=false;try {h.FillPair(l,n,.217,.31,4,true,0,&a,&a,&rejected);}catch(const std::runtime_error&){threw=true;}
  Need(threw&&!rejected,"wrong nucleus event cache fails before a cut decision");
  b.points.erase(2);threw=false;
  try {h.FillPair(l,n,.217,.31,4,true,0,&a,&b,&rejected);}catch(const std::runtime_error&){threw=true;}
  Need(threw&&!rejected,"missing saved track fails loudly instead of becoming unreachable");
  Need(CfEntries(h,"he4",true)==cf&&Qa(h)->GetEntries()==q,"cache errors cause no partial CF/cut-QA fill");
}
void He3Unchanged(FemtoLambdaLegacy& h) {
  for(int direction=0;direction<3;++direction) {
    const bool mixed=direction!=0;const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;
    FemtoCandidate l=Lambda(le),n=Nucleus(h,ne,2);
    const double physical=CfEntries(h,"he3",mixed),mapped=CfEntries(h,"he3_as_he4",mixed),q=Qa(h)->GetEntries();
    bool reason=true;
    Need(h.FillPair(l,n,.217,.31,4,mixed,0,0,0,&reason)&&!reason,"physical He3 is accepted without any pair geometry cache");
    Need(CfEntries(h,"he3",mixed)==physical+1&&CfEntries(h,"he3_as_he4",mixed)==mapped+1,
         "physical He3 and mapped He3AsHe4 keep their pair fanout");
    Need(Qa(h)->GetEntries()==q,"He3 and mapped He3AsHe4 never enter He4/proton cut counters");
  }
}
}

int test_actual_lambda_merging_paircut(const char* mainconf) {
  using namespace ActualLambdaMergingTest;
  try {
    checks=calls=0;FemtoLambdaLegacy h;
    Need(mainconf&&*mainconf&&h.Init(mainconf,std::cerr,false,false),"actual new study mainconf initializes with compact replay controls");
    Need(h.PairGeometryEnabled()&&h.IsPairNucleus(2)&&h.IsPairNucleus(3),"pair geometry and physical He3/He4 control profile enabled");
    Need(h.MassRegion(Lambda(10).reso.invMass)==1,"fixture Lambda mass lies inside configured signal window");
    h.SetPairEventPlaneBin(0);OwnEventFieldAndOwnership(h);
    for(int direction=0;direction<3;++direction) {
      const bool mixed=direction!=0;const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;
      FemtoCandidate l=Lambda(le),n=Nucleus(h,ne,3);
      EventGeometry lg=Geometry(le),ng=mixed?Geometry(ne):lg;
      Call(h,l,n,mixed,lg,ng,true); // SE and both ME directions have nonempty rejection.
      lg=Geometry(le,.20,.20,.01,.01);ng=mixed?Geometry(ne):lg;
      Call(h,l,n,mixed,lg,ng,false); // Pion close alone cannot veto the pair.
      lg=Geometry(le);ng=mixed?Geometry(ne):lg;lg.points[0].valid=false;
      if(!mixed)ng=lg;
      Call(h,l,n,mixed,lg,ng,false,false,true); // Unreachable proton stays.
      lg=Geometry(le);ng=mixed?Geometry(ne):lg;ng.points[2].valid=false;
      if(!mixed)lg=ng;
      Call(h,l,n,mixed,lg,ng,false,false,false); // Unreachable nucleus stays.
    }
    FemtoCandidate l=Lambda(10),n=Nucleus(h,10,3);
    const double edge=.05,inf=std::numeric_limits<double>::infinity();
    for(int sign=-1;sign<=1;sign+=2)for(int axis=0;axis<2;++axis) {
      EventGeometry e=Geometry(10,axis==0?sign*edge:.01,axis==1?sign*edge:.01);
      Call(h,l,n,false,e,e,false); // Exact +/- threshold is kept, strict AND.
      e.points[0]=Point(axis==0?std::nextafter(sign*edge,0.):.01,
                        axis==1?std::nextafter(sign*edge,0.):.01);
      Call(h,l,n,false,e,e,true);
      e.points[0]=Point(axis==0?std::nextafter(sign*edge,sign*inf):.01,
                        axis==1?std::nextafter(sign*edge,sign*inf):.01);
      Call(h,l,n,false,e,e,false);
    }
    // Wrapping across +/-pi must still identify a close mixed-event pair.
    EventGeometry a=Geometry(10,TMath::Pi()-.01,.01),b=Geometry(20);
    b.points[2]=Point(-TMath::Pi()+.01,0.);n=Nucleus(h,20,3);
    Call(h,l,n,true,a,b,true);
    // Identical raw track index from different events is not shared geometry.
    n=Nucleus(h,20,3,0);a=Geometry(10,.20,.01);b=Geometry(20,-.20,0.);
    Call(h,l,n,true,a,b,false);
    // A valid-marked nonfinite point is invalid geometry and follows keep policy.
    a=Geometry(10);b=Geometry(20);n=Nucleus(h,20,3);a.points[0].phi=inf;
    Call(h,l,n,true,a,b,false,false,true);
    MissingCacheRefuses(h);He3Unchanged(h);EqualAliases(h);
    std::cout<<"PASS actual Lambda merging pair cut: checks="<<checks<<" He4 calls="<<calls
             <<"; SE/ME forward/reverse, strict boundaries, wrapped phi, pion-only keep, unreachable keep,"
             <<" own-field/owned-cache, wrong-event/missing-track errors, before/after QA and He3 invariance"<<std::endl;
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"FAIL actual Lambda merging pair cut: "<<error.what()<<std::endl;return 1;
  }
}
#endif
