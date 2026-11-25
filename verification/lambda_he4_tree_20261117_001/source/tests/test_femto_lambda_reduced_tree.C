#include "include/FemtoLambdaReducedTree.h"
#include "TFileMerger.h"
#include "TH1D.h"
#include <iostream>
#include <stdexcept>

namespace FemtoReducedTreeTest {
using namespace femto_lambda_tree;
void Require(bool value,const std::string& text){if(!value)throw std::runtime_error(text);}
std::vector<std::string> Species(){std::vector<std::string> v;v.push_back("lambda");v.push_back("he3");v.push_back("he4");return v;}
std::string Path(const std::string& directory,const std::string& name){return directory+"/"+name;}
void Missing(const std::string& path){Require(gSystem->AccessPathName(path.c_str()),"refusing existing test output "+path);}
Event MakeEvent(int run,int number,int input,bool validPlane) {
  Event e;e.runId=run;e.eventId=number;e.inputIndex=input;e.cent9=6;e.cent16=12;e.nTracks=10;
  e.pvX=-.4;e.pvY=-2.2;e.pvZ=200.;e.bField=-4.98;e.refMult=123;e.fxtMult=120;e.nBTOFMatch=17;
  e.refMultCorr=119.12345;e.centWeight=1.125;e.epdValid=validPlane;
  if(validPlane){e.epdStatus=0;e.epdPsi1=.17;e.epdBin=femto_epd::XialeiBin(e.epdPsi1);e.epdRunIndex=0;e.epdCalibrationCent=6;e.epdWeight=2.;e.epdHits=3;}
  else e.epdStatus=2;
  return e;
}
FemtoCandidate Nucleus(const char* species,const Event& e,int index,double mass) {
  FemtoCandidate c;c.speciesKey=species;c.eventIndex=e.inputIndex;c.source=kFemtoCandTrack;c.charge=1;
  c.trk.trackIndex=index;c.trk.nHitsFit=29;c.trk.nSigmaHe3=.75;c.trk.nSigmaHe4=-.6;
  c.trk.nSigmaProton=8.25;c.trk.mass2=3.48123;c.trk.tofMatch=kTRUE;
  TLorentzVector p;p.SetVectM(TVector3(.6,.2,-2.),mass);c.SetP4(p);c.y+=.000001f;return c;
}
Track MakeTrack(int index,int charge) {
  Track t;t.index=index;t.id=200+index;t.charge=charge;t.nHitsFit=29;t.nHitsMax=45;t.nHitsDedx=23;
  t.pMom=TVector3(.3,.1,-1.);t.gMom=TVector3(.312345678,.098765432,-1.012345678);
  t.origin=TVector3(.1,-.2,199.3);t.dEdx=12.34567;t.nSigmaProton=2.75;t.nSigmaPion=-1.25;
  t.nSigmaKaon=3.25;t.nSigmaElectron=-4.25;return t;
}
void MakeRich(const Event& e,FemtoCandidateStore& physical,FemtoNuclearPidPairQaStore& qa,Payload& payload) {
  FemtoLambdaCandidate l;l.candidate.speciesKey="lambda";l.candidate.eventIndex=e.inputIndex;l.candidate.source=kFemtoCandResonance;
  TLorentzVector lp;lp.SetVectM(TVector3(.3,.4,-1.),1.115683);l.candidate.SetP4(lp);
  l.candidate.energy+=.00000011920928955078125f;l.candidate.y+=.000001f;
  l.candidate.reso.invMass=1.11723;l.candidate.reso.dau1EventIndex=e.inputIndex;l.candidate.reso.dau1Index=1;
  l.candidate.reso.dau2EventIndex=e.inputIndex;l.candidate.reso.dau2Index=2;
  l.protonIndex=1;l.pionIndex=2;l.protonId=201;l.pionId=202;l.pdg=3122;
  l.x=.125;l.y=-.375;l.z=199.75;l.mass=1.11723;l.massError=.00123;l.chi2Ndf=1.4567;l.topoChi2Ndf=2.1234;
  l.daughterDistance=.34567;l.distanceToPv=.1234;l.decayLength=4.567;l.decayLengthError=.4567;l.decayLengthSignificance=10.;
  l.vertexLineLength=4.578;l.vertexLineLengthError=.4578;l.vertexLineLengthSignificance=10.;l.cosPointing=.9991;
  l.protonPidPull=.4;l.pionPidPull=-.6;l.protonTofPull=1.2;l.pionTofPull=-1.3;l.protonTofM2=.88;l.pionTofM2=.019;
  l.protonDcaToPv=.234567891234;l.pionDcaToPv=.345678912345;l.protonHelixPathLength=12.3456789123;l.pionHelixPathLength=23.4567891234;
  l.imp5PathValid=true;l.protonHasTof=true;l.pionHasTof=false;
  payload.lambdas.push_back(l);physical["lambda"].push_back(l.candidate);
  physical["he3"].push_back(Nucleus("he3",e,3,2.80839));physical["he4"].push_back(Nucleus("he4",e,3,3.72742));
  Track proton=MakeTrack(1,1);proton.hasTof=kTRUE;proton.tofValid=kTRUE;proton.tofIndex=7;proton.tofCellId=333;
  proton.tofMatchFlag=1;proton.beta=.65;proton.tofTime=12.5;proton.tofHitPosition=TVector3(204.,15.,12.);
  proton.roleFlags=kProtonDaughter;payload.tracks.push_back(proton);
  Track pion=MakeTrack(2,-1);pion.roleFlags=kPionDaughter;payload.tracks.push_back(pion);
  Track helium=MakeTrack(3,1);helium.roleFlags=kHe3Track|kHe4Track;payload.tracks.push_back(helium);
  payload.hasAnyLegacyCandidate=kTRUE;payload.baselinePoolEligible=e.epdValid;payload.originalMixBin=e.epdValid?4:-1;
  payload.sourceFileIndex=3;payload.sourceEntry=12345;payload.stats.rawTracks=10;payload.stats.protonHypotheses=3;payload.stats.pionHypotheses=4;
  FemtoNuclearPidPairQaCandidate q;q.candidate=Nucleus("he3",e,4,2.80839);q.selected=false;qa["he3"].push_back(q);
}
void Write(const std::string& path,int run,bool withEmpty) {
  Missing(path);TFile f(path.c_str(),"RECREATE");Require(!f.IsZombie(),"fixture output open");
  Writer w;Require(w.Init(&f,Species(),"reconstruction","original mixing","synthetic producer",4,-1000000),w.LastError());
  Event e=MakeEvent(run,101,1,false);FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;Payload detail;MakeRich(e,p,q,detail);
  Require(w.Fill(e,p,q,detail),w.LastError());
  if(withEmpty) {
    e=MakeEvent(run,102,7,true);p.clear();q.clear();detail=Payload();
    p["deuteron"].push_back(Nucleus("deuteron",e,5,1.87561));
    detail.hasAnyLegacyCandidate=kTRUE;detail.baselinePoolEligible=kTRUE;detail.originalMixBin=4;detail.originalDeuteronCount=1;
    detail.sourceFileIndex=3;detail.sourceEntry=12351;Require(w.Fill(e,p,q,detail),w.LastError());
  }
  Require(w.Write(),w.LastError());TH1D qa("producerQA","fixed producer sample",2,0.,2.);qa.Sumw2();qa.Fill(.5,2.);qa.Write();
  TNamed completed("FemtoLambdaRunStatus","completed");completed.Write();f.Close();
}
void List(const std::string& path,const std::vector<std::string>& files) {
  Missing(path);std::ofstream out(path.c_str());for(size_t i=0;i<files.size();++i)out<<files[i]<<'\n';Require(bool(out),"fixture list write");
}
int Count(Reader& reader,int* segments=0) {
  Event e;FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;Payload detail;int n=0,s=0;
  while(reader.Next(e,p,q,detail)){++n;if(reader.NewSegment())++s;}
  Require(reader.Finished(),reader.LastError());if(segments)*segments=s;return n;
}
void Merge(const std::string& output,const std::vector<std::string>& files) {
  Missing(output);TFileMerger merger;Require(merger.OutputFile(output.c_str(),"RECREATE"),"merger output");
  for(size_t i=0;i<files.size();++i)Require(merger.AddFile(files[i].c_str()),"merger input");Require(merger.Merge(),"merger execution");
}
}

int test_femto_lambda_reduced_tree(const char* directory) {
  try {
    using namespace FemtoReducedTreeTest;
    Require(directory&&*directory,"explicit fixture directory required");
    gSystem->mkdir(directory,kTRUE);const std::string raw=gSystem->UnixPathName(directory);
    const std::string base=!raw.empty()&&raw[0]=='/'?raw:std::string(gSystem->WorkingDirectory())+"/"+raw;
    const std::string one=Path(base,"one.root"),two=Path(base,"two.root"),duplicate=Path(base,"duplicate_physics.root");
    Write(one,22001001,true);Write(two,22001002,false);Write(duplicate,22001001,false);
    Reader r;Require(r.Init(one.c_str(),Species(),"reconstruction","original mixing",true),r.LastError());
    Event e;FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;Payload detail;
    Require(r.Next(e,p,q,detail),r.LastError());Require(r.NewSegment()&&!e.epdValid&&!detail.baselinePoolEligible,"EP-invalid event retained");
    Require(detail.tracks.size()==3&&detail.lambdas.size()==1&&detail.stats.protonHypotheses==3,"rich payload lengths/stats");
    FemtoCandidateStore expected;FemtoNuclearPidPairQaStore eq;Payload ep;MakeRich(e,expected,eq,ep);
    Require(p["lambda"][0].energy==expected["lambda"][0].energy&&p["lambda"][0].y==expected["lambda"][0].y,"canonical Float values retained");
    Require(detail.lambdas[0].protonDcaToPv==ep.lambdas[0].protonDcaToPv&&detail.lambdas[0].mass==ep.lambdas[0].mass,"Lambda fit precision");
    Require(detail.tracks[0].gMom==ep.tracks[0].gMom&&detail.tracks[0].tofCellId==333&&detail.tracks[0].tofHitPosition==TVector3(204.,15.,12.),"raw TVector3/TOF");
    Require(detail.tracks[1].hasTof==false&&detail.tracks[1].beta== -999&&detail.sourceEntry==12345,"missing TOF/source provenance");
    Require(q["he3"].size()==1&&!q["he3"][0].selected&&q["he3"][0].candidate.trk.trackIndex==4,"control-only candidate retained");
    Require(r.Next(e,p,q,detail)&&!r.NewSegment()&&detail.hasAnyLegacyCandidate&&detail.originalDeuteronCount==1&&p.empty(),"d-only admission retained without d payload");
    Require(!r.Next(e,p,q,detail)&&r.Finished(),"exact EOF");
    Reader remix;Require(remix.Init(one.c_str(),Species(),"reconstruction","new mixing",false),remix.LastError());Require(Count(remix)==2,"remix input");
    Reader wrongMix;Require(!wrongMix.Init(one.c_str(),Species(),"reconstruction","new mixing",true),"strict mixing mismatch rejected");
    Reader wrongSelection;Require(!wrongSelection.Init(one.c_str(),Species(),"other reconstruction","original mixing",false),"reconstruction mismatch rejected");
    std::vector<std::string> both;both.push_back(one);both.push_back(two);
    const std::string joined=Path(base,"joined.root");Merge(joined,both);
    Reader merged;Require(merged.Init(joined.c_str(),Species(),"reconstruction","original mixing",true),merged.LastError());int segments=0;
    Require(Count(merged,&segments)==3&&segments==2,"hadd preserves original segment boundaries");
    const std::string repacked=Path(base,"repacked.root");std::vector<std::string> singleton;singleton.push_back(one);Merge(repacked,singleton);
    Reader pack;Require(pack.Init(repacked.c_str(),Species(),"reconstruction","original mixing",true),pack.LastError());Require(Count(pack)==2,"repack accepted independently of file UUID");
    std::vector<std::string> repeats;repeats.push_back(one);repeats.push_back(one);const std::string repeatedList=Path(base,"duplicate_path.txt");List(repeatedList,repeats);
    Reader repeated;Require(!repeated.Init(repeatedList.c_str(),Species(),"reconstruction","original mixing",true),"duplicate path rejected");
    repeats[1]=repacked;const std::string repeatedSegmentList=Path(base,"duplicate_segment.txt");List(repeatedSegmentList,repeats);
    Reader sameSegment;Require(!sameSegment.Init(repeatedSegmentList.c_str(),Species(),"reconstruction","original mixing",true),"repacked duplicate segment rejected");
    repeats[1]=duplicate;const std::string duplicatePhysicsList=Path(base,"duplicate_physics.txt");List(duplicatePhysicsList,repeats);
    Reader samePhysical;Require(samePhysical.Init(duplicatePhysicsList.c_str(),Species(),"reconstruction","original mixing",true),samePhysical.LastError());
    Require(samePhysical.Next(e,p,q,detail)&&samePhysical.Next(e,p,q,detail),"unique rows precede duplicate");
    Require(!samePhysical.Next(e,p,q,detail)&&!samePhysical.Finished()&&!samePhysical.LastError().empty(),"duplicate physical event rejected across independent producers");
    const std::string incomplete=Path(base,"incomplete.root");Missing(incomplete);Require(gSystem->CopyFile(one.c_str(),incomplete.c_str())==0,"copy incomplete fixture");
    {TFile f(incomplete.c_str(),"UPDATE");TNamed failed("FemtoLambdaRunStatus","failed");failed.Write("",TObject::kOverwrite);f.Close();}
    Reader bad;Require(!bad.Init(incomplete.c_str(),Species(),"reconstruction","original mixing",true),"incomplete producer rejected");
    std::cout<<"PASS reduced-tree round-trip, canonical/geometry/TOF, omitted-species admission, remix, ordered hadd/repack, duplicate and completion guards"<<std::endl;
    return 0;
  } catch(const std::exception& error) {std::cerr<<"FAIL reduced-tree test: "<<error.what()<<std::endl;return 1;}
}
