#include "../include/FemtoCandidateEventStore.h"
#include "TSystem.h"
#include "TROOT.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace FemtoCandidateEventStoreTest {
using namespace femto_candidate_store;
int checks=0;
void Require(bool condition,const std::string& text){++checks;std::cout<<"STORE_TEST_CHECK "<<checks<<" "<<text<<std::endl;if(!condition)throw std::runtime_error(text);}
std::vector<std::string> Species(){std::vector<std::string> s;s.push_back("lambda");s.push_back("he4");s.push_back("he3");return s;}
Event MakeEvent(Int_t index){Event e;e.inputIndex=index;e.runId=22123001;e.eventId=10000+index;e.cent9=7;e.cent16=14;e.refMult=43;e.fxtMult=57;e.nTracks=100;e.nBTOFMatch=39;e.pvX=-.412345678901;e.pvY=-2.312345678901;e.pvZ=200.123456789;e.bField=5.000123456;e.refMultCorr=44.987654321;e.centWeight=1.037123456;e.psi2=-1;e.eventPlaneValid=kFALSE;e.eventPlane2=.723456789;e.q2x=2.123456789;e.q2y=3.987654321;e.eventPlane2Tracks=51;e.eventPlane2Valid=kTRUE;return e;}
FemtoCandidate Track(Int_t event,const std::string& key,Int_t index){
 FemtoCandidate c;c.eventIndex=event;c.speciesKey=key;c.source=kFemtoCandTrack;c.charge=2;
 c.px=1.2345679f;c.py=-.12345679f;c.pz=3.8765432f;c.energy=5.0123458f;c.y=.31234568f;c.pt=987.125f;c.eta=-42.75f;c.phi=-.81234568f;c.qualityFlags=std::numeric_limits<UInt_t>::max();
 c.trk.nSigmaProton=1.0123457f;c.trk.nSigmaKaon=-2.1234567f;c.trk.nSigmaHe4=.9876543f;c.trk.nSigmaDeuteron=4.2345676f;c.trk.nSigmaTriton=-5.3456788f;c.trk.nSigmaHe3=6.456789f;c.trk.mass2=-999.f;c.trk.deltaOneOverBeta=999.f;c.trk.tofMatch=kTRUE;c.trk.dca=.021234568f;c.trk.nHitsFit=std::numeric_limits<Short_t>::max();c.trk.trackIndex=index;c.reso.invMass=.012345679f;c.reso.dcaDaughters=.02345679f;c.reso.betaGamma=-1.1234568f;return c;
}
FemtoCandidate Lambda(Int_t event,Int_t first){FemtoCandidate c=Track(event,"lambda",-1);c.source=kFemtoCandResonance;c.charge=0;c.reso.invMass=1.1156831f;c.reso.dau1EventIndex=event;c.reso.dau1Index=first;c.reso.dau2EventIndex=event;c.reso.dau2Index=first+1;return c;}
bool EqualCandidate(const FemtoCandidate& a,const FemtoCandidate& b){
 Buffers aa,bb;Append(aa,a,0,kPhysical,false);Append(bb,b,0,kPhysical,false);
 return a.speciesKey==b.speciesKey && aa.ints->size()==bb.ints->size() && *aa.ints==*bb.ints && *aa.quality==*bb.quality && aa.floats->size()==bb.floats->size() && std::memcmp(&(*aa.floats)[0],&(*bb.floats)[0],aa.floats->size()*sizeof(Float_t))==0;
}
bool EqualEvent(const Event& a,const Event& b){
 const Int_t ai[]={a.inputIndex,a.runId,a.eventId,a.cent9,a.cent16,a.refMult,a.fxtMult,a.nTracks,a.nBTOFMatch,a.eventPlane2Tracks,a.epdStatus,a.epdBin,a.epdRunIndex,a.epdCalibrationCent,a.epdHits};
 const Int_t bi[]={b.inputIndex,b.runId,b.eventId,b.cent9,b.cent16,b.refMult,b.fxtMult,b.nTracks,b.nBTOFMatch,b.eventPlane2Tracks,b.epdStatus,b.epdBin,b.epdRunIndex,b.epdCalibrationCent,b.epdHits};
 const Double_t ad[]={a.pvX,a.pvY,a.pvZ,a.bField,a.refMultCorr,a.centWeight,a.psi2,a.eventPlane2,a.q2x,a.q2y,a.epdPsi1,a.epdRawQx,a.epdRawQy,a.epdRecenteredQx,a.epdRecenteredQy,a.epdWeight};
 const Double_t bd[]={b.pvX,b.pvY,b.pvZ,b.bField,b.refMultCorr,b.centWeight,b.psi2,b.eventPlane2,b.q2x,b.q2y,b.epdPsi1,b.epdRawQx,b.epdRawQy,b.epdRecenteredQx,b.epdRecenteredQy,b.epdWeight};
 return std::memcmp(ai,bi,sizeof(ai))==0 && std::memcmp(ad,bd,sizeof(ad))==0 && a.eventPlaneValid==b.eventPlaneValid && a.eventPlane2Valid==b.eventPlane2Valid && a.epdValid==b.epdValid;
}
void Metadata(TFile& f,const std::string& selection="selection-fixed",const std::string& mixing="mix-original",bool epd=false){f.cd();TObjString layout(Layout(Species(),epd).c_str()),sel(selection.c_str()),mix(mixing.c_str()),uuid(f.GetUUID().AsString());Require(layout.Write(kLayoutName)>0&&sel.Write(kSelectionName)>0&&mix.Write(kMixingName)>0&&uuid.Write(kFileUuidName)>0,"metadata write");}
Event MakeEpdEvent(Int_t index) {
 Event e=MakeEvent(index);e.epdPsi1=.312345678901;e.epdRawQx=2.123456789012;e.epdRawQy=-3.987654321098;
 e.epdRecenteredQx=1.987654321012;e.epdRecenteredQy=-2.123456789012;e.epdWeight=24.987654321012;
 e.epdStatus=0;e.epdBin=4;e.epdRunIndex=21;e.epdCalibrationCent=7;e.epdHits=39;e.epdValid=kTRUE;
 if(index==11){e.epdValid=kFALSE;e.epdStatus=2;e.epdBin=-1;e.epdRunIndex=-1;e.epdCalibrationCent=-1;e.epdHits=0;e.epdPsi1=-999.;e.epdWeight=0.;}
 if(index==15){e.epdPsi1=0.;e.epdBin=3;} // Historical rounded Float_t bin edge.
 return e;
}
void EpdSchemas(const std::string& dir) {
 {TFile f((dir+"/first.root").c_str(),"READ");TTree* t=dynamic_cast<TTree*>(f.Get(kTreeName));TObjString* layout=dynamic_cast<TObjString*>(f.Get(kLayoutName));
  Require(t&&!t->GetBranch("epdPsi1")&&!t->GetBranch("epdValid")&&layout&&layout->GetString()==Layout(Species()).c_str(),"default writer keeps original schema and has no EPD branches");}
 {TFile f((dir+"/epd.root").c_str(),"CREATE");Writer w;Require(w.Init(&f,Species(),"selection-fixed","mix-original",true),w.LastError());
  const int indices[]={7,11,15};for(unsigned i=0;i<3;++i){FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;p["he4"].push_back(Track(indices[i],"he4",4));Require(w.Fill(MakeEpdEvent(indices[i]),p,q),w.LastError());}Require(w.Write(),w.LastError());f.Close();}
 Reader r;Require(r.Init((dir+"/epd.root").c_str(),Species(),"selection-fixed","mix-original",true,false,true),r.LastError());
 Event e;FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;const int indices[]={7,11,15};
 for(unsigned i=0;i<3;++i){Require(r.Next(e,p,q)&&r.HasEpdPayload(),"read schema2 event including invalid EP");Require(EqualEvent(e,MakeEpdEvent(indices[i])),"all EPD and separate TPC-Q2 fields roundtrip bitwise");Require(EqualCandidate(p["he4"][0],Track(indices[i],"he4",4)),"EPD schema preserves physical candidate payload");}
 Require(!r.Next(e,p,q)&&r.Finished(),"schema2 clean EOF");
 Reader old;Require(!old.Init((dir+"/first.root").c_str(),Species(),"selection-fixed","mix-original",true,false,true)&&old.LastError().find("EPD")!=std::string::npos,"EPD replay rejects a legacy payload");
 {std::ofstream list((dir+"/epd_then_v1.list").c_str());list<<"epd.root\nsecond.root\n";}
 Reader mixed;Require(mixed.Init((dir+"/epd_then_v1.list").c_str(),Species(),"selection-fixed","mix-original",true),mixed.LastError());
 for(unsigned i=0;i<3;++i)Require(mixed.Next(e,p,q),"schema2 precedes legacy stream");
 Require(mixed.Next(e,p,q)&&mixed.NewFile()&&!mixed.HasEpdPayload()&&EqualEvent(e,MakeEvent(7)),"legacy boundary resets every EPD field to default");
 Reader strictMixed;Require(!strictMixed.Init((dir+"/epd_then_v1.list").c_str(),Species(),"selection-fixed","mix-original",true,false,true),"EPD requirement prevalidates all streams");
 {TFile f((dir+"/epd_missing_branch.root").c_str(),"CREATE");Event event=MakeEvent(7);Buffers b;TTree t(kTreeName,"schema2 metadata on a schema1 tree");std::string error;
  Require(Bind(&t,event,b,true,error),error);Require(t.Fill()>0&&t.Write()>0,"write missing-EPD-branch fixture");Metadata(f,"selection-fixed","mix-original",true);t.SetDirectory(0);f.Close();}
 Reader missing;Require(!missing.Init((dir+"/epd_missing_branch.root").c_str(),Species(),"selection-fixed","mix-original",true,false,true),"schema2 metadata cannot mask missing EPD branches");
 for(int bad=0;bad<4;++bad){TFile f((dir+Form("/epd_bad%d.root",bad)).c_str(),"CREATE");Writer w;Require(w.Init(&f,Species(),"selection-fixed","mix-original",true),w.LastError());Event event=MakeEpdEvent(7);
  if(bad==0)event.epdBin=8;else if(bad==1)event.epdBin=3;else if(bad==2)event.epdStatus=2;else event.epdPsi1=std::numeric_limits<double>::quiet_NaN();
  p.clear();q.clear();Require(!w.Fill(event,p,q)&&!w.Write(),"schema2 rejects invalid bin, angle/bin mismatch, status or nonfinite angle");f.Close();}
}
void WriteFile(const std::string& path,int type){
 TFile f(path.c_str(),"CREATE");Require(!f.IsZombie(),"create test store");Writer w;Require(w.Init(&f,Species(),"selection-fixed","mix-original"),w.LastError());
 FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;
 if(type==1){
  Event e=MakeEvent(7);p["lambda"].push_back(Lambda(7,1));p["lambda"].push_back(Lambda(7,2));p["he4"].push_back(Track(7,"he4",4));p["he4"].push_back(Track(7,"he4",5));p["he3"];
  FemtoNuclearPidPairQaCandidate qa;qa.candidate=Track(7,"he4",6);qa.selected=false;q["he4"].push_back(qa);qa.selected=true;q["he4"].push_back(qa);q["he3"];
  Require(w.Fill(e,p,q),w.LastError());
  e=MakeEvent(11);p.clear();q.clear();Require(w.Fill(e,p,q),w.LastError());
  e=MakeEvent(15);e.eventPlane2=-1;e.eventPlane2Valid=kFALSE;e.eventPlane2Tracks=0;p["he3"].push_back(Track(15,"he3",0));qa.candidate=Track(15,"he3",0);qa.selected=false;q["he3"].push_back(qa);Require(w.Fill(e,p,q),w.LastError());
 }else if(type==2){Event e=MakeEvent(7);p["lambda"].push_back(Lambda(7,1));Require(w.Fill(e,p,q),w.LastError());}
 Require(w.Write(),w.LastError());Require(w.Entries()==(type==1?3:(type==2?1:0)),"writer entry count");f.Close();
}
void RawMalformed(const std::string& path,int bad){
 TFile f(path.c_str(),"CREATE");Require(!f.IsZombie(),"create malformed fixture");Event e=MakeEvent(7);Buffers b;TTree tree(kTreeName,"malformed fixture");tree.SetDirectory(&f);std::string error;
 if(bad==3 || bad==4){Float_t wrong=3;if(bad==3)tree.Branch("inputIndex",&e.inputIndex,"inputIndex/I");else tree.Branch("inputIndex",&wrong,"inputIndex/F");Require(tree.Fill()>0,"header fixture fill");Require(tree.Write()>0,"header fixture write");Metadata(f);tree.SetDirectory(0);f.Close();return;}
 Require(Bind(&tree,e,b,true,error),error);for(size_t i=0;i<Species().size();++i){b.physicalPresent->push_back(1);b.qaPresent->push_back(0);}Append(b,Track(7,"he4",4),1,kPhysical,false);
 if(bad==0)b.floats->pop_back();
 else if(bad==1)(*b.collection)[0]=99;
 else if(bad==2)(*b.ints)[5]=e.nTracks;

 Require(tree.Fill()>0,"malformed fixture fill");Require(tree.Write()>0,"malformed fixture write");Metadata(f,bad==5?"selection-other":"selection-fixed");tree.SetDirectory(0);f.Close();
}
void Run(const char* directory,const char* result){
 Require(directory&&*directory&&result&&*result,"explicit fixture/output paths required");std::string dir=directory;
 Require(gSystem->AccessPathName(dir.c_str()),"fixture directory must not already exist");Require(gSystem->mkdir(dir.c_str(),true)==0,"create fixture directory");
 WriteFile(dir+"/empty.root",0);WriteFile(dir+"/first.root",1);WriteFile(dir+"/second.root",2);
 std::ofstream list((dir+"/stores.list").c_str());list<<"# relative paths resolved against this list\nempty.root\nfirst.root\nsecond.root\n";list.close();
 Reader r;Require(r.Init((dir+"/stores.list").c_str(),Species(),"selection-fixed","mix-original",true),r.LastError());Require(r.EntriesTotal()==4,"all-file entry count");Event e;FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;
 Require(r.Next(e,p,q)&&r.NewFile(),"first actual file event");Require(EqualEvent(e,MakeEvent(7)),"exact event including auxiliary plane");Require(p.size()==3&&p["lambda"].size()==2&&p["he4"].size()==2&&p["he3"].empty(),"physical maps/list order");Require(EqualCandidate(p["lambda"][0],Lambda(7,1))&&EqualCandidate(p["lambda"][1],Lambda(7,2)),"all Lambda fields/order");Require(EqualCandidate(p["he4"][0],Track(7,"he4",4))&&EqualCandidate(p["he4"][1],Track(7,"he4",5)),"all nuclear fields/order/UInt max");Require(p["he4"][0].pt==987.125f&&p["he4"][0].energy==5.0123458f,"do not recalculate stored energy/derived values");Require(q.size()==2&&q["he3"].empty()&&q["he4"].size()==2&&!q["he4"][0].selected&&q["he4"][1].selected,"QA flags and empty presence");Require(EqualCandidate(q["he4"][0].candidate,Track(7,"he4",6)),"all QA candidate fields");
 Require(r.Next(e,p,q)&&!r.NewFile()&&e.inputIndex==11&&p.empty()&&q.empty(),"fully empty accepted event preserved");Require(r.Next(e,p,q)&&!r.NewFile()&&e.inputIndex==15&&p.size()==1&&p["he3"].size()==1&&!e.eventPlane2Valid,"one-species event/invalid auxiliary plane");Require(r.Next(e,p,q)&&r.NewFile()&&e.inputIndex==7,"new stream boundary/original index preserved");Require(!r.Next(e,p,q)&&r.Finished()&&r.LastError().empty()&&p.empty()&&q.empty(),"clean EOF distinguished");
 Reader wrongSelection;Require(!wrongSelection.Init((dir+"/first.root").c_str(),Species(),"different","mix-original",true)&&!wrongSelection.LastError().empty(),"reject selection provenance");Reader wrongMix;Require(!wrongMix.Init((dir+"/first.root").c_str(),Species(),"selection-fixed","different",true),"reject original mixing mismatch");Reader changedMix;Require(changedMix.Init((dir+"/first.root").c_str(),Species(),"selection-fixed","different",false)&&changedMix.Next(e,p,q),"allow explicit alternate mixing only");std::vector<std::string> reordered=Species();std::swap(reordered[0],reordered[1]);Reader wrongLayout;Require(!wrongLayout.Init((dir+"/first.root").c_str(),reordered,"selection-fixed","mix-original",true),"reject schema/species order mismatch");
 for(int bad=0;bad<6;++bad){const std::string path=dir+"/malformed"+Form("%d",bad)+".root";RawMalformed(path,bad);Reader rr;bool init=rr.Init(path.c_str(),Species(),"selection-fixed","mix-original",true);if(bad==3||bad==4||bad==5)Require(!init&&!rr.LastError().empty(),"reject header/provenance malformed store at Init");else{Require(init,rr.LastError());Require(!rr.Next(e,p,q)&&!rr.Finished()&&!rr.LastError().empty()&&p.empty()&&q.empty(),"reject malformed payload, not EOF");}}
 std::ofstream late((dir+"/latebad.list").c_str());late<<"first.root\nmalformed5.root\n";late.close();Reader allFiles;Require(!allFiles.Init((dir+"/latebad.list").c_str(),Species(),"selection-fixed","mix-original",true),"prevalidate all input files before replay");
 std::ofstream dup((dir+"/duplicate.list").c_str());dup<<"first.root\nfirst.root\n";dup.close();Reader duplicates;Require(!duplicates.Init((dir+"/duplicate.list").c_str(),Species(),"selection-fixed","mix-original",true),"reject duplicate files");
 {TFile source((dir+"/first.root").c_str(),"READ");TFile repacked((dir+"/repacked.root").c_str(),"CREATE");Require(!source.IsZombie()&&!repacked.IsZombie(),"create repacked fixture");TTree* original=dynamic_cast<TTree*>(source.Get(kTreeName));Require(original!=0,"source tree for repack");repacked.cd();TTree* clone=original->CloneTree(-1,"fast");Require(clone&&clone->Write()>0,"repack candidate tree");const char* keys[]={kLayoutName,kSelectionName,kMixingName,kFileUuidName};for(int k=0;k<4;++k){TObject* item=source.Get(keys[k]);Require(item&&item->Write(keys[k])>0,"copy repacked metadata");}repacked.Close();source.Close();}Reader repacked;Require(!repacked.Init((dir+"/repacked.root").c_str(),Species(),"selection-fixed","mix-original",true)&&repacked.LastError().find("merged/repacked")!=std::string::npos,"reject repacked/hadded original-job stores");
 Require(gSystem->CopyFile((dir+"/first.root").c_str(),(dir+"/first_copy.root").c_str(),false)==0,"bytecopy fixture");Reader copied;Require(copied.Init((dir+"/first_copy.root").c_str(),Species(),"selection-fixed","mix-original",true)&&copied.Next(e,p,q),"accept unchanged bytecopy");std::ofstream alias((dir+"/alias.list").c_str());alias<<"first.root\nfirst_copy.root\n";alias.close();Reader aliasReader;Require(!aliasReader.Init((dir+"/alias.list").c_str(),Species(),"selection-fixed","mix-original",true)&&aliasReader.LastError().find("UUID")!=std::string::npos,"reject duplicate stream under alias paths");
 Reader missingProducer;Require(!missingProducer.Init((dir+"/first.root").c_str(),Species(),"selection-fixed","mix-original",true,true)&&missingProducer.LastError().find("producer")!=std::string::npos,"require completed producer marker when requested");
 Require(gSystem->CopyFile((dir+"/first.root").c_str(),(dir+"/producer_completed.root").c_str(),false)==0,"completed producer fixture copy");{TFile f((dir+"/producer_completed.root").c_str(),"UPDATE");TNamed status(kProducerStatusName,"completed");Require(status.Write(kProducerStatusName,TObject::kOverwrite)>0,"write completed producer marker");f.Close();}Reader completedProducer;Require(completedProducer.Init((dir+"/producer_completed.root").c_str(),Species(),"selection-fixed","mix-original",true,true)&&completedProducer.Next(e,p,q),"accept explicitly completed producer");
 Require(gSystem->CopyFile((dir+"/first.root").c_str(),(dir+"/producer_incomplete.root").c_str(),false)==0,"incomplete producer fixture copy");{TFile f((dir+"/producer_incomplete.root").c_str(),"UPDATE");TNamed status(kProducerStatusName,"incomplete");Require(status.Write(kProducerStatusName,TObject::kOverwrite)>0,"write incomplete producer marker");f.Close();}Reader incompleteProducer;Require(!incompleteProducer.Init((dir+"/producer_incomplete.root").c_str(),Species(),"selection-fixed","mix-original",true,true)&&incompleteProducer.LastError().find("producer")!=std::string::npos,"reject incomplete producer with strict flag");Reader incompleteOptional;Require(!incompleteOptional.Init((dir+"/producer_incomplete.root").c_str(),Species(),"selection-fixed","mix-original",true,false),"reject existing incomplete marker even when optional");
 Require(gSystem->CopyFile((dir+"/first.root").c_str(),(dir+"/producer_wrongtype.root").c_str(),false)==0,"wrong producer marker fixture copy");{TFile f((dir+"/producer_wrongtype.root").c_str(),"UPDATE");TObjString status("completed");Require(status.Write(kProducerStatusName,TObject::kOverwrite)>0,"write malformed producer marker");f.Close();}Reader wrongProducer;Require(!wrongProducer.Init((dir+"/producer_wrongtype.root").c_str(),Species(),"selection-fixed","mix-original",true,false),"reject malformed existing producer marker");
 Reader empty;Require(empty.Init((dir+"/empty.root").c_str(),Species(),"selection-fixed","mix-original",true)&&empty.EntriesTotal()==0&&!empty.Next(e,p,q)&&empty.Finished(),"empty file EOF");
 {TFile f((dir+"/writerbad.root").c_str(),"CREATE");Writer w;Require(w.Init(&f,Species(),"selection-fixed","mix-original"),w.LastError());Event bad=MakeEvent(7);p.clear();q.clear();p["he4"].push_back(Track(8,"he4",4));Require(!w.Fill(bad,p,q)&&!w.LastError().empty()&&!w.Write(),"writer fail-closed on wrong event index");f.Close();}
 EpdSchemas(dir);
 std::ofstream report(result);Require(report.good(),"create test result");report<<"{\"phase\":\"passed\",\"checks\":"<<checks<<",\"rootVersion\":\""<<gROOT->GetVersion()<<"\",\"schema\":\""<<kSchemaVersion<<"\",\"exactFloatAndEnergyRoundtrip\":true,\"physicalAndQaOrder\":true,\"emptyAndOneSpeciesEvents\":true,\"relativeMultiFileInput\":true,\"newFileBoundaries\":true,\"allFilePrevalidation\":true,\"provenanceAndMalformedPayloadRejected\":true,\"fileUuidAndOriginalStreamBoundariesChecked\":true,\"producerCompletionGuardChecked\":true,\"optionalEpdSchema2AndLegacyBoundaryChecked\":true}\n";report.close();std::cout<<"FEMTO_CANDIDATE_EVENT_STORE_TEST_PASS checks="<<checks<<std::endl;
}
}
void testFemtoCandidateEventStore(const char* outputDirectory,const char* resultJson){try{FemtoCandidateEventStoreTest::Run(outputDirectory,resultJson);}catch(const std::exception& e){std::cerr<<"FEMTO_CANDIDATE_EVENT_STORE_TEST_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}}
