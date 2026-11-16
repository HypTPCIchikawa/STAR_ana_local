#include "../include/FemtoCandidateEventStore.h"
#include "TMD5.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SplitFemtoCandidateStore {
using namespace femto_candidate_store;
void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
std::string Quote(const std::string& s){std::string o="\"";for(size_t i=0;i<s.size();++i){if(s[i]=='\\'||s[i]=='\"')o+='\\';if(s[i]=='\n')o+="\\n";else if(s[i]=='\r')o+="\\r";else o+=s[i];}return o+'\"';}
std::string ReadText(TFile& f,const char* name){TObjString* value=dynamic_cast<TObjString*>(f.Get(name));Require(value!=0,std::string("missing source metadata ")+name);return value->GetString().Data();}
void Remap(FemtoCandidate& c,Int_t from,Int_t to){
 Require(c.eventIndex==from,"candidate inputIndex mismatch before fixture remap");c.eventIndex=to;
 if(c.reso.dau1EventIndex>=0){Require(c.reso.dau1EventIndex==from,"first daughter event mismatch");c.reso.dau1EventIndex=to;}
 if(c.reso.dau2EventIndex>=0){Require(c.reso.dau2EventIndex==from,"second daughter event mismatch");c.reso.dau2EventIndex=to;}
}
void RemapAll(FemtoCandidateStore& physical,FemtoNuclearPidPairQaStore& qa,Int_t from,Int_t to){
 for(FemtoCandidateStore::iterator s=physical.begin();s!=physical.end();++s)for(size_t j=0;j<s->second.size();++j)Remap(s->second[j],from,to);
 for(FemtoNuclearPidPairQaStore::iterator s=qa.begin();s!=qa.end();++s)for(size_t j=0;j<s->second.size();++j)Remap(s->second[j].candidate,from,to);
}
void HashFloat(TMD5& digest,const FemtoCandidate& c){Buffers b;Append(b,c,0,kPhysical,false);digest.Update(reinterpret_cast<const UChar_t*>(&(*b.floats)[0]),kFloatFields*sizeof(Float_t));}
bool SameCandidate(const FemtoCandidate& a,const FemtoCandidate& b){
 if(a.speciesKey!=b.speciesKey)return false;Buffers x,y;Append(x,a,0,kPhysical,false);Append(y,b,0,kPhysical,false);
 return *x.ints==*y.ints && *x.quality==*y.quality && std::memcmp(&(*x.floats)[0],&(*y.floats)[0],kFloatFields*sizeof(Float_t))==0;
}
bool SameEvent(const Event& a,const Event& b){
 const Int_t ai[]={a.inputIndex,a.runId,a.eventId,a.cent9,a.cent16,a.refMult,a.fxtMult,a.nTracks,a.nBTOFMatch,a.eventPlane2Tracks};
 const Int_t bi[]={b.inputIndex,b.runId,b.eventId,b.cent9,b.cent16,b.refMult,b.fxtMult,b.nTracks,b.nBTOFMatch,b.eventPlane2Tracks};
 const Double_t ad[]={a.pvX,a.pvY,a.pvZ,a.bField,a.refMultCorr,a.centWeight,a.psi2,a.eventPlane2,a.q2x,a.q2y};
 const Double_t bd[]={b.pvX,b.pvY,b.pvZ,b.bField,b.refMultCorr,b.centWeight,b.psi2,b.eventPlane2,b.q2x,b.q2y};
 return std::memcmp(ai,bi,sizeof(ai))==0 && std::memcmp(ad,bd,sizeof(ad))==0 && a.eventPlaneValid==b.eventPlaneValid && a.eventPlane2Valid==b.eventPlane2Valid;
}
void SameMaps(const FemtoCandidateStore& a,const FemtoCandidateStore& b,const FemtoNuclearPidPairQaStore& qa,const FemtoNuclearPidPairQaStore& qb,TMD5& sourceHash,TMD5& splitHash,Long64_t& physicalCount,Long64_t& qaCount){
 Require(a.size()==b.size()&&qa.size()==qb.size(),"collection map presence changed");
 FemtoCandidateStore::const_iterator y=b.begin();for(FemtoCandidateStore::const_iterator x=a.begin();x!=a.end();++x,++y){
  Require(x->first==y->first&&x->second.size()==y->second.size(),"physical species/list changed");
  for(size_t j=0;j<x->second.size();++j){Require(SameCandidate(x->second[j],y->second[j]),"physical candidate payload changed");HashFloat(sourceHash,x->second[j]);HashFloat(splitHash,y->second[j]);++physicalCount;}
 }
 FemtoNuclearPidPairQaStore::const_iterator qy=qb.begin();for(FemtoNuclearPidPairQaStore::const_iterator qx=qa.begin();qx!=qa.end();++qx,++qy){
  Require(qx->first==qy->first&&qx->second.size()==qy->second.size(),"QA species/list changed");
  for(size_t j=0;j<qx->second.size();++j){Require(qx->second[j].selected==qy->second[j].selected&&SameCandidate(qx->second[j].candidate,qy->second[j].candidate),"QA candidate/selected flag changed");HashFloat(sourceHash,qx->second[j].candidate);HashFloat(splitHash,qy->second[j].candidate);++qaCount;}
 }
}
void Run(const char* input,const char* outputDirectory){
 Require(input&&*input&&outputDirectory&&*outputDirectory,"explicit source cache and new output directory required");std::string source=input,out=outputDirectory;
 Require(gSystem->AccessPathName(out.c_str()),"fixture output directory already exists");
 std::string selection,mixing,layout,sourceUuid;std::vector<std::string> species;
 {TFile f(input,"READ");Require(!f.IsZombie(),"cannot read source cache");selection=ReadText(f,kSelectionName);mixing=ReadText(f,kMixingName);layout=ReadText(f,kLayoutName);sourceUuid=f.GetUUID().AsString();std::istringstream lines(layout);std::string line;Require(bool(std::getline(lines,line))&&line==kSchemaVersion,"unknown source schema");while(std::getline(lines,line))if(line.compare(0,8,"species ")==0)species.push_back(line.substr(8));f.Close();}
 std::string error;Require(ValidSpecies(species,error),error);Require(Layout(species)==layout,"source layout differs from supported exact schema");
 Reader reader;Require(reader.Init(input,species,selection,mixing,true,true),reader.LastError());const Long64_t total=reader.EntriesTotal();Require(total>=2&&total/2<=std::numeric_limits<Int_t>::max()&&total-total/2<=std::numeric_limits<Int_t>::max(),"source must contain at least two events and fit local Int_t indices");const Long64_t boundary=total/2;
 Require(gSystem->mkdir(out.c_str(),true)==0,"cannot create fixture directory");Long64_t written=0;std::string outputUuid[2];
 for(int part=0;part<2;++part){
  const Long64_t count=part==0?boundary:total-boundary;std::ostringstream filename;filename<<out<<"/part_"<<std::setfill('0')<<std::setw(3)<<part<<".root";
  Writer writer;TFile f(filename.str().c_str(),"CREATE");Require(!f.IsZombie(),"cannot create split cache");Require(writer.Init(&f,species,selection,mixing),writer.LastError());
  for(Long64_t row=0;row<count;++row){Event e;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;Require(reader.Next(e,physical,qa),reader.LastError().empty()?"unexpected source EOF":reader.LastError());const Int_t from=e.inputIndex,to=static_cast<Int_t>(row+1);RemapAll(physical,qa,from,to);e.inputIndex=to;Require(writer.Fill(e,physical,qa),writer.LastError());++written;}
  Require(writer.Write(),writer.LastError());f.cd();TNamed fixture("FemtoCandidateSplitFixture","closure-only split: local inputIndex restarts in each stream; all candidate payloads fixed");TObjString inputPath(input),originalUuid(sourceUuid.c_str());TParameter<Long64_t> first("sourceAcceptedEntryFirst",part==0?1:boundary+1),last("sourceAcceptedEntryLast",part==0?boundary:total);
  Require(fixture.Write()>0&&inputPath.Write("FemtoCandidateSplitFixtureSource")>0&&originalUuid.Write("FemtoCandidateSplitFixtureSourceUUID")>0&&first.Write()>0&&last.Write()>0,"fixture provenance write failed");f.Flush();Require(!f.TestBit(TFile::kWriteError),"split cache data flush failed");TNamed completed(kProducerStatusName,"completed");Require(completed.Write()>0,"fixture completed marker write failed");f.Flush();Require(!f.TestBit(TFile::kWriteError),"split cache completed flush failed");outputUuid[part]=f.GetUUID().AsString();f.Close();
 }
 Event ended;FemtoCandidateStore pEnd;FemtoNuclearPidPairQaStore qEnd;Require(written==total&&!reader.Next(ended,pEnd,qEnd)&&reader.Finished(),"source entry total/EOF mismatch");
 const std::string listPath=out+"/stores.list";{std::ofstream list(listPath.c_str());Require(list.good(),"cannot write relative store list");list<<"part_000.root\npart_001.root\n";list.close();Require(!list.fail(),"store list close failed");}
 // Independently read back source and split stores. Only deliberately changed
 // event references are restored before the full bitwise payload comparison.
 Reader original,split;Require(original.Init(input,species,selection,mixing,true,true),original.LastError());Require(split.Init(listPath.c_str(),species,selection,mixing,true,true),split.LastError());Require(original.EntriesTotal()==total&&split.EntriesTotal()==total,"verified store totals differ");TMD5 sourceFloats,splitFloats;Long64_t physicalCount=0,qaCount=0,boundaries=0;
 for(Long64_t row=0;row<total;++row){
  Event a,b;FemtoCandidateStore pa,pb;FemtoNuclearPidPairQaStore qa,qb;Require(original.Next(a,pa,qa),original.LastError());Require(split.Next(b,pb,qb),split.LastError());const Int_t expected=static_cast<Int_t>(row<boundary?row+1:row-boundary+1);Require(b.inputIndex==expected,"split local input index did not restart");Require(split.NewFile()==(row==0||row==boundary),"split file boundary flag changed");if(split.NewFile())++boundaries;RemapAll(pb,qb,b.inputIndex,a.inputIndex);b.inputIndex=a.inputIndex;Require(SameEvent(a,b),"event record changed beyond intended inputIndex");SameMaps(pa,pb,qa,qb,sourceFloats,splitFloats,physicalCount,qaCount);
 }
 Require(!original.Next(ended,pEnd,qEnd)&&original.Finished(),"original verification EOF mismatch");Require(!split.Next(ended,pEnd,qEnd)&&split.Finished(),"split verification EOF mismatch");Require(boundaries==2,"expected two independent stream boundaries");sourceFloats.Final();splitFloats.Final();const std::string sourceDigest=sourceFloats.AsString(),splitDigest=splitFloats.AsString();Require(sourceDigest==splitDigest,"Float_t payload digest changed");Require(outputUuid[0]!=sourceUuid&&outputUuid[1]!=sourceUuid&&outputUuid[0]!=outputUuid[1],"fixture streams need distinct own UUIDs");
 std::ofstream report((out+"/split_verification.json").c_str());Require(report.good(),"cannot create split verification JSON");report<<"{\n\"phase\":\"passed\",\n\"source\":"<<Quote(source)<<",\n\"sourceUUID\":"<<Quote(sourceUuid)<<",\n\"sourceEntries\":"<<total<<",\n\"partEntries\":["<<boundary<<','<<total-boundary<<"],\n\"partUUIDs\":["<<Quote(outputUuid[0])<<','<<Quote(outputUuid[1])<<"],\n\"physicalCandidatesCompared\":"<<physicalCount<<",\n\"qaCandidatesCompared\":"<<qaCount<<",\n\"floatFieldsPerCandidate\":"<<kFloatFields<<",\n\"sourceFloatPayloadMD5\":"<<Quote(sourceDigest)<<",\n\"splitFloatPayloadMD5\":"<<Quote(splitDigest)<<",\n\"allCandidateFieldsAndQAFlagsBitwisePreservedExceptEventRefs\":true,\n\"allEventFieldsBitwisePreservedExceptInputIndex\":true,\n\"localInputIndicesRestartAtOne\":true,\n\"fileBoundariesVerified\":2,\n\"list\":"<<Quote(listPath)<<"\n}\n";report.close();Require(!report.fail(),"split verification JSON close failed");std::cout<<"SPLIT_FEMTO_CANDIDATE_STORE_PASS entries="<<total<<" parts="<<boundary<<','<<total-boundary<<" physical="<<physicalCount<<" qa="<<qaCount<<" floatMD5="<<sourceDigest<<std::endl;
}
}
void splitFemtoCandidateEventStore(const char* inputCache,const char* outputDirectory){try{SplitFemtoCandidateStore::Run(inputCache,outputDirectory);}catch(const std::exception& e){std::cerr<<"SPLIT_FEMTO_CANDIDATE_STORE_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}}
