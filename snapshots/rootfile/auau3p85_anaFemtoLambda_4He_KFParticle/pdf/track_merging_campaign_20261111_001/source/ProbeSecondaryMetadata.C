#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TUUID.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

static void SecondaryString(std::ostream& o,const std::string& s){o<<'"';for(size_t i=0;i<s.size();++i){unsigned char c=s[i];if(c=='"'||c=='\\')o<<'\\'<<char(c);else if(c=='\n')o<<"\\n";else if(c=='\r')o<<"\\r";else if(c=='\t')o<<"\\t";else if(c>=32)o<<char(c);}o<<'"';}
static std::string SecondaryNamed(TFile& f,const char* key){TNamed* x=dynamic_cast<TNamed*>(f.Get(key));return x?std::string(x->GetTitle()):std::string();}
static std::string SecondaryObjectString(TFile& f,const char* key){TObjString* x=dynamic_cast<TObjString*>(f.Get(key));return x?std::string(x->GetString().Data()):std::string();}
static Long64_t SecondaryCounter(TFile& f,const char* key){TParameter<Long64_t>* x=dynamic_cast<TParameter<Long64_t>*>(f.Get(key));return x?x->GetVal():-1;}
static Long64_t SecondaryEntries(TFile& f,const char* key){TTree* x=dynamic_cast<TTree*>(f.Get(key));return x?x->GetEntries():-1;}
void ProbeSecondaryMetadata(const char* reference,const char* list,const char* output){
 try{
  const char* metadataKeys[]={"FemtoCandidateEventStoreLayout","FemtoCandidateEventStoreSelection","FemtoCandidateEventStoreMixing"};
  TFile ref(reference,"READ");if(ref.IsZombie()||ref.TestBit(TFile::kRecovered)||SecondaryNamed(ref,"FemtoLambdaRunStatus")!="completed")throw std::runtime_error("reference is not an intact completed producer");
  std::string expected[3];for(int i=0;i<3;++i){expected[i]=SecondaryObjectString(ref,metadataKeys[i]);if(expected[i].empty())throw std::runtime_error("missing reference cache metadata");}
  const std::string expectedMode=SecondaryNamed(ref,"FemtoCandidateEventStoreMode");
  if(expectedMode.empty()||SecondaryObjectString(ref,"FemtoCandidateEventStoreFileUUID")!=ref.GetUUID().AsString())throw std::runtime_error("reference cache mode/UUID invalid");
  std::ifstream in(list);if(!in)throw std::runtime_error("cannot read nine-file control list");std::ofstream out(output);if(!out)throw std::runtime_error("cannot create metadata result");
  out<<std::setprecision(17)<<"{\"rootVersion\":";SecondaryString(out,gROOT->GetVersion());out<<",\"reference\":{\"file\":";SecondaryString(out,reference);out<<",\"uuid\":";SecondaryString(out,ref.GetUUID().AsString());out<<",\"producerStatus\":\"completed\",\"mode\":";SecondaryString(out,expectedMode);out<<",\"metadata\":{";
  for(int i=0;i<3;++i){if(i)out<<',';SecondaryString(out,metadataKeys[i]);out<<':';SecondaryString(out,expected[i]);}out<<"}},\"samples\":[";
  std::string line;int count=0,valid=0;
  while(std::getline(in,line)){
   if(line.empty())continue;std::istringstream fields(line);int index;std::string path;fields>>index>>path;if(!fields)throw std::runtime_error("bad control-list line");
   TFile f(path.c_str(),"READ");const bool intact=!f.IsZombie()&&!f.TestBit(TFile::kRecovered);const std::string status=intact?SecondaryNamed(f,"FemtoLambdaRunStatus"):"";
   const std::string uuid=f.GetUUID().AsString(),storedUuid=intact?SecondaryObjectString(f,"FemtoCandidateEventStoreFileUUID"):"";
   const std::string mode=intact?SecondaryNamed(f,"FemtoCandidateEventStoreMode"):"";
   const Long64_t input=intact?SecondaryCounter(f,"inputEvents"):-1,reco=intact?SecondaryCounter(f,"reconstructedEvents"):-1,lambdas=intact?SecondaryCounter(f,"selectedLambdaCandidates"):-1;
   const Long64_t cacheRows=intact?SecondaryEntries(f,"femtoCandidateEvents"):-1,eventRows=intact?SecondaryEntries(f,"eventLedger"):-1,lambdaRows=intact?SecondaryEntries(f,"selectedLambda"):-1;
   bool matches[3];for(int i=0;i<3;++i)matches[i]=intact&&SecondaryObjectString(f,metadataKeys[i])==expected[i];
   const bool uuidValid=!storedUuid.empty()&&storedUuid==uuid;
   const bool countsValid=input>0&&reco>=0&&reco<=input&&lambdas>=0&&cacheRows==reco&&eventRows==input&&lambdaRows==lambdas;
   const bool eligible=intact&&status=="completed"&&uuidValid&&countsValid&&mode==expectedMode&&matches[0]&&matches[1]&&matches[2];
   if(count)out<<',';out<<"{\"subjob\":"<<index<<",\"cache_abs\":";SecondaryString(out,path);out<<",\"isZombie\":"<<(f.IsZombie()?"true":"false")<<",\"isRecovered\":"<<(f.TestBit(TFile::kRecovered)?"true":"false")<<",\"producerStatus\":";SecondaryString(out,status);out<<",\"fileUuid\":";SecondaryString(out,uuid);out<<",\"storedCacheUuid\":";SecondaryString(out,storedUuid);out<<",\"uuidMatches\":"<<(uuidValid?"true":"false")<<",\"mode\":";SecondaryString(out,mode);
   out<<",\"parameters\":{\"inputEvents\":"<<input<<",\"reconstructedEvents\":"<<reco<<",\"selectedLambdaCandidates\":"<<lambdas<<"},\"treeEntries\":{\"femtoCandidateEvents\":"<<cacheRows<<",\"eventLedger\":"<<eventRows<<",\"selectedLambda\":"<<lambdaRows<<"},\"countsMatch\":"<<(countsValid?"true":"false")<<",\"metadataMatchesReference\":{";
   for(int i=0;i<3;++i){if(i)out<<',';SecondaryString(out,metadataKeys[i]);out<<':'<<(matches[i]?"true":"false");}out<<"},\"strictMetadataValid\":"<<(eligible?"true":"false")<<'}';
   ++count;if(eligible)++valid;std::cout<<"SECONDARY_METADATA "<<index<<" status="<<status<<" input="<<input<<" cache="<<cacheRows<<" valid="<<eligible<<std::endl;
  }
  out<<"],\"sampleCount\":"<<count<<",\"strictMetadataValidCount\":"<<valid<<",\"candidatePayloadDecoded\":false,\"picoPayloadRead\":false}\n";out.close();if(!out)throw std::runtime_error("metadata result write failure");std::cout<<"SECONDARY_METADATA_PROBE_PASS samples="<<count<<" valid="<<valid<<std::endl;
 }catch(const std::exception& e){std::cerr<<"SECONDARY_METADATA_PROBE_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
