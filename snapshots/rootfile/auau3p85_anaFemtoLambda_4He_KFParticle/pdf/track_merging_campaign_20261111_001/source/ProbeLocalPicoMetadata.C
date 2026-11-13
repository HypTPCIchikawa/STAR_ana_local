#include "TFile.h"
#include "TTree.h"
#include "TUUID.h"
#include "TROOT.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
static void LocalPicoString(std::ostream& o,const std::string& s){o<<'"';for(size_t i=0;i<s.size();++i){if(s[i]=='"'||s[i]=='\\')o<<'\\';o<<s[i];}o<<'"';}
void ProbeLocalPicoMetadata(const char* list,const char* output){
 try{
  std::ifstream in(list);std::ofstream out(output);if(!in||!out)throw std::runtime_error("cannot open local metadata control/output");
  out<<"{\"rootVersion\":";LocalPicoString(out,gROOT->GetVersion());out<<",\"inputs\":[";std::string path;int n=0;
  while(std::getline(in,path)){
   if(path.empty())continue;TFile f(path.c_str(),"READ");const bool ok=!f.IsZombie()&&!f.TestBit(TFile::kRecovered);TTree* t=ok?dynamic_cast<TTree*>(f.Get("PicoDst")):0;
   if(n)out<<',';out<<"{\"path\":";LocalPicoString(out,path);out<<",\"nonZombie\":"<<(!f.IsZombie()?"true":"false")<<",\"nonRecovered\":"<<(!f.TestBit(TFile::kRecovered)?"true":"false")<<",\"uuid\":";LocalPicoString(out,f.GetUUID().AsString());
   out<<",\"picoTreePresent\":"<<(t?"true":"false")<<",\"entries\":"<<(t?t->GetEntries():-1)<<",\"eventBranchPresent\":"<<(t&&t->GetBranch("Event")?"true":"false")<<",\"trackBranchPresent\":"<<(t&&t->GetBranch("Track")?"true":"false")<<",\"trackCovMatrixBranchPresent\":"<<(t&&t->GetBranch("TrackCovMatrix")?"true":"false")<<'}';
   std::cout<<"LOCAL_PICO_METADATA "<<path<<" entries="<<(t?t->GetEntries():-1)<<" intact="<<ok<<std::endl;++n;
  }
  out<<"],\"count\":"<<n<<",\"eventPayloadRead\":false}\n";out.close();if(!out)throw std::runtime_error("metadata write failure");std::cout<<"LOCAL_PICO_METADATA_PROBE_COMPLETE"<<std::endl;
 }catch(const std::exception& e){std::cerr<<"LOCAL_PICO_METADATA_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
