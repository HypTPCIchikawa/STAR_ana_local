#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TUUID.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
static void PicoHeaderString(std::ostream& o,const std::string& s){o<<'"';for(size_t i=0;i<s.size();++i){if(s[i]=='"'||s[i]=='\\')o<<'\\';o<<s[i];}o<<'"';}
void ProbePicoHeaders(const char* list,const char* output){
 try{
  std::ifstream in(list);if(!in)throw std::runtime_error("cannot open input list");std::ofstream out(output);if(!out)throw std::runtime_error("cannot create header audit");
  out<<"{\"rootVersion\":";PicoHeaderString(out,gROOT->GetVersion());out<<",\"inputs\":[";std::string line;int count=0;Long64_t total=0;
  while(std::getline(in,line)){
   if(line.empty())continue;std::istringstream fields(line);std::string url;Long64_t catalog;fields>>url>>catalog;if(!fields)throw std::runtime_error("bad catalog list row");
   TFile* f=TFile::Open(url.c_str(),"READ");if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered))throw std::runtime_error("cannot read intact remote Pico");
   TTree* t=dynamic_cast<TTree*>(f->Get("PicoDst"));if(!t)throw std::runtime_error("missing PicoDst TTree");
   if(count)out<<',';out<<"{\"url\":";PicoHeaderString(out,url);out<<",\"catalogEvents\":"<<catalog<<",\"treeEntries\":"<<t->GetEntries()<<",\"uuid\":";PicoHeaderString(out,f->GetUUID().AsString());out<<'}';total+=t->GetEntries();++count;delete f;
  }
  out<<"],\"count\":"<<count<<",\"totalTreeEntries\":"<<total<<"}\n";out.close();if(!out)throw std::runtime_error("header audit write failure");std::cout<<"PICO_HEADER_PROBE_PASS inputs="<<count<<" entries="<<total<<std::endl;
 }catch(const std::exception& e){std::cerr<<"PICO_HEADER_PROBE_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
