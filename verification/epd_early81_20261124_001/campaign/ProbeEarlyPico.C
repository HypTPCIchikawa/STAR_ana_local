#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

void ProbeEarlyPico(const char* inputs,const char* output) {
  gEnv->SetValue("XNet.ConnectTimeout",15);
  gEnv->SetValue("XNet.RequestTimeout",25);
  gEnv->SetValue("XNet.MaxRedirectCount",6);
  std::ifstream list(inputs);
  std::ofstream out(output);
  if(!list || !out) { gSystem->Exit(1);return; }
  std::string line;
  int files=0;
  while(std::getline(list,line)) {
    std::istringstream fields(line);
    std::string url, extra;
    Long64_t nominal;
    if(!(fields>>url>>nominal) || (fields>>extra) || nominal!=0) {gSystem->Exit(1);return;}
    std::unique_ptr<TFile> root(TFile::Open(url.c_str(),"READ"));
    Long64_t entries=-1;
    const char* state="unreadable";
    if(root.get() && !root->IsZombie() && !root->TestBit(TFile::kRecovered)) {
      TTree* pico=dynamic_cast<TTree*>(root->Get("PicoDst"));
      if(pico && pico->GetBranch("Event") && pico->GetBranch("Track") && pico->GetBranch("EpdHit")) {
        entries=pico->GetEntries();
        state=entries>0?"has_events":"empty_tree";
      } else state="missing_tree_or_branches";
    }
    out<<url<<'\t'<<entries<<'\t'<<state<<std::endl;
    if(!out.good()) {gSystem->Exit(1);return;}
    ++files;
  }
  std::cout<<"EARLY_PICO_HEADERS_PASS files="<<files<<std::endl;
}
