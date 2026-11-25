#include "TFile.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSystem.h"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void Check(bool ok, const char* error) {
  if(!ok) throw std::runtime_error(error);
}
template<class T> T* Object(TFile& file,const char* name) {
  T* obj=dynamic_cast<T*>(file.Get(name));Check(obj,name);return obj;
}
Long64_t Count(TFile& file,const char* name) {
  return Object<TParameter<Long64_t> >(file,name)->GetVal();
}
void Inspect(TFile& file,const char* runHash) {
  Check(!file.IsZombie() && !file.TestBit(TFile::kRecovered),"corrupt output");
  Check(std::string(Object<TNamed>(file,"FemtoLambdaRunStatus")->GetTitle())=="completed","job status");
  const std::string ep=Object<TObjString>(file,"FemtoEpdEventPlaneDefinition")->GetString().Data();
  Check(ep.find(runHash)!=std::string::npos,"wrong EP run list");
  Check(Count(file,"epdEventsSeen")==Count(file,"epdEventsValid")+Count(file,"epdEventsRejected"),
        "EP status counters inconsistent");
  Check(Count(file,"epdEventsValid")>0,"no valid event plane");
}
}

void VerifyEarlyPhysicsPilot(const char* producer,const char* replay,const char* runHash,int requested) {
  try {
    TFile first(producer,"READ");Inspect(first,runHash);
    TFile second(replay,"READ");Inspect(second,runHash);
    const Long64_t tree=Object<TTree>(first,"femtoLambdaReducedEvents")->GetEntries();
    Check(Count(first,"inputEvents")==requested && tree>0 && tree<=requested,
          "producer event counts");
    Check(Count(second,"inputEvents")==tree && Count(second,"epdEventsSeen")==tree &&
          Count(first,"epdEventsValid")==Count(second,"epdEventsValid"),"replay event counts");
    Check(Object<TParameter<int> >(second,"FemtoLambdaPairEpdBinCount")->GetVal()==8,"EP8 metadata");
    Check(Object<TParameter<double> >(second,"he4RapiditySplitYLabThreshold")->GetVal()==-.5,
          "rapidity splitter metadata");
    for(const char* stem:{"hKstarSEVsCent_lambda_he4_signal","hKstarMEVsCent_lambda_he4_signal"}) {
      TH2* all=Object<TH2>(second,stem);
      TH2* below=Object<TH2>(second,(std::string(stem)+"_YLabBelow").c_str());
      TH2* above=Object<TH2>(second,(std::string(stem)+"_YLabAbove").c_str());
      for(int x=0;x<=all->GetNbinsX()+1;++x) for(int y=0;y<=all->GetNbinsY()+1;++y)
        Check(below->GetBinContent(x,y)+above->GetBinContent(x,y)<=all->GetBinContent(x,y),
              "rapidity split exceeds inclusive counts");
    }
    std::cout<<"EARLY_EP8_PHYSICS_PILOT_PASS treeEntries="<<tree
             <<" validEP="<<Count(second,"epdEventsValid")<<std::endl;
  } catch(const std::exception& e) {
    std::cerr<<"EARLY_EP8_PHYSICS_PILOT_ERROR "<<e.what()<<std::endl;
    gSystem->Exit(2);
  }
}
