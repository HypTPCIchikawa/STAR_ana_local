#include "TFile.h"
#include "TH1.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TUUID.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

static void QuickString(std::ostream& o,const std::string& s){o<<'"';for(size_t i=0;i<s.size();++i){char c=s[i];if(c=='"'||c=='\\')o<<'\\';if(c=='\n')o<<"\\n";else o<<c;}o<<'"';}
static void QuickSample(std::ostream& o,const char* path){
 TFile f(path,"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered))throw std::runtime_error("cannot read intact candidate ROOT");
 TNamed* status=dynamic_cast<TNamed*>(f.Get("FemtoLambdaRunStatus"));
 if(!status||std::string(status->GetTitle())!="completed")throw std::runtime_error("candidate producer not completed");
 o<<"{\"file\":";QuickString(o,path);o<<",\"uuid\":";QuickString(o,f.GetUUID().AsString());o<<",\"producerStatus\":";QuickString(o,status->GetTitle());
 const char* params[]={"inputEvents","reconstructedEvents","selectedLambdaCandidates","mixingVzClampedEvents"};
 o<<",\"parameters\":{";
 for(int i=0;i<4;++i){TParameter<Long64_t>* p=dynamic_cast<TParameter<Long64_t>*>(f.Get(params[i]));if(!p)throw std::runtime_error("missing producer counter");if(i)o<<',';QuickString(o,params[i]);o<<':'<<p->GetVal();}o<<'}';
 const char* trees[]={"femtoCandidateEvents","eventLedger","selectedLambda"};o<<",\"treeEntries\":{";
 for(int i=0;i<3;++i){TTree* t=dynamic_cast<TTree*>(f.Get(trees[i]));if(!t)throw std::runtime_error("missing candidate/ledger tree");if(i)o<<',';QuickString(o,trees[i]);o<<':'<<t->GetEntries();}o<<'}';
 const char* hist[]={"hKstarSE_lambda_he4_signal","hKstarSE_lambda_he3_signal","hKstarSE_lambda_he3_as_he4_signal","true/hDphiDeta_proton_4He","true/hDphiDeta_proton_4He_EventField","true/hDphiDeta_proton_3He","true/hDphiDeta_proton_3He_EventField"};
 o<<",\"histograms\":{";
 for(int i=0;i<7;++i){TH1* h=dynamic_cast<TH1*>(f.Get(hist[i]));if(!h)throw std::runtime_error("missing signal/merging histogram");if(i)o<<',';QuickString(o,hist[i]);o<<":{\"entries\":"<<h->GetEntries()<<",\"regularIntegral\":"<<h->Integral()<<",\"class\":";QuickString(o,h->ClassName());o<<'}';}
 o<<"}}";
}
void ProbeQuickInputs(const char* first,const char* second,const char* output){
 try{
  std::ofstream o(output);if(!o)throw std::runtime_error("cannot create probe JSON");o<<std::setprecision(17)<<"{\"rootVersion\":";QuickString(o,gROOT->GetVersion());o<<",\"samples\":[";
  QuickSample(o,first);o<<',';QuickSample(o,second);o<<"]}\n";o.close();if(!o)throw std::runtime_error("failed JSON write");std::cout<<"QUICK_INPUT_PROBE_PASS"<<std::endl;
 }catch(const std::exception& e){std::cerr<<"QUICK_INPUT_PROBE_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
