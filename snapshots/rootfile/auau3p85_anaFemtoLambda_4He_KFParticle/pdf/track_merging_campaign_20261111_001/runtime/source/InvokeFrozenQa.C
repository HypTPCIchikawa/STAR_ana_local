#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>
TString FrozenQaQuote(const char* p){TString s(p);s.ReplaceAll("\\","\\\\");s.ReplaceAll("\"","\\\"");return TString("\"")+s+"\"";}
void InvokeFrozenQa(const char* runtime,const char* settings){
 TString loader=TString(runtime)+"/lib/RunQuickMergingCheck_C.so";
 if(gSystem->Load(loader.Data())<0){std::cerr<<"Cannot load frozen loader"<<std::endl;gSystem->Exit(2);return;}
 TString qa=TString(runtime)+"/lib/libQuickMergingCheck.so";Int_t error=0;
 gROOT->ProcessLine(TString("RunQuickMergingCheck(")+FrozenQaQuote(qa.Data())+","+FrozenQaQuote(settings)+")",&error);
 if(error)gSystem->Exit(2);
}
