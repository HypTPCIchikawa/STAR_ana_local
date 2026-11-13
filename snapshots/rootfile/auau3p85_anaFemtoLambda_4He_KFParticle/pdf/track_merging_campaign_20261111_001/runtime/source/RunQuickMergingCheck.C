#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
typedef int (*QuickMergingEntry)(const char*);
void RunQuickMergingCheck(const char* library,const char* settings){
 if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
 Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);if(error){gSystem->Exit(2);return;}
 const char* libs[]={"StarRoot","StarClassLibrary","StPicoEvent","StPicoDstMaker","libStarAnaConfig.so"};for(unsigned i=0;i<sizeof(libs)/sizeof(libs[0]);++i)if(gSystem->Load(libs[i])<0){std::cerr<<"Missing "<<libs[i]<<std::endl;gSystem->Exit(2);return;}
 if(gSystem->Load(library)<0){gSystem->Exit(2);return;}
 QuickMergingEntry f=(QuickMergingEntry)gSystem->DynFindSymbol(library,"RunQuickMergingCheck");if(!f){gSystem->Exit(2);return;}int status=f(settings);if(status)gSystem->Exit(status);
}
