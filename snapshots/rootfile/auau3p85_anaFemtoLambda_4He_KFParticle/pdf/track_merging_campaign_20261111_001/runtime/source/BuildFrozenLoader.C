#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>
void BuildFrozenLoader(const char* runtime){TString build=TString(runtime)+"/lib";gSystem->SetBuildDir(build.Data(),kTRUE);TString source=TString(runtime)+"/source/RunQuickMergingCheck.C";if(!gSystem->CompileMacro(source.Data(),"kf")){gSystem->Exit(2);return;}}
