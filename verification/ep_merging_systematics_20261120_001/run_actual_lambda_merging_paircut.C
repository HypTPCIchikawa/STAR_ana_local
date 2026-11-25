// Thin ROOT5 loader: compile the production-type test with ACLiC, never CINT.
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"

void run_actual_lambda_merging_paircut(const char* mainconf) {
  if(!mainconf||!*mainconf){gSystem->Exit(2);return;}
  printf("ACTUAL_PAIR_CUT_ROOT %s\n",gROOT->GetVersion());
  gSystem->Exec("g++ --version");
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error){gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker","StEpdUtil"};
  for(int i=0;i<5;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(2);return;}
  const TString cwd=gSystem->WorkingDirectory();
  const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for(int i=0;i<3;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/src/third_party/yaml-cpp/include");
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString study=cwd+"/rootfile/femto_lambda_ep_merging_systematics_20261120_001";
  const TString build=TString::Format("%s/aclic_actual_paircut_%d",study.Data(),gSystem->GetPid());
  if(gSystem->mkdir(build,kTRUE)!=0){gSystem->Exit(2);return;}
  gSystem->SetBuildDir(build,kTRUE);
  if(!gSystem->CompileMacro(study+"/test_actual_lambda_merging_paircut.C","kf")){gSystem->Exit(2);return;}
  TString path(mainconf);path.ReplaceAll("\\","\\\\");path.ReplaceAll("\"","\\\"");
  const Long_t result=gROOT->ProcessLine(TString::Format("test_actual_lambda_merging_paircut(\"%s\")",path.Data()),&error);
  gSystem->Exit(error?2:int(result));
}
