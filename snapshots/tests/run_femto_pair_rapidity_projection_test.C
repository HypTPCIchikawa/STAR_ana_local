#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
void run_femto_pair_rapidity_projection_test() {
  const TString cwd=gSystem->WorkingDirectory();
  if(gSystem->Load(cwd+"/lib/libStarAnaConfig.so")<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath((cwd+"/include").Data());
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString build=TString::Format("%s/tmp/pair-rapidity-projection-test-%d",cwd.Data(),gSystem->GetPid());
  gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);
  if(!gSystem->CompileMacro((cwd+"/tests/test_femto_pair_rapidity_projection.C").Data(),"kf")){
    gSystem->Exit(2);return;
  }
  Int_t error=0;
  const Long_t result=gROOT->ProcessLine("test_femto_pair_rapidity_projection()",&error);
  gSystem->Exit(error?2:int(result));
}
