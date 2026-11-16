#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
void run_femto_epd_kstar_splits_test(const char* ep8Mainconf,const char* ep1Mainconf) {
  if(!ep8Mainconf||!*ep8Mainconf||!ep1Mainconf||!*ep1Mainconf){gSystem->Exit(2);return;}
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error){gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker","StEpdUtil"};
  for(unsigned i=0;i<5;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(2);return;}
  const TString cwd=gSystem->WorkingDirectory();
  const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so","libStNuclearIdMaker.so","libStFemtoMaker.so"};
  for(unsigned i=0;i<5;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStFemtoMaker -lStNuclearIdMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString build=TString::Format("%s/tmp/epd-kstar-splits-test-%d",cwd.Data(),gSystem->GetPid());
  gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);
  if(!gSystem->CompileMacro("tests/test_femto_epd_kstar_splits.C","kf")){gSystem->Exit(2);return;}
  TString ep8(ep8Mainconf),ep1(ep1Mainconf);
  ep8.ReplaceAll("\\","\\\\");ep8.ReplaceAll("\"","\\\"");ep1.ReplaceAll("\\","\\\\");ep1.ReplaceAll("\"","\\\"");
  const Long_t result=gROOT->ProcessLine(TString::Format("test_femto_epd_kstar_splits(\"%s\",\"%s\")",ep8.Data(),ep1.Data()),&error);
  gSystem->Exit(error?2:int(result));
}
