#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TString.h>
#include <TClass.h>
#include <TTree.h>
#include <iostream>

// The same loader/ACLiC order as tests/run_femto_lambda_reduced_tree_test.C.
// All generated libraries stay beside the study, outside shared analysis sources.
void run_grouped_tree_validation(const char* action, const char* inputs,
    const char* report, const char* merged, Long64_t maxTreeBytes,
    Int_t expectedQa, Int_t sampleEvents, Int_t maxOpenedFiles, Int_t floatUlps) {
  if (!action || !inputs || !report || !merged || maxTreeBytes <= 0 ||
      expectedQa <= 0 || sampleEvents < 1 || sampleEvents > 16 ||
      maxOpenedFiles < 3 || floatUlps < 0) { gSystem->Exit(2); return; }
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    gSystem->Exit(2); return;
  }
  Int_t error=0;
  gROOT->ProcessLine("loadSharedLibraries();", &error);
  if (error) { gSystem->Exit(2); return; }
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker","StEpdUtil"};
  for (int starIndex=0;starIndex<5;++starIndex) if (gSystem->Load(star[starIndex])<0) { gSystem->Exit(2); return; }
  const TString cwd=gSystem->WorkingDirectory();
  const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for (int libraryIndex=0;libraryIndex<3;++libraryIndex) if (gSystem->Load(cwd+"/lib/"+libs[libraryIndex])<0) { gSystem->Exit(2); return; }
  TClass* vectors=gROOT->GetClass("vector<TVector3>");
  if (!vectors || !vectors->GetCollectionProxy()) {
    std::cerr << "Missing project vector<TVector3> collection dictionary" << std::endl;
    gSystem->Exit(2); return;
  }
  gInterpreter->AddIncludePath(cwd);
  gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/src/third_party/yaml-cpp/include");
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString directory=cwd+"/rootfile/femto_lambda_tree_grouped_merge_20261118_001/validation";
  const TString build=TString::Format("%s/build_%d",directory.Data(),gSystem->GetPid());
  if (gSystem->mkdir(build,kTRUE)<0 && gSystem->AccessPathName(build)) { gSystem->Exit(2); return; }
  gSystem->SetBuildDir(build,kTRUE);
  if (!gSystem->CompileMacro(directory+"/tree_group_validation.C","kf")) { gSystem->Exit(2); return; }
  TString a(action),b(inputs),c(report),d(merged);
  TString* strings[]={&a,&b,&c,&d};
  for (int stringIndex=0;stringIndex<4;++stringIndex) { strings[stringIndex]->ReplaceAll("\\","\\\\");strings[stringIndex]->ReplaceAll("\"","\\\""); }
  const TString call=TString::Format("tree_group_validation(\"%s\",\"%s\",\"%s\",\"%s\",%lldLL,%d,%d,%d,%d)",a.Data(),b.Data(),c.Data(),d.Data(),maxTreeBytes,expectedQa,sampleEvents,maxOpenedFiles,floatUlps);
  const Long_t result=gROOT->ProcessLine(call,&error);
  gSystem->Exit(error?2:int(result));
}
