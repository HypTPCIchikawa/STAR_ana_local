// LambdaFromPico.C
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLorentzVector.h>

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

void LambdaFromPico(Int_t nEvents=100, const char* inList="test2.list",
                    const char* outFile="LambdaQA.root")
{
  // load STAR libs
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inList, "picoDst");
  StPicoTrack::Print();

}
