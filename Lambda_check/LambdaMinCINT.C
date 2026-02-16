// LambdaMinCINT.C  (ROOT5 / CINT 専用)

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
using namespace std;

// forward declaration only
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;
StPicoDstMaker *picoMaker = 0;

void LambdaMinCINT(Int_t nEvents=10,
                   const char* inList="test2.list",
                   const char* outFile="lambda_min.root")
{
  // --- load STAR libs (same as PicoQuickCheck) ---
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  // --- create chain ---
  chain = new StChain();

  // --- IMPORTANT: create PicoDstMaker via ProcessLine ---
  TString cmd;
  cmd.Form(
	   "picoMaker = new StPicoDstMaker(2, \"%s\", \"picoDst\");",
    inList
	   );
  gROOT->ProcessLine(cmd);

  if (!picoMaker) {
    cout << "ERROR: picoMaker not created" << endl;
    return;
  }

  // --- init ---
  chain->Init();

  Int_t nEntries = picoMaker->chain()->GetEntries();
  Int_t nRun = TMath::Min(nEntries, nEvents);

  cout << "[LambdaMinCINT] entries = " << nEntries
       << ", processing = " << nRun << endl;

  // --- event loop ---
  for (Int_t i=0; i<nRun; i++) {
    chain->Clear();
    chain->Make(i);

    StPicoDst *dst = picoMaker->picoDst();
    if (!dst) continue;

    StPicoEvent *evt = dst->event();
    if (!evt) continue;

    cout << "Event " << i
         << "  nTracks=" << dst->numberOfTracks()
         << "  PVz=" << evt->primaryVertex().z()
         << endl;

    // example: first track
    if (dst->numberOfTracks() > 0) {
      StPicoTrack *t = dst->track(0);
      if (t) {
        cout << "  charge=" << t->charge()
             << " pPt=" << t->pPt()
             << " nSigmaP=" << t->nSigmaProton()
             << endl;
      }
    }
  }

  chain->Finish();
}
