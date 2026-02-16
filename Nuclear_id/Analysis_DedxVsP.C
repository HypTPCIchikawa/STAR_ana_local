//============================================================
// STAR picoDst analysis macro
//  p vs dE/dx correlation
//============================================================

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

StChain *chain = 0;

void Analysis_DedxVsP(Int_t nEvents = 0,	
		      const char* inList = "test2.list",
		      //const char* inList = "test3.list",
		      TString OutputDir = "./")
{
  if (nEvents == 0) nEvents = 1000000000;

  //============================================================
  // Load STAR libraries
  //============================================================
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  //============================================================
  // Chain & PicoDstMaker
  //============================================================
  chain = new StChain();

  StPicoDstMaker *picoMaker =
    //    new StPicoDstMaker(2, inList, "picoDst");
    new StPicoDstMaker(2, inList, "PicoDst");

  if (chain->Init() == kStErr) {
    cout << "chain->Init() error!" << endl;
    return;
  }

  Int_t totalEntries = picoMaker->chain()->GetEntries();
  cout << "Total entries = " << totalEntries << endl;
  if (nEvents > totalEntries) nEvents = totalEntries;

  //============================================================
  // Histogram
  //============================================================
  TH2D *hDedxP = new TH2D(
			  "hDedxP",
			  "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			  400, 0.0, 4.0,
			  400, 0.0, 20.0);

  //============================================================
  // Event loop
  //============================================================
  for (Int_t i = 0; i < nEvents; i++) {

    if (i % 1000 == 0)
      cout << "Working on event " << i << endl;

    chain->Clear();
    Int_t iret = chain->Make(i);
    if (iret) {
      cout << "Bad return code: " << iret << endl;
      break;
    }

    StPicoDst *picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    StPicoEvent *event = picoDst->event();
    if (!event) continue;

    Int_t nTracks = picoDst->numberOfTracks();
    for (Int_t it = 0; it < nTracks; it++) {

      StPicoTrack *trk = picoDst->track(it);
      if (!trk) continue;

      // ---- basic quality cuts ----
      if (trk->nHitsDedx() < 15) continue;
      if (trk->gPt() < 0.1) continue;

      // momentum
      TVector3 pMom = trk->pMom();
      Double_t p = pMom.Mag();

      // dE/dx
      Double_t dedx = trk->dEdx();
      if (dedx <= 0) continue;

      hDedxP->Fill(p, dedx);
    }
  }

  //============================================================
  // Draw & Save
  //============================================================
  TCanvas *c1 = new TCanvas("c1", "dE/dx vs p", 800, 700);
  c1->SetLogz();
  hDedxP->Draw("colz");

  TString outname = OutputDir + "/DedxVsP.root";
  //TString outname = "DedxVsP.root";
  TFile *fout = new TFile(outname, "RECREATE");
  hDedxP->Write();
  fout->Close();

  //============================================================
  // Finish
  //============================================================
  chain->Finish();

  delete picoMaker;
  delete chain;

  cout << "******************************************" << endl;
  cout << " Analysis finished successfully" << endl;
  cout << " Total processed events = " << nEvents << endl;
  cout << " Output file: " << outname << endl;
  cout << "******************************************" << endl;
}
