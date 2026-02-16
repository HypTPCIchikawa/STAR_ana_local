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

void Analysis_DedxVsP_job(Int_t nEvents,
                          const char* inList,
                          const char* jobid)
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


  cout << "### MACRO START ###" << endl;
  cout << "STAR = " << gSystem->Getenv("STAR") << endl;
  cout << "LD_LIBRARY_PATH = " << gSystem->Getenv("LD_LIBRARY_PATH") << endl;
  cout << "FILELIST = " << inList << endl;

  cout << "inList = " << inList << endl;

  if (gSystem->AccessPathName(inList)) {
    cout << "ERROR: FILELIST does not exist!" << endl;
    gSystem->Exit(1);
  }

  //============================================================
  // Prepare output directory
  //============================================================
  const char* scratch = gSystem->Getenv("SCRATCH");
  gSystem->mkdir(Form("%s/Nuclear_id", scratch), kTRUE);
  cout << "SCRATCH = " << scratch << endl;
  gSystem->Exec("ls -l $SCRATCH");
  gSystem->Exec("ls -l $SCRATCH/Nuclear_id");


  //============================================================
  // Chain & PicoDstMaker
  //============================================================
  chain = new StChain();

  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  picoMaker->SetStatus("*",0);
  picoMaker->SetStatus("Event",1);
  picoMaker->SetStatus("Track",1);
  picoMaker->SetStatus("BTofHit",1);
  picoMaker->SetStatus("BTofPidTraits",1);
  picoMaker->SetStatus("BbcHit",1);
  picoMaker->SetStatus("EpdHit",1);
  picoMaker->SetStatus("MtdHit",1);
  picoMaker->SetStatus("BTowHit",1);
  picoMaker->SetStatus("ETofPidTraits",1);
  
  cout << "### BEFORE Init()" << endl;
  if (chain->Init() == kStErr) {
    cout << "chain->Init() error!" << endl;
    //    return;
  }
  cout << "### AFTER Init(): "<< endl;
  
  cout << "### INPUT FILE CONTENT ###" << endl;
  picoMaker->chain()->GetFile()->ls();

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
			  800, 0.0, 60.0);
  
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

  TString outname = Form("%s/Nuclear_id/DedxVsP_%s.root",
                         gSystem->Getenv("SCRATCH"), jobid);
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
