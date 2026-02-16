//============================================================
// p vs dE/dx & p vs m^2 (TOF)
//============================================================

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;

#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

StChain *chain = 0;

void Analysis_DedxVsP_M2(Int_t nEvents = 0,
                         const char* inList = "test.list",
                         TString outDir = "./")
{
  if (nEvents == 0) nEvents = 1e9;

  //=============================
  // Load STAR libraries
  //=============================
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  //=============================
  // Chain
  //=============================
  chain = new StChain();
  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  picoMaker->SetStatus("*",0);
  picoMaker->SetStatus("Event",1);
  picoMaker->SetStatus("Track",1);
  picoMaker->SetStatus("BTofPidTraits",1);

  if (chain->Init() == kStErr) {
    cout << "chain->Init() error!" << endl;
    return;
  }

  int totalEntries = picoMaker->chain()->GetEntries();
  cout << "Total entries = " << totalEntries << endl;
  if (nEvents > totalEntries) nEvents = totalEntries;

  //=============================
  // Histograms
  //=============================
  TH2D *hDedxP = new TH2D(
			  "hDedxP","TPC dE/dx vs p;p (GeV/c);dE/dx",
			  400,0,4,400,0,20);

  TH2D *hM2P = new TH2D(
			"hM2P","TOF m^{2} vs p;p (GeV/c);m^{2} (GeV^{2}/c^{4})",
			//			400,0,4,400,-0.5,3.0);
			400,0,4,1000,-0.5,25.0);

  //=============================
  // Event loop
  //=============================
  for (int i=0;i<nEvents;i++) {

    if (i%1000==0) cout<<"Event "<<i<<endl;

    chain->Clear();
    if (chain->Make(i)) break;

    StPicoDst *picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    for (UInt_t it=0; it<picoDst->numberOfTracks(); it++) {

      StPicoTrack *trk = picoDst->track(it);
      if (!trk) continue;

      // ---- track quality ----
      if (trk->nHitsDedx() < 15) continue;
      if (trk->gPt() < 0.1) continue;

      TVector3 pMom = trk->pMom();
      double p = pMom.Mag();

      // ---- dE/dx ----
      double dedx = trk->dEdx();
      if (dedx>0) hDedxP->Fill(p,dedx);

      // ---- TOF m^2 ----
      if (!trk->isTofTrack()) continue;

      int idx = trk->bTofPidTraitsIndex();
      if (idx < 0) continue;

      StPicoBTofPidTraits *btof =
        picoDst->btofPidTraits(idx);
      if (!btof) continue;

      double beta = btof->btofBeta();
      if (beta<=0) continue;

      double m2 = p*p*(1.0/(beta*beta) - 1.0);
      hM2P->Fill(p,m2);
    }
  }

  //=============================
  // Save
  //=============================
  TFile *fout =
    new TFile(outDir+"/DedxP_M2.root","RECREATE");
  hDedxP->Write();
  hM2P->Write();
  fout->Close();

  chain->Finish();
  delete picoMaker;
  delete chain;

  cout<<"=== Finished successfully ==="<<endl;
}
