// LambdaInvMassMin.C
// ROOT5 / CINT macro, no STAR headers

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

// forward declarations only
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

void LambdaInvMassMin(Int_t nEvents=100,
                      const char* inList="test2.list",
                      const char* outFile="lambda_min.root")
{
  // load STAR libs (same as PicoQuickCheck)
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  chain = new StChain();
  StPicoDstMaker* picoMaker =
    new StPicoDstMaker(2, inList, "picoDst");

  chain->Init();

  StPicoDst* picoDst = picoMaker->picoDst();

  int nEntries = picoMaker->chain()->GetEntries();
  if (nEvents > nEntries) nEvents = nEntries;

  cout << "[LambdaInvMassMin] entries=" << nEntries
       << " processing=" << nEvents << endl;

  // --- histogram ---
  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass;M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp  = 0.938272;   // proton mass [GeV]
  const double mpi = 0.139570;   // pion mass   [GeV]

  // --- event loop ---
  for (int ie=0; ie<nEvents; ie++) {
    chain->Clear();
    int iret = chain->Make(ie);
    if (iret) break;

    StPicoEvent* ev = picoDst->event();
    if (!ev) continue;

    int nTr = picoDst->numberOfTracks();

    // loop over proton candidates
    for (int i=0; i<nTr; i++) {
      StPicoTrack* trkP = picoDst->track(i);
      if (!trkP) continue;

      if (trkP->charge() <= 0) continue;        // p+
      if (fabs(trkP->nSigmaProton()) > 3.0) continue;
      if (trkP->pPt() < 0.3) continue;

      TVector3 pP = trkP->pMom();
      TLorentzVector lp;
      lp.SetXYZM(pP.x(), pP.y(), pP.z(), mp);

      // loop over pion candidates
      for (int j=0; j<nTr; j++) {
        if (j == i) continue;

        StPicoTrack* trkPi = picoDst->track(j);
        if (!trkPi) continue;

        if (trkPi->charge() >= 0) continue;     // pi-
        if (fabs(trkPi->nSigmaPion()) > 3.0) continue;
        if (trkPi->pPt() < 0.15) continue;

        TVector3 pPi = trkPi->pMom();
        TLorentzVector lpi;
        lpi.SetXYZM(pPi.x(), pPi.y(), pPi.z(), mpi);

        TLorentzVector lam = lp + lpi;
        hLam->Fill(lam.M());
      }
    }
  }

  // --- output ---
  TFile* fout = new TFile(outFile,"RECREATE");
  hLam->Write();
  fout->Close();

  TCanvas* c1 = new TCanvas("c1","Lambda invariant mass",800,600);
  hLam->Draw();

  cout << "[LambdaInvMassMin] done" << endl;
}
