// LambdaInvMass_V0.C  (ROOT5 / CINT)

#include <TROOT.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

// forward declarations
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoV0;

StChain* chain = 0;

void LambdaInvMass_V0(Int_t nEvents=1000,
                      const char* inList="test2.list")
{
  // STAR standard libs
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

  // Histogram
  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass (StPicoV0);M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    chain->Clear();
    chain->Make(iEvent);

    TClonesArray* v0Array = picoDst->v0s();
    int nV0 = v0Array->GetEntriesFast();

    for (int i = 0; i < nV0; i++) {
      StPicoV0* v0 = (StPicoV0*)v0Array->At(i);
      if (!v0) continue;

      // momenta at secondary vertex
      TVector3 pPos = v0->momPos();
      TVector3 pNeg = v0->momNeg();

      // Lambda hypothesis: p + pi-
      TLorentzVector lp, lpi;
      lp.SetVectM(pPos, mp);
      lpi.SetVectM(pNeg, mpi);

      double mass = (lp + lpi).M();
      hLam->Fill(mass);
    }
  }

  // Draw
  TCanvas* c1 = new TCanvas("c1","Lambda",800,600);
  hLam->Draw();
}
