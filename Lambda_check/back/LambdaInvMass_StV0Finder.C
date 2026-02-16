// LambdaInvMass_StV0Finder.C
// ACLiC only (ROOT5, bfc + StV0Finder)

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;

// STAR
#include "StChain.h"
#include "StMuDstMaker/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StEvent/StV0Vertex.h"

void LambdaInvMass_StV0Finder(const char* mudstList="mudst.list",
			      Int_t nEvents=1000)
{
  // ------------------------------------------------------------
  // Load bfc.C (THIS IS THE KEY)
  // ------------------------------------------------------------
  gROOT->ProcessLine(
    ".L $STAR/StRoot/StBFChain/BFC.C"
		     );

  TString chainOpt = "in,MuDst,V0Finder";
  TString cmd = Form(
		     "bfc(0,\"%s\",\"%s\")",
		     chainOpt.Data(), mudstList
		     );

  gROOT->ProcessLine(cmd);

  // ------------------------------------------------------------
  // Access chain and MuDst
  // ------------------------------------------------------------
  StChain* chain = StChain::GetChain();
  StMuDstMaker* muDstMaker =
    (StMuDstMaker*) chain->GetMaker("MuDst");

  if (!muDstMaker) {
    cout << "No MuDstMaker found" << endl;
    return;
  }

  // ------------------------------------------------------------
  // Histogram
  // ------------------------------------------------------------
  TH1D* hLam = new TH1D(
			"hLam",
			"Lambda invariant mass (StV0Finder);M_{p#pi^{-}} (GeV/c^{2})",
			200, 1.05, 1.25
			);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  // ------------------------------------------------------------
  // Event loop
  // ------------------------------------------------------------
  for (int ie=0; ie<nEvents; ie++) {
    if (chain->Make(ie)) break;

    StMuDst* muDst = muDstMaker->muDst();
    if (!muDst) continue;

    int nV0 = muDst->numberOfV0s();

    for (int iv0=0; iv0<nV0; iv0++) {
      StV0Vertex* v0 = muDst->v0(iv0);
      if (!v0) continue;

      // Lambda only
      if (v0->idTruth() != 3122) continue;

      // daughter momenta at SECONDARY vertex
      TVector3 pP  = v0->momentumPos();
      TVector3 pPi = v0->momentumNeg();

      TLorentzVector lp, lpi;
      lp.SetVectM(pP,  mp);
      lpi.SetVectM(pPi, mpi);

      hLam->Fill((lp+lpi).M());
    }
  }

  // ------------------------------------------------------------
  // Draw
  // ------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","Lambda",800,600);
  hLam->Draw();

  cout << "Done. Lambda peak should be visible." << endl;
}
