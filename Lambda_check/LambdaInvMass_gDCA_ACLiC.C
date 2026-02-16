// LambdaInvMass_gDCA_ACLiC.C
// ROOT5 + ACLiC version

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>

#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

using namespace std;

StChain* chain = 0;

void LambdaInvMass_gDCA_ACLiC(const char* inList="test2.list",
			      const char* outFile="lambda_gDCA.root",
			      Int_t nEvents=100)
{
  // --- load shared libraries (explicit for ACLiC)
  gSystem->Load("StarClassLibrary");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  chain = new StChain("StChain");

  StPicoDstMaker* picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  chain->Init();

  StPicoDst* picoDst = picoMaker->picoDst();

  int nEntries = picoMaker->chain()->GetEntries();
  if (nEvents > nEntries) nEvents = nEntries;

  cout << "[LambdaInvMass_gDCA_ACLiC] entries=" << nEntries
       << " processing=" << nEvents << endl;

  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass (gDCA cut);M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;
  const double dcaCut = 0.8; // cm

  for (int ie=0; ie<nEvents; ie++) {
    chain->Clear();
    if (chain->Make(ie)) break;

    StPicoEvent* event = picoDst->event();
    if (!event) continue;

    TVector3 pv = event->primaryVertex();
    int nTr = picoDst->numberOfTracks();

    for (int i=0; i<nTr; i++) {
      StPicoTrack* p = picoDst->track(i);
      if (!p) continue;

      if (p->charge() <= 0) continue;
      if (fabs(p->nSigmaProton()) > 3.0) continue;
      if (p->pPt() < 0.3) continue;
      if (p->gDCA(pv.x(), pv.y(), pv.z()) < dcaCut) continue;

      TVector3 pMom = p->gMom();
      TLorentzVector lp;
      lp.SetXYZM(pMom.x(), pMom.y(), pMom.z(), mp);

      for (int j=0; j<nTr; j++) {
        if (j == i) continue;

        StPicoTrack* pi = picoDst->track(j);
        if (!pi) continue;

        if (pi->charge() >= 0) continue;
        if (fabs(pi->nSigmaPion()) > 3.0) continue;
        if (pi->pPt() < 0.15) continue;
        if (pi->gDCA(pv.x(), pv.y(), pv.z()) < dcaCut) continue;

        TVector3 piMom = pi->gMom();
        TLorentzVector lpi;
        lpi.SetXYZM(piMom.x(), piMom.y(), piMom.z(), mpi);

        hLam->Fill((lp + lpi).M());
      }
    }
  }

  TFile* fout = new TFile(outFile, "RECREATE");
  hLam->Write();
  fout->Close();

  TCanvas* c1 = new TCanvas("c1","Lambda (gDCA cut)",800,600);
  hLam->Draw();

  cout << "[LambdaInvMass_gDCA_ACLiC] done" << endl;
}
