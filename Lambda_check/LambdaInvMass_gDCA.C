// LambdaInvMass_gDCA.C
// ROOT5 / CINT macro

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
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
class StPicoTrack;

StChain *chain = 0;

void LambdaInvMass_gDCA(Int_t nEvents=100,
                        const char* inList="test2.list",
                        const char* outFile="lambda_gDCA.root")
{
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

  cout << "[LambdaInvMass_gDCA] entries=" << nEntries
       << " processing=" << nEvents << endl;

  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass (gDCA cut);M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  const double dcaCut = 0.8;   // cm（まずはこのくらい）

  for (int ie=0; ie<nEvents; ie++) {
    std::cout<<"event No="<<ie<<std::endl;
    chain->Clear();
    if (chain->Make(ie)) break;

    int nTr = picoDst->numberOfTracks();

    for (int i=0; i<nTr; i++) {
      StPicoTrack* trkP = picoDst->track(i);
      if (!trkP) continue;

      // proton candidate
      if (trkP->charge() <= 0) continue;
      if (fabs(trkP->nSigmaProton()) > 3.0) continue;
      if (trkP->pPt() < 0.3) continue;
      TVector3 pv = picoDst->event()->primaryVertex();
      //if (trkP->gDCA() < dcaCut) continue;
      if (trkP->gDCA(pv.x(), pv.y(), pv.z()) < dcaCut) continue;

      //      TVector3 pP = trkP->pMom();
      TVector3 pP = trkP->gMom();
      TLorentzVector lp;
      lp.SetXYZM(pP.x(), pP.y(), pP.z(), mp);

      for (int j=0; j<nTr; j++) {
        if (j == i) continue;

        StPicoTrack* trkPi = picoDst->track(j);
        if (!trkPi) continue;

        // pion candidate
        if (trkPi->charge() >= 0) continue;
        if (fabs(trkPi->nSigmaPion()) > 3.0) continue;
        if (trkPi->pPt() < 0.15) continue;
        //if (trkPi->gDCA() < dcaCut) continue;
	if (trkPi->gDCA(pv.x(), pv.y(), pv.z()) < dcaCut) continue;

	//        TVector3 pPi = trkPi->pMom();
        TVector3 pPi = trkPi->gMom();
        TLorentzVector lpi;
        lpi.SetXYZM(pPi.x(), pPi.y(), pPi.z(), mpi);

        TLorentzVector lam = lp + lpi;
        hLam->Fill(lam.M());
      }
    }
  }

  TFile* fout = new TFile(outFile,"RECREATE");
  hLam->Write();
  fout->Close();

  TCanvas* c1 = new TCanvas("c1","Lambda (gDCA cut)",800,600);
  hLam->Draw();

  cout << "[LambdaInvMass_gDCA] done" << endl;
}
