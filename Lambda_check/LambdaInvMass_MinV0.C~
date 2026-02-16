// LambdaInvMass_MinV0.C
// picoDst + straight-line V0 finder
// ROOT5 / CINT compatible

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

// ===== utility: closest approach of two lines =====
double twoTrackDCA(const TVector3& x1, const TVector3& p1,
                   const TVector3& x2, const TVector3& p2)
{
  TVector3 w0 = x1 - x2;
  double a = p1.Dot(p1);
  double b = p1.Dot(p2);
  double c = p2.Dot(p2);
  double d = p1.Dot(w0);
  double e = p2.Dot(w0);

  double denom = a*c - b*b;
  if (fabs(denom) < 1e-5) return 999.;

  double sc = (b*e - c*d) / denom;
  double tc = (a*e - b*d) / denom;

  TVector3 dvec = w0 + sc*p1 - tc*p2;
  return dvec.Mag();
}

void LambdaInvMass_MinV0(Int_t nEvents=10,
                         const char* inList="test2.list")
{
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  chain = new StChain();
  StPicoDstMaker* picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  chain->Init();
  StPicoDst* picoDst = picoMaker->picoDst();

  int nEntries = picoMaker->chain()->GetEntries();
  if (nEvents > nEntries) nEvents = nEntries;

  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass;M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp = 0.938272;
  const double mpi = 0.139570;

  for (int ie=0; ie<nEvents; ie++) {
    std::cout<<"event:"<<ie<<std::endl;

    chain->Clear();
    if (chain->Make(ie)) break;

    StPicoEvent* evt = picoDst->event();
    TVector3 pv = evt->primaryVertex();

    int nTr = picoDst->numberOfTracks();

    for (int i=0; i<nTr; i++) {
      StPicoTrack* p = picoDst->track(i);
      if (!p) continue;
      if (p->charge() <= 0) continue;
      if (fabs(p->nSigmaProton()) > 2) continue;
      //if (p->gDCA(pv.x(),pv.y(),pv.z()) < 0.8) continue;

      //      TVector3 pMom = p->gMom();
      TVector3 pMom = p->pMom();
      TVector3 pPos = pv; // 直線近似（最小）

      for (int j=0; j<nTr; j++) {
        if (j==i) continue;
        StPicoTrack* pi = picoDst->track(j);
        if (!pi) continue;
        if (pi->charge() >= 0) continue;
        if (fabs(pi->nSigmaPion()) > 2) continue;
        if (pi->gDCA(pv.x(),pv.y(),pv.z()) < 0.8) continue;

	//        TVector3 piMom = pi->gMom();
        TVector3 piMom = pi->pMom();
        TVector3 piPos = pv;

        double dca = twoTrackDCA(pPos,pMom,piPos,piMom);
        if (dca > 1.0) continue;

        TLorentzVector lp, lpi;
        lp.SetVectM(pMom, mp);
        lpi.SetVectM(piMom, mpi);

        hLam->Fill((lp+lpi).M());
      }
    }
  }

  TCanvas* c = new TCanvas("c","Lambda",800,600);
  hLam->Draw();
}
