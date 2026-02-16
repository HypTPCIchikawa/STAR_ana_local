// LambdaInvMass_MinV0_Cuts.C
// picoDst + straight-line V0 finder with Run21-like cuts (approx w/o KF)
// ROOT5 / CINT compatible

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

// ===== utility: closest approach of two lines (x1 + s p1, x2 + t p2) =====
// returns dca and provides PCA points and mid-point (V0 vertex)
double closestApproachLines(const TVector3& x1, const TVector3& p1,
                            const TVector3& x2, const TVector3& p2,
                            TVector3& pca1, TVector3& pca2, TVector3& v0)
{
  TVector3 w0 = x1 - x2;
  double a = p1.Dot(p1);
  double b = p1.Dot(p2);
  double c = p2.Dot(p2);
  double d = p1.Dot(w0);
  double e = p2.Dot(w0);

  double denom = a*c - b*b;
  if (fabs(denom) < 1e-10) {
    pca1 = x1; pca2 = x2; v0 = 0.5*(x1+x2);
    return 999.;
  }

  double sc = (b*e - c*d) / denom;
  double tc = (a*e - b*d) / denom;

  pca1 = x1 + sc*p1;
  pca2 = x2 + tc*p2;
  v0   = 0.5*(pca1 + pca2);

  return (pca1 - pca2).Mag();
}

void LambdaInvMass_MinV0_Cuts(Int_t nEvents=30, const char* inList="test2.list")
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
  if (nEvents <= 0 || nEvents > nEntries) nEvents = nEntries;

  TH1D* hLam = new TH1D("hLam",
                        "#Lambda invariant mass;M_{p#pi^{-}} (GeV/c^{2});Counts",
                        200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  // ---- Run21 slide p.10 like cuts (those we can implement directly) ----
  const double cut_dca12_cm       = 1.0;   // dca12 < 1.0 cm
  const double cut_decayL_cm      = 5.0;   // decayLength > 5.0 cm
  const double cut_hitRatio       = 0.52;  // nHitFit/nHitPoss > 0.52

  // "chi2_topo < 5" and "chi2_primary > 10" are KFParticle-style outputs.
  // Here we keep placeholders / simple proxies.
  // You can tune these proxies once you see the Lambda peak.

  for (int ie=0; ie<nEvents; ie++) {
    std::cout<<"Events:"<<ie<<std::endl;

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

      // hit ratio
      int nFitP = p->nHitsFit();
      int nPossP = p->nHitsMax(); // in pico often "max" ~ "possible"
      if (nPossP > 0 && (double)nFitP/(double)nPossP <= cut_hitRatio) continue;

      TVector3 pMom = p->pMom();      // (proxy) momentum
      TVector3 pPos = p->origin();    // IMPORTANT: not PV

      for (int j=0; j<nTr; j++) {
        if (j==i) continue;
        StPicoTrack* pi = picoDst->track(j);
        if (!pi) continue;

        if (pi->charge() >= 0) continue;
        if (fabs(pi->nSigmaPion()) > 2) continue;

        int nFitPi = pi->nHitsFit();
        int nPossPi = pi->nHitsMax();
        if (nPossPi > 0 && (double)nFitPi/(double)nPossPi <= cut_hitRatio) continue;

        // DCA to PV proxy (your old cut). Keep for background suppression.
        if (pi->gDCA(pv.x(),pv.y(),pv.z()) < 0.8) continue;

        TVector3 piMom = pi->pMom();
        TVector3 piPos = pi->origin();

        TVector3 pca1, pca2, v0;
        double dca12 = closestApproachLines(pPos,pMom,piPos,piMom,pca1,pca2,v0);
        if (dca12 > cut_dca12_cm) continue;

        double decayL = (v0 - pv).Mag();
        if (decayL <= cut_decayL_cm) continue;

        // ---- proxies for KF chi2 cuts (optional, tune later) ----
        // pointing
        TVector3 v0Mom = pMom + piMom;
        double cosPA = (v0Mom.Dot(v0 - pv)) / (v0Mom.Mag() * (v0 - pv).Mag() + 1e-9);
        if (cosPA < 0.99) continue; // proxy for "chi2_topo < 5" (tighten/loosen later)

        // If you want to mimic chi2_primary>10, require daughters far from PV:
        // (These thresholds are *not* equivalent to KF chi2. Just a proxy.)
        if (p->gDCA(pv.x(),pv.y(),pv.z()) < 0.3) continue;
        if (pi->gDCA(pv.x(),pv.y(),pv.z()) < 0.3) continue;

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
