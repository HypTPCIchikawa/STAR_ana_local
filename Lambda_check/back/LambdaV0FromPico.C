// LambdaV0FromPico.C  (ROOT5/CINT macro)
// dE/dx nSigma only (p/pi) -> build Lambda candidates with simple V0 topology

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <vector>

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "StPhysicalHelixD.hh"

StChain* chain = 0;

// helper: distance from point X to line (P0 + t*dir)
static double distPointToLine(const TVector3& X, const TVector3& P0, const TVector3& dir)
{
  TVector3 u = dir;
  if (u.Mag() <= 0) return 1e9;
  u = u.Unit();
  return ( (X-P0).Cross(u) ).Mag();
}

struct Cand {
  int idx;
  TVector3 p;           // momentum (GeV/c)
  StPhysicalHelixD h;   // helix
  TVector3 x0;          // origin (DCA point)
};


void LambdaV0FromPico(Int_t nEvents=100,
                      const char* inList="test2.list",
                      const char* outFile="LambdaQA.root")
{
  // --- load STAR libs
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inList, "picoDst");

  // --- output
  TFile* fout = new TFile(outFile,"RECREATE");
  TH1D* hM     = new TH1D("hM", "p#pi^{-} invariant mass;M(p#pi) [GeV/c^{2}];counts", 600, 1.05, 1.25);
  TH1D* hDca12 = new TH1D("hDca12","DCA(p,#pi) at closest approach;DCA [cm];counts", 400, 0.0, 4.0);
  TH1D* hCos   = new TH1D("hCos",  "cos(pointing);cos;counts", 400, -1.0, 1.0);
  TH1D* hDL    = new TH1D("hDL",   "decay length |V0-PV|;L [cm];counts", 600, 0.0, 60.0);
  TH1D* hDcaV0 = new TH1D("hDcaV0","DCA(V0 line, PV);DCA [cm];counts", 400, 0.0, 4.0);

  // --- init chain
  chain->Init();

  const double mp  = 0.9382720813; // GeV/c^2
  const double mpi = 0.13957039;

  // ---- analysis cuts (tune later)
  const int   minNHitsFit  = 15;
  const int   minNHitsDedx = 10;
  const double minPt = 0.15;

  const double nsigP_max  = 3.0;   // |nSigmaProton|
  const double nsigPi_max = 3.0;   // |nSigmaPion|

  // V0 topology (start loose, tighten later)
  const double dca12_max = 1.0;    // cm  DCA between daughters at closest approach
  const double decayL_min = 2.0;   // cm  |V0-PV|
  const double cos_min    = 0.995; // pointing
  const double dcaV0_max  = 1.0;   // cm  DCA(V0 line, PV)

  // ---- event loop
  int nDone = 0;
  for (int i=0; i<nEvents; i++) {
    chain->Clear();
    int ierr = chain->Make(i);
    if (ierr) break;

    StPicoDst* picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    StPicoEvent* ev = picoDst->event();
    if (!ev) continue;

    // PV & B
    //StThreeVectorF pvF = ev->primaryVertex();
    //TVector3 pv(pvF.x(), pvF.y(), pvF.z());
    TVector3 pv = ev->primaryVertex();   
    const double B = ev->bField(); // STAR usually in kG

    // collect candidates first (FAST)
    std::vector<Cand> protons;
    std::vector<Cand> pions;

    const int nTrk = picoDst->numberOfTracks();
    for (int it=0; it<nTrk; it++) {
      StPicoTrack* trk = picoDst->track(it);
      if (!trk) continue;

      if (trk->nHitsFit()  < minNHitsFit)  continue;
      if (trk->nHitsDedx() < minNHitsDedx) continue;

      TVector3 p3(trk->gMom().x(), trk->gMom().y(), trk->gMom().z());
      if (p3.Pt() < minPt) continue;

      // origin (usually at DCA to PV)
      TVector3 x0(trk->origin().x(), trk->origin().y(), trk->origin().z());

      // simple DCA to PV (approx)
      double dcaPV = (x0 - pv).Mag();

      // build helix
      StThreeVectorD oD(x0.X(), x0.Y(), x0.Z());
      StThreeVectorD mD(p3.X(), p3.Y(), p3.Z());
      StPhysicalHelixD hel(mD, oD, B, trk->charge());

      // dE/dx PID only
      const double nsP  = trk->nSigmaProton();
      const double nsPi = trk->nSigmaPion();

      // ---- Lambda: p(+) + pi(-)
      if (trk->charge() > 0 && TMath::Abs(nsP) < nsigP_max) {
        // you can optionally require "not pion": |nSigmaPion|>X, etc.
        Cand c; c.idx=it; c.p=p3; c.h=hel; c.x0=x0;
        protons.push_back(c);
      }
      if (trk->charge() < 0 && TMath::Abs(nsPi) < nsigPi_max) {
        Cand c; c.idx=it; c.p=p3; c.h=hel; c.x0=x0;
        pions.push_back(c);
      }
    }

    // combine ONLY candidates (FAST)
    for (size_t ip=0; ip<protons.size(); ip++) {
      for (size_t ij=0; ij<pions.size(); ij++) {

        const Cand& P  = protons[ip];
        const Cand& Pi = pions[ij];

        // closest approach between helices
	std::pair<double,double> s = P.h.pathLengths(Pi.h);
        StThreeVectorD xP  = P.h.at(s.first);
        StThreeVectorD xPi = Pi.h.at(s.second);

        TVector3 rP(xP.x(), xP.y(), xP.z());
        TVector3 rPi(xPi.x(), xPi.y(), xPi.z());
        double dca12 = (rP - rPi).Mag();
        if (dca12 > dca12_max) continue;

        // V0 vertex: midpoint at closest approach
        TVector3 v0 = 0.5*(rP + rPi);

        // decay length
        double decayL = (v0 - pv).Mag();
        if (decayL < decayL_min) continue;

        // V0 momentum
        TVector3 pV0 = P.p + Pi.p;
        if (pV0.Mag() <= 0) continue;

        // pointing
        double cosPA = ( (v0 - pv).Dot(pV0) ) / ( (v0 - pv).Mag() * pV0.Mag() );
        if (cosPA < cos_min) continue;

        // DCA(V0 line, PV)
        double dcaV0 = distPointToLine(pv, v0, pV0);
        if (dcaV0 > dcaV0_max) continue;

        // invariant mass
        TLorentzVector lvP, lvPi;
        lvP .SetVectM(P.p , mp);
        lvPi.SetVectM(Pi.p, mpi);
        double m = (lvP + lvPi).M();

        // fill
        hM->Fill(m);
        hDca12->Fill(dca12);
        hDL->Fill(decayL);
        hCos->Fill(cosPA);
        hDcaV0->Fill(dcaV0);
      }
    }

    nDone++;
    if (nDone==1) {
      std::cout << "First event: nTracks=" << nTrk
                << " pCand=" << protons.size()
                << " piCand=" << pions.size()
                << "  PVz=" << pv.Z() << std::endl;
    }
  }

  fout->Write();
  fout->Close();
  std::cout << "Wrote: " << outFile << "  (events processed=" << nDone << ")" << std::endl;
}
