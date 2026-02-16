#include <iostream>
#include <vector>

// ROOT
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TMath.h"

// STAR
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPhysicalHelixD.h"
#include "StThreeVectorF.h"
#include "StarClassLibrary/SystemOfUnits.h"

using namespace std;

// ==========================================================
// ヘリックス生成（CINT対応・最小）
// ==========================================================
StPhysicalHelixD makeHelix(const StPicoTrack* trk,
                           const TVector3& pv,
                           double B = 0.5) // Tesla
{
  StThreeVectorF origin(trk->origin().x(),
                        trk->origin().y(),
                        trk->origin().z());

  StThreeVectorF p(trk->pMom().x(),
                   trk->pMom().y(),
                   trk->pMom().z());

  return StPhysicalHelixD(p,
                          origin,
                          B * tesla,
                          trk->charge());
}

// ==========================================================
// Helix–Helix で V0 を作る
// ==========================================================
bool MakeV0PairHelix(const StPicoTrack* t1,
                     const StPicoTrack* t2,
                     const TVector3& pv,
                     TVector3& v0,
                     double& dca12)
{
  StPhysicalHelixD h1 = makeHelix(t1, pv);
  StPhysicalHelixD h2 = makeHelix(t2, pv);

  double s1 = 0, s2 = 0;
  dca12 = h1.pathLengths(h2, s1, s2);

  if (dca12 < 0 || dca12 > 1.0) return false;

  StThreeVectorD x1 = h1.at(s1);
  StThreeVectorD x2 = h2.at(s2);

  v0.SetXYZ(0.5 * (x1.x() + x2.x()),
            0.5 * (x1.y() + x2.y()),
            0.5 * (x1.z() + x2.z()));
  return true;
}

// ==========================================================
// メイン
// ==========================================================
void LambdaFromPicoHelix(const char* infile)
{
  // --- picoDst reader
  StPicoDstReader* reader = new StPicoDstReader(infile);
  reader->Init();

  StPicoDst* picoDst = reader->picoDst();

  // --- histogram
  TH1D* hMinv = new TH1D("hMinv",
                         "#Lambda invariant mass;M(p#pi^{-}) [GeV/c^{2}];Counts",
                         200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  // ======================================================
  // event loop
  // ======================================================
  Long64_t nEvents = reader->chain()->GetEntries();
  cout << "Total events = " << nEvents << endl;

  for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
    reader->readPicoEvent(iEvent);
    if (!picoDst) continue;

    StPicoEvent* event = picoDst->event();
    if (!event) continue;

    TVector3 pv(event->primaryVertex().x(),
                event->primaryVertex().y(),
                event->primaryVertex().z());

    int nTracks = picoDst->numberOfTracks();

    // ------------------------------------------
    // track loop (p, pi-)
    // ------------------------------------------
    for (int i = 0; i < nTracks; i++) {
      StPicoTrack* pTrk = picoDst->track(i);
      if (!pTrk) continue;
      if (pTrk->charge() <= 0) continue;

      for (int j = 0; j < nTracks; j++) {
        StPicoTrack* piTrk = picoDst->track(j);
        if (!piTrk) continue;
        if (piTrk->charge() >= 0) continue;

        // --- helix V0
        TVector3 v0;
        double dca12 = -1;

        if (!MakeV0PairHelix(pTrk, piTrk, pv, v0, dca12)) continue;

        // --- DCA to PV
        StPhysicalHelixD hp  = makeHelix(pTrk, pv);
        StPhysicalHelixD hpi = makeHelix(piTrk, pv);

        double dcaP  = hp.distance(pv);
        double dcaPi = hpi.distance(pv);

        if (dcaP < 0.5)  continue;
        if (dcaPi < 0.5) continue;

        // --- invariant mass
        TVector3 pp  = pTrk->pMom();
        TVector3 ppi = piTrk->pMom();

        double Ep  = sqrt(pp.Mag2()  + mp*mp);
        double Epi = sqrt(ppi.Mag2() + mpi*mpi);

        TVector3 pTot = pp + ppi;
        double M2 = (Ep + Epi)*(Ep + Epi) - pTot.Mag2();

        if (M2 > 0) hMinv->Fill(sqrt(M2));
      }
    }
  }

  // --- output
  TFile* fout = new TFile("lambda_minv.root", "RECREATE");
  hMinv->Write();
  fout->Close();

  cout << "Output written to lambda_minv.root" << endl;
}
