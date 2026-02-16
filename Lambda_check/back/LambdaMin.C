// LambdaMin.C
// Minimal Λ(→ p + π-) reconstruction from PicoDst using StPhysicalHelixD
//
// Usage example:
// root -l -b -q 'LambdaMin.C(1000,"my.list","lambda_min.root")'
//
// Notes:
// - Uses *global* tracks (gMom, origin) and StPhysicalHelixD.
// - Assumes event->bField() is in kG (STAR convention). If your pico stores Tesla,
//   convert: 1 T = 10 kG.
// - This is a "minimum" example: you will likely tighten quality/PID/topology cuts.

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

// STAR helix classes
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

static inline double sqr(double x) { return x * x; }

static inline TLorentzVector makeLV(const StThreeVectorD& p3, double mass)
{
  TLorentzVector v;
  const double e = std::sqrt(p3.mag2() + mass * mass);
  v.SetPxPyPzE(p3.x(), p3.y(), p3.z(), e);
  return v;
}

void LambdaMin(Long64_t nEvents = 100000000,
               const char* inputList = "pico.list",
               const char* outFile   = "lambda_min.root")
{
  // masses (GeV/c^2)
  const double mP  = 0.9382720813;
  const double mPi = 0.13957039;

  // --- (very loose) analysis cuts ---
  const int    kMinNHitsFit = 15;
  const double kMaxDcaToPV  = 3.0;   // cm (loose, per daughter)
  const double kMaxPairDca  = 1.0;   // cm (DCA between daughters at closest approach)
  const double kMinV0R      = 0.5;   // cm (decay radius)
  const double kMaxV0R      = 200.0; // cm
  const double kMinCosPA    = 0.995; // pointing angle (loose-ish)
  const double kMaxNSigmaP  = 3.0;
  const double kMaxNSigmaPi = 3.0;

  // --- output ---
  TFile* fout = TFile::Open(outFile, "RECREATE");
  TH1D* hMlam = new TH1D("hMlam", "p#pi^{-} invariant mass;M_{p#pi} (GeV/c^{2});Counts", 400, 1.05, 1.25);
  TH1D* hMal  = new TH1D("hMal",  "#bar{#Lambda} (#bar{p}#pi^{+}) invariant mass;M_{#bar{p}#pi} (GeV/c^{2});Counts", 400, 1.05, 1.25);

  // --- PicoDst maker ---
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inputList);
  picoMaker->Init();

  Long64_t nProcessed = 0;

  while (true)
    {
      if (nEvents > 0 && nProcessed >= nEvents) break;

      const int iread = picoMaker->Make();
      if (iread) break; // end of chain / error

      StPicoDst* dst = picoMaker->picoDst();
      if (!dst) continue;

      StPicoEvent* ev = dst->event();
      if (!ev) continue;

      // Primary vertex
      const StThreeVectorD pv(ev->primaryVertex().x(),
			      ev->primaryVertex().y(),
			      ev->primaryVertex().z());

      // B field (STAR often uses kG; if Tesla, multiply by 10)
      const double bField = ev->bField(); // assume kG

      const int nTrk = dst->numberOfTracks();
      if (nTrk < 2) { nProcessed++; continue; }

      // Preselect indices for p/pi candidates by charge and PID
      std::vector<int> idxPPlus, idxPiMinus, idxPMinus, idxPiPlus;
      idxPPlus.reserve(nTrk);
      idxPiMinus.reserve(nTrk);
      idxPMinus.reserve(nTrk);
      idxPiPlus.reserve(nTrk);

      for (int i = 0; i < nTrk; ++i)
	{
	  StPicoTrack* trk = dst->track(i);
	  if (!trk) continue;
	  if (!trk->isPrimary()) {
	    // This flag is about "primary track container"; we still use global kinematics below.
	    // Keep it loose; many analyses use global tracks for V0 daughters.
	  }

	  if (trk->nHitsFit() < kMinNHitsFit) continue;

	  // track origin (global helix origin) and global momentum
	  const TVector3 o3 = trk->origin();
	  const TVector3 p3 = trk->gMom(pv, bField); // global momentum at DCA to PV (typical usage)
	  // If your StPicoTrack doesn't have gMom(pv,b), use: trk->gMom() or trk->gMom(ev->primaryVertex(), bField).

	  const StThreeVectorD origin(o3.X(), o3.Y(), o3.Z());
	  const StThreeVectorD mom(p3.X(), p3.Y(), p3.Z());

	  // DCA of track to PV (approx from origin - pv). Some pico formats also provide dca.
	  const double dcaPV = (origin - pv).mag();
	  if (dcaPV > kMaxDcaToPV) continue;

	  const int q = trk->charge();

	  const double nsp = trk->nSigmaProton();
	  const double nsi = trk->nSigmaPion();

	  // very loose PID buckets
	  if (std::fabs(nsp) < kMaxNSigmaP) {
	    if (q > 0) idxPPlus.push_back(i);
	    if (q < 0) idxPMinus.push_back(i);
	  }
	  if (std::fabs(nsi) < kMaxNSigmaPi) {
	    if (q < 0) idxPiMinus.push_back(i);
	    if (q > 0) idxPiPlus.push_back(i);
	  }
	}

      auto buildCandidates = [&](const std::vector<int>& idxP,
				 const std::vector<int>& idxPi,
				 TH1D* hMass,
				 bool isAnti)
	{
	  for (int ip : idxP)
	    {
	      StPicoTrack* tp = dst->track(ip);
	      if (!tp) continue;

	      const TVector3 op = tp->origin();
	      const TVector3 pp = tp->gMom(pv, bField);
	      const int qP = tp->charge();

	      const StThreeVectorD oP(op.X(), op.Y(), op.Z());
	      const StThreeVectorD pP(pp.X(), pp.Y(), pp.Z());

	      // Helix for p (or pbar)
	      StPhysicalHelixD hP(pP, oP, bField, qP);

	      for (int ii : idxPi)
		{
		  if (ii == ip) continue;
		  StPicoTrack* tpi = dst->track(ii);
		  if (!tpi) continue;

		  const TVector3 opi = tpi->origin();
		  const TVector3 ppi = tpi->gMom(pv, bField);
		  const int qPi = tpi->charge();

		  const StThreeVectorD oPi(opi.X(), opi.Y(), opi.Z());
		  const StThreeVectorD pPi(ppi.X(), ppi.Y(), ppi.Z());

		  StPhysicalHelixD hPi(pPi, oPi, bField, qPi);

		  // Closest approach between helices:
		  // pathLengths() returns the path lengths along each helix at the DCA points.
		  const std::pair<double,double> ss = hP.pathLengths(hPi);

		  const StThreeVectorD rP  = hP.at(ss.first);
		  const StThreeVectorD rPi = hPi.at(ss.second);

		  const double dcaPair = (rP - rPi).mag();
		  if (dcaPair > kMaxPairDca) continue;

		  // V0 vertex as midpoint of the two DCA points
		  const StThreeVectorD v0 = 0.5 * (rP + rPi);

		  // decay radius
		  const double rxy = std::sqrt(sqr(v0.x()) + sqr(v0.y()));
		  if (rxy < kMinV0R || rxy > kMaxV0R) continue;

		  // Daughter momenta at the vertex (at those path lengths)
		  const StThreeVectorD pP_at  = hP.momentumAt(ss.first,  bField);
		  const StThreeVectorD pPi_at = hPi.momentumAt(ss.second, bField);

		  const StThreeVectorD pV0 = pP_at + pPi_at;

		  // Pointing angle
		  const StThreeVectorD r = (v0 - pv);
		  const double cosPA = (r.mag() > 0 && pV0.mag() > 0) ? (r.dot(pV0) / (r.mag() * pV0.mag())) : -999.0;
		  if (cosPA < kMinCosPA) continue;

		  // Invariant mass
		  const TLorentzVector lvP  = makeLV(pP_at,  mP);
		  const TLorentzVector lvPi = makeLV(pPi_at, mPi);
		  const double m = (lvP + lvPi).M();

		  hMass->Fill(m);
		}
	    }
	};

      // Λ: p+ + pi-
      buildCandidates(idxPPlus, idxPiMinus, hMlam, false);
      // anti-Λ: p- + pi+
      buildCandidates(idxPMinus, idxPiPlus, hMal, true);

      nProcessed++;
    }

  fout->Write();
  fout->Close();

  // Cleanup maker
  picoMaker->Finish();
  delete picoMaker;

  printf("Done. Processed events: %lld\n", (long long)nProcessed);
  printf("Output: %s\n", outFile);
}
