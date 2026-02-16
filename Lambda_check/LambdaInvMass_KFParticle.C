// LambdaInvMass_KFParticle.C
// picoDst + KFParticle V0 (Lambda) reconstruction
// ROOT5 / CINT compatible (root4star)

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>
#include <vector>

// STAR classes: forward declarations (CINT style)
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;

StChain *chain = 0;

// ---- KFParticle headers (needed to compile in CINT) ----
#include "KFParticle.h"
#include "KFVertex.h"
#include "StRoot/KFParticle/KFParticleBase/KFPTrack.h"
#include "KFParticleDatabase.h"

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------
static inline bool GetTofBeta(const StPicoDst* picoDst, const StPicoTrack* trk, float &beta);
static inline float Mass2FromTof(float p, float beta);
static inline bool MakeKFPTrackFromPico(const StPicoTrack* trk,
                                        const TVector3& pv,
                                        float bField,
                                        KFPTrack &out);

// ------------------------------------------------------------

void LambdaInvMass_KFParticle(Int_t nEvents = 0,
                              const Char_t* inFileList = "test2.list",
                              const Char_t* outRoot = "LambdaInvMass_KFParticle.root")
{
  // Load STAR shared libraries
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("KFParticle");   // KFParticle libs

  // Histograms
  TH1D* hLam   = new TH1D("hLam",   "#Lambda #rightarrow p#pi^{-};M_{p#pi^{-}} (GeV/c^{2});counts", 400, 1.05, 1.25);
  TH1D* hALam  = new TH1D("hALam",  "#bar{#Lambda} #rightarrow #bar{p}#pi^{+};M_{#bar{p}#pi^{+}} (GeV/c^{2});counts", 400, 1.05, 1.25);

  // Makers
  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inFileList, "picoDst");
  chain->Init();

  // KF database (masses)
  KFParticleDatabase db;

  Long64_t iEvent = 0;
  while (true)
    {
      if (nEvents > 0 && iEvent >= nEvents) break;
      if (chain->Make(iEvent)) break;

      StPicoDst* picoDst = picoMaker->picoDst();
      if (!picoDst) { ++iEvent; continue; }

      StPicoEvent* ev = picoDst->event();
      if (!ev) { ++iEvent; continue; }

      // PV and B-field
      TVector3 pv = ev->primaryVertex();
      float bField = ev->bField();

      // PV as KFVertex (covariance is approximated here)
      KFVertex kfPV;
      {
	float xyz[3] = { (float)pv.X(), (float)pv.Y(), (float)pv.Z() };

	// covariance: (xx, xy, yy, xz, yz, zz) style
	// rough values are fine for first peak finding
	float cov[6] = {
	  0.01f*0.01f, 0.f,
	  0.01f*0.01f, 0.f, 0.f,
        0.05f*0.05f
	};
	kfPV = KFVertex(xyz, cov, 0, 0);
      }

      // Collect p / pi candidates to reduce combinatorics
      struct Cand {
	int idx;
	KFPTrack kfp;
	int q;
      };
      std::vector<Cand> pList, piList, apList, pipList;
      pList.reserve(2000);
      piList.reserve(2000);
      apList.reserve(2000);
      pipList.reserve(2000);

      const int nTrk = picoDst->numberOfTracks();
      for (int i = 0; i < nTrk; i++)
	{
	  StPicoTrack* trk = picoDst->track(i);
	  if (!trk) continue;

	  // Track quality (adjust as needed)
	  if (trk->nHitsFit() < 15) continue;
	  if (trk->nHitsDedx() < 10) continue;

	  // TPC PID
	  float nsP  = trk->nSigmaProton();
	  float nsPi = trk->nSigmaPion();

	  // TOF m^2 PID (optional, because you said beta exists)
	  float beta = -1.f;
	  bool hasTof = GetTofBeta(picoDst, trk, beta);

	  // Use p at PV (primary momentum) only for PID rough; KF uses its own (from gMom/origin)
	  TVector3 pMomPrim = trk->pMom();
	  float pMag = pMomPrim.Mag();
	  float m2 = hasTof ? Mass2FromTof(pMag, beta) : 1e9;

	  // Loose PID first (tune later)
	  bool isPion   = (TMath::Abs(nsPi) < 3.0) && (!hasTof || (m2 > -0.2 && m2 < 0.35));
	  bool isProton = (TMath::Abs(nsP)  < 3.0) && (!hasTof || (m2 >  0.5 && m2 < 1.35));
	  if (!isPion && !isProton) continue;

	  // Build KFPTrack (position/momentum at "track origin" + gMom at PV)
	  KFPTrack kfp;
	  if (!MakeKFPTrackFromPico(trk, pv, bField, kfp)) continue;

	  int q = trk->charge();
	  if (isProton)
	    {
	      if (q > 0) pList.push_back({i, kfp, q});
	      else       apList.push_back({i, kfp, q}); // anti-proton
	    }
	  if (isPion)
	    {
	      if (q < 0) piList.push_back({i, kfp, q});   // pi-
	      else       pipList.push_back({i, kfp, q});  // pi+
	    }
	}

      // ---- Lambda: p(+) + pi(-) ----
      for (size_t ip = 0; ip < pList.size(); ip++)
	for (size_t iPi = 0; iPi < piList.size(); iPi++)
	  {
	    const Cand& cp  = pList[ip];
	    const Cand& cpi = piList[iPi];
	    if (cp.idx == cpi.idx) continue;

	    KFParticle kfP, kfPi;
	    kfP.Create(cp.kfp,  db.GetParticleMass(2212),  2212);
	    kfPi.Create(cpi.kfp, db.GetParticleMass(-211), -211);

	    KFParticle kfLam;
	    kfLam += kfP;
	    kfLam += kfPi;

	    kfLam.SetProductionVertex(kfPV);

	    // Very loose quality (tighten later)
	    if (kfLam.GetNDF() > 0 && kfLam.GetChi2() / kfLam.GetNDF() > 20.0) continue;

	    float m = kfLam.GetMass();
	    hLam->Fill(m);
	  }

      // ---- Anti-Lambda: anti-p(-) + pi(+) ----
      for (size_t iap = 0; iap < apList.size(); iap++)
	for (size_t iPi = 0; iPi < pipList.size(); iPi++)
	  {
	    const Cand& cap  = apList[iap];
	    const Cand& cpi  = pipList[iPi];
	    if (cap.idx == cpi.idx) continue;

	    KFParticle kfAP, kfPi;
	    kfAP.Create(cap.kfp,  db.GetParticleMass(-2212), -2212);
	    kfPi.Create(cpi.kfp,  db.GetParticleMass(211),    211);

	    KFParticle kfALam;
	    kfALam += kfAP;
	    kfALam += kfPi;

	    kfALam.SetProductionVertex(kfPV);

	    if (kfALam.GetNDF() > 0 && kfALam.GetChi2() / kfALam.GetNDF() > 20.0) continue;

	    float m = kfALam.GetMass();
	    hALam->Fill(m);
	  }

      ++iEvent;
    }

  // Save
  TFile* fout = new TFile(outRoot, "RECREATE");
  hLam->Write();
  hALam->Write();
  fout->Close();

  // Quick draw
  TCanvas* c1 = new TCanvas("c1","Lambda KFParticle",1000,450);
  c1->Divide(2,1);
  c1->cd(1); hLam->Draw();
  c1->cd(2); hALam->Draw();

  std::cout << "Done. Output = " << outRoot << std::endl;

  chain->Finish();
}

// ------------------------------------------------------------
// Helper implementations (need picoDst classes at runtime via dictionary)
// ------------------------------------------------------------
static inline bool GetTofBeta(const StPicoDst* picoDst, const StPicoTrack* trk, float &beta)
{
  beta = -1.f;
  int idx = trk->bTofPidTraitsIndex();
  if (idx < 0) return false;
  const StPicoBTofPidTraits* btof = picoDst->btofPidTraits(idx);
  if (!btof) return false;
  beta = btof->btofBeta();
  return (beta > 0.f);
}

static inline float Mass2FromTof(float p, float beta)
{
  if (beta <= 0.f || beta >= 1.f) return 1e9;
  return p*p*(1.f/(beta*beta) - 1.f);
}

static inline bool MakeKFPTrackFromPico(const StPicoTrack* trk,
                                        const TVector3& pv,
                                        float bField,
                                        KFPTrack &out)
{
  // "origin" is usually the (global) DCA point used in pico track representation.
  // Momentum at PV: use gMom(pv,bField) if available; fallback to pMom()
  TVector3 r0 = trk->origin();

  TVector3 p3;
  // gMom(px,py,pz,bField) exists in many pico versions; use the vector version first
  p3 = trk->gMom(pv, bField);
  if (p3.Mag() <= 0) p3 = trk->pMom();

  out.SetX(r0.X());
  out.SetY(r0.Y());
  out.SetZ(r0.Z());
  out.SetPx(p3.X());
  out.SetPy(p3.Y());
  out.SetPz(p3.Z());
  out.SetQ(trk->charge());

  // Approximate covariance (diagonal only). Enough to see a peak first.
  float C[21] = {0};
  const float sx2 = 0.05f*0.05f; // 0.5 mm
  const float sp2 = 0.02f*0.02f; // 20 MeV/c

  C[0]  = sx2; // xx
  C[2]  = sx2; // yy
  C[5]  = sx2; // zz
  C[9]  = sp2; // pxpx
  C[14] = sp2; // pypy
  C[20] = sp2; // pzpz

  out.SetCovarianceMatrix(C);
  out.SetChi2(trk->chi2());
  out.SetNDF(TMath::Max(1, trk->nHitsFit() - 5));

  return true;
}
