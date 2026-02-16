// LambdaInvMass_KFParticle_NoKFPTrack.C
// picoDst + KFParticle V0 (Lambda) reconstruction
// ROOT5 / CINT compatible (root4star)
// NOTE: DOES NOT use KFPTrack.h

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <vector>

// forward declarations (CINT style)
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;

StChain *chain = 0;

// KFParticle headers
#include "KFParticle.h"
#include "KFVertex.h"
//#include "KFParticleDatabase.h"

// ------------------------ helpers ------------------------
static inline bool GetTofBeta(const StPicoDst* picoDst, const StPicoTrack* trk, float &beta);
static inline float Mass2FromTof(float p, float beta);

static inline bool MakeKFParticleFromPico(const StPicoTrack* trk,
                                          const TVector3& pv,
                                          float bField,
                                          int pdg,
                                          KFParticle &out);

// ------------------------ main ---------------------------
void LambdaInvMass_KFParticle_NoKFPTrack(Int_t nEvents = 0,
                                         const Char_t* inFileList = "test2.list",
                                         const Char_t* outRoot = "LambdaInvMass_KFParticle.root")
{
  // STAR shared libs
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("KFParticle");

  // hist
  TH1D* hLam  = new TH1D("hLam",  "#Lambda #rightarrow p#pi^{-};M_{p#pi^{-}} (GeV/c^{2});counts", 400, 1.05, 1.25);
  TH1D* hALam = new TH1D("hALam", "#bar{#Lambda} #rightarrow #bar{p}#pi^{+};M_{#bar{p}#pi^{+}} (GeV/c^{2});counts", 400, 1.05, 1.25);

  // chain
  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inFileList, "picoDst");
  chain->Init();

  //  KFParticleDatabase db;

  Long64_t iEvent = 0;
  while(true)
    {
      if(nEvents > 0 && iEvent >= nEvents) break;
      if(chain->Make(iEvent)) break;

      StPicoDst* picoDst = picoMaker->picoDst();
      if(!picoDst) { ++iEvent; continue; }
      StPicoEvent* ev = picoDst->event();
      if(!ev) { ++iEvent; continue; }

      TVector3 pv = ev->primaryVertex();
      float bField = ev->bField();

      // PV -> KFVertex (rough covariance)
      KFVertex kfPV;
      {
	float xyz[3] = { (float)pv.X(), (float)pv.Y(), (float)pv.Z() };
	float cov[6] = { 0.01f*0.01f, 0.f,
			 0.01f*0.01f, 0.f, 0.f,
			 0.05f*0.05f };
	kfPV = KFVertex(xyz, cov, 0, 0);
      }

      // collect candidates
      struct Idx { int i; };
      std::vector<int> pPlus, piMinus, pMinus, piPlus;
      pPlus.reserve(1000); piMinus.reserve(2000); pMinus.reserve(1000); piPlus.reserve(2000);

      const int nTrk = picoDst->numberOfTracks();
      for(int i=0;i<nTrk;i++)
	{
	  StPicoTrack* trk = picoDst->track(i);
	  if(!trk) continue;

	  if(trk->nHitsFit()  < 15) continue;
	  if(trk->nHitsDedx() < 10) continue;

	  float nsP  = trk->nSigmaProton();
	  float nsPi = trk->nSigmaPion();

	  float beta=-1.f; bool hasTof = GetTofBeta(picoDst, trk, beta);
	  TVector3 pMom = trk->pMom();
	  float pMag = pMom.Mag();
	  float m2 = hasTof ? Mass2FromTof(pMag, beta) : 1e9;

	  bool isPion   = (TMath::Abs(nsPi) < 3.0) && (!hasTof || (m2 > -0.2 && m2 < 0.35));
	  bool isProton = (TMath::Abs(nsP)  < 3.0) && (!hasTof || (m2 >  0.5 && m2 < 1.35));
	  if(!isPion && !isProton) continue;

	  int q = trk->charge();
	  if(isProton){
	    if(q>0) pPlus.push_back(i);
	    if(q<0) pMinus.push_back(i);
	  }
	  if(isPion){
	    if(q<0) piMinus.push_back(i);
	    if(q>0) piPlus.push_back(i);
	  }
	}

      // Lambda: p+ + pi-
      for(size_t ip=0; ip<pPlus.size(); ip++)
	for(size_t ii=0; ii<piMinus.size(); ii++)
	  {
	    int iP  = pPlus[ip];
	    int iPi = piMinus[ii];
	    if(iP==iPi) continue;

	    StPicoTrack* trP  = picoDst->track(iP);
	    StPicoTrack* trPi = picoDst->track(iPi);
	    if(!trP || !trPi) continue;

	    KFParticle kfP, kfPi;
	    if(!MakeKFParticleFromPico(trP,  pv, bField,  2212, kfP)) continue;
	    if(!MakeKFParticleFromPico(trPi, pv, bField,  -211, kfPi)) continue;

	    KFParticle kfLam;
	    kfLam += kfP;
	    kfLam += kfPi;
	    kfLam.SetProductionVertex(kfPV);

	    if(kfLam.GetNDF() > 0 && kfLam.GetChi2()/kfLam.GetNDF() > 20.0) continue;

	    hLam->Fill(kfLam.GetMass());
	  }

      // Anti-Lambda: p- + pi+
      for(size_t ip=0; ip<pMinus.size(); ip++)
	for(size_t ii=0; ii<piPlus.size(); ii++)
	  {
	    int iP  = pMinus[ip];
	    int iPi = piPlus[ii];
	    if(iP==iPi) continue;

	    StPicoTrack* trP  = picoDst->track(iP);
	    StPicoTrack* trPi = picoDst->track(iPi);
	    if(!trP || !trPi) continue;

	    KFParticle kfP, kfPi;
	    if(!MakeKFParticleFromPico(trP,  pv, bField, -2212, kfP)) continue;
	    if(!MakeKFParticleFromPico(trPi, pv, bField,   211, kfPi)) continue;

	    KFParticle kfALam;
	    kfALam += kfP;
	    kfALam += kfPi;
	    kfALam.SetProductionVertex(kfPV);

	    if(kfALam.GetNDF() > 0 && kfALam.GetChi2()/kfALam.GetNDF() > 20.0) continue;

	    hALam->Fill(kfALam.GetMass());
	  }

      ++iEvent;
    }

  // save
  TFile* fout = new TFile(outRoot,"RECREATE");
  hLam->Write();
  hALam->Write();
  fout->Close();

  // draw
  TCanvas* c = new TCanvas("cKF","KFParticle Lambda",1000,450);
  c->Divide(2,1);
  c->cd(1); hLam->Draw();
  c->cd(2); hALam->Draw();

  std::cout << "Done. Output: " << outRoot << std::endl;

  chain->Finish();
}

// ------------------------ helper impl ---------------------
static inline bool GetTofBeta(const StPicoDst* picoDst, const StPicoTrack* trk, float &beta)
{
  beta = -1.f;
  int idx = trk->bTofPidTraitsIndex();
  if(idx < 0) return false;
  const StPicoBTofPidTraits* btof = picoDst->btofPidTraits(idx);
  if(!btof) return false;
  beta = btof->btofBeta();
  return (beta > 0.f);
}

static inline float Mass2FromTof(float p, float beta)
{
  if(beta<=0.f || beta>=1.f) return 1e9;
  return p*p*(1.f/(beta*beta) - 1.f);
}

static inline bool MakeKFParticleFromPico(const StPicoTrack* trk,
                                          const TVector3& pv,
                                          float bField,
                                          int pdg,
                                          KFParticle &out)
{
  // position: use track origin (DCA point) if available
  TVector3 r0 = trk->origin();

  // momentum: prefer global momentum at PV if available
  TVector3 p3 = trk->gMom(pv, bField);
  if(p3.Mag() <= 0) p3 = trk->pMom();

  float par[6] = { (float)r0.X(), (float)r0.Y(), (float)r0.Z(),
                   (float)p3.X(), (float)p3.Y(), (float)p3.Z() };

  // covariance (approx diagonal) — peak確認にはこれで十分
  float cov[21]; for(int i=0;i<21;i++) cov[i]=0.f;
  const float sx2 = 0.05f*0.05f; // 0.5 mm
  const float sp2 = 0.02f*0.02f; // 20 MeV/c

  cov[0]  = sx2; // xx
  cov[2]  = sx2; // yy
  cov[5]  = sx2; // zz
  cov[9]  = sp2; // pxpx
  cov[14] = sp2; // pypy
  cov[20] = sp2; // pzpz

  int q = trk->charge();

  // Create directly from (par,cov,q,pdg)
  out.Create(par, cov, q, pdg);
  return true;
}
