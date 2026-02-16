//============================================================
//  LambdaFromPico_Linear.C
//
//  - Read picoDst with StPicoDstMaker (works in your environment)
//  - Build Lambda candidates from (p, pi-) pairs
//  - Reconstruct V0 vertex with simple 2-line PCA (linear approx)
//  - Apply loose topological + PID cuts
//  - Save QA histograms to OutputDir/lambdaQA.root
//
//  Run example:
//    root4star -l
//    .L LambdaFromPico_Linear.C
//    Analysis(2000, "200pico.list", "./", "./LocalOut")
//
//  Or take all events (dangerously long):
//    Analysis(0, "200pico.list", "./", "./LocalOut")
//
//============================================================

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TVector3.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain* chain = 0;

//-----------------------------
// helper: PCA of two lines
// line1: r1 + s*u, line2: r2 + t*v (u,v should be unit vectors)
// returns closest points p1,p2 and their distance dca
//-----------------------------
static bool LineLinePCA(const TVector3& r1, const TVector3& u,
                        const TVector3& r2, const TVector3& v,
                        TVector3& p1, TVector3& p2, double& dca)
{
  TVector3 w0 = r1 - r2;
  double a = u.Dot(u); // ~1
  double b = u.Dot(v);
  double c = v.Dot(v); // ~1
  double d = u.Dot(w0);
  double e = v.Dot(w0);
  double denom = a*c - b*b;

  if (std::fabs(denom) < 1e-12) return false; // parallel (or almost)

  double s = (b*e - c*d)/denom;
  double t = (a*e - b*d)/denom;

  p1 = r1 + s*u;
  p2 = r2 + t*v;
  dca = (p1 - p2).Mag();
  return true;
}

// distance from point p to infinite line r0 + s*u (u must be unit)
static double DcaPointToLine(const TVector3& p, const TVector3& r0, const TVector3& u)
{
  return ((p - r0).Cross(u)).Mag();
}

//============================================================
// Main
//============================================================
void Analysis(Int_t nEvents = 2000,
              TString InputFileList = "test.list",
	      //              TString /*Unused*/ = "./",
              TString OutputDir = "./LocalOut")
{
  if (nEvents == 0) nEvents = 100000000; // "all"

  // Load libraries
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  // Make sure output dir exists
  gSystem->mkdir(OutputDir, kTRUE);

  // Histograms
  TH1D* hMassLam = new TH1D("hMassLam","Lambda (p#pi^{-}) ;M_{p#pi} (GeV/c^{2});Counts", 400, 1.05, 1.25);
  TH1D* hDca12   = new TH1D("hDca12","DCA(p,#pi) (cm);DCA;Counts", 200, 0, 5);
  TH1D* hCosPA   = new TH1D("hCosPA","cos(pointing);cos;Counts", 200, 0.90, 1.0);
  TH1D* hL       = new TH1D("hL","decay length (cm);L;Counts", 300, 0, 30);
  TH1D* hNTrk    = new TH1D("hNTrk","N tracks/event;N;Events", 500, 0, 5000);

  // Chain + Pico maker
  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, InputFileList, "picoDst");

  if (chain->Init() == kStErr) {
    cout << "ERROR: chain->Init() returned kStErr" << endl;
    delete picoMaker;
    delete chain;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
  cout << "Total entries = " << total << endl;
  if (nEvents > total) nEvents = total;

  // constants
  const double mp  = 0.938272;
  const double mpi = 0.13957;

  // Loop
  for (Int_t i = 0; i < nEvents; i++) {

    //    if (i % 1000 == 0) cout << "Working on event " << i << endl;
    if (i % 10 == 0) cout << "Working on event " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);
    if (iret) { cout << "Bad return code! " << iret << endl; break; }

    StPicoDst* picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    StPicoEvent* evt = picoDst->event();
    if (!evt) continue;

    TVector3 pv = evt->primaryVertex();
    float B = evt->bField(); // kG (0でも動く)

    unsigned int nTrk = picoDst->numberOfTracks();
    hNTrk->Fill((double)nTrk);

    // track loop (p)
    for (unsigned int ip = 0; ip < nTrk; ip++) {
      StPicoTrack* tp = picoDst->track(ip);
      if (!tp) continue;

      if (tp->charge() <= 0) continue;       // p: +
      if (tp->nHitsFit() < 15) continue;
      if (tp->pMom().Mag() < 0.2) continue;

      // PID: proton
      if (std::fabs(tp->nSigmaProton()) > 3.0) continue;

      TVector3 rp = tp->origin();            // approx point
      TVector3 pp = tp->gMom(pv, B);         // approx momentum
      if (pp.Mag() < 1e-6) continue;
      TVector3 up = pp.Unit();

      // track loop (pi-)
      for (unsigned int iPi = 0; iPi < nTrk; iPi++) {
        if (iPi == ip) continue;

        StPicoTrack* tpi = picoDst->track(iPi);
        if (!tpi) continue;

        if (tpi->charge() >= 0) continue;    // pi-: -
        if (tpi->nHitsFit() < 15) continue;
        if (tpi->pMom().Mag() < 0.2) continue;

        // PID: pion
        if (std::fabs(tpi->nSigmaPion()) > 3.0) continue;

        TVector3 rpi = tpi->origin();
        TVector3 ppi = tpi->gMom(pv, B);
        if (ppi.Mag() < 1e-6) continue;
        TVector3 upi = ppi.Unit();

        // V0 vertex by 2-line PCA
        TVector3 pca1, pca2;
        double dca12 = 999;
        if (!LineLinePCA(rp, up, rpi, upi, pca1, pca2, dca12)) continue;

        TVector3 v0 = 0.5 * (pca1 + pca2);

        // topology
        double L = (v0 - pv).Mag();

        TVector3 pLam = pp + ppi;
        double cosPA = -2;
        if (pLam.Mag() > 1e-6 && L > 1e-6) {
          cosPA = pLam.Dot(v0 - pv) / (pLam.Mag() * L);
        }

        // very loose cuts (first look)
        if (dca12 > 1.0) continue;
        if (L < 2.0) continue;
        if (cosPA < 0.995) continue;

        // invariant mass
        double Ep  = std::sqrt(mp*mp   + pp.Mag2());
        double Epi = std::sqrt(mpi*mpi + ppi.Mag2());
        double m2  = (Ep + Epi)*(Ep + Epi) - (pp + ppi).Mag2();
        if (m2 < 0) continue;
        double m = std::sqrt(m2);

        // fill
        hMassLam->Fill(m);
        hDca12->Fill(dca12);
        hCosPA->Fill(cosPA);
        hL->Fill(L);
}
}
}

  // Finish chain
  chain->Finish();

  // Save histograms
  TString outName = OutputDir + "/lambdaQA.root";
  TFile* fout = new TFile(outName, "RECREATE");
  hMassLam->Write();
  hDca12->Write();
  hCosPA->Write();
  hL->Write();
  hNTrk->Write();
  fout->Close();

  cout << "==========================================" << endl;
  cout << "Done. Saved: " << outName << endl;
  cout << "==========================================" << endl;

  delete picoMaker;
  delete chain;
  chain = 0;
}
