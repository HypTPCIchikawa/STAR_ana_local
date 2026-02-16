#include <TSystem>
#include <TFile>
#include <TH1D>
#include <TVector3>
#include <TLorentzVector>
#include <TMath>
#include <iostream>

// ===== forward declarations (CINT style) =====
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

//==================================================
bool MakeV0PairLinear(const StPicoTrack* p,
                      const StPicoTrack* pi,
                      const TVector3& pv,
                      TVector3& v0,
                      double& dca12)
{
  TVector3 op = p->origin();
  TVector3 oi = pi->origin();

  TVector3 dp = p->gMom(pv).Unit();
  TVector3 di = pi->gMom(pv).Unit();

  TVector3 w0 = op - oi;
  double a = dp.Dot(dp);
  double b = dp.Dot(di);
  double c = di.Dot(di);
  double d = dp.Dot(w0);
  double e = di.Dot(w0);

  double denom = a*c - b*b;
  if (fabs(denom) < 1e-5) return false;

  double s = (b*e - c*d)/denom;
  double t = (a*e - b*d)/denom;

  TVector3 xp = op + s*dp;
  TVector3 xi = oi + t*di;

  dca12 = (xp - xi).Mag();
  if (dca12 > 1.0) return false;

  v0 = 0.5*(xp + xi);
  return true;
}

//==================================================
void Analysis(Int_t nEvents,
              TString InputFileList,
              TString OutputDir="./")
{
  if (nEvents == 0) nEvents = 100000000;

  // ===== STAR libraries =====
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  // ===== chain =====
  chain = new StChain();
  StPicoDstMaker* picoMaker =
    new StPicoDstMaker(2, InputFileList, "picoDst");

  if (chain->Init() == kStErr) {
    std::cout << "chain->Init failed" << std::endl;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
  if (nEvents > total) nEvents = total;
  std::cout << "Total entries = " << total << std::endl;

  // ===== histogram =====
  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass;M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp  = 0.938272;
  const double mpi = 0.139570;

  // ===== event loop =====
  for (int i=0;i<nEvents;i++) {
    chain->Clear();
    if (chain->Make(i)) break;

    StPicoDst* picoDst = picoMaker->picoDst();
    StPicoEvent* evt   = picoDst->event();
    if (!evt) continue;

    TVector3 pv = evt->primaryVertex();
    double pvx = pv.x(), pvy = pv.y(), pvz = pv.z();

    int nTr = picoDst->numberOfTracks();

    for (int ip=0; ip<nTr; ip++) {
      StPicoTrack* p = picoDst->track(ip);
      if (!p) continue;
      if (p->charge() <= 0) continue;
      if (fabs(p->nSigmaProton()) > 2.0) continue;

      for (int ii=0; ii<nTr; ii++) {
        if (ii == ip) continue;
        StPicoTrack* pi = picoDst->track(ii);
        if (!pi) continue;
        if (pi->charge() >= 0) continue;
        if (fabs(pi->nSigmaPion()) > 2.0) continue;
        if (pi->gDCA(pvx,pvy,pvz) < 0.8) continue;

        TVector3 v0;
        double dca12;
        if (!MakeV0PairLinear(p, pi, pv, v0, dca12)) continue;

        TVector3 momP  = p->gMom(pv);
        TVector3 momPi = pi->gMom(pv);

        TLorentzVector lp, lpi;
        lp.SetVectM(momP, mp);
        lpi.SetVectM(momPi, mpi);

        hLam->Fill((lp+lpi).M());
      }
    }
  }

  // ===== output =====
  TFile* fout = new TFile(OutputDir+"/LambdaInvMass.root","RECREATE");
  hLam->Write();
  fout->Close();
}
