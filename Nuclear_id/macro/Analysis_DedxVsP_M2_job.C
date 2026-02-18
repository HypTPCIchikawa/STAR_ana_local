//============================================================
// STAR picoDst analysis macro
//  p vs dE/dx correlation
//  (ACLiC なし .L のみで実行可能：NuclearPID をマクロ内に定義)
//============================================================

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include <iostream>
#include <cmath>

using namespace std;

//---- NuclearPID: マクロ内定義（NuclearPID.h と同期、コンパイルなし実行用） ----
struct NuclearPID {
  enum ParticleType { kDeuteron = 0, kTriton, kHe3, kHe4 };
  double fMeanParams[4][4];   // [type][p0..p3]
  double fSigmaParams[4][4];
  double fRange[4][2];       // [type][0]=pMin, [1]=pMax

  NuclearPID() {
    // Mean (p0,p1,p2,p3): deuteron, triton, He3, He4
    double mean[4][4] = {
      {10.2875, -36.6639, 67.63, -45.2483},
      {19.2177, -70.3198, 113.876, -67.0591},
      {28.8747, -68.8637, 104.455, -56.9342},
      {39.974, -103.365, 139.182, -67.1446}
    };
    // Sigma (p0,p1,p2,p3)
    double sigma[4][4] = {
      {0.771189, -2.63609, 2.8044, 0.0},
      {1.40996, -4.97147, 4.97608, 0.0},
      {1.636, -4.35588, 4.40594, 0.0},
      {2.24127, -7.08328, 8.07124, 0.0}
    };
    // Momentum range [pMin, pMax]
    double range[4][2] = {
      {0.2, 6.0},
      {0.2, 6.0},
      {0.4, 5.0},
      {0.4, 5.0}
    };
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        fMeanParams[i][j] = mean[i][j];
        fSigmaParams[i][j] = sigma[i][j];
      }
      fRange[i][0] = range[i][0];
      fRange[i][1] = range[i][1];
    }
  }

  bool IsValid(ParticleType type, double p) const {
    int t = (int)type;
    return (p >= fRange[t][0] && p <= fRange[t][1]);
  }

  double LogPoly(double p, const double* par) const {
    if (p <= 0) return 0;
    double lp = log10(p);
    return par[0] + par[1]*lp + par[2]*lp*lp + par[3]*lp*lp*lp;
  }

  double GetMean(ParticleType type, double p) const {
    if (!IsValid(type, p)) return 0;
    return LogPoly(p, fMeanParams[(int)type]);
  }

  double GetSigma(ParticleType type, double p) const {
    if (!IsValid(type, p)) return 0;
    return LogPoly(p, fSigmaParams[(int)type]);
  }

  double GetNSigma(ParticleType type, double p, double dedx) const {
    if (!IsValid(type, p)) return 999.0;
    double mean = GetMean(type, p);
    double sigma = GetSigma(type, p);
    if (sigma <= 0) return 999.0;
    return (dedx - mean) / sigma;
  }
};

StChain *chain = 0;

void Analysis_DedxVsP_M2_job(Int_t nEvents,
			     const char* inList,
			     const char* jobid)
{

  // M2 cut parameters (TOF mass squared)
  double d_mean   = 3.48096e+00;
  double d_sigma  = 1.41458e-01;
  double He3_mean = 1.92385e+00;
  double He3_sigma= 0.94e-01;
  double t_mean   = 7.76906e+00;
  double t_sigma  = 3.41755e-01;

  if (nEvents == 0) nEvents = 1000000000;

  // Nuclear PID (dE/dx nSigma for deuteron, triton, 3He, 4He)
  NuclearPID pid;

  //============================================================
  // Load STAR libraries
  //============================================================
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");


  cout << "### MACRO START ###" << endl;
  cout << "STAR = " << gSystem->Getenv("STAR") << endl;
  cout << "LD_LIBRARY_PATH = " << gSystem->Getenv("LD_LIBRARY_PATH") << endl;
  cout << "FILELIST = " << inList << endl;

  cout << "inList = " << inList << endl;

  if (gSystem->AccessPathName(inList)) {
    cout << "ERROR: FILELIST does not exist!" << endl;
    gSystem->Exit(1);
  }

  //============================================================
  // Prepare output directory
  //============================================================
  const char* scratch = gSystem->Getenv("SCRATCH");
  gSystem->mkdir(Form("%s/Nuclear_id", scratch), kTRUE);
  cout << "SCRATCH = " << scratch << endl;
  gSystem->Exec("ls -l $SCRATCH");
  gSystem->Exec("ls -l $SCRATCH/Nuclear_id");


  //============================================================
  // Chain & PicoDstMaker
  //============================================================
  chain = new StChain();

  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  picoMaker->SetStatus("*",0);
  picoMaker->SetStatus("Event",1);
  picoMaker->SetStatus("Track",1);
  picoMaker->SetStatus("BTofHit",1);
  picoMaker->SetStatus("BTofPidTraits",1);
  picoMaker->SetStatus("BbcHit",1);
  picoMaker->SetStatus("EpdHit",1);
  picoMaker->SetStatus("MtdHit",1);
  picoMaker->SetStatus("BTowHit",1);
  picoMaker->SetStatus("ETofPidTraits",1);
  
  cout << "### BEFORE Init()" << endl;
  if (chain->Init() == kStErr) {
    cout << "chain->Init() error!" << endl;
    //    return;
  }
  cout << "### AFTER Init(): "<< endl;
  
  cout << "### INPUT FILE CONTENT ###" << endl;
  picoMaker->chain()->GetFile()->ls();

  Int_t totalEntries = picoMaker->chain()->GetEntries();
  cout << "Total entries = " << totalEntries << endl;
  if (nEvents > totalEntries) nEvents = totalEntries;

  //============================================================
  // Histogram
  //============================================================
  TH2D *hDedxP = new TH2D(
			  "hDedxP",
			  "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			  600, 0.0, 6.0,
			  1200, 0.0, 120.0);

  TH2D *hDedxP_cut = new TH2D(
			      "hDedxP_cut",
			      "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			      600, 0.0, 6.0,
			      1200, 0.0, 120.0);

  TH2D *hDedxP_e = new TH2D(
			    "hDedxP_e",
			    "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			    600, 0.0, 6.0,
			    1200, 0.0, 120.0);

  TH2D *hDedxP_pi = new TH2D(
			     "hDedxP_pi",
			     "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			     600, 0.0, 6.0,
			     1200, 0.0, 120.0);
  
  TH2D *hDedxP_K = new TH2D(
			    "hDedxP_K",
			    "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			    600, 0.0, 6.0,
			    1200, 0.0, 120.0);

  TH2D *hDedxP_p = new TH2D(
			    "hDedxP_p",
			    "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			    600, 0.0, 6.0,
			    1200, 0.0, 120.0);

  TH2D *hDedxP_else = new TH2D(
			       "hDedxP_else",
			       "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			       600, 0.0, 6.0,
			       1200, 0.0, 120.0);

  TH2D *hDedxP_d = new TH2D(
			   "hDedxP_d",
			   "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			   600, 0.0, 6.0,
			   1200, 0.0, 120.0);

  TH2D *hDedxP_t = new TH2D(
			    "hDedxP_t",
			    "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			    600, 0.0, 6.0,
			    1200, 0.0, 120.0);
  
  TH2D *hDedxP_3He = new TH2D(
			      "hDedxP_3He",
			      "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			      600, 0.0, 6.0,
			      1200, 0.0, 120.0);

  TH2D *hDedxP_4He = new TH2D(
			      "hDedxP_4He",
			      "TPC dE/dx vs p; p (GeV/c); dE/dx (keV/cm)",
			      600, 0.0, 6.0,
			      1200, 0.0, 120.0);

  // 2sigma-selected: p vs PseudoRapidity, p vs Rapidity, p vs pT (d, t, 3He, 4He)
  TH2D *hPvsEta_d   = new TH2D("hPvsEta_d",   "p vs #eta (2#sigma); p (GeV/c); #eta",   600, 0., 6., 100, -3., 3.);
  TH2D *hPvsEta_t   = new TH2D("hPvsEta_t",   "p vs #eta (2#sigma); p (GeV/c); #eta",   600, 0., 6., 100, -3., 3.);
  TH2D *hPvsEta_3He = new TH2D("hPvsEta_3He", "p vs #eta (2#sigma); p (GeV/c); #eta",   600, 0., 6., 100, -3., 3.);
  TH2D *hPvsEta_4He = new TH2D("hPvsEta_4He", "p vs #eta (2#sigma); p (GeV/c); #eta",   600, 0., 6., 100, -3., 3.);
  TH2D *hPvsY_d    = new TH2D("hPvsY_d",    "p vs y (2#sigma); p (GeV/c); y",           600, 0., 6., 100, -3., 3.);
  TH2D *hPvsY_t    = new TH2D("hPvsY_t",    "p vs y (2#sigma); p (GeV/c); y",           600, 0., 6., 100, -3., 3.);
  TH2D *hPvsY_3He  = new TH2D("hPvsY_3He",  "p vs y (2#sigma); p (GeV/c); y",           600, 0., 6., 100, -3., 3.);
  TH2D *hPvsY_4He  = new TH2D("hPvsY_4He",  "p vs y (2#sigma); p (GeV/c); y",           600, 0., 6., 100, -3., 3.);
  TH2D *hPvsPt_d   = new TH2D("hPvsPt_d",   "p vs p_{T} (2#sigma); p (GeV/c); p_{T} (GeV/c)", 600, 0., 6., 600, 0., 6.);
  TH2D *hPvsPt_t   = new TH2D("hPvsPt_t",   "p vs p_{T} (2#sigma); p (GeV/c); p_{T} (GeV/c)", 600, 0., 6., 600, 0., 6.);
  TH2D *hPvsPt_3He = new TH2D("hPvsPt_3He", "p vs p_{T} (2#sigma); p (GeV/c); p_{T} (GeV/c)", 600, 0., 6., 600, 0., 6.);
  TH2D *hPvsPt_4He = new TH2D("hPvsPt_4He", "p vs p_{T} (2#sigma); p (GeV/c); p_{T} (GeV/c)", 600, 0., 6., 600, 0., 6.);

  // M2-selected (TOF mass squared cut)
  TH2D *hDedxP_d_m2 = new TH2D(
			       "hDedxP_d_m2",
			       "TPC dE/dx vs p (M2 cut); p (GeV/c); dE/dx (keV/cm)",
			       600, 0.0, 6.0,
			       1200, 0.0, 120.0);
  TH2D *hDedxP_t_m2 = new TH2D(
			       "hDedxP_t_m2",
			       "TPC dE/dx vs p (M2 cut); p (GeV/c); dE/dx (keV/cm)",
			       600, 0.0, 6.0,
			       1200, 0.0, 120.0);
  TH2D *hDedxP_3He_m2 = new TH2D(
				"hDedxP_3He_m2",
				"TPC dE/dx vs p (M2 cut); p (GeV/c); dE/dx (keV/cm)",
				600, 0.0, 6.0,
				1200, 0.0, 120.0);
  
  
  TH2D *hM2P = new TH2D(
			"hM2P","TOF m^{2} vs p;p (GeV/c);m^{2} (GeV^{2}/c^{4})",
			600,0,6,1000,-0.5,25.0);
  
  
  //============================================================
  // Event loop
  //============================================================
  for (Int_t i = 0; i < nEvents; i++) {

    if (i % 1000 == 0)
      cout << "Working on event " << i << endl;

    chain->Clear();
    Int_t iret = chain->Make(i);
    if (iret) {
      cout << "Bad return code: " << iret << endl;
      break;
    }

    StPicoDst *picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    StPicoEvent *event = picoDst->event();
    if (!event) continue;

    Int_t nTracks = picoDst->numberOfTracks();
    for (Int_t it = 0; it < nTracks; it++) {

      StPicoTrack *trk = picoDst->track(it);
      if (!trk) continue;

      // ---- basic quality cuts ----
      if (trk->nHitsDedx() < 15) continue;
      if (trk->gPt() < 0.1) continue;

      // momentum
      TVector3 pMom = trk->pMom();
      Double_t p = pMom.Mag();

      // dE/dx
      Double_t dedx = trk->dEdx();
      if (dedx <= 0) continue;

      hDedxP->Fill(p, dedx);
      
      if (fabs(trk->nSigmaElectron()) < 1.5) 
	hDedxP_e->Fill(p, dedx);
      
      if (fabs(trk->nSigmaPion()) < 1.5) 
	hDedxP_pi->Fill(p, dedx);
      
      if (fabs(trk->nSigmaKaon()) < 1.5) 
	hDedxP_K->Fill(p, dedx);

      if (fabs(trk->nSigmaProton()) < 1.5) 
	hDedxP_p->Fill(p, dedx);

      if (fabs(trk->nSigmaElectron()) >2 && 
	  fabs(trk->nSigmaPion()) >2 && 
	  fabs(trk->nSigmaKaon()) >2 && 
	  fabs(trk->nSigmaProton()) >2) 
	hDedxP_else->Fill(p, dedx);

      // dE/dx 2sigma selection (NuclearPID)
      double nSigma_d   = pid.GetNSigma(NuclearPID::kDeuteron, p, dedx);
      double nSigma_t   = pid.GetNSigma(NuclearPID::kTriton,  p, dedx);
      double nSigma_3He = pid.GetNSigma(NuclearPID::kHe3,     p, dedx);
      double nSigma_4He = pid.GetNSigma(NuclearPID::kHe4,     p, dedx);

      if (fabs(nSigma_d)   < 2.0) hDedxP_d->Fill(p, dedx);
      if (fabs(nSigma_t)   < 2.0) hDedxP_t->Fill(p, dedx);
      if (fabs(nSigma_3He) < 2.0) hDedxP_3He->Fill(p, dedx);
      if (fabs(nSigma_4He) < 2.0) hDedxP_4He->Fill(p, dedx);

      // p vs eta, p vs y, p vs pT (2sigma selection; masses in GeV/c^2)
      const double M_d   = 1.87561;
      const double M_t   = 2.80892;
      const double M_3He = 2.80839;
      const double M_4He = 3.72742;
      double eta = pMom.PseudoRapidity();
      double pT  = pMom.Perp();
      double pz  = pMom.Z();
      if (fabs(nSigma_d) < 2.0) {
        double E = TMath::Sqrt(M_d*M_d + p*p);
        double y = 0.5 * TMath::Log((E + pz) / (E - pz));
        hPvsEta_d->Fill(p, eta);
        hPvsY_d->Fill(p, y);
        hPvsPt_d->Fill(p, pT);
      }
      if (fabs(nSigma_t) < 2.0) {
        double E = TMath::Sqrt(M_t*M_t + p*p);
        double y = 0.5 * TMath::Log((E + pz) / (E - pz));
        hPvsEta_t->Fill(p, eta);
        hPvsY_t->Fill(p, y);
        hPvsPt_t->Fill(p, pT);
      }
      if (fabs(nSigma_3He) < 2.0) {
        double E = TMath::Sqrt(M_3He*M_3He + p*p);
        double y = 0.5 * TMath::Log((E + pz) / (E - pz));
        hPvsEta_3He->Fill(p, eta);
        hPvsY_3He->Fill(p, y);
        hPvsPt_3He->Fill(p, pT);
      }
      if (fabs(nSigma_4He) < 2.0) {
        double E = TMath::Sqrt(M_4He*M_4He + p*p);
        double y = 0.5 * TMath::Log((E + pz) / (E - pz));
        hPvsEta_4He->Fill(p, eta);
        hPvsY_4He->Fill(p, y);
        hPvsPt_4He->Fill(p, pT);
      }

      // ---- TOF m^2 ----
      if (!trk->isTofTrack()){
	hDedxP_cut->Fill(p, dedx);
	continue;
      }

      int idx = trk->bTofPidTraitsIndex();
      if (idx < 0) {
	hDedxP_cut->Fill(p, dedx);
	continue;
      }
      StPicoBTofPidTraits *btof =
        picoDst->btofPidTraits(idx);
      if (!btof) {
	hDedxP_cut->Fill(p, dedx);
	continue;
      }
      double beta = btof->btofBeta();
      if (beta<=0){
	hDedxP_cut->Fill(p, dedx);
	continue;
      }
      double m2 = p*p*(1.0/(beta*beta) - 1.0);
      hM2P->Fill(p,m2);


      // M2 selection (TOF mass squared)
      if (m2 > d_mean   - 1.5 * d_sigma   && m2 < d_mean   + 1.5 * d_sigma)   hDedxP_d_m2->Fill(p, dedx);
      if (m2 > t_mean   - 1.5 * t_sigma   && m2 < t_mean   + 1.5 * t_sigma)   hDedxP_t_m2->Fill(p, dedx);
      if (m2 > He3_mean - 1.5 * He3_sigma && m2 < He3_mean + 1.5 * He3_sigma) hDedxP_3He_m2->Fill(p, dedx);
      
      if(m2>0.01&&m2<10.&&p>0.2) continue;
      else hDedxP_cut->Fill(p, dedx);
      
    }
  }

  //============================================================
  // Draw & Save
  //============================================================
  TCanvas *c1 = new TCanvas("c1", "dE/dx vs p", 800, 700);
  c1->SetLogz();
  hDedxP->Draw("colz");

  TString outname = Form("%s/Nuclear_id/DedxVsP_M2_%s.root",
                         gSystem->Getenv("SCRATCH"), jobid);
  TFile *fout = new TFile(outname, "RECREATE");
  hDedxP->Write();
  hDedxP_cut->Write();
  hDedxP_e->Write();
  hDedxP_pi->Write();
  hDedxP_K->Write();
  hDedxP_p->Write();
  hDedxP_else->Write();
  hDedxP_d->Write();
  hDedxP_t->Write();
  hDedxP_3He->Write();
  hDedxP_4He->Write();
  hDedxP_d_m2->Write();
  hDedxP_t_m2->Write();
  hDedxP_3He_m2->Write();
  hPvsEta_d->Write();
  hPvsEta_t->Write();
  hPvsEta_3He->Write();
  hPvsEta_4He->Write();
  hPvsY_d->Write();
  hPvsY_t->Write();
  hPvsY_3He->Write();
  hPvsY_4He->Write();
  hPvsPt_d->Write();
  hPvsPt_t->Write();
  hPvsPt_3He->Write();
  hPvsPt_4He->Write();
  hM2P->Write();
  fout->Close();

  //============================================================
  // Finish
  //============================================================
  chain->Finish();

  delete picoMaker;
  delete chain;

  cout << "******************************************" << endl;
  cout << " Analysis finished successfully" << endl;
  cout << " Total processed events = " << nEvents << endl;
  cout << " Output file: " << outname << endl;
  cout << "******************************************" << endl;
}
