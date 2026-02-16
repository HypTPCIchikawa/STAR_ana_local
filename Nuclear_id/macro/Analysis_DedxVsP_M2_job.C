//============================================================
// STAR picoDst analysis macro
//  p vs dE/dx correlation
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

using namespace std;

StChain *chain = 0;

void Analysis_DedxVsP_M2_job(Int_t nEvents,
			     const char* inList,
			     const char* jobid)
{

  double d_mean= 3.48096e+00;
  double d_sigma= 1.41458e-01;

  double He3_mean= 1.92385e+00;
  double He3_sigma= 0.94e-01;

  double t_mean= 7.76906e+00;
  double t_sigma= 3.41755e-01;

  if (nEvents == 0) nEvents = 1000000000;

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

      if(m2>d_mean-1.5*d_sigma&&m2<d_mean+1.5*d_sigma)
	hDedxP_d->Fill(p, dedx);

      if(m2>t_mean-1.5*d_sigma&&m2<t_mean+1.5*d_sigma)
	hDedxP_t->Fill(p, dedx);

      if(m2>He3_mean-1.5*He3_sigma&&m2<He3_mean+1.5*He3_sigma)
	hDedxP_3He->Fill(p, dedx);
      
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
