//============================================================
// p vs dE/dx  and  p vs m^2 (TOF, vertex-based beta)
//============================================================

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;
class StPhysicalHelixD;
class tofPathLength;
class StPicoBTofPidTraits;

#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TFile.h"
#include "TMath.h"
#include <TSystem.h>

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"

#include <iostream>
using namespace std;

StChain *chain = 0;


//============================================================
void Analysis_DedxVsP_M2_FixedTarget(Int_t nEvents = 0,
				     const char* inList = "test.list",
				     TString outDir = "./")
{
  if (nEvents == 0) nEvents = 1e9;

  //============================
  // Load libraries
  //============================
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  //============================
  // Chain
  //============================
  chain = new StChain();

  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "PicoDst");

  picoMaker->SetStatus("*",0);
  picoMaker->SetStatus("Event",1);
  picoMaker->SetStatus("Track",1);
  picoMaker->SetStatus("BTofPidTraits",1);

  if (chain->Init() == kStErr) {
    cout << "chain->Init() error!" << endl;
    return;
  }

  int totalEntries = picoMaker->chain()->GetEntries();
  cout << "Total entries = " << totalEntries << endl;
  if (nEvents > totalEntries) nEvents = totalEntries;

  //============================
  // Histograms
  //============================
  TH2D *hDedxP = new TH2D(
			  "hDedxP","TPC dE/dx vs p;p (GeV/c);dE/dx",
			  400,0,4,400,0,20);

  TH2D *hM2P = new TH2D(
			"hM2P","TOF m^{2} vs p;p (GeV/c);m^{2} (GeV^{2}/c^{4})",
			//400,0,4,400,-0.5,.0);
			400,0,4,1000,-0.5,25.0);
  //============================
  // Event loop
  //============================
  for (int i=0;i<nEvents;i++) {

    if (i%1000==0) cout<<"Event "<<i<<endl;

    chain->Clear();
    if (chain->Make(i)) break;

    StPicoDst *picoDst = picoMaker->picoDst();
    if (!picoDst) continue;

    StPicoEvent *event = picoDst->event();
    if (!event) continue;

    TVector3 vtx = event->primaryVertex();
    float BField = event->bField();

    for (UInt_t it=0; it<picoDst->numberOfTracks(); it++) {

      StPicoTrack *trk = picoDst->track(it);
      if (!trk) continue;

      // ---- quality cuts ----
      if (trk->nHitsDedx() < 15) continue;
      if (trk->gPt() < 0.1) continue;

      TVector3 pMom = trk->pMom();
      double p = pMom.Mag();

      // ---- dE/dx ----
      double dedx = trk->dEdx();
      if (dedx>0) hDedxP->Fill(p,dedx);

      // ---- TOF ----
      double beta = -1.0;
      int idx = trk->bTofPidTraitsIndex();
      if (idx < 0) continue;

      StPicoBTofPidTraits *tofPid =
        picoDst->btofPidTraits(idx);
      if (!tofPid) continue;

      beta = tofPid->btofBeta();
      // if (beta <= 1e-4) continue;

      // ---- vertex-based beta ----
      TVector3 tofPos = tofPid->btofHitPos();
      
      StThreeVectorF vtxSt(vtx.X(), vtx.Y(), vtx.Z());
      StThreeVectorF const hitSt(tofPos.X(), tofPos.Y(), tofPos.Z());
      StThreeVectorF org(trk->origin().X(),
                         trk->origin().Y(),
                         trk->origin().Z());
      StThreeVectorF gmom(trk->gMom().X(),
                          trk->gMom().Y(),
                          trk->gMom().Z());

      StPhysicalHelixD helix(gmom, org,
                             BField*kilogauss,
			     static_cast<float>(trk->charge()) );
      //                             trk->charge());
      
      //      cout<<"vtx:"<<vtx<<", tofPos:"<<tofPos<<", roughpath="<<(vtx-tofPos).Mag()<<endl;
      
      double L   = tofPathLength(&vtxSt, &hitSt, helix.curvature());
      double roughpath   = (vtx-tofPos).Mag();
      double tof = tofPid->btof();
      //cout<<"tof="<<tof<<", L="<<L<<endl;
      //if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
      if (tof > 0) {
	beta = roughpath / (tof * (C_C_LIGHT / 1.e9));
	cout<<"beta="<<beta<<endl;
      }
      else         beta = -1.0; //std::numeric_limits<float>::quiet_NaN();
      //cout<<"beta="<<beta<<endl;
      if (beta<=0) continue;
      double m2 = p*p*(1.0/(beta*beta) - 1.0);
      
      hM2P->Fill(p,m2);
    }
  }

  //============================
  // Save
  //============================
  TFile *fout =
    new TFile(outDir+"/DedxP_M2_mod.root","RECREATE");
  hDedxP->Write();
  hM2P->Write();
  fout->Close();

  chain->Finish();
  delete picoMaker;
  delete chain;

  cout<<"=== Finished successfully ==="<<endl;
}
