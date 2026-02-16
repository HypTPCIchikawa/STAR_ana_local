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
// TOF path length calculation (STAR standard)
//============================================================
double tofPathLength(const StThreeVectorF* vtx,
                     const StThreeVectorF* hit,
                     double curvature)
{
  StThreeVectorF d = *hit - *vtx;
  double L = d.mag();

  if (curvature != 0)
    L = 2.0 * TMath::ASin(0.5 * d.mag() * curvature) / curvature;

  return fabs(L);
}

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
			400,0,4,400,-0.5,6.0);

  //============================
  // Event loop
  //============================
  for (int i=0;i<nEvents;i++) {

    //if (i%1000==0) cout<<"Event "<<i<<endl;
    cout<<"Event "<<i<<endl;

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
      int idx = trk->bTofPidTraitsIndex();
      if (idx < 0) continue;

      StPicoBTofPidTraits *tofPid =
        picoDst->btofPidTraits(idx);
      if (!tofPid) continue;

      if (tofPid->btof() <= 0) continue;

      // ---- vertex-based beta ----
      TVector3 tofPos = tofPid->btofHitPos();

      StThreeVectorF vtxSt(vtx.X(), vtx.Y(), vtx.Z());
      StThreeVectorF hitSt(tofPos.X(), tofPos.Y(), tofPos.Z());
      StThreeVectorF org(trk->origin().X(),
                         trk->origin().Y(),
                         trk->origin().Z());
      StThreeVectorF gmom(trk->gMom().X(),
                          trk->gMom().Y(),
                          trk->gMom().Z());

      StPhysicalHelixD helix(gmom, org,
                             BField*kilogauss,
                             trk->charge());

      double L = tofPathLength(&vtxSt, &hitSt,
                               helix.curvature());

      double tof = tofPid->btof(); // ns
      double beta = L / (tof * (C_C_LIGHT/1.e9));
      if (beta<=0 || beta>=1) continue;

      double m2 = p*p*(1.0/(beta*beta) - 1.0);
      hM2P->Fill(p,m2);
    }
  }

  //============================
  // Save
  //============================
  TFile *fout =
    new TFile(outDir+"/DedxP_M2.root","RECREATE");
  hDedxP->Write();
  hM2P->Write();
  fout->Close();

  chain->Finish();
  delete picoMaker;
  delete chain;

  cout<<"=== Finished successfully ==="<<endl;
}
