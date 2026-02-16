#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

//======================================================
// Helix を「安全に」作る（CINT対応）
//======================================================
StPhysicalHelixD makeHelix(const StPicoTrack* trk, double bField)
{
  StThreeVectorF p(trk->gMom().x(),
                   trk->gMom().y(),
                   trk->gMom().z());

  StThreeVectorF o(trk->origin().x(),
                   trk->origin().y(),
                   trk->origin().z());

  return StPhysicalHelixD(p, o,
                          bField * kilogauss,
                          trk->charge());
}

//======================================================
// Helix を使った「安全な」V0 位置推定
//======================================================
bool MakeV0PairHelix(const StPicoTrack* p,
                     const StPicoTrack* pi,
                     double bField,
                     TVector3& v0,
                     double& dca12)
{
  StPhysicalHelixD hp  = makeHelix(p,  bField);
  StPhysicalHelixD hpi = makeHelix(pi, bField);

  double s1 = 0, s2 = 0;
  dca12 = hp.pathLengths(hpi, s1, s2);

  if (dca12 < 0 || dca12 > 1.0) return false;

  StThreeVectorD x1 = hp.at(s1);
  StThreeVectorD x2 = hpi.at(s2);

  v0.SetXYZ(0.5*(x1.x()+x2.x()),
            0.5*(x1.y()+x2.y()),
            0.5*(x1.z()+x2.z()));

  return true;
}


//======================================================
// Main analysis
//======================================================
void LambdaInvMass_Helix(Int_t nEvents=10,
                         const char* inList="test2.list")
{
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  std::cout<<"start"<<std::endl;
  if (nEvents == 0) nEvents = 100000000;

  // ---- load STAR libraries (CINT style) ----
  gSystem->Load("StRefMultCorr");

  // ---- chain & maker ----
  chain = new StChain();
  StPicoDstMaker* picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inList, "picoDst");

  chain->Init();
  StPicoDst* picoDst = picoMaker->picoDst();

  int nEntries = picoMaker->chain()->GetEntries();
  if (nEvents > nEntries) nEvents = nEntries;

  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass;M_{p#pi^{-}} (GeV/c^{2});Counts",
			200, 1.05, 1.25);

  const double mp = 0.938272;
  const double mpi = 0.139570;

  for (int ie=0; ie<nEvents; ie++) {
    std::cout<<"event:"<<ie<<std::endl;

    chain->Clear();
    if (chain->Make(ie)) break;

    StPicoEvent* evt = picoDst->event();
    TVector3 pv = evt->primaryVertex();
    if (fabs(pv.z()) > 200) continue;  // safety



    double bField = evt->bField();
    int nTr = picoDst->numberOfTracks();
   
    for (int ip = 0; ip < nTr; ip++) {
      StPicoTrack* p = picoDst->track(ip);
      if (!p) continue;
      if (p->charge() <= 0) continue;
      if (fabs(p->nSigmaProton()) > 2) continue;
      if (p->gDCA(pv.x(), pv.y(), pv.z()) < 0.5) continue;

      for (int ii = 0; ii < nTr; ii++) {
        if (ii == ip) continue;

        StPicoTrack* pi = picoDst->track(ii);
        if (!pi) continue;

        if (pi->charge() >= 0) continue;
        if (fabs(pi->nSigmaPion()) > 2) continue;
        if (pi->gDCA(pv.x(), pv.y(), pv.z()) < 0.8) continue;

        TVector3 v0;
	double dca;
        if (!MakeV0PairHelix(p, pi, bField, v0, dca)) continue;

        TVector3 momP  = p->gMom();
        TVector3 momPi = pi->gMom();

        TLorentzVector lp, lpi;
        lp.SetVectM(momP, mp);
        lpi.SetVectM(momPi, mpi);

        hLam->Fill((lp + lpi).M());
      }
    }
  }

  // ---- save ----
  TFile* fout =
    new TFile("LambdaInvMass.root", "RECREATE");
  hLam->Write();
  fout->Close();
}
