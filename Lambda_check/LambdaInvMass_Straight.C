// LambdaInvMass_Straight.C
// ROOT5 / CINT macro (NO helix, NO StPhysicalHelix)

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

void LambdaInvMass_Straight(Int_t nEvents=100,
                            const char* inList="test2.list")
{
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  chain = new StChain();
  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(2, inList, "picoDst");

  chain->Init();
  StPicoDst *picoDst = picoMaker->picoDst();

  TH1D *hM = new TH1D("hM","Lambda invariant mass;M(p#pi^{-}) [GeV/c^{2}]",200,1.05,1.25);

  int nEntries = picoMaker->chain()->GetEntries();
  if (nEvents > nEntries) nEvents = nEntries;

  for (int ie=0; ie<nEvents; ie++) {

    std::cout<<"event:"<<ie<<std::endl;
    chain->Clear();
    chain->Make();

    StPicoEvent *ev = picoDst->event();
    TVector3 PV(ev->primaryVertex().x(),
                ev->primaryVertex().y(),
                ev->primaryVertex().z());

    int nTr = picoDst->numberOfTracks();

    for (int i=0;i<nTr;i++){
      StPicoTrack *p = picoDst->track(i);
      if (!p) continue;

      if (p->charge()<=0) continue;
      if (fabs(p->nSigmaProton())>2.0) continue;
      if (p->gDCA(PV.x(),PV.y(),PV.z()) < 1.0) continue;

      TVector3 pMom = p->gMom();

      for (int j=0;j<nTr;j++){
        StPicoTrack *pi = picoDst->track(j);
        if (!pi) continue;

        if (pi->charge()>=0) continue;
        if (fabs(pi->nSigmaPion())>2.0) continue;
        if (pi->gDCA(PV.x(),PV.y(),PV.z()) < 1.0) continue;

        TVector3 piMom = pi->gMom();

        TLorentzVector lp, lpi;
        lp.SetVectM(pMom, 0.938272);
        lpi.SetVectM(piMom, 0.139570);

        double mass = (lp+lpi).M();
        hM->Fill(mass);
      }
    }
  }

  TCanvas *c = new TCanvas("c","Lambda",800,600);
  hM->Draw();
}
