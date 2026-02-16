#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>

// 絶対パスincludeでもOK（あなたの環境に合わせた）
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoDstReader.h"
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoDst.h"
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoEvent.h"
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoV0.h"   // ★これが必要

void testLambdaFromPico(Int_t nMaxEvents = 2000)
{
  // ★PicoDstReader用：これだけでOK（StPicoDstMakerは不要）
  int rc = gSystem->Load("libStPicoEvent");
  if (rc < 0) {
    std::cout << "ERROR: cannot load libStPicoEvent\n";
    return;
  }

  const char* inFile =
    "/star/data107/reco/production_7p3GeV_fixedTarget_2020/"
    "ReversedFullField/P24ia/2020/035/21035003/"
    "st_physics_21035003_raw_1000019.picoDst.root";

  StPicoDstReader* reader = new StPicoDstReader(inFile);
  reader->Init();

  if (!reader->chain()) {
    std::cout << "ERROR: Cannot initialize PicoDstReader\n";
    return;
  }

  TH1D* hMassLam = new TH1D("hMassLam", "#Lambda mass;M_{p#pi} (GeV/c^{2});Counts", 200, 1.08, 1.15);
  TH1D* hNV0     = new TH1D("hNV0",     "N_{V0} per event;N_{V0};Events",           50,  0,   50);

  int nEvents = 0;

  while (reader->Make()) {
    if (nMaxEvents > 0 && nEvents >= nMaxEvents) break;

    StPicoDst* picoDst = reader->picoDst();
    if (!picoDst) continue;

    StPicoEvent* event = picoDst->event();
    if (!event) continue;

    nEvents++;

    unsigned int nV0 = picoDst->numberOfV0s();
    hNV0->Fill(nV0);

    for (unsigned int i=0; i<nV0; i++) {
      StPicoV0* v0 = picoDst->v0(i);
      if (!v0) continue;

      // loose topo cuts for existence check
      if (v0->decayLength() < 1.0) continue;
      if (v0->cosTheta() < 0.99) continue;

      hMassLam->Fill(v0->massLambda());
    }
  }

  std::cout << "Processed events = " << nEvents << std::endl;

  TCanvas* c1 = new TCanvas("c1","Lambda test",1200,500);
  c1->Divide(2,1);
  c1->cd(1); hNV0->Draw();
  c1->cd(2); hMassLam->Draw();
}
