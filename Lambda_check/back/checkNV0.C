#include <iostream>
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoDstReader.h"
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoDst.h"
#include "/afs/rhic.bnl.gov/star/packages/SL24a/StRoot/StPicoEvent/StPicoEvent.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
//#include "StPicoV0.h"


void checkNV0(const char* inFile =
	      "/star/data107/reco/production_7p3GeV_fixedTarget_2020/ReversedFullField/P24ia/2020/035/21035003/st_physics_21035003_raw_1000019.picoDst.root",
	      Long64_t nMax = 2000)
{
  // 依存ロード（あなたの手順通り）
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("libStPicoEvent");

  StPicoDstReader* r = new StPicoDstReader(inFile);
  r->Init();

  TChain* ch = r->chain();
  if (!ch) { std::cout << "ERROR: chain() is null\n"; return; }

  Long64_t nEntries = ch->GetEntries();
  Long64_t nToRun = nEntries;
  if (nMax > 0 && nToRun > nMax) nToRun = nMax;

  long long sumV0 = 0;
  Long64_t nEvt = 0;

  for (Long64_t i = 0; i < nToRun; i++) {
    if (!r->readPicoEvent(i)) continue;   // ★引数付き

    StPicoDst* dst = r->picoDst();
    if (!dst) continue;

    sumV0 += (long long)dst->numberOfV0s();
    nEvt++;
  }

  std::cout << "entries=" << nEntries
            << " processed=" << nEvt
            << " <N_V0>=" << (nEvt ? (double)sumV0/nEvt : 0.0)
            << std::endl;
}
