// PicoQuickCheck2.C  (ROOT5/CINT macro: no STAR headers)

#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <cstdio>

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

StChain *chain = 0;

void PicoQuickCheck2(Int_t nEvents=1, const char* inList="test2.list")
{
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inList, "picoDst");

  chain->Init();

  for (Int_t i=0; i<nEvents; i++) {
    Int_t ierr = chain->Make(i);
    if (ierr) break;

    StPicoDst* dst = picoMaker->picoDst();
    StPicoEvent* ev = dst ? dst->event() : 0;
    if (!ev) continue;

    printf("RunId=%d EventId=%d PVz=%g nTracks=%d\n",
           ev->runId(), ev->eventId(), ev->primaryVertex().Z(),
           dst->numberOfTracks());

    if (dst->numberOfTracks() > 0) {
      StPicoTrack* tr = dst->track(0);
      if (tr) {
        TVector3 p = tr->gMom();   // ←あなたの環境ではこれが使える
        printf("track0 id=%d q=%d p=(%g,%g,%g) nSigPi=%g nSigP=%g\n",
               tr->id(), tr->charge(), p.X(), p.Y(), p.Z(),
               tr->nSigmaPion(), tr->nSigmaProton());
      }
    }
  }

  chain->Finish();
}
