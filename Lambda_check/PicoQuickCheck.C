// PicoQuickCheck.C  (ROOT5/CINT macro)
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
using namespace std;

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;

StChain *chain = 0;

void PicoQuickCheck(Int_t nEvents=1, const char* inList="test2.list")
{
  // --- load STAR libs
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr"); // 無くても動く場合あり

  // --- chain & maker
  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2, inList, "picoDst");

  // --- init
  chain->Init();

  // --- event loop (min)
  for(int ie=0; ie<nEvents; ++ie){
    int iret = chain->Make(ie);
    if(iret){
      cout << "chain->Make returned " << iret << " at event " << ie << endl;
      break;
    }

    StPicoDst* dst = picoMaker->picoDst();
    if(!dst){
      cout << "picoDst is null" << endl;
      continue;
    }

    StPicoEvent* ev = dst->event();
    if(ev){
      cout << "RunId=" << ev->runId()
           << " EventId=" << ev->eventId()
           << " PVz=" << ev->primaryVertex().z()
           << " nTracks=" << dst->numberOfTracks()
           << endl;
    } else {
      cout << "event() is null" << endl;
    }

    if(dst->numberOfTracks()<=0){
      cout << "No tracks" << endl;
      continue;
    }

    StPicoTrack* trk = dst->track(0);
    if(!trk){
      cout << "track(0) is null" << endl;
      continue;
    }

    // --- try to print basic PID-ish quantities
    cout << "==== track(0) ====" << endl;
    cout << "id=" << trk->id() << endl;

    // charge(): 実装によっては sign() / charge() のどちらか
    // まずは charge() を試す
    cout << "charge=" << trk->charge() << endl;

    // nSigma (TPC PID)
    cout << "nSigmaPion="   << trk->nSigmaPion()
         << " nSigmaKaon="  << trk->nSigmaKaon()
         << " nSigmaProton="<< trk->nSigmaProton()
         << " nSigmaElectron="<< trk->nSigmaElectron()
         << endl;

    // dE/dx
    cout << "dEdx=" << trk->dEdx()
         << " nHitsFit=" << trk->nHitsFit()
         << " nHitsDedx=" << trk->nHitsDedx()
         << endl;

    // --- TOF beta: BTofPidTraits に紐づいていることが多い
    int ibtof = trk->bTofPidTraitsIndex();
    if(ibtof>=0){
      StPicoBTofPidTraits* btof = dst->btofPidTraits(ibtof);
      if(btof){
        cout << "btofMatchFlag=" << (int)btof->btofMatchFlag()
             << " beta=" << btof->btofBeta()
             << " tof="  << btof->btof()
             << endl;
      } else {
        cout << "btofPidTraits("<<ibtof<<") is null" << endl;
      }
    } else {
      cout << "No BTofPidTraitsIndex" << endl;
    }
  }

  chain->Finish();
}
