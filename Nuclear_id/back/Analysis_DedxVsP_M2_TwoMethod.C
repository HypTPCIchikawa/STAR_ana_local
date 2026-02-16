//============================================================
// picoDst analysis
//  p vs m^2 (two methods)
//   1) m2_simple : using btofBeta()
//   2) m2_path   : using recalculated TOF path length
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
#include "TMath.h"


#include <iostream>
using namespace std;

StChain *chain = 0;

// ------------------------------------------------------------
// simple helix path-length approximation
// ------------------------------------------------------------
float tofPathLength(const StThreeVectorF* vtx,
                    const StThreeVectorF* hit,
                    float curvature)
{
  float dx = hit->x() - vtx->x();
  float dy = hit->y() - vtx->y();
  float dz = hit->z() - vtx->z();
  float dxy = ::sqrt(dx*dx + dy*dy);

  if (fabs(curvature) < 1e-5) {
    return ::sqrt(dxy*dxy + dz*dz);
  }

  float R = 1.0 / fabs(curvature);
  float phi = dxy / R;
  float Lxy = R * phi;

  return ::sqrt(Lxy*Lxy + dz*dz);
}

// ------------------------------------------------------------
void Analysis_DedxVsP_M2_TwoMethod(Int_t nEvents = 0,
				   const char* inList = "test.list",
				   const char* outFile = "m2_compare.root")
{
  if (nEvents == 0) nEvents = 1000000000;

  // Load STAR libraries
  gROOT->LoadMacro(
		   "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");

  chain = new StChain();
  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(2, inList, "PicoDst");

  picoMaker->SetStatus("*",0);
  picoMaker->SetStatus("Event",1);
  picoMaker->SetStatus("Track",1);
  picoMaker->SetStatus("BTofPidTraits",1);

  if (chain->Init() == kStErr) {
    cout << "chain->Init() failed" << endl;
    return;
  }
  


  int total = picoMaker->chain()->GetEntries();
  if (nEvents > total) nEvents = total;
  cout << "Total entries = " << total << endl;

  // Histograms
  TH2D *hP_M2_simple =
    new TH2D("hP_M2_simple",
             "p vs m^{2} (btofBeta);p (GeV/c);m^{2} (GeV^{2}/c^{4})",
             400,0,4, 400,-0.5,5);

  TH2D *hP_M2_path =
    new TH2D("hP_M2_path",
             "p vs m^{2} (path-length corrected);p (GeV/c);m^{2} (GeV^{2}/c^{4})",
             400,0,4, 400,-0.5,5);

  const float c_light = 29.9792458; // cm/ns

  // Event loop
  for (int i=0;i<nEvents;i++) {

    if (i%1000==0) cout << "Event " << i << endl;

    chain->Clear();
    if (chain->Make(i)) break;

    StPicoDst *dst = picoMaker->picoDst();
    if (!dst) continue;

    StPicoEvent *event = dst->event();
    if (!event) continue;

    TVector3 vtx = event->primaryVertex();
    float BField = event->bField();

    for (int it=0; it<dst->numberOfTracks(); it++) {

      StPicoTrack *trk = dst->track(it);
      if (!trk) continue;

      if (trk->nHitsDedx() < 15) continue;
      if (trk->gPt() < 0.2) continue;

      TVector3 pvec = trk->pMom();
      float p = pvec.Mag();

      int idx = trk->btofPidTraitsIndex();
      if (idx < 0) continue;

      StPicoBTofPidTraits *tofPid = dst->btofPidTraits(idx);
      if (!tofPid) continue;

      // ----------------------------
      // m2_simple (collider style)
      // ----------------------------
      float beta1 = tofPid->btofBeta();
      if (beta1 > 0)
        hP_M2_simple->Fill(p, p*p*(1.0/(beta1*beta1)-1.0));

      // ----------------------------
      // m2_path (vertex-aware)
      // ----------------------------
      TVector3 tofPos = tofPid->btofHitPos();
      StThreeVectorF vtxSt(vtx.X(),vtx.Y(),vtx.Z());
      StThreeVectorF hitSt(tofPos.X(),tofPos.Y(),tofPos.Z());
      StThreeVectorF gmom(pvec.X(),pvec.Y(),pvec.Z());
      StThreeVectorF org(trk->origin().X(),
                         trk->origin().Y(),
                         trk->origin().Z());

      StPhysicalHelixD helix(gmom, org,
                             BField*kilogauss,
                             trk->charge());

      float L = tofPathLength(&vtxSt, &hitSt, helix.curvature());
      float tof = tofPid->btof();

      if (tof > 0) {
        float beta2 = L / (tof * c_light);
        if (beta2 > 0)
          hP_M2_path->Fill(p, p*p*(1.0/(beta2*beta2)-1.0));
      }
    }
  }

  // Save
  TFile *fout = new TFile(outFile,"RECREATE");
  hP_M2_simple->Write();
  hP_M2_path->Write();
  fout->Close();

  chain->Finish();
  delete picoMaker;
  delete chain;

  cout << "Done. Output = " << outFile << endl;
}
