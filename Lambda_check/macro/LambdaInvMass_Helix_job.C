////////////////////////////////////////////////////////////
//
//  LambdaInvMass_Helix.C
//  - CINT compatible
//  - ROOT5 / root4star
//  - Helix-based full Lambda reconstruction
//  - starver SL24y
////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPhysicalHelixD;

StChain* chain = 0;

//======================================================
// Helix 作成（CINT安全）
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
// Helix V0 再構成 + 運動量取得
//======================================================
bool MakeLambdaHelix(const StPicoTrack* p,
                     const StPicoTrack* pi,
                     double bField,
                     TVector3& v0,
                     TVector3& momP,
                     TVector3& momPi,
                     double& dca12)
{
  StPhysicalHelixD hp  = makeHelix(p,  bField);
  StPhysicalHelixD hpi = makeHelix(pi, bField);

  pair<double,double> s=hp.pathLengths(hpi);
  if (fabs(s.first) > 100 || fabs(s.second) > 100) return false;
  
  StThreeVectorD dcaA = hp.at(s.first);
  StThreeVectorD dcaB = hpi.at(s.second);

  StThreeVectorD v0_  = (dcaA+dcaB)*0.5;
  float dca12 = (dcaA-dcaB).mag();

  if (dca12 < 0 || dca12 > 1.0) return false;


  // ---- momentum at V0 ----
  StThreeVectorD pp  = hp.momentumAt(s.first,  bField * kilogauss);
  StThreeVectorD ppi = hpi.momentumAt(s.second, bField * kilogauss);

  momP.SetXYZ(pp.x(), pp.y(), pp.z());
  momPi.SetXYZ(ppi.x(), ppi.y(), ppi.z());

  v0.SetXYZ(v0_.x(), v0_.y(), v0_.z());

  return true;
}

//======================================================
// Main analysis
//======================================================
// void LambdaInvMass_Helix(Int_t nEvents = 1000,
//                          //const char* inList = "test2.list")
// 			 const char* inList = "test.list")

void LambdaInvMass_Helix_job(Int_t nEvents = 1000,
			     const char* inList,
			     const char* jobid)
{
  if (nEvents == 0) nEvents = 1000000000;
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
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

  chain = new StChain();
  StPicoDstMaker* picoMaker =
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

  TH1D* hLam = new TH1D("hLam",
			"Lambda invariant mass (Helix);M_{p#pi^{-}} (GeV/c^{2});Counts",
			1000, 1.05, 1.45);

  const double mp  = 0.938272;
  const double mpi = 0.139570;


  for (Int_t i = 0; i < nEvents; i++) {

    //    if (i % 1000 == 0)
    cout << "Working on event " << i << endl;
    
    chain->Clear();
    Int_t iret = chain->Make(i);
    if (iret) {
      cout << "Bad return code: " << iret << endl;
      break;
    }
    
    StPicoDst *picoDst = picoMaker->picoDst();
    if (!picoDst) continue;
    
    StPicoEvent *evt = picoDst->event();
    if (!evt) continue;
    
      
    TVector3 pv = evt->primaryVertex();
    //if (fabs(pv.z()) > 200) continue;
      
    double bField = evt->bField();
    int nTr = picoDst->numberOfTracks();
    std::cout<<"nTr="<<nTr<<std::endl;
    if (nTr > 300) continue;
      
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
	
	TVector3 v0, momP, momPi;
	double dca12;
	
	if (!MakeLambdaHelix(p, pi, bField,
			     v0, momP, momPi, dca12))
	  continue;
	
	TVector3 pLam = momP + momPi;
	double pLamMag = pLam.Mag();
	if (pLamMag < 1e-5) continue;
	
	TVector3 pLamUnit = pLam * (1.0 / pLamMag);
	
	// DCA of Lambda to PV
	TVector3 diff = pv - v0;
	double dcaV0 = (diff.Cross(pLamUnit)).Mag();
	if (dcaV0 > 1.0) continue;   // 典型値: 0.5–1.0 cm
	
	TVector3 flight = v0 - pv;
	double cosPoint =
	  flight.Dot(pLam) / (flight.Mag() * pLam.Mag() + 1e-10);
	
	if (cosPoint < 0.995) continue;   // BES energies

	TLorentzVector lp, lpi;
	lp.SetVectM(momP,  mp);
	lpi.SetVectM(momPi, mpi);
	
	hLam->Fill((lp + lpi).M());
		
      }
    }
  }
  
  TString outname = Form("%s/Nuclear_id/LambdaInvMass_Helix_%s.root",
                         gSystem->Getenv("SCRATCH"), jobid);
  TFile *fout = new TFile(outname, "RECREATE");
  hLam->Write();
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
