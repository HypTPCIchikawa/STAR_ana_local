#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>

#include "../include/HistConfig.h"
#include "../include/CutConfig.h"

class StChain;
class StPicoDstMaker;

StChain *chain = 0;

void readPicoDst(const Char_t *inputFile="test.list",
                 const Char_t *outputFile="test.root")
{
  TStopwatch timer;
  timer.Start();

  {
    TString out(outputFile);
    TString outDir = gSystem->DirName(out);
    if (outDir.Length() > 0 && gSystem->AccessPathName(outDir)) {
      gSystem->mkdir(outDir, kTRUE);
    }
  }

  TFile *fout = TFile::Open(outputFile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    ::Error("readPicoDst", "Cannot open output file: %s", outputFile);
    return;
  }

  // -----------------------------------------------------------------------
  // Histogram Definitions (Using Config namespace)
  // -----------------------------------------------------------------------
  
  // 1. Event Level
  TH1F *hVz = new TH1F("hVz", "Primary Vertex Z;V_{z} [cm];Counts", Config::nBinsVz, Config::minVz, Config::maxVz);
  TH2F *hVxVy = new TH2F("hVxVy", "Primary Vertex X vs Y;V_{x} [cm];V_{y} [cm]", Config::nBinsVxy, Config::minVxy, Config::maxVxy, Config::nBinsVxy, Config::minVxy, Config::maxVxy);
  TH1F *hRefMult = new TH1F("hRefMult", "RefMult Distribution;RefMult;Counts", Config::nBinsRefMult, Config::minRefMult, Config::maxRefMult);
  TH2F *hVzVsRun = new TH2F("hVzVsRun", "V_{z} vs RunID;RunID;V_{z} [cm]", Config::nBinsRun, Config::minRunId, Config::maxRunId, Config::nBinsVz, Config::minVz, Config::maxVz);
  TH2F *hRefMultVsVz = new TH2F("hRefMultVsVz", "RefMult vs V_{z};V_{z} [cm];RefMult", Config::nBinsVz, Config::minVz, Config::maxVz, Config::nBinsRefMult, Config::minRefMult, Config::maxRefMult);
  // Vz(TPC) - Vz(VPD) : Pileup rejection check
  TH1F *hVzDiff = new TH1F("hVzDiff", "V_{z}^{TPC} - V_{z}^{VPD};#Delta V_{z} [cm];Counts", 200, -10, 10); 
  // Trigger ID : To check physics triggers
  TH1F *hTriggerIds = new TH1F("hTriggerIds", "Trigger IDs;Trigger ID;Counts", 100, 0, 1000000); // Bin range depends on Run year

  // 2. Track Level
  TH1F *hPt = new TH1F("hPt", "Transverse Momentum p_{T};p_{T} [GeV/c];Counts", Config::nBinsPt, Config::minPt, Config::maxPt);
  TH1F *hEta = new TH1F("hEta", "Pseudorapidity #eta;#eta;Counts", Config::nBinsEta, Config::minEta, Config::maxEta);
  TH1F *hPhi = new TH1F("hPhi", "Azimuthal Angle #phi;#phi [rad];Counts", Config::nBinsPhi, Config::minPhi, Config::maxPhi);
  TH1F *hNHitsFit = new TH1F("hNHitsFit", "nHitsFit;nHitsFit;Counts", Config::nBinsHits, Config::minHits, Config::maxHits);
  TH1F *hNHitsRatio = new TH1F("hNHitsRatio", "nHitsFit / nHitsMax;Ratio;Counts", 100, 0, 1.05);
  TH1F *hDCA = new TH1F("hDCA", "Global DCA to PV;DCA [cm];Counts", Config::nBinsDCA, Config::minDCA, Config::maxDCA);
  TH1F *hCharge = new TH1F("hCharge", "Track Charge;Charge;Counts", 5, -2.5, 2.5);
  TH1F *hChi2 = new TH1F("hChi2", "Track #chi^{2};#chi^{2};Counts", 100, 0, 10);

  // 3. PID
  TH2F *hDedxVsP = new TH2F("hDedxVsP", "dE/dx vs Momentum;p [GeV/c];dE/dx [GeV/cm]", Config::nBinsPt, Config::minPt, Config::maxPt, Config::nBinsDedx, Config::minDedx, Config::maxDedx);
  TH2F *hNSigmaPionVsP = new TH2F("hNSigmaPionVsP", "n#sigma_{#pi} vs p;p [GeV/c];n#sigma_{#pi}", Config::nBinsPt, Config::minPt, Config::maxPt, 400, -20, 20);
  TH2F *hNSigmaKaonVsP = new TH2F("hNSigmaKaonVsP", "n#sigma_{K} vs p;p [GeV/c];n#sigma_{K}", Config::nBinsPt, Config::minPt, Config::maxPt, 400, -20, 20);
  TH2F *hNSigmaProtonVsP = new TH2F("hNSigmaProtonVsP", "n#sigma_{p} vs p;p [GeV/c];n#sigma_{p}", Config::nBinsPt, Config::minPt, Config::maxPt, 400, -20, 20);
  
  TH2F *hBetaVsP = new TH2F("hBetaVsP", "1/#beta vs p;p [GeV/c];1/#beta", Config::nBinsPt, Config::minPt, Config::maxPt, Config::nBinsBeta, Config::minBeta, Config::maxBeta); 
  TH2F *hMass2VsP = new TH2F("hMass2VsP", "m^{2} vs p;p [GeV/c];m^{2} [(GeV/c^{2})^{2}]", Config::nBinsPt, Config::minPt, Config::maxPt, Config::nBinsMass2, Config::minMass2, Config::maxMass2); 

  // 4. Event Plane / Flow (Raw check)
  TH1F *hPsi2 = new TH1F("hPsi2", "Raw Event Plane #Psi_{2};#Psi_{2} [rad];Counts", Config::nBinsPsi, Config::minPsi, Config::maxPsi);
  TH2F *hQxQy = new TH2F("hQxQy", "Q-Vector (Raw);Q_{x};Q_{y}", 200, -100, 100, 200, -100, 100);

  // 5. Detector Sanity
  TH1F *hBbcAdcSum = new TH1F("hBbcAdcSum", "BBC ADC Sum (East+West);ADC Sum;Counts", 1000, 0, 60000);
  TH1F *hTofMatchMult = new TH1F("hTofMatchMult", "Number of TOF Matched Tracks;N_{TOF};Counts", 500, 0, 500);

  TH1I *hN = new TH1I("hN", "N processed events;dummy;count", 1, 0, 1);

  // -----------------------------------------------------------------------
  // TTree Definitions for Resonance Reconstruction Analysis
  // -----------------------------------------------------------------------
  
  // Event Tree
  TTree *eventTree = new TTree("EventTree", "Event information for resonance analysis");
  Float_t ev_Vz, ev_Vx, ev_Vy, ev_Vr;
  Int_t ev_refMult, ev_runId, ev_eventId;
  Float_t ev_centrality;  // Will be calculated later, set to -1 for now
  Float_t ev_Qx, ev_Qy, ev_psi2;
  Int_t ev_nTracks;
  Float_t ev_vzVpd;
  
  eventTree->Branch("Vz", &ev_Vz, "Vz/F");
  eventTree->Branch("Vx", &ev_Vx, "Vx/F");
  eventTree->Branch("Vy", &ev_Vy, "Vy/F");
  eventTree->Branch("Vr", &ev_Vr, "Vr/F");
  eventTree->Branch("refMult", &ev_refMult, "refMult/I");
  eventTree->Branch("runId", &ev_runId, "runId/I");
  eventTree->Branch("eventId", &ev_eventId, "eventId/I");
  eventTree->Branch("centrality", &ev_centrality, "centrality/F");
  eventTree->Branch("Qx", &ev_Qx, "Qx/F");
  eventTree->Branch("Qy", &ev_Qy, "Qy/F");
  eventTree->Branch("psi2", &ev_psi2, "psi2/F");
  eventTree->Branch("nTracks", &ev_nTracks, "nTracks/I");
  eventTree->Branch("vzVpd", &ev_vzVpd, "vzVpd/F");
  
  // Track Tree
  TTree *trackTree = new TTree("TrackTree", "Track information for resonance analysis");
  Int_t tr_eventIndex;
  Float_t tr_pT, tr_eta, tr_phi;
  Short_t tr_charge;
  Short_t tr_nHitsFit, tr_nHitsMax, tr_nHitsDedx;
  Float_t tr_DCA, tr_chi2;
  Float_t tr_nSigmaPion, tr_nSigmaKaon, tr_nSigmaProton;
  Float_t tr_beta, tr_mass2;
  Bool_t tr_tofMatch;
  
  trackTree->Branch("eventIndex", &tr_eventIndex, "eventIndex/I");
  trackTree->Branch("pT", &tr_pT, "pT/F");
  trackTree->Branch("eta", &tr_eta, "eta/F");
  trackTree->Branch("phi", &tr_phi, "phi/F");
  trackTree->Branch("charge", &tr_charge, "charge/S");
  trackTree->Branch("nHitsFit", &tr_nHitsFit, "nHitsFit/S");
  trackTree->Branch("nHitsMax", &tr_nHitsMax, "nHitsMax/S");
  trackTree->Branch("nHitsDedx", &tr_nHitsDedx, "nHitsDedx/S");
  trackTree->Branch("DCA", &tr_DCA, "DCA/F");
  trackTree->Branch("chi2", &tr_chi2, "chi2/F");
  trackTree->Branch("nSigmaPion", &tr_nSigmaPion, "nSigmaPion/F");
  trackTree->Branch("nSigmaKaon", &tr_nSigmaKaon, "nSigmaKaon/F");
  trackTree->Branch("nSigmaProton", &tr_nSigmaProton, "nSigmaProton/F");
  trackTree->Branch("beta", &tr_beta, "beta/F");
  trackTree->Branch("mass2", &tr_mass2, "mass2/F");
  trackTree->Branch("tofMatch", &tr_tofMatch, "tofMatch/O");
  
  // Set compression for trees
  eventTree->SetAutoSave(1000000);  // Auto-save every 1M entries
  trackTree->SetAutoSave(10000000); // Auto-save every 10M entries

  Long64_t nEvents = 10000000;
  Int_t currentEventIndex = 0;  // Counter for event indexing

  // ---- Load STAR shared libraries (root4star) ----
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // PicoDst
  if (gSystem->Load("StPicoEvent") < 0) {
    std::cerr << "ERROR: failed to load StPicoEvent" << std::endl;
    return;
  }
  if (gSystem->Load("StPicoDstMaker") < 0) {
    if (gSystem->Load("StPicoDstMaker.so") < 0) {
      std::cerr << "ERROR: failed to load StPicoDstMaker" << std::endl;
      return;
    }
  }

  // ---- Build chain ----
  chain = new StChain("StChain");

  StPicoDstMaker *picoMaker =
    new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");

  picoMaker->SetStatus("*", 0);
  picoMaker->SetStatus("Event", 1);
  picoMaker->SetStatus("Track", 1);
  picoMaker->SetStatus("BTofHit", 1);
  picoMaker->SetStatus("BTofPidTraits", 1);
  picoMaker->SetStatus("BbcHit", 1);
  picoMaker->SetStatus("EpdHit", 1);
  picoMaker->SetStatus("MtdHit", 1);
  picoMaker->SetStatus("BTowHit", 1);
  picoMaker->SetStatus("ETofPidTraits", 1);

  // ---- Init ----
  if (chain->Init() == kStErr) {
    std::cerr << "ERROR: chain->Init() returned kStErr" << std::endl;
    fout->Close();
    delete picoMaker;
    delete chain;
    chain = 0;
    return;
  }

  // ---- Determine entries ----
  Long64_t totalEntries = picoMaker->chain() ? picoMaker->chain()->GetEntries() : 0;
  std::cout << "Total entries = " << totalEntries << std::endl;

  if (totalEntries <= 0) {
    std::cerr << "ERROR: no entries found. Check inputFile list/root files." << std::endl;
    chain->Finish();
    delete picoMaker;
    delete chain;
    chain = 0;
    return;
  }

  if (nEvents > totalEntries) nEvents = totalEntries;

  // ---- Event loop ----
  for (Long64_t i = 0; i < nEvents; ++i) {
    if (i % 1000 == 0) std::cout << "Working on event " << i << std::endl;

    chain->Clear();
    Int_t iret = chain->Make(i);
    if (iret) {
      std::cerr << "Bad return code: " << iret << " at event " << i << std::endl;
      break;
    }

    StPicoDst *dst = picoMaker->picoDst();
    if (!dst) {
      continue;
    }

    StPicoEvent *event = dst->event();
    if (!event) {
      continue;
    }

    // --- 1. Event Level Fill ---
    TVector3 pVtx = event->primaryVertex();    
    // Quick Event Cuts (Optional: Add loose cuts if necessary, e.g. |Vz| < 100)
    // if (fabs(pVtx.z()) > 100) continue;
    
    hVz->Fill(pVtx.z());
    hVxVy->Fill(pVtx.x(), pVtx.y());
    hRefMult->Fill(event->refMult());
    
    // Use float runID for plotting if convenient, or just check binning
    hVzVsRun->Fill(event->runId(), pVtx.z());
    hRefMultVsVz->Fill(pVtx.z(), event->refMult());

    Float_t vzVpd = event->vzVpd(); // Available from tree branch 43    
    // Check if VPD Vz is available (usually -999 or similar if undefined)
    if (fabs(vzVpd) < 200) { 
      hVzDiff->Fill(pVtx.z() - vzVpd);
    }
    
    // Trigger ID Fill
    std::vector<unsigned int> triggerIds = event->triggerIds();
    for(UInt_t itrig = 0; itrig < triggerIds.size(); itrig++) {
      hTriggerIds->Fill(triggerIds[itrig]);
    }
    
    // --- Prepare Event Tree variables ---
    ev_Vz = pVtx.z();
    ev_Vx = pVtx.x();
    ev_Vy = pVtx.y();
    ev_Vr = TMath::Sqrt(pVtx.x()*pVtx.x() + pVtx.y()*pVtx.y());
    ev_refMult = event->refMult();
    ev_runId = event->runId();
    ev_eventId = event->eventId();
    ev_centrality = -1.0;  // To be calculated later
    ev_vzVpd = vzVpd;
    ev_nTracks = 0;  // Will be filled after track loop

    // --- 5. Detector Sanity (Event level part) ---
    // BBC Sum - tree_print shows mBbcAdcEast/West [cite: 71, 72]
    // StPicoEvent usually has accessor like bbcTriggerStat() or similar, 
    // but often users sum ADCs manually if not directly available in aggregate.
    // Here assuming typical user check might just use existing RefMult or look deeper.
    // For simple check, we skip raw BBC loop to save time unless requested explicitly.
    
    // --- Track Loop Variables for Event Plane ---
    Double_t Qx = 0.0;
    Double_t Qy = 0.0;
    Int_t nTofMatch = 0;

    Int_t nTracks = dst->numberOfTracks();
    Int_t nTracksSaved = 0;  // Counter for tracks saved to tree

    // --- 2. Track Level Loop ---
    for (Int_t itrk = 0; itrk < nTracks; itrk++) {
      StPicoTrack *trk = dst->track(itrk);
      if (!trk) continue;

      // Primary Momentum check:
      // PicoDst default stores global tracks. 
      // Need to ensure we use primary momentum for physics if available, 
      // or global momentum for QA. Let's use Global for raw QA as it's always there.
      TVector3 mom = trk->gMom(pVtx, event->bField()); 
      // Note: tree_print [cite: 85] shows mPMomentumX exists, so primary tracks are stored.
      // If you want primary track explicitly:
      TVector3 pMom = trk->pMom(); 

      // Basic Cuts for QA (visualize everything, or cut duplicates?)
      // Let's visualize tracks that are at least decent
      if (trk->nHitsFit() < 10) continue; // Loose cut

      Float_t pt = pMom.Perp();
      Float_t eta = pMom.PseudoRapidity();
      Float_t phi = pMom.Phi();

      // Skip non-primary tracks (pMom is zero for non-primaries usually in newer PicoDst)
      if (fabs(pMom.Mag()) < 1e-4) continue; 

      hPt->Fill(pt);
      hEta->Fill(eta);
      hPhi->Fill(phi);
      hNHitsFit->Fill(trk->nHitsFit());
      hNHitsRatio->Fill((Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax());
      hDCA->Fill(trk->gDCA(pVtx).Mag());
      hCharge->Fill(trk->charge());
      hChi2->Fill(trk->chi2());
      
      // --- 3. PID (TPC) ---
      hDedxVsP->Fill(pMom.Mag(), trk->dEdx());
      hNSigmaPionVsP->Fill(pMom.Mag(), trk->nSigmaPion());
      hNSigmaKaonVsP->Fill(pMom.Mag(), trk->nSigmaKaon());
      hNSigmaProtonVsP->Fill(pMom.Mag(), trk->nSigmaProton());

      // --- 4. Event Plane Accumulation (Raw, no weights) ---
      // Usually exclude high pt or specific eta region
      if (pt > 0.15 && pt < 2.0 && fabs(eta) < 1.0) {
	Qx += cos(2.0 * phi);
	Qy += sin(2.0 * phi);
      }

      // --- Fill Track Tree (save all tracks, cuts applied later in analysis) ---
      tr_eventIndex = currentEventIndex;
      tr_pT = pt;
      tr_eta = eta;
      tr_phi = phi;
      tr_charge = trk->charge();
      tr_nHitsFit = trk->nHitsFit();
      tr_nHitsMax = trk->nHitsMax();
      tr_nHitsDedx = trk->nHitsDedx();
      tr_DCA = trk->gDCA(pVtx).Mag();
      tr_chi2 = trk->chi2();
      tr_nSigmaPion = trk->nSigmaPion();
      tr_nSigmaKaon = trk->nSigmaKaon();
      tr_nSigmaProton = trk->nSigmaProton();
      
      // TOF information
      Int_t btofIndex = trk->bTofPidTraitsIndex();
      if (btofIndex >= 0) {
	StPicoBTofPidTraits *tof = dst->btofPidTraits(btofIndex);
	if (tof) {
	  nTofMatch++;
	  Double_t beta = tof->btofBeta();
	  if (beta > 1e-4) {
	    Double_t oneOverBeta = 1.0 / beta;
	    Double_t mass2 = pMom.Mag2() * (oneOverBeta*oneOverBeta - 1.0);
	    hBetaVsP->Fill(pMom.Mag(), oneOverBeta);
	    hMass2VsP->Fill(pMom.Mag(), mass2);
	    
	    tr_beta = beta;
	    tr_mass2 = mass2;
	    tr_tofMatch = kTRUE;
	  } else {
	    tr_beta = -999.0;
	    tr_mass2 = -999.0;
	    tr_tofMatch = kFALSE;
	} 
	} else {
	  tr_beta = -999.0;
	  tr_mass2 = -999.0;
	  tr_tofMatch = kFALSE;
	}
      } else {
	tr_beta = -999.0;
	tr_mass2 = -999.0;
	tr_tofMatch = kFALSE;
      }
      
      trackTree->Fill();
      nTracksSaved++;
    } // End Track Loop

    // --- 4. Event Plane Fill ---
    TVector2 Q(Qx, Qy);
    hQxQy->Fill(Qx, Qy);
    if (Q.Mod() > 0) {
      Double_t psi2 = 0.5 * TMath::ATan2(Qy, Qx);
      // Shift to 0 - pi
      if (psi2 < 0) psi2 += TMath::Pi();
      hPsi2->Fill(psi2);
      ev_psi2 = psi2;
    } else {
      ev_psi2 = -999.0;
    }
    ev_Qx = Qx;
    ev_Qy = Qy;
    ev_nTracks = nTracksSaved;
	
    // --- 5. Detector Sanity Fill ---
    hTofMatchMult->Fill(nTofMatch);
    hN->Fill(0);
    
    // --- Fill Event Tree ---
    eventTree->Fill();
    currentEventIndex++;
  }

  // ---- Finish ----
  chain->Finish();

  fout->cd();
  
  // Write Trees
  eventTree->Write();
  trackTree->Write();
  
  // Write Histograms
  hN->Write();
  hVz->Write();
  hVxVy->Write();
  hRefMult->Write();
  hVzVsRun->Write();
  hRefMultVsVz->Write();
  hVzDiff->Write();
  hTriggerIds->Write();

  hPt->Write();
  hEta->Write();
  hPhi->Write();
  hNHitsFit->Write();
  hNHitsRatio->Write();
  hDCA->Write();
  hCharge->Write();
  hChi2->Write();

  hDedxVsP->Write();
  hNSigmaPionVsP->Write();
  hNSigmaKaonVsP->Write();
  hNSigmaProtonVsP->Write();
  hBetaVsP->Write();
  hMass2VsP->Write();

  hPsi2->Write();
  hQxQy->Write();  
  hTofMatchMult->Write();

  fout->Write();
  fout->Close();

  timer.Stop();
  std::cout << "Processed events: " << nEvents << std::endl;
  std::cout << "RealTime: " << timer.RealTime()
            << " CpuTime: " << timer.CpuTime() << std::endl;

  delete picoMaker;
  delete chain;
  chain = 0;
}
