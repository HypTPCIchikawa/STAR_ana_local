//============================================================
// Multi-PID log-polynomial fit with
// TSpectrum peak search + multi-Gaussian
// π / K / p / d / t / 4He
//============================================================

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;

//--------------------------------------------------
// log-polynomial model
//--------------------------------------------------
double LogPoly(double *x, double *par)
{
  double lp = log10(x[0]);
  return par[0]
       + par[1]*lp
       + par[2]*lp*lp
    + par[3]*lp*lp*lp;
}

//--------------------------------------------------
void Fit_PID_LogPoly_PeakSearch()
{
  //==============================
  // Input
  //==============================
  TFile *fin = new TFile(
			 //"DedxVsP_M2_production_3p85GeV_fixedTarget_2019.root");
			 "DedxVsP_M2_production_13p5GeV_2020_tmp.root");
  //  TH2D *h2 = (TH2D*)fin->Get("hDedxP_cut");
  TH2D *h2 = (TH2D*)fin->Get("hDedxP_d");

  if (!h2) {
    cout << "ERROR: hDedxP_cut not found" << endl;
    return;
  }

  //==============================
  // Particle definitions
  //==============================
  const int NPID = 6;
  const char* pidName[NPID] = {"pi","K","p","d","t","He4"};
  int color[NPID]={kBlue,kGreen+2,kRed,kMagenta,kOrange+1,kBlack};

  TGraphErrors *gr[NPID];
  for (int i=0;i<NPID;i++) gr[i] = new TGraphErrors();

  //==============================
  // Momentum slicing
  //==============================
  int    nBinsP = 60;
  double pMin  = 0.2;
  double pMax  = 6.0;
  double dp    = (pMax-pMin)/nBinsP;

  int ipoint[NPID]={0};

  TCanvas *c0 = new TCanvas("c0","Slice hist",900,700);


  //==============================
  // Loop over p slices
  //==============================
  for (int ib=0; ib<nBinsP; ib++) {
    std::cout<<"ib="<<ib<<std::endl;
    double p1 = pMin + ib*dp;
    double p2 = p1 + dp;
    int bin1 = h2->GetXaxis()->FindBin(p1);
    int bin2 = h2->GetXaxis()->FindBin(p2);

    TH1D *hproj = h2->ProjectionY(
				  Form("hproj_%d",ib), bin1, bin2);

    if (hproj->GetEntries() < 100) continue;

    double pcenter = 0.5*(p1+p2);

    //==============================
    // Peak search
    //==============================
    TSpectrum spec(10);
    int nPeaks = spec.Search(hproj, 2, "", 0.05);

    if (nPeaks < 1) continue;

    double *xpeaks = spec.GetPositionX();
    vector<double> peaks;
    
    for (int i=0;i<nPeaks;i++) peaks.push_back(xpeaks[i]);
    sort(peaks.begin(), peaks.end()); // dE/dx order

    //==============================
    // Limit to expected PID count
    //==============================
    //int useN = min((int)peaks.size(), NPID);
    int useN = (int)peaks.size();

    //==============================
    // Build multi-Gaussian
    //==============================
    TString formula;
    for (int i=0;i<useN;i++) {
      formula += Form("gaus(%d)",3*i);
      if (i!=useN-1) formula += "+";
    }

    TF1 *f = new TF1("f_multi", formula,
                     0.1, 150.0);

    
    // Initial parameters from peaks
    for (int i=0;i<useN;i++) {
    

      double peak = peaks[i];
      double amp  = hproj->GetBinContent(
					 hproj->FindBin(peak));
      std::cout<<"amp="<<amp<<", peak="<<peak<<std::endl;
      f->SetParameter(3*i,   amp);
      f->SetParameter(3*i+1, peak);
      f->SetParameter(3*i+2, 0.1*peak);
    }

    //hproj->Fit(f,"QNR");
    hproj->Fit(f,"","",peaks[0]-0.2, peaks[useN-1]+0.2);
    c0->Update();
    getchar();
    

    //==============================
    // Store PID points
    //==============================
    for (int ip=0; ip<useN; ip++) {
      double mean  = f->GetParameter(3*ip+1);
      double sigma = f->GetParameter(3*ip+2);

      gr[ip]->SetPoint(ipoint[ip], pcenter, mean);
      gr[ip]->SetPointError(ipoint[ip], dp/2, sigma);
      ipoint[ip]++;
    }
  }

  //==============================
  // Draw & log-polynomial fit
  //==============================
  TCanvas *c = new TCanvas("c","PID log-poly",900,700);
  h2->Draw("colz");

  for (int ip=0; ip<NPID; ip++) {

    TF1 *flp = new TF1(
		       Form("f_%s",pidName[ip]),
		       LogPoly, 0.2, 4.0, 4);

    flp->SetParameters(1,-1,0.2,0.0);
    flp->SetLineColor(color[ip]);
    flp->SetLineWidth(2);

    gr[ip]->Fit(flp,"QR");
    flp->Draw("same");

    gr[ip]->SetMarkerStyle(20);
    gr[ip]->SetMarkerColor(color[ip]);
    gr[ip]->Draw("P same");
  }

  c->SaveAs("PID_logpoly_peaksearch_pi_K_p_d_t_He4.pdf");

  cout << "=== Peak-search PID fit finished ===" << endl;
}
