//============================================================
// Multi-PID log-polynomial fit for
// π / K / p / d / t / 4He
//============================================================

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
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
void Fit_PID_LogPoly()
{
  //==============================
  // Input
  //==============================
  TFile *fin = new TFile("DedxVsP_M2_production_3p85GeV_fixedTarget_2019.root");
  //TH2D  *h2   = (TH2D*)fin->Get("hDedxP");
  TH2D  *h2   = (TH2D*)fin->Get("hDedxP_cut");

  if (!h2) {
    cout << "ERROR: hDedxP not found" << endl;
    return;
  }

  //==============================
  // Particle definitions
  //==============================
  const int NPID = 6;
  const char* pidName[NPID] = {"pi","K","p","d","t","He4"};

  // rough dE/dx windows (keV/cm) — tune if needed
  double dedxMin[NPID] = {1, 1, 1, 1, 1, 1};
  double dedxMax[NPID] = {6, 30, 70, 80, 100, 120};

  TGraphErrors *gr[NPID];

  for (int i=0;i<NPID;i++)
    gr[i] = new TGraphErrors();

  //==============================
  // Momentum slicing
  //==============================
  int nBinsP = 60;
  double pMin = 0.2;
  double pMax = 6.0;
  double dp   = (pMax-pMin)/nBinsP;

  int ipoint[NPID] = {0};

  for (int ib=0; ib<nBinsP; ib++) {

    double p1 = pMin + ib*dp;
    double p2 = p1 + dp;
    int bin1 = h2->GetXaxis()->FindBin(p1);
    int bin2 = h2->GetXaxis()->FindBin(p2);

    TH1D *hproj = h2->ProjectionY(
				  Form("hproj_%d",ib), bin1, bin2);

    double pcenter = 0.5*(p1+p2);
    if (hproj->GetEntries()<200) continue;

    //==============================
    // Loop over PID species
    //==============================
    for (int ip=0; ip<NPID; ip++) {

      int y1 = h2->GetYaxis()->FindBin(dedxMin[ip]);
      int y2 = h2->GetYaxis()->FindBin(dedxMax[ip]);

      TH1D *hs = (TH1D*)hproj->Clone("hs");
      hs->GetXaxis()->SetRange(y1,y2);

      if (hs->GetEntries()<50) continue;

      TF1 *g = new TF1("g","gaus",
                       dedxMin[ip],dedxMax[ip]);
      hs->Fit(g,"QNR");

      double mean  = g->GetParameter(1);
      double sigma = g->GetParameter(2);

      gr[ip]->SetPoint(
		       ipoint[ip], pcenter, mean);
      gr[ip]->SetPointError(
			    ipoint[ip], dp/2, sigma);

      ipoint[ip]++;
    }
  }

  //==============================
  // Fit log-polynomial
  //==============================
  TCanvas *c = new TCanvas("c","PID log-polynomial",900,700);
  h2->GetXaxis()->SetRangeUser(0.01, 6.);
  h2->Draw("colz");

  int color[NPID]={kBlue,kGreen+2,kRed,kMagenta,kOrange+1,kBlack};

  for (int ip=0; ip<NPID; ip++) {

    TF1 *f = new TF1(
		     Form("f_%s",pidName[ip]),
		     LogPoly, 0.2, 4.0, 4);

    f->SetParameters(1, -1, 0.2, 0.0);
    f->SetLineColor(color[ip]);
    f->SetLineWidth(2);

    gr[ip]->Fit(f,"QR");
    f->Draw("same");

    gr[ip]->SetMarkerStyle(20);
    gr[ip]->SetMarkerColor(color[ip]);
    gr[ip]->Draw("P same");
  }

  c->SaveAs("PID_logpoly_pi_K_p_d_t_He4.pdf");

  cout << "=== Fit finished successfully ===" << endl;
}
