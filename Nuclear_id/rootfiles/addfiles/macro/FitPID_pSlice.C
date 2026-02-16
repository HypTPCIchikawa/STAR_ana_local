
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

void FitPID_pSlice()
{
  
  TFile *fin = new TFile("DedxVsP_M2_production_3p85GeV_fixedTarget_2019.root");
  TH2D  *h2  = (TH2D*)fin->Get("hDedxP");
  //TH2D  *h2   = (TH2D*)fin->Get("hDedxP_cut");

  if (!h2) {
    cout << "ERROR: hDedxP not found" << endl;
    return;
  }


  const int NPID = 6;
  const char* pname[NPID] = {"pi","K","p","d","t","He4"};

  const double pmin = 0.2;
  const double pmax = 6.0;
  const double dp   = 0.05;

  TGraphErrors* gr[NPID];
  for(int i=0;i<NPID;i++)
    gr[i] = new TGraphErrors();

  int ipoint[NPID] = {0};

  // ---- p-slice loop ----
  for(double p = pmin; p < pmax; p += dp) {
    std::cout<<"p="<<p<<std::endl;
    int bin1 = h2->GetXaxis()->FindBin(p);
    int bin2 = h2->GetXaxis()->FindBin(p+dp);

    TH1D* h = h2->ProjectionY(
			      Form("hdedx_p_%.2f",p), bin1, bin2);
    
    if(h->GetEntries() < 100) continue;

    // log(dE/dx)
    for(int ib=1; ib<=h->GetNbinsX(); ib++){
      double y = h->GetBinCenter(ib);
      if(y>0)
        h->SetBinContent(ib, h->GetBinContent(ib));
      else
        h->SetBinContent(ib,0);
    }

    // ---- multi-Gaussian ----
    TF1* f = new TF1("f",
		     "gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",
		     log(2), log(200));

    // 初期値（重要）
    double mean0[NPID] = {
      log(2.5),  // pi
      log(3.5),  // K
      log(5.0),  // p
      log(10.0), // d
      log(18.0), // t
      log(30.0)  // He4
    };

    for(int i=0;i<NPID;i++){
      f->SetParameter(3*i+0, h->GetMaximum()/NPID);
      f->SetParameter(3*i+1, mean0[i]);
      f->SetParameter(3*i+2, 0.15);

      f->SetParLimits(3*i+2, 0.05, 0.3);
    }

    h->Fit(f,"QNR");

    // ---- store result ----
    for(int i=0;i<NPID;i++){
      double mu  = f->GetParameter(3*i+1);
      double emu = f->GetParError (3*i+1);

      gr[i]->SetPoint     (ipoint[i], p+dp/2., mu);
      gr[i]->SetPointError(ipoint[i], dp/2., emu);
      ipoint[i]++;
    }
  }

  // ---- global log-polynomial fit ----
  TCanvas* c = new TCanvas("cPID","PID",1200,800);
  c->Divide(3,2);

  for(int i=0;i<NPID;i++){
    c->cd(i+1);
    gr[i]->SetTitle(Form("%s",pname[i]));
    gr[i]->SetMarkerStyle(20);

    TF1* fbb = new TF1(
		       Form("f_%s",pname[i]),
		       "[0]+[1]*log(x)+[2]*log(x)*log(x)",
		       pmin,pmax);

    fbb->SetParameters(1,1,0);
    gr[i]->Fit(fbb,"R");
    gr[i]->Draw("AP");
  }
}
