/////////////////////////////////////////////////////////
// Extract d and 4He dE/dx band from TPC dE/dx vs p
/////////////////////////////////////////////////////////

void Fit_PID_Set()
{
  TFile *f = new TFile("DedxVsP_M2_production_13p5GeV_2020_tmp.root");
  TH2D *h = (TH2D*)f->Get("hDedxP_d");

  const int binStep = 50;
  int nPx = h->GetNbinsX();

  // =======================
  // Seed table from image digitization
  // =======================
  const int Nseed = 50;
  double p_seed[Nseed];
  double d_seed[Nseed];
  double he_seed[Nseed];

  for(int i=0;i<Nseed;i++){
    p_seed[i] = (i+0.5)*0.1;
    d_seed[i]  = 45./pow(p_seed[i],1.8) + 12.;
    he_seed[i] = 140./pow(p_seed[i],1.8) + 12.;
  }

  vector<double> p_center;
  vector<double> dedx_d;
  vector<double> dedx_he;

  for(int ibin=1; ibin<nPx; ibin+=binStep){

    int binLow = ibin;
    int binHigh = ibin + binStep - 1;
    if(binHigh > nPx) binHigh = nPx;

    // p center
    double p1 = h->GetXaxis()->GetBinCenter(binLow);
    double p2 = h->GetXaxis()->GetBinCenter(binHigh);
    double pcent = 0.5*(p1+p2);

    // find nearest seed index
    int iseed = (int)(pcent/0.1);
    if(iseed<0) iseed=0;
    if(iseed>=Nseed) iseed=Nseed-1;

    double d_guess  = d_seed[iseed];
    double he_guess = he_seed[iseed];

    // projection
    TH1D *proj = h->ProjectionY(Form("proj_%d",ibin), binLow, binHigh);
    if(proj->GetEntries() < 200) continue;

    // Fit function
    TF1 *f2 = new TF1(Form("f2_%d",ibin),
		      "gaus(0)+gaus(3)", 5, 200);

    // Initial parameters
    f2->SetParameters(proj->GetMaximum(), d_guess, 3,
		      proj->GetMaximum()/10., he_guess, 6);

    // Parameter constraints
    f2->SetParLimits(1, d_guess-8, d_guess+8);
    f2->SetParLimits(2, 1, 10);

    f2->SetParLimits(4, he_guess-20, he_guess+20);
    f2->SetParLimits(5, 2, 20);

    // Fit
    proj->Fit(f2,"RQ0");

    double d_mean  = f2->GetParameter(1);
    double he_mean = f2->GetParameter(4);

    p_center.push_back(pcent);
    dedx_d.push_back(d_mean);
    dedx_he.push_back(he_mean);
  }

  // =======================
  // Draw result
  // =======================
  TCanvas *c = new TCanvas("c","",1000,800);
  h->Draw("colz");

  TGraph *gD  = new TGraph(p_center.size());
  TGraph *gHe = new TGraph(p_center.size());

  for(int i=0;i<p_center.size();i++){
    gD ->SetPoint(i, p_center[i], dedx_d[i]);
    gHe->SetPoint(i, p_center[i], dedx_he[i]);
  }

  gD ->SetMarkerStyle(20);
  gD ->SetMarkerColor(kRed);
  gD ->SetLineColor(kRed);

  gHe->SetMarkerStyle(21);
  gHe->SetMarkerColor(kBlue);
  gHe->SetLineColor(kBlue);

  gD ->Draw("P same");
  gHe->Draw("P same");

  c->SaveAs("dedx_band_fit.png");
}
