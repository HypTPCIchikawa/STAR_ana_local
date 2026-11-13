#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TBox.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TUUID.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace SEME {
void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error(message);
}
int edge(const TAxis *axis, double value) {
  for (int i=1; i<=axis->GetNbins()+1; ++i)
    if (std::fabs(axis->GetBinLowEdge(i)-value)<1e-10) return i;
  throw std::runtime_error("Requested boundary is not a histogram edge");
}
TH2 *get2(TFile &f, const char *name) {
  TH2 *h=dynamic_cast<TH2*>(f.Get(name));
  require(h!=0, std::string("Missing TH2: ")+name);
  h->GetXaxis()->UnZoom(); h->GetYaxis()->UnZoom();
  return h;
}
void checkAxes(const TH1 *a, const TH1 *b) {
  require(a->GetNbinsX()==b->GetNbinsX(),"Bin count mismatch");
  for (int i=1; i<=a->GetNbinsX()+1; ++i)
    require(a->GetXaxis()->GetBinLowEdge(i)==b->GetXaxis()->GetBinLowEdge(i),"Bin edge mismatch");
}
void exportBins(std::ofstream &stream, const std::string &id, const char *kind,
                const TH1 *se, const TH1 *me, const TH1 *scaled) {
  for (int i=0; i<=se->GetNbinsX()+1; ++i)
    stream << id << '\t' << kind << '\t' << i << '\t'
      << se->GetXaxis()->GetBinLowEdge(i) << '\t' << se->GetXaxis()->GetBinUpEdge(i) << '\t'
      << se->GetBinContent(i) << '\t' << me->GetBinContent(i) << '\t' << scaled->GetBinContent(i) << '\t'
      << se->GetBinError(i) << '\t' << me->GetBinError(i) << '\t' << scaled->GetBinError(i) << '\n';
}
void panel(TH1D *se, TH1D *me, double lo, double hi, double nlo, double nhi, double factor) {
  gPad->SetLeftMargin(.15); gPad->SetRightMargin(.035);
  gPad->SetTopMargin(.035); gPad->SetBottomMargin(.125);
  TH1D *s=static_cast<TH1D*>(se->Clone());
  TH1D *m=static_cast<TH1D*>(me->Clone());
  s->SetDirectory(0); m->SetDirectory(0);
  int first=edge(s->GetXaxis(),lo), last=edge(s->GetXaxis(),hi)-1;
  double ymax=0.;
  for (int i=first; i<=last; ++i)
    ymax=std::max(ymax,std::max(s->GetBinContent(i)+s->GetBinError(i),m->GetBinContent(i)));
  require(ymax>0.,"Empty display range");
  s->SetTitle(""); s->SetStats(0); s->SetMinimum(0.); s->SetMaximum(ymax*1.34);
  s->GetXaxis()->SetRange(first,last);
  s->GetXaxis()->SetTitle("k^{*} (GeV/c)");
  s->GetYaxis()->SetTitle(Form("Pairs / %.3f GeV/c",s->GetBinWidth(1)));
  s->GetXaxis()->SetTitleSize(.047); s->GetYaxis()->SetTitleSize(.047);
  s->GetXaxis()->SetLabelSize(.04); s->GetYaxis()->SetLabelSize(.04);
  s->GetYaxis()->SetTitleOffset(1.55); s->GetXaxis()->SetTitleOffset(1.1);
  s->SetMarkerStyle(20); s->SetMarkerSize(.65); s->SetMarkerColor(kBlack); s->SetLineColor(kBlack);
  m->SetLineColor(kRed+1); m->SetLineWidth(2); m->SetStats(0);
  // Establish axes, then shade the normalization interval behind both spectra.
  s->Draw("AXIS");
  TBox *box=new TBox(nlo,0.,nhi,ymax*1.34); box->SetFillColor(19); box->SetLineColor(19); box->Draw("SAME");
  m->Draw("HIST SAME"); s->Draw("E1 SAME");
  TLegend *legend=new TLegend(.19,.755,.75,.95);
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(.038);
  legend->AddEntry(s,"SE: raw pair counts","lep");
  legend->AddEntry(m,Form("ME #times %.8f",factor),"l");
  legend->AddEntry(box,Form("Match integrals: %.1f < k^{*} < %.1f",nlo,nhi),"f");
  legend->Draw();
  gPad->RedrawAxis();
}
}

void DrawSEME(const char *input, const char *outputDir, const char *groupFile,
              const char *seName, const char *meName, const char *seQA, const char *meQA,
              const char *expectedUUID, double normLo, double normHi, int rebin,
              double fullLo, double fullHi, double zoomLo, double zoomHi) {
  try {
    gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
    gStyle->SetCanvasColor(0); gStyle->SetPadColor(0); gStyle->SetFrameFillColor(0);
    TH1::AddDirectory(kFALSE);
    TFile inputFile(input,"READ");
    SEME::require(!inputFile.IsZombie()&&!inputFile.TestBit(TFile::kRecovered),"Invalid or recovered input ROOT");
    SEME::require(std::string(inputFile.GetUUID().AsString())==expectedUUID,"Input UUID changed");
    TH2 *sourceSE=SEME::get2(inputFile,seName), *sourceME=SEME::get2(inputFile,meName);
    SEME::checkAxes(sourceSE,sourceME);
    TFile products(Form("%s/overlay_histograms.root",outputDir),"RECREATE");
    SEME::require(!products.IsZombie(),"Cannot create output ROOT");
    std::ofstream norms(Form("%s/normalization.tsv",outputDir));
    std::ofstream bins(Form("%s/spectra.tsv",outputDir));
    norms << std::setprecision(17); bins << std::setprecision(17);
    norms << "group\tlabel\tSE_norm\tME_norm\tME_scale\tscaled_ME_norm\tSE_regular\tME_regular\tQA_max_abs_difference\n";
    bins << "group\tbinning\tbin\tlow\thigh\tSE\tME\tscaled_ME\tSE_error\tME_error\tscaled_ME_error\n";
    std::ifstream groups(groupFile); SEME::require(groups.good(),"Cannot read group configuration");
    std::string id,label; int centFirst,centLast; int pages=0;
    std::string allPDF=std::string(outputDir)+"/se_me_overlay.pdf";
    TCanvas pdfBook("se_me_pdf_book","",1600,850); pdfBook.Print((allPDF+"[").c_str());
    while (groups>>id>>label>>centFirst>>centLast) {
      int yFirst=sourceSE->GetYaxis()->FindBin(centFirst), yLast=sourceSE->GetYaxis()->FindBin(centLast);
      SEME::require(sourceME->GetYaxis()->FindBin(centFirst)==yFirst && sourceME->GetYaxis()->FindBin(centLast)==yLast,"Centrality axis mismatch");
      TH1D *se=sourceSE->ProjectionX((id+"_SE_native").c_str(),yFirst,yLast,"e");
      TH1D *me=sourceME->ProjectionX((id+"_ME_raw_native").c_str(),yFirst,yLast,"e");
      se->SetDirectory(0); me->SetDirectory(0); SEME::checkAxes(se,me);
      std::vector<double> qSE(se->GetNbinsX()+2,0.), qME(me->GetNbinsX()+2,0.);
      for (int cent=centFirst; cent<=centLast; ++cent) {
        TH2 *qaSE=SEME::get2(inputFile,Form(seQA,cent));
        TH2 *qaME=SEME::get2(inputFile,Form(meQA,cent));
        SEME::checkAxes(se,qaSE); SEME::checkAxes(me,qaME);
        for (int x=0; x<=se->GetNbinsX()+1; ++x)
          for (int y=0; y<=qaSE->GetNbinsY()+1; ++y) qSE[x]+=qaSE->GetBinContent(x,y);
        for (int x=0; x<=me->GetNbinsX()+1; ++x)
          for (int y=0; y<=qaME->GetNbinsY()+1; ++y) qME[x]+=qaME->GetBinContent(x,y);
      }
      double maxDiff=0.;
      for (int x=1; x<=se->GetNbinsX(); ++x) {
        maxDiff=std::max(maxDiff,std::fabs(qSE[x]-se->GetBinContent(x)));
        maxDiff=std::max(maxDiff,std::fabs(qME[x]-me->GetBinContent(x)));
      }
      SEME::require(maxDiff==0.,"Physical histogram differs from prior CF projection");
      int nf=SEME::edge(se->GetXaxis(),normLo), nl=SEME::edge(se->GetXaxis(),normHi)-1;
      double sn=se->Integral(nf,nl), mn=me->Integral(nf,nl);
      SEME::require(sn>0.&&mn>0.,"Nonpositive normalization integral");
      double scale=sn/mn;
      TH1D *scaled=static_cast<TH1D*>(me->Clone((id+"_ME_scaled_native").c_str())); scaled->Scale(scale);
      SEME::exportBins(bins,id,"native",se,me,scaled);
      products.cd(); se->Write(); me->Write(); scaled->Write();
      SEME::require(rebin>0 && se->GetNbinsX()%rebin==0,"Invalid display rebin factor");
      TH1D *s=static_cast<TH1D*>(se->Rebin(rebin,(id+"_SE_display").c_str()));
      TH1D *m=static_cast<TH1D*>(me->Rebin(rebin,(id+"_ME_raw_display").c_str()));
      TH1D *ms=static_cast<TH1D*>(scaled->Rebin(rebin,(id+"_ME_scaled_display").c_str()));
      int df=SEME::edge(s->GetXaxis(),normLo), dl=SEME::edge(s->GetXaxis(),normHi)-1;
      SEME::require(s->Integral(df,dl)==sn && m->Integral(df,dl)==mn,"Rebin changed normalization counts");
      SEME::require(std::fabs(ms->Integral(df,dl)-sn)<sn*1e-12,"Scaled normalization mismatch");
      norms << id << '\t' << label << '\t' << sn << '\t' << mn << '\t' << scale << '\t'
        << ms->Integral(df,dl) << '\t' << se->Integral() << '\t' << me->Integral() << '\t' << maxDiff << '\n';
      SEME::exportBins(bins,id,"display",s,m,ms);
      products.cd(); s->Write(); m->Write(); ms->Write();
      TParameter<double>((id+"_ME_scale").c_str(),scale).Write();
      TCanvas *c=new TCanvas(("se_me_"+id).c_str(),label.c_str(),1600,850);
      TPad *left=new TPad((id+"_full").c_str(),"",0.,.075,.5,.92);
      TPad *right=new TPad((id+"_zoom").c_str(),"",.5,.075,1.,.92);
      left->Draw(); right->Draw();
      left->cd(); SEME::panel(s,ms,fullLo,fullHi,normLo,normHi,scale);
      right->cd(); SEME::panel(s,ms,zoomLo,zoomHi,normLo,normHi,scale);
      c->cd(); TLatex title; title.SetNDC(); title.SetTextFont(42); title.SetTextSize(.035);
      title.DrawLatex(.085,.955,Form("#Lambda - ^{4}He  |  Centrality %s  |  k^{*} spectra before CF division",label.c_str()));
      title.SetTextSize(.023);
      title.DrawLatex(.085,.040,Form("#Lambda signal gate; selected ^{4}He; no background subtraction. Norm: SE %.0f / ME %.0f",sn,mn));
      title.DrawLatex(.085,.010,"2021 P24iy  |  1CE production (26945 jobs)  |  Left: full range; right: high-k^{*} detail  |  SE kept at raw counts");
      c->SaveAs(Form("%s/se_me_%s.png",outputDir,id.c_str()));
      c->SaveAs(Form("%s/se_me_%s.pdf",outputDir,id.c_str()));
      c->Print(allPDF.c_str()); products.cd(); c->Write(); ++pages;
      std::cout << "SEME_GROUP " << id << " SE_norm=" << sn << " ME_norm=" << mn << " scale=" << std::setprecision(17) << scale << std::endl;
    }
    SEME::require(pages>0,"No configured groups");
    pdfBook.Print((allPDF+"]").c_str());
    norms.close(); bins.close(); products.Close(); inputFile.Close();
    std::cout << "SEME_SUCCESS pages=" << pages << " root=" << gROOT->GetVersion() << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "SEME_FAILURE: " << e.what() << std::endl;
    gSystem->Exit(2);
  }
}
