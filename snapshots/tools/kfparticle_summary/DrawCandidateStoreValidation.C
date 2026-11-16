#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <iostream>
#include <stdexcept>
#include <vector>
static TH1* StoreValidationHist(TFile& f,const char* key,const char* name,int rebin){TH1* src=dynamic_cast<TH1*>(f.Get(key));if(!src)throw std::runtime_error(key);TH1* h=static_cast<TH1*>(src->Clone(name));h->SetDirectory(0);if(rebin<1||h->GetNbinsX()%rebin)throw std::runtime_error("invalid plotting rebin");h->Rebin(rebin);h->SetTitle("");h->GetXaxis()->SetTitle("k* (GeV/c)");h->GetYaxis()->SetTitle("Pairs / displayed bin");return h;}
static void StoreValidationAxes(TH1* h){h->GetXaxis()->SetTitleSize(.05);h->GetYaxis()->SetTitleSize(.05);h->GetXaxis()->SetLabelSize(.043);h->GetYaxis()->SetLabelSize(.043);h->GetYaxis()->SetTitleOffset(1.25);}
void DrawCandidateStoreValidation(const char* baseline,const char* replay,const char* variant,const char* output,int rebin){
 try{if(!gSystem->AccessPathName(output))throw std::runtime_error("refusing existing PDF");TFile a(baseline,"READ"),b(replay,"READ"),c(variant,"READ");if(a.IsZombie()||b.IsZombie()||c.IsZombie())throw std::runtime_error("invalid ROOT input");gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(26,20);TCanvas canvas("candidate_cache_validation","Candidate cache validation",1400,1000);canvas.Divide(2,2,.003,.003);
  const char* keys[]={"hKstarSE_lambda_he4_signal","hKstarME_lambda_he4_signal"};
  for(int i=0;i<3;++i){canvas.cd(i+1);gPad->SetTopMargin(.13);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);TH1* x=StoreValidationHist(a,keys[i==0?0:1],Form("cv_a_%d",i),rebin);TH1* y=StoreValidationHist(i<2?b:c,keys[i==0?0:1],Form("cv_b_%d",i),rebin);StoreValidationAxes(x);x->SetMinimum(0);x->SetMaximum(1.5*std::max(x->GetMaximum(),y->GetMaximum())+.5);x->SetLineColor(kBlack);x->SetMarkerColor(kBlack);x->SetMarkerStyle(20);y->SetMarkerStyle(24);y->SetMarkerSize(1.1);y->SetLineColor(i<2?kBlue+1:kRed+1);y->SetMarkerColor(i<2?kBlue+1:kRed+1);x->Draw("E");y->Draw("E SAME");TLegend* l=new TLegend(.43,.66,.94,.84);l->SetFillStyle(0);l->SetBorderSize(0);l->SetTextSize(.036);l->AddEntry(x,Form("Direct: %.0f entries",x->GetEntries()),"lp");l->AddEntry(y,Form(i<2?"Cached replay: %.0f":"Fine Vz replay: %.0f",y->GetEntries()),"lp");l->Draw();TLatex t;t.SetNDC();t.SetTextSize(.047);t.DrawLatex(.14,.94,i==0?"SE: direct vs cached replay":(i==1?"ME: direct vs cached replay":"ME: Vz bins 10 vs 20"));}
  canvas.cd(4);gPad->SetTopMargin(.13);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);TH1* x=StoreValidationHist(a,keys[1],"cv_ref_difference",rebin);TH1* y=StoreValidationHist(c,keys[1],"cv_difference",rebin);for(int i=0;i<=x->GetNbinsX()+1;++i){y->SetBinContent(i,y->GetBinContent(i)-x->GetBinContent(i));y->SetBinError(i,0);}StoreValidationAxes(y);y->GetYaxis()->SetTitle("ME(fine Vz) - ME(original)");y->SetMinimum(std::min(-1.,1.3*y->GetMinimum()));y->SetMaximum(std::max(1.,1.3*y->GetMaximum()));y->SetLineColor(kRed+1);y->SetFillColor(kRed-9);y->Draw("HIST");TLine zero(y->GetXaxis()->GetXmin(),0,y->GetXaxis()->GetXmax(),0);zero.SetLineStyle(2);zero.Draw();TLatex t;t.SetNDC();t.SetTextSize(.045);t.DrawLatex(.14,.94,"Fixed candidates; only ME changes");t.SetTextSize(.035);t.DrawLatex(.16,.82,"Native-bin closure is checked separately.");t.DrawLatex(.16,.77,"Display rebin only; no physics normalization.");canvas.Print(output);
 }catch(const std::exception& e){std::cerr<<"CANDIDATE_STORE_PLOT_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
