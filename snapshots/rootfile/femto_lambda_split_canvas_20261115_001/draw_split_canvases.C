void draw_split_canvases(const char* inputRoot, const char* planTsv,
                         const char* datasetsCsv, const char* outputPdf,
                         const char* outputRoot);

#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace split_canvas {
void Need(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
std::vector<std::string> Split(const std::string& text, char sep) {
  std::vector<std::string> out; std::istringstream in(text); std::string field;
  while (std::getline(in, field, sep)) out.push_back(field); return out;
}
double Parameter(TDirectory* directory, const char* key) {
  TParameter<double>* value = dynamic_cast<TParameter<double>*>(directory->Get(key));
  Need(value != 0, std::string("Missing parameter: ") + key); return value->GetVal();
}
void Label(double x, double y, const std::string& text, double size) {
  TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(size);
  label.DrawLatex(x, y, text.c_str());
}
std::string Number(double value) { std::ostringstream out; out.precision(5); out << value; return out.str(); }
void Style() {
  gPad->SetLeftMargin(.13); gPad->SetRightMargin(.045);
  gPad->SetTopMargin(.18); gPad->SetBottomMargin(.15); gPad->SetTicks(1,1);
}
void Axes(TH1* histogram) {
  histogram->SetStats(kFALSE);
  histogram->GetXaxis()->SetLabelFont(42); histogram->GetYaxis()->SetLabelFont(42);
  histogram->GetXaxis()->SetTitleFont(42); histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelSize(.038); histogram->GetYaxis()->SetLabelSize(.038);
  histogram->GetXaxis()->SetTitleSize(.044); histogram->GetYaxis()->SetTitleSize(.044);
  histogram->GetXaxis()->SetTitleOffset(1.1); histogram->GetYaxis()->SetTitleOffset(1.35);
}
struct Data {
  std::string dataset, note; TH1D *se, *me; TGraphErrors* cf;
  double normLow, normHigh, displayLow, displayHigh, highLow, highHigh, scale;
  Data(TFile& input, const std::string& ds, const std::string& cent, const std::string& selection) : dataset(ds), se(0), me(0), cf(0), scale(0.) {
    TDirectory* directory = input.GetDirectory((ds + "/" + cent + "/" + selection).c_str());
    Need(directory != 0, "Missing selection directory");
    TH1D* original = dynamic_cast<TH1D*>(directory->Get("SE"));
    Need(original != 0, "Missing SE"); se = static_cast<TH1D*>(original->Clone()); se->SetDirectory(0);
    original = dynamic_cast<TH1D*>(directory->Get("ME_scaled"));
    if (original) { me = static_cast<TH1D*>(original->Clone()); me->SetDirectory(0); scale=Parameter(directory,"normA"); }
    TGraphErrors* graph = dynamic_cast<TGraphErrors*>(directory->Get("CF_graph"));
    Need(graph != 0, "Missing CF_graph"); cf = static_cast<TGraphErrors*>(graph->Clone());
    normLow=Parameter(directory,"normLow"); normHigh=Parameter(directory,"normHigh");
    displayLow=Parameter(directory,"displayLow"); displayHigh=Parameter(directory,"displayHigh");
    highLow=Parameter(directory,"highLow"); highHigh=Parameter(directory,"highHigh");
  }
  ~Data() { delete se; delete me; delete cf; }
};
struct Range {
  double low, high; Range() : low(1.), high(1.) {}
  void Add(const Data& data, bool zoom) {
    const double xmin=zoom?data.highLow:data.displayLow, xmax=zoom?data.highHigh:data.displayHigh;
    for (int i=0; i<data.cf->GetN(); ++i) {
      double x,y; data.cf->GetPoint(i,x,y); if (x<xmin || x>=xmax) continue;
      low=std::min(low,y-data.cf->GetErrorY(i)); high=std::max(high,y+data.cf->GetErrorY(i));
    }
  }
  void Pad() { double width=high-low; if (!(width>0.)) width=.2; low-=.12*width; high+=.18*width; }
};
void Heading(const Data& data, const std::string& title, const std::string& kind) {
  Label(.13,.956,"#Lambda-^{4}He signal; " + title + "; " + data.dataset,.034);
  Label(.13,.89,kind + "; norm: " + Number(data.normLow) + " #leq k^{*} < " + Number(data.normHigh) + " GeV/c",.030);
}
void Spectra(Data& data, const std::string& title) {
  Style(); TH1D* se=data.se; TH1D* me=data.me;
  double ymax=0.;
  for (int i=1; i<=se->GetNbinsX(); ++i) {
    ymax=std::max(ymax,se->GetBinContent(i)+se->GetBinError(i));
    if (me) ymax=std::max(ymax,me->GetBinContent(i));
  }
  se->SetTitle(";k^{*} (GeV/c);Pair counts / bin"); Axes(se);
  se->GetXaxis()->SetRangeUser(data.displayLow,data.displayHigh);
  se->SetMinimum(0.); se->SetMaximum(ymax>0.?ymax*1.25:1.);
  se->SetLineColor(kBlack); se->SetMarkerColor(kBlack); se->SetMarkerStyle(20); se->SetMarkerSize(.55);
  se->Draw("E1");
  TLegend legend(.59,.66,.95,.81); legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextFont(42); legend.SetTextSize(.036);
  legend.AddEntry(se,"Same event","pe");
  if (me) { me->SetLineColor(kRed+1); me->SetLineWidth(2); me->Draw("HIST same"); legend.AddEntry(me,"Mixed event (scaled)","l"); }
  legend.DrawClone(); Heading(data,title,"SE + scaled ME");
  if (me) Label(.16,.75,"ME scale = " + Number(data.scale),.029);
  else Label(.16,.75,"No SE normalization counts; scale unavailable",.029);
  if (data.note!="-") Label(.16,.67,data.note,.029);
  gPad->RedrawAxis();
}
void Correlation(Data& data, const std::string& title, bool zoom, const Range& range, std::vector<TObject*>& owned) {
  Style(); double xmin=zoom?data.highLow:data.displayLow, xmax=zoom?data.highHigh:data.displayHigh;
  TH1D* frame = new TH1D("frame",";k^{*} (GeV/c);C(k^{*})",1,xmin,xmax); frame->SetDirectory(0); owned.push_back(frame);
  Axes(frame); frame->SetMinimum(range.low); frame->SetMaximum(range.high); frame->Draw("AXIS");
  double a=std::max(xmin,data.normLow), b=std::min(xmax,data.normHigh);
  if (a<b) { TBox box(a,range.low,b,range.high); box.SetFillColor(kGray); box.SetFillStyle(3001); box.DrawClone("same"); }
  TLine line(xmin,1.,xmax,1.); line.SetLineColor(kGray+2); line.SetLineStyle(2); line.DrawClone("same");
  data.cf->SetMarkerStyle(20); data.cf->SetMarkerSize(.65); data.cf->SetMarkerColor(kBlue+1); data.cf->SetLineColor(kBlue+1);
  if (data.cf->GetN()) data.cf->Draw("P same");
  else Label(.19,.53,"CF unavailable: no normalization SE counts",.036);
  Heading(data,title,zoom?"CF: high k^{*}":"CF: full range");
  Label(.16,.75,"Raw CF; shared-normalization statistical errors",.028); gPad->RedrawAxis();
  if (data.note!="-") Label(.16,.68,data.note,.029);
}
void Run(const char* inputRoot, const char* planTsv, const char* datasetsCsv, const char* outputPdf, const char* outputRoot) {
  Need(gSystem->AccessPathName(outputPdf) && gSystem->AccessPathName(outputRoot),"Output already exists");
  gROOT->SetBatch(kTRUE); TH1::AddDirectory(kFALSE); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
  TFile input(inputRoot,"READ"); Need(!input.IsZombie(),"Cannot open input ROOT");
  TFile output(outputRoot,"CREATE"); Need(!output.IsZombie(),"Cannot create output ROOT");
  std::vector<std::string> datasets=Split(datasetsCsv,','); Need(!datasets.empty(),"No datasets supplied");
  std::ifstream plan(planTsv); Need(bool(plan),"Cannot read page plan"); std::string line; int pages=0;
  TCanvas canvas("split_canvas","SE/ME and CF",1600,1450); canvas.Print((std::string(outputPdf)+"[").c_str());
  while (std::getline(plan,line)) {
    if (line.empty()) continue; std::vector<std::string> fields=Split(line,'\t'); Need(fields.size()==4,"Invalid page row");
    std::vector<Data*> data; std::vector<TObject*> owned; Range full, high;
    for (size_t i=0; i<datasets.size(); ++i) { data.push_back(new Data(input,datasets[i],fields[0],fields[1])); data.back()->note=fields[3]; full.Add(*data.back(),false); high.Add(*data.back(),true); }
    full.Pad(); high.Pad(); canvas.Clear(); canvas.Divide(datasets.size(),3,.002,.002);
    for (size_t i=0; i<data.size(); ++i) {
      canvas.cd(i+1); Spectra(*data[i],fields[2]);
      canvas.cd(datasets.size()+i+1); Correlation(*data[i],fields[2],false,full,owned);
      canvas.cd(2*datasets.size()+i+1); Correlation(*data[i],fields[2],true,high,owned);
    }
    canvas.Print(outputPdf); output.cd(); Need(canvas.Write(fields[1].c_str())>0,"Cannot save Canvas"); ++pages;
    canvas.Clear(); for (size_t i=0;i<owned.size();++i) delete owned[i]; for (size_t i=0;i<data.size();++i) delete data[i];
  }
  Need(pages>0 && !plan.bad(),"Page plan is empty/unreadable"); canvas.Print((std::string(outputPdf)+"]").c_str());
  output.Close(); input.Close(); std::cout<<"DRAW_SPLIT_CANVASES_OK pages="<<pages<<std::endl;
}
}
void draw_split_canvases(const char* inputRoot, const char* planTsv, const char* datasetsCsv, const char* outputPdf, const char* outputRoot) {
  try { split_canvas::Run(inputRoot,planTsv,datasetsCsv,outputPdf,outputRoot); }
  catch (const std::exception& error) { std::cerr<<"DRAW_SPLIT_CANVASES_FAILED: "<<error.what()<<std::endl; gSystem->Exit(2); }
}
#endif
