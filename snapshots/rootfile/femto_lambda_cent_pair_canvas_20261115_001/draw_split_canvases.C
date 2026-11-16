void draw_split_canvases(const char* inputRoot, const char* planTsv,
                         const char* datasetsCsv, const char* outputPdf,
                         const char* outputRoot);

#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace split_canvas {
bool rootErrorSeen = false;
std::string rootErrorText;
void RootErrorHandler(Int_t level, Bool_t, const char* where, const char* message) {
  if (level >= kError) {
    rootErrorSeen = true;
    rootErrorText = std::string(where ? where : "ROOT") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, kFALSE, where, message);
}
struct RootErrorGuard {
  ErrorHandlerFunc_t previous;
  RootErrorGuard() : previous(SetErrorHandler(RootErrorHandler)) { rootErrorSeen=false; rootErrorText.clear(); }
  ~RootErrorGuard() { SetErrorHandler(previous); }
};
void Need(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
void CheckRoot() { Need(!rootErrorSeen, "ROOT error: " + rootErrorText); }
bool Finite(double value) { return TMath::Finite(value) != 0; }
std::vector<std::string> Split(const std::string& text, char sep) {
  std::vector<std::string> out; std::istringstream in(text); std::string field;
  while (std::getline(in, field, sep)) out.push_back(field); return out;
}
double Parameter(TDirectory* directory, const char* key) {
  TParameter<double>* value=dynamic_cast<TParameter<double>*>(directory->Get(key)); CheckRoot();
  Need(value != 0 && Finite(value->GetVal()), std::string("Missing/nonfinite parameter: ")+key);
  return value->GetVal();
}
void Label(double x, double y, const std::string& text, double size) {
  TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(size);
  label.DrawLatex(x,y,text.c_str());
}
void TitleLabel(double x, double y, const std::string& text) {
  const double nominal=.033;
  const double size=text.size()>72 ? nominal*72./text.size() : nominal;
  Label(x,y,text,std::max(.021,size));
}
std::string Number(double value) { std::ostringstream out; out.precision(6); out<<value; return out.str(); }
std::string CountNumber(double value) {
  std::ostringstream out;
  if (value==std::floor(value)) out<<std::fixed<<std::setprecision(0)<<value;
  else out<<std::setprecision(12)<<value;
  return out.str();
}
std::string IntegerText(int value) { std::ostringstream out; out<<value; return out.str(); }
void Style() {
  gPad->SetLeftMargin(.13); gPad->SetRightMargin(.045);
  gPad->SetTopMargin(.21); gPad->SetBottomMargin(.15); gPad->SetTicks(1,1);
}
void Axes(TH1* histogram) {
  histogram->SetStats(kFALSE);
  histogram->GetXaxis()->SetLabelFont(42); histogram->GetYaxis()->SetLabelFont(42);
  histogram->GetXaxis()->SetTitleFont(42); histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelSize(.038); histogram->GetYaxis()->SetLabelSize(.038);
  histogram->GetXaxis()->SetTitleSize(.044); histogram->GetYaxis()->SetTitleSize(.044);
  histogram->GetXaxis()->SetTitleOffset(1.1); histogram->GetYaxis()->SetTitleOffset(1.35);
}
void SameAxis(const TH1D& a,const TH1D& b) {
  Need(a.GetNbinsX()==b.GetNbinsX(),"Histogram x-bin counts differ");
  for (int i=1;i<=a.GetNbinsX()+1;++i)
    Need(a.GetXaxis()->GetBinLowEdge(i)==b.GetXaxis()->GetBinLowEdge(i),"Histogram x-bin edges differ");
}
TH1D* Histogram(TDirectory* directory,const char* key,bool required) {
  TH1D* original=dynamic_cast<TH1D*>(directory->Get(key)); CheckRoot();
  Need(original != 0 || !required,std::string("Missing histogram: ")+key);
  if (!original) return 0;
  TH1D* result=static_cast<TH1D*>(original->Clone()); result->SetDirectory(0); CheckRoot();
  for (int i=0;i<=result->GetNbinsX()+1;++i)
    Need(Finite(result->GetBinContent(i)) && Finite(result->GetBinError(i)),"Nonfinite histogram bin");
  return result;
}
void WindowBins(TH1D& histogram,double low,double high,int& first,int& last) {
  Need(Finite(low) && Finite(high) && low<high,"Invalid integration interval");
  const double magnitude=std::max(1.,std::max(std::fabs(low),std::fabs(high)));
  const double epsilon=64.*std::numeric_limits<double>::epsilon()*magnitude;
  first=std::max(1,histogram.GetXaxis()->FindBin(low+epsilon));
  last=std::min(histogram.GetNbinsX(),histogram.GetXaxis()->FindBin(high-epsilon));
  // All displayed counts are sums of existing bins; no fractional-bin weights
  // and no GetEntries/event-count interpretation are introduced here.
}
double WindowIntegral(TH1D& histogram,double low,double high) {
  int first=0,last=-1; WindowBins(histogram,low,high,first,last);
  const double count=first<=last ? histogram.Integral(first,last) : 0.;
  Need(Finite(count),"Nonfinite histogram integral"); return count;
}
int EmptyBins(TH1D& histogram,double low,double high) {
  int first=0,last=-1,count=0; WindowBins(histogram,low,high,first,last);
  for (int i=first;i<=last;++i) if (!(histogram.GetBinContent(i)>0.)) ++count;
  return count;
}
struct Data {
  std::string dataset;
  TH1D *se,*rawMe,*me;
  TGraphErrors* cf;
  double normLow,normHigh,displayLow,displayHigh,highLow,highHigh,scale;
  double totalSe,totalMe,normSe,normMe;
  bool normalized;
  Data(TFile& input,const std::string& ds,const std::string& cent,const std::string& selection)
      : dataset(ds),se(0),rawMe(0),me(0),cf(0),scale(0.),totalSe(0.),totalMe(0.),normSe(0.),normMe(0.),normalized(false) {
    TDirectory* directory=input.GetDirectory((ds+"/"+cent+"/"+selection).c_str()); CheckRoot();
    Need(directory != 0,"Missing selection directory: "+ds+"/"+cent+"/"+selection);
    se=Histogram(directory,"SE",true); rawMe=Histogram(directory,"ME",true); me=Histogram(directory,"ME_scaled",false);
    SameAxis(*se,*rawMe); if (me) SameAxis(*se,*me);
    TParameter<double>* norm=dynamic_cast<TParameter<double>*>(directory->Get("normA")); CheckRoot();
    if (norm && Finite(norm->GetVal())) scale=norm->GetVal();
    TGraphErrors* graph=dynamic_cast<TGraphErrors*>(directory->Get("CF_graph")); CheckRoot();
    Need(graph != 0,"Missing CF_graph"); cf=static_cast<TGraphErrors*>(graph->Clone());
    normLow=Parameter(directory,"normLow"); normHigh=Parameter(directory,"normHigh");
    displayLow=Parameter(directory,"displayLow"); displayHigh=Parameter(directory,"displayHigh");
    highLow=Parameter(directory,"highLow"); highHigh=Parameter(directory,"highHigh");
    Need(normLow<normHigh && displayLow<displayHigh && highLow<highHigh,"Reversed display/normalization interval");
    totalSe=se->Integral(1,se->GetNbinsX()); totalMe=rawMe->Integral(1,rawMe->GetNbinsX());
    normSe=WindowIntegral(*se,normLow,normHigh); normMe=WindowIntegral(*rawMe,normLow,normHigh);
    Need(Finite(totalSe) && Finite(totalMe),"Nonfinite normal-bin pair count");
    normalized=normSe>0. && normMe>0. && scale>0. && me!=0;
    for (int i=0;i<cf->GetN();++i) {
      double x=0.,y=0.; cf->GetPoint(i,x,y);
      Need(Finite(x) && Finite(y) && Finite(cf->GetErrorX(i)) && Finite(cf->GetErrorY(i)) &&
           cf->GetErrorX(i)>=0. && cf->GetErrorY(i)>=0.,"Invalid CF point/error");
    }
  }
  ~Data() { delete se; delete rawMe; delete me; delete cf; }
  std::string NormalizationStatus() const {
    if (!(normSe>0.) && !(normMe>0.)) return "CF unavailable: S_{H}=0 and M_{H}=0";
    if (!(normSe>0.)) return "CF unavailable: S_{H}=0";
    if (!(normMe>0.)) return "CF unavailable: M_{H}=0";
    if (!normalized) return "CF unavailable: scale/ME_{scaled} missing";
    return "";
  }
};
int ValidPoints(const Data& data,bool zoom) {
  if (!data.normalized) return 0;
  const double low=zoom?data.highLow:data.displayLow,high=zoom?data.highHigh:data.displayHigh;
  int count=0;
  for (int i=0;i<data.cf->GetN();++i) { double x=0.,y=0.; data.cf->GetPoint(i,x,y); if (x>=low && x<high) ++count; }
  return count;
}
int VisibleBins(Data& data,bool zoom) {
  int first=0,last=-1; WindowBins(*data.se,zoom?data.highLow:data.displayLow,zoom?data.highHigh:data.displayHigh,first,last);
  return std::max(0,last-first+1);
}
struct Range {
  double low,high; Range() : low(1.),high(1.) {}
  void Add(const Data& data,bool zoom) {
    if (!data.normalized) return;
    const double xmin=zoom?data.highLow:data.displayLow,xmax=zoom?data.highHigh:data.displayHigh;
    for (int i=0;i<data.cf->GetN();++i) {
      double x=0.,y=0.; data.cf->GetPoint(i,x,y); if (x<xmin || x>=xmax) continue;
      low=std::min(low,y-data.cf->GetErrorY(i)); high=std::max(high,y+data.cf->GetErrorY(i));
    }
  }
  void Pad() { double width=high-low; if (!(width>0.)) width=.2; low-=.12*width; high+=.25*width; }
};
void Heading(const Data& data,const std::string& title,const std::string& kind) {
  Label(.13,.960,"#Lambda-^{4}He signal; "+data.dataset,.034);
  TitleLabel(.13,.913,title);
  Label(.13,.855,kind+"; H: "+Number(data.normLow)+" #leq k^{*} < "+Number(data.normHigh)+" GeV/c",.027);
}
void Spectra(Data& data,const std::string& title) {
  Style(); TH1D* se=data.se; TH1D* me=data.normalized?data.me:0;
  double ymax=0.;
  for (int i=1;i<=se->GetNbinsX();++i) {
    if (se->GetXaxis()->GetBinUpEdge(i)<=data.displayLow || se->GetXaxis()->GetBinLowEdge(i)>=data.displayHigh) continue;
    ymax=std::max(ymax,se->GetBinContent(i)+se->GetBinError(i)); if (me) ymax=std::max(ymax,me->GetBinContent(i));
  }
  se->SetTitle(";k^{*} (GeV/c);Pair counts / bin"); Axes(se);
  se->GetXaxis()->SetRangeUser(data.displayLow,data.displayHigh);
  se->SetMinimum(0.); se->SetMaximum(ymax>0.?ymax*1.55:1.);
  se->SetLineColor(kBlack); se->SetMarkerColor(kBlack); se->SetMarkerStyle(20); se->SetMarkerSize(.55); se->Draw("E1");
  TLegend legend(.62,.625,.95,.785); legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextFont(42); legend.SetTextSize(.033);
  legend.AddEntry(se,"Same event","pe");
  if (me) { me->SetLineColor(kRed+1); me->SetLineWidth(2); me->Draw("HIST same"); legend.AddEntry(me,"Mixed event (scaled)","l"); }
  legend.DrawClone(); Heading(data,title,"SE + scaled ME");
  Label(.16,.755,"Pairs: S="+CountNumber(data.totalSe)+", M="+CountNumber(data.totalMe),.027);
  Label(.16,.700,"H counts: S_{H}="+CountNumber(data.normSe)+", M_{H}="+CountNumber(data.normMe),.027);
  if (me) Label(.16,.645,"ME scale A="+Number(data.scale),.027);
  else Label(.16,.645,data.NormalizationStatus(),.025);
  Label(.16,.590,"Pair-bin integrals; flow bins excluded",.024);
  gPad->RedrawAxis();
}
void Correlation(Data& data,const std::string& title,bool zoom,const Range& range,std::vector<TObject*>& owned) {
  Style(); const double xmin=zoom?data.highLow:data.displayLow,xmax=zoom?data.highHigh:data.displayHigh;
  TH1D* frame=new TH1D("frame",";k^{*} (GeV/c);C(k^{*})",1,xmin,xmax); frame->SetDirectory(0); owned.push_back(frame);
  Axes(frame); frame->SetMinimum(range.low); frame->SetMaximum(range.high); frame->Draw("AXIS");
  const double a=std::max(xmin,data.normLow),b=std::min(xmax,data.normHigh);
  if (a<b) { TBox box(a,range.low,b,range.high); box.SetFillColor(kGray); box.SetFillStyle(3001); box.DrawClone("same"); }
  TLine line(xmin,1.,xmax,1.); line.SetLineColor(kGray+2); line.SetLineStyle(2); line.DrawClone("same");
  data.cf->SetMarkerStyle(20); data.cf->SetMarkerSize(.65); data.cf->SetMarkerColor(kBlue+1); data.cf->SetLineColor(kBlue+1);
  const int valid=ValidPoints(data,zoom),visible=VisibleBins(data,zoom);
  if (valid) data.cf->Draw("P same");
  else Label(.17,.50,data.normalized?"No valid CF points in this range":data.NormalizationStatus(),.032);
  Heading(data,title,zoom?"CF: high k^{*}":"CF: full range");
  Label(.16,.755,"Valid CF bins: "+IntegerText(valid)+" / "+IntegerText(visible),.027);
  const double low=zoom?data.highLow:data.displayLow,high=zoom?data.highHigh:data.displayHigh;
  Label(.16,.704,"Empty raw bins: S="+IntegerText(EmptyBins(*data.se,low,high))+", M="+IntegerText(EmptyBins(*data.rawMe,low,high)),.027);
  gPad->RedrawAxis();
}
void Run(const char* inputRoot,const char* planTsv,const char* datasetsCsv,const char* outputPdf,const char* outputRoot) {
  RootErrorGuard errors;
  Need(gSystem->AccessPathName(outputPdf) && gSystem->AccessPathName(outputRoot),"Output already exists");
  gROOT->SetBatch(kTRUE); const Bool_t oldAddDirectory=TH1::AddDirectoryStatus(); TH1::AddDirectory(kFALSE);
  struct AddDirectoryGuard { Bool_t previous; explicit AddDirectoryGuard(Bool_t value):previous(value){} ~AddDirectoryGuard(){TH1::AddDirectory(previous);} } addDirectoryGuard(oldAddDirectory);
  gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
  TFile input(inputRoot,"READ"); CheckRoot(); Need(!input.IsZombie() && !input.TestBit(TFile::kRecovered),"Cannot open clean input ROOT");
  TFile output(outputRoot,"CREATE"); CheckRoot(); Need(!output.IsZombie() && output.IsWritable(),"Cannot create output ROOT");
  const std::vector<std::string> datasets=Split(datasetsCsv,','); Need(!datasets.empty(),"No datasets supplied");
  std::ifstream plan(planTsv); Need(bool(plan),"Cannot read page plan"); std::string line; int pages=0; std::set<std::string> canvasKeys;
  TCanvas canvas("split_canvas","SE/ME and CF",1600,1450); canvas.Print((std::string(outputPdf)+"[").c_str()); CheckRoot();
  while (std::getline(plan,line)) {
    if (!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
    if (line.empty()) continue;
    const std::vector<std::string> fields=Split(line,'\t'); Need(fields.size()==4,"Invalid page row");
    const std::string canvasKey=fields[0]+"__"+fields[1]; Need(canvasKeys.insert(canvasKey).second,"Duplicate cent__selection Canvas key");
    std::vector<Data*> data; std::vector<TObject*> owned; Range full,high;
    for (size_t i=0;i<datasets.size();++i) {
      data.push_back(new Data(input,datasets[i],fields[0],fields[1]));
      full.Add(*data.back(),false); high.Add(*data.back(),true);
      std::cout<<"COLUMN_COUNTS\t"<<fields[0]<<'\t'<<fields[1]<<'\t'<<datasets[i]
               <<"\tSE="<<CountNumber(data.back()->totalSe)<<"\tME="<<CountNumber(data.back()->totalMe)
               <<"\tSH="<<CountNumber(data.back()->normSe)<<"\tMH="<<CountNumber(data.back()->normMe)
               <<"\tnorm="<<(data.back()->normalized?"valid":"unavailable")
               <<"\tvalidFull="<<ValidPoints(*data.back(),false)<<"\tvalidHigh="<<ValidPoints(*data.back(),true)<<std::endl;
    }
    full.Pad(); high.Pad(); canvas.Clear(); canvas.SetName(canvasKey.c_str()); canvas.Divide(datasets.size(),3,.002,.002);
    for (size_t i=0;i<data.size();++i) {
      canvas.cd(i+1); Spectra(*data[i],fields[2]);
      canvas.cd(datasets.size()+i+1); Correlation(*data[i],fields[2],false,full,owned);
      canvas.cd(2*datasets.size()+i+1); Correlation(*data[i],fields[2],true,high,owned);
    }
    const std::string bookmark="Title:"+fields[0]+" / "+fields[1];
    canvas.Print(outputPdf,bookmark.c_str()); CheckRoot();
    output.cd(); Need(canvas.Write(canvasKey.c_str())>0,"Cannot save Canvas"); CheckRoot(); ++pages;
    canvas.Clear(); for (size_t i=0;i<owned.size();++i) delete owned[i]; for (size_t i=0;i<data.size();++i) delete data[i];
  }
  Need(pages>0 && !plan.bad(),"Page plan is empty/unreadable"); canvas.Print((std::string(outputPdf)+"]").c_str()); CheckRoot();
  output.Close(); input.Close(); CheckRoot(); std::cout<<"DRAW_SPLIT_CANVASES_OK pages="<<pages<<std::endl;
}
}
void draw_split_canvases(const char* inputRoot,const char* planTsv,const char* datasetsCsv,const char* outputPdf,const char* outputRoot) {
  try { split_canvas::Run(inputRoot,planTsv,datasetsCsv,outputPdf,outputRoot); }
  catch (const std::exception& error) { std::cerr<<"DRAW_SPLIT_CANVASES_FAILED: "<<error.what()<<std::endl; gSystem->Exit(2); }
  catch (...) { std::cerr<<"DRAW_SPLIT_CANVASES_FAILED: unknown exception"<<std::endl; gSystem->Exit(2); }
}
#endif
