void DrawHe4Y(const char* metadataPath,const char* binsPath,const char* settingsPath,const char* outputDirectory);
#ifndef __CINT__
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TNamed.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TError.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PeriodCFDraw {
bool rootError=false;
void Handler(Int_t level,Bool_t abort,const char* where,const char* text) {if(level>=kError)rootError=true;DefaultErrorHandler(level,abort,where,text);}
void Need(bool ok,const std::string& message) {if(!ok)throw std::runtime_error(message);}
std::vector<std::string> Fields(const std::string& text) {
  std::vector<std::string> fields;size_t start=0;
  for(;;){size_t end=text.find('\t',start);fields.push_back(text.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}return fields;
}
double Number(const std::string& value) {std::istringstream s(value);double n=0;Need(bool(s>>n)&&TMath::Finite(n),"invalid numeric table cell");s>>std::ws;Need(s.eof(),"trailing numeric table cell");return n;}
std::string Key(const std::string& id,const std::string& cent,int plane) {std::ostringstream s;s<<id<<'/'<<cent<<'/'<<plane;return s.str();}
struct Bin {double low,high,se,sv,me,mv,cf,error;bool valid;};
struct Curve {
  std::string id,cent,label;int plane,ep,cut;double alpha,sh,mh;bool valid;
  std::vector<Bin> bins;TGraphErrors* cf;TGraphErrors* se;TGraphErrors* me;
  Curve():plane(-2),ep(0),cut(0),alpha(0),sh(0),mh(0),valid(false),cf(0),se(0),me(0){}
};
typedef std::map<std::string,Curve> Curves;
struct Settings {double low,high,cfLow,cfHigh,normLow,normHigh,highKLow,highKHigh,highKCFMin,highKCFMax;bool fixedCFRange;std::string mode;std::vector<std::string> centralities,diagnostics;Settings():low(0),high(0),cfLow(0),cfHigh(0),normLow(0),normHigh(0),highKLow(0),highKHigh(0),highKCFMin(0),highKCFMax(0),fixedCFRange(false){};};
Curve& Find(Curves& curves,const std::string& cent,int ep,int cut,int plane) {
  for(Curves::iterator it=curves.begin();it!=curves.end();++it)if(it->second.cent==cent&&it->second.ep==ep&&it->second.cut==cut&&it->second.plane==plane)return it->second;
  throw std::runtime_error("missing requested plot curve");
}
void Text(double x,double y,const std::string& value,double size=.034) {TLatex* text=new TLatex(x,y,value.c_str());text->SetNDC();text->SetTextFont(42);text->SetTextSize(size);text->Draw();}
double Maximum(const TGraphErrors* graph,const Settings& settings) {
  double result=0;for(int i=0;i<graph->GetN();++i){double x=0,y=0;graph->GetPoint(i,x,y);if(x>=settings.low&&x<=settings.high)result=std::max(result,y+graph->GetErrorY(i));}return result;
}
void Frame(const std::string& title,const Settings& settings,double low,double high,const char* yTitle) {
  gPad->SetLeftMargin(.13);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.16);
  TH1* frame=gPad->DrawFrame(settings.low,low,settings.high,high);frame->SetTitle(title.c_str());frame->GetXaxis()->SetTitle("k* [GeV/c]");frame->GetYaxis()->SetTitle(yTitle);frame->GetXaxis()->SetTitleSize(.045);frame->GetYaxis()->SetTitleSize(.045);frame->GetXaxis()->SetLabelSize(.04);frame->GetYaxis()->SetLabelSize(.04);frame->GetYaxis()->SetTitleOffset(1.2);
}
TLegend* Legend() {TLegend* legend=new TLegend(.57,.71,.95,.84);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.032);return legend;}
void CFPanel(const std::string& title,const std::vector<Curve*>& sources,const Settings& settings,bool plane) {
  double high=settings.cfHigh;if(!settings.fixedCFRange)for(size_t i=0;i<sources.size();++i)high=std::max(high,1.1*Maximum(sources[i]->cf,settings));
  Frame(title,settings,settings.cfLow,high,"C(k*) raw signal");
  TLine* line=new TLine(settings.low,1,settings.high,1);line->SetLineColor(kGray+2);line->SetLineStyle(2);line->Draw();
  TLegend* legend=Legend();bool available=false;
  for(size_t i=0;i<sources.size();++i) {
    Curve& c=*sources[i];TGraphErrors* graph=c.cf;
    graph->SetLineColor(plane?(i==0?kBlue+1:i==1?kBlack:kGreen+2):(c.cut?kRed+1:kBlue+1));graph->SetMarkerColor(graph->GetLineColor());graph->SetLineStyle(plane&&i?(i==1?2:3):1);graph->SetMarkerStyle(plane&&i?24:(c.cut?24:20));graph->SetMarkerSize(.45);if(graph->GetN()>0)graph->Draw("P E1 SAME");available=available||graph->GetN()>0;
    std::string label=plane?(i?(i==1?"all Y inclusive":"sum in-range Y bins"):Form("He4 Y bin %d",c.plane)):(c.plane==-1?"all Y inclusive":c.plane==-2?"sum in-range Y bins":c.plane==-3?"Y underflow control":"Y overflow control");legend->AddEntry(graph,label.c_str(),"lep");
  }
  legend->Draw();if(!sources[0]->valid)Text(.18,.48,"Y slice: normalization unavailable",.038);else if(!available)Text(.22,.48,"No valid normalization/count bins",.038);
  Text(.14,.89,Form("EP8, cut Off; He4 y_{lab}; H=[%.2f,%.2f); fixed-alpha errors",settings.normLow,settings.normHigh),.028);
  Text(.14,.02,"Raw counts summed before normalization; no background subtraction",.024);
}
void RawPanel(const std::string& title,const std::vector<Curve*>& sources,const Settings& settings) {
  double high=0;for(size_t i=0;i<sources.size();++i){high=std::max(high,Maximum(sources[i]->se,settings));high=std::max(high,Maximum(sources[i]->me,settings));}
  Frame(title,settings,0,high>0?high*1.35:1,"Pairs/bin divided by H integral");
  TLegend* legend=new TLegend(.52,.65,.96,.84);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.03);
  for(size_t i=0;i<sources.size();++i) {
    Curve& c=*sources[i];int color=c.cut?kRed+1:kBlue+1;
    c.se->SetLineColor(color);c.se->SetMarkerColor(color);c.se->SetMarkerStyle(c.cut?24:20);c.se->SetMarkerSize(.4);c.se->SetLineStyle(1);if(c.se->GetN()>0)c.se->Draw("LP SAME");
    c.me->SetLineColor(kRed+1);c.me->SetMarkerColor(kRed+1);c.me->SetLineStyle(2);if(c.me->GetN()>0)c.me->Draw("L SAME");
    legend->AddEntry(c.se,(c.cut?"SE / SE(H)":"SE / SE(H)"),"lp");legend->AddEntry(c.me,(c.cut?"ME / ME(H)":"ME / ME(H)"),"l");
  }
  legend->Draw();if(sources[0]->sh<=0||sources[0]->mh<=0)Text(.18,.48,"Y slice: normalization unavailable",.038);Text(.14,.89,Form("EP8, cut Off; He4 y_{lab}; H=[%.2f,%.2f); shape scales",settings.normLow,settings.normHigh),.028);Text(.14,.02,"Unscaled SE/ME counts and Y-flow populations are in he4_y_cf_metrics.json",.024);
}
void Pages(Curves& curves,const Settings& settings,const std::string& output) {
 TCanvas canvas("he4_y_canvas","Single He4 lab rapidity",1600,1100);
 const int groups[]={4,8};
 for(int ig=0;ig<2;++ig){int n=groups[ig];std::string id=Form("he4Y%d",n),cf=output+Form("/he4_y%d_CF.pdf",n),raw=output+Form("/he4_y%d_SE_ME_scaled.pdf",n),highk=output+Form("/he4_y%d_CF_highk.pdf",n);canvas.Print((cf+"[").c_str());canvas.Print((raw+"[").c_str());canvas.Print((highk+"[").c_str());
  for(size_t ic=0;ic<settings.centralities.size();++ic){const std::string cent=settings.centralities[ic];
   canvas.Clear();canvas.Divide(n==4?2:4,2);for(int k=0;k<n;++k){canvas.cd(k+1);Curve& c=Find(curves,cent,n,0,k);std::vector<Curve*> sources;sources.push_back(&c);sources.push_back(&Find(curves,cent,n,0,-1));sources.push_back(&Find(curves,cent,n,0,-2));CFPanel(c.label,sources,settings,true);}canvas.Print(cf.c_str());
   Settings zoom=settings;zoom.low=settings.highKLow;zoom.high=settings.highKHigh;zoom.cfLow=settings.highKCFMin;zoom.cfHigh=settings.highKCFMax;zoom.fixedCFRange=true;canvas.Clear();canvas.Divide(n==4?2:4,2);for(int k=0;k<n;++k){canvas.cd(k+1);Curve& c=Find(curves,cent,n,0,k);std::vector<Curve*> sources;sources.push_back(&c);sources.push_back(&Find(curves,cent,n,0,-1));sources.push_back(&Find(curves,cent,n,0,-2));CFPanel(c.label,sources,zoom,true);}canvas.Print(highk.c_str());
   canvas.Clear();canvas.Divide(n==4?2:4,2);for(int k=0;k<n;++k){canvas.cd(k+1);Curve& c=Find(curves,cent,n,0,k);std::vector<Curve*> one(1,&c);RawPanel(c.label,one,settings);}canvas.Print(raw.c_str());
  }canvas.Print((cf+"]").c_str());canvas.Print((raw+"]").c_str());canvas.Print((highk+"]").c_str());
 }
 std::string flow=output+"/he4_y_flow_controls.pdf";canvas.Print((flow+"[").c_str());
 for(size_t ic=0;ic<settings.diagnostics.size();++ic){canvas.Clear();canvas.Divide(2,2);for(int k=0;k<4;++k){canvas.cd(k+1);Curve& c=Find(curves,settings.diagnostics[ic],8,0,-1-k);std::vector<Curve*> one(1,&c);CFPanel(c.label,one,settings,false);}canvas.Print(flow.c_str());canvas.Clear();canvas.Divide(2,2);for(int k=0;k<4;++k){canvas.cd(k+1);Curve& c=Find(curves,settings.diagnostics[ic],8,0,-1-k);std::vector<Curve*> one(1,&c);RawPanel(c.label,one,settings);}canvas.Print(flow.c_str());}canvas.Print((flow+"]").c_str());
}
} // namespace PeriodCFDraw

void DrawHe4Y(const char* metadataPath,const char* binsPath,const char* settingsPath,const char* outputDirectory) {
  using namespace PeriodCFDraw;
  ErrorHandlerFunc_t previous=SetErrorHandler(Handler);rootError=false;TH1::AddDirectory(kFALSE);
  try {
    Need(metadataPath&&binsPath&&settingsPath&&outputDirectory,"explicit plot arguments required");
    Curves curves;std::ifstream metadata(metadataPath);Need(bool(metadata),"cannot read curve metadata");std::string line;
    while(std::getline(metadata,line)) {
      std::vector<std::string> f=Fields(line);Need(f.size()==9,"metadata requires nine columns");Curve c;c.id=f[0];c.cent=f[1];c.label=f[2];c.plane=int(Number(f[3]));c.ep=int(Number(f[4]));c.cut=int(Number(f[5]));c.valid=f[6]!="NA";if(c.valid)c.alpha=Number(f[6]);c.sh=Number(f[7]);c.mh=Number(f[8]);Need(curves.insert(std::make_pair(Key(c.id,c.cent,c.plane),c)).second,"duplicate curve metadata");
    }
    Need(metadata.eof()&&!curves.empty(),"empty/failed metadata");
    std::ifstream bins(binsPath);Need(bool(bins),"cannot read curve bins");
    while(std::getline(bins,line)) {
      std::vector<std::string> f=Fields(line);Need(f.size()==12,"bin requires twelve columns");Curves::iterator it=curves.find(Key(f[0],f[1],int(Number(f[2]))));Need(it!=curves.end(),"bin has no metadata");Need(int(Number(f[3]))==int(it->second.bins.size()),"bin order mismatch");Bin b;b.low=Number(f[4]);b.high=Number(f[5]);b.se=Number(f[6]);b.sv=Number(f[7]);b.me=Number(f[8]);b.mv=Number(f[9]);b.valid=f[10]!="NA";Need(b.valid==(f[11]!="NA"),"inconsistent CF/error availability");b.cf=b.error=0;if(b.valid){b.cf=Number(f[10]);b.error=Number(f[11]);}Need(b.high>b.low&&b.se>=0&&b.me>=0&&b.sv>=0&&b.mv>=0&&b.error>=0,"invalid count/bin");if(!it->second.bins.empty())Need(it->second.bins.back().high==b.low,"noncontiguous bins");it->second.bins.push_back(b);
    }
    Need(bins.eof(),"bin input read failure");Settings settings;std::ifstream inputSettings(settingsPath);Need(bool(inputSettings),"cannot read settings");std::set<std::string> scalarKeys;
    while(std::getline(inputSettings,line)) {
      std::vector<std::string> f=Fields(line);Need(f.size()>=2,"invalid plot settings");
      if(f[0]=="centrality")settings.centralities.push_back(f[1]);
      else if(f[0]=="diagnostic")settings.diagnostics.push_back(f[1]);
      else {Need(scalarKeys.insert(f[0]).second,"duplicate plot scalar");if(f[0]=="mode"){Need(f.size()==2,"bad mode row");settings.mode=f[1];}else {Need(f.size()==3,"invalid range row");double a=Number(f[1]),b=Number(f[2]);Need(a<b,"invalid plot range");if(f[0]=="kRange"){settings.low=a;settings.high=b;}else if(f[0]=="cfRange"){settings.cfLow=a;settings.cfHigh=b;}else if(f[0]=="normRange"){settings.normLow=a;settings.normHigh=b;}else if(f[0]=="highKRange"){settings.highKLow=a;settings.highKHigh=b;}else if(f[0]=="highKCFRange"){settings.highKCFMin=a;settings.highKCFMax=b;}else Need(false,"unknown plot setting");}}
    }
    Need(inputSettings.eof()&&scalarKeys.size()==6&&!settings.centralities.empty()&&!settings.diagnostics.empty()&&settings.mode=="he4Y","incomplete plot settings");
    std::string rootPath=std::string(outputDirectory)+"/cf_histograms.root";Need(gSystem->AccessPathName(rootPath.c_str()),"refusing existing companion ROOT");TFile* output=TFile::Open(rootPath.c_str(),"CREATE");Need(output&&!output->IsZombie(),"cannot create companion ROOT");
    for(Curves::iterator it=curves.begin();it!=curves.end();++it) {
      Curve& c=it->second;Need(!c.bins.empty(),"curve has no bins");std::vector<double> edges;for(size_t i=0;i<c.bins.size();++i)edges.push_back(c.bins[i].low);edges.push_back(c.bins.back().high);
      output->cd();TDirectory* caseDir=output->GetDirectory(c.id.c_str());if(!caseDir)caseDir=output->mkdir(c.id.c_str());Need(caseDir!=0,"cannot create case directory");TDirectory* centDir=caseDir->GetDirectory(c.cent.c_str());if(!centDir)centDir=caseDir->mkdir(c.cent.c_str());Need(centDir!=0,"cannot create centrality directory");TDirectory* dir=centDir->mkdir(c.plane==-1?"inclusive_all_y":c.plane==-2?"inclusive_in_range":c.plane==-3?"rapidity_underflow":c.plane==-4?"rapidity_overflow":Form("rapidity_bin_%02d",c.plane));Need(dir!=0,"cannot create plane directory");dir->cd();
      c.cf=new TGraphErrors;c.cf->SetName("CF_raw_signal");c.se=new TGraphErrors;c.se->SetName("SE_divided_by_norm_integral");c.me=new TGraphErrors;c.me->SetName("ME_divided_by_norm_integral");
      TH1D rawSE("SE_raw","Raw signal SE;k* [GeV/c];Pairs",int(c.bins.size()),&edges[0]);TH1D rawME("ME_raw","Raw signal ME;k* [GeV/c];Pairs",int(c.bins.size()),&edges[0]);rawSE.Sumw2();rawME.Sumw2();
      for(size_t i=0;i<c.bins.size();++i) {const Bin& b=c.bins[i];double x=(b.low+b.high)/2,ex=(b.high-b.low)/2;rawSE.SetBinContent(i+1,b.se);rawSE.SetBinError(i+1,std::sqrt(b.sv));rawME.SetBinContent(i+1,b.me);rawME.SetBinError(i+1,std::sqrt(b.mv));if(b.valid){int point=c.cf->GetN();c.cf->SetPoint(point,x,b.cf);c.cf->SetPointError(point,ex,b.error);}if(c.sh>0){c.se->SetPoint(i,x,b.se/c.sh);c.se->SetPointError(i,ex,std::sqrt(b.sv)/c.sh);}if(c.mh>0){c.me->SetPoint(i,x,b.me/c.mh);c.me->SetPointError(i,ex,std::sqrt(b.mv)/c.mh);}}
      Need(rawSE.Write()>0&&rawME.Write()>0&&c.cf->Write()>0&&c.se->Write()>0&&c.me->Write()>0,"companion object write failed");
      TNamed alpha("normalization",c.valid?Form("alpha=%.17g; SE(H)=%.17g; ME(H)=%.17g; CF=SE/(alpha ME); errors hold alpha fixed",c.alpha,c.sh,c.mh):"unavailable: empty normalization window");Need(alpha.Write()>0,"normalization write failed");
    }
    output->Close();Need(!output->TestBit(TFile::kWriteError),"companion write error");delete output;
    gStyle->SetOptStat(0);gStyle->SetTitleFontSize(.044);gStyle->SetPaperSize(28,21);Pages(curves,settings,outputDirectory);Need(!rootError,"ROOT error while drawing");
    for(Curves::iterator it=curves.begin();it!=curves.end();++it){delete it->second.cf;delete it->second.se;delete it->second.me;}
    SetErrorHandler(previous);std::cout<<"CF_CANVAS_PASS curves="<<curves.size()<<std::endl;
  } catch(const std::exception& e){SetErrorHandler(previous);std::cerr<<"CF_CANVAS_FAILED "<<e.what()<<std::endl;gSystem->Exit(1);}
}
#endif
