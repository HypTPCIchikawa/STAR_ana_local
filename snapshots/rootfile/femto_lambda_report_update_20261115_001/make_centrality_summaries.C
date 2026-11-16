void make_centrality_summaries(const char* metricsTsv,const char* parametersTsv,
                               const char* plotDir,const char* outputRoot);
#ifndef __CINT__
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace centrality_summary {
bool rootError=false; std::string rootMessage;
void ErrorHandler(Int_t level,Bool_t,const char* where,const char* message) {
  if(level>=kError){rootError=true;rootMessage=std::string(where?where:"ROOT")+": "+(message?message:"");}
  DefaultErrorHandler(level,kFALSE,where,message);
}
struct ErrorGuard {
  ErrorHandlerFunc_t previous;
  ErrorGuard():previous(SetErrorHandler(ErrorHandler)){rootError=false;rootMessage.clear();}
  ~ErrorGuard(){SetErrorHandler(previous);}
};
void Need(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
void Check(){Need(!rootError,"ROOT error: "+rootMessage);}
bool Finite(double value){return TMath::Finite(value)!=0;}
double Number(const std::string& text){std::istringstream in(text);double value=0.;Need(bool(in>>value)&&Finite(value),"Invalid number: "+text);in>>std::ws;Need(in.eof(),"Trailing numeric text");return value;}
std::string Format(double value,int digits=3){std::ostringstream out;out<<std::setprecision(digits)<<value;return out.str();}
std::vector<std::string> Fields(const std::string& line){std::vector<std::string> result;std::istringstream in(line);std::string field;while(std::getline(in,field,'\t'))result.push_back(field);return result;}
struct Params {
  std::map<std::string,std::string> values;
  explicit Params(const char* path){std::ifstream in(path);Need(bool(in),"Cannot read parameters");std::string line;while(std::getline(in,line)){if(line.empty())continue;std::vector<std::string> f=Fields(line);Need(f.size()==2,"Invalid parameter row");Need(values.insert(std::make_pair(f[0],f[1])).second,"Duplicate parameter");}Need(!in.bad(),"Parameter read failed");}
  std::string Text(const std::string& key)const{std::map<std::string,std::string>::const_iterator it=values.find(key);Need(it!=values.end(),"Missing parameter: "+key);return it->second;}
  double Value(const std::string& key)const{return Number(Text(key));}
  int Integer(const std::string& key)const{double value=Value(key);Need(value==std::floor(value),"Noninteger parameter");return int(value);}
};
struct Cent {std::string id,label;int cent9;double low,high;};
bool CentLess(const Cent& a,const Cent& b){return a.low<b.low;}
struct Row {
  std::string dataset,selection,metric;Cent cent;bool valid;double value,error,normSE,normME,totalSE,totalME;
};
typedef std::vector<Row> Rows;
void ReadRows(const char* path,Rows& rows,std::vector<Cent>& cents){
  std::ifstream in(path);Need(bool(in),"Cannot read metric TSV");std::string line;std::map<std::string,Cent> unique;std::set<std::string> seen;
  while(std::getline(in,line)){
    if(line.empty())continue;std::vector<std::string> f=Fields(line);Need(f.size()==14,"Invalid metric TSV row");Row r;
    r.dataset=f[0];r.selection=f[1];r.cent.id=f[2];r.cent.label=f[3];r.cent.cent9=int(Number(f[4]));r.cent.low=Number(f[5]);r.cent.high=Number(f[6]);r.metric=f[7];
    Need((r.dataset=="EP1"||r.dataset=="EP8")&&(r.metric=="tail"||r.metric=="shape")&&r.cent.low<r.cent.high,"Invalid metric identifiers");
    Need((f[8]=="NA")== (f[9]=="NA"),"Metric/error NA masks differ");r.valid=f[8]!="NA";r.value=0.;r.error=0.;
    if(r.valid){r.value=Number(f[8]);r.error=Number(f[9]);Need(r.error>=0.,"Negative metric error");}
    r.normSE=Number(f[10]);r.normME=Number(f[11]);r.totalSE=Number(f[12]);r.totalME=Number(f[13]);
    Need(seen.insert(r.dataset+"/"+r.cent.id+"/"+r.selection+"/"+r.metric).second,"Duplicate metric cell");
    unique[r.cent.id]=r.cent;rows.push_back(r);
  }
  Need(!in.bad()&&!rows.empty(),"Empty/failed metric TSV");
  for(std::map<std::string,Cent>::const_iterator it=unique.begin();it!=unique.end();++it)cents.push_back(it->second);
  std::sort(cents.begin(),cents.end(),CentLess);
}
const Row& Find(const Rows& rows,const std::string& dataset,const std::string& selection,const std::string& metric,const std::string& cent){
  for(size_t i=0;i<rows.size();++i)if(rows[i].dataset==dataset&&rows[i].selection==selection&&rows[i].metric==metric&&rows[i].cent.id==cent)return rows[i];
  throw std::runtime_error("Missing metric cell: "+dataset+"/"+cent+"/"+selection+"/"+metric);
}
void Label(double x,double y,const std::string& text,double size=.032){TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(size);t.DrawLatex(x,y,text.c_str());}
void Style(TPad* pad){pad->SetLeftMargin(.13);pad->SetRightMargin(.04);pad->SetTopMargin(.19);pad->SetBottomMargin(.14);pad->SetTicks(1,1);}
void Axes(TH1& h){h.SetStats(kFALSE);h.GetXaxis()->SetLabelFont(42);h.GetYaxis()->SetLabelFont(42);h.GetXaxis()->SetTitleFont(42);h.GetYaxis()->SetTitleFont(42);h.GetXaxis()->SetTitleSize(.044);h.GetYaxis()->SetTitleSize(.044);h.GetXaxis()->SetLabelSize(.034);h.GetYaxis()->SetLabelSize(.034);h.GetYaxis()->SetTitleOffset(1.3);}
std::string Window(const Params& p,const std::string& kind){return Format(p.Value(kind+".low"))+" #leq k^{*} < "+Format(p.Value(kind+".high"));}
std::string Metric(const std::string& metric){return metric=="tail"?"C_{tail}":"T";}
std::string MetricNote(const Params& p,const std::string& metric){return metric=="tail"?"Tail: "+Window(p,"tail")+" GeV/c; H: "+Format(p.Value("norm.low"))+"-"+Format(p.Value("norm.high")):"T=(S/M)_{"+Format(p.Value("tail.low"))+"-"+Format(p.Value("tail.high"))+"}/(S/M)_{"+Format(p.Value("reference.low"))+"-"+Format(p.Value("reference.high"))+"}; normalization cancels";}
void GraphStyle(TGraphErrors& graph,const Params& p,const std::string& dataset){int color=TColor::GetColor(p.Text("color."+dataset).c_str());graph.SetMarkerColor(color);graph.SetLineColor(color);graph.SetMarkerStyle(dataset=="EP1"?20:21);graph.SetMarkerSize(.95);graph.SetLineWidth(2);}
struct Range {double low,high;Range():low(1.),high(1.){}void Add(double y,double e){low=std::min(low,y-e);high=std::max(high,y+e);}void Pad(){double span=high-low;if(!(span>0.))span=.2;low-=.12*span;high+=.24*span;Need(Finite(low)&&Finite(high),"Invalid display range");}};
void AddRange(const Rows& rows,const std::vector<Cent>& cents,const std::string& selection,const std::string& metric,Range& range){
  const char* datasets[]={"EP1","EP8"};for(int d=0;d<2;++d)for(size_t i=0;i<cents.size();++i){const Row& r=Find(rows,datasets[d],selection,metric,cents[i].id);if(r.valid)range.Add(r.value,r.error);}
}
void SummaryPad(const Rows& rows,const std::vector<Cent>& cents,const Params& p,const std::string& selection,const std::string& metric,const std::string& title,const Range& range,std::vector<TObject*>& owned){
  Style(static_cast<TPad*>(gPad));const double xmin=cents.front().low,xmax=cents.back().high;
  TH1D* frame=new TH1D(("frame_"+selection+"_"+metric).c_str(),(";Centrality (%);"+Metric(metric)).c_str(),1,xmin,xmax);frame->SetDirectory(0);owned.push_back(frame);Axes(*frame);frame->SetMinimum(range.low);frame->SetMaximum(range.high);frame->Draw("AXIS");
  TLine unity(xmin,1.,xmax,1.);unity.SetLineColor(kGray+2);unity.SetLineStyle(2);unity.DrawClone("same");
  TLegend legend(.67,.70,.94,.81);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.030);
  const char* datasets[]={"EP1","EP8"};
  for(int d=0;d<2;++d){
    TGraphErrors* graph=new TGraphErrors();owned.push_back(graph);graph->SetName((selection+"_"+metric+"_"+datasets[d]).c_str());GraphStyle(*graph,p,datasets[d]);
    for(size_t i=0;i<cents.size();++i){const Row& r=Find(rows,datasets[d],selection,metric,cents[i].id);if(!r.valid)continue;int n=graph->GetN();graph->SetPoint(n,(r.cent.low+r.cent.high)/2.,r.value);graph->SetPointError(n,(r.cent.high-r.cent.low)/2.,r.error);}
    if(graph->GetN())graph->Draw("P same");
    const std::string label=std::string(datasets[d])+" ("+Format(graph->GetN())+"/"+Format(cents.size())+")";legend.AddEntry(graph,label.c_str(),"pe");
  }
  legend.DrawClone();Label(.13,.954,p.Text("header"),.030);Label(.13,.902,title,.035);Label(.13,.851,MetricNote(p,metric),.026);
  Label(.15,.755,"Single cent9 bins; all available errors shown",.027);Label(.15,.704,"Missing values are omitted",.026);gPad->RedrawAxis();
}
void Save(TCanvas& canvas,TFile& output,const std::string& directory,const std::string& stem){
  const std::string path=directory+"/"+stem+".pdf";canvas.Print(path.c_str());Check();output.cd();Need(canvas.Write(stem.c_str())>0,"Cannot write Canvas");Check();std::ifstream file(path.c_str(),std::ios::binary|std::ios::ate);Need(bool(file)&&file.tellg()>0,"Missing/empty PDF");std::cout<<"PLOT_OK\t"<<path<<std::endl;
}
void Summary(const Rows& rows,const std::vector<Cent>& cents,const Params& p,const std::string& selection,const std::string& metric,const std::string& title,const std::string& directory,const std::string& stem,TFile& output){
  Range range;AddRange(rows,cents,selection,metric,range);range.Pad();TCanvas canvas("summary_canvas","Centrality summary",1100,620);std::vector<TObject*> owned;
  SummaryPad(rows,cents,p,selection,metric,title,range,owned);Save(canvas,output,directory,stem);canvas.Clear();for(size_t i=0;i<owned.size();++i)delete owned[i];
}
void He4Summary(const Rows& rows,const std::vector<Cent>& cents,const Params& p,const std::string& metric,const std::string& directory,const std::string& stem,TFile& output){
  Range range;AddRange(rows,cents,"he4_below",metric,range);AddRange(rows,cents,"he4_above",metric,range);range.Pad();TCanvas canvas("he4_summary_canvas","He4 rapidity and centrality",1400,620);canvas.Divide(2,1,.002,.002);std::vector<TObject*> owned;
  canvas.cd(1);SummaryPad(rows,cents,p,"he4_below",metric,"^{4}He y_{lab} < "+Format(p.Value("he4.threshold")),range,owned);
  canvas.cd(2);SummaryPad(rows,cents,p,"he4_above",metric,"^{4}He y_{lab} > "+Format(p.Value("he4.threshold")),range,owned);
  Save(canvas,output,directory,stem);canvas.Clear();for(size_t i=0;i<owned.size();++i)delete owned[i];
}
struct Palette {
  int first,steps;double low,high;
  Palette(const Params& p,const std::string& metric):first(0),steps(p.Integer("map.color_steps")),low(p.Value("map."+metric+".min")),high(p.Value("map."+metric+".max")){
    Need(steps>1&&low<1.&&high>1.,"Invalid colour range");Double_t stops[3]={0.,(1.-low)/(high-low),1.};Double_t red[3]={.14,.98,.78};Double_t green[3]={.31,.98,.32};Double_t blue[3]={.49,.96,.12};first=TColor::CreateGradientColorTable(3,stops,red,green,blue,steps);Need(first>=0,"Cannot create palette");
  }
  int Color(double value)const{double f=(value-low)/(high-low);f=std::max(0.,std::min(1.,f));return first+int(f*(steps-1));}
};
void CellText(double x,double y,const std::string& text,double size,int color=kBlack){TLatex label;label.SetTextAlign(22);label.SetTextFont(42);label.SetTextSize(size);label.SetTextColor(color);label.DrawLatex(x,y,text.c_str());}
void Map(const Rows& rows,const std::vector<Cent>& cents,const Params& p,const std::string& dataset,const std::string& metric,const std::string& directory,const std::string& stem,TFile& output){
  const int nPair=p.Integer("pair.bins"),nCent=int(cents.size());const double xmin=p.Value("pair.min"),xmax=p.Value("pair.max"),width=(xmax-xmin)/nPair;Palette palette(p,metric);
  TCanvas canvas("map_canvas","Pair rapidity and centrality",1400,950);
  TPad* mapPad=new TPad("map_pad","",0.,.045,.88,1.);mapPad->Draw();mapPad->cd();mapPad->SetLeftMargin(.15);mapPad->SetRightMargin(.025);mapPad->SetBottomMargin(.12);mapPad->SetTopMargin(.18);mapPad->SetTicks(1,1);
  TH2D frame("map_axis_frame",";y_{pair}^{lab};Centrality (%)",nPair,xmin,xmax,nCent,0.,nCent);frame.SetDirectory(0);Axes(frame);frame.GetXaxis()->SetNdivisions(510);frame.GetYaxis()->SetTickLength(0.);
  for(int c=0;c<nCent;++c)frame.GetYaxis()->SetBinLabel(nCent-c,cents[c].label.c_str());frame.Draw("AXIS");
  int missing=0,clippedLow=0,clippedHigh=0,wideError=0,finite=0;TGraphErrors values;values.SetName("valid_map_cells");
  for(int c=0;c<nCent;++c)for(int b=0;b<nPair;++b){
    std::ostringstream selection;selection<<"pair_bin"<<std::setw(2)<<std::setfill('0')<<b+1;
    const Row& r=Find(rows,dataset,selection.str(),metric,cents[c].id);const double xlow=xmin+b*width,xhigh=xlow+width,ylow=nCent-c-1.,yhigh=ylow+1.,x=(xlow+xhigh)/2.,y=ylow+.5;
    TBox cell(xlow,ylow,xhigh,yhigh);cell.SetLineColor(TColor::GetColor("#AAAAAA"));cell.SetLineWidth(1);cell.SetFillColor(r.valid?palette.Color(r.value):TColor::GetColor("#EEEEEE"));cell.DrawClone("same");
    if(!r.valid){++missing;CellText(x,y,"NA",.026,kGray+2);continue;}
    ++finite;const bool below=r.value<palette.low,above=r.value>palette.high,wide=r.error>=palette.high-palette.low;
    if(below)++clippedLow;if(above)++clippedHigh;if(wide)++wideError;
    if(below||above){TBox hatch(xlow,ylow,xhigh,yhigh);hatch.SetFillColor(kBlack);hatch.SetFillStyle(3004);hatch.SetLineColor(kBlack);hatch.DrawClone("same");}
    const bool dark=std::fabs(r.value-1.)>.32*(palette.high-palette.low);const int textColor=dark?kWhite:kBlack;
    CellText(x,y+.16,Format(r.value)+(wide?" !":""),.024,textColor);CellText(x,y-.17,"#pm "+Format(r.error),.020,textColor);
    int point=values.GetN();values.SetPoint(point,c*nPair+b,r.value);values.SetPointError(point,0.,r.error);
  }
  frame.Draw("AXIS same");Label(.15,.957,p.Text("header")+"; "+dataset,.030);Label(.15,.909,Metric(metric)+" by pair rapidity and single cent9",.035);Label(.15,.859,MetricNote(p,metric),.025);
  canvas.cd();TPad* legend=new TPad("color_pad","",.88,.045,1.,1.);legend->Draw();legend->cd();legend->Range(0.,0.,1.,1.);
  const int blocks=palette.steps;for(int i=0;i<blocks;++i){double a=.26+.50*i/blocks,b=.26+.50*(i+1)/blocks;TBox box(.08,a,.30,b);box.SetLineColor(palette.first+i);box.SetFillColor(palette.first+i);box.DrawClone();}
  // Text size uses the smaller pad dimension; this narrow pad needs a
  // larger relative size than the main data pad to preserve readable labels.
  Label(.04,.89,"Colour",.12);Label(.04,.845,Metric(metric),.14);Label(.35,.755,"#geq "+Format(palette.high),.11);Label(.35,.255,"#leq "+Format(palette.low),.11);
  double one=.26+.50*(1.-palette.low)/(palette.high-palette.low);TLine marker(.05,one,.33,one);marker.SetLineColor(kBlack);marker.DrawClone();Label(.35,one-.01,"1",.11);
  Label(.04,.20,"NA: "+Format(missing),.105);Label(.04,.16,"Low: "+Format(clippedLow),.105);Label(.04,.12,"High: "+Format(clippedHigh),.105);Label(.04,.08,"!: "+Format(wideError),.105);
  canvas.cd();Label(.13,.023,"Gray: NA; hatched: outside colour range; !: #sigma #geq colour span. Values/errors unchanged.",.020);
  Save(canvas,output,directory,stem);output.mkdir((stem+"_data").c_str());output.cd((stem+"_data").c_str());Need(values.Write()>0,"Cannot write valid map points");Check();
  std::cout<<"MAP_COUNTS\t"<<dataset<<'\t'<<metric<<"\tvalid="<<finite<<"\tNA="<<missing<<"\tclipLow="<<clippedLow<<"\tclipHigh="<<clippedHigh<<"\twideError="<<wideError<<std::endl;
  canvas.Clear();
}
void Run(const char* metricsTsv,const char* parametersTsv,const char* plotDir,const char* outputRoot){
  ErrorGuard errors;gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);gStyle->SetOptStat(0);gStyle->SetEndErrorSize(3);
  Params p(parametersTsv);Rows rows;std::vector<Cent> cents;ReadRows(metricsTsv,rows,cents);Need(int(cents.size())==p.Integer("cent.count"),"Centrality count mismatch");
  const char* names[]={"inclusive_ctail_vs_centrality","inclusive_shape_vs_centrality","EP1_pairY_ctail_map","EP8_pairY_ctail_map","EP1_pairY_shape_map","EP8_pairY_shape_map","fixedmerge02_ctail_vs_centrality","fixedmerge02_shape_vs_centrality","he4_splits_ctail_vs_centrality","he4_splits_shape_vs_centrality"};
  Need(gSystem->AccessPathName(outputRoot),"Output ROOT exists");for(int i=0;i<10;++i)Need(gSystem->AccessPathName((std::string(plotDir)+"/"+names[i]+".pdf").c_str()),"Output PDF exists");
  if(gSystem->AccessPathName(plotDir))Need(gSystem->mkdir(plotDir,kTRUE)==0,"Cannot create plot directory");
  TFile output(outputRoot,"CREATE");Check();Need(!output.IsZombie()&&output.IsWritable(),"Cannot create output ROOT");
  Summary(rows,cents,p,"inclusive","tail","Inclusive pairs",plotDir,names[0],output);Summary(rows,cents,p,"inclusive","shape","Inclusive pairs",plotDir,names[1],output);
  Map(rows,cents,p,"EP1","tail",plotDir,names[2],output);Map(rows,cents,p,"EP8","tail",plotDir,names[3],output);Map(rows,cents,p,"EP1","shape",plotDir,names[4],output);Map(rows,cents,p,"EP8","shape",plotDir,names[5],output);
  const std::string merge=Format(p.Value("merge02.min"))+" #leq y_{pair}^{lab} < "+Format(p.Value("merge02.max"));
  Summary(rows,cents,p,"pair_merge02","tail",merge,plotDir,names[6],output);Summary(rows,cents,p,"pair_merge02","shape",merge,plotDir,names[7],output);
  He4Summary(rows,cents,p,"tail",plotDir,names[8],output);He4Summary(rows,cents,p,"shape",plotDir,names[9],output);
  output.Close();Check();std::cout<<"MAKE_CENTRALITY_SUMMARIES_OK plots=10 rows="<<rows.size()<<std::endl;
}
}
void make_centrality_summaries(const char* metricsTsv,const char* parametersTsv,const char* plotDir,const char* outputRoot){
  try{centrality_summary::Run(metricsTsv,parametersTsv,plotDir,outputRoot);}
  catch(const std::exception& e){std::cerr<<"MAKE_CENTRALITY_SUMMARIES_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
  catch(...){std::cerr<<"MAKE_CENTRALITY_SUMMARIES_FAILED: unknown exception"<<std::endl;gSystem->Exit(2);}
}
#endif
