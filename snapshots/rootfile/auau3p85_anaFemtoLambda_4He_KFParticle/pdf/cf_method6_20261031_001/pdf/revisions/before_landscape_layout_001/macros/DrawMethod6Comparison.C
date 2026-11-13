// Read-only rendering of explicitly computed CF and purity points.
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Method6Plot {
struct Point { double x,y,e; };
typedef std::map<std::string,std::vector<Method6Plot::Point> > Table;
std::vector<std::string> Split(const std::string &s) {
 std::vector<std::string> out;std::string::size_type a=0,b;
 while((b=s.find('\t',a))!=std::string::npos){out.push_back(s.substr(a,b-a));a=b+1;}
 out.push_back(s.substr(a));return out;
}
double Number(const std::string &s) {
 char* end=0;const double v=std::strtod(s.c_str(),&end);
 if(!end || *end)throw std::runtime_error("invalid number: "+s);return v;
}
void Read(const char* path,bool purity,Table &data,std::vector<std::string> &groups,
          std::map<std::string,std::string> &labels) {
 std::ifstream f(path);if(!f.good())throw std::runtime_error("cannot read points");
 std::string line;std::getline(f,line);
 const std::string expected=purity?"method\tgroup\tlabel\tside\tbin\tlow\thigh\tvalue\terror":
  "method\tgroup\tlabel\tbin\tlow\thigh\tvalue\terror";
 if(line!=expected)throw std::runtime_error("unexpected TSV header");
 while(std::getline(f,line)) {
  if(line.empty())continue;const std::vector<std::string> c=Split(line);
  if(c.size()!=(purity?9u:8u))throw std::runtime_error("point column count differs");
  if(labels.find(c[1])==labels.end()){groups.push_back(c[1]);labels[c[1]]=c[2];}
  if(labels[c[1]]!=c[2])throw std::runtime_error("group label mismatch");
  const int offset=purity?1:0;Point p;
  p.x=(Number(c[4+offset])+Number(c[5+offset]))/2;
  p.y=Number(c[6+offset]);p.e=Number(c[7+offset]);
  data[c[0]+"/"+c[1]+(purity?"/"+c[3]:"")].push_back(p);
 }
}
TGraphErrors* Graph(const std::vector<Method6Plot::Point>& points,double lo,double hi,int color,int marker) {
 TGraphErrors* g=new TGraphErrors();g->SetMarkerColor(color);g->SetLineColor(color);
 g->SetMarkerStyle(marker);g->SetMarkerSize(.85);
 for(unsigned i=0;i<points.size();++i) {
  const Point &p=points[i];if(p.x<lo || p.x>=hi || !TMath::Finite(p.y) || !TMath::Finite(p.e) || p.e<0)continue;
  const int n=g->GetN();g->SetPoint(n,p.x,p.y);g->SetPointError(n,0,p.e);
 }
 return g;
}
void Draw(const Table &data,const std::vector<std::string> &groups,
 const std::map<std::string,std::string> &labels,const std::string &output,
 bool purity,double xlo,double xhi,double normlo,double normhi) {
 const char* methods[]={"method3","method5","method6"};
 const char* names[]={"Method 3","Method 5","Method 6: shared ME"};
 const int colors[]={kBlack,kRed+1,kBlue+1};
 TCanvas c(purity?"m6_purity":"m6_cf","",1500,650);c.Print((output+"[").c_str());
 for(int set=0;set<2;++set) {
  c.Clear();const int first=set?2:0,n=set?3:2;c.Divide(n,1,.008,.04);
  for(int j=0;j<n;++j) {
   c.cd(j+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.13);
   const std::string id=groups[first+j];std::vector<TGraphErrors*> graphs;std::vector<std::string> titles;
   for(int m=0;m<(purity?2:3);++m)for(int s=0;s<(purity?2:1);++s) {
    const std::string key=std::string(methods[m])+"/"+id+(purity?(s?"/ME":"/SE"):"");
    Table::const_iterator it=data.find(key);if(it==data.end())throw std::runtime_error("missing series: "+key);
    graphs.push_back(Graph(it->second,xlo,xhi,colors[m],purity?(s?24:20):20+m));
    titles.push_back(std::string(names[m])+(purity?(s?" ME":" SE"):""));
   }
   double ymin=purity?.65:.7,ymax=purity?1.02:1.3;
   for(unsigned m=0;m<graphs.size();++m)for(int b=0;b<graphs[m]->GetN();++b){double x,y;graphs[m]->GetPoint(b,x,y);ymin=std::min(ymin,y-graphs[m]->GetErrorY(b));ymax=std::max(ymax,y+graphs[m]->GetErrorY(b));}
   const double padding=.08*(ymax-ymin);ymin-=padding;ymax+=padding;
   TH1D* frame=new TH1D(Form("m6_frame_%d_%d_%d",purity,set,j),"",100,xlo,xhi);frame->SetDirectory(0);
   frame->SetMinimum(ymin);frame->SetMaximum(ymax);frame->GetXaxis()->SetTitle("k* (GeV/c)");frame->GetYaxis()->SetTitle(purity?"#Lambda purity":"C(k*)");
   frame->GetXaxis()->SetTitleSize(.05);frame->GetYaxis()->SetTitleSize(.05);frame->GetXaxis()->SetLabelSize(.043);frame->GetYaxis()->SetLabelSize(.043);frame->GetYaxis()->SetTitleOffset(1.25);frame->Draw("AXIS");
   TLine* unity=new TLine(xlo,1,xhi,1);unity->SetLineStyle(2);unity->Draw();
   TLegend* legend=new TLegend(.43,.68,.93,.86);legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(.035);
   for(unsigned m=0;m<graphs.size();++m){graphs[m]->Draw("P SAME");legend->AddEntry(graphs[m],(titles[m]+(graphs[m]->GetN()?"":" (NA)")).c_str(),"lp");}
   legend->Draw();TPaveText* heading=new TPaveText(.10,.89,.98,.99,"NDC");heading->SetFillStyle(0);heading->SetBorderSize(0);heading->SetTextSize(.04);heading->AddText(("Centrality "+labels.find(id)->second).c_str());heading->Draw();
   if(purity){TLatex text;text.SetNDC();text.SetTextSize(.028);text.DrawLatex(.18,.19,"Method 6 uses the Method 5 ME purity and error");}
  }
  c.cd();TLatex text;text.SetNDC();text.SetTextSize(.022);text.DrawLatex(.02,.973,Form("#Lambda-^{4}He; normalization %.2f-%.2f GeV/c; %s",normlo,normhi,purity?"SE/ME purity fits":"correlated shared-P and normalization errors"));
  c.Print(output.c_str());
 }
 c.Print((output+"]").c_str());
}
}
void DrawMethod6Comparison(const char* cfPath,const char* purityPath,const char* outputDir,
 double xMin,double xMax,double normLo,double normHi,double purityMax) {
 try {
  if(!cfPath || !*cfPath || !purityPath || !*purityPath || !outputDir || !*outputDir || !(xMin<xMax && normLo<normHi && xMin<purityMax))throw std::runtime_error("explicit paths/ranges required");
  const std::string out(outputDir),cf=out+"/cf_methods_3_5_6.pdf",purity=out+"/purity_se_me.pdf";
  if(!gSystem->AccessPathName(cf.c_str()) || !gSystem->AccessPathName(purity.c_str()))throw std::runtime_error("refusing existing comparison PDFs");
  gStyle->SetOptStat(0);gStyle->SetOptTitle(0);Method6Plot::Table cfs,ps;std::vector<std::string> groups,pg;std::map<std::string,std::string> labels,pl;
  Method6Plot::Read(cfPath,false,cfs,groups,labels);Method6Plot::Read(purityPath,true,ps,pg,pl);
  if(groups.size()!=5 || groups!=pg || labels!=pl)throw std::runtime_error("two/three-group layout requires five matching ordered groups");
  Method6Plot::Draw(cfs,groups,labels,cf,false,xMin,xMax,normLo,normHi);
  Method6Plot::Draw(ps,groups,labels,purity,true,xMin,purityMax,normLo,normHi);
 }catch(const std::exception &e){std::cerr<<"METHOD6_PLOT_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
