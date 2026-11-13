#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
struct BaselineRow{std::string kind,group,label,series;double lo,hi,value,error;};
static std::vector<std::string> BaselineSplit(const std::string& s){std::vector<std::string> v;std::istringstream q(s);std::string x;while(std::getline(q,x,'\t'))v.push_back(x);return v;}
static void BaselinePage(const std::vector<BaselineRow>& rows,const std::vector<std::string>& groups,const std::string& kind,const std::vector<std::string>& series,const std::string& title,const std::string& output,double xmin,double xmax,double ymin,double ymax){
 const int cols=groups.size()>5?3:3,rs=(groups.size()+cols-1)/cols;TCanvas c("baseline_canvas","Mixing baseline audit",1500,rs*440+75);c.Divide(cols,rs,.002,.002);
 for(size_t ip=0;ip<groups.size();++ip){TPad* p=static_cast<TPad*>(c.GetPad(ip+1));p->SetPad(p->GetXlowNDC(),p->GetYlowNDC()*.94,p->GetXlowNDC()+p->GetWNDC(),(p->GetYlowNDC()+p->GetHNDC())*.94);}
 const int colors[]={kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1};
 const char* descriptions[]={"Raw mass gate","Method3 (L)","Method5 (L)","Method6 (shared ME, L)","Per-cent exposure matched"};
 for(size_t j=0;j<groups.size();++j){c.cd(j+1);gPad->SetTopMargin(.16);gPad->SetBottomMargin(.15);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.035);
  TH1D* f=new TH1D(Form("bl_frame_%s_%lu",kind.c_str(),(unsigned long)j),"",1,xmin,xmax);f->SetDirectory(0);f->SetMinimum(ymin);f->SetMaximum(ymax);f->GetXaxis()->SetTitle("k* (GeV/c)");f->GetYaxis()->SetTitle("SE/ME, scaled in 0.40-0.50");f->GetXaxis()->SetTitleSize(.047);f->GetYaxis()->SetTitleSize(.047);f->GetXaxis()->SetLabelSize(.043);f->GetYaxis()->SetLabelSize(.043);f->GetYaxis()->SetTitleOffset(1.45);f->Draw("AXIS");
  TLine* unity=new TLine(xmin,1,xmax,1);unity->SetLineStyle(2);unity->Draw();TLine* lo=new TLine(.4,ymin,.4,ymax);TLine* hi=new TLine(.5,ymin,.5,ymax);lo->SetLineStyle(3);hi->SetLineStyle(3);lo->SetLineColor(kGray+1);hi->SetLineColor(kGray+1);lo->Draw();hi->Draw();
  TLegend* l=new TLegend(.16,.63,.90,.83);l->SetFillStyle(0);l->SetBorderSize(0);l->SetTextSize(.032);
  std::string label=groups[j];
  for(size_t is=0;is<series.size();++is){int style=0;if(series[is]=="he3raw")style=1;else if(series[is]=="he3MappedRaw")style=2;else if(series[is]=="method3")style=1;else if(series[is]=="method5")style=2;else if(series[is]=="method6")style=3;else if(series[is]=="centExposureMatched")style=4;
   TGraphErrors* g=new TGraphErrors;g->SetMarkerColor(colors[style]);g->SetLineColor(colors[style]);g->SetMarkerStyle(20+style);g->SetMarkerSize(.7);
   for(size_t n=0;n<rows.size();++n){const BaselineRow& r=rows[n];if(r.kind!=kind||r.group!=groups[j]||r.series!=series[is])continue;label=r.label;if(!std::isfinite(r.value)||!std::isfinite(r.error)||r.lo<xmin-1.e-8||r.hi>xmax+1.e-8)continue;int i=g->GetN();g->SetPoint(i,(r.lo+r.hi)/2,r.value);g->SetPointError(i,0,r.error);}
   g->Draw("P SAME");const char* description=series[is]=="he4raw"?"#Lambda-^{4}He raw":(series[is]=="he3raw"?"#Lambda-^{3}He raw":(series[is]=="he3MappedRaw"?"^{3}He under ^{4}He mass hypothesis":descriptions[style]));l->AddEntry(g,description,"lp");
  }l->Draw();TLatex* h=new TLatex(.16,.91,("Centrality "+label).c_str());h->SetNDC();h->SetTextSize(.047);h->Draw();
 }
 c.cd();TLatex t;t.SetNDC();t.SetTextSize(.021);t.DrawLatex(.025,.987,title.c_str());t.SetTextSize(.016);t.DrawLatex(.025,.965,"0.40-0.50 is a reference scale, not an established physical plateau. Count-only statistical errors.");c.Print(output.c_str());
}
void DrawMixingBaseline(const char* input,const char* outDir,double xmin,double xmax,double ymin,double ymax,double fineYmin,double fineYmax){
 try{gStyle->SetOptStat(0);gStyle->SetOptTitle(0);std::ifstream in(input);if(!in)throw std::runtime_error("missing baseline TSV");std::string s;std::getline(in,s);if(s!="kind\tgroup\tlabel\tseries\tlow\thigh\tvalue\terror")throw std::runtime_error("baseline TSV header mismatch");std::vector<BaselineRow> rows;
  while(std::getline(in,s)){if(s.empty())continue;std::vector<std::string> v=BaselineSplit(s);if(v.size()!=8)throw std::runtime_error("invalid row");BaselineRow r;r.kind=v[0];r.group=v[1];r.label=v[2];r.series=v[3];r.lo=std::strtod(v[4].c_str(),0);r.hi=std::strtod(v[5].c_str(),0);r.value=std::strtod(v[6].c_str(),0);r.error=std::strtod(v[7].c_str(),0);rows.push_back(r);}
  std::vector<std::string> groups;const char* ids[]={"cent20_80","cent0_20","cent0_10","cent10_20","cent20_60"};for(int i=0;i<5;++i)groups.push_back(ids[i]);std::vector<std::string> fine;for(int i=0;i<9;++i)fine.push_back(Form("cent%d",i));std::string out=outDir;gSystem->mkdir(outDir,true);
  if(!rows.empty()&&rows[0].kind=="control"){std::vector<std::string> control;control.push_back("he4raw");control.push_back("he3raw");control.push_back("he3MappedRaw");BaselinePage(rows,groups,"control",control,"Species/mass-hypothesis control: raw mass-gate ratios, no purity correction",out+"/baseline_species_control.pdf",xmin,xmax,ymin,ymax);return;}
  std::vector<std::string> ss;ss.push_back("raw");ss.push_back("method3");ss.push_back("method5");ss.push_back("method6");BaselinePage(rows,groups,"group",ss,"Grouped high-k*: raw and adopted L-derived CFs",out+"/baseline_group_L.pdf",xmin,xmax,ymin,ymax);
  ss.resize(2);BaselinePage(rows,fine,"fine",ss,"Fine centralities: raw and Method3 Poisson-L CF",out+"/baseline_fine_L.pdf",xmin,xmax,fineYmin,fineYmax);
  ss[1]="centExposureMatched";BaselinePage(rows,groups,"group",ss,"Centrality-composition diagnostic only: per-cent ME exposure matching",out+"/baseline_cent_exposure.pdf",xmin,xmax,ymin,ymax);
 }catch(const std::exception& e){std::cerr<<"BASELINE_DRAW_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
