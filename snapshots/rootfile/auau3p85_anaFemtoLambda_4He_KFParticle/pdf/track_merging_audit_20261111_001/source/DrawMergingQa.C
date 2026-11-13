#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
void DrawMergingQa(const char*,const char*,const char*,double,double,double,double,double);
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ function DrawMergingQa(const char*,const char*,const char*,double,double,double,double,double);
#endif
static std::vector<std::string> MQSplit(const std::string& line) {
 std::vector<std::string> cols;std::istringstream s(line);std::string x;while(std::getline(s,x,'\t'))cols.push_back(x);return cols;
}
static void MQFrame(TPad* p) {p->SetLeftMargin(.14);p->SetRightMargin(.17);p->SetTopMargin(.20);p->SetBottomMargin(.14);}
static void MQAxes(TH1* h) {
 h->GetXaxis()->SetTitleSize(.047);h->GetYaxis()->SetTitleSize(.047);h->GetXaxis()->SetLabelSize(.042);h->GetYaxis()->SetLabelSize(.042);h->GetYaxis()->SetTitleOffset(1.35);
}
static void MQText(double x,double y,const std::string& text,double size=.035) {TLatex t;t.SetNDC();t.SetTextSize(size);t.DrawLatex(x,y,text.c_str());}
static void MQZero(double low,double high) {TLine* x=new TLine(low,0,high,0);TLine* y=new TLine(0,low,0,high);x->SetLineStyle(3);y->SetLineStyle(3);x->SetLineColor(kGray+2);y->SetLineColor(kGray+2);x->Draw();y->Draw();}
void DrawMergingQa(const char* metaFile,const char* cellFile,const char* output,double globalLow,double globalHigh,double zoomLow,double zoomHigh,double projectionHalfWidth) {
 try {
  if(!gSystem->AccessPathName(output))throw std::runtime_error("output exists");
  gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPalette(1);gStyle->SetNumberContours(60);
  std::ifstream meta(metaFile),cells(cellFile);if(!meta||!cells)throw std::runtime_error("missing TSV input");
  std::string line;std::getline(meta,line);if(line!="id\tspecies\tvariant\tname\tnx\txmin\txmax\tny\tymin\tymax\tentries\tregular\tregular_fraction\tsignal_entries\tfill_fraction\torigin005\torigin010\torigin020")throw std::runtime_error("metadata header mismatch");
  TH2D* hist[6]={0,0,0,0,0,0};std::string species[6],variant[6];double entries[6],regular[6],regfraction[6],fillfraction[6];int found=0;
  while(std::getline(meta,line)) {
   std::vector<std::string> v=MQSplit(line);if(v.size()!=18)throw std::runtime_error("bad metadata row");int id=std::atoi(v[0].c_str());if(id<0||id>=6||hist[id])throw std::runtime_error("invalid histogram id");
   hist[id]=new TH2D(Form("merging_%d",id),"",std::atoi(v[4].c_str()),std::atof(v[5].c_str()),std::atof(v[6].c_str()),std::atoi(v[7].c_str()),std::atof(v[8].c_str()),std::atof(v[9].c_str()));hist[id]->SetDirectory(0);hist[id]->GetXaxis()->SetTitle("#Delta#varphi* (rad)");hist[id]->GetYaxis()->SetTitle("#Delta#eta");hist[id]->GetZaxis()->SetTitle("SE pairs / bin");
   species[id]=v[1];variant[id]=v[2];entries[id]=std::atof(v[10].c_str());regular[id]=std::atof(v[11].c_str());regfraction[id]=std::atof(v[12].c_str());fillfraction[id]=std::atof(v[14].c_str());++found;
  }
  if(found!=6)throw std::runtime_error("six histograms required");
  std::getline(cells,line);if(line!="id\tix\tiy\tcontent\terror2")throw std::runtime_error("cells header mismatch");long rowCount=0;
  int id,ix,iy;double value,variance;
  while(cells>>id>>ix>>iy>>value>>variance) {if(id<0||id>=6||ix<0||ix>hist[id]->GetNbinsX()+1||iy<0||iy>hist[id]->GetNbinsY()+1||!TMath::Finite(value)||!TMath::Finite(variance)||variance<0)throw std::runtime_error("invalid cell");hist[id]->SetBinContent(ix,iy,value);hist[id]->SetBinError(ix,iy,std::sqrt(variance));++rowCount;}
  if(!cells.eof())throw std::runtime_error("truncated cells TSV");long expected=0;for(int j=0;j<6;++j){expected+=(hist[j]->GetNbinsX()+2)*(hist[j]->GetNbinsY()+2);hist[j]->SetEntries(entries[j]);}if(rowCount!=expected)throw std::runtime_error("cell count mismatch");
  TCanvas cover("merging_begin","",1300,1050);cover.Print((std::string(output)+"[").c_str());unsigned pages=0;
  for(int sp=0;sp<3;++sp) {
   const int first=sp*2;const std::string pretty=sp==0?"^{4}He":(sp==1?"^{3}He":"^{3}He under ^{4}He mass hypothesis");
   for(int zoom=0;zoom<2;++zoom) {
    const double lo=zoom?zoomLow:globalLow,hi=zoom?zoomHigh:globalHigh;TCanvas c(Form("merging_page_%d_%d",sp,zoom),"",1300,1050);c.Divide(2,2,.002,.002);for(int k=1;k<=4;++k){TPad* q=static_cast<TPad*>(c.GetPad(k));q->SetPad(q->GetXlowNDC(),.045+.91*q->GetYlowNDC(),q->GetXlowNDC()+q->GetWNDC(),.045+.91*(q->GetYlowNDC()+q->GetHNDC()));}
    double maximum=0;
    for(int j=first;j<first+2;++j)for(int x=1;x<=hist[j]->GetNbinsX();++x)for(int y=1;y<=hist[j]->GetNbinsY();++y)if(hist[j]->GetXaxis()->GetBinCenter(x)>=lo&&hist[j]->GetXaxis()->GetBinCenter(x)<=hi&&hist[j]->GetYaxis()->GetBinCenter(y)>=lo&&hist[j]->GetYaxis()->GetBinCenter(y)<=hi)maximum=TMath::Max(maximum,hist[j]->GetBinContent(x,y));
    for(int panel=0;panel<4;++panel) {c.cd(panel+1);TPad* pad=static_cast<TPad*>(gPad);MQFrame(pad);const bool logz=panel%2;int j=first+panel/2;pad->SetLogz(logz);hist[j]->GetXaxis()->SetRangeUser(lo,hi);hist[j]->GetYaxis()->SetRangeUser(lo,hi);hist[j]->SetMinimum(logz?.5:0.);hist[j]->SetMaximum(TMath::Max(logz?1.:1.,maximum));MQAxes(hist[j]);hist[j]->Draw("COLZ");MQZero(lo,hi);
     MQText(.14,.94,pretty+" : "+(variant[j]=="legacy"?"legacy fixed field":"event field, physical charge/p_{T}")+(logz?" (log z)":" (linear z)"),.034);
     MQText(.15,.90,Form("Entries %.0f; regular %.0f (%.2f%%)",entries[j],regular[j],100*regfraction[j]),.027);MQText(.15,.85,Form("QA fill / all-k* signal pairs = %.2f%%",100*fillfraction[j]),.027);
    }
    c.cd();MQText(.04,.014,"Raw signal-window SE only; no ME angular reference. Lines mark origin, not an applied cut.",.018);c.Print(output);++pages;
   }
   TCanvas c(Form("merging_projection_%d",sp),"",1300,1050);c.Divide(2,2,.002,.002);for(int k=1;k<=4;++k){TPad* q=static_cast<TPad*>(c.GetPad(k));q->SetPad(q->GetXlowNDC(),.045+.91*q->GetYlowNDC(),q->GetXlowNDC()+q->GetWNDC(),.045+.91*(q->GetYlowNDC()+q->GetHNDC()));}
   for(int panel=0;panel<4;++panel) {
    c.cd(panel+1);MQFrame(static_cast<TPad*>(gPad));gPad->SetRightMargin(.06);const bool yaxis=panel>=2,zoom=panel%2;TH1D* proj[2];double ymax=0;
    for(int v=0;v<2;++v) {TH2D* h=hist[first+v];h->GetXaxis()->SetRange(0,0);h->GetYaxis()->SetRange(0,0);TAxis* other=yaxis?h->GetXaxis():h->GetYaxis();const int a=other->FindBin(-projectionHalfWidth+1.e-9),b=other->FindBin(projectionHalfWidth-1.e-9);proj[v]=yaxis?h->ProjectionY(Form("mq_proj_%d_%d_%d",sp,panel,v),a,b,"e"):h->ProjectionX(Form("mq_proj_%d_%d_%d",sp,panel,v),a,b,"e");proj[v]->SetDirectory(0);proj[v]->GetXaxis()->SetTitle(yaxis?"#Delta#eta":"#Delta#varphi* (rad)");proj[v]->GetYaxis()->SetTitle("SE pairs / bin");proj[v]->GetXaxis()->SetRangeUser(zoom?zoomLow:globalLow,zoom?zoomHigh:globalHigh);ymax=TMath::Max(ymax,proj[v]->GetMaximum());proj[v]->SetLineColor(v?kRed+1:kBlue+1);proj[v]->SetMarkerColor(v?kRed+1:kBlue+1);proj[v]->SetMarkerStyle(v?24:20);proj[v]->SetMarkerSize(.5);MQAxes(proj[v]);}
    proj[0]->SetMinimum(0.);proj[0]->SetMaximum(ymax*1.4+1);proj[0]->Draw("E1");proj[1]->Draw("E1 SAME");TLegend* leg=new TLegend(.15,.74,.65,.87);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->AddEntry(proj[0],"Legacy fixed field","lp");leg->AddEntry(proj[1],"Event-field QA","lp");leg->Draw();MQText(.14,.94,pretty+" : "+(yaxis?"#Delta#eta projection":"#Delta#varphi* projection"),.034);MQText(.15,.67,Form("Other axis restricted to #pm%.3f",projectionHalfWidth),.030);
   }
   c.cd();MQText(.04,.014,"Raw SE projections with count-only errors; finite angular window. A flat raw map does not establish two-track efficiency.",.018);c.Print(output);++pages;
  }
  cover.Print((std::string(output)+"]").c_str());std::cout<<"MERGING_PLOT_PASS pages="<<pages<<" histograms="<<found<<" cells="<<rowCount<<std::endl;
 }catch(const std::exception& e){std::cerr<<"MERGING_PLOT_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
