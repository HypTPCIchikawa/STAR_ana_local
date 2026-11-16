// Draw saved Lambda-mass fit diagnostics without refitting or opening ROOT data.
// ROOT 5 / C++98; all plotting ranges arrive through explicit arguments.
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PurityFitReviewDraw {
void Require(bool ok,const std::string &message) {
  if(!ok)throw std::runtime_error(message);
}
bool Finite(double x) {return TMath::Finite(x);}
bool Close(double a,double b) {
  return std::fabs(a-b)<=1.e-12*std::max(1.,std::max(std::fabs(a),std::fabs(b)));
}
std::vector<std::string> Fields(const std::string &line) {
  std::vector<std::string> result;std::string::size_type first=0;
  while(true) {
    const std::string::size_type next=line.find('\t',first);
    if(next==std::string::npos) {result.push_back(line.substr(first));return result;}
    result.push_back(line.substr(first,next-first));first=next+1;
  }
}
double Number(const std::string &text,bool optional=false) {
  if(optional && (text.empty() || text=="NA" || text=="null" ||
                  text=="nan" || text=="NaN"))
    return std::numeric_limits<double>::quiet_NaN();
  char *after=0;errno=0;const double value=std::strtod(text.c_str(),&after);
  Require(!text.empty() && after && *after=='\0' && after!=text.c_str() &&
          errno!=ERANGE && Finite(value),"invalid numeric TSV field: "+text);
  return value;
}
int Integer(const std::string &text) {
  const double value=Number(text);
  Require(value>=0 && value<=std::numeric_limits<int>::max() && std::floor(value)==value,
          "invalid nonnegative integer TSV field: "+text);
  return static_cast<int>(value);
}
std::string Display(std::string text) {
  const char *from[]={"\342\200\223","\342\200\224","\342\210\222"};
  for(int i=0;i<3;++i) {
    std::string::size_type at=0;
    while((at=text.find(from[i],at))!=std::string::npos) {text.replace(at,3,"-");++at;}
  }
  return text;
}
std::string Format(double value,int precision=4) {
  if(!Finite(value))return "NA";
  std::ostringstream out;out.precision(precision);out<<value;return out.str();
}
struct Row {
  int page,bin;std::string group,label,method,side,reason;
  double low,high,windowLow,windowHigh,status,covariance,ndf,chi2,probability;
  double purity,purityError,observations,parameters[6];
};
struct Cell {
  int index;double low,high,count,error2,centralCount,centralError2;
};
struct Page {
  std::map<std::string,Row> rows;
  std::map<std::string,std::vector<PurityFitReviewDraw::Cell> > bins;
};
typedef std::map<int,Page> Pages;

bool NextLine(std::ifstream &in,std::string &line) {
  while(std::getline(in,line)) {
    if(!line.empty() && line[line.size()-1]=='\r')line.resize(line.size()-1);
    if(!line.empty() && line[0]!='#')return true;
  }
  return false;
}
void ReadMeta(const char *path,Pages &pages) {
  std::ifstream in(path);Require(in.good(),std::string("cannot read metadata TSV: ")+path);
  std::string line;Require(NextLine(in,line),"empty metadata TSV");
  const std::string header="page\tgroup\tlabel\tmethod\tbin\tlow\thigh\twindowLow\twindowHigh\tside\tstatus\tcovariance\tndf\tchi2\tprobability\tpurity\tpurityError\tobservations\tp0\tp1\tp2\tp3\tp4\tp5\treason";
  Require(line==header,"unexpected metadata TSV header");
  while(NextLine(in,line)) {
    const std::vector<std::string> v=Fields(line);Require(v.size()==25,"metadata row must have 25 columns");
    Row r;r.page=Integer(v[0]);r.group=v[1];r.label=Display(v[2]);r.method=v[3];r.bin=Integer(v[4]);
    r.low=Number(v[5]);r.high=Number(v[6]);r.windowLow=Number(v[7]);r.windowHigh=Number(v[8]);r.side=v[9];
    r.status=Number(v[10],true);r.covariance=Number(v[11],true);r.ndf=Number(v[12],true);
    r.chi2=Number(v[13],true);r.probability=Number(v[14],true);r.purity=Number(v[15],true);
    r.purityError=Number(v[16],true);r.observations=Number(v[17],true);
    for(int p=0;p<6;++p)r.parameters[p]=Number(v[p+18],true);r.reason=Display(v[24]);
    Require(r.page>0 && r.bin>0,"page and grid bin must be one-based");
    Require(!r.group.empty() && r.group.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos,
            "unsafe group identifier");
    Require(r.method=="method3" || r.method=="method5","unknown fit method");
    Require(r.side=="SE" || r.side=="ME","unknown sample side");
    Require(r.high>r.low && r.windowHigh>r.windowLow &&
            r.windowLow<=r.low && r.windowHigh>=r.high,"invalid kstar cell/window bounds");
    Require(!pages[r.page].rows.count(r.side),"duplicate metadata page/side");
    pages[r.page].rows[r.side]=r;
  }
  Require(!pages.empty(),"no metadata rows");
}
void ReadBins(const char *path,Pages &pages) {
  std::ifstream in(path);Require(in.good(),std::string("cannot read mass-bin TSV: ")+path);
  std::string line;Require(NextLine(in,line),"empty mass-bin TSV");
  Require(line=="page\tside\tindex\tlow\thigh\tcount\terror2\tcentralCount\tcentralError2",
          "unexpected mass-bin TSV header");
  while(NextLine(in,line)) {
    const std::vector<std::string> v=Fields(line);Require(v.size()==9,"mass-bin row must have 9 columns");
    const int page=Integer(v[0]);const std::string side=v[1];
    Require(pages.count(page) && pages[page].rows.count(side),"mass bins reference an unknown page/side");
    Cell c;c.index=Integer(v[2]);c.low=Number(v[3]);c.high=Number(v[4]);
    c.count=Number(v[5]);c.error2=Number(v[6]);c.centralCount=Number(v[7]);c.centralError2=Number(v[8]);
    std::vector<PurityFitReviewDraw::Cell> &cells=pages[page].bins[side];
    Require(c.index==static_cast<int>(cells.size())+1,"mass bins must be ordered, contiguous and one-based");
    Require(c.high>c.low && c.count>=0 && c.error2>=0 && c.centralCount>=0 && c.centralError2>=0,
            "invalid mass-bin bounds/count/variance");
    if(!cells.empty())Require(Close(c.low,cells.back().high),"noncontiguous mass edges");
    cells.push_back(c);
  }
}
void Validate(const Pages &pages,double fitLo,double fitHi) {
  std::set<std::string> seen;
  for(Pages::const_iterator it=pages.begin();it!=pages.end();++it) {
    const Page &p=it->second;
    Require(p.rows.size()==2 && p.rows.count("SE") && p.rows.count("ME"),"each page requires SE and ME metadata");
    const Row &a=p.rows.find("SE")->second,&b=p.rows.find("ME")->second;
    Require(a.group==b.group && a.label==b.label && a.method==b.method && a.bin==b.bin &&
            Close(a.low,b.low) && Close(a.high,b.high) && Close(a.windowLow,b.windowLow) &&
            Close(a.windowHigh,b.windowHigh),"SE/ME page metadata disagree");
    Require(seen.insert(a.method+"/"+a.group+"/"+Format(a.bin,12)).second,"duplicate group/method/grid bin");
    Require(p.bins.count("SE") && p.bins.count("ME"),"each page requires SE and ME mass bins");
    const std::vector<PurityFitReviewDraw::Cell> &se=p.bins.find("SE")->second,&me=p.bins.find("ME")->second;
    Require(!se.empty() && se.size()==me.size(),"different or empty SE/ME mass axes");
    for(unsigned i=0;i<se.size();++i)
      Require(Close(se[i].low,me[i].low) && Close(se[i].high,me[i].high),"different SE/ME mass edges");
    Require((fitLo>=se.front().low || Close(fitLo,se.front().low)) &&
            (fitHi<=se.back().high || Close(fitHi,se.back().high)),"display fit range outside saved mass bins");
  }
}

struct Drawing {
  TCanvas *canvas;std::vector<TObject*> owned;
  Drawing(int page):canvas(new TCanvas(Form("pfr_page_%d",page),"Purity fit review",1500,900)) {
    canvas->Divide(3,2,.004,.004);
  }
  ~Drawing() {
    delete canvas;
    for(unsigned i=0;i<owned.size();++i)delete owned[i];
  }
  void Keep(TObject *object) {owned.push_back(object);}
};
void Text(Drawing &d,double x,double y,const std::string &text,double size=.031) {
  TLatex *t=new TLatex(x,y,text.c_str());t->SetNDC();t->SetTextFont(42);t->SetTextSize(size);t->Draw();d.Keep(t);
}
void Pad(Drawing &d,int number,bool log=false) {
  d.canvas->cd(number);gPad->SetLeftMargin(.13);gPad->SetRightMargin(.035);
  gPad->SetBottomMargin(.14);gPad->SetTopMargin(.265);gPad->SetLogy(log);
  gPad->SetTicks(1,1);
}
void Axes(TH1D *h,const char *yTitle,double fitLo,double fitHi) {
  h->SetStats(kFALSE);h->SetTitle("");h->GetXaxis()->SetTitle("M(p#pi^{-}) [GeV/c^{2}]");
  h->GetYaxis()->SetTitle(yTitle);h->GetXaxis()->SetRangeUser(fitLo,fitHi);
  h->GetXaxis()->SetTitleSize(.047);h->GetYaxis()->SetTitleSize(.047);
  h->GetXaxis()->SetLabelSize(.04);h->GetYaxis()->SetLabelSize(.04);
  h->GetXaxis()->SetTitleOffset(1.2);h->GetYaxis()->SetTitleOffset(1.3);
}
TH1D *Histogram(Drawing &d,const std::vector<PurityFitReviewDraw::Cell> &cells,const std::string &name,int color) {
  std::vector<double> edges;edges.push_back(cells.front().low);
  for(unsigned i=0;i<cells.size();++i)edges.push_back(cells[i].high);
  TH1D *h=new TH1D(name.c_str(),"",cells.size(),&edges[0]);h->SetDirectory(0);h->Sumw2();
  for(unsigned i=0;i<cells.size();++i) {h->SetBinContent(i+1,cells[i].count);h->SetBinError(i+1,std::sqrt(cells[i].error2));}
  h->SetLineColor(color);h->SetMarkerColor(color);h->SetMarkerStyle(20);h->SetMarkerSize(.45);d.Keep(h);return h;
}
TF1 *Model(Drawing &d,const Row &r,double fitLo,double fitHi,bool background) {
  for(int i=background?3:0;i<6;++i)if(!Finite(r.parameters[i]))return 0;
  if(!background && r.parameters[2]<=0)return 0;
  TF1 *f=new TF1(Form("pfr_%s_%s_%d_%u",background?"bg":"model",r.side.c_str(),r.page,static_cast<unsigned>(d.owned.size())),
                background?"pol2(0)":"gaus(0)+pol2(3)",fitLo,fitHi);
  for(int i=0;i<(background?3:6);++i)f->SetParameter(i,r.parameters[i+(background?3:0)]);
  f->SetLineColor(background?kGreen+2:kMagenta+1);f->SetLineStyle(background?2:1);
  f->SetLineWidth(2);f->SetNpx(500);d.Keep(f);return f;
}
void SignalLines(Drawing &d,double signalLo,double signalHi,double ymin,double ymax) {
  const double positions[]={signalLo,signalHi};
  for(int i=0;i<2;++i) {TLine *line=new TLine(positions[i],ymin,positions[i],ymax);
    line->SetLineColor(kGray+1);line->SetLineStyle(3);line->Draw();d.Keep(line);}
}
std::string WindowLabel(const Row &r) {
  return "k^{*}: ["+Format(r.low)+", "+Format(r.high)+
    "); fit window: ["+Format(r.windowLow)+", "+Format(r.windowHigh)+") GeV/c";
}
double CentralTotal(const std::vector<PurityFitReviewDraw::Cell> &cells) {
  double total=0;for(unsigned i=0;i<cells.size();++i)total+=cells[i].centralCount;return total;
}
void MassPanel(Drawing &d,int pad,const Row &r,const std::vector<PurityFitReviewDraw::Cell> &cells,
               double fitLo,double fitHi,double signalLo,double signalHi) {
  Pad(d,pad);TH1D *h=Histogram(d,cells,Form("pfr_mass_%s_%d",r.side.c_str(),r.page),r.side=="SE"?kRed+1:kBlue+1);
  TF1 *model=Model(d,r,fitLo,fitHi,false),*background=Model(d,r,fitLo,fitHi,true);
  double ymin=0,ymax=0;
  for(unsigned i=0;i<cells.size();++i) {
    const double x=(cells[i].low+cells[i].high)/2.;if(x<fitLo || x>fitHi)continue;
    ymax=std::max(ymax,cells[i].count+std::sqrt(cells[i].error2));
    TF1 *functions[]={model,background};
    for(int f=0;f<2;++f)if(functions[f]) {const double y=functions[f]->Eval(x);
      if(Finite(y)) {ymin=std::min(ymin,y);ymax=std::max(ymax,y);}}
  }
  ymax=ymax>0?ymax*1.15:1.;ymin=ymin<0?ymin*1.1:0;
  Axes(h,"Window counts / bin",fitLo,fitHi);h->SetMinimum(ymin);h->SetMaximum(ymax);h->Draw("E1");
  if(model)model->Draw("same");if(background)background->Draw("same");
  SignalLines(d,signalLo,signalHi,ymin,ymax);
  Text(d,.13,.963,r.label+" | "+r.method+" | "+r.side,.04);
  Text(d,.13,.917,WindowLabel(r),.029);
  Text(d,.13,.873,"status="+Format(r.status)+", cov="+Format(r.covariance)+
       ", #chi^{2}/ndf="+Format(r.chi2)+"/"+Format(r.ndf)+", prob="+Format(r.probability),.028);
  const std::string purity=Finite(r.purity) && Finite(r.purityError)?
    Format(r.purity)+" #pm "+Format(r.purityError):"NA";
  Text(d,.13,.829,"P="+purity+"; N(window)="+Format(r.observations,6)+
       ", N(central)="+Format(CentralTotal(cells),6),.028);
  Text(d,.13,.786,r.reason+"; window plotted, central N only",.022);
  TLegend *legend=new TLegend(.60,.60,.96,.73);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.028);
  legend->AddEntry(h,r.side.c_str(),"lep");if(model)legend->AddEntry(model,"Saved Gaussian + pol2","l");
  if(background)legend->AddEntry(background,"Saved background","l");legend->Draw();d.Keep(legend);
  if(!model)Text(d,.19,.66,"Saved model: NA",.04);
}
void PullPanel(Drawing &d,int pad,const Row &r,const std::vector<PurityFitReviewDraw::Cell> &cells,double fitLo,double fitHi) {
  Pad(d,pad);TF1 *model=Model(d,r,fitLo,fitHi,false);TGraph *points=new TGraph();d.Keep(points);
  double extent=4.;
  if(model)for(unsigned i=0;i<cells.size();++i) {
    const double x=(cells[i].low+cells[i].high)/2.;
    if(x<fitLo || x>fitHi || cells[i].error2<=0)continue;
    const double prediction=model->Eval(x);if(!Finite(prediction))continue;
    const double pull=(cells[i].count-prediction)/std::sqrt(cells[i].error2);
    if(!Finite(pull))continue;points->SetPoint(points->GetN(),x,pull);extent=std::max(extent,std::fabs(pull)*1.15);
  }
  TH1D *frame=new TH1D(Form("pfr_pull_frame_%s_%d",r.side.c_str(),r.page),"",1,fitLo,fitHi);
  frame->SetDirectory(0);d.Keep(frame);Axes(frame,"(Data - saved fit) / error",fitLo,fitHi);
  frame->SetMinimum(-extent);frame->SetMaximum(extent);frame->Draw("axis");
  TLine *zero=new TLine(fitLo,0,fitHi,0);zero->SetLineStyle(2);zero->SetLineColor(kGray+1);zero->Draw();d.Keep(zero);
  points->SetMarkerStyle(20);points->SetMarkerSize(.5);points->SetMarkerColor(r.side=="SE"?kRed+1:kBlue+1);points->Draw("P same");
  Text(d,.13,.963,r.side+" pull from saved fit",.04);Text(d,.13,.917,WindowLabel(r),.029);
  Text(d,.13,.873,"Zero-error bins are missing; no refit.",.03);
  Text(d,.13,.829,"Available pull bins: "+Format(points->GetN(),8),.03);
  if(!points->GetN())Text(d,.2,.58,"Pull: NA",.045);
}
double Normalize(TH1D *h,double fitLo,double fitHi) {
  double total=0;
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=fitLo && h->GetBinCenter(b)<=fitHi)total+=h->GetBinContent(b);
  if(total>0)h->Scale(1./total);return total;
}
void OverlayPanel(Drawing &d,int pad,const Row &r,const std::vector<PurityFitReviewDraw::Cell> &se,
                  const std::vector<PurityFitReviewDraw::Cell> &me,double fitLo,double fitHi,
                  double signalLo,double signalHi,bool log) {
  Pad(d,pad,log);TH1D *a=Histogram(d,se,Form("pfr_overlay_se_%d_%d",r.page,pad),kRed+1);
  TH1D *b=Histogram(d,me,Form("pfr_overlay_me_%d_%d",r.page,pad),kBlue+1);
  const double nse=Normalize(a,fitLo,fitHi),nme=Normalize(b,fitLo,fitHi);
  double ymax=0,positive=std::numeric_limits<double>::max();TH1D *histograms[]={a,b};
  for(int s=0;s<2;++s)for(int bin=1;bin<=histograms[s]->GetNbinsX();++bin) {
    const double x=histograms[s]->GetBinCenter(bin);if(x<fitLo || x>fitHi)continue;
    const double value=histograms[s]->GetBinContent(bin);
    ymax=std::max(ymax,value+histograms[s]->GetBinError(bin));if(value>0)positive=std::min(positive,value);
  }
  const double ymin=log?(positive<std::numeric_limits<double>::max()?positive*.4:1.e-5):0;
  ymax=ymax>0?ymax*(log?3.:1.2):1.;
  Axes(a,"Normalized window counts / bin",fitLo,fitHi);a->SetMinimum(ymin);a->SetMaximum(ymax);
  a->Draw("E1");b->Draw("E1 same");SignalLines(d,signalLo,signalHi,ymin,ymax);
  Text(d,.13,.963,std::string("SE / ME mass shape | ")+(log?"log":"linear"),.04);
  Text(d,.13,.917,WindowLabel(r),.029);
  Text(d,.13,.873,"Unit area in displayed fit range (window counts).",.028);
  Text(d,.13,.829,"SE normalization N="+(nse>0?Format(nse,6):"NA")+
       "; ME normalization N="+(nme>0?Format(nme,6):"NA"),.028);
  Text(d,.13,.786,"Central-cell counts are not drawn.",.026);
  TLegend *legend=new TLegend(.72,.63,.95,.73);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.032);
  legend->AddEntry(a,"SE","lep");legend->AddEntry(b,"ME","lep");legend->Draw();d.Keep(legend);
}
void DrawPage(Drawing &d,const Page &page,double fitLo,double fitHi,double signalLo,double signalHi) {
  const Row &se=page.rows.find("SE")->second,&me=page.rows.find("ME")->second;
  const std::vector<PurityFitReviewDraw::Cell> &a=page.bins.find("SE")->second,&b=page.bins.find("ME")->second;
  MassPanel(d,1,se,a,fitLo,fitHi,signalLo,signalHi);MassPanel(d,2,me,b,fitLo,fitHi,signalLo,signalHi);
  OverlayPanel(d,3,se,a,b,fitLo,fitHi,signalLo,signalHi,false);
  PullPanel(d,4,se,a,fitLo,fitHi);PullPanel(d,5,me,b,fitLo,fitHi);
  OverlayPanel(d,6,se,a,b,fitLo,fitHi,signalLo,signalHi,true);d.canvas->Update();
}
std::string SinglePath(const std::string &out,const Row &r) {
  return out+"/"+Form("fit_%s_%s_%03d.pdf",r.method.c_str(),r.group.c_str(),r.bin);
}
void Run(const char *metaPath,const char *binsPath,const char *outputDir,
         double fitLo,double fitHi,double signalLo,double signalHi) {
  Require(metaPath && *metaPath && binsPath && *binsPath && outputDir && *outputDir,
          "metaPath, binsPath, outputDir and all display/signal bounds are required");
  Require(Finite(fitLo) && Finite(fitHi) && Finite(signalLo) && Finite(signalHi) &&
          fitHi>fitLo && signalHi>signalLo && signalLo>=fitLo && signalHi<=fitHi,
          "invalid explicit display/signal bounds");
  Pages pages;ReadMeta(metaPath,pages);ReadBins(binsPath,pages);Validate(pages,fitLo,fitHi);
  const std::string out(outputDir);const char *methods[]={"method3","method5"};
  std::map<std::string,int> counts;
  for(Pages::const_iterator it=pages.begin();it!=pages.end();++it) {
    const Row &r=it->second.rows.find("SE")->second;++counts[r.method];
    Require(gSystem->AccessPathName(SinglePath(out,r).c_str()),"refusing to overwrite existing single-page PDF");
  }
  for(int m=0;m<2;++m) {
    Require(counts[methods[m]]>0,std::string("no pages for ")+methods[m]);
    Require(gSystem->AccessPathName((out+"/atlas_"+methods[m]+".pdf").c_str()),"refusing to overwrite existing atlas PDF");
  }
  if(gSystem->AccessPathName(out.c_str()))Require(gSystem->mkdir(out.c_str(),kTRUE)==0,"cannot create output directory");
  for(int m=0;m<2;++m) {
    const std::string atlas=out+"/atlas_"+methods[m]+".pdf";
    TCanvas book(Form("pfr_book_%d",m),"Atlas",1500,900);book.Print((atlas+"[").c_str());
    for(Pages::const_iterator it=pages.begin();it!=pages.end();++it) {
      const Row &r=it->second.rows.find("SE")->second;if(r.method!=methods[m])continue;
      Drawing drawing(it->first);DrawPage(drawing,it->second,fitLo,fitHi,signalLo,signalHi);
      drawing.canvas->Print(atlas.c_str());drawing.canvas->Print(SinglePath(out,r).c_str());
    }
    book.Print((atlas+"]").c_str());
    Require(!gSystem->AccessPathName(atlas.c_str()),"atlas PDF was not created");
    std::cout<<"PURITY_FIT_ATLAS_COMPLETE "<<atlas<<" pages="<<counts[methods[m]]<<std::endl;
  }
}
} // namespace PurityFitReviewDraw

void DrawPurityFitReview(const char *metaPath,const char *binsPath,const char *outputDir,
                         double fitLo,double fitHi,double signalLo,double signalHi) {
  try {PurityFitReviewDraw::Run(metaPath,binsPath,outputDir,fitLo,fitHi,signalLo,signalHi);}
  catch(const std::exception &e) {std::cerr<<"DrawPurityFitReview: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
