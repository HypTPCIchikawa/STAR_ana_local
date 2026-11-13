#include <TCanvas.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace HighKPlot {
struct Point { double low, high, value, error; };
struct Series {
  std::string id, group, label, family, species, stage, mass, selection;
  std::vector<HighKPlot::Point> bins;
};
struct Curve {
  std::string family, species, stage, mass, selection, legend;
  Curve() {}
  Curve(const char* f, const char* s, const char* st, const char* m,
        const char* cut, const char* text):
    family(f),species(s),stage(st),mass(m),selection(cut),legend(text) {}
};
struct Page {
  std::string title, subtitle;
  std::vector<std::string> groups;
  std::vector<HighKPlot::Curve> curves;
};
struct PointOrder { bool operator()(const Point& a,const Point& b) const { return a.low<b.low; } };

bool Fail(const std::string& text) {
  std::cerr << "DrawHighK ERROR: " << text << std::endl;
  return false;
}
std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> cells;
  std::string::size_type begin=0,pos;
  while ((pos=line.find('\t',begin))!=std::string::npos) {
    cells.push_back(line.substr(begin,pos-begin)); begin=pos+1;
  }
  cells.push_back(line.substr(begin));
  return cells;
}
bool Number(const std::string& text,double& value,bool allowNaN=false) {
  if (text.empty()) return false;
  char* end=0; errno=0;
  value=std::strtod(text.c_str(),&end);
  return end==text.c_str()+text.size() && errno!=ERANGE && (TMath::Finite(value) || (allowNaN && value!=value));
}
bool Read(const char* path,std::map<std::string,HighKPlot::Series>& series) {
  std::ifstream file(path);
  if (!file) return Fail(std::string("cannot open ")+path);
  std::string line;
  if (!std::getline(file,line)) return Fail("empty input");
  if (!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
  const std::string expected="series_id\tgroup\tlabel\tfamily\tspecies\tstage\tmass\tselection\tlow\thigh\tvalue\terror";
  if (line!=expected) return Fail("points.tsv header does not match the required 12-column schema");
  unsigned long row=1,points=0,missing=0;
  while (std::getline(file,line)) {
    ++row;
    if (!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
    std::ostringstream context; context << "row " << row << ": ";
    if (line.empty()) return Fail(context.str()+"empty row");
    std::vector<std::string> f=Split(line);
    if (f.size()!=12) return Fail(context.str()+"expected 12 tab-separated fields");
    for (unsigned j=0;j<8;++j) {
      if (f[j].empty()) return Fail(context.str()+"empty metadata field");
      for (unsigned k=0;k<f[j].size();++k)
        if ((unsigned char)f[j][k]<32 || (unsigned char)f[j][k]>126)
          return Fail(context.str()+"metadata must contain printable ASCII");
    }
    Point p;
    if (!Number(f[8],p.low) || !Number(f[9],p.high) ||
        !Number(f[10],p.value,true) || !Number(f[11],p.error,true) ||
        p.low<0 || p.high<=p.low || p.error<0)
      return Fail(context.str()+"invalid bin edge, value, or uncertainty");
    Series value;
    value.id=f[0];value.group=f[1];value.label=f[2];value.family=f[3];
    value.species=f[4];value.stage=f[5];value.mass=f[6];value.selection=f[7];
    std::map<std::string,HighKPlot::Series>::iterator found=series.find(value.id);
    if (found==series.end()) series[value.id]=value;
    else {
      Series& prev=found->second;
      if (prev.group!=value.group || prev.label!=value.label || prev.family!=value.family ||
          prev.species!=value.species || prev.stage!=value.stage || prev.mass!=value.mass ||
          prev.selection!=value.selection)
        return Fail(context.str()+"series metadata changed within one series_id");
    }
    if (p.value!=p.value || p.error!=p.error) { ++missing; continue; }
    series[value.id].bins.push_back(p); ++points;
  }
  if (file.bad() || !points) return Fail("input read failed or has no points");
  for (std::map<std::string,HighKPlot::Series>::iterator i=series.begin();i!=series.end();++i) {
    std::vector<HighKPlot::Point>& bins=i->second.bins;
    std::sort(bins.begin(),bins.end(),PointOrder());
    for (unsigned j=1;j<bins.size();++j)
      if (bins[j].low<bins[j-1].high-1.e-10)
        return Fail("duplicate or overlapping bins in series "+i->first);
  }
  std::cout << "DrawHighK: validated " << points << " finite bins in " << series.size() << " series; skipped " << missing << " NA bins (not drawn, not replaced by zero)" << std::endl;
  return true;
}
const Series* Find(const std::map<std::string,HighKPlot::Series>& all,const std::string& group,const Curve& c) {
  const Series* result=0;
  for (std::map<std::string,HighKPlot::Series>::const_iterator i=all.begin();i!=all.end();++i) {
    const Series& s=i->second;
    if (s.group==group && s.family==c.family && s.species==c.species &&
        s.stage==c.stage && s.mass==c.mass && s.selection==c.selection) {
      if (result) { Fail("ambiguous series selector for "+group+" "+c.family+" "+c.selection); gSystem->Exit(2); }
      result=&s;
    }
  }
  return result;
}
std::string GroupLabel(const std::string& group) {
  if (group=="cent0_20") return "0-20%";
  if (group=="cent20_80") return "20-80%";
  if (group=="cent10_20") return "10-20%";
  const char* labels[]={"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  for (int i=0;i<9;++i) if (group==Form("cent%d",i)) return labels[i];
  return group;
}
std::vector<std::string> Major() {
  std::vector<std::string> result;
  result.push_back("cent0_20"); result.push_back("cent20_80"); result.push_back("cent10_20");
  return result;
}
std::vector<std::string> Fine() {
  std::vector<std::string> result;
  for (int i=8;i>=0;--i) result.push_back(Form("cent%d",i));
  return result;
}
Page Make(const char* title,const char* subtitle,bool fine) {
  Page p;p.title=title;p.subtitle=subtitle;p.groups=fine?Fine():Major();return p;
}
void Pair(std::vector<HighKPlot::Page>& pages,const Page& coarse,const char* fineTitle) {
  pages.push_back(coarse); Page fine=coarse;fine.groups=Fine();fine.title=fineTitle;pages.push_back(fine);
}
std::vector<HighKPlot::Page> Pages() {
  std::vector<HighKPlot::Page> pages;
  Page p=Make("Raw mass-gated CF: species comparison","Selected PID; Lambda signal mass gate",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","#Lambda-^{4}He"));
  p.curves.push_back(Curve("baseline","3He","selected","signal","all","#Lambda-^{3}He"));
  p.curves.push_back(Curve("baseline","3HeAs4He","selected","signal","all","#Lambda-^{3}He (^{4}He mass)"));
  Pair(pages,p,"Raw mass-gated CF: species by cent9");
  p=Make("Lambda mass-gate dependence: selected 4He","Signal halves and sidebands are raw mass-gated CFs",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","signal"));
  p.curves.push_back(Curve("baseline","4He","selected","signalLow","all","signal low half"));
  p.curves.push_back(Curve("baseline","4He","selected","signalHigh","all","signal high half"));
  p.curves.push_back(Curve("baseline","4He","selected","SBNeg","all","Lambda left sideband"));
  p.curves.push_back(Curve("baseline","4He","selected","SBPos","all","Lambda right sideband"));
  Pair(pages,p,"Lambda mass gates by cent9: selected 4He");
  p=Make("Selected 4He: dE/dx nSigma dependence","Lambda signal; one-dimensional nSigma selection",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","selected: all"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","abs06","|n#sigma| < 0.6"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","abs10","|n#sigma| < 1.0"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","abs14","|n#sigma| < 1.4"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","abs20","|n#sigma| < 2.0"));
  Pair(pages,p,"Selected 4He nSigma cuts by cent9");
  p=Make("Selected 4He: nSigma sign comparison","Lambda signal; one-dimensional nSigma selection",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","selected: all"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","negative","-2 #leq n#sigma < 0"));
  p.curves.push_back(Curve("nsigma","4He","selected","signal","positive","0 #leq n#sigma < 2"));
  pages.push_back(p);
  p=Make("Selected 4He: rigidity and TOF boundary","Lambda signal; rigidity R=|p/q| [GeV/c]",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","selected: all"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","Rlow","R < 2: no TOF required"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","Rhigh","2 #leq R < 5: TOF + m^{2}"));
  Pair(pages,p,"Selected 4He: TOF boundary by cent9");
  p=Make("Selected 4He: finer rigidity intervals","Lambda signal; independent one-dimensional R selections",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","selected: all"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","Rlt10","R < 1.0"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","R10_15","1.0 #leq R < 1.5"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","R15_20","1.5 #leq R < 2.0"));
  p.curves.push_back(Curve("rigidity","4He","selected","signal","R20_50","2.0 #leq R < 5.0"));
  pages.push_back(p);
  p=Make("Selected 4He: pair transverse-mass dependence","Lambda signal; independent one-dimensional pair mT selections [GeV/c^{2}]",false);
  p.curves.push_back(Curve("baseline","4He","selected","signal","all","selected: all m_{T}"));
  p.curves.push_back(Curve("mt","4He","selected","signal","mt24_26","2.4 #leq m_{T} < 2.6"));
  p.curves.push_back(Curve("mt","4He","selected","signal","mt26_28","2.6 #leq m_{T} < 2.8"));
  p.curves.push_back(Curve("mt","4He","selected","signal","mt28_30","2.8 #leq m_{T} < 3.0"));
  p.curves.push_back(Curve("mt","4He","selected","signal","mt30_32","3.0 #leq m_{T} < 3.2"));
  Pair(pages,p,"Selected 4He: pair mT dependence by cent9");
  p=Make("4He hypothesis before final PID: nSigma dependence","PrePID is a mass hypothesis control; not an identified 4He sample",false);
  p.curves.push_back(Curve("prepid","4He","selected","signal","all","selected: all"));
  p.curves.push_back(Curve("prepid","4He","prePID","signal","abs06","prePID |n#sigma| < 0.6"));
  p.curves.push_back(Curve("prepid","4He","prePID","signal","abs10","prePID |n#sigma| < 1.0"));
  p.curves.push_back(Curve("prepid","4He","prePID","signal","abs14","prePID |n#sigma| < 1.4"));
  p.curves.push_back(Curve("prepid","4He","prePID","signal","abs20","prePID |n#sigma| < 2.0"));
  Pair(pages,p,"4He hypothesis prePID controls by cent9");
  p=Make("Normalization-window dependence: selected 4He","Lambda signal; each curve uses the displayed normalization interval",false);
  p.curves.push_back(Curve("normalization","4He","selected","signal","norm30_40","norm 0.30-0.40"));
  p.curves.push_back(Curve("normalization","4He","selected","signal","norm40_50","norm 0.40-0.50"));
  p.curves.push_back(Curve("normalization","4He","selected","signal","norm50_60","norm 0.50-0.60"));
  p.curves.push_back(Curve("normalization","4He","selected","signal","norm50_80","norm 0.50-0.80"));
  Pair(pages,p,"Normalization-window dependence by cent9");
  p=Make("Production comparison: selected 4He","Lambda signal; samples have different small excluded-input sets",false);
  p.curves.push_back(Curve("production","4He","selected","signal","old04BF","old 04BF"));
  p.curves.push_back(Curve("production","4He","selected","signal","new1CE","new 1CE hist-only"));
  Pair(pages,p,"Production comparison by cent9");
  return pages;
}
double Quantile(std::vector<double> values,double q) {
  std::sort(values.begin(),values.end());
  if (values.empty()) return 1.;
  return values[(unsigned)(q*(values.size()-1))];
}
void DrawPanel(const Page& page,const std::string& group,const std::map<std::string,HighKPlot::Series>& all,int unique) {
  gPad->SetLeftMargin(.12);gPad->SetRightMargin(.025);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.09);
  std::vector<const Series*> inputs;std::vector<double> lows,highs;
  for (unsigned c=0;c<page.curves.size();++c) {
    const Series* s=Find(all,group,page.curves[c]);inputs.push_back(s);
    if (!s) continue;
    for (unsigned b=0;b<s->bins.size();++b) {
      const Point& x=s->bins[b];
      if (x.high>.5 && x.low<1.) { lows.push_back(x.value-x.error);highs.push_back(x.value+x.error); }
    }
  }
  double lower=.65,upper=1.35;
  if (!lows.empty()) {
    lower=std::min(.85,Quantile(lows,.03));upper=std::max(1.15,Quantile(highs,.97));
    double margin=.12*(upper-lower);lower-=margin;upper+=margin;
  }
  TH1D* frame=new TH1D(Form("highk_frame_%d",unique),Form("%s;k* [GeV/c];raw CF",GroupLabel(group).c_str()),100,0,1.);
  frame->SetDirectory(0);frame->SetStats(0);frame->SetMinimum(lower);frame->SetMaximum(upper);
  frame->GetXaxis()->SetTitleSize(.047);frame->GetYaxis()->SetTitleSize(.047);frame->GetXaxis()->SetLabelSize(.04);frame->GetYaxis()->SetLabelSize(.04);
  frame->GetYaxis()->SetTitleOffset(1.15);frame->GetXaxis()->CenterTitle();frame->Draw();
  TBox* band=new TBox(.4,lower,.5,upper);band->SetFillColor(kGray);band->SetFillStyle(3004);band->SetLineColor(kGray);band->Draw("same");
  TLine* unity=new TLine(0,1,1,1);unity->SetLineStyle(2);unity->SetLineColor(kGray+2);unity->Draw();
  const int colors[]={kBlack,kBlue+1,kRed+1,kGreen+2,kMagenta+1,kOrange+7};
  const int markers[]={20,21,22,23,24,25};
  TLegend* legend=new TLegend(.15,.73,.96,.90);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.032);legend->SetNColumns(page.curves.size()>3?2:1);
  unsigned shown=0,clipped=0;
  for (unsigned c=0;c<inputs.size();++c) {
    const Series* s=inputs[c];if (!s || s->bins.empty()) continue;
    TGraphErrors* graph=new TGraphErrors((int)s->bins.size());
    graph->SetName(Form("highk_graph_%d_%u",unique,c));
    for (unsigned b=0;b<s->bins.size();++b) {
      const Point& x=s->bins[b]; graph->SetPoint(b,.5*(x.low+x.high),x.value);graph->SetPointError(b,0,x.error);
      if (x.low<1 && (x.value<lower || x.value>upper)) ++clipped;
    }
    graph->SetLineColor(colors[c%6]);graph->SetMarkerColor(colors[c%6]);graph->SetMarkerStyle(markers[c%6]);graph->SetMarkerSize(.55);graph->SetLineWidth(1);
    graph->Draw("PZ same");legend->AddEntry(graph,page.curves[c].legend.c_str(),"pe");++shown;
  }
  legend->Draw();
  TLatex text;text.SetNDC();text.SetTextSize(.028);text.SetTextColor(kGray+2);
  if (!shown) { text.SetTextSize(.05);text.DrawLatex(.25,.5,"No available series"); }
  else if (clipped) text.DrawLatex(.15,.015,"Frame follows high-k range; out-of-frame points omitted");
  gPad->RedrawAxis();
}
} // namespace HighKPlot

void DrawHighK(const char* points,const char* outPdf) {
  using namespace HighKPlot;
  if (!points || !outPdf || !*points || !*outPdf) { Fail("explicit input and output paths required");gSystem->Exit(2);return; }
  std::string output(outPdf);
  if (output.size()<4 || output.substr(output.size()-4)!=".pdf") { Fail("output must end in .pdf");gSystem->Exit(2);return; }
  if (!gSystem->AccessPathName(outPdf)) { Fail("refusing to overwrite an existing output PDF");gSystem->Exit(2);return; }
  std::map<std::string,HighKPlot::Series> all;
  if (!Read(points,all)) { gSystem->Exit(2);return; }
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetTitleFontSize(.045);
  std::vector<HighKPlot::Page> pages=Pages();unsigned drawn=0;
  TCanvas canvas("highk_systematics","High-k CF diagnostic atlas",1500,1050);
  canvas.Print((output+"[").c_str());
  for (unsigned i=0;i<pages.size();++i) {
    Page& p=pages[i];unsigned available=0;
    for (unsigned g=0;g<p.groups.size();++g) for (unsigned c=0;c<p.curves.size();++c) if (Find(all,p.groups[g],p.curves[c])) ++available;
    if (!available) { std::cout << "DrawHighK: skip unavailable page " << p.title << std::endl;continue; }
    canvas.Clear();canvas.cd();
    TLatex title;title.SetNDC();title.SetTextFont(42);title.SetTextSize(.026);title.DrawLatex(.035,.976,p.title.c_str());
    title.SetTextSize(.017);title.DrawLatex(.035,.947,p.subtitle.c_str());
    TPad* panels=new TPad(Form("highk_panels_%u",i),"panels",0,.087,1,.925);panels->Draw();panels->Divide(3,p.groups.size()>3?3:1,.002,.002);
    for (unsigned g=0;g<p.groups.size();++g) { panels->cd(g+1);DrawPanel(p,p.groups[g],all,(int)(100*i+g)); }
    canvas.cd();title.SetTextSize(.016);title.SetTextColor(kGray+2);
    title.DrawLatex(.035,.064,"Raw mass-gate ratios; count-only errors. Each nSigma, R, or mT cut is applied separately in one dimension.");
    title.DrawLatex(.035,.039,"Reference normalization: 0.40-0.50 GeV/c (shaded); no high-k plateau is assumed. Normalization scan: see legend.");
    title.DrawLatex(.035,.014,"Common normalization and overlapping selections correlate points and curves; error bars omit these covariances. Empty bins are not drawn.");
    canvas.Print(outPdf);++drawn;
    std::cout << "DrawHighK PAGE " << drawn << ": " << p.title << "; available series=" << available << std::endl;
  }
  canvas.Print((output+"]").c_str());
  if (!drawn) { Fail("none of the requested page selectors matched the input");gSystem->Exit(2);return; }
  std::cout << "DrawHighK SUCCESS: " << drawn << " pages written to " << outPdf << std::endl;
}

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ function DrawHighK(const char*,const char*);
#endif
