// Render saved SE/ME mass projections with identical area normalization.
#include "DrawPurityFitReview.C"
#include "TStyle.h"

namespace NormalizedPurityCases {
using namespace PurityFitReviewDraw;

void Overlay(Drawing &d, int pad, const Row &r,
             const std::vector<PurityFitReviewDraw::Cell> &se,
             const std::vector<PurityFitReviewDraw::Cell> &me,
             double massLo, double massHi, double signalLo, double signalHi,
             bool central, bool logarithmic) {
  Pad(d,pad,logarithmic);
  std::vector<PurityFitReviewDraw::Cell> a=se,b=me;
  if(central)for(unsigned i=0;i<a.size();++i) {
    a[i].count=a[i].centralCount;a[i].error2=a[i].centralError2;
    b[i].count=b[i].centralCount;b[i].error2=b[i].centralError2;
  }
  TH1D *hs=Histogram(d,a,Form("npc_se_%d_%d",r.page,pad),kRed+1);
  TH1D *hm=Histogram(d,b,Form("npc_me_%d_%d",r.page,pad),kBlue+1);
  // Select the same physical mass bins before ROOT computes the normalization.
  Axes(hs,"Fraction / mass bin",massLo,massHi);
  Axes(hm,"Fraction / mass bin",massLo,massHi);
  int first=0,last=0;
  for(unsigned i=0;i<a.size();++i) {
    const double x=(a[i].low+a[i].high)/2.;
    if(massLo<=x && x<massHi){if(!first)first=i+1;last=i+1;}
  }
  Require(first && last>=first,"empty normalized mass range");
  hs->GetXaxis()->SetRange(first,last);hm->GetXaxis()->SetRange(first,last);
  const double ns=hs->Integral(hs->GetXaxis()->GetFirst(),hs->GetXaxis()->GetLast());
  const double nm=hm->Integral(hm->GetXaxis()->GetFirst(),hm->GetXaxis()->GetLast());
  double ymax=0,ymin=1.;
  for(unsigned i=0;i<a.size();++i) {
    const double x=(a[i].low+a[i].high)/2.;if(x<massLo || x>massHi)continue;
    const double values[]={ns>0?a[i].count/ns:0.,nm>0?b[i].count/nm:0.};
    const double errors[]={ns>0?std::sqrt(a[i].error2)/ns:0.,nm>0?std::sqrt(b[i].error2)/nm:0.};
    for(int s=0;s<2;++s) {
      ymax=std::max(ymax,values[s]+errors[s]);
      if(values[s]>0)ymin=std::min(ymin,values[s]);
    }
  }
  if(ymax<=0)ymax=1.;
  const double bottom=logarithmic?std::max(ymin*.35,1.e-9):0.;
  const double top=ymax*(logarithmic?3.:1.15);
  TH1D *frame=static_cast<TH1D*>(hs->Clone(Form("npc_frame_%d_%d",r.page,pad)));
  frame->Reset();frame->SetDirectory(0);d.Keep(frame);
  frame->SetMinimum(bottom);frame->SetMaximum(top);frame->Draw("AXIS");
  TH1 *shownSE=ns>0?hs->DrawNormalized("E1 SAME",1.):0;
  TH1 *shownME=nm>0?hm->DrawNormalized("E1 SAME",1.):0;
  SignalLines(d,signalLo,signalHi,bottom,top);
  Text(d,.13,.963,r.label+" | "+r.method+" | "+(central?"central cell":"fit window"),.04);
  Text(d,.13,.917,WindowLabel(r),.029);
  Text(d,.13,.873,"ROOT DrawNormalized: unit area in displayed mass range",.027);
  Text(d,.13,.829,"SE N="+Format(ns,7)+"; ME N="+Format(nm,7),.031);
  Text(d,.13,.786,"Normalized statistical errors; common area is treated as fixed",.024);
  TLegend *legend=new TLegend(.63,.59,.96,.73);
  legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.032);
  if(shownSE)legend->AddEntry(shownSE,"SE (true)","lep");
  if(shownME)legend->AddEntry(shownME,"ME (mix)","lep");
  legend->Draw();d.Keep(legend);
  if(!shownSE || !shownME)Text(d,.17,.68,"A projection has zero integral",.035);
}
}

void DrawNormalizedPurityCases(const char *metaPath,const char *binsPath,
                               const char *outputDir,double massLo,double massHi,
                               double signalLo,double signalHi) {
  using namespace PurityFitReviewDraw;
  try {
    Require(metaPath && binsPath && outputDir && *metaPath && *binsPath && *outputDir,
            "explicit paths required");
    Require(massLo<massHi && massLo<signalLo && signalLo<signalHi && signalHi<massHi,
            "invalid mass limits");
    Pages pages;ReadMeta(metaPath,pages);ReadBins(binsPath,pages);Validate(pages,massLo,massHi);
    const std::string book=std::string(outputDir)+"/normalized_purity_cases.pdf";
    Require(gSystem->AccessPathName(book.c_str()),"refusing existing overlay PDF");
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(20,26);
    bool opened=false;
    for(Pages::const_iterator it=pages.begin();it!=pages.end();++it) {
      const Page &p=it->second;
      const Row &a=p.rows.find("SE")->second,&b=p.rows.find("ME")->second;
      const std::vector<PurityFitReviewDraw::Cell> &se=p.bins.find("SE")->second;
      const std::vector<PurityFitReviewDraw::Cell> &me=p.bins.find("ME")->second;
      Drawing d(a.page);
      MassPanel(d,1,a,se,massLo,massHi,signalLo,signalHi);
      MassPanel(d,2,b,me,massLo,massHi,signalLo,signalHi);
      NormalizedPurityCases::Overlay(d,3,a,se,me,massLo,massHi,signalLo,signalHi,false,false);
      PullPanel(d,4,a,se,massLo,massHi);
      PullPanel(d,5,b,me,massLo,massHi);
      NormalizedPurityCases::Overlay(d,6,a,se,me,massLo,massHi,signalLo,signalHi,true,true);
      const std::string one=std::string(outputDir)+Form("/normalized_case_%03d.pdf",a.page);
      Require(gSystem->AccessPathName(one.c_str()),"refusing existing case PDF");
      if(!opened){d.canvas->Print((book+"[").c_str());opened=true;}
      d.canvas->Print(one.c_str());d.canvas->Print(book.c_str());
    }
    TCanvas close("npc_close","",20,20);if(opened)close.Print((book+"]").c_str());
  }catch(const std::exception &e) {
    std::cerr<<"NORMALIZED_PURITY_CASES_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);
  }
}
