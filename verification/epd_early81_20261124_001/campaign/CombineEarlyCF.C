#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace epdEarlyCF {
void Check(bool valid,const std::string& why) {
  if(!valid)throw std::runtime_error(why);
}
template<class T> T* Object(TFile& f,const char* key) {
  T* out=dynamic_cast<T*>(f.Get(key));Check(out!=0,key);return out;
}
std::string Text(TFile& f,const char* key) {
  TObject* value=f.Get(key);
  if(TObjString* text=dynamic_cast<TObjString*>(value)) return text->GetString().Data();
  if(TNamed* name=dynamic_cast<TNamed*>(value)) return name->GetTitle();
  throw std::runtime_error(key);
}
std::string Selection(const std::string& stored) {
  const std::string prefix="stored EPD definition/calibration\n";
  const size_t pos=stored.find(prefix);
  Check(pos!=std::string::npos,"missing EPD section");
  return stored.substr(0,pos);
}
std::string EP(const std::string& stored) {
  const std::set<std::string> paths={"recenterFile","recenterSha256","shiftFile","shiftSha256",
                                     "runListFile","runListSha256"};
  std::ostringstream result;
  std::istringstream text(stored);
  std::string line;
  while(std::getline(text,line)) if(!paths.count(line.substr(0,line.find(':'))))result<<line<<'\n';
  return result.str();
}
struct Scope {
  std::string selection,mixing,reduction,y,ep;
  std::set<std::string> inputPaths;
  std::set<std::string> baselineEP;
};

void Read(const char* list, bool early, Scope& scope,TH2D& se,TH2D& me,
          long long& saved,long long& epValid) {
  std::ifstream inputs(list);
  std::string path;
  int files=0;
  while(std::getline(inputs,path)) {
    Check(scope.inputPaths.insert(path).second,"duplicate replay ROOT "+path);
    TDirectory::TContext current(gDirectory);
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
    Check(f.get() && !f->IsZombie() && !f->TestBit(TFile::kRecovered),"invalid ROOT "+path);
    Check(Text(*f,"FemtoLambdaRunStatus")=="completed","incomplete replay "+path);
    const std::string s=Selection(Text(*f,"FemtoLambdaReducedReplaySelection"));
    const std::string mix=Text(*f,"FemtoLambdaReducedReplayMixing");
    const std::string reduced=Text(*f,"FemtoLambdaReducedTreeSettings");
    const std::string strict=Text(*f,"FemtoHe4RapiditySplitDefinition");
    const std::string ep=Text(*f,"FemtoEpdEventPlaneDefinition");
    Check(mix.find("bufferSize: 100")!=std::string::npos && mix.find("nEventPlaneBins: 8")!=std::string::npos &&
          reduced.find("resetPoolAtSegmentBoundary: false")!=std::string::npos,"wrong mixing");
    Check(Object<TParameter<Int_t> >(*f,"FemtoLambdaPairEpdBinCount")->GetVal()==8 &&
          Object<TParameter<Double_t> >(*f,"he4RapiditySplitYLabThreshold")->GetVal()==-.5,
          "wrong EP/Y controls");
    Check(!f->GetKey("FemtoLambdaMergingPairCutEnabled"),"unexpected pair cut");
    if(scope.selection.empty()) {
      scope.selection=s;scope.mixing=mix;scope.reduction=reduced;scope.y=strict;scope.ep=EP(ep);
    }
    Check(s==scope.selection&&mix==scope.mixing&&reduced==scope.reduction&&strict==scope.y&&EP(ep)==scope.ep,
          "selections differ across cohorts");
    if(!early) scope.baselineEP.insert(ep);
    if(early) Check(!scope.baselineEP.count(ep),"early output has late EP provenance");
    saved+=Object<TParameter<Long64_t> >(*f,"epdEventsSeen")->GetVal();
    epValid+=Object<TParameter<Long64_t> >(*f,"epdEventsValid")->GetVal();
    const char* keys[]={"hKstarSEVsCent_lambda_he4_signal_YLabBelow","hKstarMEVsCent_lambda_he4_signal_YLabBelow"};
    TH2D* destination[]={&se,&me};
    for(int k=0;k<2;++k) {
      TH2* raw=Object<TH2>(*f,keys[k]);
      TH2D* target=destination[k];
      Check(raw->GetNbinsX()==200 && raw->GetNbinsY()==9,"axis counts");
      for(int i=1;i<=201;++i) Check(std::fabs(raw->GetXaxis()->GetBinLowEdge(i)-target->GetXaxis()->GetBinLowEdge(i))<1e-12,"kstar edges");
      for(int i=1;i<=10;++i) Check(raw->GetYaxis()->GetBinLowEdge(i)==target->GetYaxis()->GetBinLowEdge(i),"cent edges");
      Check(target->Add(raw),"adding pair hist failed");
    }
    ++files;
  }
  Check(inputs.eof() && (early ? files>0 : files==20),"period coverage mismatch");
}
struct PairCounts {
  double se,me,sv,mv;
  PairCounts():se(0),me(0),sv(0),mv(0){}
};
PairCounts Sum(TH2D& se,TH2D& me,int c0,int c1,int x0,int x1) {
  PairCounts out;
  for(int cent=c0+1;cent<=c1+1;++cent)for(int x=x0;x<=x1;++x) {
    out.se+=se.GetBinContent(x,cent);
    out.me+=me.GetBinContent(x,cent);
    out.sv+=std::pow(se.GetBinError(x,cent),2);
    out.mv+=std::pow(me.GetBinError(x,cent),2);
  }
  return out;
}
double CF(const PairCounts& w, double alpha) {return w.se/(alpha*w.me);}
void Draw(TCanvas& page,const char* name, double end, const PairCounts* norm,
          TH2D** sh,TH2D** mh,int c0,int c1,const char* file) {
  TGraphErrors lines[3];
  const int color[]={kGray+2,kBlue+1,kRed+1};
  const char* label[]={"Previous","New early runs","Combined"};
  double minimum=.5,maximum=1.35;
  for(int j=0;j<3;++j) {
    if(norm[j].se<=0 || norm[j].me<=0)continue;
    const double a=norm[j].se/norm[j].me;
    lines[j].SetMarkerColor(color[j]);lines[j].SetLineColor(color[j]);
    lines[j].SetMarkerStyle(20+j);lines[j].SetMarkerSize(.65);
    for(int x=1;x+3<=200;x+=4) {
      const double lo=(x-1)*.005,hi=(x+3)*.005;
      if(hi>end+1e-12)break;
      const PairCounts c=Sum(*sh[j],*mh[j],c0,c1,x,x+3);
      if(c.me<=0)continue;
      const double v=CF(c,a);
      const double dv=std::sqrt(c.sv/(c.me*c.me)+c.se*c.se*c.mv/std::pow(c.me,4))/a;
      int i=lines[j].GetN();lines[j].SetPoint(i,(lo+hi)/2,v);lines[j].SetPointError(i,.01,dv);
      minimum=std::min(minimum,v-1.05*dv);
      maximum=std::max(maximum,v+1.05*dv);
    }
  }
  page.Clear();
  TH1* frame=page.DrawFrame(0,minimum,end,maximum);
  frame->SetTitle(Form("%s; k* [GeV/c]; C(k*) raw signal",name));
  TLine one(0,1,end,1);one.SetLineStyle(2);one.Draw("same");
  TLegend legend(.57,.70,.9,.89);legend.SetFillStyle(0);legend.SetBorderSize(0);
  for(int j=0;j<3;++j) {
    if(!lines[j].GetN())continue;
    lines[j].Draw("P E1 SAME");legend.AddEntry(&lines[j],label[j],"lp");
  }
  legend.Draw();
  page.Print(file);
  page.Clear();
}
}

void CombineEarlyCF(const char* oldList,const char* newList,const char* directory) {
  using namespace epdEarlyCF;
  try {
    const std::string dir(directory);
    TH2D oldSE("old_SE","previous SE",200,0.,1.,9,-.5,8.5);
    TH2D oldME("old_ME","previous ME",200,0.,1.,9,-.5,8.5);
    TH2D newSE("early_SE","early SE",200,0.,1.,9,-.5,8.5);
    TH2D newME("early_ME","early ME",200,0.,1.,9,-.5,8.5);
    TH2D allSE("combined_SE","combined SE",200,0.,1.,9,-.5,8.5);
    TH2D allME("combined_ME","combined ME",200,0.,1.,9,-.5,8.5);
    TH2D* se[]={&oldSE,&newSE,&allSE},*me[]={&oldME,&newME,&allME};
    for(int k=0;k<3;++k){se[k]->Sumw2();me[k]->Sumw2();}
    Scope scope;long long oldEvents=0,newEvents=0,oldValid=0,newValid=0;
    Read(oldList,false,scope,oldSE,oldME,oldEvents,oldValid);
    Read(newList,true,scope,newSE,newME,newEvents,newValid);
    allSE.Add(&oldSE,&newSE);allME.Add(&oldME,&newME);
    const PairCounts norm=Sum(oldSE,oldME,2,8,81,100);
    const PairCounts window=Sum(oldSE,oldME,2,8,121,200);
    Check(oldEvents==1819447887LL && oldValid==1736875306LL && norm.se==46968. &&
          std::fabs(CF(window,norm.se/norm.me)-.9777243956034068)<1e-12,
          "published old low-rapidity result was not reproduced");
    TFile root((dir+"/pair_counts.root").c_str(),"CREATE");
    Check(!root.IsZombie(),"pair counts output exists or is unwritable");
    for(int k=0;k<3;++k) Check(se[k]->Write()>0&&me[k]->Write()>0,"cannot write pair counts");
    root.Close();
    std::ofstream out((dir+"/integral_cf.tsv").c_str());
    out<<std::setprecision(12)<<"centrality\tcohort\tsavedEventsAllCent\tepValidEventsAllCent"
        "\tSEAllK\tMEAllK\tSEH\tMEH\tCF_060_100\n";
    const char* centers[]={"0-60%","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%"};
    const int low[]={2,7,6,5,4,3,2},high[]={8,8,6,5,4,3,2};
    const long long saved[]={oldEvents,newEvents,oldEvents+newEvents};
    const long long valid[]={oldValid,newValid,oldValid+newValid};
    const std::string pdf=dir+"/ep8_low_y_cf_before_after.pdf";
    TCanvas page("early81Combined","EP8 low Y",1050,760);
    page.Print((pdf+"[").c_str());
    for(int c=0;c<7;++c) {
      PairCounts norms[3];
      for(int j=0;j<3;++j) {
        const PairCounts all=Sum(*se[j],*me[j],low[c],high[c],0,201);
        norms[j]=Sum(*se[j],*me[j],low[c],high[c],81,100);
        const PairCounts w=Sum(*se[j],*me[j],low[c],high[c],121,200);
        out<<centers[c]<<'\t'<<(j==0?"previous":j==1?"early":"combined")<<'\t'<<saved[j]<<'\t'
           <<valid[j]<<'\t'<<all.se<<'\t'<<all.me<<'\t'<<norms[j].se<<'\t'<<norms[j].me<<'\t';
        if(norms[j].se>0&&norms[j].me>0&&w.me>0) out<<CF(w,norms[j].se/norms[j].me);else out<<"NA";
        out<<'\n';
      }
      const TString title=Form("%s | EP8, He4 y_{lab} < -0.5 | H=[0.4,0.5)",centers[c]);
      Draw(page,title,1.,norms,se,me,low[c],high[c],pdf.c_str());
      Draw(page,title,.3,norms,se,me,low[c],high[c],pdf.c_str());
    }
    page.Print((pdf+"]").c_str());
    Check(out.good(),"writing numerical comparison failed");
    std::cout<<"EARLY_WITH_BASELINE_CF_PASS previous="<<oldEvents<<" added="<<newEvents<<std::endl;
  } catch(const std::exception& ex) { std::cerr<<"EARLY_CF_ERROR "<<ex.what()<<std::endl; gSystem->Exit(1); }
}
