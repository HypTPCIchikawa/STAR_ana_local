#include "../femto_lambda_period_mix100_20261119_001/ComparePeriodMix.C"
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TROOT.h>
#include <cstdlib>

namespace MergingComparison {
using namespace PeriodMixComparison;
struct BoxStats {
  double count,variance,area;
  BoxStats():count(0),variance(0),area(0){}
};
struct RatioStats { BoxStats inner,ring;double ratio,error,covariance; };
double Overlap(double a,double b,double lo,double hi) {
  return std::max(0.,std::min(b,hi)-std::max(a,lo));
}
RatioStats Density(const TH2* h,double dx,double dy,double scale) {
  RatioStats s;s.ratio=s.error=s.covariance=0;
  for(int x=1;x<=h->GetNbinsX();++x)for(int y=1;y<=h->GetNbinsY();++y) {
    const double x0=h->GetXaxis()->GetBinLowEdge(x),x1=h->GetXaxis()->GetBinUpEdge(x);
    const double y0=h->GetYaxis()->GetBinLowEdge(y),y1=h->GetYaxis()->GetBinUpEdge(y);
    const double area=(x1-x0)*(y1-y0);
    const double ai=Overlap(x0,x1,-dx,dx)*Overlap(y0,y1,-dy,dy);
    const double ao=Overlap(x0,x1,-scale*dx,scale*dx)*Overlap(y0,y1,-scale*dy,scale*dy);
    const double wi=ai/area,wr=(ao-ai)/area,c=h->GetBinContent(x,y),e=h->GetBinError(x,y);
    Need(TMath::Finite(c)&&TMath::Finite(e)&&c>=0&&e>=0,"invalid angular-map content/error");
    s.inner.count+=wi*c;s.inner.variance+=wi*wi*e*e;s.inner.area+=ai;
    s.ring.count+=wr*c;s.ring.variance+=wr*wr*e*e;s.ring.area+=ao-ai;
    s.covariance+=wi*wr*e*e;
  }
  Need(s.inner.area>0&&s.ring.area>0&&s.inner.count>0&&s.ring.count>0,"empty density comparison");
  s.ratio=s.inner.count/s.inner.area/(s.ring.count/s.ring.area);
  const double relative=s.inner.variance/(s.inner.count*s.inner.count)+s.ring.variance/(s.ring.count*s.ring.count)-2*s.covariance/(s.inner.count*s.ring.count);
  Need(relative>=-1.e-12,"invalid density variance");s.error=s.ratio*std::sqrt(std::max(0.,relative));return s;
}
void RatioJson(std::ostream& out,const RatioStats& s) {
  out<<"{\"innerCount\":"<<s.inner.count<<",\"innerArea\":"<<s.inner.area
     <<",\"ringCount\":"<<s.ring.count<<",\"ringArea\":"<<s.ring.area
     <<",\"areaNormalizedInnerToRingRatio\":"<<s.ratio<<",\"countStatisticalError\":"<<s.error<<'}';
}
void Exact(const TH2* a,const TH2* b) {
  Need(std::string(a->ClassName())==b->ClassName()&&a->GetEntries()==b->GetEntries()&&Axis(a->GetXaxis(),b->GetXaxis())&&Axis(a->GetYaxis(),b->GetYaxis()),"bare/EventField aliases differ in class, axis or Entries");
  for(int x=0;x<=a->GetNbinsX()+1;++x)for(int y=0;y<=a->GetNbinsY()+1;++y)
    Need(a->GetBinContent(x,y)==b->GetBinContent(x,y)&&a->GetBinError(x,y)==b->GetBinError(x,y),"bare/EventField aliases differ in content/error");
}
TH2* CloneMap(TFile* f,const std::string& key,const std::string& name) {
  TH2* source=dynamic_cast<TH2*>(f->Get(key.c_str()));Need(source!=0,"missing angular histogram "+key);
  TH2* h=dynamic_cast<TH2*>(source->Clone(name.c_str()));Need(h!=0,"histogram clone failed");h->SetDirectory(0);return h;
}
double RangeMax(const TH2* h,double limit) {
  double result=0;
  for(int x=1;x<=h->GetNbinsX();++x)for(int y=1;y<=h->GetNbinsY();++y)
    if(std::fabs(h->GetXaxis()->GetBinCenter(x))<limit&&std::fabs(h->GetYaxis()->GetBinCenter(y))<limit)result=std::max(result,h->GetBinContent(x,y));
  return std::max(result,1.e-15);
}
void PrepareCanvas(TCanvas& canvas,const char* title,const std::string& line1,const std::string& line2) {
  canvas.Divide(2,2,.004,.004);
  for(int p=1;p<=4;++p) {
    canvas.cd(p);const double x=gPad->GetXlowNDC(),y=gPad->GetYlowNDC(),w=gPad->GetWNDC(),h=gPad->GetHNDC();
    gPad->SetPad(x,.13+.75*y,x+w,.13+.75*(y+h));gPad->SetLeftMargin(.14);gPad->SetRightMargin(.18);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.22);
  }
  canvas.cd();TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.028);text.DrawLatex(.025,.977,title);
  text.SetTextSize(.019);text.DrawLatex(.025,.942,line1.c_str());text.DrawLatex(.025,.914,line2.c_str());
  text.SetTextSize(.018);text.DrawLatex(.025,.095,"Rectangle: |#Delta#varphi*| < 0.05 rad AND |#Delta#eta| < 0.05. Target: proton-^{4}He pair only.");
  text.DrawLatex(.025,.066,"Pion rectangle is a reference, not the requested cut. Raw signal SE only; no ME or background subtraction.");
  text.DrawLatex(.025,.037,"Map integrals are bin-level reference estimates, not exact rejected-pair counts. A dip is not an efficiency measurement.");
}
void Panel(TH2* h,int panel,const char* title,double limit,double zmax,double phiCut,double etaCut,bool proton,bool normalized,const RatioStats& metric) {
  gPad->GetCanvas()->cd(panel);
  h->SetTitle("");h->GetXaxis()->SetTitle("#Delta#varphi* [rad]");h->GetYaxis()->SetTitle("#Delta#eta_{global momentum}");
  h->GetZaxis()->SetTitle(normalized?"Fraction / display bin":"Raw SE / display bin");
  h->GetXaxis()->SetRangeUser(-limit,limit);h->GetYaxis()->SetRangeUser(-limit,limit);
  h->GetXaxis()->SetTitleSize(.045);h->GetYaxis()->SetTitleSize(.045);h->GetZaxis()->SetTitleSize(.035);
  h->GetXaxis()->SetLabelSize(.038);h->GetYaxis()->SetLabelSize(.038);h->GetZaxis()->SetLabelSize(.030);h->GetYaxis()->SetTitleOffset(1.25);
  h->SetMinimum(0);h->SetMaximum(zmax*1.02);h->Draw("COLZ");
  TLine* a=new TLine(-limit,0,limit,0);TLine* b=new TLine(0,-limit,0,limit);a->SetLineColor(kGray+2);b->SetLineColor(kGray+2);a->SetLineStyle(3);b->SetLineStyle(3);a->Draw();b->Draw();
  TBox* box=new TBox(-phiCut,-etaCut,phiCut,etaCut);box->SetFillStyle(0);box->SetLineWidth(2);box->SetLineColor(proton?kRed+1:kGray+3);box->SetLineStyle(proton?1:2);box->Draw();
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.035);text.DrawLatex(.14,.95,title);
  text.SetTextSize(.028);text.DrawLatex(.14,.882,Form("Display %d x %d bins; R_{inner/ring}(0.05) = %.3f #pm %.3f",h->GetNbinsX(),h->GetNbinsY(),metric.ratio,metric.error));
}
} // namespace MergingComparison

// Metadata rows: currentInput currentSaved oldInput oldCache oldSelected oldJobs oldExcluded.
void DrawMergingCurrentVsOld(const char* inputList,const char* oldRoot,const char* metadataPath,
    const char* outputPrefix,const char* jsonPath,const char* aggregateRoot,
    double phiCut,double etaCut,double ringScale,double fullRange,double zoomRange,
    int rawRebin,int commonBins) {
  using namespace MergingComparison;
  TH1::AddDirectory(kFALSE);rootError=false;SetErrorHandler(ErrorHandler);
  try {
    Need(inputList&&oldRoot&&metadataPath&&outputPrefix&&jsonPath&&aggregateRoot&&phiCut>0&&etaCut>0&&ringScale>1&&fullRange>zoomRange&&zoomRange>=ringScale*std::max(phiCut,etaCut)&&rawRebin>0&&commonBins>0,"invalid explicit plot arguments");
    const TString pdf=TString(outputPrefix)+".pdf";
    Need(gSystem->AccessPathName(pdf)&&gSystem->AccessPathName(jsonPath)&&gSystem->AccessPathName(aggregateRoot),"refusing existing output");
    std::ifstream metadata(metadataPath);Long64_t currentInput=0,currentSaved=0,oldInput=0,oldCache=0,oldSelected=0;int oldJobs=0,oldExcluded=0;
    Need(bool(metadata>>currentInput>>currentSaved>>oldInput>>oldCache>>oldSelected>>oldJobs>>oldExcluded)&&currentInput>=currentSaved&&currentSaved>0&&oldInput>=oldCache&&oldSelected>0&&oldJobs>0,"invalid metadata values");
    std::ifstream list(inputList);Need(bool(list),"cannot open 20-output list");std::string path;TH2* current[2]={0,0};TH2* aliases[2]={0,0};const char* daughters[]={"proton","pion"};
    int sources=0;Long64_t saved=0;double selected=0;
    while(list>>path) {
      TFile* f=TFile::Open(path.c_str(),"READ");Need(f&&!f->IsZombie()&&!f->TestBit(TFile::kRecovered),"bad current input ROOT");
      Need(Object(f,"FemtoLambdaRunStatus")=="completed","current producer incomplete");
      saved+=Parameter(f,"reconstructedEvents");TH1* signal=dynamic_cast<TH1*>(f->Get("hKstarSE_lambda_he4_signal"));Need(signal!=0,"missing signal pair histogram");selected+=signal->GetEntries();
      for(int d=0;d<2;++d) {
        const std::string key=std::string("true/hDphiDeta_")+daughters[d]+"_4He";
        TH2* h=CloneMap(f,key,Form("part_%d_%d",sources,d));TH2* alias=CloneMap(f,key+"_EventField",Form("alias_%d_%d",sources,d));Exact(h,alias);
        if(!sources){current[d]=dynamic_cast<TH2*>(h->Clone(Form("current_%s",daughters[d])));current[d]->SetDirectory(0);aliases[d]=dynamic_cast<TH2*>(alias->Clone(Form("current_%s_EventField",daughters[d])));aliases[d]->SetDirectory(0);}
        else {Need(current[d]->Add(h)&&aliases[d]->Add(alias),"map sum failed");}
        delete h;delete alias;
      }
      f->Close();delete f;++sources;
    }
    Need(list.eof()&&sources==20&&saved==currentSaved&&!rootError,"current source coverage mismatch");
    TFile* old=TFile::Open(oldRoot,"READ");Need(old&&!old->IsZombie()&&!old->TestBit(TFile::kRecovered),"bad old QA ROOT");
    TH2* oldMaps[2];RatioStats metric[2][2];std::ostringstream body;body<<std::setprecision(17)<<"\"currentSources\":"<<sources<<",\"currentSourceInputEvents\":"<<currentInput<<",\"currentSavedEvents\":"<<saved<<",\"currentHe4SignalPairs\":"<<selected<<",\"oldSourceInputEvents\":"<<oldInput<<",\"oldCacheEvents\":"<<oldCache<<",\"oldHe4SignalPairs\":"<<oldSelected<<",\"oldSuccessfulJobs\":"<<oldJobs<<",\"oldExcludedJobs\":"<<oldExcluded<<",\"histograms\":[";
    for(int d=0;d<2;++d) {
      Exact(current[d],aliases[d]);oldMaps[d]=CloneMap(old,std::string("hDphiDeta_selected_4He_")+daughters[d]+"_HelixFixedGlobalEta",std::string("old_")+daughters[d]);
      TH2* h[2]={oldMaps[d],current[d]};
      for(int sample=0;sample<2;++sample) {
        TH2* map=h[sample];Need(map->GetEntries()>0,"empty angular map");
        if(d||sample)body<<',';
        const double regular=map->Integral(1,map->GetNbinsX(),1,map->GetNbinsY()),all=map->Integral(0,map->GetNbinsX()+1,0,map->GetNbinsY()+1);
        Need(all==map->GetEntries(),"raw map Entries differ from all-bin integral");
        body<<"{\"sample\":"<<Json(sample?"current":"old")<<",\"daughter\":"<<Json(daughters[d])<<",\"entries\":"<<map->GetEntries()<<",\"regularIntegral\":"<<regular<<",\"flowIntegral\":"<<all-regular<<",\"nativeBinsX\":"<<map->GetNbinsX()<<",\"nativeBinsY\":"<<map->GetNbinsY()<<",\"densityRatios\":[";
        const double widths[]={phiCut/4,phiCut/2,phiCut};
        for(int w=0;w<3;++w){if(w)body<<',';body<<"{\"phiHalfWidth\":"<<widths[w]<<",\"etaHalfWidth\":"<<etaCut*widths[w]/phiCut<<",\"measurement\":";RatioStats result=Density(map,widths[w],etaCut*widths[w]/phiCut,ringScale);RatioJson(body,result);body<<'}';if(w==2)metric[sample][d]=result;}
        body<<"]}";
      }
    }
    old->Close();delete old;
    TFile* aggregate=TFile::Open(aggregateRoot,"CREATE");Need(aggregate&&!aggregate->IsZombie(),"cannot create map aggregate");
    for(int d=0;d<2;++d){Need(current[d]->Write()>0&&aliases[d]->Write()>0&&oldMaps[d]->Write()>0,"cannot save summed map");}aggregate->Close();delete aggregate;
    gStyle->SetOptStat(0);gStyle->SetPalette(1);gStyle->SetPaintTextFormat(".2g");
    const std::string scope1=Form("Current: 20 EP8/depth100 period outputs, %lld saved events, %.0f He4 signal SE pairs.",currentSaved,selected);
    const std::string scope2=Form("Old QA: %d successful jobs (%d excluded), %lld cache frames, %lld signal pairs; cohorts/selections differ.",oldJobs,oldExcluded,oldCache,oldSelected);
    TCanvas cover("merging_cover","Merging comparison",1200,1000);cover.Print(pdf+"[");
    std::vector<TH2*> drawn;
    TCanvas raw("current_raw","Current uncut raw SE maps",1200,1000);PrepareCanvas(raw,"Current uncut #Lambda-^{4}He angular maps; fixed r = 1.4 m",scope1,"Global rigidity helix and each event B; daughter-position phi difference; global-momentum eta difference.");
    for(int d=0;d<2;++d)for(int zoom=0;zoom<2;++zoom){TH2* h=dynamic_cast<TH2*>(current[d]->Clone(Form("raw_%d_%d",d,zoom)));h->SetDirectory(0);Need(h->GetNbinsX()%rawRebin==0&&h->GetNbinsY()%rawRebin==0,"raw rebin does not divide axes");h->Rebin2D(rawRebin,rawRebin);drawn.push_back(h);raw.cd();const double limit=zoom?zoomRange:fullRange;Panel(h,1+d*2+zoom,Form("%s - ^{4}He, %s",d?"#pi^{-}":"p",zoom?"origin zoom":"full native range"),limit,RangeMax(h,limit),phiCut,etaCut,d==0,false,metric[1][d]);}
    raw.cd();raw.Print(pdf);raw.Print(TString(outputPrefix)+"_current_raw.png");
    TH2* normalized[2][2];
    for(int sample=0;sample<2;++sample)for(int d=0;d<2;++d){TH2* source=sample?current[d]:oldMaps[d];TH2* h=dynamic_cast<TH2*>(source->Clone(Form("shape_%d_%d",sample,d)));h->SetDirectory(0);Need(h->GetNbinsX()%commonBins==0&&h->GetNbinsY()%commonBins==0,"common display bins do not divide native axes");h->Rebin2D(h->GetNbinsX()/commonBins,h->GetNbinsY()/commonBins);h->Scale(1./source->GetEntries());normalized[sample][d]=h;drawn.push_back(h);}
    for(int d=0;d<2;++d)Need(Axis(normalized[0][d]->GetXaxis(),normalized[1][d]->GetXaxis())&&Axis(normalized[0][d]->GetYaxis(),normalized[1][d]->GetYaxis()),"old/new angular ranges disagree");
    for(int zoom=0;zoom<2;++zoom){TCanvas comparison(Form("comparison_%d",zoom),"Old/current angular shape comparison",1200,1000);PrepareCanvas(comparison,zoom?"Old/current shape comparison: origin zoom":"Old/current shape comparison: full range",scope1,scope2);const double limit=zoom?zoomRange:fullRange;
      for(int d=0;d<2;++d){double maximum=std::max(RangeMax(normalized[0][d],limit),RangeMax(normalized[1][d],limit));for(int sample=0;sample<2;++sample){comparison.cd();Panel(normalized[sample][d],1+d*2+sample,Form("%s: %s - ^{4}He",sample?"Current EP8":"Old QA",d?"#pi^{-}":"p"),limit,maximum,phiCut,etaCut,d==0,true,metric[sample][d]);}}
      comparison.cd();comparison.Print(pdf);comparison.Print(TString(outputPrefix)+(zoom?"_comparison_zoom.png":"_comparison_full.png"));}
    cover.Print(pdf+"]");
    std::ofstream report(jsonPath);Need(bool(report),"cannot create metrics report");
    report<<"{\"passed\":true,"<<body.str()<<"],\"bareEventFieldAliasesExactInEverySourceAndSum\":true,\"cutTarget\":\"proton-4He only; pion rectangle is reference\",\"cutPhiHalfWidth\":"<<phiCut<<",\"cutEtaHalfWidth\":"<<etaCut<<",\"adjacentRingScale\":"<<ringScale<<",\"currentRawDisplayRebin\":"<<rawRebin<<",\"commonDisplayBins\":"<<commonBins<<",\"densityMethod\":\"bin-overlap-weighted rectangle and adjacent rectangular ring; count-only error with shared-bin covariance; not exact event-cut counts\",\"comparisonScope\":\"raw signal SE, all kstar and centrality; current EP-valid cohort differs from old isolated QA cohort; no matched-cohort closure or efficiency inference\",\"geometry\":\"STAR global raw-rigidity helix at r=1.4m; own event magnetic field; first forward intersection inside |z|<=210cm; delta phi from positions, delta eta from global momentum\",\"normalizedDisplay\":\"divide by native map Entries including angular flow, then common 40-bin rebin; same z scale for old/current per daughter\",\"pages\":3,\"fullEventReplayPerformed\":false,\"pdf\":"<<Json(pdf.Data())<<",\"aggregateRoot\":"<<Json(aggregateRoot)<<"}\n";
    report.close();Need(bool(report)&&!rootError,"plot/report ROOT error");
    for(size_t i=0;i<drawn.size();++i)delete drawn[i];for(int d=0;d<2;++d){delete current[d];delete aliases[d];delete oldMaps[d];}
    SetErrorHandler(DefaultErrorHandler);std::cout<<"MERGING_CURRENT_OLD_PLOT_PASS pages=3 sources="<<sources<<std::endl;
  } catch(const std::exception& e){SetErrorHandler(DefaultErrorHandler);std::cerr<<"MERGING_CURRENT_OLD_PLOT_FAILED "<<e.what()<<std::endl;gSystem->Exit(1);}
}
