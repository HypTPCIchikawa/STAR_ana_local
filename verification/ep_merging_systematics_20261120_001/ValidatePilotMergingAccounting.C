// Read eight existing real-data pilots; no event replay or shared build.
int ValidatePilotMergingAccounting(const char* study);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TClass.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PilotMergingAccounting {
unsigned checks=0;
void Need(bool ok,const std::string& why){++checks;if(!ok)throw std::runtime_error(why);}
std::string Number(int n){std::ostringstream s;s<<n;return s.str();}
TH1* Hist(TFile& f,const std::string& k){
  TH1* h=dynamic_cast<TH1*>(f.Get(k.c_str()));
  if(!h)h=dynamic_cast<TH1*>(f.Get(("true/"+k).c_str()));
  if(!h)h=dynamic_cast<TH1*>(f.Get(("mix/"+k).c_str()));
  Need(h!=0,std::string(f.GetName())+": missing "+k);return h;
}
int Cells(const TH1& h){return (h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1)*(h.GetDimension()>2?h.GetNbinsZ()+2:1);}
double Sum(const TH1& h){double n=0.;for(int i=0;i<Cells(h);++i)n+=h.GetBinContent(i);return n;}
void Unit(const TH1& h){
  Need(Sum(h)==h.GetEntries(),std::string(h.GetName())+": total including flows equals entries");
  for(int i=0;i<Cells(h);++i){const double n=h.GetBinContent(i),v=std::pow(h.GetBinError(i),2);
    Need(std::isfinite(n)&&n>=0.&&n==std::floor(n),"nonnegative integer counts");
    Need(std::fabs(v-n)<1.e-7*(1.+n),"unit-count variances");}
}
void Same(const TH1& a,const TH1& b,const std::string& why){
  Need(a.GetDimension()==b.GetDimension()&&Cells(a)==Cells(b)&&a.GetEntries()==b.GetEntries(),why+": shape/entries");
  for(int i=0;i<Cells(a);++i)Need(a.GetBinContent(i)==b.GetBinContent(i)&&a.GetBinError(i)==b.GetBinError(i),why+": exact bins/errors incl flow");
}
double Q(TH2& h,int x,int mode){return h.GetBinContent(h.FindBin(x,mode));}
double Flow(const TH2& h){return Sum(h)-h.Integral(1,h.GetNbinsX(),1,h.GetNbinsY());}
double Interior(const TH2& h){double n=0.;
  // Boundary-touching bins are excluded: finite bins cannot count a strict cut exactly.
  for(int x=1;x<=h.GetNbinsX();++x)for(int y=1;y<=h.GetNbinsY();++y)
    if(h.GetXaxis()->GetBinLowEdge(x)>-.05&&h.GetXaxis()->GetBinUpEdge(x)<.05&&
       h.GetYaxis()->GetBinLowEdge(y)>-.05&&h.GetYaxis()->GetBinUpEdge(y)<.05)n+=h.GetBinContent(x,y);
  return n;}
std::string MapKey(const char* d,bool mixed,bool after){return std::string("hDphiDetaPair_")+d+"_"+(mixed?"Mixed_":"")+"4He_"+(after?"after":"before");}
void He3(TDirectory& d,const std::string& prefix,std::map<std::string,TH1*>& out){
  TIter next(d.GetListOfKeys());TKey* k=0;
  while((k=static_cast<TKey*>(next()))){const std::string path=prefix+k->GetName();TClass* type=TClass::GetClass(k->GetClassName());
    if(type&&type->InheritsFrom("TDirectory")){TDirectory* sub=dynamic_cast<TDirectory*>(d.Get(k->GetName()));Need(sub!=0,"read directory");He3(*sub,path+"/",out);}
    else if(type&&type->InheritsFrom("TH1")&&(path.find("he3")!=std::string::npos||path.find("3He")!=std::string::npos)){
      TH1* h=dynamic_cast<TH1*>(d.Get(k->GetName()));Need(h!=0,"read He3 histogram");out[path]=h;}}
}
void Check(TFile& f,bool enabled,std::ostringstream& r){
  TH2* qa=dynamic_cast<TH2*>(Hist(f,"hDphiDetaPairCutQA"));Need(qa&&qa->GetNbinsX()==12&&qa->GetNbinsY()==2,"12x2 QA shape");Unit(*qa);Need(Flow(*qa)==0.,"QA has no flow");r<<"[";
  for(int m=0;m<2;++m){const bool mixed=m!=0;
    const double all=Q(*qa,0,m),pass=Q(*qa,1,m),rej=Q(*qa,2,m),vp=Q(*qa,3,m),ip=Q(*qa,4,m),vpi=Q(*qa,5,m),ipi=Q(*qa,6,m),
      active=Q(*qa,7,m),sig=Q(*qa,8,m),sp=Q(*qa,9,m),sr=Q(*qa,10,m),si=Q(*qa,11,m);
    Need(all>0.&&all==pass+rej&&all==vp+ip&&all==vpi+ipi,"all-pair acceptance/validity accounting");
    Need(sig==sp+sr&&sig<=all&&sr<=rej&&si<=sp,"signal/unreachable-keep accounting");
    Need(active==(enabled?all:0.),"enabled counter matches case");Need(enabled||(rej==0.&&sr==0.),"off never rejects");
    const std::string stem=std::string("hKstar")+(mixed?"ME":"SE")+"_lambda_he4";
    Need(Hist(f,stem)->GetEntries()==pass&&Hist(f,stem+"_signal")->GetEntries()==sp,"QA pass equals physical kstar entries");
    if(m)r<<",";r<<"{\"mode\":\""<<(mixed?"ME":"SE")<<"\",\"considered\":"<<all<<",\"passed\":"<<pass<<",\"rejected\":"<<rej
      <<",\"signalConsidered\":"<<sig<<",\"signalPassed\":"<<sp<<",\"signalRejected\":"<<sr<<",\"signalProtonUnreachableKept\":"<<si<<",\"maps\":[";
    const char* daughters[]={"proton","pion"};
    for(int d=0;d<2;++d){TH2* b=dynamic_cast<TH2*>(Hist(f,MapKey(daughters[d],mixed,false))),*a=dynamic_cast<TH2*>(Hist(f,MapKey(daughters[d],mixed,true)));
      Need(b&&a&&Cells(*b)==Cells(*a),"map shapes");Unit(*b);Unit(*a);const double loss=b->GetEntries()-a->GetEntries();Need(loss>=0.&&loss<=sr,"map loss bounded by signal rejection");
      if(d==0){Need(b->GetEntries()==sig-si&&a->GetEntries()==sp-si&&loss==sr,"proton map/invalid-kept/rejection accounting");
        Need(Flow(*b)==Flow(*a),"proton cut does not remove map flow");Need(!enabled||Interior(*a)==0.,"cut-on proton after has zero strict interior bins");}
      if(!enabled)Same(*b,*a,"off before/after invariance");
      for(int i=0;i<Cells(*b);++i)Need(a->GetBinContent(i)<=b->GetBinContent(i),"after is subset of before");
      if(d)r<<",";r<<"{\"daughter\":\""<<daughters[d]<<"\",\"beforeEntries\":"<<b->GetEntries()<<",\"afterEntries\":"<<a->GetEntries()
        <<",\"beforeFlow\":"<<Flow(*b)<<",\"afterFlow\":"<<Flow(*a)<<",\"afterStrictInteriorBins\":"<<Interior(*a)<<"}";}
    r<<"]}";}r<<"]";
}
}

int ValidatePilotMergingAccounting(const char* study){using namespace PilotMergingAccounting;
  try{checks=0;Need(study&&*study,"explicit study");const std::string base=study,result=base+"/pilot_merging_accounting.json";
    Need(gSystem->AccessPathName(result.c_str()),"refuse gate overwrite");const int bins[]={1,8,16,32};std::ostringstream body;body.precision(17);body<<"[";
    unsigned he3=0;bool nonemptySE=false,nonemptyME=false;
    for(int e=0;e<4;++e){const std::string stem=base+"/pilots/ep"+Number(bins[e])+"_merging";
      TFile off((stem+"Off.root").c_str(),"READ"),on((stem+"On.root").c_str(),"READ");Need(!off.IsZombie()&&!on.IsZombie(),"pilot ROOTs open");
      std::ostringstream oldReport,newReport;oldReport.precision(17);newReport.precision(17);Check(off,false,oldReport);Check(on,true,newReport);
      TH2* oq=dynamic_cast<TH2*>(Hist(off,"hDphiDetaPairCutQA")),*nq=dynamic_cast<TH2*>(Hist(on,"hDphiDetaPairCutQA"));
      for(int m=0;m<2;++m){const bool mixed=m!=0;
        for(int x=0;x<12;++x)if(x==0||x==3||x==4||x==5||x==6||x==8||x==11)Need(Q(*oq,x,m)==Q(*nq,x,m),"on/off pre-cut/invalid accounting unchanged");
        const double rej=Q(*nq,10,m);const std::string key=std::string("hKstar")+(mixed?"ME":"SE")+"_lambda_he4_signal";
        Need(Hist(off,key)->GetEntries()-Hist(on,key)->GetEntries()==rej,"off-on signal CF difference equals signalRejected");
        if(m)nonemptyME=nonemptyME||rej>0.;else nonemptySE=nonemptySE||rej>0.;
        const char* ds[]={"proton","pion"};for(int d=0;d<2;++d)Same(*Hist(off,MapKey(ds[d],mixed,false)),*Hist(on,MapKey(ds[d],mixed,false)),"on/off before map invariance");}
      std::map<std::string,TH1*> oldHe3,newHe3;He3(off,"",oldHe3);He3(on,"",newHe3);Need(!oldHe3.empty()&&oldHe3.size()==newHe3.size(),"matching nonempty He3 hist set");
      for(std::map<std::string,TH1*>::const_iterator h=oldHe3.begin();h!=oldHe3.end();++h){Need(newHe3.count(h->first)==1,"matching He3 key");Same(*h->second,*newHe3[h->first],"He3 control invariance");++he3;}
      if(e)body<<",";body<<"{\"epBins\":"<<bins[e]<<",\"offPath\":\""<<stem<<"Off.root\",\"onPath\":\""<<stem<<"On.root\",\"off\":"<<oldReport.str()<<",\"on\":"<<newReport.str()<<"}";}
    body<<"]";std::ofstream output(result.c_str());Need(bool(output),"create accounting gate");
    output<<"{\"passed\":true,\"checks\":"<<checks<<",\"he3HistogramsCompared\":"<<he3<<",\"nonemptySignalRejectionSE\":"<<(nonemptySE?"true":"false")
      <<",\"nonemptySignalRejectionME\":"<<(nonemptyME?"true":"false")<<",\"rectangleNote\":\"Boundary-touching bins are excluded from the interior check. Bin integrals are not exact strict-cut counts; signalRejected is authoritative. Real-pilot zero rejection is allowed; actual-Legacy synthetic tests separately verify nonempty SE and both ME directions.\",\"cases\":"<<body.str()<<"}\n";
    output.close();Need(bool(output),"write accounting gate");std::cout<<"PASS pilot merging accounting checks="<<checks<<" He3 histograms="<<he3<<" nonempty signal rejection SE="<<nonemptySE<<" ME="<<nonemptyME<<std::endl;return 0;
  }catch(const std::exception& error){std::cerr<<"FAIL pilot merging accounting: "<<error.what()<<std::endl;gSystem->Exit(1);return 1;}}
#endif
