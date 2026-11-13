#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TArrayD.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <set>
#include <vector>
#include <cmath>

namespace HighKExport {
void String(std::ostream& o,const std::string& s){
 o<<'"';for(std::string::const_iterator i=s.begin();i!=s.end();++i){
  const unsigned char c=*i;
  if(c=='"')o<<"\\\"";else if(c=='\\')o<<"\\\\";
  else if(c=='\n')o<<"\\n";else if(c=='\r')o<<"\\r";else if(c=='\t')o<<"\\t";
  else if(c<32){std::ostringstream t;t<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c);o<<t.str();}
  else o<<*i;
 }o<<'"';
}
void Number(std::ostream& o,double v){if(!TMath::Finite(v))throw std::runtime_error("nonfinite histogram value");o<<v;}
void Edges(std::ostream& o,const TAxis* a){o<<'[';for(int i=1;i<=a->GetNbins()+1;++i){if(i>1)o<<',';Number(o,a->GetBinLowEdge(i));}o<<']';}
void Metadata(std::ostream& o,TH1* h,const std::string& path){
 o<<"{\"name\":";String(o,h->GetName());o<<",\"path\":";String(o,path);o<<",\"title\":";String(o,h->GetTitle());
 o<<",\"class\":";String(o,h->ClassName());o<<",\"dimension\":"<<h->GetDimension()<<",\"entries\":";Number(o,h->GetEntries());
 o<<",\"sumw2N\":"<<h->GetSumw2N()<<",\"xEdges\":";Edges(o,h->GetXaxis());o<<",\"xTitle\":";String(o,h->GetXaxis()->GetTitle());
 if(h->GetDimension()>=2){o<<",\"yEdges\":";Edges(o,h->GetYaxis());o<<",\"yTitle\":";String(o,h->GetYaxis()->GetTitle());}
 if(h->GetDimension()>=3){o<<",\"zEdges\":";Edges(o,h->GetZaxis());o<<",\"zTitle\":";String(o,h->GetZaxis()->GetTitle());}
}
void Payload(std::ostream& o,TH1* h,const std::string& path){
 Metadata(o,h,path);
 const int nx=h->GetNbinsX(),ny=h->GetDimension()==2?h->GetNbinsY():-1;
 o<<",\"shape\":[";if(ny>=0)o<<ny+2<<',';o<<nx+2<<"],\"content\":[";
 const int n=(nx+2)*(ny>=0?ny+2:1);
 for(int i=0;i<n;++i){if(i)o<<',';Number(o,h->GetBinContent(i));}
 o<<"],\"error2\":[";
 for(int i=0;i<n;++i){if(i)o<<',';const double e=h->GetSumw2N()?h->GetSumw2()->At(i):h->GetBinError(i)*h->GetBinError(i);if(e<0)throw std::runtime_error("negative bin variance");Number(o,e);}
 o<<"]}";
}
void Walk(TDirectory* d,const std::string& prefix,std::ostream& o,bool& comma,std::vector<std::string>& meta,int& count){
 std::set<std::string> seen;TIter next(d->GetListOfKeys());TKey* k=0;
 while((k=static_cast<TKey*>(next()))){
  const std::string name=k->GetName(),path=prefix+name;
  if(seen.count(name))continue;seen.insert(name);
  TClass* cl=gROOT->GetClass(k->GetClassName());if(!cl)continue;
  if(cl->InheritsFrom(TDirectory::Class())){TDirectory* sub=dynamic_cast<TDirectory*>(d->Get(name.c_str()));if(!sub)throw std::runtime_error("cannot read subdirectory");Walk(sub,path+"/",o,comma,meta,count);continue;}
  if(name.find("hKstar")!=0||!cl->InheritsFrom(TH1::Class()))continue;
  TH1* h=dynamic_cast<TH1*>(d->Get(name.c_str()));if(!h)throw std::runtime_error("cannot read histogram "+path);
  if(h->GetDimension()>2){std::ostringstream m;m<<std::setprecision(17);Metadata(m,h,path);m<<'}';meta.push_back(m.str());delete h;continue;}
  if(comma)o<<',';comma=true;String(o,path);o<<':';Payload(o,h,path);++count;
  if(count%100==0)std::cout<<"exported "<<count<<" histograms"<<std::endl;
  delete h;
 }
}
}
void ExportHighKHistograms(const char* input,const char* output){
 try{
  TFile f(input,"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered))throw std::runtime_error("invalid or recovered source ROOT");
  std::ofstream o(output);if(!o)throw std::runtime_error("cannot create output");o<<std::setprecision(17);
  o<<"{\"schemaVersion\":1,\"input\":";HighKExport::String(o,input);o<<",\"rootVersion\":";HighKExport::String(o,gROOT->GetVersion());
  o<<",\"uuid\":";HighKExport::String(o,f.GetUUID().AsString());o<<",\"storageOrder\":\"ROOT global bins, x fastest, flow included\",\"histograms\":{";
  bool comma=false;int count=0;std::vector<std::string> meta;HighKExport::Walk(&f,"",o,comma,meta,count);
  if(count==0)throw std::runtime_error("no selected histograms");
  o<<"},\"metadataOnly\":[";for(unsigned i=0;i<meta.size();++i){if(i)o<<',';o<<meta[i];}o<<"],\"histogramCount\":"<<count<<"}\n";
  o.close();if(!o)throw std::runtime_error("output write failed");
  std::cout<<"HIGHK_EXPORT_OK histograms="<<count<<" metadataOnly="<<meta.size()<<std::endl;
 }catch(const std::exception& e){std::cerr<<"HIGHK_EXPORT_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
