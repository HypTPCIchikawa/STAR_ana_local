#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TParameter.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>

static void MixingQaHist(std::ostream& o,TH1* h){
 if(!h){o<<"null";return;}
 o<<"{\"class\":\""<<h->ClassName()<<"\",\"entries\":"<<h->GetEntries()<<",\"nx\":"<<h->GetNbinsX()<<",\"ny\":"<<h->GetNbinsY()<<",\"xMin\":"<<h->GetXaxis()->GetXmin()<<",\"xMax\":"<<h->GetXaxis()->GetXmax()<<",\"bins\":[";
 bool comma=false;
 for(int ix=0;ix<=h->GetNbinsX()+1;++ix){
  const int firstY=h->GetDimension()==2?0:1,lastY=h->GetDimension()==2?h->GetNbinsY()+1:1;
  for(int iy=firstY;iy<=lastY;++iy){
   const double n=h->GetDimension()==2?h->GetBinContent(ix,iy):h->GetBinContent(ix);
   if(n==0)continue;
   if(comma)o<<',';comma=true;
   o<<'['<<ix<<','<<iy<<','<<h->GetXaxis()->GetBinLowEdge(ix)<<','<<h->GetXaxis()->GetBinUpEdge(ix)<<','<<n<<']';
  }
 }
 o<<"]}";
}

void ProbeMixingQa(const char* input,const char* output){
 try{
  TFile f(input,"READ");if(f.IsZombie())throw std::runtime_error("cannot open source ROOT");
  std::ofstream o(output);if(!o)throw std::runtime_error("cannot create QA output");o<<std::setprecision(17);
  o<<"{\"histograms\":{";const char* names[]={"hVz","hCentralityVsVz","hMixSamplerQA","hN","hCentrality","hKstarSEVsCent_lambda_he4_signal","hKstarMEVsCent_lambda_he4_signal","hKstarSEVsCent_lambda_he3_signal","hKstarMEVsCent_lambda_he3_signal","hKstarSEVsCent_lambda_he3_as_he4_signal","hKstarMEVsCent_lambda_he3_as_he4_signal"};
  for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i){if(i)o<<',';o<<'\"'<<names[i]<<"\":";MixingQaHist(o,dynamic_cast<TH1*>(f.Get(names[i])));}
  o<<"},\"parameters\":{";bool comma=false;TIter next(f.GetListOfKeys());TKey* key=0;
  while((key=static_cast<TKey*>(next()))){
   const std::string name=key->GetName(),cl=key->GetClassName();
   if(name.find("Vz")!=std::string::npos||name.find("Clamp")!=std::string::npos||name.find("clamp")!=std::string::npos||cl.find("TParameter")!=std::string::npos)std::cout<<name<<" "<<cl<<std::endl;
   if(cl.find("TParameter")!=0)continue;
   TObject* obj=f.Get(name.c_str());TParameter<Long64_t>* l=dynamic_cast<TParameter<Long64_t>*>(obj);TParameter<double>* d=dynamic_cast<TParameter<double>*>(obj);
   if(!l&&!d)continue;if(comma)o<<',';comma=true;o<<'\"'<<name<<"\":";if(l)o<<l->GetVal();else o<<d->GetVal();
  }
  o<<"}}\n";o.close();if(!o)throw std::runtime_error("QA output write failure");
 }catch(const std::exception& e){std::cerr<<"MIXING_QA_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
