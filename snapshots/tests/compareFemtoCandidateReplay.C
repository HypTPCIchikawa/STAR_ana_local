#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "TTree.h"
#include "TH1.h"
#include "THnBase.h"
#include "TAxis.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FemtoReplayComparison {
void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
std::string Number(Long64_t value){std::ostringstream s;s<<value;return s.str();}
std::string Quote(const std::string& s){std::string out="\"";for(size_t i=0;i<s.size();++i){const char c=s[i];if(c=='\\'||c=='\"')out+='\\';if(c=='\n')out+="\\n";else if(c=='\r')out+="\\r";else if(c=='\t')out+="\\t";else out+=c;}return out+'\"';}
bool Mixed(const std::string& path,const std::string& leaf){return path.compare(0,4,"mix/")==0||path.find("/mix/")!=std::string::npos||leaf.compare(0,8,"hKstarME")==0||leaf.find("Mixed_")!=std::string::npos;}
bool Selected(const std::string& path,bool includePrePid,bool sameEventOnly){
 const size_t slash=path.rfind('/');const std::string leaf=slash==std::string::npos?path:path.substr(slash+1);
 if(leaf=="hNuclearPidPairJoint")return includePrePid&&!sameEventOnly;
 if(leaf=="hMixSamplerQA")return !sameEventOnly;
 if(leaf.compare(0,6,"hKstar")!=0&&leaf.compare(0,5,"hQlab")!=0)return false;
 if(!includePrePid&&leaf.find("PrePid")!=std::string::npos)return false;
 return !sameEventOnly||!Mixed(path,leaf);
}
void Keys(TDirectory& d,const std::string& prefix,bool includePrePid,bool sameEventOnly,std::set<std::string>& names){
 TIter next(d.GetListOfKeys());TKey* key=0;
 while((key=static_cast<TKey*>(next()))){const std::string name=key->GetName(),cl=key->GetClassName(),path=prefix+name;
  if(cl=="TDirectoryFile"||cl=="TDirectory"){TDirectory* sub=d.GetDirectory(name.c_str());Require(sub!=0,"cannot read directory "+path);Keys(*sub,path+"/",includePrePid,sameEventOnly,names);}
  else if(Selected(path,includePrePid,sameEventOnly))Require(names.insert(path).second,"duplicate key cycle: "+path);
 }
}
void Axis(const TAxis& a,const TAxis& b,const std::string& path){
 Require(a.GetNbins()==b.GetNbins()&&a.GetXmin()==b.GetXmin()&&a.GetXmax()==b.GetXmax(),"axis dimensions/range: "+path);
 for(int i=1;i<=a.GetNbins()+1;++i)Require(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),"axis edge: "+path+" bin "+Number(i));
}
Int_t DenseCells(const TH1& h){
 Require(h.GetDimension()>=1&&h.GetDimension()<=3,"dense histogram dimension must be 1-3");
 const Long64_t maximum=std::numeric_limits<Int_t>::max();
 Long64_t cells=static_cast<Long64_t>(h.GetNbinsX())+2;
 Require(cells>0&&cells<=maximum,"dense X cell count exceeds Int_t indexing");
 if(h.GetDimension()>1){const Long64_t factor=static_cast<Long64_t>(h.GetNbinsY())+2;Require(factor>0&&cells<=maximum/factor,"dense XY cell count exceeds Int_t indexing");cells*=factor;}
 if(h.GetDimension()>2){const Long64_t factor=static_cast<Long64_t>(h.GetNbinsZ())+2;Require(factor>0&&cells<=maximum/factor,"dense XYZ cell count exceeds Int_t indexing");cells*=factor;}
 return static_cast<Int_t>(cells);
}
Long64_t Dense(TH1& a,TH1& b,const std::string& path){
 Require(a.GetDimension()==b.GetDimension(),"dimension: "+path);
 Axis(*a.GetXaxis(),*b.GetXaxis(),path+" X");if(a.GetDimension()>1)Axis(*a.GetYaxis(),*b.GetYaxis(),path+" Y");if(a.GetDimension()>2)Axis(*a.GetZaxis(),*b.GetZaxis(),path+" Z");
 Require(a.GetEntries()==b.GetEntries(),"Entries: "+path);
 const Int_t cells=DenseCells(a);Require(cells==DenseCells(b)&&a.GetSumw2N()==b.GetSumw2N(),"cells/Sumw2 shape: "+path);
 for(int i=0;i<cells;++i){Require(a.GetBinContent(i)==b.GetBinContent(i),"bin content: "+path+" cell "+Number(i));Require(a.GetBinError(i)==b.GetBinError(i),"bin error: "+path+" cell "+Number(i));}
 for(int i=0;i<a.GetSumw2N();++i)Require(a.GetSumw2()->At(i)==b.GetSumw2()->At(i),"Sumw2: "+path+" cell "+Number(i));
 Double_t sa[13]={0},sb[13]={0};a.GetStats(sa);b.GetStats(sb);for(int i=0;i<13;++i)Require(sa[i]==sb[i],"statistics: "+path+" field "+Number(i));
 return cells;
}
Long64_t Sparse(THnBase& a,THnBase& b,const std::string& path){
 Require(a.GetNdimensions()==b.GetNdimensions()&&a.GetEntries()==b.GetEntries()&&a.GetNbins()==b.GetNbins(),"sparse shape/Entries/cells: "+path);
 for(int i=0;i<a.GetNdimensions();++i)Axis(*a.GetAxis(i),*b.GetAxis(i),path+" axis "+Number(i));
 std::vector<Int_t> coordinate(a.GetNdimensions());
 for(Long64_t i=0;i<a.GetNbins();++i){const double value=a.GetBinContent(i,&coordinate[0]);const Long64_t j=b.GetBin(&coordinate[0],kFALSE);Require(j>=0,"sparse coordinate missing: "+path);Require(value==b.GetBinContent(j),"sparse content: "+path+" cell "+Number(i));Require(a.GetBinError2(i)==b.GetBinError2(j),"sparse error2: "+path+" cell "+Number(i));}
 return a.GetNbins();
}
}

// Reconstruction and helix/merging QA require PicoDst and are outside replay closure.
// Core pair histograms include signal/sidebands, fine centrality, mT and mass mappings.
void compareFemtoCandidateReplay(const char* reference,const char* replay,const char* report,
                                Bool_t includePrePidQa,const char* comparisonMode){
 using namespace FemtoReplayComparison;
 try{
  Require(reference&&*reference&&replay&&*replay&&report&&*report&&comparisonMode,"explicit arguments required");
  const std::string mode=comparisonMode;Require(mode=="closure"||mode=="sameEventOnly","comparisonMode must be closure or sameEventOnly");const bool seOnly=mode=="sameEventOnly";
  Require(gSystem->AccessPathName(report),"refusing existing comparison report");
  TFile a(reference,"READ"),b(replay,"READ");Require(!a.IsZombie()&&!b.IsZombie()&&!a.TestBit(TFile::kRecovered)&&!b.TestBit(TFile::kRecovered),"valid closed ROOT inputs required");
  TNamed* ca=dynamic_cast<TNamed*>(a.Get("FemtoLambdaRunStatus"));TNamed* cb=dynamic_cast<TNamed*>(b.Get("FemtoLambdaRunStatus"));Require(ca&&cb&&std::string(ca->GetTitle())=="completed"&&std::string(cb->GetTitle())=="completed","both runs must be completed");
  std::set<std::string> ka,kb;Keys(a,"",includePrePidQa,seOnly,ka);Keys(b,"",includePrePidQa,seOnly,kb);Require(!ka.empty()&&ka==kb,"selected physics key sets must match");
  TH1* hse=dynamic_cast<TH1*>(a.Get("hKstarSE_lambda_he4_signal"));Require(hse&&hse->GetEntries()>0,"pilot must exercise real He4 signal SE pairs");const double seEntries=hse->GetEntries();double meEntries=0;
  if(!seOnly){TH1* hme=dynamic_cast<TH1*>(a.Get("hKstarME_lambda_he4_signal"));Require(hme&&hme->GetEntries()>0,"pilot must exercise real He4 signal ME pairs");meEntries=hme->GetEntries();}
  TTree* storeA=dynamic_cast<TTree*>(a.Get("femtoCandidateEvents"));TTree* storeB=dynamic_cast<TTree*>(b.Get("femtoCandidateEvents"));
  const Long64_t rowsA=storeA?storeA->GetEntries():-1,rowsB=storeB?storeB->GetEntries():-1;
  bool meChanged=false;Long64_t differentMeCells=0;double alternateMeEntries=0;
  if(seOnly){TH1* ar=dynamic_cast<TH1*>(a.Get("hKstarME_lambda_he4_signal"));TH1* br=dynamic_cast<TH1*>(b.Get("hKstarME_lambda_he4_signal"));Require(ar&&br,"ME probe histograms required");const Int_t meCells=DenseCells(*ar);Require(meCells==DenseCells(*br),"ME probe cell counts");meEntries=ar->GetEntries();alternateMeEntries=br->GetEntries();meChanged=meEntries!=alternateMeEntries;for(int i=0;i<meCells;++i)if(ar->GetBinContent(i)!=br->GetBinContent(i)||ar->GetBinError(i)!=br->GetBinError(i)){++differentMeCells;meChanged=true;}}
  Long64_t denseCells=0,sparseCells=0;unsigned denseHistograms=0,sparseHistograms=0;
  for(std::set<std::string>::const_iterator i=ka.begin();i!=ka.end();++i){TObject* oa=a.Get(i->c_str());TObject* ob=b.Get(i->c_str());Require(oa&&ob&&std::string(oa->ClassName())==ob->ClassName(),"object type: "+*i);
   TH1* ha=dynamic_cast<TH1*>(oa);TH1* hb=dynamic_cast<TH1*>(ob);THnBase* sa=dynamic_cast<THnBase*>(oa);THnBase* sb=dynamic_cast<THnBase*>(ob);
   if(ha&&hb){denseCells+=Dense(*ha,*hb,*i);++denseHistograms;}else if(sa&&sb){sparseCells+=Sparse(*sa,*sb,*i);++sparseHistograms;}else throw std::runtime_error("selected key is not histogram: "+*i);
   delete oa;delete ob;
  }
  std::ofstream out(report);Require(bool(out),"cannot create comparison report");out<<std::setprecision(17)<<"{\n  \"phase\": \"passed\",\n  \"comparisonMode\": "<<Quote(mode)<<",\n  \"reference\": "<<Quote(reference)<<",\n  \"replay\": "<<Quote(replay)<<",\n  \"includePrePidQa\": "<<(includePrePidQa?"true":"false")<<",\n  \"denseHistograms\": "<<denseHistograms<<",\n  \"sparseHistograms\": "<<sparseHistograms<<",\n  \"denseCellsIncludingFlows\": "<<denseCells<<",\n  \"sparseCells\": "<<sparseCells<<",\n  \"exactContentErrorEntriesAndStatistics\": true,\n  \"referenceHe4SignalSEEntries\": "<<seEntries<<",\n  \"referenceHe4SignalMEEntries\": "<<meEntries<<",\n  \"referenceCandidateStoreRows\": "<<rowsA<<",\n  \"replayCandidateStoreRows\": "<<rowsB<<",\n  \"probeMixedEventChanged\": "<<(meChanged?"true":"false")<<",\n  \"probeMixedEventDifferentCells\": "<<differentMeCells<<",\n  \"probeMixedEventReplayEntries\": "<<alternateMeEntries<<",\n  \"excluded\": \"reconstruction, track/helix merging QA and metadata not used for physical pair closure\"\n}\n";Require(bool(out),"comparison report write failed");
  std::cout<<"PASS candidate replay "<<mode<<": "<<denseHistograms<<" dense, "<<sparseHistograms<<" sparse histograms; exact contents/errors/Entries/statistics"<<std::endl;
 }catch(const std::exception& e){std::cerr<<"CANDIDATE_REPLAY_COMPARISON_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
