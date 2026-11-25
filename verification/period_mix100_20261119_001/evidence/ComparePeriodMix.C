#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <stdexcept>

// Pair/Mixed classification and flow-inclusive comparison follow the existing
// rootfile/femto_lambda_tree_implementation_20261117_001/verify_outputs_v3.C.
// Configuration/job-id objects are reported separately from physical SE histograms.
namespace PeriodMixComparison {
bool rootError=false;
void ErrorHandler(Int_t level,Bool_t abort,const char* where,const char* message) {
  if(level>=kError)rootError=true;DefaultErrorHandler(level,abort,where,message);
}
void Need(bool ok,const std::string& message) { if(!ok)throw std::runtime_error(message); }
std::string Json(const std::string& value) {
  std::ostringstream out;out<<'"';
  for(size_t i=0;i<value.size();++i) {
    unsigned char c=value[i];
    if(c=='"'||c=='\\')out<<'\\'<<char(c);
    else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
    else out<<char(c);
  }
  return out.str()+'"';
}
std::string Base(const std::string& name) { size_t p=name.find_last_of('/');return p==std::string::npos?name:name.substr(p+1); }
bool Starts(const std::string& text,const char* prefix) { return text.compare(0,std::string(prefix).size(),prefix)==0; }
bool Pair(const std::string& path) {
  std::string n=Base(path);return Starts(n,"hKstar")||Starts(n,"hQlab")||Starts(n,"hDphi")||Starts(n,"hNuclearPidPair");
}
bool Mixed(const std::string& path) { return Starts(path,"mix/")||Starts(Base(path),"hKstarME")||Base(path).find("_Mixed_")!=std::string::npos; }
void Inventory(TDirectory* directory,const std::string& prefix,std::map<std::string,TKey*>& se) {
  TIter next(directory->GetListOfKeys());TKey* key=0;std::set<std::string> seen;
  while((key=dynamic_cast<TKey*>(next()))) {
    const std::string name=key->GetName();if(!seen.insert(name).second)continue;
    key=directory->GetKey(name.c_str());Need(key!=0,"missing latest key");
    TClass* type=TClass::GetClass(key->GetClassName());Need(type!=0,"unknown ROOT class");
    const std::string path=prefix+name;
    if(type->InheritsFrom(TDirectory::Class())) {
      TDirectory* child=directory->GetDirectory(name.c_str());Need(child!=0,"cannot read ROOT directory");Inventory(child,path+"/",se);
    } else if(Pair(path)&&!Mixed(path)) {
      Need(type->InheritsFrom(TH1::Class()),"unsupported non-histogram SE object: "+path);
      Need(se.insert(std::make_pair(path,key)).second,"duplicate SE path");
    }
  }
}
bool Axis(const TAxis* a,const TAxis* b) {
  if(a->GetNbins()!=b->GetNbins())return false;
  for(int i=1;i<=a->GetNbins()+1;++i)if(!TMath::Finite(a->GetBinLowEdge(i))||a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;
  for(int i=1;i<=a->GetNbins();++i)if(std::string(a->GetBinLabel(i))!=b->GetBinLabel(i))return false;
  return true;
}
Long64_t Cells(const TH1* h) { Long64_t n=h->GetNbinsX()+2;if(h->GetDimension()>1)n*=h->GetNbinsY()+2;if(h->GetDimension()>2)n*=h->GetNbinsZ()+2;Need(n>0&&n<=2147483647LL,"unsupported histogram shape");return n; }
void Equal(TKey* first,TKey* second,const std::string& name,Long64_t& bins) {
  TH1* a=dynamic_cast<TH1*>(first->ReadObj());TH1* b=dynamic_cast<TH1*>(second->ReadObj());
  Need(a&&b&&!rootError,"cannot read SE histogram: "+name);a->SetDirectory(0);b->SetDirectory(0);
  const int cells=int(Cells(a));bool equal=std::string(a->ClassName())==b->ClassName()&&a->GetDimension()==b->GetDimension()&&cells==Cells(b)&&a->GetEntries()==b->GetEntries()&&TMath::Finite(a->GetEntries())&&a->GetSumw2N()==b->GetSumw2N();
  equal=equal&&Axis(a->GetXaxis(),b->GetXaxis())&&(a->GetDimension()<2||Axis(a->GetYaxis(),b->GetYaxis()))&&(a->GetDimension()<3||Axis(a->GetZaxis(),b->GetZaxis()));
  if(equal)for(int k=0;k<cells;++k) {
    const double x=a->GetBinContent(k),y=b->GetBinContent(k),e=a->GetBinError(k),f=b->GetBinError(k);
    if(!TMath::Finite(x)||!TMath::Finite(y)||!TMath::Finite(e)||!TMath::Finite(f)||e<0||f<0||x!=y||e!=f){equal=false;break;}
    if(a->GetSumw2N()) { const double v=a->GetSumw2()->At(k),w=b->GetSumw2()->At(k);if(!TMath::Finite(v)||!TMath::Finite(w)||v<0||w<0||v!=w){equal=false;break;} }
  }
  delete a;delete b;Need(equal,"SE content/error/Entries/axis mismatch: "+name);bins+=cells;
}
std::string Object(TFile* file,const char* name) {
  TObject* value=file->Get(name);Need(value!=0,std::string("missing metadata: ")+name);
  TObjString* text=dynamic_cast<TObjString*>(value);if(text)return text->GetString().Data();
  TNamed* named=dynamic_cast<TNamed*>(value);Need(named!=0,"unsupported metadata object");return named->GetTitle();
}
Long64_t Parameter(TFile* file,const char* name) { TParameter<Long64_t>* p=dynamic_cast<TParameter<Long64_t>*>(file->Get(name));Need(p!=0,std::string("missing counter: ")+name);return p->GetVal(); }
std::string Trim(const std::string& text) { size_t a=text.find_first_not_of(" \t\r");if(a==std::string::npos)return "";size_t b=text.find_last_not_of(" \t\r");return text.substr(a,b-a+1); }
std::map<std::string,std::string> Scalars(const std::string& text) {
  std::map<std::string,std::string> out;std::istringstream input(text);std::string line;
  while(std::getline(input,line)) { size_t p=line.find(':');if(p==std::string::npos)continue;std::string key=Trim(line.substr(0,p));std::string value=Trim(line.substr(p+1));if(key.empty()||key[0]=='#'||key[0]=='-')continue;Need(out.insert(std::make_pair(key,value)).second,"duplicate scalar in metadata"); }
  return out;
}
void SameLedger(TFile** files,Long64_t expected,std::vector<Long64_t>& starts) {
  const char* fields[]={"inputIndex","runId","eventId","sourceInputIndex","candidateStoreStream","status","cent9","selectedLambda"};
  TTree* trees[3];Int_t values[3][8];int previous=-1;
  for(int f=0;f<3;++f) {
    trees[f]=dynamic_cast<TTree*>(files[f]->Get("eventLedger"));Need(trees[f]&&trees[f]->GetEntries()==expected,"eventLedger prefix count mismatch");trees[f]->SetBranchStatus("*",0);
    for(int j=0;j<8;++j){Need(trees[f]->GetBranch(fields[j])!=0,"missing ledger identity/segment branch");trees[f]->SetBranchStatus(fields[j],1);Need(trees[f]->SetBranchAddress(fields[j],&values[f][j])>=0,"cannot bind ledger");}
  }
  for(Long64_t row=0;row<expected;++row) {
    for(int f=0;f<3;++f)Need(trees[f]->GetEntry(row)>0&&!rootError,"ledger read failed");
    for(int f=1;f<3;++f)for(int j=0;j<8;++j)Need(values[0][j]==values[f][j],"inputs differ: eventLedger row identity/payload/segment mismatch");
    Need(values[0][0]==row+1&&values[0][5]==4,"ledger index/completion status mismatch");
    if(values[0][4]!=previous){Need(values[0][4]==previous+1,"segment stream is not contiguous");starts.push_back(row+1);previous=values[0][4];}
  }
  for(int f=0;f<3;++f){trees[f]->ResetBranchAddresses();trees[f]->SetBranchStatus("*",1);}
  Need(starts.size()>=3,"pilot has not crossed the first two segment boundaries");
}
} // namespace PeriodMixComparison

// Input order: B=continuous100, A=reset100, C=continuous20.
void ComparePeriodMix(const char* pathB,const char* pathA,const char* pathC,
                      const char* reportPath,Long64_t expectedEvents=400000) {
  using namespace PeriodMixComparison;
  rootError=false;TH1::AddDirectory(kFALSE);SetErrorHandler(ErrorHandler);
  if(!reportPath||!*reportPath||!gSystem->AccessPathName(reportPath)){std::cerr<<"Refusing missing/existing report path"<<std::endl;gSystem->Exit(2);return;}
  std::ofstream report(reportPath);TFile* files[3]={0,0,0};std::ostringstream body;body<<std::setprecision(17);
  try {
    Need(bool(report)&&pathB&&pathA&&pathC&&expectedEvents>0,"invalid explicit arguments");
    const char* paths[]={pathB,pathA,pathC};const char* labels[]={"B_continuous100","A_reset100","C_continuous20"};
    std::map<std::string,TKey*> se[3];std::string selection,input,definitions;
    std::map<std::string,std::string> baselineMix;
    body<<"\"expectedEvents\":"<<expectedEvents<<",\"inputs\":[";
    for(int f=0;f<3;++f) {
      files[f]=TFile::Open(paths[f],"READ");Need(files[f]&&!files[f]->IsZombie()&&!files[f]->TestBit(TFile::kRecovered)&&!rootError,"invalid pilot ROOT");
      Need(Object(files[f],"FemtoLambdaRunStatus")=="completed"&&Object(files[f],"FemtoLambdaReducedTreeMode")=="replay","pilot producer status/mode mismatch");
      Need(Parameter(files[f],"inputEvents")==expectedEvents&&Parameter(files[f],"reconstructedEvents")==expectedEvents,"pilot event count mismatch");
      const std::string settings=Object(files[f],"FemtoLambdaReducedTreeSettings"),mixing=Object(files[f],"FemtoLambdaReducedReplayMixing");
      std::map<std::string,std::string> s=Scalars(settings),m=Scalars(mixing);
      Need(s["mode"]=="replay"&&s["requireOriginalMixing"]=="false"&&s["resetPoolAtSegmentBoundary"]==(f==1?"true":"false"),"B/A/C reduced-tree settings mislabeled");
      Need(m["bufferSize"]==(f==2?"20":"100")&&m["maxMixEvents"]==(f==2?"20":"100")&&m["mixingMode"]=="bufferAll"&&m["mixBothDirections"]=="true","B/A/C buffer/mixing settings mislabeled");
      m.erase("bufferSize");m.erase("maxMixEvents");
      if(!f){baselineMix=m;selection=Object(files[f],"FemtoLambdaReducedReplaySelection");input=Object(files[f],"FemtoLambdaReducedTreeInput");definitions=Object(files[f],"FemtoLambdaReducedQAHistogramDefinitions");}
      Need(m==baselineMix&&Object(files[f],"FemtoLambdaReducedReplaySelection")==selection&&Object(files[f],"FemtoLambdaReducedTreeInput")==input&&Object(files[f],"FemtoLambdaReducedQAHistogramDefinitions")==definitions,"B/A/C physics scope or input differs beyond intended settings");
      Inventory(files[f],"",se[f]);Need(!se[f].empty(),"no SE histograms");
      if(f)body<<',';body<<"{\"label\":"<<Json(labels[f])<<",\"path\":"<<Json(paths[f])<<",\"mainconf\":"<<Json(Object(files[f],"FemtoLambdaMainconf"))<<",\"jobId\":"<<Json(Object(files[f],"FemtoLambdaReducedTreeJobId"))<<",\"reducedSettings\":"<<Json(settings)<<",\"mixingProvenance\":"<<Json(mixing)<<",\"parameters\":{";
      TIter keys(files[f]->GetListOfKeys());TKey* key=0;bool first=true;std::set<std::string> seen;
      while((key=dynamic_cast<TKey*>(keys()))) {
        if(!seen.insert(key->GetName()).second)continue;
        if(!Starts(key->GetClassName(),"TParameter<"))continue;
        TParameter<Long64_t>* value=dynamic_cast<TParameter<Long64_t>*>(files[f]->Get(key->GetName()));
        if(value){if(!first)body<<',';first=false;body<<Json(key->GetName())<<':'<<value->GetVal();}
      }
      body<<"}}";
    }
    body<<"],\"samePrefixVerifiedByFullScalarLedger\":true";
    std::vector<Long64_t> starts;SameLedger(files,expectedEvents,starts);body<<",\"segmentsObserved\":"<<starts.size()<<",\"segmentStartEvents\":[";
    for(size_t i=0;i<starts.size();++i){if(i)body<<',';body<<starts[i];}body<<']';
    Long64_t bins=0;int nonempty=0;
    Need(se[0].size()==se[1].size()&&se[0].size()==se[2].size(),"SE histogram union differs");
    for(std::map<std::string,TKey*>::const_iterator i=se[0].begin();i!=se[0].end();++i) {
      for(int f=1;f<3;++f){Need(se[f].count(i->first)!=0,"missing SE histogram");Equal(i->second,se[f][i->first],i->first,bins);}
      TH1* h=dynamic_cast<TH1*>(i->second->ReadObj());if(h&&h->GetEntries()>0)++nonempty;delete h;
    }
    Need(nonempty>0,"empty SE comparison is inconclusive");
    body<<",\"se\":{\"histograms\":"<<se[0].size()<<",\"nonemptyHistograms\":"<<nonempty<<",\"flowInclusiveBinComparisons\":"<<bins<<",\"allContentErrorEntriesAxesExact\":true}";
    const char* qaNames[]={"attempted","eligible","filled","skipped","capSkipped","pairCutSkipped","eligibleDirections","emptyDirections","selectedEmptyDirections","filledForward","filledReverse","eligibleForward","eligibleReverse","overlapSkipped","signalWindowSkipped"};
    std::vector<std::vector<double> > totals(3,std::vector<double>(15,0.));
    body<<",\"meTelemetry\":[";
    for(int f=0;f<3;++f) {
      TH2* h=dynamic_cast<TH2*>(files[f]->Get("hMixSamplerQA"));Need(h&&h->GetNbinsX()==15,"missing/wrong sampler QA");
      if(f)body<<',';body<<"{\"label\":"<<Json(labels[f])<<",\"channels\":[";
      for(int y=1;y<=h->GetNbinsY();++y) {
        if(y>1)body<<',';body<<"{\"channelIndex\":"<<h->GetYaxis()->GetBinCenter(y);
        for(int q=0;q<15;++q){double value=h->GetBinContent(h->GetXaxis()->FindBin(q),y);Need(TMath::Finite(value)&&value>=0,"invalid sampler counter");totals[f][q]+=value;body<<','<<Json(qaNames[q])<<':'<<value;}
        body<<'}';
      }
      Need(totals[f][0]==totals[f][1]&&totals[f][4]==0&&totals[f][2]==totals[f][9]+totals[f][10]&&totals[f][1]==totals[f][11]+totals[f][12],"sampler accounting/cap failure");
      body<<"],\"totals\":{";
      for(int q=0;q<15;++q){if(q)body<<',';body<<Json(qaNames[q])<<':'<<totals[f][q];}body<<"}}";
    }
    body<<"],\"BminusA\":{\"eligible\":"<<totals[0][1]-totals[1][1]<<",\"filled\":"<<totals[0][2]-totals[1][2]<<"},\"BminusC\":{\"eligible\":"<<totals[0][1]-totals[2][1]<<",\"filled\":"<<totals[0][2]-totals[2][2]<<'}';
    Need(totals[0][1]>totals[1][1]&&totals[0][1]>totals[2][1]&&totals[0][2]>totals[1][2]&&totals[0][2]>totals[2][2],"ME increase is absent/inconclusive");
    Need(totals[0][9]>0&&totals[0][10]>0,"both ME directions must be observed");
    TH1* signal=dynamic_cast<TH1*>(files[0]->Get("hKstarME_lambda_he4_signal"));Need(signal&&signal->GetEntries()>0,"nonempty He4 signal ME is required");
    body<<",\"Bhe4SignalMeEntries\":"<<signal->GetEntries()<<",\"noPayloadRescan\":true";
    Need(!rootError,"ROOT error during comparison");
    report<<"{\"passed\":true,"<<body.str()<<"}\n";report.close();Need(bool(report),"report write failed");
    for(int f=0;f<3;++f){files[f]->Close();delete files[f];}SetErrorHandler(DefaultErrorHandler);
    std::cout<<"PERIOD_MIX_COMPARISON_PASS "<<reportPath<<std::endl;
  } catch(const std::exception& error) {
    if(report){report<<"{\"passed\":false,\"error\":"<<Json(error.what())<<"}\n";report.close();}
    std::cerr<<"PERIOD_MIX_COMPARISON_FAILED "<<error.what()<<std::endl;SetErrorHandler(DefaultErrorHandler);gSystem->Exit(1);
  }
}
