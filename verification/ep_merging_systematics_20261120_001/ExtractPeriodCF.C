void ExtractPeriodCF(const char* planPath, const char* requirementsPath, const char* rawPath);
#ifndef __CINT__
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TNamed.h>
#include <TMath.h>
#include <TError.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PeriodCFExtract {
bool rootError=false;
void Handler(Int_t level,Bool_t abort,const char* where,const char* message) {
  if(level>=kError)rootError=true;
  DefaultErrorHandler(level,abort,where,message);
}
void Need(bool ok,const std::string& text) {if(!ok)throw std::runtime_error(text);}
std::string Json(const std::string& value) {
  std::ostringstream out;out<<'"';
  for(size_t i=0;i<value.size();++i) {
    const unsigned char c=value[i];
    if(c=='"'||c=='\\')out<<'\\'<<char(c);
    else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
    else out<<char(c);
  }
  return out.str()+'"';
}
std::string Trim(const std::string& text) {
  const size_t a=text.find_first_not_of(" \t\r");
  return a==std::string::npos?"":text.substr(a,text.find_last_not_of(" \t\r")-a+1);
}
std::vector<std::string> Fields(const std::string& text) {
  std::vector<std::string> result;size_t start=0;
  for(;;) {size_t end=text.find('\t',start);result.push_back(text.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}
  return result;
}
Long64_t Integer(const std::string& text) {
  std::istringstream s(text);Long64_t n=0;Need(bool(s>>n),"invalid integer: "+text);s>>std::ws;Need(s.eof(),"trailing integer text");return n;
}
double Number(const std::string& text) {
  std::istringstream s(text);double n=0;Need(bool(s>>n)&&TMath::Finite(n),"invalid number: "+text);s>>std::ws;Need(s.eof(),"trailing numeric text");return n;
}
double DefinedNumber(const std::string& text,const std::string& prefix,const std::string& suffix) {
  const size_t a=text.find(prefix);Need(a!=std::string::npos,"missing required cut definition field: "+prefix);
  const size_t first=a+prefix.size(),last=text.find(suffix,first);Need(last!=std::string::npos,"unterminated cut definition field");
  return Number(text.substr(first,last-first));
}
std::string Object(TFile* file,const char* name) {
  TObject* object=file->Get(name);Need(object!=0,std::string("missing metadata: ")+name);
  TObjString* text=dynamic_cast<TObjString*>(object);if(text)return text->GetString().Data();
  TNamed* named=dynamic_cast<TNamed*>(object);Need(named!=0,"metadata is not text");return named->GetTitle();
}
Long64_t Counter(TFile* file,const char* name) {
  TParameter<Long64_t>* value=dynamic_cast<TParameter<Long64_t>*>(file->Get(name));Need(value!=0,std::string("missing counter: ")+name);return value->GetVal();
}
std::map<std::string,std::string> Scalars(const std::string& text) {
  std::map<std::string,std::string> values;std::istringstream s(text);std::string line;
  while(std::getline(s,line)) {size_t p=line.find(':');if(p==std::string::npos)continue;std::string k=Trim(line.substr(0,p));if(k.empty()||k[0]=='-'||k[0]=='#')continue;Need(values.insert(std::make_pair(k,Trim(line.substr(p+1)))).second,"duplicate metadata scalar");}
  return values;
}
bool Axis(const TAxis* a,const TAxis* b) {
  if(a->GetNbins()!=b->GetNbins())return false;
  for(int i=1;i<=a->GetNbins()+1;++i)if(a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;
  return true;
}
void IntegerAxis(const TAxis* a,int n) {
  Need(a->GetNbins()==n,"integer axis size mismatch");
  for(int i=1;i<=n+1;++i)Need(a->GetBinLowEdge(i)==double(i)-1.5,"integer axis edge mismatch");
}
struct Row {std::string id,path,source;int ep,cut,baseline,period;Long64_t events;};
struct Data {
  int ep,cut,baseline,nx,sources;Long64_t events;
  std::vector<double> values,variances;
  Data():ep(0),cut(0),baseline(0),nx(0),sources(0),events(0){}
  size_t Index(int mixed,int plane,int cent,int x)const {return (((mixed*(ep+1)+(plane+1))*9+cent)*(nx+2)+x);}
};
void Array(std::ostream& out,const std::vector<double>& values,size_t first,size_t count) {
  out<<'[';for(size_t i=0;i<count;++i){if(i)out<<',';out<<values[first+i];}out<<']';
}
} // namespace PeriodCFExtract

// Read four signal histograms and small provenance objects, never TTree entries.
void ExtractPeriodCF(const char* planPath,const char* requirementsPath,const char* rawPath) {
  using namespace PeriodCFExtract;
  ErrorHandlerFunc_t previous=SetErrorHandler(Handler);rootError=false;TH1::AddDirectory(kFALSE);
  try {
    Need(planPath&&requirementsPath&&rawPath&&gSystem->AccessPathName(rawPath),"missing arguments or existing raw output");
    std::map<std::string,std::string> requirements;std::ifstream req(requirementsPath);Need(bool(req),"cannot open requirements");std::string line;
    while(std::getline(req,line)){std::vector<std::string> f=Fields(line);Need(f.size()==2&&requirements.insert(std::make_pair(f[0],f[1])).second,"invalid requirements row");}
    Need(req.eof(),"requirements read failed");
    const char* needed[]={"SECent","MECent","SEPlaneCent","MEPlaneCent","absoluteTolerance","relativeTolerance","baselineJobid","cut.phiMax","cut.etaMax","cut.radiusMeters"};
    for(size_t i=0;i<sizeof(needed)/sizeof(*needed);++i)Need(requirements.count(needed[i])!=0,"missing requirement");
    double atol=Number(requirements["absoluteTolerance"]),rtol=Number(requirements["relativeTolerance"]);Need(atol>=0&&rtol>=0,"negative tolerance");
    std::ifstream plan(planPath);Need(bool(plan),"cannot open plan");std::vector<Row> rows;std::set<std::string> paths,identities;
    while(std::getline(plan,line)) {
      const std::vector<std::string> f=Fields(line);Need(f.size()==8,"plan row must have eight columns");
      Row r;r.id=f[0];r.ep=int(Integer(f[1]));r.cut=int(Integer(f[2]));r.baseline=int(Integer(f[3]));r.period=int(Integer(f[4]));r.events=Integer(f[5]);r.source=f[6];r.path=f[7];
      Need((r.ep==1||r.ep==8||r.ep==16||r.ep==32)&&(r.cut==0||r.cut==1)&&(r.baseline==0||r.baseline==1)&&r.period>=0&&r.events>0,"invalid plan values");
      Need(paths.insert(r.path).second&&identities.insert(r.id+"/"+f[4]).second,"duplicate ROOT path or case/period");rows.push_back(r);
    }
    Need(plan.eof()&&!rows.empty(),"empty/failed plan");
    std::map<std::string,Data> data;std::vector<double> edges;
    std::string selection,epDefinition,reducedSettings;std::map<std::string,std::string> mixingCommon;
    std::map<int,std::string> periodEventCounters;std::string commonCutDefinition;std::ostringstream sourceJson;sourceJson<<std::setprecision(17)<<'[';
    for(size_t row=0;row<rows.size();++row) {
      const Row& r=rows[row];
      TFile* file=TFile::Open(r.path.c_str(),"READ");Need(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered),"invalid ROOT: "+r.path);
      Need(Object(file,"FemtoLambdaRunStatus")=="completed"&&Object(file,"FemtoLambdaReducedTreeMode")=="replay","incomplete/non-replay ROOT: "+r.path);
      const Long64_t read=Counter(file,"inputEvents"),done=Counter(file,"reconstructedEvents");Need(read==r.events&&done==r.events,"source event count mismatch: "+r.path);
      const std::string source=Object(file,"FemtoLambdaReducedSourceInput");std::istringstream sourceLines(source);std::string sourceLine;int sourceCount=0;
      while(std::getline(sourceLines,sourceLine)) {
        sourceLine=Trim(sourceLine);if(sourceLine.empty()||sourceLine.compare(0,16,"source argument:")==0)continue;
        const size_t slash=sourceLine.find_last_of('/');Need(sourceLine.substr(slash==std::string::npos?0:slash+1)==r.source,"source period identity/order mismatch: "+r.path);++sourceCount;
      }
      Need(sourceCount==1,"output is not from exactly one expected period");
      const std::string settingsText=Object(file,"FemtoLambdaReducedTreeSettings"),selectionText=Object(file,"FemtoLambdaReducedReplaySelection"),epText=Object(file,"FemtoEpdEventPlaneDefinition");
      std::map<std::string,std::string> settings=Scalars(settingsText),mixing=Scalars(Object(file,"FemtoLambdaReducedReplayMixing"));
      for(std::map<std::string,std::string>::const_iterator q=requirements.begin();q!=requirements.end();++q) {
        if(q->first.compare(0,4,"mix.")==0)Need(mixing[q->first.substr(4)]==q->second,"mixing setting mismatch: "+q->first);
        if(q->first.compare(0,8,"reduced.")==0)Need(settings[q->first.substr(8)]==q->second,"reduced setting mismatch: "+q->first);
      }
      Need(Integer(mixing["nEventPlaneBins"])==r.ep,"mixing EP bin count mismatch");mixing.erase("nEventPlaneBins");
      TParameter<Int_t>* cut=dynamic_cast<TParameter<Int_t>*>(file->Get("FemtoLambdaMergingPairCutEnabled"));
      TParameter<Int_t>* pairBins=dynamic_cast<TParameter<Int_t>*>(file->Get("FemtoLambdaPairEpdBinCount"));
      if(r.baseline) {
        Need(r.ep==8&&r.cut==0&&Object(file,"FemtoLambdaReducedTreeJobId").compare(0,requirements["baselineJobid"].size(),requirements["baselineJobid"])==0,"legacy exception is not the pinned baseline");
        if(cut)Need(cut->GetVal()==0,"baseline cut unexpectedly enabled");if(pairBins)Need(pairBins->GetVal()==r.ep,"baseline pair EP count mismatch");
      } else {
        Need(cut&&cut->GetVal()==r.cut&&pairBins&&pairBins->GetVal()==r.ep,"missing/mismatched explicit cut or pair EP metadata");
        std::string definition=Object(file,"FemtoLambdaMergingPairCutDefinition");
        Need(definition.compare(0,10,r.cut?"enabled=1;":"enabled=0;")==0,"cut definition enabled field disagrees");
        Need(DefinedNumber(definition,"abs(wrapped deltaPhi)<"," AND")==Number(requirements["cut.phiMax"])&&
             DefinedNumber(definition,"abs(deltaEta)<",";")==Number(requirements["cut.etaMax"])&&
             DefinedNumber(definition,"radiusMeters=",";")==Number(requirements["cut.radiusMeters"]),"cut threshold or radius mismatch");
        Need(definition.find("physical He4/proton;")!=std::string::npos&&definition.find("same SE/both ME directions")!=std::string::npos&&
             definition.find("unreachable retained and counted")!=std::string::npos&&definition.find("missing track/event cache is an error")!=std::string::npos&&
             definition.find("all Lambda mass views share cut; pion diagnostic only; He3 unchanged")!=std::string::npos,"pair-cut acceptance semantics differ");
        definition.erase(0,10);if(commonCutDefinition.empty())commonCutDefinition=definition;
        Need(definition==commonCutDefinition,"pair-cut definition differs beyond enabled flag");
      }
      if(row==0){selection=selectionText;epDefinition=epText;reducedSettings=settingsText;mixingCommon=mixing;}
      Need(selectionText==selection&&epText==epDefinition&&settingsText==reducedSettings&&mixing==mixingCommon,"common selection/EP calibration/reduced settings/mixing differ beyond N");
      const Long64_t epSeen=Counter(file,"epdEventsSeen"),epValid=Counter(file,"epdEventsValid"),epRejected=Counter(file,"epdEventsRejected");
      const Long64_t selected=Counter(file,"selectedLambdaCandidates");
      std::ostringstream counters;counters<<read<<','<<done<<','<<epSeen<<','<<epValid<<','<<epRejected<<','<<selected;
      if(!periodEventCounters.count(r.period))periodEventCounters[r.period]=counters.str();
      Need(periodEventCounters[r.period]==counters.str(),"EP-valid event population differs between cases");
      TH2* h2[2];TH3* h3[2];
      h2[0]=dynamic_cast<TH2*>(file->Get(requirements["SECent"].c_str()));h2[1]=dynamic_cast<TH2*>(file->Get(requirements["MECent"].c_str()));
      h3[0]=dynamic_cast<TH3*>(file->Get(requirements["SEPlaneCent"].c_str()));h3[1]=dynamic_cast<TH3*>(file->Get(requirements["MEPlaneCent"].c_str()));
      Need(h2[0]&&h2[1]&&h3[0]&&h3[1],"missing signal histogram");
      const int nx=h2[0]->GetNbinsX();Need(nx>0&&h2[0]->GetDimension()==2&&h2[1]->GetDimension()==2,"invalid kstar/cent histograms");
      for(int m=0;m<2;++m) {
        Need(Axis(h2[0]->GetXaxis(),h2[m]->GetXaxis())&&Axis(h2[0]->GetXaxis(),h3[m]->GetXaxis()),"kstar axes differ");
        IntegerAxis(h2[m]->GetYaxis(),9);IntegerAxis(h3[m]->GetZaxis(),9);IntegerAxis(h3[m]->GetYaxis(),r.ep);
      }
      if(edges.empty())for(int x=1;x<=nx+1;++x)edges.push_back(h2[0]->GetXaxis()->GetBinLowEdge(x));
      Need(edges.size()==size_t(nx+1),"native kstar bin count differs");
      for(int x=1;x<=nx+1;++x)Need(TMath::Finite(edges[x-1])&&edges[x-1]==h2[0]->GetXaxis()->GetBinLowEdge(x)&&(x==1||edges[x-1]>edges[x-2]),"native kstar edges differ/invalid");
      Data& d=data[r.id];
      if(d.sources==0){d.ep=r.ep;d.cut=r.cut;d.baseline=r.baseline;d.nx=nx;d.values.assign(2*(r.ep+1)*9*(nx+2),0.);d.variances=d.values;}
      Need(d.ep==r.ep&&d.cut==r.cut&&d.baseline==r.baseline&&d.nx==nx,"case plan inconsistent");
      for(int m=0;m<2;++m)for(int c=0;c<9;++c)for(int x=0;x<=nx+1;++x) {
        const double count=h2[m]->GetBinContent(x,c+1),error=h2[m]->GetBinError(x,c+1);
        Need(TMath::Finite(count)&&TMath::Finite(error)&&count>=0&&error>=0,"invalid raw count/error");
        const size_t inclusive=d.Index(m,-1,c,x);d.values[inclusive]+=count;d.variances[inclusive]+=error*error;
        double planeSum=0;
        for(int p=0;p<r.ep;++p) {
          const double v=h3[m]->GetBinContent(x,p+1,c+1),e=h3[m]->GetBinError(x,p+1,c+1);
          Need(TMath::Finite(v)&&TMath::Finite(e)&&v>=0&&e>=0,"invalid plane raw count/error");
          const size_t i=d.Index(m,p,c,x);d.values[i]+=v;d.variances[i]+=e*e;planeSum+=v;
        }
        Need(std::fabs(planeSum-count)<=atol+rtol*std::max(planeSum,count),"EP projection does not close to inclusive signal histogram");
      }
      // Centrality/EP flow would fall outside every declared analysis selection.
      for(int m=0;m<2;++m)for(int x=0;x<=nx+1;++x) {
        Need(h2[m]->GetBinContent(x,0)==0&&h2[m]->GetBinContent(x,10)==0,"nonzero centrality flow");
        for(int c=0;c<=10;++c)Need(h3[m]->GetBinContent(x,0,c)==0&&h3[m]->GetBinContent(x,r.ep+1,c)==0,"nonzero plane flow");
        for(int p=0;p<=r.ep+1;++p)Need(h3[m]->GetBinContent(x,p,0)==0&&h3[m]->GetBinContent(x,p,10)==0,"nonzero plane histogram centrality flow");
      }
      if(row)sourceJson<<',';
      sourceJson<<"{\"caseId\":"<<Json(r.id)<<",\"period\":"<<r.period<<",\"path\":"<<Json(r.path)<<",\"sourcePeriod\":"<<Json(r.source)<<",\"inputEvents\":"<<read<<",\"completedEvents\":"<<done<<",\"epdEventsSeen\":"<<epSeen<<",\"epdEventsValid\":"<<epValid<<",\"epdEventsRejected\":"<<epRejected<<",\"selectedLambdaCandidates\":"<<selected<<'}';
      ++d.sources;d.events+=done;file->Close();delete file;Need(!rootError,"ROOT error during extraction");
      std::cout<<"CF_HIST_SOURCE_PASS "<<r.id<<" period="<<r.period<<std::endl;
    }
    sourceJson<<']';
    std::ofstream output(rawPath);Need(bool(output),"cannot create raw JSON");output<<std::setprecision(17);
    output<<"{\"schemaVersion\":1,\"passed\":true,\"treeEntriesRead\":0,\"selection\":"<<Json(selection)<<",\"epDefinition\":"<<Json(epDefinition)<<",\"reducedSettings\":"<<Json(reducedSettings)<<",\"nativeEdges\":";Array(output,edges,0,edges.size());
    output<<",\"sources\":"<<sourceJson.str()<<",\"cases\":[";bool firstCase=true;
    for(std::map<std::string,Data>::const_iterator it=data.begin();it!=data.end();++it) {
      if(!firstCase)output<<',';firstCase=false;const Data& d=it->second;
      output<<"{\"caseId\":"<<Json(it->first)<<",\"epBins\":"<<d.ep<<",\"cutEnabled\":"<<(d.cut?"true":"false")<<",\"sources\":"<<d.sources<<",\"completedEvents\":"<<d.events<<",\"slices\":[";bool firstSlice=true;
      for(int p=-1;p<d.ep;++p)for(int c=0;c<9;++c) {
        if(!firstSlice)output<<',';firstSlice=false;output<<"{\"plane\":"<<p<<",\"cent9\":"<<c;
        for(int m=0;m<2;++m){output<<','<<Json(m?"ME":"SE")<<':';Array(output,d.values,d.Index(m,p,c,0),d.nx+2);output<<','<<Json(m?"MEVariance":"SEVariance")<<':';Array(output,d.variances,d.Index(m,p,c,0),d.nx+2);}
        output<<'}';
      }
      output<<"]}";
    }
    output<<"]}\n";output.close();Need(bool(output)&&!rootError,"raw output failed");SetErrorHandler(previous);std::cout<<"CF_HIST_EXTRACTION_PASS sources="<<rows.size()<<std::endl;
  } catch(const std::exception& e) {SetErrorHandler(previous);std::cerr<<"CF_HIST_EXTRACTION_FAILED "<<e.what()<<std::endl;gSystem->Exit(1);}
}
#endif
