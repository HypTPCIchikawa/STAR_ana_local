#include "include/FemtoLambdaReducedTree.h"
#include <TFileMerger.h>
#include <TH1.h>
#include <TAxis.h>
#include <TClass.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TKey.h>
#include <TParameter.h>
#include <TArrayD.h>
#include <TError.h>
#include <TSystem.h>
#include <TROOT.h>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>

namespace GroupedTreeValidation {
bool readError=false;
void ErrorHandler(Int_t level,Bool_t abort,const char* place,const char* message) {
  if(level>=kError) readError=true;
  DefaultErrorHandler(level,abort,place,message);
}
void Need(bool ok,const std::string& message) { if(!ok) throw std::runtime_error(message); }
std::string Json(const std::string& text) {
  std::ostringstream out;out<<'"';
  for(size_t i=0;i<text.size();++i) {
    const unsigned char c=text[i];
    if(c=='"'||c=='\\') out<<'\\'<<char(c);
    else if(c=='\n') out<<"\\n";
    else if(c=='\r') out<<"\\r";
    else if(c=='\t') out<<"\\t";
    else if(c<32) out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c)<<std::dec;
    else out<<char(c);
  }
  return out.str()+'"';
}
std::string Hash(const std::string& text) {
  ULong64_t value=14695981039346656037ULL;
  for(size_t i=0;i<text.size();++i) { value^=static_cast<unsigned char>(text[i]);value*=1099511628211ULL; }
  std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<value;return out.str();
}
Long64_t Size(const std::string& path) {
  struct stat s;Need(lstat(path.c_str(),&s)==0 && S_ISREG(s.st_mode),"not a regular nonsymlink file: "+path);return s.st_size;
}
std::vector<std::string> Paths(const char* list) {
  std::ifstream in(list);Need(bool(in),"cannot read input list");std::vector<std::string> paths;std::set<std::string> seen;std::string line;
  while(std::getline(in,line)) {
    Need(!line.empty() && line[0]=='/' && line.find_first_of("\r\t\n")==std::string::npos,"list requires canonical absolute paths");
    Need(seen.insert(line).second,"duplicate input path");Size(line);paths.push_back(line);
  }
  Need(in.eof()&&!paths.empty(),"empty or incomplete input list");return paths;
}
TFile* Open(const std::string& path) {
  TFile* file=TFile::Open(path.c_str(),"READ");
  Need(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered)&&!readError,"unreadable/unclosed ROOT: "+path);return file;
}
std::vector<std::string> Species() { std::vector<std::string> x;x.push_back("lambda");x.push_back("he3");x.push_back("he4");return x; }
void Branches(TObjArray* branches,Long64_t entries,std::ostringstream& schema) {
  for(int i=0;i<branches->GetEntriesFast();++i) {
    TBranch* b=dynamic_cast<TBranch*>(branches->At(i));Need(b!=0,"not a branch");
    Need(b->GetEntries()==entries,"branch entry mismatch: "+std::string(b->GetName()));
    schema<<b->GetName()<<'\t'<<b->GetClassName()<<'\t'<<b->GetTitle()<<'\n';
    TObjArray* leaves=b->GetListOfLeaves();
    for(int j=0;j<leaves->GetEntriesFast();++j) {
      TLeaf* l=dynamic_cast<TLeaf*>(leaves->At(j));Need(l!=0,"not a leaf");
      schema<<l->GetName()<<':'<<l->GetTypeName()<<':'<<l->GetLenStatic()<<':'<<(l->GetLeafCount()?l->GetLeafCount()->GetName():"")<<'\n';
    }
    Branches(b->GetListOfBranches(),entries,schema);
  }
}
std::string TreeSchema(TTree* t) { std::ostringstream out;Branches(t->GetListOfBranches(),t->GetEntries(),out);return out.str(); }
std::string HistSchema(const TH1* h) {
  std::ostringstream out;out<<std::setprecision(17)<<h->ClassName()<<':'<<h->GetDimension();
  const TAxis* axes[]={h->GetXaxis(),h->GetYaxis(),h->GetZaxis()};
  for(int a=0;a<h->GetDimension();++a) {
    const TAxis* axis=axes[a];out<<':'<<axis->GetNbins()<<':'<<axis->GetXmin()<<':'<<axis->GetXmax();
    const TArrayD* bins=axis->GetXbins();
    for(int b=0;b<bins->GetSize();++b) out<<':'<<bins->At(b);
    for(int b=1;b<=axis->GetNbins();++b) if(*axis->GetBinLabel(b)) out<<':'<<b<<'='<<axis->GetBinLabel(b);
  }
  return out.str();
}
Long64_t Cells(const TH1* h) {
  Long64_t n=h->GetNbinsX()+2;
  if(h->GetDimension()>1)n*=h->GetNbinsY()+2;
  if(h->GetDimension()>2)n*=h->GetNbinsZ()+2;
  Need(n>0&&n<=std::numeric_limits<int>::max(),"unsupported QA bin count");return n;
}
struct Segment {
  ULong64_t high,low;Long64_t entries;std::string producer,layout,selection,mixing;
  bool operator==(const Segment& s)const { return high==s.high&&low==s.low&&entries==s.entries&&producer==s.producer&&layout==s.layout&&selection==s.selection&&mixing==s.mixing; }
};
struct Header {
  Long64_t entries,inputEvents,reconstructedEvents;int compression;
  std::string eventSchema,metadataSchema;
  std::vector<GroupedTreeValidation::Segment> segments;
  std::map<std::string,TH1*> hist;
  std::map<std::string,std::pair<Long64_t,std::string> > trees;
};
Long64_t Counter(TFile* f,const char* name) {
  TParameter<Long64_t>* x=dynamic_cast<TParameter<Long64_t>*>(f->Get(name));
  Need(x!=0,std::string("missing Long64 counter: ")+name);return x->GetVal();
}
void Walk(TDirectory* directory,const std::string& prefix,Header& h) {
  TIter next(directory->GetListOfKeys());TKey* key=0;
  while((key=dynamic_cast<TKey*>(next()))) {
    TClass* type=gROOT->GetClass(key->GetClassName());Need(type!=0,"unknown persistent class");
    if(!type->InheritsFrom(TH1::Class())&&!type->InheritsFrom(TTree::Class())&&!type->InheritsFrom(TDirectory::Class()))continue;
    TObject* value=key->ReadObj();Need(value!=0&&!readError,"object read failed");
    const std::string name=prefix+key->GetName();
    if(value->InheritsFrom(TDirectory::Class()))Walk(static_cast<TDirectory*>(value),name+"/",h);
    else if(value->InheritsFrom(TH1::Class())) {
      TH1* hist=static_cast<TH1*>(value);hist->SetDirectory(0);
      Need(h.hist.insert(std::make_pair(name,hist)).second,"multiple QA key cycles: "+name);
    } else {
      TTree* tree=static_cast<TTree*>(value);
      Need(h.trees.insert(std::make_pair(name,std::make_pair(tree->GetEntries(),TreeSchema(tree)))).second,"multiple tree key cycles: "+name);
      // TDirectory owns tree objects read from this file.
    }
  }
}
void Endpoints(TTree* tree,const std::vector<GroupedTreeValidation::Segment>& segments) {
  ULong64_t high=0,low=0;Long64_t ordinal=-1;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("segmentHigh",1);tree->SetBranchStatus("segmentLow",1);tree->SetBranchStatus("ordinal",1);
  Need(tree->SetBranchAddress("segmentHigh",&high)>=0&&tree->SetBranchAddress("segmentLow",&low)>=0&&tree->SetBranchAddress("ordinal",&ordinal)>=0,"cannot bind endpoint identity");
  Long64_t offset=0;
  for(size_t i=0;i<segments.size();++i) {
    const Segment& s=segments[i];
    if(s.entries) {
      const Long64_t positions[]={0,s.entries/2,s.entries-1};
      for(int j=0;j<3;++j) {
        Need(tree->GetEntry(offset+positions[j])>0&&!readError,"segment endpoint read failed");
        Need(high==s.high&&low==s.low&&ordinal==positions[j],"segment endpoint identity/order mismatch");
      }
    }
    offset+=s.entries;
  }
  tree->ResetBranchAddresses();tree->SetBranchStatus("*",1);
  Need(offset==tree->GetEntries(),"segment endpoint extent mismatch");
}
void ReadHeader(TFile* file,Header& h,int expectedQa) {
  TNamed* status=dynamic_cast<TNamed*>(file->Get(femto_candidate_store::kProducerStatusName));
  Need(status&&std::string(status->GetTitle())=="completed","producer status is not completed");
  TTree* tree=dynamic_cast<TTree*>(file->Get(femto_lambda_tree::kTreeName));
  TTree* meta=dynamic_cast<TTree*>(file->Get(femto_lambda_tree::kMetadataName));
  Need(tree&&meta&&meta->GetEntries()>0,"missing event/producer metadata tree");
  h.entries=tree->GetEntries();h.compression=file->GetCompressionLevel();h.inputEvents=Counter(file,"inputEvents");h.reconstructedEvents=Counter(file,"reconstructedEvents");
  Need(h.entries>=0&&h.reconstructedEvents==h.entries&&h.inputEvents>=h.entries,"event/counter mismatch");
  Need(tree->GetFileNumber()==0,"automatic TTree file splitting detected");
  h.eventSchema=TreeSchema(tree);h.metadataSchema=TreeSchema(meta);
  femto_lambda_tree::Metadata row;std::string error;
  Need(femto_lambda_tree::BindMetadata(meta,row,false,error),error);
  std::set<femto_lambda_tree::SegmentKey> seen;Long64_t count=0;
  for(Long64_t i=0;i<meta->GetEntries();++i) {
    Need(meta->GetEntry(i)>0&&!readError,"cannot read producer metadata");
    Need(row.completed&&row.entries>=0&&(row.segmentHigh||row.segmentLow)&&row.producerId&&row.layout&&row.selection&&row.mixing&&!row.producerId->empty()&&!row.selection->empty()&&!row.mixing->empty(),"incomplete producer metadata");
    Need(*row.layout==femto_lambda_tree::Layout(Species()),"reduced-tree layout mismatch");
    Need(seen.insert(femto_lambda_tree::SegmentKey(row.segmentHigh,row.segmentLow)).second,"duplicate segment key in file");
    Need(row.entries<=std::numeric_limits<Long64_t>::max()-count,"segment count overflow");count+=row.entries;
    Segment s;s.high=row.segmentHigh;s.low=row.segmentLow;s.entries=row.entries;s.producer=*row.producerId;s.layout=*row.layout;s.selection=*row.selection;s.mixing=*row.mixing;
    if(!h.segments.empty())Need(s.selection==h.segments[0].selection&&s.mixing==h.segments[0].mixing,"inconsistent segment provenance");
    h.segments.push_back(s);
  }
  meta->ResetBranchAddresses();Need(count==h.entries,"metadata entries do not cover whole event tree");
  Endpoints(tree,h.segments);Walk(file,"",h);
  Need(int(h.hist.size())==expectedQa,"producer QA histogram count mismatch");
  for(std::map<std::string,TH1*>::const_iterator q=h.hist.begin();q!=h.hist.end();++q)
    Need(TMath::Finite(q->second->GetEntries())&&q->second->GetEntries()>=0,"invalid QA Entries");
}
void Release(Header& h) { for(std::map<std::string,TH1*>::iterator i=h.hist.begin();i!=h.hist.end();++i)delete i->second;h.hist.clear(); }
int ReaderSample(const std::string& path,const Segment& segment,int requested,Long64_t entries) {
  femto_lambda_tree::Reader reader;
  Need(reader.Init(path.c_str(),Species(),segment.selection,segment.mixing,true),reader.LastError());
  Need(reader.EntriesTotal()==entries,"Reader preflight tree count mismatch");
  const Long64_t n=std::min(entries,Long64_t(requested));
  femto_lambda_tree::Event event;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;femto_lambda_tree::Payload payload;
  for(Long64_t i=0;i<n;++i)Need(reader.Next(event,physical,qa,payload),reader.LastError());
  return int(n);
}
void WriteFile(std::ostream& out,const std::string& path,const Header& h,int samples) {
  out<<"{\"path\":"<<Json(path)<<",\"entries\":"<<h.entries<<",\"inputEvents\":"<<h.inputEvents<<",\"reconstructedEvents\":"<<h.reconstructedEvents<<",\"compression\":"<<h.compression<<",\"eventSchemaFnv1a64\":"<<Json(Hash(h.eventSchema))<<",\"metadataSchemaFnv1a64\":"<<Json(Hash(h.metadataSchema))<<",\"readerSampleEvents\":"<<samples<<",\"segments\":[";
  for(size_t j=0;j<h.segments.size();++j) {
    const Segment& s=h.segments[j];if(j)out<<',';
    out<<"{\"high\":"<<Json(Form("%016llx",s.high))<<",\"low\":"<<Json(Form("%016llx",s.low))<<",\"entries\":"<<s.entries<<",\"producerId\":"<<Json(s.producer)<<",\"selectionFnv1a64\":"<<Json(Hash(s.selection))<<",\"mixingFnv1a64\":"<<Json(Hash(s.mixing))<<"}";
  }
  out<<"],\"qa\":[";bool first=true;
  for(std::map<std::string,TH1*>::const_iterator i=h.hist.begin();i!=h.hist.end();++i) {
    if(!first)out<<',';first=false;out<<"{\"name\":"<<Json(i->first)<<",\"schemaFnv1a64\":"<<Json(Hash(HistSchema(i->second)))<<",\"entries\":"<<i->second->GetEntries()<<",\"sumw2Size\":"<<i->second->GetSumw2N()<<"}";
  }
  out<<"]}";
}
ULong64_t FloatDistance(double a,double b) {
  union Number { Float_t f;UInt_t bits; } x,y;x.f=Float_t(a);y.f=Float_t(b);
  const ULong64_t mid=0x80000000ULL;
  const ULong64_t ix=(x.bits&0x80000000U)?mid-(x.bits&0x7fffffffU):mid+x.bits;
  const ULong64_t iy=(y.bits&0x80000000U)?mid-(y.bits&0x7fffffffU):mid+y.bits;
  return ix>iy?ix-iy:iy-ix;
}
void ValidErrors(const TH1* h,const std::string& name) {
  const int cells=int(Cells(h));const int n=h->GetSumw2N();
  Need(n==0||n==cells,"malformed Sumw2: "+name);
  for(int b=0;b<cells;++b) {
    Need(TMath::Finite(h->GetBinContent(b))&&TMath::Finite(h->GetBinError(b))&&h->GetBinError(b)>=0,"nonfinite QA content/error: "+name);
    if(n)Need(TMath::Finite(h->GetSumw2()->At(b))&&h->GetSumw2()->At(b)>=0,"negative/nonfinite stored Sumw2: "+name);
  }
}
} // namespace GroupedTreeValidation

int tree_group_validation(const char* action,const char* inputList,const char* reportPath,
    const char* mergedPath,Long64_t maxTreeBytes,int expectedQa,int sampleEvents,
    int maxOpenedFiles,int floatUlps) {
  using namespace GroupedTreeValidation;
  TH1::AddDirectory(kFALSE);readError=false;SetErrorHandler(ErrorHandler);
  try {
    Need(action&&inputList&&reportPath&&mergedPath&&maxTreeBytes>0&&expectedQa>0&&sampleEvents>=1&&sampleEvents<=16&&maxOpenedFiles>=3&&floatUlps>=0,"invalid explicit validation arguments");
    Need(gSystem->AccessPathName(reportPath),"refusing existing report");
    const Long64_t originalMax=TTree::GetMaxTreeSize();TTree::SetMaxTreeSize(maxTreeBytes);
    Need(TTree::GetMaxTreeSize()==maxTreeBytes,"MaxTreeSize setting failed");
    TClass* vectors=gROOT->GetClass("vector<TVector3>");Need(vectors&&vectors->GetCollectionProxy(),"missing vector<TVector3> dictionary");
    std::ofstream report(reportPath);Need(bool(report),"cannot create report");report<<std::setprecision(17);
    if(std::string(action)=="probe") {
      TFileMerger merger(kFALSE,kTRUE);merger.SetFastMethod(kTRUE);merger.SetMaxOpenedFiles(maxOpenedFiles);
      report<<"{\"passed\":true,\"rootVersion\":"<<Json(gROOT->GetVersion())<<",\"originalMaxTreeSize\":"<<originalMax<<",\"configuredMaxTreeSize\":"<<TTree::GetMaxTreeSize()<<",\"maxOpenedFiles\":"<<merger.GetMaxOpenedFies()<<",\"vectorTVector3Dictionary\":true,\"fastMethodRequested\":true,\"localInputCopies\":false}\n";
      report.close();Need(bool(report),"probe report write failed");return 0;
    }
    const std::vector<std::string> paths=Paths(inputList);
    if(std::string(action)=="merge") {
      Need(*mergedPath&&gSystem->AccessPathName(mergedPath),"refusing existing merge output");
      Long64_t bytes=0;int compression=-1;
      for(size_t i=0;i<paths.size();++i) {
        const Long64_t size=Size(paths[i]);Need(size<=std::numeric_limits<Long64_t>::max()-bytes,"input byte overflow");bytes+=size;
        TFile* file=Open(paths[i]);if(compression<0)compression=file->GetCompressionLevel();
        Need(compression==file->GetCompressionLevel(),"fast merge requires identical input compression");file->Close();delete file;
      }
      Need(maxTreeBytes>bytes+bytes/10+1000000000LL,"MaxTreeSize lacks margin above whole group bytes");
      {
        TFileMerger merger(kFALSE,kTRUE);merger.SetFastMethod(kTRUE);merger.SetMaxOpenedFiles(maxOpenedFiles);
        Need(merger.GetMaxOpenedFies()>=int(paths.size())+2,"open-file limit would introduce hidden merge batches");
        Need(merger.OutputFile(mergedPath,"CREATE",compression),"cannot create fresh staged output");
        for(size_t i=0;i<paths.size();++i)Need(merger.AddFile(paths[i].c_str()),"cannot add merge input");
        Need(!merger.HasCompressionChange(),"fast merge would recompress inputs");
        Need(merger.Merge()&&!readError,"fast TFileMerger failed");
      }
      TFile* file=Open(mergedPath);TTree* events=dynamic_cast<TTree*>(file->Get(femto_lambda_tree::kTreeName));
      Need(events&&events->GetFileNumber()==0,"merge split/lost event tree");file->Close();delete file;
      const std::string out(mergedPath);const std::string split=out.substr(0,out.size()-5)+"_1.root";
      Need(gSystem->AccessPathName(split.c_str()),"unexpected automatic split sibling");
      report<<"{\"passed\":true,\"mode\":\"merge_only_requires_independent_check\",\"inputFiles\":"<<paths.size()<<",\"sourceBytes\":"<<bytes<<",\"output\":"<<Json(mergedPath)<<",\"outputBytes\":"<<Size(mergedPath)<<",\"maxTreeSize\":"<<maxTreeBytes<<",\"compression\":"<<compression<<",\"fast\":true,\"localCopies\":false,\"published\":false}\n";
      report.close();Need(bool(report),"merge report write failed");return 0;
    }
    Need(std::string(action)=="catalog"||std::string(action)=="check","unknown action");
    const bool checking=std::string(action)=="check";
    Header total;total.entries=total.inputEvents=total.reconstructedEvents=0;total.compression=-1;
    std::set<femto_lambda_tree::SegmentKey> segments;
    std::map<std::string,std::string> qaSchemas;
    report<<"{\"scope\":\"all headers/metadata/QA; three identity samples per segment and bounded official Reader samples; no full event replay\",\"action\":"<<Json(action)<<",\"files\":[";
    for(size_t k=0;k<paths.size();++k) {
      TFile* file=Open(paths[k]);Header one;ReadHeader(file,one,expectedQa);
      if(!k){total.eventSchema=one.eventSchema;total.metadataSchema=one.metadataSchema;total.compression=one.compression;}
      Need(one.eventSchema==total.eventSchema&&one.metadataSchema==total.metadataSchema,"source branch schemas differ");
      Need(one.compression==total.compression,"source compression differs");
      for(std::map<std::string,TH1*>::const_iterator q=one.hist.begin();q!=one.hist.end();++q) {
        if(!k)qaSchemas[q->first]=HistSchema(q->second);
        else Need(qaSchemas.count(q->first)&&qaSchemas[q->first]==HistSchema(q->second),"source QA schema differs");
      }
      for(size_t j=0;j<one.segments.size();++j) {
        const Segment& s=one.segments[j];Need(segments.insert(femto_lambda_tree::SegmentKey(s.high,s.low)).second,"duplicate source producer segment");
        if(!total.segments.empty())Need(s.layout==total.segments[0].layout&&s.selection==total.segments[0].selection&&s.mixing==total.segments[0].mixing,"source provenance differs");
        total.segments.push_back(s);
      }
      total.entries+=one.entries;total.inputEvents+=one.inputEvents;total.reconstructedEvents+=one.reconstructedEvents;
      for(std::map<std::string,std::pair<Long64_t,std::string> >::const_iterator t=one.trees.begin();t!=one.trees.end();++t) {
        if(!k)total.trees[t->first]=t->second;
        else {Need(total.trees.count(t->first)&&total.trees[t->first].second==t->second.second,"other tree schema differs");total.trees[t->first].first+=t->second.first;}
      }
      if(k)Need(total.trees.size()==one.trees.size(),"tree name union differs");
      if(checking)for(std::map<std::string,TH1*>::iterator h=one.hist.begin();h!=one.hist.end();++h) {
        ValidErrors(h->second,h->first);
        if(!k){TH1* clone=dynamic_cast<TH1*>(h->second->Clone());Need(clone!=0,"QA clone failed");clone->SetDirectory(0);total.hist[h->first]=clone;}
        else {Need(total.hist.count(h->first)&&HistSchema(total.hist[h->first])==HistSchema(h->second),"QA schema differs");Need(total.hist[h->first]->Add(h->second),"native QA Add failed");}
      }
      const int sampled=ReaderSample(paths[k],one.segments[0],sampleEvents,one.entries);
      if(k)report<<',';WriteFile(report,paths[k],one,sampled);report.flush();
      Release(one);file->Close();delete file;
      std::cout<<"SOURCE_HEADER_PASS "<<k+1<<" "<<paths[k]<<std::endl;
    }
    report<<"],\"sourceEvents\":"<<total.entries<<",\"sourceInputEvents\":"<<total.inputEvents<<",\"producerSegments\":"<<total.segments.size();
    if(checking) {
      Need(*mergedPath,"check requires merged path");TFile* file=Open(mergedPath);Header actual;ReadHeader(file,actual,expectedQa);
      Need(actual.entries==total.entries&&actual.inputEvents==total.inputEvents&&actual.reconstructedEvents==total.reconstructedEvents,"merged counters/entries differ from source sum");
      Need(actual.segments==total.segments,"metadata producer segment values/order changed");
      Need(actual.eventSchema==total.eventSchema&&actual.metadataSchema==total.metadataSchema&&actual.trees==total.trees,"merged tree schema/count changed");
      Need(actual.hist.size()==total.hist.size(),"merged QA union differs");
      ULong64_t worstUlp=0;Long64_t nonidentical=0;
      report<<",\"qaComparison\":[";bool first=true;
      for(std::map<std::string,TH1*>::iterator i=total.hist.begin();i!=total.hist.end();++i) {
        Need(actual.hist.count(i->first),"missing merged QA");TH1* a=i->second;TH1* b=actual.hist[i->first];
        Need(HistSchema(a)==HistSchema(b)&&a->GetEntries()==b->GetEntries(),"merged QA schema/Entries mismatch");ValidErrors(b,i->first);
        Need(a->GetSumw2N()==b->GetSumw2N(),"merged QA Sumw2 presence mismatch");
        const std::string type=a->ClassName();const bool floats=type=="TH1F"||type=="TH2F"||type=="TH3F";
        ULong64_t maxUlp=0;Long64_t changes=0;const int cells=int(Cells(a));
        for(int bin=0;bin<cells;++bin) {
          const double x=a->GetBinContent(bin),y=b->GetBinContent(bin);
          Need(TMath::Finite(x)&&TMath::Finite(y),"nonfinite QA sum");
          if(x!=y) {
            ++changes;
            if(floats){const ULong64_t distance=FloatDistance(x,y);maxUlp=std::max(maxUlp,distance);Need(distance<=ULong64_t(floatUlps),"QA float ULP bound exceeded: "+i->first);}
            else Need(std::fabs(x-y)<=(2.*paths.size()+8.)*DBL_EPSILON*(1.+std::fabs(x)+std::fabs(y)),"QA double sum bound exceeded: "+i->first);
          }
          if(a->GetSumw2N()) {
            const double x2=a->GetSumw2()->At(bin),y2=b->GetSumw2()->At(bin);
            Need(TMath::Finite(x2)&&TMath::Finite(y2)&&x2>=0&&y2>=0,"invalid QA sum variance");
            Need(std::fabs(x2-y2)<=(2.*paths.size()+8.)*DBL_EPSILON*(1.+std::fabs(x2)+std::fabs(y2)),"QA variance sum bound exceeded: "+i->first);
          }
        }
        worstUlp=std::max(worstUlp,maxUlp);nonidentical+=changes;
        if(!first)report<<',';first=false;report<<"{\"name\":"<<Json(i->first)<<",\"nonidenticalBins\":"<<changes<<",\"maxFloatUlps\":"<<maxUlp<<"}";
      }
      const int samples=ReaderSample(mergedPath,actual.segments[0],sampleEvents,actual.entries);
      report<<"],\"merged\":";WriteFile(report,mergedPath,actual,samples);
      report<<",\"qaNonidenticalBins\":"<<nonidentical<<",\"maximumFloatUlps\":"<<worstUlp<<",\"configuredFloatUlps\":"<<floatUlps;
      Release(actual);file->Close();delete file;
    }
    Release(total);Need(!readError,"ROOT emitted a read/write error");
    report<<",\"passed\":true,\"filesChecked\":"<<paths.size()<<",\"fullEventReplayPerformed\":false}\n";report.close();Need(bool(report),"report write failed");
    SetErrorHandler(DefaultErrorHandler);return 0;
  } catch(const std::exception& error) { std::cerr<<"GROUPED_TREE_VALIDATION_FAILED "<<error.what()<<std::endl;SetErrorHandler(DefaultErrorHandler);return 1; }
}
