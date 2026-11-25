// Read-only sample closure and storage audit. Evidence paths must be new.
#include "../../include/FemtoLambdaReducedTree.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"
#include "TAxis.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TArrayD.h"
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace LambdaTreeVerification {
void Runtime() {
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0)throw std::runtime_error("cannot load STAR shared-library loader");
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error)throw std::runtime_error("STAR shared-library loader failed");
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker","StEpdUtil"};
  for(unsigned i=0;i<sizeof(star)/sizeof(star[0]);++i)if(gSystem->Load(star[i])<0)throw std::runtime_error(std::string("cannot load ")+star[i]);
  const char* project[]={"libStarAnaConfig.so","libStRefMultCorr.so","libKFParticle.so","libStKfParticleCommon.so","libStCommon.so","libStFemtoMaker.so"};
  const TString cwd=gSystem->WorkingDirectory();
  for(unsigned i=0;i<sizeof(project)/sizeof(project[0]);++i)if(gSystem->Load(cwd+"/lib/"+project[i])<0)throw std::runtime_error(std::string("cannot load ")+project[i]);
}
std::string Json(const std::string& text) {
  std::ostringstream out;out<<'"';
  for(size_t i=0;i<text.size();++i) {
    const unsigned char c=static_cast<unsigned char>(text[i]);
    if(c=='"'||c=='\\')out<<'\\'<<char(c);
    else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
    else out<<char(c);
  }
  out<<'"';return out.str();
}
std::string Base(const std::string& s) {const size_t p=s.find_last_of('/');return p==std::string::npos?s:s.substr(p+1);}
bool Starts(const std::string& s,const char* p){return s.compare(0,std::string(p).size(),p)==0;}
bool Pair(const std::string& path) {
  const std::string n=Base(path);return Starts(n,"hKstar")||Starts(n,"hQlab")||Starts(n,"hDphi")||Starts(n,"hNuclearPidPair");
}
bool PrePidPair(const std::string& path){const std::string n=Base(path);return n.find("NSigmaPrePid")!=std::string::npos||n.find("POverQPrePid")!=std::string::npos;}
bool HasPrePidControl(TFile* file) {
  const char* marker="prePID control bank omitted; nonselected prePID pair QA is outside the stored envelope";
  TObjString* selection=dynamic_cast<TObjString*>(file->Get("FemtoLambdaReducedReplaySelection"));
  const bool omitted=selection&&std::string(selection->GetString().Data()).find(marker)!=std::string::npos;
  TParameter<Int_t>* explicitValue=dynamic_cast<TParameter<Int_t>*>(file->Get("FemtoLambdaReducedHasPrePidControl"));
  if(!explicitValue)return !omitted;
  if(explicitValue->GetVal()!=0&&explicitValue->GetVal()!=1)throw std::runtime_error("invalid explicit prePID capability");
  const bool enabled=explicitValue->GetVal()!=0;
  if(enabled==omitted)throw std::runtime_error("prePID capability disagrees with selection fingerprint");
  return enabled;
}
bool Mixed(const std::string& p){return Starts(p,"mix/")||Starts(Base(p),"hKstarME")||Base(p).find("_Mixed_")!=std::string::npos;}
std::string Category(const std::string& path) {
  const std::string n=Base(path);
  if(n.find("PairY")!=std::string::npos)return "pairY";
  if(n.find("YLab")!=std::string::npos)return "he4YLab";
  if(n.find("Mt")!=std::string::npos)return "mT";
  if(n.find("NSigma")!=std::string::npos||n.find("POverQ")!=std::string::npos||n.find("PrePid")!=std::string::npos)return "PID_pair";
  if(Starts(n,"hDphi"))return "merging";if(Starts(n,"hQlab"))return "Qlab";
  if(n.find("Mass")!=std::string::npos)return "kstar_mass";return "kstar";
}
struct Inventory {
  std::map<std::string,TKey*> hist;std::vector<std::string> unsupportedPair;
  Long64_t keys,zip,total,keyBytes,histBytes,logical,allocated;
  Inventory():keys(0),zip(0),total(0),keyBytes(0),histBytes(0),logical(0),allocated(-1){}
};
void Inspect(TDirectory* directory,const std::string& prefix,const char* label,Inventory& data,std::ofstream& sizes) {
  TIter next(directory->GetListOfKeys());std::set<std::string> names;TKey* listed=0;
  while((listed=dynamic_cast<TKey*>(next()))) {
    const std::string name=listed->GetName();if(!names.insert(name).second)continue;
    TKey* key=directory->GetKey(name.c_str());if(!key)throw std::runtime_error("missing latest ROOT key "+name);
    const std::string path=prefix+name;TClass* type=TClass::GetClass(key->GetClassName());
    ++data.keys;data.keyBytes+=key->GetNbytes();Long64_t zip=0,total=0;
    if(type&&type->InheritsFrom(TDirectory::Class())) {
      TDirectory* child=directory->GetDirectory(name.c_str());if(!child)throw std::runtime_error("missing directory "+path);
      sizes<<label<<'\t'<<path<<'\t'<<key->GetClassName()<<'\t'<<key->GetNbytes()<<'\t'<<key->GetObjlen()<<"\t0\t0\n";
      Inspect(child,path+"/",label,data,sizes);continue;
    }
    if(type&&type->InheritsFrom(TH1::Class())){data.hist[path]=key;data.histBytes+=key->GetNbytes();}
    else if(Pair(path))data.unsupportedPair.push_back(path+" ["+key->GetClassName()+"]");
    if(type&&type->InheritsFrom(TTree::Class())) {
      TTree* tree=dynamic_cast<TTree*>(key->ReadObj());if(!tree)throw std::runtime_error("cannot read tree "+path);
      zip=tree->GetZipBytes();total=tree->GetTotBytes();data.zip+=zip;data.total+=total;
      tree->SetDirectory(0);delete tree;
    }
    sizes<<label<<'\t'<<path<<'\t'<<key->GetClassName()<<'\t'<<key->GetNbytes()<<'\t'<<key->GetObjlen()<<'\t'<<zip<<'\t'<<total<<'\n';
  }
}
struct Counter {
  Long64_t compared,failed,nonempty,content,error,entries,axis,stats;
  Double_t maxContent,maxError;
  Counter():compared(0),failed(0),nonempty(0),content(0),error(0),entries(0),axis(0),stats(0),maxContent(0),maxError(0){}
  void Add(const Counter& c){compared+=c.compared;failed+=c.failed;nonempty+=c.nonempty;content+=c.content;error+=c.error;entries+=c.entries;axis+=c.axis;stats+=c.stats;maxContent=TMath::Max(maxContent,c.maxContent);maxError=TMath::Max(maxError,c.maxError);}
};
struct Checks {
  std::vector<std::string> errors;std::ofstream* detail;Long64_t rows;
  Checks():detail(0),rows(0){}
  void Fail(const std::string& text){errors.push_back(text);std::cerr<<"FAIL "<<text<<std::endl;}
  void Row(const std::string& group,const std::string& key,const char* field,Long64_t bin,Double_t a,Double_t b){if(detail&&rows<5000){*detail<<group<<'\t'<<key<<'\t'<<field<<'\t'<<bin<<'\t'<<std::setprecision(17)<<a<<'\t'<<b<<'\n';++rows;}}
};
bool Different(Double_t a,Double_t b,Double_t tolerance){return !TMath::Finite(a)||!TMath::Finite(b)||TMath::Abs(a-b)>tolerance;}
bool SameAxis(const TAxis* a,const TAxis* b){if(a->GetNbins()!=b->GetNbins())return false;for(int i=1;i<=a->GetNbins()+1;++i)if(a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;return true;}
Int_t Cells(const TH1* histogram){Int_t n=histogram->GetNbinsX()+2;if(histogram->GetDimension()>1)n*=histogram->GetNbinsY()+2;if(histogram->GetDimension()>2)n*=histogram->GetNbinsZ()+2;return n;}
Counter Compare(TKey* left,TKey* right,const std::string& name,const char* group,Checks& check,Double_t tolerance) {
  Counter c;++c.compared;if(!left||!right){++c.failed;check.Fail(std::string(group)+": missing "+name);return c;}
  TH1* a=dynamic_cast<TH1*>(left->ReadObj());TH1* b=dynamic_cast<TH1*>(right->ReadObj());
  if(!a||!b){delete a;delete b;throw std::runtime_error("cannot read hist "+name);}
  bool failed=false;if(a->GetEntries()>0)++c.nonempty;
  if(a->GetDimension()!=b->GetDimension()||Cells(a)!=Cells(b)||!SameAxis(a->GetXaxis(),b->GetXaxis())||(a->GetDimension()>1&&!SameAxis(a->GetYaxis(),b->GetYaxis()))||(a->GetDimension()>2&&!SameAxis(a->GetZaxis(),b->GetZaxis()))){++c.axis;failed=true;check.Row(group,name,"axis",-1,Cells(a),Cells(b));}
  else for(Int_t bin=0;bin<Cells(a);++bin) {
    const Double_t ac=a->GetBinContent(bin),bc=b->GetBinContent(bin),ae=a->GetBinError(bin),be=b->GetBinError(bin);
    if(Different(ac,bc,tolerance)){++c.content;failed=true;check.Row(group,name,"content",bin,ac,bc);}
    if(Different(ae,be,tolerance)){++c.error;failed=true;check.Row(group,name,"error",bin,ae,be);}
    if(TMath::Finite(ac)&&TMath::Finite(bc))c.maxContent=TMath::Max(c.maxContent,TMath::Abs(ac-bc));
    if(TMath::Finite(ae)&&TMath::Finite(be))c.maxError=TMath::Max(c.maxError,TMath::Abs(ae-be));
  }
  if(Different(a->GetEntries(),b->GetEntries(),tolerance)){++c.entries;failed=true;check.Row(group,name,"entries",-1,a->GetEntries(),b->GetEntries());}
  Double_t av[13]={0},bv[13]={0};a->GetStats(av);b->GetStats(bv);
  for(int k=0;k<13;++k)if(Different(av[k],bv[k],tolerance)){++c.stats;check.Row(group,name,"stats_advisory",k,av[k],bv[k]);}
  if(failed){++c.failed;check.Fail(std::string(group)+": differs "+name);}delete a;delete b;return c;
}
TKey* Key(const Inventory& data,const std::string& path){std::map<std::string,TKey*>::const_iterator it=data.hist.find(path);return it==data.hist.end()?0:it->second;}
Double_t Entries(const Inventory& data,const char* path){TKey* key=Key(data,path);if(!key)return -1;TH1* h=dynamic_cast<TH1*>(key->ReadObj());if(!h)return -1;const Double_t n=h->GetEntries();delete h;return n;}
Long64_t Parameter(TFile* file,const char* name){TParameter<Long64_t>* p=dynamic_cast<TParameter<Long64_t>*>(file->Get(name));return p?p->GetVal():-1;}
int MixingBins(TFile* file){TObjString* text=dynamic_cast<TObjString*>(file->Get("FemtoLambdaReducedReplayMixing"));if(!text)return -1;const std::string s=text->GetString().Data();const size_t p=s.find("nEventPlaneBins");if(p==std::string::npos)return -1;const size_t colon=s.find(':',p);if(colon==std::string::npos)return -1;std::istringstream in(s.substr(colon+1));int bins=-1;in>>bins;return in?bins:-1;}
Double_t Variance(TH1* histogram,Int_t bin){return histogram->GetSumw2N()?histogram->GetSumw2()->At(bin):TMath::Abs(histogram->GetBinContent(bin));}
Counter PlaneSums(const Inventory& data,const char* label,Checks& check,Double_t tolerance,bool required,bool& available) {
  Counter total;available=false;
  const char* species[]={"he4","he3"};const char* tags[]={"","_signal","_leftSB","_rightSB"};const char* event[]={"SE","ME"};
  for(int s=0;s<2;++s)for(int t=0;t<4;++t)for(int e=0;e<2;++e)if(Key(data,std::string("hKstar")+event[e]+"VsEpdBinVsCent_lambda_"+species[s]+tags[t]))available=true;
  if(!available&&!required)return total;
  for(int s=0;s<2;++s)for(int t=0;t<4;++t)for(int e=0;e<2;++e) {
    Counter result;++result.compared;const std::string suffix=std::string("lambda_")+species[s]+tags[t];
    const std::string name=std::string("hKstar")+event[e]+"VsEpdBinVsCent_"+suffix;
    const std::string base=std::string("hKstar")+event[e]+"VsCent_"+suffix;
    TKey* sourceKey=Key(data,name);TKey* targetKey=Key(data,base);
    if(!sourceKey||!targetKey){++result.failed;check.Fail(std::string(label)+": missing plane/base hist "+name);total.Add(result);continue;}
    TH1* source=dynamic_cast<TH1*>(sourceKey->ReadObj());TH1* target=dynamic_cast<TH1*>(targetKey->ReadObj());
    if(!source||!target){delete source;delete target;throw std::runtime_error("cannot read plane/base hist "+name);}
    bool failed=false;if(source->GetEntries()>0)++result.nonempty;
    if(source->GetDimension()!=3||target->GetDimension()!=2||!SameAxis(source->GetXaxis(),target->GetXaxis())||!SameAxis(source->GetZaxis(),target->GetYaxis())||source->GetYaxis()->GetNbins()!=8||source->GetYaxis()->GetXmin()!=-0.5||source->GetYaxis()->GetXmax()!=7.5){++result.axis;failed=true;}
    else for(Int_t x=0;x<=source->GetNbinsX()+1;++x)for(Int_t z=0;z<=source->GetNbinsZ()+1;++z) {
      Double_t content=0,variance=0;for(Int_t y=0;y<=source->GetNbinsY()+1;++y){const Int_t bin=source->GetBin(x,y,z);content+=source->GetBinContent(bin);variance+=Variance(source,bin);}
      const Int_t targetBin=target->GetBin(x,z);const Double_t reference=target->GetBinContent(targetBin),error=TMath::Sqrt(variance),referenceError=target->GetBinError(targetBin);
      if(Different(reference,content,tolerance)){++result.content;failed=true;check.Row(label,name,"EP_sum_content",targetBin,reference,content);}
      if(Different(referenceError,error,tolerance)){++result.error;failed=true;check.Row(label,name,"EP_sum_error",targetBin,referenceError,error);}
      result.maxContent=TMath::Max(result.maxContent,TMath::Abs(reference-content));result.maxError=TMath::Max(result.maxError,TMath::Abs(referenceError-error));
    }
    if(Different(source->GetEntries(),target->GetEntries(),tolerance)){++result.entries;failed=true;check.Row(label,name,"EP_sum_entries",-1,target->GetEntries(),source->GetEntries());}
    if(failed){++result.failed;check.Fail(std::string(label)+": EP sum differs "+name);}total.Add(result);delete source;delete target;
  }
  return total;
}
struct TreeCounts {
  Long64_t entries,declared,tracks,lambda,he3,he4,controls,invalidEP,emptyTarget,dtOnly,pool,noTof,slots,extra,physical,admission,source;
  TreeCounts():entries(0),declared(0),tracks(0),lambda(0),he3(0),he4(0),controls(0),invalidEP(0),emptyTarget(0),dtOnly(0),pool(0),noTof(0),slots(0),extra(0),physical(0),admission(0),source(0){}
};
TreeCounts InspectTree(TFile* producer,const char* path,Checks& check) {
  TreeCounts c;TTree* tree=dynamic_cast<TTree*>(producer->Get(femto_lambda_tree::kTreeName));TTree* metadata=dynamic_cast<TTree*>(producer->Get(femto_lambda_tree::kMetadataName));
  if(!tree||!metadata)throw std::runtime_error("producer tree or metadata missing");
  const char* required[]={"epdPsi1","epdValid","baselinePoolEligible","hasAnyLegacyCandidate","sourceEntry","sourceFileIndex","segmentHigh","segmentLow","ordinal","trackPrimaryMomentum","trackGlobalMomentum","trackOrigin","trackFloat","trackDouble","trackFlags","lambdaP","he3P","he4P","lambdaFloat","lambdaDouble","lambdaProtonTrackSlot","lambdaPionTrackSlot","he3TrackSlot","he4TrackSlot","candidateFloat","candidateCollection"};
  for(unsigned i=0;i<sizeof(required)/sizeof(required[0]);++i)if(!tree->GetBranch(required[i]))check.Fail(std::string("missing branch ")+required[i]);
  femto_lambda_tree::Metadata m;std::string error,selection,mixing;
  if(!femto_lambda_tree::BindMetadata(metadata,m,false,error))throw std::runtime_error(error);
  for(Long64_t i=0;i<metadata->GetEntries();++i){if(metadata->GetEntry(i)<=0||!m.completed)throw std::runtime_error("incomplete producer metadata");c.declared+=m.entries;if(!i){selection=*m.selection;mixing=*m.mixing;}else if(selection!=*m.selection||mixing!=*m.mixing)throw std::runtime_error("inconsistent segment definitions");}
  metadata->ResetBranchAddresses();std::vector<std::string> species;species.push_back("lambda");species.push_back("he3");species.push_back("he4");
  femto_lambda_tree::Reader reader;if(!reader.Init(path,species,selection,mixing,true))throw std::runtime_error(reader.LastError());
  femto_lambda_tree::Event event;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;femto_lambda_tree::Payload payload;
  while(reader.Next(event,physical,qa,payload)) {
    ++c.entries;c.lambda+=physical["lambda"].size();c.he3+=physical["he3"].size();c.he4+=physical["he4"].size();c.controls+=qa["he3"].size()+qa["he4"].size();c.tracks+=payload.tracks.size();
    if(!event.epdValid)++c.invalidEP;if(payload.baselinePoolEligible)++c.pool;
    const bool empty=physical["lambda"].empty()&&physical["he3"].empty()&&physical["he4"].empty();if(empty)++c.emptyTarget;
    if(empty&&payload.baselinePoolEligible&&(payload.originalDeuteronCount||payload.originalTritonCount))++c.dtOnly;
    const bool any=!empty||payload.originalDeuteronCount>0||payload.originalTritonCount>0;
    if(any!=bool(payload.hasAnyLegacyCandidate)||bool(payload.baselinePoolEligible)!=(any&&bool(event.epdValid)))++c.admission;
    if(payload.sourceFileIndex<0||payload.sourceEntry<0)++c.source;
    std::map<int,unsigned> expected;
    for(size_t j=0;j<payload.lambdas.size();++j){expected[payload.lambdas[j].protonIndex]|=femto_lambda_tree::kProtonDaughter;expected[payload.lambdas[j].pionIndex]|=femto_lambda_tree::kPionDaughter;}
    for(int s=0;s<2;++s){const char* name=s?"he4":"he3";const unsigned role=s?femto_lambda_tree::kHe4Track:femto_lambda_tree::kHe3Track;for(size_t j=0;j<physical[name].size();++j)expected[physical[name][j].trk.trackIndex]|=role;}
    if(expected.size()!=payload.tracks.size())++c.extra;
    std::map<int,const femto_lambda_tree::Track*> actual;
    for(size_t j=0;j<payload.tracks.size();++j){const femto_lambda_tree::Track& t=payload.tracks[j];if(!actual.insert(std::make_pair(t.index,&t)).second||!expected.count(t.index)||expected[t.index]!=t.roleFlags)++c.extra;if(!t.hasTof)++c.noTof;}
    for(std::map<int,unsigned>::const_iterator it=expected.begin();it!=expected.end();++it)if(!actual.count(it->first))++c.slots;
    for(int s=0;s<2;++s){const char* name=s?"he4":"he3";for(size_t j=0;j<physical[name].size();++j){const FemtoCandidate& candidate=physical[name][j];if(!actual.count(candidate.trk.trackIndex)){++c.slots;continue;}const TVector3& p=actual[candidate.trk.trackIndex]->pMom;if(candidate.px!=Float_t(2*p.X())||candidate.py!=Float_t(2*p.Y())||candidate.pz!=Float_t(2*p.Z()))++c.physical;}}
  }
  if(!reader.Finished())throw std::runtime_error(reader.LastError());
  if(c.entries!=tree->GetEntries()||c.entries!=c.declared||c.entries!=Parameter(producer,"reconstructedEvents"))check.Fail("producer tree/metadata/reconstructed counts differ");
  if(!c.entries)check.Fail("producer tree is empty");
  if(c.slots||c.extra||c.physical||c.admission||c.source)check.Fail("target-track/admission/provenance validation failed");return c;
}
void PrintCounter(std::ostream& out,const Counter& c){out<<"{\"compared\":"<<c.compared<<",\"failed\":"<<c.failed<<",\"nonemptyReferenceHistograms\":"<<c.nonempty<<",\"contentBinMismatches\":"<<c.content<<",\"errorBinMismatches\":"<<c.error<<",\"entriesMismatches\":"<<c.entries<<",\"axisMismatches\":"<<c.axis<<",\"statsComponentsDifferingAdvisory\":"<<c.stats<<",\"maxContentDifference\":"<<c.maxContent<<",\"maxErrorDifference\":"<<c.maxError<<'}';}
} // namespace LambdaTreeVerification

void verify_outputs_v3(const char* baselinePath,const char* producerPath,const char* replayPath,const char* ep1Path,const char* summaryPath,Double_t absoluteTolerance=0.,Bool_t requirePlaneHistograms=kFALSE) {
  using namespace LambdaTreeVerification;
  if(!baselinePath||!*baselinePath||!producerPath||!*producerPath||!replayPath||!*replayPath||!ep1Path||!*ep1Path||!summaryPath||!*summaryPath||absoluteTolerance<0){std::cerr<<"ERROR: baseline/producer/replay/EP1/new-summary paths required"<<std::endl;gSystem->Exit(2);return;}
  const std::string output(summaryPath),keyPath=output+".keys.tsv",diffPath=output+".differences.tsv";
  if(!gSystem->AccessPathName(output.c_str())||!gSystem->AccessPathName(keyPath.c_str())||!gSystem->AccessPathName(diffPath.c_str())){std::cerr<<"ERROR: refusing existing verification evidence"<<std::endl;gSystem->Exit(2);return;}
  const TString dir=gSystem->DirName(summaryPath);if(gSystem->mkdir(dir,kTRUE)&&gSystem->AccessPathName(dir)){gSystem->Exit(2);return;}
  std::ofstream summary(output.c_str()),keys(keyPath.c_str()),differences(diffPath.c_str());if(!summary||!keys||!differences){std::cerr<<"ERROR: cannot create evidence"<<std::endl;gSystem->Exit(2);return;}
  keys<<"file\tkey\tclass\tkey_bytes\tserialized_object_bytes\ttree_zip_bytes\ttree_total_bytes\n";differences<<"group\tkey\tfield\tglobal_bin_including_flow\treference\tactual\n";
  Checks check;check.detail=&differences;Counter qa,pairs,same,planes;std::map<std::string,Counter> groups;Inventory data[4];TFile* files[4]={0,0,0,0};bool replayPlanes=false,ep1Planes=false;
  bool hasPrePidControl=true;std::vector<std::string> omittedPrePidHistograms;
  const char* paths[]={baselinePath,producerPath,replayPath,ep1Path};const char* labels[]={"baseline","producer","replay","ep1"};TreeCounts tree;Double_t counts[6]={-1,-1,-1,-1,-1,-1};bool meaningful=false;int replayBins=-1,ep1Bins=-1;Long64_t changedMixed=0;
  const Bool_t addDirectory=TH1::AddDirectoryStatus();TH1::AddDirectory(kFALSE);
  try {
    Runtime();
    for(int i=0;i<4;++i){files[i]=TFile::Open(paths[i],"READ");if(!files[i]||files[i]->IsZombie())throw std::runtime_error(std::string("cannot open ")+labels[i]);TNamed* status=dynamic_cast<TNamed*>(files[i]->Get("FemtoLambdaRunStatus"));if(!status||std::string(status->GetTitle())!="completed")check.Fail(std::string(labels[i])+" run is not completed");data[i].logical=files[i]->GetSize();struct stat info;if(stat(paths[i],&info)==0)data[i].allocated=Long64_t(info.st_blocks)*512;Inspect(files[i],"",labels[i],data[i],keys);for(size_t j=0;j<data[i].unsupportedPair.size();++j)check.Fail(std::string(labels[i])+": unsupported pair object "+data[i].unsupportedPair[j]);}
    hasPrePidControl=HasPrePidControl(files[1]);
    if(HasPrePidControl(files[2])!=hasPrePidControl||HasPrePidControl(files[3])!=hasPrePidControl)check.Fail("producer/replay prePID capabilities differ");
    if(!hasPrePidControl)for(int file=1;file<4;++file)for(std::map<std::string,TKey*>::const_iterator key=data[file].hist.begin();key!=data[file].hist.end();++key)if(PrePidPair(key->first))check.Fail(std::string(labels[file])+": compact output unexpectedly books nonselected prePID pair QA "+key->first);
    for(std::map<std::string,TKey*>::const_iterator it=data[1].hist.begin();it!=data[1].hist.end();++it){if(Pair(it->first))check.Fail("producer contains pair hist "+it->first);qa.Add(Compare(Key(data[0],it->first),it->second,it->first,"producer_QA",check,absoluteTolerance));}
    if(qa.compared!=108)check.Fail("producer QA count is not 108");
    for(std::map<std::string,TKey*>::const_iterator it=data[0].hist.begin();it!=data[0].hist.end();++it)if(Pair(it->first)) {
      if(!hasPrePidControl&&PrePidPair(it->first)){omittedPrePidHistograms.push_back(it->first);continue;}
      TKey* replay=Key(data[2],it->first);TKey* remix=Key(data[3],it->first);const Counter result=Compare(it->second,replay,it->first,"pair_closure",check,absoluteTolerance);pairs.Add(result);groups[Category(it->first)].Add(result);
      if(!Mixed(it->first))same.Add(Compare(replay,remix,it->first,"ep1_SE_invariance",check,absoluteTolerance));
      else if(replay&&remix){TH1* a=dynamic_cast<TH1*>(replay->ReadObj());TH1* b=dynamic_cast<TH1*>(remix->ReadObj());if(a&&b&&Cells(a)==Cells(b))for(Int_t bin=0;bin<Cells(a);++bin)if(Different(a->GetBinContent(bin),b->GetBinContent(bin),absoluteTolerance))++changedMixed;delete a;delete b;}
      else check.Fail("EP1 mixed hist missing "+it->first);
    }
    if(!pairs.compared)check.Fail("no pair histograms compared");tree=InspectTree(files[1],producerPath,check);
    if(!hasPrePidControl&&tree.controls!=0)check.Fail("compact tree contains nonselected prePID control rows");
    planes.Add(PlaneSums(data[2],"replay_EP_sum",check,absoluteTolerance,requirePlaneHistograms,replayPlanes));
    planes.Add(PlaneSums(data[3],"ep1_EP_sum",check,absoluteTolerance,requirePlaneHistograms,ep1Planes));
    if(Parameter(files[0],"inputEvents")!=Parameter(files[1],"inputEvents")||Parameter(files[0],"reconstructedEvents")!=tree.entries)check.Fail("baseline/producer input scopes differ");
    if(Parameter(files[2],"inputEvents")!=tree.entries||Parameter(files[3],"inputEvents")!=tree.entries)check.Fail("replay did not cover all stored events");
    replayBins=MixingBins(files[2]);ep1Bins=MixingBins(files[3]);if(replayBins!=8||ep1Bins!=1)check.Fail("EP8 replay/EP1 remix provenance not confirmed");
    counts[0]=Entries(data[0],"hKstarSE_lambda_he4_signal");counts[1]=Entries(data[0],"hKstarME_lambda_he4_signal");counts[2]=Entries(data[2],"hKstarSE_lambda_he4_signal");counts[3]=Entries(data[2],"hKstarME_lambda_he4_signal");counts[4]=Entries(data[3],"hKstarSE_lambda_he4_signal");counts[5]=Entries(data[3],"hKstarME_lambda_he4_signal");meaningful=true;for(int i=0;i<6;++i)meaningful=meaningful&&counts[i]>0;
    if(!meaningful)check.Fail("INSUFFICIENT_STATISTICS: nonempty He4 signal SE/ME required; empty closure is not meaningful");
  }catch(const std::exception& e){check.Fail(e.what());}
  summary<<std::setprecision(17)<<"{\n\"passed\":"<<(check.errors.empty()?"true":"false")<<",\n\"scope\":\"sample closure; not a full-production quota guarantee\",\n\"absoluteTolerance\":"<<absoluteTolerance<<",\n\"includesAllUnderflowOverflowCells\":true,\n\"meaningfulNonemptyHe4PairValidation\":"<<(meaningful?"true":"false")<<",\n\"producerQA\":";PrintCounter(summary,qa);summary<<",\n\"pairs\":";PrintCounter(summary,pairs);summary<<",\n\"ep1SameEventInvariance\":";PrintCounter(summary,same);
  summary<<",\n\"hasPrePidControl\":"<<(hasPrePidControl?"true":"false")<<",\n\"nonselectedPrePidPairComparison\":"<<Json(hasPrePidControl?"required full-control closure":"outside compact stored envelope; exact PrePid group omitted")<<",\n\"selectedPidRigidityAndOrdinaryPairsStillCompared\":true,\n\"omittedPrePidPairHistograms\":[";
  for(size_t i=0;i<omittedPrePidHistograms.size();++i){if(i)summary<<',';summary<<Json(omittedPrePidHistograms[i]);}summary<<"]";
  summary<<",\n\"planeHistogramsRequired\":"<<(requirePlaneHistograms?"true":"false")<<",\n\"planeHistogramsAvailableInReplay\":"<<(replayPlanes?"true":"false")<<",\n\"planeHistogramsAvailableInEP1\":"<<(ep1Planes?"true":"false")<<",\n\"planeSumClosure\":";PrintCounter(summary,planes);
  summary<<",\n\"pairCategories\":{";bool first=true;for(std::map<std::string,Counter>::const_iterator it=groups.begin();it!=groups.end();++it){if(!first)summary<<',';first=false;summary<<Json(it->first)<<':';PrintCounter(summary,it->second);}summary<<"},\n\"he4SignalPairEntries\":{\"baselineSE\":"<<counts[0]<<",\"baselineME\":"<<counts[1]<<",\"replaySE\":"<<counts[2]<<",\"replayME\":"<<counts[3]<<",\"ep1SE\":"<<counts[4]<<",\"ep1ME\":"<<counts[5]<<"},\n\"replayEPBins\":"<<replayBins<<",\n\"ep1EPBins\":"<<ep1Bins<<",\n\"ep1MixedContentBinsChanged\":"<<changedMixed;
  summary<<",\n\"tree\":{\"entries\":"<<tree.entries<<",\"declaredMetadataEntries\":"<<tree.declared<<",\"savedRawTracks\":"<<tree.tracks<<",\"lambda\":"<<tree.lambda<<",\"he3\":"<<tree.he3<<",\"he4\":"<<tree.he4<<",\"prePidControls\":"<<tree.controls<<",\"invalidEPEvents\":"<<tree.invalidEP<<",\"emptyTargetEvents\":"<<tree.emptyTarget<<",\"dtOnlyPoolAdmissionEvents\":"<<tree.dtOnly<<",\"poolEligibleEvents\":"<<tree.pool<<",\"savedTracksWithoutTof\":"<<tree.noTof<<",\"slotErrors\":"<<tree.slots<<",\"extraTrackOrRoleErrors\":"<<tree.extra<<",\"physicalVectorErrors\":"<<tree.physical<<",\"admissionErrors\":"<<tree.admission<<",\"sourceProvenanceErrors\":"<<tree.source<<",\"unobservedEdgeCasesAreNotValidated\":true},\n\"files\":{\n";
  for(int i=0;i<4;++i){if(i)summary<<",\n";summary<<Json(labels[i])<<":{\"path\":"<<Json(paths[i])<<",\"logicalBytes\":"<<data[i].logical<<",\"allocatedBytes\":"<<data[i].allocated<<",\"keyCount\":"<<data[i].keys<<",\"histogramCount\":"<<data[i].hist.size()<<",\"summedKeyBytes\":"<<data[i].keyBytes<<",\"histogramKeyBytes\":"<<data[i].histBytes<<",\"treeZipBytes\":"<<data[i].zip<<",\"treeTotalBytes\":"<<data[i].total<<'}';}
  summary<<"\n},\n\"keySizesTsv\":"<<Json(keyPath)<<",\n\"differencesTsv\":"<<Json(diffPath)<<",\n\"differenceDetailRows\":"<<check.rows<<",\n\"errors\":[";for(size_t i=0;i<check.errors.size();++i){if(i)summary<<',';summary<<Json(check.errors[i]);}summary<<"]\n}\n";summary.close();keys.close();differences.close();
  for(int i=0;i<4;++i)if(files[i]){files[i]->Close();delete files[i];}TH1::AddDirectory(addDirectory);const bool ok=check.errors.empty();std::cout<<(ok?"PASS":"FAIL")<<" Lambda reduced-tree verification: "<<summaryPath<<std::endl;gSystem->Exit(ok?0:1);
}
