#include "../femto_lambda_period_mix100_20261119_001/ComparePeriodMix.C"
#include <TH3.h>
namespace EpdMatrixComparison {
using namespace PeriodMixComparison;
void AllHistograms(TDirectory* directory,const std::string& prefix,std::map<std::string,TKey*>& out) {
 TIter next(directory->GetListOfKeys());TKey* key=0;std::set<std::string> seen;
 while((key=dynamic_cast<TKey*>(next()))) {
  const std::string name=key->GetName();if(!seen.insert(name).second)continue;
  key=directory->GetKey(name.c_str());TClass* type=TClass::GetClass(key->GetClassName());Need(type!=0,"unknown class");
  if(type->InheritsFrom(TDirectory::Class()))AllHistograms(directory->GetDirectory(name.c_str()),prefix+name+"/",out);
  else if(type->InheritsFrom(TH1::Class()))Need(out.insert(std::make_pair(prefix+name,key)).second,"duplicate histogram");
 }
}
int IntParameter(TFile* f,const char* key) {TParameter<Int_t>* p=dynamic_cast<TParameter<Int_t>*>(f->Get(key));Need(p!=0,std::string("missing int parameter ")+key);return p->GetVal();}
void SameLedgerAll(TFile** files,int count,Long64_t expected,std::vector<Long64_t>& starts) {
 const char* fields[]={"inputIndex","runId","eventId","sourceInputIndex","candidateStoreStream","status","cent9","selectedLambda"};
 TTree* trees[9];Int_t values[9][8];int previous=-1;
 for(int f=0;f<count;++f){trees[f]=dynamic_cast<TTree*>(files[f]->Get("eventLedger"));Need(trees[f]&&trees[f]->GetEntries()==expected,"ledger count mismatch");trees[f]->SetBranchStatus("*",0);for(int j=0;j<8;++j){Need(trees[f]->GetBranch(fields[j])!=0,"missing ledger field");trees[f]->SetBranchStatus(fields[j],1);Need(trees[f]->SetBranchAddress(fields[j],&values[f][j])>=0,"bind ledger failed");}}
 for(Long64_t row=0;row<expected;++row){
  for(int f=0;f<count;++f)Need(trees[f]->GetEntry(row)>0&&!rootError,"ledger read failed");
  for(int f=1;f<count;++f)for(int j=0;j<8;++j)Need(values[0][j]==values[f][j],"full scalar ledger identity differs");
  Need(values[0][0]==row+1&&values[0][5]==4,"ledger completion/index differs");
  if(values[0][4]!=previous){Need(values[0][4]==previous+1,"noncontiguous candidate stream");starts.push_back(row+1);previous=values[0][4];}
 }
 for(int f=0;f<count;++f){trees[f]->ResetBranchAddresses();trees[f]->SetBranchStatus("*",1);}
 Need(starts.size()==3&&starts[0]==1&&starts[1]==182434&&starts[2]==357248,"expected first two segment boundaries absent");
}
double Variance(TH1* h,int bin) {return h->GetSumw2N()?h->GetSumw2()->At(bin):h->GetBinContent(bin);}
Long64_t Closure(TFile* f,int n) {
 const char* mode[]={"SE","ME"};const char* species[]={"he4","he3"};const char* suffix[]={"","_signal","_leftSB","_rightSB"};Long64_t bins=0;
 for(int m=0;m<2;++m)for(int s=0;s<2;++s)for(int t=0;t<4;++t){
  std::string pre=std::string("hKstar")+mode[m];std::string tail=std::string("_lambda_")+species[s]+suffix[t];
  std::string planeName=pre+"VsEpdBinVsCent"+tail, inclusiveName=pre+"VsCent"+tail;
  TH3* a=dynamic_cast<TH3*>(f->Get(planeName.c_str()));TH2* b=dynamic_cast<TH2*>(f->Get(inclusiveName.c_str()));
  Need(a&&b,"missing plane/inclusive histogram "+planeName);Need(a->GetNbinsY()==n&&a->GetYaxis()->GetXmin()==-.5&&a->GetYaxis()->GetXmax()==n-.5,"plane axis wrong");
  Need(Axis(a->GetXaxis(),b->GetXaxis())&&Axis(a->GetZaxis(),b->GetYaxis())&&a->GetEntries()==b->GetEntries(),"closure axes/entries mismatch");
  for(int x=0;x<=a->GetNbinsX()+1;++x)for(int z=0;z<=a->GetNbinsZ()+1;++z){double sum=0,var=0;for(int y=0;y<=n+1;++y){int bin=a->GetBin(x,y,z);double v=a->GetBinContent(bin);Need(TMath::Finite(v)&&v>=0,"bad plane content");if(y==0||y==n+1)Need(v==0&&Variance(a,bin)==0,"plane flow nonzero");sum+=v;var+=Variance(a,bin);}int bin=b->GetBin(x,z);Need(sum==b->GetBinContent(bin)&&var==Variance(b,bin),"plane content/variance closure mismatch "+planeName);++bins;}
 }
 return bins;
}
void PoolQA(TFile* f,int n,Long64_t valid) {
 TH2* a=dynamic_cast<TH2*>(f->Get("hEpdPoolBinVsCent"));TH2* b=dynamic_cast<TH2*>(f->Get("hEpdMixingBinVsCent"));Need(a&&b,"missing EP pool/canonical QA");
 Need(a->GetNbinsX()==n&&a->GetNbinsY()==9&&b->GetNbinsX()==8&&b->GetNbinsY()==9,"pool/canonical axes invalid");
 Need(a->GetEntries()==valid&&b->GetEntries()==valid,"pool/canonical EPvalid population differs");
 for(int z=0;z<11;++z){if(n==1){double sum=0;for(int k=1;k<=8;++k)sum+=b->GetBinContent(k,z);Need(a->GetBinContent(1,z)==sum,"EP1 valid selection differs");}
  else for(int k=0;k<8;++k){double sum=0;for(int j=0;j<n/8;++j)sum+=a->GetBinContent(k*(n/8)+j+1,z);Need(sum==b->GetBinContent(k+1,z),"finer EP pools not nested in canonical8");}}
}
} // namespace
void CompareEpdMatrix(const char* study="rootfile/femto_lambda_ep_merging_systematics_20261120_001",const char* reportPath="rootfile/femto_lambda_ep_merging_systematics_20261120_001/pilot_ep_comparison.json") {
 using namespace EpdMatrixComparison;using namespace PeriodMixComparison;rootError=false;TH1::AddDirectory(kFALSE);SetErrorHandler(ErrorHandler);
 if(!gSystem->AccessPathName(reportPath)){std::cerr<<"Refuse existing report"<<std::endl;gSystem->Exit(2);return;}
 std::ofstream report(reportPath);TFile* files[9]={0};std::map<std::string,TKey*> hist[9];std::ostringstream body;body<<std::setprecision(17);
 try {
  const int ns[]={1,1,8,8,16,16,32,32};const char* ids[]={"ep1_mergingOff","ep1_mergingOn","ep8_mergingOff","ep8_mergingOn","ep16_mergingOff","ep16_mergingOn","ep32_mergingOff","ep32_mergingOn"};
  std::string selection,input,defs;std::map<std::string,std::string> mixingBase;Long64_t epvalid=-1,seen=-1,rejected=-1,lambda=-1;std::vector<Long64_t> closureCounts;
  body<<"\"expectedEvents\":400000,\"cases\":[";
  for(int i=0;i<9;++i){std::string path=i==8?"rootfile/femto_lambda_period_mix100_20261119_001/pilot_B_period100.root":std::string(study)+"/pilots/"+ids[i]+".root";files[i]=TFile::Open(path.c_str(),"READ");Need(files[i]&&!files[i]->IsZombie()&&!files[i]->TestBit(TFile::kRecovered)&&!rootError,"bad ROOT file");
   Need(Object(files[i],"FemtoLambdaRunStatus")=="completed"&&Object(files[i],"FemtoLambdaReducedTreeMode")=="replay","producer status/mode wrong");
   Need(Parameter(files[i],"inputEvents")==400000&&Parameter(files[i],"reconstructedEvents")==400000,"event count wrong");
   if(!i){epvalid=Parameter(files[i],"epdEventsValid");seen=Parameter(files[i],"epdEventsSeen");rejected=Parameter(files[i],"epdEventsRejected");lambda=Parameter(files[i],"selectedLambdaCandidates");selection=Object(files[i],"FemtoLambdaReducedReplaySelection");input=Object(files[i],"FemtoLambdaReducedTreeInput");}
   Need(Parameter(files[i],"epdEventsValid")==epvalid&&Parameter(files[i],"epdEventsSeen")==seen&&Parameter(files[i],"epdEventsRejected")==rejected&&Parameter(files[i],"selectedLambdaCandidates")==lambda,"same EP-valid population/candidates not retained");
   Need(Object(files[i],"FemtoLambdaReducedReplaySelection")==selection&&Object(files[i],"FemtoLambdaReducedTreeInput")==input,"selection/input differs");
   std::map<std::string,std::string> s=Scalars(Object(files[i],"FemtoLambdaReducedTreeSettings")),m=Scalars(Object(files[i],"FemtoLambdaReducedReplayMixing"));
   Need(s["mode"]=="replay"&&s["requireOriginalMixing"]=="false"&&s["resetPoolAtSegmentBoundary"]=="false","replay boundary setting wrong");
   Need(m["bufferSize"]=="100"&&m["maxMixEvents"]=="100"&&m["mixingMode"]=="bufferAll"&&m["mixBothDirections"]=="true","pool settings wrong");
   std::ostringstream nvalue;nvalue<<(i==8?8:ns[i]);Need(m["nEventPlaneBins"]==nvalue.str(),"EP mixing count metadata wrong");m.erase("nEventPlaneBins");if(!i)mixingBase=m;Need(m==mixingBase,"mixing settings differ beyond N");
   AllHistograms(files[i],"",hist[i]);Need(!hist[i].empty(),"empty histogram inventory");
   if(i<8){Need(IntParameter(files[i],"FemtoLambdaPairEpdBinCount")==ns[i]&&IntParameter(files[i],"FemtoLambdaMergingPairCutEnabled")==i%2,"new metadata mismatch");PoolQA(files[i],ns[i],epvalid);closureCounts.push_back(Closure(files[i],ns[i]));if(i)body<<',';body<<"{\"caseId\":"<<Json(ids[i])<<",\"epBins\":"<<ns[i]<<",\"histograms\":"<<hist[i].size()<<",\"planeClosureBins\":"<<closureCounts.back()<<"}";}
  }
  body<<"],\"epdEventsSeen\":"<<seen<<",\"epdEventsValid\":"<<epvalid<<",\"epdEventsRejected\":"<<rejected;
  std::vector<Long64_t> starts;SameLedgerAll(files,9,400000,starts);body<<",\"fullScalarLedgerMatchesAllNineFiles\":true,\"segmentStartEvents\":[1,182434,357248]";
  Long64_t oldbins=0;int oldcount=0,oldnonempty=0;
  for(std::map<std::string,TKey*>::const_iterator h=hist[8].begin();h!=hist[8].end();++h){Need(hist[2].count(h->first)!=0,"old histogram missing from EP8Off: "+h->first);Equal(h->second,hist[2][h->first],h->first,oldbins);++oldcount;TH1* v=dynamic_cast<TH1*>(h->second->ReadObj());if(v->GetEntries()>0)++oldnonempty;delete v;}
  Need(oldnonempty>0,"old baseline comparison empty");body<<",\"oldEP8Off\":{\"allOldHistogramsExact\":true,\"histograms\":"<<oldcount<<",\"nonempty\":"<<oldnonempty<<",\"flowInclusiveBins\":"<<oldbins<<"}";
  Long64_t sebins=0;int sec=0;
  for(int cut=0;cut<2;++cut)for(std::map<std::string,TKey*>::const_iterator h=hist[cut].begin();h!=hist[cut].end();++h){const std::string base=Base(h->first);bool physical=Starts(base,"hKstar")||Starts(base,"hQlab")||Starts(base,"hNuclearPidPair");if(!physical||Mixed(h->first)||base.find("VsEpdBinVsCent")!=std::string::npos)continue;for(int p=1;p<4;++p){int i=2*p+cut;Need(hist[i].count(h->first)!=0,"missing invariant SE hist");Equal(h->second,hist[i][h->first],h->first,sebins);}++sec;}
  Need(sec>0,"cross EP SE check empty");body<<",\"crossEpSameCutSE\":{\"exact\":true,\"histogramsPerCutSummed\":"<<sec<<",\"flowInclusiveBinComparisons\":"<<sebins<<"},\"sampler\":[";
  for(int i=0;i<8;++i){TH2* h=dynamic_cast<TH2*>(files[i]->Get("hMixSamplerQA"));Need(h&&h->GetNbinsX()==15,"samplerQA missing");double q[15]={0};for(int j=0;j<15;++j)for(int y=1;y<=h->GetNbinsY();++y)q[j]+=h->GetBinContent(h->GetXaxis()->FindBin(j),y);Need(q[0]==q[1]&&q[4]==0&&q[2]==q[9]+q[10]&&q[1]==q[11]+q[12],"sampler accounting failed");Need(q[9]>0&&q[10]>0,"empty ME direction");if(i)body<<',';body<<"{\"caseId\":"<<Json(ids[i])<<",\"attempted\":"<<q[0]<<",\"eligible\":"<<q[1]<<",\"filled\":"<<q[2]<<",\"pairCutSkipped\":"<<q[5]<<",\"capSkipped\":"<<q[4]<<",\"filledForward\":"<<q[9]<<",\"filledReverse\":"<<q[10]<<"}";}
  body<<"]";Need(!rootError,"ROOT error during comparison");report<<"{\"passed\":true,"<<body.str()<<"}\n";report.close();Need(bool(report),"cannot write result");for(int i=0;i<9;++i){files[i]->Close();delete files[i];}SetErrorHandler(DefaultErrorHandler);std::cout<<"EP_MATRIX_COMPARISON_PASS "<<reportPath<<std::endl;
 }catch(const std::exception& e){report<<"{\"passed\":false,\"error\":"<<Json(e.what())<<"}\n";report.close();std::cerr<<"EP_MATRIX_COMPARISON_FAILED "<<e.what()<<std::endl;SetErrorHandler(DefaultErrorHandler);gSystem->Exit(1);}
}
