// Independent replay control: ten pair-Y bins merged by two must equal five bins.
#include "verify_outputs.C"

namespace LambdaPairYVerification {
LambdaTreeVerification::Counter Rebin(TKey* tenKey,TKey* fiveKey,const std::string& name,
    LambdaTreeVerification::Checks& check,Double_t tolerance) {
  using namespace LambdaTreeVerification;Counter c;++c.compared;
  if(!tenKey||!fiveKey){++c.failed;check.Fail("pairY rebin missing "+name);return c;}
  TH1* source=dynamic_cast<TH1*>(tenKey->ReadObj());TH1* target=dynamic_cast<TH1*>(fiveKey->ReadObj());
  if(!source||!target){delete source;delete target;throw std::runtime_error("cannot read pairY hist "+name);}
  bool failed=false;if(source->GetEntries()>0)++c.nonempty;
  if(source->GetDimension()!=3||target->GetDimension()!=3||!SameAxis(source->GetXaxis(),target->GetXaxis())||!SameAxis(source->GetZaxis(),target->GetZaxis())||source->GetNbinsY()!=10||target->GetNbinsY()!=5||source->GetYaxis()->GetXmin()!=target->GetYaxis()->GetXmin()||source->GetYaxis()->GetXmax()!=target->GetYaxis()->GetXmax()) {++c.axis;failed=true;}
  else {
    for(int bin=1;bin<=6;++bin)if(source->GetYaxis()->GetBinLowEdge(2*bin-1)!=target->GetYaxis()->GetBinLowEdge(bin)){++c.axis;failed=true;}
    if(!failed)for(Int_t x=0;x<=source->GetNbinsX()+1;++x)for(Int_t z=0;z<=source->GetNbinsZ()+1;++z)for(Int_t y=0;y<=target->GetNbinsY()+1;++y) {
      const Int_t first=y==0?0:(y==target->GetNbinsY()+1?source->GetNbinsY()+1:2*y-1);
      const Int_t last=y==0||y==target->GetNbinsY()+1?first:first+1;
      Double_t content=0,variance=0;
      for(Int_t fine=first;fine<=last;++fine){const Int_t bin=source->GetBin(x,fine,z);content+=source->GetBinContent(bin);variance+=Variance(source,bin);}
      const Int_t targetBin=target->GetBin(x,y,z);const Double_t observed=target->GetBinContent(targetBin),observedError=target->GetBinError(targetBin),error=TMath::Sqrt(variance);
      if(Different(content,observed,tolerance)){++c.content;failed=true;check.Row("pairY_10to5",name,"content",targetBin,content,observed);}
      if(Different(error,observedError,tolerance)){++c.error;failed=true;check.Row("pairY_10to5",name,"error",targetBin,error,observedError);}
      c.maxContent=TMath::Max(c.maxContent,TMath::Abs(content-observed));c.maxError=TMath::Max(c.maxError,TMath::Abs(error-observedError));
    }
  }
  if(Different(source->GetEntries(),target->GetEntries(),tolerance)){++c.entries;failed=true;check.Row("pairY_10to5",name,"entries",-1,source->GetEntries(),target->GetEntries());}
  if(failed){++c.failed;check.Fail("pairY 10-to-5 rebin differs "+name);}delete source;delete target;return c;
}
}

void verify_pair_rapidity_rebin(const char* tenPath,const char* fivePath,const char* summaryPath,Double_t tolerance=0.) {
  using namespace LambdaTreeVerification;
  if(!tenPath||!*tenPath||!fivePath||!*fivePath||!summaryPath||!*summaryPath||!TMath::Finite(tolerance)||tolerance<0){std::cerr<<"ERROR: ten-bin/five-bin/new-summary paths are required"<<std::endl;gSystem->Exit(2);return;}
  const std::string output(summaryPath),keyPath=output+".keys.tsv",diffPath=output+".differences.tsv";
  if(!gSystem->AccessPathName(output.c_str())||!gSystem->AccessPathName(keyPath.c_str())||!gSystem->AccessPathName(diffPath.c_str())){std::cerr<<"ERROR: refusing existing rebin evidence"<<std::endl;gSystem->Exit(2);return;}
  const TString directory=gSystem->DirName(summaryPath);if(gSystem->mkdir(directory,kTRUE)&&gSystem->AccessPathName(directory)){gSystem->Exit(2);return;}
  std::ofstream summary(output.c_str()),keys(keyPath.c_str()),differences(diffPath.c_str());if(!summary||!keys||!differences){std::cerr<<"ERROR: cannot create rebin evidence"<<std::endl;gSystem->Exit(2);return;}
  keys<<"file\tkey\tclass\tkey_bytes\tserialized_object_bytes\ttree_zip_bytes\ttree_total_bytes\n";differences<<"group\tkey\tfield\tglobal_bin_including_flow\treference\tactual\n";
  Checks check;check.detail=&differences;Inventory data[2];TFile* files[2]={0,0};const char* paths[]={tenPath,fivePath};const char* labels[]={"pairY10","pairY5"};Counter rebin,inclusive;Double_t counts[4]={-1,-1,-1,-1};bool meaningful=false,sourceEqual=false;
  const Bool_t addDirectory=TH1::AddDirectoryStatus();TH1::AddDirectory(kFALSE);
  try {
    Runtime();
    for(int i=0;i<2;++i){files[i]=TFile::Open(paths[i],"READ");if(!files[i]||files[i]->IsZombie())throw std::runtime_error(std::string("cannot open ")+labels[i]);TNamed* status=dynamic_cast<TNamed*>(files[i]->Get("FemtoLambdaRunStatus"));if(!status||std::string(status->GetTitle())!="completed")check.Fail(std::string(labels[i])+" run not completed");if(MixingBins(files[i])!=8)check.Fail(std::string(labels[i])+" is not EP8 mixing");data[i].logical=files[i]->GetSize();struct stat info;if(stat(paths[i],&info)==0)data[i].allocated=Long64_t(info.st_blocks)*512;Inspect(files[i],"",labels[i],data[i],keys);}
    TObjString* source0=dynamic_cast<TObjString*>(files[0]->Get("FemtoLambdaReducedTreeInput"));TObjString* source1=dynamic_cast<TObjString*>(files[1]->Get("FemtoLambdaReducedTreeInput"));
    sourceEqual=source0&&source1&&source0->GetString()==source1->GetString();if(!sourceEqual)check.Fail("pairY replays do not identify the same stored input");
    if(Parameter(files[0],"inputEvents")!=Parameter(files[1],"inputEvents")||Parameter(files[0],"inputEvents")<=0)check.Fail("pairY replays have different or empty event scopes");
    for(std::map<std::string,TKey*>::const_iterator it=data[0].hist.begin();it!=data[0].hist.end();++it)if(Pair(it->first)) {
      if(it->first.find("VsPairYVsCent")!=std::string::npos)rebin.Add(LambdaPairYVerification::Rebin(it->second,Key(data[1],it->first),it->first,check,tolerance));
      else inclusive.Add(Compare(it->second,Key(data[1],it->first),it->first,"inclusive_pair_invariance",check,tolerance));
    }
    if(rebin.compared!=8)check.Fail("expected eight pairY SE/ME/mass-class histograms");
    if(!inclusive.compared)check.Fail("no inclusive pair invariance compared");
    counts[0]=Entries(data[0],"hKstarSE_lambda_he4_signal");counts[1]=Entries(data[0],"hKstarME_lambda_he4_signal");counts[2]=Entries(data[1],"hKstarSE_lambda_he4_signal");counts[3]=Entries(data[1],"hKstarME_lambda_he4_signal");meaningful=true;for(int i=0;i<4;++i)meaningful=meaningful&&counts[i]>0;
    if(!meaningful)check.Fail("INSUFFICIENT_STATISTICS: empty He4 signal SE/ME cannot establish pairY control");
  }catch(const std::exception& e){check.Fail(e.what());}
  summary<<std::setprecision(17)<<"{\n\"passed\":"<<(check.errors.empty()?"true":"false")<<",\n\"scope\":\"same stored sample, EP8 mixing, pairY10 merged by two versus pairY5\",\n\"absoluteTolerance\":"<<tolerance<<",\n\"includesAllUnderflowOverflowCells\":true,\n\"sameStoredInput\":"<<(sourceEqual?"true":"false")<<",\n\"meaningfulNonemptyHe4Signal\":"<<(meaningful?"true":"false")<<",\n\"pairRapidityRebin\":";PrintCounter(summary,rebin);summary<<",\n\"inclusivePairInvariance\":";PrintCounter(summary,inclusive);
  summary<<",\n\"he4SignalEntries\":{\"pairY10SE\":"<<counts[0]<<",\"pairY10ME\":"<<counts[1]<<",\"pairY5SE\":"<<counts[2]<<",\"pairY5ME\":"<<counts[3]<<"},\n\"files\":{";
  for(int i=0;i<2;++i){if(i)summary<<',';summary<<Json(labels[i])<<":{\"path\":"<<Json(paths[i])<<",\"logicalBytes\":"<<data[i].logical<<",\"allocatedBytes\":"<<data[i].allocated<<",\"histogramKeyBytes\":"<<data[i].histBytes<<",\"treeZipBytes\":"<<data[i].zip<<'}';}
  summary<<"},\n\"keySizesTsv\":"<<Json(keyPath)<<",\n\"differencesTsv\":"<<Json(diffPath)<<",\n\"errors\":[";for(size_t i=0;i<check.errors.size();++i){if(i)summary<<',';summary<<Json(check.errors[i]);}summary<<"]\n}\n";summary.close();keys.close();differences.close();
  for(int i=0;i<2;++i)if(files[i]){files[i]->Close();delete files[i];}TH1::AddDirectory(addDirectory);const bool ok=check.errors.empty();std::cout<<(ok?"PASS":"FAIL")<<" pairY10-to-5 verification: "<<summaryPath<<std::endl;gSystem->Exit(ok?0:1);
}
