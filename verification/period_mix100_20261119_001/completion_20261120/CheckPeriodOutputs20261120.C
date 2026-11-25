#include "ComparePeriodMix.C"

// Check run completion and stored settings with two scalar-ledger reads per file.
void CheckPeriodOutputs20261120(const char* manifestPath,const char* reportPath) {
  using namespace PeriodMixComparison;
  rootError=false;SetErrorHandler(ErrorHandler);TH1::AddDirectory(kFALSE);
  try {
    Need(manifestPath&&reportPath&&gSystem->AccessPathName(reportPath),"explicit fresh report required");
    std::ifstream input(manifestPath);Need(bool(input),"cannot open output manifest");
    std::ostringstream body;body<<std::setprecision(17)<<"\"files\":[";
    std::string path;Long64_t expected=0,total=0;int segments=0,index=-1,count=0;
    double he4se=0,he4me=0,he3se=0,he3me=0,eligible=0,filled=0;
    while(input>>index>>expected>>segments>>path) {
      Need(index==count&&expected>0&&segments>0,"invalid manifest order/count");
      TFile* file=TFile::Open(path.c_str(),"READ");
      Need(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered)&&!rootError,"invalid/recovered output ROOT");
      Need(Object(file,"FemtoLambdaRunStatus")=="completed"&&Object(file,"FemtoLambdaReducedTreeMode")=="replay","incomplete replay status");
      Need(Parameter(file,"inputEvents")==expected&&Parameter(file,"reconstructedEvents")==expected,"output event counters differ from source period");
      const std::map<std::string,std::string> setting=Scalars(Object(file,"FemtoLambdaReducedTreeSettings"));
      const std::map<std::string,std::string> mix=Scalars(Object(file,"FemtoLambdaReducedReplayMixing"));
      Need(setting.find("resetPoolAtSegmentBoundary")!=setting.end()&&setting.find("requireOriginalMixing")!=setting.end(),"missing replay policy");
      Need(setting.find("resetPoolAtSegmentBoundary")->second=="false"&&setting.find("requireOriginalMixing")->second=="false","period pool continuation not enabled");
      Need(mix.find("bufferSize")!=mix.end()&&mix.find("maxMixEvents")!=mix.end()&&mix.find("mixingMode")!=mix.end()&&mix.find("mixBothDirections")!=mix.end(),"missing mixing policy");
      Need(mix.find("bufferSize")->second=="100"&&mix.find("maxMixEvents")->second=="100"&&mix.find("mixingMode")->second=="bufferAll"&&mix.find("mixBothDirections")->second=="true","unexpected depth/mixing policy");
      TTree* ledger=dynamic_cast<TTree*>(file->Get("eventLedger"));Need(ledger&&ledger->GetEntries()==expected,"ledger event extent mismatch");
      Int_t ordinal=0,stream=-1,status=-1;
      ledger->SetBranchStatus("*",0);
      const char* names[]={"inputIndex","candidateStoreStream","status"};Int_t* addresses[]={&ordinal,&stream,&status};
      for(int i=0;i<3;++i){ledger->SetBranchStatus(names[i],1);Need(ledger->SetBranchAddress(names[i],addresses[i])>=0,"cannot bind scalar ledger endpoint");}
      Need(ledger->GetEntry(0)>0&&ordinal==1&&stream==0&&status==4,"first event endpoint mismatch");
      Need(ledger->GetEntry(expected-1)>0&&ordinal==expected&&stream==segments-1&&status==4,"last event/segment endpoint mismatch");
      ledger->ResetBranchAddresses();
      const char* histogramNames[]={"hKstarSE_lambda_he4_signal","hKstarME_lambda_he4_signal","hKstarSE_lambda_he3_signal","hKstarME_lambda_he3_signal"};
      double pairs[4]={0,0,0,0};
      for(int h=0;h<4;++h){TH1* value=dynamic_cast<TH1*>(file->Get(histogramNames[h]));Need(value&&TMath::Finite(value->GetEntries())&&value->GetEntries()>0,"missing/empty signal pair histogram");pairs[h]=value->GetEntries();}
      TH2* qa=dynamic_cast<TH2*>(file->Get("hMixSamplerQA"));Need(qa&&qa->GetNbinsX()==15,"missing sampler QA");
      double counters[15]={0};
      for(int y=1;y<=qa->GetNbinsY();++y)for(int q=0;q<15;++q){double value=qa->GetBinContent(qa->GetXaxis()->FindBin(q),y);Need(TMath::Finite(value)&&value>=0,"invalid sampler QA value");counters[q]+=value;}
      Need(counters[4]==0&&counters[0]==counters[1]&&counters[2]>0&&counters[9]>0&&counters[10]>0,"sampler/cap/direction mismatch");
      if(count)body<<',';
      body<<"{\"subjob\":"<<index<<",\"path\":"<<Json(path)<<",\"events\":"<<expected<<",\"producerSegments\":"<<segments<<",\"he4SignalSE\":"<<pairs[0]<<",\"he4SignalME\":"<<pairs[1]<<",\"he3SignalSE\":"<<pairs[2]<<",\"he3SignalME\":"<<pairs[3]<<",\"meEligible\":"<<counters[1]<<",\"meFilled\":"<<counters[2]<<",\"capSkipped\":"<<counters[4]<<",\"periodContinuation\":true,\"bothDepthLimits\":100}";
      total+=expected;he4se+=pairs[0];he4me+=pairs[1];he3se+=pairs[2];he3me+=pairs[3];eligible+=counters[1];filled+=counters[2];
      file->Close();delete file;++count;
    }
    Need(input.eof()&&count==20&&!rootError,"incomplete output manifest/check");
    std::ofstream report(reportPath);Need(bool(report),"cannot create output check report");
    report<<std::setprecision(17)<<"{\"passed\":true,"<<body.str()<<"],\"sourcePeriods\":"<<count<<",\"storedEventsProcessed\":"<<total<<",\"he4SignalSE\":"<<he4se<<",\"he4SignalME\":"<<he4me<<",\"he3SignalSE\":"<<he3se<<",\"he3SignalME\":"<<he3me<<",\"allChannelMeEligible\":"<<eligible<<",\"allChannelMeFilled\":"<<filled<<",\"fullEventReplayPerformed\":false,\"ledgerPayloadRowsRead\":"<<2*count<<"}\n";
    report.close();Need(bool(report),"output check report write failed");SetErrorHandler(DefaultErrorHandler);
    std::cout<<"PRODUCTION_OUTPUT_CHECK_PASS files="<<count<<" events="<<total<<std::endl;
  } catch(const std::exception& error) {
    SetErrorHandler(DefaultErrorHandler);std::cerr<<"PRODUCTION_OUTPUT_CHECK_FAILED "<<error.what()<<std::endl;gSystem->Exit(1);
  }
}
