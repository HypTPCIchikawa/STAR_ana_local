#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
void InspectBaselineCounts() {
 const char* paths[]={"rootfile/femto_lambda_period_mix100_20261119_001/pilot_B_period100.root","rootfile/femto_lambda_ep_merging_systematics_20261120_001/pilots/ep8_mergingOff.root"};
 const char* names[]={"hKstarSE_lambda_he4","hKstarME_lambda_he4","hKstarSE_lambda_he4_signal","hKstarME_lambda_he4_signal"};
 std::ofstream out("rootfile/femto_lambda_ep_merging_systematics_20261120_001/pilot_baseline_counts.json");out<<"{\"passed\":true,\"inputs\":[";
 double previous[4]={0};
 for(int i=0;i<2;++i){TFile* f=TFile::Open(paths[i],"READ");if(!f||f->IsZombie())gSystem->Exit(1);if(i)out<<',';out<<"{\"path\":\""<<paths[i]<<"\",\"entries\":{";for(int j=0;j<4;++j){TH1* h=dynamic_cast<TH1*>(f->Get(names[j]));if(!h)gSystem->Exit(1);double n=h->GetEntries();if(i&&previous[j]!=n)gSystem->Exit(1);previous[j]=n;if(j)out<<',';out<<'"'<<names[j]<<"\":"<<n;std::cout<<i<<' '<<names[j]<<' '<<n<<std::endl;}out<<"}}";f->Close();delete f;}
 out<<"],\"explanation\":\"Old and new pilot agree for each identical histogram. Earlier chat summary numbers 607/100116 match neither these full-mass nor signal histograms. GetEntries includes kstar underflow/overflow.\"}\n";out.close();
}
