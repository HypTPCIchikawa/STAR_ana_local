#include "TH1D.h"
#include "TSystem.h"
#include "TMath.h"
#include <limits>
#include <stdexcept>
#include "/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/tools/kfparticle_summary/KFSummarySupport.h"
bool CheckNormalization(TH1D* hSigTrue, TH1D* hSigMix) {
 TH1D cf("cf_guard","",10,0,1); TH1D* hCF=&cf;
 int bin1=4,bin2=7;
 double intTrue=hSigTrue->Integral(bin1,bin2),intMix=hSigMix->Integral(bin1,bin2);

bool kfNormValid = intTrue>0 && intMix>0 && TMath::Finite(intTrue) && TMath::Finite(intMix);
for(int kfBin=bin1; kfBin<=bin2; ++kfBin) {
  const double kfSE=hSigTrue->GetBinContent(kfBin), kfME=hSigMix->GetBinContent(kfBin);
  const double kfSEError=hSigTrue->GetBinError(kfBin), kfMEError=hSigMix->GetBinError(kfBin);
  if (!(kfSE>0 && kfME>0 && kfSEError>0 && kfMEError>0 &&
        TMath::Finite(kfSE) && TMath::Finite(kfME) &&
        TMath::Finite(kfSEError) && TMath::Finite(kfMEError))) kfNormValid=false;
}
KFReport::Record("normalization",hCF->GetName(),kfNormValid?0:-1,-1,false,
  kfNormValid?"complete_normalization_interval":"invalid_normalization_bin",
  intTrue,intMix,bin2-bin1+1);

 return kfNormValid;
}
void test_cf_normalization_20261030(const char* auditPath) {
 const double cfg[9]={1.08,1.15,1.11,1.12,.001,.008,.003,10,2};
 KFReport::Init(auditPath,1.11596,.00296,cfg);
 TH1D se("se","",10,0,1),me("me","",10,0,1);se.Sumw2();me.Sumw2();
 for(int b=1;b<=10;++b){se.SetBinContent(b,100);se.SetBinError(b,10);me.SetBinContent(b,1000);me.SetBinError(b,32);}
 if(!CheckNormalization(&se,&me))gSystem->Exit(1);
 // A failed normalization fit is zero/error-zero in the generated yields.
 me.SetBinContent(5,0);me.SetBinError(5,0);
 if(CheckNormalization(&se,&me))gSystem->Exit(2);
 me.SetBinContent(5,1000);me.SetBinError(5,32);
 // Missing fits outside the normalization range must not veto other points.
 me.SetBinContent(2,0);me.SetBinError(2,0);
 if(!CheckNormalization(&se,&me))gSystem->Exit(3);
 se.SetBinContent(6,-1);
 if(CheckNormalization(&se,&me))gSystem->Exit(4);
 se.SetBinContent(6,100);se.SetBinError(6,std::numeric_limits<double>::quiet_NaN());
 if(CheckNormalization(&se,&me))gSystem->Exit(5);
 std::cout<<"COMPLETE_NORMALIZATION_GUARD_PASS"<<std::endl;
}
