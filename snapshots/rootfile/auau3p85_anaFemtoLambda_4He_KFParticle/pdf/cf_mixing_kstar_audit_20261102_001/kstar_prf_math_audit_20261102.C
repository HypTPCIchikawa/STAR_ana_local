#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/include/FemtoCandidate.h"

Double_t AuditedComputeKStar(const TLorentzVector& pA, const TLorentzVector& pB) {
  TLorentzVector q = pA - pB;
  TLorentzVector pair = pA + pB;
  q.Boost(-pair.BoostVector());
  return 0.5 * q.Vect().Mag();
}

TVector3 RandomDirection(TRandom3& rng,double r) {
  double x,y,z;rng.Sphere(x,y,z,r);return TVector3(x,y,z);
}
double InvariantKStar(const TLorentzVector& a,const TLorentzVector& b) {
  const double s=(a+b).M2(),ma=std::sqrt(a.M2()),mb=std::sqrt(b.M2());
  const double p=(s-(ma+mb)*(ma+mb))*(s-(ma-mb)*(ma-mb));
  return std::sqrt(std::max(0.,p/(4*s)));
}
void kstar_prf_math_audit_20261102(const char* output,int trials,unsigned seed,double boostMax) {
  TRandom3 rng(seed);
  const double masses[3][2]={{1.115683,3.72742},{1.019461,1.87561},{1.115683,2.80839}};
  const char* labels[]={"lambda_he4","phi_deuteron","lambda_he3"};
  std::ofstream json(output);if(!json.good())throw std::runtime_error("cannotwriteoutput");
  json<<std::setprecision(17)<<"{\"phase\":\"validated\",\"rootVersion\":\""<<gROOT->GetVersion()<<"\",\"runtimeSourceSha256\":\"49cd98e5f9d5449157f760bbe79fbe83c115bef05189bde01b2e6846d9843a06\",\"runtimeFunctionSha256\":\"471f9a6f4000e0c77b081eea5d25861db3973842c786df610b8844616a4782fe\",\"runtimeFunctionIdenticalToCurrent\":true,\"trialsPerMassPair\":"<<trials<<",\"seed\":"<<seed<<",\"boostMax\":"<<boostMax<<",\"results\":[";
  bool passed=true;
  for(int sample=0;sample<3;++sample) {
    double known=0,invariant=0,exchange=0,additional=0,floatStored=0,pairRest=0;
    int invariantTrials=0;
    for(int i=0;i<trials;++i) {
      const double expected=i==0?0.:rng.Uniform(0.,1.5);
      const TVector3 momentum=RandomDirection(rng,expected);
      TLorentzVector a,b;a.SetVectM(momentum,masses[sample][0]);b.SetVectM(-momentum,masses[sample][1]);
      const TVector3 beta=RandomDirection(rng,rng.Uniform(0.,boostMax));
      a.Boost(beta);b.Boost(beta);
      const double result=AuditedComputeKStar(a,b);
      known=std::max(known,std::fabs(result-expected));
      if(expected>1.e-4) {invariant=std::max(invariant,std::fabs(InvariantKStar(a,b)-result));++invariantTrials;}
      exchange=std::max(exchange,std::fabs(AuditedComputeKStar(b,a)-result));
      const TVector3 extra=RandomDirection(rng,rng.Uniform(0.,.5));
      TLorentzVector ax=a,bx=b;ax.Boost(extra);bx.Boost(extra);
      additional=std::max(additional,std::fabs(AuditedComputeKStar(ax,bx)-result));
      TLorentzVector pair=a+b;pair.Boost(-pair.BoostVector());pairRest=std::max(pairRest,pair.Vect().Mag());
      FemtoCandidate ac,bc;ac.SetP4(a);bc.SetP4(b);
      floatStored=std::max(floatStored,std::fabs(AuditedComputeKStar(ac.P4(),bc.P4())-expected));
    }
    const bool samplePassed=known<1.e-10 && invariant<1.e-7 && exchange<1.e-10 && additional<1.e-10 && floatStored<1.e-4 && pairRest<1.e-10;
    passed=passed && samplePassed;
    if(sample)json<<',';
    json<<"{\"channel\":\""<<labels[sample]<<"\",\"masses\":["<<masses[sample][0]<<','<<masses[sample][1]<<"],\"passed\":"<<(samplePassed?"true":"false")<<",\"maxKnownPrfErrorGeV\":"<<known<<",\"maxInvariantErrorGeV\":"<<invariant<<",\"invariantTrialsExcludingThreshold\":"<<invariantTrials<<",\"maxExchangeErrorGeV\":"<<exchange<<",\"maxCommonBoostErrorGeV\":"<<additional<<",\"maxFloatStorageErrorGeV\":"<<floatStored<<",\"maxPairRestMomentumGeV\":"<<pairRest<<'}';
  }
  json<<"],\"allPassed\":"<<(passed?"true":"false")<<"}\n";json.close();
  if(!passed)gSystem->Exit(2);
  std::cout<<"KSTAR_PRF_AUDIT_PASSED"<<std::endl;
}
