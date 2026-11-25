#include "CombineEarlyCF.C"

void CheckBaselineLowY(const char* list) {
  using namespace epdEarlyCF;
  try {
    TH2D se("baselineSE","",200,0.,1.,9,-.5,8.5);
    TH2D me("baselineME","",200,0.,1.,9,-.5,8.5);
    se.Sumw2();me.Sumw2();Scope scope;
    long long saved=0,valid=0;
    Read(list,false,scope,se,me,saved,valid);
    const PairCounts norm=Sum(se,me,2,8,81,100);
    const PairCounts high=Sum(se,me,2,8,121,200);
    const double c=CF(high,norm.se/norm.me);
    Check(saved==1819447887LL && valid==1736875306LL && norm.se==46968. &&
          std::fabs(c-.9777243956034068)<1e-12,"previous result did not reproduce");
    std::cout<<std::setprecision(14)<<"EARLY_COMBINATION_BASELINE_PASS saved="<<saved
             <<" valid="<<valid<<" normSE="<<norm.se<<" highCF="<<c<<std::endl;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<std::endl;gSystem->Exit(2);
  }
}
