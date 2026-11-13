// Read-only histogram/counter inspection. No event loop, fits or CF calculation.
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
bool CfPrepAxis(const TAxis* a,int n,double lo,double hi){return a->GetNbins()==n&&std::fabs(a->GetXmin()-lo)<1.e-10&&std::fabs(a->GetXmax()-hi)<1.e-10;}
void check_cf_input_20261030(const char* input,const char* output,Long64_t expectedInput,Long64_t expectedReco,int expectedHist){
 if(!input||!*input||!output||!*output||!gSystem->AccessPathName(output)){std::cerr<<"Explicit input and NEW output required"<<std::endl;gSystem->Exit(2);return;}
 TFile f(input,"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered)){gSystem->Exit(2);return;}
 std::vector<std::string> names;names.push_back("hLambda_InvMass");names.push_back("hRefMult");names.push_back("hRefMultCorr");
 const char* regions[]={"","_SBNeg","_SBPos","_signalLow","_signalHigh"};const char* common[]={"signal","leftSB","rightSB","signalLow","signalHigh"};
 for(int c=0;c<9;++c)for(int m=0;m<2;++m){const char* d=m?"mix":"true",*t=m?"Mixed_":"";names.push_back(Form("%s/hKstarMass_%s4He_CentBin%d",d,t,c));for(int r=0;r<5;++r){if(r<3)names.push_back(Form("%s/hKstar_%s4He%s_CentBin%d",d,t,regions[r],c));names.push_back(Form("%s/hKstarMt_%s4He%s_CentBin%d",d,t,regions[r],c));}}
 for(int m=0;m<2;++m)for(int r=0;r<5;++r)names.push_back(Form("hKstar%sVsCent_lambda_he4_%s",m?"ME":"SE",common[r]));
 bool pass=true;int ordinary=0;const char* dirs[]={"","true","mix"};
 for(int i=0;i<3;++i){TDirectory* scanDirectory=i?f.GetDirectory(dirs[i]):&f;if(!scanDirectory){pass=false;continue;}TIter next(scanDirectory->GetListOfKeys());TKey* key=0;while((key=static_cast<TKey*>(next()))){TClass* cls=TClass::GetClass(key->GetClassName());if(cls&&cls->InheritsFrom(TH1::Class()))++ordinary;}}
 std::ofstream o(output);if(!o.good()){gSystem->Exit(2);return;}o<<std::setprecision(17)<<"{\n\"inputRoot\":\""<<input<<"\",\n\"histograms\":[\n";
 for(unsigned i=0;i<names.size();++i){if(i)o<<",\n";TH1* h=dynamic_cast<TH1*>(f.Get(names[i].c_str()));o<<"{\"name\":\""<<names[i]<<"\",\"found\":"<<(h?"true":"false");if(!h){pass=false;o<<"}";continue;}
 const bool mass=names[i].find("hKstarMass_")!=std::string::npos,kstar=names[i].find("/hKstar_")!=std::string::npos,mt=names[i].find("hKstarMt_")!=std::string::npos,cent=names[i].find("VsCent_lambda_he4_")!=std::string::npos;
 bool ok=TMath::Finite(h->GetEntries())&&h->GetEntries()>0.;if(mass||kstar||mt||cent)ok=ok&&CfPrepAxis(h->GetXaxis(),200,0.,1.);if(mass)ok=ok&&h->GetDimension()==2&&CfPrepAxis(h->GetYaxis(),200,1.05,1.25);if(kstar)ok=ok&&h->GetDimension()==1;if(mt)ok=ok&&h->GetDimension()==2&&CfPrepAxis(h->GetYaxis(),200,0.,10.);if(cent)ok=ok&&h->GetDimension()==2&&CfPrepAxis(h->GetYaxis(),9,-.5,8.5);pass=pass&&ok;
 o<<",\"class\":\""<<h->ClassName()<<"\",\"dimension\":"<<h->GetDimension()<<",\"entries\":"<<h->GetEntries()<<",\"sumw2Cells\":"<<h->GetSumw2N()<<",\"contractMatches\":"<<(ok?"true":"false")<<",\"axes\":[";const TAxis* axes[]={h->GetXaxis(),h->GetYaxis()};for(int a=0;a<h->GetDimension()&&a<2;++a){if(a)o<<",";o<<"{\"bins\":"<<axes[a]->GetNbins()<<",\"min\":"<<axes[a]->GetXmin()<<",\"max\":"<<axes[a]->GetXmax()<<"}";}o<<"]}";}
 TH1* inputs=dynamic_cast<TH1*>(f.Get("hRefMult"));TH1* recos=dynamic_cast<TH1*>(f.Get("hRefMultCorr"));TNamed* status=dynamic_cast<TNamed*>(f.Get("FemtoLambdaRunStatus"));TNamed* calibration=dynamic_cast<TNamed*>(f.Get("FemtoNuclearDedxCalibration"));TNamed* domain=dynamic_cast<TNamed*>(f.Get("FemtoNuclearPidDomain"));
 const bool counts=inputs&&recos&&inputs->GetEntries()==expectedInput&&recos->GetEntries()==expectedReco,completed=status&&std::string(status->GetTitle())=="completed",accepted=calibration&&std::string(calibration->GetTitle()).find("he4_accepted_20261022")!=std::string::npos,split=f.Get("NuclearPidPairJointSplitManifest")!=0,noJoint=f.GetKey("hNuclearPidPairJoint")==0;
 pass=pass&&ordinary==expectedHist&&counts&&completed&&accepted&&domain&&split&&noJoint;
 o<<"\n],\n\"ordinaryHistogramCount\":"<<ordinary<<",\"checkedHistogramCount\":"<<names.size()<<",\"inputEvents\":"<<(inputs?inputs->GetEntries():-1)<<",\"reconstructedEvents\":"<<(recos?recos->GetEntries():-1)<<",\"counterAgreement\":"<<(counts?"true":"false")<<",\"runCompleted\":"<<(completed?"true":"false")<<",\"acceptedHe4CalibrationRecorded\":"<<(accepted?"true":"false")<<",\"pidDomainRecorded\":"<<(domain?"true":"false")<<",\"splitManifestPresent\":"<<(split?"true":"false")<<",\"giantJointAbsentFromNormalRoot\":"<<(noJoint?"true":"false")<<",\"fitPerformed\":false,\"cfCalculated\":false,\"pass\":"<<(pass?"true":"false")<<"}\n";o.close();f.Close();std::cout<<(pass?"PASS":"FAIL")<<" CF input preparation: "<<names.size()<<" histogram contracts, "<<ordinary<<" ordinary histograms"<<std::endl;if(!pass)gSystem->Exit(1);
}
