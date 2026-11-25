// Production-helper integration test with synthetic pair fixtures, not events.
// The caller supplies the compact replay profile; both helpers omit pre-PID controls.
int test_he4_ylab_binning(const char* mainconf);
#ifndef __CINT__
#include "FemtoCandidate.h"
#include <iosfwd>
#include <string>
#include <vector>
#define private public
#include "StMaker/common/FemtoLambdaLegacy.h"
#undef private
#include "HistManager.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMemFile.h"
#include "TFile.h"
#include "TFileMerger.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

namespace He4YLabTest {
unsigned long long checks=0;unsigned calls=0;
void Need(bool pass,const char* message){++checks;if(!pass)throw std::runtime_error(message);}
int Cells(const TH1& h){return (h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1)*(h.GetDimension()>2?h.GetNbinsZ()+2:1);}
double Var(const TH1& h,int bin){return h.GetSumw2N()?h.GetSumw2()->At(bin):h.GetBinContent(bin);}
std::string Key(int mixed,int region){const char* s[]={"","_signal","_leftSB","_rightSB"};return std::string("hKstar")+(mixed?"ME":"SE")+"VsHe4YLabVsCent_lambda_he4"+s[region];}
std::string Inclusive(int mixed,int region){const char* s[]={"","_signal","_leftSB","_rightSB"};return std::string("hKstar")+(mixed?"ME":"SE")+"_lambda_he4"+s[region];}
void Equal(const TH1* a,const TH1* b){Need(a&&b&&Cells(*a)==Cells(*b),"histogram shapes match");Need(a->GetEntries()==b->GetEntries(),"entries match");for(int i=0;i<Cells(*a);++i){Need(a->GetBinContent(i)==b->GetBinContent(i),"all contents including flows match");Need(a->GetBinError(i)==b->GetBinError(i),"all errors including flows match");}}
void EqualManager(HistManager* a,HistManager* b,const std::vector<std::string>& keys,bool omitNew){if(!a&&!b)return;Need(a&&b,"same optional histogram bank");for(size_t i=0;i<keys.size();++i){if(omitNew&&keys[i].find("VsHe4YLabVsCent")!=std::string::npos)continue;Equal(a->Get(keys[i].c_str()),b->Get(keys[i].c_str()));}}
FemtoCandidate Candidate(const char* species,int event,int track,double y,double mass,double pt){FemtoCandidate c;c.speciesKey=species;c.eventIndex=event;c.source=c.speciesKey=="lambda"?kFemtoCandResonance:kFemtoCandTrack;c.charge=c.source==kFemtoCandTrack?2:0;c.trk.trackIndex=track;const double mt=std::sqrt(mass*mass+pt*pt);TLorentzVector p(pt,0.,mt*std::sinh(y),mt*std::cosh(y));c.SetP4(p);c.reso.dau1EventIndex=c.reso.dau2EventIndex=event;c.reso.dau1Index=0;c.reso.dau2Index=1;return c;}
double Y(const FemtoCandidate& c){const double e=c.energy,z=c.pz;return .5*std::log((e+z)/(e-z));}
void Projection(FemtoLambdaLegacy& h){for(int mode=0;mode<2;++mode)for(int region=0;region<4;++region){TH3* s=dynamic_cast<TH3*>(h.RootHistograms()->Get(Key(mode,region).c_str()));TH1* incl=h.RootHistograms()->Get(Inclusive(mode,region).c_str());Need(s&&incl,"inclusive and new histograms exist");double under=0,over=0;for(int x=0;x<=s->GetNbinsX()+1;++x){double all=0,var=0,regular=0,coarse=0;for(int z=0;z<=s->GetNbinsZ()+1;++z)for(int y=0;y<=s->GetNbinsY()+1;++y){const int c=s->GetBin(x,y,z);const double v=s->GetBinContent(c);all+=v;var+=Var(*s,c);if(y==0)under+=v;else if(y==s->GetNbinsY()+1)over+=v;else regular+=v;}for(int g=0;g<s->GetNbinsY()/2;++g)for(int y=2*g+1;y<=2*g+2;++y)for(int z=0;z<=s->GetNbinsZ()+1;++z)coarse+=s->GetBinContent(s->GetBin(x,y,z));Need(all==incl->GetBinContent(x),"all Y bins and flows reproduce inclusive kstar");Need(std::fabs(var-Var(*incl,x))<1.e-9,"Y projection preserves variance");Need(coarse==regular,"8-to-4 raw merging preserves all regular counts");}Need(under>0&&over>0,"both rapidity flows exercised in every mode/mass sample");for(int y=1;y<=s->GetNbinsY();++y)Need(s->Integral(0,s->GetNbinsX()+1,y,y,0,s->GetNbinsZ()+1)>0,"all regular rapidity bins exercised");}}
void Metadata(FemtoLambdaLegacy& h){TMemFile output("he4YLabMetadata.root","RECREATE");Need(h.Write(&output),"production Write succeeds");const char* names[]={"he4RapidityLabNBins","he4RapidityLabMin","he4RapidityLabMax"};TParameter<int>* n=dynamic_cast<TParameter<int>*>(output.Get(names[0]));TParameter<double>* lo=dynamic_cast<TParameter<double>*>(output.Get(names[1]));TParameter<double>* hi=dynamic_cast<TParameter<double>*>(output.Get(names[2]));TNamed* def=dynamic_cast<TNamed*>(output.Get("FemtoHe4RapidityBinningDefinition"));Need(n&&n->GetVal()==h.mHe4RapidityNBins&&n->IsConst(),"bin metadata typed and first-value");Need(lo&&lo->GetVal()==h.mHe4RapidityMin&&lo->IsConst()&&hi&&hi->GetVal()==h.mHe4RapidityMax&&hi->IsConst(),"range metadata typed and first-value");Need(def&&std::string(def->GetTitle()).find("He4YLab")!=std::string::npos&&std::string(def->GetTitle()).find("nucleus candidate's own P4")!=std::string::npos,"metadata records own-nucleus definition");const TString dir=TString::Format("%s/rootfile/lambda_he4_ylab_split_20261121_001/fixture_metadata_%d",gSystem->WorkingDirectory(),gSystem->GetPid());Need(gSystem->mkdir(dir,kTRUE)==0,"fresh metadata fixture directory");for(int i=0;i<2;++i){const TString path=TString::Format("%s/worker%d.root",dir.Data(),i);TFile f(path,"CREATE");Need(!f.IsZombie(),"metadata worker created");for(int j=0;j<3;++j){f.cd();Need(output.Get(names[j])->Write()>0,"parameter copied");}f.cd();Need(output.Get(Key(0,1).c_str())->Write()>0,"populated production histogram copied");f.Close();}const TString path=dir+"/merged.root";{TFileMerger merger;Need(merger.OutputFile(path,"RECREATE")&&merger.AddFile(dir+"/worker0.root")&&merger.AddFile(dir+"/worker1.root")&&merger.Merge(),"metadata workers merge");}TFile merged(path,"READ");n=dynamic_cast<TParameter<int>*>(merged.Get(names[0]));lo=dynamic_cast<TParameter<double>*>(merged.Get(names[1]));hi=dynamic_cast<TParameter<double>*>(merged.Get(names[2]));Need(n&&n->GetVal()==h.mHe4RapidityNBins&&n->IsConst()&&lo&&lo->GetVal()==h.mHe4RapidityMin&&lo->IsConst()&&hi&&hi->GetVal()==h.mHe4RapidityMax&&hi->IsConst(),"hadd preserves bin/range constants");TH1* before=dynamic_cast<TH1*>(output.Get(Key(0,1).c_str()));TH1* after=dynamic_cast<TH1*>(merged.Get(Key(0,1).c_str()));Need(before&&after&&after->GetEntries()==2*before->GetEntries(),"hadd adds histogram entries");for(int i=0;i<Cells(*before);++i)Need(after->GetBinContent(i)==2*before->GetBinContent(i),"hadd adds double-precision cells including flow");}
void Run(const char* mainconf){FemtoLambdaLegacy on,off;Need(on.Init(mainconf,std::cerr,false,false)&&off.Init(mainconf,std::cerr,false,false),"helpers initialize with explicit mainconf");Need(on.mHe4RapidityBinningEnabled&&on.mHe4RapidityNBins==8,"requested eight-bin profile enabled");off.mHe4RapidityBinningEnabled=false;on.mPairGeometryConfigured=off.mPairGeometryConfigured=false; // Isolate the additive diagnostic from separately tested track geometry.
std::map<std::string,TH3*> expected;for(int m=0;m<2;++m)for(int r=0;r<4;++r){TH3* a=dynamic_cast<TH3*>(on.RootHistograms()->Get(Key(m,r).c_str()));Need(a&&a->IsA()==TH3D::Class(),"all eight histograms double precision");expected[Key(m,r)]=static_cast<TH3*>(a->Clone(("expected_"+Key(m,r)).c_str()));expected[Key(m,r)]->SetDirectory(0);expected[Key(m,r)]->Reset();}
const double width=(on.mHe4RapidityMax-on.mHe4RapidityMin)/on.mHe4RapidityNBins;std::vector<double> ys;for(int i=0;i<on.mHe4RapidityNBins;++i)ys.push_back(on.mHe4RapidityMin+(i+.5)*width);ys.push_back(on.mHe4RapidityMin-width);ys.push_back(on.mHe4RapidityMax+width);ys.push_back(0.);const double masses[]={on.mMean+(on.mOuterFactor+2)*on.mWindow,on.mMean,on.mMean-.5*(1+on.mOuterFactor)*on.mWindow,on.mMean+.5*(1+on.mOuterFactor)*on.mWindow};const int cent[]={0,4,8};TH1* xr=on.RootHistograms()->Get(Inclusive(0,0).c_str());const double xlo=xr->GetXaxis()->GetXmin(),xhi=xr->GetXaxis()->GetXmax(),dx=xr->GetXaxis()->GetBinWidth(1);const double ks[]={xlo-dx,xlo+.237*(xhi-xlo),xhi+dx};const double hm=on.NuclearP4(TVector3(),3).M();
for(int direction=0;direction<3;++direction)for(size_t yi=0;yi<ys.size();++yi)for(int region=0;region<4;++region)for(int ci=0;ci<3;++ci)for(int xi=0;xi<3;++xi){const bool mixed=direction!=0;const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;FemtoCandidate l=Candidate("lambda",le,-1,3.,on.mLambdaPairMass,1.2);l.reso.invMass=masses[region];FemtoCandidate n=Candidate("he4",ne,mixed?0:2,ys[yi],hm,1.7);n.y=yi%2?999.f:-999.f;Need(on.MassRegion(l.reso.invMass)==region,"fixture mass falls in intended window");Need(on.FillPair(l,n,ks[xi],.2,cent[ci],mixed)&&off.FillPair(l,n,ks[xi],.2,cent[ci],mixed),"SE and both ME directions keep identical pairs");++calls;const double y=Y(n);Need(std::fabs(y-n.P4().Rapidity())<1.e-12,"independent log-energy rapidity equals stored P4");Need(y!=n.y,"cached y sentinel ignored");for(int r=0;r<4;++r)if(r==0||r==region)expected[Key(mixed,r)]->Fill(ks[xi],y,cent[ci]);if(yi<8&&xi==1)Need((l.P4()+n.P4()).Rapidity()>on.mHe4RapidityMax,"fixture pair rapidity differs from single He4 rapidity");}
for(std::map<std::string,TH3*>::iterator i=expected.begin();i!=expected.end();++i){Equal(on.RootHistograms()->Get(i->first.c_str()),i->second);delete i->second;}
Projection(on);
FemtoCandidate l=Candidate("lambda",30,-1,1.,on.mLambdaPairMass,1.);l.reso.invMass=on.mMean;FemtoCandidate n=Candidate("he4",30,2,on.mHe4RapidityMin+.5*width,hm,1.);const double before=on.RootHistograms()->Get(Key(0,1).c_str())->GetEntries();n.energy=std::fabs(n.pz);Need(on.FillPair(l,n,.37,.2,4,false)&&off.FillPair(l,n,.37,.2,4,false),"invalid rapidity never vetoes physical pair");Need(on.RootHistograms()->Get(Key(0,1).c_str())->GetEntries()==before,"invalid rapidity skips only new diagnostic");
if(on.IsPairNucleus(2)){n=Candidate("he3",30,2,-.7,on.NuclearP4(TVector3(),2).M(),1.);Need(on.FillPair(l,n,.37,.2,4,false)&&off.FillPair(l,n,.37,.2,4,false),"He3 and mapped control retained");Need(on.RootHistograms()->Get(Key(0,1).c_str())->GetEntries()==before,"He3AsHe4 does not enter physical He4 view");}
EqualManager(on.RootHistograms(),off.RootHistograms(),on.mRootKeys,true);EqualManager(on.NuclearHistograms(),off.NuclearHistograms(),on.mNuclearKeys,false);EqualManager(on.NuclearPidPairQaHistograms(),off.NuclearPidPairQaHistograms(),on.mNuclearPidPairQaKeys,false);
TH3* edges=dynamic_cast<TH3*>(on.RootHistograms()->Get(Key(0,1).c_str()));const double inf=std::numeric_limits<double>::infinity();Need(edges->GetYaxis()->FindFixBin(on.mHe4RapidityMin)==1&&edges->GetYaxis()->FindFixBin(on.mHe4RapidityMax)==on.mHe4RapidityNBins+1,"native half-open global endpoints");for(int edge=0;edge<=on.mHe4RapidityNBins;++edge){const double boundary=on.mHe4RapidityMin+edge*width;Need(edges->GetYaxis()->FindFixBin(boundary)==edge+1,"exact axis boundary enters upper native interval");for(int side=0;side<2;++side){FemtoCandidate t=Candidate("he4",31,2,boundary,hm,1.);if(boundary==0.)t.pz=(side?1.f:-1.f)*float(64.*std::numeric_limits<double>::epsilon()*std::max(1.,double(t.energy)));for(unsigned i=0;i<1024;++i){const double y=Y(t);const int bin=edges->GetYaxis()->FindFixBin(y);if((side?y>boundary:y<boundary)&&bin==(side?edge+1:edge))break;t.pz=::nextafterf(t.pz,side?float(inf):float(-inf));t.energy=std::sqrt(hm*hm+t.px*t.px+t.py*t.py+t.pz*t.pz);Need(i<1023,"representable candidate found near boundary");}const double y=Y(t);edges->Reset();on.FillHe4RapidityBinning(t,.237,4,false,"_signal");Need(edges->GetEntries()==1&&edges->GetBinContent(edges->FindFixBin(.237,y,4))==1,"actual helper follows each boundary from either side");}}
Need(on.He4LabRapidityRegion(on.mHe4RapidityThreshold)==0,"legacy strict equality behavior preserved");Need(edges->GetYaxis()->FindFixBin(on.mHe4RapidityThreshold)>0,"new native threshold remains a defined bin");
FemtoCandidate invalid=Candidate("he4",31,2,-.7,hm,1.);invalid.energy=std::numeric_limits<float>::infinity();const double cnt=edges->GetEntries();on.FillHe4RapidityBinning(invalid,.237,4,false,"_signal");invalid.energy=std::fabs(invalid.pz);on.FillHe4RapidityBinning(invalid,.237,4,false,"_signal");Need(edges->GetEntries()==cnt,"nonfinite and nonphysical rapidity never become flow");edges->Fill(.237,on.mHe4RapidityMin+.5*width,4,16777216.);edges->Fill(.237,on.mHe4RapidityMin+.5*width,4);Need(edges->GetBinContent(edges->FindFixBin(.237,on.mHe4RapidityMin+.5*width,4))==16777217.,"counts above 2^24 retained by double histogram");Metadata(on);
std::cout<<"HE4_YLAB_BINNING_TEST_PASS checks="<<checks<<" calls="<<calls<<std::endl;}
}
#endif
int test_he4_ylab_binning(const char* mainconf){
#ifndef __CINT__
try{He4YLabTest::Run(mainconf);return 0;}catch(const std::exception& e){std::cerr<<"HE4_YLAB_BINNING_TEST_FAIL: "<<e.what()<<std::endl;return 1;}
#else
return 1;
#endif
}
