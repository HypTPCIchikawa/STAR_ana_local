// ROOT5/STAR synthetic regression. Numerical momenta are discriminating fixtures,
// not recommended cuts. The production histogram/configuration comes from caller.
int femto_lambda_pair_diagnostics(const char* mainconf);
#ifndef __CINT__
#include "FemtoCandidate.h"
#include <iosfwd>
#include <string>
#include <vector>
#define private public
#include "StMaker/common/FemtoLambdaLegacy.h"
#undef private
#include "HistManager.h"
#include "YamlParser.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMemFile.h"
#include "TFileMerger.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TVector2.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace PairDiagnosticsTest {
unsigned checks=0;
unsigned pairRapidityChecks=0;
void Need(bool pass,const char* message) { ++checks;if(!pass)throw std::runtime_error(message); }
int Cells(const TH1& h) { return (h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1)*(h.GetDimension()>2?h.GetNbinsZ()+2:1); }
void Equal(const TH1* actual,const TH1* expected) {
  Need(actual&&expected&&Cells(*actual)==Cells(*expected),"matching diagnostic histogram shapes");
  Need(actual->GetEntries()==expected->GetEntries(),"exact histogram entry count");
  for(int i=0;i<Cells(*actual);++i) {
    Need(actual->GetBinContent(i)==expected->GetBinContent(i),"exact bin contents including flows");
    Need(actual->GetBinError(i)==expected->GetBinError(i),"exact Sumw2 including flows");
  }
}
FemtoCandidate Candidate(const char* species,int event,int track,double y,double mass,double pt) {
  FemtoCandidate c;c.speciesKey=species;c.eventIndex=event;
  c.source=std::string(species)=="lambda"?kFemtoCandResonance:kFemtoCandTrack;
  c.trk.trackIndex=track;c.charge=c.source==kFemtoCandTrack?2:0;
  TLorentzVector p4;const double mt=std::sqrt(mass*mass+pt*pt);
  p4.SetPxPyPzE(pt,0.,mt*std::sinh(y),mt*std::cosh(y));c.SetP4(p4);
  c.reso.invMass=1.11596f;c.reso.dau1EventIndex=event;c.reso.dau2EventIndex=event;
  c.reso.dau1Index=0;c.reso.dau2Index=1;return c;
}
std::string Key(int mode,int cent,int region,int side) {
  const char* suffix[]={"","_signal","_leftSB","_rightSB"};
  return std::string("hKstar")+(mode?"ME":"SE")+(cent?"VsCent":"")+
    "_lambda_he4"+suffix[region]+(side?"_YLabAbove":"_YLabBelow");
}
void Rapidity(const char* mainconf) {
  FemtoLambdaLegacy on,off;
  Need(on.Init(mainconf,std::cerr)&&off.Init(mainconf,std::cerr),"rapidity helpers initialize");
  Need(on.mHe4RapiditySplitEnabled&&on.mHe4RapidityThreshold==-.5,"explicit requested He4 threshold");
  off.mHe4RapiditySplitEnabled=false; // Diagnostic bypass, identical physical selection/config.
  const double edge=on.mHe4RapidityThreshold,inf=std::numeric_limits<double>::infinity();
  Need(on.He4LabRapidityRegion(edge)==0,"exact equality enters neither side");
  Need(on.He4LabRapidityRegion(std::nextafter(edge,-inf))==-1&&
       on.He4LabRapidityRegion(std::nextafter(edge,inf))==1,"strict adjacent double boundaries");
  Need(on.He4LabRapidityRegion(inf)==0&&on.He4LabRapidityRegion(-inf)==0&&
       on.He4LabRapidityRegion(std::numeric_limits<double>::quiet_NaN())==0,"nonfinite rapidity rejected only by diagnostic");
  std::map<std::string,TH1*> expected;
  for(int mode=0;mode<2;++mode)for(int cent=0;cent<2;++cent)
    for(int region=0;region<4;++region)for(int side=0;side<2;++side) {
      const std::string key=Key(mode,cent,region,side);TH1* h=on.RootHistograms()->Get(key.c_str());
      Need(h!=0,"all 32 rapidity diagnostics exist");
      expected[key]=static_cast<TH1*>(h->Clone(("expected_"+key).c_str()));
      expected[key]->SetDirectory(0);expected[key]->Reset();
    }
  // Full mass, both signal halves, both sidebands and outside the windows.
  const float masses[]={1.112f,1.119f,1.095f,1.14f,1.20f};
  const int regions[]={1,1,2,3,0},centralities[]={0,4,8};
  const double kstars[]={-.1,.223,1.1};
  const double nuclearMass=on.NuclearP4(TVector3(),3).M();
  for(int direction=0;direction<3;++direction)for(int side=0;side<2;++side)
    for(unsigned m=0;m<sizeof(masses)/sizeof(masses[0]);++m)
      for(unsigned c=0;c<sizeof(centralities)/sizeof(centralities[0]);++c)
        for(unsigned k=0;k<sizeof(kstars)/sizeof(kstars[0]);++k) {
          const bool mixed=direction!=0;const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;
          FemtoCandidate l=Candidate("lambda",le,-1,side?-4.:3.,1.115683,1.1);
          FemtoCandidate n=Candidate("he4",ne,mixed?0:2,edge+(side?.2:-.2),nuclearMass,1.7);
          l.reso.invMass=masses[m];n.y=side?-999.f:999.f; // Must use physical P4, never cached y.
          Need((l.P4()+n.P4()).Rapidity()*(side?1.:-1.)<0.,"pair rapidity deliberately has opposite side");
          Need(on.FillPair(l,n,kstars[k],.2,centralities[c],mixed)&&
               off.FillPair(l,n,kstars[k],.2,centralities[c],mixed),"both ME directions and SE retain pairs");
          for(int region=0;region<4;++region)if(region==0||region==regions[m]) {
            expected[Key(mixed,0,region,side)]->Fill(kstars[k]);
            static_cast<TH2*>(expected[Key(mixed,1,region,side)])->Fill(kstars[k],centralities[c]);
          }
        }
  // An invalid diagnostic P4 does not remove an otherwise admitted physical pair.
  FemtoCandidate l=Candidate("lambda",30,-1,.2,1.115683,.8);
  FemtoCandidate n=Candidate("he4",30,2,-.7,nuclearMass,1.2);n.energy=std::fabs(n.pz);
  Need(on.FillPair(l,n,.37,.2,4,false)&&off.FillPair(l,n,.37,.2,4,false),"invalid rapidity only skips the new diagnostic");
  // Exact equality is representable at pz=0 using an artificial zero threshold.
  on.mHe4RapidityThreshold=0.;n=Candidate("he4",30,2,0.,nuclearMass,1.2);
  Need(n.P4().Rapidity()==0.&&on.FillPair(l,n,.41,.2,4,false)&&off.FillPair(l,n,.41,.2,4,false),
       "exact equality skips diagnostics and keeps inclusive CF");
  on.mHe4RapidityThreshold=edge;
  if(on.IsPairNucleus(2)) {
    n=Candidate("he3",30,2,-.8,on.NuclearP4(TVector3(),2).M(),1.2);
    Need(on.FillPair(l,n,.47,.2,4,false)&&off.FillPair(l,n,.47,.2,4,false),"He3 and mapped He3 control accepted");
  }
  for(std::map<std::string,TH1*>::iterator i=expected.begin();i!=expected.end();++i) {
    Equal(on.RootHistograms()->Get(i->first.c_str()),i->second);delete i->second;
  }
  const char* suffix[]={"","_signal","_leftSB","_rightSB","_signalLow","_signalHigh"};
  for(int mode=0;mode<2;++mode)for(int cent=0;cent<2;++cent)for(int region=0;region<6;++region) {
    const std::string key=std::string("hKstar")+(mode?"ME":"SE")+(cent?"VsCent":"")+"_lambda_he4"+suffix[region];
    Equal(on.RootHistograms()->Get(key.c_str()),off.RootHistograms()->Get(key.c_str()));
  }
  Need(on.RootHistograms()->Get("hKstarSE_lambda_he4_signal")->GetEntries()>0&&
       on.RootHistograms()->Get("hKstarME_lambda_he4_signal")->GetEntries()>0,"nonempty inclusive invariance controls");
}

std::string PairYKey(int mode,int region) {
  const char* suffix[]={"","_signal","_leftSB","_rightSB"};
  return std::string("hKstar")+(mode?"ME":"SE")+"VsPairYVsCent_lambda_he4"+suffix[region];
}
std::string PairTestDirectory() {
  const TString path=TString::Format("%s/tmp/pair-rapidity-test-%d",gSystem->WorkingDirectory(),gSystem->GetPid());
  Need(gSystem->mkdir(path,kTRUE)==0||!gSystem->AccessPathName(path),"pair test artifact directory");
  return path.Data();
}
void PairMetadataMerge(TMemFile& output) {
  const char* names[]={"he4PairRapidityLabNBins","he4PairRapidityLabMin",
                       "he4PairRapidityLabMax","he4RapiditySplitYLabThreshold"};
  const std::string directory=PairTestDirectory();
  for(int worker=0;worker<2;++worker) {
    const TString path=TString::Format("%s/worker%d.root",directory.c_str(),worker);
    Need(gSystem->AccessPathName(path),"fresh metadata worker artifact");
    TFile file(path,"CREATE");Need(!file.IsZombie(),"metadata worker writable");
    for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i) {
      TObject* parameter=output.Get(names[i]);
      file.cd();
      Need(parameter&&parameter->Write()>0,"production parameter copied into synthetic worker");
    }
    TH3* counts=dynamic_cast<TH3*>(output.Get(PairYKey(0,1).c_str()));
    file.cd();
    Need(counts&&counts->GetEntries()>0&&counts->Write()>0,"populated production pair histogram copied");
    file.Close();
  }
  const std::string mergedPath=directory+"/merged.root";
  Need(gSystem->AccessPathName(mergedPath.c_str()),"fresh merged metadata artifact");
  {
    TFileMerger merger;
    Need(merger.OutputFile(mergedPath.c_str(),"RECREATE"),"merge output configured");
    Need(merger.AddFile((directory+"/worker0.root").c_str())&&
         merger.AddFile((directory+"/worker1.root").c_str())&&merger.Merge(),"two-worker ROOT merge succeeds");
  }
  TFile merged(mergedPath.c_str(),"READ");Need(!merged.IsZombie(),"merged metadata reopens");
  TParameter<int>* bins=dynamic_cast<TParameter<int>*>(merged.Get(names[0]));
  TParameter<double>* low=dynamic_cast<TParameter<double>*>(merged.Get(names[1]));
  TParameter<double>* high=dynamic_cast<TParameter<double>*>(merged.Get(names[2]));
  TParameter<double>* cut=dynamic_cast<TParameter<double>*>(merged.Get(names[3]));
  Need(bins&&bins->GetVal()==10&&bins->IsConst(),"hadd preserves pair-Y bin count rather than summing it");
  Need(low&&low->GetVal()==-1.5&&low->IsConst()&&high&&high->GetVal()==.5&&high->IsConst(),
       "hadd preserves both pair-Y range constants");
  Need(cut&&cut->GetVal()==-.5&&cut->IsConst(),"hadd preserves the existing He4 rapidity threshold");
  TH3* before=dynamic_cast<TH3*>(output.Get(PairYKey(0,1).c_str()));
  TH3* after=dynamic_cast<TH3*>(merged.Get(PairYKey(0,1).c_str()));
  Need(before&&after&&after->GetEntries()==2*before->GetEntries(),"hadd doubles raw pair counts");
  for(int cell=0;cell<Cells(*before);++cell)
    Need(after->GetBinContent(cell)==2*before->GetBinContent(cell),"hadd preserves double-precision 3D counts and flows");
}
double IndependentPairY(const FemtoCandidate& a,const FemtoCandidate& b) {
  const double energy=double(a.energy)+double(b.energy);
  const double pz=double(a.pz)+double(b.pz);
  if(!std::isfinite(energy)||!std::isfinite(pz)||!(energy>std::fabs(pz)))
    return std::numeric_limits<double>::quiet_NaN();
  return .5*std::log((energy+pz)/(energy-pz));
}
double Variance(const TH1& h,int bin) {
  return h.GetSumw2N()?h.GetSumw2()->At(bin):h.GetBinContent(bin);
}
void PairProjectionClosure(FemtoLambdaLegacy& helper) {
  const char* suffix[]={"","_signal","_leftSB","_rightSB"};
  for(int mode=0;mode<2;++mode)for(int region=0;region<4;++region) {
    TH3* source=dynamic_cast<TH3*>(helper.RootHistograms()->Get(PairYKey(mode,region).c_str()));
    const std::string old=std::string("hKstar")+(mode?"ME":"SE")+"_lambda_he4"+suffix[region];
    TH1* inclusive=helper.RootHistograms()->Get(old.c_str());
    Need(source&&inclusive,"inclusive and pair-rapidity histograms exist");
    double normal=0.,under=0.,over=0.;
    for(int x=0;x<=source->GetNbinsX()+1;++x) {
      double all=0.,variance=0.,regular=0.,adjacent=0.;
      for(int z=0;z<=source->GetNbinsZ()+1;++z)
        for(int y=0;y<=source->GetNbinsY()+1;++y) {
          const int cell=source->GetBin(x,y,z);const double value=source->GetBinContent(cell);
          all+=value;variance+=Variance(*source,cell);
          if(y==0)under+=value;
          else if(y==source->GetNbinsY()+1)over+=value;
          else {normal+=value;regular+=value;}
        }
      // Merge neighboring rapidity intervals into five intervals. Counts must
      // agree with the ten-bin projection without using the producer's helper.
      for(int group=0;group<5;++group)
        for(int y=2*group+1;y<=2*group+2;++y)
          for(int z=0;z<=source->GetNbinsZ()+1;++z)
            adjacent+=source->GetBinContent(source->GetBin(x,y,z));
      Need(all==inclusive->GetBinContent(x),"all pair-Y bins plus flows reproduce inclusive kstar");
      Need(variance==Variance(*inclusive,x),"pair-Y projection preserves counting variance");
      Need(adjacent==regular,"adjacent pair-Y merging preserves regular-bin counts");
    }
    Need(normal>0.&&under>0.&&over>0.,"every mass sample exercises regular pair-Y and both flows");
    for(int y=1;y<=source->GetNbinsY();++y)
      Need(source->Integral(0,source->GetNbinsX()+1,y,y,0,source->GetNbinsZ()+1)>0.,
           "all ten pair-rapidity intervals are populated");
  }
}
void PairRapidity(const char* mainconf) {
  const unsigned before=checks;
  FemtoLambdaLegacy on,off;
  Need(on.Init(mainconf,std::cerr)&&off.Init(mainconf,std::cerr),"pair-rapidity helpers initialize");
  if(!on.mPairRapiditySplitEnabled) {
    std::cout<<"PAIR_RAPIDITY_TEST_SKIPPED disabled profile"<<std::endl;return;
  }
  Need(on.mPairRapidityNBins==10&&on.mPairRapidityMin==-1.5&&on.mPairRapidityMax==.5,
       "explicit lab pair-rapidity ten-bin configuration");
  off.mPairRapiditySplitEnabled=false;
  std::map<std::string,TH3*> expected;
  for(int mode=0;mode<2;++mode)for(int region=0;region<4;++region) {
    const std::string key=PairYKey(mode,region);
    TH3D* source=dynamic_cast<TH3D*>(on.RootHistograms()->Get(key.c_str()));
    Need(source!=0,"all eight pair-rapidity histograms have double precision");
    Need(source->GetNbinsY()==10&&source->GetYaxis()->GetXmin()==-1.5&&
         source->GetYaxis()->GetXmax()==.5&&source->GetNbinsZ()==9,"pair-Y and centrality axes");
    expected[key]=static_cast<TH3*>(source->Clone(("expected_"+key).c_str()));
    expected[key]->SetDirectory(0);expected[key]->Reset();
  }
  const float masses[]={1.112f,1.119f,1.095f,1.14f,1.20f};
  const int regions[]={1,1,2,3,0},centralities[]={0,4,8};
  const double kstars[]={-.1,.223,1.1};
  const double mass=on.NuclearP4(TVector3(),3).M();
  const double width=(on.mPairRapidityMax-on.mPairRapidityMin)/on.mPairRapidityNBins;
  for(int direction=0;direction<3;++direction)for(int interval=-1;interval<=10;++interval)
    for(unsigned m=0;m<sizeof(masses)/sizeof(masses[0]);++m)
      for(unsigned c=0;c<sizeof(centralities)/sizeof(centralities[0]);++c)
        for(unsigned k=0;k<sizeof(kstars)/sizeof(kstars[0]);++k) {
          const bool mixed=direction!=0;
          const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;
          const double y=on.mPairRapidityMin+(interval+.5)*width;
          FemtoCandidate l=Candidate("lambda",le,-1,y,1.115683,1.1);
          FemtoCandidate n=Candidate("he4",ne,mixed?0:2,y,mass,1.7);
          l.reso.invMass=masses[m];l.y=999.f;n.y=-999.f;
          const double pairY=IndependentPairY(l,n);
          Need(std::isfinite(pairY),"independent sum-P4 rapidity is finite");
          Need(on.FillPair(l,n,kstars[k],.2,centralities[c],mixed)&&
               off.FillPair(l,n,kstars[k],.2,centralities[c],mixed),"SE and both ME directions retain pair");
          for(int region=0;region<4;++region)if(region==0||region==regions[m])
            expected[PairYKey(mixed,region)]->Fill(kstars[k],pairY,centralities[c]);
        }
  // This pair lies in different bins for pair Y, nucleus Y and Lambda Y.
  // Stale cached y and a CM boost would also classify it incorrectly.
  for(int direction=0;direction<3;++direction) {
    const bool mixed=direction!=0;const int le=direction==2?20:10,ne=mixed?(direction==2?10:20):10;
    FemtoCandidate l=Candidate("lambda",le,-1,-1.8,1.115683,1.1);
    FemtoCandidate n=Candidate("he4",ne,mixed?0:2,.1,mass,1.7);
    l.y=999.f;n.y=-999.f;const double pairY=IndependentPairY(l,n);
    TH3* axis=expected[PairYKey(mixed,0)];
    const int pairBin=axis->GetYaxis()->FindFixBin(pairY);
    Need(pairBin!=axis->GetYaxis()->FindFixBin(n.P4().Rapidity())&&
         pairBin!=axis->GetYaxis()->FindFixBin(l.P4().Rapidity()),"pair rapidity differs from either particle");
    Need(on.FillPair(l,n,.317,.2,4,mixed)&&off.FillPair(l,n,.317,.2,4,mixed),"discriminating pair accepted");
    expected[PairYKey(mixed,0)]->Fill(.317,pairY,4);
    expected[PairYKey(mixed,1)]->Fill(.317,pairY,4);
  }
  for(std::map<std::string,TH3*>::iterator i=expected.begin();i!=expected.end();++i)
    Equal(on.RootHistograms()->Get(i->first.c_str()),i->second);
  PairProjectionClosure(on);
  // Invalid pair energy must not silently lose an otherwise accepted physical
  // pair. Only the additional rapidity diagnostic is unavailable.
  for(int mode=0;mode<2;++mode) {
    FemtoCandidate l=Candidate("lambda",10,-1,-.3,1.115683,.9);
    FemtoCandidate n=Candidate("he4",mode?20:10,2,-.6,mass,1.2);
    n.energy=std::numeric_limits<float>::quiet_NaN();
    Need(on.FillPair(l,n,.419,.2,4,mode)&&off.FillPair(l,n,.419,.2,4,mode),"nonfinite pair-Y skips diagnostic only");
  }
  if(on.IsPairNucleus(2)) {
    FemtoCandidate l=Candidate("lambda",30,-1,-.2,1.115683,.9);
    FemtoCandidate n=Candidate("he3",30,2,-.7,on.NuclearP4(TVector3(),2).M(),1.2);
    Need(on.FillPair(l,n,.471,.2,4,false)&&off.FillPair(l,n,.471,.2,4,false),"physical and mass-mapped He3 controls accepted");
  }
  for(std::map<std::string,TH3*>::iterator i=expected.begin();i!=expected.end();++i) {
    Equal(on.RootHistograms()->Get(i->first.c_str()),i->second);delete i->second;
  }
  for(size_t i=0;i<on.mRootKeys.size();++i) {
    const std::string& key=on.mRootKeys[i];
    if(key.find("hKstar")==0&&key.find("VsPairYVsCent")==std::string::npos)
      Equal(on.RootHistograms()->Get(key.c_str()),off.RootHistograms()->Get(key.c_str()));
  }
  {
    TMemFile output("pairRapidityMetadata.root","RECREATE");
    Need(on.Write(&output),"pair diagnostic output metadata is writable");
    TParameter<int>* bins=dynamic_cast<TParameter<int>*>(output.Get("he4PairRapidityLabNBins"));
    TParameter<double>* low=dynamic_cast<TParameter<double>*>(output.Get("he4PairRapidityLabMin"));
    TParameter<double>* high=dynamic_cast<TParameter<double>*>(output.Get("he4PairRapidityLabMax"));
    TNamed* definition=dynamic_cast<TNamed*>(output.Get("FemtoHe4PairRapiditySplitDefinition"));
    Need(bins&&bins->GetVal()==10&&low&&low->GetVal()==-1.5&&high&&high->GetVal()==.5,
         "typed pair-Y axis metadata exactly matches configuration");
    Need(definition&&std::string(definition->GetTitle()).find("lab")!=std::string::npos,
         "pair-Y definition records the laboratory frame");
    PairMetadataMerge(output);
  }
  // Exercise every axis edge from each side using continuous double-P4 input.
  // Stored Float candidates have also been tested independently above.
  TH3* edgeHistogram=static_cast<TH3*>(on.RootHistograms()->Get(PairYKey(0,1).c_str()));
  Need(edgeHistogram->GetYaxis()->FindFixBin(on.mPairRapidityMin)==1&&
       edgeHistogram->GetYaxis()->FindFixBin(on.mPairRapidityMax)==11,"ROOT half-open global interval");
  for(int edge=0;edge<=10;++edge)for(int side=0;side<2;++side) {
    const double boundary=on.mPairRapidityMin+edge*width;
    const int yBin=side?edge+1:edge;
    double pz=8.*std::tanh(boundary),y=0.;unsigned steps=0;
    do {
      pz=std::nextafter(pz,side?std::numeric_limits<double>::infinity():-std::numeric_limits<double>::infinity());
      y=.5*std::log((8.+pz)/(8.-pz));++steps;
    // ROOT's subtraction and scaling can round a nominal-edge neighbor back
    // across that edge. Select the first nearby representable value whose
    // independent rapidity also belongs to the intended native ROOT cell.
    }while(((side?!(y>boundary):!(y<boundary))||
            edgeHistogram->GetYaxis()->FindFixBin(y)!=yBin)&&steps<1024);
    Need(steps<1024,"construct representable pair just across rapidity edge");
    TLorentzVector pair(0.,0.,pz,8.);edgeHistogram->Reset();
    on.FillHe4PairRapidity(pair,.237,4,false,"_signal");
    const int cell=edgeHistogram->GetBin(edgeHistogram->GetXaxis()->FindFixBin(.237),yBin,5);
    Need(edgeHistogram->GetEntries()==1.&&edgeHistogram->GetBinContent(cell)==1.,
         "all ten pair-Y intervals and under/overflow use native edge convention");
  }
  const TLorentzVector invalid[]={TLorentzVector(0.,0.,1.,1.),TLorentzVector(0.,0.,2.,1.),
    TLorentzVector(0.,0.,0.,-1.),TLorentzVector(0.,0.,0.,std::numeric_limits<double>::infinity())};
  edgeHistogram->Reset();
  for(unsigned i=0;i<sizeof(invalid)/sizeof(invalid[0]);++i)
    on.FillHe4PairRapidity(invalid[i],.237,4,false,"_signal");
  Need(edgeHistogram->GetEntries()==0.,"unphysical or nonfinite pair rapidity never becomes a flow entry");
  // An individual rapidity cell may accumulate more than 2^24 unit pairs.
  // Float histograms lose the next count; the declared TH3D must retain it.
  edgeHistogram->Fill(.237,0.,4.,16777216.);edgeHistogram->Fill(.237,0.,4.);
  Need(edgeHistogram->GetBinContent(edgeHistogram->FindFixBin(.237,0.,4.))==16777217.,
       "pair-Y double histogram preserves a count beyond Float integer precision");
  pairRapidityChecks+=checks-before;
}

std::string ReadText(const std::string& path) {
  std::ifstream file(path.c_str());Need(bool(file),"fixture source readable");
  std::ostringstream text;text<<file.rdbuf();return text.str();
}
void WriteText(const std::string& path,const std::string& text) {
  Need(gSystem->AccessPathName(path.c_str()),"fixture never replaces existing data");
  std::ofstream file(path.c_str());file<<text;Need(bool(file),"fixture write succeeds");
}
std::string ReplacePairBlock(const std::string& input,const std::string& block) {
  size_t start=input.find("\npairRapiditySplit:");
  if(start!=std::string::npos)++start;
  else if(input.find("pairRapiditySplit:")==0)start=0;
  Need(start!=std::string::npos,"valid maker has explicit pair-Y block");
  size_t end=input.find('\n',start);
  while(end!=std::string::npos&&end+1<input.size()) {
    const size_t next=end+1;
    if(input[next]!=' '&&input[next]!='\t'&&input[next]!='\r'&&input[next]!='\n')break;
    end=input.find('\n',next);
  }
  if(end==std::string::npos)end=input.size();else ++end;
  return input.substr(0,start)+block+input.substr(end);
}
void PairConfiguration(const char* mainconf) {
  const unsigned before=checks;
  std::string path=mainconf;
  if(path.empty()||path[0]!='/')path=std::string(gSystem->WorkingDirectory())+"/"+path;
  const size_t mark=path.rfind("/config/");Need(mark!=std::string::npos,"explicit config-root entry");
  const std::string configRoot=path.substr(0,mark)+"/config/";
  std::map<std::string,std::string> main;
  Need(YamlParser::ParseFile(path.c_str(),main)&&main.count("maker")&&main.count("femtoHist"),"scalar mainconf parsed");
  for(std::map<std::string,std::string>::iterator i=main.begin();i!=main.end();++i)
    if(!i->second.empty()&&i->second[0]!='/')i->second=configRoot+i->second;
  const std::string maker=ReadText(main["maker"]);
  if(maker.find("pairRapiditySplit:")==std::string::npos)return;
  const std::string directory=PairTestDirectory()+"/config";
  Need(gSystem->mkdir((directory+"/mainconf").c_str(),kTRUE)==0||
       !gSystem->AccessPathName((directory+"/mainconf").c_str()),"isolated config-root fixture");
  Need(gSystem->Symlink((configRoot+"calibration").c_str(),(directory+"/calibration").c_str())==0,
       "nested detector calibration remains anchored to the original read-only directory");
  const char* invalid[]={
    "pairRapiditySplit:\n  enabled: true\n  frame: cm\n  nBins: 10\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 0\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: -2\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10.5\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  min: 0.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  min: 1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  min: .nan\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  min: -1.5\n  max: .inf\n",
    "pairRapiditySplit:\n  frame: lab\n  nBins: 10\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  nBins: 10\n  min: -1.5\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  max: 0.5\n",
    "pairRapiditySplit:\n  enabled: true\n  frame: lab\n  nBins: 10\n  min: -1.5\n  max: 0.5\n  extra: 7\n",
    "pairRapiditySplit:\n  enabled: false\n  frame: lab\n  nBins: 10\n  min: -1.5\n  max: 0.5\n"
  };
  for(unsigned test=0;test<=sizeof(invalid)/sizeof(invalid[0]);++test) {
    const TString stem=TString::Format("%s/case%u",directory.c_str(),test);
    const std::string makerPath=std::string(stem.Data())+"_maker.yaml";
    WriteText(makerPath,test?ReplacePairBlock(maker,invalid[test-1]):maker);
    std::ostringstream text;
    for(std::map<std::string,std::string>::const_iterator i=main.begin();i!=main.end();++i)
      text<<i->first<<": "<<(i->first=="maker"?makerPath:i->second)<<"\n";
    const std::string mainPath=std::string(stem.Data())+"_main.yaml";
    WriteText(mainPath,text.str());FemtoLambdaLegacy helper;std::ostringstream diagnostics;
    const bool valid=helper.Init(mainPath.c_str(),diagnostics);
    Need(test?!valid:valid,test?"invalid pair-Y config rejected":"valid relocated config control succeeds");
    if(test)Need(!diagnostics.str().empty(),"invalid pair-Y config reports its failure");
  }
  const std::string histogram=ReadText(main["femtoHist"]);
  const std::string key=PairYKey(0,0);
  const size_t keyAt=histogram.find(key+":");Need(keyAt!=std::string::npos,"pair-Y contract entry in histogram profile");
  for(int test=0;test<2;++test) {
    std::string changed=histogram;
    if(test==0)changed.replace(keyAt,key.size(),key+"_missing");
    else {
      const size_t typeAt=changed.find("type: TH3D",keyAt);
      Need(typeAt!=std::string::npos,"explicit TH3D definition in diagnostic profile");
      changed.replace(typeAt,10,"type: TH3F");
    }
    const TString stem=TString::Format("%s/histcase%d",directory.c_str(),test);
    const std::string histPath=std::string(stem.Data())+"_hist.yaml";
    WriteText(histPath,changed);std::ostringstream text;
    for(std::map<std::string,std::string>::const_iterator i=main.begin();i!=main.end();++i)
      text<<i->first<<": "<<(i->first=="femtoHist"?histPath:i->second)<<"\n";
    const std::string mainPath=std::string(stem.Data())+"_main.yaml";
    WriteText(mainPath,text.str());FemtoLambdaLegacy helper;std::ostringstream diagnostics;
    Need(!helper.Init(mainPath.c_str(),diagnostics),"missing or reduced-precision pair histogram rejected");
    Need(!diagnostics.str().empty(),"histogram contract rejection is reported");
  }
  pairRapidityChecks+=checks-before;
}
struct PicoFixture {
  TClonesArray events,tracks;TClonesArray* arrays[StPicoArrays::NAllPicoArrays];StPicoDst dst;
  PicoFixture():events("StPicoEvent"),tracks("StPicoTrack") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;
    new(events[0])StPicoEvent();StPicoDst::set(arrays);
  }
  ~PicoFixture(){StPicoDst::unset();}
  StPicoTrack* Add(int index,double pt,double phi,double pz,int charge) {
    StPicoTrack* t=new(tracks[index])StPicoTrack();t->setId(100+index);t->setNHitsFit(25*charge);
    t->setGlobalMomentum(pt*std::cos(phi),pt*std::sin(phi),pz);t->setOrigin(0.,0.,0.);return t;
  }
};
double PhiAtRadius(const StPicoTrack& t,double kg,double radius) {
  // Independent origin-at-beamline circle geometry; no STAR helix calls.
  return t.gMom().Phi()-std::asin(.299792458*(.1*kg)*t.charge()*radius/(2*t.gMom().Pt()));
}
void Merging(const char* mainconf) {
  FemtoLambdaLegacy helper;Need(helper.Init(mainconf,std::cerr),"merging helper initializes");
  Need(helper.mMergingIncludePion,"test profile explicitly includes pion diagnostics");
  const double r=helper.mRadiusMeters,z=helper.mMergingTpcHalfLengthCm,tol=helper.mMergingRadiusToleranceCm;
  PicoFixture f;StPicoTrack* proton=f.Add(0,.8,.22,-.12,1);
  StPicoTrack* pion=f.Add(1,.6,-.15,.08,-1);StPicoTrack* nucleus=f.Add(2,1.2,.1,-.3,1);
  Need(proton->pMom().Mag2()==0&&pion->pMom().Mag2()==0,"both Lambda daughter primary fits are absent");
  FemtoCandidate l=Candidate("lambda",0,-1,.1,1.115683,1.);
  FemtoCandidate n=Candidate("he4",0,2,-.7,helper.NuclearP4(TVector3(),3).M(),1.2);
  n.SetP4(helper.NuclearP4(nucleus->gMom(),3)); // q=2, physical momentum=2*rigidity.
  const char* daughters[]={"proton","pion"};StPicoTrack* track[]={proton,pion};
  TH2* expected[2];
  for(int i=0;i<2;++i) {
    TH2* h=static_cast<TH2*>(helper.NuclearHistograms()->Get((std::string("hDphiDeta_")+daughters[i]+"_4He").c_str()));
    Need(h!=0,"both daughter merger maps exist");expected[i]=static_cast<TH2*>(h->Clone((std::string("expected_merging_")+daughters[i]).c_str()));
    expected[i]->SetDirectory(0);expected[i]->Reset();
  }
  const double fields[]={5.,-5.,0.};
  for(unsigned b=0;b<sizeof(fields)/sizeof(fields[0]);++b) {
    f.dst.event()->setBField(fields[b]);TVector3 np;
    Need(FemtoLambdaLegacy::MergingTrackPositionAtRadius(nucleus,fields[b],r,z,tol,np),"reachable nuclear helix");
    Need(std::fabs(TVector2::Phi_mpi_pi(np.Phi()-PhiAtRadius(*nucleus,fields[b],r)))<1.e-7,"raw He rigidity and STAR curvature sign");
    for(int i=0;i<2;++i) {
      TVector3 p;Need(FemtoLambdaLegacy::MergingTrackPositionAtRadius(track[i],fields[b],r,z,tol,p),"secondary global helix reaches cylinder");
      const double phi=PhiAtRadius(*track[i],fields[b],r);
      Need(std::fabs(TVector2::Phi_mpi_pi(p.Phi()-phi))<1.e-7,"both charge signs agree with independent circle geometry");
      const double dphi=TVector2::Phi_mpi_pi(phi-PhiAtRadius(*nucleus,fields[b],r));
      expected[i]->Fill(dphi,track[i]->gMom().Eta()-nucleus->gMom().Eta());
      if(b==0) {
        const double oldDphi=TVector2::Phi_mpi_pi(2*track[i]->gMom().Phi()-phi-
          (2*nucleus->gMom().Phi()-PhiAtRadius(*nucleus,fields[b],r)));
        Need(expected[i]->GetXaxis()->FindBin(dphi)!=expected[i]->GetXaxis()->FindBin(oldDphi),"fixture distinguishes historical plus-asin bug");
      }
    }
    Need(helper.FillPair(l,n,.2,.3,4,false,&f.dst),"corrected merger QA retains signal SE");
  }
  // Unreachable low-pT helix and a global-z crossing outside the TPC must not
  // become clamped asin entries, and must not impose any physical pair cut.
  f.dst.event()->setBField(5.);TVector3 p;
  nucleus->setGlobalMomentum(.02,0.,0.);
  Need(!FemtoLambdaLegacy::MergingTrackPositionAtRadius(nucleus,5.,r,z,tol,p),"unreachable helix rejected");
  Need(helper.FillPair(l,n,.2,.3,4,false,&f.dst),"unreachable diagnostic keeps physical SE");
  nucleus->setGlobalMomentum(1.2,0.,0.);nucleus->setOrigin(0.,0.,z+10.);
  Need(!FemtoLambdaLegacy::MergingTrackPositionAtRadius(nucleus,5.,r,z,tol,p),"global TPC z reach enforced");
  Need(helper.FillPair(l,n,.2,.3,4,false,&f.dst),"out-of-TPC diagnostic keeps physical SE");
  for(int i=0;i<2;++i) {
    const std::string key=std::string("hDphiDeta_")+daughters[i]+"_4He";
    Equal(helper.NuclearHistograms()->Get(key.c_str()),expected[i]);
    Equal(helper.NuclearHistograms()->Get((key+"_EventField").c_str()),expected[i]);delete expected[i];
  }
  Need(helper.RootHistograms()->Get("hKstarSE_lambda_he4_signal")->GetEntries()==5.,"geometry diagnostics never drop accepted physical pairs");
}
}
#endif
int femto_lambda_pair_diagnostics(const char* mainconf) {
#ifndef __CINT__
  try { PairDiagnosticsTest::checks=0;PairDiagnosticsTest::pairRapidityChecks=0;
    PairDiagnosticsTest::Rapidity(mainconf);PairDiagnosticsTest::Merging(mainconf);
    PairDiagnosticsTest::PairRapidity(mainconf);PairDiagnosticsTest::PairConfiguration(mainconf);
    std::cout<<"FEMTO_LAMBDA_PAIR_DIAGNOSTICS_PASS checks="<<PairDiagnosticsTest::checks
             <<" pairRapidityChecks="<<PairDiagnosticsTest::pairRapidityChecks<<std::endl;return 0;
  } catch(const std::exception& e) { std::cerr<<"FEMTO_LAMBDA_PAIR_DIAGNOSTICS_FAIL: "<<e.what()<<std::endl;return 1; }
#else
  return 1;
#endif
}
