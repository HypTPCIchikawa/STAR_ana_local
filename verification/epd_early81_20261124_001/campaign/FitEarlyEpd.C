#include "TChain.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TNamed.h"
#include "TSystem.h"
#include "/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/include/FemtoEpdEventPlane.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace earlyepd {
const int nQ=36, terms=20, abX=16;
const double pi=TMath::Pi();

struct Estimate {
  Long64_t count;
  double sum[nQ], square[nQ];
  Estimate():count(0) {
    std::fill(sum,sum+nQ,0.); std::fill(square,square+nQ,0.);
  }
  void Add(const double* values) {
    ++count;
    for(int i=0;i<nQ;++i) {
      if(!std::isfinite(values[i])) throw std::runtime_error("nonfinite raw Q");
      sum[i]+=values[i]; square[i]+=values[i]*values[i];
    }
  }
  double Mean(int axis) const { return sum[axis]/count; }
  double Error(int axis) const {
    const double mean=Mean(axis);
    return std::sqrt(std::max(0.,square[axis]/count-mean*mean)/count);
  }
};

struct Fourier {
  Long64_t count;
  double sine[terms], cosine[terms];
  Fourier():count(0) {std::fill(sine,sine+terms,0.); std::fill(cosine,cosine+terms,0.);}
  void Add(double angle) {
    ++count;
    for(int i=0;i<terms;++i) {
      sine[i]+=std::sin((i+1)*angle);cosine[i]+=std::cos((i+1)*angle);
    }
  }
  double Correct(double angle) const {
    double corrected=angle;
    for(int i=0;i<terms;++i) corrected+=2./(i+1)*
      (-sine[i]/count*std::cos((i+1)*angle)+cosine[i]/count*std::sin((i+1)*angle));
    if(corrected < -pi) corrected+=2*pi;
    if(corrected > pi) corrected-=2*pi;
    return corrected;
  }
};

struct Cell {
  Estimate raw, residual;
  Fourier recentered, shifted;
  Long64_t validationRaw, invalid;
  Long64_t bins[8];
  int model;
  Cell():validationRaw(0),invalid(0),model(0) {std::fill(bins,bins+8,0);}
};
typedef std::pair<int,int> Key;
typedef std::map<Key,Cell> Cells;
struct Entry {
  Int_t run, cent, fold;
  Double_t q[nQ];
};
int Day(int run) {return (run/1000)%1000;}
void Require(bool condition,const char* message) {
  if(!condition) throw std::runtime_error(message);
}
bool Center(const Entry& row,const Estimate& mean,double* centered) {
  for(int a=0;a<nQ;++a) centered[a]=row.q[a]-mean.Mean(a);
  for(int a=0;a<nQ;++a) if(!std::isfinite(centered[a]) || centered[a]==0.) return false;
  return true;
}
void WriteProfiles(const std::string& out,const std::vector<int>& runs,
                   const Cells& perRun,const Cells& perDay) {
  TFile rec((out+"/candidate_recenterpar.root").c_str(),"CREATE");
  Require(!rec.IsZombie(),"could not create recenter ROOT");
  for(int i=0;i<nQ;++i) {
    const int j=(i<8?i:i-8);
    const int h=j/(i<8?4:14), group=(j/2)%(i<8?2:7);
    const char axis=i%2?'y':'x';
    const TString name=i<8?Form("TPCQ%c_%c_recen%s",axis,'A'+group,h?"_psi2":""):
      Form("EPDQ%c_ring%d_recen%s",axis,group,h?"_psi2":"");
    TProfile2D profile(name,name,9,-.5,8.5,runs.size(),-.5,runs.size()-.5);
    // Store the measured mean and the count that contributed to it.
    // Profile errors are not used. Held-out errors are recorded in the QA.
    for(size_t r=0;r<runs.size();++r) for(int cent=0;cent<9;++cent) {
      const Cell& cell=perRun.at(Key(runs[r],cent));
      const Cell& model=cell.model==1?cell:perDay.at(Key(Day(runs[r]),cent));
      profile.Fill(cent,r,model.raw.Mean(i),static_cast<double>(model.raw.count));
    }
    Require(profile.Write()>0,"recenter profile write failed");
  }
  rec.Close();
  Require(!rec.TestBit(TFile::kWriteError),"recenter close failed");
  TFile shift((out+"/candidate_shiftpar.root").c_str(),"CREATE");
  Require(!shift.IsZombie(),"could not create shift ROOT");
  for(int component=0;component<2;++component) {
    const char* name=component?"EPDshiftpar_cos_ring4":"EPDshiftpar_sin_ring4";
    TProfile3D profile(name,name,9,-.5,8.5,terms,-.5,terms-.5,runs.size(),-.5,runs.size()-.5);
    for(size_t r=0;r<runs.size();++r) for(int cent=0;cent<9;++cent) {
      const Cell& cell=perRun.at(Key(runs[r],cent));
      const Cell& model=cell.model==1?cell:perDay.at(Key(Day(runs[r]),cent));
      for(int i=0;i<terms;++i) {
        const double total=component?model.recentered.cosine[i]:model.recentered.sine[i];
        profile.Fill(cent,i,r,total/model.recentered.count,model.recentered.count);
      }
    }
    Require(profile.Write()>0,"shift profile write failed");
  }
  shift.Close();
  Require(!shift.TestBit(TFile::kWriteError),"shift close failed");
}

bool Pass(const Cell& cell,Long64_t minValid,double qTol,double fourierTol,
          double countAllowance,double binTol) {
  if(cell.shifted.count<minValid || cell.invalid*100>cell.validationRaw) return false;
  for(int a=0;a<nQ;++a)
    if(std::fabs(cell.residual.Mean(a)) > std::max(qTol,countAllowance*cell.residual.Error(a))) return false;
  const double n=cell.shifted.count;
  const double fLimit=std::max(fourierTol,countAllowance/std::sqrt(2*n));
  for(int h=0;h<terms;++h)
    if(std::max(std::fabs(cell.shifted.sine[h]),std::fabs(cell.shifted.cosine[h]))/n>fLimit) return false;
  const double p=1./8.;
  const double bLimit=std::max(binTol,countAllowance*std::sqrt(p*(1-p)/n));
  for(int b=0;b<8;++b) if(std::fabs(cell.bins[b]/n-p)>bLimit) return false;
  return true;
}
} // namespace earlyepd

void FitEarlyEpd(const char* list,const char* out,int minFit,int minValid,
                 double qTol,double fourierTol,double countAllowance,double binTol) {
  using namespace earlyepd;
  try {
    Require(minFit>0 && minValid>0 && qTol>0 && fourierTol>0 && countAllowance>0 && binTol>0,
            "invalid validation thresholds");
    if(gSystem->AccessPathName(out) || !gSystem->AccessPathName(Form("%s/qa.tsv",out)))
      throw std::runtime_error("fit output must be a fresh existing directory");
    TChain chain("epdRaw");
    std::ifstream paths(list);
    std::string path;
    while(std::getline(paths,path)) Require(chain.Add(path.c_str())==1,"bad cache path");
    const Long64_t entries=chain.GetEntries();
    Require(bool(paths.eof()) && entries>0,"empty cache");
    Entry row;
    Require(chain.SetBranchAddress("run",&row.run)>=0 &&
            chain.SetBranchAddress("calibrationCent",&row.cent)>=0 &&
            chain.SetBranchAddress("fold",&row.fold)>=0 &&
            chain.SetBranchAddress("q",row.q)>=0,"cache schema mismatch");
    chain.SetBranchStatus("*",0);
    for(const char* name : {"run","calibrationCent","fold","q"}) chain.SetBranchStatus(name,1);
    Cells perRun,perDay,dayQA;
    std::set<int> runs;
    for(Long64_t i=0;i<entries;++i) {
      Require(chain.GetEntry(i)>0,"cache I/O");
      Require(Day(row.run)>=121 && Day(row.run)<=125 && row.cent>=0 && row.cent<9 &&
              (row.fold==0 || row.fold==1),"unexpected cache data");
      runs.insert(row.run);
      if(row.fold==0) {
        perRun[Key(row.run,row.cent)].raw.Add(row.q);
        perDay[Key(Day(row.run),row.cent)].raw.Add(row.q);
      }
    }
    for(int run:runs) for(int cent=0;cent<9;++cent) {
      Cell& cell=perRun[Key(run,cent)];
      cell.model=cell.raw.count>=minFit?1:
                 (perDay[Key(Day(run),cent)].raw.count>=minFit?2:0);
    }
    for(Long64_t i=0;i<entries;++i) {
      Require(chain.GetEntry(i)>0,"cache I/O");
      if(row.fold) continue;
      Cell& run=perRun.at(Key(row.run,row.cent));
      Cell& day=perDay.at(Key(Day(row.run),row.cent));
      double centered[nQ];
      if(run.model==1 && Center(row,run.raw,centered))
        run.recentered.Add(std::atan2(centered[abX+1],centered[abX]));
      if(day.raw.count>=minFit && Center(row,day.raw,centered))
        day.recentered.Add(std::atan2(centered[abX+1],centered[abX]));
    }
    for(Long64_t i=0;i<entries;++i) {
      Require(chain.GetEntry(i)>0,"cache I/O");
      if(!row.fold) continue;
      Cell& run=perRun.at(Key(row.run,row.cent));
      const Cell& day=perDay.at(Key(Day(row.run),row.cent));
      Cell* targets[]={&run,&dayQA[Key(Day(row.run),row.cent)]};
      const Cell* models[]={run.model==1?&run:&day,&day};
      for(int which=0;which<2;++which) {
        Cell& target=*targets[which];
        const Cell& model=*models[which];
        ++target.validationRaw;
        double centered[nQ];
        if(model.recentered.count<minFit || !Center(row,model.raw,centered)) {
          ++target.invalid;continue;
        }
        const double flat=model.recentered.Correct(std::atan2(centered[abX+1],centered[abX]));
        const int bin=femto_epd::XialeiBin(flat);
        if(bin<0) {++target.invalid;continue;}
        target.residual.Add(centered);
        target.shifted.Add(flat);
        ++target.bins[bin];
      }
    }
    const std::string dir(out);
    std::ofstream qa((dir+"/qa.tsv").c_str());
    qa<<"run\tcalCent\tmodel\tfitRun\tfitDay\tfitShift\tvalidRaw\tvalidFinal\tinvalid\tmaxAbsMeanQ\tmaxAbsFourier\tmaxBinDev\tpass";
    for(int b=0;b<8;++b) qa<<"\tbin"<<b;
    qa<<'\n'<<std::setprecision(10);
    std::vector<int> accepted;
    for(int run:runs) {
      bool good=true, evaluated=false;
      for(int cent=0;cent<9;++cent) {
        Cell& cell=perRun.at(Key(run,cent));
        const Key dayKey(Day(run),cent);
        const Cell& day=perDay.at(dayKey);
        const Cell& model=cell.model==1?cell:day;
        bool ok=cell.model!=0 && model.recentered.count>=minFit;
        if(cell.shifted.count>=minValid) {
          ok=ok && Pass(cell,minValid,qTol,fourierTol,countAllowance,binTol);
          evaluated=true;
        } else {
          // Sparse cells are explicitly validated on the same-day independent sample.
          // Run-specific constants are not promoted without their own validation.
          ok=ok && cell.model==2 && Pass(dayQA[dayKey],minValid,qTol,fourierTol,countAllowance,binTol);
        }
        if(cell.model==2)
          ok=ok && Pass(dayQA[dayKey],minValid,qTol,fourierTol,countAllowance,binTol);
        if(!ok) good=false;
        double qMax=0., fMax=0., bMax=0.;
        if(cell.shifted.count) {
          for(int a=0;a<nQ;++a) qMax=std::max(qMax,std::fabs(cell.residual.Mean(a)));
          const double n=cell.shifted.count;
          for(int h=0;h<terms;++h) fMax=std::max(fMax,
            std::max(std::fabs(cell.shifted.sine[h]),std::fabs(cell.shifted.cosine[h]))/n);
          for(int b=0;b<8;++b) bMax=std::max(bMax,std::fabs(cell.bins[b]/n-1./8.));
        }
        qa<<run<<'\t'<<cent<<'\t'<<cell.model<<'\t'<<cell.raw.count<<'\t'
          <<day.raw.count<<'\t'<<model.recentered.count<<'\t'<<cell.validationRaw<<'\t'
          <<cell.shifted.count<<'\t'<<cell.invalid<<'\t'<<qMax<<'\t'<<fMax<<'\t'<<bMax<<'\t'<<ok;
        for(int b=0;b<8;++b) qa<<'\t'<<cell.bins[b];
        qa<<'\n';
      }
      if(good && evaluated) accepted.push_back(run);
    }
    qa.close();Require(qa.good(),"QA write failed");
    std::ofstream acceptedFile((dir+"/candidate_validated_runs.txt").c_str());
    for(int run:accepted) acceptedFile<<run<<'\n';
    acceptedFile.close();
    Require(acceptedFile.good(),"run list write failed");
    if(!accepted.empty()) WriteProfiles(dir,accepted,perRun,perDay);
    std::cout<<"EPD_CALIBRATION_FIT_PASS rawEntries="<<entries<<" runs="<<runs.size()
             <<" validatedRuns="<<accepted.size()<<std::endl;
  } catch(const std::exception& e) {
    std::cerr<<"EPD_CALIBRATION_FIT_ERROR: "<<e.what()<<std::endl;
    gSystem->Exit(2);
  }
}
