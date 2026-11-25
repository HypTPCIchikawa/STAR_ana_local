int test_femto_reduced_geometry(const char*,const char*);
#ifndef __CINT__
#include "test_femto_nuclear_pid_domain.C"
#include "FemtoLambdaProvider.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"

namespace {
std::string ReducedWriteYaml(const YAML::Node& value,const std::string& path) {
  std::ofstream out(path.c_str());out<<value;out.close();
  DomainRequire(out.good(),"write reduced-geometry fixture YAML");return path;
}
std::string ReducedFixture(const char* source,const std::string& directory,bool producer,bool compact=false) {
  std::string full(source?source:"");
  DomainRequire(!full.empty(),"explicit source mainconf");
  if(full[0]!='/')full=std::string(gSystem->WorkingDirectory())+"/"+full;
  const size_t marker=full.rfind("/config/");
  DomainRequire(marker!=std::string::npos,"source mainconf below config");
  const std::string base=full.substr(0,marker)+"/config/";
  YAML::Node main=YAML::LoadFile(full);
  for(YAML::const_iterator it=main.begin();it!=main.end();++it) {
    const std::string key=it->first.as<std::string>();
    const std::string path=it->second.as<std::string>();
    if(!path.empty() && path[0]!='/')main[key]=base+path;
  }
  const std::string prefix=directory+"/config/"+(producer?"producer_":(compact?"compact_":"pair_"));
  const char* profiles[]={"femtoHist","nuclearHist","nuclearPidPairQaHist"};
  for(unsigned p=0;p<3;++p) {
    if(!main[profiles[p]])continue;
    YAML::Node profile=YAML::LoadFile(main[profiles[p]].as<std::string>());
    YAML::Node hist=profile["histograms"];
    std::vector<std::string> remove;
    for(YAML::const_iterator it=hist.begin();it!=hist.end();++it) {
      const std::string key=it->first.as<std::string>();
      if(producer && p<2 && (key.find("hKstar")==0 || key.find("hQlab")==0 ||
         key.find("hDphiDeta")==0 || key=="hMixSamplerQA"))remove.push_back(key);
      if(compact && p==2 && key.find("PrePid")!=std::string::npos)remove.push_back(key);
      const char* axes[]={"axis","xAxis","yAxis","zAxis"};
      for(unsigned a=0;a<4;++a) {
        YAML::Node axis=hist[key][axes[a]];
        if(axis.IsMap() && axis["nBins"] && axis["nBins"].as<int>()>16)axis["nBins"]=16;
      }
      if(hist[key]["nBins"] && hist[key]["nBins"].as<int>()>16)hist[key]["nBins"]=16;
    }
    for(size_t i=0;i<remove.size();++i)hist.remove(remove[i]);
    if(!producer && p==0) {
      const char* species[]={"he4","he3"};
      const char* suffix[]={"","_signal","_leftSB","_rightSB"};
      for(int mode=0;mode<2;++mode)for(int s=0;s<2;++s)for(int region=0;region<4;++region) {
        const std::string stem="hKstar"+std::string(mode?"ME":"SE");
        const std::string channel="_lambda_"+std::string(species[s])+suffix[region];
        const std::string key=stem+"VsEpdBinVsCent"+channel;
        YAML::Node item;item["type"]="TH3D";
        item["title"]="synthetic current-event EPD-bin QA;k*;current-event EPD bin;cent9";
        item["xAxis"]=hist[stem+channel]["axis"];
        if(!item["xAxis"].IsMap())item["xAxis"]=hist[stem+channel]["xAxis"];
        item["yAxis"]["nBins"]=8;item["yAxis"]["min"]=-.5;item["yAxis"]["max"]=7.5;
        item["zAxis"]["nBins"]=9;item["zAxis"]["min"]=-.5;item["zAxis"]["max"]=8.5;
        hist[key]=item;
      }
    }
    // The fixture tests geometry and dense QA, not sparse storage.
    if(p==2)profile["jointSparse"]=false;
    main[profiles[p]]=ReducedWriteYaml(profile,prefix+profiles[p]+".yaml");
  }
  YAML::Node nuclear=YAML::LoadFile(main["nuclearid"].as<std::string>());
  if(nuclear["dedxCalibration"] && nuclear["dedxCalibration"]["calibrationFile"]) {
    const std::string path=nuclear["dedxCalibration"]["calibrationFile"].as<std::string>();
    if(!path.empty() && path[0]!='/')nuclear["dedxCalibration"]["calibrationFile"]=base+path;
  }
  main["nuclearid"]=ReducedWriteYaml(nuclear,prefix+"nuclear.yaml");
  return ReducedWriteYaml(main,prefix+"main.yaml");
}
FemtoLambdaLegacy::TrackGeometry ReducedGeometry(const StPicoTrack& track,int index) {
  FemtoLambdaLegacy::TrackGeometry saved;saved.trackIndex=index;
  saved.globalMomentum=track.gMom();saved.origin=track.origin();saved.charge=track.charge();
  return saved;
}
void ReducedSameHistogram(const TH1* pico,const TH1* replay) {
  DomainRequire(pico && replay && pico->GetDimension()==replay->GetDimension() &&
      pico->GetNbinsX()==replay->GetNbinsX() && pico->GetNbinsY()==replay->GetNbinsY() &&
      pico->GetNbinsZ()==replay->GetNbinsZ(),"matching histogram inventory");
  DomainRequire(pico->GetEntries()==replay->GetEntries(),"matching histogram entries");
  DomainRequire(pico->GetSumw2N()==replay->GetSumw2N(),"matching error-storage policy");
  const int cells=(pico->GetNbinsX()+2)*(pico->GetDimension()>1?pico->GetNbinsY()+2:1)*
      (pico->GetDimension()>2?pico->GetNbinsZ()+2:1);
  for(int bin=0;bin<cells;++bin) {
    DomainRequire(pico->GetBinContent(bin)==replay->GetBinContent(bin),"matching contents including flow");
    DomainRequire(pico->GetBinError(bin)==replay->GetBinError(bin),"matching bin errors including flow");
  }
}
unsigned ReducedSameQaDirectory(TDirectory* full,TDirectory* compact) {
  DomainRequire(full && compact,"matching QA output directories");
  unsigned count=0;TIter next(full->GetListOfKeys());
  while(TKey* key=static_cast<TKey*>(next())) {
    const TH1* hist=dynamic_cast<TH1*>(full->Get(key->GetName()));
    if(hist) {
      ReducedSameHistogram(hist,dynamic_cast<TH1*>(compact->Get(key->GetName())));++count;
    } else {
      TDirectory* child=full->GetDirectory(key->GetName());
      if(child)count+=ReducedSameQaDirectory(child,compact->GetDirectory(key->GetName()));
    }
  }
  return count;
}
}
#endif

int test_femto_reduced_geometry(const char* source,const char* output) {
#ifndef __CINT__
  try {
    DomainChecks=0;
    DomainRequire(output && *output,"explicit output directory required");
    std::string directory(output);
    if(directory[0]!='/')directory=std::string(gSystem->WorkingDirectory())+"/"+directory;
    DomainRequire(gSystem->AccessPathName(directory.c_str()),"new output directory required");
    DomainRequire(gSystem->mkdir((directory+"/config").c_str(),kTRUE)==0,"create test output");
    const std::string pairConfig=ReducedFixture(source,directory,false);
    const std::string producerConfig=ReducedFixture(source,directory,true);
    const std::string compactConfig=ReducedFixture(source,directory,false,true);
    FemtoLambdaLegacy pico,replay,producer,compact,compactProducer;
    DomainRequire(pico.Init(pairConfig.c_str(),std::cerr),"Pico helper Init");
    DomainRequire(replay.Init(pairConfig.c_str(),std::cerr),"replay helper Init");
    DomainRequire(producer.Init(producerConfig.c_str(),std::cerr,true),"non-pair producer QA Init");
    DomainRequire(compact.Init(compactConfig.c_str(),std::cerr,false,false),"compact selected-pair QA Init");
    DomainRequire(compactProducer.Init(producerConfig.c_str(),std::cerr,true,false),"compact producer QA Init");
    DomainRequire(!producer.RootHistograms()->Get("hKstarSE_lambda_he4") &&
        !producer.NuclearPidPairQaHistograms() && !producer.NuclearPidPairSparse(),
        "producer does not book pair or sparse histograms");
    DomainRequire(producer.NuclearPidPairQaEnabled(),"producer retains prePID candidate collection");
    DomainRequire(!compact.NuclearPidPairQaEnabled() && compact.NuclearPidPairQaHistograms() &&
        !compact.NuclearPidPairQaHistograms()->Get("hKstarNSigmaPrePid_4He_CentBin8"),
        "compact retains selected-pair QA while omitting prePID histograms and controls");
    DomainRequire(!compactProducer.NuclearPidPairQaEnabled() && !compactProducer.NuclearPidPairQaHistograms(),
        "compact producer omits prePID collection without booking pair QA");
    YAML::Node main=YAML::LoadFile(pairConfig);
    YAML::Node maker=YAML::LoadFile(main["maker"].as<std::string>());
    const double radius=maker["mergingRadiusMeters"].as<double>();
    const double halfLength=maker["mergingGeometry"]["tpcHalfLengthCm"].as<double>();
    const double tolerance=maker["mergingGeometry"]["radiusToleranceCm"].as<double>();
    const double fields[]={0.,5.000000123,-5.000000123};
    const TVector3 momenta[]={TVector3(.85,.24,.13),TVector3(.45,.18,.08),
        TVector3(-.26,.35,-.1),TVector3(.0001,0.,0.),TVector3(.5,.2,20.)};
    unsigned reached=0,missed=0;
    for(int charge=-1;charge<=1;charge+=2)for(unsigned f=0;f<3;++f)for(unsigned m=0;m<5;++m) {
      StPicoTrack track;track.setNHitsFit(charge*25);track.setGlobalMomentum(momenta[m]);
      track.setOrigin(TVector3(.06,-.04,.1));
      const FemtoLambdaLegacy::TrackGeometry saved=ReducedGeometry(track,3);
      TVector3 one,two;
      const bool original=FemtoLambdaLegacy::MergingTrackPositionAtRadius(&track,fields[f],radius,halfLength,tolerance,one);
      const bool stored=FemtoLambdaLegacy::MergingTrackPositionAtRadius(saved,fields[f],radius,halfLength,tolerance,two);
      DomainRequire(original==stored,"same helix reachability for charges and field signs");
      if(original) {
        ++reached;
        DomainRequire(one.X()==two.X() && one.Y()==two.Y() && one.Z()==two.Z(),
            "bitwise-identical cylinder intersection including Float_t field conversion");
      } else ++missed;
    }
    DomainRequire(reached>0 && missed>0,"geometry test includes actual intersections and failures");
    DomainPicoFixture fixture;fixture.dst.event()->setBField(5.);
    const TVector3 global[]={momenta[0],momenta[1],momenta[2]};
    std::vector<FemtoLambdaLegacy::TrackGeometry> bank;
    for(int i=0;i<3;++i) {
      StPicoTrack* track=new(fixture.tracks[i])StPicoTrack();
      track->setNHitsFit(i==2?-25:25);track->setNHitsDedx(20);
      track->setGlobalMomentum(global[i]);track->setPrimaryMomentum(TVector3(1.1,-.15,.4));
      track->setOrigin(TVector3(.06,-.04,.1));bank.push_back(ReducedGeometry(*track,i));
    }
    replay.SetReplayTrackGeometry(bank,fixture.dst.event()->bField());
    compact.SetReplayTrackGeometry(bank,fixture.dst.event()->bField());
    FemtoLambdaCandidate lambda;lambda.protonIndex=1;lambda.pionIndex=2;
    lambda.mass=maker["lambdaSignalMean"].as<float>();
    lambda.candidate.speciesKey="lambda";lambda.candidate.eventIndex=47;
    lambda.candidate.source=kFemtoCandResonance;lambda.candidate.reso.invMass=lambda.mass;
    lambda.candidate.reso.dau1Index=1;lambda.candidate.reso.dau2Index=2;
    lambda.candidate.reso.dau1EventIndex=lambda.candidate.reso.dau2EventIndex=47;
    FemtoCandidate nucleus;nucleus.speciesKey="he4";nucleus.eventIndex=47;
    nucleus.source=kFemtoCandTrack;nucleus.trk.trackIndex=0;nucleus.charge=2;
    nucleus.trk.nSigmaHe4=.3;nucleus.trk.nSigmaHe3=.5;
    nucleus.SetP4(pico.NuclearP4(fixture.dst.track(0)->pMom(),3));
    for(int i=0;i<2;++i) {
      TLorentzVector v;v.SetXYZM(.4,.2,i?-.6:.6,maker["lambdaPairMass"].as<double>());
      lambda.candidate.SetP4(v);
      DomainRequire(pico.FillLambdaAcceptance(lambda,&fixture.dst),"Pico daughter acceptance");
      DomainRequire(replay.FillLambdaAcceptanceFromReplay(lambda),"saved daughter acceptance");
      DomainRequire(compact.FillLambdaAcceptanceFromReplay(lambda),"compact saved daughter acceptance");
      replay.SetPairEventPlaneBin(i?7:0);
      compact.SetPairEventPlaneBin(i?7:0);
      DomainRequire(pico.FillPair(lambda.candidate,nucleus,.15,.3,8,false,&fixture.dst),"Pico signal pair");
      DomainRequire(replay.FillPair(lambda.candidate,nucleus,.15,.3,8,false,0),"saved signal pair");
      DomainRequire(compact.FillPair(lambda.candidate,nucleus,.15,.3,8,false,0),"compact signal pair");
      FemtoCandidate he3=nucleus;he3.speciesKey="he3";
      he3.SetP4(pico.NuclearP4(fixture.dst.track(0)->pMom(),2));
      DomainRequire(pico.FillPair(lambda.candidate,he3,.15,.3,8,false,&fixture.dst),"Pico He3 control pair");
      DomainRequire(replay.FillPair(lambda.candidate,he3,.15,.3,8,false,0),"saved He3 control pair");
      DomainRequire(compact.FillPair(lambda.candidate,he3,.15,.3,8,false,0),"compact He3 control pair");
      const FemtoCandidate nuclei[]={nucleus,he3};
      for(int s=0;s<2;++s) {
        FemtoCandidate mixed=nuclei[s];mixed.eventIndex=48;
        DomainRequire(pico.FillPair(lambda.candidate,mixed,.15,.3,8,true,0),"Pico helper mixed pair");
        DomainRequire(replay.FillPair(lambda.candidate,mixed,.15,.3,8,true,0),"saved helper mixed pair");
        DomainRequire(compact.FillPair(lambda.candidate,mixed,.15,.3,8,true,0),"compact helper mixed pair");
      }
    }
    const char* acceptance[]={"hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal", "hLambdaPion_PtVsYLab_signal"};
    for(unsigned i=0;i<3;++i) {
      ReducedSameHistogram(pico.RootHistograms()->Get(acceptance[i]),replay.RootHistograms()->Get(acceptance[i]));
      ReducedSameHistogram(pico.RootHistograms()->Get(acceptance[i]),compact.RootHistograms()->Get(acceptance[i]));
      DomainRequire(replay.RootHistograms()->Get(acceptance[i])->GetEntries()==2,
          "one acceptance fill per Lambda candidate despite shared daughters");
    }
    const char* daughters[]={"proton","pion"};
    const int nDaughters=maker["mergingGeometry"]["includePion"].as<bool>()?2:1;
    for(int d=0;d<nDaughters;++d)for(int alias=0;alias<2;++alias) {
      const std::string key=std::string("hDphiDeta_")+daughters[d]+"_4He"+(alias?"_EventField":"");
      const TH1* original=pico.NuclearHistograms()->Get(key.c_str());
      const TH1* stored=replay.NuclearHistograms()->Get(key.c_str());
      ReducedSameHistogram(original,stored);
      ReducedSameHistogram(original,compact.NuclearHistograms()->Get(key.c_str()));
      DomainRequire(stored->GetEntries()==2,"merging comparison is nonempty for both accepted signal pairs");
    }
    const char* species[]={"he4","he3"};
    const char* suffix[]={"","_signal","_leftSB","_rightSB"};
    for(int mode=0;mode<2;++mode)for(int s=0;s<2;++s)for(int region=0;region<4;++region) {
      const std::string stem="hKstar"+std::string(mode?"ME":"SE");
      const std::string channel="_lambda_"+std::string(species[s])+suffix[region];
      const std::string key=stem+"VsEpdBinVsCent"+channel;
      TH3* h=dynamic_cast<TH3*>(replay.RootHistograms()->Get(key.c_str()));
      const TH1* plain=replay.RootHistograms()->Get((stem+channel).c_str());
      const TH1* old=pico.RootHistograms()->Get(key.c_str());
      DomainRequire(h && plain && old && old->GetEntries()==0,"EP diagnostic remains disabled without explicit event label");
      ReducedSameHistogram(h,compact.RootHistograms()->Get(key.c_str()));
      DomainRequire(h->GetEntries()==(region<2?2:0),"EP diagnostic fills physical helium SE/ME mass class once");
      if(region<2)DomainRequire(h->GetBinContent(h->FindBin(.15,0,8))==1 &&
          h->GetBinContent(h->FindBin(.15,7,8))==1,"first and last EPD bins receive the supplied current-event label");
      for(int x=0;x<=h->GetNbinsX()+1;++x) {
        double sum=0;
        for(int y=0;y<=h->GetNbinsY()+1;++y)for(int z=0;z<=h->GetNbinsZ()+1;++z)
          sum+=h->GetBinContent(h->GetBin(x,y,z));
        DomainRequire(sum==plain->GetBinContent(x),"sum of EP-bin histogram including flow equals inclusive kstar");
      }
    }
    const char* representations[]={"4He","3He","3HeAs4He"};
    const char* variables[]={"NSigma","POverQ"};
    const std::string signalHalf=replay.SignalHalf(lambda.mass)==1?"_signalLow":"_signalHigh";
    const std::string massClasses[]={"",signalHalf};
    for(int mode=0;mode<2;++mode)for(int r=0;r<3;++r)for(int v=0;v<2;++v)for(int m=0;m<2;++m) {
      const std::string key="hKstar"+std::string(variables[v])+"_"+(mode?"Mixed_":"")+
          representations[r]+massClasses[m]+"_CentBin8";
      const TH1* selected=compact.NuclearPidPairQaHistograms()->Get(key.c_str());
      ReducedSameHistogram(replay.NuclearPidPairQaHistograms()->Get(key.c_str()),selected);
      DomainRequire(selected->GetEntries()==2,"compact selected-pair PID/rigidity QA remains nonempty");
    }
    DomainRequire(!compact.FillPrePidNuclearPairQa(lambda.candidate,nucleus,.15,8,false),
        "compact explicitly disables direct prePID pair QA calls");
    fixture.Add(producer,.9,false);
    FemtoCandidateStore physical;FemtoNuclearPidPairQaStore pidQa;
    producer.CollectNuclei(&fixture.dst,47,8,physical,&pidQa);
    DomainRequire(!pidQa["he4"].empty(),"prePID payload exists without booked pair histograms");
    FemtoCandidateStore compactPhysical;FemtoNuclearPidPairQaStore compactControls;
    compactProducer.CollectNuclei(&fixture.dst,47,8,compactPhysical,&compactControls);
    DomainRequire(compactControls.empty(),"compact producer persists no prePID bank");
    const char* physicalSpecies[]={"deuteron","triton","he3","he4"};
    for(unsigned s=0;s<4;++s) {
      const std::vector<FemtoCandidate>& a=physical[physicalSpecies[s]];
      const std::vector<FemtoCandidate>& b=compactPhysical[physicalSpecies[s]];
      DomainRequire(a.size()==b.size(),"compact keeps the same selected physical multiplicities");
      for(size_t i=0;i<a.size();++i)DomainRequire(a[i].px==b[i].px && a[i].py==b[i].py &&
          a[i].pz==b[i].pz && a[i].energy==b[i].energy,"compact keeps canonical physical P4");
    }
    const std::string root=directory+"/producer_qa.root";
    {TFile file(root.c_str(),"CREATE");DomainRequire(!file.IsZombie() && producer.Write(&file),"write compact producer QA");}
    {TFile file(root.c_str(),"READ");DomainRequire(file.Get("FemtoLambdaProducerQaScope") && !file.GetDirectory("mix"),
        "producer QA is marked and has no mixed-event directory");}
    const std::string compactRoot=directory+"/producer_compact_qa.root";
    {TFile file(compactRoot.c_str(),"CREATE");DomainRequire(!file.IsZombie() && compactProducer.Write(&file),"write compact-profile producer QA");}
    {TFile a(root.c_str(),"READ"),b(compactRoot.c_str(),"READ");
      DomainRequire(ReducedSameQaDirectory(&a,&b)>0,"full and compact producer QA histograms agree exactly");}
    replay.ClearReplayTrackGeometry();
    DomainRequire(!replay.FillLambdaAcceptanceFromReplay(lambda),"clearing geometry prevents stale daughter reuse");
    std::cout<<"PASS reduced geometry: "<<DomainChecks<<" checks, "<<reached<<" intersections, "<<missed<<" misses"<<std::endl;
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"FAIL reduced geometry after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;
    return 1;
  }
#else
  return 1;
#endif
}
