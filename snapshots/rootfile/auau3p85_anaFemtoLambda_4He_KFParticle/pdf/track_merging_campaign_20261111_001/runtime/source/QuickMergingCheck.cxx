#include "StChain.h"
#include "StMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "FemtoCandidateEventStore.h"
#include "FemtoCandidateEventStoreConfig.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TVector2.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace quick_merging {
struct Config {
 std::string sourceMain,cache,picoList,output,json;int bins,centMin,centMax;double halfWidth,radiusM,tpcZ,legacyB,mean,sigma,nSigma,he4Mass;bool prePid;std::vector<double> radii,kstarMax;
 explicit Config(const char* path){
  YAML::Node c=YAML::LoadFile(path);sourceMain=c["sourceMainconf"].as<std::string>();cache=c["candidateRoot"].as<std::string>();picoList=c["picoList"].as<std::string>();output=c["outputRoot"].as<std::string>();json=c["outputJson"].as<std::string>();bins=c["angularBins"].as<int>();halfWidth=c["angularHalfWidth"].as<double>();radiusM=c["radiusM"].as<double>();tpcZ=c["tpcHalfLengthCm"].as<double>();prePid=c["prePidControl"].as<bool>();centMin=c["cent9Min"].as<int>();centMax=c["cent9Max"].as<int>();
  for(size_t i=0;i<c["commonRadiiM"].size();++i)radii.push_back(c["commonRadiiM"][i].as<double>());
  for(size_t i=0;i<c["kstarMax"].size();++i)kstarMax.push_back(c["kstarMax"][i].as<double>());
  if(bins<1||halfWidth<=0||radiusM<=0||radii.empty())throw std::runtime_error("invalid QA axes/radii");
  YAML::Node m=YAML::LoadFile(sourceMain);YAML::Node maker=YAML::LoadFile(femto_candidate_store::ConfigReference(sourceMain,m["maker"].as<std::string>()));
  mean=maker["lambdaSignalMean"].as<double>();sigma=maker["lambdaSignalSigma"].as<double>();nSigma=maker["lambdaSignalNSigma"].as<double>();legacyB=maker["mergingFieldTesla"].as<double>();YAML::Node nuclear=YAML::LoadFile(femto_candidate_store::ConfigReference(sourceMain,m["nuclearid"].as<std::string>()));he4Mass=nuclear["nuclearMassHe4"].as<double>();
 }
};
class Maker:public StMaker {
 Config cfg;femto_candidate_store::Settings store;femto_candidate_store::Reader reader;TChain* pico;TClonesArray* arrays[StPicoArrays::NAllPicoArrays];std::map<std::string,TH2D*> hist;std::map<std::string,Long64_t> count;Long64_t expected;bool ended;std::string failure;
public:
 Maker(const Config& c):StMaker("mergingQuick"),cfg(c),pico(0),expected(0),ended(false){std::fill(arrays,arrays+StPicoArrays::NAllPicoArrays,(TClonesArray*)0);}
 ~Maker(){StPicoDst::unset();delete pico;for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)delete arrays[i];for(std::map<std::string,TH2D*>::iterator i=hist.begin();i!=hist.end();++i)delete i->second;}
 Long64_t Expected()const{return expected;}bool Ended()const{return ended;}const std::string& Failure()const{return failure;}
 Int_t Init(){try{
  {TFile integrity(cfg.cache.c_str(),"READ");if(integrity.IsZombie()||integrity.TestBit(TFile::kRecovered))throw std::runtime_error("candidate ROOT is zombie or recovered");}store.Load(cfg.sourceMain);if(!store.enabled)throw std::runtime_error("source mainconf does not describe original candidate store");
  if(!reader.Init(cfg.cache.c_str(),store.species,store.selectionProvenance,store.mixingProvenance,true,true))throw std::runtime_error(reader.LastError());expected=reader.EntriesTotal();if(expected<1)throw std::runtime_error("empty candidate store");
  pico=new TChain("PicoDst");std::ifstream in(cfg.picoList.c_str());if(!in)throw std::runtime_error("missing explicit Pico list");std::string line;int files=0;
  while(std::getline(in,line)){std::istringstream s(line);std::string path,extra;if(!(s>>path)||path[0]=='#')continue;if(s>>extra)throw std::runtime_error("Pico list requires one path per line");if(pico->Add(path.c_str(),0)!=1)throw std::runtime_error("cannot add Pico file "+path);++files;}
  if(files<1||pico->GetEntries()<1||pico->GetListOfFiles()->GetEntries()!=files)throw std::runtime_error("empty or incomplete Pico chain");{TFile coverage(cfg.cache.c_str(),"READ");TParameter<Long64_t>* inputs=dynamic_cast<TParameter<Long64_t>*>(coverage.Get("inputEvents"));if(!inputs||inputs->GetVal()!=pico->GetEntries())throw std::runtime_error("original Pico input coverage differs from completed cache inputEvents");}pico->SetBranchStatus("*",0);
  arrays[StPicoArrays::Event]=new TClonesArray("StPicoEvent");arrays[StPicoArrays::Track]=new TClonesArray("StPicoTrack");
  const char* names[]={"Event","Track"};const int ids[]={StPicoArrays::Event,StPicoArrays::Track};
  for(int i=0;i<2;++i){if(!pico->GetBranch(names[i]))throw std::runtime_error("missing Pico branch");pico->SetBranchStatus((std::string(names[i])+"*").c_str(),1);if(pico->SetBranchAddress(names[i],&arrays[ids[i]])<0)throw std::runtime_error("cannot bind Pico branch");}
  pico->SetBranchStatus("TrackCovMatrix*",0);StPicoDst::set(arrays);
  const char* species[]={"4He","3He","3HeAs4He"};
  for(int i=0;i<3;++i){H(std::string("hDphiDeta_proton_")+species[i],true);H(std::string("hDphiDeta_proton_")+species[i]+"_EventField",true);}
  return kStOK;
 }catch(const std::exception& e){failure=e.what();std::cerr<<"QA_INIT_FAILED "<<failure<<std::endl;return kStErr;}}
 TH2D* H(const std::string& key,bool old=false){if(hist.count(key))return hist[key];int n=old?200:cfg.bins;double w=old?.25:cfg.halfWidth;TH2D* h=new TH2D(key.c_str(),(key+";#Delta#phi* (rad);#Delta#eta").c_str(),n,-w,w,n,-w,w);h->SetDirectory(0);hist[key]=h;return h;}
 bool Signal(const FemtoCandidate& l){return std::abs(double(l.reso.invMass)-cfg.mean)<=cfg.nSigma*cfg.sigma;}
 double Kstar(const FemtoCandidate& a,const FemtoCandidate& b){TLorentzVector p=a.P4(),q=b.P4(),sum=p+q;p.Boost(-sum.BoostVector());q.Boost(-sum.BoostVector());return .5*(p.Vect()-q.Vect()).Mag();}
 static double Angle(double arg){return std::asin(std::max(-1.,std::min(1.,arg)));}
 void Legacy(const StPicoTrack* p,const StPicoTrack* n,const std::string& species,double field){
  if(p->pPt()<=0||n->pPt()<=0){++count["old_primary_gate_failed_"+species];return;}TVector3 a=p->pMom(),b=n->pMom();double de=a.Eta()-b.Eta();
  double d=a.Phi()-b.Phi()+Angle(.3*cfg.legacyB*cfg.radiusM/(2*a.Pt()))-Angle(.3*2*cfg.legacyB*cfg.radiusM/(2*b.Pt()));H("hDphiDeta_proton_"+species,true)->Fill(TVector2::Phi_mpi_pi(d),de);
  d=a.Phi()-b.Phi()+Angle(.3*p->charge()*field*cfg.radiusM/(2*a.Pt()))-Angle(.3*n->charge()*2*field*cfg.radiusM/(2*b.Pt()*2));H("hDphiDeta_proton_"+species+"_EventField",true)->Fill(TVector2::Phi_mpi_pi(d),de);
 }
 bool Position(const StPicoTrack* t,double bKg,double radius,TVector3& pos){
  StPicoPhysicalHelix h=t->helix(bKg);std::pair<double,double> ss=h.pathLength(radius*100.);double best=1.e99;
  if(std::isfinite(ss.first)&&ss.first>=0&&ss.first<best)best=ss.first;if(std::isfinite(ss.second)&&ss.second>=0&&ss.second<best)best=ss.second;
  if(best>=1.e98)return false;pos=h.at(best);return std::isfinite(pos.X())&&std::isfinite(pos.Y())&&std::isfinite(pos.Z())&&std::abs(pos.Perp()-100*radius)<1.e-3&&std::abs(pos.Z())<=cfg.tpcZ;
 }
 void NewMaps(const StPicoTrack* a,const StPicoTrack* b,const std::string& prefix,double bKg,double pvZ,double k){
  TVector3 pa=a->gMom(),pb=b->gMom();if(pa.Pt()<=0||pb.Pt()<=0){++count["invalid_global_"+prefix];return;}double de=pa.Eta()-pb.Eta();double aa=.299792458*a->charge()*(.1*bKg)*cfg.radiusM/(2*pa.Pt()),ab=.299792458*b->charge()*(.1*bKg)*cfg.radiusM/(2*pb.Pt());
  std::vector<std::string> slices(1,"");for(size_t i=0;i<cfg.kstarMax.size();++i)if(k<cfg.kstarMax[i]){std::ostringstream s;s<<"_Klt"<<int(cfg.kstarMax[i]*1000+.5);slices.push_back(s.str());}
  if(std::abs(aa)<=1&&std::abs(ab)<=1){double d=TVector2::Phi_mpi_pi(pa.Phi()-pb.Phi()-std::asin(aa)+std::asin(ab));for(size_t i=0;i<slices.size();++i)H(prefix+"_GlobalMinusAsin"+slices[i])->Fill(d,de);++count["global_formula_valid_"+prefix];}else ++count["global_unreachable_radius_"+prefix];
  TVector3 posa,posb;if(Position(a,bKg,cfg.radiusM,posa)&&Position(b,bKg,cfg.radiusM,posb)){
   double d=TVector2::Phi_mpi_pi(posa.Phi()-posb.Phi());for(size_t i=0;i<slices.size();++i){H(prefix+"_HelixFixedGlobalEta"+slices[i])->Fill(d,de);double da=std::asinh((posa.Z()-pvZ)/(100*cfg.radiusM))-std::asinh((posb.Z()-pvZ)/(100*cfg.radiusM));H(prefix+"_HelixFixedSpatialEta"+slices[i])->Fill(d,da);}++count["helix_fixed_valid_"+prefix];
  }else ++count["helix_fixed_unreachable_"+prefix];
  double best=1.e99,bd=0,be=0,beg=de,br=0;
  for(size_t i=0;i<cfg.radii.size();++i)if(Position(a,bKg,cfg.radii[i],posa)&&Position(b,bKg,cfg.radii[i],posb)){
   double d=TVector2::Phi_mpi_pi(posa.Phi()-posb.Phi()),dz=posa.Z()-posb.Z(),distance=std::sqrt(std::pow(100*cfg.radii[i]*d,2)+dz*dz);
   if(distance<best){best=distance;bd=d;be=std::asinh((posa.Z()-pvZ)/(100*cfg.radii[i]))-std::asinh((posb.Z()-pvZ)/(100*cfg.radii[i]));br=cfg.radii[i];}
  }
  if(br>0){for(size_t i=0;i<slices.size();++i){H(prefix+"_HelixClosestGlobalEta"+slices[i])->Fill(bd,beg);H(prefix+"_HelixClosestSpatialEta"+slices[i])->Fill(bd,be);}++count["helix_common_valid_"+prefix];}else ++count["helix_no_common_radius_"+prefix];
 }
 void Pair(const FemtoCandidate& l,const FemtoCandidate& n,const std::string& species,const std::string& stage,const femto_candidate_store::Event& e){
  if(!Signal(l)||FemtoCandidatesShareTrack(l,n))return;FemtoCandidate mapped=n;if(species=="3HeAs4He"){TLorentzVector v=n.P4();v.SetVectM(v.Vect(),cfg.he4Mass);mapped.SetP4(v);}double k=Kstar(l,mapped);if(!std::isfinite(k))throw std::runtime_error("nonfinite cached pair kstar");++count["signal_pairs_"+stage+"_"+species];
  if(l.reso.dau1Index<0||l.reso.dau2Index<0||n.trk.trackIndex<0||l.reso.dau1Index>=e.nTracks||l.reso.dau2Index>=e.nTracks||n.trk.trackIndex>=e.nTracks)throw std::runtime_error("bad track reference");
  StPicoTrack* p=StPicoDst::track(l.reso.dau1Index);StPicoTrack* pi=StPicoDst::track(l.reso.dau2Index);StPicoTrack* nuc=StPicoDst::track(n.trk.trackIndex);if(!p||!pi||!nuc||p->charge()==0||pi->charge()==0||nuc->charge()==0)throw std::runtime_error("invalid daughter/nuclear track identity");if(p->charge()<0||pi->charge()>0)++count["unexpected_daughter_charge_"+stage];if(nuc->charge()<0)++count["negative_nuclear_charge_"+stage];
  if(stage=="selected")Legacy(p,nuc,species,.1*e.bField);
  if(e.cent9<cfg.centMin||e.cent9>cfg.centMax)return;
  std::string prefix="hDphiDeta_"+stage+"_"+species;NewMaps(p,nuc,prefix+"_proton",e.bField,e.pvZ,k);NewMaps(pi,nuc,prefix+"_pion",e.bField,e.pvZ,k);
 }
 Int_t Make(){try{
  femto_candidate_store::Event e;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;if(!reader.Next(e,physical,qa)){if(reader.Finished()){ended=true;return kStEOF;}throw std::runtime_error(reader.LastError());}++count["cache_events"];if(physical.count("lambda"))count["cache_lambda_candidates"]+=physical["lambda"].size();
  FemtoCandidateStore::const_iterator li=physical.find("lambda");if(li==physical.end()||li->second.empty())return kStOK;bool anyLambda=false;for(size_t i=0;i<li->second.size();++i)if(Signal(li->second[i]))anyLambda=true;if(!anyLambda)return kStOK;
  bool nucleus=physical.count("he4")&&!physical["he4"].empty()||physical.count("he3")&&!physical["he3"].empty()||cfg.prePid&&qa.count("he4")&&!qa["he4"].empty();if(!nucleus)return kStOK;
  if(e.inputIndex<1||e.inputIndex>pico->GetEntries())throw std::runtime_error("cache inputIndex outside original Pico chain");if(pico->GetEntry(e.inputIndex-1)<=0)throw std::runtime_error("Pico event read failed");++count["pico_events_read"];
  StPicoEvent* pe=StPicoDst::event();if(!pe||pe->runId()!=e.runId||pe->eventId()!=e.eventId||StPicoDst::numberOfTracks()!=static_cast<unsigned>(e.nTracks))throw std::runtime_error("Pico/cache identity mismatch");TVector3 pv=pe->primaryVertex();if(std::abs(pv.X()-e.pvX)>1.e-6||std::abs(pv.Y()-e.pvY)>1.e-6||std::abs(pv.Z()-e.pvZ)>1.e-6||std::abs(pe->bField()-e.bField)>1.e-6)throw std::runtime_error("Pico/cache vertex or field mismatch");
  const char* keys[]={"he4","he3"};const char* labels[]={"4He","3He"};for(int sp=0;sp<2;++sp)if(physical.count(keys[sp]))for(size_t i=0;i<li->second.size();++i)for(size_t j=0;j<physical[keys[sp]].size();++j){Pair(li->second[i],physical[keys[sp]][j],labels[sp],"selected",e);if(sp==1)Pair(li->second[i],physical[keys[sp]][j],"3HeAs4He","selected",e);}
  if(cfg.prePid&&qa.count("he4"))for(size_t i=0;i<li->second.size();++i)for(size_t j=0;j<qa["he4"].size();++j)Pair(li->second[i],qa["he4"][j].candidate,"4He","prePIDControl",e);
  return kStOK;
 }catch(const std::exception& e){failure=e.what();std::cerr<<"QA_MAKE_FAILED "<<failure<<std::endl;return kStErr;}}
 Int_t Finish(){try{
  if(!failure.empty())throw std::runtime_error(failure);if(count["cache_events"]!=expected)throw std::runtime_error("did not process declared cache frames");
  TFile reference(cfg.cache.c_str(),"READ");TTree* eventLedger=dynamic_cast<TTree*>(reference.Get("eventLedger"));TTree* selectedLambda=dynamic_cast<TTree*>(reference.Get("selectedLambda"));TParameter<Long64_t>* inputCount=dynamic_cast<TParameter<Long64_t>*>(reference.Get("inputEvents"));TParameter<Long64_t>* reconstructed=dynamic_cast<TParameter<Long64_t>*>(reference.Get("reconstructedEvents"));TParameter<Long64_t>* selectedParameter=dynamic_cast<TParameter<Long64_t>*>(reference.Get("selectedLambdaCandidates"));if(!eventLedger||!selectedLambda||!inputCount||!reconstructed||!selectedParameter||eventLedger->GetEntries()!=inputCount->GetVal()||reconstructed->GetVal()!=expected||selectedLambda->GetEntries()!=count["cache_lambda_candidates"]||selectedLambda->GetEntries()!=selectedParameter->GetVal())throw std::runtime_error("candidate/ledger/event count integrity failure");const char* labels[]={"4He","3He","3HeAs4He"};const char* canonical[]={"he4","he3","he3_as_he4"};Long64_t checks=0;
  for(int sp=0;sp<3;++sp){TH1* signal=dynamic_cast<TH1*>(reference.Get((std::string("hKstarSE_lambda_")+canonical[sp]+"_signal").c_str()));if(!signal||signal->GetEntries()!=count["signal_pairs_selected_"+std::string(labels[sp])])throw std::runtime_error("physical pair closure failed");++checks;
   for(int mode=0;mode<2;++mode){std::string key="hDphiDeta_proton_"+std::string(labels[sp])+(mode?"_EventField":"");TH2* old=dynamic_cast<TH2*>(reference.Get(("true/"+key).c_str()));TH2D* made=H(key,true);if(!old||old->GetEntries()!=made->GetEntries())throw std::runtime_error("old QA Entries closure failed "+key);for(int ix=0;ix<=201;++ix)for(int iy=0;iy<=201;++iy){if(old->GetBinContent(ix,iy)!=made->GetBinContent(ix,iy))throw std::runtime_error("old QA bin closure failed "+key);++checks;}}
  }
  TFile result(cfg.output.c_str(),"CREATE");if(result.IsZombie())throw std::runtime_error("cannot create new QA ROOT");for(std::map<std::string,TH2D*>::iterator i=hist.begin();i!=hist.end();++i)i->second->Write();for(std::map<std::string,Long64_t>::const_iterator i=count.begin();i!=count.end();++i){TParameter<Long64_t> c(i->first.c_str(),i->second);c.Write();}TNamed status("QuickMergingStatus","completed");status.Write();result.Close();
  std::ofstream o(cfg.json.c_str());o<<std::setprecision(17)<<"{\"passed\":true,\"cacheRoot\":\""<<cfg.cache<<"\",\"outputRoot\":\""<<cfg.output<<"\",\"closureChecks\":"<<checks<<",\"histogramCount\":"<<hist.size()<<",\"counters\":{";bool comma=false;for(std::map<std::string,Long64_t>::const_iterator i=count.begin();i!=count.end();++i){if(comma)o<<',';comma=true;o<<'\"'<<i->first<<"\":"<<i->second;}o<<"},\"scope\":\"Exact frozen selected candidates, raw SE only; no KF/PID/event-selection rerun. prePIDControl excludes TPC own-pull and best species only.\"}\n";o.close();if(!o)throw std::runtime_error("cannot write summary JSON");std::cout<<"QUICK_MERGING_PASS cache="<<count["cache_events"]<<" picoRead="<<count["pico_events_read"]<<" He4="<<count["signal_pairs_selected_4He"]<<" hist="<<hist.size()<<" closure="<<checks<<std::endl;return kStOK;
 }catch(const std::exception& e){failure=e.what();std::cerr<<"QA_FINISH_FAILED "<<failure<<std::endl;return kStErr;}}
};
}
extern "C" int RunQuickMergingCheck(const char* settings){try{
 quick_merging::Config cfg(settings);if(!gSystem->AccessPathName(cfg.output.c_str())||!gSystem->AccessPathName(cfg.json.c_str()))throw std::runtime_error("refuse existing output");TStopwatch timer;timer.Start();StChain chain("quickMergingChain");quick_merging::Maker* maker=new quick_merging::Maker(cfg);if(chain.Init()!=kStOK||!maker->Failure().empty())return 2;Long64_t n=maker->Expected();for(Long64_t i=0;i<n;++i){chain.Clear();if(chain.Make(i)!=kStOK||!maker->Failure().empty())return 3;}if(chain.Finish()!=kStOK||!maker->Failure().empty())return 4;timer.Stop();std::cout<<"QUICK_MERGING_TIME real="<<timer.RealTime()<<" cpu="<<timer.CpuTime()<<std::endl;return 0;
 }catch(const std::exception& e){std::cerr<<"QUICK_MERGING_FAILED "<<e.what()<<std::endl;return 1;}}
