#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TArrayD.h>
#include <TSystem.h>
#include <TMath.h>
#include <TError.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
namespace MatchedClosure {
void Need(bool good,const std::string& text){if(!good)throw std::runtime_error(text);}
std::string J(const std::string& text){std::ostringstream s;s<<'"';for(size_t i=0;i<text.size();++i){unsigned char c=text[i];if(c=='"'||c=='\\')s<<'\\'<<char(c);else if(c=='\n')s<<"\\n";else if(c=='\r')s<<"\\r";else if(c=='\t')s<<"\\t";else if(c<32)s<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else s<<char(c);}return s.str()+'"';}
bool Starts(const std::string& s,const char* prefix){return s.compare(0,std::string(prefix).size(),prefix)==0;}
std::string Base(const std::string& s){size_t p=s.find_last_of('/');return p==std::string::npos?s:s.substr(p+1);}
bool Pair(const std::string& path){std::string s=Base(path);return Starts(s,"hKstar")||Starts(s,"hQlab")||Starts(s,"hDphi")||Starts(s,"hNuclearPidPair");}
bool PrePid(const std::string& path){return path.find("NSigmaPrePid")!=std::string::npos||path.find("POverQPrePid")!=std::string::npos;}
void Inventory(TDirectory* d,const std::string& prefix,std::map<std::string,TKey*>& out){TIter next(d->GetListOfKeys());TKey* key=0;std::set<std::string> seen;while((key=dynamic_cast<TKey*>(next()))){std::string n=key->GetName();if(!seen.insert(n).second)continue;key=d->GetKey(n.c_str());TClass* c=TClass::GetClass(key->GetClassName());Need(c!=0,"unknown class");if(c->InheritsFrom(TDirectory::Class()))Inventory(d->GetDirectory(n.c_str()),prefix+n+"/",out);else if(c->InheritsFrom(TH1::Class()))out[prefix+n]=key;}}
bool Axis(TAxis* a,TAxis* b){if(a->GetNbins()!=b->GetNbins())return false;for(int i=1;i<=a->GetNbins()+1;++i)if(a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;for(int i=1;i<=a->GetNbins();++i)if(std::string(a->GetBinLabel(i))!=b->GetBinLabel(i))return false;return true;}
int Cells(TH1* h){int n=h->GetNbinsX()+2;if(h->GetDimension()>1)n*=h->GetNbinsY()+2;if(h->GetDimension()>2)n*=h->GetNbinsZ()+2;return n;}
void Equal(TKey* ak,TKey* bk,const std::string& name,Long64_t& bins){TH1* a=dynamic_cast<TH1*>(ak->ReadObj());TH1* b=dynamic_cast<TH1*>(bk->ReadObj());Need(a&&b,"cannot read "+name);a->SetDirectory(0);b->SetDirectory(0);Need(std::string(a->ClassName())==b->ClassName()&&a->GetDimension()==b->GetDimension()&&Cells(a)==Cells(b)&&a->GetEntries()==b->GetEntries()&&a->GetSumw2N()==b->GetSumw2N(),"shape/entries/sumw2 mismatch "+name);Need(Axis(a->GetXaxis(),b->GetXaxis())&&(a->GetDimension()<2||Axis(a->GetYaxis(),b->GetYaxis()))&&(a->GetDimension()<3||Axis(a->GetZaxis(),b->GetZaxis())),"axis mismatch "+name);for(int k=0;k<Cells(a);++k){Need(a->GetBinContent(k)==b->GetBinContent(k)&&a->GetBinError(k)==b->GetBinError(k)&&TMath::Finite(a->GetBinContent(k))&&TMath::Finite(a->GetBinError(k)),"bin content/error mismatch "+name);if(a->GetSumw2N())Need(a->GetSumw2()->At(k)==b->GetSumw2()->At(k),"raw variance mismatch "+name);}bins+=Cells(a);delete a;delete b;}
void Metadata(TFile* f,std::ostream& out){out<<'{';TIter next(f->GetListOfKeys());TKey* k=0;std::set<std::string> seen;bool first=true;while((k=dynamic_cast<TKey*>(next()))){std::string name=k->GetName();if(!seen.insert(name).second)continue;TClass* c=TClass::GetClass(k->GetClassName());if(!c||c->InheritsFrom(TH1::Class())||c->InheritsFrom(TDirectory::Class()))continue;TObject* o=f->Get(name.c_str());TParameter<Long64_t>* ll=dynamic_cast<TParameter<Long64_t>*>(o);TParameter<Int_t>* ii=dynamic_cast<TParameter<Int_t>*>(o);TTree* t=dynamic_cast<TTree*>(o);TObjString* os=dynamic_cast<TObjString*>(o);TNamed* named=dynamic_cast<TNamed*>(o);if(!(ll||ii||t||os||named))continue;if(name.find("HistogramDefinitions")!=std::string::npos)continue;if(!first)out<<',';first=false;out<<J(name)<<':';if(ll)out<<ll->GetVal();else if(ii)out<<ii->GetVal();else if(t)out<<"{\"class\":\"TTree\",\"entries\":"<<t->GetEntries()<<",\"payloadRead\":false}";else if(os)out<<J(os->GetString().Data());else out<<J(named->GetTitle());}out<<'}';}
void Raw(TH1* h,std::ostream& out){Need(h&&h->GetDimension()==1,"expected raw1D");out<<"{\"entries\":"<<h->GetEntries()<<",\"axisTitle\":"<<J(h->GetXaxis()->GetTitle())<<",\"bins\":"<<h->GetNbinsX()<<",\"edges\":[";for(int i=1;i<=h->GetNbinsX()+1;++i){if(i>1)out<<',';out<<h->GetXaxis()->GetBinLowEdge(i);}out<<"],\"contentIncludingFlow\":[";for(int i=0;i<=h->GetNbinsX()+1;++i){if(i)out<<',';out<<h->GetBinContent(i);}out<<"],\"errorIncludingFlow\":[";for(int i=0;i<=h->GetNbinsX()+1;++i){if(i)out<<',';out<<h->GetBinError(i);}out<<"]}";}
}
void audit_matched_closure(){using namespace MatchedClosure;TH1::AddDirectory(kFALSE);const char* dest="rootfile/lambda_ep8_cf_closure_report_20261120_001/closure_root_audit.json";Need(gSystem->AccessPathName(dest),"refuse existing audit");std::ofstream out(dest);out<<std::setprecision(17);const std::string d="rootfile/femto_lambda_tree_implementation_20261117_001/";const char* names[]={"baseline_10k.root","producer_compact_10k.root","compact_replay_10k.root"};TFile* f[3]={0};std::map<std::string,TKey*> hist[3];try{for(int i=0;i<3;++i){f[i]=TFile::Open((d+names[i]).c_str(),"READ");Need(f[i]&&!f[i]->IsZombie()&&!f[i]->TestBit(TFile::kRecovered),"bad input ROOT");Inventory(f[i],"",hist[i]);}Long64_t bins=0;int qa=0,pairs=0,nonempty=0,omitted=0;std::vector<std::string> pairNames;
for(std::map<std::string,TKey*>::const_iterator it=hist[1].begin();it!=hist[1].end();++it){Need(hist[0].count(it->first)!=0,"producerQA absent baseline");Equal(hist[0][it->first],it->second,it->first,bins);++qa;}
for(std::map<std::string,TKey*>::const_iterator it=hist[0].begin();it!=hist[0].end();++it){if(!Pair(it->first))continue;if(PrePid(it->first)){++omitted;continue;}Need(hist[2].count(it->first)!=0,"pair absent replay "+it->first);Equal(it->second,hist[2][it->first],it->first,bins);pairNames.push_back(it->first);++pairs;TH1* h=dynamic_cast<TH1*>(it->second->ReadObj());if(h->GetEntries()>0)++nonempty;delete h;}
Need(qa==108&&pairs==1186&&omitted==540&&nonempty>0,"unexpected compact histogram scope");out<<"{\"passed\":true,\"scope\":\"saved matched 10000-Pico-event direct-Pico versus compact tree writer/replay; no event replay\",\"producerQACompared\":"<<qa<<",\"pairHistogramsCompared\":"<<pairs<<",\"nonemptyPairHistograms\":"<<nonempty<<",\"omittedPrePidHistograms\":"<<omitted<<",\"flowInclusiveBinsCompared\":"<<bins<<",\"allContentErrorEntriesAxesAndRawSumw2Exact\":true,\"files\":[";for(int i=0;i<3;++i){if(i)out<<',';out<<"{\"path\":"<<J(d+names[i])<<",\"histograms\":"<<hist[i].size()<<",\"metadata\":";Metadata(f[i],out);out<<"}";}out<<"],\"histograms\":{";
const char* species[]={"he4","he3"};const char* tags[]={"","_signal","_leftSB","_rightSB"};bool first=true;
for(int s=0;s<2;++s)for(int t=0;t<4;++t){std::string suffix=std::string("_lambda_")+species[s]+tags[t];for(int e=0;e<2;++e){std::string name=std::string("hKstar")+(e?"ME":"SE")+suffix;if(!first)out<<',';first=false;out<<J(name)<<":{\"old\":";Raw(dynamic_cast<TH1*>(f[0]->Get(name.c_str())),out);out<<",\"treeReplay\":";Raw(dynamic_cast<TH1*>(f[2]->Get(name.c_str())),out);out<<"}";}}
out<<"},\"comparedPairNames\":[";for(size_t j=0;j<pairNames.size();++j){if(j)out<<',';out<<J(pairNames[j]);}out<<"],\"eventTreePayloadRead\":false}\n";out.close();for(int i=0;i<3;++i){f[i]->Close();delete f[i];}std::cout<<"MATCHED_CLOSURE_AUDIT_PASS qa="<<qa<<" pairs="<<pairs<<" bins="<<bins<<std::endl;}catch(const std::exception& err){out.close();std::cerr<<"MATCHED_CLOSURE_AUDIT_FAIL "<<err.what()<<std::endl;gSystem->Exit(1);}}
