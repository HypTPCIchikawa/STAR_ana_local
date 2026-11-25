void ExtractHe4Y(const char* planPath,const char* requirementsPath,const char* outputPath);
#ifndef __CINT__
#include "../femto_lambda_ep_merging_systematics_20261120_001/ExtractPeriodCF.C"
namespace He4YExtract {
using namespace PeriodCFExtract;
struct Input {int period;Long64_t events;std::string source,path,reference;};
struct Bank {
 int nx,ny; std::vector<double> values,variances;
 Bank():nx(0),ny(0){}
 size_t At(int m,int y,int c,int x)const{return (((m*(ny+3)+(y+1))*9+c)*(nx+2)+x);}
};
void Exact(const TH1* a,const TH1* b,const std::string& key) {
 Need(a&&b&&a->ClassName()==std::string(b->ClassName())&&a->GetDimension()==b->GetDimension(),"reference type mismatch: "+key);
 Need(Axis(a->GetXaxis(),b->GetXaxis()),"reference X axis mismatch: "+key);
 if(a->GetDimension()>1)Need(Axis(a->GetYaxis(),b->GetYaxis()),"reference Y axis mismatch: "+key);
 const int cells=(a->GetNbinsX()+2)*(a->GetDimension()>1?a->GetNbinsY()+2:1)*(a->GetDimension()>2?a->GetNbinsZ()+2:1);
 Need(a->GetEntries()==b->GetEntries(),"reference entries mismatch: "+key);
 for(int i=0;i<cells;++i)Need(a->GetBinContent(i)==b->GetBinContent(i)&&a->GetBinError(i)==b->GetBinError(i),"reference count/error mismatch: "+key);
}
bool Close(double a,double b,double at,double rt){return std::fabs(a-b)<=at+rt*std::max(std::fabs(a),std::fabs(b));}
void SourceIdentity(TFile* file,const Input& row) {
 Need(Object(file,"FemtoLambdaRunStatus")=="completed"&&Object(file,"FemtoLambdaReducedTreeMode")=="replay","incomplete replay ROOT");
 Need(Counter(file,"inputEvents")==row.events&&Counter(file,"reconstructedEvents")==row.events,"processed event count mismatch");
 std::istringstream in(Object(file,"FemtoLambdaReducedSourceInput"));std::string line,argument;int count=0;
 while(std::getline(in,line)){line=Trim(line);if(line.empty())continue;if(line.compare(0,16,"source argument:")==0){Need(argument.empty(),"duplicate source argument");argument=Trim(line.substr(16));continue;}size_t slash=line.find_last_of('/');Need(line.substr(slash==std::string::npos?0:slash+1)==row.source,"source period identity differs");++count;}
 Need(!argument.empty(),"missing source argument");
 const bool direct=argument.size()>=5&&argument.compare(argument.size()-5,5,".root")==0;
 if(direct){const size_t slash=argument.find_last_of('/');Need(argument.substr(slash==std::string::npos?0:slash+1)==row.source&&count==0,"direct source period identity differs");}
 else Need(count==1,"not exactly one period source");
}
}
void ExtractHe4Y(const char* planPath,const char* requirementsPath,const char* outputPath) {
 using namespace He4YExtract;ErrorHandlerFunc_t previous=SetErrorHandler(Handler);rootError=false;TH1::AddDirectory(kFALSE);
 try {
  Need(planPath&&requirementsPath&&outputPath&&gSystem->AccessPathName(outputPath),"explicit arguments and fresh output required");
  std::map<std::string,std::string> req;std::string line;std::ifstream rs(requirementsPath);Need(bool(rs),"requirements missing");
  while(std::getline(rs,line)){std::vector<std::string> f=Fields(line);Need(f.size()==2&&req.insert(std::make_pair(f[0],f[1])).second,"bad requirements row");}Need(rs.eof(),"requirements read error");
  const char* keys[]={"nRapidityBins","rapidityMin","rapidityMax","absoluteTolerance","relativeTolerance","metadataBins","metadataMin","metadataMax","metadataDefinition"};
  for(size_t i=0;i<sizeof(keys)/sizeof(*keys);++i)Need(req.count(keys[i]),"missing extraction requirement");
  const int ny=int(Integer(req["nRapidityBins"]));const double ymin=Number(req["rapidityMin"]),ymax=Number(req["rapidityMax"]),at=Number(req["absoluteTolerance"]),rt=Number(req["relativeTolerance"]);
  Need(ny>0&&ny%2==0&&ymin<ymax&&at>=0&&rt>=0,"invalid extraction requirements");
  std::ifstream ps(planPath);Need(bool(ps),"extraction plan missing");std::vector<Input> rows;std::set<std::string> paths;
  while(std::getline(ps,line)){std::vector<std::string> f=Fields(line);Need(f.size()==5,"plan requires period/events/source/output/reference");Input row;row.period=int(Integer(f[0]));row.events=Integer(f[1]);row.source=f[2];row.path=f[3];row.reference=f[4];Need(row.period>=0&&row.events>0&&paths.insert(row.path).second,"invalid duplicate plan");rows.push_back(row);}Need(ps.eof()&&!rows.empty(),"empty/failed plan");
  Bank bank;bank.ny=ny;std::vector<double> xedges,yedges;std::string selection,epDefinition,reducedDefinition,yDefinition;std::map<std::string,std::string> mixCommon;std::ostringstream sources;sources<<'[';Long64_t events=0,epSeen=0,epValid=0,epRejected=0;Long64_t checks=0;double maxDifference=0;
  const char* mass[]={"","_signal","_leftSB","_rightSB"};
  for(size_t ir=0;ir<rows.size();++ir){const Input& row=rows[ir];TFile* file=TFile::Open(row.path.c_str(),"READ");Need(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered),"bad output ROOT: "+row.path);SourceIdentity(file,row);
   TFile* ref=0;if(row.reference!="-"){ref=TFile::Open(row.reference.c_str(),"READ");Need(ref&&!ref->IsZombie()&&!ref->TestBit(TFile::kRecovered),"bad reference ROOT");SourceIdentity(ref,row);}
   const std::string sel=Object(file,"FemtoLambdaReducedReplaySelection"),ep=Object(file,"FemtoEpdEventPlaneDefinition"),red=Object(file,"FemtoLambdaReducedTreeSettings"),yd=Object(file,req["metadataDefinition"].c_str());
   std::map<std::string,std::string> mix=Scalars(Object(file,"FemtoLambdaReducedReplayMixing")),reduced=Scalars(red);
   for(std::map<std::string,std::string>::const_iterator it=req.begin();it!=req.end();++it){if(it->first.compare(0,4,"mix.")==0)Need(mix[it->first.substr(4)]==it->second,"mixing mismatch "+it->first);if(it->first.compare(0,8,"reduced.")==0)Need(reduced[it->first.substr(8)]==it->second,"reduced mismatch "+it->first);if(it->first.compare(0,11,"definition.")==0)Need(yd.find(it->second)!=std::string::npos,"rapidity definition mismatch "+it->first);}
   TParameter<Int_t>* yn=dynamic_cast<TParameter<Int_t>*>(file->Get(req["metadataBins"].c_str()));TParameter<Double_t>* yl=dynamic_cast<TParameter<Double_t>*>(file->Get(req["metadataMin"].c_str()));TParameter<Double_t>* yh=dynamic_cast<TParameter<Double_t>*>(file->Get(req["metadataMax"].c_str()));Need(yn&&yl&&yh&&yn->GetVal()==ny&&yl->GetVal()==ymin&&yh->GetVal()==ymax,"rapidity metadata differs");
   TParameter<Int_t>* cut=dynamic_cast<TParameter<Int_t>*>(file->Get("FemtoLambdaMergingPairCutEnabled"));Need(!cut||cut->GetVal()==0,"reference analysis must be pair-cut Off");
   if(ir==0){selection=sel;epDefinition=ep;reducedDefinition=red;yDefinition=yd;mixCommon=mix;}Need(sel==selection&&ep==epDefinition&&red==reducedDefinition&&yd==yDefinition&&mix==mixCommon,"selection/EP/reduced/binning differ across periods");
   if(ref){Need(Object(ref,"FemtoLambdaReducedReplaySelection")==sel&&Object(ref,"FemtoEpdEventPlaneDefinition")==ep&&Object(ref,"FemtoLambdaReducedTreeSettings")==red&&Scalars(Object(ref,"FemtoLambdaReducedReplayMixing"))==mix,"reference reconstruction/mixing differs");const char* counters[]={"epdEventsSeen","epdEventsValid","epdEventsRejected","selectedLambdaCandidates"};for(size_t i=0;i<4;++i)Need(Counter(ref,counters[i])==Counter(file,counters[i]),"reference event population differs");}
   for(int mode=0;mode<2;++mode)for(int ma=0;ma<4;++ma){const std::string stem="hKstar"+std::string(mode?"ME":"SE"),channel="_lambda_he4"+std::string(mass[ma]),k1=stem+channel,k2=stem+"VsCent"+channel,k3=stem+"VsHe4YLabVsCent"+channel;TH1* h1=dynamic_cast<TH1*>(file->Get(k1.c_str()));TH2* h2=dynamic_cast<TH2*>(file->Get(k2.c_str()));TH3D* h3=dynamic_cast<TH3D*>(file->Get(k3.c_str()));Need(h1&&h2&&h3&&Axis(h1->GetXaxis(),h2->GetXaxis())&&Axis(h2->GetXaxis(),h3->GetXaxis())&&Axis(h2->GetYaxis(),h3->GetZaxis()),"missing/mismatched He4Y histograms");IntegerAxis(h2->GetYaxis(),9);IntegerAxis(h3->GetZaxis(),9);Need(h3->GetYaxis()->GetNbins()==ny&&h3->GetYaxis()->GetXmin()==ymin&&h3->GetYaxis()->GetXmax()==ymax,"He4Y axis differs");for(int y=1;y<=ny+1;++y)Need(Close(h3->GetYaxis()->GetBinLowEdge(y),ymin+(y-1)*(ymax-ymin)/ny,1.e-14,0),"nonuniform Y axis");
    if(ref){Exact(h1,dynamic_cast<TH1*>(ref->Get(k1.c_str())),k1);Exact(h2,dynamic_cast<TH1*>(ref->Get(k2.c_str())),k2);}
    const int nx=h2->GetNbinsX();if(bank.nx==0){bank.nx=nx;bank.values.assign(2*(ny+3)*9*(nx+2),0.);bank.variances=bank.values;for(int x=1;x<=nx+1;++x)xedges.push_back(h2->GetXaxis()->GetBinLowEdge(x));for(int y=1;y<=ny+1;++y)yedges.push_back(h3->GetYaxis()->GetBinLowEdge(y));}Need(nx==bank.nx,"kstar bins differ");for(int x=1;x<=nx+1;++x)Need(h2->GetXaxis()->GetBinLowEdge(x)==xedges[x-1],"kstar edges differ");
    for(int x=0;x<=nx+1;++x)for(int cent=0;cent<9;++cent){double sum=0,var=0;for(int y=0;y<=ny+1;++y){const double count=h3->GetBinContent(x,y,cent+1),er=h3->GetBinError(x,y,cent+1);Need(TMath::Finite(count)&&TMath::Finite(er)&&count>=0&&er>=0,"invalid raw Y cell");sum+=count;var+=er*er;if(ma==1){size_t atIndex=bank.At(mode,y,cent,x);bank.values[atIndex]+=count;bank.variances[atIndex]+=er*er;}}
     const double count=h2->GetBinContent(x,cent+1),er=h2->GetBinError(x,cent+1);maxDifference=std::max(maxDifference,std::fabs(sum-count));Need(Close(sum,count,at,rt)&&Close(var,er*er,at,rt),"He4Y sum including flow does not close to inclusive");++checks;if(ma==1){size_t atIndex=bank.At(mode,-1,cent,x);bank.values[atIndex]+=count;bank.variances[atIndex]+=er*er;}}
    for(int x=0;x<=nx+1;++x){Need(h2->GetBinContent(x,0)==0&&h2->GetBinContent(x,10)==0,"inclusive centrality flow nonzero");for(int y=0;y<=ny+1;++y)Need(h3->GetBinContent(x,y,0)==0&&h3->GetBinContent(x,y,10)==0,"He4Y centrality flow nonzero");}
   }
   if(ir)sources<<',';sources<<"{\"period\":"<<row.period<<",\"path\":"<<Json(row.path)<<",\"reference\":"<<Json(row.reference)<<",\"source\":"<<Json(row.source)<<",\"inputEvents\":"<<row.events<<",\"epdEventsSeen\":"<<Counter(file,"epdEventsSeen")<<",\"epdEventsValid\":"<<Counter(file,"epdEventsValid")<<",\"epdEventsRejected\":"<<Counter(file,"epdEventsRejected")<<'}';events+=row.events;epSeen+=Counter(file,"epdEventsSeen");epValid+=Counter(file,"epdEventsValid");epRejected+=Counter(file,"epdEventsRejected");if(ref){ref->Close();delete ref;}file->Close();delete file;Need(!rootError,"ROOT error during Y extraction");std::cout<<"HE4Y_SOURCE_PASS period="<<row.period<<std::endl;
  }
  sources<<']';std::ofstream out(outputPath);Need(bool(out),"cannot create raw output");out<<std::setprecision(17)<<"{\"schemaVersion\":1,\"passed\":true,\"treeEntriesRead\":0,\"sources\":"<<sources.str()<<",\"events\":"<<events<<",\"epdEventsSeen\":"<<epSeen<<",\"epdEventsValid\":"<<epValid<<",\"epdEventsRejected\":"<<epRejected<<",\"checks\":"<<checks<<",\"maxProjectionCountDifference\":"<<maxDifference<<",\"selection\":"<<Json(selection)<<",\"epDefinition\":"<<Json(epDefinition)<<",\"reducedSettings\":"<<Json(reducedDefinition)<<",\"rapidityDefinition\":"<<Json(yDefinition)<<",\"nativeEdges\":";Array(out,xedges,0,xedges.size());out<<",\"rapidityEdges\":";Array(out,yedges,0,yedges.size());out<<",\"slices\":[";bool first=true;
  for(int y=-1;y<=ny+1;++y)for(int cent=0;cent<9;++cent){if(!first)out<<',';first=false;out<<"{\"rootYBin\":"<<y<<",\"cent9\":"<<cent;for(int mode=0;mode<2;++mode){out<<','<<Json(mode?"ME":"SE")<<':';Array(out,bank.values,bank.At(mode,y,cent,0),bank.nx+2);out<<','<<Json(mode?"MEVariance":"SEVariance")<<':';Array(out,bank.variances,bank.At(mode,y,cent,0),bank.nx+2);}out<<'}';}out<<"]}\n";out.close();Need(bool(out)&&!rootError,"raw output failed");SetErrorHandler(previous);std::cout<<"HE4Y_EXTRACTION_PASS sources="<<rows.size()<<" checks="<<checks<<std::endl;
 }catch(const std::exception& ex){SetErrorHandler(previous);std::cerr<<"HE4Y_EXTRACTION_FAILED "<<ex.what()<<std::endl;gSystem->Exit(1);}
}
#endif
