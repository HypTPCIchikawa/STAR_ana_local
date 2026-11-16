// Read-only Lambda-mass projections for inspecting the existing purity fits.
// The included extractor supplies the fit, covariance and acceptance rules.
#include "ExtractMassHalfCF.C"

namespace PurityFitReview {
using namespace MassHalfExtract;

void WriteMassSpectrum(std::ostream &o,const TH1D *h,bool error2) {
  std::vector<double> values;
  for(int b=1;b<=h->GetNbinsX();++b)
    values.push_back(error2?Variance(h,b):h->GetBinContent(b));
  Array(o,values);
}

void WriteReviewBin(std::ostream &o,TH2D *h,int bx,int halfWidth,
                    int originalRebin,bool isPurityGrid,
                    const Settings &c,const Partition &p) {
  // Serialize the existing row once so its single fit and all validity rules
  // remain identical to ExtractMassHalfCF. The additional projections only
  // record the distributions; they do not run another fit.
  std::ostringstream saved;
  WriteFitBin(saved,h,bx,halfWidth,originalRebin,isPurityGrid,c,p);
  const std::string row=saved.str();
  Require(!row.empty() && row[row.size()-1]=='}',"invalid serialized fit row");
  const int first=std::max(1,bx-halfWidth);
  const int last=std::min(h->GetNbinsX(),bx+halfWidth);
  TH1D *window=h->ProjectionY("pfr_window_projection",first,last,"e");
  TH1D *center=h->ProjectionY("pfr_central_projection",bx,bx,"e");
  Require(window && center,"cannot project review mass distributions");
  window->SetDirectory(0);center->SetDirectory(0);
  o<<row.substr(0,row.size()-1)<<",\"massCounts\":";
  WriteMassSpectrum(o,window,false);
  o<<",\"massError2\":";WriteMassSpectrum(o,window,true);
  o<<",\"centralMassCounts\":";WriteMassSpectrum(o,center,false);
  o<<",\"centralMassError2\":";WriteMassSpectrum(o,center,true);
  o<<'}';
  delete center;delete window;
}

void Run(const char *inputRoot,const char *outputDir,
         const char *settingsFile,double kstarMax) {
  Require(inputRoot && *inputRoot && outputDir && *outputDir &&
          settingsFile && *settingsFile,
          "inputRoot, outputDir, settingsFile and kstarMax are all required");
  Require(Finite(kstarMax) && kstarMax>0,"kstarMax must be finite and positive");
  const Settings c=ReadSettings(settingsFile);
  const std::string out(outputDir);
  const std::string target=out+"/purity_fit_review.json";
  const std::string temporary=target+".tmp";
  Require(gSystem->AccessPathName(target.c_str()) &&
          gSystem->AccessPathName(temporary.c_str()),
          "refusing to overwrite existing purity_fit_review.json or .tmp");
  if(gSystem->AccessPathName(out.c_str()))
    Require(gSystem->mkdir(out.c_str(),kTRUE)==0,"cannot create output directory");
  TFile *f=TFile::Open(inputRoot,"READ");
  Require(f && !f->IsZombie(),"cannot open input ROOT read-only");
  TH1 *inclusive=Histogram(f,"hLambda_InvMass",1);
  const Partition p=MakePartition(inclusive->GetXaxis(),c);
  Require(c.fitLo>=inclusive->GetXaxis()->GetXmin() &&
          c.fitHi<=inclusive->GetXaxis()->GetXmax(),
          "fit interval outside stored mass axis");
  const std::string referenceName=Form("true/hKstarMass_4He_CentBin%d",c.groups[0].lo);
  TH2 *reference=static_cast<TH2*>(Histogram(f,referenceName,2));
  SameAxis(inclusive->GetXaxis(),reference->GetYaxis(),referenceName);
  UniformAxis(reference->GetXaxis(),"kstar");
  Require(kstarMax>reference->GetXaxis()->GetXmin() &&
          kstarMax<=reference->GetXaxis()->GetXmax(),
          "kstarMax outside stored kstar axis");
  const std::vector<double> massEdges=RebinEdges(inclusive->GetXaxis(),1);
  const std::vector<double> originalEdges=RebinEdges(reference->GetXaxis(),1);
  const std::vector<double> fineEdges=RebinEdges(reference->GetXaxis(),c.rebin);
  const std::vector<double> coarseEdges=RebinEdges(reference->GetXaxis(),c.purityRebin);
  std::ofstream o(temporary.c_str());Require(o.good(),"cannot open review output JSON");
  o<<"{\n\"schemaVersion\":1,\"inputRoot\":"<<Quote(inputRoot)
   <<",\"settingsFile\":"<<Quote(settingsFile)<<",\"kstarMax\":";
  Num(o,kstarMax);o<<",\"settings\":";WriteSettings(o,c);
  o<<",\"massPartition\":";WritePartition(o,p,c);
  o<<",\"massEdges\":";Array(o,massEdges);
  o<<",\"originalKstarEdges\":";Array(o,originalEdges);
  o<<",\"kstarEdges\":";Array(o,fineEdges);
  o<<",\"method5PurityEdges\":";Array(o,coarseEdges);
  o<<",\"caveats\":["
   "\"Only rows with grid center x < kstarMax are exported. Sliding fit windows retain the configured neighboring bins and are clipped only at the stored histogram edges.\","
   "\"massCounts/massError2 describe the fitted window; centralMassCounts/centralMassError2 describe its central grid cell. All regular stored mass bins are retained without a display-range cut.\","
   "\"Each row reuses one WriteFitBin call and therefore one mass fit. No fit threshold, parameter limit or background rule is changed.\","
   "\"Centralities are summed after conversion to TH2D, preserving the stored Sumw2 or the existing unit-fill Poisson variance policy.\","
   "\"Observed-count/background covariance and covariance between overlapping fit windows are not supplied by the existing extractor.\"],\n\"groups\":[\n";
  const char *modes[]={"SE","ME"};
  for(unsigned gi=0;gi<c.groups.size();++gi) {
    const Group &g=c.groups[gi];
    if(gi)o<<",\n";
    std::cout<<"Purity-fit review group "<<g.id<<std::endl;
    TH2D *original[2],*fine[2],*coarse[2];
    for(int m=0;m<2;++m) {
      original[m]=GroupedMass(f,g,m==1,1,inclusive->GetXaxis(),reference->GetXaxis());
      fine[m]=static_cast<TH2D*>(original[m]->Clone(Form("pfr_fine_%s_%d",g.id.c_str(),m)));
      coarse[m]=static_cast<TH2D*>(original[m]->Clone(Form("pfr_coarse_%s_%d",g.id.c_str(),m)));
      fine[m]->SetDirectory(0);coarse[m]->SetDirectory(0);
      fine[m]->RebinX(c.rebin);coarse[m]->RebinX(c.purityRebin);
    }
    o<<"{\"id\":"<<Quote(g.id)<<",\"centMin\":"<<g.lo
     <<",\"centMax\":"<<g.hi<<",\"kstarEdges\":";Array(o,fineEdges);
    o<<",\"method5PurityEdges\":";Array(o,coarseEdges);
    o<<",\"originalKstarEdges\":";Array(o,originalEdges);
    o<<",\"method5Base\":{";
    for(int m=0;m<2;++m) {if(m)o<<',';o<<Quote(modes[m])<<':';WriteMethod5Base(o,original[m],p);}
    o<<"},\"methods\":{";
    for(int method=0;method<2;++method) {
      if(method)o<<',';
      o<<Quote(method?"method5":"method3")<<":{";
      for(int m=0;m<2;++m) {
        if(m)o<<',';
        o<<Quote(modes[m])<<":[";
        TH2D *h=method?coarse[m]:fine[m];bool first=true;
        for(int b=1;b<=h->GetNbinsX();++b) {
          if(h->GetXaxis()->GetBinCenter(b)>=kstarMax)continue;
          if(!first)o<<',';first=false;
          WriteReviewBin(o,h,b,method?c.halfWidth:0,
                         method?c.purityRebin:c.rebin,method!=0,c,p);
        }
        o<<']';
      }
      o<<'}';
    }
    for(int m=0;m<2;++m) {delete fine[m];delete coarse[m];delete original[m];}
    o<<"}}";o.flush();Require(o.good(),"review JSON write failed");
  }
  o<<"\n]}\n";o.close();Require(!o.fail(),"cannot close review output JSON");
  f->Close();delete f;
  Require(gSystem->Rename(temporary.c_str(),target.c_str())==0,
          "cannot finalize review output JSON");
  std::cout<<"EXTRACT_PURITY_FIT_REVIEW_COMPLETE "<<target<<std::endl;
}
} // namespace PurityFitReview

void ExtractPurityFitReview(const char *inputRoot,const char *outputDir,
                           const char *settingsFile,double kstarMax) {
  try { PurityFitReview::Run(inputRoot,outputDir,settingsFile,kstarMax); }
  catch(const std::exception &e) {
    std::cerr<<"ExtractPurityFitReview: "<<e.what()<<std::endl;
    gSystem->Exit(2);
  }
}
