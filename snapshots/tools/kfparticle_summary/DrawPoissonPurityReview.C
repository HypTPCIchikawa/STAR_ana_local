// Render saved Poisson fits and residuals, retaining empty mass bins.
#include "DrawNormalizedPurityCases.C"
namespace PoissonReview {
using namespace PurityFitReviewDraw;
typedef std::vector<PurityFitReviewDraw::Cell> Cells;
std::map<std::string,std::vector<double> > bases;
void ReadBasis(const char *path){
 std::ifstream f(path);Require(f.good(),"cannot read centered-basis TSV");std::string line;std::getline(f,line);
 Require(line=="page\tside\tmidpoint\tscale\tb0\tb1\tb2\tconstantFallback","unexpected centered-basis header");
 while(std::getline(f,line)){if(line.empty())continue;std::vector<std::string> v=Fields(line);Require(v.size()==8,"basis row requires eight columns");
  std::string key=v[0]+"/"+v[1];Require(!bases.count(key),"duplicate centered-basis row");
  std::vector<double> values;for(int i=2;i<8;++i)values.push_back(Number(v[i]));Require(values[1]>0,"nonpositive basis scale");bases[key]=values;}
}
TF1 *StableModel(Drawing &d,const Row &r,double lo,double hi,bool background){
 for(int i=0;i<6;++i)if(!Finite(r.parameters[i]))return 0;
 if(r.parameters[2]<=0)return 0;
 const std::string key=Form("%d/%s",r.page,r.side.c_str());Require(bases.count(key),"missing centered-basis parameters");
 const std::vector<double> &b=bases[key];
 const std::string poly=Form("((x-(%.17g))/(%.17g))",b[0],b[1]);
 const std::string formula=background?"[0]+[1]*"+poly+"+[2]*"+poly+"^2":"gaus(0)+[3]+[4]*"+poly+"+[5]*"+poly+"^2";
 TF1 *f=new TF1(Form("pl_stable_%d_%s_%u",r.page,r.side.c_str(),static_cast<unsigned>(d.owned.size())),formula.c_str(),lo,hi);
 if(!background)for(int i=0;i<3;++i)f->SetParameter(i,r.parameters[i]);
 for(int i=0;i<3;++i)f->SetParameter(i+(background?0:3),b[i+2]);
 f->SetLineColor(background?kGreen+2:kMagenta+1);f->SetLineStyle(background?2:1);f->SetLineWidth(2);f->SetNpx(500);d.Keep(f);return f;
}
bool In(const Cell &c,double lo,double hi){double x=(c.low+c.high)/2.;return lo<=x&&x<hi;}
bool Residual(double n,double mu,double &r,double &v){
 if(!Finite(n)||n<0||!Finite(mu)||mu<0||(n>0&&mu<=0))return false;
 v=n==0?2*mu:2*(mu-n+n*std::log(n/mu));v=std::max(0.,v);
 r=(n<mu?-1.:n>mu?1.:0.)*std::sqrt(v);return Finite(r)&&Finite(v);
}
void Summary(TF1 *f,const Cells &cells,double lo,double hi,int &bins,int &zeros,int &bad,double &n,double &deviance){
 bins=zeros=bad=0;n=deviance=0;
 for(unsigned i=0;i<cells.size();++i){if(!In(cells[i],lo,hi))continue;++bins;n+=cells[i].count;if(cells[i].count==0)++zeros;
  double r=0,v=0;if(!f||!Residual(cells[i].count,f->Eval((cells[i].low+cells[i].high)/2.),r,v)){++bad;continue;}deviance+=v;}
}
void Mass(Drawing &d,int pad,const Row &r,const Cells &cells,double lo,double hi,double slo,double shi){
 const std::string basisKey=Form("%d/%s",r.page,r.side.c_str());
 Require(bases.count(basisKey) && bases[basisKey].size()==6,"missing background model metadata");
 const bool constant=bases[basisKey][5]!=0;
 Pad(d,pad);TH1D *h=Histogram(d,cells,Form("pl_mass_%s_%d",r.side.c_str(),r.page),r.side=="SE"?kRed+1:kBlue+1);
 TF1 *f=StableModel(d,r,lo,hi,false),*b=StableModel(d,r,lo,hi,true);double ymin=0,ymax=0;
 for(unsigned i=0;i<cells.size();++i){if(!In(cells[i],lo,hi))continue;double x=(cells[i].low+cells[i].high)/2.;ymax=std::max(ymax,cells[i].count+std::sqrt(cells[i].error2));
  TF1 *fs[]={f,b};for(int j=0;j<2;++j)if(fs[j]){double y=fs[j]->Eval(x);if(Finite(y)){ymin=std::min(ymin,y);ymax=std::max(ymax,y);}}}
 ymin=ymin<0?ymin*1.1:0;ymax=ymax>0?ymax*1.15:1;
 Axes(h,"Window counts / bin",lo,hi);h->SetMinimum(ymin);h->SetMaximum(ymax);h->Draw("E1");if(f)f->Draw("same");if(b)b->Draw("same");SignalLines(d,slo,shi,ymin,ymax);
 int bins,zeros,bad;double n,D;Summary(f,cells,lo,hi,bins,zeros,bad,n,D);
 if(!bad && Finite(r.chi2))Require(std::fabs(D-r.chi2)<1.e-5+1.e-7*std::fabs(r.chi2),"plotted and saved Poisson deviance disagree");
 Text(d,.13,.963,r.label+" | "+r.method+" | "+r.side+" | L | "+(constant?"constant BG":"quadratic BG"),.035);
 Text(d,.13,.917,WindowLabel(r),.029);
 Text(d,.13,.873,"status="+Format(r.status)+", cov="+Format(r.covariance)+"; D/ndf="+(bad?"NA":Format(D)+"/"+Format(r.ndf)),.028);
 const std::string purity=Finite(r.purity)&&Finite(r.purityError)?Format(r.purity)+" #pm "+Format(r.purityError):"NA";
 Text(d,.13,.829,"P="+purity+"; N(fit)="+Format(n,7)+"; zero="+Format(zeros,8)+"/"+Format(bins,8),.028);
 Text(d,.13,.786,r.reason+"; saved curve, no refit",.022);
 TLegend *legend=new TLegend(.61,.60,.96,.73);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.028);legend->AddEntry(h,r.side.c_str(),"lep");if(f)legend->AddEntry(f,constant?"Gaussian + b (b #geq 0)":"Gaussian + pol2 (L)","l");if(b)legend->AddEntry(b,constant?"Slope terms fixed to zero":"Saved background","l");legend->Draw();d.Keep(legend);
 if(!f)Text(d,.18,.65,"Saved model: NA",.04);else if(bad)Text(d,.18,.65,"Nonphysical prediction in fit range",.03);
}
void Pull(Drawing &d,int pad,const Row &r,const Cells &cells,double lo,double hi){
 Pad(d,pad);TF1 *f=StableModel(d,r,lo,hi,false);TGraph *g=new TGraph();d.Keep(g);double extent=4.;
 int bins,zeros,bad;double n,D;Summary(f,cells,lo,hi,bins,zeros,bad,n,D);
 if(f)for(unsigned i=0;i<cells.size();++i){if(!In(cells[i],lo,hi))continue;double x=(cells[i].low+cells[i].high)/2.,res=0,v=0;
  if(!Residual(cells[i].count,f->Eval(x),res,v))continue;g->SetPoint(g->GetN(),x,res);extent=std::max(extent,std::fabs(res)*1.15);}
 TH1D *frame=new TH1D(Form("pl_pull_%s_%d",r.side.c_str(),r.page),"",1,lo,hi);frame->SetDirectory(0);d.Keep(frame);Axes(frame,"Signed Poisson deviance residual",lo,hi);frame->SetMinimum(-extent);frame->SetMaximum(extent);frame->Draw("AXIS");
 TLine *line=new TLine(lo,0,hi,0);line->SetLineStyle(2);line->SetLineColor(kGray+1);line->Draw();d.Keep(line);
 g->SetMarkerStyle(20);g->SetMarkerSize(.5);g->SetMarkerColor(r.side=="SE"?kRed+1:kBlue+1);g->Draw("P SAME");
 Text(d,.13,.963,r.side+" deviance residual | Poisson L",.038);Text(d,.13,.917,WindowLabel(r),.029);
 Text(d,.13,.873,"Empty bins retained: r(n=0) = -#sqrt{2#mu}",.029);
 Text(d,.13,.829,"Points="+Format(g->GetN(),8)+"; zero="+Format(zeros,8)+"; invalid #mu="+Format(bad,8),.029);
 Text(d,.13,.786,"D recomputed from saved model; no refit",.025);if(!f)Text(d,.2,.59,"Residuals: NA",.04);
}
void PageDraw(Drawing &d,const Page &p,double lo,double hi,double slo,double shi){
 const Row &a=p.rows.find("SE")->second,&b=p.rows.find("ME")->second;
 const Cells &se=p.bins.find("SE")->second,&me=p.bins.find("ME")->second;
 Mass(d,1,a,se,lo,hi,slo,shi);Mass(d,2,b,me,lo,hi,slo,shi);
 NormalizedPurityCases::Overlay(d,3,a,se,me,lo,hi,slo,shi,false,false);
 Pull(d,4,a,se,lo,hi);Pull(d,5,b,me,lo,hi);
 NormalizedPurityCases::Overlay(d,6,a,se,me,lo,hi,slo,shi,true,true);d.canvas->Update();
}
void Run(const char *meta,const char *bins,const char *basis,const char *directory,double lo,double hi,double slo,double shi){
 Require(meta&&bins&&directory&&*meta&&*bins&&*directory,"explicit paths required");Require(Finite(lo)&&Finite(hi)&&Finite(slo)&&Finite(shi)&&lo<hi&&lo<=slo&&slo<shi&&shi<=hi,"invalid explicit ranges");
 Pages pages;ReadMeta(meta,pages);ReadBins(bins,pages);Validate(pages,lo,hi);bases.clear();Require(basis && *basis,"explicit basis file required");ReadBasis(basis);std::string out(directory);std::map<std::string,int> counts;
 if(gSystem->AccessPathName(out.c_str()))Require(gSystem->mkdir(out.c_str(),kTRUE)==0,"cannot create output directory");
 for(Pages::const_iterator i=pages.begin();i!=pages.end();++i){const Row &r=i->second.rows.find("SE")->second;++counts[r.method];Require(gSystem->AccessPathName(SinglePath(out,r).c_str()),"refusing existing fit PDF");}
 gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(20,26);const char *methods[]={"method3","method5"};
 for(int m=0;m<2;++m){if(!counts[methods[m]])continue;std::string atlas=out+"/atlas_"+methods[m]+".pdf";Require(gSystem->AccessPathName(atlas.c_str()),"refusing existing atlas PDF");
  TCanvas book(Form("pl_book_%d",m),"Poisson fit atlas",1500,900);book.Print((atlas+"[").c_str());
  for(Pages::const_iterator i=pages.begin();i!=pages.end();++i){const Row &r=i->second.rows.find("SE")->second;if(r.method!=methods[m])continue;Drawing d(i->first);PageDraw(d,i->second,lo,hi,slo,shi);d.canvas->Print(atlas.c_str());d.canvas->Print(SinglePath(out,r).c_str());}
  book.Print((atlas+"]").c_str());std::cout<<"POISSON_ATLAS_COMPLETE "<<atlas<<" pages="<<counts[methods[m]]<<std::endl;}
}
}
void DrawPoissonPurityReview(const char *meta,const char *bins,const char *basis,const char *out,double lo,double hi,double slo,double shi){
 try{PoissonReview::Run(meta,bins,basis,out,lo,hi,slo,shi);}catch(const std::exception &e){std::cerr<<"POISSON_PURITY_REVIEW_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
