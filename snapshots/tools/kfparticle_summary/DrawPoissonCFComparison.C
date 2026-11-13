// Compare explicitly computed CFs under identical normalization and selections.
#include "DrawMethod6Comparison.C"
namespace PoissonCFCompare {
void Run(const char* path,const char* dir,double lo,double hi,double nlo,double nhi) {
 using namespace Method6Plot;
 std::ifstream f(path);if(!f.good())throw std::runtime_error("cannot read comparison TSV");
 std::string line;std::getline(f,line);
 if(line!="objective\tmethod\tgroup\tlabel\tbin\tlow\thigh\tvalue\terror")throw std::runtime_error("bad TSV header");
 Table data;std::vector<std::string> groups;std::map<std::string,std::string> labels;
 while(std::getline(f,line)){
  if(line.empty())continue;const std::vector<std::string> v=Split(line);
  if(v.size()!=9)throw std::runtime_error("nine columns required");
  if(v[0]!="poisson" && v[0]!="chiSquare")throw std::runtime_error("unknown objective");
  if(!labels.count(v[2])){groups.push_back(v[2]);labels[v[2]]=v[3];}
  if(labels[v[2]]!=v[3])throw std::runtime_error("label mismatch");
  Point p;p.x=(Number(v[5])+Number(v[6]))/2.;p.y=Number(v[7]);p.e=Number(v[8]);
  data[v[0]+"/"+v[1]+"/"+v[2]].push_back(p);
 }
 const char* required[]={"cent20_80","cent0_20","cent0_10","cent10_20","cent20_60"};
 if(groups.size()!=5)throw std::runtime_error("five groups required");
 for(unsigned i=0;i<5;++i)if(groups[i]!=required[i])throw std::runtime_error("group order differs");
 const char* methods[]={"method3","method5","method6"};const char* names[]={"M3","M5","M6 shared ME"};
 const char* objectives[]={"poisson","chiSquare"};const int colors[]={kBlack,kRed+1,kBlue+1};
 const std::string output=std::string(dir)+"/cf_objective_comparison.pdf";
 if(!gSystem->AccessPathName(output.c_str()))throw std::runtime_error("refusing existing PDF");
 gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(20,26);
 TCanvas c("poisson_cf_compare","",1500,650);c.Print((output+"[").c_str());
 for(int set=0;set<2;++set){
  c.Clear();const int first=set?2:0,n=set?3:2;c.Divide(n,1,.008,.04);
  for(int j=0;j<n;++j){
   c.cd(j+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.31);gPad->SetTicks(1,1);
   const std::string group=groups[first+j];std::vector<TGraphErrors*> graphs;std::vector<std::string> titles;double ymin=.7,ymax=1.3;
   for(int m=0;m<3;++m)for(int o=0;o<2;++o){
    const std::string key=std::string(objectives[o])+"/"+methods[m]+"/"+group;
    if(!data.count(key))throw std::runtime_error("missing series: "+key);
    TGraphErrors* g=Graph(data.find(key)->second,lo,hi,colors[m],(o?24:20)+m);g->SetMarkerSize(o?.95:.8);g->SetLineStyle(o?2:1);graphs.push_back(g);titles.push_back(std::string(names[m])+(o?" #chi^{2}":" L"));
    for(int b=0;b<g->GetN();++b){double x,y;g->GetPoint(b,x,y);ymin=std::min(ymin,y-g->GetErrorY(b));ymax=std::max(ymax,y+g->GetErrorY(b));}
   }
   const double pad=.08*(ymax-ymin);ymin-=pad;ymax+=pad;
   TH1D* frame=new TH1D(Form("poisson_cf_frame_%d_%d",set,j),"",100,lo,hi);frame->SetDirectory(0);frame->SetMinimum(ymin);frame->SetMaximum(ymax);frame->GetXaxis()->SetTitle("k* (GeV/c)");frame->GetYaxis()->SetTitle("C(k*)");frame->GetXaxis()->SetTitleSize(.05);frame->GetYaxis()->SetTitleSize(.05);frame->GetXaxis()->SetLabelSize(.043);frame->GetYaxis()->SetLabelSize(.043);frame->GetYaxis()->SetTitleOffset(1.25);frame->Draw("AXIS");
   TLine* unity=new TLine(lo,1,hi,1);unity->SetLineStyle(2);unity->Draw();
   TLegend* legend=new TLegend(.16,.72,.96,.905);legend->SetNColumns(2);legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(.030);
   for(int o=1;o>=0;--o)for(int m=0;m<3;++m)graphs[m*2+o]->Draw("P SAME");
   for(unsigned i=0;i<graphs.size();++i)legend->AddEntry(graphs[i],(titles[i]+(graphs[i]->GetN()?"":" (NA)")).c_str(),"lp");legend->Draw();
   TLatex* heading=new TLatex(.14,.95,("Centrality "+labels[group]).c_str());heading->SetNDC();heading->SetTextFont(42);heading->SetTextSize(.044);heading->Draw();
  }
  c.cd();TLatex title;title.SetNDC();title.SetTextFont(42);title.SetTextSize(.022);title.DrawLatex(.02,.976,Form("#Lambda-^{4}He CF: Poisson L (filled) versus #chi^{2} (open); common normalization %.2f-%.2f GeV/c",nlo,nhi));c.Print(output.c_str());
 }
 c.Print((output+"]").c_str());
}
}
void DrawPoissonCFComparison(const char* path,const char* dir,double lo,double hi,double nlo,double nhi){
 try{if(!path || !*path || !dir || !*dir || !(lo<hi) || !(nlo<nhi))throw std::runtime_error("explicit valid paths/ranges required");PoissonCFCompare::Run(path,dir,lo,hi,nlo,nhi);}catch(const std::exception& e){std::cerr<<"POISSON_CF_COMPARISON_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
