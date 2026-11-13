// Plot saved fine-centrality purity estimates without rerunning mass fits.
#include "DrawMethod6Comparison.C"
namespace PoissonFine {
using namespace Method6Plot;
void Draw(const Table &data,const std::map<std::string,std::string> &labels,const std::string &method,const std::string &path,bool compare,double xlo,double xhi,double ylo,double yhi){
 gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(20,26);
 TCanvas c(Form("fine_%s_%d",method.c_str(),compare),"Fine centrality purity",1500,1050);c.Divide(3,3,.006,.01);
 for(int p=1;p<=9;++p){TPad *pad=static_cast<TPad*>(c.GetPad(p));const double x=pad->GetXlowNDC(),y=pad->GetYlowNDC(),w=pad->GetWNDC(),h=pad->GetHNDC();pad->SetPad(x,y*.95,x+w,(y+h)*.95);}
 const char *objectives[]={"poisson","chiSquare"};const char *sides[]={"SE","ME"};
 for(int j=0;j<9;++j){
  c.cd(j+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.30);gPad->SetTicks(1,1);
  const std::string group=Form("cent%d",j);if(!labels.count(group))throw std::runtime_error("Missing fine centrality: "+group);
  TH1D *frame=new TH1D(Form("fine_frame_%s_%d_%d",method.c_str(),compare,j),"",100,xlo,xhi);frame->SetDirectory(0);frame->SetMinimum(ylo);frame->SetMaximum(yhi);
  frame->GetXaxis()->SetTitle("k* (GeV/c)");frame->GetYaxis()->SetTitle("#Lambda purity");frame->GetXaxis()->SetTitleSize(.05);frame->GetYaxis()->SetTitleSize(.05);frame->GetXaxis()->SetLabelSize(.043);frame->GetYaxis()->SetLabelSize(.043);frame->GetYaxis()->SetTitleOffset(1.18);frame->Draw("AXIS");
  TLine *unity=new TLine(xlo,1,xhi,1);unity->SetLineStyle(3);unity->SetLineColor(kGray+1);unity->Draw();
  TLegend *legend=new TLegend(.14,.72,.97,.89);legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(.034);legend->SetNColumns(2);
  for(int pass=0;pass<(compare?2:1);++pass){int objective=compare?1-pass:0;
   for(int side=0;side<2;++side){std::string key=std::string(objectives[objective])+"/"+method+"/"+group+"/"+sides[side];
    Table::const_iterator it=data.find(key);if(it==data.end())throw std::runtime_error("Missing purity series: "+key);
    TGraphErrors *g=Graph(it->second,xlo,xhi,side?kBlue+1:kRed+1,objective?(side?25:24):(side?21:20));g->SetMarkerSize(.72);g->SetLineStyle(objective?2:1);g->SetLineWidth(1);g->Draw("P SAME");
    std::string name=std::string(objective?"#chi^{2} ":"L ")+sides[side]+(g->GetN()?"":" (NA)");legend->AddEntry(g,name.c_str(),"lp");
   }
  }
  legend->Draw();TLatex *title=new TLatex(.14,.95,("Centrality "+labels.find(group)->second).c_str());title->SetNDC();title->SetTextSize(.046);title->Draw();
 }
 c.cd();TLatex heading;heading.SetNDC();heading.SetTextSize(.019);heading.DrawLatex(.025,.975,(method+": "+(compare?"Poisson L vs preceding #chi^{2}":"Poisson L")+"; red SE (true), blue ME (mix); invalid estimates omitted").c_str());
 c.Print(path.c_str());
}
}
void DrawPoissonFinePurity(const char *tablePath,const char *outputDir,double xlo,double xhi,double ylo,double yhi){
 using namespace Method6Plot;
 try{
  if(!tablePath||!outputDir||!*tablePath||!*outputDir||!(xlo<xhi&&ylo<yhi))throw std::runtime_error("Explicit paths and valid plot ranges required");
  std::ifstream f(tablePath);if(!f.good())throw std::runtime_error("Cannot read fine-purity TSV");std::string line;std::getline(f,line);
  if(line!="objective\tmethod\tgroup\tlabel\tside\tbin\tlow\thigh\tvalue\terror")throw std::runtime_error("Unexpected fine-purity header");
  Table data;std::map<std::string,std::string> labels;
  while(std::getline(f,line)){
   if(line.empty())continue;const std::vector<std::string> v=Split(line);if(v.size()!=10)throw std::runtime_error("Fine-purity row must have 10 columns");
   if(v[0]!="poisson"&&v[0]!="chiSquare")throw std::runtime_error("Unknown fit objective");if(v[1]!="method3"&&v[1]!="method5")throw std::runtime_error("Unknown purity method");if(v[4]!="SE"&&v[4]!="ME")throw std::runtime_error("Unknown side");
   if(labels.count(v[2])&&labels[v[2]]!=v[3])throw std::runtime_error("Centrality label differs between series");labels[v[2]]=v[3];
   Point p;p.x=(Number(v[6])+Number(v[7]))/2.;p.y=Number(v[8]);p.e=Number(v[9]);data[v[0]+"/"+v[1]+"/"+v[2]+"/"+v[4]].push_back(p);
  }
  if(labels.size()!=9)throw std::runtime_error("Expected exactly 9 centralities");
  const std::string out(outputDir);const char *methods[]={"method3","method5"};
  for(int m=0;m<2;++m)for(int comparison=0;comparison<2;++comparison){const std::string path=out+"/purity_fine_"+methods[m]+(comparison?"_comparison.pdf":"_L.pdf");if(!gSystem->AccessPathName(path.c_str()))throw std::runtime_error("Refusing existing fine-purity PDF");}
  if(gSystem->AccessPathName(out.c_str())&&gSystem->mkdir(out.c_str(),kTRUE)!=0)throw std::runtime_error("Cannot create output directory");
  for(int m=0;m<2;++m)for(int comparison=0;comparison<2;++comparison)PoissonFine::Draw(data,labels,methods[m],out+"/purity_fine_"+methods[m]+(comparison?"_comparison.pdf":"_L.pdf"),comparison,xlo,xhi,ylo,yhi);
 }catch(const std::exception &e){std::cerr<<"POISSON_FINE_PURITY_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
