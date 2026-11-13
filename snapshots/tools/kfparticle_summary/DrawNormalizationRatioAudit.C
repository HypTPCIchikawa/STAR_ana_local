// Draw saved count-ratio diagnostics without opening or changing the input ROOT.
#include "DrawMethod6Comparison.C"
#include "TBox.h"

void DrawNormalizationRatioAudit(const char *tablePath,const char *outputPath,
                                double xlo,double xhi,double normlo,double normhi,
                                double oldlo,double oldhi) {
  using namespace Method6Plot;
  try {
    if(!(tablePath && outputPath && xlo<xhi && normlo<normhi && oldlo<oldhi))
      throw std::runtime_error("Explicit valid paths and ranges required");
    if(!gSystem->AccessPathName(outputPath))throw std::runtime_error("Refusing existing ratio PDF");
    std::ifstream f(tablePath);if(!f.good())throw std::runtime_error("Cannot read ratio TSV");
    std::string line;std::getline(f,line);
    if(line!="group\tbin\tlow\thigh\tSE\tME\tratio\terror\tcf\tcfError")
      throw std::runtime_error("Unexpected ratio TSV header");
    Table data;std::vector<std::string> groups;
    while(std::getline(f,line)) {
      if(line.empty())continue;const std::vector<std::string> c=Split(line);
      if(c.size()!=10)throw std::runtime_error("Ratio TSV column count differs");
      if(!data.count(c[0]))groups.push_back(c[0]);
      Point p;p.x=(Number(c[2])+Number(c[3]))/2.;p.y=Number(c[8]);p.e=Number(c[9]);
      data[c[0]].push_back(p);
    }
    if(groups.size()!=5)throw std::runtime_error("Expected five centrality groups");
    const char *labels[]={"20-80%","0-20%","0-10%","10-20%","20-60%"};
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPaperSize(20,26);
    TCanvas c("norm_ratio_audit","",1500,900);c.Divide(3,2,.01,.015);
    std::vector<TObject*> owned;
    for(unsigned j=0;j<groups.size();++j) {
      c.cd(j+1);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.035);
      gPad->SetBottomMargin(.16);gPad->SetTopMargin(.17);gPad->SetTicks(1,1);
      TGraphErrors *g=Graph(data[groups[j]],xlo,xhi,kBlack,20);owned.push_back(g);
      double ymin=.8,ymax=1.15;
      for(int b=0;b<g->GetN();++b){double x,y;g->GetPoint(b,x,y);ymin=std::min(ymin,y-g->GetErrorY(b));ymax=std::max(ymax,y+g->GetErrorY(b));}
      const double pad=.05*(ymax-ymin);ymin-=pad;ymax+=pad;
      TH1D *frame=new TH1D(Form("norm_ratio_frame_%u",j),"",100,xlo,xhi);owned.push_back(frame);frame->SetDirectory(0);
      frame->SetMinimum(ymin);frame->SetMaximum(ymax);
      frame->GetXaxis()->SetTitle("k* (GeV/c)");frame->GetYaxis()->SetTitle("Raw window SE/ME, normalized");
      frame->GetXaxis()->SetTitleSize(.046);frame->GetYaxis()->SetTitleSize(.046);
      frame->GetXaxis()->SetLabelSize(.04);frame->GetYaxis()->SetLabelSize(.04);
      frame->GetYaxis()->SetTitleOffset(1.5);frame->Draw("AXIS");
      TBox *box=new TBox(normlo,ymin,normhi,ymax);owned.push_back(box);
      box->SetFillColor(kGreen-9);box->SetFillStyle(3001);box->SetLineColor(kGreen+2);box->Draw();
      TLine *unity=new TLine(xlo,1,xhi,1);owned.push_back(unity);unity->SetLineStyle(2);unity->Draw();
      const double old[]={oldlo,oldhi};
      for(int k=0;k<2;++k){TLine *l=new TLine(old[k],ymin,old[k],ymax);owned.push_back(l);l->SetLineColor(kGray+2);l->SetLineStyle(3);l->Draw();}
      g->Draw("P SAME");
      TLatex *t=new TLatex(.15,.91,Form("Centrality %s",labels[j]));owned.push_back(t);t->SetNDC();t->SetTextSize(.047);t->Draw();
    }
    c.cd(6);TLatex t;t.SetNDC();t.SetTextSize(.039);
    t.DrawLatex(.08,.84,"Raw counts in the fixed #Lambda mass gate");
    t.DrawLatex(.08,.74,Form("Green: selected norm %.2f-%.2f GeV/c",normlo,normhi));
    t.DrawLatex(.08,.64,Form("Dotted: previous norm %.2f-%.2f",oldlo,oldhi));
    t.DrawLatex(.08,.54,"The raw points do not use fitted purity.");
    t.DrawLatex(.08,.44,"A normalization change does not remove a trend.");
    t.DrawLatex(.08,.34,"Point errors include count/norm covariance;");
    t.DrawLatex(.08,.25,"event and mixed-pair correlations are absent.");
    c.Print(outputPath);
    // Canvas pads own their drawing lists for the duration of this batch process.
  }catch(const std::exception &e){std::cerr<<"NORMALIZATION_RATIO_AUDIT_FAILED: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
