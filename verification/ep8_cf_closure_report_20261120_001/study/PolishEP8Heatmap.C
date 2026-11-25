void PolishEP8Heatmap(const char* tables,const char* figures);
#ifndef __CINT__
#include <TCanvas.h>
#include <TH2D.h>
#include <TBox.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>
#include <algorithm>
namespace EP8HeatmapPolish {
void Need(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
std::vector<std::string> Fields(const std::string& s){std::vector<std::string> f;size_t a=0;for(;;){size_t b=s.find('\t',a);f.push_back(s.substr(a,b==std::string::npos?b:b-a));if(b==std::string::npos)break;a=b+1;}return f;}
double Num(const std::string& s){std::istringstream in(s);double v=0;Need(bool(in>>v)&&TMath::Finite(v),"Invalid numeric field");return v;}
std::string Key(const std::string& c,int p){std::ostringstream s;s<<c<<'/'<<p;return s.str();}
struct Cell{double value;bool valid;Cell():value(0),valid(false){}};
void Note(double x,double y,const char* s,double size){TLatex* t=new TLatex(x,y,s);t->SetNDC();t->SetTextFont(42);t->SetTextSize(size);t->Draw();}
}
void PolishEP8Heatmap(const char* tables,const char* figures){using namespace EP8HeatmapPolish;try{
  Need(tables&&figures,"Explicit paths required");std::string base=tables,line;std::map<std::string,std::string> labels;std::map<std::string,Cell> cells;
  std::ifstream meta((base+"/ep8_curve_metadata.tsv").c_str());Need(bool(meta),"Missing metadata");while(std::getline(meta,line)){std::vector<std::string> f=Fields(line);Need(f.size()==7,"Bad metadata row");labels[f[0]]=f[1];}Need(meta.eof(),"Metadata read failed");
  std::ifstream in((base+"/ep8_plot_windows.tsv").c_str());Need(bool(in),"Missing window table");while(std::getline(in,line)){std::vector<std::string> f=Fields(line);Need(f.size()==9,"Bad window row");if(f[2]!="k060_100")continue;Cell c;c.valid=f[7]!="NA";if(c.valid)c.value=Num(f[7]);Need(cells.insert(std::make_pair(Key(f[0],int(Num(f[1]))),c)).second,"Duplicate window cell");}Need(in.eof(),"Window read failed");
  std::vector<std::string> rows;rows.push_back("c00_60");for(int c=8;c>=0;--c)rows.push_back(Form("cent9_%d",c));
  const std::string pdf=std::string(figures)+"/ep8_highk_primary_heatmap_v2.pdf",png=std::string(figures)+"/ep8_highk_primary_heatmap_v2.png";Need(gSystem->AccessPathName(pdf.c_str())&&gSystem->AccessPathName(png.c_str()),"Refusing existing polished figure");
  TH1::AddDirectory(kFALSE);gStyle->SetOptStat(0);gStyle->SetPalette(1);gStyle->SetPaperSize(28,21);TCanvas canvas("ep8_heatmap_polish","EP8 integrated CF",1600,1100);canvas.SetLeftMargin(.21);canvas.SetRightMargin(.15);canvas.SetTopMargin(.16);canvas.SetBottomMargin(.16);
  TH2D h("ep8_heatmap_readable","",8,-.5,7.5,10,0,10);h.SetDirectory(0);double low=.85,high=1.05;int valid=0,missing=0,zeros=0;
  for(size_t r=0;r<rows.size();++r){h.GetYaxis()->SetBinLabel(10-r,labels[rows[r]].c_str());for(int p=0;p<8;++p){std::string k=Key(rows[r],p);Need(cells.count(k),"Missing requested cell");const Cell& c=cells[k];if(c.valid){h.SetBinContent(p+1,10-r,c.value);low=std::min(low,c.value);high=std::max(high,c.value);++valid;if(c.value==0)++zeros;}else ++missing;}}
  h.SetMinimum(low-.015);h.SetMaximum(high+.015);h.GetXaxis()->SetTitle("EP bin");h.GetXaxis()->SetTitleSize(.038);h.GetXaxis()->SetLabelSize(.035);h.GetYaxis()->SetTitle("");h.GetYaxis()->SetLabelSize(.033);h.GetZaxis()->SetTitle("Integrated raw-signal CF");h.GetZaxis()->SetTitleSize(.032);h.GetZaxis()->SetLabelSize(.030);h.GetZaxis()->SetTitleOffset(1.05);h.Draw("COLZ");
  for(size_t r=0;r<rows.size();++r)for(int p=0;p<8;++p){const Cell& c=cells[Key(rows[r],p)];double y=9.5-r;if(!c.valid){TBox* box=new TBox(p-.495,y-.495,p+.495,y+.495);box->SetFillColor(kWhite);box->SetLineColor(kGray);box->Draw();}TLatex* t=new TLatex(p,y,c.valid?Form("%.3f",c.value):"NA");t->SetTextFont(42);t->SetTextAlign(22);t->SetTextSize(.033);t->SetTextColor(c.valid&&c.value<.4?kWhite:kBlack);t->Draw();}
  Note(.08,.955,"EP8 cut Off: window-integrated C(k*) for 0.6 #leq k* < 1.0 GeV/c",.031);Note(.08,.913,"C_{W}=SE(W)/[#alpha ME(W)],  #alpha=SE(H)/ME(H),  H=[0.4,0.5) GeV/c",.025);
  Note(.08,.077,"NA: no normalization-window SE/ME; 0.000: valid H normalization but zero SE in the high-k window.",.020);Note(.08,.038,"Peripheral bins have small H counts. Shared normalization/event correlations are not a significance test.",.020);
  canvas.Print(pdf.c_str());canvas.Print(png.c_str());std::cout<<"EP8_HEATMAP_POLISH_PASS valid="<<valid<<" missing="<<missing<<" zero="<<zeros<<std::endl;
}catch(const std::exception& e){std::cerr<<"EP8_HEATMAP_POLISH_FAILED "<<e.what()<<std::endl;gSystem->Exit(1);}}
#endif
