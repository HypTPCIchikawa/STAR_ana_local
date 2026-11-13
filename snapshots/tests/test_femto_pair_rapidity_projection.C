// ROOT5 regression for the actual postprocessing projection helper.
int test_femto_pair_rapidity_projection();
#ifndef __CINT__
#include "../common/macro/checkHistAnaFemtoLambda.C"
#include "TH3D.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace PairProjectionTest {
unsigned checks=0;
void Need(bool condition,const char* message) {
  ++checks;if(!condition)throw std::runtime_error(message);
}
void Near(double actual,double expected,const char* message) {
  Need(std::fabs(actual-expected)<1.e-12*(1.+std::fabs(expected)),message);
}
TH1* Projection(TH3& source,const char* name,int firstY,int lastY,int firstCent=1,int lastCent=1) {
  TH1* result=FemtoLambdaQa::ProjectPairRapidity(source,name,firstY,lastY,firstCent,lastCent);
  Need(result&&result->GetDirectory()==0,"projection exists and belongs to its caller");
  if(result->GetNbinsX()!=source.GetNbinsX())
    std::cerr<<"PROJECTION_AXIS_MISMATCH name="<<name<<" sourceBins="<<source.GetNbinsX()
             <<" projectedBins="<<result->GetNbinsX()<<" sourceRange="<<source.GetXaxis()->GetFirst()
             <<","<<source.GetXaxis()->GetLast()<<std::endl;
  Need(result->GetNbinsX()==source.GetNbinsX(),"complete original kstar axis retained");
  for(int bin=1;bin<=source.GetNbinsX()+1;++bin)
    Near(result->GetXaxis()->GetBinLowEdge(bin),source.GetXaxis()->GetBinLowEdge(bin),"kstar bin edges retained");
  return result;
}
void Cell(const TH1& histogram,int bin,double count,double variance) {
  Near(histogram.GetBinContent(bin),count,"exact projected content");
  Near(histogram.GetBinError(bin)*histogram.GetBinError(bin),variance,"projected counting variance");
}
void WeightedAndCentrality() {
  TH3D source("weightedPairProjection","weightedPairProjection",2,0.,2.,2,0.,2.,9,-.5,8.5);
  source.SetDirectory(0);source.Sumw2();
  // The selected actual centrality index is 1, corresponding to ROOT Z bin 2.
  source.Fill(.25,-.25,1.,1.);source.Fill(.25,.25,1.,2.);
  source.Fill(.25,1.25,1.,4.);source.Fill(.25,2.25,1.,8.);
  // Distinct counts in other physical centralities and Z flow bins detect a
  // mistaken centrality index/bin conversion or the ROOT (0,0) default range.
  source.Fill(.25,-.25,0.,64.);source.Fill(.25,-.25,2.,128.);
  source.Fill(.25,-.25,-1.,256.);source.Fill(.25,-.25,9.,512.);
  source.GetXaxis()->SetRange(2,2);source.GetYaxis()->SetRange(2,2);source.GetZaxis()->SetRange(5,5);
  TH1* under=Projection(source,"qaUnderOnly",0,0);
  TH1* regular=Projection(source,"qaAdjacentRegular",1,2);
  TH1* over=Projection(source,"qaOverOnly",3,3);
  TH1* all=Projection(source,"qaAllYWithFlows",0,3);
  Cell(*under,1,1.,1.);Cell(*regular,1,6.,20.);Cell(*over,1,8.,64.);Cell(*all,1,15.,85.);
  Near(under->Integral(0,3),1.,"underflow-only call excludes all regular and overflow cells");
  Near(regular->Integral(0,3),6.,"two neighboring Y intervals merge by count");
  Near(over->Integral(0,3),8.,"overflow-only call remains distinct");
  Near(all->Integral(0,3),15.,"all Y bins plus flows close to the expected total");
  for(int x=0;x<=3;++x) {
    Near(all->GetBinContent(x),under->GetBinContent(x)+regular->GetBinContent(x)+over->GetBinContent(x),
         "Y partition count closure including X flows");
    Near(all->GetBinError(x)*all->GetBinError(x),
         under->GetBinError(x)*under->GetBinError(x)+regular->GetBinError(x)*regular->GetBinError(x)+
         over->GetBinError(x)*over->GetBinError(x),"Y partition variance closure including X flows");
  }
  TH1* first=Projection(source,"qaFirstRegular",1,1);
  TH1* second=Projection(source,"qaSecondRegular",2,2);
  Cell(*first,1,2.,4.);Cell(*second,1,4.,16.);
  TH1* widerCent=Projection(source,"qaActualCentZeroToTwo",0,3,0,2);
  Cell(*widerCent,1,207.,20565.);
  Need(source.GetXaxis()->GetFirst()==2&&source.GetXaxis()->GetLast()==2&&
       source.GetYaxis()->GetFirst()==2&&source.GetYaxis()->GetLast()==2&&
       source.GetZaxis()->GetFirst()==5&&source.GetZaxis()->GetLast()==5,
       "source display ranges are not used or changed by the explicit projection");
  Near(source.GetBinContent(source.GetBin(1,0,2)),1.,"source underflow remains unchanged");
  Near(source.GetSumw2()->At(source.GetBin(1,0,2)),1.,"source variance remains unchanged");
  delete under;delete regular;delete over;delete all;delete first;delete second;delete widerCent;
  source.Fill(-.25,.25,1.,16.);source.Fill(2.25,1.25,1.,32.);
  TH1* xflows=Projection(source,"qaKstarFlows",0,3);
  Cell(*xflows,0,16.,256.);Cell(*xflows,1,15.,85.);Cell(*xflows,3,32.,1024.);
  Near(xflows->Integral(0,3),63.,"kstar underflow and overflow are retained");
  delete xflows;
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badYNegative",-1,1,1,1),"negative Y bin rejected");
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badYOrder",2,1,1,1),"reversed Y range rejected");
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badYOverflow",1,4,1,1),"Y range outside stored flow rejected");
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badCentNegative",1,2,-1,1),"negative physical centrality rejected");
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badCentOrder",1,2,2,1),"reversed physical centralities rejected");
  Need(!FemtoLambdaQa::ProjectPairRapidity(source,"badCentOverflow",1,2,1,9),"physical centrality above eight rejected");
}
void ImplicitPoisson() {
  TH3D source("poissonPairProjection","poissonPairProjection",2,0.,2.,2,0.,2.,9,-.5,8.5);
  source.SetDirectory(0);
  source.SetBinContent(source.GetBin(1,0,2),3.);source.SetBinContent(source.GetBin(1,1,2),5.);
  source.SetBinContent(source.GetBin(1,2,2),7.);source.SetBinContent(source.GetBin(1,3,2),11.);
  Need(source.GetSumw2N()==0,"fixture represents unit-count histograms without explicit Sumw2");
  TH1* under=Projection(source,"qaPoissonUnder",0,0);
  TH1* regular=Projection(source,"qaPoissonRegular",1,2);
  TH1* over=Projection(source,"qaPoissonOver",3,3);
  TH1* all=Projection(source,"qaPoissonAll",0,3);
  Cell(*under,1,3.,3.);Cell(*regular,1,12.,12.);Cell(*over,1,11.,11.);Cell(*all,1,26.,26.);
  delete under;delete regular;delete over;delete all;
}
void VariableAxes() {
  const double x[]={0.,.1,.5,2.},y[]={0.,1.,2.},z[]={-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
  TH3D source("variablePairProjection","variablePairProjection",3,x,2,y,9,z);
  source.SetDirectory(0);source.Sumw2();source.GetXaxis()->SetTitle("kstar fixture");
  source.Fill(.05,.25,1.,3.);source.Fill(.3,-.25,1.,5.);source.Fill(1.,2.25,1.,7.);
  source.Fill(-.1,1.25,1.,11.);source.Fill(3.,.25,1.,13.);
  source.GetXaxis()->SetRange(2,2);source.GetYaxis()->SetRange(2,2);source.GetZaxis()->SetRange(1,1);
  TH1* result=Projection(source,"qaVariableXEdges",0,3);
  Need(result->GetXaxis()->GetXbins()->GetSize()==4,"variable kstar edges remain explicitly stored");
  Need(std::string(result->GetXaxis()->GetTitle())=="kstar fixture","projected kstar axis title retained");
  Cell(*result,0,11.,121.);Cell(*result,1,3.,9.);Cell(*result,2,5.,25.);
  Cell(*result,3,7.,49.);Cell(*result,4,13.,169.);
  Need(source.GetXaxis()->GetFirst()==2&&source.GetXaxis()->GetLast()==2&&
       source.GetYaxis()->GetFirst()==2&&source.GetYaxis()->GetLast()==2&&
       source.GetZaxis()->GetFirst()==1&&source.GetZaxis()->GetLast()==1,
       "variable-axis projection leaves source display ranges unchanged");
  delete result;
}
}
#endif
int test_femto_pair_rapidity_projection() {
#ifndef __CINT__
  try {
    PairProjectionTest::checks=0;
    PairProjectionTest::WeightedAndCentrality();PairProjectionTest::ImplicitPoisson();PairProjectionTest::VariableAxes();
    std::cout<<"FEMTO_PAIR_RAPIDITY_PROJECTION_PASS checks="<<PairProjectionTest::checks<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"FEMTO_PAIR_RAPIDITY_PROJECTION_FAIL "<<e.what()<<std::endl;return 1;}
#else
  return 1;
#endif
}
