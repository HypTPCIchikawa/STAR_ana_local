#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TSystem.h"
#include <fstream>
#include <cmath>
void CheckPoissonEmptyBins(const char *out){
 TH1D h("empty_toy","",10,0,10);h.SetBinContent(5,1);h.SetBinError(5,1);
 TF1 old("old_toy","pol0",0,10);old.SetParameter(0,.3);
 TFitResultPtr c=h.Fit(&old,"SQNR0");
 TF1 poisson("poisson_toy","pol0",0,10);poisson.SetParameter(0,.3);poisson.SetParLimits(0,1.e-8,10);
 TFitResultPtr l=h.Fit(&poisson,"SLQNR0");
 const double b=poisson.GetParameter(0),d=2*(10*b-1+std::log(1./b));
 std::ofstream f(out);f.precision(17);
 f<<"{\"chiSquareEstimate\":"<<old.GetParameter(0)<<",\"poissonEstimate\":"<<b
  <<",\"poissonDeviance\":"<<d<<",\"chiSquareStatus\":"<<int(c)<<",\"poissonStatus\":"<<int(l)<<"}\n";
 f.close();if(std::fabs(old.GetParameter(0)-1)>1.e-3 || std::fabs(b-.1)>1.e-3 || int(l)!=0)gSystem->Exit(2);
}
