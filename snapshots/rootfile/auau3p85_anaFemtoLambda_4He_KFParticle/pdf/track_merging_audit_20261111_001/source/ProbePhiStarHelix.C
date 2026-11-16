#include "/star/nfs4/AFS/star/packages/SL24y/StRoot/StarClassLibrary/StPhysicalHelixD.hh"
#include "/star/nfs4/AFS/star/packages/SL24y/StRoot/StarClassLibrary/SystemOfUnits.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <utility>

static double WrapPhiHelix(double x) { return std::atan2(std::sin(x),std::cos(x)); }
void ProbePhiStarHelix(const char* settings,const char* output) {
 try {
  std::ifstream in(settings);if(!in)throw std::runtime_error("Cannot read explicit settings");
  std::ofstream out(output);if(!out)throw std::runtime_error("Cannot open output");
  out<<std::setprecision(17);
  out<<"fieldT\tcharge\tptGeV\tphi0\tpzGeV\tradiusM\tspatialPhi\tplusAsin03\tminusAsin03\tminusAsinExactC\tplusErrorRad\tminusErrorRad\texactErrorRad\tradiusErrorM\ttangentPhi\n";
  std::string line;std::getline(in,line);unsigned n=0;double maxPlus=0,maxMinus=0,maxExact=0,maxR=0;
  while(std::getline(in,line)) {
   if(line.empty())continue;
   double field,q,pt,phi0,pz,r;std::istringstream s(line);
   if(!(s>>field>>q>>pt>>phi0>>pz>>r))throw std::runtime_error("Invalid settings row");
   StPhysicalHelixD h(StThreeVectorD(pt*std::cos(phi0)*GeV,pt*std::sin(phi0)*GeV,pz*GeV),StThreeVectorD(0,0,0),field*tesla,q);
   std::pair<double,double> paths=h.pathLength(r*meter);
   double path=paths.first>=0?paths.first:paths.second;
   if(paths.first>=0 && paths.second>=0)path=std::min(paths.first,paths.second);
   if(path<0||path>1.e20)throw std::runtime_error("No outgoing cylinder crossing");
   const StThreeVectorD position=h.at(path),momentum=h.momentumAt(path,field*tesla);
   const double actual=std::atan2(position.y(),position.x());
   const double arg=.3*q*field*r/(2*pt),exactArg=.299792458*q*field*r/(2*pt);
   if(std::fabs(arg)>=1)throw std::runtime_error("Synthetic track cannot reach radius");
   const double plus=WrapPhiHelix(phi0+std::asin(arg)),minus=WrapPhiHelix(phi0-std::asin(arg));
   const double exact=WrapPhiHelix(phi0-std::asin(exactArg));
   const double ep=std::fabs(WrapPhiHelix(plus-actual)),em=std::fabs(WrapPhiHelix(minus-actual)),ee=std::fabs(WrapPhiHelix(exact-actual));
   const double er=std::fabs(position.perp()/meter-r);
   maxPlus=std::max(maxPlus,ep);maxMinus=std::max(maxMinus,em);maxExact=std::max(maxExact,ee);maxR=std::max(maxR,er);
   out<<field<<'\t'<<q<<'\t'<<pt<<'\t'<<phi0<<'\t'<<pz<<'\t'<<r<<'\t'<<actual<<'\t'<<plus<<'\t'<<minus<<'\t'<<exact<<'\t'<<ep<<'\t'<<em<<'\t'<<ee<<'\t'<<er<<'\t'<<std::atan2(momentum.y(),momentum.x())<<'\n';++n;
  }
  out.close();if(!out)throw std::runtime_error("Write failure");
  std::cout<<std::setprecision(17)<<"HELIX_PHI_TEST rows="<<n<<" maxPlusError="<<maxPlus<<" maxMinusError="<<maxMinus<<" maxExactError="<<maxExact<<" maxRadiusErrorM="<<maxR<<std::endl;
  if(n==0||maxExact>1.e-10||maxR>1.e-9)throw std::runtime_error("Helix oracle consistency failure");
 }catch(const std::exception& e){std::cerr<<"HELIX_PHI_TEST_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
