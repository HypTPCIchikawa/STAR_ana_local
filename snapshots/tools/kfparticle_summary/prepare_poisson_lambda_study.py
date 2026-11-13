#!/usr/bin/env python3
"""Create isolated, explicitly configured Poisson-L copies of the mass extractor."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_once(text, old, new):
    if text.count(old)!=1:
        raise ValueError('Expected one source match: '+old[:80])
    return text.replace(old,new)


FIT = r'''FitDiagnostics FitMass(TH1D *h,const Settings &c) {
  FitDiagnostics d;d.fitBins=0;d.emptyBins=0;d.minExpected=0;
  d.rootStatistic=0;d.expectationValid=false;
  if(h->Integral()<c.minEntries){d.reason="insufficient_entries";return d;}
  TF1 fit("mh_total_fit","gaus(0)+pol2(3)",c.fitLo,c.fitHi);
  fit.SetParameters(h->GetMaximum(),c.mean,c.sigmaInitial,0.,0.,0.);
  fit.SetParLimits(1,c.meanLo,c.meanHi);fit.SetParLimits(2,c.sigmaLo,c.sigmaHi);
  TFitResultPtr result=h->Fit(&fit,c.fitOptions.c_str(),"",c.fitLo,c.fitHi);
  d.status=int(result);d.covStatus=result.Get()?result->CovMatrixStatus():-1;
  d.isValid=result.Get() && result->IsValid();d.boundary=Boundary(fit);
  d.rootStatistic=fit.GetChisquare();d.expectationValid=true;
  double deviance=0;d.minExpected=std::numeric_limits<double>::infinity();
  for(int b=1;b<=h->GetNbinsX();++b){
    const double x=h->GetBinCenter(b);if(x<c.fitLo || x>c.fitHi)continue;
    ++d.fitBins;const double n=h->GetBinContent(b),mu=fit.Eval(x);
    if(n==0)++d.emptyBins;
    if(Finite(mu))d.minExpected=std::min(d.minExpected,mu);
    if(!Finite(mu) || mu<0 || (n>0 && mu<=0)){d.expectationValid=false;continue;}
    const double term=n>0?mu-n+n*std::log(n/mu):mu;
    deviance+=2*std::max(0.,term);
  }
  d.chi2=d.expectationValid?deviance:std::numeric_limits<double>::quiet_NaN();
  d.ndf=std::max(0,d.fitBins-fit.GetNpar());
  d.probability=d.ndf>0 && Finite(d.chi2)?TMath::Prob(d.chi2,d.ndf):0;
  bool finite=true;for(int i=0;i<6;++i){d.parameters[i]=fit.GetParameter(i);finite=finite && Finite(d.parameters[i]);}
  if(result.Get()){
    TMatrixDSym cov=result->GetCovarianceMatrix();d.finiteCov=cov.GetNrows()==6;
    for(int i=0;i<6 && cov.GetNrows()==6;++i)for(int j=0;j<6;++j){
      d.covariance[i][j]=cov(i,j);d.finiteCov=d.finiteCov && Finite(cov(i,j));
    }
  }
  d.numericallyValid=d.status==0 && d.isValid && d.covStatus>=c.minCov && finite && d.finiteCov && d.expectationValid;
  d.reason=!d.expectationValid?"nonphysical_poisson_expectation":d.numericallyValid?
    (d.boundary?"numerically_valid_boundary_warning":"numerically_valid"):"invalid_fit_or_covariance";
  return d;
}
void WriteFit(std::ostream &o,const FitDiagnostics &d,double observations){
  o<<"{\"status\":"<<d.status<<",\"covarianceStatus\":"<<d.covStatus<<",\"isValid\":"<<Bool(d.isValid)
   <<",\"numericallyValid\":"<<Bool(d.numericallyValid)<<",\"boundary\":"<<Bool(d.boundary)
   <<",\"reason\":"<<Quote(d.reason)<<",\"chi2\":";
  if(d.status!=-1)Num(o,d.chi2);else o<<"null";
  o<<",\"ndf\":"<<d.ndf<<",\"probability\":";if(d.status!=-1)Num(o,d.probability);else o<<"null";
  o<<",\"observationsInStoredMassRange\":";Num(o,observations);o<<",\"parameters\":[";
  for(int i=0;i<6;++i){if(i)o<<',';if(d.status!=-1)Num(o,d.parameters[i]);else o<<"null";}
  o<<"],\"covariance\":[";
  for(int i=0;i<6;++i){if(i)o<<',';o<<'[';for(int j=0;j<6;++j){if(j)o<<',';
    if(d.finiteCov)Num(o,d.covariance[i][j]);else o<<"null";}o<<']';}
  o<<"],\"objective\":\"poisson\",\"fitOptions\":\"SLQNR0\",\"deviance\":";
  if(d.status!=-1)Num(o,d.chi2);else o<<"null";
  o<<",\"rootReportedStatistic\":";if(d.status!=-1)Num(o,d.rootStatistic);else o<<"null";
  o<<",\"usedFitBins\":"<<d.fitBins<<",\"emptyFitBins\":"<<d.emptyBins
   <<",\"expectationValid\":"<<Bool(d.expectationValid)<<",\"minimumExpectedCount\":";
  if(d.status!=-1)Num(o,d.minExpected);else o<<"null";
  o<<",\"probabilityCaveat\":\"Asymptotic reference only; not a calibrated sparse-count goodness-of-fit p-value.\"}";
}
'''


TOY = r'''#include "TH1D.h"
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
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--group-reference',type=Path,required=True)
    parser.add_argument('--fine-reference',type=Path,required=True)
    parser.add_argument('--mainconf',type=Path,required=True)
    parser.add_argument('--target-report',type=Path,required=True)
    parser.add_argument('--center-polynomial',action='store_true',
                        help='Fit an exactly equivalent centered/scaled quadratic basis.')
    parser.add_argument('--background-seed-floor',type=float)
    parser.add_argument('--maximum-function-calls',type=int)
    parser.add_argument('--negative-background-constant',action='store_true')
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[2];source=repo/'tools/kfparticle_summary'
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'poisson_context.json').exists():raise FileExistsError('Existing prepared study')
    macros=out/'macros';macros.mkdir(exist_ok=True)
    original=(source/'ExtractMassHalfCF.C').read_text();code=original
    code=replace_once(code,'#include <algorithm>','#include <algorithm>\n#include <limits>')
    code=replace_once(code,'int rebin,purityRebin,halfWidth,minCov;','int rebin,purityRebin,halfWidth,minCov;\n  std::string fitOptions;')
    code=replace_once(code,'else if(key=="fitRange")','else if(key=="fitOptions") s>>c.fitOptions;\n    else if(key=="fitRange")')
    code=replace_once(code,'"fitMeanLimits","fitSigmaLimits","fitSigmaInitial","minimumEntries","minimumCovarianceStatus"};',
                          '"fitMeanLimits","fitSigmaLimits","fitSigmaInitial","minimumEntries","minimumCovarianceStatus","fitOptions"};')
    code=replace_once(code,'Require(!c.groups.empty(),"no explicit centrality groups");',
      'Require(c.fitOptions=="SLQNR0","Poisson study requires explicit SLQNR0");\n  Require(!c.groups.empty(),"no explicit centrality groups");')
    code=replace_once(code,'double chi2,probability,parameters[6],covariance[6][6];std::string reason;',
      'double chi2,probability,parameters[6],covariance[6][6];std::string reason;\n  int fitBins,emptyBins;double rootStatistic,minExpected;bool expectationValid;')
    a=code.index('FitDiagnostics FitMass(');b=code.index('void WriteFitBin(',a)
    code=code[:a]+FIT+'\n'+code[b:]
    if args.center_polynomial:
        code=replace_once(code,'int fitBins,emptyBins;double rootStatistic,minExpected;bool expectationValid;',
          'int fitBins,emptyBins;double rootStatistic,minExpected;bool expectationValid;\n  double basisParameters[3],basisCovariance[3][3],basisMidpoint,basisScale;')
        code=replace_once(code,'FitDiagnostics d;d.fitBins=0;d.emptyBins=0;d.minExpected=0;',
          'FitDiagnostics d;d.fitBins=0;d.emptyBins=0;d.minExpected=0;\n  d.basisMidpoint=(c.fitLo+c.fitHi)/2.;d.basisScale=(c.fitHi-c.fitLo)/2.;\n  for(int a=0;a<3;++a){d.basisParameters[a]=0;for(int b=0;b<3;++b)d.basisCovariance[a][b]=0;}')
        code=replace_once(code,'TF1 fit("mh_total_fit","gaus(0)+pol2(3)",c.fitLo,c.fitHi);',r'''const double midpoint=(c.fitLo+c.fitHi)/2.,scale=(c.fitHi-c.fitLo)/2.;
  TF1 fit("mh_total_fit",Form("gaus(0)+[3]+[4]*((x-(%.17g))/(%.17g))+[5]*((x-(%.17g))/(%.17g))^2",midpoint,scale,midpoint,scale),c.fitLo,c.fitHi);''')
        marker='  d.numericallyValid=d.status==0'
        conversion=r'''  // This invertible change of coordinates preserves the full quadratic family.
  const double original[3]={fit.GetParameter(3),fit.GetParameter(4),fit.GetParameter(5)};
  for(int a=0;a<3;++a){d.basisParameters[a]=original[a];for(int b=0;b<3;++b)d.basisCovariance[a][b]=d.covariance[a+3][b+3];}
  double jacobian[6][6]={{0}};
  for(int i=0;i<3;++i)jacobian[i][i]=1.;
  jacobian[3][3]=1.;jacobian[3][4]=-midpoint/scale;jacobian[3][5]=midpoint*midpoint/(scale*scale);
  jacobian[4][4]=1./scale;jacobian[4][5]=-2*midpoint/(scale*scale);jacobian[5][5]=1./(scale*scale);
  d.parameters[3]=original[0]-midpoint/scale*original[1]+midpoint*midpoint/(scale*scale)*original[2];
  d.parameters[4]=original[1]/scale-2*midpoint/(scale*scale)*original[2];
  d.parameters[5]=original[2]/(scale*scale);
  if(d.finiteCov){
    double transformed[6][6];
    for(int i=0;i<6;++i)for(int j=0;j<6;++j){
      long double value=0;
      for(int a=0;a<6;++a)for(int b=0;b<6;++b)
        value+=static_cast<long double>(jacobian[i][a])*d.covariance[a][b]*jacobian[j][b];
      transformed[i][j]=static_cast<double>(value);
    }
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)d.covariance[i][j]=transformed[i][j];
  }
'''
        code=replace_once(code,marker,conversion+marker)
        code=replace_once(code,'gradient[r][1]=gradient[r][0]*(hi+lo)/2.L;\n    gradient[r][2]=gradient[r][0]*(hi*hi+hi*lo+lo*lo)/3.L;',
          'const long double a=(lo-fit.basisMidpoint)/fit.basisScale,b=(hi-fit.basisMidpoint)/fit.basisScale;\n    gradient[r][1]=gradient[r][0]*(a+b)/2.L;\n    gradient[r][2]=gradient[r][0]*(a*a+a*b+b*b)/3.L;')
        code=replace_once(code,'static_cast<long double>(fit.parameters[a+3])','static_cast<long double>(fit.basisParameters[a])')
        code=replace_once(code,'static_cast<long double>(fit.covariance[a+3][b+3])','static_cast<long double>(fit.basisCovariance[a][b])')
        endpoint=code.index('probabilityCaveat')
        start=code.rfind('  o<<',0,endpoint)
        prefix=r'''  o<<",\"polynomialBasis\":\"centered_scaled\",\"basisMidpoint\":";Num(o,d.basisMidpoint);
  o<<",\"basisScale\":";Num(o,d.basisScale);o<<",\"basisParameters\":[";
  for(int a=0;a<3;++a){if(a)o<<',';Num(o,d.basisParameters[a]);}
  o<<"],\"basisCovariance\":[";
  for(int a=0;a<3;++a){if(a)o<<',';o<<'[';for(int b=0;b<3;++b){if(b)o<<',';Num(o,d.basisCovariance[a][b]);}o<<']';}
  o<<']';
'''
        code=code[:start]+prefix+code[start:]
    seeded=args.background_seed_floor is not None or args.maximum_function_calls is not None
    if seeded:
        if not args.center_polynomial or args.background_seed_floor is None or args.background_seed_floor<=0 or not args.maximum_function_calls or args.maximum_function_calls<=0:
            raise ValueError('Positive initialization needs centered basis, positive seed floor and explicit call budget')
        code=replace_once(code,'#include <limits>','#include <limits>\n#include "Math/MinimizerOptions.h"')
        code=replace_once(code,'std::string fitOptions;','std::string fitOptions;\n  double backgroundSeedFloor;int maximumFunctionCalls;')
        code=replace_once(code,'else if(key=="fitOptions")',
            'else if(key=="backgroundSeedFloor") s>>c.backgroundSeedFloor;\n    else if(key=="maximumFunctionCalls") s>>c.maximumFunctionCalls;\n    else if(key=="fitOptions")')
        code=replace_once(code,'"minimumCovarianceStatus","fitOptions"};',
            '"minimumCovarianceStatus","fitOptions","backgroundSeedFloor","maximumFunctionCalls"};')
        code=replace_once(code,'fit.SetParameters(h->GetMaximum(),c.mean,c.sigmaInitial,0.,0.,0.);',r'''double sidebandCount=0;int sidebandBins=0;
  for(int b=1;b<=h->GetNbinsX();++b){const double x=h->GetBinCenter(b);
    if(x>=c.fitLo && x<=c.fitHi && std::fabs(x-c.mean)>c.nsigma*c.sigma){sidebandCount+=h->GetBinContent(b);++sidebandBins;}}
  const double backgroundSeed=std::max(c.backgroundSeedFloor,sidebandBins?sidebandCount/sidebandBins:0.);
  fit.SetParameters(h->GetMaximum(),c.mean,c.sigmaInitial,backgroundSeed,0.,0.);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(c.maximumFunctionCalls);''')
    if args.negative_background_constant:
        if not seeded:raise ValueError('Constant fallback requires explicit positive seed and call budget')
        code=replace_once(code,'std::string fitOptions;', 'std::string fitOptions;\n  int constantFallbackPolicy;')
        code=replace_once(code,'else if(key=="fitOptions")',
          'else if(key=="constantFallbackPolicy") s>>c.constantFallbackPolicy;\n    else if(key=="fitOptions")')
        code=replace_once(code,'"maximumFunctionCalls"};','"maximumFunctionCalls","constantFallbackPolicy"};')
        code=replace_once(code,'Require(c.fitOptions=="SLQNR0",',
          'Require(c.constantFallbackPolicy==1,"constant fallback policy must be explicit");\n  Require(c.fitOptions=="SLQNR0",')
        code=replace_once(code,'bool expectationValid;',
          'bool expectationValid;\n  bool constantFallback;double primaryBackgroundMinimum,backgroundMinimum,constantUpperBound;int primaryStatus,freeParameters;')
        code=replace_once(code,'d.rootStatistic=0;d.expectationValid=false;',
          'd.rootStatistic=0;d.expectationValid=false;d.constantFallback=false;d.primaryBackgroundMinimum=0;d.backgroundMinimum=0;d.constantUpperBound=0;d.primaryStatus=-1;d.freeParameters=6;')
        code=replace_once(code,'  d.status=int(result);',r'''  d.primaryStatus=int(result);
  const double b0=fit.GetParameter(3),b1=fit.GetParameter(4),b2=fit.GetParameter(5);
  d.primaryBackgroundMinimum=std::min(b0-b1+b2,b0+b1+b2);
  if(b2>0){const double u=-b1/(2*b2);if(u>-1 && u<1)d.primaryBackgroundMinimum=std::min(d.primaryBackgroundMinimum,b0+b1*u+b2*u*u);}
  if(Finite(d.primaryBackgroundMinimum) && d.primaryBackgroundMinimum<0){
    d.constantFallback=true;d.freeParameters=4;
    double totalFit=0;for(int b=1;b<=h->GetNbinsX();++b){const double x=h->GetBinCenter(b);if(x>=c.fitLo && x<=c.fitHi)totalFit+=h->GetBinContent(b);}
    d.constantUpperBound=std::max(1.,totalFit);
    fit.FixParameter(4,0.);fit.FixParameter(5,0.);
    fit.SetParameter(3,std::min(backgroundSeed,d.constantUpperBound));
    fit.SetParLimits(3,0.,d.constantUpperBound);
    result=h->Fit(&fit,c.fitOptions.c_str(),"",c.fitLo,c.fitHi);
  }
  d.backgroundMinimum=d.constantFallback?fit.GetParameter(3):d.primaryBackgroundMinimum;
  d.status=int(result);''')
        code=replace_once(code,'d.ndf=std::max(0,d.fitBins-fit.GetNpar());','d.ndf=std::max(0,d.fitBins-d.freeParameters);')
        endpoint=code.index('probabilityCaveat');start=code.rfind('  o<<',0,endpoint)
        extra=r'''  o<<",\"constantBackgroundFallback\":"<<Bool(d.constantFallback)
   <<",\"backgroundModel\":"<<Quote(d.constantFallback?"constant_nonnegative":"quadratic")
   <<",\"primaryStatus\":"<<d.primaryStatus<<",\"primaryBackgroundMinimum\":";Num(o,d.primaryBackgroundMinimum);
  o<<",\"backgroundMinimum\":";Num(o,d.backgroundMinimum);
  o<<",\"constantUpperBound\":";Num(o,d.constantUpperBound);o<<",\"freeParameters\":"<<d.freeParameters;
'''
        code=code[:start]+extra+code[start:]
    (macros/'ExtractMassHalfCF.C').write_text(code)
    review=(source/'ExtractPurityFitReview.C').read_text()
    review=replace_once(review,
       'Each row reuses one WriteFitBin call and therefore one mass fit. No fit threshold, parameter limit or background rule is changed.',
       'Each spectrum cell records a primary fit and any configured refit. Model selection and acceptance diagnostics accompany the saved fit.')
    review=replace_once(review,'void ExtractPurityFitReview(', 'void ExtractPurityFitReviewPoisson(')
    review=replace_once(review,'WriteSettings(o,c);',
       'WriteSettings(o,c);o<<",\\\"fitOptions\\\":"<<Quote(c.fitOptions)<<",\\\"fitObjective\\\":\\\"poisson\\\"";')
    (macros/'ExtractPurityFitReviewPoisson.C').write_text(review)
    (macros/'CheckPoissonEmptyBins.C').write_text(TOY)
    reference_paths=dict(group=args.group_reference.resolve(strict=True),fine=args.fine_reference.resolve(strict=True))
    refs={k:json.loads(p.read_text()) for k,p in reference_paths.items()}
    root=Path(refs['group']['inputRoot']).resolve(strict=True)
    if Path(refs['fine']['inputRoot']).resolve()!=root:raise ValueError('Reference sample mismatch')
    commands={}
    for name,reference in refs.items():
        folder=out/(name+'_extract');folder.mkdir(exist_ok=True)
        text=Path(reference['settingsFile']).read_text()
        settings=folder/'extract_settings.txt'
        extra='\nfitOptions SLQNR0\n'
        if seeded:extra+='backgroundSeedFloor '+str(args.background_seed_floor)+'\nmaximumFunctionCalls '+str(args.maximum_function_calls)+'\n'
        if args.negative_background_constant:extra+='constantFallbackPolicy 1\n'
        settings.write_text(text+extra)
        values=[str(root),str(folder),str(settings)]
        call=str(macros/'ExtractPurityFitReviewPoisson.C')+'+('+','.join(json.dumps(x) for x in values)+',1.0)'
        commands[name]=['bash','script/singularity_root_batch.sh',str(args.mainconf.resolve()),'root4star','-l','-b','-q',call]
    toy=str(macros/'CheckPoissonEmptyBins.C')+'+('+json.dumps(str(out/'empty_bin_toy.json'))+')'
    commands['toy']=['bash','script/singularity_root_batch.sh',str(args.mainconf.resolve()),'root4star','-l','-b','-q',toy]
    target=args.target_report.resolve(strict=True)
    context=dict(phase='prepared',preparedAtJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        fitOptions='SLQNR0',objective='poisson',input=dict(path=str(root),size=root.stat().st_size,mtime_ns=root.stat().st_mtime_ns),
        polynomialBasis='centered_scaled' if args.center_polynomial else 'raw_mass',
        backgroundSeedFloor=args.background_seed_floor,maximumFunctionCalls=args.maximum_function_calls,
        negativeBackgroundConstant=args.negative_background_constant,
        references={k:dict(path=str(p),sha256=sha(p)) for k,p in reference_paths.items()},
        targetReport=dict(path=str(target),sha256=sha(target)),mainconf=str(args.mainconf.resolve()),commands=commands,
        sourceCopies={str(p):sha(p) for p in macros.glob('*.C')},
        changes=['Fit objective L; function, Gaussian seed/bounds, count gate and covariance threshold retained.',
                 'Background starts at the observed sideband mean with the explicit floor; optimizer call budget is explicit.' if seeded else 'The background coefficient seed is zero, as in the chi-square reference.',
                 'A negative quadratic background anywhere in the closed fit range triggers a new L fit with both slope terms fixed to zero and a nonnegative constant; the data-count upper guard is recorded.' if args.negative_background_constant else 'No constant-background fallback is used.',
                 'Poisson-domain validation requires nonnegative expected counts, strictly positive where observed counts are nonzero.',
                 'Deviance is computed independently including empty bins; probability is an asymptotic reference only.'])
    (out/'poisson_context.json').write_text(json.dumps(context,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),commands=list(commands))))


if __name__=='__main__':main()
