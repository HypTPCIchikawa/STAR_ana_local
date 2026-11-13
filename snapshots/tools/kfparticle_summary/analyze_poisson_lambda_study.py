#!/usr/bin/env python3
"""Compare matching saved chi-square/Poisson fits and prepare fixed diagnostics."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

import yaml

from select_purity_diagnostic_cases import candidates, write_fit_tables


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def write_basis(cases,folder):
    lines=['page\tside\tmidpoint\tscale\tb0\tb1\tb2\tconstantFallback']
    for case in cases:
        for side in ('SE','ME'):
            fit=case['_rows'][side]['fit']
            assert fit['polynomialBasis']=='centered_scaled'
            lines.append('\t'.join(map(str,[case['page'],side,fit['basisMidpoint'],fit['basisScale'],*fit['basisParameters'],int(fit.get('constantBackgroundFallback',False))])))
    (folder/'basis_fit_parameters.tsv').write_text('\n'.join(lines)+'\n')


def omit_uncalibrated_plot_probabilities(folder):
    path=folder/'cases_fit_meta.tsv'
    with path.open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');columns=reader.fieldnames;rows=list(reader)
    for row in rows:row['probability']='nan'
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def diagnostic(row,edges,limits):
    fit=row['fit'];p=fit['parameters']
    indices=[i for i,(a,b) in enumerate(zip(edges,edges[1:])) if limits[0]<=(a+b)/2<limits[1]]
    n=sum(row['massCounts'][i] for i in indices)
    zeros=sum(row['massCounts'][i]==0 for i in indices)
    out=dict(observations=n,zeroBins=zeros,totalBins=len(indices),deviance=None,
             minimumExpected=None,domainValid=False)
    if any(x is None or not math.isfinite(x) for x in p) or p[2]<=0:return out
    means=[];deviance=0
    for i in indices:
        x=(edges[i]+edges[i+1])/2
        if fit.get('polynomialBasis')=='centered_scaled':
            u=(x-fit['basisMidpoint'])/fit['basisScale'];b=fit['basisParameters']
            background=b[0]+b[1]*u+b[2]*u*u
        else:
            background=p[3]+p[4]*x+p[5]*x*x
        mu=p[0]*math.exp(-.5*((x-p[1])/p[2])**2)+background
        means.append(mu);count=row['massCounts'][i]
        if not math.isfinite(mu) or mu<0 or (count>0 and mu<=0):
            out['minimumExpected']=min(means);return out
        deviance+=2*max(0.,mu-count+count*math.log(count/mu) if count else mu)
    out.update(domainValid=True,deviance=deviance,minimumExpected=min(means))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--previous-cases',type=Path,required=True)
    parser.add_argument('--case-settings',type=Path,required=True)
    parser.add_argument('--inspection-max',type=float,required=True)
    parser.add_argument('--raw-likelihood-dir',type=Path,required=True)
    parser.add_argument('--zero-seed-centered-dir',type=Path,required=True)
    parser.add_argument('--before-fallback-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    context=json.loads((out/'poisson_context.json').read_text())
    cfg=yaml.safe_load(args.case_settings.read_text())
    data={key:json.loads((out/(key+'_extract/purity_fit_review.json')).read_text()) for key in ('group','fine')}
    old={key:json.loads(Path(item['path']).read_text()) for key,item in context['references'].items()}
    compared=[];summaries=[];deviance_checks=0
    for key,current in data.items():
        assert current['fitOptions']=='SLQNR0' and current['fitObjective']=='poisson'
        assert current['settings']==old[key]['settings']
        for group,previous in zip(current['groups'],old[key]['groups']):
            assert group['id']==previous['id']
            for method in ('method3','method5'):
                for side in ('SE','ME'):
                    before=previous['methods'][method][side];after=group['methods'][method][side]
                    assert len(before)==len(after)
                    for left,right in zip(before,after):
                        for field in ('massCounts','massError2','centralMassCounts','centralMassError2','windowNtot','centralNtot'):
                            assert left[field]==right[field],(key,group['id'],method,side,left['bin'],field)
                        first=diagnostic(left,current['massEdges'],current['settings']['fitRange'])
                        second=diagnostic(right,current['massEdges'],current['settings']['fitRange'])
                        fit=right['fit']
                        if fit['status']!=-1:
                            assert fit['usedFitBins']==second['totalBins']
                            assert fit['emptyFitBins']==second['zeroBins']
                            assert fit['expectationValid']==second['domainValid']
                            if second['domainValid']:
                                assert math.isclose(fit['deviance'],second['deviance'],rel_tol=2.e-7,abs_tol=2.e-6)
                                deviance_checks+=1
                        compared.append(dict(source=key,group=group['id'],method=method,side=side,bin=left['bin'],x=left['x'],
                            chiSquare=first,poisson=second,oldValid=left['full']['valid'],newValid=right['full']['valid'],
                            oldReason=left['fit']['reason'],newReason=right['fit']['reason']))
                    inspected=[(a,b) for a,b in zip(before,after) if b['x']<args.inspection_max]
                    summaries.append(dict(source=key,group=group['id'],method=method,side=side,total=len(inspected),
                        oldValid=sum(bool(a['full']['valid']) for a,b in inspected),
                        newValid=sum(bool(b['full']['valid']) for a,b in inspected),
                        newReasons=dict(collections.Counter(b['fit']['reason'] for a,b in inspected))))
    fixed=json.loads(args.previous_cases.read_text())
    all_cases=candidates(data['fine'],'fine',cfg)+candidates(data['group'],'grouped',cfg)
    alternatives={}
    for name,folder in [('rawLikelihoodVersion',args.raw_likelihood_dir),('zeroSeedCenteredVersion',args.zero_seed_centered_dir),('beforeConstantFallbackVersion',args.before_fallback_dir)]:
        fine=json.loads((folder/'fine_extract/purity_fit_review.json').read_text())
        grouped=json.loads((folder/'group_extract/purity_fit_review.json').read_text())
        alternatives[name]=candidates(fine,'fine',cfg)+candidates(grouped,'grouped',cfg)
    selected=[]
    for previous in fixed['cases']:
        matches=[c for c in all_cases if c['sourceInput']==previous['sourceInput'] and c['id']==previous['id']]
        assert len(matches)==1
        case=matches[0];case['page']=previous['page'];case['category']=previous['category']
        case['chiSquareVersion']=previous
        for name,other in alternatives.items():
            matched=next(c for c in other if c['sourceInput']==previous['sourceInput'] and c['id']==previous['id'])
            case[name]={k:v for k,v in matched.items() if not k.startswith('_')}
        case['poissonDiagnostic']={s:diagnostic(case['_rows'][s],case['_massEdges'],cfg['mass_display_range']) for s in ('SE','ME')}
        selected.append(case)
    write_fit_tables(selected,out)
    omit_uncalibrated_plot_probabilities(out)
    write_basis(selected,out)
    save(out/'fixed_cases_comparison.json',dict(cases=[{k:v for k,v in c.items() if not k.startswith('_')} for c in selected],
        policy='Same nine case identities as the preceding chi-square report; no reselection by improved appearance.'))
    atlas=[]
    for case in candidates(data['group'],'grouped',cfg):
        if case['low']>=args.inspection_max:continue
        case['page']=len(atlas)+1;atlas.append(case)
    folder=out/'atlas_tables';folder.mkdir(exist_ok=True);write_fit_tables(atlas,folder)
    omit_uncalibrated_plot_probabilities(folder)
    write_basis(atlas,folder)
    fine_rows=['objective\tmethod\tgroup\tlabel\tside\tbin\tlow\thigh\tvalue\terror']
    for objective,source in [('poisson',data['fine']),('chiSquare',old['fine'])]:
        for case in candidates(source,'fine',cfg):
            for side in ('SE','ME'):
                value=case['purity'+side]
                fields=[objective,case['method'],case['group'],case['centLabel'],side,case['bin'],case['low'],case['high'],
                        value['value'] if value else 'nan',value['error'] if value else 'nan']
                fine_rows.append('\t'.join(map(str,fields)))
    (out/'fine_purity_comparison.tsv').write_text('\n'.join(fine_rows)+'\n')
    save(out/'fit_objective_comparison.json',dict(rows=compared,summary=summaries,
        matchedCountRows=len(compared),independentDevianceChecks=deviance_checks,
        inspectionRange=[0,args.inspection_max],atlasPages=len(atlas)))
    norm=yaml.safe_load((Path(context['references']['fine']['path']).parents[1]/'normalization_settings.yaml').read_text())
    norm.update(input=str(out/'group_extract/purity_fit_review.json'),output_dir=str(out),
                candidates=[{'id':'norm_0p4_0p5','range':[.4,.5]},{'id':'norm_0p3_0p4','range':[.3,.4]}],selected='norm_0p4_0p5')
    (out/'normalization_settings.yaml').write_text(yaml.safe_dump(norm,sort_keys=False))
    save(out/'analysis_validation.json',dict(matchedCountRows=len(compared),independentDevianceChecks=deviance_checks,
        fixedCaseIdentities=[c['id'] for c in selected],atlasPages=len(atlas),inputBinningGateAndGaussianLimitsUnchanged=True,
        negativeBackgroundConstant=context.get('negativeBackgroundConstant'),
        polynomialBasis=context.get('polynomialBasis'),backgroundSeedFloor=context.get('backgroundSeedFloor'),
        maximumFunctionCalls=context.get('maximumFunctionCalls')))
    print(json.dumps(dict(compared=len(compared),devianceChecks=deviance_checks,cases=len(selected),atlasPages=len(atlas))))


if __name__=='__main__':main()
