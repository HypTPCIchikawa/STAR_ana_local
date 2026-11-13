"""Select reproducible mass-fit examples from saved fine-centrality fits.

The saved fits are never rerun. All selection thresholds and requested example
counts arrive in an explicit YAML profile. Shape comparisons are descriptive:
event and mixed-pair correlations are not recoverable from these projections.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import yaml


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bytes=path.stat().st_size)


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def mass_indices(edges, interval):
    return [i for i, (lo, hi) in enumerate(zip(edges, edges[1:]))
            if interval[0] <= (lo + hi) / 2 < interval[1]]


def count(values, indices):
    return sum(values[i] for i in indices)


def estimate(row):
    full = row['full']
    total = full['windowNtot']
    background = full['background']
    if not full['valid'] or not finite(total) or total <= 0 or not finite(background):
        return None
    purity = 1 - background / total
    if not 0 < purity <= 1:
        return None
    variance = (full['backgroundVariance'] / total ** 2 +
                background ** 2 * full['windowError2'] / total ** 4)
    require(finite(variance) and variance >= 0, 'invalid purity variance')
    return dict(value=purity, variance=variance, error=math.sqrt(variance))


def shape_comparison(first, second, indices):
    """Compare normalized shapes with diagonal-only descriptive chi square."""
    a, va, b, vb = first[0], first[1], second[0], second[1]
    na, nb = count(a, indices), count(b, indices)
    if na <= 0 or nb <= 0:
        return dict(available=False, reason='nonpositive mass-range integral')
    terms = []
    for i in indices:
        pa, pb = a[i] / na, b[i] / nb
        # This is the same fixed-scale diagonal variance as TH1::DrawNormalized.
        variance = va[i] / na ** 2 + vb[i] / nb ** 2
        if variance > 0:
            terms.append((pa - pb) ** 2 / variance)
    return dict(available=True, normalizationSE=na, normalizationME=nb,
                diagonalChi2=sum(terms), comparedBins=len(terms),
                descriptiveNdf=max(0, len(terms) - 1),
                diagonalChi2PerNdf=(sum(terms)/(len(terms)-1) if len(terms)>1 else None),
                totalVariation=0.5*sum(abs(a[i]/na-b[i]/nb) for i in indices),
                maximumAbsoluteBinDifference=max(abs(a[i]/na-b[i]/nb) for i in indices),
                caveat='Fixed-scale diagonal diagnostic only; normalization covariance and event/mixed-pair correlations are omitted. No formal p value is assigned.')


def side_summary(row, indices, warning):
    fit = dict(row['fit'])
    ratio = fit['chi2']/fit['ndf'] if finite(fit['chi2']) and fit['ndf']>0 else None
    fit['chi2PerNdf'] = ratio
    fit['modelWarning'] = ratio is not None and ratio > warning
    return dict(storedMassObservations=sum(row['massCounts']),
                fitRangeObservations=count(row['massCounts'], indices),
                centralStoredMassObservations=sum(row['centralMassCounts']),
                centralFitRangeObservations=count(row['centralMassCounts'], indices),
                signalWindowCounts=row['full']['windowNtot'],
                centralSignalWindowCounts=row['full']['centralNtot'],
                purity=estimate(row), fit=fit,
                validityReason=row['fit']['reason'],
                fullValidity=row['full']['valid'])


def candidates(data, source_key, cfg):
    indices = mass_indices(data['massEdges'], cfg['mass_display_range'])
    require(indices, 'mass display range contains no bins')
    labels = {int(k): v for k, v in cfg['fine_centrality_labels'].items()}
    output = []
    for group in data['groups']:
        fine = group['centMin'] == group['centMax']
        label = labels[group['centMin']] if fine else cfg['group_labels'][group['id']]
        for method, grid in [('method3', 'kstarEdges'), ('method5', 'method5PurityEdges')]:
            edges = group[grid]
            se_rows, me_rows = group['methods'][method]['SE'], group['methods'][method]['ME']
            require(len(se_rows)==len(me_rows), 'SE/ME row count differs')
            for index, (se, me) in enumerate(zip(se_rows, me_rows)):
                low, high = edges[index:index+2]
                summary_se = side_summary(se, indices, cfg['model_warning_chi2_ndf'])
                summary_me = side_summary(me, indices, cfg['model_warning_chi2_ndf'])
                ps, pm = summary_se['purity'], summary_me['purity']
                delta = ps['value']-pm['value'] if ps is not None and pm is not None else None
                native = group['originalKstarEdges']
                first, last = se['windowOriginalFirstBin']-1, se['windowOriginalLastBin']
                require(first==me['windowOriginalFirstBin']-1 and last==me['windowOriginalLastBin'],
                        'SE/ME fit window differs')
                output.append(dict(id='%s_%s_%03d' % (method, group['id'], se['bin']),
                    method=method, group=group['id'], centLabel=label,
                    centMin=group['centMin'], centMax=group['centMax'], bin=se['bin'],
                    low=low, high=high, windowLow=native[first], windowHigh=native[last],
                    sourceInput=source_key, sourceFigure=(cfg['source_figures'][method] if fine else None),
                    puritySE=ps, purityME=pm, purityDifference=delta,
                    absolutePurityDifference=abs(delta) if delta is not None else None,
                    SE=summary_se, ME=summary_me,
                    shapeComparison=dict(
                        fitWindow=shape_comparison((se['massCounts'],se['massError2']),
                                                   (me['massCounts'],me['massError2']), indices),
                        centralCell=shape_comparison((se['centralMassCounts'],se['centralMassError2']),
                                                     (me['centralMassCounts'],me['centralMassError2']), indices)),
                    _rows=dict(SE=se, ME=me), _massEdges=data['massEdges']))
    return output


def interpretation(case):
    messages = []
    if case['puritySE'] is None:
        messages.append('The SE purity is unavailable; inspect sparse counts and the recorded fit rejection rather than interpreting it as zero signal.')
    if case['purityME'] is None:
        messages.append('The ME purity is unavailable. If this cell contributes to normalization, complete-normalization rules mask the entire CF.')
    if case['purityDifference'] is not None:
        messages.append('The displayed SE-minus-ME purity difference is %.6g.' % case['purityDifference'])
    for side in ('SE','ME'):
        s = case[side]
        if s['fit']['modelWarning']:
            messages.append('%s has chi2/ndf=%.4g: numerical convergence or purity agreement alone does not establish a suitable mass model.' % (side,s['fit']['chi2PerNdf']))
        if s['fit']['boundary']:
            messages.append('%s has a fitted parameter at a configured bound.' % side)
    if case['method']=='method5':
        messages.append('Method5 fits the neighboring-cell window, not only the central 50-MeV cell; assess its fit-window statistics separately from central-cell statistics.')
    messages.append('Normalized-shape agreement and fitted-purity agreement are different diagnostics; this example alone does not establish a unique cause or a total purity uncertainty.')
    return ' '.join(messages)


def match(case, rule):
    n = case['SE']['fitRangeObservations']
    delta = case['absolutePurityDifference']
    if case['method'] != rule['method']:
        return False
    if rule['require_valid_se'] != (case['puritySE'] is not None):
        return False
    if rule['require_valid_me'] != (case['purityME'] is not None):
        return False
    if rule['minimum_se_fit_count'] is not None and n < rule['minimum_se_fit_count']:
        return False
    if rule['maximum_se_fit_count'] is not None and n > rule['maximum_se_fit_count']:
        return False
    if rule['minimum_abs_delta'] is not None and (delta is None or delta < rule['minimum_abs_delta']):
        return False
    if rule['maximum_abs_delta'] is not None and (delta is None or delta > rule['maximum_abs_delta']):
        return False
    return True


def ranking(case, order):
    delta = case['absolutePurityDifference']
    if order == 'largest_difference':
        return (-delta, -case['SE']['fitRangeObservations'],case['centMin'],case['bin'])
    if order == 'smallest_difference':
        return (delta, -case['SE']['fitRangeObservations'],case['centMin'],case['bin'])
    if order == 'largest_se_count':
        return (-case['SE']['fitRangeObservations'],case['centMin'],case['bin'])
    raise ValueError('Unknown ranking order: '+order)


def clean(text):
    return str(text).replace('\t',' ').replace('\n',' ')


def write_fit_tables(selected, target):
    meta_columns = ['page','group','label','method','bin','low','high','windowLow','windowHigh',
                    'side','status','covariance','ndf','chi2','probability','purity','purityError',
                    'observations','p0','p1','p2','p3','p4','p5','reason']
    bin_columns = ['page','side','index','low','high','count','error2','centralCount','centralError2']
    with (target/'cases_fit_meta.tsv').open('w') as mf, (target/'cases_fit_bins.tsv').open('w') as bf:
        mw=csv.writer(mf,delimiter='\t',lineterminator='\n')
        bw=csv.writer(bf,delimiter='\t',lineterminator='\n')
        mw.writerow(meta_columns);bw.writerow(bin_columns)
        for case in selected:
            for side in ('SE','ME'):
                row=case['_rows'][side];fit=row['fit'];purity=case[side]['purity']
                values=[case['page'],case['group'],case['centLabel'],case['method'],case['bin'],
                        case['low'],case['high'],case['windowLow'],case['windowHigh'],side,
                        fit['status'],fit['covarianceStatus'],fit['ndf'],fit['chi2'],fit['probability'],
                        purity['value'] if purity else None,purity['error'] if purity else None,
                        fit['observationsInStoredMassRange']]+fit['parameters']+[fit['reason']]
                mw.writerow(['nan' if v is None else clean(v) for v in values])
                edges=case['_massEdges']
                for i,(lo,hi) in enumerate(zip(edges,edges[1:])):
                    bw.writerow([case['page'],side,i+1,lo,hi,row['massCounts'][i],row['massError2'][i],
                                 row['centralMassCounts'][i],row['centralMassError2'][i]])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fine-input',required=True,type=Path)
    parser.add_argument('--grouped-input',required=True,type=Path)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    cfg=yaml.safe_load(args.settings.read_text())
    fine=json.loads(args.fine_input.read_text());grouped=json.loads(args.grouped_input.read_text())
    require(Path(fine['inputRoot']).resolve()==Path(grouped['inputRoot']).resolve(),
            'fine and grouped inputs refer to different ROOT paths')
    require(fine['settings']['fitRange']==cfg['mass_display_range'],'configured mass range differs from fits')
    require(all(g['centMin']==g['centMax'] for g in fine['groups']),'fine input contains grouped centralities')
    require({g['centMin'] for g in fine['groups']}==set(range(9)),'fine input must contain all 9 centralities')
    all_fine=candidates(fine,'fine',cfg);all_grouped=candidates(grouped,'grouped',cfg)
    all_cases=all_fine+all_grouped
    selected=[];used=set();audit=[]
    for explicit in cfg['explicit_cases']:
        matches=[c for c in all_cases if c['sourceInput']==explicit['source'] and
                 c['method']==explicit['method'] and c['group']==explicit['group'] and c['bin']==explicit['bin']]
        require(len(matches)==1,'Explicit example is missing or ambiguous: '+repr(explicit))
        case=matches[0];case['category']=explicit['category'];selected.append(case)
        used.add((case['sourceInput'],case['id']))
    lo,hi=cfg['inspection_range']
    for rule in cfg['selection_rules']:
        pool=[c for c in all_fine if lo<=c['low'] and c['high']<=hi+1.e-10 and
              ('fine',c['id']) not in used and match(c,rule)]
        pool.sort(key=lambda c:ranking(c,rule['ranking']))
        taken=[];centralities=set()
        for case in pool:
            if rule['distinct_centralities'] and case['centMin'] in centralities:
                continue
            case['category']=rule['category'];selected.append(case);taken.append(case['id'])
            centralities.add(case['centMin']);used.add(('fine',case['id']))
            if len(taken)==rule['count']:
                break
        audit.append(dict(rule=rule,eligibleCandidates=len(pool),selected=taken,
                          requested=rule['count'],actual=len(taken),
                          shortage=len(taken)<rule['count']))
    require(selected,'No cases were selected')
    target=args.output_dir;target.mkdir(parents=True,exist_ok=True)
    products=['selected_purity_cases.json','all_purity_case_candidates.tsv','cases_fit_meta.tsv','cases_fit_bins.tsv']
    require(all(not (target/p).exists() for p in products),'Refusing to overwrite case products')
    for page,case in enumerate(selected,1):
        case['page']=page;case['interpretation']=interpretation(case)
    write_fit_tables(selected,target)
    with (target/'all_purity_case_candidates.tsv').open('w') as f:
        columns=['source','id','method','centrality','bin','low','high','windowLow','windowHigh','P_SE','P_ME',
                 'deltaP','SEfitCount','MEfitCount','SEcentralFitCount','MEcentralFitCount','SEchi2Ndf','MEchi2Ndf',
                 'SEstatus','MEstatus','SEcov','MEcov']
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(columns)
        for c in all_cases:
            se,me=c['SE'],c['ME']
            w.writerow([c['sourceInput'],c['id'],c['method'],c['centLabel'],c['bin'],c['low'],c['high'],
                        c['windowLow'],c['windowHigh'],c['puritySE']['value'] if c['puritySE'] else 'NA',
                        c['purityME']['value'] if c['purityME'] else 'NA',c['purityDifference'],
                        se['fitRangeObservations'],me['fitRangeObservations'],
                        se['centralFitRangeObservations'],me['centralFitRangeObservations'],
                        se['fit']['chi2PerNdf'],me['fit']['chi2PerNdf'],se['fit']['status'],me['fit']['status'],
                        se['fit']['covarianceStatus'],me['fit']['covarianceStatus']])
    serial=[]
    for c in selected:
        serial.append({key:value for key,value in c.items() if not key.startswith('_')})
    payload=dict(schemaVersion=1,sources=dict(fine=fingerprint(args.fine_input),
                 grouped=fingerprint(args.grouped_input),settings=fingerprint(args.settings)),
                 settings=cfg,figureReferences=cfg['source_figures'],selectionAudit=audit,cases=serial,
                 limitations=[
                    'Purity statistical errors omit fitted-background/count covariance.',
                    'Normalized-shape chi square is a fixed-scale diagonal descriptive diagnostic, not a formal hypothesis-test statistic.',
                    'No event-level or mixed-pair covariance, mass-model systematic, or purity-transfer systematic is included.',
                    'Agreement of fitted purities does not establish adequacy of either fit model; disagreement alone does not prove low statistics.'
                 ])
    (target/'selected_purity_cases.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(selected),selectionShortages=[a['rule']['category'] for a in audit if a['shortage']],
                         output=str(target/'selected_purity_cases.json'))))


if __name__=='__main__':
    main()
