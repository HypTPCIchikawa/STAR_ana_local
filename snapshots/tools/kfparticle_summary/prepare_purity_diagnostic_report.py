#!/usr/bin/env python3
"""Prepare report prose, diagnostic tables and labeled source-figure copies."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def num(value, digits=4):
    return format(value, '.%dg'%digits) if isinstance(value,(int,float)) else 'NA'


def measured(value):
    return num(value['value'],5)+' +/- '+num(value['error'],3) if value else 'NA'


def interval(values):
    return num(values[0],3)+'-'+num(values[1],3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve(strict=True)
    repo=Path(__file__).resolve().parents[2]
    audit=json.loads((out/'normalization_audit.json').read_text())
    examples=json.loads((out/'selected_purity_cases.json').read_text())
    cited=json.loads((out/'cited_sources.json').read_text())
    selected=next(c for c in audit['candidates'] if c['id']==audit['settings']['selected'])
    reference={g['id']:g for g in selected['groups']}
    sensitivity=[]
    for candidate in audit['candidates']:
        for group in candidate['groups']:
            for method,result in group['methods'].items():
                first=result['normalization']['factor']
                second=reference[group['id']]['methods'][method]['normalization']['factor']
                shift=first['value']/second['value']-1 if first and second else None
                result['relativeScaleChangeFromSelected']=shift
                sensitivity.append([candidate['id'],group['id'],method,shift])
    audit['recommendation'].update(
        status='applied as an explicit common normalization reference',
        normalization=selected['normalization'],
        reason='The nearest disjoint 100-MeV interval excludes the failed purity cell and mixed CF bins; all 15 grouped normalizations are valid. This is a reference convention, not a proof of a physics plateau in all centralities.')
    save(out/'normalization_audit.json',audit)
    with (out/'normalization_sensitivity.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['candidate','group','method','relativeScaleChangeFromSelected'])
        writer.writerows(sensitivity)
    raw=out/'quoted_figures_raw';fig=out/'quoted_figures'
    raw.mkdir(exist_ok=True);fig.mkdir(exist_ok=True)
    for era in ('historical','updated'):
        for method in ('method3','method5'):
            source=Path(cited['files'][era+'/kstar_cent_CFsub_'+method+'_4He.pdf']['path'])
            target=raw/(era+'_'+method+'.pdf')
            if target.exists():
                raise FileExistsError(target)
            shutil.copy2(str(source),str(target))
    subprocess.run(['python3',str(repo/'tools/kfparticle_summary/crop_method6_figures.py'),
                    '--input-dir',str(raw),'--output-dir',str(fig),'--margin-points','6'],check=True)
    failed=audit['failedFits'][0]
    sections=[
      dict(title='Scope, updated result and sample',paragraphs=[
        'The common normalization reference for Method3, Method5 and Method6 is changed from 0.30-0.40 to 0.40-0.50 GeV/c. The same mass model and fit-validity gates are retained. Comparisons use 20-80%/0-20% and 0-10%/10-20%/20-60%.',
        'The input is the 2021 P24iy ROOT from SUMS 04BF2DAE47164B6F040E3A460E539289: 26,955/26,994 subjobs, 39 excluded, 2,100,726,543 input events and 2,067,647,474 reconstructed events. PID, Maker selection and mixing are unchanged. No production resubmission is made.',
        'The top-level PDF cited in the request is the historical October 3 report (old 0BEF sample). The October 30 report being updated uses the 04BF sample. Both versions of Figs.13/14 are quoted with explicit labels. New mass diagnostics correspond to the updated 04BF points, not a re-fit of the historical sample.']),
      dict(title='Why Method5 20-80% could not be normalized',newpage=True,paragraphs=[
        'The failed ME purity cell is 0.30-0.35 GeV/c, fitted in the sliding k* window 0.20-0.45. Purity is applied to native 5-MeV counts before combining to 20-MeV CF bins. Thus bins 0.30-0.32, 0.32-0.34 and 0.34-0.36 are unavailable. The old normalization interval 0.30-0.40 contains them. Requiring all normalization contributions masks the entire CF.',
        'This is not a sparse ME spectrum: '+format(failed['observationsInFitRange'],',')+' counts lie in the fitted mass range 1.08-1.15 GeV/c^2, versus '+format(failed['storedMassObservations'],',')+' over the full stored mass axis. Status=4 and isValid=false. Covariance status=2 alone cannot rescue an invalid fit; background integral/covariance is unavailable.',
        'The unsuccessful saved Gaussian plus quadratic-background curve is drawn only for diagnosis. Large structured peak/tail pulls are visible. The precise contributions from peak-response shape, background shape and minimization cannot be separated by this one fit. High statistics, numerical convergence and model adequacy are distinct issues.'],
        tables=[dict(columns=['ME diagnostic','Value'],rows=[
          ['Purity cell [GeV/c]','0.30-0.35'],['Sliding k* window','0.20-0.45'],
          ['Counts in fit mass range',format(failed['observationsInFitRange'],',')],
          ['Counts over stored mass axis',format(failed['storedMassObservations'],',')],
          ['Status / covariance status','4 / 2'],
          ['chi2 / NDF',num(failed['fit']['chi2'],8)+' / '+str(failed['fit']['ndf'])],
          ['Unaccepted Gaussian mean',num(failed['fit']['parameters'][1],10)],
          ['Unaccepted Gaussian sigma',num(failed['fit']['parameters'][2],10)]])]),
      dict(title='Normalization reference and range sensitivity',newpage=True,paragraphs=[
        'Use 0.40-0.50 GeV/c as the common reference: the nearest disjoint interval above the old range, with the same 100-MeV width. All 15 grouped normalizations become valid. This does not recover the original failed CF bins or Method3 low-k* missing points.',
        'The raw SE/ME ratio in the fixed Lambda mass gate provides a check independent of fitted purity. For 20-80%, its constant-fit chi2/NDF is 3.155/4 over 0.40-0.50, with slope -0.002608 +/- 0.007277 per GeV/c. For 0-20%, the corresponding values are 10.017/4 and -0.007428 +/- 0.002817. Count-only diagnostics do not certify a flat physics plateau for every centrality.',
        'At larger k* the raw ratio falls. Widening the interval therefore changes CF scale, not only statistical precision. The table shows common multiplicative shifts relative to 0.40-0.50. They are range sensitivities, not independent statistical errors or a completed systematic uncertainty, and are not included in plotted error bars.',
        'Integrated CF means the normalized ratio of corrected yields integrated over k*=0-0.20 GeV/c, not an arithmetic average of CF points. Method3 is unavailable there if low-k* fits are missing.'],
        equations=[r"C_i=\frac{S_i}{M_i}\frac{\sum_{j\in N}M_j}{\sum_{j\in N}S_j},\qquad C_i^{(N')}/C_i^{(N)}=A_{N'}/A_N."],
        tables=[dict(columns=['Norm [GeV/c]','20-80 M5','20-80 M6','0-20 M5'],rows=[]),
                dict(columns=['Centrality','M5 integrated CF','M6 integrated CF'],rows=[])]),
      dict(title='Figs.13 and 14: contrasting mass spectra',newpage=True,paragraphs=[
        'Figs.13/14 contain nine fine centralities, 70-80% through 0-5%. All nine were extracted from the same updated ROOT with unchanged fit settings: 1,260 fits (900 Method3, 360 Method5) over k*=0-1 GeV/c.',
        'Explicit recorded criteria select multiple discrepant points and agreement controls. Four large-difference examples have only 27-52 SE counts in the actual fit mass range, with purity errors 0.159-0.423. These demonstrate unstable low-statistics estimation, not a precise physical SE/ME purity difference.',
        'Method3 fits its central 20-MeV cell. Method5 fits a neighboring-cell window, up to 250 MeV, distinct from its central 50-MeV purity cell. Each page shows raw fit-window spectra, saved fits and pulls, ROOT DrawNormalized overlays of fit-window spectra, and normalized central-cell spectra. Unit area is defined in 1.08-1.15 GeV/c^2; red is SE, blue ME.',
        'Only one Method5 agreement control satisfies the declared distinct-centrality criteria. The missing second control is recorded rather than obtained by relaxing thresholds.'],
        tables=[dict(columns=['Case','M','Cent.','k* cell','P(SE)','P(ME)','Nfit SE/ME'],rows=[])]),
      dict(title='Interpretation: statistics and model adequacy',newpage=True,paragraphs=[
        'Cases2/3 (Method3) and cases6/7 (Method5) have sparse SE spectra, large purity errors and several parameter-bound warnings. Case5 is an adjacent sparse SE fit rejected altogether. Separated purity points are not precise conflicting measurements. Small chi2/NDF in sparse spectra is not proof of a good background estimate: zero-error bins are omitted from the saved chi-square fit and pull.',
        'Method3 agreement case4 (0-5%, 0.06-0.08) has 524 SE fit-range counts, versus 27-28 in the discrepant examples. SE and ME purity agree; SE chi2/NDF is 0.705. This is a useful contrast, though ME chi2/NDF=3.77 still indicates model concerns.',
        'Method5 agreement case8 (0-5%, 0.25-0.30) has 89,303 SE and 3,307,029 ME fit-range counts, close normalized mass shapes, and purity difference 0.000705. Yet fit chi2/NDF is 18.44/599.79. Similar model bias can coexist with matching purity. Agreement alone does not validate the estimate or its uncertainty.',
        'Method5 case9 (40-50%, 0.20-0.25) has 1,083/48,935 SE/ME fit-range counts, purity difference 0.06315, different normalized shapes and ME chi2/NDF=7.83. This is not simply a few-dozen-count SE fit. Mass-model inadequacy and different signal/background weights within the sliding window remain possible; the origin is not uniquely identified.',
        'SE/ME purity need not be identical if signal and background correlations differ. For ME purity P, relative SE signal/background weights Cs/Cb give the relation below. Method6 is the requested ME-to-SE transfer hypothesis; its copied statistical error does not include transfer or model systematics.'],
        equations=[r'P_{\rm SE}=C_sP_{\rm ME}/[C_sP_{\rm ME}+C_b(1-P_{\rm ME})].']),
      dict(title='Errors and reproducibility',paragraphs=[
        'Method6 shares one Method5 ME purity estimator and its error across both SE and ME. Weighting is done on native 5-MeV bins before 20-MeV rebinning. Common-purity and normalization covariance are propagated, without independently counting the same purity error twice. Twelve existing shared-purity and normalization tests pass.',
        'Uncertainties omit event/mixed-pair covariance absent from these histograms, fitted-background/count covariance, covariance between distinct background fits, and mass-model/ME-to-SE transfer systematics. The normalized-shape chi-square is a fixed-area diagonal descriptive diagnostic with no formal p-value.',
        'The original 0.30-0.40 CF/LL-fit figures remain reference results. The new 0.40-0.50 comparisons take precedence for the stated comparison. LL parameters are not refitted. Input ROOT, original PDFs, commands, settings, full fits, case selection and hashes are retained.'])]
    for candidate in audit['candidates']:
        groups={g['id']:g for g in candidate['groups']}
        row=[interval(candidate['normalization'])]
        for gid,method in [('cent20_80','method5'),('cent20_80','method6'),('cent0_20','method5')]:
            shift=groups[gid]['methods'][method]['relativeScaleChangeFromSelected']
            row.append('%+.2f%%'%(100*shift) if shift is not None else 'NA')
        sections[2]['tables'][0]['rows'].append(row)
    for group in selected['groups']:
        sections[2]['tables'][1]['rows'].append([group['label'],
            measured(group['methods']['method5']['integrated']['cf']),
            measured(group['methods']['method6']['integrated']['cf'])])
    for case in examples['cases']:
        sections[3]['tables'][0]['rows'].append([str(case['page']),case['method'][-1],case['centLabel'],
            interval([case['low'],case['high']]),measured(case['puritySE']),measured(case['purityME']),
            str(case['SE']['fitRangeObservations'])+'/'+str(case['ME']['fitRangeObservations'])])
    figures=[dict(path=str(out/'figures/normalized_case_001.pdf'),
        caption='Case 1: the original normalization failure',paragraphs=[
        '20-80%, Method5, purity cell 0.30-0.35 and sliding window 0.20-0.45 GeV/c. The ME curve fails with '+
        format(failed['observationsInFitRange'],',')+' fit-range counts and large structured pulls. No ME purity is assigned.']),
        dict(path=str(out/'figures/raw_ratio_normalization.pdf'),caption='Raw SE/ME ratio and reference interval',paragraphs=[
        'Green: selected 0.40-0.50; dotted: preceding 0.30-0.40. The trend at larger k* shows why normalization-range dependence must be retained. Raw mass-gate count ratios do not use fitted purity; event/mixed-pair covariance is absent.']),
        dict(path=str(out/'figures/cf_methods_3_5_6.pdf')),
        dict(path=str(out/'figures/purity_se_me.pdf'))]
    for era in ('historical','updated'):
        for method,number in [('method3',13),('method5',14)]:
            label='Historical October3, old 0BEF sample' if era=='historical' else 'Updated October30, current 04BF sample'
            figures.append(dict(path=str(fig/(era+'_'+method+'.pdf')),page=2,
                caption='Quoted Fig.'+str(number)+': '+label,
                paragraphs=['Unchanged purity panels from '+method+' source PDF, page2. The following new mass cases diagnose the updated 04BF points.']))
    for case in examples['cases'][1:]:
        figures.append(dict(path=str(out/'figures'/('normalized_case_%03d.pdf'%case['page'])),
            caption='Case '+str(case['page'])+': Fig.'+str(case['sourceFigure']['figure'])+'; '+case['centLabel']+'; '+interval([case['low'],case['high']])+' GeV/c',
            paragraphs=['Purity cell '+interval([case['low'],case['high']])+', actual k* fit window '+interval([case['windowLow'],case['windowHigh']])+' GeV/c. P(SE)='+measured(case['puritySE'])+'; P(ME)='+measured(case['purityME'])+'. Fit-range counts SE/ME='+str(case['SE']['fitRangeObservations'])+'/'+str(case['ME']['fitRangeObservations'])+'.',case['interpretation']]))
    save(out/'report_content.json',dict(stem='summary_4He_KFParticle_3p85_updated_20261102',
        title='mass-fit diagnostics and normalization update',cover=[
          'This edition adds a diagnostic study of Figs.13/14 and changes the common Method3/5/6 normalization reference to 0.40-0.50 GeV/c in the addendum. Original pages and figure numbers remain references. The addendum supersedes the corresponding 0.30-0.40 CF comparison; original LL parameters are not refitted.',
          'The 20-80% Method5/6 CF can now be computed. The failed purity cell, missing CF bins, mass-model problems and normalization-range sensitivity remain explicitly documented.',
          'Several discrepant fine-centrality points are unstable low-statistics SE estimates. Agreement controls and a higher-statistics mismatch show why purity agreement alone cannot validate the fits.',
          'Historical top-level Figs.13/14 (October3, 0BEF) and updated Figs.13/14 (October30, 04BF) are quoted separately. All new mass diagnostics use the updated 2021 P24iy 04BF ROOT.'],sections=sections,figures=figures))
    files={}
    for name in ('DrawNormalizedPurityCases.C','DrawNormalizationRatioAudit.C','audit_normalization_ranges.py',
                 'select_purity_diagnostic_cases.py','purity_diagnostic_cases.yaml','build_purity_diagnostic_update.py',
                 'prepare_purity_diagnostic_report.py','validate_purity_diagnostic_update.py',
                 'publish_purity_diagnostic_update.py','method6_statistics.py'):
        path=repo/'tools/kfparticle_summary'/name
        files[name]=dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    save(out/'diagnostic_source_code.json',files)
    print(json.dumps(dict(sections=len(sections),figures=len(figures),cases=len(examples['cases']))))


if __name__=='__main__':
    main()
