#!/usr/bin/env python3
"""Prepare a report of Poisson fits with a nonnegative constant fallback."""
import argparse
import json
from pathlib import Path


def number(x,n=4):return format(x,'.%dg'%n) if isinstance(x,(int,float)) else 'NA'
def measured(p):return number(p['value'],5)+' +/- '+number(p['error'],3) if p else 'NA'
def interval(a,b):return number(a,3)+'-'+number(b,3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    context=json.loads((out/'poisson_context.json').read_text())
    examples=json.loads((out/'fixed_cases_comparison.json').read_text())['cases']
    comparison=json.loads((out/'fit_objective_comparison.json').read_text())
    audit=json.loads((out/'normalization_audit.json').read_text())
    results=json.loads((out/'methods_comparison_selected.json').read_text())
    toy=json.loads((out/'empty_bin_toy.json').read_text())
    sources={k:json.loads((out/(k+'_extract/purity_fit_review.json')).read_text()) for k in ('group','fine')}
    fallback_counts=[]
    zero_boundary=[]
    for name,source in sources.items():
        for method in ('method3','method5'):
            for side in ('SE','ME'):
                rows=[r for g in source['groups'] for r in g['methods'][method][side]]
                fallback_counts.append([name,method[-1],side,str(sum(r['fit']['constantBackgroundFallback'] for r in rows)),str(len(rows)),str(sum(r['full']['valid'] for r in rows))])
                for g in source['groups']:
                    for r in g['methods'][method][side]:
                        if r['fit']['constantBackgroundFallback'] and r['fit']['basisParameters'][0]<1.e-8:
                            zero_boundary.append([name,g['id'],method[-1],side,number(r['x'],4)])
    former=examples[0]
    sections=[dict(title='Updated analysis and settings',paragraphs=[
        'The Lambda invariant-mass fits are re-evaluated with ROOT Poisson-likelihood option L, using SLQNR0 on unscaled, unsubtracted count histograms. Empty mass bins participate in the objective. Gaussian initialization/limits for the primary fit, mass range 1.08-1.15 GeV/c^2, count/covariance thresholds, signal mass gate and 20/50-MeV k* grids are retained. The Gaussian parameters remain free in the fallback, starting from the primary solution.',
        'The sample is the same 2021 P24iy 04BF2DAE47164B6F040E3A460E539289 merged ROOT: 26,955/26,994 subjobs, 39 excluded, 2,100,726,543 input events and 2,067,647,474 reconstructed events. No PID, Maker, mixing or production-job change is made.',
        'All 700 grouped and 1,260 fine-centrality fits are repeated. The same nine mass-spectrum examples from the preceding chi-square report are fixed in advance; examples are not reselected for agreement.',
        'For comparison with the preceding CF result, the primary normalization remains 0.40-0.50 GeV/c for all five selections and Method3/5/6. The old 0.30-0.40 interval is also recomputed to test whether its failed normalization is recovered.']),
      dict(title='Empty-bin effect and nonnegative-background rule',newpage=True,paragraphs=[
        'The default observed-error chi-square fit omits bins with zero counts and zero error. Those bins constrain the background in a Poisson fit: their deviance contribution is 2*mu. A constant-background toy with nine zero bins and one bin containing one count returns the values below in the actual ROOT5 runtime.',
        'The requested background rule is implemented as a refit, not as an unrecorded alteration of a fitted curve. First fit Gaussian plus quadratic background. Evaluate the background minimum over the closed mass-fit range, including any quadratic vertex. If that minimum is negative, refit with both the linear and quadratic background terms fixed to zero and the remaining constant b constrained to b>=0.',
        'The fallback has four free parameters (three Gaussian parameters and b), versus six for the quadratic fit. Its residual degrees of freedom are adjusted accordingly. b is fitted, so it may be positive or reach zero. The identity and diagnostics of every fallback are recorded separately for SE and ME.',
        'For ROOT finite-bound handling, the upper guard on b is max(1,Nfit), where Nfit is the total count within the fitted mass range. This is recorded and checked to be inactive; it is not a lower background floor. The numerical starting value uses the observed sideband mean with the explicit 0.1-count/bin seed floor, which does not constrain the final fitted b.'],
        equations=[r'D=2\sum_i[\mu_i-n_i+n_i\log(n_i/\mu_i)],\quad D_i(n_i=0)=2\mu_i.',
                   r'B_{\rm primary}(u)=b_0+b_1u+b_2u^2;\quad \min B<0\Rightarrow B_{\rm refit}(u)=b\geq0.'],
        tables=[dict(columns=['Toy diagnostic','Value'],rows=[['Chi-square background [count/bin]',number(toy['chiSquareEstimate'],8)],['Poisson background [count/bin]',number(toy['poissonEstimate'],8)],['Poisson deviance',number(toy['poissonDeviance'],8)]]),
                dict(columns=['Fits','M','Side','Constant refits','Total','Usable purity'],rows=fallback_counts)]),
      dict(title='Numerical stability and error interpretation',newpage=True,paragraphs=[
        'The full quadratic family uses an invertible coordinate: the mass offset from the fit midpoint divided by the half width of the fixed fit range. The curve remains Gaussian plus pol2. Coefficients and the full 6x6 covariance are converted back for compatibility; background integrals and their covariance are evaluated directly in the stable centered basis.',
        'Raw-coordinate L, centered-coordinate L with the old zero background seed, and the positive-seed version are retained as numerical trials. The final version uses an explicit 50,000-function-call budget. Therefore changes are not attributed to the letter L alone: initialization, numerical coordinates and the requested background fallback are all documented.',
        'Deviance and signed deviance residuals are independently recomputed from saved predictions, including empty bins. These residuals replace the preceding observed-error pull that dropped empty bins. The reported D/NDF is a diagnostic, not a calibrated sparse-count p-value.',
        'The quoted fit errors use the local covariance. In particular, b=0 is a boundary estimate, not proof that the background is exactly absent. Asymmetric profile intervals and model-selection uncertainty from choosing the fallback are not evaluated. Event/mixed-pair correlations and fitted-background/count covariance not available in the stored histograms are also omitted.'],
        tables=[dict(columns=['b near zero','Group','M','Side','k* center'],rows=zero_boundary)] if zero_boundary else []),
      dict(title='The same nine examples: chi-square versus L',newpage=True,paragraphs=[
        'The table reports accepted full-window purities. NA means a failed fit, unusable covariance or nonphysical result; it does not mean zero purity. The last column identifies the final SE/ME background models (Q: quadratic; C: nonnegative constant fallback). Counts and mass projections are unchanged.',
        'Fine-centrality comparisons reproduce the Method3/Method5 selections underlying the preceding Figs.13/14. The updated red/blue purity plots and the old chi-square estimates are displayed under the same k* and centrality definitions.'],
        tables=[dict(columns=['Case','M','Cent. / k*','Old P(SE)','L P(SE)','L P(ME)','BG S/M'],rows=[]),
                dict(columns=['Case','SE count','SE zero bins','D(SE)/NDF','D(ME)/NDF'],rows=[])]),
      dict(title='Normalization recovery and CF comparison',newpage=True,paragraphs=[
        'Method3 fits each 20-MeV CF cell. Method5 applies its 50-MeV sliding-window purity to native 5-MeV counts before 20-MeV aggregation. Method6 applies the same Method5 ME purity and error to both SE and ME as a shared estimator, retaining common-purity and normalization covariance.',
        'The previously failed 20-80% ME purity cell at 0.30-0.35 has L purity '+measured(former['purityME'])+'. Its status is '+str(former['ME']['fit']['status'])+', covariance status '+str(former['ME']['fit']['covarianceStatus'])+'. Numerical recovery does not guarantee the Gaussian/background model is adequate: the independently evaluated deviance remains '+number(former['ME']['fit']['deviance'],7)+'.',
        'The common 0.40-0.50 comparison isolates the change from the preceding chi-square result at the same normalization. The old 0.30-0.40 calculation is shown separately. Missing fit contributions continue to mask their CF bins; invalid fits are never silently replaced with zero background or an interpolated purity.',
        'Integrated CF values are ratios of corrected yields summed over k*=0-0.20 GeV/c with the stated normalization, not arithmetic averages of plotted CF points. No LL/source-radius fit is repeated in this mass-purity/CF update.'],
        tables=[dict(columns=['Norm','Centrality','Valid methods'],rows=[]),
                dict(columns=['Centrality','M3 points','M5 points','M6 points','Integrated M5','Integrated M6'],rows=[])]),
      dict(title='What can be concluded',paragraphs=[
        'The empty-bin mechanism is real and demonstrated by the ROOT5 toy. Sparse background bins must be included in the Poisson objective. With the requested fallback, a negative polynomial background is replaced through a new fit by a physically nonnegative constant.',
        'Agreement of SE/ME purity is assessed rather than imposed. Some points can remain different because of statistical limitations, different signal/background weights, parameter-bound behavior and model inadequacy. High-statistics structured residuals are not cured merely by choosing a Poisson objective.',
        'The report preserves all numerical trials, the previous chi-square PDF, the fixed example identities and all fit/selection metadata. The new report and new CF figures supersede the corresponding chi-square estimates; the old QA/LL material remains in the preserved reference PDF.'])]
    for c in examples:
        bg=lambda s:'C' if c[s]['fit'].get('constantBackgroundFallback') else 'Q'
        sections[3]['tables'][0]['rows'].append([str(c['page']),c['method'][-1],c['centLabel']+' / '+interval(c['low'],c['high']),
            measured(c['chiSquareVersion']['puritySE']),measured(c['puritySE']),measured(c['purityME']),bg('SE')+'/'+bg('ME')])
        def dndf(side):
            f=c[side]['fit'];return number(f['deviance'],5)+'/'+str(f['ndf']) if f['deviance'] is not None else 'NA'
        sections[3]['tables'][1]['rows'].append([str(c['page']),str(c['SE']['fitRangeObservations']),str(c['poissonDiagnostic']['SE']['zeroBins']),dndf('SE'),dndf('ME')])
    for candidate in audit['candidates']:
        for g in candidate['groups']:
            valid=[m[-1] for m,r in g['methods'].items() if r['normalization']['valid']]
            sections[4]['tables'][0]['rows'].append([interval(*candidate['normalization']),g['label'],', '.join(valid) or 'NA'])
    for g in results['groups']:
        m=g['methods'];sections[4]['tables'][1]['rows'].append([g['label'],str(m['method3']['valid_plot_bins']),str(m['method5']['valid_plot_bins']),str(m['method6']['valid_plot_bins']),measured(m['method5']['integrated']['cf']),measured(m['method6']['integrated']['cf'])])
    figures=[]
    for method in ('method3','method5'):
        figures.append(dict(path=str(out/'comparison'/('purity_fine_'+method+'_L.pdf')),caption='Updated '+method+' purity: Poisson L',paragraphs=['Red: SE (true); blue: ME (mix). Negative primary backgrounds trigger a nonnegative constant refit. Invalid estimates are omitted; no equality of SE/ME is imposed.']))
        figures.append(dict(path=str(out/'comparison'/('purity_fine_'+method+'_comparison.pdf')),caption=method+': fixed-bin chi-square/L comparison',paragraphs=['Filled: final L; open: preceding chi-square. The displayed physical-purity range is common across centralities; exact errors and rejected-fit diagnostics are available in the tables and JSON.']))
    figures.extend([dict(path=str(out/'comparison/cf_methods_3_5_6.pdf')),
                    dict(path=str(out/'comparison/cf_objective_comparison.pdf')),
                    dict(path=str(out/'comparison/purity_se_me.pdf'))])
    for c in examples:
        path=out/'cases'/('fit_%s_%s_%03d.pdf'%(c['method'],c['group'],c['bin']))
        types=['constant b>=0' if c[s]['fit'].get('constantBackgroundFallback') else 'quadratic' for s in ('SE','ME')]
        figures.append(dict(path=str(path),caption='Fixed case '+str(c['page'])+': '+c['method']+', '+c['centLabel']+', '+interval(c['low'],c['high'])+' GeV/c',paragraphs=[
            'Actual k* fit window '+interval(c['windowLow'],c['windowHigh'])+'. Background SE/ME: '+' / '.join(types)+'. P(SE)='+measured(c['puritySE'])+'; P(ME)='+measured(c['purityME'])+'.',
            'The lower-left and lower-middle panels include empty bins as signed Poisson deviance residuals. The normalized mass shapes are unchanged input data; the fitted background, purity and residuals are re-evaluated.']))
    content=dict(stem='summary_4He_KFParticle_3p85_PoissonL_20261102',sections=sections,figures=figures)
    (out/'report_content.json').write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(sections=len(sections),figureRecords=len(figures),cases=len(examples))))


if __name__=='__main__':main()
