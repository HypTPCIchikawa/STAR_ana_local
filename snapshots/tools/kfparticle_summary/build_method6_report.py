#!/usr/bin/env python3
"""Build the Method3/5/6 comparison report from audited saved figures and JSON."""
import argparse
import collections
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess

import yaml

from sample_provenance import resolve_sample

OWNER = 'tools/kfparticle_summary/build_method6_report.py'
METHOD_NAMES = {'method3': 'Method 3', 'method5': 'Method 5', 'method6': 'Method 6'}


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, data):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def esc(value):
    chars = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
             '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}',
             '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    return ''.join(chars.get(c, c) for c in str(value))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def number(value, digits=5):
    return format(value, '.%dg' % digits) if finite(value) else 'NA'


def measured(value, tex=False):
    if not value or not finite(value.get('value')) or not finite(value.get('error')):
        return 'NA'
    if tex:
        return '$%s\\pm%s$' % (number(value['value']), number(value['error'], 3))
    return '%s ± %s' % (number(value['value']), number(value['error'], 3))


def pdf_pages(path):
    if path.read_bytes()[:5] != b'%PDF-':
        raise ValueError('Invalid or incomplete PDF: '+str(path))
    gs = shutil.which('gs')
    if not gs:
        raise RuntimeError('Ghostscript is required for PDF page verification.')
    result = subprocess.run([gs, '-q', '-dNODISPLAY', '-dSAFER',
        '--permit-file-read='+str(path), '-sPDFFile='+str(path), '-c',
        'PDFFile (r) file runpdfbegin pdfpagecount = quit'],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
        timeout=60, check=True)
    return int(result.stdout.strip().splitlines()[-1])


def record(path, pages=None):
    result = {'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size}
    if pages is not None:
        actual = pdf_pages(path)
        if actual != pages:
            raise ValueError('%s: expected %d pages, found %d.' % (path.name, pages, actual))
        result['pages'] = actual
    return result


def inspected_rows(group, method, side, limits):
    return [row for row in group['methods'][method][side]
            if limits[0] <= row['x'] < limits[1]]


def fit_diagnostics(review, cfg):
    labels = {g['id']: g['label'] for g in cfg['groups']}
    output = []
    for group in review['groups']:
        if group['id'] not in labels:
            continue
        for method in ('method3', 'method5'):
            for side in ('SE', 'ME'):
                rows = inspected_rows(group, method, side, cfg['inspection_range'])
                chi = [r['fit']['chi2']/r['fit']['ndf'] for r in rows
                       if finite(r['fit'].get('chi2')) and r['fit'].get('ndf', 0)>0]
                obs = [r['fit']['observationsInStoredMassRange'] for r in rows]
                counts = collections.Counter(r['fit']['reason'] for r in rows)
                output.append({'group': group['id'], 'label': labels[group['id']],
                    'method': method, 'side': side, 'total': len(rows),
                    'fullValid': sum(bool(r['full']['valid']) for r in rows),
                    'numericallyValid': sum(bool(r['fit']['numericallyValid']) for r in rows),
                    'boundary': sum(bool(r['fit']['boundary']) for r in rows),
                    'observationsMin': min(obs) if obs else None,
                    'observationsMedian': statistics.median(obs) if obs else None,
                    'medianChi2Ndf': statistics.median(chi) if chi else None,
                    'maxChi2Ndf': max(chi) if chi else None,
                    'fitReasons': dict(counts)})
    return output


def purity(row):
    full = row['full']
    n, b, vb = full['windowNtot'], full.get('background'), full.get('backgroundVariance')
    if not full['valid'] or not finite(n) or n<=0 or not finite(b) or not finite(vb):
        return None
    variance = vb/(n*n)+b*b*full['windowError2']/(n*n*n*n)
    if variance<0 or not finite(variance):
        return None
    return {'value': 1-b/n, 'error': math.sqrt(variance)}


def representative_inputs(out, review, cfg):
    groups = {g['id']: g for g in review['groups']}
    selected = []
    for method in ('method3', 'method5'):
        for group_id in cfg['group_sets'][0]:
            group = groups[group_id]
            for index in cfg['representative_bins'][method]:
                if isinstance(index, bool) or not isinstance(index, int) or index<1:
                    raise ValueError('Representative bins must be positive integers.')
                rows = {side: next(r for r in group['methods'][method][side] if r['bin']==index)
                        for side in ('SE', 'ME')}
                name = 'fit_%s_%s_%03d.pdf' % (method, group_id, index)
                selected.append({'name': name, 'method': method, 'group': group_id,
                    'bin': index, 'x': rows['SE']['x'], 'rows': rows})
    return selected


def include_pdf(name, pages='-'):
    return r'\includepdf[pages={'+str(pages)+r'},fitpaper=true,pagecommand={}]{'+name+'}\n'


def tex_report(sample, review, data, cfg, base, diagnostics, representatives):
    labels = {g['id']: g['label'] for g in cfg['groups']}
    groups = {g['id']: g for g in data['groups']}
    date = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).strftime('%Y-%m-%d %H:%M JST')
    native_width = review['originalKstarEdges'][1]-review['originalKstarEdges'][0]
    fit = review['settings']
    n_events = sample['eventCounters']['inputEventsFromRefMultEntries']
    n_reco = sample['eventCounters']['reconstructedEventsFromRefMultCorrEntries']
    parts = [r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry,graphicx,amsmath,booktabs,longtable,pdfpages,hyperref}
\geometry{margin=1.8cm}
\setlength{\parindent}{0pt}\setlength{\parskip}{0.4em}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\title{$\Lambda$--$^{4}\mathrm{He}$: Purity checks and Method 3/5/6 comparison}
\author{STAR analysis working report}
\date{'''+esc(date)+r'''}
\begin{document}\maketitle
\section{Sample and calculation}
''',
        'Au+Au fixed-target beam-energy label '+esc(sample['beam'])+' GeV; '+str(sample['year'])+
        r' P24iy. The beam label is not a statement about $\sqrt{s_{NN}}$. '+
        'The input contains '+format(sample['included'], ',')+'/'+format(sample['expected'], ',')+
        ' submitted subjobs; '+str(sample['expected']-sample['included'])+' are excluded.\n',
        r'\textbf{SUMS:} \texttt{'+esc(sample['jobid'])+r'}. '+
        r'\textbf{ROOT:} \path{'+Path(sample['inputRoot']).name+'}.\n',
        'Measured input/reconstructed events: '+format(n_events, ',.0f')+' / '+format(n_reco, ',.0f')+
        '. Existing ROOTs, PID coefficients, Maker selection and mixing are unchanged.\n',
        r'\begin{center}\begin{tabular}{lll}\toprule Method & Purity source & Application\\\midrule'+'\n',
        r'3 & SE and ME separately & Per-CF-bin mass-fit signal $N-B$\\'+'\n',
        r'5 & SE and ME separately & Coarse sliding-window purity on native counts\\'+'\n',
        r'6 & ME only & Same ME purity and error on both SE and ME\\\bottomrule\end{tabular}\end{center}'+'\n',
        'Native bins: '+number(native_width*1000)+' MeV/$c$; CF bins: '+number(native_width*fit['rebin']*1000)+
        ' MeV/$c$; Method 5/6 purity cells: '+number(native_width*fit['purityRebin']*1000)+
        ' MeV/$c$. The sliding window uses the center and '+str(fit['slidingHalfWidth'])+
        ' neighbors on either side, clipped at the stored endpoints.\n',
        'CF display: '+number(base['plot_x_range'][0])+'--'+number(base['plot_x_range'][1])+
        ' GeV/$c$; normalization: '+number(base['normalization'][0])+'--'+number(base['normalization'][1])+
        ' GeV/$c$. The normalization interval remains in use outside the displayed CF range. '+
        'The purity plots and fit inspection include '+number(cfg['inspection_range'][0])+'--'+number(cfg['inspection_range'][1])+' GeV/$c$.\n',
        r'The mass model is Gaussian plus quadratic background for the \emph{$\Lambda$ invariant mass}. '+
        'This does not modify the accepted helium dE/dx PID model. No mass-model shape or acceptance threshold is relaxed to obtain successful fits.\n',
        r'\subsection*{Centrality selections}'+'\n']
    for index, ids in enumerate(cfg['group_sets']):
        parts.append(str(len(ids))+' selections: '+', '.join(esc(labels[g]) for g in ids)+'.\n')
    parts += [r'\clearpage\section{Mass-fit reliability checks}'+'\n',
        'The table inspects fit cells with centers in the stated inspection range, including normalization. '+
        'Accepted means the full-window background is physically usable under the numerical checks; it does not establish model adequacy. '+
        'A boundary warning is recorded separately. The minimum and median observation counts refer to the entire fitted mass projection, not only its signal window or central purity cell.\n',
        r'{\scriptsize\setlength{\tabcolsep}{3pt}\begin{longtable}{llcrrrrrr}\toprule Centrality & Method & Side & Accepted/total & Boundary & Min obs. & Median obs. & Median $\chi^2/\nu$ & Max $\chi^2/\nu$\\\midrule\endhead'+'\n']
    for d in diagnostics:
        parts.append(' & '.join([esc(d['label']), d['method'][-1], d['side'],
            '%d/%d'%(d['fullValid'], d['total']), str(d['boundary']), number(d['observationsMin'], 4),
            number(d['observationsMedian'], 4), number(d['medianChi2Ndf'], 4), number(d['maxChi2Ndf'], 4)])+r'\\'+'\n')
    parts += [r'\bottomrule\end{longtable}}'+'\n',
        'Sparse Method 3 SE spectra can have poorly determined backgrounds even when the minimizer returns status zero. '+
        'For Method 5, a low-count central cell is fitted using a broader sliding projection: the central-cell count alone is not the fit sample size. '+
        r'Large $\chi^2/\nu$ or small fit probabilities remain a separate problem at high statistics. '+
        'The atlas retains all projections, including failures and boundary solutions, with model curves, pulls and normalized SE/ME mass shapes.\n',
        'Purity differences alone do not identify the cause. The figures permit checking low statistics, nonphysical fitted background, covariance failure, and persistent mass-shape/model discrepancy separately. '+
        'No numerical goodness threshold is silently used as a new selection cut.\n',
        r'\subsection*{Representative low-$k^*$ and normalization fits}'+'\n',
        r'{\scriptsize\begin{tabular}{llcrrrr}\toprule Centrality & Method & Cell & $k^*$ & $P_{\rm SE}$ & $P_{\rm ME}$ & SE/ME fitted observations\\\midrule'+'\n']
    for rep in representatives:
        se, me = rep['rows']['SE'], rep['rows']['ME']
        parts.append(' & '.join([esc(labels[rep['group']]), rep['method'][-1], str(rep['bin']),
            number(rep['x']), measured(purity(se), True), measured(purity(me), True),
            number(se['fit']['observationsInStoredMassRange'],4)+' / '+number(me['fit']['observationsInStoredMassRange'],4)])+r'\\'+'\n')
    parts += [r'\bottomrule\end{tabular}}'+'\n',
        r'\clearpage\section{Meaning and uncertainty of Method 6}'+'\n',
        r'SE and ME purity equality is an assumption to test, not a requirement. Let the mixed-event signal/background be $S_m,B_m$, and let the same-event components be proportional to $C_sS_m,C_bB_m$ with a common scale. Then'+'\n',
        r'\[P_{\rm SE}=\frac{C_sP_{\rm ME}}{C_sP_{\rm ME}+C_b(1-P_{\rm ME})}.\]'+'\n',
        'Thus different signal and background correlations can produce different purities even with arbitrarily large statistics. '+
        'Method 6 is the requested common-ME-purity transfer comparison; the agreement imposed by its construction is not evidence that the transfer is unbiased.\n',
        r'For each native cell, $A=P_{\rm ME}N_{\rm SE}$ and $B=P_{\rm ME}N_{\rm ME}$ use the same random purity estimate, including the same ME-derived error. '+
        r'For a ratio wholly inside one purity cell, $A/B=N_{\rm SE}/N_{\rm ME}$ and the local purity uncertainty cancels. '+
        'An output bin crossing two purity cells, and the weighted normalization integrals, can retain purity dependence. '+
        'Purity is applied before rebinning; assigning one purity using the output-bin center would give a different answer.\n',
        'The calculation retains the shared-purity SE--ME covariance, normalization covariance, and covariance caused by reuse of observed counts in overlapping purity windows. '+
        'The saved full covariance matrices should be used instead of treating plotted CF points as independent. '+
        'Errors are first-order estimates using the supplied histogram variances and background-fit covariance.\n',
        r'\textbf{Not included:} background-fit/count covariance, cross-fit background-parameter covariance, missing event/track/mixed-pair reuse correlations, cross-method covariance, and the model/transfer systematic uncertainty of applying ME purity to SE. '+
        'Near zero observed counts, plug-in Poisson variances and Gaussian error propagation have limited coverage. '+
        'A smaller Method 6 error bar does not establish that its transfer assumption is correct. No independent-error significance is assigned to differences between methods.\n',
        r'\subsection*{Missing fits and normalization}'+'\n',
        'All normalization bins require usable positive SE and ME signal yields. If any contributing fit is unavailable, the entire group CF is NA; no smaller normalization interval is substituted. '+
        'A failed ME purity masks both Method 6 sides. Sharing ME purity cannot repair a failed ME normalization fit. '+
        'Missing values are not physical zero measurements.\n',
        r'\clearpage\section{CF comparison and figure guide}'+'\n',
        'The integrated CF below is the ratio of summed yields over '+number(cfg['comparison_range'][0])+'--'+number(cfg['comparison_range'][1])+
        ' GeV/$c$, multiplied by the separately evaluated normalization. It is not an arithmetic mean of plotted CF points. '+
        'The numbers use the present shared-input covariance propagation; their errors should not be read as independent errors across methods.\n',
        r'\begin{center}\begin{tabular}{lrrr}\toprule Centrality & Method 3 & Method 5 & Method 6\\\midrule'+'\n']
    for group in cfg['groups']:
        g = groups[group['id']]
        parts.append(esc(group['label'])+' & '+' & '.join(measured(g['methods'][m]['integrated']['cf'], True)
            for m in ('method3','method5','method6'))+r'\\'+'\n')
    parts += [r'\bottomrule\end{tabular}\end{center}'+'\n',r'\subsection*{Groups unavailable because of normalization fits}'+'\n']
    missing = []
    for group in cfg['groups']:
        for method in ('method3', 'method5', 'method6'):
            norm = groups[group['id']]['methods'][method]['normalization']
            if not norm['valid']:
                missing.append(esc(group['label'])+' / '+METHOD_NAMES[method]+': '+esc(norm['reason']))
    if missing:
        parts.append(r'{\small\begin{itemize}'+'\n'+''.join(r'\item '+line+'\n' for line in missing)+r'\end{itemize}}'+'\n')
    else:
        parts.append('All requested group/method normalization intervals are available in this calculation.\n')
    parts += [r'\subsection*{Pages that follow}'+'\n',
        'The next two original pages compare CFs for the two-way and three-way centrality selections; the next two compare SE/ME purities over the inspection range. '+
        'Method 6 uses identical ME purity/error on both sides by construction. '+
        'The final representative pages show the configured cells for both groups of the two-way selection, first Method 3 and then Method 5. '+
        'Full atlases remain separate so that every inspected fit is available without expanding this report to hundreds of pages.\n',
        r'\path{methods_comparison.json} contains CF values, covariance, purity and masks; '+
        r'\path{purity_fit_review.json} contains the fitted mass projections and fit parameters; '+
        r'\path{report_metadata.json} records dependencies and hashes. '+
        'The comparison does not refit or choose final interaction parameters.\n',
        include_pdf('cf_methods_3_5_6.pdf'), include_pdf('purity_se_me.pdf')]
    parts.extend(include_pdf(rep['name'], '1') for rep in representatives)
    parts.append(r'\end{document}'+'\n')
    return ''.join(parts)


def readme(sample, data, cfg, diagnostics, report_name, representatives, out, settings):
    labels = {g['id']:g['label'] for g in cfg['groups']}
    lines = ['# Method3・5・6：Λ純度fit点検とCF比較', '',
        '- [図付きレポート](%s)' % report_name,
        '- [CF比較：中央度2分割／3分割](cf_methods_3_5_6.pdf)',
        '- [SE/ME純度比較：規格化域を含む](purity_se_me.pdf)',
        '- [Method3 全fit atlas](atlas_method3.pdf)／[Method5 全fit atlas](atlas_method5.pdf)',
        '- 入力ROOT: [%s](%s)' % (Path(sample['inputRoot']).name, os.path.relpath(sample['inputRoot'], str(out))),
        '- 2021 P24iy、SUMS `%s`、%s / %s subjobs（除外%s件）。' %
            (sample['jobid'], format(sample['included'],','), format(sample['expected'],','), sample['expected']-sample['included']),
        '- 計算設定: `%s`。' % settings, '', '## 定義と誤差', '',
        '- Method3：CF binごとのSE/ME Λ質量fitから別々に信号N−Bを求める。',
        '- Method5：SE/MEそれぞれの粗格子・滑動窓fitの純度を元binへ適用してからCF binへ集約。',
        '- Method6：Method5のME純度とその誤差を、同一推定量としてSE/ME両方へ共通適用。',
        '- 同じ純度を独立な2つの誤差として加算しない。単一純度cell内の比では純度誤差が相殺し、格子跨ぎや規格化では依存が残り得る。',
        '- SE/ME同純度は検証すべき仮定。signal/BGの相関や質量形状が異なれば十分な統計でも純度差が生じる。Method6による一致は構成上のもので、仮定の検証にはならない。',
        '- ME→SEへの移行・モデル系統はME fitの統計誤差に含まれない。規格化や共有純度の共分散は保持するが、fit/count、別fit背景間、event/track/mixed-pair再利用、cross-methodの未取得共分散は未反映。',
        '- 規格化域の全binが必要。不成立fitがあればCF全体をNAとし、残りbinだけで規格化しない。ME規格化fitの失敗はMethod6へ移行しても回復しない。',
        '- Method3/5は前回と同じ定義・入力で比較するが、今回の誤差は規格化と共有入力の共分散を扱う近似。同じerror barの再現や方式差の独立有意度は主張しない。', '',
        '## fit診断（inspection range内）', '',
        '| 中央度 | Method | sample | usable/total | boundary | min / median observations | median / max χ²/NDF |',
        '|---|---|---|---:|---:|---:|---:|']
    for d in diagnostics:
        lines.append('| %s | %s | %s | %d/%d | %d | %s / %s | %s / %s |' %
            (d['label'], d['method'][-1], d['side'], d['fullValid'],d['total'],d['boundary'],
             number(d['observationsMin']),number(d['observationsMedian']),number(d['medianChi2Ndf']),number(d['maxChi2Ndf'])))
    lines += ['', '数値収束とモデル適合の良さは別の診断。Method3の低統計fit、Method5の実fit窓の統計量、両sampleのχ²・pull・規格化質量形状をatlasで確認する。中心cellの少なさだけでMethod5のfit全体が低統計とは断定しない。', '',
        '## 積分CF', '', '| 中央度 | Method3 | Method5 | Method6 |', '|---|---:|---:|---:|']
    groups = {g['id']:g for g in data['groups']}
    for group in cfg['groups']:
        values=[measured(groups[group['id']]['methods'][method]['integrated']['cf']) for method in ('method3','method5','method6')]
        lines.append('| '+group['label']+' | '+' | '.join(values)+' |')
    lines += ['', '## 規格化不成立', '']
    missing=[]
    for group in cfg['groups']:
        for method in ('method3','method5','method6'):
            norm=groups[group['id']]['methods'][method]['normalization']
            if not norm['valid']:missing.append('- %s / %s: %s' % (group['label'],METHOD_NAMES[method],norm['reason']))
    lines += missing or ['規格化区間は全群・全方式で成立。']
    lines += ['', '## レポートに収録した代表fit', '']
    for rep in representatives:
        lines.append('- [%s / %s / bin %d](%s)' %
            (labels[rep['group']], METHOD_NAMES[rep['method']], rep['bin'], rep['name']))
    lines += ['', '## 再生成', '', '```bash',
        'python3 tools/kfparticle_summary/build_method6_report.py \\',
        '  --output-dir '+str(out)+' \\', '  --settings '+str(settings), '```', '',
        'JSON・PDF・設定・生成codeのhashと頁数はreport_metadata.jsonへ保存。元ROOTとproduction条件は変更しない。', '']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--settings',type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve(strict=True);settings=args.settings.resolve(strict=True)
    cfg=yaml.safe_load(settings.read_text());base_path=(settings.parent/cfg['summary_settings']).resolve(strict=True)
    base=yaml.safe_load(base_path.read_text())
    paths=[out/'methods_comparison.json',out/'purity_fit_review.json']
    data,review=[json.loads(p.read_text()) for p in paths]
    sample=data['sample'];proof=sample['provenance']
    verified=resolve_sample('3p85',sample['inputRoot'],proof['outputManifest']['path'],proof['runmeta']['path'])
    if sample!=verified or Path(review['inputRoot']).resolve()!=Path(sample['inputRoot']).resolve():
        raise ValueError('Sample/input mismatch between comparison, mass review and submission evidence.')
    if data['group_sets']!=cfg['group_sets']:
        raise ValueError('Comparison centrality group sets differ from explicit settings.')
    ids=[g['id'] for g in cfg['groups']]
    if set(ids)!={g['id'] for g in data['groups']} or set(ids)!={g['id'] for g in review['groups']}:
        raise ValueError('Requested groups differ between settings, CF data and mass review.')
    for group in data['groups']:
        if group['normalization']!=base['normalization'] or group['plot_range']!=base['plot_x_range'] or group['comparison_range']!=cfg['comparison_range']:
            raise ValueError('Group ranges differ from the explicit calculation settings.')
    stem=cfg.get('report_stem','method6_cf_purity_report')
    if not isinstance(stem,str) or not stem or Path(stem).name!=stem or not all(c.isalnum() or c in '_-' for c in stem):
        raise ValueError('report_stem must be a safe plain basename.')
    name=stem+'.pdf';representatives=representative_inputs(out,review,cfg)
    atlas_pages={method:sum(len(inspected_rows(g,method,'SE',cfg['inspection_range'])) for g in review['groups'])
                 for method in ('method3','method5')}
    figures={}
    for file,pages in [('cf_methods_3_5_6.pdf',2),('purity_se_me.pdf',2),
            ('atlas_method3.pdf',atlas_pages['method3']),('atlas_method5.pdf',atlas_pages['method5'])]+[(rep['name'],1) for rep in representatives]:
        figures[file]=record(out/file,pages)
    sources={str(p):record(p) for p in paths+[settings,base_path,Path(__file__).resolve(),Path(__file__).with_name('sample_provenance.py').resolve()]}
    quality=out/'fit_quality.json'
    if quality.exists():sources[str(quality)]=record(quality)
    signature={'owner':OWNER,'inputRoot':sample['inputRoot'],'sample':sample,'settingsSha256':sha(settings)}
    metadata_path=out/'report_metadata.json'
    if metadata_path.exists():
        old=json.loads(metadata_path.read_text())
        if any(old.get(k)!=v for k,v in signature.items()):
            raise ValueError('Existing report belongs to another owner, sample or settings.')
    elif any((out/n).exists() for n in (name,stem+'.tex','README.md')):
        raise ValueError('Refusing to overwrite unowned report products.')
    diagnostics=fit_diagnostics(review,cfg)
    metadata=dict(signature,schemaVersion=1,phase='building',figures=figures,sources=sources,fitDiagnostics=diagnostics)
    save(metadata_path,metadata)
    try:
        tex=out/(stem+'.tex');tex.write_text(tex_report(sample,review,data,cfg,base,diagnostics,representatives))
        (out/'README.md').write_text(readme(sample,data,cfg,diagnostics,name,representatives,out,settings))
        logs=out/'logs';logs.mkdir(exist_ok=True)
        latex=shutil.which('pdflatex')
        if not latex:raise RuntimeError('pdflatex is required.')
        for iteration in (1,2):
            with (logs/('method6_pdflatex_%d.log'%iteration)).open('w') as handle:
                subprocess.run([latex,'-no-shell-escape','-interaction=nonstopmode','-halt-on-error','-file-line-error',tex.name],cwd=str(out),stdout=handle,stderr=subprocess.STDOUT,timeout=120,check=True)
        result=out/name
        pages=pdf_pages(result)
        expected=4+4+len(representatives)
        if pages<expected:
            raise RuntimeError('Report lost required pages: expected at least %d, got %d.'%(expected,pages))
        for path,rec in sources.items():
            if sha(Path(path))!=rec['sha256']:raise RuntimeError('Report source changed during generation: '+path)
        for file,rec in figures.items():
            if sha(out/file)!=rec['sha256']:raise RuntimeError('Report figure changed during generation: '+file)
        metadata.update(phase='complete',report=record(result,pages),texSha256=sha(tex),readmeSha256=sha(out/'README.md'),
            completedAtJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat())
        save(metadata_path,metadata)
        print('Created %s (%d pages)'%(result,pages))
    except Exception as exc:
        metadata.update(phase='failed',error=str(exc));save(metadata_path,metadata);raise


if __name__=='__main__':
    main()
