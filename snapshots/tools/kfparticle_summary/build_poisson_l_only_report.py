#!/usr/bin/env python3
"""Build a current-results report containing only audited Poisson-L figures."""
import argparse
import copy
import datetime
import hashlib
import json
from pathlib import Path

from build_method6_report import esc, pdf_pages
from build_purity_diagnostic_update import compile_pdf


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--target-report', type=Path, required=True)
    parser.add_argument('--report-stem', required=True)
    args = parser.parse_args()
    study = args.study_dir.resolve(strict=True)
    out = args.output_dir.resolve()
    if out.exists():
        raise FileExistsError(out)
    if not args.report_stem or Path(args.report_stem).name != args.report_stem:
        raise ValueError('A simple report stem is required.')
    target = args.target_report.resolve(strict=True)
    previous = json.loads((study/'poisson_report_metadata.json').read_text())
    scientific = json.loads((study/'validation.json').read_text())
    content = copy.deepcopy(json.loads((study/'report_content.json').read_text()))
    assert scientific['phase'] == 'validated'
    assert scientific['reportSha256'] == previous['report']['sha256'] == sha(target)
    source_hashes = {}
    for name in ['report_content.json', 'validation.json', 'poisson_report_metadata.json',
                 'cf_points.tsv', 'purity_points.tsv']:
        path = study/name
        source_hashes[str(path)] = sha(path)
    fit_records = {}
    for kind in ['group', 'fine']:
        path = study/(kind+'_extract/purity_fit_review.json')
        data = json.loads(path.read_text())
        assert data['fitOptions'] == 'SLQNR0'
        settings = study/(kind+'_extract/extract_settings.txt')
        assert 'fitOptions SLQNR0' in settings.read_text().splitlines()
        fit_records[kind] = dict(path=str(path), sha256=sha(path), fitOptions=data['fitOptions'],
                                 settingsPath=str(settings), settingsSha256=sha(settings))
    sections = content['sections']
    sections[0]['title'] = 'Current analysis: Poisson likelihood L'
    sections[0]['paragraphs'][2] = (
        'All 700 grouped and 1,260 fine-centrality spectrum cells were re-evaluated with L. '
        'Every fitted curve, purity point and correlation function shown in this report comes '
        'from this adopted L analysis. The nine diagnostic examples use the same fixed '
        'centralities and k* cells as the prior inspection.')
    sections[0]['paragraphs'][3] = (
        'The common CF normalization is 0.40-0.50 GeV/c for all five centrality selections '
        'and Method3/5/6. The alternative 0.30-0.40 interval is evaluated with the same L fits.')
    sections[1]['paragraphs'][0] = (
        'Empty mass bins participate in the Poisson objective. Their deviance contribution '
        'is 2*mu. A ROOT5 constant-background check using nine empty bins and one bin with '
        'one count gives the Poisson estimate and deviance below.')
    sections[1]['tables'][0]['rows'] = [
        row for row in sections[1]['tables'][0]['rows'] if not row[0].startswith('Chi-square')]
    sections[2]['paragraphs'][1] = (
        'The adopted fit uses an explicit 50,000-function-call budget and a positive '
        'sideband-based starting background. Primary and conditional constant-background '
        'fits both use SLQNR0. All displayed results use this same adopted numerical setup.')
    sections[2]['paragraphs'][2] = (
        'Deviance and signed deviance residuals are independently recomputed from saved '
        'predictions, including empty bins. The reported D/NDF is a diagnostic, not a '
        'calibrated sparse-count p-value.')
    cases = sections[3]
    cases['title'] = 'Nine diagnostic examples: L fit results'
    cases['paragraphs'][0] = (
        'All numbers in these tables use ROOT SLQNR0. NA denotes a failed or unusable '
        'estimate, not zero purity. Q denotes a quadratic background and C a nonnegative '
        'constant refit. Counts and mass projections are unchanged.')
    cases['paragraphs'][1] = (
        'The following Method3/Method5 purity figures show the L estimates alone: red '
        'for SE (true), blue for ME (mix). The representative mass spectra also show only '
        'L-fitted models. No agreement between SE and ME is imposed.')
    old_column = cases['tables'][0]['columns'].index('Old P(SE)')
    cases['tables'][0]['columns'].pop(old_column)
    for row in cases['tables'][0]['rows']:
        row.pop(old_column)
    sections[4]['title'] = 'L-derived normalization and correlation functions'
    sections[4]['paragraphs'][2] = (
        'All CF points and both normalization checks use the adopted L-derived purities. '
        'The plotted CFs use 0.40-0.50 GeV/c. Missing fit contributions mask their CF bins; '
        'invalid fits are never replaced with zero background or interpolated purity.')
    sections[5]['paragraphs'][-1] = (
        'The report presents only the adopted Poisson-L results. Prior reports, numerical '
        'trials and objective-comparison plots remain archived separately. The fitted model '
        'and fit quality are recorded for each SE and ME spectrum.')
    excluded = []
    figures = []
    for item in content['figures']:
        source = Path(item['path']).resolve(strict=True)
        if source.name.endswith('_comparison.pdf'):
            excluded.append(dict(path=str(source), sha256=sha(source), reason='contains non-L results'))
            continue
        if source.parent == study/'comparison':
            if source.name.startswith('purity_fine_method'):
                label = 'Method '+('3' if 'method3' in source.name else '5')+' purity, nine centralities'
                captions = [label]
                paragraphs = ['Red: SE (true); blue: ME (mix). All points use L fits. '
                              'Invalid estimates are omitted. Negative primary backgrounds '
                              'trigger a nonnegative constant-background L refit.']
            elif source.name == 'cf_methods_3_5_6.pdf':
                captions = ['Method3/5/6 CF: two centrality selections',
                            'Method3/5/6 CF: three centrality selections']
                paragraphs = ['All CF points use the adopted L-derived purities; normalization '
                              '0.40-0.50 GeV/c. Method6 shares Method5 ME purity and its error '
                              'between SE and ME.']
            elif source.name == 'purity_se_me.pdf':
                captions = ['SE/ME purity: two centrality selections',
                            'SE/ME purity: three centrality selections']
                paragraphs = ['All purity points use L fits. Method3 and Method5 fit SE and ME '
                              'separately. Method6 applies the same Method5 ME purity and error '
                              'to both.']
            else:
                raise ValueError('Unreviewed comparison figure: '+str(source))
        elif source.parent == study/'cases':
            captions = [item['caption'].replace('Fixed case ', 'Mass-spectrum case ')]
            paragraphs = item.get('paragraphs', [])
        else:
            raise ValueError('Unreviewed figure: '+str(source))
        assert len(captions) == pdf_pages(source)
        for page, caption in enumerate(captions, 1):
            figures.append(dict(path=str(source), page=page, fitOptions='SLQNR0',
                                caption=caption+' | Poisson L (SLQNR0)', paragraphs=paragraphs,
                                sha256=sha(source)))
    assert len(excluded) == 3 and len(figures) == 15
    content.update(stem=args.report_stem, sections=sections, figures=figures,
                   resultPolicy='Poisson L only', fitOptions='SLQNR0')
    out.mkdir(parents=True)
    save(out/'report_content.json', content)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9)))
    text = [r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}
\usepackage{geometry,graphicx,amsmath,booktabs,longtable,hyperref,fancyhdr}
\geometry{margin=1.7cm,headheight=14pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{0.4em}\emergencystretch=3em
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small $\Lambda$--$^4$He: adopted results}
\fancyhead[R]{\small Poisson likelihood L (SLQNR0)}\fancyfoot[C]{\thepage}
\begin{document}
\begin{center}{\Large $\Lambda$ purity and $\Lambda$--$^4$He CF}\\[0.4em]
{\Large\bfseries Poisson likelihood L only (SLQNR0)}\end{center}
''', esc(now.strftime('%Y-%m-%d %H:%M JST'))+'\n\n',
        r'\fbox{\parbox{0.95\linewidth}{\textbf{Adopted results only.} All mass fits use ROOT '
        r'\texttt{SLQNR0}, including option \texttt{L}. All purity and CF points come from these '
        r'fits. No non-L estimates or non-L comparison overlays are included.}}'+'\n\n']
    for section in sections:
        if section.get('newpage'):
            text.append('\\clearpage\n')
        text.append('\\section*{'+esc(section['title'])+'}\n')
        for paragraph in section.get('paragraphs', []):
            text.append(esc(paragraph)+'\n\n')
        for equation in section.get('equations', []):
            text.append('\\['+equation+'\\]\n')
        for table in section.get('tables', []):
            columns = table['columns']
            text.append(r'{\small\setlength{\tabcolsep}{3pt}\begin{longtable}{'+'l'*len(columns)+r'}\toprule'+'\n')
            text.append(' & '.join(esc(x) for x in columns)+r'\\\midrule\endhead'+'\n')
            for row in table['rows']:
                assert len(row) == len(columns)
                text.append(' & '.join(esc(x) for x in row)+r'\\'+'\n')
            text.append(r'\bottomrule\end{longtable}}'+'\n')
    text.append('\\clearpage\n')
    for index, item in enumerate(figures, 1):
        text.append('\\section*{'+esc('Figure '+str(index)+'. '+item['caption'])+'}\n')
        for paragraph in item['paragraphs']:
            text.append(esc(paragraph)+'\n\n')
        text.append('\\begin{center}\\includegraphics[page='+str(item['page'])+
                    ',width=\\textwidth,height=.70\\textheight,keepaspectratio]{'+item['path']+
                    '}\\end{center}\\clearpage\n')
    text.append('\\end{document}\n')
    tex = out/(args.report_stem+'.tex')
    tex.write_text(''.join(text))
    report = compile_pdf(tex, out)
    for path, digest in source_hashes.items():
        assert sha(path) == digest
    for item in figures:
        assert sha(item['path']) == item['sha256']
        assert item['fitOptions'] == 'SLQNR0' and '_comparison.pdf' not in item['path']
    metadata = dict(phase='staged', layout='poisson_L_only', builtAtJst=now.isoformat(),
                    report=report, originalReport=dict(path=str(target), sha256=sha(target)),
                    figures=figures, excludedFigures=excluded, fitOptions='SLQNR0',
                    content=dict(path=str(out/'report_content.json'), sha256=sha(out/'report_content.json')),
                    tex=dict(path=str(tex), sha256=sha(tex)),
                    sourceCode=dict(path=str(Path(__file__).resolve()), sha256=sha(__file__)),
                    scientificStudy=dict(path=str(study), sources=source_hashes, fits=fit_records))
    save(out/'poisson_report_metadata.json', metadata)
    validation = dict(phase='validated', reportSha256=report['sha256'], reportPages=report['pages'],
                      fitOptions='SLQNR0', figurePages=len(figures), nonLFiguresIncluded=0,
                      excludedComparisonFiles=len(excluded), scientificResultsUnchanged=True,
                      sourceHashesVerified=True, scientificValidation=scientific,
                      scientificValidationPath=str(study/'validation.json'),
                      scientificValidationSha256=sha(study/'validation.json'))
    save(out/'validation.json', validation)
    print(json.dumps(dict(report=report, figurePages=len(figures), nonLFiguresIncluded=0)))


if __name__ == '__main__':
    main()
