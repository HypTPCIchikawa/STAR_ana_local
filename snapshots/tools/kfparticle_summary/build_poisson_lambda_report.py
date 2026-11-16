#!/usr/bin/env python3
"""Build a standalone Poisson-L purity/CF report from explicit saved content."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from build_method6_report import esc, pdf_pages
from build_purity_diagnostic_update import compile_pdf, include


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    content=json.loads((out/'report_content.json').read_text())
    context=json.loads((out/'poisson_context.json').read_text())
    stem=content['stem'];target=out/(stem+'.pdf')
    if target.exists():raise FileExistsError(target)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9)))
    text=[r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}
\usepackage{geometry,graphicx,amsmath,booktabs,longtable,pdfpages,hyperref}
\geometry{margin=1.7cm}\setlength{\parindent}{0pt}\setlength{\parskip}{0.4em}
\emergencystretch=3em
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\begin{document}
\begin{center}{\Large $\Lambda$--$^{4}$He: Poisson-likelihood purity and CF update}\end{center}
''',esc(now.strftime('%Y-%m-%d %H:%M JST'))+'\n\n']
    for section in content['sections']:
        if section.get('newpage'):text.append('\\clearpage\n')
        text.append('\\section*{'+esc(section['title'])+'}\n')
        for p in section.get('paragraphs',[]):text.append(esc(p)+'\n\n')
        for equation in section.get('equations',[]):text.append('\\['+equation+'\\]\n')
        for table in section.get('tables',[]):
            columns=table['columns']
            text.append(r'{\small\setlength{\tabcolsep}{3pt}\begin{longtable}{'+'l'*len(columns)+r'}\toprule'+'\n')
            text.append(' & '.join(esc(x) for x in columns)+r'\\\midrule\endhead'+'\n')
            for row in table['rows']:
                assert len(row)==len(columns)
                text.append(' & '.join(esc(x) for x in row)+r'\\'+'\n')
            text.append(r'\bottomrule\end{longtable}}'+'\n')
    text.append('\\clearpage\n');figures=[]
    for item in content['figures']:
        path=Path(item['path']).resolve(strict=True)
        figures.append(dict(path=str(path),sha256=sha(path),pages=pdf_pages(path)))
        if item.get('caption'):
            text.append('\\section*{'+esc(item['caption'])+'}\n')
            for p in item.get('paragraphs',[]):text.append(esc(p)+'\n\n')
            text.append('\\begin{center}\\includegraphics[page='+str(item.get('page',1))+
                ',width=\\textwidth,height=.72\\textheight,keepaspectratio]{'+str(path)+
                '}\\end{center}\\clearpage\n')
        else:text.append(include(path,item.get('pages','-')))
    text.append('\\end{document}\n')
    tex=out/(stem+'.tex');tex.write_text(''.join(text))
    report=compile_pdf(tex,out)
    metadata=dict(phase='staged',layout='standalone_poisson',builtAtJst=now.isoformat(),
        report=report,originalReport=context['targetReport'],figures=figures,
        content=dict(path=str(out/'report_content.json'),sha256=sha(out/'report_content.json')),
        tex=dict(path=str(tex),sha256=sha(tex)),
        sourceCode=dict(path=str(Path(__file__).resolve()),sha256=sha(__file__)))
    (out/'poisson_report_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':main()
