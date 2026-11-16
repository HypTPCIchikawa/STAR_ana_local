#!/usr/bin/env python3
"""Assemble an auditable addendum and a staged update of a cited CF report."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from build_method6_report import esc, number, pdf_pages


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def include(path, pages='-'):
    return '\\includepdf[pages={'+str(pages)+'},fitpaper=true,pagecommand={}]{'+str(path)+'}\n'


def compile_pdf(tex, output):
    for unused in range(2):
        result = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                                 '-output-directory', str(output), str(tex)],
                                cwd=str(output), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True)
        (output/(tex.stem+'.build.log')).write_text(result.stdout)
        if result.returncode:
            raise RuntimeError('PDF compilation failed: '+tex.name)
    pdf = output/(tex.stem+'.pdf')
    return dict(path=str(pdf), sha256=digest(pdf), bytes=pdf.stat().st_size,
                pages=pdf_pages(pdf))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--content-json', type=Path, required=True,
                        help='Explicit prose, table and ordered figure specification.')
    parser.add_argument('--sources', type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve(strict=True)
    content = json.loads(args.content_json.read_text())
    sources = json.loads(args.sources.read_text())
    for source in sources['files'].values():
        if digest(source['path']) != source['sha256']:
            raise ValueError('Cited source changed: '+source['path'])
    original = Path(sources['files']['updated/summary_4He_KFParticle_3p85.pdf']['path'])
    original_pages = pdf_pages(original)
    figures = []
    for item in content['figures']:
        path = Path(item['path']).resolve(strict=True)
        figures.append(dict(path=str(path), sha256=digest(path), pages=pdf_pages(path)))
    stem = content['stem']
    if not stem or Path(stem).name != stem:
        raise ValueError('Report stem must be a simple filename stem.')
    if (out/(stem+'.pdf')).exists():
        raise ValueError('Refusing an existing staged report.')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9)))
    preamble = r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry,graphicx,amsmath,booktabs,longtable,pdfpages,hyperref}
\geometry{margin=1.7cm}
\setlength{\parindent}{0pt}\setlength{\parskip}{0.45em}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\begin{document}
'''
    header = r'\begin{center}{\Large $\Lambda$--$^{4}$He CF: '+esc(content['title'])+r'}\end{center}'+'\n'
    header += r'\textbf{Updated:} '+esc(now.strftime('%Y-%m-%d %H:%M JST'))+r'\quad\textbf{Sample:} 2021 P24iy.'+'\n'
    header += r'\textbf{SUMS:} \texttt{04BF2DAE47164B6F040E3A460E539289}.'+'\n'
    text = [preamble, header]
    for section in content['sections']:
        if section.get('newpage'):
            text.append('\\clearpage\n')
        text.append('\\section*{'+esc(section['title'])+'}\n')
        for paragraph in section.get('paragraphs', []):
            text.append(esc(paragraph)+'\n\n')
        for equation in section.get('equations', []):
            text.append('\\['+equation+'\\]\n')
        for table in section.get('tables', []):
            columns = table['columns']
            text.append(r'{\small\setlength{\tabcolsep}{3pt}\begin{longtable}{'+
                        'l'*len(columns)+r'}\toprule'+'\n')
            text.append(' & '.join(esc(v) for v in columns)+r'\\\midrule\endhead'+'\n')
            for row in table['rows']:
                if len(row) != len(columns):
                    raise ValueError('Table row width differs.')
                text.append(' & '.join(esc(v) for v in row)+r'\\'+'\n')
            text.append(r'\bottomrule\end{longtable}}'+'\n')
    text.append('\\clearpage\n')
    for item in content['figures']:
        if item.get('caption'):
            text.append('\\section*{'+esc(item['caption'])+'}\n')
            for paragraph in item.get('paragraphs', []):
                text.append(esc(paragraph)+'\n\n')
            # Keep the caption and the scientific figure on the same sheet.
            text.append('\\begin{center}\\includegraphics[page='+str(item.get('page',1))+
                        ',width=\\textwidth,height=.73\\textheight,keepaspectratio]{'+
                        str(Path(item['path']).resolve())+'}\\end{center}\\clearpage\n')
        else:
            text.append(include(Path(item['path']).resolve(), item.get('pages','-')))
    text.append('\\end{document}\n')
    addendum_tex = out/(stem+'_addendum.tex')
    addendum_tex.write_text(''.join(text))
    addendum = compile_pdf(addendum_tex, out)
    # The cover explicitly identifies which reference results the addendum replaces.
    full = [preamble,header,r'\section*{Reading this updated report}'+'\n']
    for paragraph in content['cover']:
        full.append(esc(paragraph)+'\n\n')
    full += ['The preserved reference report occupies the next '+str(original_pages)+
             ' pages. The diagnostic addendum begins on PDF page '+str(original_pages+2)+'.\n',
             '\\clearpage\n',include(original),include(Path(addendum['path'])),
             '\\end{document}\n']
    tex = out/(stem+'.tex')
    tex.write_text(''.join(full))
    assembled = compile_pdf(tex, out)
    if assembled['pages'] != 1+original_pages+addendum['pages']:
        raise ValueError('Unexpected combined report page count.')
    for source in sources['files'].values():
        if digest(source['path']) != source['sha256']:
            raise ValueError('Source changed during assembly.')
    save(out/'diagnostic_report_metadata.json',dict(
        phase='staged', builtAtJst=now.isoformat(), originalPages=original_pages,
        originalReport=sources['files']['updated/summary_4He_KFParticle_3p85.pdf'],
        addendum=addendum, report=assembled, figures=figures,
        sourceManifest=dict(path=str(args.sources.resolve()),sha256=digest(args.sources)),
        content=dict(path=str(args.content_json.resolve()),sha256=digest(args.content_json)),
        sourceCode=dict(path=str(Path(__file__).resolve()),sha256=digest(__file__))))
    print(json.dumps(dict(report=assembled,addendum=addendum)))


if __name__ == '__main__':
    main()
