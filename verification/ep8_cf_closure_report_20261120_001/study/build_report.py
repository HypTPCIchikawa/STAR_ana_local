#!/usr/bin/env python3
"""Build and check the Japanese report from verified local figure artifacts."""
import argparse,json,hashlib,subprocess,shutil,re
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--figure-directory',type=Path,required=True);a=ap.parse_args()
    f=a.figure_directory.resolve();src=f/'lambda_he4_ep8_cf_and_tree_closure_ja.tex'
    required=['computed_summary.tex','report_validation_note.tex','ep8_cf_c00_60_panels.pdf','ep8_cf_c00_60_highk_panels.pdf','ep8_seme_c00_60_panels.pdf','ep8_highk_c00_60_table.pdf','ep8_highk_primary_heatmap_v2.pdf','ep8_cf_six_centralities_overlay.pdf','closure_matched_seme.pdf','ep8_cf_nine_cent9_overlay.pdf']
    for name in required: assert (f/name).is_file() and (f/name).stat().st_size>0,name
    exe=shutil.which('xelatex');assert exe,'xelatex missing'
    for i in [1,2]:
        with (f/('report_xelatex_pass'+str(i)+'.log')).open('w') as out:
            subprocess.run([exe,'-no-shell-escape','-interaction=nonstopmode','-halt-on-error',src.name],cwd=str(f),stdout=out,stderr=subprocess.STDOUT,check=True)
    pdf=f/src.with_suffix('.pdf').name
    gs=shutil.which('gs');assert gs,'Ghostscript missing'
    log=src.with_suffix('.log').read_text(errors='replace')
    match=re.search(r'Output written on .*?\((\d+) pages?',log)
    assert match,'Missing XeLaTeX page count'
    pages=int(match.group(1));assert pages==12,('Unexpected report page count',pages)
    info='Pages: '+str(pages)+' (XeLaTeX log; independently rendered below)\n'
    (f/'report_pdfinfo.txt').write_text(info)
    subprocess.run([gs,'-q','-dBATCH','-dNOPAUSE','-sDEVICE=txtwrite',
        '-sOutputFile='+str(f/'report_extracted_text.txt'),str(pdf)],check=True)
    previews=f/'report_previews';previews.mkdir(exist_ok=True)
    subprocess.run([gs,'-q','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-r90',
        '-sOutputFile='+str(previews/'page_%02d.png'),str(pdf)],check=True)
    assert len(list(previews.glob('page_*.png')))==pages,'Rendered page count differs'

    text=(f/'report_extracted_text.txt').read_text()
    assert all(s in text for s in ['高','tree','0.899','1186','108']),'Missing report text'
    log=src.with_suffix('.log').read_text(errors='replace')
    assert 'Missing character:' not in log,'Missing glyphs'
    assert 'Overfull \\hbox' not in log and 'Overfull \\vbox' not in log,'Report content overflows'
    result=dict(generated=True,pdf=str(pdf),bytes=pdf.stat().st_size,sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),pdfinfo=info,figureFiles=required,visualReviewPending=True)
    (f/'report_build.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
