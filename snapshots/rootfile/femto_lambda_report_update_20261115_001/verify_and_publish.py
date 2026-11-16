#!/usr/bin/env python3
"""Check the staged Japanese report and replace the requested PDF atomically."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import re
import shutil
import subprocess
import zlib
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
STAGE=HERE/'pdf_build'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path))


def main():
    cfg=yaml.safe_load((HERE/'report_settings.yaml').read_text())
    pdf=STAGE/'cf_ep_rapidity_report_ja.pdf';target=ROOT/cfg['target_pdf'];out=target.parent
    count=int(subprocess.check_output(['gs','-q','-dNODISPLAY','-dNOSAFER',
        '-sPDFname='+str(pdf),'-c','PDFname (r) file runpdfbegin pdfpagecount = quit'],text=True).strip())
    assert count==33
    assert len(list(STAGE.glob('page-[0-9][0-9].png')))==count
    log=(STAGE/'build-2.log').read_text()
    warnings={word:log.count(word) for word in ('Missing character','Overfull','Underfull','Warning')}
    assert not any(warnings.values()),warnings
    text=STAGE/'cf_ep_rapidity_report_ja.txt'
    subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=txtwrite',
        '-sOutputFile='+str(text),str(pdf)],check=True)
    words=text.read_text()
    assert all(word in words for word in ['中央度','低統計','0.906','0.899','0.943'])
    raw=pdf.read_bytes();decoded=bytearray(raw)
    for match in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream',raw,re.S):
        try:decoded.extend(zlib.decompress(match.group(1)))
        except zlib.error:pass
    fonts={name:name.encode() in decoded for name in ('DroidSansFallback','/FontFile','/ToUnicode')}
    assert all(fonts.values()),fonts
    science=ROOT/cfg['centrality_study']
    identities={}
    for ds,entry in json.loads((science/'input_manifest.json').read_text())['inputs'].items():
        st=Path(entry['path']).stat();actual=dict(size=st.st_size,mtimeNs=st.st_mtime_ns,inode=st.st_ino,device=st.st_dev)
        assert all(actual[k]==entry[k] for k in actual)
        identities[ds]=dict(path=entry['path'],identity=actual,unchanged=True,
            check='size, mtime_ns, inode, device; full-file SHA not computed')
    independent=json.loads((science/'independent_validation.json').read_text());assert independent['status']=='PASS'
    plots=ROOT/cfg['plots_directory'];assert len(list(plots.glob('*.pdf')))==10
    plot_manifest=HERE/'summary_plot_manifest.json';assert plot_manifest.exists()
    archive=ROOT/cfg['archive_directory']
    archive_info=json.loads((archive/'archive_manifest.json').read_text())
    previous=next(x for x in archive_info['files'] if x['name']==target.name)
    assert sha(archive/target.name)==previous['sha256']
    assert sha(target)==previous['sha256'],'Requested target changed since the original version was archived'
    paths=[HERE/x for x in ('report_settings.yaml','centrality_pages.py','build_updated_report.py',
        'verify_and_publish.py','build_summary.json','centrality_scientific_summary.json',
        'centrality_scientific_summary.md','summary_display.yaml','summary_plots.log',
        'plot_input_provenance.json','summary_plot_manifest.json')]
    paths += [science/x for x in ('study.yaml','analysis_results.json','metrics.csv',
        'independent_validation.json','closure_checks.json','input_manifest.json')]
    paths += [ROOT/cfg['original_study']/'build_report.py']
    paths += [STAGE/'cf_ep_rapidity_report_ja.tex',STAGE/'build-2.log',text]+list(sorted(plots.glob('*.pdf')))
    manifest=dict(updatedAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        report=dict(path=str(target.relative_to(ROOT)),bytes=pdf.stat().st_size,sha256=sha(pdf)),pages=count,
        previous_version=dict(directory=cfg['archive_directory'],pages=22,pdf_sha256=previous['sha256']),
        settings=cfg,inputIdentityChecks=identities,
        checks=dict(pdf_pages=count,rendered_pages=count,typesetting_warning_counts=warnings,
            fonts=fonts,japanese_text_extraction=True,source_independent_validation=independent['status'],
            scientific_review='centrality_scientific_summary.json and .md',new_vector_plots=10,
            root_plot_input_validation='see summary_plot_manifest.json',
            visual_review='All 33 pages in contact sheet; individual pages 20 and 26 inspected; no clipping observed.'),
        production_code_update='Paused by user; this task changes report postprocessing only.',
        limitations=['Fixed comparison normalization is not an established physical plateau.',
            'kstar range ends at 1 GeV/c; no asymptotic convergence claim.',
            'Rapidity selections do not split ME pools; EP-plane-resolved kstar data are unavailable.',
            'Event reuse and EP1/EP8 covariance are not available; sparse Gaussian errors are diagnostic.'],
        artifacts=[record(path) for path in paths])
    # The archived report remains immutable; only validated derived outputs are published.
    temporary=out/'.cf_ep_rapidity_report_ja.updating.pdf';shutil.copy2(pdf,temporary)
    assert sha(temporary)==manifest['report']['sha256'];temporary.replace(target)
    for name in ('cf_ep_rapidity_report_ja.tex','cf_ep_rapidity_report_ja.txt','build-1.log','build-2.log','contact_sheet.png'):
        shutil.copy2(STAGE/name,out/name)
    for path in STAGE.glob('page-[0-9][0-9].png'):shutil.copy2(path,out/path.name)
    (out/'report_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (HERE/'publication_result.json').write_text(json.dumps(dict(report=manifest['report'],pages=count,archive=cfg['archive_directory']),indent=2)+'\n')
    print('JAPANESE_REPORT_UPDATED pages=33 sha256='+manifest['report']['sha256'])


if __name__=='__main__':main()
