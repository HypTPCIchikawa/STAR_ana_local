#!/usr/bin/env python3
"""Validate ROOT-generated PDFs, input identities and their reproducibility records."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import re
import subprocess
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pages(path):
    return int(subprocess.check_output(['gs','-q','-dNODISPLAY','-dNOSAFER',
        '-sPDFname='+str(path),'-c','PDFname (r) file runpdfbegin pdfpagecount = quit'],text=True).strip())


def main():
    cfg=yaml.safe_load((HERE/'study.yaml').read_text());out=ROOT/cfg['canvas']['output_directory']
    index=json.loads((out/'pdf_index.json').read_text())
    files=[ROOT/index['combined_pdf']]+[ROOT/c['pdf'] for c in index['centralities']]
    counts=[pages(path) for path in files]
    assert counts==[121]+[11]*11,counts
    log=(HERE/'draw.log').read_text()
    assert 'RETURN_CODE 0' in log and 'DRAW_CANVAS_BATCH_OK files=12' in log
    assert not re.search(r'(^|\n)(Error in |SysError in |Fatal in |.*_FAILED:)',log)
    source=json.loads((HERE/'input_manifest.json').read_text());identities={}
    for ds,entry in source['inputs'].items():
        st=Path(entry['path']).stat()
        actual=dict(size=st.st_size,mtimeNs=st.st_mtime_ns,inode=st.st_ino,device=st.st_dev)
        assert all(actual[k]==entry[k] for k in actual)
        identities[ds]=dict(path=entry['path'],identity=actual,unchanged=True,
                           check='size, mtime_ns, inode, device; full-file SHA not computed')
    closure=json.loads((HERE/'closure_checks.json').read_text())
    assert len(closure)==396 and all(c['passed_regular'] for c in closure)
    independent=json.loads((HERE/'independent_validation.json').read_text())
    assert independent['status']=='PASS' and independent['inputsUnchanged']
    results=json.loads((HERE/'analysis_results.json').read_text())
    seen={}
    for line in log.splitlines():
        if not line.startswith('COLUMN_COUNTS\t'):continue
        fields=line.split('\t');cent,selection,dataset=fields[1:4]
        values=dict(field.split('=',1) for field in fields[4:])
        key='__'.join((dataset,cent,selection));result=results['curves'][key]
        for label,expected in [('SE','totalSE'),('ME','totalME'),('SH','normSE'),('MH','normME')]:
            assert float(values[label])==result[expected],(key,label)
        full=sum(p['value'] is not None for p in result['points'])
        high=sum(p['value'] is not None and cfg['display']['high_k_range'][0]<=(p['low']+p['high'])/2<cfg['display']['high_k_range'][1]
                 for p in result['points'])
        assert int(values['validFull'])==full and int(values['validHigh'])==high,key
        seen[key]=seen.get(key,0)+1
    assert len(seen)==242 and set(seen.values())=={2}
    text=out/'centrality_pair_rapidity10_SE_ME_CF_root.txt'
    subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=txtwrite',
        '-sOutputFile='+str(text),str(files[0])],check=True)
    extracted=text.read_text()
    assert 'EP1' in extracted and 'EP8' in extracted and 'cent9' in extracted
    rendered=list((out/'previews').glob('page-[0-9][0-9][0-9].png'))
    assert len(rendered)==121
    artifact_paths=files+[HERE/x for x in ('study.yaml','input_manifest.json','projection_plan.tsv',
        'prepare_study.py','extract_cf_inputs.C','extract.log','raw_projections.tsv',
        'raw_projections.root','source_axes.tsv','analyze_cf.py','analyze.log',
        'analysis_results.json','metrics.csv','cf_covariances.npz','closure_checks.json',
        'write_cf_objects.C','write_cf_objects.log','cf_histograms.root',
        'draw_split_canvases.C','draw_canvas_batch.C','run_data.py','run_canvases.py','render_previews.py','canvases_all.root',
        'canvas_jobs.tsv','all_pages.tsv','draw.log','verify_outputs.py',
        'independent_validate.py','independent_validation.json','source_copy_manifest.json')]
    manifest=dict(finalizedAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        settings=cfg,input_root_identities=identities,
        checks=dict(pdf_page_counts=counts,rendered_pages=len(rendered),
            root_exit_code=0,regular_projection_closures=len(closure),
            maximum_regular_count_difference=max(c['max_regular_count_difference'] for c in closure),
            text_extraction=True,drawn_column_counts_and_valid_points=242,
            main_and_part_counts_match=True,independent_validation=independent),
        production_code_update='Paused by explicit user instruction; no production source edits, Pico reads or jobs in this task.',
        artifacts=[dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path)) for path in artifact_paths])
    (out/'canvas_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('CENTRALITY_PDF_VALIDATION_OK files=12 pages=121 previews=121 input_roots_unchanged=YES')


if __name__=='__main__':main()
