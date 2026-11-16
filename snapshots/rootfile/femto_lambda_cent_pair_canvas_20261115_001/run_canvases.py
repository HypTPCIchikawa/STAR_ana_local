#!/usr/bin/env python3
"""Write the complete and per-centrality PDFs with ROOT Canvas Divide."""
from pathlib import Path
import json
import subprocess
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    cfg=yaml.safe_load((HERE/'study.yaml').read_text())
    out=ROOT/cfg['canvas']['output_directory'];out.mkdir(parents=True,exist_ok=True)
    parts=out/'by_centrality';parts.mkdir(exist_ok=True)
    rapidity=cfg['pair_rapidity'];step=(rapidity['maximum']-rapidity['minimum'])/rapidity['n_bins']
    pages=[];jobs=[];index=[]
    input_root=str(HERE/'cf_histograms.root')
    for cent in cfg['centralities']:
        block=[]
        for sel in cfg['canvas']['selections']:
            if sel=='inclusive':
                label='all pair rapidities'
            else:
                i=int(sel.rsplit('bin',1)[1]);low=rapidity['minimum']+(i-1)*step;high=low+step
                label='%.1f #leq y_{pair,Lab} < %.1f'%(low,high)
            title=cent['label']+' (cent9 '+str(cent['min'])
            if cent['min']!=cent['max']:title+='-'+str(cent['max'])
            title+='); '+label
            block.append('\t'.join((cent['id'],sel,title,'-')))
        plan=HERE/(cent['id']+'_pages.tsv');plan.write_text('\n'.join(block)+'\n')
        pdf=parts/('lambda_he4_'+cent['id']+'_SE_ME_CF_root.pdf')
        canvas_root=HERE/('canvases_'+cent['id']+'.root')
        jobs.append('\t'.join((input_root,str(plan),','.join(cfg['canvas']['columns']),str(pdf),str(canvas_root))))
        index.append(dict(centrality=cent,pdf=str(pdf.relative_to(ROOT)),pages=len(block),
                          combined_first_page=len(pages)+1,combined_last_page=len(pages)+len(block)))
        pages.extend(block)
    plan=HERE/'all_pages.tsv';plan.write_text('\n'.join(pages)+'\n')
    combined=out/'centrality_pair_rapidity10_SE_ME_CF_root.pdf'
    jobs.insert(0,'\t'.join((input_root,str(plan),','.join(cfg['canvas']['columns']),str(combined),str(HERE/'canvases_all.root'))))
    jobfile=HERE/'canvas_jobs.tsv';jobfile.write_text('\n'.join(jobs)+'\n')
    (out/'pdf_index.json').write_text(json.dumps(dict(combined_pdf=str(combined.relative_to(ROOT)),
        combined_pages=len(pages),centralities=index),ensure_ascii=False,indent=2)+'\n')
    invocation=str(HERE/'draw_canvas_batch.C')+'+('+json.dumps(str(jobfile))+')'
    command=['bash',str(ROOT/'script/singularity_root_batch.sh'),str(ROOT/cfg['inputs']['EP8']['mainconf']),
             'root4star','-b','-q',invocation]
    with (HERE/'draw.log').open('w') as log:
        log.write('COMMAND '+json.dumps(command)+'\n');log.flush()
        result=subprocess.run(command,cwd=str(ROOT),stdout=log,stderr=subprocess.STDOUT)
        log.write('\nRETURN_CODE %s\n'%result.returncode)
    if result.returncode:raise RuntimeError('ROOT Canvas output failed; see draw.log')
    assert 'DRAW_CANVAS_BATCH_OK files=12' in (HERE/'draw.log').read_text()
    print('CENTRALITY_PDFS_OK combined_pages=%d pdf_files=%d'%(len(pages),len(jobs)))


if __name__=='__main__':main()
