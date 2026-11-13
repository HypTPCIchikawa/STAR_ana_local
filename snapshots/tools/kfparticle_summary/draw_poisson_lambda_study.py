#!/usr/bin/env python3
"""Draw saved Poisson fits, fixed examples and matched-objective CF comparisons."""
import argparse
import json
from pathlib import Path
import subprocess

import yaml


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--plot-settings',type=Path,required=True)
    parser.add_argument('--chi-square-comparison',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    context=json.loads((out/'poisson_context.json').read_text())
    cfg=yaml.safe_load(args.plot_settings.read_text())
    data=json.loads((out/'group_extract/purity_fit_review.json').read_text())
    current=json.loads((out/'methods_comparison_selected.json').read_text())
    old=json.loads(args.chi_square_comparison.read_text())
    assert [g['id'] for g in current['groups']]==[g['id'] for g in old['groups']]
    for a,b in zip(current['groups'],old['groups']):
        for method in ('method3','method5','method6'):
            assert a['methods'][method]['normalization']['bins']==b['methods'][method]['normalization']['bins']
    rows=['objective\tmethod\tgroup\tlabel\tbin\tlow\thigh\tvalue\terror']
    for objective,doc in [('poisson',current),('chiSquare',old)]:
        for g in doc['groups']:
            for method,result in g['methods'].items():
                for point in result['points']:
                    value=point['cf'];row=[objective,method,g['id'],g['label'],point['bin'],point['low'],point['high'],
                                          value['value'] if value else 'nan',value['error'] if value else 'nan']
                    rows.append('\t'.join(map(str,row)))
    (out/'cf_objective_comparison.tsv').write_text('\n'.join(rows)+'\n')
    def literal(x):return json.dumps(str(x)) if isinstance(x,(str,Path)) else format(x,'.17g')
    def run(name,macro,values):
        call='tools/kfparticle_summary/'+macro+'+('+','.join(literal(v) for v in values)+')'
        cmd=['bash','script/singularity_root_batch.sh',context['mainconf'],'root4star','-l','-b','-q',call]
        (out/(name+'_command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
        with (out/(name+'.log')).open('w') as f:result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
        text=(out/(name+'.log')).read_text()
        if result.returncode or 'Compilation failed' in text or '_FAILED:' in text:
            print(text[-7000:],flush=True);raise RuntimeError('ROOT drawing failed: '+name)
        print('Completed '+name,flush=True)
    s=data['settings'];bounds=[*s['fitRange'],s['mean']-s['nSigma']*s['sigma'],s['mean']+s['nSigma']*s['sigma']]
    for name,table in [('cases',out),('atlas',out/'atlas_tables')]:
        raw=out/(name+'_raw');raw.mkdir(exist_ok=True)
        with (table/'cases_fit_meta.tsv').open() as handle:expected=(sum(1 for line in handle)-1)//2+2
        log=out/('draw_'+name+'.log')
        if args.resume and len(list(raw.glob('*.pdf')))==expected and log.exists() and log.read_text().count('POISSON_ATLAS_COMPLETE')==2:
            print('Retained complete '+name+' drawings',flush=True);continue
        if args.resume and log.exists():
            archive=out/'revisions';archive.mkdir(exist_ok=True)
            saved=archive/(log.name+'.previous')
            if not saved.exists():saved.write_bytes(log.read_bytes())
        run('draw_'+name,'DrawPoissonPurityReview.C',[table/'cases_fit_meta.tsv',table/'cases_fit_bins.tsv',table/'basis_fit_parameters.tsv',raw,*bounds])
    raw=out/'comparison_raw';raw.mkdir(exist_ok=True)
    run('draw_fine_purity','DrawPoissonFinePurity.C',[out/'fine_purity_comparison.tsv',raw,*cfg['fine_x_range'],*cfg['fine_y_range']])
    run('draw_cf','DrawMethod6Comparison.C',[out/'cf_points.tsv',out/'purity_points.tsv',raw,*cfg['cf_range'],*cfg['normalization'],cfg['group_purity_max']])
    run('draw_cf_objectives','DrawPoissonCFComparison.C',[out/'cf_objective_comparison.tsv',raw,*cfg['cf_range'],*cfg['normalization']])
    for name in ('cases','atlas','comparison'):
        target=out/name;target.mkdir(exist_ok=True)
        result=subprocess.run(['python3','tools/kfparticle_summary/crop_method6_figures.py','--input-dir',str(out/(name+'_raw')),
            '--output-dir',str(target),'--margin-points',str(cfg['margin_points'])],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        print(name,result.stdout,flush=True)
        if result.returncode:raise RuntimeError('Cropping failed: '+name)


if __name__=='__main__':main()
