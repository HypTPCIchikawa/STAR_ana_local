#!/usr/bin/env python3
"""Convert numerical results into explicit, ROOT5-readable plot tables."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def token(value):
    return 'NA' if value is None else repr(value)


def main():
    data=json.loads((HERE/'analysis_results.json').read_text());cfg=data['config'];curves=data['curves']
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    conditions={c['id']:c for c in manifest['conditions']}
    cent_labels={x['id']:x['label'] for x in cfg['centralities']}
    y=cfg['pair_rapidity'];step=(y['maximum']-y['minimum'])/y['n_bins']
    def bounds(condition):
        return y['minimum']+(condition['first']-1)*step,y['minimum']+condition['last']*step
    def label(selection):
        if selection=='inclusive':return 'inclusive'
        if selection=='he4_below':return 'y_{^{4}He,lab} < %.1f'%cfg['he4_rapidity_threshold']
        if selection=='he4_above':return 'y_{^{4}He,lab} > %.1f'%cfg['he4_rapidity_threshold']
        if selection.startswith(('pair_bin','pair_merge')):
            low,high=bounds(conditions[selection]);return '%.1f #leq y_{pair,lab} < %.1f'%(low,high)
        return selection
    with (HERE/'plot_metadata.tsv').open('w') as f:
        for r in curves.values():
            fields=[r['dataset'],r['centrality'],r['selection'],token(r['MEscaleToSE'])]
            fields+=list(map(token,cfg['normalization']+cfg['display']['k_range']+cfg['display']['high_k_range']))
            f.write('\t'.join(fields)+'\n')
    selected=['inclusive','he4_below','he4_above']+['pair_bin%02d'%i for i in range(1,11)]+['pair_merge%02d'%i for i in range(1,6)]
    with (HERE/'plot_panels.tsv').open('w') as f:
        for c in cfg['centralities']:
            for selection in selected:
                f.write('\t'.join([c['id'],selection,'#Lambda-^{4}He signal, '+cent_labels[c['id']]+', '+label(selection)])+'\n')
    with (HERE/'plot_summary.tsv').open('w') as f:
        for c in cfg['centralities']:
            for family,names in [('pair10',['pair_bin%02d'%i for i in range(1,11)]),
                                 ('pair5',['pair_merge%02d'%i for i in range(1,6)]),
                                 ('he4',['he4_below','he4_above'])]:
                for index,selection in enumerate(names):
                    if family=='he4':x,ex=float(index),0.
                    else:
                        low,high=bounds(conditions[selection]);x,ex=(low+high)/2,(high-low)/2
                    for ds in cfg['inputs']:
                        r=curves['__'.join([ds,c['id'],selection])];t=r['windows']['k080_100'];s=r['shape_ratio']
                        fields=[ds,c['id'],family,token(x),token(ex),label(selection),token(t['value']),token(t['error']),token(s['value']),token(s['error'])]
                        f.write('\t'.join(fields)+'\n')
    print('Prepared',len(selected)*len(cfg['centralities']),'plot panels and summary series')


if __name__=='__main__':main()
