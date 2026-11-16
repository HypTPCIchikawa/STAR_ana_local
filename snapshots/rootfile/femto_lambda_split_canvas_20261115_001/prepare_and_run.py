#!/usr/bin/env python3
"""Make ROOT Divide canvases from existing normalized projected histograms."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'rootfile/femto_lambda_ep_rapidity_cf_20261115_001'
OUT = ROOT / 'share/figure/20261115_lambda_split_canvases'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    config = yaml.safe_load((SOURCE / 'study.yaml').read_text())
    data = json.loads((SOURCE / 'analysis_results.json').read_text())
    input_root = SOURCE / 'cf_histograms.root'
    original_identity = {key: input_root.stat().__getattribute__(key)
                         for key in ('st_size', 'st_mtime_ns', 'st_ino', 'st_dev')}
    cent = next(c for c in config['centralities'] if c['id'] == 'c0_20')
    rapidity = config['pair_rapidity']
    step = (rapidity['maximum'] - rapidity['minimum']) / rapidity['n_bins']
    rows = []
    for i in range(1, rapidity['n_bins'] + 1):
        selection = 'pair_bin%02d' % i
        low = rapidity['minimum'] + (i - 1) * step
        high = rapidity['minimum'] + i * step
        title = '%s; %.1f #leq y_{pair,Lab} < %.1f' % (cent['label'], low, high)
        counts = [data['curves']['__'.join((ep, cent['id'], selection))]['normSE']
                  for ep in config['inputs']]
        if all(x == 0 for x in counts):
            note = 'No SE normalization counts'
        elif all(x == 10 for x in counts):
            note = 'Low statistics: S_{H} = 10'
        else:
            note = '-'
        rows.append('\t'.join((cent['id'], selection, title, note)))
    plan = HERE / 'pair_rapidity_pages.tsv'
    plan.write_text('\n'.join(rows) + '\n')
    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / 'pair_rapidity10_SE_ME_CF_root.pdf'
    canvases = HERE / 'pair_rapidity10_canvases.root'
    args = [str(input_root), str(plan), ','.join(config['inputs']), str(pdf), str(canvases)]
    invocation = str(HERE / 'draw_split_canvases.C') + '+(' + ','.join(json.dumps(a) for a in args) + ')'
    command = ['bash', str(ROOT / 'script/singularity_root_batch.sh'),
               str(ROOT / config['inputs']['EP8']['mainconf']), 'root4star', '-b', '-q', invocation]
    settings = dict(source_study=str(SOURCE / 'study.yaml'), input_root=str(input_root),
                    datasets=list(config['inputs']), centrality=cent, pair_rapidity=rapidity,
                    normalization=config['normalization'], rebin_factor=config['rebin_factor'],
                    display=config['display'], canvas_columns=2, canvas_rows=3,
                    rows=['SE + scaled ME', 'CF full', 'CF high'],
                    statistical_assumptions=config['statistical_model'],
                    skipped_ep_planes_reason='No kstar by event-plane joint histogram is saved in the production ROOT.',
                    command=command)
    (HERE / 'canvas_settings.yaml').write_text(yaml.safe_dump(settings, sort_keys=False))
    with (HERE / 'draw.log').open('w') as log:
        log.write('COMMAND ' + json.dumps(command) + '\n'); log.flush()
        result = subprocess.run(command, cwd=str(ROOT), stdout=log, stderr=subprocess.STDOUT)
        log.write('\nRETURN_CODE %s\n' % result.returncode)
    if result.returncode:
        raise RuntimeError('ROOT drawing failed; see draw.log')
    assert 'DRAW_SPLIT_CANVASES_OK pages=10' in (HERE / 'draw.log').read_text()
    assert all(getattr(input_root.stat(), k) == v for k, v in original_identity.items())
    manifest = dict(createdAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
                    settings=settings, input_root_unchanged=True, artifacts=[])
    for path in [input_root, SOURCE / 'analysis_results.json', SOURCE / 'study.yaml',
                 HERE / 'draw_split_canvases.C', Path(__file__), plan,
                 HERE / 'canvas_settings.yaml', HERE / 'draw.log', pdf, canvases]:
        manifest['artifacts'].append(dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path)))
    (OUT / 'canvas_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print('PAIR_RAPIDITY_CANVASES_OK', pdf)


if __name__ == '__main__':
    main()
