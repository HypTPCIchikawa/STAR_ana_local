#!/usr/bin/env python3
"""Freeze centrality-resolved projections from the verified production ROOTs."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import shutil
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'rootfile/femto_lambda_ep_rapidity_cf_20261115_001'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cfg = copy.deepcopy(yaml.safe_load((SOURCE / 'study.yaml').read_text()))
    cfg['centralities'] = [
        dict(id='c0_20', label='0-20%', min=6, max=8),
        dict(id='c0_10', label='0-10%', min=7, max=8),
        dict(id='c0_5', label='0-5%', min=8, max=8),
        dict(id='c5_10', label='5-10%', min=7, max=7),
        dict(id='c10_20', label='10-20%', min=6, max=6),
        dict(id='c20_30', label='20-30%', min=5, max=5),
        dict(id='c30_40', label='30-40%', min=4, max=4),
        dict(id='c40_50', label='40-50%', min=3, max=3),
        dict(id='c50_60', label='50-60%', min=2, max=2),
        dict(id='c60_70', label='60-70%', min=1, max=1),
        dict(id='c70_80', label='70-80%', min=0, max=0)]
    cfg['canvas'] = dict(columns=['EP1', 'EP8'],
        rows=['SE + scaled ME', 'CF full', 'CF high'],
        selections=['inclusive'] + ['pair_bin%02d' % i for i in range(1, cfg['pair_rapidity']['n_bins'] + 1)],
        ordering='centrality then inclusive/pair rapidity',
        normalization='Separate SE_H / ME_H for each EP, centrality and pair rapidity',
        counts='Regular kstar range only; pair counts, not event counts',
        output_directory='share/figure/20261115_lambda_centrality_pair_canvases')
    (HERE / 'study.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    copies = {}
    for name in ('extract_cf_inputs.C', 'analyze_cf.py', 'prepare_plots.py'):
        destination = HERE / name
        assert not destination.exists()
        shutil.copy2(SOURCE / name, destination)
        copies[name] = sha(destination)
        assert copies[name] == sha(SOURCE / name)
    spec = importlib.util.spec_from_file_location('original_prepare_study', SOURCE / 'prepare_study.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.HERE = HERE; module.ROOT = ROOT; module.main()
    original = json.loads((SOURCE / 'input_manifest.json').read_text())['inputs']
    current = json.loads((HERE / 'input_manifest.json').read_text())['inputs']
    assert original == current
    (HERE / 'source_copy_manifest.json').write_text(json.dumps(dict(
        original_study=str(SOURCE), source_copies_sha256=copies,
        prepare_source_sha256=sha(SOURCE / 'prepare_study.py'),
        input_identities_identical=True, production_code_update='paused by user'), indent=2) + '\n')
    print('CENTRALITY_STUDY_PREPARED centralities=11 displayed_pages=121')


if __name__ == '__main__':
    main()
