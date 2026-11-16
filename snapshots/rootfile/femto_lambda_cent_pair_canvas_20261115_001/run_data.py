#!/usr/bin/env python3
"""Execute projections, normalized CF calculations and ROOT object export."""
from pathlib import Path
import json
import subprocess
import sys
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def logged(command, name):
    with (HERE / name).open('w') as log:
        log.write('COMMAND ' + json.dumps(command) + '\n'); log.flush()
        result = subprocess.run(command, cwd=str(ROOT), stdout=log, stderr=subprocess.STDOUT)
        log.write('\nRETURN_CODE %s\n' % result.returncode)
    if result.returncode:
        raise RuntimeError('Command failed: ' + name)


def root_command(macro, args, mainconf):
    invocation = str(HERE / macro) + '+(' + ','.join(json.dumps(str(x)) for x in args) + ')'
    return ['bash', str(ROOT / 'script/singularity_root_batch.sh'),
            str(ROOT / mainconf), 'root4star', '-b', '-q', invocation]


def main():
    cfg = yaml.safe_load((HERE / 'study.yaml').read_text())
    conf = cfg['inputs']['EP8']['mainconf']
    logged(root_command('extract_cf_inputs.C', [HERE / 'projection_plan.tsv',
        HERE / 'raw_projections.tsv', HERE / 'raw_projections.root', HERE / 'source_axes.tsv'], conf), 'extract.log')
    assert 'EXTRACT_CF_INPUTS_OK requests=484' in (HERE / 'extract.log').read_text()
    logged([sys.executable, '-I', str(HERE / 'analyze_cf.py')], 'analyze.log')
    logged([sys.executable, '-I', str(HERE / 'prepare_plots.py')], 'prepare_plots.log')
    logged(root_command('write_cf_objects.C', [HERE / 'curve_bins.tsv',
        HERE / 'plot_metadata.tsv', HERE / 'cf_histograms.root'], conf), 'write_cf_objects.log')
    assert 'WRITE_CF_OBJECTS_OK selections=484' in (HERE / 'write_cf_objects.log').read_text()
    print('CENTRALITY_DATA_OK series=484')


if __name__ == '__main__':
    main()
