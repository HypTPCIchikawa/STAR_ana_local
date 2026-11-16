#!/usr/bin/env python3
"""Prepare an opt-in EP-class histogram configuration from the covered-run input."""
from pathlib import Path
import copy
import hashlib
import json
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10'
NAME = 'auau3p85_anaFemtoLambda_4He_KFParticle_epd8_kstarSplit_pairY10'
OLD_MAIN = ROOT / ('config/mainconf/main_' + BASE + '_covered25.yaml')
NEW_MAIN = ROOT / ('config/mainconf/main_' + NAME + '_covered25.yaml')
NEW_HIST = ROOT / 'config/hist/hist_auau3p85_femtoLambda_he4_he3_control_epd8KstarPairY10_20261115_kf.yaml'
NEW_INFO = ROOT / ('config/analysis/analysis_info_' + NAME + '_covered25.yaml')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, content):
    with path.open('x') as stream:
        stream.write(content)


def main():
    main = yaml.safe_load(OLD_MAIN.read_text())
    old_hist = ROOT / 'config' / main['femtoHist']
    old_info = ROOT / 'config' / main['analysis']
    mix = yaml.safe_load((ROOT / 'config' / main['mixing']).read_text())
    hist = yaml.safe_load(old_hist.read_text())
    previous = copy.deepcopy(hist)
    class_axis = copy.deepcopy(hist['histograms']['hEpdMixingBinVsCent']['xAxis'])
    assert class_axis == dict(nBins=mix['nEventPlaneBins'], min=-0.5,
                              max=mix['nEventPlaneBins']-0.5)
    added = []
    for side in ('SE', 'ME'):
        for suffix in ('', '_signal', '_leftSB', '_rightSB'):
            name = 'hKstar' + side + 'VsEpBinVsCent_lambda_he4' + suffix
            reference = hist['histograms']['hKstar' + side + 'VsPairYVsCent_lambda_he4' + suffix]
            hist['histograms'][name] = dict(type='TH3D',
                xAxis=copy.deepcopy(reference['xAxis']), yAxis=copy.deepcopy(class_axis),
                zAxis=copy.deepcopy(reference['zAxis']),
                title=side + ' Lambda-He4' + suffix + ', mixing EP class;k* [GeV/c];mixing EP class;cent9')
            added.append(name)
    for key, value in previous['histograms'].items():
        assert hist['histograms'][key] == value
    info = yaml.safe_load(old_info.read_text())
    info['analysis'] = {key: value.replace(BASE, NAME) if isinstance(value, str) else value
                        for key, value in info['analysis'].items()}
    assert info['dataset'] == yaml.safe_load(old_info.read_text())['dataset']
    main['analysis'] = str(NEW_INFO.relative_to(ROOT / 'config'))
    main['femtoHist'] = str(NEW_HIST.relative_to(ROOT / 'config'))
    assert all(not p.exists() for p in (NEW_MAIN, NEW_HIST, NEW_INFO))
    sources = {str(path.relative_to(ROOT)): sha(path) for path in
               (OLD_MAIN, old_hist, old_info, ROOT / 'config' / main['mixing'])}
    write_new(NEW_HIST, yaml.safe_dump(hist, sort_keys=False))
    write_new(NEW_INFO, yaml.safe_dump(info, sort_keys=False))
    write_new(NEW_MAIN, yaml.safe_dump(main, sort_keys=False))
    manifest = dict(originals_sha256=sources, added_histograms=added,
        mixing_ep_bins=mix['nEventPlaneBins'], mainconf=str(NEW_MAIN.relative_to(ROOT)),
        original_histogram_definitions_unchanged=True, dataset_unchanged=True,
        generated={str(p.relative_to(ROOT)): sha(p) for p in (NEW_MAIN, NEW_HIST, NEW_INFO)})
    (HERE / 'config_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('EP_KSTAR_CONFIGS_OK', NEW_MAIN)


if __name__ == '__main__':
    main()
