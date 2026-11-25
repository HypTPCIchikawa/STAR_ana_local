#!/usr/bin/env python3
"""Create dedicated EP/merging variants without changing the reference configs."""
import copy
import datetime
import hashlib
import json
from pathlib import Path
import yaml


ROOT = Path(__file__).resolve().parents[2]
STUDY = Path(__file__).resolve().parent
BASE_STUDY = ROOT / 'rootfile/femto_lambda_period_mix100_20261119_001'
BASE_MAIN = ROOT / 'config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_yaml(path):
    with Path(path).open() as source:
        value = yaml.safe_load(source)
    assert isinstance(value, dict), 'Expected a YAML mapping: ' + str(path)
    return value


def fresh_bytes(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as output:
        output.write(data)


def dump_yaml(value):
    return yaml.safe_dump(value, sort_keys=False, allow_unicode=True).encode()


def main():
    assert not (STUDY / 'cases.json').exists(), 'Existing cases require inspection before regeneration'
    main = load_yaml(BASE_MAIN)
    refs = {key: ROOT / 'config' / main[key]
            for key in ['analysis', 'maker', 'femtoHist', 'nuclearHist', 'mixing', 'reducedTree', 'eventPlane']}
    originals = {str(path): sha(path) for path in [BASE_MAIN] + list(refs.values())}
    analysis, maker = load_yaml(refs['analysis']), load_yaml(refs['maker'])
    root_hist, nuclear_hist = load_yaml(refs['femtoHist']), load_yaml(refs['nuclearHist'])
    mixing, reduced, plane = load_yaml(refs['mixing']), load_yaml(refs['reducedTree']), load_yaml(refs['eventPlane'])
    assert reduced['mode'] == 'replay' and not reduced['resetPoolAtSegmentBoundary'] and not reduced['requireOriginalMixing']
    assert plane['bins'] == 8 and plane['harmonic'] == 1
    assert mixing['bufferSize'] == mixing['maxMixEvents'] == 100
    assert mixing['mixingMode'] == 'bufferAll' and mixing['mixBothDirections'] is True
    assert 'mergingPairCut' not in maker, 'Reference maker unexpectedly contains a pair-cut block'
    assert 'hEpdPoolBinVsCent' not in root_hist['histograms']
    reference_map = nuclear_hist['histograms']['hDphiDeta_proton_4He']
    assert reference_map['xAxis']['nBins'] > 0 and reference_map['yAxis']['nBins'] > 0

    period_list = BASE_STUDY / 'periods.list'
    periods = [line for line in period_list.read_text().splitlines() if line.strip()]
    assert len(periods) == len(set(periods)) == 20
    source_manifest_path = BASE_STUDY / 'input_manifest.json'
    source_manifest = json.loads(source_manifest_path.read_text())
    assert periods == [row['path'] for row in source_manifest['sourceFiles']]
    for row in source_manifest['sourceFiles']:
        info = Path(row['path']).stat()
        assert (info.st_size, info.st_mtime_ns, info.st_ino, info.st_dev) == (
            row['bytes'], row['mtimeNs'], row['inode'], row['device'])
    baseline = json.loads((BASE_STUDY / 'completion_verified_20261120.json').read_text())
    baseline_outputs = json.loads((BASE_STUDY / 'output_snapshot_20261120.json').read_text())
    assert baseline['passed'] and baseline['outputRoots'] == 20
    baseline_paths = [row['path'] for row in baseline_outputs['outputs']]
    assert len(baseline_paths) == len(set(baseline_paths)) == 20
    for row in baseline_outputs['outputs']:
        info = Path(row['path']).stat()
        assert (info.st_size, info.st_mtime_ns, info.st_ino, info.st_dev) == (
            row['bytes'], row['mtimeNs'], row['inode'], row['device'])

    pending, cases = [], []
    for bins in [1, 8, 16, 32]:
        for enabled in [False, True]:
            cut_name = 'On' if enabled else 'Off'
            case_id = 'ep' + str(bins) + '_merging' + cut_name
            ana = 'auau3p85_anaFemtoLambda_4He_tree_' + case_id + '_mix100'
            paths = {
                'analysis': Path('analysis') / ('analysis_info_' + ana + '.yaml'),
                'maker': Path('maker') / ('maker_' + ana + '.yaml'),
                'femtoHist': Path('hist') / ('hist_' + ana + '_kf.yaml'),
                'nuclearHist': Path('hist') / ('hist_' + ana + '_nuclear.yaml'),
                'mixing': Path('cuts/mixing') / ('mixing_' + ana + '.yaml')}
            variants = {key: copy.deepcopy(value) for key, value in [
                ('analysis', analysis), ('maker', maker), ('femtoHist', root_hist),
                ('nuclearHist', nuclear_hist), ('mixing', mixing)]}
            info = variants['analysis']['analysis']
            for key in ['anaName', 'jobName', 'name', 'outputFileStem', 'scratchSubdir']:
                info[key] = ana
            info.update(filesPerJob=1, maxEvents=-1,
                        logDir=str(ROOT / 'log' / ana), errDir=str(ROOT / 'err' / ana), workDir=str(ROOT))
            variants['maker']['mergingPairCut'] = dict(enabled=enabled, nucleus='he4', daughter='proton',
                phiMax=0.05, etaMax=0.05, unreachablePolicy='keep', writeQA=True)
            variants['mixing']['nEventPlaneBins'] = bins
            ep_axis = dict(nBins=bins, min=-0.5, max=bins-0.5)
            hist = variants['femtoHist']
            hist['axes']['EpdBin'] = ep_axis
            hist['axes']['EpdPoolBin'] = ep_axis
            plane_keys = []
            for key, definition in hist['histograms'].items():
                if 'VsEpdBinVsCent_lambda_' in key:
                    definition['yAxis'] = ep_axis
                    definition['title'] = ('Accepted ' + key + ', EP' + str(bins) +
                        ' current-event labels;k* [GeV/c];current-event EPD pairing bin;cent9')
                    plane_keys.append(key)
            assert len(plane_keys) == 16
            hist['histograms']['hEpdPoolBinVsCent'] = dict(type='TH2D',
                title='Accepted EPD events classified for EP' + str(bins) +
                      ' mixing;current-event EPD pool bin;cent9',
                xAxis=ep_axis, yAxis=hist['axes']['Cent9'])
            nuclear = variants['nuclearHist']['histograms']
            assert 'hDphiDetaPairCutQA' not in nuclear
            nuclear['hDphiDetaPairCutQA'] = dict(type='TH2D',
                title='He4 merging pair-cut decisions;decision code;pair type (0=SE, 1=ME)',
                xAxis=dict(nBins=12, min=-0.5, max=11.5), yAxis=dict(nBins=2, min=-0.5, max=1.5))
            for daughter in ['proton', 'pion']:
                for mixed in [False, True]:
                    for stage in ['before', 'after']:
                        name = 'hDphiDetaPair_' + daughter + '_' + ('Mixed_' if mixed else '') + '4He_' + stage
                        assert name not in nuclear
                        nuclear[name] = dict(type='TH2D',
                            title=('ME' if mixed else 'SE') + ' signal ' + daughter + '-4He ' + stage +
                                  ' pair cut;#Delta#phi* [rad];#Delta#eta',
                            xAxis=copy.deepcopy(reference_map['xAxis']), yAxis=copy.deepcopy(reference_map['yAxis']))
            variant_main = copy.deepcopy(main)
            for key, relative in paths.items():
                variant_main[key] = relative.as_posix()
                pending.append((ROOT / 'config' / relative, dump_yaml(variants[key])))
            main_relative = Path('config/mainconf') / ('main_' + ana + '.yaml')
            pending.append((ROOT / main_relative, dump_yaml(variant_main)))
            reuse = bins == 8 and not enabled
            cases.append(dict(caseId=case_id, epBins=bins, cutEnabled=enabled, anaName=ana,
                mainconf=main_relative.as_posix(), outputFileStem=ana, filesPerJob=1,
                outputDirectory=str(ROOT / 'rootfile' / ana / 'production20261120'),
                pilotRequired=True, productionRequired=not reuse, baselineReuse=reuse,
                existingBaselineOutputDirectory=str(Path(baseline_paths[0]).parent) if reuse else None,
                existingBaselineOutputs=baseline_paths if reuse else [],
                configFiles={key: ('config/' + value.as_posix()) for key, value in paths.items()}))

    assert len(pending) == len({path for path, _ in pending}) == 48
    for path, _ in pending:
        assert not path.exists(), 'Refusing an existing config: ' + str(path)
    assert not (STUDY / 'periods.list').exists(), 'Refusing an existing study period list'
    assert all(sha(path) == expected for path, expected in originals.items()), 'Reference config changed'
    for path, data in pending:
        fresh_bytes(path, data)
    fresh_bytes(STUDY / 'periods.list', period_list.read_bytes())
    records = [dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest()) for path, data in pending]
    result = dict(schemaVersion=1, createdJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        projectRoot=str(ROOT), studyDirectory=str(STUDY), inputList=str(STUDY / 'periods.list'),
        inputListSha256=sha(STUDY / 'periods.list'), inputManifest=str(source_manifest_path),
        inputManifestSha256=sha(source_manifest_path), expectedPeriods=20,
        storedEvents=source_manifest['storedEvents'], storedCanonicalEpdBins=8,
        calibrationConfig=main['eventPlane'], reducedTreeConfig=main['reducedTree'],
        baselineJobid=baseline['jobid'], baselineCompletion=str(BASE_STUDY / 'completion_verified_20261120.json'),
        baselineReuseCaseId='ep8_mergingOff', pilotCases=[case['caseId'] for case in cases],
        productionCases=[case['caseId'] for case in cases if case['productionRequired']],
        plannedProductionJobs=sum(20 for case in cases if case['productionRequired']), cases=cases,
        referenceConfigSha256=originals, generatedConfigFiles=records,
        deferredMissingProducers=34, deferredUncalibratedRuns=81)
    fresh_bytes(STUDY / 'cases.json', (json.dumps(result, indent=2) + '\n').encode())
    assert all(sha(path) == expected for path, expected in originals.items()), 'Reference config changed after generation'
    for case in cases:
        current = load_yaml(ROOT / case['mainconf'])
        assert current['eventPlane'] == main['eventPlane'] and current['reducedTree'] == main['reducedTree']
        for key in main:
            if key not in paths:
                assert current[key] == main[key]
        observed_maker = load_yaml(ROOT / case['configFiles']['maker'])
        cut = observed_maker.pop('mergingPairCut')
        assert observed_maker == maker and cut['enabled'] == case['cutEnabled']
        observed_mix = load_yaml(ROOT / case['configFiles']['mixing'])
        assert observed_mix.pop('nEventPlaneBins') == case['epBins']
        assert observed_mix == {key: value for key, value in mixing.items() if key != 'nEventPlaneBins'}
        observed_root = load_yaml(ROOT / case['configFiles']['femtoHist'])
        assert observed_root['axes']['EpdMixingBin'] == root_hist['axes']['EpdMixingBin']
        assert observed_root['histograms'].pop('hEpdPoolBinVsCent')['xAxis']['nBins'] == case['epBins']
        observed_root['axes'].pop('EpdPoolBin')
        observed_root['axes']['EpdBin'] = copy.deepcopy(root_hist['axes']['EpdBin'])
        for key in plane_keys:
            assert observed_root['histograms'][key]['yAxis']['nBins'] == case['epBins']
            observed_root['histograms'][key] = copy.deepcopy(root_hist['histograms'][key])
        assert observed_root == root_hist, 'Unintended root-hist changes'
        observed_nuclear = load_yaml(ROOT / case['configFiles']['nuclearHist'])
        expected_new_keys = {'hDphiDetaPairCutQA'}
        expected_new_keys.update('hDphiDetaPair_' + daughter + '_' + mode + '4He_' + stage
            for daughter in ['proton', 'pion'] for mode in ['', 'Mixed_'] for stage in ['before', 'after'])
        assert set(observed_nuclear['histograms']) - set(nuclear_hist['histograms']) == expected_new_keys
        for key in expected_new_keys:
            definition = observed_nuclear['histograms'].pop(key)
            assert definition['type'] == 'TH2D'
            if key != 'hDphiDetaPairCutQA':
                assert definition['xAxis'] == reference_map['xAxis'] and definition['yAxis'] == reference_map['yAxis']
        assert observed_nuclear == nuclear_hist, 'Unintended nuclear-hist changes'
    validation = dict(passed=True, cases=len(cases), sourcePeriods=20,
        allReferenceConfigsUnchanged=True, allNonEPHistogramDefinitionsPreserved=True,
        allOriginalNuclearMapsPreserved=True, canonicalEpdBins=plane['bins'],
        allReducedModesFalseFalse=True, productionCases=len(result['productionCases']),
        casesManifestSha256=sha(STUDY / 'cases.json'))
    fresh_bytes(STUDY / 'config_validation.json', (json.dumps(validation, indent=2) + '\n').encode())
    print(json.dumps(dict(cases=len(cases), generatedConfigs=len(pending),
        productionCases=len(result['productionCases']), plannedProductionJobs=result['plannedProductionJobs'],
        baselineReuseCaseId=result['baselineReuseCaseId'], manifest=str(STUDY / 'cases.json'))))


if __name__ == '__main__':
    main()
