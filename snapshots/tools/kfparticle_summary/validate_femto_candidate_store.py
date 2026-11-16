#!/usr/bin/env python3
"""Run a bounded real-data candidate-store/replay comparison after a STAR build."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import yaml


REPO = Path(__file__).resolve().parents[2]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def config_reference(mainconf, reference):
    if not isinstance(reference, str) or not reference:
        raise ValueError('A nonempty config reference is required')
    if Path(reference).is_absolute():
        return Path(reference).resolve(strict=True)
    text = str(mainconf)
    index = text.rfind('/config/')
    if index < 0:
        raise ValueError('The mainconf must reside below config/')
    return (Path(text[:index])/'config'/reference).resolve(strict=True)


def read_main(path, mode):
    path = path.resolve(strict=True)
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise ValueError('mainconf must be a mapping: '+str(path))
    feature = None
    if 'candidateStore' in data:
        feature = yaml.safe_load(config_reference(path, data['candidateStore']).read_text())
    if mode is None:
        if feature is not None:
            raise ValueError('The baseline must have candidateStore absent')
    elif not isinstance(feature, dict) or feature.get('mode') != mode:
        raise ValueError('Wrong candidateStore mode in '+str(path))
    return path, data, feature


def source_hashes(configs):
    paths = set()
    for path, data, unused in configs:
        paths.add(path)
        for value in data.values():
            if isinstance(value, str) and value.endswith(('.yaml', '.yml')):
                paths.add(config_reference(path, value))
    for name in ['libStarAnaConfig.so', 'libStRefMultCorr.so', 'libKFParticle.so',
                 'libStKfParticleCommon.so', 'libStCommon.so', 'libStFemtoMaker.so']:
        paths.add((REPO/'lib'/name).resolve(strict=True))
    for name in ['analysis/anaFemtoLambda_4He.C', 'analysis/run_anaFemtoLambda_4He.C',
                 'analysis/FemtoLambdaRunSupport.h', 'tests/compareFemtoCandidateReplay.C']:
        paths.add(REPO/name)
    return {str(path): {'sha256': sha(path), 'bytes': path.stat().st_size} for path in sorted(paths)}


def cpp_string(value):
    return json.dumps(str(value), ensure_ascii=True)


def root_fingerprint(path):
    path = path.resolve(strict=True)
    with path.open('rb') as source:
        if source.read(4) != b'root':
            raise ValueError('Invalid cached ROOT header: '+str(path))
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtimeNs=stat.st_mtime_ns, sha256=sha(path))


def reusable_reconstruction(record_path, provenance, requested, current_hashes, configs):
    path = record_path.resolve(strict=True)
    old = json.loads(path.read_text())
    if old.get('input') != provenance or old.get('requestedPicoEvents') != requested:
        raise ValueError('Reused reconstruction has different frozen input or requested events')
    previous = old.get('configAndLibraryHashesBefore')
    if not isinstance(previous, dict) or set(previous) != set(current_hashes):
        raise ValueError('Reused reconstruction source/config/library inventory differs')
    checker = str(REPO/'tests/compareFemtoCandidateReplay.C')
    if checker not in previous:
        raise ValueError('Producer record does not fingerprint the comparison macro')
    for name, value in current_hashes.items():
        if name != checker and value != previous[name]:
            raise ValueError('Reused reconstruction source/config/library changed: '+name)
    commands = old.get('commands', [])
    roots = {}
    counters = {}
    for index, stage in enumerate(['baseline', 'write']):
        source = path.parent/(stage+'.root')
        matches = [c for c in commands if c.get('stage') == stage]
        if len(matches) != 1 or matches[0].get('returnCode') != 0:
            raise ValueError('Producer record lacks one successful '+stage+' command')
        argv = matches[0].get('argv', [])
        if (len(argv) != 7 or argv[0] != 'bash'
                or Path(argv[1]).resolve() != REPO/'script/singularity_run_anaFemtoLambda_4He.sh'
                or Path(argv[2]).resolve() != configs[index][0]
                or Path(argv[3]).resolve() != Path(provenance['path']).resolve()
                or Path(argv[4]).resolve() != source.resolve()
                or argv[6] != str(requested)):
            raise ValueError('Producer '+stage+' command does not match the explicit validation inputs')
        count = old.get(stage+'Counters')
        if (not isinstance(count, dict) or count.get('requested', 0) <= 0
                or count.get('requested') != count.get('completed')
                or count.get('reconstructed', 0) <= 0):
            raise ValueError('Producer '+stage+' completion counters are not valid')
        if requested != -1 and count['requested'] != requested:
            raise ValueError('Producer '+stage+' counters differ from requested limit')
        roots[stage] = root_fingerprint(source)
        counters[stage] = count
    if counters['baseline'] != counters['write']:
        raise ValueError('Producer baseline/write completion counters differ')
    return dict(producerRecord=str(path), producerRecordSha256=sha(path),
                producerPhase=old.get('phase'),
                producerCommands=[c for c in commands if c.get('stage') in ['baseline', 'write']],
                producerCounters=counters, originalRootsBefore=roots,
                strictMatchingExceptComparisonMacro=True,
                comparisonMacro=dict(path=checker, before=previous[checker], after=current_hashes[checker],
                                     changed=previous[checker] != current_hashes[checker]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['baseline-mainconf', 'write-mainconf', 'replay-mainconf', 'input',
                 'input-provenance', 'output-dir']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--variant-mainconf', type=Path)
    parser.add_argument('--reuse-reconstruction-record', type=Path,
                        help='Reuse successful baseline/write ROOTs only after strict source/input checks.')
    parser.add_argument('--events', type=int, required=True)
    parser.add_argument('--job-prefix', required=True)
    parser.add_argument('--pre-pid-qa', choices=['include', 'exclude'], required=True)
    args = parser.parse_args()
    if args.events == 0 or args.events < -1 or not args.job_prefix:
        raise ValueError('Explicit positive events or -1 and a job prefix are required')
    configurations = [read_main(args.baseline_mainconf, None),
                      read_main(args.write_mainconf, 'write'),
                      read_main(args.replay_mainconf, 'replay')]
    baseline = configurations[0][1]
    physical = {key: value for key, value in baseline.items() if key != 'candidateStore'}
    for path, cfg, unused in configurations[1:]:
        if {k: v for k, v in cfg.items() if k != 'candidateStore'} != physical:
            raise ValueError('Writer/replay must use the same production references: '+str(path))
    if args.variant_mainconf:
        variant = read_main(args.variant_mainconf, 'replay')
        if {k: v for k, v in variant[1].items() if k not in ['candidateStore', 'mixing']} != {
                k: v for k, v in physical.items() if k != 'mixing'}:
            raise ValueError('The remix variant may change only mixing and candidateStore settings')
        configurations.append(variant)
    source = args.input.resolve(strict=True)
    provenance = json.loads(args.input_provenance.resolve(strict=True).read_text())
    if (provenance.get('year') != 2021 or provenance.get('production') != 'P24iy'
            or Path(provenance['path']).resolve(strict=True) != source
            or source.stat().st_size != provenance['size'] or sha(source) != provenance['sha256']):
        raise ValueError('The explicit pilot input must match its frozen2021 P24iy provenance')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    before = source_hashes(configurations)
    reuse = (reusable_reconstruction(args.reuse_reconstruction_record, provenance, args.events, before,
                                     configurations) if args.reuse_reconstruction_record else None)
    status = dict(phase='running', startedAtJst=datetime.datetime.now(
        datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        input=provenance, requestedPicoEvents=args.events, commands=[],
        configAndLibraryHashesBefore=before, comparisons={},
        runnerSource=dict(path=str(Path(__file__).resolve()), sha256=sha(Path(__file__))))
    if reuse:
        status['reusedReconstruction'] = reuse
        status['reconstructionReused'] = True
        status['baselineCounters'] = reuse['producerCounters']['baseline']
        status['writeCounters'] = reuse['producerCounters']['write']
    save(out/'validation.json', status)

    def execute(label, command):
        logfile = out/(label+'.log')
        start = time.monotonic()
        with logfile.open('x') as log:
            result = subprocess.run(command, cwd=str(REPO), stdout=log,
                                    stderr=subprocess.STDOUT, universal_newlines=True)
        status['commands'].append(dict(stage=label, argv=command, returnCode=result.returncode,
                                       wallSeconds=time.monotonic()-start, log=str(logfile)))
        save(out/'validation.json', status)
        if result.returncode:
            raise RuntimeError(label+' failed: '+str(logfile))
        return logfile.read_text()

    def reconstruct(label, mainconf, input_path, events):
        output = out/(label+'.root')
        command = ['bash', str(REPO/'script/singularity_run_anaFemtoLambda_4He.sh'),
                   str(mainconf), str(input_path), str(output), args.job_prefix+'_'+label, str(events)]
        text = execute(label, command)
        if not output.is_file():
            raise RuntimeError(label+' did not create its ROOT output')
        if label in ['replay', 'variant']:
            if '[candidate-store-input] picoDstReader=false lambdaProvider=false' not in text:
                raise RuntimeError('Replay did not confirm that PicoDst and Lambda provider are absent')
            match = re.search(r'FemtoLambda replay job=.*? requested=(\d+) completed=(\d+) available=(\d+) status=(\d+)', text)
            if not match or match.group(1) != match.group(2) or int(match.group(4)):
                raise RuntimeError('Invalid replay completion counters')
            status[label+'Counters'] = dict(requested=int(match.group(1)), completed=int(match.group(2)),
                                            availableStoreRows=int(match.group(3)), picoDstReader=False,
                                            lambdaProvider=False)
        else:
            match = re.search(r'FemtoLambda job=.*? requested=(\d+) completed=(\d+) reconstructed=(\d+) status=(\d+)', text)
            if not match or match.group(1) != match.group(2) or int(match.group(4)):
                raise RuntimeError('Invalid PicoDst completion counters')
            status[label+'Counters'] = dict(requested=int(match.group(1)), completed=int(match.group(2)),
                                            reconstructed=int(match.group(3)))
        save(out/'validation.json', status)
        return output

    def compare(label, reference, other, mode):
        report = out/(label+'.json')
        expression = 'tests/compareFemtoCandidateReplay.C+('+','.join([
            cpp_string(reference), cpp_string(other), cpp_string(report),
            'true' if args.pre_pid_qa == 'include' else 'false', cpp_string(mode)])+')'
        command = ['bash', str(REPO/'script/singularity_root_batch.sh'),
                   str(configurations[0][0]), 'root4star', '-l', '-b', '-q', expression]
        execute(label, command)
        record = json.loads(report.read_text())
        if record.get('phase') != 'passed':
            raise RuntimeError('Comparison did not pass: '+str(report))
        status['comparisons'][label] = record
        save(out/'validation.json', status)
        return record

    try:
        if reuse:
            direct = Path(reuse['originalRootsBefore']['baseline']['path'])
            writer = Path(reuse['originalRootsBefore']['write']['path'])
        else:
            direct = reconstruct('baseline', configurations[0][0], source, args.events)
            writer = reconstruct('write', configurations[1][0], source, args.events)
        writer_comparison = compare('baseline_write', direct, writer, 'closure')
        if writer_comparison['replayCandidateStoreRows'] <= 0:
            raise RuntimeError('Writer did not persist candidate event rows')
        replay = reconstruct('replay', configurations[2][0], writer, -1)
        compare('baseline_replay', direct, replay, 'closure')
        if status['replayCounters']['availableStoreRows'] != writer_comparison['replayCandidateStoreRows']:
            raise RuntimeError('Writer tree row count differs from replay declared entries')
        if args.variant_mainconf:
            variant = reconstruct('variant', configurations[3][0], writer, -1)
            variation = compare('baseline_variant_same_event', direct, variant, 'sameEventOnly')
            if not variation['probeMixedEventChanged']:
                raise RuntimeError('Variant did not change the observed He4 ME signal histogram')
        after = source_hashes(configurations)
        if before != after or sha(source) != provenance['sha256']:
            raise RuntimeError('Input, configuration, macro or library changed during validation')
        if reuse:
            fingerprints = {name: root_fingerprint(Path(value['path']))
                            for name, value in reuse['originalRootsBefore'].items()}
            if fingerprints != reuse['originalRootsBefore']:
                raise RuntimeError('Original reconstructed ROOTs changed during reuse validation')
            status['reusedReconstruction']['originalRootsAfter'] = fingerprints
            status['reusedReconstruction']['originalRootsUnchanged'] = True
            status['reusedReconstruction']['producerRecordUnchanged'] = (
                sha(Path(reuse['producerRecord'])) == reuse['producerRecordSha256'])
            if not status['reusedReconstruction']['producerRecordUnchanged']:
                raise RuntimeError('Original producer validation record changed during reuse')
        status.update(phase='passed', configAndLibraryHashesAfter=after, sourceInputUnchanged=True,
                      exactlyReproducedCorePairs=True, candidateRows=writer_comparison['replayCandidateStoreRows'],
                      variantSameEventUnchanged=bool(args.variant_mainconf),
                      variantMixedEventChanged=bool(args.variant_mainconf))
    except Exception as exc:
        status.update(phase='failed', error=str(exc))
        if reuse:
            try:
                fingerprints = {name: root_fingerprint(Path(value['path']))
                                for name, value in reuse['originalRootsBefore'].items()}
                status['reusedReconstruction']['originalRootsAfter'] = fingerprints
                status['reusedReconstruction']['originalRootsUnchanged'] = (
                    fingerprints == reuse['originalRootsBefore'])
                status['reusedReconstruction']['producerRecordUnchanged'] = (
                    sha(Path(reuse['producerRecord'])) == reuse['producerRecordSha256'])
            except Exception as fingerprint_error:
                status['reusedReconstruction']['fingerprintError'] = str(fingerprint_error)
        save(out/'validation.json', status)
        raise
    status['finishedAtJst'] = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    save(out/'validation.json', status)
    print(json.dumps(dict(phase=status['phase'], output=str(out),
                         candidateRows=status['candidateRows'],
                         variantMixedEventChanged=status['variantMixedEventChanged'])))


if __name__ == '__main__':
    main()
