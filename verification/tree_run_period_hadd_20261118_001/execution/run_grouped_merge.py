#!/usr/bin/env python3
"""Validate and concatenate complete reduced-tree producer files in fixed groups."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import socket
import stat
import subprocess
import sys


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fresh_json(path, value):
    with Path(path).open('x') as output:
        json.dump(value, output, indent=2, ensure_ascii=False)
        output.write('\n')


def state_json(path, value):
    temporary = Path(str(path) + '.next')
    with temporary.open('w') as output:
        json.dump(value, output, indent=2, ensure_ascii=False)
        output.write('\n')
        output.flush()
        os.fsync(output.fileno())
    os.replace(str(temporary), str(path))


def fingerprint(row):
    path = Path(row['path'])
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode), 'Source is not a regular file: ' + str(path))
    actual = dict(bytes=info.st_size, inode=info.st_ino, device=info.st_dev,
                  mtimeNs=info.st_mtime_ns)
    require(all(actual[k] == row[k] for k in actual), 'Source changed: ' + str(path))


def run(command, cwd, log):
    print(now(), 'START', log.name, flush=True)
    with log.open('x') as output:
        output.write(json.dumps(command) + '\n')
        output.flush()
        result = subprocess.run(command, cwd=str(cwd), stdout=output, stderr=subprocess.STDOUT)
    require(result.returncode == 0, 'Command failed (%d): %s' % (result.returncode, log))
    print(now(), 'PASS', log.name, flush=True)


def quota(config, evidence, remaining):
    result = subprocess.run(config['quotaCommand'], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True, timeout=45)
    with evidence.with_suffix('.log').open('x') as output:
        output.write(result.stdout + '\nSTDERR\n' + result.stderr)
    require(result.returncode == 0, 'Quota query failed')
    head = next(x for x in result.stdout.splitlines() if x.startswith('mmlsquota:user:HEADER:'))
    row = next(x for x in result.stdout.splitlines() if x.startswith('mmlsquota:user:0:'))
    values = dict(zip(head.split(':')[3:], row.split(':')[3:]))
    require(values['name'] == config['quotaUser'], 'Quota user mismatch')
    free = (int(values['blockQuota']) - int(values['blockUsage']) - int(values['blockInDoubt'])) * 1024
    reserve = int(remaining * config['outputExpansionAllowance']) + config['safetyReserveBytes']
    proof = dict(checkedJst=now(), fields=values, softHeadroomBytesAfterInDoubt=free,
                 remainingSourceBytes=remaining, requiredReserveBytes=reserve, passed=free >= reserve)
    fresh_json(evidence, proof)
    require(proof['passed'], 'Quota reserve insufficient; no next merge started')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root, study = Path(config['projectRoot']), Path(config['studyDirectory'])
    lock = (study / 'merge.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    status_path = study / 'merge_status.json'
    require(not status_path.exists(), 'Prior execution exists; inspect before restarting')
    status = dict(phase='preflight', startedJst=now(), pid=os.getpid(), host=socket.getfqdn(),
                  completedGroups=[], totalGroups=config['expectedGroups'], originalSourcesRetained=True,
                  stage2Started=False, configSha256=sha(args.config), currentGroup=None)
    state_json(status_path, status)
    try:
        require(socket.getfqdn() == config['executionHost'], 'Unexpected execution host')
        for name, expected in config['fingerprints'].items():
            require(sha(name) == expected, 'Pinned artifact changed: ' + name)
        snapshot = json.loads(Path(config['sourceSnapshot']).read_text())
        plan = json.loads(Path(config['assignment']).read_text())
        sources = {x['path']: x for x in snapshot['paths']}
        require(len(sources) == snapshot['frozenSourceCount'] == config['expectedSources'], 'Snapshot count differs')
        groups = plan['groups']
        require(len(groups) == config['expectedGroups'] and plan['wholeProducersPreserved'], 'Unexpected grouping')
        flat = [p for group in groups for p in group['inputs']]
        require(len(flat) == len(set(flat)) == len(sources) and set(flat) == set(sources), 'Grouping lost/repeated sources')
        for row in snapshot['paths']:
            fingerprint(row)
        output_dir = Path(config['outputDirectory'])
        require(not output_dir.exists(), 'Refusing existing output directory')
        output_dir.mkdir()
        operations = study / 'operations'
        operations.mkdir()
        needed_open = max(len(g['inputs']) for g in groups) + 8
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        require(hard == resource.RLIM_INFINITY or hard >= needed_open, 'Open-file hard limit too small')
        if soft < needed_open:
            resource.setrlimit(resource.RLIMIT_NOFILE, (needed_open, hard))
        source_records, seen_segments, schema, provenance = {}, set(), None, None
        remaining = sum(r['bytes'] for r in sources.values())
        for group in groups:
            name = group['group']
            require(name.replace('_', '').isalnum(), 'Unsafe group label')
            status.update(phase='source_validation', currentGroup=name, updatedJst=now())
            state_json(status_path, status)
            directory = operations / name
            directory.mkdir()
            rows = [sources[p] for p in group['inputs']]
            for row in rows:
                fingerprint(row)
            input_list = directory / 'inputs.list'
            with input_list.open('x') as output:
                output.write(''.join(r['path'] + '\n' for r in rows))
            subset = dict(snapshot)
            subset.update(paths=rows, frozenSourceCount=len(rows),
                          inputEventsFromCompletionLogs=sum(r['inputEventsFromCompletionLog'] for r in rows),
                          acceptedEventsFromCompletionLogs=sum(r['acceptedEventsFromCompletionLog'] for r in rows))
            fresh_json(directory / 'source_snapshot.json', subset)
            fresh_json(directory / 'assignment.json', group)
            final = output_dir / (name + '_runs_%s_%s.root' % (group['runMin'], group['runMax']))
            partial = output_dir / (final.stem + '.partial.root')
            require(not os.path.lexists(str(final)) and not os.path.lexists(str(partial)), 'Output already exists')

            def root_action(action):
                expression = '%s("%s","%s","%s","%s",%dLL,%d,%d,%d,%d)' % (
                    config['rootMacro'], action, input_list, directory / (action + '.json'), partial,
                    config['maxTreeBytes'], config['expectedQa'], config['sampleEvents'],
                    needed_open, config['floatUlps'])
                command = ['bash', config['batchWrapper'], config['mainconf'],
                           'root4star', '-l', '-b', '-q', expression]
                run(command, root, directory / (action + '.log'))
                report = json.loads((directory / (action + '.json')).read_text())
                require(report.get('passed') is True, 'Incomplete ROOT report: ' + action)
                return report

            catalog = root_action('catalog')
            run([sys.executable, config['catalogChecker'], '--snapshot', str(directory / 'source_snapshot.json'),
                 '--catalog', str(directory / 'catalog.json'), '--output', str(directory / 'catalog_validation.json'),
                 '--expected-qa', str(config['expectedQa'])], root, directory / 'catalog_validation.log')
            for record in catalog['files']:
                current_schema = (record['eventSchemaFnv1a64'], record['metadataSchemaFnv1a64'],
                                  tuple((q['name'], q['schemaFnv1a64']) for q in record['qa']), record['compression'])
                if schema is None:
                    schema = current_schema
                require(schema == current_schema, 'Schema/compression differs between periods')
                require(record['path'] not in source_records, 'Repeated source catalog')
                for segment in record['segments']:
                    current_provenance = (segment['selectionFnv1a64'], segment['mixingFnv1a64'])
                    if provenance is None:
                        provenance = current_provenance
                    require(provenance == current_provenance, 'Producer provenance differs between periods')
                    key = (segment['high'], segment['low'])
                    require(key not in seen_segments, 'Producer segment repeated between periods')
                    seen_segments.add(key)
                source_records[record['path']] = record
            quota(config, directory / 'quota_before.json', remaining)
            status.update(phase='merging', updatedJst=now())
            state_json(status_path, status)
            merge = root_action('merge')
            status.update(phase='merged_output_validation', updatedJst=now())
            state_json(status_path, status)
            check = root_action('check')
            require(check['filesChecked'] == len(rows), 'Output check source count differs')
            require(check['sourceInputEvents'] == subset['inputEventsFromCompletionLogs'] and
                    check['sourceEvents'] == subset['acceptedEventsFromCompletionLogs'], 'Output event counters differ')
            for row in rows:
                fingerprint(row)
            # Linking refuses a pre-existing destination and uses no second data copy.
            os.link(str(partial), str(final))
            os.unlink(str(partial))
            info = final.stat()
            published = dict(group=name, path=str(final), sourceFiles=len(rows),
                             bytes=info.st_size, allocatedBytes=info.st_blocks * 512,
                             inode=info.st_ino, device=info.st_dev, mtimeNs=info.st_mtime_ns,
                             inputEvents=check['sourceInputEvents'], acceptedEvents=check['sourceEvents'],
                             producerSegments=check['producerSegments'], publishedJst=now(),
                             inputsListSha256=sha(input_list), checkSha256=sha(directory / 'check.json'),
                             qaNonidenticalBins=check['qaNonidenticalBins'])
            fresh_json(directory / 'published.json', published)
            status['completedGroups'].append(published)
            remaining -= sum(r['bytes'] for r in rows)
            status.update(phase='group_completed', updatedJst=now())
            state_json(status_path, status)
            print(now(), 'PUBLISHED', final, flush=True)
        ordered_records = [source_records[row['path']] for row in snapshot['paths']]
        global_catalog = dict(passed=True, fullEventReplayPerformed=False,
                              files=ordered_records, filesChecked=len(ordered_records),
                              sourceEvents=snapshot['acceptedEventsFromCompletionLogs'],
                              sourceInputEvents=snapshot['inputEventsFromCompletionLogs'])
        fresh_json(study / 'all_sources_catalog.json', global_catalog)
        run([sys.executable, config['catalogChecker'], '--snapshot', config['sourceSnapshot'],
             '--catalog', str(study / 'all_sources_catalog.json'), '--output', str(study / 'all_sources_validation.json'),
             '--expected-qa', str(config['expectedQa'])], root, study / 'all_sources_validation.log')
        for row in snapshot['paths']:
            fingerprint(row)
        for row in status['completedGroups']:
            fingerprint(row)
        require(sum(x['sourceFiles'] for x in status['completedGroups']) == len(sources), 'Published source coverage differs')
        require(sum(x['acceptedEvents'] for x in status['completedGroups']) == snapshot['acceptedEventsFromCompletionLogs'], 'Published event sum differs')
        with (output_dir / 'grouped_roots.list').open('x') as output:
            output.write(''.join(x['path'] + '\n' for x in status['completedGroups']))
        status.update(phase='complete', completedJst=now(), currentGroup=None,
                      producerSegments=len(seen_segments), fullEventReplayPerformed=False)
        state_json(status_path, status)
        fresh_json(study / 'completed.json', status)
        return 0
    except BaseException as error:
        status.update(phase='attention_required', error=str(error), stoppedJst=now())
        state_json(status_path, status)
        print(now(), 'STOPPED', repr(error), flush=True)
        raise


if __name__ == '__main__':
    sys.exit(main())
