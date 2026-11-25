#!/usr/bin/env python3
"""Verify seven submitted period cases, postprocess once, notify once."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import email.message
import email.policy
import email.utils
import fcntl
import getpass
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time
import types
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def load(path):
    with Path(path).open() as stream:
        return json.load(stream)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def resolve(value, base):
    path = Path(value)
    return (path if path.is_absolute() else base / path).resolve()


def inside(path, directory):
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def stat_record(path, allow_empty=False):
    st = path.stat()
    need(path.is_file() and (allow_empty or st.st_size > 0), 'Missing/empty file: ' + str(path))
    return dict(path=str(path), bytes=st.st_size, device=st.st_dev, inode=st.st_ino, mtimeNs=st.st_mtime_ns)


def same_stat(record):
    current = stat_record(Path(record['path']), allow_empty=record['bytes'] == 0)
    need(all(current[k] == record[k] for k in ['bytes', 'device', 'inode', 'mtimeNs']),
         'Pinned input stat changed: ' + record['path'])


def atomic_json(path, value):
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    with tmp.open('x') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(tmp), str(path))


def exclusive_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def pin(path):
    result = stat_record(path)
    result['sha256'] = sha(path)
    return result


def check_pins(pins, source_stats):
    for record in pins:
        same_stat(record)
        need(sha(record['path']) == record['sha256'], 'Pinned source/config SHA changed: ' + record['path'])
    for record in source_stats:
        same_stat(record)


def list_rows(path):
    return [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def native_inputs(meta, sources):
    prefix = meta['jobPrefix']
    need(re.fullmatch(r'[A-Za-z0-9_]+', prefix) is not None, 'Unsafe native job prefix')
    paths = list((PROJECT / 'job/run').glob(prefix + '_*.list'))
    found, records = {}, []
    for path in paths:
        match = re.fullmatch(re.escape(prefix) + r'_(0|[1-9][0-9]*)\.list', path.name)
        need(match is not None, 'Noncanonical native list filename')
        index = int(match.group(1))
        need(index not in found and 0 <= index < len(sources), 'Duplicate/out-of-range native period index')
        rows = list_rows(path)
        need(len(rows) == 1 and len(rows[0]) in (1, 2), 'Native worker does not have exactly one input')
        if len(rows[0]) == 2:
            need(rows[0][1].isdigit(), 'Invalid native event count field')
        source = Path(rows[0][0]).resolve()
        need(str(source) == str(Path(sources[index]['path']).resolve()), 'Native input order/source path mismatch')
        found[index] = source
        record = pin(path)
        record.update(period=index, source=str(source))
        records.append(record)
    need(set(found) == set(range(len(sources))), 'Native input coverage is not exactly periods 0..19')
    return sorted(records, key=lambda row: row['period'])


def uri_path(value):
    need(value.startswith('file:'), 'Only local file: paths are allowed')
    return Path(value[5:])


def build_run(case, directory, config_path, sources, input_hash, baseline=False):
    state_path = directory / 'launch_state.json'
    state, config = load(state_path), load(config_path)
    accepted = ['submitted_verified'] if not baseline else ['submitted_verified', 'monitor_started', 'notification_completed']
    need(state['phase'] in accepted and state['actualJobs'] == len(sources), 'Launch is not a verified 20-worker submission: ' + case['caseId'])
    need(state.get('submitReturnCode', 0) == 0, 'Nonzero submit return code')
    jobid = state['jobid']
    need(re.fullmatch(r'[A-Fa-f0-9]{32}', jobid) is not None, 'Invalid SUMS job ID')
    need(config['inputListSha256'] == input_hash and sha(config['inputList']) == input_hash, 'Launch input-list pin differs')
    configured_sources = config['sources']
    need(len(configured_sources) == len(sources), 'Launch source count differs')
    need([(str(Path(s['path']).resolve()), s['acceptedEvents']) for s in configured_sources] ==
         [(str(Path(s['path']).resolve()), s['acceptedEvents']) for s in sources], 'Launch source order/events differ')
    if baseline:
        expected_dir = Path(case['existingBaselineOutputDirectory']).resolve()
    else:
        expected_dir = Path(case['outputDirectory']).resolve()
        need(config['anaName'] == case['anaName'] and resolve(config['mainconf'], PROJECT) == resolve(case['mainconf'], PROJECT),
             'Case config anaName/mainconf differs')
    need(Path(config['outputDirectory']).resolve() == expected_dir, 'Output directory differs from case contract')
    runmeta_path = resolve(state['runmeta'], PROJECT)
    meta = load(runmeta_path)
    need(meta['jobid'] == jobid and Path(meta['projectRoot']).resolve() == PROJECT, 'Runmeta job/project mismatch')
    snapshot = meta['submit']['joblistSnapshot']
    xml_path = Path(snapshot['path'])
    need(sha(xml_path) == snapshot['sha256'] and xml_path.stat().st_size == snapshot['sizeBytes'], 'Submitted XML snapshot pin mismatch')
    xml = ET.parse(str(xml_path)).getroot()
    stdout = uri_path(xml.find('stdout').attrib['URL'])
    stderr = uri_path(xml.find('stderr').attrib['URL'])
    outputs = [node for node in xml.findall('output') if node.attrib['fromScratch'].endswith('.root')]
    need(len(outputs) == 1, 'Expected one ROOT output element')
    output = outputs[0]
    need(uri_path(output.attrib['toURL']).resolve() == expected_dir, 'XML output directory differs')
    template = Path(output.attrib['fromScratch']).name
    need('$JOBID' in template and '$JOBID' in str(stdout) and '$JOBID' in str(stderr), 'Unrecognized output/log filename template')
    expected_template = config['outputFileStem'] + '_$JOBID.root'
    need(template == expected_template, 'XML ROOT stem differs from launch config')
    native = native_inputs(meta, sources)
    roots = [expected_dir / template.replace('$JOBID', jobid + '_' + str(i)) for i in range(len(sources))]
    if baseline:
        need([str(p) for p in roots] == case['existingBaselineOutputs'], 'Baseline output order/job differs from pinned result')
    return dict(caseId=case['caseId'], baselineReuse=baseline, directory=str(directory), configPath=str(config_path),
                jobid=jobid, runmeta=str(runmeta_path), eventLog=str(resolve(state['eventLog'], PROJECT)),
                submitHost=config['submitHost'], nativeInputs=native, roots=[str(p) for p in roots],
                stdoutTemplate=str(stdout), stderrTemplate=str(stderr),
                staticPins=[pin(p) for p in [state_path, config_path, runmeta_path, xml_path]])


def event_parser(path):
    module = types.ModuleType('matrix_shared_event_parser')
    module.__file__ = str(path)
    exec(compile(path.read_text(), str(path), 'exec'), module.__dict__)
    return module.parse_events


def queue_snapshot(runs, timeout):
    hosts = {run['submitHost'] for run in runs}
    need(len(hosts) == 1, 'All cases must use one submit host')
    ids = [run['jobid'] for run in runs]
    expression = 'regexp("' + '|'.join(ids) + '",Out)'
    result = subprocess.run(['condor_q', '-name', next(iter(hosts)), '-constraint', expression,
                             '-af', 'ClusterId', 'ProcId', 'JobStatus', 'Out'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                            check=True, timeout=timeout)
    queues = {jobid: {} for jobid in ids}
    for line in result.stdout.splitlines():
        fields = line.split(None, 3)
        need(len(fields) == 4 and all(field.isdigit() for field in fields[:3]), 'Unrecognized condor_q row')
        matches = [jobid for jobid in ids if jobid in fields[3]]
        need(len(matches) == 1, 'Queue row does not map uniquely to a SUMS id')
        key = str(int(fields[0])) + '.' + str(int(fields[1]))
        need(key not in queues[matches[0]], 'Duplicate queue identity')
        queues[matches[0]][key] = int(fields[2])
    return queues


def poll(runs, parse_events, timeout):
    queues = queue_snapshot(runs, timeout)
    result = []
    for run in runs:
        event = Path(run['eventLog'])
        need(event.is_file(), 'Missing Condor event log: ' + str(event))
        submitted, terminal = parse_events(event.read_text(errors='replace'))
        need(len(submitted) <= len(run['roots']) and set(terminal) <= submitted, 'Unexpected Condor event identities/count')
        queue = queues[run['jobid']]
        need(len(queue) <= len(run['roots']), 'More queue workers than expected')
        if len(submitted) == len(run['roots']):
            need(set(queue) <= submitted, 'Queue identities absent from complete submitted-event proof')
        need(not any(state == 5 for state in queue.values()), 'Held worker in ' + run['caseId'])
        bad = {key: value for key, value in terminal.items() if value != 'exit_0'}
        need(not bad, 'Nonzero/removed/unknown worker termination in ' + run['caseId'] + ': ' + repr(bad))
        ready = len(submitted) == len(run['roots']) and set(terminal) == submitted and not queue
        result.append(dict(caseId=run['caseId'], jobid=run['jobid'], submitted=len(submitted), terminal=len(terminal),
                           outcomes=dict(Counter(terminal.values())), queueStates=dict(Counter(queue.values())),
                           allTerminal=ready, eventLog=pin(event)))
    return dict(checkedAtUtc=now(), allTerminal=all(row['allTerminal'] for row in result), cases=result)


REPLAY = re.compile(r'^FemtoLambda tree replay job=(\S+) requested=(\d+) completed=(\d+) available=(\d+) status=(-?\d+) RealTime=([0-9.eE+\-]+) CpuTime=([0-9.eE+\-]+)\s*$')
MAKER = re.compile(r'^\[StFemtoMaker Lambda\] inputEvents=(\d+) reconstructedEvents=(\d+) selectedLambda=(\d+) mixingVzClamped=(\d+) status=(\w+)\s*$')


def scan_log(path, fatal_patterns, stdout=False):
    before = stat_record(path, allow_empty=not stdout)
    h = hashlib.sha256()
    result = dict(path=str(path), bytes=before['bytes'], fatalLines=[], replay=[], maker=[], normalNativeExit=False,
                  chainCounts=[], endRaw=None)
    chain_next = False
    with path.open('rb') as stream:
        for block in stream:
            h.update(block)
            line = block.decode('utf-8', errors='replace').strip()
            if any(pattern.search(line) for pattern in fatal_patterns) and len(result['fatalLines']) < 12:
                result['fatalLines'].append(line[:1200])
            if not stdout:
                continue
            match = REPLAY.fullmatch(line)
            if match:
                job, requested, completed, available, status, wall, cpu = match.groups()
                result['replay'].append(dict(job=job, requested=int(requested), completed=int(completed),
                                             available=int(available), status=int(status), wallSeconds=float(wall), cpuSeconds=float(cpu)))
            match = MAKER.fullmatch(line)
            if match:
                read, done, selected, clamped, status = match.groups()
                result['maker'].append(dict(inputEvents=int(read), reconstructedEvents=int(done), selectedLambda=int(selected),
                                            mixingVzClamped=int(clamped), status=status))
            if 'nStOK' in line and 'nStWarn' in line and 'nStEOF' in line and 'nStErr' in line and 'nStFatal' in line:
                chain_next = True
            elif chain_next:
                match = re.fullmatch(r'StInfo:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
                if match:
                    result['chainCounts'].append(dict(zip(['nStOK', 'nStWarn', 'nStEOF', 'nStErr', 'nStFatal'], map(int, match.groups()))))
                chain_next = False
            if 'has ended, exiting normally' in line:
                result['normalNativeExit'] = True
            match = re.search(r'We END on \[(.*)\]', line)
            if match:
                result['endRaw'] = match.group(1)
    same_stat(before)
    result['sha256'] = h.hexdigest()
    result['stat'] = before
    return result


def verify_outputs(run, sources, fatal_patterns):
    workers, output_stats = [], []
    expected_indices = set(range(len(sources)))
    root_directory = Path(run['roots'][0]).parent
    # Enumerate only this exact SUMS stem, never unrelated historical outputs.
    stem = Path(run['roots'][0]).name.rsplit(run['jobid'] + '_', 1)[0]
    found = set()
    for path in root_directory.glob(stem + run['jobid'] + '_*.root'):
        match = re.fullmatch(re.escape(stem + run['jobid']) + r'_(0|[1-9][0-9]*)\.root', path.name)
        need(match is not None, 'Noncanonical case output filename')
        found.add(int(match.group(1)))
    need(found == expected_indices, 'Case output ROOT indices are not exactly 0..19: ' + run['caseId'])
    for i, source in enumerate(sources):
        subjob = run['jobid'] + '_' + str(i)
        stdout = Path(run['stdoutTemplate'].replace('$JOBID', subjob))
        stderr = Path(run['stderrTemplate'].replace('$JOBID', subjob))
        out = scan_log(stdout, fatal_patterns, stdout=True)
        err = scan_log(stderr, fatal_patterns)
        need(not out['fatalLines'] and not err['fatalLines'], 'Fatal/error signature in native logs: ' + subjob)
        need(len(out['replay']) == len(out['maker']) == len(out['chainCounts']) == 1,
             'Missing/duplicate completion markers: ' + subjob)
        marker, maker, counts = out['replay'][0], out['maker'][0], out['chainCounts'][0]
        expected = source['acceptedEvents']
        need(marker['job'] == subjob and marker['requested'] == marker['completed'] == marker['available'] == expected and marker['status'] == 0,
             'Replay marker counts/status differ from stored period: ' + subjob)
        need(maker['inputEvents'] == maker['reconstructedEvents'] == expected and maker['status'] == 'completed', 'Maker incomplete/count mismatch: ' + subjob)
        need(counts == dict(nStOK=expected, nStWarn=0, nStEOF=0, nStErr=0, nStFatal=0), 'Nonzero chain warning/error/EOF or incomplete nStOK: ' + subjob)
        need(out['normalNativeExit'] and out['endRaw'] is not None, 'Missing normal wrapper exit/END: ' + subjob)
        root = Path(run['roots'][i])
        stat = stat_record(root)
        need(stat['bytes'] > 4, 'Empty ROOT output')
        with root.open('rb') as stream:
            need(stream.read(4) == b'root', 'Invalid ROOT header: ' + str(root))
        same_stat(stat)
        stat.update(caseId=run['caseId'], period=i)
        output_stats.append(stat)
        workers.append(dict(period=i, source=source['path'], expectedStoredEvents=expected, root=stat,
                            stdout=out, stderr=err, completion=marker))
    return dict(caseId=run['caseId'], jobid=run['jobid'], passed=True, workers=workers,
                completedEvents=sum(row['completion']['completed'] for row in workers), outputStats=output_stats,
                totalOutputBytes=sum(row['bytes'] for row in output_stats))


def send_once(state, settings, successful):
    attempt, sent = state / 'delivery_attempt.json', state / 'sent.json'
    need(not attempt.exists() and not sent.exists(), 'Delivery already attempted; automatic resend disabled')
    message = email.message.EmailMessage(policy=email.policy.SMTP)
    message['From'] = '{}@{}'.format(getpass.getuser(), socket.getfqdn())
    message['To'] = settings['recipient']
    message['Date'] = email.utils.format_datetime(datetime.now().astimezone())
    message['Subject'] = '[STAR解析 完了]' if successful else '[STAR解析 要確認]'
    message['Message-ID'] = email.utils.make_msgid()
    message.set_content(('ご依頼の解析jobの終了・出力回収と、条件比較図の作成を確認しました。\n' if successful else
                         'ご依頼の解析または条件比較の確認で、確認が必要な点があります。\n') +
                        '詳細はこのチャットをご確認ください。\n')
    data = message.as_bytes()
    with (state / 'notification.eml').open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    exclusive_json(attempt, dict(attemptedAtUtc=now(), recipient=settings['recipient'],
                                messageId=str(message['Message-ID']), messageSha256=hashlib.sha256(data).hexdigest(),
                                kind='success' if successful else 'attention', automaticResend=False))
    try:
        subprocess.run([settings['sendmail'], '-i', '-t'], input=data, check=True, timeout=settings['commandTimeoutSeconds'])
    except Exception as error:
        atomic_json(state / 'delivery_uncertain.json', dict(atUtc=now(), error=str(error), automaticResend=False))
        raise RuntimeError('Local sendmail outcome uncertain; automatic resend disabled')
    exclusive_json(sent, dict(sentAtUtc=now(), recipient=settings['recipient'], mailStatus='accepted_by_local_sendmail',
                             messageId=str(message['Message-ID']), kind='success' if successful else 'attention'))


def validate_settings(args):
    """Pre-submission checks without launch states, Condor, outputs or mail."""
    need(args.matrix is not None, '--validate-settings requires --matrix')
    sp, mp = args.settings.resolve(), args.matrix.resolve()
    settings, matrix = load(sp), load(mp)
    need(settings['schemaVersion'] == matrix['schemaVersion'] == 1 and settings['pollSeconds'] == 300 and settings['settleSeconds'] >= 0, 'Invalid schema/poll/settle settings')
    need(re.fullmatch(r'[^\s<>@]+@[^\s<>@]+', settings['recipient']) is not None and Path(settings['sendmail']).is_absolute(), 'Invalid explicit recipient/sendmail')
    need(settings['commandTimeoutSeconds'] > 0 and settings['postprocessTimeoutSeconds'] > 0, 'Invalid timeouts')
    for pattern in settings['fatalPatterns']:
        re.compile(pattern)
    cp = resolve(settings['cases'], sp.parent)
    cases = load(cp)
    by_id = {c['caseId']: c for c in cases['cases']}
    need(cases['schemaVersion'] == 1 and cases['expectedPeriods'] == 20 and len(by_id) == 8 and len(cases['productionCases']) == 7, 'Invalid case matrix')
    need({(c['epBins'], c['cutEnabled']) for c in cases['cases']} == {(n, cut) for n in [1, 8, 16, 32] for cut in [False, True]}, 'Invalid EP/cut combinations')
    need(Path(matrix['projectRoot']).resolve() == PROJECT and Path(matrix['studyDirectory']).resolve() == HERE and resolve(matrix['casesManifest'], mp.parent) == cp and matrix['casesManifestSha256'] == sha(cp), 'Matrix project/study/case-manifest pin differs')
    rows = matrix['cases']
    need(len(rows) == len({r['caseId'] for r in rows}) == 7 and {r['caseId'] for r in rows} == set(cases['productionCases']), 'Matrix must contain exactly seven new cases')
    ip = resolve(cases['inputManifest'], cp.parent)
    lp = resolve(cases['inputList'], cp.parent)
    need(sha(ip) == cases['inputManifestSha256'] and sha(lp) == cases['inputListSha256'], 'Input manifest/list SHA differs')
    sources = load(ip)['sourceFiles']
    need(len(sources) == 20 and sum(s['acceptedEvents'] for s in sources) == cases['storedEvents'], 'Source event coverage differs')
    need([str(Path(r[0]).resolve()) for r in list_rows(lp)] == [str(Path(s['path']).resolve()) for s in sources], 'Source order differs')
    source_stats = [dict(path=s['path'], **{k: s[k] for k in ['bytes', 'device', 'inode', 'mtimeNs']}) for s in sources]
    for record in source_stats:
        same_stat(record)
    paths = [Path(__file__).resolve(), sp, mp, cp, ip, lp, resolve(settings['sharedEventParser'], sp.parent), resolve(settings['postprocessSettings'], sp.parent), HERE/'postprocess_cf.py', HERE/'ExtractPeriodCF.C', HERE/'DrawPeriodCF.C']
    paths += [resolve(v, sp.parent) for v in settings['additionalPinFiles']]
    pins = {str(p): pin(p) for p in paths}
    for row in rows:
        case = by_id[row['caseId']]
        directory, config_path = resolve(row['launchDirectory'], mp.parent), resolve(row['config'], mp.parent)
        need(inside(directory, HERE/'production') and inside(config_path, directory), 'Launch/config outside case directory')
        config = load(config_path)
        need(config['schemaVersion'] == 1 and config['anaName'] == case['anaName'] and resolve(config['mainconf'], PROJECT) == resolve(case['mainconf'], PROJECT) and Path(config['outputDirectory']).resolve() == Path(case['outputDirectory']).resolve(), 'Launch config identity differs')
        need(config['inputListSha256'] == cases['inputListSha256'] and sha(config['inputList']) == cases['inputListSha256'], 'Launch input-list SHA differs')
        need([(str(Path(s['path']).resolve()), s['acceptedEvents']) for s in config['sources']] == [(str(Path(s['path']).resolve()), s['acceptedEvents']) for s in sources], 'Launch source order/events differ')
        need(config['fingerprints'], 'Launch config has no fingerprints')
        for value, expected in config['fingerprints'].items():
            p = resolve(value, PROJECT)
            if str(p) not in pins:
                pins[str(p)] = pin(p)
            need(pins[str(p)]['sha256'] == expected, 'Launch fingerprint differs: '+str(p))
        pins[str(config_path)] = pin(config_path)
    post_settings = load(resolve(settings['postprocessSettings'], sp.parent))
    centrality_ids = {c['id'] for c in post_settings['centralities']}
    need(post_settings['requiredValidCentralities'] and set(post_settings['requiredValidCentralities']) <= centrality_ids,
         'Missing/unknown required main CF centralities')
    return dict(passed=True, validateSettingsOnly=True, launchStatesRead=False, schedulerQueried=False, outputsRead=False, mailSent=False, postprocessRun=False, analysisSubmitted=False, newCases=7, newWorkers=140, matrixSha256=sha(mp), pins=list(pins.values()), sourceStats=source_stats)


def prepare(args):
    settings_path = args.settings.resolve()
    settings = load(settings_path)
    need(settings['schemaVersion'] == 1 and settings['pollSeconds'] == 300 and settings['settleSeconds'] >= 0,
         'Unsupported settings/poll interval')
    need(re.fullmatch(r'[^\s<>@]+@[^\s<>@]+', settings['recipient']) is not None and
         Path(settings['sendmail']).is_absolute(), 'Invalid explicit recipient/sendmail path')
    cases_path = resolve(settings['cases'], settings_path.parent)
    cases = load(cases_path)
    need(cases['schemaVersion'] == 1 and cases['expectedPeriods'] == 20, 'Unsupported case schema/period count')
    by_id = {case['caseId']: case for case in cases['cases']}
    need(len(by_id) == 8 and len(cases['productionCases']) == 7, 'Expected eight cases / seven new productions')
    manifest_path = resolve(cases['inputManifest'], cases_path.parent)
    need(sha(manifest_path) == cases['inputManifestSha256'], 'Input manifest SHA mismatch')
    manifest = load(manifest_path)
    sources = manifest['sourceFiles']
    need(len(sources) == cases['expectedPeriods'] and sum(source['acceptedEvents'] for source in sources) == cases['storedEvents'], 'Source event coverage mismatch')
    source_stats = [dict(path=source['path'], **{key: source[key] for key in ['bytes', 'device', 'inode', 'mtimeNs']}) for source in sources]
    for record in source_stats:
        same_stat(record)
    parser_path = resolve(settings['sharedEventParser'], settings_path.parent)
    post_settings = resolve(settings['postprocessSettings'], settings_path.parent)
    pins = [pin(path) for path in [Path(__file__).resolve(), settings_path, cases_path, manifest_path, parser_path, post_settings,
                                  HERE / 'postprocess_cf.py', HERE / 'ExtractPeriodCF.C', HERE / 'DrawPeriodCF.C']]
    pins += [pin(resolve(value, settings_path.parent)) for value in settings['additionalPinFiles']]
    baseline = by_id[cases['baselineReuseCaseId']]
    need(baseline['baselineReuse'] and baseline['epBins'] == 8 and baseline['cutEnabled'] is False, 'Baseline contract differs')
    baseline_dir = resolve(settings['baselineLaunchDirectory'], settings_path.parent)
    baseline_config = resolve(settings['baselineLaunchConfig'], settings_path.parent)
    baseline_run = build_run(baseline, baseline_dir, baseline_config, sources, cases['inputListSha256'], baseline=True)
    need(baseline_run['jobid'] == cases['baselineJobid'], 'Pinned baseline SUMS id differs')
    baseline_proof_path = resolve(settings['baselineInputProvenance'], settings_path.parent)
    baseline_completed_path = resolve(settings['baselineCompleted'], settings_path.parent)
    baseline_proof, baseline_completed = load(baseline_proof_path), load(baseline_completed_path)
    need(baseline_completed['passed'] and baseline_proof['cases'] == [baseline['caseId']], 'Baseline postprocess proof missing/mismatched')
    baseline_stats = baseline_proof['sourceStats']
    need([row['path'] for row in baseline_stats] == baseline_run['roots'], 'Baseline stat proof output order differs')
    for record in baseline_stats:
        same_stat(record)
    source_stats += baseline_stats
    pins += [pin(baseline_proof_path), pin(baseline_completed_path)]
    if args.check_baseline:
        runs = [baseline_run]
    else:
        need(args.matrix is not None, '--matrix is required except --check --check-baseline')
        matrix_path = args.matrix.resolve()
        matrix = load(matrix_path)
        need(matrix['schemaVersion'] == 1, 'Unsupported matrix schema')
        need(Path(matrix['projectRoot']).resolve() == PROJECT and Path(matrix['studyDirectory']).resolve() == HERE and resolve(matrix['casesManifest'], matrix_path.parent) == cases_path and matrix['casesManifestSha256'] == sha(cases_path), 'Matrix project/study/case-manifest pin differs')
        rows = matrix['cases']
        need(len(rows) == 7 and len({row['caseId'] for row in rows}) == 7 and
             {row['caseId'] for row in rows} == set(cases['productionCases']), 'Matrix must contain exactly the seven new cases')
        pins.append(pin(matrix_path))
        runs = []
        for row in rows:
            directory = resolve(row['launchDirectory'], matrix_path.parent)
            config_path = resolve(row['config'], matrix_path.parent)
            need(inside(directory, HERE / 'production') and inside(config_path, directory), 'Production launch/config is outside its case directory')
            runs.append(build_run(by_id[row['caseId']], directory, config_path, sources, cases['inputListSha256']))
        need(len({run['jobid'] for run in runs}) == 7 and cases['baselineJobid'] not in {run['jobid'] for run in runs}, 'Duplicate/new-baseline SUMS ids')
    for run in runs + ([] if args.check_baseline else [baseline_run]):
        pins += run['staticPins'] + run['nativeInputs']
    unique = {}
    for record in pins:
        if record['path'] in unique:
            need(record['sha256'] == unique[record['path']]['sha256'], 'Inconsistent duplicate pin')
        unique[record['path']] = record
    return settings, cases, sources, list(unique.values()), source_stats, runs, baseline_run, event_parser(parser_path), post_settings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--settings', type=Path, required=True)
    ap.add_argument('--matrix', type=Path)
    ap.add_argument('--state-dir', type=Path)
    ap.add_argument('--check', action='store_true', help='One read-only pass; no files/mail/postprocess are created')
    ap.add_argument('--validate-settings', action='store_true', help='Pre-submit config/pin/source-stat validation; no launch states/queue/outputs/mail')
    ap.add_argument('--check-baseline', action='store_true', help='With --check, independently verify only the prior actual EP8 baseline run')
    args = ap.parse_args()
    need(not args.validate_settings or not (args.check or args.check_baseline), 'Validation cannot be combined with check modes')
    if args.validate_settings:
        print(json.dumps(validate_settings(args), ensure_ascii=False, indent=2, allow_nan=False))
        return
    need(not args.check_baseline or args.check, '--check-baseline is allowed only with --check')
    if args.check:
        settings, cases, sources, pins, source_stats, runs, baseline, parser, post_settings = prepare(args)
        status = poll(runs, parser, settings['commandTimeoutSeconds'])
        status.update(checkOnly=True, mailSent=False, postprocessRun=False, analysisSubmitted=False)
        if status['allTerminal']:
            patterns = [re.compile(pattern) for pattern in settings['fatalPatterns']]
            status['completion'] = [verify_outputs(run, sources, patterns) for run in runs]
        check_pins(pins, source_stats)
        print(json.dumps(status, ensure_ascii=False, indent=2, allow_nan=False))
        return
    need(args.state_dir is not None, 'Normal watch mode requires --state-dir')
    state = args.state_dir.resolve()
    need(inside(state, HERE) and state != HERE, 'State directory must be a subdirectory of this study')
    state.mkdir(parents=True, exist_ok=True)
    lock = (state / 'watcher.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    settings = load(args.settings.resolve())
    need(not (state / 'sent.json').exists() and not (state / 'delivery_attempt.json').exists(), 'Notification already sent/attempted; no automatic restart')
    need(not (state / 'attention.json').exists(), 'Attention state requires manual review; automatic restart disabled')
    try:
        settings, cases, sources, pins, source_stats, runs, baseline, parser, post_settings = prepare(args)
        startup_path = state / 'startup_pins.json'
        if startup_path.exists():
            previous = load(startup_path)
            need({p['path']: p['sha256'] for p in previous['pins']} == {p['path']: p['sha256'] for p in pins} and
                 previous['sourceStats'] == source_stats, 'Watcher restart input pins differ')
        else:
            exclusive_json(startup_path, dict(startedAtUtc=now(), pins=pins, sourceStats=source_stats,
                                              runs=[{key: run[key] for key in ['caseId', 'jobid', 'runmeta', 'eventLog']} for run in runs],
                                              baselineJobid=baseline['jobid'], automaticResubmit=False, automaticHadd=False,
                                              automaticRetry=False))
        atomic_json(state / 'pid.json', dict(pid=os.getpid(), host=socket.getfqdn(), startedAtUtc=now()))
        patterns = [re.compile(pattern) for pattern in settings['fatalPatterns']]
        first_terminal = None
        while True:
            check_pins(pins, source_stats)
            status = poll(runs, parser, settings['commandTimeoutSeconds'])
            status.update(phase='waiting', expectedNewCases=7, expectedNewWorkers=140, reusedBaselineWorkers=20,
                          automaticRetry=False, automaticResubmit=False, automaticHadd=False)
            if status['allTerminal']:
                if first_terminal is None:
                    first_terminal = time.monotonic()
                status['phase'] = 'all_terminal_settling'
                status['settleElapsedSeconds'] = time.monotonic() - first_terminal
                if status['settleElapsedSeconds'] >= settings['settleSeconds']:
                    break
            else:
                first_terminal = None
            atomic_json(state / 'status.json', status)
            print(json.dumps(status, ensure_ascii=False), flush=True)
            time.sleep(settings['pollSeconds'])
        reports = [verify_outputs(run, sources, patterns) for run in runs + [baseline]]
        need(all(report['completedEvents'] == cases['storedEvents'] for report in reports), 'Per-case total stored events differ')
        check_pins(pins, source_stats)
        for report in reports:
            for record in report['outputStats']:
                same_stat(record)
        exclusive_json(state / 'completion_verified.json', dict(passed=True, verifiedAtUtc=now(), cases=reports,
                                                               eventsPerCase=cases['storedEvents'], treeEntriesRead=0))
        list_directory = state / 'output_lists'
        list_directory.mkdir()
        mapping = dict(schemaVersion=1, cases={})
        for run in runs + [baseline]:
            path = list_directory / (run['caseId'] + '.list')
            with path.open('x') as stream:
                stream.write('\n'.join(run['roots']) + '\n')
            mapping['cases'][run['caseId']] = dict(outputList=str(path))
        mapping_path = state / 'case_outputs.json'
        exclusive_json(mapping_path, mapping)
        result_dir = state / 'cf_comparison'
        attempt = state / 'postprocess_attempt.json'
        need(not attempt.exists(), 'Postprocess already attempted; automatic retry disabled')
        command = [sys.executable, str(HERE / 'postprocess_cf.py'), '--cases', str(resolve(settings['cases'], args.settings.resolve().parent)),
                   '--settings', str(post_settings), '--outputs', str(mapping_path), '--mode', 'compare', '--output-dir', str(result_dir)]
        exclusive_json(attempt, dict(startedAtUtc=now(), command=command, automaticRetry=False))
        atomic_json(state / 'status.json', dict(phase='postprocessing', atUtc=now(), allTerminal=True, completionVerified=True))
        with (state / 'postprocess.log').open('x') as stream:
            subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT, check=True,
                           timeout=settings['postprocessTimeoutSeconds'])
        completed = load(result_dir / 'completed.json')
        need(completed['passed'] and completed['mode'] == 'compare' and set(completed['cases']) == {case['caseId'] for case in cases['cases']} and completed['sources'] == 160,
             'Postprocess did not verify all eight cases')
        comparison = load(result_dir / 'cf_comparison.json')
        post_config = load(post_settings)
        for case in cases['cases']:
            for cent_id in post_config['requiredValidCentralities']:
                main_curves = [c for c in comparison['curves'] if c['caseId'] == case['caseId'] and c['centralityId'] == cent_id and c['plane'] == -1]
                need(len(main_curves) == 1 and main_curves[0]['available'] and main_curves[0]['normSE'] > 0 and main_curves[0]['normME'] > 0,
                     'Required main CF missing/unnormalized: ' + case['caseId'] + '/' + cent_id)
        check_pins(pins, source_stats)
        for report in reports:
            for record in report['outputStats']:
                same_stat(record)
        exclusive_json(state / 'postprocess_verified.json', dict(passed=True, finishedAtUtc=now(), completion=pin(result_dir / 'completed.json'),
                                                                comparison=pin(result_dir / 'cf_comparison.json')))
        atomic_json(state / 'status.json', dict(phase='completed_verified', atUtc=now(), newCases=7, newWorkers=140,
                                               baselineReused=True, comparisonDirectory=str(result_dir), physicalSignificanceAssessed=False))
        send_once(state, settings, True)
        print('MATRIX_COMPLETION_AND_COMPARISON_PASS mail=accepted_by_local_sendmail', flush=True)
    except Exception as error:
        attention = dict(phase='attention', atUtc=now(), error=str(error), automaticRetry=False,
                         automaticResubmit=False, automaticHadd=False, physicsCompletionClaim=False)
        atomic_json(state / 'attention.json', attention)
        atomic_json(state / 'status.json', attention)
        if not (state / 'delivery_attempt.json').exists() and not (state / 'sent.json').exists():
            try:
                send_once(state, settings, False)
            except Exception as delivery_error:
                atomic_json(state / 'notification_error.json', dict(atUtc=now(), error=str(delivery_error), automaticResend=False))
        raise


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('MATRIX_WATCH_FAILED: ' + str(error), file=sys.stderr, flush=True)
        sys.exit(1)
