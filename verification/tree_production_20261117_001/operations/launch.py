#!/usr/bin/env python3
"""Submit each approved campaign once and persist its completion-mail monitor."""
import argparse
from collections import Counter
from datetime import datetime, timezone, timedelta
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time


def now():
    return datetime.now(timezone(timedelta(hours=9))).isoformat()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(config):
    if socket.getfqdn() != config['submitHost']:
        raise RuntimeError('wrong submission/notification host')
    for path, expected in config['fingerprints'].items():
        if sha(path) != expected:
            raise RuntimeError('validated file changed: ' + path)


def coverage(root, entry, jobid, config):
    pattern = re.compile(re.escape(entry['anaName'] + jobid) + r'_(\d+)\.list$')
    paths = [f for f in (root / 'job/run').glob(entry['anaName'] + jobid + '_*.list') if pattern.fullmatch(f.name)]
    if not paths:
        raise RuntimeError('no generated worker input lists')
    paths.sort(key=lambda f: int(pattern.fullmatch(f.name).group(1)))
    expected = [Path(line.split()[0]).name for line in Path(config['inputList']).read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    groups = []
    for f in paths:
        names = [Path(line.split()[0]).name for line in f.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
        if not 1 <= len(names) <= 25:
            raise RuntimeError('invalid worker file count: ' + str(f))
        groups.append(names)
    actual = [name for group in groups for name in group]
    if len(expected) != config['expectedFiles'] or Counter(actual) != Counter(expected):
        raise RuntimeError('submitted inputs do not match the frozen covered-run list')
    signatures = [hashlib.sha256(json.dumps(group, separators=(',', ':')).encode()).hexdigest() for group in groups]
    return dict(workerJobs=len(groups), inputFiles=len(actual), orderedGroupSignatures=signatures,
                lists=[str(f) for f in paths], maxFiles=max(map(len, groups)), groupingOrder='within-worker order preserved')


def event_log(root, entry, jobid):
    values = set()
    files = list((root / 'job/run').glob(entry['anaName'] + jobid + '*.condor'))
    for f in files[:10]:
        for line in f.read_text(errors='replace').splitlines():
            match = re.match(r'\s*log\s*=\s*(.+?)\s*$', line, re.I)
            if match:
                values.add(match.group(1).strip().strip('"'))
    if len(values) == 1:
        value = values.pop()
        if jobid not in value or '$' in value:
            raise RuntimeError('unresolved Condor event-log path')
        return value
    if values:
        raise RuntimeError('more than one Condor event log; manual audit required')
    candidate = Path('/tmp/yichikawa') / ('sched' + jobid + '.condor.log')
    if not candidate.exists():
        raise RuntimeError('cannot locate actual Condor event log')
    return str(candidate)


def one(config, key):
    entry = next(e for e in config['entries'] if e['key'] == key)
    root, directory = Path(config['projectRoot']), Path(entry['stateDirectory'])
    lock = (directory / 'launcher.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    state_path = directory / 'launch_state.json'
    state = dict(phase='preflight', startedAtJst=now(), pid=os.getpid(), host=socket.getfqdn(),
                 automaticResubmit=False, automaticMerge=False, mainconf=entry['mainconf'], anaName=entry['anaName'])
    save(state_path, state)
    try:
        verify(config)
        with (directory / 'submit_attempt.json').open('x') as f:
            json.dump(dict(attemptedAtJst=now(), joblist=entry['joblist'], mainconf=entry['mainconf']), f, indent=2)
        state.update(phase='submitting', submitStartedAtJst=now()); save(state_path, state)
        with (directory / 'submit.log').open('w') as log:
            result = subprocess.run(['/bin/bash', str(root / 'job/run/submit.sh'), entry['joblist']], cwd=str(root / 'job/run'), stdout=log, stderr=subprocess.STDOUT)
        text = (directory / 'submit.log').read_text(errors='replace')
        ids = set(re.findall(re.escape(entry['anaName']) + r'([0-9A-Fa-f]{32})\.report', text))
        state.update(submitReturnCode=result.returncode, submitReturnedAtJst=now(), possibleReportIds=sorted(ids)); save(state_path, state)
        if result.returncode != 0:
            raise RuntimeError('submission returned nonzero; no automatic retry')
        if len(ids) != 1:
            raise RuntimeError('one unambiguous SUMS ID required; no automatic retry')
        jobid = ids.pop()
        runmeta = root / 'job/run/runmeta' / ('runmeta_' + entry['anaName'] + '_' + jobid + '.json')
        meta = json.loads(runmeta.read_text())
        if meta['jobid'] != jobid or meta['submit']['mainconf'] != entry['mainconf'] or meta['submit']['host'] != config['submitHost']:
            raise RuntimeError('runmeta identity/mainconf/host mismatch')
        proof = coverage(root, entry, jobid, config); save(directory / 'input_coverage.json', proof)
        event = event_log(root, entry, jobid)
        notify = directory / 'notify'
        watcher = [config['python'], '-u', str(root / 'script/watch_job_notify.py'), '--runmeta', str(runmeta),
                   '--event-log', event, '--email', config['email'], '--email-content', config['emailContent'],
                   '--poll-sec', str(config['pollSeconds']), '--settle-sec', str(config['settleSeconds']),
                   '--state-dir', str(notify), '--watch']
        state.update(phase='monitor_starting', jobid=jobid, runmeta=str(runmeta), eventLog=event,
                     actualJobs=proof['workerJobs'], inputFiles=proof['inputFiles'], monitorCommand=watcher)
        save(state_path, state)
        with (directory / 'completion_email_watch.log').open('w') as log:
            monitor = subprocess.Popen(watcher, cwd=str(root), stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        state.update(phase='monitor_started', monitorPid=monitor.pid, monitorStartedAtJst=now()); save(state_path, state)
        code = monitor.wait()
        if code != 0 or not (notify / 'sent.json').is_file():
            raise RuntimeError('monitor stopped without confirmed local-mail acceptance')
        state.update(phase='notification_completed', completedAtJst=now(), sentRecord=str(notify / 'sent.json')); save(state_path, state)
        return 0
    except Exception as error:
        state.update(phase='attention_required', error=str(error), failedAtJst=now()); save(state_path, state)
        print(key + ': ' + str(error), flush=True)
        return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--entry')
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if args.entry:
        return one(config, args.entry)
    root, directory = Path(config['projectRoot']), Path(config['directory'])
    lock = (directory / 'launcher.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    verify(config)
    for entry in config['entries']:
        if (Path(entry['stateDirectory']) / 'submit_attempt.json').exists():
            raise RuntimeError('prior submission intent exists; inspect before any restart')
        if list((root / 'job/run/runmeta').glob('runmeta_' + entry['anaName'] + '_*.json')) or list((root / 'job/run').glob(entry['anaName'] + '*.report')):
            raise RuntimeError('matching submission artifacts already exist')
    quota = subprocess.run(config['quotaCommand'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    (directory / 'quota_before_launch.log').write_text(quota.stdout + '\nSTDERR\n' + quota.stderr)
    header = next(line for line in quota.stdout.splitlines() if line.startswith('mmlsquota:user:HEADER:'))
    row = next(line for line in quota.stdout.splitlines() if line.startswith('mmlsquota:user:0:'))
    fields = dict(zip(header.split(':')[3:], row.split(':')[3:]))
    free = (int(fields['blockQuota']) - int(fields['blockUsage']) - int(fields['blockInDoubt'])) * 1024
    if quota.returncode or free < config['minimumSoftHeadroomBytes']:
        raise RuntimeError('quota preflight does not meet the approved two-campaign reserve')
    state = dict(phase='launching', startedAtJst=now(), pid=os.getpid(), host=socket.getfqdn(), softHeadroomBytes=free,
                 automaticResubmit=False, automaticMerge=False, entries={}, submitted=False)
    save(directory / 'launch_state.json', state)
    processes = {}
    for entry in config['entries']:
        log = (Path(entry['stateDirectory']) / 'launcher.log').open('w')
        processes[entry['key']] = subprocess.Popen([config['python'], '-u', str(Path(__file__).resolve()), '--config', str(args.config.resolve()), '--entry', entry['key']], cwd=str(root), stdout=log, stderr=subprocess.STDOUT)
        log.close()
        state['entries'][entry['key']] = dict(pid=processes[entry['key']].pid)
    save(directory / 'launch_state.json', state)
    audited = False
    while True:
        entries = {}
        for entry in config['entries']:
            path = Path(entry['stateDirectory']) / 'launch_state.json'
            if path.exists():
                entries[entry['key']] = json.loads(path.read_text())
        state['entries'] = entries
        if len(entries) == len(config['entries']) and all(e.get('jobid') for e in entries.values()):
            state['submitted'] = True; state['phase'] = 'submitted_and_monitors_running'
            if not audited:
                proofs = [json.loads((Path(e['stateDirectory']) / 'input_coverage.json').read_text()) for e in config['entries']]
                equal = Counter(proofs[0]['orderedGroupSignatures']) == Counter(proofs[1]['orderedGroupSignatures'])
                state['sameOrderedWorkerGroups'] = equal
                save(directory / 'grouping_comparison.json', dict(passed=equal, atJst=now(), workerJobs=[p['workerJobs'] for p in proofs], inputFiles=[p['inputFiles'] for p in proofs]))
                audited = True
        if state.get('sameOrderedWorkerGroups') is False:
            state['phase'] = 'submitted_grouping_attention'
        state['updatedAtJst'] = now(); save(directory / 'launch_state.json', state)
        if all(process.poll() is not None for process in processes.values()):
            state['phase'] = 'notification_completed' if all(process.returncode == 0 for process in processes.values()) else 'attention_required'
            state['finishedAtJst'] = now(); save(directory / 'launch_state.json', state)
            return 0 if state['phase'] == 'notification_completed' else 1
        time.sleep(30)


if __name__ == '__main__':
    sys.exit(main())
