#!/usr/bin/env python3
"""Submit a validated matrix once and remain in its aggregate watcher."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone
E = Path(__file__).resolve().parent
ROOT = E.parents[1]

def now():
    return datetime.now(timezone.utc).isoformat()

def save(name, value):
    temp = E / ('.'+name+'.'+str(os.getpid())+'.tmp')
    with temp.open('x') as f:
        json.dump(value,f,indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    temp.replace(E/name)

def main():
    lock = (E/'production_matrix.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (E/'matrix_launch_attempt.json').exists(), 'No automatic repeat'
    matrix_path = E/'production_matrix.json'
    matrix = json.loads(matrix_path.read_text())
    assert len(matrix['cases']) == 7 and matrix['plannedProductionJobs'] == 140
    for item in matrix['cases']:
        d = Path(item['launchDirectory'])
        assert d.parent == E/'production'
        with (d/'preflight_stdout.json').open('x') as out:
            subprocess.run([sys.executable,str(d/'launch_once.py'),'--config',item['config'],'--check'],
                cwd=str(ROOT),stdout=out,check=True)
    assert (E/'watch_matrix_settings.json').is_file()
    subprocess.run([sys.executable,str(E/'watch_production_matrix.py'),'--matrix',str(matrix_path),
        '--settings',str(E/'watch_matrix_settings.json'),'--state-dir',str(E/'aggregate_watch_001'),
        '--validate-settings'],cwd=str(ROOT),check=True)
    with (E/'matrix_launch_attempt.json').open('x') as out:
        json.dump(dict(startedAtUTC=now(),matrixSha256=hashlib.sha256(matrix_path.read_bytes()).hexdigest(),
                       plannedCases=7,plannedJobs=140,automaticResubmit=False),out,indent=2)
        out.write('\n')
        out.flush()
        os.fsync(out.fileno())
    state = dict(phase='submitting',startedAtUTC=now(),cases=[],automaticResubmit=False)
    save('matrix_launch_state.json',state)
    for item in matrix['cases']:
        d = Path(item['launchDirectory'])
        with (d/'launcher_stdout.log').open('x') as out:
            result = subprocess.run([sys.executable,'-u',str(d/'launch_once.py'),'--config',item['config'],'--submit'],
                cwd=str(ROOT),stdout=out,stderr=subprocess.STDOUT)
        if result.returncode:
            state.update(phase='attention_required',caseId=item['caseId'],returnCode=result.returncode)
            save('matrix_launch_state.json',state)
            raise RuntimeError('Submission failed/ambiguous: '+item['caseId']+'; no retry')
        case_state = json.loads((d/'launch_state.json').read_text())
        assert case_state['phase'] == 'submitted_verified' and case_state['actualJobs'] == 20
        state['cases'].append(dict(caseId=item['caseId'],jobid=case_state['jobid'],
                                   runmeta=case_state['runmeta'],eventLog=case_state['eventLog']))
        save('matrix_launch_state.json',state)
        print(json.dumps(state['cases'][-1]),flush=True)
    state.update(phase='watching',actualJobs=140,submittedAtUTC=now())
    save('matrix_launch_state.json',state)
    command = [sys.executable,'-u',str(E/'watch_production_matrix.py'),'--matrix',str(matrix_path),
        '--settings',str(E/'watch_matrix_settings.json'),'--state-dir',str(E/'aggregate_watch_001')]
    with (E/'production_matrix_watch.log').open('x') as out:
        result = subprocess.run(command,cwd=str(ROOT),stdout=out,stderr=subprocess.STDOUT)
    state.update(phase='finished' if result.returncode == 0 else 'attention_required',
                 watcherReturnCode=result.returncode,endedAtUTC=now())
    save('matrix_launch_state.json',state)
    return result.returncode

if __name__ == '__main__':
    sys.exit(main())
