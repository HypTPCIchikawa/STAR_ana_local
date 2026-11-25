#!/usr/bin/env python3
"""Run fixed-prefix real-tree controls without changing production inputs."""
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def save(path, value):
    temporary = Path(str(path) + '.next')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(str(temporary), str(path))


def main():
    if len(sys.argv) != 4:
        raise RuntimeError('Arguments: PROJECT STUDY NUMBER_OF_STORED_EVENTS')
    root, study, count = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    assert count > 357247
    manifest = json.loads((study / 'input_manifest.json').read_text())
    first = manifest['sourceFiles'][0]
    source = Path(first['path'])
    state_path = study / 'pilot_status.json'
    assert not state_path.exists()
    start_stat = source.stat()
    assert (start_stat.st_size, start_stat.st_ino, start_stat.st_dev, start_stat.st_mtime_ns) == (
        first['bytes'], first['inode'], first['device'], first['mtimeNs'])
    labels = [('B_period100', 'periodMix100'), ('A_reset100', 'periodMix100_resetCheck'),
              ('C_period20', 'periodMix20_check')]
    state = dict(startedJst=now(), storedEventsPerPilot=count, input=str(source), pilots={}, rssSamples=[])
    processes = {}
    for label, suffix in labels:
        output = study / ('pilot_' + label + '.root')
        assert not output.exists()
        conf = 'config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_' + suffix + '.yaml'
        command = ['/usr/bin/time', '-v', '-o', str(study / ('pilot_' + label + '.time')),
                   'bash', 'script/singularity_run_anaFemtoLambda_4He_tree.sh', conf,
                   str(source), str(output), 'periodMix100Check_' + label, str(count)]
        log = (study / ('pilot_' + label + '.log')).open('x')
        process = subprocess.Popen(command, cwd=str(root), stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        log.close()
        processes[label] = process
        state['pilots'][label] = dict(command=command, pid=process.pid, output=str(output))
    save(state_path, state)
    while any(process.poll() is None for process in processes.values()):
        readings = subprocess.check_output(['ps', '-eo', 'pgid=,rss=,comm='], universal_newlines=True)
        sample = dict(checkedJst=now(), groups={label: 0 for label in processes})
        for line in readings.splitlines():
            fields = line.split(None, 2)
            if len(fields) < 3:
                continue
            for label, process in processes.items():
                if int(fields[0]) == process.pid:
                    sample['groups'][label] += int(fields[1])
        state['rssSamples'].append(sample)
        for label, process in processes.items():
            if process.poll() is not None:
                state['pilots'][label]['returnCode'] = process.returncode
        save(state_path, state)
        time.sleep(15)
    final_stat = source.stat()
    assert (start_stat.st_size, start_stat.st_ino, start_stat.st_dev, start_stat.st_mtime_ns) == (
        final_stat.st_size, final_stat.st_ino, final_stat.st_dev, final_stat.st_mtime_ns)
    state['sourceUnchanged'] = True
    for label, process in processes.items():
        record = state['pilots'][label]
        record['returnCode'] = process.returncode
        lines = (study / ('pilot_' + label + '.log')).read_text().splitlines()
        record['completionMarkers'] = [line for line in lines if line.startswith('FemtoLambda tree replay job=')]
        record['passed'] = (process.returncode == 0 and len(record['completionMarkers']) == 1 and
                            'completed=' + str(count) + ' ' in record['completionMarkers'][0] and
                            'status=0 ' in record['completionMarkers'][0])
        output = Path(record['output'])
        record['outputBytes'] = output.stat().st_size if output.exists() else None
    state['finishedJst'] = now()
    state['passed'] = all(record['passed'] for record in state['pilots'].values())
    save(state_path, state)
    print(json.dumps({key: value for key, value in state.items() if key != 'rssSamples'}), flush=True)
    return 0 if state['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
