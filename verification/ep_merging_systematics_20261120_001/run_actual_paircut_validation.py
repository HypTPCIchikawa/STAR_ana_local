#!/usr/bin/env python3
"""Run one small production-helper test, preserve logs and gate on real PASS."""
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

project = Path.cwd()
study = Path('rootfile/femto_lambda_ep_merging_systematics_20261120_001')
mainconf = sys.argv[1] if len(sys.argv) > 1 else 'config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_ep8_mergingOn_mix100.yaml'
result_path = study / 'test_actual_lambda_merging_paircut_result.json'
if result_path.exists():
    raise SystemExit('Refusing to overwrite existing actual pair-cut validation result')

sources = [mainconf, str(study/'test_actual_lambda_merging_paircut.C'), str(study/'run_actual_lambda_merging_paircut.C'),
           'include/FemtoLambdaPairGeometry.h', 'StMaker/common/FemtoLambdaLegacy.h',
           'StMaker/common/FemtoLambdaLegacy.cxx', 'StMaker/StFemtoMaker/StFemtoMaker.h',
           'StMaker/StFemtoMaker/StFemtoMaker.cxx', 'script/singularity_root_batch.sh',
           'lib/libStCommon.so', 'lib/libStarAnaConfig.so', 'lib/libStRefMultCorr.so']

def fingerprint(name):
    path = Path(name)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return {'path': name, 'size': path.stat().st_size, 'sha256': digest.hexdigest()}

before = [fingerprint(name) for name in sources]
started = datetime.datetime.now(datetime.timezone.utc)
stamp = started.strftime('%Y%m%dT%H%M%S%fZ')
log = study / ('test_actual_lambda_merging_paircut_' + stamp + '.log')
macro = str(study/'run_actual_lambda_merging_paircut.C') + '(' + json.dumps(mainconf) + ')'
command = ['bash', 'script/singularity_root_batch.sh', mainconf, 'root4star', '-l', '-b', '-q', macro]
with log.open('x') as handle:
    proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
text = log.read_text(errors='replace')
print(text, end='')
matched = re.search(r'PASS actual Lambda merging pair cut: checks=(\d+) He4 calls=(\d+)', text)
if proc.returncode or not matched or 'FAIL actual Lambda merging pair cut:' in text:
    print('Actual pair-cut validation FAILED; result gate not written; log=' + str(log), file=sys.stderr)
    raise SystemExit(proc.returncode or 1)
after = [fingerprint(name) for name in sources]
if before != after:
    raise SystemExit('Source/library fingerprints changed during actual pair-cut validation; gate not written')
toolchain = [line for line in text.splitlines() if line.startswith('ACTUAL_PAIR_CUT_ROOT ') or line.startswith('g++ (')]
result = {'passed': True, 'returnCode': proc.returncode, 'mainconf': mainconf,
          'checks': int(matched.group(1)), 'he4Calls': int(matched.group(2)),
          'log': str(log), 'logSha256': fingerprint(str(log))['sha256'],
          'toolchain': toolchain, 'command': command,
          'startedUtc': started.isoformat(), 'finishedUtc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'inputFingerprintsUnchanged': True, 'sourceAndLibraryFingerprints': before,
          'coverage': ['SE and both ME directions with nonempty rejection', 'strict +/-0.05 boundaries and adjacent doubles',
                       'wrapped phi', 'pion-only close pair kept', 'unreachable proton/nucleus kept',
                       'own-event magnetic field and owned cache', 'different-event identical track index',
                       'wrong-event/missing-track errors', 'nonfinite-point keep policy',
                       'before/after QA and physical CF counters', 'legacy raw aliases before rejection',
                       'He3 and He3AsHe4 invariance'],
          'scope': 'Small explicit synthetic inputs against the actual compiled Legacy helper; no data event loop or batch job'}
with result_path.open('x') as handle:
    json.dump(result, handle, ensure_ascii=False, indent=2)
    handle.write('\n')
print('Saved PASS gate: ' + str(result_path))
