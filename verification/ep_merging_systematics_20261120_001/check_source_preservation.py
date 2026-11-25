#!/usr/bin/env python3
"""Check the frozen initial 38-file source snapshot without modifying source."""
import datetime
import hashlib
import json
from pathlib import Path

study = Path(__file__).resolve().parent
project = study.parents[1]
manifest_path = study / 'source_before_manifest.json'
manifest = json.loads(manifest_path.read_text())
allowed = {
    'Makefile', 'include/FemtoEpdEventPlane.h',
    'StMaker/StFemtoMaker/StFemtoMaker.cxx', 'StMaker/StFemtoMaker/StFemtoMaker.h',
    'StMaker/common/FemtoLambdaLegacy.cxx', 'StMaker/common/FemtoLambdaLegacy.h',
    'StMaker/StFemtoMaker/README.md',
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

allowed_changed, unchanged, unexpected = [], [], []
for row in manifest['files']:
    name = row['path']
    current = project / name
    saved = Path(row['saved'])
    item = dict(path=name, beforeSha256=row['sha256'],
                currentSha256=sha(current) if current.is_file() else None,
                snapshotSha256=sha(saved) if saved.is_file() else None,
                originalTrackedChange=row.get('originalTrackedChange'))
    if item['snapshotSha256'] != row['sha256']:
        item['reason'] = 'initial snapshot missing or altered'
        unexpected.append(item)
    elif item['currentSha256'] == row['sha256']:
        unchanged.append(item)
    elif current.is_file() and name in allowed:
        allowed_changed.append(item)
    else:
        item['reason'] = 'existing source missing or changed outside allowlist'
        unexpected.append(item)

cases_path = study / 'cases.json'
cases = json.loads(cases_path.read_text())
expected_new = []
for row in cases.get('generatedConfigFiles', []):
    path = Path(row['path'])
    expected_new.append(dict(path=str(path.relative_to(project)), exists=path.is_file(),
                             sha256=sha(path) if path.is_file() else None))
new_geometry = project/'include/FemtoLambdaPairGeometry.h'
expected_new.append(dict(path=str(new_geometry.relative_to(project)), exists=new_geometry.is_file(),
                         sha256=sha(new_geometry) if new_geometry.is_file() else None))
for pattern in ('*pair_geometry*.C', '*ep*partition*.C', '*ep*mix*bin*.C'):
    for path in sorted((project/'tests').glob(pattern)):
        item = dict(path=str(path.relative_to(project)), exists=True, sha256=sha(path))
        if item not in expected_new:
            expected_new.append(item)

result = dict(passed=not unexpected, checkedAtUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              sourceManifest=str(manifest_path), sourceManifestSha256=sha(manifest_path),
              baselineFilesChecked=len(manifest['files']), allowedChanged=allowed_changed,
              unchanged=unchanged, unexpected=unexpected, expectedNewFilesObserved=expected_new,
              scope='All existing files in the initial source_before_manifest are compared by SHA256. '
                    'Every preserved source snapshot is also verified. New case configurations and '
                    'pair-geometry/test additions are listed separately; source files are never edited.')
destination = study/'source_preservation_check.json'
with destination.open('x') as output:
    json.dump(result, output, ensure_ascii=False, indent=2)
    output.write('\n')
print(json.dumps(dict(passed=result['passed'], baselineFilesChecked=len(manifest['files']),
                     allowedChanged=[r['path'] for r in allowed_changed], unchanged=len(unchanged),
                     unexpected=unexpected, output=str(destination)), indent=2))
raise SystemExit(0 if result['passed'] else 1)
