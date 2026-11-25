#!/usr/bin/env python3
"""Read-only verification against campaign snapshots; writes only its JSON report."""
import datetime
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
D = Path(__file__).resolve().parent
REPORT = D / 'source_preservation.json'
ALLOWED = {
    'StMaker/common/FemtoLambdaLegacy.cxx',
    'StMaker/common/FemtoLambdaLegacy.h',
    'StMaker/StFemtoMaker/README.md',
}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def main():
    require(not REPORT.exists(), 'Refusing to overwrite existing source preservation report')
    campaign = json.loads((D / 'analysis_campaign.json').read_text())
    snapshot = json.loads((D / 'source_before_manifest.json').read_text())
    status = json.loads((D / 'source_implementation_status.json').read_text())
    require(len(snapshot['files']) == 24, 'Expected the complete 24-file source snapshot')
    require(set(status['publicFilesChanged']) == ALLOWED, 'Implementation source list differs')
    files = []
    changed = set()
    for record in snapshot['files']:
        path = record['path']
        backup = D / 'source_before' / path
        require(backup.is_file(), 'Missing snapshot backup: ' + path)
        require(backup.stat().st_size == record['bytes'], 'Snapshot size mismatch: ' + path)
        require(sha(backup) == record['sha256'], 'Snapshot fingerprint mismatch: ' + path)
        current = REPO / path
        actual = sha(current)
        is_changed = actual != record['sha256']
        if is_changed:
            changed.add(path)
        require(not is_changed or path in ALLOWED, 'Unrelated source changed: ' + path)
        if path in ALLOWED:
            require(actual == status['sourceSha256'][path], 'Implementation fingerprint differs: ' + path)
        files.append({'path': path, 'snapshotSha256': record['sha256'], 'currentSha256': actual,
                      'changed': is_changed, 'backupVerified': True})
    require(changed == ALLOWED, 'Expected exactly the three implementation source changes')
    configs = []
    for label in ['originalReferenceConfigSha256', 'generatedFiles']:
        for path, expected in campaign[label].items():
            actual = sha(path)
            require(actual == expected, label + ' fingerprint mismatch: ' + path)
            configs.append({'category': label, 'path': path, 'sha256': actual})
    provenance = []
    for name in ['sourceManifest', 'inputList']:
        path = campaign[name]
        actual = sha(path)
        require(actual == campaign[name + 'Sha256'], 'Input provenance mismatch: ' + path)
        provenance.append({'kind': name, 'path': path, 'sha256': actual})
    inputs = []
    for record in campaign['sources']:
        p = Path(record['path'])
        st = p.stat()
        actual = {'bytes': st.st_size, 'inode': st.st_ino, 'device': st.st_dev, 'mtimeNs': st.st_mtime_ns}
        require(all(actual[k] == record[k] for k in actual), 'Rich period ROOT changed: ' + str(p))
        inputs.append({'path': str(p), **actual})
    require(len(inputs) == campaign['expectedPeriods'] == 20, 'Expected 20 unchanged rich period ROOTs')
    before_legacy = (D / 'source_before/StMaker/common/FemtoLambdaLegacy.cxx').read_text()
    require('FemtoLambdaPairEpdBinCount' in before_legacy, 'Extra EP metadata must predate this campaign')
    result = {
        'passed': True,
        'verifiedUTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'snapshotManifest': str(D / 'source_before_manifest.json'),
        'snapshotManifestSha256': sha(D / 'source_before_manifest.json'),
        'snapshotFiles': 24,
        'allowedChangedFiles': sorted(changed),
        'unrelatedSnapshotFilesUnchanged': 21,
        'makefileUnchanged': 'Makefile' not in changed,
        'files': files,
        'originalReferenceConfigsUnchanged': len(campaign['originalReferenceConfigSha256']),
        'generatedConfigsUnchanged': len(campaign['generatedFiles']),
        'configFingerprints': configs,
        'inputProvenance': provenance,
        'richPeriodROOTsUnchangedByStatFingerprint': inputs,
        'preexistingEPBinMetadataConfirmedInSourceSnapshot': True,
        'sourceImplementationStatus': str(D / 'source_implementation_status.json'),
        'sourceImplementationStatusHistoricalBuildFlags': 'Recorded before parent build and fixture; not used as execution-result authority.',
        'noSourceLibraryConfigurationOrInputMutationsByThisAudit': True,
    }
    REPORT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print('SOURCE_PRESERVATION_PASS', REPORT)

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        if not REPORT.exists():
            REPORT.write_text(json.dumps({'passed': False, 'error': str(exc)}, ensure_ascii=False, indent=2) + '\n')
        raise
