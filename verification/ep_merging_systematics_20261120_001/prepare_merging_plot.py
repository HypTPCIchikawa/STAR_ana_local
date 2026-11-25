#!/usr/bin/env python3
"""Pin existing histograms and cohort counts for the merging-map comparison."""
from pathlib import Path
import datetime
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
STUDY = Path(__file__).resolve().parent
OLD = ROOT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/completed.json'
PERIOD = ROOT / 'rootfile/femto_lambda_period_mix100_20261119_001'


def fresh(path, data):
    if path.exists():
        assert path.read_bytes() == data, 'Existing study artifact differs: ' + str(path)
    else:
        with path.open('xb') as output:
            output.write(data)


def main():
    (ROOT / 'share/figure/20261120_lambda_ep_merging_systematics').mkdir(parents=True, exist_ok=True)
    current = json.loads((PERIOD / 'input_manifest.json').read_text())
    snapshot = json.loads((PERIOD / 'output_snapshot_20261120.json').read_text())
    old = json.loads(OLD.read_text())
    old_root = Path(old['sourceRoot'])
    digest = hashlib.sha256()
    with old_root.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    assert digest.hexdigest() == old['sourceRootSha256']
    rows = snapshot['outputs']
    assert len(rows) == 20
    for row in rows:
        info = Path(row['path']).stat()
        assert (info.st_size, info.st_ino, info.st_dev, info.st_mtime_ns) == (
            row['bytes'], row['inode'], row['device'], row['mtimeNs'])
    fresh(STUDY / 'merging_current_outputs.list', ''.join(x['path'] + '\n' for x in rows).encode())
    values = [current['inputEvents'], current['storedEvents'], old['sourceInputEventsRepresented'],
              old['acceptedCacheFrames'], old['selectedHe4Pairs'], old['successfulSourceJobs'], old['excludedSourceJobs']]
    fresh(STUDY / 'merging_plot_metadata.txt', (' '.join(map(str, values)) + '\n').encode())
    info = old_root.stat()
    proof = dict(createdJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
                 currentOutputs=rows,
                 oldQaRoot=dict(path=str(old_root), bytes=info.st_size, mtimeNs=info.st_mtime_ns,
                                inode=info.st_ino, device=info.st_dev, sha256=digest.hexdigest()),
                 oldMetadata=str(OLD), oldMetadataSha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),
                 metadataColumns=['currentInput', 'currentSaved', 'oldInput', 'oldCache', 'oldSelectedHe4Pairs',
                                  'oldSuccessfulJobs', 'oldExcludedJobs'], metadataValues=values,
                 requestedCutTarget='proton-4He only; pion QA/reference, no pion cut',
                 closePairCutAppliedInInputs=False,
                 newSourceMainconf='config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml')
    destination = STUDY / 'merging_plot_inputs.json'
    if destination.exists():
        previous = json.loads(destination.read_text())
        proof['createdJst'] = previous['createdJst']
    fresh(destination, (json.dumps(proof, indent=2) + '\n').encode())
    print(json.dumps(dict(currentSources=len(rows), metadata=values, oldRootSha256=digest.hexdigest())))


if __name__ == '__main__':
    main()
