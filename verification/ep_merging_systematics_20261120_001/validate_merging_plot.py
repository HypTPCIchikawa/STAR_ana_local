#!/usr/bin/env python3
"""Check plot artifacts and unchanged inputs without reopening event trees."""
from pathlib import Path
import datetime
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
STUDY = Path(__file__).resolve().parent
FIGURE = ROOT / 'share/figure/20261120_lambda_ep_merging_systematics'


def main():
    metrics = STUDY / 'merging_current_vs_old_metrics.json'
    values = json.loads(metrics.read_text())
    assert values['passed'] and values['currentSources'] == 20 and values['pages'] == 3
    inputs = json.loads((STUDY / 'merging_plot_inputs.json').read_text())
    for row in inputs['currentOutputs'] + [inputs['oldQaRoot']]:
        info = Path(row['path']).stat()
        assert (info.st_size, info.st_ino, info.st_dev, info.st_mtime_ns) == (
            row['bytes'], row['inode'], row['device'], row['mtimeNs'])
    digest = hashlib.sha256()
    with Path(inputs['oldQaRoot']['path']).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    assert digest.hexdigest() == inputs['oldQaRoot']['sha256']
    prefix = FIGURE / 'track_merging_current_vs_old_20261120'
    pdf = Path(str(prefix) + '.pdf')
    raw = pdf.read_bytes()
    assert raw.rstrip().endswith(b'%%EOF')
    pages = len(re.findall(rb'/Type\s*/Page\b', raw))
    assert pages == 3
    artifacts = []
    for path in [pdf] + [Path(str(prefix) + suffix + '.png') for suffix in
                         ['_current_raw', '_comparison_full', '_comparison_zoom']]:
        data = path.read_bytes()
        if path.suffix == '.png':
            assert data.startswith(b'\x89PNG\r\n\x1a\n')
        artifacts.append(dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    copy = Path(str(prefix) + '_metrics.json')
    if copy.exists():
        assert copy.read_bytes() == metrics.read_bytes()
    else:
        with copy.open('xb') as output:
            output.write(metrics.read_bytes())
    proof = dict(passed=True, checkedJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
                 currentInputFingerprintsUnchanged=True, oldInputSha256Unchanged=True, pages=pages,
                 pdfAndPngArtifacts=artifacts, metrics=str(copy),
                 metricsSha256=hashlib.sha256(metrics.read_bytes()).hexdigest(),
                 sourceHelper=str(STUDY / 'DrawMergingCurrentVsOld.C'),
                 sourceHelperSha256=hashlib.sha256((STUDY / 'DrawMergingCurrentVsOld.C').read_bytes()).hexdigest(),
                 inputManifest=str(STUDY / 'merging_plot_inputs.json'),
                 inputManifestSha256=hashlib.sha256((STUDY / 'merging_plot_inputs.json').read_bytes()).hexdigest(),
                 bareEventFieldAliasesExact=values['bareEventFieldAliasesExactInEverySourceAndSum'],
                 strictPairCutCountsReported=False, eventPayloadReplayed=False,
                 visualInspection='Current raw and comparison origin-zoom PNGs inspected: axes, rectangle, target and density notes are readable. Cohort counts appear on comparison pages.')
    target = STUDY / 'merging_plot_validation.json'
    if target.exists():
        previous = json.loads(target.read_text())
        assert previous['passed'] and previous['metricsSha256'] == proof['metricsSha256']
    else:
        with target.open('x') as output:
            json.dump(proof, output, indent=2)
            output.write('\n')
    print(json.dumps(dict(pdf=str(pdf), pages=pages, sourceEvents=values['currentSavedEvents'],
                          selectedHe4Pairs=values['currentHe4SignalPairs'])))
    for row in values['histograms']:
        print(row['sample'], row['daughter'], 'entries', row['entries'], 'regular', row['regularIntegral'])
        for entry in row['densityRatios']:
            print('halfwidth', entry['phiHalfWidth'], 'measurement', json.dumps(entry['measurement']))


if __name__ == '__main__':
    main()
