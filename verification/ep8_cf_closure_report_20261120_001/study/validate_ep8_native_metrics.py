#!/usr/bin/env python3
"""Independent arithmetic from frozen exported histogram bins; no ROOT/event IO."""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

destination = Path(__file__).resolve().parent
project = destination.parents[1]
source = project/'rootfile/femto_lambda_ep_merging_systematics_20261120_001/cf_baseline_ep8Off_002'
raw_path, cf_path = source/'raw_histograms.json', source/'cf_comparison.json'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

source_hashes = {str(p): digest(p) for p in (raw_path, cf_path)}
raw, reference = json.loads(raw_path.read_text()), json.loads(cf_path.read_text())
checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)

def close(actual, expected, message, tolerance=2.e-12):
    require(math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance), message)

require(raw['passed'] is True and reference['passed'] is True, 'input exports passed')
require(raw['treeEntriesRead'] == reference['treeEntriesRead'] == 0, 'histogram-only exports')
require(len(raw['cases']) == 1, 'one completed baseline case')
case = raw['cases'][0]
require(case['caseId'] == 'ep8_mergingOff' and case['epBins'] == 8 and not case['cutEnabled'], 'EP8 off baseline')
edges = raw['nativeEdges']
require(len(edges) == 201 and edges[0] == 0 and edges[-1] == 1, 'native 200-bin [0,1] axis')
require(all(a < b for a, b in zip(edges, edges[1:])), 'strictly ordered native axis')
size = len(edges) + 1
indexed = {(row['plane'], row['cent9']): row for row in case['slices']}
require(len(indexed) == 81 and len(case['slices']) == 81, 'nine plane labels by nine centrality labels')
fields = ('SE', 'ME', 'SEVariance', 'MEVariance')
for plane in range(-1, 8):
    for cent in range(9):
        row = indexed[plane, cent]
        for name in fields:
            require(len(row[name]) == size, 'explicit native bins including both X flow bins')
            require(all(math.isfinite(x) and x >= 0 for x in row[name]), 'finite nonnegative bins')
        for name in ('SE', 'ME'):
            require(all(x == int(x) for x in row[name]), 'raw pair counts are integers')

# Closure is tested at every native X bin, explicitly including bin 0 and N+1.
closure = []
for cent in range(9):
    for name in fields:
        summed = [sum(indexed[p, cent][name][i] for p in range(8)) for i in range(size)]
        inclusive = indexed[-1, cent][name]
        for actual, expected in zip(summed, inclusive):
            if name in ('SE', 'ME'):
                require(actual == expected, 'eight-plane count sum equals inclusive including flow')
            else:
                close(actual, expected, 'eight-plane variance sum closure', 1.e-10)
        closure.append(dict(cent9=cent, field=name, binsChecked=size,
                            underflow=summed[0], regular=sum(summed[1:-1]), overflow=summed[-1],
                            maxAbsoluteDifference=max(abs(a-b) for a, b in zip(summed, inclusive))))

def edge_index(value):
    index = min(range(len(edges)), key=lambda i: abs(edges[i]-value))
    require(abs(edges[index]-value) < 1.e-12, 'window bound exactly matches a native edge')
    return index

def window(values, low, high):
    first, stop = edge_index(low), edge_index(high)
    require(first < stop, 'nonempty half-open native-bin interval')
    return sum(values[first+1:stop+1])

windows = {'k050_060': [0.5, 0.6], 'k060_080': [0.6, 0.8],
           'k080_100': [0.8, 1.0], 'k060_100': [0.6, 1.0]}
curves = []
for plane in range(-1, 8):
    values = {name: [sum(indexed[plane, c][name][i] for c in range(2, 9))
                     for i in range(size)] for name in fields}
    ref = next(c for c in reference['curves'] if c['centralityId'] == 'c00_60' and c['plane'] == plane)
    require(ref['cent9'] == [2, 8] and ref['normalizationGeVc'] == [0.4, 0.5], 'reference centrality/norm condition')
    norm_se, norm_me = window(values['SE'], .4, .5), window(values['ME'], .4, .5)
    require(norm_se > 0 and norm_me > 0, 'positive independent normalization counts')
    alpha = norm_se/norm_me
    require(norm_se == ref['normSE'] and norm_me == ref['normME'], 'reference normalization raw counts reproduced')
    close(alpha, ref['normAlpha'], 'reference normalization alpha reproduced')
    unity = norm_se/(alpha*norm_me)
    close(unity, 1., 'normalization-window integrated CF equals unity by construction')
    totals = {}
    for name in ('SE', 'ME'):
        totals[name] = dict(underflow=values[name][0], regular=sum(values[name][1:-1]),
                            overflow=values[name][-1], allIncludingFlow=sum(values[name]))
        for label in ('Underflow', 'Regular', 'Overflow'):
            require(totals[name][label.lower()] == ref[name+label], 'reference X range/flow totals reproduced')
    high = {}
    for key, bounds in windows.items():
        se, me = window(values['SE'], *bounds), window(values['ME'], *bounds)
        require(me > 0, 'positive denominator in each high-k window')
        high[key] = dict(windowGeVc=bounds, SE=se, ME=me, alphaScaledME=alpha*me,
                         integralCF=se/(alpha*me))
    original_high = ref['highK']['k060_100']
    require(high['k060_100']['SE'] == original_high['SE'] and
            high['k060_100']['ME'] == original_high['ME'], 'reference high-k raw counts reproduced')
    close(high['k060_100']['integralCF'], original_high['integralCF'], 'reference high-k integral CF reproduced')
    # Rebin only for comparison to the existing figure; high-k metrics above use native bins.
    for display in ref['bins']:
        se, me = window(values['SE'], display['low'], display['high']), window(values['ME'], display['low'], display['high'])
        require(se == display['SE'] and me == display['ME'], 'every display bin reproduced from native counts')
        if me:
            close(se/(alpha*me), display['CF'], 'every existing display CF reproduced')
    curves.append(dict(plane=plane, centrality='0-60%', cent9=[2, 8], nativeBinWidthGeVc=.005,
                       normalizationGeVc=[.4, .5], normSE=norm_se, normME=norm_me,
                       normAlpha=alpha, normalizationWindowIntegralCF=unity,
                       countsIncludingFlow=totals, highK=high))

summary = {}
for key in windows:
    values = [(c['plane'], c['highK'][key]['integralCF']) for c in curves if c['plane'] >= 0]
    summary[key] = dict(windowGeVc=windows[key], minimum=min(v for _, v in values),
                        maximum=max(v for _, v in values),
                        minimumPlane=min(values, key=lambda p:p[1])[0], maximumPlane=max(values, key=lambda p:p[1])[0],
                        allEightBelowUnity=all(v < 1. for _, v in values),
                        allEightAboveUnity=all(v > 1. for _, v in values))
require(source_hashes == {str(p): digest(p) for p in (raw_path, cf_path)}, 'input JSON fingerprints unchanged')
result = dict(passed=True, schemaVersion=1, createdUtc=datetime.now(timezone.utc).isoformat(), checks=checks,
              sourceSha256=source_hashes, helperSha256=digest(Path(__file__)), treeEntriesRead=0,
              centralityCombination='Sum raw native-bin SE and ME over cent9 2-8, then derive one alpha for each plane.',
              normalization='alpha=SE([.4,.5))/ME([.4,.5)); integralCF(W)=SE(W)/(alpha*ME(W)).',
              intervalConvention='Half-open intervals; all bounds align with native X edges. X under/overflow are included in closure and listed separately, not folded into high-k windows.',
              closure=closure, curves=curves, planeRanges=summary,
              limitations=['Descriptive ratios only; no significance test or fitted convergence claim.',
                           'Normalization-window unity is imposed by construction.',
                           'All windows for one curve share its estimated alpha, so their values are correlated.',
                           'Pairs reuse events and tracks in SE/ME; event/shared-ME and normalization covariance are not evaluated.',
                           'Native k* axis ends at 1 GeV/c. Overflow preserves total counts but cannot recover the shape above 1 GeV/c.',
                           'Raw Lambda-He4 signal before background subtraction/purity correction; this independent arithmetic does not establish original direct-Pico equivalence.'])
output = destination/'independent_metrics_validation.json'
with output.open('x') as handle:
    json.dump(result, handle, ensure_ascii=False, indent=2)
    handle.write('\n')
print(json.dumps(dict(passed=True, checks=checks, output=str(output), planeRanges=summary,
                     inclusive=curves[0]['highK']), indent=2))
