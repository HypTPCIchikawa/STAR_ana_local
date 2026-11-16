#!/usr/bin/env python3
"""Read-only numerical cross-check of saved projections and CF products.

No analysis module is imported. Native-bin sums, ratio derivatives, covariance,
and high-k GLS are calculated independently from the exported count arrays.
"""
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import traceback
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'femto_lambda_ep_rapidity_cf_20261115_001'
REPORT = HERE / 'independent_validation.json'


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def fingerprint(path):
    st = path.stat()
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return dict(path=str(path), size=st.st_size, mtimeNs=st.st_mtime_ns,
                inode=st.st_ino, device=st.st_dev, sha256=sha.hexdigest())


def raw(path):
    result = {}
    with path.open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if not row or row[0] != 'BIN':
                continue
            require(len(row) == 10, 'native TSV row width')
            ds, selection, cent, side = row[1:5]
            key = '__'.join((ds, cent, selection))
            bins = result.setdefault(key, {}).setdefault(side, {})
            index = int(row[5])
            require(index not in bins, 'duplicate native bin: ' + key)
            bins[index] = tuple(map(float, row[6:]))
    for key, sides in result.items():
        require(set(sides) == {'SE', 'ME'}, 'missing native side: ' + key)
        for side, bins in sides.items():
            require(set(bins) == set(range(202)), 'missing native/flow bin: ' + key)
            sides[side] = np.array([bins[i] for i in range(202)])
            require(np.isfinite(sides[side]).all(), 'nonfinite native array')
            require((sides[side][:, 2:] >= 0).all(), 'negative native count/variance')
        require(np.array_equal(sides['SE'][:, :2], sides['ME'][:, :2]), 'SE/ME edges differ')
    return result


def window(edges, bounds):
    found = [np.flatnonzero(np.abs(edges - x) < 1.e-9) for x in bounds]
    require(all(len(x) == 1 for x in found), 'window not aligned to display edges')
    return np.arange(int(found[0][0]), int(found[1][0]))


def ratio(s, m, vs, vm, target, norm):
    sw, mw = float(s[target].sum()), float(m[target].sum())
    sn, mn = float(s[norm].sum()), float(m[norm].sum())
    if min(sw, mw, sn, mn) <= 0:
        return None, None
    value = sw * mn / (mw * sn)
    ds, dm = np.zeros(len(s)), np.zeros(len(s))
    ds[target] += value / sw
    ds[norm] -= value / sn
    dm[norm] += value / mn
    dm[target] -= value / mw
    return value, math.sqrt(float(np.dot(vs, ds * ds) + np.dot(vm, dm * dm)))


def validate(report):
    new = json.loads((HERE / 'analysis_results.json').read_text())
    old = json.loads((OLD / 'analysis_results.json').read_text())
    cfg, curves = new['config'], new['curves']
    native, previous = raw(HERE / 'raw_projections.tsv'), raw(OLD / 'raw_projections.tsv')
    require(cfg['normalization'] == [.4, .5] and cfg['rebin_factor'] == 4, 'fixed comparison conditions')
    canonical = {'c70_80': 0, 'c60_70': 1, 'c50_60': 2, 'c40_50': 3,
                 'c30_40': 4, 'c20_30': 5, 'c10_20': 6, 'c5_10': 7, 'c0_5': 8}
    centcfg = {c['id']: c for c in cfg['centralities']}
    require(set(centcfg) == set(canonical) | {'c0_10', 'c0_20'}, 'eleven centrality slices')
    for name, index in canonical.items():
        require(centcfg[name]['min'] == centcfg[name]['max'] == index, 'canonical cent9: ' + name)
    require((centcfg['c0_10']['min'], centcfg['c0_10']['max']) == (7, 8), '0-10 definition')
    require((centcfg['c0_20']['min'], centcfg['c0_20']['max']) == (6, 8), '0-20 definition')
    require(len(curves) == len(native) == 484, '484 CF selections / 968 native histograms')
    selections = sorted({r['selection'] for r in curves.values()})
    require(len(selections) == 22, '22 selections')
    aggregates, totals = 0, []
    for ds in cfg['inputs']:
        for sel in selections:
            get = lambda cent, side: native['__'.join((ds, cent, sel))][side]
            for group, parts in [('c0_10', ['c0_5', 'c5_10']),
                                 ('c0_20', ['c0_5', 'c5_10', 'c10_20'])]:
                for side in ('SE', 'ME'):
                    arrays = [get(c, side) for c in parts]
                    reference = get(group, side)
                    require(all(np.array_equal(x[:, :2], reference[:, :2]) for x in arrays), 'aggregate axis agreement')
                    require(np.array_equal(sum(x[:, 2:] for x in arrays), reference[:, 2:]),
                            'exact centrality count/variance closure including flows: ' + str((ds, sel, group, side)))
                    aggregates += 1
            total = dict(dataset=ds, selection=sel)
            for side in ('SE', 'ME'):
                arrays = [get(c, side) for c in canonical]
                require(all(np.array_equal(x[:, :2], arrays[0][:, :2]) for x in arrays), 'nine-single axes agree')
                combined = sum(x[:, 2:] for x in arrays)
                total[side] = dict(regularCount=float(combined[1:-1, 0].sum()),
                                  regularVariance=float(combined[1:-1, 1].sum()),
                                  flowCount=combined[[0, -1], 0].tolist(),
                                  flowVariance=combined[[0, -1], 1].tolist(),
                                  allBinsCount=float(combined[:, 0].sum()),
                                  allBinsVariance=float(combined[:, 1].sum()))
            totals.append(total)
    report['centrality'] = dict(canonicalSingleBins=canonical, exactAggregateSideComparisons=aggregates,
        binsPerSide=202, nineSingleCoverage='PASS: exactly cent9 0..8 for each dataset/selection',
        allNineTotals=totals, independent0to80Reference='NOT_RUN: no separate 0-80 projection was produced')
    for key, sides in previous.items():
        require(key in native, 'old raw selection retained')
        for side, values in sides.items():
            require(np.array_equal(values, native[key][side]), 'old native content/variance/edges exact: ' + key + side)
    for key, result in old['curves'].items():
        require(curves[key] == result, 'old CF series including points/errors/windows/fits exact: ' + key)
    report['previousStudyComparison'] = dict(rawHistograms=len(previous) * 2, rawNativeBinsPerHistogram=202,
        rawEdgesCountsVariances='EXACT_EQUAL', cfSeries=len(old['curves']), cfFullResultObjects='EXACT_EQUAL')

    max_value = max_error = max_cov = max_fit = max_fit_error = max_chi2 = 0.
    points_checked = windows_checked = fits_checked = empty_points = 0
    saved_cov = np.load(str(HERE / 'cf_covariances.npz'))
    require(set(saved_cov.files) == set(curves), 'all covariance matrices present')
    for key, r in curves.items():
        sides = native[key]
        edges = np.r_[sides['SE'][1:-1, 0], sides['SE'][-2, 1]][::4]
        data = {side: sides[side][1:-1, 2:].reshape(50, 4, 2).sum(axis=1) for side in ('SE', 'ME')}
        s, vs = data['SE'].T
        m, vm = data['ME'].T
        require(np.array_equal(edges, r['edges']), 'rebinned edges exact')
        for name, values in [('SE', s), ('ME', m), ('SEvariance', vs), ('MEvariance', vm)]:
            require(np.array_equal(values, r[name]), 'native->display count/variance exact: ' + key + name)
        for side in ('SE', 'ME'):
            require(r['kstarFlow' + side] == sides[side][[0, -1], 2].tolist(), 'flow metadata exact')
        norm = window(edges, cfg['normalization'])
        require(float(s[norm].sum()) == r['normSE'] and float(m[norm].sum()) == r['normME'], 'normalization counts')

        def compare(target, reference, output):
            nonlocal max_value, max_error
            value, error = ratio(s, m, vs, vm, target, reference)
            require((value is None) == (output['value'] is None) and
                    (error is None) == (output['error'] is None), 'NA mask agrees with positive counts')
            if value is not None:
                dv, de = abs(value - output['value']), abs(error - output['error'])
                max_value, max_error = max(max_value, dv), max(max_error, de)
                require(dv <= 1.e-12 * max(1., abs(value)) and de <= 1.e-12 * max(1., error), 'independent ratio/error')

        for i, point in enumerate(r['points']):
            require(point['low'] == edges[i] and point['high'] == edges[i + 1], 'point edges')
            compare(np.array([i]), norm, point)
            points_checked += 1
            empty_points += point['value'] is None
        for name, bounds in cfg['high_windows'].items():
            target = window(edges, bounds)
            compare(target, norm, r['windows'][name])
            require(r['windows'][name]['SE'] == float(s[target].sum()) and
                    r['windows'][name]['ME'] == float(m[target].sum()), 'window count integral')
            windows_checked += 1
        compare(window(edges, [.8, 1.]), window(edges, [.5, .7]), r['shape_ratio'])
        if r['normSE'] > 0 and r['normME'] > 0:
            value, error = ratio(s, m, vs, vm, norm, norm)
            require(value == 1. and error == 0., 'normalization self-ratio closure')
        ids = [i for i, x in enumerate(r['points']) if x['value'] is not None]
        require(ids == r['valid_covariance_bins'], 'covariance mask/order')
        if not ids:
            require(saved_cov[key].shape == (0, 0), 'empty covariance')
            continue
        values = np.array([r['points'][i]['value'] for i in ids])
        sn, mn = r['normSE'], r['normME']
        beta = vs[norm].sum() / sn ** 2 + vm[norm].sum() / mn ** 2
        overlap = np.array([float(i in norm) * (vs[i] / (s[i] * sn) + vm[i] / (m[i] * mn)) for i in ids])
        diagonal = np.array([values[j] ** 2 * (vs[i] / s[i] ** 2 + vm[i] / m[i] ** 2) for j, i in enumerate(ids)])
        cov = np.diag(diagonal) + np.outer(values, values) * (beta - overlap[:, None] - overlap[None, :])
        max_cov = max(max_cov, float(abs(cov - saved_cov[key]).max()))
        require(np.allclose(cov, saved_cov[key], rtol=1.e-11, atol=1.e-11), 'independent shared-normalization covariance')
        require(np.linalg.eigvalsh(cov).min() >= -1.e-10 * max(1., float(np.diag(cov).max())), 'covariance PSD')
        for name, bounds in cfg['fit_ranges'].items():
            fit = r['fits'][name]
            ii = [i for i in ids if edges[i] >= bounds[0] - 1.e-9 and edges[i + 1] <= bounds[1] + 1.e-9]
            if len(ii) < 3:
                require(fit['slope'] is None, 'too few fit bins')
                continue
            require(not set(ii) & set(norm), 'Sherman-Morrison fit uses disjoint normalization')
            y = np.array([r['points'][i]['value'] for i in ii])
            x = np.array([(edges[i] + edges[i + 1]) / 2 for i in ii]); x -= x.mean()
            d = np.array([y[j] ** 2 * (vs[i] / s[i] ** 2 + vm[i] / m[i] ** 2) for j, i in enumerate(ii)])
            require((d > 0).all(), 'positive independent diagonal variance')
            invd = 1 / d; v = invd * y
            inv = np.diag(invd) - beta * np.outer(v, v) / (1 + beta * np.dot(y, v))
            design = np.column_stack((np.ones(len(x)), x)); normal = design.T.dot(inv).dot(design)
            pcov = np.linalg.inv(normal); pars = pcov.dot(design.T).dot(inv).dot(y)
            require(fit['slope'] is not None and fit['rank'] == len(ii), 'full-rank high-k fit')
            max_fit = max(max_fit, abs(float(pars[0]) - fit['intercept']), abs(float(pars[1]) - fit['slope']))
            max_fit_error = max(max_fit_error, abs(math.sqrt(float(pcov[1, 1])) - fit['slope_error']))
            residual = y - design.dot(pars)
            max_chi2 = max(max_chi2, abs(float(residual.dot(inv).dot(residual)) - fit['linear_chi2']),
                           abs(float((y - 1).dot(inv).dot(y - 1)) - fit['unity_chi2']))
            fits_checked += 1
    require(max_fit < 1.e-9 and max_fit_error < 1.e-9 and max_chi2 < 1.e-7, 'independent GLS tolerance')
    report['numerical'] = dict(series=len(curves), pointRatiosAndErrors=points_checked, pointsMarkedNA=empty_points,
        integratedWindowRatios=windows_checked, shapeRatios=len(curves), sharedCovarianceMatrices=len(curves),
        independentGLSFits=fits_checked, nativeRebin='EXACT_EQUAL', maxAbsRatioDifference=max_value,
        maxAbsErrorDifference=max_error, maxAbsCovarianceDifference=max_cov, maxAbsFitParameterDifference=max_fit,
        maxAbsSlopeErrorDifference=max_fit_error, maxAbsChi2Difference=max_chi2)
    closure = json.loads((HERE / 'closure_checks.json').read_text())
    require(all(x['passed_regular'] for x in closure), 'producer regular-bin closure records')
    report['reportedProjectionClosure'] = dict(checks=len(closure), allRegularPassed=True,
        maxRegularCountDifference=max(x['max_regular_count_difference'] for x in closure),
        nonzeroFlowChecks=[x for x in closure if any(x['kstar_flow_count_difference'])])
    with (HERE / 'metrics.csv').open() as stream:
        metrics = list(csv.DictReader(stream))
    require(len(metrics) == len(curves), 'CSV metric rows')
    for row in metrics:
        r = curves['__'.join((row['dataset'], row['centrality'], row['selection']))]
        expected = dict(SEtotal=r['totalSE'], MEtotal=r['totalME'], SEnorm=r['normSE'], MEnorm=r['normME'],
            CFtail=r['windows']['k080_100']['value'], CFtailError=r['windows']['k080_100']['error'],
            shapeRatio=r['shape_ratio']['value'], shapeRatioError=r['shape_ratio']['error'],
            slope=r['fits']['high'].get('slope'), slopeError=r['fits']['high'].get('slope_error'),
            unityChi2=r['fits']['high'].get('unity_chi2'), unityNdf=r['fits']['high'].get('unity_ndf'))
        for name, value in expected.items():
            require((None if row[name] == '' else float(row[name])) == value, 'CSV/JSON metric agreement')
    report['metricsRowsCompared'] = len(metrics)
    report['limitations'] = [
        'No independent 0-80 projection exists; nine-single sums are recorded, not compared to an external reference.',
        'Count model includes common normalization covariance but omits event/pair reuse and EP1/EP8 cross-covariance.',
        'Zero-count points are omitted; sparse-window Gaussian errors/GLS are diagnostic only.',
        'Input ROOT files and ROOT drawing objects are not reopened by this numerical validation.'
    ]


def main():
    paths = [HERE / x for x in ['raw_projections.tsv', 'analysis_results.json', 'cf_covariances.npz',
             'metrics.csv', 'closure_checks.json', 'cf_histograms.root']]
    paths += [OLD / 'raw_projections.tsv', OLD / 'analysis_results.json', Path(__file__).resolve()]
    report = dict(startedAt=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  inputBefore=[fingerprint(p) for p in paths], status='RUNNING')
    try:
        validate(report)
        report['status'] = 'PASS'
    except Exception:
        report['status'] = 'FAIL'
        report['exception'] = traceback.format_exc()
    report['inputAfter'] = [fingerprint(p) for p in paths]
    report['inputsUnchanged'] = report['inputBefore'] == report['inputAfter']
    if not report['inputsUnchanged']:
        report['status'] = 'FAIL'
    report['completedAt'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: report[k] for k in ['status', 'inputsUnchanged']}, ensure_ascii=False))
    if 'numerical' in report:
        print(json.dumps(report['numerical'], ensure_ascii=False))
    if 'exception' in report:
        print(report['exception'])
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
