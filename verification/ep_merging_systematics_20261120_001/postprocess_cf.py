#!/usr/bin/env python3
"""Histogram-only raw-signal CF comparison with explicit case/output coverage."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    with Path(path).open() as stream:
        return json.load(stream)


def save(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def resolve(path, base):
    p = Path(path)
    return (p if p.is_absolute() else base / p).resolve()


def stat_record(path):
    st = path.stat()
    need(path.is_file() and st.st_size > 0, 'Missing/empty regular file: ' + str(path))
    return dict(path=str(path), bytes=st.st_size, device=st.st_dev,
                inode=st.st_ino, mtimeNs=st.st_mtime_ns)


def lines(path):
    result = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        need(line == line.strip() and not any(c.isspace() for c in line),
             'Each list row must be one whitespace-free path: ' + str(path))
        result.append(resolve(line, path.parent))
    need(result and len(result) == len(set(result)), 'Empty/duplicate input list: ' + str(path))
    return result


def root_run(mainconf, macro, arguments, log):
    invocation = str(macro) + '+(' + ','.join(json.dumps(str(x)) for x in arguments) + ')'
    command = ['bash', str(PROJECT / 'script/singularity_root_batch.sh'),
               str(mainconf), 'root4star', '-l', '-b', '-q', invocation]
    with log.open('x') as stream:
        subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT, check=True)


def edge_index(edges, value):
    found = [i for i, edge in enumerate(edges) if math.isclose(edge, value, rel_tol=0, abs_tol=1.e-10)]
    need(len(found) == 1, 'Requested window edge is not a native histogram edge: ' + repr(value))
    return found[0]


def integral(values, edges, interval):
    lo, hi = interval
    if lo < edges[0] - 1.e-10 or hi > edges[-1] + 1.e-10:
        return None
    a, b = edge_index(edges, lo), edge_index(edges, hi)
    need(a < b, 'Reversed/empty interval')
    return sum(values[a:b])


def curve(case, cent, plane, slices, edges, settings):
    native = len(edges) - 1
    summed = {key: [sum(s[key][x] for s in slices) for x in range(native + 2)]
              for key in ['SE', 'ME', 'SEVariance', 'MEVariance']}
    se, me, sv, mv = (summed[key][1:-1] for key in ['SE', 'ME', 'SEVariance', 'MEVariance'])
    norm = settings['normalizationGeVc']
    sh, mh = integral(se, edges, norm), integral(me, edges, norm)
    need(sh is not None and mh is not None, 'Normalization is outside native histogram range')
    alpha = sh / mh if sh > 0 and mh > 0 else None
    result = dict(caseId=case['caseId'], epBins=case['epBins'], cutEnabled=case['cutEnabled'],
                  centralityId=cent['id'], centralityLabel=cent['label'], cent9=cent['cent9'],
                  plane=plane, normalizationGeVc=norm, normSE=sh, normME=mh, normAlpha=alpha,
                  available=alpha is not None, unavailableReason=None if alpha else 'zero SE or ME in normalization window',
                  SERegular=sum(se), MERegular=sum(me), SEUnderflow=summed['SE'][0],
                  SEOverflow=summed['SE'][-1], MEUnderflow=summed['ME'][0], MEOverflow=summed['ME'][-1])
    step = settings['displayRebin']
    need(native % step == 0, 'displayRebin does not exactly divide native bins')
    bins = []
    for first in range(0, native, step):
        s, m, vs, vm = (sum(v[first:first + step]) for v in [se, me, sv, mv])
        value = s / (alpha * m) if alpha is not None and m > 0 else None
        variance = (vs / (alpha * alpha * m * m) + s * s * vm / (alpha * alpha * m ** 4)
                    if value is not None else None)
        bins.append(dict(low=edges[first], high=edges[first + step], SE=s, ME=m,
                         SEVariance=vs, MEVariance=vm, CF=value,
                         fixedAlphaCountError=math.sqrt(max(0, variance)) if variance is not None else None))
    result['bins'] = bins
    result['highK'] = {}
    for name, window in settings['highKWindowsGeVc'].items():
        s, m = integral(se, edges, window), integral(me, edges, window)
        complete = s is not None and m is not None
        points = [b['CF'] for b in bins if b['low'] >= window[0] - 1.e-10
                  and b['high'] <= window[1] + 1.e-10 and b['CF'] is not None]
        result['highK'][name] = dict(windowGeVc=window, available=bool(complete and alpha is not None and m > 0),
                                    unavailableReason=('requested interval outside native axis' if not complete else
                                                       ('zero normalization or ME window' if alpha is None or m <= 0 else None)),
                                    SE=s, ME=m, integralCF=(s / (alpha * m) if complete and alpha is not None and m > 0 else None),
                                    arithmeticMeanOfValidDisplayBins=(sum(points) / len(points) if complete and points else None),
                                    validDisplayBins=len(points) if complete else 0)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cases', type=Path, required=True)
    ap.add_argument('--settings', type=Path, required=True)
    ap.add_argument('--outputs', type=Path, help='Explicit caseId -> outputList manifest; required in compare mode')
    ap.add_argument('--mode', choices=['baseline', 'compare'], required=True)
    ap.add_argument('--output-dir', type=Path, required=True, help='New directory; existing directories are refused')
    args = ap.parse_args()
    case_path, settings_path = args.cases.resolve(), args.settings.resolve()
    cases, settings = load(case_path), load(settings_path)
    need(cases['schemaVersion'] == settings['schemaVersion'] == 1, 'Unsupported schema')
    required = [(n, cut) for n in [1, 8, 16, 32] for cut in [False, True]]
    case_by_id = {c['caseId']: c for c in cases['cases']}
    need(len(case_by_id) == len(cases['cases']) == 8 and
         sorted((c['epBins'], c['cutEnabled']) for c in cases['cases']) == sorted(required),
         'cases.json must define exactly EP1/8/16/32 x cutOff/On')
    baseline = case_by_id[cases['baselineReuseCaseId']]
    need(baseline['epBins'] == 8 and baseline['cutEnabled'] is False and baseline['baselineReuse'] is True,
         'Pinned baseline must be EP8 cutOff')
    need([c['caseId'] for c in cases['cases'] if c['baselineReuse']] == [baseline['caseId']],
         'Only the pinned baseline may use legacy metadata exceptions')
    need(settings['baselineReuseJobid'] == cases['baselineJobid'], 'Baseline job id mismatch')
    selected = [baseline] if args.mode == 'baseline' else sorted(cases['cases'], key=lambda c: (c['epBins'], c['cutEnabled']))
    need(settings['normalizationInterval'] == 'half-open [low, high)', 'Unknown normalization convention')
    need(isinstance(settings['displayRebin'], int) and settings['displayRebin'] > 0, 'Invalid display rebin')
    for interval in [settings['normalizationGeVc'], settings['displayKGeVc'], settings['displayCFRange']] + list(settings['highKWindowsGeVc'].values()):
        need(len(interval) == 2 and all(math.isfinite(x) for x in interval) and interval[0] < interval[1], 'Invalid interval')
    centrals = settings['centralities']
    need(len({c['id'] for c in centrals}) == len(centrals), 'Duplicate centrality id')
    for c in centrals:
        need(0 <= c['cent9'][0] <= c['cent9'][1] <= 8 and all(isinstance(x, int) for x in c['cent9']), 'Invalid cent9 bounds')
        need(all(ch.isalnum() or ch == '_' for ch in c['id']) and '\t' not in c['label'] and '\n' not in c['label'], 'Unsafe centrality text')
    need(set(settings['diagnosticCentralities']) <= {c['id'] for c in centrals}, 'Unknown diagnostic centrality')
    need(settings['requiredValidCentralities'] and set(settings['requiredValidCentralities']) <= {c['id'] for c in centrals},
         'Missing/unknown required main centrality')
    expected = cases['expectedPeriods']
    need(expected == 20, 'This study requires all 20 periods')
    input_manifest_path = resolve(cases['inputManifest'], case_path.parent)
    period_list_path = resolve(cases['inputList'], case_path.parent)
    need(digest(input_manifest_path) == cases['inputManifestSha256'] and digest(period_list_path) == cases['inputListSha256'], 'Input manifest/list hash mismatch')
    manifest = load(input_manifest_path)
    source_periods = manifest['sourceFiles']
    period_paths = lines(period_list_path)
    need(len(source_periods) == len(period_paths) == expected, 'Incomplete period source manifest')
    need([str(p) for p in period_paths] == [str(Path(s['path']).resolve()) for s in source_periods], 'Period manifest order mismatch')
    need(sum(s['acceptedEvents'] for s in source_periods) == cases['storedEvents'] == manifest['storedEvents'], 'Stored event total mismatch')
    pins = [{'path': str(p), 'sha256': digest(p)} for p in [case_path, settings_path, input_manifest_path, period_list_path,
                                                          HERE / 'postprocess_cf.py', HERE / 'ExtractPeriodCF.C', HERE / 'DrawPeriodCF.C']]
    mapping = None
    if args.outputs:
        output_manifest_path = args.outputs.resolve()
        output_manifest = load(output_manifest_path)
        need(output_manifest['schemaVersion'] == 1, 'Unsupported output manifest schema')
        mapping = output_manifest['cases']
        need(not set(mapping) - set(case_by_id), 'Unknown case in output manifest')
        if args.mode == 'compare':
            need(set(mapping) == set(case_by_id), 'Compare mode requires explicit output lists for all eight cases')
        pins.append({'path': str(output_manifest_path), 'sha256': digest(output_manifest_path)})
    need(args.mode != 'compare' or mapping is not None, 'Compare mode requires --outputs')
    plans, stats, identities = [], [], set()
    for c in selected:
        if mapping is not None and c['caseId'] in mapping:
            list_path = resolve(mapping[c['caseId']]['outputList'], args.outputs.resolve().parent)
            outputs = lines(list_path)
            pins.append({'path': str(list_path), 'sha256': digest(list_path)})
        else:
            need(c['baselineReuse'], 'Missing output list for ' + c['caseId'])
            outputs = [Path(p).resolve() for p in c['existingBaselineOutputs']]
        need(len(outputs) == expected and len(set(outputs)) == expected, 'Each case must contain exactly 20 unique outputs')
        if c['baselineReuse']:
            need([str(p) for p in outputs] == c['existingBaselineOutputs'], 'Baseline output list differs from pinned existing outputs/order')
        for i, (path, src) in enumerate(zip(outputs, source_periods)):
            st = stat_record(path)
            identity = (st['device'], st['inode'])
            need(identity not in identities, 'Duplicate/hardlinked ROOT input across cases: ' + str(path))
            identities.add(identity)
            st.update(caseId=c['caseId'], period=i)
            stats.append(st)
            fields = [c['caseId'], str(c['epBins']), str(int(c['cutEnabled'])), str(int(c['baselineReuse'])),
                      str(i), str(src['acceptedEvents']), Path(src['path']).name, str(path)]
            need(not any('\t' in x or '\n' in x for x in fields), 'Unsafe TSV field')
            plans.append(fields)
    out = args.output_dir.resolve()
    need(not out.exists(), 'Refusing existing output directory: ' + str(out))
    out.mkdir(parents=True)
    save(out / 'input_provenance.json', dict(mode=args.mode, cases=[c['caseId'] for c in selected], pins=pins,
                                           sourceStats=stats, requestedStoredEvents=cases['storedEvents'],
                                           treeEntriesWillBeRead=0))
    with (out / 'extraction_plan.tsv').open('x') as stream:
        for row in plans:
            stream.write('\t'.join(row) + '\n')
    requirements = dict(settings['histograms'])
    requirements.update({'absoluteTolerance': settings['pairProjectionAbsoluteTolerance'],
                         'relativeTolerance': settings['pairProjectionRelativeTolerance'], 'baselineJobid': cases['baselineJobid']})
    requirements.update({'mix.' + k: v for k, v in settings['requiredMixing'].items()})
    requirements.update({'reduced.' + k: v for k, v in settings['requiredReducedSettings'].items()})
    requirements.update({'cut.' + k: v for k, v in settings['requiredPairCut'].items()})
    with (out / 'requirements.tsv').open('x') as stream:
        for k, v in sorted(requirements.items()):
            stream.write(k + '\t' + str(v) + '\n')
    mainconf = resolve(settings['mainconfForRoot'], PROJECT)
    try:
        root_run(mainconf, HERE / 'ExtractPeriodCF.C', [out / 'extraction_plan.tsv', out / 'requirements.tsv', out / 'raw_histograms.json'], out / 'extract.log')
        raw = load(out / 'raw_histograms.json')
        need(raw['passed'] and raw['treeEntriesRead'] == 0, 'Extractor did not report histogram-only success')
        need({c['caseId'] for c in raw['cases']} == {c['caseId'] for c in selected}, 'Extracted case coverage mismatch')
        if args.mode == 'compare':
            reference_se = {}
            for c in raw['cases']:
                for s in c['slices']:
                    if s['plane'] != -1:
                        continue
                    key = (c['cutEnabled'], s['cent9'])
                    if key not in reference_se:
                        reference_se[key] = s['SE']
                    need(all(abs(a - b) <= settings['pairProjectionAbsoluteTolerance'] +
                             settings['pairProjectionRelativeTolerance'] * max(a, b)
                             for a, b in zip(s['SE'], reference_se[key])),
                         'SE signal counts differ across EP bin counts for the same cut and cent9')
        curves = []
        for c in raw['cases']:
            need(c['sources'] == expected and c['completedEvents'] == cases['storedEvents'], 'Extracted event/source total mismatch')
            lookup = {(s['plane'], s['cent9']): s for s in c['slices']}
            need(len(lookup) == (c['epBins'] + 1) * 9, 'Duplicate/missing native projection')
            for cent in centrals:
                planes = [-1] + (list(range(c['epBins'])) if cent['id'] in settings['diagnosticCentralities'] else [])
                for plane in planes:
                    slices = [lookup[(plane, ci)] for ci in range(cent['cent9'][0], cent['cent9'][1] + 1)]
                    curves.append(curve(c, cent, plane, slices, raw['nativeEdges'], settings))
        lookup = {(c['caseId'], c['centralityId'], c['plane']): c for c in curves}
        for selected_case in selected:
            for cent_id in settings['requiredValidCentralities']:
                main_curve = lookup[(selected_case['caseId'], cent_id, -1)]
                need(main_curve['available'] and main_curve['normSE'] > 0 and main_curve['normME'] > 0,
                     'Required main CF has no valid normalization: ' + selected_case['caseId'] + '/' + cent_id)
        comparisons, ep_comparisons = [], []
        if args.mode == 'compare':
            for cent in centrals:
                for n in [1, 8, 16, 32]:
                    pair = [next(c for c in selected if c['epBins'] == n and c['cutEnabled'] == on) for on in [False, True]]
                    off, on = [lookup[(c['caseId'], cent['id'], -1)] for c in pair]
                    row = dict(epBins=n, centralityId=cent['id'], cutOffCase=off['caseId'], cutOnCase=on['caseId'],
                               SERawOnOverOff=on['SERegular'] / off['SERegular'] if off['SERegular'] > 0 else None,
                               MERawOnOverOff=on['MERegular'] / off['MERegular'] if off['MERegular'] > 0 else None,
                               highK={})
                    for name in settings['highKWindowsGeVc']:
                        a, b = off['highK'][name], on['highK'][name]
                        row['highK'][name] = dict(cutOffIntegralCF=a['integralCF'], cutOnIntegralCF=b['integralCF'],
                                                  cutOnOverOffIntegralCF=(b['integralCF'] / a['integralCF']
                                                                         if a['integralCF'] is not None and a['integralCF'] > 0 and b['integralCF'] is not None else None))
                    comparisons.append(row)
                for on in [False, True]:
                    reference_case = next(c for c in selected if c['epBins'] == 8 and c['cutEnabled'] == on)
                    reference = lookup[(reference_case['caseId'], cent['id'], -1)]
                    for n in [1, 16, 32]:
                        candidate_case = next(c for c in selected if c['epBins'] == n and c['cutEnabled'] == on)
                        candidate = lookup[(candidate_case['caseId'], cent['id'], -1)]
                        row = dict(epBins=n, cutEnabled=on, centralityId=cent['id'],
                                   referenceCase=reference['caseId'], candidateCase=candidate['caseId'], highK={})
                        for name in settings['highKWindowsGeVc']:
                            a, b = reference['highK'][name], candidate['highK'][name]
                            row['highK'][name] = dict(referenceIntegralCF=a['integralCF'], candidateIntegralCF=b['integralCF'],
                                                      candidateOverEP8IntegralCF=(b['integralCF'] / a['integralCF']
                                                                                 if a['integralCF'] is not None and a['integralCF'] > 0 and b['integralCF'] is not None else None))
                        ep_comparisons.append(row)
        metrics = dict(schemaVersion=1, passed=True, mode=args.mode,
                       createdUtc=datetime.now(timezone.utc).isoformat(), cases=[c['caseId'] for c in selected],
                       normalization='alpha=SE(H)/ME(H); CF=SE/(alpha*ME); H is half-open',
                       centralityCombination=settings['centralityCombination'], signal='raw Lambda-He4 signal before sideband subtraction or purity correction',
                       eventPopulation='same stored EP-valid period sample, including EP1',
                       statistics='Display errors use raw histogram Sumw2 with alpha held fixed. Normalization covariance, event/shared-ME correlations and cross-case covariance are not included. No significance test is performed.',
                       highKMeaning='integralCF is the ratio of window counts using the full-curve normAlpha; arithmeticMeanOfValidDisplayBins is a separate unweighted descriptive mean. Outside-axis windows are unavailable, with no clipping or extrapolation.',
                       planeMeaning='plane=-1 is raw plane-count inclusive, normalized once; diagnostic planes are separately normalized. ME label is the currently processed event in both directions.',
                       inputEventsPerCase=cases['storedEvents'], treeEntriesRead=0,
                       nativeKRangeGeVc=[raw['nativeEdges'][0], raw['nativeEdges'][-1]],
                       settings=settings, curves=curves, cutComparisons=comparisons, epComparisons=ep_comparisons)
        save(out / 'cf_comparison.json', metrics)
        with (out / 'curve_metadata.tsv').open('x') as stream:
            for c in curves:
                values = [c['caseId'], c['centralityId'], c['centralityLabel'], c['plane'], c['epBins'], int(c['cutEnabled']),
                          'NA' if c['normAlpha'] is None else repr(c['normAlpha']), c['normSE'], c['normME']]
                stream.write('\t'.join(str(x) for x in values) + '\n')
        with (out / 'curve_bins.tsv').open('x') as stream:
            for c in curves:
                for i, b in enumerate(c['bins']):
                    values = [c['caseId'], c['centralityId'], c['plane'], i, b['low'], b['high'], b['SE'], b['SEVariance'], b['ME'], b['MEVariance'],
                              'NA' if b['CF'] is None else repr(b['CF']), 'NA' if b['fixedAlphaCountError'] is None else repr(b['fixedAlphaCountError'])]
                    stream.write('\t'.join(str(x) for x in values) + '\n')
        with (out / 'plot_settings.tsv').open('x') as stream:
            for name, vals in [('kRange', settings['displayKGeVc']), ('cfRange', settings['displayCFRange']), ('normRange', settings['normalizationGeVc'])]:
                stream.write(name + '\t' + '\t'.join(map(str, vals)) + '\n')
            stream.write('mode\t' + args.mode + '\n')
            for cent in centrals:
                stream.write('centrality\t' + cent['id'] + '\n')
            for cent in settings['diagnosticCentralities']:
                stream.write('diagnostic\t' + cent + '\n')
        root_run(mainconf, HERE / 'DrawPeriodCF.C', [out / 'curve_metadata.tsv', out / 'curve_bins.tsv', out / 'plot_settings.tsv', out], out / 'draw.log')
        for st in stats:
            after = stat_record(Path(st['path']))
            need(all(after[k] == st[k] for k in ['bytes', 'device', 'inode', 'mtimeNs']), 'Input ROOT changed during postprocessing')
        for pin in pins:
            need(digest(pin['path']) == pin['sha256'], 'Pinned helper/config/list changed during postprocessing')
        artifacts = []
        for name in ['cf_comparison.json', 'cf_comparison.pdf', 'cf_plane_diagnostic.pdf', 'raw_se_me_scaled.pdf', 'cf_histograms.root',
                     'raw_histograms.json', 'curve_metadata.tsv', 'curve_bins.tsv']:
            p = out / name
            st = stat_record(p)
            st['sha256'] = digest(p)
            artifacts.append(st)
        save(out / 'completed.json', dict(passed=True, mode=args.mode, cases=[c['caseId'] for c in selected],
                                          sources=len(plans), curves=len(curves), inputStatsUnchanged=True,
                                          pinsUnchanged=True, treeEntriesRead=0, artifacts=artifacts))
        print('CF_POSTPROCESS_PASS mode=%s cases=%d sources=%d curves=%d output=%s' %
              (args.mode, len(selected), len(plans), len(curves), out))
    except Exception as error:
        save(out / 'failed.json', dict(passed=False, error=str(error), mode=args.mode))
        raise


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('CF_POSTPROCESS_FAILED: ' + str(error), file=sys.stderr)
        sys.exit(1)
