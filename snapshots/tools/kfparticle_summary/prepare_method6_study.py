#!/usr/bin/env python3
"""Prepare shared-ME-purity CF statistics and fit-atlas tables from saved fits."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import sys

import yaml

from method6_statistics import analyze_group, METHODS
from sample_provenance import resolve_sample


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OWNER = 'tools/kfparticle_summary/prepare_method6_study.py'
META_COLUMNS = ('page group label method bin low high windowLow windowHigh side status '
                'covariance ndf chi2 probability purity purityError observations '
                'p0 p1 p2 p3 p4 p5 reason').split()
BIN_COLUMNS = 'page side index low high count error2 centralCount centralError2'.split()
VALUE_REL_TOLERANCE = 1.e-8
VALUE_ABS_TOLERANCE = 1.e-10


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), sha256=digest(path))


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns, inode=stat.st_ino)


def read_json(path):
    with Path(path).open(encoding='utf-8') as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError('Expected a JSON object: %s' % path)
    return value


def read_yaml(path):
    with Path(path).open(encoding='utf-8') as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError('Expected a YAML mapping: %s' % path)
    return value


def save_json(path, value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n',
                         encoding='utf-8')
    temporary.replace(path)


def text(value):
    return str(value).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')


def numeric(value):
    if value is None:
        return 'nan'
    value = float(value)
    return format(value, '.17g') if math.isfinite(value) else 'nan'


def save_tsv(path, columns, rows):
    lines = ['\t'.join(columns)]
    for row in rows:
        if len(row) != len(columns):
            raise ValueError('Wrong TSV row width for %s.' % path)
        lines.append('\t'.join(text(value) for value in row))
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    temporary.replace(path)


def config_records(mainconf):
    """Record mainconf and recursively reachable local YAML/calibration inputs."""
    output, visited = [], set()

    def resolve(value, parent):
        candidate = Path(value)
        candidates = ([candidate] if candidate.is_absolute() else
                      [PROJECT/'config'/candidate, PROJECT/candidate, parent/candidate])
        return next((p.resolve() for p in candidates if p.is_file()), None)

    def inspect(path):
        path = path.resolve(strict=True)
        if path in visited:
            return
        visited.add(path)
        output.append(record(path))
        if path.suffix.lower() not in ('.yaml', '.yml'):
            return
        doc = read_yaml(path)

        def walk(value, key=''):
            if isinstance(value, dict):
                for name, child in value.items():
                    walk(child, str(name))
            elif isinstance(value, list):
                for child in value:
                    walk(child, key)
            elif isinstance(value, str):
                suffix = Path(value).suffix.lower()
                if suffix in ('.yaml', '.yml') or key == 'calibrationFile':
                    target = resolve(value, path.parent)
                    if target is None:
                        raise ValueError('Missing referenced configuration %s in %s.' % (value, path))
                    inspect(target)
        walk(doc)
    inspect(Path(mainconf))
    return sorted(output, key=lambda item: item['path'])


def validate_settings(settings, summary, extract, groups):
    if settings.get('method6_purity_source') != 'method5_ME' or settings.get('shared_purity_error') is not True:
        raise ValueError('Method6 requires the shared Method5 ME purity estimator.')
    if settings.get('fit_model') != 'gaus(0)+pol2(3)':
        raise ValueError('Unexpected mass-fit model.')
    actual = extract['settings']
    checks = [('rebin', summary['rebin_group']), ('purityRebin', summary['purity_rebin']),
              ('slidingHalfWidth', summary['sliding_half_width'])]
    for key, expected in checks:
        if actual.get(key) != expected:
            raise ValueError('Extraction settings.%s differs from the summary YAML.' % key)
    # Additional names depend on the existing WriteSettings schema; verify those
    # present rather than silently replacing the extractor's recorded settings.
    fit = summary['mass_fit']
    aliases = [('fitRange', fit['range']), ('fitMeanLimits', fit['mean_limits']),
               ('fitSigmaLimits', fit['sigma_limits']), ('fitSigmaInitial', fit['sigma_initial']),
               ('minimumEntries', fit['minimum_entries']), ('minimumCovarianceStatus', fit['minimum_covariance_status']),
               ('normalization', summary['normalization'])]
    for key, expected in aliases:
        if key not in actual or actual[key] != expected:
            raise ValueError('Extraction settings.%s differs from the summary YAML.' % key)
    if not math.isclose(float(extract['kstarMax']), float(settings['extract_kstar_max']), abs_tol=1.e-12):
        raise ValueError('Extracted kstar maximum differs from the study settings.')
    known = set(groups)
    for group_set in settings['group_sets']:
        if not group_set or any(group_id not in known for group_id in group_set):
            raise ValueError('group_sets references an absent group.')


def compare_reference(reference, results, reference_dir, current_fingerprint):
    context_path = reference_dir/'mass_halves'/'context.json'
    old_context = read_json(context_path)
    old_input = old_context.get('input')
    if old_input != current_fingerprint:
        raise ValueError('Reference-study ROOT fingerprint differs from the current input.')
    if Path(reference['inputRoot']).resolve(strict=True) != Path(current_fingerprint['path']):
        raise ValueError('Reference comparison used a different ROOT input.')
    old_groups = {group['id']: group for group in reference['groups']}
    summary = dict(compatible=True, valueRelativeTolerance=VALUE_REL_TOLERANCE,
                   valueAbsoluteTolerance=VALUE_ABS_TOLERANCE, comparedBins=0,
                   matchedNA=0, errorChanges=0, maxAbsoluteValueDifference=0.,
                   maxAbsoluteErrorDifference=0., mismatches=[])
    for group in results:
        if group['id'] not in old_groups:
            raise ValueError('Missing reference group: '+group['id'])
        previous = old_groups[group['id']]
        for method in ('method3', 'method5'):
            old = previous['methods'][method]['curves']['full']
            new = group['methods'][method]['points']
            if len(old) != len(new):
                raise ValueError('Reference CF-bin count differs for %s/%s.' % (group['id'], method))
            for index, (before, after) in enumerate(zip(old, new)):
                before, after = before['cf'], after['cf']
                summary['comparedBins'] += 1
                identity = dict(group=group['id'], method=method, bin=index+1)
                if before is None or after is None:
                    if before is None and after is None:
                        summary['matchedNA'] += 1
                    else:
                        summary['mismatches'].append(dict(identity, reason='NA status differs',
                                                          reference=before, current=after))
                    continue
                difference = abs(before['value']-after['value'])
                summary['maxAbsoluteValueDifference'] = max(summary['maxAbsoluteValueDifference'], difference)
                if not math.isclose(before['value'], after['value'], rel_tol=VALUE_REL_TOLERANCE,
                                    abs_tol=VALUE_ABS_TOLERANCE):
                    summary['mismatches'].append(dict(identity, reason='CF value differs',
                                                      reference=before['value'], current=after['value']))
                error_difference = abs(before['error']-after['error'])
                summary['maxAbsoluteErrorDifference'] = max(summary['maxAbsoluteErrorDifference'], error_difference)
                if error_difference > 1.e-14:
                    summary['errorChanges'] += 1
    summary['compatible'] = not summary['mismatches']
    summary['referenceComparison'] = record(reference_dir/'mass_halves'/'comparison.json')
    summary['referenceContext'] = record(context_path)
    return summary


def sideband_proxy(extract, row):
    try:
        bounds = extract['massPartition']['binned']
        edges, counts, error2 = extract['massEdges'], row['massCounts'], row['massError2']
        estimates = {}
        for region in ('full', 'leftSB', 'rightSB'):
            low, high = bounds[region]['lower'], bounds[region]['upper']
            indices = [i for i in range(len(counts)) if low <= (edges[i]+edges[i+1])/2 < high]
            estimates[region] = dict(count=sum(counts[i] for i in indices),
                                     variance=sum(error2[i] for i in indices), width=high-low)
        signal, left, right = (estimates[key] for key in ('full', 'leftSB', 'rightSB'))
        if signal['count'] <= 0. or left['width']+right['width'] <= 0.:
            return dict(value=None, error=None, reason='nonpositive signal count or sideband width')
        alpha = signal['width']/(left['width']+right['width'])
        side_count, side_variance = left['count']+right['count'], left['variance']+right['variance']
        value = 1-alpha*side_count/signal['count']
        variance = alpha**2*(side_variance/signal['count']**2+
                              side_count**2*signal['variance']/signal['count']**4)
        return dict(value=value, error=math.sqrt(variance), variance=variance, alpha=alpha,
                    signal=signal, left=left, right=right,
                    reason='flat sideband diagnostic; disjoint mass-bin errors only')
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        return dict(value=None, error=None, reason='sideband proxy unavailable: '+str(exc))


def fit_side(extract, row, purity, window):
    fit = row['fit']
    count = sum(row['windowNtot']) if row.get('windowNtot') is not None else None
    central = sum(row['centralNtot']) if row.get('centralNtot') is not None else None
    chi2, ndf = fit.get('chi2'), fit.get('ndf')
    return dict(purity=purity, status=fit.get('status'), covarianceStatus=fit.get('covarianceStatus'),
                ndf=ndf, chi2=chi2, chi2Ndf=chi2/ndf if chi2 is not None and ndf and ndf > 0 else None,
                probability=fit.get('probability'), observations=fit.get('observationsInStoredMassRange'),
                numericallyValid=fit.get('numericallyValid'), isValid=fit.get('isValid'),
                boundary=fit.get('boundary'), reason=fit.get('reason'),
                window=list(window), signalWindowCount=count, centralSignalWindowCount=central,
                sidebandProxy=sideband_proxy(extract, row),
                parameters=fit.get('parameters'), parameterCovariance=fit.get('covariance'))


def fit_tables(extract, settings, group_definitions, results):
    mass_edges = extract['massEdges']
    meta_rows, bin_rows, pairs = [], [], []
    result_by_id = {group['id']: group for group in results}
    extracted = {group['id']: group for group in extract['groups']}
    page = 0
    inspection = settings['inspection_range']
    for method in ('method3', 'method5'):
        for definition in group_definitions:
            group_id, label = definition['id'], definition['label']
            group, calculation = extracted[group_id], result_by_id[group_id]['methods'][method]
            edges = group['kstarEdges'] if method == 'method3' else group['method5PurityEdges']
            for index, (low, high) in enumerate(zip(edges, edges[1:])):
                center = (low+high)/2
                if not inspection[0] <= center < inspection[1]:
                    continue
                page += 1
                pair = dict(page=page, group=group_id, label=label, method=method, bin=index+1,
                            low=low, high=high, x=center, SE=None, ME=None)
                for side in ('SE', 'ME'):
                    row = group['methods'][method][side][index]
                    if method == 'method5':
                        native = group['originalKstarEdges']
                        window = (native[row['windowOriginalFirstBin']-1], native[row['windowOriginalLastBin']])
                    else:
                        window = (edges[row['windowFirstBin']-1], edges[row['windowLastBin']])
                    purity = calculation['purity'][side][index]
                    diagnostic = fit_side(extract, row, purity, window)
                    pair[side] = diagnostic
                    fit = row['fit']
                    parameters = fit.get('parameters') or [None]*6
                    if len(parameters) != 6:
                        raise ValueError('A Gaussian+pol2 fit must have six parameters.')
                    reason = fit.get('reason') or row.get('full', {}).get('reason') or ''
                    meta_rows.append([page, group_id, label, method, index+1, numeric(low), numeric(high),
                                      numeric(window[0]), numeric(window[1]), side,
                                      numeric(fit.get('status')), numeric(fit.get('covarianceStatus')),
                                      numeric(fit.get('ndf')), numeric(fit.get('chi2')),
                                      numeric(fit.get('probability')),
                                      numeric(purity['value'] if purity else None),
                                      numeric(purity['error'] if purity else None),
                                      numeric(fit.get('observationsInStoredMassRange'))]+
                                     [numeric(value) for value in parameters]+[reason])
                    arrays = [row[key] for key in ('massCounts', 'massError2',
                                                   'centralMassCounts', 'centralMassError2')]
                    if any(len(values) != len(mass_edges)-1 for values in arrays):
                        raise ValueError('Mass-bin arrays disagree with massEdges.')
                    for mass_index in range(len(mass_edges)-1):
                        values = [values[mass_index] for values in arrays]
                        if any(not math.isfinite(float(value)) for value in values):
                            raise ValueError('Nonfinite stored mass-bin measurement.')
                        bin_rows.append([page, side, mass_index+1, numeric(mass_edges[mass_index]),
                                         numeric(mass_edges[mass_index+1])]+[numeric(value) for value in values])
                se, me = pair['SE']['purity'], pair['ME']['purity']
                difference, difference_error, pull = None, None, None
                if se is not None and me is not None:
                    cross = calculation['purity']['covariance']['SE_ME'][index][index]
                    variance = se['variance']+me['variance']-2*cross
                    difference = se['value']-me['value']
                    difference_error = math.sqrt(max(variance, 0.))
                    pull = difference/difference_error if difference_error > 0. else None
                pair.update(purityDifference=difference, purityDifferenceError=difference_error,
                            descriptivePull=pull, pullIsDescriptive=True)
                pairs.append(pair)
    return meta_rows, bin_rows, pairs


def quality_tsv(pairs):
    columns = ('page group label method bin low high x PSE errorPSE PME errorPME '
               'difference differenceError descriptivePull').split()
    side_columns = ('status covarianceStatus ndf chi2 chi2Ndf probability observations '
                    'windowLow windowHigh signalWindowCount centralSignalWindowCount sidebandPurity sidebandPurityError').split()
    columns += [side+'_'+key for side in ('SE', 'ME') for key in side_columns]
    rows = []
    for pair in pairs:
        se, me = pair['SE']['purity'], pair['ME']['purity']
        row = [pair[key] for key in ('page', 'group', 'label', 'method', 'bin')]
        row += [numeric(pair[key]) for key in ('low', 'high', 'x')]
        row += [numeric(se['value'] if se else None), numeric(se['error'] if se else None),
                numeric(me['value'] if me else None), numeric(me['error'] if me else None)]
        row += [numeric(pair[key]) for key in ('purityDifference', 'purityDifferenceError', 'descriptivePull')]
        for side in ('SE', 'ME'):
            diagnostic = pair[side]
            row += [numeric(diagnostic.get(key)) for key in side_columns[:7]]
            row += [numeric(diagnostic['window'][0]), numeric(diagnostic['window'][1]),
                    numeric(diagnostic['signalWindowCount']), numeric(diagnostic['centralSignalWindowCount']),
                    numeric(diagnostic['sidebandProxy']['value']), numeric(diagnostic['sidebandProxy']['error'])]
        rows.append(row)
    return columns, rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-root', 'mainconf', 'output-manifest', 'runmeta', 'settings',
                 'output-dir', 'reference-study'):
        parser.add_argument('--'+name, required=True, type=Path)
    args = parser.parse_args(argv)
    out, settings_path = args.output_dir.resolve(strict=True), args.settings.resolve(strict=True)
    input_path, mainconf = args.input_root.resolve(strict=True), args.mainconf.resolve(strict=True)
    reference_dir = args.reference_study.resolve(strict=True)
    initial_fingerprint = fingerprint(input_path)
    sample = resolve_sample('3p85', input_path, args.output_manifest, args.runmeta)
    settings = read_yaml(settings_path)
    summary_path = (settings_path.parent/settings['summary_settings']).resolve(strict=True)
    summary = read_yaml(summary_path)
    configurations = config_records(mainconf)
    extract_path = out/'purity_fit_review.json'
    extract = read_json(extract_path)
    if Path(extract['inputRoot']).resolve(strict=True) != input_path:
        raise ValueError('Purity-fit extraction used a different ROOT input.')
    from compare_mass_halves import settings_text
    main_values = read_yaml(mainconf)
    maker_values = read_yaml(PROJECT/'config'/main_values['maker'])
    extract_settings_path = Path(extract['settingsFile']).resolve(strict=True)
    if extract_settings_path != (out/'extract_settings.txt').resolve(strict=True):
        raise ValueError('Purity-fit extraction used a different explicit settings file.')
    if extract_settings_path.read_text() != settings_text(maker_values, settings, summary):
        raise ValueError('Extracted flat settings do not match the Maker/study/summary YAMLs.')
    for actual, key in [('mean','lambdaSignalMean'),('sigma','lambdaSignalSigma'),
                        ('nSigma','lambdaSignalNSigma'),('sidebandOuterFactor','lambdaSidebandOuterFactor')]:
        if extract['settings'][actual] != maker_values[key]:
            raise ValueError('Extracted mass gate differs from Maker YAML: '+actual)
    copied_sources = {}
    for filename in ('ExtractPurityFitReview.C', 'ExtractMassHalfCF.C'):
        source, copied = HERE/filename, out/'macros'/filename
        if digest(source) != digest(copied):
            raise ValueError('Copied extractor source differs from the current source: '+filename)
        copied_sources[filename] = dict(source=record(source), copied=record(copied))
    definitions = settings['groups']
    if len(definitions) != 5 or len({group['id'] for group in definitions}) != 5:
        raise ValueError('This study requires five uniquely named centrality groups.')
    extracted_groups = {group['id']: group for group in extract['groups']}
    if len(extracted_groups) != len(extract['groups']) or set(extracted_groups) != {group['id'] for group in definitions}:
        raise ValueError('Extracted centrality groups differ from the study settings.')
    validate_settings(settings, summary, extract, extracted_groups)
    results = []
    for definition in definitions:
        group = extracted_groups[definition['id']]
        if group['centMin'] != definition['min'] or group['centMax'] != definition['max']:
            raise ValueError('Extracted centrality bounds differ: '+definition['id'])
        calculation = analyze_group(group, extract['settings'], summary['normalization'],
                                    summary['plot_x_range'], settings['comparison_range'])
        calculation['label'] = definition['label']
        results.append(calculation)
    reference_path = reference_dir/'mass_halves'/'comparison.json'
    reference_check = compare_reference(read_json(reference_path), results, reference_dir, initial_fingerprint)
    sources = {name: record(HERE/name) for name in (
        'prepare_method6_study.py', 'method6_statistics.py', 'compare_mass_halves.py',
        'mass_half_statistics.py', 'sample_provenance.py')}
    context = dict(owner=OWNER, schemaVersion=1, phase='prepared',
                   preparedAt=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   input=initial_fingerprint, sample=sample,
                   mainconf=record(mainconf), configurationFiles=configurations,
                   settings=record(settings_path), summarySettings=record(summary_path),
                   extract=record(extract_path), extractSettings=record(extract_settings_path), copiedExtractors=copied_sources,
                   sourceFiles=sources, referenceStudy=str(reference_dir), referenceComparison=reference_check)
    context_path = out/'context.json'
    if context_path.exists():
        previous = read_json(context_path)
        if previous.get('owner') != OWNER or previous.get('input') != initial_fingerprint:
            raise ValueError('Refusing to replace an unrelated study context.')
    if not reference_check['compatible']:
        context['phase'] = 'reference_mismatch'
        save_json(context_path, context)
        raise ValueError('Fresh Method3/5 CF values or NA status differ from the reference; see context.json.')
    meta_rows, bin_rows, pairs = fit_tables(extract, settings, definitions, results)
    if len(pairs) != 140:
        raise ValueError('Expected 140 inspection pages; generated %d.' % len(pairs))
    if fingerprint(input_path) != initial_fingerprint:
        raise RuntimeError('Input ROOT changed while preparing the study.')
    comparison = dict(schemaVersion=1, inputRoot=str(input_path), sample=sample,
                      settings=settings, summarySettings=summary, group_sets=settings['group_sets'],
                      normalization=summary['normalization'], plot_range=summary['plot_x_range'],
                      comparison_range=settings['comparison_range'], inspection_range=settings['inspection_range'],
                      referenceComparison=reference_check, groups=results)
    cf_rows, purity_rows = [], []
    for group in results:
        for method in METHODS:
            calculation = group['methods'][method]
            for point in calculation['points']:
                estimate = point['cf']
                cf_rows.append([method, group['id'], group['label'], point['bin'], numeric(point['low']),
                                numeric(point['high']), numeric(estimate['value'] if estimate else None),
                                numeric(estimate['error'] if estimate else None)])
            purity = calculation['purity']
            for side in ('SE', 'ME'):
                for index, estimate in enumerate(purity[side]):
                    purity_rows.append([method, group['id'], group['label'], side, index+1,
                                        numeric(purity['edges'][index]), numeric(purity['edges'][index+1]),
                                        numeric(estimate['value'] if estimate else None),
                                        numeric(estimate['error'] if estimate else None)])
    save_json(out/'methods_comparison.json', comparison)
    save_tsv(out/'cf_points.tsv', 'method group label bin low high value error'.split(), cf_rows)
    save_tsv(out/'purity_points.tsv', 'method group label side bin low high value error'.split(), purity_rows)
    save_tsv(out/'fit_meta.tsv', META_COLUMNS, meta_rows)
    save_tsv(out/'fit_bins.tsv', BIN_COLUMNS, bin_rows)
    quality = dict(schemaVersion=1, inputRoot=str(input_path), sample=sample,
                   inspection_range=settings['inspection_range'], totalPages=len(pairs),
                   method6AgreementIsByConstruction=True, pairs=pairs,
                   caveats=['Purity difference pulls are descriptive; no formal agreement p-value is estimated.',
                            'Mass-fit chi2/probability diagnose the mass-shape model, not SE/ME purity agreement.',
                            'The flat sideband proxy is a diagnostic and does not replace fitted purity.'])
    save_json(out/'fit_quality.json', quality)
    quality_columns, quality_rows = quality_tsv(pairs)
    save_tsv(out/'fit_quality.tsv', quality_columns, quality_rows)
    context['atlasPages'] = len(pairs)
    context['products'] = {name: record(out/name) for name in (
        'methods_comparison.json', 'cf_points.tsv', 'purity_points.tsv',
        'fit_meta.tsv', 'fit_bins.tsv', 'fit_quality.json', 'fit_quality.tsv')}
    save_json(context_path, context)
    print(json.dumps(dict(outputDir=str(out), groups=len(results), atlasPages=len(pairs),
                          referenceCompatible=reference_check['compatible']), sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
