"""Validate saved Method6 provenance, shared errors and numerical identities."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from build_report import pdf_pages
from sample_provenance import resolve_sample


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same(a, b):
    return math.isclose(a, b, rel_tol=1.e-9, abs_tol=1.e-12)


def check_record(record):
    if digest(record['path']) != record['sha256']:
        raise ValueError('Changed dependency: '+record['path'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve(strict=True)
    context = read(out/'context.json')
    results = read(out/'methods_comparison.json')
    review = read(out/'purity_fit_review.json')
    sample = results['sample']
    proof = sample['provenance']
    if sample != resolve_sample('3p85', results['inputRoot'], proof['outputManifest']['path'], proof['runmeta']['path']):
        raise ValueError('Sample provenance changed.')
    original = context['input']
    stat = Path(original['path']).stat()
    current = dict(path=original['path'], size=stat.st_size, mtime_ns=stat.st_mtime_ns, inode=stat.st_ino)
    if current != original or context['sample'] != sample:
        raise ValueError('Input fingerprint or sample changed.')
    for key in ('mainconf', 'settings', 'summarySettings', 'extract', 'extractSettings'):
        check_record(context[key])
    for record in context['configurationFiles']:
        check_record(record)
    for record in context['sourceFiles'].values():
        check_record(record)
    for pair in context['copiedExtractors'].values():
        check_record(pair['source'])
        check_record(pair['copied'])
    for record in context['products'].values():
        check_record(record)
    reference = results['referenceComparison']
    if not reference['compatible'] or reference['comparedBins'] != 500:
        raise ValueError('Method3/5 reference comparison was not completed.')
    check_record(reference['referenceComparison'])
    check_record(reference['referenceContext'])
    sources = {g['id']: g for g in review['groups']}
    checked_ratio, checked_points, unavailable = 0, 0, []
    settings = review['settings']
    for group in results['groups']:
        original_group = sources[group['id']]
        m5, m6 = group['methods']['method5'], group['methods']['method6']
        p = m6['purity']
        if not p['shared'] or not p['identicalSharedP'] or p['source'] != 'ME':
            raise ValueError('Method6 did not declare shared ME purity.')
        if p['SE'] != p['ME'] or p['SE'] != m5['purity']['ME']:
            raise ValueError('Method6 purity/error were not copied unchanged from Method5 ME.')
        if not (p['covariance']['SE'] == p['covariance']['ME'] == p['covariance']['SE_ME']):
            raise ValueError('Copied purity uncertainties were treated as independent.')
        for method, result in group['methods'].items():
            points, covariance = result['points'], result['covariance']
            if not result['normalization']['valid']:
                unavailable.append(dict(group=group['id'], method=method, reason=result['normalization']['reason']))
                if any(point['cf'] is not None for point in points) or any(value is not None for row in covariance for value in row):
                    raise ValueError('Incomplete normalization did not mask the entire CF.')
            for i, point in enumerate(points):
                estimate = point['cf']
                if estimate is None:
                    if any(value is not None for value in covariance[i]):
                        raise ValueError('Invalid CF has non-null covariance.')
                    continue
                checked_points += 1
                if not same(estimate['variance'], estimate['error']**2) or not same(estimate['variance'], covariance[i][i]):
                    raise ValueError('CF error and covariance diagonal differ.')
                for j, value in enumerate(covariance[i]):
                    if value is not None and not same(value, covariance[j][i]):
                        raise ValueError('CF covariance is not symmetric.')
        for i, point in enumerate(m6['points']):
            first, last = i*settings['rebin'], (i+1)*settings['rebin']
            cells = {index//settings['purityRebin'] for index in range(first, last)}
            if any(p['ME'][cell] is None for cell in cells):
                if m6['yields']['SE'][i] is not None or m6['yields']['ME'][i] is not None:
                    raise ValueError('Failed ME purity was not shared as a mask.')
                continue
            if len(cells) != 1 or point['yieldRatio'] is None:
                continue
            totals, variances = {}, {}
            for side in ('SE', 'ME'):
                native = original_group['method5Base'][side][first:last]
                totals[side] = sum(sum(row['counts']) for row in native)
                variances[side] = sum(sum(row['error2']) for row in native)
            if totals['ME'] <= 0:
                continue
            raw = totals['SE']/totals['ME']
            variance = variances['SE']/totals['ME']**2 + totals['SE']**2*variances['ME']/totals['ME']**4
            if not same(point['yieldRatio']['value'], raw) or not same(point['yieldRatio']['variance'], variance):
                raise ValueError('Common-P cancellation failed within one purity cell.')
            checked_ratio += 1
    report = read(out/'report_metadata.json')
    if report['phase'] != 'complete':
        raise ValueError('Report is not complete.')
    check_record(report['report'])
    for record in report['sources'].values():
        check_record(record)
    for record in report['figures'].values():
        check_record(record)
    report_path = Path(report['report']['path'])
    if digest(report_path.with_suffix('.tex')) != report['texSha256'] or digest(out/'README.md') != report['readmeSha256']:
        raise ValueError('Report text or index changed after assembly.')
    cropping = read(out/'crop_manifest.json')
    if cropping['phase'] != 'complete' or len(cropping['files']) != 144:
        raise ValueError('Vector plot geometry processing is incomplete.')
    for name, record in cropping['files'].items():
        if digest(out/name) != record['outputSha256'] or digest(Path(record['source'])) != record['sourceSha256']:
            raise ValueError('Raw/final vector figure changed: '+name)
    figures = {}
    for name, pages in [('cf_methods_3_5_6.pdf', 2), ('purity_se_me.pdf', 2), ('atlas_method3.pdf', 100), ('atlas_method5.pdf', 40)]:
        found = pdf_pages(out/name)
        if found != pages:
            raise ValueError('PDF page count differs: '+name)
        figures[name] = dict(pages=found, sha256=digest(out/name), bytes=(out/name).stat().st_size)
    singles = sorted(out.glob('fit_method[35]_*.pdf'))
    if len(singles) != 140:
        raise ValueError('Expected 140 individual mass-fit pages.')
    result = dict(phase='validated', input_unchanged=True, groups=len(results['groups']),
                  referenceComparedBins=500, checkedCfPoints=checked_points,
                  sharedPurityCancellationBins=checked_ratio, normalizationUnavailable=unavailable,
                  plots=figures, individualFitPdfs=len(singles))
    (out/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
