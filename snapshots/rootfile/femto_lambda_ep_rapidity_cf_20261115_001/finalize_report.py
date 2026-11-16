#!/usr/bin/env python3
"""Verify the final diagnostic report and record its artifact provenance."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import hashlib
import json
import re
import subprocess
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / 'share/figure/20261115_lambda_he4_ep_rapidity_cf'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def describe(path):
    return {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size,
            'sha256': sha(path)}


def main():
    original = json.loads((HERE / 'input_manifest.json').read_text())
    results = json.loads((HERE / 'analysis_results.json').read_text())
    closures = json.loads((HERE / 'closure_checks.json').read_text())
    identities = {}
    for dataset, entry in original['inputs'].items():
        st = Path(entry['path']).stat()
        actual = dict(size=st.st_size, mtimeNs=st.st_mtime_ns,
                      inode=st.st_ino, device=st.st_dev)
        unchanged = all(actual[k] == entry[k] for k in actual)
        assert unchanged, 'Input identity changed: ' + dataset
        identities[dataset] = dict(path=entry['path'], identity=actual,
                                   unchanged=True,
                                   check='size, mtime_ns, inode, device; no full-file SHA computed')
    assert len(closures) == 108
    assert all(c['passed_regular'] and c['max_regular_count_difference'] == 0
               for c in closures)
    pdf = OUT / 'cf_ep_rapidity_report_ja.pdf'
    page_count = int(subprocess.check_output([
        'gs', '-q', '-dNODISPLAY', '-dNOSAFER', '-sPDFname=' + str(pdf),
        '-c', 'PDFname (r) file runpdfbegin pdfpagecount = quit'], text=True).strip())
    assert page_count == 22
    assert len(list(OUT.glob('page-[0-9][0-9].png'))) == 22
    text_path = OUT / 'cf_ep_rapidity_report_ja.txt'
    subprocess.run(['gs', '-q', '-dBATCH', '-dNOPAUSE', '-sDEVICE=txtwrite',
                    '-sOutputFile=' + str(text_path), str(pdf)], check=True)
    text = text_path.read_text()
    assert 'ラピディティ' in text and all(x in text for x in ('0.866', '0.906', '0.968'))
    build_log = (OUT / 'build-2.log').read_text()
    warnings = {x: build_log.count(x) for x in
                ('Missing character', 'Overfull', 'Underfull', 'Warning')}
    assert not any(warnings.values()), warnings
    decoded = bytearray(pdf.read_bytes())
    for match in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', bytes(decoded), re.S):
        try:
            decoded.extend(zlib.decompress(match.group(1)))
        except zlib.error:
            pass
    font_checks = {name: name.encode() in decoded for name in
                   ('DroidSansFallback', '/FontFile', '/ToUnicode')}
    assert all(font_checks.values()), font_checks
    crop_manifest = OUT / 'plots_cropped/crop_manifest.json'
    assert crop_manifest.is_file()
    assert len(list((OUT / 'plots').glob('*.pdf'))) == 180
    assert len(list((OUT / 'plots_cropped').glob('*.pdf'))) == 180
    files = [HERE / x for x in (
        'study.yaml', 'input_manifest.json', 'projection_plan.tsv',
        'prepare_study.py', 'extract_cf_inputs.C', 'extract_retry.log',
        'raw_projections.tsv', 'raw_projections.root', 'source_axes.tsv',
        'analyze_cf.py', 'analysis_results.json', 'metrics.csv',
        'cf_covariances.npz', 'closure_checks.json', 'curve_bins.tsv',
        'prepare_plots.py', 'plot_metadata.tsv', 'plot_panels.tsv',
        'plot_summary.tsv', 'make_cf_plots.C', 'plots.log',
        'cf_histograms.root', 'build_report.py', 'finalize_report.py')]
    files += [pdf, OUT / 'cf_ep_rapidity_report_ja.tex', text_path,
              OUT / 'build-2.log', crop_manifest]
    reference = Path(results['ratio_source'])
    if not reference.is_absolute():
        reference = ROOT / reference
    assert sha(reference) == results['ratio_source_sha256']
    manifest = {
        'finalizedAtJst': datetime.now(timezone(timedelta(hours=9))).isoformat(),
        'report': describe(pdf), 'pages': page_count,
        'dataset': 'STAR Au+Au 3p85GeV fixed target / 2021 P24iy / Lambda-He4 signal',
        'scope': dict(common_success_groups=original['commonJobs'],
                      pico_input_references=72424, excluded_groups=35,
                      calibrated_run_catalogue=388, uncalibrated_runs_on_hold=81),
        'inputIdentityChecks': identities,
        'analysisSettings': results['config'],
        'statisticalAssumptions': results['statistical_assumption'],
        'referenceRatioSource': describe(reference),
        'checks': dict(regular_projection_closures=len(closures),
                       maximum_regular_count_difference=0,
                       numerical_independent_review='132 curves and 184 GLS recomputations passed before final typesetting',
                       final_pdf_independent_agent_review=False,
                       rendered_pages=22, japanese_text_extraction=True,
                       fonts=font_checks, typesetting_warning_counts=warnings,
                       visual_review='Contact sheet of all 22 pages and individual pages 1,3,4,7,8,13,16,20; root-agent inspection',
                       vector_plots=180, cropped_vector_plots=180),
        'limitations': [
            'Fixed normalization [0.4,0.5) GeV/c is a comparison convention, not a demonstrated uncorrelated plateau.',
            'Evaluation ends at k*=1 GeV/c; asymptotic convergence is not established.',
            'Rapidity cuts are common SE/ME pair selections; only EP changes the mixing pool.',
            'He4 rapidity and pair rapidity are not stored jointly.',
            'Event/pair reuse and shared EP1/EP8 covariance are unavailable.',
            'k* overflow ME residuals EP1=3 and EP8=2 counts lie outside the CF range; cause not directly established.'
        ],
        'artifacts': [describe(path) for path in files],
    }
    target = OUT / 'report_manifest.json'
    target.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(dict(manifest=str(target), pages=page_count,
                          pdf_sha256=sha(pdf), input_identity_unchanged=True,
                          regular_closures=108, warnings=warnings), ensure_ascii=False))


if __name__ == '__main__':
    main()
