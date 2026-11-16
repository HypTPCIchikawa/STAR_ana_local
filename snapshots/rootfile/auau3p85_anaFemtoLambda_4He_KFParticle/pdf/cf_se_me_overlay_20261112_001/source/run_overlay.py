#!/usr/bin/env python3
"""Run ROOT plotting from explicit study settings and verify normalized spectra."""
import argparse
import csv
import datetime
import hashlib
import json
import math
import re
from types import SimpleNamespace
from pathlib import Path
import subprocess


def fingerprint(path):
    s = path.stat()
    return dict(size=s.st_size, mtimeNs=s.st_mtime_ns, ino=s.st_ino, device=s.st_dev)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run(settings_path, validate_existing=False):
    settings_path = settings_path.resolve()
    out = settings_path.parent
    cfg = json.loads(settings_path.read_text())
    project = Path(cfg['source_cf_study']).parents[3]
    input_root = Path(cfg['input_root'])
    before = fingerprint(input_root)
    assert before == cfg['input_fingerprint'], (before, cfg['input_fingerprint'])
    if validate_existing:
        assert (out / 'overlay_histograms.root').exists()
        # A prior successful ROOT invocation is also required by its explicit log marker below.
        proc = SimpleNamespace(returncode=0)
    else:
        assert not (out / 'overlay_histograms.root').exists(), 'Output ROOT already exists'
        group_file = out / 'groups.tsv'
        lines = []
        for group in cfg['groups']:
            cents = group['cent9']
            assert cents == list(range(min(cents), max(cents) + 1))
            lines.append('{}\t{}\t{}\t{}'.format(group['id'], group['label'], min(cents), max(cents)))
        group_file.write_text('\n'.join(lines) + '\n')
        strings = [str(input_root), str(out), str(group_file), cfg['se_histogram'], cfg['me_histogram'],
                   cfg['se_qa_template'], cfg['me_qa_template'], cfg['input_uuid']]
        numbers = cfg['normalization'] + [cfg['rebin_factor']] + sum(cfg['display_ranges'], [])
        call = 'DrawSEME({});'.format(', '.join([json.dumps(s) for s in strings] + [repr(x) for x in numbers]))
        commands = '.L {}+\n{}\n.q\n'.format(out / 'source' / 'DrawSEME.C', call)
        (out / 'root_commands.txt').write_text(commands)
        argv = ['bash', 'script/singularity_root_batch.sh', cfg['mainconf'], 'root4star', '-l', '-b']
        with (out / 'root.log').open('w') as log:
            proc = subprocess.run(argv, input=commands, universal_newlines=True, stdout=log,
                                  stderr=subprocess.STDOUT, cwd=str(project))
    after = fingerprint(input_root)
    assert after == before, 'Input ROOT fingerprint changed'
    assert proc.returncode == 0, 'ROOT failed; see root.log'
    log_text = (out / 'root.log').read_text()
    assert 'SEME_SUCCESS pages={}'.format(len(cfg['groups'])) in log_text, 'Missing success marker'
    assert 'SEME_FAILURE' not in log_text
    with (out / 'normalization.tsv').open() as f:
        norms = list(csv.DictReader(f, delimiter='\t'))
    with (out / 'spectra.tsv').open() as f:
        spectra = list(csv.DictReader(f, delimiter='\t'))
    assert len(norms) == len(cfg['groups'])
    checks = []
    for row in norms:
        gid = row['group']
        se, me = cfg['reference_norm_counts'][gid]
        assert float(row['SE_norm']) == se and float(row['ME_norm']) == me
        scale = se / me
        assert math.isclose(float(row['ME_scale']), scale, rel_tol=cfg['numeric_rtol'])
        assert math.isclose(float(row['scaled_ME_norm']), se, rel_tol=cfg['numeric_rtol'])
        assert float(row['QA_max_abs_difference']) == cfg['regular_bin_closure_atol']
        native = [r for r in spectra if r['group'] == gid and r['binning'] == 'native']
        display = [r for r in spectra if r['group'] == gid and r['binning'] == 'display']
        n, d = len(native) - 2, len(display) - 2
        assert n == d * cfg['rebin_factor']
        for kind in [native, display]:
            for r in kind:
                assert math.isclose(float(r['scaled_ME']), float(r['ME']) * scale,
                                    rel_tol=cfg['numeric_rtol'], abs_tol=1e-10)
                assert math.isclose(float(r['scaled_ME_error']), float(r['ME_error']) * scale,
                                    rel_tol=cfg['numeric_rtol'], abs_tol=1e-10)
            band = [r for r in kind[1:-1] if float(r['low']) >= cfg['normalization'][0]-1e-10
                    and float(r['high']) <= cfg['normalization'][1]+1e-10]
            assert math.fsum(float(r['SE']) for r in band) == se
            assert math.fsum(float(r['ME']) for r in band) == me
            assert math.isclose(math.fsum(float(r['scaled_ME']) for r in band), se,
                                rel_tol=cfg['numeric_rtol'])
        for key in ['SE', 'ME']:
            assert math.fsum(float(r[key]) for r in native[1:-1]) == math.fsum(float(r[key]) for r in display[1:-1])
            for i, dr in enumerate(display[1:-1]):
                nr = native[1 + i * cfg['rebin_factor']:1 + (i + 1) * cfg['rebin_factor']]
                assert math.fsum(float(r[key]) for r in nr) == float(dr[key])
                assert math.isclose(math.fsum(float(r[key + '_error']) ** 2 for r in nr),
                                    float(dr[key + '_error']) ** 2, rel_tol=cfg['numeric_rtol'])
        checks.append(dict(group=gid, SE_norm=se, ME_norm=me, ME_scale=scale,
                           native_bins=n, display_bins=d, normalization_closure='PASS',
                           rebin_count_and_error_closure='PASS', raw_vs_prior_cf_projection='EXACT'))
        for suffix in ['.png', '.pdf']:
            file = out / ('se_me_' + gid + suffix)
            assert file.stat().st_size > 1000
            assert file.read_bytes().startswith(b'\x89PNG' if suffix == '.png' else b'%PDF-')
    pdf = out / 'se_me_overlay.pdf'
    pdf_data = pdf.read_bytes()
    pages = len(re.findall(rb'/Type\s*/Page(?:\s|/|>)', pdf_data))
    assert pages == len(cfg['groups'])
    assert b'%%EOF' in pdf_data[-1024:]
    for file in [pdf] + [out / ('se_me_' + g['id'] + '.pdf') for g in cfg['groups']]:
        result = subprocess.run(['/usr/bin/gs', '-dSAFER', '-dBATCH', '-dNOPAUSE',
                                 '-sDEVICE=nullpage', str(file)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        (out / (file.stem + '_pdf_validation.log')).write_text(result.stdout)
        assert result.returncode == 0, 'Ghostscript could not read PDF: ' + str(file)
        expected_pages = pages if file == pdf else 1
        assert len(re.findall(rb'/Type\s*/Page(?:\s|/|>)', file.read_bytes())) == expected_pages
    assert (out / 'overlay_histograms.root').read_bytes()[:4] == b'root'
    validation = dict(status='PASS', created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      input_before=before, input_after=after, input_unchanged=True,
                      root_success_marker_verified=True, validation_existing=validate_existing,
                      groups=checks, pdf_pages=pages, pdf_parser='Ghostscript nullpage + PDF Page object count',
                      production_sources_modified=False, jobs_submitted=False, mail_sent=False)
    (out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    outputs = [out / 'settings.json', out / 'normalization.tsv', out / 'spectra.tsv',
               out / 'overlay_histograms.root', out / 'se_me_overlay.pdf',
               out / 'source' / 'DrawSEME.C', out / 'source' / 'run_overlay.py']
    outputs += [out / ('se_me_' + g['id'] + s) for g in cfg['groups'] for s in ['.png', '.pdf']]
    validation['sha256'] = {str(f.relative_to(out)): sha256(f) for f in outputs}
    (out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('settings', type=Path)
    parser.add_argument('--validate-existing', action='store_true')
    args = parser.parse_args()
    run(args.settings, args.validate_existing)
