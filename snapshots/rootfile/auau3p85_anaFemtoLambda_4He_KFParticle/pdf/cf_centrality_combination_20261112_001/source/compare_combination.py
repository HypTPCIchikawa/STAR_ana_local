#!/usr/bin/env python3
"""Compare explicit CF centrality estimands using saved native counts only."""
import argparse
import csv
import datetime
import hashlib
import json
import math
import mmap
from pathlib import Path


def fp(path):
    st = path.stat()
    return dict(size=st.st_size, mtime_ns=st.st_mtime_ns, inode=st.st_ino, device=st.st_dev)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def extract_histogram(mm, name):
    marker = json.dumps(name).encode()
    offset = 0
    while True:
        found = mm.find(marker, offset)
        if found < 0:
            raise ValueError('Missing histogram: ' + name)
        q = found + len(marker)
        while mm[q:q+1].isspace():
            q += 1
        if mm[q:q+1] == b':':
            break
        offset = q
    q += 1
    while mm[q:q+1].isspace():
        q += 1
    assert mm[q:q+1] == b'{'
    start, depth, quoted, escape = q, 0, False, False
    for j in range(start, len(mm)):
        c = mm[j]
        if quoted:
            if escape:
                escape = False
            elif c == 92:
                escape = True
            elif c == 34:
                quoted = False
        else:
            if c == 34:
                quoted = True
            elif c == 123:
                depth += 1
            elif c == 125:
                depth -= 1
                if depth == 0:
                    return json.loads(mm[start:j+1])
    raise ValueError('Unterminated histogram object')


def close(a, b, cfg):
    return math.isclose(a, b, rel_tol=cfg['numeric_rtol'], abs_tol=cfg['numeric_atol'])


def edge_index(edges, value, cfg):
    i = min(range(len(edges)), key=lambda k: abs(edges[k] - value))
    assert close(edges[i], value, cfg), 'Unaligned k* boundary'
    return i


def indices(edges, bounds, cfg):
    return list(range(edge_index(edges, bounds[0], cfg), edge_index(edges, bounds[1], cfg)))


def total(data, side, cents, bins, key='content'):
    return math.fsum(data[c][side][key][i] for c in cents for i in bins)


def blank_gradient(data):
    return {side: {c: [0.] * len(data[c][side]['content']) for c in data} for side in ['SE', 'ME']}


def variance(gradient, data):
    return math.fsum(gradient[side][c][i] ** 2 * data[c][side]['error2'][i]
                     for side in gradient for c in data for i in range(len(gradient[side][c])))


def aggregate(data, groups, norm, target, fixed_weights=False):
    """Return inclusive ratio or fixed-normalization-SE-weight average and Jacobian."""
    cents = list(data)
    sn = total(data, 'SE', cents, norm)
    sw = total(data, 'SE', cents, target)
    pieces = []
    for cs in groups:
        sh, mh = total(data, 'SE', cs, norm), total(data, 'ME', cs, norm)
        sx, mx = total(data, 'SE', cs, target), total(data, 'ME', cs, target)
        assert sh > 0 and mh > 0
        pieces.append((cs, sh, mh, sx, mx))
    gradient = blank_gradient(data)
    if fixed_weights:
        if any(sx <= 0 or mx <= 0 for _, _, _, sx, mx in pieces):
            return dict(value=None, error=None, reason='undefined individual CF'), None
        z = math.fsum(sx * mh / mx for _, _, mh, sx, mx in pieces)
        value = z / sn
        for cs, sh, mh, sx, mx in pieces:
            for c in cs:
                for i in target:
                    gradient['SE'][c][i] += mh / (mx * sn)
                    gradient['ME'][c][i] -= sx * mh / (mx * mx * sn)
                for i in norm:
                    gradient['SE'][c][i] -= z / (sn * sn)
                    gradient['ME'][c][i] += sx / (mx * sn)
        weights = [{'cent9': cs, 'weight': sh / sn} for cs, sh, _, _, _ in pieces]
        result = dict(value=value, reference_denominator=None, numerator=sw,
                      kind='fixed normalization-SE-pair-weight mean of individual normalized CFs', weights=weights)
    else:
        d = math.fsum(sh / mh * mx for _, sh, mh, _, mx in pieces)
        if sw <= 0 or d <= 0:
            return dict(value=None, error=None, reason='nonpositive target count'), None
        value = sw / d
        for cs, sh, mh, sx, mx in pieces:
            for c in cs:
                for i in target:
                    gradient['SE'][c][i] += 1. / d
                    gradient['ME'][c][i] -= sw * sh / (mh * d * d)
                for i in norm:
                    gradient['SE'][c][i] -= sw * mx / (mh * d * d)
                    gradient['ME'][c][i] += sw * sh * mx / (mh * mh * d * d)
        weights = [{'cent9': cs, 'weight': sh / mh * mx / d} for cs, sh, mh, _, mx in pieces]
        result = dict(value=value, reference_denominator=d, numerator=sw,
                      kind='sum SE / sum individually normalization-scaled ME', weights=weights)
    result['error'] = math.sqrt(max(0., variance(gradient, data)))
    return result, gradient


def difference_error(a, b, data):
    if a is None or b is None:
        return None
    diff = {side: {c: [x-y for x, y in zip(a[side][c], b[side][c])]
                   for c in data} for side in a}
    return math.sqrt(max(0., variance(diff, data)))


def evaluate(data, cfg, norm, target):
    output, gradients = {}, {}
    for method, groups in cfg['partitions'].items():
        output[method], gradients[method] = aggregate(data, groups, norm, target,
                                                    method in cfg['fixed_weight_methods'])
    for method in output:
        if output[method]['value'] is not None and output['global']['value'] is not None:
            output[method]['minus_global'] = output[method]['value'] - output['global']['value']
            output[method]['difference_count_model_error'] = difference_error(gradients[method], gradients['global'], data)
    return output


def write_table(path, rows, fields):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main(settings_path):
    settings_path = settings_path.resolve()
    cfg = json.loads(settings_path.read_text())
    assert cfg['fixed_weight_extra_global_renormalization'] is False
    out = settings_path.parent
    assert not (out / 'results.json').exists(), 'Refuse existing results'
    input_path = Path(cfg['input_export'])
    before = fp(input_path)
    input_sha = digest(input_path)
    assert input_sha == cfg['expected_input_sha256']
    data, edges = {}, None
    with input_path.open('rb') as stream:
        with mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            for cent in cfg['fine_cent9']:
                data[cent] = {}
                for side, template in [('SE', cfg['se_template']), ('ME', cfg['me_template'])]:
                    h = extract_histogram(mm, template % cent)
                    ny, nx = h['shape']
                    assert nx == cfg['native_axis']['bins'] + 2
                    if edges is None:
                        edges = h['xEdges']
                    assert edges == h['xEdges']
                    data[cent][side] = {key: [math.fsum(h[key][y*nx+x] for y in range(ny))
                                                    for x in range(1, nx-1)]
                                        for key in ['content', 'error2']}
    assert close(edges[0], cfg['native_axis']['range'][0], cfg)
    assert close(edges[-1], cfg['native_axis']['range'][1], cfg)
    assert all(close(edges[i+1]-edges[i], cfg['native_axis']['width'], cfg) for i in range(len(edges)-1))
    norm = indices(edges, cfg['normalization'], cfg)
    counts, fine_norms = [], []
    for cent in cfg['fine_cent9']:
        sh = total(data, 'SE', [cent], norm)
        mh = total(data, 'ME', [cent], norm)
        assert [sh, mh] == cfg['expected_norm_counts'][str(cent)]
        fine_norms.append(dict(cent9=cent, label=cfg['fine_labels'][str(cent)], SE_norm=sh, ME_norm=mh,
                               alpha=sh/mh, ME_SE_exposure=mh/sh))
        for i in range(len(edges)-1):
            counts.append(dict(cent9=cent, bin=i+1, low=edges[i], high=edges[i+1],
                               SE=data[cent]['SE']['content'][i], ME=data[cent]['ME']['content'][i],
                               SE_variance=data[cent]['SE']['error2'][i], ME_variance=data[cent]['ME']['error2'][i]))
    windows, window_rows = {}, []
    for name, bounds in cfg['windows'].items():
        answer = evaluate(data, cfg, norm, indices(edges, bounds, cfg))
        assert close(answer['global']['value'], cfg['expected_global_windows'][name], cfg)
        windows[name] = dict(bounds=bounds, methods=answer)
        for method, values in answer.items():
            window_rows.append(dict(window=name, low=bounds[0], high=bounds[1], method=method, **values))
    factor = cfg['rebin_factor']
    assert (len(edges)-1) % factor == 0
    assert edge_index(edges, cfg['normalization'][0], cfg) % factor == 0
    assert edge_index(edges, cfg['normalization'][1], cfg) % factor == 0
    points, point_rows = [], []
    for i in range(0, len(edges)-1, factor):
        answer = evaluate(data, cfg, norm, list(range(i, i+factor)))
        points.append(dict(low=edges[i], high=edges[i+factor], methods=answer))
        for method, values in answer.items():
            point_rows.append(dict(low=edges[i], high=edges[i+factor], method=method, **values))
    closure = evaluate(data, cfg, norm, norm)
    assert all(close(row['value'], 1., cfg) and close(row['error'], 0., cfg) for row in closure.values())
    common = json.loads(json.dumps(data))
    common = {int(c): sides for c, sides in common.items()}
    for c in common:
        scale = total(data, 'SE', [c], norm) / total(data, 'ME', [c], norm)
        common[c]['ME']['content'] = [x*scale for x in data[c]['ME']['content']]
        common[c]['ME']['error2'] = [x*scale*scale for x in data[c]['ME']['error2']]
    common_checks = []
    for name, bounds in cfg['windows'].items():
        answer = evaluate(common, cfg, norm, indices(edges, bounds, cfg))
        values = [answer[m]['value'] for m in cfg['inclusive_methods']]
        assert all(close(v, values[0], cfg) for v in values)
        common_checks.append(dict(window=name, inclusive_values=values))
    toy = cfg['toy']
    toy_data = {int(c): {'SE': {'content': s, 'error2': s},
                         'ME': {'content': [x*toy['ME_exposures'][c] for x in s],
                                'error2': [x*toy['ME_exposures'][c] for x in s]}}
                for c, s in toy['SE'].items()}
    toy_results = []
    for target in toy['targets']:
        global_toy, _ = aggregate(toy_data, [list(toy_data)], toy['normalization_indices'], target)
        balanced_toy, _ = aggregate(toy_data, [[c] for c in toy_data], toy['normalization_indices'], target)
        fixed_toy, _ = aggregate(toy_data, [[c] for c in toy_data], toy['normalization_indices'], target, True)
        assert close(balanced_toy['value'], 1., cfg) and close(fixed_toy['value'], 1., cfg)
        assert not close(global_toy['value'], 1., cfg)
        toy_results.append(dict(target=target, global_ratio=global_toy['value'], balanced=balanced_toy['value'], fixed_weight_mean=fixed_toy['value']))
    after = fp(input_path)
    assert before == after
    results = dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   fine_norms=fine_norms, windows=windows, points=points,
                   interpretation={'global':'raw counts summed across cent9, one common normalization, then ratio',
                                   'balanced':'each partition ME scaled to its SE normalization yield, then sum counts and divide; k-dependent scaled-ME weights of individual CFs',
                                   'fixed_norm_SE_weights':'weighted average of individual normalized CFs using constant SE normalization-yield fractions; a different estimand from inclusive sum-ratio'},
                   uncertainty=cfg['error_model'])
    results['fixed_weight_extra_global_renormalization'] = False
    norm_map = {row['cent9']: row for row in fine_norms}
    results['partition_norms'] = {
        method: [dict(cent9=g,
                      SE_norm=sum(norm_map[c]['SE_norm'] for c in g),
                      ME_norm=sum(norm_map[c]['ME_norm'] for c in g),
                      alpha=sum(norm_map[c]['SE_norm'] for c in g) / sum(norm_map[c]['ME_norm'] for c in g))
                 for g in groups]
        for method, groups in cfg['partitions'].items()}
    (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    write_table(out/'counts.tsv', counts, ['cent9','bin','low','high','SE','ME','SE_variance','ME_variance'])
    write_table(out/'normalization.tsv', fine_norms, ['cent9','label','SE_norm','ME_norm','alpha','ME_SE_exposure'])
    fields = ['method','value','error','minus_global','difference_count_model_error','numerator','reference_denominator']
    write_table(out/'window_results.tsv', window_rows, ['window','low','high']+fields)
    write_table(out/'display_points.tsv', point_rows, ['low','high']+fields)
    validation = dict(status='PASS', input_before=before, input_after=after, input_unchanged=True,
                      input_sha256=input_sha, selected_histograms=2*len(data), native_bins_per_slice=len(edges)-1,
                      displayed_bins=len(points), normalization_closure='PASS for H treated as a single window target in all five methods; not a constraint on fixed-weight point means under the inclusive D(k) measure',
                      fixed_weight_extra_global_renormalization=False,
                      common_alpha_identity='PASS for the three inclusive sum-ratio methods; fixed-weight averages are a different observable',
                      common_alpha_checks=common_checks, toy_individual_CF_unity_but_raw_combination_distorted=toy_results,
                      window_global_matches_previous_study=True, covariance_of_differences='Shared input samples included through difference of Jacobians',
                      source_files_modified=False, root_opened=False, jobs_submitted=False)
    validation['fixed_weight_point_norm_checks'] = {}
    for fixed, balanced in [('fixed_norm_SE_weights2', 'balanced2'), ('fixed_norm_SE_weights3', 'balanced3')]:
        band = [row for row in points
                if row['low'] >= cfg['normalization'][0]-1e-10
                and row['high'] <= cfg['normalization'][1]+1e-10]
        mean = math.fsum(row['methods'][fixed]['value'] * row['methods'][balanced]['reference_denominator'] for row in band) / math.fsum(row['methods'][balanced]['reference_denominator'] for row in band)
        validation['fixed_weight_point_norm_checks'][fixed] = dict(
            D_weighted_point_mean_in_H=mean, renormalized=False,
            meaning='Reported without requiring unity; not the same as first forming one H-window CF in each class')
    output_paths = [settings_path, Path(__file__).resolve(), out/'results.json', out/'counts.tsv',
                    out/'normalization.tsv', out/'window_results.tsv', out/'display_points.tsv']
    validation['sha256'] = {str(p.relative_to(out)): digest(p) for p in output_paths}
    (out/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
    print(json.dumps(dict(study=str(out),validation='PASS',fine_norms=fine_norms,windows=windows),indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('settings', type=Path)
    main(parser.parse_args().settings)
