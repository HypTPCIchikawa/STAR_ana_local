"""Full-mass-window Method3/5/6 CFs with shared-input delta-method errors.

Method6 uses the Method5 *ME* purity as one shared random quantity on both
SE and ME.  Purity is applied to native bins before CF rebinning.  No ROOT
objects are read here and no analysis settings have implicit defaults.
"""
import math

from compare_mass_halves import fit_inputs, selected_bins
from mass_half_statistics import ErrorModel, Quantity, integral


METHODS = ('method3', 'method5', 'method6')
SIDES = ('SE', 'ME')
LIMITATIONS = [
    'First-order delta-method errors use supplied Sumw2 and background-fit covariance.',
    'Fit-background/count covariance and covariance between distinct background fits are unavailable and omitted.',
    'Overlapping purity-window count covariance and shared SE/ME Method6 purity are retained.',
    'Event-level and mixed-pair correlations absent from the input histograms cannot be recovered.',
    'Cross-method covariance and significance of method differences are not estimated.',
    'ME-to-SE purity-transfer and mass-model systematic uncertainty are not included in the copied ME statistical error.',
]


def _estimate(model, quantity):
    return model.estimate(quantity) if quantity is not None else None


def _matrix(model, first, second=None):
    second = first if second is None else second
    return [[model.covariance(a, b) if a is not None and b is not None else None
             for b in second] for a in first]


def _reuse(model, key, value):
    """Refer to an input already registered by fit_inputs, never duplicate it."""
    if key not in model.variances:
        raise ValueError('Missing registered fit input: %r' % (key,))
    return Quantity(model, value, {key: 1.0})


def _reason(row):
    return (row.get('full', {}).get('reason') or row.get('reason') or
            '; '.join(str(reason) for reason in row.get('regionReason', []) if reason) or
            'full-window mass fit unavailable')


def _physical_purity(purity):
    return purity is not None and 0.0 < purity.value <= 1.0


def _validate_edges(group, settings):
    edges = group['kstarEdges']
    original = group['originalKstarEdges']
    coarse = group['method5PurityEdges']
    for name, values in [('CF', edges), ('native', original), ('purity', coarse)]:
        if (len(values) < 2 or any(not math.isfinite(float(value)) for value in values) or
                any(b <= a for a, b in zip(values, values[1:]))):
            raise ValueError('%s edges must be finite and strictly increasing.' % name)
    for key in ('rebin', 'purityRebin'):
        if isinstance(settings[key], bool) or not isinstance(settings[key], int) or settings[key] <= 0:
            raise ValueError('%s must be a positive integer.' % key)
    cf_rebin, purity_rebin = settings['rebin'], settings['purityRebin']
    n_native = len(original) - 1
    if n_native % cf_rebin or n_native % purity_rebin:
        raise ValueError('CF/purity binning must divide the native bins.')
    for name, actual, expected in [('CF', edges, original[::cf_rebin]),
                                   ('purity', coarse, original[::purity_rebin])]:
        if len(actual) != len(expected) or any(abs(a-b) > 1.e-9 for a, b in zip(actual, expected)):
            raise ValueError('%s edges disagree with the configured native-bin grouping.' % name)
    for side in SIDES:
        if len(group['method5Base'][side]) != n_native:
            raise ValueError('Method5 native-count size differs from its edges.')
        if len(group['methods']['method3'][side]) != len(edges)-1:
            raise ValueError('Method3 fit-row size differs from CF edges.')
        rows = group['methods']['method5'][side]
        if len(rows) != len(coarse)-1:
            raise ValueError('Method5 fit-row size differs from purity edges.')
        for row in rows:
            first, last = row['windowOriginalFirstBin'], row['windowOriginalLastBin']
            if (isinstance(first, bool) or isinstance(last, bool) or
                    not isinstance(first, int) or not isinstance(last, int) or
                    not 1 <= first <= last <= n_native):
                raise ValueError('Method5 purity window is outside the native-bin range.')


def _method3_purities(model, group, source):
    output, reasons = {}, []
    for side in SIDES:
        output[side] = []
        for index, row in enumerate(group['methods']['method3'][side]):
            count = sum((_reuse(model, ('method3', side, 'count', index, half),
                                row['centralNtot'][half]) for half in range(2)),
                        model.constant(0))
            signal = source[side][index]
            purity = signal/count if signal is not None and count.value > 0.0 else None
            if not _physical_purity(purity):
                source[side][index] = None
                purity = None
                reasons.append(dict(sample=side, bin=index+1, reason=_reason(row)))
            output[side].append(purity)
    return output, reasons


def _method5_purities(model, group):
    output, reasons = {}, []
    for side in SIDES:
        cells = [[_reuse(model, ('method5', side, 'base', index, half), row['counts'][half])
                  for half in range(2)]
                 for index, row in enumerate(group['method5Base'][side])]
        output[side] = []
        for index, row in enumerate(group['methods']['method5'][side]):
            prefix = ('method5', side, 'background', index)
            purity = None
            if (row['full']['valid'] and row['background'] is not None and
                    all((prefix, half) in model.variances for half in range(2))):
                background = sum((_reuse(model, (prefix, half), row['background'][half])
                                  for half in range(2)), model.constant(0))
                first, last = row['windowOriginalFirstBin']-1, row['windowOriginalLastBin']
                total = sum((sum(cells[i], model.constant(0)) for i in range(first, last)),
                            model.constant(0))
                if total.value > 0.0:
                    purity = 1-background/total
            if not _physical_purity(purity):
                purity = None
                reasons.append(dict(sample=side, purity_bin=index+1, reason=_reason(row)))
            output[side].append(purity)
    return output, reasons


def _method6_inputs(model, group, settings):
    """Build one ME-derived purity, with native ME-count derivative information."""
    cells = {}
    for side in SIDES:
        cells[side] = [[model.variable(('method6', side, 'base', index, half),
                                     row['counts'][half], row['error2'][half])
                        for half in range(2)]
                       for index, row in enumerate(group['method5Base'][side])]
    purities, reasons = [], []
    for index, row in enumerate(group['methods']['method5']['ME']):
        purity, reason = None, _reason(row)
        if row['full']['valid'] and row['background'] is not None and row['backgroundCovariance'] is not None:
            try:
                background = model.correlated_block(('method6', 'ME', 'background', index),
                                                    row['background'], row['backgroundCovariance'])
                first, last = row['windowOriginalFirstBin']-1, row['windowOriginalLastBin']
                total = sum((sum(cells['ME'][i], model.constant(0)) for i in range(first, last)),
                            model.constant(0))
                if total.value > 0.0:
                    purity = 1-sum(background, model.constant(0))/total
                else:
                    reason = 'nonpositive native ME count in purity window'
            except (TypeError, ValueError) as exc:
                reason = 'invalid ME background fit/covariance: ' + str(exc)
        if not _physical_purity(purity):
            if purity is not None:
                reason = 'ME purity is outside (0, 1]'
            purity = None
            reasons.append(dict(sample='ME', appliesTo=['SE', 'ME'], purity_bin=index+1, reason=reason))
        purities.append(purity)
    source = {side: [] for side in SIDES}
    cf_rebin, purity_rebin = settings['rebin'], settings['purityRebin']
    for first in range(0, len(cells['ME']), cf_rebin):
        indices = range(first, first+cf_rebin)
        # A failed native contribution masks the entire output bin, even when
        # the surviving contribution across a 50-MeV boundary would be positive.
        if any(purities[index//purity_rebin] is None for index in indices):
            for side in SIDES:
                source[side].append(None)
            continue
        for side in SIDES:
            source[side].append(sum((sum(cells[side][index], model.constant(0)) *
                                     purities[index//purity_rebin] for index in indices),
                                    model.constant(0)))
    # The same objects, including their ME-input gradients, are intentionally reused.
    return source, {'SE': purities, 'ME': purities}, reasons


def _mask_method5(source, purities, settings):
    """Also reject nonphysical coarse P and preserve partial-bin invalidity."""
    for side in SIDES:
        for index in range(len(source[side])):
            first = index*settings['rebin']
            if any(purities[side][native//settings['purityRebin']] is None
                   for native in range(first, first+settings['rebin'])):
                source[side][index] = None


def _normalization(model, source, indices):
    problems = []
    for index in indices:
        for side in SIDES:
            quantity = source[side][index]
            if quantity is None:
                problems.append('%s bin %d: invalid contributing mass fit' % (side, index+1))
            elif quantity.value <= 0.0:
                problems.append('%s bin %d: nonpositive corrected yield' % (side, index+1))
            else:
                variance = model.variance(quantity)
                if not math.isfinite(variance) or variance <= 0.0:
                    problems.append('%s bin %d: nonfinite or nonpositive corrected-yield variance'
                                    % (side, index+1))
    if problems:
        return None, dict(valid=False, reason='; '.join(problems), problems=problems,
                          factor=None, numerator=None, denominator=None, bins=[i+1 for i in indices])
    numerator, denominator = (integral(source[side], indices) for side in SIDES)
    if numerator.value <= 0.0 or denominator.value <= 0.0:
        return None, dict(valid=False, reason='nonpositive normalization integral', problems=[],
                          factor=None, numerator=_estimate(model, numerator),
                          denominator=_estimate(model, denominator), bins=[i+1 for i in indices])
    factor = denominator/numerator
    return factor, dict(valid=True, reason=None, problems=[], factor=model.estimate(factor),
                        numerator=model.estimate(numerator), denominator=model.estimate(denominator),
                        bins=[i+1 for i in indices])


def _analyze_method(model, method, source, purities, purity_edges, reasons,
                    edges, norm_indices, plot_indices, comparison_indices):
    factor, normalization = _normalization(model, source, norm_indices)
    quantities, points = [], []
    for index, (same, mixed) in enumerate(zip(source['SE'], source['ME'])):
        quantity, ratio, reason = None, None, None
        if same is None or mixed is None:
            reason = 'invalid contributing mass fit or purity'
        elif mixed.value <= 0.0:
            reason = 'nonpositive mixed-event denominator'
        else:
            ratio = same/mixed
        if factor is None:
            reason = 'normalization unavailable: ' + normalization['reason']
        elif ratio is not None:
            quantity = factor*ratio
        quantities.append(quantity)
        points.append(dict(bin=index+1, low=edges[index], high=edges[index+1],
                           cf=_estimate(model, quantity), yieldRatio=_estimate(model, ratio), reason=reason))
    integrated, integrated_reason = None, None
    if factor is None:
        integrated_reason = 'normalization unavailable: ' + normalization['reason']
    elif any(source[side][i] is None for side in SIDES for i in comparison_indices):
        integrated_reason = 'invalid contributing mass fit or purity in comparison interval'
    else:
        same, mixed = (integral(source[side], comparison_indices) for side in SIDES)
        if mixed.value > 0.0:
            integrated = factor*same/mixed
        else:
            integrated_reason = 'nonpositive integrated mixed-event denominator'
    return dict(
        points=points, covariance=_matrix(model, quantities), normalization=normalization,
        yields=dict(SE=[_estimate(model, q) for q in source['SE']],
                    ME=[_estimate(model, q) for q in source['ME']],
                    covariance={'SE': _matrix(model, source['SE']),
                                'ME': _matrix(model, source['ME']),
                                'SE_ME': _matrix(model, source['SE'], source['ME'])}),
        purity=dict(edges=list(purity_edges), SE=[_estimate(model, p) for p in purities['SE']],
                    ME=[_estimate(model, p) for p in purities['ME']],
                    shared=(method == 'method6'), identicalSharedP=(method == 'method6'),
                    source='ME' if method == 'method6' else 'independent SE and ME fits',
                    covariance={'SE': _matrix(model, purities['SE']),
                                'ME': _matrix(model, purities['ME']),
                                'SE_ME': _matrix(model, purities['SE'], purities['ME'])}),
        integrated=dict(cf=_estimate(model, integrated), reason=integrated_reason),
        fit_problems=reasons,
        valid_bins=sum(q is not None for q in quantities),
        valid_plot_bins=sum(quantities[i] is not None for i in plot_indices),
        total_plot_bins=len(plot_indices))


def analyze_group(group, settings, normalization, plot_range, comparison_range):
    """Analyze one saved mass_half_extract group, full mass window only.

    ``settings`` is the extraction JSON's ``settings`` object (not the plotting
    YAML).  All requested limits must coincide with CF-bin edges.  Returned
    indices and point bins are one-based; arrays retain all saved CF bins.
    Missing CFs/covariances are JSON null, never a numerical zero or NaN.
    """
    _validate_edges(group, settings)
    edges = group['kstarEdges']
    ranges = [normalization, plot_range, comparison_range]
    for limits in ranges:
        if (len(limits) != 2 or not all(math.isfinite(float(x)) for x in limits) or
                limits[0] >= limits[1]):
            raise ValueError('Every analysis range must contain two finite increasing limits.')
    norm_indices, plot_indices, comparison_indices = [selected_bins(edges, limits) for limits in ranges]
    result = dict(schemaVersion=1, id=group['id'], label=group.get('label', group['id']),
                  kstarEdges=list(edges), originalKstarEdges=list(group['originalKstarEdges']),
                  method5PurityEdges=list(group['method5PurityEdges']),
                  normalization=list(normalization), plot_range=list(plot_range),
                  comparison_range=list(comparison_range),
                  normalization_bins=[i+1 for i in norm_indices],
                  plot_bins=[i+1 for i in plot_indices],
                  comparison_bins=[i+1 for i in comparison_indices],
                  limitations=list(LIMITATIONS), methods={})
    for method in METHODS:
        model = ErrorModel()
        if method == 'method6':
            source, purities, reasons = _method6_inputs(model, group, settings)
        else:
            fitted, reasons = fit_inputs(model, group, method, settings)
            source = {side: fitted[side]['full'] for side in SIDES}
            if method == 'method3':
                purities, extra_reasons = _method3_purities(model, group, source)
            else:
                purities, extra_reasons = _method5_purities(model, group)
                _mask_method5(source, purities, settings)
            reasons += extra_reasons
        result['methods'][method] = _analyze_method(
            model, method, source, purities,
            edges if method == 'method3' else group['method5PurityEdges'], reasons,
            edges, norm_indices, plot_indices, comparison_indices)
    return result
