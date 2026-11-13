"""Statistical contracts for shared-ME-purity Method6; no ROOT execution."""
import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from method6_statistics import analyze_group


SETTINGS = {'rebin': 4, 'purityRebin': 10}


def make_group(se_counts=None, me_counts=None, me_purities=None, background_variance=4.0):
    se_counts = se_counts or [20.0]*40
    me_counts = me_counts or [10.0]*40
    me_purities = me_purities or [.8]*4
    native = [index*.005 for index in range(41)]
    group = dict(id='toy', label='Toy full window', originalKstarEdges=native,
                 kstarEdges=native[::4], method5PurityEdges=native[::10],
                 method5Base={}, methods={'method3': {}, 'method5': {}})
    covariance = [[background_variance, .2*background_variance],
                  [.2*background_variance, background_variance]]
    for side, values in [('SE', se_counts), ('ME', me_counts)]:
        group['method5Base'][side] = [dict(counts=[.4*value, .6*value],
                                               error2=[.4*value, .6*value]) for value in values]
        group['methods']['method3'][side] = []
        group['methods']['method5'][side] = []
        for first in range(0, len(values), 4):
            total = sum(values[first:first+4])
            background = .2*total
            group['methods']['method3'][side].append(dict(
                centralNtot=[.4*total, .6*total], centralError2=[.4*total, .6*total],
                background=[.4*background, .6*background],
                backgroundCovariance=copy.deepcopy(covariance), valid=[True, True],
                regionReason=['ok', 'ok'], full={'valid': True, 'reason': 'ok'}))
        for index, first in enumerate(range(0, len(values), 10)):
            total = sum(values[first:first+10])
            purity = me_purities[index] if side == 'ME' else .6
            background = total*(1-purity)
            group['methods']['method5'][side].append(dict(
                windowOriginalFirstBin=first+1, windowOriginalLastBin=first+10,
                background=[.4*background, .6*background],
                backgroundCovariance=copy.deepcopy(covariance), valid=[True, True],
                regionReason=['ok', 'ok'], full={'valid': True, 'reason': 'ok'}))
    return group


def analyze(group, normalization=(.12, .16), comparison=(0., .04)):
    return analyze_group(group, SETTINGS, normalization, (0., .2), comparison)


def invalidate(group, side, coarse_bin):
    row = group['methods']['method5'][side][coarse_bin]
    row['full'] = {'valid': False, 'reason': 'invalid_fit_or_covariance'}
    row['valid'] = [False, False]
    row['regionReason'] = ['invalid_fit_or_covariance']*2
    row['background'] = None
    row['backgroundCovariance'] = None


class Method6StatisticsTests(unittest.TestCase):
    def test_shared_purity_is_identical_on_both_sides_with_shared_covariance(self):
        group = make_group()
        result = analyze(group)['methods']
        purity = result['method6']['purity']
        self.assertTrue(purity['shared'])
        self.assertTrue(purity['identicalSharedP'])
        self.assertEqual(purity['source'], 'ME')
        self.assertEqual(purity['SE'], purity['ME'])
        self.assertEqual(purity['covariance']['SE'], purity['covariance']['SE_ME'])
        self.assertAlmostEqual(purity['ME'][0]['value'], .8)
        # P=1-B/N: VarP=VarB/N^2+B^2 VarN/N^4, using native ME counts.
        n, b = 100.0, 20.0
        background_variance = 4.0+4.0+2*.8
        expected = background_variance/n**2 + b**2*n/n**4
        self.assertAlmostEqual(purity['ME'][0]['variance'], expected)
        self.assertAlmostEqual(result['method5']['purity']['ME'][0]['variance'], expected)
        self.assertFalse(result['method5']['purity']['shared'])
        self.assertAlmostEqual(result['method5']['purity']['SE'][0]['value'], .6)

    def test_common_purity_error_cancels_in_cell_ratio_and_common_normalization(self):
        first = analyze(make_group(background_variance=1.0))['methods']['method6']
        second = analyze(make_group(background_variance=1.e6))['methods']['method6']
        for one, two in zip(first['points'], second['points']):
            self.assertAlmostEqual(one['yieldRatio']['value'], 2.0)
            self.assertAlmostEqual(one['yieldRatio']['variance'], two['yieldRatio']['variance'], places=10)
            self.assertAlmostEqual(one['cf']['value'], 1.0)
            self.assertAlmostEqual(one['cf']['variance'], two['cf']['variance'], places=10)

    def test_fifty_mev_boundary_is_weighted_on_native_counts_before_twenty_mev_sum(self):
        se = [5.0]*10 + [20.0]*30
        me = [10.0]*40
        result = analyze(make_group(se, me, [.8, .6, .7, .9]))['methods']['method6']
        # Bin .04-.06 has two native cells on either side of the .05 boundary.
        self.assertAlmostEqual(result['yields']['SE'][2]['value'], 2*5*.8 + 2*20*.6)
        self.assertAlmostEqual(result['yields']['ME'][2]['value'], 2*10*.8 + 2*10*.6)
        self.assertAlmostEqual(result['points'][2]['yieldRatio']['value'], 32.0/28.0)
        self.assertNotAlmostEqual(result['points'][2]['yieldRatio']['value'], 50.0/40.0)

    def test_me_count_and_purity_covariance_recovers_n_minus_background_variance(self):
        result = analyze(make_group())['methods']['method6']
        # The first five CF bins span two complete coarse purity cells.
        value = sum(row['value'] for row in result['yields']['ME'][:5])
        covariance = result['yields']['covariance']['ME']
        variance = sum(covariance[i][j] for i in range(5) for j in range(5))
        self.assertAlmostEqual(value, 200.0-40.0)
        self.assertAlmostEqual(variance, 200.0+2*(4.0+4.0+2*.8))
        self.assertGreater(result['yields']['covariance']['SE_ME'][0][0], 0.0)

    def test_bad_se_fit_does_not_invalidate_me_derived_method6(self):
        group = make_group()
        baseline = analyze(group)['methods']['method6']
        for index in range(4):
            invalidate(group, 'SE', index)
        changed = analyze(group)['methods']
        self.assertEqual(changed['method6'], baseline)
        self.assertFalse(changed['method5']['normalization']['valid'])

    def test_failed_me_purity_masks_both_sides_including_partial_boundary_bin(self):
        group = make_group()
        invalidate(group, 'ME', 1)
        result = analyze(group, normalization=(.16, .2), comparison=(0., .08))['methods']['method6']
        self.assertTrue(result['normalization']['valid'])
        for index in (2, 3, 4):
            self.assertIsNone(result['yields']['SE'][index])
            self.assertIsNone(result['yields']['ME'][index])
            self.assertIsNone(result['points'][index]['cf'])
        self.assertIsNotNone(result['points'][1]['cf'])
        self.assertIsNone(result['integrated']['cf'])
        self.assertIn('comparison interval', result['integrated']['reason'])
        self.assertIsNone(result['purity']['SE'][1])
        self.assertIsNone(result['purity']['ME'][1])

    def test_all_normalization_bins_required_not_only_surviving_positive_integral(self):
        group = make_group()
        invalidate(group, 'ME', 2)
        result = analyze(group, normalization=(.12, .2))['methods']
        for method in ('method5', 'method6'):
            curve = result[method]
            self.assertFalse(curve['normalization']['valid'])
            self.assertTrue(all(point['cf'] is None for point in curve['points']))
            self.assertTrue(all(value is None for row in curve['covariance'] for value in row))
            self.assertTrue(any(value is not None for value in curve['yields']['ME']))
        self.assertTrue(result['method3']['normalization']['valid'])
        json.dumps(result, allow_nan=False)

    def test_normalization_covariance_and_exact_self_normalized_bin_for_all_methods(self):
        result = analyze(make_group(), normalization=(.12, .14))['methods']
        for method in ('method3', 'method5', 'method6'):
            curve = result[method]
            self.assertAlmostEqual(curve['points'][6]['cf']['value'], 1.0)
            self.assertAlmostEqual(curve['points'][6]['cf']['variance'], 0.0, places=20)
            for index, point in enumerate(curve['points']):
                self.assertAlmostEqual(curve['covariance'][index][index], point['cf']['variance'])
                self.assertAlmostEqual(curve['covariance'][index][6], 0.0, places=15)
            self.assertGreater(curve['covariance'][0][1], 0.0)

    def test_method3_uses_registered_count_background_errors(self):
        result = analyze(make_group())['methods']['method3']
        # SE central N=80, B=16; Var(N-B)=80+Var(B), not fixed-P count error.
        self.assertAlmostEqual(result['yields']['SE'][0]['value'], 64.0)
        self.assertAlmostEqual(result['yields']['SE'][0]['variance'], 80.0+9.6)
        self.assertAlmostEqual(result['purity']['SE'][0]['variance'], 9.6/80**2+16**2*80/80**4)

    def test_overlap_of_me_purity_windows_retains_shared_native_count_covariance(self):
        group = make_group()
        first = group['methods']['method5']['ME'][0]
        second = group['methods']['method5']['ME'][1]
        first['windowOriginalLastBin'] = 20
        second['windowOriginalFirstBin'] = 1
        second['windowOriginalLastBin'] = 20
        result = analyze(group)['methods']['method6']['purity']
        # Both have B=20 and the same N=200. Their fitted B's are independent,
        # but common count uncertainty gives Cov(P0,P1)=20^2*200/200^4.
        self.assertAlmostEqual(result['covariance']['ME'][0][1], 20.0**2*200.0/200.0**4)

    def test_malformed_edges_and_ranges_are_rejected_without_mutating_group(self):
        group = make_group()
        before = copy.deepcopy(group)
        analyze(group)
        self.assertEqual(group, before)
        bad = copy.deepcopy(group)
        bad['kstarEdges'][1] += .001
        with self.assertRaises(ValueError):
            analyze(bad)
        with self.assertRaises(ValueError):
            analyze(group, normalization=(.121, .16))
        with self.assertRaises(ValueError):
            analyze(group, comparison=(.04, .0))

    def test_nonpositive_normalization_variance_is_not_accepted(self):
        group = make_group(background_variance=0.0)
        for row in group['method5Base']['SE']+group['method5Base']['ME']:
            row['error2'] = [0.0, 0.0]
        result = analyze(group)['methods']['method6']
        self.assertFalse(result['normalization']['valid'])
        self.assertIn('variance', result['normalization']['reason'])


if __name__ == '__main__':
    unittest.main()
