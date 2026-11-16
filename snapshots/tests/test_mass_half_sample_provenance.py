"""Mass-half provenance checks without ROOT execution or production writes."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools/kfparticle_summary'))
import compare_mass_halves as compare
import report_mass_halves as report


def current_sample(root):
    return dict(year=2021, included=26955, expected=26994, jobid='A' * 32,
                beam='3.85', stem='summary_4He_KFParticle_3p85', inputRoot=str(root),
                provenance={
                    'outputManifest': dict(path='/tmp/output_manifest.json', sha256='a' * 64),
                    'runmeta': dict(path='/tmp/runmeta.json', sha256='b' * 64)})


class MassHalfSampleTests(unittest.TestCase):
    def test_saved_sample_is_verified_without_replacement(self):
        root = '/tmp/current.root'
        sample = current_sample(root)
        context = dict(energy='3p85', sample=sample)
        before = copy.deepcopy(context)
        with mock.patch.object(report, 'resolve_sample', return_value=copy.deepcopy(sample)) as resolve:
            self.assertIs(report.validated_sample(context, root), sample)
        resolve.assert_called_once_with('3p85', root,
            output_manifest='/tmp/output_manifest.json', runmeta='/tmp/runmeta.json')
        self.assertEqual(context, before)

    def test_changed_evidence_is_rejected(self):
        sample = current_sample('/tmp/current.root')
        fresh = copy.deepcopy(sample)
        fresh['provenance']['runmeta']['sha256'] = 'c' * 64
        with mock.patch.object(report, 'resolve_sample', return_value=fresh):
            with self.assertRaisesRegex(ValueError, 'differs from the saved'):
                report.validated_sample(dict(energy='3p85', sample=sample), sample['inputRoot'])

    def test_legacy_context_keeps_historical_coverage(self):
        root = '/tmp/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root'
        with mock.patch.object(report, 'resolve_sample') as resolve:
            sample = report.validated_sample(dict(energy='3p85'), root)
        resolve.assert_not_called()
        self.assertEqual(sample['included'], 26956)
        self.assertEqual(sample['jobid'], '0BEFC64E62F19BBA534A129CABC9B106')

    def test_legacy_context_rejects_unrecorded_current_input(self):
        with self.assertRaisesRegex(ValueError, 'Coverage is unrecorded'):
            report.validated_sample(dict(energy='3p85'), '/tmp/current.root')

    def test_report_uses_current_coverage_and_relative_root_link(self):
        root = Path('/tmp/results/current.root')
        sample = current_sample(root)
        context = dict(energy='3p85', sample=sample,
                       settings=dict(plot_x_max=.3),
                       base=dict(rebin_group=4, purity_rebin=10, sliding_half_width=2))
        estimate = dict(value=1., error=.01)
        method = dict(integrated={r: dict(cf=estimate, vsFull=dict(ratio=estimate))
                                  for r in report.REGIONS},
                      valid_plot_bins={r: 1 for r in report.REGIONS}, total_plot_bins=1,
                      invalid={r: '' for r in report.REGIONS})
        result = dict(inputRoot=str(root), comparison_range=[0., .2], normalization=[.3, .4],
                      purity={r: dict(value=.97, error=.001, width_scale=1./3)
                              for r in report.REGIONS},
                      groups=[dict(id='all_0_80', label='0--80%',
                                   methods={m: copy.deepcopy(method) for m in report.METHODS})])
        nominal = dict(low=1.10, high=1.13, mean=1.115, sideLow=1.08, sideHigh=1.15)
        extract = dict(massPartition=dict(nominal=nominal,
                          binned={r: dict(lower=1.10, upper=1.13)
                                  for r in ('full', 'low', 'high', 'leftSB', 'rightSB')}),
                       groups=[dict(kstarEdges=[0., .02],
                                    methods=dict(method3=dict(SE=[dict(full=dict(background=None))])))])
        with mock.patch.object(report, 'resolve_sample', return_value=copy.deepcopy(sample)):
            readme = report.build_readme(Path('/tmp/results/pdf/current/mass_halves'),
                                         context, result, extract, [])
            tex = report.build_tex(context, result, extract)
        self.assertIn('[current.root](../../../current.root)', readme)
        for text in (readme, tex):
            self.assertIn('26,955', text)
            self.assertNotIn('26,956', text)
            self.assertIn(sample['jobid'], text)
        self.assertIn('outputManifest', readme)


class BeforeRootExecution(Exception):
    pass


class MassHalfRunnerTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='mass-half-provenance-')
        self.addCleanup(self.directory.cleanup)
        self.project = Path(self.directory.name)
        (self.project / 'config').mkdir()
        self.root = self.project / 'current.root'
        self.root.write_bytes(b'root-test-placeholder')
        self.mainconf = self.project / 'main.yaml'
        self.mainconf.write_text(yaml.safe_dump(dict(maker='maker.yaml')))
        maker = dict(lambdaSignalMean=1.11596, lambdaSignalSigma=.00296,
                     lambdaSignalNSigma=3., lambdaSidebandOuterFactor=4.)
        (self.project / 'config/maker.yaml').write_text(yaml.safe_dump(maker))
        base = dict(rebin_group=4, purity_rebin=10, sliding_half_width=2,
                    normalization=[.3, .4], mass_fit=dict(range=[1.08, 1.15],
                    mean_limits=[1.11, 1.12], sigma_limits=[.001, .005],
                    sigma_initial=.003, minimum_entries=20, minimum_covariance_status=2))
        (self.project / 'base.yaml').write_text(yaml.safe_dump(base))
        self.settings = self.project / 'settings.yaml'
        self.settings.write_text(yaml.safe_dump(dict(summary_settings='base.yaml',
            groups=[dict(id='all_0_80', min=0, max=8)])))
        self.out = self.project / 'output'
        self.sample = current_sample(self.root)
        self.arguments = ['--mainconf', str(self.mainconf), '--input-root', str(self.root),
                          '--energy', '3p85', '--settings', str(self.settings),
                          '--output-dir', str(self.out), '--output-manifest', '/tmp/output_manifest.json',
                          '--runmeta', '/tmp/runmeta.json']
        self.project_patch = mock.patch.object(compare, 'PROJECT', self.project)
        self.project_patch.start()
        self.addCleanup(self.project_patch.stop)

    def run_before_root(self, sample, resume=False):
        with mock.patch.object(compare, 'resolve_sample', return_value=copy.deepcopy(sample)) as resolve:
            with mock.patch.object(compare, 'command', side_effect=BeforeRootExecution):
                with self.assertRaises(BeforeRootExecution):
                    compare.main(self.arguments + (['--resume'] if resume else []))
        return resolve

    def test_explicit_evidence_is_recorded_and_unchanged_resume_is_allowed(self):
        resolve = self.run_before_root(self.sample)
        resolve.assert_called_once_with('3p85', self.root,
            output_manifest=Path('/tmp/output_manifest.json'), runmeta=Path('/tmp/runmeta.json'))
        saved = json.loads((self.out / 'context.json').read_text())
        self.assertEqual(saved['sample'], self.sample)
        self.run_before_root(self.sample, resume=True)

    def test_resume_rejects_modified_manifest_or_runmeta(self):
        self.run_before_root(self.sample)
        for name in ('outputManifest', 'runmeta'):
            with self.subTest(evidence=name):
                changed = copy.deepcopy(self.sample)
                changed['provenance'][name]['sha256'] = 'c' * 64
                with mock.patch.object(compare, 'resolve_sample', return_value=changed):
                    with mock.patch.object(compare, 'command') as command:
                        with self.assertRaisesRegex(ValueError, 'Existing output/config differs'):
                            compare.main(self.arguments + ['--resume'])
                    command.assert_not_called()


if __name__ == '__main__':
    unittest.main()
