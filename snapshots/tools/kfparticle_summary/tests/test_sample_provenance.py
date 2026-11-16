"""Input-provenance gates: these tests do not run ROOT or a CF calculation."""
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / 'sample_provenance.py'
SPEC = importlib.util.spec_from_file_location('sample_provenance_under_test', str(MODULE))
provenance = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(provenance)


class SampleProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base / 'current_merge.root'
        self.root.write_bytes(b'root-test-placeholder')
        self.manifest_path = self.base / 'outputs.json'
        self.runmeta_path = self.base / 'runmeta.json'
        self.bundle_path = self.base / 'runtime.tar.gz'
        self.member = 'config/analysis/analysis_info_current.yaml'
        self.yaml = ('starTag:\n'
                     '  year: year2021\n'
                     '  productionTag: P24iy\n'
                     '  triggerSets: production_3p85GeV_fixedTarget_2021\n'
                     '  libraryTag: SL24y\n').encode('utf-8')
        self.write_bundle()
        self.manifest = {
            'normalRoot': str(self.root), 'jobid': 'A' * 32,
            'anaName': 'auau3p85_anaFemtoLambda_4He_KFParticle',
            'inputCount': 2, 'expectedSubjobs': 3, 'excludedIndices': [1],
            'parts': [{'index': 0, 'file': 'part000.root', 'subjobCount': 2}],
            'eventCounters': {'inputEventsFromRefMultEntries': 1000.0,
                              'reconstructedEventsFromRefMultCorrEntries': 950.0}}
        self.runmeta = {
            'jobid': 'A' * 32, 'anaName': self.manifest['anaName'],
            'submit': {'analysisInfo': 'analysis/analysis_info_current.yaml'},
            'runtime': {'runtimeBundle': self.bundle_record()}}
        self.write_jsons()

    def write_bundle(self, member=None, kind=None):
        with tarfile.open(str(self.bundle_path), 'w:gz') as archive:
            entry = tarfile.TarInfo(member or self.member)
            if kind == 'symlink':
                entry.type = tarfile.SYMTYPE
                entry.linkname = '/outside/analysis.yaml'
                archive.addfile(entry)
            else:
                entry.size = len(self.yaml)
                archive.addfile(entry, io.BytesIO(self.yaml))

    def bundle_record(self):
        return {'path': str(self.bundle_path), 'sha256': provenance.sha256(self.bundle_path),
                'sizeBytes': self.bundle_path.stat().st_size}

    def write_jsons(self):
        self.manifest_path.write_text(json.dumps(self.manifest), encoding='utf-8')
        self.runmeta_path.write_text(json.dumps(self.runmeta), encoding='utf-8')

    def resolve(self):
        self.write_jsons()
        return provenance.resolve_sample('3p85', self.root, self.manifest_path,
                                         self.runmeta_path)

    def refresh_bundle(self):
        self.runmeta['runtime']['runtimeBundle'] = self.bundle_record()

    def test_valid_sample_records_real_paths_and_hashes_without_modifying_inputs(self):
        inputs = [self.root, self.bundle_path]
        before = [p.read_bytes() for p in inputs]
        sample = self.resolve()
        self.assertEqual(sample['included'], 2)
        self.assertEqual(sample['expected'], 3)
        self.assertEqual(sample['year'], 2021)
        self.assertEqual(sample['beam'], '3.85')
        self.assertEqual(sample['stem'], provenance.SAMPLES['3p85']['stem'])
        self.assertEqual(sample['jobid'], self.manifest['jobid'])
        self.assertEqual(sample['inputRoot'], str(self.root.resolve()))
        for key, path in [('outputManifest', self.manifest_path),
                          ('runmeta', self.runmeta_path), ('runtimeBundle', self.bundle_path)]:
            self.assertEqual(sample['provenance'][key]['path'], str(path.resolve()))
            self.assertEqual(sample['provenance'][key]['sha256'], provenance.sha256(path))
        self.assertEqual(sample['submittedAnalysis']['member'], self.member)
        self.assertEqual(sample['submittedAnalysis']['sha256'], hashlib.sha256(self.yaml).hexdigest())
        self.assertEqual(before, [p.read_bytes() for p in inputs])

    def test_historical_samples_keep_counts_and_require_exact_documented_filename(self):
        for energy, historical in provenance.SAMPLES.items():
            with self.subTest(energy=energy):
                name = ('auau%s_anaFemtoLambda_4He_KFParticle_%dP24iy_partial_merge.root'
                        % (energy, historical['year']))
                path = self.base / name
                path.write_bytes(b'root')
                sample = provenance.resolve_sample(energy, path)
                for key, value in historical.items():
                    self.assertEqual(sample[key], value)
                original_included = historical['included']
                sample['included'] = -1
                self.assertEqual(provenance.SAMPLES[energy]['included'], original_included)
        with self.assertRaises(ValueError):
            provenance.resolve_sample('3p85', self.root)

    def test_unknown_energy_and_implicit_runmeta_are_rejected(self):
        with self.assertRaises(ValueError):
            provenance.resolve_sample('2019', self.root)
        with self.assertRaises(ValueError):
            provenance.resolve_sample('3p85', self.root, self.manifest_path)
        with self.assertRaises(ValueError):
            provenance.resolve_sample('3p85', self.root, runmeta=self.runmeta_path)

    def test_input_real_path_must_match_manifest(self):
        other = self.base / 'different.root'
        other.write_bytes(b'root')
        self.manifest['normalRoot'] = str(other)
        with self.assertRaisesRegex(ValueError, 'normalRoot'):
            self.resolve()
        self.manifest['normalRoot'] = str(self.root)
        link = self.base / 'same.root'
        link.symlink_to(self.root)
        self.write_jsons()
        sample = provenance.resolve_sample('3p85', link, self.manifest_path, self.runmeta_path)
        self.assertEqual(sample['inputRoot'], str(self.root.resolve()))

    def test_manifest_and_submission_identity_must_match(self):
        for document, key, value in [('manifest', 'jobid', 'not-a-sums-id'),
                                     ('manifest', 'anaName', 'auau2019_wrong_analysis'),
                                     ('runmeta', 'jobid', 'B' * 32),
                                     ('runmeta', 'anaName', 'another-analysis')]:
            with self.subTest(document=document, key=key):
                target = getattr(self, document)
                original = target[key]
                target[key] = value
                with self.assertRaises(ValueError):
                    self.resolve()
                target[key] = original

    def test_exclusions_and_count_arithmetic_are_checked(self):
        invalid = [('inputCount', 0), ('inputCount', 1), ('inputCount', True),
                   ('expectedSubjobs', 2.0), ('excludedIndices', [1, 1]),
                   ('excludedIndices', [-1]), ('excludedIndices', [3]),
                   ('excludedIndices', [True]), ('excludedIndices', '1')]
        for key, value in invalid:
            with self.subTest(key=key, value=value):
                original = self.manifest[key]
                self.manifest[key] = value
                with self.assertRaises(ValueError):
                    self.resolve()
                self.manifest[key] = original

    def test_part_counts_and_unique_indices_are_required(self):
        invalid = [[], [{'index': 0, 'subjobCount': 1}],
                   [{'index': 0, 'subjobCount': 1}, {'index': 0, 'subjobCount': 1}],
                   [{'index': 0, 'subjobCount': True}], [None]]
        for parts in invalid:
            with self.subTest(parts=parts):
                self.manifest['parts'] = parts
                with self.assertRaises(ValueError):
                    self.resolve()

    def test_bundle_content_hash_and_size_are_checked(self):
        for key, value in [('sha256', '0' * 64), ('sha256', 'bad-hash'), ('sizeBytes', 1)]:
            with self.subTest(key=key):
                self.refresh_bundle()
                self.runmeta['runtime']['runtimeBundle'][key] = value
                with self.assertRaises(ValueError):
                    self.resolve()

    def test_bundled_analysis_must_be_2021_p24iy_correct_production(self):
        original = self.yaml
        for old, new in [(b'year2021', b'year2019'), (b'P24iy', b'P21ib'),
                         (b'production_3p85GeV_fixedTarget_2021', b'production_wrong_2021')]:
            with self.subTest(replacement=new):
                self.yaml = original.replace(old, new)
                self.write_bundle()
                self.refresh_bundle()
                with self.assertRaisesRegex(ValueError, 'starTag'):
                    self.resolve()

    def test_archive_yaml_is_not_followed_or_extracted(self):
        self.write_bundle(kind='symlink')
        self.refresh_bundle()
        with self.assertRaisesRegex(ValueError, 'regular'):
            self.resolve()
        self.write_bundle()
        self.refresh_bundle()
        self.runmeta['submit']['analysisInfo'] = '../outside.yaml'
        with self.assertRaisesRegex(ValueError, 'Unsafe'):
            self.resolve()
        self.runmeta['submit']['analysisInfo'] = '/outside.yaml'
        with self.assertRaisesRegex(ValueError, 'Unsafe'):
            self.resolve()

    def test_event_counters_match_inspected_hist_entries_and_ignore_parameters(self):
        sample = self.resolve()
        info = {'eventCounters': {'inputEvents': 1000, 'reconstructedEvents': 950.0},
                'eventCounterParameters': {'inputEvents': 1, 'reconstructedEvents': 2}}
        self.assertIsNone(provenance.validate_sample_counters(sample, info))
        for value in [949, None, float('nan'), 950.5, True]:
            with self.subTest(value=value):
                bad = copy.deepcopy(info)
                bad['eventCounters']['reconstructedEvents'] = value
                with self.assertRaises(ValueError):
                    provenance.validate_sample_counters(sample, bad)
        with self.assertRaises(ValueError):
            provenance.validate_sample_counters(sample, {})
        self.assertIsNone(provenance.validate_sample_counters({'year': 2021}, {}))

    def test_nonfinite_or_missing_manifest_counters_are_rejected(self):
        for value in [None, -1, float('inf'), True, '1000']:
            with self.subTest(value=value):
                self.manifest['eventCounters']['inputEventsFromRefMultEntries'] = value
                with self.assertRaises(ValueError):
                    self.resolve()

    def test_relative_json_paths_are_resolved_against_their_metadata(self):
        self.manifest['normalRoot'] = self.root.name
        self.runmeta['runtime']['runtimeBundle']['path'] = self.bundle_path.name
        sample = self.resolve()
        self.assertEqual(sample['inputRoot'], str(self.root.resolve()))
        self.assertEqual(sample['provenance']['runtimeBundle']['path'], str(self.bundle_path.resolve()))


if __name__ == '__main__':
    unittest.main()
