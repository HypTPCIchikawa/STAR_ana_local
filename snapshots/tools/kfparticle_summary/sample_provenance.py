"""Resolve historical or explicitly documented KFParticle input samples.

No ROOT import is required.  Current inputs need both a merge-output manifest
and the submission runmeta, including the hash-verified runtime bundle.  Bundle
members are read in place; this module never extracts an archive.
"""
import copy
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import tarfile


SAMPLES = {
    '13p5': {'year': 2020, 'included': 542, 'expected': 545,
             'jobid': '1213058CE62BCF113B2E3C75ACBD7053',
             'beam': '13.5', 'stem': 'summary_4He_KFParticle'},
    '3p85': {'year': 2021, 'included': 26956, 'expected': 26994,
             'jobid': '0BEFC64E62F19BBA534A129CABC9B106',
             'beam': '3.85', 'stem': 'summary_4He_KFParticle_3p85'},
}


def sha256(path):
    """Hash a file without loading a runtime bundle into memory."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _file(path, description, base=None):
    if not isinstance(path, (str, Path)) or not str(path):
        raise ValueError('%s must be a file path.' % description)
    path = Path(path)
    if not path.is_absolute() and base is not None:
        path = base / path
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError('%s is not a regular file: %s' % (description, path))
    return path


def _json(path, description):
    with path.open(encoding='utf-8') as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError('%s must contain a JSON object.' % description)
    return value


def _integer(value, description, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError('%s must be an integer >= %d.' % (description, minimum))
    return value


def _counter(value, description):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or int(value) != value):
        raise ValueError('%s must be a finite, nonnegative integer count.' % description)
    return int(value)


def _record(path):
    return {'path': str(path), 'sha256': sha256(path)}


def _identity(document, jobid, ana_name, description):
    if document.get('jobid') != jobid or document.get('anaName') != ana_name:
        raise ValueError('%s jobid/anaName does not match the output manifest.' % description)


def _bundle_analysis(bundle_path, analysis_info, energy, year):
    """Validate the production recorded inside the submitted runtime bundle."""
    if not isinstance(analysis_info, str):
        raise ValueError('runmeta submit.analysisInfo must be a relative YAML path.')
    relative = PurePosixPath(analysis_info)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Unsafe submit.analysisInfo member path.')
    member_name = str(relative if relative.parts[0] == 'config'
                      else PurePosixPath('config') / relative)
    try:
        import yaml
    except ImportError:
        raise RuntimeError('PyYAML is required to validate submitted analysis metadata.')
    with tarfile.open(str(bundle_path), mode='r:*') as archive:
        matches = [member for member in archive.getmembers() if member.name == member_name]
        if len(matches) != 1 or not matches[0].isfile():
            raise ValueError('Runtime bundle must contain one regular %s.' % member_name)
        member = matches[0]
        if member.size > 2 * 1024 * 1024:
            raise ValueError('Submitted analysis YAML exceeds the metadata size limit.')
        with archive.extractfile(member) as handle:
            payload = handle.read()
        analysis = yaml.safe_load(payload)
    if not isinstance(analysis, dict) or not isinstance(analysis.get('starTag'), dict):
        raise ValueError('Submitted analysis YAML is missing starTag metadata.')
    tags = analysis['starTag']
    expected = {'year': 'year%d' % year, 'productionTag': 'P24iy',
                'triggerSets': 'production_%sGeV_fixedTarget_%d' % (energy, year)}
    for key, value in expected.items():
        if tags.get(key) != value:
            raise ValueError('Submitted analysis starTag.%s must be %s.' % (key, value))
    return {'member': member_name, 'sha256': hashlib.sha256(payload).hexdigest(),
            'starTag': copy.deepcopy(tags)}


def resolve_sample(energy, input_root, output_manifest=None, runmeta=None):
    """Return sample metadata after verifying its declared input provenance.

    With no manifest, retain the historical production filename check and sample
    counts.  A manifest always requires an explicit runmeta; none is discovered
    from a similarly named file.  Relative paths inside JSON are resolved against
    that JSON's directory.  ROOT contents are checked separately by inspection.
    """
    if energy not in SAMPLES:
        raise ValueError('Unknown sample energy: %s' % energy)
    input_path = _file(input_root, 'Input ROOT')
    sample = copy.deepcopy(SAMPLES[energy])
    sample['inputRoot'] = str(input_path)
    sample['provenance'] = {}
    if output_manifest is None:
        if runmeta is not None:
            raise ValueError('runmeta requires an explicit output manifest.')
        expected_name = ('auau%s_anaFemtoLambda_4He_KFParticle_%dP24iy_partial_merge.root'
                         % (energy, sample['year']))
        if input_path.name != expected_name:
            raise ValueError('Historical sample requires input named %s.' % expected_name)
        return sample
    if runmeta is None:
        raise ValueError('An explicit runmeta is required with an output manifest.')
    manifest_path = _file(output_manifest, 'Output manifest')
    runmeta_path = _file(runmeta, 'Submission runmeta')
    manifest = _json(manifest_path, 'Output manifest')
    metadata = _json(runmeta_path, 'Submission runmeta')
    normal_root = _file(manifest.get('normalRoot'), 'Manifest normalRoot', manifest_path.parent)
    if normal_root != input_path:
        raise ValueError('Manifest normalRoot does not match the input ROOT real path.')
    jobid = manifest.get('jobid')
    if not isinstance(jobid, str) or not re.match(r'^[0-9A-Fa-f]{32}$', jobid):
        raise ValueError('Manifest jobid must be a 32-character hexadecimal SUMS ID.')
    ana_name = 'auau%s_anaFemtoLambda_4He_KFParticle' % energy
    if manifest.get('anaName') != ana_name:
        raise ValueError('Manifest anaName does not match the selected energy/analysis.')
    _identity(metadata, jobid, ana_name, 'Submission runmeta')
    included = _integer(manifest.get('inputCount'), 'Manifest inputCount', 1)
    expected = _integer(manifest.get('expectedSubjobs'), 'Manifest expectedSubjobs', 1)
    excluded = manifest.get('excludedIndices')
    if not isinstance(excluded, list):
        raise ValueError('Manifest excludedIndices must be a list.')
    for index in excluded:
        _integer(index, 'Excluded subjob index')
        if index >= expected:
            raise ValueError('Excluded subjob index is outside the zero-based subjob range.')
    if len(set(excluded)) != len(excluded):
        raise ValueError('Manifest excludedIndices contains duplicate subjobs.')
    if included + len(excluded) != expected:
        raise ValueError('Manifest inputCount + excluded count does not equal expectedSubjobs.')
    parts = manifest.get('parts')
    if not isinstance(parts, list) or not parts:
        raise ValueError('Manifest parts must be a nonempty list.')
    part_indices = set()
    part_count = 0
    for part in parts:
        if not isinstance(part, dict):
            raise ValueError('Each manifest part must be an object.')
        index = _integer(part.get('index'), 'Manifest part index')
        if index in part_indices:
            raise ValueError('Manifest parts contain duplicate indices.')
        part_indices.add(index)
        part_count += _integer(part.get('subjobCount'), 'Manifest part subjobCount', 1)
    if part_count != included:
        raise ValueError('Manifest part subjobCount sum does not match inputCount.')
    counters = manifest.get('eventCounters')
    if not isinstance(counters, dict):
        raise ValueError('Manifest eventCounters must be an object.')
    for key in ('inputEventsFromRefMultEntries', 'reconstructedEventsFromRefMultCorrEntries'):
        _counter(counters.get(key), 'Manifest eventCounters.' + key)
    runtime = metadata.get('runtime')
    bundle = runtime.get('runtimeBundle') if isinstance(runtime, dict) else None
    if not isinstance(bundle, dict):
        raise ValueError('runmeta runtime.runtimeBundle metadata is required.')
    bundle_path = _file(bundle.get('path'), 'Runtime bundle', runmeta_path.parent)
    bundle_record = _record(bundle_path)
    declared_hash = bundle.get('sha256')
    if (not isinstance(declared_hash, str) or
            not re.match(r'^[0-9A-Fa-f]{64}$', declared_hash) or
            bundle_record['sha256'] != declared_hash.lower()):
        raise ValueError('Runtime bundle SHA256 does not match the submission runmeta.')
    bundle_record['sizeBytes'] = bundle_path.stat().st_size
    if ('sizeBytes' in bundle and
            _integer(bundle['sizeBytes'], 'Runtime bundle sizeBytes') != bundle_record['sizeBytes']):
        raise ValueError('Runtime bundle sizeBytes does not match the submission runmeta.')
    submit = metadata.get('submit')
    if not isinstance(submit, dict):
        raise ValueError('runmeta submit metadata is required.')
    analysis = _bundle_analysis(bundle_path, submit.get('analysisInfo'), energy, sample['year'])
    sample.update(included=included, expected=expected, jobid=jobid, anaName=ana_name,
                  excludedIndices=sorted(excluded), eventCounters=copy.deepcopy(counters),
                  submittedAnalysis=analysis,
                  provenance={'outputManifest': _record(manifest_path),
                              'runmeta': _record(runmeta_path), 'runtimeBundle': bundle_record})
    return sample


def validate_sample_counters(sample, info):
    """Require inspected histogram-entry event counts to match the manifest.

    Historical samples have no manifest counters and retain their old behavior.
    TParameter counters may also be retained by the inspector as diagnostics,
    but are deliberately not used for this comparison.
    """
    expected = sample.get('eventCounters')
    if expected is None:
        return
    actual = info.get('eventCounters')
    if not isinstance(actual, dict):
        raise ValueError('Inspection eventCounters are required for the declared sample.')
    for inspected, declared in (
            ('inputEvents', 'inputEventsFromRefMultEntries'),
            ('reconstructedEvents', 'reconstructedEventsFromRefMultCorrEntries')):
        measured = _counter(actual.get(inspected), 'Inspection eventCounters.' + inspected)
        reference = _counter(expected.get(declared), 'Manifest eventCounters.' + declared)
        if measured != reference:
            raise ValueError('Inspection %s (%s) does not match the output manifest (%s).'
                             % (inspected, measured, reference))
