#!/usr/bin/env python3
"""Validate and aggregate finalized local-preview QA sources, separately from farm QA."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic(path, value):
    temp = path.with_name(path.name + '.partial')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(str(temp), str(path))


def exports(path):
    value = load(path)
    hist = {item['name']: item for item in value['histograms'].values()}
    if len(hist) != value['histogramCount']:
        raise ValueError('Histogram identity/count mismatch')
    for name, item in hist.items():
        if any(not math.isfinite(x) or x < 0 for x in item['content'] + item['error2']):
            raise ValueError('Invalid raw bin in ' + name)
        if abs(sum(item['content']) - item['entries']) > 1.e-6:
            raise ValueError('Raw bins/Entries mismatch in ' + name)
    return hist


def exact_sum(parts, merged):
    names = set().union(*(set(part) for part in parts))
    if set(merged) != names:
        raise ValueError('Merged histogram membership differs from source union')
    bins_checked = 0
    error2_rounded_values = 0
    max_error2_difference = 0.0
    exact_sumw2_histograms = 0
    fallback_error_histograms = 0
    for name in sorted(names):
        present = [part[name] for part in parts if name in part]
        result = merged[name]
        for item in present:
            if any(item[key] != result[key] for key in ['class', 'dimension', 'xEdges', 'yEdges', 'shape']):
                raise ValueError('Histogram schema differs in ' + name)
        if result['entries'] != sum(item['entries'] for item in present):
            raise ValueError('Merged/source Entries differ in ' + name)
        exact_error2 = result['sumw2N'] > 0 and all(item['sumw2N'] > 0 for item in present)
        if exact_error2:
            exact_sumw2_histograms += 1
        else:
            fallback_error_histograms += 1
        for field in ['content', 'error2']:
            target = result[field]
            if any(len(item[field]) != len(target) for item in present):
                raise ValueError('Histogram array size differs in ' + name)
            expected = [sum(values) for values in zip(*(item[field] for item in present))]
            if field == 'content' and target != expected:
                raise ValueError('Merged raw content differs from source sum in ' + name)
            if field == 'error2' and exact_error2 and target != expected:
                raise ValueError('Merged stored Sumw2 differs from source sum in ' + name)
            if field == 'error2' and not exact_error2:
                # Without Sumw2 the exporter squares GetBinError(), which can round.
                for actual, total in zip(target, expected):
                    if actual == total:
                        continue
                    scale = max(actual, total)
                    exponent = math.frexp(scale)[1]
                    tolerance = math.ldexp(1.0, exponent - 53) * 8 if scale else 0.0
                    difference = abs(actual - total)
                    if difference > tolerance:
                        raise ValueError('Merged error2 exceeds eight double ULP in ' + name)
                    error2_rounded_values += 1
                    max_error2_difference = max(max_error2_difference, difference)
            bins_checked += len(target)
    return dict(passed=True, histogramCount=len(names), rawArrayValuesChecked=bins_checked,
                includesUnderflowOverflow=True, includesErrorSquared=True, contentAndEntriesExact=True,
                error2UlpTolerance=8, error2RoundedValues=error2_rounded_values,
                maximumError2Difference=max_error2_difference,
                exactStoredSumw2Histograms=exact_sumw2_histograms,
                fallbackErrorHistograms=fallback_error_histograms)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'qa_input_manifest.json')
    parser.add_argument('--once', action='store_true', required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve()
    local = source.parents[1]
    campaign = local.parent
    manifest_path = args.input_manifest.resolve()
    if manifest_path.parent != local:
        raise ValueError('Input manifest must belong to the local extension')
    if not manifest_path.is_file():
        print('LOCAL_AGGREGATE_WAIT finalized QA input manifest is not available')
        return 0
    manifest_raw = manifest_path.read_bytes()
    manifest = json.loads(manifest_raw)
    if manifest.get('status') != 'final':
        print('LOCAL_AGGREGATE_WAIT QA sources have not all reached a final outcome')
        return 0
    if manifest.get('schemaVersion') != 1:
        raise ValueError('Unexpected local QA manifest schema')
    rows = {str(row['id']): row for row in manifest['sources']}
    if set(rows) != {'0', '1', '2', '3'} or len(rows) != len(manifest['sources']):
        raise ValueError('Expected exactly local source identities 0..3')
    if rows['0']['kind'] != 'baseline' or any(rows[s]['kind'] != 'additional' for s in ['1', '2', '3']):
        raise ValueError('Baseline/additional source identity differs')
    runtime = Path(manifest['runtime']).resolve(strict=True)
    policy = Path(manifest['policy']).resolve(strict=True)
    validation_path = Path(manifest['runtimeValidation']).resolve(strict=True)
    validation = load(validation_path)
    if runtime != local / 'runtime' or policy.parent != runtime or validation.get('status') != 'passed':
        raise ValueError('Local frozen runtime has not passed validation')
    runtime_hash, policy_hash = sha(runtime / 'manifest.json'), sha(policy)
    if Path(validation['runtime']).resolve() != runtime or validation['manifestSha256'] != runtime_hash:
        raise ValueError('Frozen local runtime validation identity mismatch')
    frozen = load(runtime / 'manifest.json')
    for relative, entry in frozen['files'].items():
        if sha(runtime / relative) != entry['sha256']:
            raise ValueError('Frozen local runtime file changed: ' + relative)
    provenance_path = local / 'input_provenance.json'
    provenance = load(provenance_path)
    baseline = load(campaign / 'validation_v2/frozen_smoke_001/receipt.json')
    if baseline.get('status') != 'completed' or not baseline.get('qaSummary', {}).get('passed'):
        raise ValueError('Known baseline receipt has not passed')
    proven = {str(i): row for i, row in enumerate(provenance['inputs'], 1)}
    settings_path = campaign / 'monitor_settings.json'
    cfg = load(settings_path)
    inventory_path = Path(cfg['inventory']).resolve(strict=True)
    if sha(inventory_path) != cfg['inventory_sha256']:
        raise ValueError('Farm dispatch inventory hash changed')
    with inventory_path.open() as handle:
        farm_lfns = {url.rsplit('/', 1)[-1] for row in csv.DictReader(handle, delimiter='\t') for url in json.loads(row['urls_json'])}
    local_lfns = [Path(original['path']).name for original in provenance['inputs']]
    local_lfns += [Path(item['source']).name for item in baseline['inputs']]
    if len(local_lfns) != len(set(local_lfns)) or set(local_lfns) & farm_lfns:
        raise ValueError('Local inputs duplicate each other or overlap the farm dispatch inventory')
    helper_path = campaign / 'source/watch_campaign.py'
    spec = importlib.util.spec_from_file_location('local_aggregate_tools', str(helper_path))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    aggregate = local / 'aggregate'
    aggregate.mkdir(exist_ok=True)
    with (aggregate / 'aggregate.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verified, excluded = {}, {}
        for sid, row in sorted(rows.items(), key=lambda item: int(item[0])):
            if row['qaOutcome'] == 'failed':
                if not row.get('exclusionReason'):
                    raise ValueError('Failed local source has no exclusion reason')
                excluded[sid] = row['exclusionReason']
                continue
            if row['qaOutcome'] != 'completed':
                raise ValueError('Final manifest contains a nonfinal source')
            directory = local / 'qa' / sid
            q, receipt, marker = (load(directory / name) for name in ['summary.json', 'receipt.json', 'worker_result.json'])
            root = directory / 'qa.root'
            if not q.get('passed') or q.get('closureChecks') != 244827 or receipt != marker or receipt.get('status') != 'completed' or receipt.get('qaSummary') != q:
                raise ValueError('QA commit/full legacy closure mismatch for ' + sid)
            cache, listing = Path(row['cacheRoot']).resolve(strict=True), Path(row['picoList']).resolve(strict=True)
            if Path(receipt['sourceCache']).resolve() != cache or Path(receipt['sourcePicoList']).resolve() != listing:
                raise ValueError('QA input identity mismatch for ' + sid)
            if receipt['runtimeManifestSha256'] != runtime_hash or receipt['policySha256'] != policy_hash or sha(listing) != receipt['inputListSha256']:
                raise ValueError('QA frozen runtime/policy/list mismatch for ' + sid)
            if Path(receipt['runtime']).resolve() != runtime or sha(root) != receipt['outputRootSha256'] or sha(cache) != receipt['sourceCacheSha256']:
                raise ValueError('QA ROOT/cache integrity mismatch for ' + sid)
            fp = cache.stat()
            expected_fp = receipt['sourceCacheFingerprint']
            if (fp.st_size, fp.st_mtime_ns, fp.st_ino) != (expected_fp['size'], expected_fp['mtimeNs'], expected_fp['ino']):
                raise ValueError('Candidate fingerprint mismatch for ' + sid)
            if q['counters']['cache_events'] != row['acceptedCacheEvents'] or not 0 <= row['acceptedCacheEvents'] <= row['sourceInputEventsMeasured']:
                raise ValueError('Accepted/source-event coverage mismatch for ' + sid)
            metadata_path = Path(row['rootMetadata']).resolve(strict=True)
            metadata = load(metadata_path)
            if metadata_path != directory / 'root_metadata.json' or sha(metadata_path) != row['rootMetadataSha256']:
                raise ValueError('Strict ROOT metadata identity differs for ' + sid)
            if metadata.get('schemaVersion') != 1 or metadata.get('passed') is not True or metadata.get('runStatus') != 'completed' or metadata.get('oldHistogramClosureChecks') != 244827:
                raise ValueError('Strict ROOT metadata has not passed for ' + sid)
            if metadata['sourceCacheSha256'] != receipt['sourceCacheSha256'] or not metadata['uuid']:
                raise ValueError('Strict ROOT cache identity differs for ' + sid)
            if any(metadata[key] != row['sourceInputEventsMeasured'] for key in ['inputEvents', 'eventLedgerEntries', 'picoEntries']):
                raise ValueError('Strict full source-event coverage differs for ' + sid)
            if any(metadata[key] != row['acceptedCacheEvents'] for key in ['reconstructedEvents', 'cacheEntries']):
                raise ValueError('Strict accepted-cache coverage differs for ' + sid)
            if any(metadata[key] != q['counters']['cache_lambda_candidates'] for key in ['selectedLambdaCandidates', 'selectedLambdaEntries', 'cachedLambdaCandidates']):
                raise ValueError('Strict Lambda candidate coverage differs for ' + sid)
            if sid == '0':
                if cache != Path(baseline['sourceCache']).resolve() or receipt['inputListSha256'] != baseline['inputListSha256'] or receipt['sourceCacheSha256'] != baseline['sourceCacheSha256']:
                    raise ValueError('Known baseline candidate/input-list identity differs')
                if row['sourceInputEventsMeasured'] != 4027 or row['acceptedCacheEvents'] != 3858:
                    raise ValueError('Known baseline coverage differs')
            else:
                original = proven[sid]
                producer_path = Path(row['producerResult']).resolve(strict=True)
                producer = load(producer_path)
                if producer_path != local / 'producer' / sid / 'result.json' or producer['returnCode'] != 0 or not producer['inputUnchanged'] or not producer['producerLibrariesUnchanged']:
                    raise ValueError('Individual producer has not passed for ' + sid)
                if Path(producer['candidateRoot']).resolve() != cache or Path(producer['picoList']).resolve() != listing or row['sourceInputEventsMeasured'] != original['sourceInputEventsMeasured']:
                    raise ValueError('Producer/provenance input identity mismatch for ' + sid)
                if original['overlapsCombinedCampaignLfn'] or Path(producer['input']['path']).resolve() != Path(original['path']).resolve():
                    raise ValueError('Local source overlap/provenance mismatch for ' + sid)
                lines = [line.strip() for line in listing.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
                if len(lines) != 1 or Path(lines[0]).resolve() != Path(original['path']).resolve():
                    raise ValueError('Local producer input-list content differs for ' + sid)
            verified[sid] = dict(row=row, root=root, summary=q, receipt=receipt, rootMetadata=metadata)
        if not verified:
            raise ValueError('No successful local QA sources')
        jobs = sorted(verified, key=int)
        identity = hashlib.sha256((hashlib.sha256(manifest_raw).hexdigest() + ''.join(verified[s]['receipt']['outputRootSha256'] for s in jobs)).encode()).hexdigest()[:12]
        stem = 'qa_' + str(len(jobs)).zfill(5) + '_' + identity
        root = aggregate / (stem + '.root')
        extraction = aggregate / (stem + '.histograms.json')
        previous = aggregate / (stem + '.json')
        if previous.exists():
            result = load(previous)
            if not result.get('passed') or result['rootSha256'] != sha(root) or sha(extraction) != result['histogramExportSha256']:
                raise ValueError('Existing local aggregate integrity mismatch')
            atomic(aggregate / 'latest.json', result)
        else:
            counters = Counter()
            parts = []
            source_validation = {}
            for sid in jobs:
                item = verified[sid]
                counts = item['summary']['counters']
                counters.update(counts)
                exported = aggregate / ('source_' + sid + '_' + item['receipt']['outputRootSha256'][:12] + '.histograms.json')
                source_validation[sid] = helper.inspect(cfg, item['root'], exported, counts)
                if source_validation[sid]['histogramCount'] != item['summary']['histogramCount']:
                    raise ValueError('Source exported histogram count differs from QA summary for ' + sid)
                parts.append(exports(exported))
            temporary = aggregate / (stem + '.pending.root')
            command = helper.utility(cfg, 'hadd', ['-f', str(temporary)] + [str(verified[s]['root']) for s in jobs])
            with (aggregate / (stem + '.hadd.log')).open('w') as log:
                done = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if done.returncode or not temporary.is_file():
                raise RuntimeError('Local QA hadd failed')
            os.replace(str(temporary), str(root))
            check = helper.inspect(cfg, root, extraction, dict(counters))
            check['fullRawBinSum'] = exact_sum(parts, exports(extraction))
            if any(sha(verified[s]['root']) != verified[s]['receipt']['outputRootSha256'] for s in jobs) or manifest_path.read_bytes() != manifest_raw:
                raise ValueError('Committed local QA inputs changed during aggregation')
            coverage = {}
            for kind in ['baseline', 'additional']:
                subset = [verified[s] for s in jobs if verified[s]['row']['kind'] == kind]
                coverage[kind] = dict(sourceCount=len(subset), sourceInputEventsMeasured=sum(x['row']['sourceInputEventsMeasured'] for x in subset), acceptedCacheFrames=sum(x['row']['acceptedCacheEvents'] for x in subset), picoEventsActuallyRead=sum(x['summary']['counters']['pico_events_read'] for x in subset), selectedHe4Pairs=sum(x['summary']['counters'].get('signal_pairs_selected_4He', 0) for x in subset))
            result = dict(passed=True, createdAtUtc=datetime.now(timezone.utc).isoformat(), scope='local preview only; kept separate from farm campaign', dataset=provenance['dataset'], farmDispatchedJobsAdded=0, farmDispatchInputOverlap=0, farmInventorySha256=sha(inventory_path), expectedLocalSources=len(rows), excludedSourceCount=len(excluded), root=str(root), rootSha256=sha(root), sourceJobs=jobs, sourceCount=len(jobs), counters=dict(counters), validation=check, inputQaRoots=[str(verified[s]['root']) for s in jobs], sourceInputEventsRepresented=sum(verified[s]['row']['sourceInputEventsMeasured'] for s in jobs), sourceAcceptedEventsRepresented=sum(verified[s]['row']['acceptedCacheEvents'] for s in jobs), picoEventsActuallyRead=counters['pico_events_read'], coverageByKind=coverage, excludedSources=excluded, histogramExport=str(extraction), histogramExportSha256=sha(extraction), sourceManifest=str(manifest_path), sourceManifestSha256=hashlib.sha256(manifest_raw).hexdigest(), inputProvenanceSha256=sha(provenance_path), localRuntimeManifestSha256=runtime_hash, localPolicySha256=policy_hash, sourceValidation=source_validation, sourceQaHashes={s: verified[s]['receipt']['outputRootSha256'] for s in jobs}, rootMetadata={s: verified[s]['rootMetadata'] for s in jobs}, aggregateSourceSha256=sha(source), helperSourceSha256=sha(helper_path), monitorSettingsSha256=sha(settings_path))
            atomic(previous, result)
            atomic(aggregate / 'latest.json', result)
        plotter = campaign / 'source/plot_campaign.py'
        command = ['python3', str(plotter), '--campaign', str(local), '--output-dir', str(local / 'plots'), '--once']
        with (aggregate / 'plot.log').open('w') as log:
            done = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if done.returncode:
            raise RuntimeError('Local aggregate verified, but drawing failed; see aggregate/plot.log')
        print('LOCAL_AGGREGATE_PASS sources', len(jobs), 'He4', result['counters'].get('signal_pairs_selected_4He', 0), 'excluded', excluded, aggregate / 'latest.json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
