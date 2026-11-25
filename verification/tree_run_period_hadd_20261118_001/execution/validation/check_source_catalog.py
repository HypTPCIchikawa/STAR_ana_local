#!/usr/bin/env python3
"""Check the frozen producer snapshot against bounded ROOT metadata/header QA."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-qa', type=int, required=True)
    args = parser.parse_args()
    require(not os.path.lexists(str(args.output)), 'Refusing an existing report')
    snapshot = json.loads(args.snapshot.read_text())
    catalog = json.loads(args.catalog.read_text())
    sources = snapshot['paths']
    records = catalog['files']
    require(catalog.get('passed') is True, 'ROOT catalog did not finish successfully')
    require(catalog.get('fullEventReplayPerformed') is False, 'Unexpected validation scope')
    require(len(records) == len(sources) == snapshot['frozenSourceCount'] == catalog['filesChecked'],
            'Source/catalog count mismatch')
    require([x['path'] for x in records] == [x['path'] for x in sources],
            'Catalog input order differs from frozen source order')
    seen_paths, seen_segments = set(), set()
    total_entries = total_input = samples = 0
    reference = None
    for original, record in zip(sources, records):
        path = Path(original['path'])
        require(str(path) not in seen_paths, 'Duplicate source path')
        seen_paths.add(str(path))
        info = path.lstat()
        require(stat.S_ISREG(info.st_mode), 'Source is not a regular nonsymlink file: ' + str(path))
        actual = dict(bytes=info.st_size, inode=info.st_ino, device=info.st_dev, mtimeNs=info.st_mtime_ns)
        require(all(actual[k] == original[k] for k in actual), 'Frozen source fingerprint changed: ' + str(path))
        require(record['entries'] == record['reconstructedEvents'] == original['acceptedEventsFromCompletionLog'],
                'Accepted event count differs from completion log: ' + str(path))
        require(record['inputEvents'] == original['inputEventsFromCompletionLog'],
                'Input event count differs from completion log: ' + str(path))
        require(len(record['qa']) == args.expected_qa, 'QA count mismatch')
        schema = (record['eventSchemaFnv1a64'], record['metadataSchemaFnv1a64'],
                  tuple((q['name'], q['schemaFnv1a64']) for q in record['qa']), record['compression'])
        if reference is None:
            reference = schema
        require(schema == reference, 'Source event/metadata/QA schema or compression changed')
        declared = 0
        for segment in record['segments']:
            key = (segment['high'], segment['low'])
            require(key != ('0000000000000000', '0000000000000000') and key not in seen_segments,
                    'Invalid or repeated producer segment')
            seen_segments.add(key)
            require(segment['producerId'] and segment['entries'] >= 0, 'Invalid producer metadata')
            declared += segment['entries']
        require(declared == record['entries'], 'Metadata segment extent mismatch')
        total_entries += record['entries']
        total_input += record['inputEvents']
        samples += record['readerSampleEvents']
    require(total_entries == snapshot['acceptedEventsFromCompletionLogs'] == catalog['sourceEvents'],
            'Total accepted entries mismatch')
    require(total_input == snapshot['inputEventsFromCompletionLogs'] == catalog['sourceInputEvents'],
            'Total input entries mismatch')
    result = dict(passed=True, sourceFiles=len(records), producerSegments=len(seen_segments),
                  inputEvents=total_input, acceptedEvents=total_entries,
                  officialReaderSampleEvents=samples, expectedQaPerSource=args.expected_qa,
                  frozenSourceFingerprintsUnchanged=True, fullEventReplayPerformed=False,
                  snapshotSha256=digest(args.snapshot), catalogSha256=digest(args.catalog),
                  scope='Metadata/all tree schemas and counts/QA headers; bounded endpoint and Reader samples. No full event replay.')
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
