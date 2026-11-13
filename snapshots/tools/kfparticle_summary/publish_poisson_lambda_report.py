#!/usr/bin/env python3
"""Publish a verified standalone Poisson report and preserve the previous PDF."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil

from build_method6_report import pdf_pages


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--backup-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    metadata=json.loads((out/'poisson_report_metadata.json').read_text())
    validation=json.loads((out/'validation.json').read_text())
    assert metadata['phase']=='staged' and validation['phase']=='validated'
    source=Path(metadata['report']['path']);target=Path(metadata['originalReport']['path'])
    assert sha(source)==metadata['report']['sha256']==validation['reportSha256']
    assert sha(target)==metadata['originalReport']['sha256']
    assert pdf_pages(source)==metadata['report']['pages']
    backup=args.backup_dir.resolve()
    if backup.exists():raise FileExistsError(backup)
    backup.mkdir(parents=True);preserved=backup/target.name;shutil.copy2(str(target),str(preserved))
    assert sha(preserved)==metadata['originalReport']['sha256']
    temporary=target.with_name(target.name+'.poisson.tmp')
    if temporary.exists():raise FileExistsError(temporary)
    shutil.copy2(str(source),str(temporary));assert sha(temporary)==sha(source)
    os.replace(str(temporary),str(target));assert sha(target)==sha(source)
    record=dict(phase='published',publishedAtJst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        target=str(target),sha256=sha(target),pages=metadata['report']['pages'],bytes=target.stat().st_size,
        backup=str(preserved),backupSha256=sha(preserved),stagedSource=str(source),validationSha256=sha(out/'validation.json'))
    for path in (out/'publication.json',backup/'preservation.json'):
        path.write_text(json.dumps(record,indent=2)+'\n')
    target.with_name('summary_poisson_update.md').write_text('# Poisson likelihood update\n\n'
        'The current summary is the standalone Poisson-L purity/CF report. Negative fitted backgrounds '
        'trigger a nonnegative constant-background refit, as specified in the study settings. '
        'The preceding chi-square report is preserved, with all prior figures unchanged.\n\n'
        '- [Study, settings, source and validation records]('+str(out/'README.md')+')\n'
        '- [Preserved preceding report]('+str(preserved)+')\n'
        '- [Publication record]('+str(out/'publication.json')+')\n\n'
        'The original summary TeX and older metadata describe previous reference versions; '
        'the current report source is recorded in poisson_report_metadata.json.\n')
    print(json.dumps(record))


if __name__=='__main__':main()
