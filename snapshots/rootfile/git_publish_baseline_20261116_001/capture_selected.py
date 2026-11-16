#!/usr/bin/env python3
"""Refresh explicit snapshots while preserving the audited historical inventory."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import shutil
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIVATE=ROOT/'star_local'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    spec=importlib.util.spec_from_file_location('workspace_sync',PRIVATE/'scripts/sync_workspace.py')
    sync=importlib.util.module_from_spec(spec);spec.loader.exec_module(sync);cfg=sync.policy()
    metadata=PRIVATE/'metadata/snapshots.json';old=json.loads(metadata.read_text());records=dict(old['files'])
    selected=json.loads((HERE/'private_capture_paths.json').read_text())
    pending=[]
    for relative in selected:
        sync.safe_relative(relative);source=ROOT/relative
        ok,reason=sync.permitted(source,relative,cfg,False)
        if not ok:raise ValueError('Snapshot rejected: '+relative+' '+reason)
        digest=sha(source);dest=PRIVATE/'snapshots'/relative
        assert not dest.is_symlink()
        if dest.exists() and relative in records:
            assert sha(dest) in {records[relative]['sha256'],digest},'Snapshot independently edited: '+relative
        record=dict(sha256=digest,bytes=source.stat().st_size,
                    kind=records.get(relative,{}).get('kind','personal_document'))
        pending.append((source,dest,relative,record))
    copied=0
    for source,dest,relative,record in pending:
        if not dest.exists() or sha(dest)!=record['sha256']:
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest);copied+=1
        assert sha(dest)==record['sha256'] and sha(source)==record['sha256']
        records[relative]=record
    old['files']=dict(sorted(records.items()))
    old['parentBaseCommit']=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    metadata.write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
    result=dict(capturedAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        method='Explicit incremental refresh; preserve all previously audited snapshot records',
        selected=len(selected),copied=copied,total_records=len(records),
        historical_snapshot_deletions=0,all_selected_source_and_destination_sha256_match=True)
    (HERE/'private_capture_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PRIVATE_CAPTURE_OK selected=%d copied=%d total=%d'%(len(selected),copied,len(records)))


if __name__=='__main__':main()
