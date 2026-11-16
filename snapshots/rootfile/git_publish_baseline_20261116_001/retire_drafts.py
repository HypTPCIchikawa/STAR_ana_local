#!/usr/bin/env python3
"""Retire an unimplemented opt-in study without altering the production baseline."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIVATE=ROOT/'star_local'
DEST=PRIVATE/'drafts/ep_kstar_20261116_001'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert not DEST.exists()
    state=json.loads((HERE/'snapshot_before_draft_archive.json').read_text())
    entries=state['draft_files'];assert len(entries)==9
    for entry in entries:
        p=ROOT/entry['path'];assert p.is_file() and not p.is_symlink()
        assert sha(p)==entry['sha256'] and p.stat().st_size==entry['bytes']
    for name,where in [('main',ROOT),('private',PRIVATE)]:
        diff=subprocess.check_output(['git','-C',str(where),'diff','HEAD','--binary'])
        with (HERE/('before_'+name+'.patch')).open('xb') as stream:stream.write(diff)
    DEST.mkdir(parents=True)
    moved=[]
    for entry in entries:
        source=ROOT/entry['path'];destination=DEST/entry['path']
        destination.parent.mkdir(parents=True,exist_ok=True)
        assert not destination.exists();source.rename(destination)
        assert not source.exists() and sha(destination)==entry['sha256']
        moved.append(dict(entry,retired_path=str(destination.relative_to(PRIVATE))))
    old=ROOT/'rootfile/femto_lambda_ep_kstar_validation_20261115_001'
    if old.exists() and not any(old.iterdir()):old.rmdir()
    manifest=dict(retiredAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        reason='User requested the November 13 production baseline before changing the analysis approach.',
        state='Unimplemented, uncompiled and unexecuted draft; inactive archive only.',
        files=moved,production_source_changed=False,production_configuration_changed=False)
    (DEST/'retirement_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (DEST/'README.md').write_text('# EP区分別k*保存：未実装draftの退避\n\n'
        '2026-11-16、ユーザー指示により11/13本番投入時の解析機能を維持し、この下書きをactiveなconfig/tests/rootfileから退避しました。\n\n'
        '9ファイルの元パス・SHAはretirement_manifest.jsonに記録しています。C++本体は未実装で、テストも未コンパイル・未実行です。'
        '復元してそのまま本番投入できる状態ではありません。既存のEP8混合・pairY10・He4 YLab分割・merging修正は維持しています。\n\n'
        '中断した実装を再開する場合は、別途明示指示に従って実装・検証してください。\n')
    (HERE/'draft_retirement.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('DRAFT_RETIREMENT_OK files=9 production_source_unchanged=YES')


if __name__=='__main__':main()
