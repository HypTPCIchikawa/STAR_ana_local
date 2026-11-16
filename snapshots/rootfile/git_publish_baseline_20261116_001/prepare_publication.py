#!/usr/bin/env python3
"""Prepare explicit public paths and private postprocessing snapshot additions."""
from pathlib import Path
import importlib.util
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIVATE=ROOT/'star_local'


def names(*args):
    return [x.decode() for x in subprocess.check_output(['git','-C',str(ROOT),*args]).split(b'\0') if x]


def main():
    changed=names('diff','--name-only','-z','HEAD')
    public=[x for x in changed if x in {'AGENTS.md','Makefile','docs/REFERENCE.md',
        'docs/ai/AGENT_RULES.md','job/run/README.md','src/HistManager.cpp',
        'script/archive_merged_subjobs.py','script/watch_job_merge_notify.py','script/watch_job_notify.py',
        'tests/test_archive_merged_subjobs.py','tests/test_watch_job_notify.py'}
        or x.startswith(('StMaker/','include/','analysis/'))]
    audit=json.loads((HERE/'baseline_comparison.json').read_text())
    tracked=set(names('ls-files','-z'))
    for group in ('StMaker','include','analysis'):
        for item in audit['groups'][group]['files']:
            path=item['path']
            if path not in tracked and Path(path).suffix.lower() in {'.c','.cpp','.cxx','.h','.hpp','.md'}:
                public.append(path)
    public+=['docs/ai/skills/approval-recovery.md','.cursor/skills/approval-recovery/SKILL.md',
             'tests/test_merge_job_selection.py']
    public=sorted(set(public))
    for path in public:
        assert (ROOT/path).is_file()
        assert not path.startswith(('config/','common/macro/','tools/','rootfile/','share/','job/run/auau'))
        assert Path(path).suffix.lower() not in {'.root','.so','.o','.d','.log','.pdf','.zip','.gz'}
    (HERE/'public_paths.json').write_text(json.dumps(public,indent=2)+'\n')
    (HERE/'public_paths.txt').write_text('\n'.join(public)+'\n')
    spec=importlib.util.spec_from_file_location('workspace_sync',PRIVATE/'scripts/sync_workspace.py')
    sync=importlib.util.module_from_spec(spec);spec.loader.exec_module(sync)
    cfg=sync.policy();candidates=set(changed)-{'star_local'}
    candidates.update(names('ls-files','--others','--exclude-standard','-z',
        'StMaker','include','analysis','common/macro','tools','script','tests','docs/ai'))
    candidates.add('tests/test_merge_job_selection.py')
    study_dirs=[
        'rootfile/femto_lambda_ep_rapidity_cf_20261115_001',
        'rootfile/femto_lambda_split_canvas_20261115_001',
        'rootfile/femto_lambda_cent_pair_canvas_20261115_001',
        'rootfile/femto_lambda_report_update_20261115_001',
        'rootfile/git_publish_baseline_20261116_001']
    for directory in study_dirs:
        for path in (ROOT/directory).iterdir():
            if path.is_file():candidates.add(str(path.relative_to(ROOT)))
    hadd=ROOT/'rootfile/femto_lambda_pair_rapidity_hadd_20261114_001'
    candidates.update(str(p.relative_to(ROOT)) for p in hadd.rglob('merge_notify_audit.C'))
    figures=ROOT/'share/figure/20261115_lambda_he4_ep_rapidity_cf'
    for path in figures.iterdir():
        if path.is_file():candidates.add(str(path.relative_to(ROOT)))
    for folder in ('plots','plots_cropped','centrality_update_plots','archive/20261115_initial22p'):
        for path in (figures/folder).iterdir():
            if path.is_file():candidates.add(str(path.relative_to(ROOT)))
    for folder in ('share/figure/20261115_lambda_split_canvases',
                   'share/figure/20261115_lambda_centrality_pair_canvases'):
        for path in (ROOT/folder).iterdir():
            if path.is_file():candidates.add(str(path.relative_to(ROOT)))
    for path in (ROOT/'share/figure/20261115_lambda_centrality_pair_canvases/by_centrality').glob('*.pdf'):
        candidates.add(str(path.relative_to(ROOT)))
    # Retain this one complete, small scientific report without broadening the general limit.
    detail='share/figure/20261115_lambda_centrality_pair_canvases/centrality_pair_rapidity10_SE_ME_CF_root.pdf'
    assert (ROOT/detail).read_bytes().startswith(b'%PDF-')
    cfg['extraLargeReports'][detail]=4*1024*1024
    accepted=[];rejected=[]
    for relative in sorted(candidates):
        path=ROOT/relative
        if not path.is_file():continue
        ok,reason=sync.permitted(path,relative,cfg,False)
        if ok:accepted.append(relative)
        else:rejected.append(dict(path=relative,reason=reason,bytes=path.stat().st_size))
    cfg['extraPaths']=sorted(set(cfg['extraPaths'])|set(accepted))
    (PRIVATE/'policy.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    (HERE/'private_capture_paths.json').write_text(json.dumps(accepted,indent=2)+'\n')
    (HERE/'publication_plan.json').write_text(json.dumps(dict(public_paths=public,
        private_refresh_paths=accepted,excluded_generated_or_large=rejected,
        specific_report_limit={detail:cfg['extraLargeReports'][detail]},
        active_ep_drafts_removed=True,physics_baseline_preserved=True,
        public_exclusions=['config/ individual settings','common/macro/ postprocessing','tools/',
                           'rootfile/','share/figure/','job runtime packages','all logs/ROOT/build products']),
        ensure_ascii=False,indent=2)+'\n')
    print('PUBLICATION_PLAN_READY public=%d private_refresh=%d excluded=%d'%(len(public),len(accepted),len(rejected)))


if __name__=='__main__':main()
