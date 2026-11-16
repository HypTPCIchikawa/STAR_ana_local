#!/usr/bin/env python3
"""Validate the frozen QA cohort and prepare Condor files without submitting jobs."""
import argparse
import csv
import getpass
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

FIELDS = ['subjob','cache_abs','original_list_abs','size','mtime_ns','inode','device','uuid','producer_log_abs','input_events','accepted_events','selected_lambdas','input_count','list_sha256','producer_root_metadata_verified','requires_worker_root_metadata_verification','urls_json']
SAFE_PATH = re.compile(r'^[A-Za-z0-9_./:+@%=-]+$')
SUCCESS = {'complete', 'completed', 'passed', 'success'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def checked_path(value, must_exist=True):
    p=Path(value)
    if not p.is_absolute() or not SAFE_PATH.fullmatch(str(p)):
        raise ValueError('Unsafe or nonabsolute path: '+str(value))
    p=p.resolve(strict=must_exist)
    if not SAFE_PATH.fullmatch(str(p)):
        raise ValueError('Unsafe canonical path: '+str(value))
    return p


def fingerprint_ok(row):
    st=Path(row['cache_abs']).stat()
    return (st.st_size,st.st_mtime_ns,st.st_ino,st.st_dev)==tuple(int(row[k]) for k in ['size','mtime_ns','inode','device'])


def active_jobs(campaign):
    constraint='Owner == '+json.dumps(getpass.getuser())+' && MergingCampaign == '+json.dumps(campaign)
    cmd=['condor_q','-json','-constraint',constraint,'-attributes','ClusterId,ProcId,JobStatus,MergingSubjob,MergingStage']
    r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if r.returncode:
        raise RuntimeError('Cannot verify campaign queue: '+r.stderr.strip())
    data=json.loads(r.stdout) if r.stdout.strip() else []
    result={}
    for job in data:
        if 'MergingSubjob' not in job:raise ValueError('Campaign queue job missing MergingSubjob')
        ident=int(job['MergingSubjob']);result.setdefault(ident,[]).append(job)
    return result


def completion(row, target, policy_hash, runtime_hash=None):
    """Accept a completed output only when the receipt binds it to this frozen input/policy."""
    if not target.exists():return None
    receipt=target/'receipt.json';summary=target/'summary.json';root=target/'qa.root'
    if not all(p.is_file() for p in [receipt,summary,root]):
        return {'state':'existing_incomplete_or_failed','directory':str(target)}
    r=json.loads(receipt.read_text());s=json.loads(summary.read_text())
    if r.get('status') not in SUCCESS or s.get('passed') is not True or s.get('closureChecks')!=244827 or root.stat().st_size<=0:
        return {'state':'existing_incomplete_or_failed','directory':str(target)}
    if r.get('outputRootSha256')!=sha(root):
        return {'state':'existing_output_checksum_mismatch','directory':str(target)}
    expected_fp={'size':int(row['size']),'mtimeNs':int(row['mtime_ns']),'ino':int(row['inode']),'device':int(row['device'])}
    if r.get('sourceCacheFingerprint')!=expected_fp or s.get('counters',{}).get('cache_events')!=int(row['accepted_events']):
        return {'state':'receipt_candidate_fingerprint_or_count_mismatch','directory':str(target)}
    # The worker receipt's names are explicit rather than inferred from an output filename.
    if r.get('sourceCache') != row['cache_abs'] or r.get('sourceList') != row['original_list_abs']:
        return {'state':'receipt_input_mismatch','directory':str(target)}
    if r.get('policySha256') != policy_hash or r.get('inputListSha256') != row['list_sha256']:
        return {'state':'receipt_policy_or_list_mismatch','directory':str(target)}
    if runtime_hash is not None and r.get('runtimeManifestSha256')!=runtime_hash:
        return {'state':'receipt_runtime_mismatch','directory':str(target)}
    return {'state':'verified_complete','directory':str(target),'receiptSha256':sha(receipt),'summarySha256':sha(summary),'rootSize':root.stat().st_size}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--worker',type=Path,required=True)
    p.add_argument('--policy',type=Path,required=True)
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--stage',choices=['pilot','main'],required=True)
    p.add_argument('--inventory-index',type=Path)
    p.add_argument('--inventory-validation',type=Path)
    p.add_argument('--runtime-dir',type=Path)
    p.add_argument('--runtime-validation',type=Path)
    p.add_argument('--results-dir',type=Path)
    p.add_argument('--submission-label')
    p.add_argument('--require-successful-pilots',help='Comma-separated cohort IDs required to pass before actual submission')
    p.add_argument('--ids',help='Explicit comma-separated subjob IDs; required for pilot')
    p.add_argument('--max-materialize',type=int,default=32)
    p.add_argument('--max-idle',type=int,default=16)
    p.add_argument('--wall-seconds',type=int,default=7200)
    a=p.parse_args()
    if a.max_materialize<1 or a.max_idle<1 or a.max_idle>a.max_materialize or a.wall_seconds<1:
        raise ValueError('Invalid operational limits')
    if a.stage=='pilot' and not a.ids:raise ValueError('Pilot requires explicit --ids')
    campaign=checked_path(str(a.campaign));inventory=checked_path(str(a.inventory));worker=checked_path(str(a.worker));policy=checked_path(str(a.policy))
    if inventory.parent != campaign:raise ValueError('Inventory must be inside the selected campaign')
    runtime=checked_path(str(a.runtime_dir)) if a.runtime_dir else worker.parent
    if runtime.parent!=campaign or not re.fullmatch(r'runtime(?:_v[0-9]+)?',runtime.name):raise ValueError('Runtime must be a versioned frozen campaign directory')
    if worker.parent!=runtime or worker.name!='worker.sh':raise ValueError('Worker must be selected runtime/worker.sh')
    if policy.parent!=runtime:raise ValueError('Policy must be inside the selected frozen runtime')
    results=checked_path(str(a.results_dir),False) if a.results_dir else campaign/'results'
    if results.parent!=campaign:raise ValueError('Results directory must be an immediate campaign child')
    label=a.submission_label or a.stage
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',label):raise ValueError('Unsafe submission label')
    runtime_manifest=runtime/'manifest.json'
    runtime_doc=json.loads(runtime_manifest.read_text())
    for rel,item in runtime_doc['files'].items():
        if sha(runtime/rel)!=item['sha256']:raise ValueError('Frozen runtime file changed: '+rel)
    runtime_hash=sha(runtime_manifest)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',campaign.name):raise ValueError('Unsafe campaign tag')
    meta_path=checked_path(str(a.inventory_index)) if a.inventory_index else campaign/'input_inventory.json'
    if meta_path.parent!=campaign:raise ValueError('Inventory index must be in selected campaign')
    meta=json.loads(meta_path.read_text())
    if meta.get('status') not in ['inventory_complete_primary_preflight_cohort','dispatch_inventory_complete_preflight_cohort']:raise ValueError('Inventory cohort preflight has not passed')
    if Path(meta['candidates_tsv']).resolve()!=inventory or sha(inventory)!=meta['candidates_tsv_sha256']:
        raise ValueError('Inventory manifest path/hash mismatch')
    if meta.get('dataset')!='production_3p85GeV_fixedTarget_2021 / P24iy':raise ValueError('Wrong dataset')
    independent_path=checked_path(str(a.inventory_validation)) if a.inventory_validation else campaign/'input_inventory_validation.json'
    if independent_path.parent!=campaign:raise ValueError('Inventory validation must be in selected campaign')
    independent=json.loads(independent_path.read_text())
    if independent.get('passed') is not True:raise ValueError('Independent inventory validation has not passed')
    for path in [meta_path,inventory]:
        if sha(path)!=independent['hashes'][path.name]:raise ValueError('Independent inventory hash mismatch')
    explicit_ids=set(int(x) for x in a.ids.split(',')) if a.ids else set()
    with inventory.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=FIELDS:raise ValueError('Inventory TSV header mismatch')
        rows=list(reader)
    if len(rows)!=meta['summary']['preflight_eligible_jobs']:raise ValueError('Inventory row count mismatch')
    if not rows:raise ValueError('Empty eligible cohort')
    byid={};urls=set()
    for row in rows:
        ident=int(row['subjob'])
        if str(ident)!=row['subjob'] or ident in byid:raise ValueError('Invalid or duplicate subjob ID')
        byid[ident]=row
        if int(row['input_events'])<=0 or int(row['accepted_events'])<=0 or int(row['accepted_events'])>int(row['input_events']):raise ValueError('Invalid measured event counts')
        if row['requires_worker_root_metadata_verification']!='True':raise ValueError('Worker metadata verification must be required')
        if row['producer_root_metadata_verified'] not in ['True','False']:raise ValueError('Invalid ROOT-verification state')
        for key in ['cache_abs','original_list_abs','producer_log_abs']:
            if not Path(row[key]).is_absolute() or not SAFE_PATH.fullmatch(row[key]):raise ValueError('Unsafe inventory path')
        # The frozen full-cohort inventory was already independently checked.
        # Reopen only explicit pilot/selected IDs here; every worker repeats ROOT/layout/input checks.
        if ident in explicit_ids:
            checked_path(row['cache_abs']);listing=checked_path(row['original_list_abs']);checked_path(row['producer_log_abs'])
            if not fingerprint_ok(row):raise ValueError('Candidate file fingerprint changed: '+str(ident))
            if sha(listing)!=row['list_sha256']:raise ValueError('Input-list fingerprint changed: '+str(ident))
        raw_urls=json.loads(row['urls_json'])
        if len(raw_urls)!=int(row['input_count']) or len(set(raw_urls))!=len(raw_urls):raise ValueError('Bad URL inventory')
        if ident in explicit_ids:
            actual=[line.split()[0] for line in listing.read_text().splitlines() if line.strip()]
            if actual!=raw_urls:raise ValueError('Original input order changed')
        for url in raw_urls:
            if not (url.startswith('root://xrdstar.rcf.bnl.gov:1095//') and '/production_3p85GeV_fixedTarget_2021/' in url and '/P24iy/2021/' in url and url.endswith('.picoDst.root')):
                raise ValueError('Wrong PicoDst dataset in cohort')
            if url in urls:raise ValueError('PicoDst overlap in cohort')
            urls.add(url)
    requested=None
    if a.ids:
        requested=[int(x) for x in a.ids.split(',')]
        if len(set(requested))!=len(requested) or any(x not in byid for x in requested):raise ValueError('Requested IDs not unique eligible members')
    selected=[byid[x] for x in sorted(requested)] if requested is not None else sorted(rows,key=lambda r:int(r['subjob']))
    active=active_jobs(campaign.name);policy_hash=sha(policy)
    include=[];excluded=[]
    for row in selected:
        ident=int(row['subjob']);dest=results/str(ident)
        existing=completion(row,dest,policy_hash,runtime_hash)
        if existing:
            excluded.append(dict(subjob=ident,**existing));continue
        if ident in active:
            excluded.append({'subjob':ident,'state':'already_queued','jobs':active[ident]});continue
        include.append(row)
    if not include:raise ValueError('No new jobs remain; refusing empty or duplicate submission plan')
    submit=campaign/(label+'.submit');items=campaign/(label+'.items.tsv');manifest=campaign/(label+'_submit_manifest.json');ads=campaign/(label+'_dry_run.ads');drylog=campaign/(label+'_dry_run.log')
    for path in [submit,items,manifest,ads,drylog]:
        if path.exists():raise FileExistsError('Preserve existing preparation: '+str(path))
    (campaign/'logs').mkdir(exist_ok=True);results.mkdir(exist_ok=True)
    with items.open('w') as f:
        for row in include:
            f.write('\t'.join([row['subjob'],row['cache_abs'],row['original_list_abs'],str(results/row['subjob']),str(policy)])+'\n')
    content='''Universe = vanilla
Executable = /bin/bash
Arguments = "'''+str(worker)+''' $(cache_root) $(input_list) $(out_dir) $(policy_file)"
Initialdir = '''+str(campaign)+'''
GetEnv = False
Notification = never
should_transfer_files = NO
request_cpus = 1
request_memory = 2GB
request_disk = 8GB
Requirements = (TARGET.HasSingularity =?= true)
kill_sig = SIGINT
PeriodicRemove = (JobStatus == 2) && ((CurrentTime - JobCurrentStartDate) > '''+str(a.wall_seconds)+''')
Output = '''+str(campaign)+'''/logs/$(subjob).$(Cluster).$(Process).out
Error = '''+str(campaign)+'''/logs/$(subjob).$(Cluster).$(Process).err
Log = '''+str(campaign)+'''/logs/'''+label+'''.condor.log
+MergingCampaign = "'''+campaign.name+'''"
+MergingStage = "'''+a.stage+'''"
+MergingSubjob = $(subjob)
+MergingSubmissionLabel = "'''+label+'''"
'''
    if a.stage=='main':content+='max_materialize = '+str(a.max_materialize)+'\nmax_idle = '+str(a.max_idle)+'\n'
    content+='queue subjob,cache_root,input_list,out_dir,policy_file from '+str(items)+'\n'
    submit.write_text(content)
    runtime_sources=[path for path in sorted(runtime.rglob('*')) if path.is_file() and path.suffix in ['.sh','.py','.json','.yaml','.yml','.cxx','.C','.h','.so']]
    runtime_sources=[path for path in runtime_sources if '/validation/' not in str(path)]
    doc={'status':'prepared_not_submitted','createdAtUtc':datetime.now(timezone.utc).isoformat(),'stage':a.stage,'submissionLabel':label,'runtimeDirectory':str(runtime),'runtimeManifestSha256':runtime_hash,'resultsDirectory':str(results),'campaign':campaign.name,'dataset':meta['dataset'],'sourceJobid':meta['source_jobid'],'eligibleCohortJobs':len(rows),'plannedJobs':len(include),'plannedSubjobs':[int(row['subjob']) for row in include],'excluded':excluded,'measuredInputEvents':sum(int(row['input_events']) for row in include),'measuredAcceptedEvents':sum(int(row['accepted_events']) for row in include),'counterProvenance':meta.get('scope'),'sourceRootMetadataMayNeedWorkerVerification':True,'preparationFilesystemRechecks':'Explicit --ids only; full inventory structure/hash/URL uniqueness independently validated; every worker repeats ROOT and original-input checks','independentInventoryValidationSha256':sha(independent_path),'requiredWorkerChecks':meta['requirements_before_dispatch'],'inputs':{'inventory':str(inventory),'inventorySha256':sha(inventory),'inventoryManifest':str(meta_path),'inventoryManifestSha256':sha(meta_path),'worker':str(worker),'workerSha256':sha(worker),'policy':str(policy),'policySha256':policy_hash},'submit':str(submit),'submitSha256':sha(submit),'items':str(items),'itemsSha256':sha(items),'sourceSha256':sha(Path(__file__)),'runtimeFileSha256':{str(path.relative_to(campaign)):sha(path) for path in runtime_sources},'sampleJobs':[{'subjob':int(row['subjob']),'cache':row['cache_abs'],'inputList':row['original_list_abs'],'outputDirectory':str(results/row['subjob'])} for row in include[:3]],'limits':{'cpus':1,'memoryMB':2048,'diskKB':8388608,'wallSeconds':a.wall_seconds,'maxMaterialize':a.max_materialize if a.stage=='main' else None,'maxIdle':a.max_idle if a.stage=='main' else None},'queueSnapshot':{str(k):v for k,v in active.items()},'executionBlockedUntilFrozenRuntimePass':True}
    if a.runtime_validation:
        rv=checked_path(str(a.runtime_validation));validation=json.loads(rv.read_text())
        if validation.get('status')!='passed':raise ValueError('Frozen runtime validation has not passed')
        if validation.get('manifestSha256')!=runtime_hash or Path(validation.get('runtime','')).resolve()!=runtime:
            raise ValueError('Runtime PASS does not match the selected frozen runtime')
        doc.update(executionBlockedUntilFrozenRuntimePass=False,runtimeValidation=str(rv),runtimeValidationSha256=sha(rv))
    cmd=['condor_submit','-dry-run',str(ads),'-single-cluster','-maxjobs',str(len(include)),str(submit)]
    r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    drylog.write_text(r.stdout);doc['dryRunCommand']=cmd;doc['dryRunReturnCode']=r.returncode
    if r.returncode:
        doc['status']='dry_run_failed_not_submitted';write_json(manifest,doc);raise RuntimeError('Condor dry-run failed; see '+str(drylog))
    text=ads.read_text();ad_jobs=len(re.findall(r'^ProcId\s*=',text,re.MULTILINE))
    for key,value in [('RequestCpus',1),('RequestMemory',2048),('RequestDisk',8388608)]:
        found=re.findall(r'^'+key+r'\s*=\s*(\d+)\s*$',text,re.MULTILINE)
        if not found or any(int(x)!=value for x in found):raise ValueError('Unexpected dry-run '+key)
    if str(worker) not in text or 'MergingCampaign' not in text or 'HasSingularity' not in text:raise ValueError('Dry-run arguments or custom attributes missing')
    doc.update(status='dry_run_passed_not_submitted',dryRunAdsSha256=sha(ads),dryRunMaterializedAdCount=ad_jobs)
    required=[int(x) for x in a.require_successful_pilots.split(',')] if a.require_successful_pilots else []
    if len(set(required))!=len(required) or any(x not in byid for x in required):raise ValueError('Required pilots are not distinct cohort members')
    states={str(x):completion(byid[x],results/str(x),policy_hash,runtime_hash) for x in required}
    passed=all(x is not None and x.get('state')=='verified_complete' for x in states.values())
    doc.update(requiresSuccessfulPilotSubjobs=required,executionBlockedUntilPilotPass=bool(required) and not passed,pilotGateSnapshot=states)
    write_json(manifest,doc)
    print('SUBMIT_PREPARATION_PASS',a.stage,'planned',len(include),'excluded',len(excluded),'NO_ACTUAL_SUBMISSION')
    print(submit)


if __name__=='__main__':main()
