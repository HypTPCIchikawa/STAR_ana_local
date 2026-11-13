#!/usr/bin/env python3
"""Submit the frozen main factory once, only after both real-data pilots validate."""
import argparse,csv,fcntl,hashlib,json,os,re,subprocess,time,traceback
from datetime import datetime,timezone
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def load(p):return json.loads(p.read_text())
def atomic(p,d):
 t=p.with_name(p.name+'.partial');t.write_text(json.dumps(d,indent=2)+'\n');os.replace(str(t),str(p))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--campaign',type=Path,required=True);ap.add_argument('--check-only',action='store_true');ap.add_argument('--poll-seconds',type=int,required=True);a=ap.parse_args()
 c=a.campaign.resolve(strict=True)
 with (c/'advance.lock').open('a+') as lock:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:raise RuntimeError('Factory advancer already active')
  pilot=load(c/'pilot_v2_submission.json');pilot_ids={str(x) for x in pilot['subjobs']}
  with (c/'dispatch_candidates.tsv').open() as f:rows={r['subjob']:r for r in csv.DictReader(f,delimiter='\t')}
  while True:
   done=[];pending=[];errors=[]
   for sid in sorted(pilot_ids,key=int):
    out=c/'results_v2'/sid;marker=out/'worker_result.json'
    if not marker.exists():pending.append(sid);continue
    try:
     r=load(marker)
     if r.get('status')=='failed':raise RuntimeError('Pilot worker failed: '+str(r.get('error',r)))
     if r.get('status')!='completed':pending.append(sid);continue
     s=load(out/'summary.json');row=rows[sid];fp=r['sourceCacheFingerprint']
     if not s.get('passed') or s.get('closureChecks')!=244827:raise ValueError('Pilot full closure not passed')
     if r['runtimeManifestSha256']!=pilot['runtimeManifestSha256'] or r['policySha256']!=pilot['policySha256']:raise ValueError('Pilot used a different runtime/policy')
     if r['sourceCache']!=row['cache_abs'] or r['sourcePicoList']!=row['original_list_abs'] or r['inputListSha256']!=row['list_sha256']:raise ValueError('Pilot input identity mismatch')
     if tuple(fp[k] for k in ['size','mtimeNs','ino'])!=tuple(int(row[k]) for k in ['size','mtime_ns','inode']):raise ValueError('Pilot candidate fingerprint mismatch')
     if s['counters']['cache_events']!=int(row['accepted_events']):raise ValueError('Pilot accepted-frame coverage mismatch')
     if sha(out/'qa.root')!=r['outputRootSha256']:raise ValueError('Pilot QA checksum mismatch')
     done.append(sid)
    except Exception as e:errors.append(dict(subjob=sid,error=str(e)))
   state={'updatedAtUtc':now(),'phase':'waiting_for_pilots','pilotSubjobs':sorted(pilot_ids,key=int),'validatedPilots':done,'pendingPilots':pending,'errors':errors}
   if errors:
    state['phase']='stopped_pilot_failure';atomic(c/'advance_status.json',state);print(json.dumps(state),flush=True);return 2
   if set(done)!=pilot_ids:
    atomic(c/'advance_status.json',state);print(now(),'waiting pilots',pending,flush=True)
    if a.check_only:return 0
    time.sleep(a.poll_seconds);continue
   mpath=c/'main_v2_submit_manifest.json'
   if not mpath.exists():
    state['phase']='waiting_for_main_preparation';atomic(c/'advance_status.json',state)
    if a.check_only:return 0
    time.sleep(a.poll_seconds);continue
   m=load(mpath);sub=c/'main_v2.submit';items=c/'main_v2.items.tsv'
   if m.get('status')!='dry_run_passed_not_submitted':raise ValueError('Main factory dry-run not passed')
   if sha(sub)!=m['submitSha256'] or sha(items)!=m['itemsSha256']:raise ValueError('Prepared factory changed')
   if sha(c/'dispatch_candidates.tsv')!=m['inputs']['inventorySha256']:raise ValueError('Dispatch cohort changed')
   if set(str(x) for x in m['plannedSubjobs'])!=set(rows)-pilot_ids:raise ValueError('Factory does not cover exactly the remaining cohort')
   if m['limits']['maxMaterialize']!=32 or m['limits']['maxIdle']!=16:raise ValueError('Factory concurrency differs from reviewed limits')
   if sha(c/'runtime_v2/manifest.json')!=pilot['runtimeManifestSha256'] or sha(c/'runtime_v2/policy.json')!=pilot['policySha256']:raise ValueError('Frozen runtime changed after pilot')
   frozen=load(c/'runtime_v2/manifest.json')
   for rel,d in frozen['files'].items():
    if sha(c/'runtime_v2'/rel)!=d['sha256']:raise ValueError('Frozen runtime member changed: '+rel)
   state.update(phase='pilots_validated_factory_ready',plannedFactoryJobs=m['plannedJobs'],plannedFactorySubjobs=m['plannedSubjobs'],factorySubmit=str(sub))
   atomic(c/'advance_status.json',state)
   gate={'passed':True,'checkedAtUtc':now(),'pilotSubjobs':sorted(pilot_ids,key=int),'runtimeManifestSha256':pilot['runtimeManifestSha256'],'policySha256':pilot['policySha256'],'checks':'Both real-data pilot completion receipts, exact old QA closure, source/list/fingerprints/coverage and QA SHA verified'}
   atomic(c/'main_v2_pilot_gate_validation.json',gate)
   if a.check_only:print('FACTORY_READY_NO_SUBMISSION',m['plannedJobs'],flush=True);return 0
   receipt=c/'main_v2_submission.json';intent=c/'main_v2_submission.intent'
   if receipt.exists():
    state['phase']='factory_already_submitted';atomic(c/'advance_status.json',state);return 0
   with intent.open('x') as f:f.write(now()+'\n')
   if m.get('executionBlockedUntilPilotPass'):
    preserved=c/'main_v2_submit_manifest.pre_pilot.json'
    if not preserved.exists():atomic(preserved,m)
    m['executionBlockedUntilPilotPass']=False;m['pilotGateValidation']=str(c/'main_v2_pilot_gate_validation.json');m['pilotGateVerifiedAtUtc']=gate['checkedAtUtc'];atomic(mpath,m)
   cmd=['condor_submit','-terse',str(sub)];r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
   (c/'main_v2_submission.stdout').write_text(r.stdout);(c/'main_v2_submission.stderr').write_text(r.stderr)
   if r.returncode:raise RuntimeError('Factory submission failed; reconcile preserved intent, do not resubmit blindly')
   matched=re.match(r'\s*(\d+)\.(\d+)',r.stdout)
   if not matched:raise RuntimeError('Factory submission returned without parseable ID; reconcile queue')
   d={'status':'submitted','submittedAtUtc':now(),'cluster':int(matched.group(1)),'jobCount':m['plannedJobs'],'subjobs':m['plannedSubjobs'],'stdout':r.stdout,'stderr':r.stderr,'submitSha256':m['submitSha256'],'itemsSha256':m['itemsSha256'],'runtimeManifestSha256':pilot['runtimeManifestSha256'],'policySha256':pilot['policySha256'],'validatedPilotSubjobs':sorted(pilot_ids,key=int),'concurrency':{'maxMaterialize':32,'maxIdle':16},'command':cmd}
   atomic(receipt,d);state.update(phase='factory_submitted',factoryCluster=d['cluster']);atomic(c/'advance_status.json',state)
   print('FACTORY_SUBMITTED',d['cluster'],m['plannedJobs'],flush=True);return 0

if __name__=='__main__':
 try:raise SystemExit(main())
 except Exception as e:
  traceback.print_exc()
  try:
   import sys
   c=Path(sys.argv[sys.argv.index('--campaign')+1]).resolve();atomic(c/'advance_failure.json',{'phase':'stopped_validation_or_submission_error','timeUtc':now(),'error':str(e)})
  except Exception:pass
  raise SystemExit(3)
