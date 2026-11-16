#!/usr/bin/env python3
"""Watch committed QA results and merge each successful source job exactly once."""
import argparse,csv,fcntl,hashlib,json,math,os,subprocess,time,traceback
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path


def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def atomic(p,obj):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.partial')
 tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(str(tmp),str(p))
def read_if(p):
 try:return load(p)
 except (FileNotFoundError,json.JSONDecodeError):return None

def utility(cfg,command,args):
 return ['bash',cfg['root_batch_wrapper'],cfg['frozen_mainconf'],command]+list(args)

def inspect(cfg,root,jsonfile,counters):
 expression=cfg['existing_histogram_exporter']+'+("'+str(root)+'","'+str(jsonfile)+'")'
 cmd=utility(cfg,'root4star',['-l','-b','-q',expression])
 with jsonfile.with_suffix('.log').open('w') as log:
  r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError('Merged histogram extraction failed')
 d=load(jsonfile);hs={h['name']:h for h in d['histograms'].values()};metrics=[]
 for h in hs.values():
  if not all(math.isfinite(v) and v>=0 for v in h['content']+h['error2']):raise ValueError('Invalid merged bin')
  if abs(sum(h['content'])-h['entries'])>1.e-6:raise ValueError('Merged Entries/bin-sum mismatch '+h['name'])
 for name,template in cfg['expected_counter_templates'].items():
  for species in cfg['species']:
   for daughter in cfg['daughters']:
    key='hDphiDeta_selected_'+species+'_'+daughter+'_'+name
    counter=template.format(species=species,daughter=daughter)
    expected=counters.get(counter,0);h=hs.get(key)
    if h is None:
     if expected:raise ValueError('Missing expected nonempty histogram '+key)
     continue
    if h['entries']!=expected:raise ValueError('Merged counter/Entries mismatch '+key)
    nx,ny=h['shape'][1],h['shape'][0];a=h['content']
    regular=sum(a[y*nx+x] for y in range(1,ny-1) for x in range(1,nx-1));boxes={}
    for w in cfg['origin_halfwidths']:
     def inds(edges):
      ix=[i for i,v in enumerate(edges) if abs(v+w)<1.e-9];jx=[i for i,v in enumerate(edges) if abs(v-w)<1.e-9]
      if len(ix)!=1 or len(jx)!=1:raise ValueError('Origin window not aligned to histogram edges')
      return range(ix[0]+1,jx[0]+1)
     ix,iy=inds(h['xEdges']),inds(h['yEdges']);boxes[str(w)]=sum(a[y*nx+x] for y in iy for x in ix)
    metrics.append(dict(name=key,entries=h['entries'],regular=regular,flow=h['entries']-regular,originBoxes=boxes))
 return dict(passed=True,histogramCount=len(hs),countersChecked=True,histograms=metrics,command=cmd)

def main():
 p=argparse.ArgumentParser();p.add_argument('--settings',type=Path,required=True);p.add_argument('--once',action='store_true');a=p.parse_args()
 cfg=load(a.settings);campaign=Path(cfg['campaign']);project=Path(cfg['project']);os.chdir(str(project))
 with (campaign/'monitor.lock').open('a+') as lock:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:raise RuntimeError('Campaign monitor is already active')
  inventory=Path(cfg['inventory']);assert sha(inventory)==cfg['inventory_sha256']
  with inventory.open() as f:rows={str(r['subjob']):r for r in csv.DictReader(f,delimiter='\t')}
  previous=read_if(campaign/'aggregate/latest.json');merged=set(previous['sourceJobs']) if previous else set()
  verified={};last_merge=time.monotonic();terminal=False
  while not terminal:
   errors=[];failed={};pending=[]
   results_dir=Path(cfg['results_dir'])
   live_dirs={e.name for e in os.scandir(str(results_dir)) if e.is_dir(follow_symlinks=False)} if results_dir.exists() else set()
   for sid,row in rows.items():
    if sid not in live_dirs:
     pending.append(sid);continue
    d=results_dir/sid;marker=read_if(d/'worker_result.json')
    if not marker:pending.append(sid);continue
    if marker.get('status')=='failed':failed[sid]=marker;continue
    if marker.get('status')!='completed':pending.append(sid);continue
    try:
     qa=d/'qa.root';summary=d/'summary.json';receipt=d/'receipt.json'
     st=qa.stat();finger=(st.st_size,st.st_mtime_ns,st.st_ino)
     if sid in verified:
      if verified[sid]['finger']!=finger:raise ValueError('Committed QA ROOT fingerprint changed')
      continue
     q=load(summary);r=load(receipt)
     if not q.get('passed') or r.get('status')!='completed':raise ValueError('Incomplete worker transaction')
     if str(Path(r['sourceCache']).resolve())!=str(Path(row['cache_abs']).resolve()):raise ValueError('Worker cache identity mismatch')
     if r['inputListSha256']!=row['list_sha256']:raise ValueError('Worker input-list identity mismatch')
     fp=r['sourceCacheFingerprint']
     if (fp['size'],fp['mtimeNs'],fp['ino'])!=(int(row['size']),int(row['mtime_ns']),int(row['inode'])):raise ValueError('Source candidate fingerprint differs from frozen inventory')
     if r['runtimeManifestSha256']!=cfg['runtime_manifest_sha256'] or r['policySha256']!=cfg['policy_sha256']:raise ValueError('Worker runtime/policy differs from frozen campaign')
     if sha(qa)!=r['outputRootSha256']:raise ValueError('Committed QA ROOT hash mismatch')
     if q['counters']['cache_events']!=int(row['accepted_events']):raise ValueError('Worker accepted-event count mismatch')
     verified[sid]=dict(summary=q,receipt=r,finger=finger)
    except Exception as e:errors.append(dict(subjob=sid,error=str(e)))
   counters=Counter();stage_bytes=0;stage_seconds=0.;qa_seconds=0.
   for x in verified.values():
    counters.update(x['summary']['counters']);r=x['receipt'];stage_bytes+=r['picoStageBytes'];stage_seconds+=r['stageSeconds'];qa_seconds+=r['qaProcessSeconds']
   queue_ok=False;queue=[];queue_error=None
   constraint='Owner == "'+cfg['owner']+'" && MergingCampaign == "'+cfg['campaign_tag']+'"'
   qr=subprocess.run(['condor_q','-constraint',constraint,'-json'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
   if qr.returncode==0:
    try:queue=json.loads(qr.stdout) if qr.stdout.strip() else [];queue_ok=True
    except Exception as e:queue_error=str(e)
   else:queue_error=qr.stderr.strip()
   counts=Counter(str(x.get('JobStatus')) for x in queue)
   status=dict(updatedAtUtc=now(),phase='monitoring',expectedSourceJobs=len(rows),completedJobs=len(verified),failedJobs=len(failed),pendingJobs=len(rows)-len(verified)-len(failed),queueReadOk=queue_ok,queueByStatus=dict(counts),queueError=queue_error,counters=dict(counters),sourceInputEventsRepresented=sum(int(rows[s]['input_events']) for s in verified),sourceAcceptedEventsRepresented=sum(int(rows[s]['accepted_events']) for s in verified),picoFilesStaged=sum(x['receipt']['picoFileCount'] for x in verified.values()),picoStageBytes=stage_bytes,summedStageSeconds=stage_seconds,summedQaProcessSeconds=qa_seconds,integrityErrors=errors,latestAggregate=str(campaign/'aggregate/latest.json') if previous else None)
   status['dispatchStage']='main_v2' if (campaign/'main_v2_submission.json').exists() else 'pilot_v2'
   if errors:status['phase']='needs_review_integrity'
   new=set(verified)-merged
   if new and not errors and (not previous or len(new)>=cfg['merge_min_new_jobs'] or time.monotonic()-last_merge>=cfg['merge_interval_seconds'] or len(verified)+len(failed)==len(rows)):
    try:
     aggregate=campaign/'aggregate';aggregate.mkdir(exist_ok=True)
     union=set(merged);current=Path(previous['root']) if previous else None
     ordered=sorted(new,key=int)
     for pos in range(0,len(ordered),cfg['merge_batch_size']):
      batch=ordered[pos:pos+cfg['merge_batch_size']];union.update(batch)
      stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');stem='qa_'+str(len(union)).zfill(5)+'_'+stamp
      target=aggregate/(stem+'.root');temp=aggregate/(stem+'.pending.root');inputs=([str(current)] if current else [])+[str(Path(cfg['results_dir'])/s/'qa.root') for s in batch]
      cmd=utility(cfg,'hadd',['-f',str(temp)]+inputs)
      with (aggregate/(stem+'.hadd.log')).open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
      if r.returncode or not temp.is_file():raise RuntimeError('QA-only hadd failed')
      os.replace(str(temp),str(target));current=target
      partc=Counter()
      for sid in union:partc.update(verified[sid]['summary']['counters'])
      extraction=aggregate/(stem+'.histograms.json');check=inspect(cfg,current,extraction,dict(partc))
      manifest=dict(passed=True,createdAtUtc=now(),root=str(current),rootSha256=sha(current),sourceJobs=sorted(union,key=int),sourceCount=len(union),counters=dict(partc),validation=check,inputQaRoots=[str(Path(cfg['results_dir'])/s/'qa.root') for s in sorted(union,key=int)],sourceInputEventsRepresented=sum(int(rows[s]['input_events']) for s in union),sourceAcceptedEventsRepresented=sum(int(rows[s]['accepted_events']) for s in union),histogramExport=str(extraction))
      atomic(aggregate/(stem+'.json'),manifest);atomic(aggregate/'latest.json',manifest);previous=manifest;merged=set(union)
     last_merge=time.monotonic();status['latestAggregate']=str(aggregate/'latest.json')
    except Exception as e:
     status['mergeError']=str(e);status['phase']='needs_review_merge';traceback.print_exc()
   atomic(campaign/'failed_jobs.json',failed)
   if len(verified)+len(failed)==len(rows) and not errors and merged==set(verified) and queue_ok and not queue:
    status['phase']='complete' if not failed else 'complete_with_failures';terminal=True
   atomic(campaign/'status.json',status)
   print(now(),status['phase'],'complete',len(verified),'failed',len(failed),'remaining',status['pendingJobs'],'queue',dict(counts),'He4',counters.get('signal_pairs_selected_4He',0),flush=True)
   if a.once or terminal:break
   time.sleep(cfg['poll_seconds'])

if __name__=='__main__':main()
