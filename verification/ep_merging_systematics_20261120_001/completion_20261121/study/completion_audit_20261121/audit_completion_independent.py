#!/usr/bin/env python3
"""Independently recheck stored successful completions without replay or writes to inputs."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,re,struct,subprocess
ROOT=Path(__file__).resolve().parents[3]
E=Path(__file__).resolve().parent.parent
W=E/'aggregate_watch_001'
load=lambda p:json.loads(Path(p).read_text())
checks=0

def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def statmatch(x):
 s=Path(x['path']).stat()
 for k,v in [('bytes',s.st_size),('device',s.st_dev),('inode',s.st_ino),('mtimeNs',s.st_mtime_ns)]:need(x[k]==v,x['path']+' '+k)
 return s

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def main():
 c=load(W/'completion_verified.json');spec=load(E/'cases.json');matrix=load(E/'production_matrix.json');manifest=load(spec['inputManifest']);sources=manifest['sourceFiles'];wanted=[str(Path(x['path']).resolve()) for x in sources]
 need(c['passed'] and c['eventsPerCase']==1819447887,'completion gate/events')
 need(len(c['cases'])==8 and len(sources)==20 and len(matrix['cases'])==7,'case/input dimensions')
 for x in sources:statmatch(x)
 post={x['path']:x for x in load(W/'cf_comparison/input_provenance.json')['sourceStats']};need(len(post)==160,'postprocess160')
 reports=[];ends=[];logs=0
 for case in c['cases']:
  need(case['passed'] and len(case['workers'])==20,'case20passed')
  need(case['completedEvents']==1819447887,'caseevents')
  need([x['period'] for x in case['workers']]==list(range(20)),'periodorder')
  need([str(Path(x['source']).resolve()) for x in case['workers']]==wanted,'sourceorder')
  total=0;times=[]
  for x,source in zip(case['workers'],sources):
   n=source['acceptedEvents'];need(x['expectedStoredEvents']==n,'expectedcounts')
   q=x['completion'];need(q['requested']==q['completed']==q['available']==n and q['status']==0,'allprocessed')
   out=x['stdout'];err=x['stderr'];need(out['normalNativeExit'] and not out['fatalLines'] and not err['fatalLines'],'normal/fatal')
   need(len(out['maker'])==1 and out['maker'][0]['status']=='completed','Makercompleted')
   need(out['maker'][0]['inputEvents']==out['maker'][0]['reconstructedEvents']==n,'Makercounts')
   need(len(out['chainCounts'])==1,'chainrecord');ch=out['chainCounts'][0]
   need(ch['nStOK']==n and all(ch[k]==0 for k in ['nStWarn','nStEOF','nStErr','nStFatal']),'nochainerrors')
   for z in [out,err]:statmatch(z['stat']);need(sha(z['path'])==z['sha256'],'logSHA');logs+=1
   text=Path(out['path']).read_text(errors='replace');m=re.findall(r'FemtoLambda tree replay job=(\S+) requested=(\d+) completed=(\d+) available=(\d+) status=(\d+)',text)
   need(len(m)==1 and m[0][0]==q['job'] and tuple(map(int,m[0][1:]))==(n,n,n,0),'livereplaymarker')
   rr=x['root'];ss=statmatch(rr);statmatch(post[rr['path']]);total+=ss.st_size
   with Path(rr['path']).open('rb') as f:h=f.read(32)
   need(h[:4]==b'root','rootmagic');ver=struct.unpack('>I',h[4:8])[0];begin=struct.unpack('>I',h[8:12])[0];end=struct.unpack('>Q' if ver>=1000000 else '>I',h[12:20] if ver>=1000000 else h[12:16])[0];need(100<=begin<=end<=ss.st_size,'headerbounds')
   dt=datetime.strptime(out['endRaw'],'%a %b %d %I:%M:%S %p EST %Y').replace(tzinfo=timezone(timedelta(hours=-5)));times.append(dt);ends.append((dt,case['caseId'],x['period']))
  need(total==case['totalOutputBytes'],'casebytes')
  reports.append(dict(caseId=case['caseId'],workers=20,events=case['completedEvents'],rootBytes=total,lastFinishedAmericaNewYork=max(times).isoformat()))
 for item in matrix['cases']:
  coverage=load(Path(item['launchDirectory'])/'input_coverage.json');need(coverage['passed'] and coverage['actualJobs']==20,'coverage20')
  need([str(Path(x['input']).resolve()) for x in coverage['orderedInputs']]==wanted,'nativeorder')
  for row in coverage['orderedInputs']:
   lines=[z.split() for z in Path(row['nativeList']).read_text().splitlines() if z.strip() and not z.lstrip().startswith('#')]
   need(len(lines)==1 and str(Path(lines[0][0]).resolve())==wanted[row['subjob']],'currentnativelist')
 ids=[x['jobid'] for x in load(E/'matrix_launch_state.json')['cases']];host=load(matrix['cases'][0]['config'])['submitHost'];constraint=' || '.join('regexp("'+x+'",Out)' for x in ids);queues=[]
 for fmt in [['-json','-attributes','ClusterId,ProcId,JobStatus,Out'],['-autoformat','ClusterId','ProcId','JobStatus']]:
  command=['condor_q','-name',host,'-constraint',constraint]+fmt;rr=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=45,check=True)
  need(not rr.stderr.strip(),'queue stderr')
  if fmt[0]=='-json':need(not rr.stdout.strip() or json.loads(rr.stdout)==[],'JSON queueempty')
  else:need(not rr.stdout.strip(),'autoformatqueueempty')
  queues.append(dict(command=command,returnCode=rr.returncode,stdout=rr.stdout,stderr=rr.stderr))
 sent=load(W/'sent.json');need(sent['kind']=='success' and sent['mailStatus']=='accepted_by_local_sendmail','successmail')
 last=max(ends)
 result=dict(passed=True,checkedAtUTC=datetime.now(timezone.utc).isoformat(),checks=checks,newWorkers=140,reusedBaselineWorkers=20,totalRoots=160,normalWorkers=160,missingRoots=0,badRootHeaders=0,eventsPerCase=1819447887,sourceRichRootsStatChecked=20,rootStatsAndHeadersChecked=160,logSHA256Checked=logs,liveReplayCompletionMarkersChecked=160,currentNativeListsChecked=140,freshTargetedQueueCount=0,queueChecks=queues,cases=reports,totalOutputBytes=sum(x['rootBytes'] for x in reports),newOutputBytes=sum(x['rootBytes'] for x in reports if x['caseId']!='ep8_mergingOff'),latestNativeEnd=dict(americaNewYork=last[0].isoformat(),caseId=last[1],period=last[2]),mailAcceptedAmericaNewYork=datetime.fromisoformat(sent['sentAtUtc']).astimezone(timezone(timedelta(hours=-5))).isoformat(),mailResent=False,jobResubmitted=False,haddExecuted=False)
 dest=Path(__file__).parent/'independent_completion_agent.json'
 with dest.open('x') as f:json.dump(result,f,indent=2)
 print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
