#!/usr/bin/env python3
"""Read-only independent audit of once-submitted matrix; never submits or retries."""
import pathlib,json,re,hashlib,zipfile,time,datetime,os,sys
E=pathlib.Path(__file__).resolve().parent;R=E.parents[1]
M=json.loads((E/'production_matrix.json').read_text());verified={};digest_cache={}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):
 p=pathlib.Path(path);s=p.stat();key=(str(p),s.st_size,s.st_mtime_ns)
 if key not in digest_cache:
  h=hashlib.sha256()
  with p.open('rb') as f:
   while True:
    data=f.read(1024*1024)
    if not data:break
    h.update(data)
  digest_cache[key]=h.hexdigest()
 return digest_cache[key]
def load(p):return json.loads(pathlib.Path(p).read_text())
def save(v):
 p=E/'independent_submission_audit.json';t=pathlib.Path(str(p)+'.next');t.write_text(json.dumps(v,indent=2)+'\n');os.replace(t,p)
def check(case,state):
 d=pathlib.Path(case['launchDirectory']);c=load(case['config']);jobid=state['jobid'];assert re.fullmatch('[0-9A-Fa-f]{32}',jobid)
 assert state['phase']=='submitted_verified' and state['actualJobs']==state['sourceCount']==20
 assert state['submitReturnCode']==0 and state['possibleReportIds']==[jobid]
 assert state['anaName']==c['anaName'] and state['submitName']==c['submitName']
 text=(d/'submit.log').read_text();ids=set(re.findall(re.escape(c['submitName'])+'([0-9A-Fa-f]{32})'+r'\.report',text));assert ids=={jobid} and 'Scheduling successful' in text
 assert sha(c['inputList'])==c['inputListSha256'];inputs=[str(pathlib.Path(x).resolve()) for x in pathlib.Path(c['inputList']).read_text().splitlines() if x.strip()]
 assert len(inputs)==len(set(inputs))==20
 for row in c['sources']:
  s=pathlib.Path(row['path']).stat();assert s.st_size==row['bytes'] and s.st_ino==row['inode'] and s.st_dev==row['device'] and s.st_mtime_ns==row['mtimeNs']
 for path,expected in c['fingerprints'].items():assert sha(path)==expected,('pin',path)
 coverage=load(d/'input_coverage.json');assert coverage['passed'] and coverage['actualJobs']==coverage['expectedJobs']==20 and coverage['filesPerWorker']==1
 assert coverage['inputListSha256']==c['inputListSha256']
 assert [row['subjob'] for row in coverage['orderedInputs']]==list(range(20))
 for row,expected in zip(coverage['orderedInputs'],inputs):
  lines=[x.split() for x in pathlib.Path(row['nativeList']).read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
  assert len(lines)==1 and len(lines[0]) in (1,2) and str(pathlib.Path(lines[0][0]).resolve())==row['input']==expected
 artifacts=load(d/'submission_artifacts.json')
 for kind,record in artifacts.items():
  p=pathlib.Path(record['path']);assert p.stat().st_size==record['bytes'] and sha(p)==record['sha256'],('artifact',kind)
 meta=load(state['runmeta']);assert meta['jobid']==jobid and meta['anaName']==c['submitName'] and meta['jobPrefix']==c['submitName']+jobid
 assert meta['submit']['mainconf']==c['mainconf'] and str(pathlib.Path(meta['projectRoot']).resolve())==str(R)
 assert artifacts['runmeta']['path']==state['runmeta']
 assert artifacts['joblistlog']['sha256']==c['fingerprints'][c['joblist']]
 shipped={'analysis','lib','config','include','StMaker','StRoot'};count=0
 with zipfile.ZipFile(artifacts['runtimeZip']['path']) as frozen:
  for path,expected in c['fingerprints'].items():
   relative=pathlib.Path(path).relative_to(R)
   if relative.parts[0] in shipped:assert hashlib.sha256(frozen.read(relative.as_posix())).hexdigest()==expected;count+=1
 assert count==artifacts['runtimeZip']['pinnedShippedMembersVerified']
 native=load(d/'native_resources.json');assert native['descriptions']
 for desc in native['descriptions']:
  fields={}
  for line in pathlib.Path(desc['path']).read_text().splitlines():
   m=re.match(r'\s*([+A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$',line)
   if m:fields[m.group(1).lower()]=m.group(2)
  for key,expected in c['nativeResources'].items():assert float(fields[key])==expected==float(desc['resources'][key])
  assert fields['log'].strip('"')==state['eventLog']
 receipt=native['queue'];ads=receipt['ads'];assert receipt['queueJobCount']==len(ads)==20
 seen=set()
 for ad in ads:
  m=re.search(re.escape(jobid)+r'_(\d+)\.',ad['Out']);assert m
  index=int(m.group(1));assert index not in seen;seen.add(index)
  assert ad['RequestMemory']==4096 and ad['RequestDisk']==65536000
 assert seen==set(range(20))
 return dict(caseId=case['caseId'],jobid=jobid,actualJobs=20,sourceCount=20,pinnedCurrentFiles=len(c['fingerprints']),pinnedShippedMembers=count,artifacts=len(artifacts),nativeDescriptions=len(native['descriptions']),queueReceiptJobs=20,queueReceiptCheckedAtUTC=receipt['checkedAtUTC'],resources=c['nativeResources'],runmeta=state['runmeta'],eventLog=state['eventLog'],checkedAtUTC=now(),passed=True)
start=time.monotonic();last=None
while time.monotonic()-start<1800:
 status=load(E/'matrix_launch_state.json')
 try:
  for case in M['cases']:
   if case['caseId'] in verified:continue
   path=pathlib.Path(case['launchDirectory'])/'launch_state.json'
   if not path.exists():continue
   state=load(path)
   if state['phase']=='attention_required':raise RuntimeError((case['caseId'],state))
   if state['phase']=='submitted_verified':
    verified[case['caseId']]=check(case,state);print(now(),'VERIFIED',case['caseId'],verified[case['caseId']]['jobid'],flush=True)
  if status['phase']=='attention_required':raise RuntimeError(status)
  complete=len(verified)==7
  out=dict(passed=complete,phase='complete' if complete else 'monitoring',checkedAtUTC=now(),matrixPhase=status['phase'],verifiedCases=len(verified),actualJobs=sum(x['actualJobs'] for x in verified.values()),cases=verified,readOnly=True,noSubmissionOrRetry=True)
  if complete:
   assert len(set(x['jobid'] for x in verified.values()))==7
   assert status['phase'] in ['watching','finished'] and status['actualJobs']==140
   assert {r['caseId']:r['jobid'] for r in status['cases']}=={k:v['jobid'] for k,v in verified.items()}
   settings=load(E/'watch_matrix_settings.json');out['watchSettings']=dict(path=str(E/'watch_matrix_settings.json'),sha256=sha(E/'watch_matrix_settings.json'),settings=settings)
   aggregate=E/'aggregate_watch_001';out['aggregateWatchDirectoryExists']=aggregate.is_dir();out['aggregateWatchFiles']=sorted(p.name for p in aggregate.iterdir()) if aggregate.is_dir() else []
   save(out);print(now(),'ALL_7_VERIFIED actualJobs=140',flush=True);sys.exit(0)
  save(out)
  label=(len(verified),status['phase'])
  if label!=last:print(now(),'PROGRESS',label,flush=True);last=label
 except Exception as e:
  save(dict(passed=False,phase='attention_required',checkedAtUTC=now(),verifiedCases=len(verified),cases=verified,error=str(e),readOnly=True,noSubmissionOrRetry=True));raise
 time.sleep(30)
save(dict(passed=False,phase='monitor_timeout',checkedAtUTC=now(),cases=verified,readOnly=True,noSubmissionOrRetry=True));sys.exit(2)
