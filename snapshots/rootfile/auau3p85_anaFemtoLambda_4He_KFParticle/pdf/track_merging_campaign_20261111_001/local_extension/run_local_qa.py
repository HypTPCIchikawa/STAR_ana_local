#!/usr/bin/env python3
"""Independent, immutable local-only QA copy; never submits a farm job."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,shutil,subprocess,time

C=Path(__file__).resolve().parents[1]
L=C/'local_extension'
R=L/'runtime'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(p)
def read(p):return json.loads(Path(p).read_text())

def probe(source,receipt,summary):
 # ROOT built-in metadata reads only; no event loop or physics calculation.
 # The worker has already removed its private scratch. Use the durable
 # local-extension directory and omit that now-expired scratch bind.
 raw=receipt['qaCommand']
 idx=raw.index('-l');cmd=raw[:idx]+['-l','-b']
 cmd[cmd.index('--pwd')+1]=str(L)
 old_bind=receipt['scratch']+':'+receipt['scratch']
 for j in range(len(cmd)-1,0,-1):
  if cmd[j]==old_bind and cmd[j-1]=='-B':del cmd[j-1:j+1]
 cache=source['cacheRoot'];pico=source['inputPath']
 q=lambda x:json.dumps(str(x))
 script='\n'.join([
  '#include "TFile.h"','#include "TTree.h"','#include "TParameter.h"','#include "TObjString.h"','#include "TNamed.h"',
  'TFile *f=TFile::Open('+q(cache)+',"READ");',
  'if(!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) {gSystem->Exit(11);}',
  'TNamed *status=dynamic_cast<TNamed*>(f->Get("FemtoLambdaRunStatus"));',
  'TObjString *uuid=dynamic_cast<TObjString*>(f->Get("FemtoCandidateEventStoreFileUUID"));',
  'TTree *e=dynamic_cast<TTree*>(f->Get("eventLedger"));',
  'TTree *c=dynamic_cast<TTree*>(f->Get("femtoCandidateEvents"));',
  'TTree *l=dynamic_cast<TTree*>(f->Get("selectedLambda"));',
  'TParameter<Long64_t> *ni=dynamic_cast<TParameter<Long64_t>*>(f->Get("inputEvents"));',
  'TParameter<Long64_t> *nr=dynamic_cast<TParameter<Long64_t>*>(f->Get("reconstructedEvents"));',
  'TParameter<Long64_t> *nl=dynamic_cast<TParameter<Long64_t>*>(f->Get("selectedLambdaCandidates"));',
  'if(!status || !uuid || !e || !c || !l || !ni || !nr || !nl) {gSystem->Exit(12);}',
  'TFile *p=TFile::Open('+q(pico)+',"READ");',
  'if(!p || p->IsZombie() || p->TestBit(TFile::kRecovered)) {gSystem->Exit(13);}',
  'TTree *pt=dynamic_cast<TTree*>(p->Get("PicoDst"));',
  'if(!pt) {gSystem->Exit(14);}',
  'std::cout<<"LOCAL_ROOT_METADATA|"<<f->GetUUID().AsString()<<"|"<<uuid->GetString().Data()<<"|"<<status->GetTitle()<<"|"<<ni->GetVal()<<"|"<<e->GetEntries()<<"|"<<pt->GetEntries()<<"|"<<nr->GetVal()<<"|"<<c->GetEntries()<<"|"<<nl->GetVal()<<"|"<<l->GetEntries()<<std::endl;',
  'gSystem->Exit(0);',''])
 done=subprocess.run(cmd,input=script,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
 out=L/'qa'/source['id']/'root_metadata_probe.log';out.write_text(done.stdout)
 match=re.search(r'^LOCAL_ROOT_METADATA\|([^\n\r]+)$',done.stdout,re.M)
 if done.returncode or not match:raise RuntimeError('Direct ROOT metadata probe failed; see '+str(out))
 vals=match.group(1).split('|')
 if len(vals)!=10:raise RuntimeError('Malformed ROOT metadata probe output: '+repr(vals))
 meta={'schemaVersion':1,'passed':False,'sourceCache':cache,'sourceCacheSha256':receipt['sourceCacheSha256'],'uuid':vals[0],'storedFileUuid':vals[1],'runStatus':vals[2],'nonZombie':True,'nonRecovered':True,'picoNonZombie':True,'picoNonRecovered':True,'inputEvents':int(vals[3]),'eventLedgerEntries':int(vals[4]),'picoEntries':int(vals[5]),'reconstructedEvents':int(vals[6]),'cacheEntries':int(vals[7]),'selectedLambdaCandidates':int(vals[8]),'selectedLambdaEntries':int(vals[9]),'cachedLambdaCandidates':summary['counters'].get('cache_lambda_candidates',0),'oldHistogramClosureChecks':summary['closureChecks'],'evidence':'Direct ROOT metadata reads plus strict frozen worker PASS, complete original input Entries, identity checks for all QA-read events, full cached-frame decode and exact old histogram closure','probeLog':str(out),'probeLogSha256':sha(out),'verifiedAtUtc':now()}
 assert meta['uuid']==meta['storedFileUuid'] and meta['runStatus']=='completed'
 assert meta['inputEvents']==meta['eventLedgerEntries']==meta['picoEntries']==source['sourceInputEventsMeasured']
 assert meta['cacheEntries']==meta['reconstructedEvents']==source['acceptedCacheEvents']==summary['counters']['cache_events']
 assert meta['selectedLambdaCandidates']==meta['selectedLambdaEntries']==meta['cachedLambdaCandidates']==source['expectedLambda']
 assert summary['passed'] is True and meta['oldHistogramClosureChecks']==244827
 meta['passed']=True
 return meta

assert not R.exists(),'Do not replace an existing frozen local runtime'
status=read(L/'producer_status.json')
assert status['phase']=='producer_stage_finished_pending_strict_QA' and status['completedProcessCount']==3
baseline=read(C/'validation_v2/frozen_smoke_001/receipt.json')
baseline_list=L/'input_lists/0.list';baseline_list.parent.mkdir(exist_ok=True)
baseline_list.write_text(baseline['inputs'][0]['source']+'\n')
sources=[{'id':'0','cacheRoot':baseline['sourceCache'],'picoList':str(baseline_list),'inputPath':baseline['inputs'][0]['source'],'sourceInputEventsMeasured':4027,'acceptedCacheEvents':3858,'expectedLambda':62,'kind':'baseline','producerResult':None,'qaOutcome':'pending','exclusionReason':None}]
allow=read(C/'runtime_v2/validation/local_allowlist.json')
producer_checks=[]
for i,expected in ((1,34603),(2,34841),(3,13140)):
 result_path=L/f'producer/{i}/result.json';r=read(result_path)
 log=L/f'producer/{i}/producer.log';s=log.read_text(errors='replace')
 final=re.findall(r'FemtoLambda job=\S+ requested=(\d+) completed=(\d+) reconstructed=(\d+) status=(\d+)',s)
 row=re.findall(r'\[StFemtoMaker Lambda\] inputEvents=(\d+) reconstructedEvents=(\d+) selectedLambda=(\d+) mixingVzClamped=\d+ status=(\S+)',s)
 assert len(final)==len(row)==1,('Ambiguous producer completion',i)
 requested,completed,reconstructed,code=map(int,final[0]);inputs,accepted,lam=map(int,row[0][:3])
 passed=r['phase']=='producer_returned_zero_pending_strict_QA' and r['returnCode']==0 and r['inputUnchanged'] is True and r['producerLibrariesUnchanged'] is True and requested==completed==inputs==expected and reconstructed==accepted and code==0 and row[0][3]=='completed'
 check={'id':str(i),'passed':passed,'requested':requested,'completed':completed,'reconstructedEvents':reconstructed,'selectedLambda':lam,'statusCode':code,'runStatus':row[0][3],'result':str(result_path),'resultSha256':sha(result_path),'producerLog':str(log),'producerLogSha256':sha(log),'wallSeconds':r['wallSeconds']}
 producer_checks.append(check)
 source={'id':str(i),'cacheRoot':r['candidateRoot'],'picoList':r['picoList'],'inputPath':r['input']['path'],'sourceInputEventsMeasured':expected,'acceptedCacheEvents':accepted,'expectedLambda':lam,'kind':'additional','producerResult':str(result_path),'qaOutcome':'pending' if passed else 'failed','exclusionReason':None if passed else 'Producer failed full input completion/integrity gate'}
 sources.append(source)
 if passed:
  assert sha(source['inputPath'])==r['input']['sha256'],'Producer input changed'
  allow[source['inputPath']]=r['input']['sha256'];allow[source['cacheRoot']]=sha(source['cacheRoot'])
write(L/'producer_completion_validation.json',{'status':'passed' if all(x['passed'] for x in producer_checks) else 'partial','checks':producer_checks,'checkedAtUtc':now()})
shutil.copytree(C/'runtime_v2',R,ignore=shutil.ignore_patterns('__pycache__','*.pyc','build.log'))
policy=read(R/'policy.json');policy['localValidation']=True
policy['scope']='Independent local preview only; four exact local sources, no farm cohort, no remote transfer'
write(R/'policy.json',policy);write(R/'validation/local_allowlist.json',allow)
write(R/'local_copy_revision.json',{'sourceRuntime':str(C/'runtime_v2'),'sourceManifestSha256':sha(C/'runtime_v2/manifest.json'),'qaCodeConfigLibrariesUnchanged':True,'copiedAtUtc':now(),'change':'Local validation policy and exact SHA256 allowlist only. Analysis source/config/libraries unchanged. No Condor submission.'})
manifest=read(C/'runtime_v2/manifest.json');manifest['scope']='Independent local extension preview';manifest['sourceRuntimeManifestSha256']=sha(C/'runtime_v2/manifest.json')
manifest['files']={str(p.relative_to(R)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(R.rglob('*')) if p.is_file() and p.name!='manifest.json'}
write(R/'manifest.json',manifest)
for rel,old in read(C/'runtime_v2/manifest.json')['files'].items():
 if rel.startswith(('source/','include/','config/','lib/')):assert sha(R/rel)==old['sha256'],rel
runtime_sha=sha(R/'manifest.json');policy_sha=sha(R/'policy.json')
manifest_path=L/'qa_input_manifest.json';validation_path=L/'runtime_validation.json'
input_manifest={'schemaVersion':1,'status':'running','runtimeValidation':str(validation_path),'runtime':str(R),'policy':str(R/'policy.json'),'runtimeManifestSha256':runtime_sha,'policySha256':policy_sha,'sourceCoverageIsNotQaPicoReadCount':True,'farmCohortOverlapLfnCount':0,'sources':sources}
write(manifest_path,input_manifest)
write(validation_path,{'status':'running','runtime':str(R),'manifestSha256':runtime_sha,'policy':str(R/'policy.json'),'policySha256':policy_sha,'qaCodeConfigLibrariesUnchanged':True,'sourceRuntimeValidation':str(C/'runtime_validation_v2.json'),'producerCompletionValidation':str(L/'producer_completion_validation.json')})
print('LOCAL_RUNTIME_FROZEN',R,runtime_sha,flush=True)
started=time.time()
for source in sources:
 if source['qaOutcome']=='failed':continue
 i=source['id'];out=L/'qa'/i;out.mkdir(parents=True,exist_ok=False)
 cmd=['bash',str(R/'worker.sh'),source['cacheRoot'],source['picoList'],str(out),str(R/'policy.json')]
 try:
  with (out/'worker_console.log').open('w') as log:
   done=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=600)
  assert done.returncode==0,'Worker return code '+str(done.returncode)
  receipt=read(out/'receipt.json');marker=read(out/'worker_result.json');summary=read(out/'summary.json')
  assert receipt==marker and receipt['status']=='completed' and summary==receipt['qaSummary'] and summary['passed'] is True
  assert receipt['runtimeManifestSha256']==runtime_sha and receipt['policySha256']==policy_sha
  assert receipt['sourceCacheSha256']==allow[source['cacheRoot']] and receipt['inputListSha256']==sha(source['picoList'])
  assert sha(out/'qa.root')==receipt['outputRootSha256']
  meta=probe(source,receipt,summary);write(out/'root_metadata.json',meta)
  source.update(qaOutcome='completed',rootMetadata=str(out/'root_metadata.json'),rootMetadataSha256=sha(out/'root_metadata.json'),receipt=str(out/'receipt.json'),receiptSha256=sha(out/'receipt.json'),qaPicoEventsRead=summary['counters']['pico_events_read'],selectedHe4Pairs=summary['counters'].get('signal_pairs_selected_4He',0),selectedHe3Pairs=summary['counters'].get('signal_pairs_selected_3He',0))
  print('LOCAL_QA_PASS',i,'input',meta['inputEvents'],'accepted',meta['cacheEntries'],'picoRead',source['qaPicoEventsRead'],'He4',source['selectedHe4Pairs'],'He3',source['selectedHe3Pairs'],'wall',receipt['wallSeconds'],flush=True)
 except Exception as e:
  source.update(qaOutcome='failed',exclusionReason=str(e));write(out/'local_validation_failure.json',{'passed':False,'error':repr(e),'atUtc':now()})
  print('LOCAL_QA_FAILED',i,repr(e),flush=True)
 write(manifest_path,input_manifest)
for rel,item in manifest['files'].items():assert sha(R/rel)==item['sha256'],'Runtime changed during local QA '+rel
passed=[s for s in sources if s['qaOutcome']=='completed']
validation={'status':'passed' if len(passed)==len(sources) else 'partial','runtime':str(R),'manifestSha256':runtime_sha,'policy':str(R/'policy.json'),'policySha256':policy_sha,'sourceRuntimeManifestSha256':manifest['sourceRuntimeManifestSha256'],'sourceRuntimeValidation':str(C/'runtime_validation_v2.json'),'qaCodeConfigLibrariesUnchanged':True,'runtimeUnchangedAfterQa':True,'producerCompletionValidation':str(L/'producer_completion_validation.json'),'qaSourceIdsPassed':[s['id'] for s in passed],'qaSourceIdsExcluded':[s['id'] for s in sources if s not in passed],'oldHistogramClosureChecksPerSource':244827,'sourceInputEventsRepresented':sum(s['sourceInputEventsMeasured'] for s in passed),'acceptedCacheEvents':sum(s['acceptedCacheEvents'] for s in passed),'qaPicoEventsRead':sum(s['qaPicoEventsRead'] for s in passed),'selectedHe4Pairs':sum(s['selectedHe4Pairs'] for s in passed),'selectedHe3Pairs':sum(s['selectedHe3Pairs'] for s in passed),'wallSeconds':time.time()-started,'completedAtUtc':now(),'condorSubmitted':False,'remoteFullstreamTransferred':False}
write(validation_path,validation)
input_manifest.update(status='final',completedAtUtc=now(),runtimeValidationSha256=sha(validation_path));write(manifest_path,input_manifest)
print('LOCAL_QA_FINAL',json.dumps(validation),flush=True)
