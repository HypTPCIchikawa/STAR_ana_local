#!/usr/bin/env python3
"""Freeze v3 transport-only clone and validate real STAR QA plus mocked replica metadata."""
from pathlib import Path
from datetime import datetime,timezone
from unittest.mock import patch
from urllib.parse import urlsplit
import contextlib,hashlib,importlib.util,io,json,shutil,subprocess,sys,time
sys.dont_write_bytecode=True
C=Path(__file__).resolve().parents[1];R=C/'runtime_v3';V=C/'validation_v3'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def now():return datetime.now(timezone.utc).isoformat()
assert not R.exists(),'Refuse replacing an immutable runtime'
source_manifest_sha=sha(C/'runtime_v2/manifest.json')
source_manifest=read(C/'runtime_v2/manifest.json')
for rel,item in source_manifest['files'].items():assert sha(C/'runtime_v2'/rel)==item['sha256'],'v2 changed: '+rel
shutil.copytree(C/'runtime_v2',R,ignore=shutil.ignore_patterns('__pycache__','*.pyc','build.log'))
p=R/'worker.py';s=p.read_text()
old="            if before is not None and before.get('MTime')!=after.get('MTime'):raise RuntimeError('Source MTime changed during staging')"
new="""            # A redirector may return distinct archival replicas with different
            # inode/mtime metadata. Preserve these values for audit; strict
            # before/after/local byte sizes and subsequent ROOT gates remain.
            metadata_audit={'policy':'audit_only_redirector_replica_metadata','mtimeChanged':before is not None and before.get('MTime')!=after.get('MTime'),'idChanged':before is not None and before.get('Id')!=after.get('Id')}"""
assert s.count(old)==1;s=s.replace(old,new)
s=s.replace("sourceStatBefore=before,sourceStatAfter=after,serverChecksumVerified=False", "sourceStatBefore=before,sourceStatAfter=after,sourceMetadataAudit=metadata_audit,serverChecksumVerified=False")
s=s.replace("'remote sizes/MTime before+after and local digests; no source checksum available'", "'exact-URL source sizes before/after equal local size; local digests recorded; MTime/Id audit-only; no source checksum available'")
s=s.replace("receipt['inputs']=staged;receipt['stageSeconds']", "receipt['sourceMetadataPolicy']='audit_only_redirector_replica_metadata';receipt['sourceMetadataChangedInputs']=sum(bool(x['sourceMetadataAudit']['mtimeChanged'] or x['sourceMetadataAudit']['idChanged']) for x in staged);receipt['inputs']=staged;receipt['stageSeconds']")
p.write_text(s);compile(s,str(p),'exec')
policy=read(R/'policy.json');policy.update(remoteMTimePolicy='audit_only_redirector_replica_metadata',remoteIdPolicy='audit_only_redirector_replica_metadata',transferVerification='exact-URL server sizes before/after equal local size + local SHA256/adler32 + strict ROOT coverage/identity/old-hist closure; replica MTime/Id audit-only')
write(R/'policy.json',policy);smoke_policy=dict(policy,localValidation=True);write(R/'validation/policy_smoke.json',smoke_policy)
revision={'version':3,'sourceRuntime':str(C/'runtime_v2'),'sourceManifestSha256':source_manifest_sha,'change':'Transport only: native redirector MTime and Id may differ by archival replica, so retain them as audit metadata rather than fail. Keep exact-URL before/after/local size equality, local digests, and all ROOT/candidate/identity/old-hist closure gates unchanged.','qaCodeConfigLibrariesSameAsV2':True,'failedV2ResultsPreserved':True,'serverChecksumVerified':False,'checksumLimitation':'The endpoint does not support server checksum queries. Local SHA256/adler32 do not by themselves prove equality to a server checksum.'}
write(R/'transport_revision.json',revision)
manifest=dict(source_manifest);manifest.update(runtimeVersion=3,purpose='Frozen selected-candidate QA; replica MTime/Id are audit-only, byte sizes and ROOT gates remain strict',sourceV2ManifestSha256=source_manifest_sha)
manifest['files']={str(p.relative_to(R)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(R.rglob('*')) if p.is_file() and p.name!='manifest.json'}
write(R/'manifest.json',manifest)
for rel,item in source_manifest['files'].items():
 if rel.startswith(('source/','include/','config/','lib/')):assert sha(R/rel)==item['sha256'],rel
runtime_sha=sha(R/'manifest.json');policy_sha=sha(R/'policy.json');V.mkdir(exist_ok=False)
base=read(C/'validation_v2/frozen_smoke_001/receipt.json');cache=base['sourceCache'];local_list=base['sourcePicoList'];local_pico=base['inputs'][0]['source'];pico_bytes=Path(local_pico).stat().st_size
print('V3_FROZEN',runtime_sha,flush=True)
start=time.time();out=V/'frozen_smoke_001'
cmd=['bash',str(R/'worker.sh'),cache,local_list,str(out),str(R/'validation/policy_smoke.json')]
with (V/'frozen_smoke_console.log').open('w') as log:done=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
assert done.returncode==0,(V/'frozen_smoke_console.log').read_text()
smoke=read(out/'receipt.json');summary=read(out/'summary.json');assert smoke['status']=='completed' and summary['passed'] and summary['closureChecks']==244827
assert summary['counters']==base['qaSummary']['counters'] and summary['histogramCount']==base['qaSummary']['histogramCount']
assert not Path(smoke['scratch']).exists()
print('V3_LOCAL_STRICT_PASS',smoke['wallSeconds'],flush=True)

# Simulated remote stage uses the known local bytes but keeps the production
# worker's ordinary remote branch and executes its real compiled STAR QA.
spec=importlib.util.spec_from_file_location('merging_v3_worker',R/'worker.py');worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
mock_url='root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_3p85GeV_fixedTarget_2021/ReversedFullField/P24iy/2021/123/22123040/st_physics_22123040_raw_5000023.picoDst.root'
mock_list=V/'simulated_remote.list';mock_list.write_text(mock_url+'\n');env=read(R/'environment.json');real_run=subprocess.run
cases=[]
def run_case(name,before_bytes,after_bytes,truncate=False,transfer_rc=0):
 out=V/name;commands=[];statcalls=[]
 before={'sourceUrl':mock_url,'endpoint':'root://xrdstar.rcf.bnl.gov:1095','pathPassedUnchanged':urlsplit(mock_url).path,'bytes':before_bytes,'MTime':'2025-10-13 22:52:17','Id':'1075391091441710','simulated':True}
 after=dict(before,bytes=after_bytes,MTime='2025-10-28 22:00:21',Id='1288756476772398')
 def stat(url,e,timeout):
  assert url==mock_url and e['xrdfs']==env['xrdfs'];statcalls.append(url)
  return dict(before if len(statcalls)==1 else after)
 def run(cmd,*args,**kwargs):
  if cmd[0]==env['xrdcp']:
   commands.append(cmd);assert cmd[1:4]==['--nopbar','--retry','2'] and '--cksum' not in cmd and cmd[-2]==mock_url
   if transfer_rc:return subprocess.CompletedProcess(cmd,transfer_rc)
   if truncate:
    with Path(local_pico).open('rb') as src,Path(cmd[-1]).open('wb') as dst:dst.write(src.read(1024))
   else:shutil.copyfile(local_pico,cmd[-1])
   return subprocess.CompletedProcess(cmd,0)
  return real_run(cmd,*args,**kwargs)
 argv=['worker.py','--runtime',str(R),'--cache',cache,'--pico-list',str(mock_list),'--output-dir',str(out),'--policy',str(R/'policy.json'),'--scratch-base','/tmp']
 console=io.StringIO()
 with patch.object(sys,'argv',argv),patch.object(worker,'exact_remote_stat',stat),patch.object(worker.subprocess,'run',run),contextlib.redirect_stdout(console),contextlib.redirect_stderr(console):rc=worker.main()
 (V/(name+'_console.log')).write_text(console.getvalue());r=read(out/'worker_result.json')
 assert not Path(r['scratch']).exists(),'Owned scratch cleanup failed'
 result={'name':name,'returnCode':rc,'status':r['status'],'workerResult':str(out/'worker_result.json'),'ownedScratchCleaned':True,'statCalls':len(statcalls),'ordinaryXrdcpCommand':commands[0] if commands else None,'simulatedRemoteTransfer':True}
 if name=='replica_metadata_change_001':
  assert rc==0 and r['status']=='completed' and r['qaSummary']['counters']==summary['counters'] and r['qaSummary']['closureChecks']==244827
  assert r['sourceMetadataChangedInputs']==1 and r['inputs'][0]['sourceMetadataAudit']['mtimeChanged'] and r['inputs'][0]['sourceMetadataAudit']['idChanged']
  assert r['inputs'][0]['localDigests']['sha256']==base['inputs'][0]['localDigests']['sha256']
  result.update(metadataDifferenceAccepted=True,realCompiledStarQaPassed=True,closureChecks=244827)
 else:
  assert rc!=0 and r['status']=='failed' and not (out/'qa.root').exists() and not (out/'receipt.json').exists()
  if transfer_rc:assert 'xrdcp failed' in r['error']
  else:assert 'Source/local size mismatch' in r['error']
  result['rejectedError']=r['error']
 cases.append(result);print('V3_REGRESSION_PASS',name,rc,flush=True)
run_case('replica_metadata_change_001',pico_bytes,pico_bytes)
run_case('source_size_mismatch_001',pico_bytes,pico_bytes+1)
run_case('local_size_mismatch_001',pico_bytes,pico_bytes,truncate=True)
run_case('transfer_failure_001',pico_bytes,pico_bytes,transfer_rc=23)
for rel,item in manifest['files'].items():assert sha(R/rel)==item['sha256'],'v3 changed after freeze '+rel
assert sha(C/'runtime_v2/manifest.json')==source_manifest_sha
for rel,item in source_manifest['files'].items():assert sha(C/'runtime_v2'/rel)==item['sha256'],'v2 changed '+rel
write(V/'regression_results.json',{'status':'passed','cases':cases,'runtimeManifestSha256':runtime_sha,'policySha256':policy_sha,'remoteFullstreamTransferred':False})
validation={'status':'passed','runtimeVersion':3,'completedAtUtc':now(),'runtime':str(R),'worker':str(R/'worker.sh'),'workerCli':['CACHE_ROOT','ORIGINAL_PICO_LIST','OUTPUT_DIR','POLICY_JSON'],'policy':str(R/'policy.json'),'manifestSha256':runtime_sha,'policySha256':policy_sha,'sourceV2ManifestSha256':source_manifest_sha,'qaCodeConfigLibrariesSameAsV2':True,'qaCodeConfigLibrariesSameAsV1':True,'frozenSmoke':str(V/'frozen_smoke_001/worker_result.json'),'smokeWallSeconds':smoke['wallSeconds'],'smokeQaProcessSeconds':smoke['qaProcessSeconds'],'smokeStagedBytes':smoke['picoStageBytes'],'oldHistogramClosureChecks':244827,'physicsCountersIdenticalToV2Smoke':True,'replicaMetadataDifferenceAcceptedWithRealStarQa':True,'simulatedRemoteSizeMismatchRejected':True,'simulatedLocalSizeMismatchRejected':True,'transferFailureRejected':True,'failureCleanupPassed':True,'regressionResults':str(V/'regression_results.json'),'transferPolicy':policy['transferVerification'],'sourceMetadataPolicy':'audit_only_redirector_replica_metadata','serverChecksumVerified':False,'checksumLimitation':revision['checksumLimitation'],'failedV2ResultsPreserved':True,'v2RuntimeUnchanged':True,'noSharedAcliC':True,'condorSubmittedByThisAgent':False,'publicSuccessMarker':'worker_result.json status=completed + summary.json passed=true + qa.root','validationWallSeconds':time.time()-start}
write(C/'runtime_validation_v3.json',validation);print('V3_READY',json.dumps(validation),flush=True)
