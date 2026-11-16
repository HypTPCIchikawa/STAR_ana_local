#!/usr/bin/env python3
"""Complete independent direct metadata validation after all four workers returned.
The initial metadata-only failure was a cleaned scratch bind; no QA rerun occurs.
"""
from pathlib import Path
import ast,json,time
p=Path(__file__).with_name('run_local_qa.py').resolve()
tree=ast.parse(p.read_text())
nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef)) or isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('C','L','R') for t in n.targets)]
ns={'__file__':str(p)};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(p),'exec'),ns)
L,R,C,read,write,sha=map(ns.__getitem__,['L','R','C','read','write','sha'])
manifest_path=L/'qa_input_manifest.json';validation_path=L/'runtime_validation.json'
m=read(manifest_path);old_validation=read(validation_path)
assert all((L/'qa'/s['id']/'worker_result.json').exists() for s in m['sources'])
for path in (manifest_path,validation_path):
 backup=path.with_name(path.stem+'_initial_metadata_failure.json')
 if not backup.exists():write(backup,read(path))
m['status']='metadata_verifying';write(manifest_path,m)
runtime=read(R/'manifest.json');runtime_sha=sha(R/'manifest.json');policy_sha=sha(R/'policy.json')
assert runtime_sha==m['runtimeManifestSha256'] and policy_sha==m['policySha256']
start=time.time();results=[]
for s in m['sources']:
 out=L/'qa'/s['id'];r=read(out/'receipt.json');summary=read(out/'summary.json');marker=read(out/'worker_result.json')
 assert r==marker and r['status']=='completed' and summary==r['qaSummary'] and summary['passed'] is True
 assert r['runtimeManifestSha256']==runtime_sha and r['policySha256']==policy_sha
 assert r['sourceCache']==s['cacheRoot'] and r['inputListSha256']==sha(s['picoList']) and sha(s['cacheRoot'])==r['sourceCacheSha256']
 assert sha(out/'qa.root')==r['outputRootSha256'] and summary['closureChecks']==244827
 meta_path=out/'root_metadata.json'
 if not meta_path.exists():
  log=out/'root_metadata_probe.log'
  if log.exists():log.rename(out/'root_metadata_probe_initial_failure.log')
  meta=ns['probe'](s,r,summary);write(meta_path,meta)
 meta=read(meta_path)
 assert meta['passed'] is True and meta['runStatus']=='completed' and meta['uuid']==meta['storedFileUuid']
 assert meta['sourceCacheSha256']==r['sourceCacheSha256']
 assert meta['inputEvents']==meta['eventLedgerEntries']==meta['picoEntries']==s['sourceInputEventsMeasured']
 assert meta['reconstructedEvents']==meta['cacheEntries']==s['acceptedCacheEvents']==summary['counters']['cache_events']
 assert meta['selectedLambdaCandidates']==meta['selectedLambdaEntries']==meta['cachedLambdaCandidates']==summary['counters']['cache_lambda_candidates']==s['expectedLambda']
 s.update(qaOutcome='completed',exclusionReason=None,rootMetadata=str(meta_path),rootMetadataSha256=sha(meta_path),receipt=str(out/'receipt.json'),receiptSha256=sha(out/'receipt.json'),qaPicoEventsRead=summary['counters']['pico_events_read'],selectedHe4Pairs=summary['counters'].get('signal_pairs_selected_4He',0),selectedHe3Pairs=summary['counters'].get('signal_pairs_selected_3He',0))
 results.append({'id':s['id'],'inputEvents':meta['inputEvents'],'accepted':meta['cacheEntries'],'lambdaCandidates':meta['selectedLambdaCandidates'],'picoRead':s['qaPicoEventsRead'],'He4':s['selectedHe4Pairs'],'He3':s['selectedHe3Pairs'],'wallSeconds':r['wallSeconds'],'qaProcessSeconds':r['qaProcessSeconds'],'receipt':s['receipt'],'rootMetadata':str(meta_path)})
 print('LOCAL_STRICT_PASS',json.dumps(results[-1]),flush=True)
 write(manifest_path,m)
for rel,item in runtime['files'].items():assert sha(R/rel)==item['sha256'],'Local runtime changed: '+rel
original=read(C/'runtime_v2/manifest.json');assert sha(C/'runtime_v2/manifest.json')==runtime['sourceRuntimeManifestSha256']
for rel,item in original['files'].items():
 if rel.startswith(('source/','include/','config/','lib/')):assert sha(R/rel)==item['sha256'],'QA implementation differs: '+rel
v=dict(old_validation)
v.update(status='passed',qaSourceIdsPassed=[s['id'] for s in m['sources']],qaSourceIdsExcluded=[],sourceInputEventsRepresented=sum(s['sourceInputEventsMeasured'] for s in m['sources']),acceptedCacheEvents=sum(s['acceptedCacheEvents'] for s in m['sources']),qaPicoEventsRead=sum(s['qaPicoEventsRead'] for s in m['sources']),selectedHe4Pairs=sum(s['selectedHe4Pairs'] for s in m['sources']),selectedHe3Pairs=sum(s['selectedHe3Pairs'] for s in m['sources']),completedAtUtc=ns['now'](),metadataProbeRecovery={'cause':'Initial post-worker probe reused already-cleaned private scratch bind','fix':'Metadata-only probe uses durable local-extension directory and omits expired scratch bind','qaWorkersRerun':False,'qaCodeOrFrozenRuntimeChanged':False,'recoverySeconds':time.time()-start},sourceResults=results,qaWorkerWallSecondsSum=sum(x['wallSeconds'] for x in results),qaProcessSecondsSum=sum(x['qaProcessSeconds'] for x in results),runtimeUnchangedAfterQa=True,qaCodeConfigLibrariesUnchanged=True)
write(validation_path,v);m.update(status='final',completedAtUtc=ns['now'](),runtimeValidationSha256=sha(validation_path));write(manifest_path,m)
print('LOCAL_FINAL_PASS',json.dumps({k:v[k] for k in ['status','sourceInputEventsRepresented','acceptedCacheEvents','qaPicoEventsRead','selectedHe4Pairs','selectedHe3Pairs','qaWorkerWallSecondsSum','qaProcessSecondsSum']}),flush=True)
