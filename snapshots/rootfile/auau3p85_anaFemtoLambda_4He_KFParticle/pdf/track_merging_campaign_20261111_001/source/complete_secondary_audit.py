#!/usr/bin/env python3
"""Join strict ROOT metadata with the separate nine-file inventory; leave primary inputs frozen."""
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();out=a.study.resolve()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
primary=json.loads((out/'input_inventory.json').read_text());original_validation=json.loads((out/'input_inventory_validation.json').read_text());original_secondary=json.loads((out/'secondary_header_consistent_unconfirmed.json').read_text());rootmeta=json.loads((out/'secondary_root_metadata.json').read_text())
for name in ['input_inventory.json','candidates.tsv','excluded_candidates.tsv','secondary_header_consistent_unconfirmed.json']:
 assert sha(out/name)==original_validation['hashes'][name],('Frozen primary artifact changed',name)
expected=[15642,16558,16648,16655,16732,16746,16753,16768,16988]
assert [r['subjob'] for r in original_secondary['rows']]==expected
assert [r['subjob'] for r in rootmeta['samples']]==expected and rootmeta['sampleCount']==rootmeta['strictMetadataValidCount']==9
original_by_id={r['subjob']:r for r in original_secondary['rows']};primary_rows=[r for r in primary['rows'] if r['preflight_eligible']]
primary_urls={i['url'] for r in primary_rows for i in r['inputs']};primary_lfns={u.rsplit('/',1)[-1] for u in primary_urls};seen_urls=set();seen_lfns=set();records=[]
for m in rootmeta['samples']:
 r=dict(original_by_id[m['subjob']]);fp=Path(r['cache_abs']).stat();assert (fp.st_size,fp.st_mtime_ns,fp.st_ino,fp.st_dev)==(r['size'],r['mtime_ns'],r['inode'],r['device'])
 assert r['header_status']=='consistent_fend' and r['root_header_end']==r['size']
 assert m['strictMetadataValid'] and not m['isZombie'] and not m['isRecovered'] and m['producerStatus']=='completed' and m['uuidMatches'] and m['countsMatch']
 assert all(m['metadataMatchesReference'].values()) and m['mode']==rootmeta['reference']['mode']
 q=m['parameters'];t=m['treeEntries'];assert q['inputEvents']==t['eventLedger'] and q['reconstructedEvents']==t['femtoCandidateEvents'] and q['selectedLambdaCandidates']==t['selectedLambda']
 assert sha(Path(r['original_list_abs']))==r['list_sha256'] and r['dataset_verified'] and r['input_count']==3
 urls=[i['url'] for i in r['inputs']];lfns=[u.rsplit('/',1)[-1] for u in urls]
 assert all('/production_3p85GeV_fixedTarget_2021/' in u and '/P24iy/2021/' in u for u in urls)
 assert len(set(urls))==len(set(lfns))==3 and not (set(urls)&(primary_urls|seen_urls)) and not (set(lfns)&(primary_lfns|seen_lfns))
 seen_urls.update(urls);seen_lfns.update(lfns)
 r.update(uuid=m['fileUuid'],producer_root_metadata_verified=True,root_metadata_verified=True,requires_worker_root_metadata_verification=True,input_events=q['inputEvents'],accepted_events=q['reconstructedEvents'],selected_lambdas=q['selectedLambdaCandidates'],counter_source='ROOT TParameter and matching TTree metadata',no_stdoutinfo=True,no_final_stdout_completion_marker=True,stdout_status='final_completion_marker_missing',strict_root_metadata=m,conditional_admission_basis='Completed original ROOT producer, nonRecovered, exact UUID/provenance and tree/counter agreement; must still pass full campaign worker validation.',primary_dispatch_member=False,potentially_admissible=True)
 records.append(r)
assert len({r['uuid'] for r in records})==9 and rootmeta['reference']['uuid'] not in {r['uuid'] for r in records}
with (out/'candidates.tsv').open() as f:fields=next(csv.reader(f,delimiter='\t'))
with (out/'secondary_candidates.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
 for r in records:
  d=dict(r,urls_json=json.dumps([i['url'] for i in r['inputs']],separators=(',',':')));w.writerow([d.get(k,'') for k in fields])
summary={'sample_count':9,'strict_root_metadata_valid_count':9,'additional_unique_pico_urls':len(seen_urls),'additional_unique_lfns':len(seen_lfns),'url_overlap_with_primary_or_self':0,'lfn_overlap_with_primary_or_self':0,'input_events_from_root_metadata':sum(r['input_events'] for r in records),'accepted_events_from_root_metadata':sum(r['accepted_events'] for r in records),'selected_lambdas_from_root_metadata':sum(r['selected_lambdas'] for r in records),'combined_jobs_if_root_admits_secondary':len(primary_rows)+len(records),'combined_unique_pico_files_if_admitted':len(primary_urls)+len(seen_urls),'combined_measured_input_events_if_admitted':primary['summary']['eligible_measured_input_events_from_producer_logs']+sum(r['input_events'] for r in records),'combined_accepted_events_if_admitted':primary['summary']['eligible_accepted_events_from_producer_logs']+sum(r['accepted_events'] for r in records)}
source_hashes={name:sha(out/name) for name in ['secondary_root_metadata.json','secondary_header_consistent_unconfirmed.json','secondary_metadata_paths.tsv','secondary_metadata_probe_command.json','secondary_metadata_probe.log','source/ProbeSecondaryMetadata.C','source/complete_secondary_audit.py','secondary_candidates.tsv']}
reference_metadata_hashes={k:hashlib.sha256(v.encode()).hexdigest() for k,v in rootmeta['reference']['metadata'].items()}
result={'status':'potentially_admissible_completed_producers_pending_full_worker_validation','createdAtUtc':datetime.now(timezone.utc).isoformat(),'dataset':'production_3p85GeV_fixedTarget_2021 / P24iy','producer_jobid':'2AF731F2C75544890B4E858354DDE43D','rootVersion':rootmeta['rootVersion'],'summary':summary,'reference':{'file':rootmeta['reference']['file'],'uuid':rootmeta['reference']['uuid'],'mode':rootmeta['reference']['mode'],'metadata_sha256':reference_metadata_hashes},'primary_inventory_modified':False,'primary_candidates_tsv_sha256_unchanged':sha(out/'candidates.tsv'),'no_dispatch_performed':True,'candidatePayloadDecoded':False,'picoPayloadRead':False,'workerRequirements':['Open source READ and recheck nonZombie/nonRecovered, producer completed, UUID, schema, species, selection/mixing provenance and all scalar/tree-count invariants.','Decode every candidate event and match event/run/input/track references against original ordered2021P24iyPicoDsts.','Read original PicoDst Entries and require exact coverage sum=inputEvents; never substitute catalog event counts.','Require exact old merging-histogram closure and all campaign completion/output checks before publishing.','These nine have no final stdout completion marker; counter source is ROOT metadata, so they must remain explicitly distinguished in campaign provenance.'],'sourceHashes':source_hashes,'rows':records}
(out/'secondary_metadata_probe.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
with (out/'secondary_candidates.tsv').open() as f:readback=list(csv.DictReader(f,delimiter='\t'))
assert list(readback[0])==fields and [int(r['subjob']) for r in readback]==expected
assert all(json.loads(a['urls_json'])==[i['url'] for i in b['inputs']] for a,b in zip(readback,records))
validation={'passed':True,'createdAtUtc':datetime.now(timezone.utc).isoformat(),'samples':9,'allStrictRootMetadataValid':True,'allInputListsExactThree2021P24iy':True,'allTreeCountsMatchProducerCounters':True,'allCacheUuidsOriginalAndDistinct':True,'allProvenanceMatchesValidatedReference':True,'allInputFingerprintsUnchanged':True,'allPrimaryArtifactsUnchanged':True,'primaryAndSecondaryTsvColumnsIdentical':True,'secondaryTsvReadbackMatches':True,'primaryAndSelfUrlOverlap':0,'primaryAndSelfLfnOverlap':0,'fullWorkerValidationStillRequired':True,'hashes':{name:sha(out/name) for name in ['secondary_metadata_probe.json','secondary_candidates.tsv','secondary_root_metadata.json','source/ProbeSecondaryMetadata.C','source/complete_secondary_audit.py']}}
(out/'secondary_metadata_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print(json.dumps(summary,indent=2));print('SECONDARY_METADATA_VALIDATION_PASS')
