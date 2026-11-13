#!/usr/bin/env python3
"""Inventory preserved candidate jobs without decoding candidate baskets or copying PicoDsts."""
import argparse,csv,hashlib,json,os,re,struct,tarfile,time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--project',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--ana-name',required=True);p.add_argument('--jobid',required=True);p.add_argument('--threads',type=int,required=True);p.add_argument('--stdout-tail-bytes',type=int,required=True);p.add_argument('--expected-files-per-job',type=int,required=True);a=p.parse_args()
 project=a.project.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
 rootdir=project/'rootfile'/a.ana_name;rootprefix=a.ana_name+'_'+a.jobid+'_';listprefix=a.ana_name+a.jobid+'_';logdir=project/'log'/a.ana_name
 rootpat=re.compile(re.escape(rootprefix)+r'(\d+)\.root$');rows=[]
 with os.scandir(rootdir) as entries:
  for entry in entries:
   match=rootpat.fullmatch(entry.name)
   if not match or not entry.is_file():continue
   st=entry.stat();rows.append({'subjob':int(match.group(1)),'cache_abs':str(Path(entry.path).resolve()),'size':st.st_size,'mtime_ns':st.st_mtime_ns,'inode':st.st_ino,'device':st.st_dev})
 rows.sort(key=lambda row:row['subjob'])
 marker=re.compile(r'FemtoLambda job='+re.escape(a.jobid)+r'_(\d+) requested=(\d+) completed=(\d+) reconstructed=(\d+) status=(\d+) RealTime=([0-9.]+) CpuTime=([0-9.]+)')
 counters=re.compile(r'\[StFemtoMaker Lambda\] inputEvents=(\d+) reconstructedEvents=(\d+) selectedLambda=(\d+) mixingVzClamped=(\d+) status=(\w+)')
 def inspect(row):
  row=dict(row);reasons=[];row['original_list_abs']=str(project/'job/run'/(listprefix+str(row['subjob'])+'.list'));row['producer_log_abs']=str(logdir/('stdout.'+a.jobid+'_'+str(row['subjob'])+'.out'));row['uuid']=None;row['producer_root_metadata_verified']=False
  if row['size']<16:row['header_status']='empty_or_short';reasons.append('root_empty_or_short')
  else:
   try:
    with open(row['cache_abs'],'rb') as f:header=f.read(24)
    magic=header[:4];version=struct.unpack('>I',header[4:8])[0];begin=struct.unpack('>I',header[8:12])[0]
    end=struct.unpack('>Q' if version>=1000000 else '>I',header[12:20] if version>=1000000 else header[12:16])[0]
    row.update(root_header_version=version,root_header_begin=begin,root_header_end=end)
    good=magic==b'root' and begin>=64 and end==row['size'];row['header_status']='consistent_fend' if good else 'inconsistent_magic_or_fend'
    if not good:reasons.append('root_header_inconsistent')
   except Exception as e:row['header_status']='read_error';row['header_error']=str(e);reasons.append('root_header_read_error')
  try:
   data=Path(row['original_list_abs']).read_bytes();row['list_sha256']=hashlib.sha256(data).hexdigest();row['list_bytes']=len(data);inputs=[]
   for number,line in enumerate(data.decode().splitlines(),1):
    if not line.strip():continue
    fields=line.split()
    if len(fields) not in (1,2):raise ValueError('bad row '+str(number))
    url=fields[0]
    if not (url.startswith('root://xrdstar.rcf.bnl.gov:1095//') and '/production_3p85GeV_fixedTarget_2021/' in url and '/P24iy/2021/' in url and url.endswith('.picoDst.root')):raise ValueError('wrong dataset/url row '+str(number))
    ev=int(fields[1]) if len(fields)==2 else None
    if ev is not None and ev<0:raise ValueError('negative catalog count')
    inputs.append({'url':url,'catalog_events':ev})
   if len(inputs)!=a.expected_files_per_job:reasons.append('unexpected_input_count')
   if len({x['url'] for x in inputs})!=len(inputs):reasons.append('duplicate_url_within_job')
   row.update(inputs=inputs,input_count=len(inputs),dataset_verified=True,catalog_events=sum(x['catalog_events'] for x in inputs) if all(x['catalog_events'] is not None for x in inputs) else None)
  except Exception as e:row['inputs']=[];row['input_count']=0;row['dataset_verified']=False;row['input_list_error']=str(e);reasons.append('input_list_invalid_or_missing')
  try:
   with open(row['producer_log_abs'],'rb') as f:
    f.seek(0,2);length=f.tell();f.seek(max(0,length-a.stdout_tail_bytes));tail=f.read()
   text=tail.decode(errors='replace');ms=list(marker.finditer(text));cs=list(counters.finditer(text));row['stdout_tail_sha256']=hashlib.sha256(tail).hexdigest();row['stdout_size']=length
   if not ms:reasons.append('no_final_completion_marker');row['stdout_completed_zero']=False
   else:
    m=ms[-1];row.update(log_subjob=int(m.group(1)),input_events=int(m.group(2)),completed_events=int(m.group(3)),accepted_events=int(m.group(4)),exit_status=int(m.group(5)),historical_real_seconds=float(m.group(6)),historical_cpu_seconds=float(m.group(7)),final_completion_marker=m.group(0))
    valid=int(m.group(1))==row['subjob'] and int(m.group(2))==int(m.group(3)) and int(m.group(2))>0 and int(m.group(5))==0;row['stdout_completed_zero']=valid
    if not valid:reasons.append('completion_marker_not_success')
    if not cs:reasons.append('no_producer_counter_summary')
    else:
     c=cs[-1];row.update(producer_input_events=int(c.group(1)),producer_accepted_events=int(c.group(2)),selected_lambdas=int(c.group(3)),mixing_vz_clamped=int(c.group(4)),producer_log_status=c.group(5))
     if (int(c.group(1))!=row.get('input_events') or int(c.group(2))!=row.get('accepted_events') or c.group(5)!='completed'):reasons.append('producer_log_counter_disagreement')
  except Exception as e:row['stdout_completed_zero']=False;row['stdout_error']=str(e);reasons.append('producer_stdout_missing_or_unreadable')
  row['exclusion_reasons']=reasons;row['preflight_eligible']=not reasons;row['requires_worker_root_metadata_verification']=True
  return row
 began=time.time()
 with ThreadPoolExecutor(max_workers=a.threads) as pool:rows=list(pool.map(inspect,rows))
 # A PicoDst may contribute only once. Reject the entire later job if its input list overlaps.
 owners={};duplicates=[]
 for row in rows:
  if not row['preflight_eligible']:continue
  overlap=[{'url':i['url'],'first_subjob':owners[i['url']],'later_subjob':row['subjob']} for i in row['inputs'] if i['url'] in owners]
  if overlap:
   row['preflight_eligible']=False;row['exclusion_reasons'].append('pico_url_overlap_with_earlier_eligible_job');duplicates.extend(overlap)
  else:
   for i in row['inputs']:owners[i['url']]=row['subjob']
 # Reuse the prior small, explicitly validated ROOT-metadata probe only when the file fingerprint still matches.
 prior=project/'rootfile'/a.ana_name/'pdf/track_merging_quickcheck_20261111_001/input_audit.json'
 if prior.exists():
  known=json.loads(prior.read_text())
  for sample in known['samples']:
   for row in rows:
    if row['subjob']!=sample['subjob']:continue
    fp=sample['fingerprint']
    if (row['cache_abs'],row['size'],row['mtime_ns'],row['inode'],row['device'])==(fp['path'],fp['sizeBytes'],fp['mtimeNs'],fp['inode'],fp['device']):
     row['uuid']=sample['uuid'];row['producer_root_metadata_verified']=sample['producerStatus']=='completed';row['known_root_metadata_probe']=str(prior);row['known_signal_pairs_all_kstar']={k:v['entries'] for k,v in sample['histograms'].items() if k.startswith('hKstarSE')}
 eligible=[r for r in rows if r['preflight_eligible']];needs=[r for r in rows if r['size']>0 and not r['preflight_eligible']]
 fields=['subjob','cache_abs','original_list_abs','size','mtime_ns','inode','device','uuid','producer_log_abs','input_events','accepted_events','selected_lambdas','input_count','list_sha256','producer_root_metadata_verified','requires_worker_root_metadata_verification','urls_json']
 with (out/'candidates.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
  for row in eligible:
   d=dict(row,urls_json=json.dumps([i['url'] for i in row['inputs']],separators=(',',':')));w.writerow([d.get(k,'') for k in fields])
 with (out/'excluded_candidates.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['subjob','cache_abs','size','header_status','reasons'])
  for row in rows:
   if not row['preflight_eligible']:w.writerow([row[k] for k in ['subjob','cache_abs','size','header_status']]+[';'.join(row['exclusion_reasons'])])
 summary={'root_files':len(rows),'positive_root_files':sum(r['size']>0 for r in rows),'root_bytes':sum(r['size'] for r in rows),'header_consistent_files':sum(r['header_status']=='consistent_fend' for r in rows),'preflight_eligible_jobs':len(eligible),'excluded_jobs':len(rows)-len(eligible),'positive_but_not_primary_eligible':len(needs),'unique_eligible_pico_files':len(owners),'eligible_measured_input_events_from_producer_logs':sum(r['input_events'] for r in eligible),'eligible_accepted_events_from_producer_logs':sum(r['accepted_events'] for r in eligible),'eligible_selected_lambdas_from_logs':sum(r['selected_lambdas'] for r in eligible),'eligible_candidate_root_bytes':sum(r['size'] for r in eligible),'catalog_events_eligible_inputs':sum(r['catalog_events'] or 0 for r in eligible),'duplicate_pico_overlap_jobs_rejected':len({d['later_subjob'] for d in duplicates}),'exclusion_reasons':dict(Counter(x for r in rows for x in r['exclusion_reasons'])),'root_producer_metadata_already_verified_jobs':[r['subjob'] for r in eligible if r['producer_root_metadata_verified']],'elapsed_seconds':time.time()-began}
 prior_samples=None
 archive=project/'rootfile'/a.ana_name/'tests/storage_cleanup_20261104_001/cleanup_audit_records.tar.gz'
 if archive.exists():
  with tarfile.open(archive,'r:gz') as t:
   target=[m for m in t.getmembers() if m.name=='current_job_verified_size_samples.json']
   if len(target)==1:prior_samples=json.loads(t.extractfile(target[0]).read())
  (out/'historical_verified_size_samples.json').write_text(json.dumps(prior_samples,indent=2)+'\n')
 doc={'status':'inventory_complete_primary_preflight_cohort','created_at_utc':datetime.now(timezone.utc).isoformat(),'source_jobid':a.jobid,'ana_name':a.ana_name,'dataset':'production_3p85GeV_fixedTarget_2021 / P24iy','project':str(project),'scope':'Nonrecursive preserved candidate ROOT inventory; header consistency, matching input-list dataset, exact stdout completion/counters. No candidate basket or Pico payload scan.','requirements_before_dispatch':['Worker must open original candidate ROOT READ and require nonZombie/nonRecovered, FemtoLambdaRunStatus=completed, expected candidate/ledger trees and matching measured counts.','Worker must validate original cache UUID/selection layout and track/event alignment against the listed 2021 P24iy PicoDst inputs.','Do not merge partial/failed output or count catalog events as measured coverage.'],'settings':{'threads':a.threads,'stdout_tail_bytes':a.stdout_tail_bytes,'expected_files_per_job':a.expected_files_per_job},'summary':summary,'duplicate_pico_overlaps':duplicates,'candidates_tsv':str(out/'candidates.tsv'),'candidates_tsv_sha256':hashlib.sha256((out/'candidates.tsv').read_bytes()).hexdigest(),'rows':rows}
 (out/'input_inventory.json').write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n')
 print(json.dumps(summary,indent=2),flush=True);print('INPUT_INVENTORY_PASS',flush=True)
if __name__=='__main__':main()
