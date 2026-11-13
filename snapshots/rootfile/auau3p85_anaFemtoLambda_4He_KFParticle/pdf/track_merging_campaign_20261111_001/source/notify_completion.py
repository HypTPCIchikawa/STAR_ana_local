#!/usr/bin/env python3
"""Send one user-requested email after final QA aggregation and figures validate."""
import argparse,csv,datetime,email.message,email.policy,email.utils,fcntl,getpass,hashlib,importlib.util,json,os,socket,subprocess,time
from pathlib import Path


def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,d):
 p=Path(p);tmp=p.with_name(p.name+'.partial')
 with tmp.open('w') as f:
  json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(str(tmp),str(p))
def wait(reason,status=None):
 return dict(ready=False,checkedAtUtc=now(),reason=reason,completedJobs=(status or {}).get('completedJobs'),failedJobs=(status or {}).get('failedJobs'),pendingJobs=(status or {}).get('pendingJobs'))

def queue_records(monitor):
 constraint='Owner == "'+monitor['owner']+'" && MergingCampaign == "'+monitor['campaign_tag']+'"'
 r=subprocess.run(['condor_q','-constraint',constraint,'-json'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=45)
 if r.returncode:raise RuntimeError('Cannot confirm empty queue: '+r.stderr.strip())
 return json.loads(r.stdout) if r.stdout.strip() else []

def inspect_completion(cfg,queue_checker=queue_records):
 c=Path(cfg['campaign']);s=read(c/'status.json')
 if s['phase'] not in ['complete','complete_with_failures']:return wait('campaign_not_terminal',s)
 if s['pendingJobs'] or s['completedJobs']+s['failedJobs']!=s['expectedSourceJobs']:return wait('incomplete_source_partition',s)
 if not s.get('queueReadOk') or s.get('queueByStatus') or s.get('integrityErrors') or s.get('mergeError'):return wait('campaign_checks_not_clear',s)
 monitor=read(cfg['monitor_settings']);inventory=Path(monitor['inventory'])
 if sha(inventory)!=monitor['inventory_sha256']:raise ValueError('Frozen inventory changed')
 with inventory.open() as f:rows={r['subjob']:r for r in csv.DictReader(f,delimiter='\t')}
 if len(rows)!=s['expectedSourceJobs']:raise ValueError('Expected source membership differs')
 success=set();failure=set();routes=monitor.get('source_result_routes',{})
 for sid,row in rows.items():
  spec=routes.get(sid,dict(results_dir=monitor['results_dir'],runtime_manifest_sha256=monitor['runtime_manifest_sha256'],policy_sha256=monitor['policy_sha256']))
  out=Path(spec['results_dir'])/sid;marker=out/'worker_result.json'
  if not marker.exists():return wait('missing_current_worker_marker:'+sid,s)
  r=read(marker)
  if r.get('status')=='failed':failure.add(sid);continue
  if r.get('status')!='completed':return wait('nonterminal_worker:'+sid,s)
  for name in ['summary.json','receipt.json','qa.root']:
   if not (out/name).is_file():return wait('missing_committed_output:'+sid+'/'+name,s)
  q=read(out/'summary.json');receipt=read(out/'receipt.json')
  if receipt!=r or q!=r.get('qaSummary') or not q.get('passed') or q.get('closureChecks')!=244827:raise ValueError('Committed source proof differs: '+sid)
  if r['runtimeManifestSha256']!=spec['runtime_manifest_sha256'] or r['policySha256']!=spec['policy_sha256']:raise ValueError('Source runtime route differs: '+sid)
  if Path(r['sourceCache']).resolve()!=Path(row['cache_abs']).resolve() or r['inputListSha256']!=row['list_sha256']:raise ValueError('Source identity differs: '+sid)
  if q['counters']['cache_events']!=int(row['accepted_events']) or (out/'qa.root').stat().st_size!=r['outputRootBytes']:raise ValueError('Source count/output size differs: '+sid)
  success.add(sid)
 if len(success)!=s['completedJobs'] or len(failure)!=s['failedJobs'] or success&failure or success|failure!=set(rows):raise ValueError('Current result partition differs')
 recorded_failed=read(c/'failed_jobs.json')
 if set(recorded_failed)!=failure:raise ValueError('Recorded failures differ from current markers')
 if queue_checker(monitor):return wait('queue_is_not_empty',s)
 aggregate=None;plots=None;pdf=None
 if success:
  aggregate=read(c/'aggregate/latest.json')
  if not aggregate.get('passed') or not aggregate.get('validation',{}).get('passed') or not aggregate['validation'].get('countersChecked'):raise ValueError('Final histogram validation not passed')
  if aggregate['sourceCount']!=len(success) or len(aggregate['sourceJobs'])!=len(success) or set(aggregate['sourceJobs'])!=success: return wait('aggregate_is_not_final',s)
  if aggregate['counters']!=s['counters']:raise ValueError('Final counters differ')
  if sha(aggregate['root'])!=aggregate['rootSha256']:raise ValueError('Final aggregate ROOT hash differs')
  if not (c/'plots/latest.json').exists():return wait('final_figures_pending',s)
  plots=read(c/'plots/latest.json')
  if plots.get('status')!='passed' or not plots.get('inputUnchanged') or plots['sourceRootSha256']!=aggregate['rootSha256'] or plots['sourceCount']!=len(success) or plots['sourceJobs']!=aggregate['sourceJobs']:return wait('final_figures_pending',s)
  if 'merging_maps.pdf' not in plots['artifacts']:raise ValueError('Final PDF has no hash proof')
  folder=Path(plots['manifest']).parent
  for name,artifact in plots['artifacts'].items():
   path=folder/name
   if not path.is_file():return wait('missing_figure:'+name,s)
   if path.stat().st_size!=artifact['bytes'] or sha(path)!=artifact['sha256']:raise ValueError('Final figure hash differs: '+name)
  pdf=str(folder/'merging_maps.pdf')
  if not Path(pdf).is_file():return wait('missing_final_pdf',s)
  if read(c/'aggregate/latest.json')!=aggregate or read(c/'plots/latest.json')!=plots:return wait('final_snapshot_advanced',s)
 if read(c/'status.json')!=s:return wait('campaign_status_advanced',s)
 if queue_checker(monitor):return wait('queue_became_nonempty',s)
 return dict(ready=True,checkedAtUtc=now(),completedJobs=len(success),failedJobs=len(failure),expectedSourceJobs=len(rows),failedSourceJobs=sorted(failure,key=int),status=s,aggregate=aggregate,plots=plots,pdf=pdf,verification='current route markers/receipt/summary, full source partition, empty queue, final ROOT SHA/counters and matching PDF/PNG hashes')

def make_message(cfg,proof):
 s=proof['status'];counts=s.get('counters',{});label='終了・除外あり' if proof['failedJobs'] else '完了'
 msg=email.message.EmailMessage(policy=email.policy.SMTP)
 msg['From']=getpass.getuser()+'@'+socket.getfqdn();msg['To']=cfg['recipient'];msg['Date']=email.utils.format_datetime(datetime.datetime.now().astimezone())
 msg['Subject']=cfg['subject_prefix']+' '+label
 msg['Message-ID']='<'+Path(cfg['campaign']).name+'-completion@'+socket.getfqdn()+'>'
 jst=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).strftime('%Y-%m-%d %H:%M:%S JST')
 lines=['3.85 GeV Λ–⁴He track-merging QAのjobがすべて終了しました。', '成功分のhadd・検証・最終図更新まで確認しました。' if proof['completedJobs'] else '成功したQA出力はありません。失敗記録を確認してください。', '', '確認時刻: '+jst, '成功: '+str(proof['completedJobs'])+' / 対象: '+str(proof['expectedSourceJobs']), '失敗・除外: '+str(proof['failedJobs']), '成功分の元入力イベント数: '+format(s.get('sourceInputEventsRepresented',0),','), '⁴He–Λ信号SEペア: '+format(counts.get('signal_pairs_selected_4He',0),','), '³He–Λ信号SEペア: '+format(counts.get('signal_pairs_selected_3He',0),','), '', '最終確認図PDF: '+str(proof['pdf'] or 'なし'), '結果・除外一覧: '+str(Path(cfg['campaign'])/'status.json'), '監視記録: '+str(cfg['state_dir']), '', 'この通知はjob・成功分の統合と図の完成を知らせるものです。原点ピークや効率の物理的判断は別途確認します。']
 msg.set_content('\n'.join(lines)+'\n');return msg

def load_sender(cfg):
 spec=importlib.util.spec_from_file_location('star_existing_completion_sender',cfg['sender_module']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--settings',type=Path,required=True);mode=ap.add_mutually_exclusive_group(required=True);mode.add_argument('--check',action='store_true');mode.add_argument('--watch',action='store_true');a=ap.parse_args();cfg=read(a.settings)
 if '\n' in cfg['recipient'] or email.utils.parseaddr(cfg['recipient'])[1]!=cfg['recipient']:raise ValueError('Invalid recipient')
 if cfg['poll_seconds']<30:raise ValueError('Polling interval too small')
 state=Path(cfg['state_dir']);state.mkdir(parents=True,exist_ok=True)
 if a.check:print(json.dumps(inspect_completion(cfg),ensure_ascii=False));return 0
 with (state/'watch.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (state/'sent.json').exists():print('MAIL_ALREADY_ACCEPTED');return 0
  if (state/'delivery_attempt.json').exists():raise RuntimeError('Previous delivery attempted; verify delivery before retry')
  while True:
   try:
    proof=inspect_completion(cfg)
   except Exception as e:proof=dict(ready=False,checkedAtUtc=now(),reason='verification_error',error=str(e))
   atomic(state/'watch_status.json',proof)
   if proof['ready']:
    message=make_message(cfg,proof);(state/'completion_message.eml').write_bytes(message.as_bytes());atomic(state/'completion_verification.json',proof)
    sender=load_sender(cfg)
    try:sender.send_once(state,message,proof,cfg['recipient'])
    except SystemExit as e:
     atomic(state/'delivery_status.json',dict(status='outcome_uncertain_no_automatic_retry',timeUtc=now(),reason=str(e)));raise
    atomic(state/'delivery_status.json',dict(status='accepted_by_local_sendmail',timeUtc=now(),recipient=cfg['recipient'],messageId=str(message['Message-ID'])));print('COMPLETION_MAIL_ACCEPTED',flush=True);return 0
   print(now(),'waiting',proof.get('reason'),'completed',proof.get('completedJobs'),'remaining',proof.get('pendingJobs'),flush=True)
   time.sleep(cfg['poll_seconds'])

if __name__=='__main__':raise SystemExit(main())
