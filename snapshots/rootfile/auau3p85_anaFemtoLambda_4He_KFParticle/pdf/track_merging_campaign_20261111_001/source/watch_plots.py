#!/usr/bin/env python3
"""Refresh verified raw-SE figures when the aggregate changes."""
import argparse,fcntl,json,os,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path

def read(p):
 try:return json.loads(p.read_text())
 except FileNotFoundError:return None

def atomic(p,d):
 t=p.with_name(p.name+'.partial');t.write_text(json.dumps(d,indent=2)+'\n');os.replace(str(t),str(p))


def terminal_snapshot(c):
 """Exit only when fresh final status, aggregate, and plots agree."""
 status=read(c/'status.json') or {}
 if status.get('phase') not in ['complete','complete_with_failures']:
  return False,{}
 # Read after terminal status: the aggregate used before drawing may be old.
 latest=read(c/'aggregate/latest.json')
 plot=read(c/'plots/latest.json')
 def wait(reason):return False,{'finalReadinessError':reason}
 completed=status.get('completedJobs');failed=status.get('failedJobs');expected=status.get('expectedSourceJobs')
 if not all(isinstance(v,int) and v>=0 for v in [completed,failed,expected]):
  return wait('Invalid final job counts')
 if status.get('pendingJobs')!=0 or completed+failed!=expected:
  return wait('Final job outcomes are incomplete')
 if status.get('queueReadOk') is not True or status.get('queueByStatus') or status.get('integrityErrors') or status.get('mergeError'):
  return wait('Final campaign has queued jobs or unresolved errors')
 if completed==0:
  return True,dict(phase='complete',sourceCount=0,sourceRootSha256=None,plotManifest=None,noSuccessfulSources=True)
 if not latest or latest.get('passed') is not True or latest.get('validation',{}).get('passed') is not True:
  return wait('Final verified aggregate is missing')
 jobs=latest.get('sourceJobs',[])
 if latest.get('sourceCount')!=completed or len(jobs)!=completed or len(set(jobs))!=completed:
  return wait('Final aggregate membership/count differs from completed jobs')
 if not plot or plot.get('status')!='passed' or plot.get('inputUnchanged') is not True:
  return wait('Final verified plots are missing')
 if not latest.get('rootSha256') or plot.get('sourceRootSha256')!=latest['rootSha256'] or plot.get('sourceCount')!=completed or plot.get('sourceJobs')!=jobs:
  return wait('Plots do not represent the final aggregate')
 return True,dict(phase='complete',sourceCount=completed,sourceRootSha256=latest['rootSha256'],plotManifest=plot.get('manifest'),noSuccessfulSources=False)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--campaign',type=Path,required=True);ap.add_argument('--poll-seconds',type=int,required=True);ap.add_argument('--once',action='store_true');a=ap.parse_args();c=a.campaign.resolve(strict=True)
 with (c/'plot_monitor.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  while True:
   state={'updatedAtUtc':datetime.now(timezone.utc).isoformat(),'phase':'waiting_for_verified_aggregate'}
   latest=read(c/'aggregate/latest.json')
   if latest:
    try:
     r=subprocess.run([sys.executable,str(Path(__file__).resolve().with_name('plot_campaign.py')),'--campaign',str(c),'--output-dir',str(c/'plots'),'--once'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=300)
     if r.returncode:raise RuntimeError(r.stdout[-4000:])
     plot=read(c/'plots/latest.json')
     state.update(phase='updated' if plot and plot.get('sourceRootSha256')==latest['rootSha256'] else 'aggregate_advanced_retry_next_poll',sourceCount=latest['sourceCount'],sourceRootSha256=latest['rootSha256'],plotManifest=plot.get('manifest') if plot else None)
     if r.stdout.strip():print(r.stdout.strip(),flush=True)
    except Exception as e:state.update(phase='needs_review_plot',error=str(e));print(state,flush=True)
   terminal,final_state=terminal_snapshot(c)
   if terminal:state.update(final_state)
   elif final_state:state.update(phase='waiting_for_final_plot',**final_state)
   atomic(c/'plot_status.json',state)
   if a.once or terminal:return
   time.sleep(a.poll_seconds)

if __name__=='__main__':main()
