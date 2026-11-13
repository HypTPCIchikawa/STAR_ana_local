#!/usr/bin/env python3
"""Run each proven local input separately; candidate acceptance requires later QA."""
import datetime,hashlib,json,os,subprocess,time
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic(p,d):
 t=p.with_name(p.name+'.partial');t.write_text(json.dumps(d,indent=2)+'\n');os.replace(str(t),str(p))

c=Path(__file__).resolve().parent
project=c.parents[4]
assert (project/'StMaker').is_dir(),str(project)
provenance=c/'input_provenance.json';cfg=json.loads(provenance.read_text());inputs=cfg['inputs'];assert len(inputs)==3
mainconf=project/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml'
ref=project/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/built_libraries.json'
libs=json.loads(ref.read_text());assert all(sha(Path(p))==v['sha256'] for p,v in libs.items())
record=[]
for i,x in enumerate(inputs,1):
 p=Path(x['path']).resolve(strict=True);url=x['originalUrl'];assert url.rsplit('/',1)[-1]==p.name and '/production_3p85GeV_fixedTarget_2021/' in url and '/P24iy/2021/' in url
 st=p.stat();record.append(dict(id=i,path=str(p),originalUrl=url,size=st.st_size,mtimeNs=st.st_mtime_ns,inode=st.st_ino,sha256=sha(p)))
assert len(set(x['path'] for x in record))==3
with (c/'producer_started.intent').open('x') as f:f.write(now()+'\n')
atomic(c/'producer_inputs_frozen.json',dict(timeUtc=now(),inputProvenanceSha256=sha(provenance),inputs=record,mainconf=str(mainconf),mainconfSha256=sha(mainconf),libraries=libs))
results=[]
for x in record:
 sid=str(x['id']);out=c/'producer'/sid;out.mkdir(parents=True,exist_ok=False);listing=out/'pico.list';listing.write_text(x['path']+'\n');target=out/'candidates.root'
 cmd=['nice','-n','10','timeout','--signal=INT','--kill-after=30s','1200s','bash',str(project/'script/singularity_run_anaFemtoLambda_4He.sh'),str(mainconf),str(listing),str(target),'merging_local_'+sid+'_20261111','-1']
 atomic(c/'producer_status.json',dict(updatedAtUtc=now(),phase='running_individual_producer',currentInput=x,completedProcessCount=len(results),results=results))
 start=time.time()
 with (out/'producer.log').open('w') as log:r=subprocess.run(cmd,cwd=str(project),stdout=log,stderr=subprocess.STDOUT)
 original=Path(x['path']);st=original.stat();unchanged=(st.st_size,st.st_mtime_ns,st.st_ino)==(x['size'],x['mtimeNs'],x['inode']) and sha(original)==x['sha256']
 same_libs=all(sha(Path(p))==v['sha256'] for p,v in libs.items())
 result=dict(id=x['id'],finishedAtUtc=now(),phase='producer_returned_zero_pending_strict_QA' if r.returncode==0 and unchanged and same_libs else 'producer_failed_pending_strict_QA',returnCode=r.returncode,wallSeconds=time.time()-start,input=x,inputUnchanged=unchanged,producerLibrariesUnchanged=same_libs,command=cmd,candidateRoot=str(target),picoList=str(listing),mainconf=str(mainconf))
 atomic(out/'result.json',result);results.append(result)
 if not unchanged or not same_libs:raise RuntimeError('Original input or producer libraries changed')
atomic(c/'producer_status.json',dict(updatedAtUtc=now(),phase='producer_stage_finished_pending_strict_QA',completedProcessCount=len(results),results=results))
