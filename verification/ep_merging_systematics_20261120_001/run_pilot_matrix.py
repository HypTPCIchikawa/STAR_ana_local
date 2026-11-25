#!/usr/bin/env python3
"""Run the approved fixed-prefix replay controls, at most four at a time."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
STUDY=Path(__file__).resolve().parent
EVENTS=400000
CONCURRENCY=4

def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fresh(path,value):
    with Path(path).open('x') as output:
        json.dump(value,output,indent=2);output.write('\n')

def save(path,value):
    temporary=Path(str(path)+'.next')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    os.replace(str(temporary),str(path))

def stat_record(path):
    path=Path(path);s=path.stat()
    return dict(path=str(path),bytes=s.st_size,inode=s.st_ino,device=s.st_dev,mtimeNs=s.st_mtime_ns)

def timing(path):
    result={}
    for line in Path(path).read_text().splitlines():
        for label,key in [('User time (seconds)','userSeconds'),('System time (seconds)','systemSeconds'),
                          ('Maximum resident set size (kbytes)','maxRssKiB'),('Exit status','timeExitStatus')]:
            if line.strip().startswith(label+':'):result[key]=float(line.split(':',1)[1].strip())
    return result

def main():
    matrix=json.loads((STUDY/'cases.json').read_text());assert len(matrix['cases'])==8
    source_info=json.loads(Path(matrix['inputManifest']).read_text())['sourceFiles'][0]
    source=Path(source_info['path']);source_stat=stat_record(source)
    assert all(source_stat[k]==source_info[k] for k in ['bytes','inode','device','mtimeNs'])
    assert 357247<EVENTS<source_info['acceptedEvents']
    reference=ROOT/'rootfile/femto_lambda_period_mix100_20261119_001/pilot_B_period100.root'
    directory=STUDY/'pilots';directory.mkdir(exist_ok=False)
    paths={Path(row['path']) for row in matrix['generatedConfigFiles']}
    before=json.loads((STUDY/'source_before_manifest.json').read_text())['files']
    for row in before:
        path=Path(row['path']);paths.add(path if path.is_absolute() else ROOT/path)
    paths.update((ROOT/'include').glob('FemtoLambda*.h'))
    paths.update(ROOT/'lib'/name for name in ['libStarAnaConfig.so','libStRefMultCorr.so','libKFParticle.so',
                  'libStKfParticleCommon.so','libStCommon.so','libStFemtoMaker.so'])
    paths.update([ROOT/'include/FemtoEpdEventPlane.h',ROOT/'include/FemtoCandidate.h',Path(__file__).resolve(),
                  ROOT/'script/singularity_run_anaFemtoLambda_4He_tree.sh',ROOT/'script/run_femtoLambda_tree_common.sh'])
    pins={str(path):sha(path) for path in sorted(paths)}
    for row in matrix['generatedConfigFiles']:assert pins[row['path']]==row['sha256']
    reference_sha=sha(reference)
    fresh(STUDY/'pilot_matrix_pins.json',dict(createdJst=now(),source=source_stat,
          referencePilot=dict(**stat_record(reference),sha256=reference_sha),sourceConfigLibSha256=pins,
          casesManifestSha256=sha(STUDY/'cases.json'),buildLogSha256=sha(STUDY/'pilot_incremental_build.log')))
    state=dict(startedJst=now(),phase='running',storedEventsPerCase=EVENTS,maximumConcurrentProcesses=CONCURRENCY,
               source=str(source),cases={},rssSamples=[])
    pending=list(matrix['cases']);processes={}
    while pending or any('returnCode' not in state['cases'][name] for name in processes):
        live=sum(p.poll() is None for p in processes.values())
        while pending and live<CONCURRENCY:
            case=pending.pop(0);name=case['caseId'];output=directory/(name+'.root')
            command=['/usr/bin/time','-v','-o',str(directory/(name+'.time')),'bash',
                'script/singularity_run_anaFemtoLambda_4He_tree.sh',case['mainconf'],str(source),str(output),
                'epMergingPilot_'+name,str(EVENTS)]
            with (directory/(name+'.log')).open('x') as log:
                process=subprocess.Popen(command,cwd=str(ROOT),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            processes[name]=process
            state['cases'][name]=dict(pid=process.pid,command=command,startedJst=now(),output=str(output),
                mainconf=case['mainconf'],epBins=case['epBins'],cutEnabled=case['cutEnabled'])
            print(now(),'START',name,flush=True);live+=1
        readings=subprocess.check_output(['ps','-eo','pgid=,rss='],universal_newlines=True)
        rss={name:0 for name in processes}
        for line in readings.splitlines():
            group,usage=map(int,line.split())
            for name,process in processes.items():
                if group==process.pid:rss[name]+=usage
        state['rssSamples'].append(dict(checkedJst=now(),processGroupRssKiB=rss))
        for name,process in processes.items():
            if process.poll() is None or 'returnCode' in state['cases'][name]:continue
            row=state['cases'][name];row.update(returnCode=process.returncode,finishedJst=now())
            text=(directory/(name+'.log')).read_text(errors='replace')
            markers=[line for line in text.splitlines() if line.startswith('FemtoLambda tree replay job=')]
            row['completionMarkers']=markers
            match=re.search(r'requested=(\d+) completed=(\d+) available=(\d+) status=(\d+) RealTime=([0-9.eE+-]+) CpuTime=([0-9.eE+-]+)',markers[0]) if len(markers)==1 else None
            row['passed']=(process.returncode==0 and match is not None and int(match.group(1))==int(match.group(2))==EVENTS
                and int(match.group(3))==source_info['acceptedEvents'] and int(match.group(4))==0)
            if row['passed']:row.update(eventLoopWallSeconds=float(match.group(5)),eventLoopCpuSeconds=float(match.group(6)))
            row.update(timing(directory/(name+'.time')))
            output=Path(row['output']);row['outputStat']=stat_record(output) if output.exists() else None
            fresh(directory/(name+'.result.json'),row)
            print(now(),'FINISHED',name,'passed='+str(row['passed']),flush=True)
        save(STUDY/'pilot_matrix_status.json',state)
        if pending or any(p.poll() is None for p in processes.values()):time.sleep(10)
    assert stat_record(source)==source_stat and sha(reference)==reference_sha
    assert all(sha(path)==expected for path,expected in pins.items()),'Code/config/lib changed during pilot'
    state.update(finishedJst=now(),phase='complete',sourceUnchanged=True,sourceConfigLibPinsUnchanged=True,
                 passed=len(state['cases'])==8 and all(row['passed'] for row in state['cases'].values()))
    save(STUDY/'pilot_matrix_status.json',state)
    fresh(STUDY/'pilot_matrix_results.json',{k:v for k,v in state.items() if k!='rssSamples'})
    fresh(directory/'resource.json',dict(passed=state['passed'],maximumConcurrentProcesses=CONCURRENCY,
          cases={name:{key:row.get(key) for key in ['eventLoopWallSeconds','eventLoopCpuSeconds','userSeconds',
            'systemSeconds','maxRssKiB','outputStat']} for name,row in state['cases'].items()},rssSamples=state['rssSamples']))
    return 0 if state['passed'] else 1

if __name__=='__main__':sys.exit(main())
