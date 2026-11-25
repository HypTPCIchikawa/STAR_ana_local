#!/usr/bin/env python3
"""Verify one 20-period replay, produce four/eight Y CFs, notify once."""
import argparse,fcntl,hashlib,json,os,re,subprocess,sys,time,types
from pathlib import Path
import xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]

def load(path):
    with Path(path).open() as f:return json.load(f)
def need(ok,message):
    if not ok:raise RuntimeError(message)
def resolve(value,base):
    p=Path(value);return (p if p.is_absolute() else base/p).resolve()
def module(path,name):
    m=types.ModuleType(name);m.__file__=str(path);exec(compile(path.read_text(),str(path),'exec'),m.__dict__);return m

def prepare(args,launched):
    sp=args.settings.resolve();s=load(sp);need(s['schemaVersion']==1 and s['pollSeconds']>=30 and s['settleSeconds']>=0,'invalid explicit watcher settings')
    need(s['commandTimeoutSeconds']>0 and s['postprocessTimeoutSeconds']>0 and re.fullmatch(r'[^\s<>@]+@[^\s<>@]+',s['recipient']),'invalid timeouts/recipient')
    for p in s['fatalPatterns']:re.compile(p)
    helper_path=resolve(s['sharedCompletionHelper'],sp.parent);m=module(helper_path,'he4_y_completion_helpers')
    campaign_path=resolve(s['analysisCampaign'],sp.parent);campaign=load(campaign_path);need(campaign['schemaVersion']==1 and campaign['expectedPeriods']==20,'campaign schema/coverage')
    source_manifest=Path(campaign['sourceManifest']);source_list=Path(campaign['inputList']);need(m.sha(source_manifest)==campaign['sourceManifestSha256'] and m.sha(source_list)==campaign['inputListSha256'],'source manifest/list pin mismatch')
    sources=load(source_manifest)['sourceFiles'];need(len(sources)==20 and sum(x['acceptedEvents'] for x in sources)==campaign['storedEvents'],'source event coverage mismatch')
    need([str(Path(x[0]).resolve()) for x in m.list_rows(source_list)]==[str(Path(x['path']).resolve()) for x in sources],'source order mismatch')
    source_stats=[dict(path=x['path'],**{k:x[k] for k in ['bytes','device','inode','mtimeNs']}) for x in sources]
    for st in source_stats:m.same_stat(st)
    post_path=resolve(s['postprocessSettings'],sp.parent);post=module(HERE/'postprocess_he4_y.py','he4_y_postprocess_validation');post_settings=post.validate_settings(post_path)
    need(post_settings['mainconfForRoot']==campaign['mainconf'] and post_settings['rapidityBinning']=={k:campaign['rapidity'][k] for k in ['frame','nBins','min','max']},'postprocess campaign binning/config mismatch')
    cp=args.launch_config.resolve();config=load(cp);need(config['anaName']==campaign['anaName'] and config['mainconf']==campaign['mainconf'] and config['outputFileStem']==campaign['outputFileStem'] and Path(config['outputDirectory']).resolve()==Path(campaign['outputDirectory']).resolve(),'launch config identity mismatch')
    need(config['inputListSha256']==campaign['inputListSha256'] and m.sha(config['inputList'])==campaign['inputListSha256'],'launch source list differs')
    need([(x['path'],x['acceptedEvents']) for x in config['sources']]==[(x['path'],x['acceptedEvents']) for x in sources],'launch source records differ')
    paths=[Path(__file__).resolve(),sp,campaign_path,source_manifest,source_list,cp,post_path,HERE/'postprocess_he4_y.py',HERE/'ExtractHe4Y.C',HERE/'DrawHe4Y.C',helper_path,resolve(s['sharedEventParser'],sp.parent),resolve(post_settings['sharedCFMath'],post_path.parent),resolve(post_settings['sharedExtractorUtilities'],post_path.parent),resolve(campaign['mainconf'],PROJECT)]
    paths += [Path(x) for x in campaign['analysisFiles'].values()]
    pins=[m.pin(p) for p in dict.fromkeys(paths)]
    if not launched:return s,campaign,sources,m,pins,source_stats,None,None,post_path
    state=load(args.launch_state.resolve());need(state['phase'] in ['submitted_verified','monitor_starting','monitor_started'],'launch is not verified')
    need(state['actualJobs']==20 and state.get('submitReturnCode',0)==0 and re.fullmatch(r'[0-9A-Fa-f]{32}',state['jobid']),'invalid submitted job identity/count')
    meta_path=resolve(state['runmeta'],PROJECT);meta=load(meta_path);need(meta['jobid']==state['jobid'] and Path(meta['projectRoot']).resolve()==PROJECT,'runmeta identity differs')
    snapshot=meta['submit']['joblistSnapshot'];xp=Path(snapshot['path']);need(m.sha(xp)==snapshot['sha256'] and xp.stat().st_size==snapshot['sizeBytes'],'submitted XML pin differs')
    xml=ET.parse(str(xp)).getroot();stdout=m.uri_path(xml.find('stdout').get('URL'));stderr=m.uri_path(xml.find('stderr').get('URL'));outputs=[x for x in xml.findall('output') if x.get('fromScratch').endswith('.root')];need(len(outputs)==1,'one ROOT output required');out=outputs[0];od=Path(config['outputDirectory']).resolve();need(m.uri_path(out.get('toURL')).resolve()==od and Path(out.get('fromScratch')).name==config['outputFileStem']+'_$JOBID.root','submitted ROOT directory/stem differs')
    need('$JOBID' in str(stdout) and '$JOBID' in str(stderr),'unrecognized native log templates')
    native=m.native_inputs(meta,sources);run=dict(caseId='he4Y8',jobid=state['jobid'],runmeta=str(meta_path),eventLog=str(resolve(state['eventLog'],PROJECT)),submitHost=config['submitHost'],nativeInputs=native,roots=[str(od/(config['outputFileStem']+'_'+state['jobid']+'_'+str(i)+'.root')) for i in range(20)],stdoutTemplate=str(stdout),stderrTemplate=str(stderr))
    pins += [m.pin(meta_path),m.pin(xp)]
    parser=m.event_parser(resolve(s['sharedEventParser'],sp.parent))
    return s,campaign,sources,m,pins,source_stats,run,parser,post_path

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--launch-config',type=Path,required=True);ap.add_argument('--launch-state',type=Path,required=True);ap.add_argument('--settings',type=Path,required=True);ap.add_argument('--state-dir',type=Path);mode=ap.add_mutually_exclusive_group();mode.add_argument('--validate-settings',action='store_true');mode.add_argument('--check',action='store_true');a=ap.parse_args()
    if a.validate_settings:
        s,c,sources,m,pins,stats,_,_,_=prepare(a,False);m.check_pins(pins,stats);print(json.dumps(dict(passed=True,newWorkers=20,launchStatesRead=False,schedulerQueried=False,postprocessRun=False,mailSent=False,pins=pins,sourceStats=stats),indent=2));return
    s,c,sources,m,pins,stats,run,parser,post=prepare(a,True)
    if a.check:
        result=m.poll([run],parser,s['commandTimeoutSeconds']);result.update(checkOnly=True,mailSent=False,postprocessRun=False)
        if result['allTerminal']:result['completion']=m.verify_outputs(run,sources,[re.compile(x) for x in s['fatalPatterns']])
        m.check_pins(pins,stats);print(json.dumps(result,indent=2));return
    need(a.state_dir is not None,'normal watch requires state-dir');state=a.state_dir.resolve();need(m.inside(state,HERE) and state!=HERE,'state must be under this study');state.mkdir(exist_ok=True);lock=(state/'watcher.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    need(not any((state/x).exists() for x in ['delivery_attempt.json','sent.json','attention.json']),'prior attention/mail attempt requires manual review')
    try:
        startup=dict(pins=pins,sourceStats=stats,jobid=run['jobid'],runmeta=run['runmeta'],eventLog=run['eventLog'],automaticRetry=False,automaticResubmit=False,automaticHadd=False)
        if (state/'startup_pins.json').exists():
            old=load(state/'startup_pins.json');need(old['pins']==pins and old['sourceStats']==stats and old['jobid']==run['jobid'],'restart pins differ')
        else:m.exclusive_json(state/'startup_pins.json',startup)
        import socket
        m.atomic_json(state/'pid.json',dict(pid=os.getpid(),host=socket.getfqdn(),startedUTC=m.now()))
        first=None
        while True:
            m.check_pins(pins,stats);status=m.poll([run],parser,s['commandTimeoutSeconds']);status.update(phase='waiting',expectedWorkers=20)
            if status['allTerminal']:
                if first is None:first=time.monotonic()
                status['phase']='all_terminal_settling';status['settleElapsedSeconds']=time.monotonic()-first
                if status['settleElapsedSeconds']>=s['settleSeconds']:break
            else:first=None
            m.atomic_json(state/'status.json',status);print(json.dumps(status),flush=True);time.sleep(s['pollSeconds'])
        completion=m.verify_outputs(run,sources,[re.compile(x) for x in s['fatalPatterns']]);need(completion['completedEvents']==c['storedEvents'],'completed event total differs');m.check_pins(pins,stats);m.exclusive_json(state/'completion_verified.json',completion)
        output_list=state/'output_roots.list'
        with output_list.open('x') as f:f.write('\n'.join(run['roots'])+'\n')
        result=state/'he4_y_results';command=[sys.executable,str(HERE/'postprocess_he4_y.py'),'--settings',str(post),'--outputs',str(output_list),'--output-dir',str(result)]
        m.exclusive_json(state/'postprocess_attempt.json',dict(command=command,startedUTC=m.now(),automaticRetry=False));m.atomic_json(state/'status.json',dict(phase='postprocessing',allTerminal=True))
        with (state/'postprocess.log').open('x') as f:subprocess.run(command,cwd=str(PROJECT),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=s['postprocessTimeoutSeconds'])
        done=load(result/'completed.json');need(done['passed'] and not done['pilot'] and done['sources']==20 and done['events']==c['storedEvents'],'full postprocess incomplete');m.check_pins(pins,stats)
        for st in completion['outputStats']:m.same_stat(st)
        m.exclusive_json(state/'postprocess_verified.json',dict(passed=True,completed=m.pin(result/'completed.json'),metrics=m.pin(result/'he4_y_cf_metrics.json')));m.atomic_json(state/'status.json',dict(phase='completed_verified',workers=20,events=c['storedEvents'],comparisonDirectory=str(result),physicsSignificanceAssessed=False));m.send_once(state,s,True);print('HE4_Y_COMPLETION_PASS mail=accepted_by_local_sendmail',flush=True)
    except Exception as ex:
        m.atomic_json(state/'attention.json',dict(error=str(ex),atUTC=m.now(),automaticRetry=False,automaticResubmit=False));m.atomic_json(state/'status.json',dict(phase='attention',error=str(ex)))
        if not (state/'delivery_attempt.json').exists() and not (state/'sent.json').exists():
            try:m.send_once(state,s,False)
            except Exception as delivery:m.atomic_json(state/'notification_error.json',dict(error=str(delivery),automaticResend=False))
        raise
if __name__=='__main__':
    try:main()
    except Exception as ex:print('HE4_Y_WATCH_FAILED: '+str(ex),file=sys.stderr);sys.exit(1)
