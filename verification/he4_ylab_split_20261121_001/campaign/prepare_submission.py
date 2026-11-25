#!/usr/bin/env python3
"""Freeze a gated single-He4-rapidity campaign; this does not submit jobs."""
import hashlib,json,socket,subprocess,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import yaml
D=Path(__file__).resolve().parent;ROOT=D.parents[1];P=ROOT/'rootfile/femto_lambda_period_mix100_20261119_001'
LABEL='he4YLab820261121'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def once(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b,1)
def save(p,x):
    with p.open('x') as out:json.dump(x,out,indent=2);out.write('\n')
def main():
    assert not (D/'launch_config.json').exists(),'Inspect existing campaign before repeating'
    for n in ['pilot_run_result.json','pilot_comparison.json','fixture_result.json','source_preservation.json']:
        assert json.loads((D/n).read_text())['passed'] is True,n
    post=json.loads((D/'pilot_postprocess_003/completed.json').read_text())
    assert post['passed'] is True and post['pilot'] is True and post['events']==400000 and post['curves']==320
    assert post['inputStatsUnchanged'] is True and post['pinsUnchanged'] is True and post['treeEntriesRead']==0
    c=json.loads((D/'analysis_campaign.json').read_text());conf=c['mainconf'];ana=c['anaName'];output=Path(c['outputDirectory']);output.mkdir(parents=True,exist_ok=True);assert not list(output.glob('*.root'))
    text=(P/'launch_once.py').read_text()
    text=once(text,'STUDY.name == "femto_lambda_period_mix100_20261119_001"','STUDY.name == "lambda_he4_ylab_split_20261121_001"')
    text=once(text,'config["anaName"].endswith("_tree_periodMix100")','config["anaName"] == "'+ana+'"')
    text=once(text,'config["label"] == "periodMix10020261119"','config["label"] == "'+LABEL+'"')
    text=once(text,'float(xml.get("filesPerHour", "0")) == 0.025','float(xml.get("filesPerHour", "0")) == 1.0 / 24.0')
    text=once(text,'reviewed 40-hour-per-period scheduler estimate','reviewed 24-hour-per-period scheduler estimate')
    begin=text.index('        monitor = start_monitor(config, meta, runmeta, event, state)');end=text.index('    except Exception as error:',begin)
    text=text[:begin]+'''        command=[config["python"],"-u",str(STUDY/"watch_he4_y.py"),"--launch-config",str(args.config.resolve()),
                 "--launch-state",str(STUDY/"launch_state.json"),"--settings",str(STUDY/"he4_y_watch_settings.json"),
                 "--state-dir",str(STUDY/"completion_watch_001")]
        with (STUDY/"completion_watch.log").open("x") as output:
            monitor=subprocess.Popen(command,cwd=str(ROOT),stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
        save("watch_started.json",dict(pid=monitor.pid,command=command,startedUTC=now(),jobid=jobid))
        print(json.dumps(dict(phase="submitted_verified",jobid=jobid,actualJobs=20,watchPid=monitor.pid,automaticResubmit=False)),flush=True)
        code=monitor.wait()
        require(code==0 and (STUDY/"completion_watch_001/sent.json").is_file(),"Watcher stopped before confirmed completion/mail; no automatic restart")
        sent=load(STUDY/"completion_watch_001/sent.json")
        require(sent.get("mailStatus")=="accepted_by_local_sendmail" and sent.get("kind")=="success","Watcher did not record successful comparison mail")
        state.update(phase="notification_completed",completedAtUTC=now(),watcherReturnCode=code)
        save("launch_state.json",state)
        return 0
'''+text[end:]
    text=once(text,'The existing completion-only notifier sends one minimal email;\nthis launcher does not implement mail delivery or edit shared job tooling.',
        'The dedicated He4-rapidity watcher verifies outputs and creates CF plots\nbefore sending one minimal completion email. Shared tooling is not edited.')
    launcher=D/'launch_once.py';launcher.write_text(text)
    subprocess.run([sys.executable,str(ROOT/'script/generate_joblist_from_list.py'),'--mainconf',conf,'--input',str(D/'periods.list'),'--files-per-subjob','1','--label',LABEL,'--output-dir',str(output)],cwd=str(ROOT),check=True)
    xml_path=ROOT/'job/joblist'/('joblist_'+ana+'_'+LABEL+'.xml');xml=ET.parse(xml_path);job=xml.getroot();assert job.get('maxFilesPerProcess')=='1' and job.find('ResourceUsage') is None
    job.set('minFilesPerProcess','1');job.set('filesPerHour',str(1.0/24.0));resource=ET.fromstring('<ResourceUsage><Memory><MinMemory>4096</MinMemory><MaxMemory>4096</MaxMemory></Memory><StorageSpace><MinStorage>65536</MinStorage><MaxStorage>65536</MaxStorage></StorageSpace></ResourceUsage>');job.insert(list(job).index(job.find('command'))+1,resource);xml.write(str(xml_path),encoding='utf-8',xml_declaration=True)
    base=json.loads((P/'launch_config.json').read_text());fixed=[Path(p) for p in base['fingerprints'] if not p.startswith(str(P)) and '/job/joblist/' not in p and '/config/' not in p]
    fixed += [Path(__file__),launcher,xml_path,ROOT/conf,D/'analysis_campaign.json',D/'rapidity_binning.json',D/'fixture_result.json',D/'source_preservation.json']
    fixed += [ROOT/'config'/v for v in yaml.safe_load((ROOT/conf).read_text()).values() if isinstance(v,str) and v.endswith('.yaml')]
    fixed += [Path(p) for p in base['fingerprints'] if '/config/event_plane/xialei_epd_2021P24iy/' in p or '/config/calibration/' in p]
    fixed += [ROOT/'include/FemtoLambdaPairGeometry.h', ROOT/'include/FemtoEpdEventPlane.h', ROOT/'include/FemtoLambdaProvider.h', ROOT/'Makefile']
    previous=ROOT/'rootfile/femto_lambda_ep_merging_systematics_20261120_001'
    fixed += [previous/n for n in ['watch_production_matrix.py','postprocess_cf.py','ExtractPeriodCF.C']]
    fixed += [D/n for n in ['pilot_run_result.json','pilot_postprocess_003/completed.json']]
    fixed += [D/n for n in ['postprocess_he4_y.py','ExtractHe4Y.C','DrawHe4Y.C','he4_y_postprocess_settings.json','watch_he4_y.py','he4_y_watch_settings.json']]
    launch=dict(base);launch.update(directory=str(D),python=sys.executable,submitHost=socket.getfqdn(),anaName=ana,label=LABEL,submitName=ana+'_'+LABEL,outputFileStem=c['outputFileStem'],mainconf=conf,joblist=str(xml_path),outputDirectory=str(output),inputList=str(D/'periods.list'),inputListSha256=sha(D/'periods.list'),sources=c['sources'],fingerprints={str(p.resolve()):sha(p) for p in fixed},pilotComparison=dict(path=str(D/'pilot_comparison.json'),sha256=sha(D/'pilot_comparison.json')),resourceUsageXml=ET.tostring(resource,encoding='unicode'),estimatedAdditionalOutputBytes=15000000000)
    save(D/'launch_config.json',launch)
    print(json.dumps(dict(prepared=True,jobs=20,joblist=str(xml_path),pinned=len(launch['fingerprints']),submitted=False)))
if __name__=='__main__':main()
