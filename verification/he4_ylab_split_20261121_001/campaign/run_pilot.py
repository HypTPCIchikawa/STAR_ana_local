#!/usr/bin/env python3
"""Run one explicitly limited period-prefix pilot and pin its physics inputs."""
import argparse,hashlib,json,re,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
D=Path(__file__).resolve().parent;ROOT=D.parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stored-events',type=int,required=True);a=ap.parse_args();assert a.stored_events>357247
    c=json.loads((D/'analysis_campaign.json').read_text());row=c['sources'][0];src=Path(row['path']);s=src.stat();assert (s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev)==(row['bytes'],row['mtimeNs'],row['inode'],row['device'])
    out=D/('pilot_'+str(a.stored_events)+'.root');assert not out.exists()
    fixed=[ROOT/p for p in ['StMaker/common/FemtoLambdaLegacy.cxx','StMaker/common/FemtoLambdaLegacy.h','StMaker/StFemtoMaker/StFemtoMaker.cxx','StMaker/StFemtoMaker/StFemtoMaker.h','analysis/anaFemtoLambda_4He_tree.C','analysis/run_anaFemtoLambda_4He_tree.C','include/FemtoLambdaReducedTree.h','include/FemtoLambdaReducedTreeConfig.h','include/FemtoLambdaTreeTypes.h','include/FemtoLambdaTreeLinkDef.h']]
    fixed += [Path(p) for p in c['generatedFiles']]+[Path(p) for p in c['originalReferenceConfigSha256']]
    fixed += [ROOT/'lib'/n for n in ['libStarAnaConfig.so','libStRefMultCorr.so','libKFParticle.so','libStKfParticleCommon.so','libStCommon.so','libStFemtoMaker.so']]
    pins={str(p):sha(p) for p in fixed};cmd=['/usr/bin/time','-v','-o',str(D/'pilot.time'),'bash','script/singularity_run_anaFemtoLambda_4He_tree.sh',c['mainconf'],str(src),str(out),'he4YLab8Pilot_'+str(a.stored_events),str(a.stored_events)]
    with (D/'pilot_command.json').open('x') as f:json.dump(dict(command=cmd,source=row,pins=pins,startedUTC=datetime.now(timezone.utc).isoformat()),f,indent=2)
    with (D/'pilot.log').open('x') as f:result=subprocess.run(cmd,cwd=str(ROOT),stdout=f,stderr=subprocess.STDOUT)
    text=(D/'pilot.log').read_text(errors='replace');markers=[x for x in text.splitlines() if x.startswith('FemtoLambda tree replay job=')]
    passed=result.returncode==0 and len(markers)==1 and ('requested='+str(a.stored_events)+' ') in markers[0] and ('completed='+str(a.stored_events)+' ') in markers[0] and 'status=0 ' in markers[0]
    assert all(sha(p)==v for p,v in pins.items()),'Physics code/config/lib changed during pilot'
    t=src.stat();assert (t.st_size,t.st_mtime_ns,t.st_ino,t.st_dev)==(row['bytes'],row['mtimeNs'],row['inode'],row['device'])
    (D/'pilot.list').write_text(str(out)+'\n')
    record=dict(passed=passed,returnCode=result.returncode,completedUTC=datetime.now(timezone.utc).isoformat(),events=a.stored_events,markers=markers,output=str(out),outputBytes=out.stat().st_size if out.exists() else None,physicsPinsUnchanged=True,sourceUnchanged=True)
    (D/'pilot_run_result.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record));return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
