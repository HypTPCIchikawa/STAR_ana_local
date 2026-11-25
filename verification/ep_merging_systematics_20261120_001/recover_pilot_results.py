#!/usr/bin/env python3
"""Harvest existing completed pilots after driver final-harvest race; never starts jobs."""
import json, pathlib, re, datetime, hashlib, importlib.util
E=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pilot_driver',E/'run_pilot_matrix.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
pins=json.loads((E/'pilot_matrix_pins.json').read_text());s=json.loads((E/'pilot_matrix_status.json').read_text())
assert len(s['cases'])==8
matrix=json.loads((E/'cases.json').read_text())
source_info=json.loads(pathlib.Path(matrix['inputManifest']).read_text())['sourceFiles'][0]
driver=str(E/'run_pilot_matrix.py')
changes=[]
for path,expected in pins['sourceConfigLibSha256'].items():
 actual=d.sha(path)
 if actual!=expected:
  assert pathlib.Path(path)==E/'run_pilot_matrix.py'
  assert d.sha(E/'run_pilot_matrix_before_harvest_fix.py')==expected
  changes.append(dict(path=path,before=expected,after=actual,reason='Final harvest loop must include completed but unharvested children; no physics code or configuration change.'))
assert len(changes)==1
assert d.stat_record(pins['source']['path'])==pins['source']
ref=pins['referencePilot'];assert d.sha(ref['path'])==ref['sha256']
assert d.sha(E/'cases.json')==pins['casesManifestSha256']
assert d.sha(E/'pilot_incremental_build.log')==pins['buildLogSha256']
for name,row in s['cases'].items():
 log=E/'pilots'/(name+'.log');time=E/'pilots'/(name+'.time');output=pathlib.Path(row['output'])
 markers=[x for x in log.read_text(errors='replace').splitlines() if x.startswith('FemtoLambda tree replay job=')]
 assert len(markers)==1,(name,markers)
 match=re.search(r'requested=(\d+) completed=(\d+) available=(\d+) status=(\d+) RealTime=([0-9.eE+-]+) CpuTime=([0-9.eE+-]+)',markers[0]);assert match
 assert int(match.group(1))==int(match.group(2))==400000
 assert int(match.group(3))==source_info['acceptedEvents'] and int(match.group(4))==0
 resource=d.timing(time);assert resource['timeExitStatus']==0
 if 'returnCode' not in row:
  row.update(returnCode=0,finishedJst=datetime.datetime.fromtimestamp(time.stat().st_mtime,datetime.timezone(datetime.timedelta(hours=9))).isoformat(),completionMarkers=markers,passed=True,eventLoopWallSeconds=float(match.group(5)),eventLoopCpuSeconds=float(match.group(6)),outputStat=d.stat_record(output),**resource)
  row['harvestRecoveredFromCompletedLog']=True
  d.fresh(E/'pilots'/(name+'.result.json'),row)
 else:
  assert row['passed'] and row['returnCode']==0
  assert row['completionMarkers']==markers and row['outputStat']==d.stat_record(output)
  for k,v in resource.items():assert row[k]==v
 assert output.stat().st_size>0
s.update(finishedJst=d.now(),phase='complete',sourceUnchanged=True,sourceConfigLibPinsUnchanged=True,pilotDriverOnlyHarvestFix=changes,passed=True,recovery='Existing completed logs/time/output harvested; no ROOT jobs rerun.')
d.save(E/'pilot_matrix_status.json',s)
d.fresh(E/'pilot_matrix_results.json',{k:v for k,v in s.items() if k!='rssSamples'})
d.fresh(E/'pilots/resource.json',dict(passed=True,maximumConcurrentProcesses=4,cases={name:{key:row.get(key) for key in ['eventLoopWallSeconds','eventLoopCpuSeconds','userSeconds','systemSeconds','maxRssKiB','outputStat']} for name,row in s['cases'].items()},rssSamples=s['rssSamples']))
d.fresh(E/'pilot_recovery.json',dict(passed=True,recoveredJst=d.now(),noRootJobsRerun=True,runtimeSourceConfigLibUnchanged=True,driverChanges=changes,recoveredCases=[name for name,row in s['cases'].items() if row.get('harvestRecoveredFromCompletedLog')]))
print(json.dumps(dict(passed=True,cases={name:{k:row[k] for k in ['eventLoopWallSeconds','maxRssKiB']} for name,row in s['cases'].items()}),indent=2))
