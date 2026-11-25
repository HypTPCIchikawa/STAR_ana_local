#!/usr/bin/env python3
"""Assemble the production gate from completed, independent replay checks."""
import json, pathlib, hashlib, datetime
E=pathlib.Path(__file__).resolve().parent
R=E.parents[1]
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,v):
 with pathlib.Path(p).open('x') as f:json.dump(v,f,indent=2);f.write('\n')
def record(p):
 p=pathlib.Path(p);s=p.stat();return dict(path=str(p),bytes=s.st_size,mtimeNs=s.st_mtime_ns,sha256=sha(p))
names=['pilot_matrix_results.json','pilot_recovery.json','pilot_ep_comparison.json','pilot_merging_accounting.json','test_actual_lambda_merging_paircut_result.json','pilot_baseline_counts.json','config_validation.json']
reports={name:json.loads((E/name).read_text()) for name in names}
assert all(row['passed'] is True for row in reports.values())
pins=json.loads((E/'pilot_matrix_pins.json').read_text());recovery=reports['pilot_recovery.json']
changes={row['path']:row for row in recovery['driverChanges']}
readme=str(R/'StMaker/StFemtoMaker/README.md')
changes[readme]=dict(path=readme,before=pins['sourceConfigLibSha256'][readme],after=sha(readme),reason='Parent updated explanatory EP1/8/16/32 and cut-QA documentation after pilots; this Markdown file is not part of runtime code, configuration, or libraries.')
final=[]
for path,expected in pins['sourceConfigLibSha256'].items():
 current=sha(path)
 if current!=expected:
  assert path in changes and changes[path]['before']==expected and changes[path]['after']==current
  assert pathlib.Path(path) in [E/'run_pilot_matrix.py',R/'StMaker/StFemtoMaker/README.md']
 final.append(dict(path=path,sha256=current,pilotSha256=expected,unchanged=current==expected))
ref=pins['referencePilot'];assert sha(ref['path'])==ref['sha256']
src=pathlib.Path(pins['source']['path']);s=src.stat()
assert all(pins['source'][k]==v for k,v in [('bytes',s.st_size),('inode',s.st_ino),('device',s.st_dev),('mtimeNs',s.st_mtime_ns)])
assert sha(E/'cases.json')==pins['casesManifestSha256']
for row in reports['test_actual_lambda_merging_paircut_result.json']['sourceAndLibraryFingerprints']:
 assert sha(R/row['path'])==row['sha256'],row['path']
outputs=[record(pathlib.Path(row['output'])) for row in reports['pilot_matrix_results.json']['cases'].values()]
ep_log=(E/'test_ep_bins.log').read_text()
assert 'PASS' in ep_log and 'Error:' not in ep_log and 'FAIL' not in ep_log
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
finalpins=dict(passed=True,createdUtc=now,runtimeSourceConfigLibrariesUnchanged=True,source=pins['source'],referencePilot=pins['referencePilot'],sourceConfigLibSha256={row['path']:row['sha256'] for row in final},pilotRuntimeSha256={row['path']:row['pilotSha256'] for row in final},explicitNonRuntimeChanges=list(changes.values()),casesManifestSha256=sha(E/'cases.json'),reportPins=[record(E/name) for name in names],pilotOutputs=outputs)
save(E/'pilot_final_pins.json',finalpins)
matrix=json.loads((E/'cases.json').read_text());manifest=json.loads(pathlib.Path(matrix['inputManifest']).read_text())
maxperiod=max(manifest['sourceFiles'],key=lambda x:x['acceptedEvents']);events=maxperiod['acceptedEvents'];scale=events/400000.
resources={}
for case,row in reports['pilot_matrix_results.json']['cases'].items():
 resources[case]=dict(pilotLoopCpuSeconds=row['eventLoopCpuSeconds'],pilotLoopWallSeconds=row['eventLoopWallSeconds'],pilotProcessCpuSeconds=row['userSeconds']+row['systemSeconds'],maxRssKiB=row['maxRssKiB'],linearLongestPeriodLoopCpuHours=row['eventLoopCpuSeconds']*scale/3600.,linearLongestPeriodWallHours=row['eventLoopWallSeconds']*scale/3600.,linearLongestPeriodProcessCpuHours=(row['userSeconds']+row['systemSeconds'])*scale/3600.)
save(E/'pilot_resource_projection.json',dict(passed=True,pilotEvents=400000,longestPeriod=maxperiod['path'],longestPeriodStoredEvents=events,cases=resources,recommendedSchedulerEstimateHours=48,existingDefaultHardLimitHours=72,caveats='A 400000-event prefix includes cold pools, spans three candidate streams, and ran at four-way concurrency. Linear scaling is approximate; period candidate density and CPU contention may differ. No hard-limit override is inferred.'))
ep=reports['pilot_ep_comparison.json'];cut=reports['pilot_merging_accounting.json']
gate=dict(passed=True,createdUtc=now,pilotCases=8,numericalEpTest=dict(passed=True,log=record(E/'test_ep_bins.log'),coverage='N=1/8/16/32, invalid values, exact Float_t hierarchy, all32 boundaries with nextafter, 20000-point deterministic grid; existing original assertions retained.'),storedEventsPerCase=400000,maximumConcurrentProcesses=4,runtimeSourceConfigLibrariesUnchanged=True,noRootJobsRerunForDriverRecovery=True,allOldEp8OffHistogramsExact=ep['oldEP8Off'],crossEpSameCutSE=ep['crossEpSameCutSE'],sameFullScalarLedgerAcrossEightCasesAndOldBaseline=True,segmentStartEvents=ep['segmentStartEvents'],epdEventsSeen=ep['epdEventsSeen'],epdEventsValid=ep['epdEventsValid'],epdEventsRejected=ep['epdEventsRejected'],allEightPlaneClosureAndCanonicalPoolHierarchy=True,cutAccountingChecks=cut['checks'],realPilotNonemptySignalRejectionSE=cut['nonemptySignalRejectionSE'],realPilotNonemptySignalRejectionME=cut['nonemptySignalRejectionME'],syntheticActualLegacyNonemptySEAndBothMEDirections=True,unreachablePolicy='keep',notes=['The 400k real-data SE signal has zero strict-box rejections; this is permitted. Actual-Legacy synthetic integration separately exercises nonempty SE and both ME directions.','Previous chat-summary 100116 is unsubstantiated by prior study JSON/log. Previous chat-summary 607 was the number of compared SE histograms, not SE pair entries. Direct old/new full He4 SE=560 ME=116183, signal SE=256 ME=51337.','The pilot-driver harvest race affected final bookkeeping only; completed output/log/time files were recovered without replay. Driver before/after hashes are recorded; runtime physics code/config/libs are unchanged.'],reports={name:record(E/name) for name in names},finalPins=record(E/'pilot_final_pins.json'),resources=record(E/'pilots/resource.json'),resourceProjection=record(E/'pilot_resource_projection.json'))
save(E/'pilot_comparison.json',gate)
print(json.dumps(dict(passed=True,gate=str(E/'pilot_comparison.json'),runtimePins=len(final),oldHistogramComparisons=ep['oldEP8Off']['histograms'],largestPeriodN32WallHours=max(x['linearLongestPeriodWallHours'] for x in resources.values())),indent=2))
