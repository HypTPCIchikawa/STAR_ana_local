#!/usr/bin/env python3
import pathlib,json,hashlib,re,datetime
E=pathlib.Path(__file__).resolve().parent;P=E.parent/'femto_lambda_period_mix100_20261119_001'
a=json.loads((P/'completion_log_audit_20261120.json').read_text());assert a['passed'] and len(a['workers'])==20
new=json.loads((E/'pilot_matrix_results.json').read_text())
old=json.loads((P/'pilot_status.json').read_text())
marker=old['pilots']['B_period100']['completionMarkers'][0]
old_cpu=float(re.search(r'CpuTime=([0-9.]+)',marker).group(1));old_wall=float(re.search(r'RealTime=([0-9.]+)',marker).group(1))
checks=[]
for w in a['workers']:
 p=pathlib.Path(w['log']);blob=p.read_bytes();assert hashlib.sha256(blob).hexdigest()==w['logSha256']
 lines=[s for s in blob.decode(errors='replace').splitlines() if s.startswith('FemtoLambda tree replay job=')];assert len(lines)==1
 m=re.search(r'requested=(\d+) completed=(\d+) available=(\d+) status=(\d+) RealTime=([0-9.eE+-]+) CpuTime=([0-9.eE+-]+)',lines[0]);assert m
 assert int(m.group(1))==int(m.group(2))==int(m.group(3))==w['expectedStoredEvents']
 assert int(m.group(4))==0 and float(m.group(5))==w['wallSeconds'] and float(m.group(6))==w['cpuSeconds']
 checks.append(dict(group=w['group'],cpuSeconds=w['cpuSeconds'],wallSeconds=w['wallSeconds'],storedEvents=w['completed'],log=str(p),logSha256=w['logSha256']))
long_cpu=max(checks,key=lambda r:r['cpuSeconds']);long_wall=max(checks,key=lambda r:r['wallSeconds'])
ratios={}
for name,row in new['cases'].items():
 ratios[name]=dict(pilotCpuRatioToOldEp8=row['eventLoopCpuSeconds']/old_cpu,pilotWallRatioToOldEp8=row['eventLoopWallSeconds']/old_wall,predictedMaxPeriodCpuHours=long_cpu['cpuSeconds']*row['eventLoopCpuSeconds']/old_cpu/3600.,predictedMaxPeriodWallHours=long_wall['wallSeconds']*row['eventLoopWallSeconds']/old_wall/3600.)
max_cpu=max(v['predictedMaxPeriodCpuHours'] for v in ratios.values());max_wall=max(v['predictedMaxPeriodWallHours'] for v in ratios.values())
# Additional conservative comparison: 400k prefix direct per-event extrapolation.
prefix=json.loads((E/'pilot_resource_projection.json').read_text())
max_wall=max(max_wall,max(v['linearLongestPeriodWallHours'] for v in prefix['cases'].values()))
assert max_cpu<48 and max_wall<48
record=dict(passed=True,createdUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),baselineActualLogsChecked=20,baselineAuditSha256=hashlib.sha256((P/'completion_log_audit_20261120.json').read_bytes()).hexdigest(),baselineMaxCpuPeriod=long_cpu,baselineMaxWallPeriod=long_wall,oldEp8PilotCpuSeconds=old_cpu,oldEp8PilotWallSeconds=old_wall,cases=ratios,conservativePredictedMaxCpuHours=max_cpu,conservativePredictedMaxWallHours=max_wall,schedulerEstimateHours=48,defaultHardLimitHours=72,hardLimitToPredictedCpuRatio=72/max_cpu,notes=['All twenty old full-period log hashes and exact completion markers match the prior audit.','Pilot ratios include current QA overhead and EP32 cost; four-way pilot CPU contention and period candidate density make projections approximate.','48h scheduler estimate is about 1.54 times the max CPU prediction; the existing default72h is about2.3 times. No hard-limit override is requested.'])
with (E/'pilot_completed_baseline_resource_projection.json').open('x') as out:json.dump(record,out,indent=2);out.write('\n')
print(json.dumps({k:record[k] for k in ['passed','baselineActualLogsChecked','conservativePredictedMaxCpuHours','conservativePredictedMaxWallHours','hardLimitToPredictedCpuRatio']},indent=2))
