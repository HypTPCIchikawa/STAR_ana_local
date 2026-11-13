#!/usr/bin/env python3
"""Measure high-k* baselines and centrality-composition effects without refitting."""
import argparse, hashlib, json, math
from pathlib import Path
import yaml
from audit_normalization_ranges import raw_ratio, raw_cells, plateau, aligned_indices
from method6_statistics import analyze_group


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--settings',type=Path,required=True);a=p.parse_args()
 cfg=yaml.safe_load(a.settings.read_text());out=Path(cfg['output_dir']).resolve(strict=True)
 if (out/'baseline_audit.json').exists():raise FileExistsError(out/'baseline_audit.json')
 data={k:json.loads(Path(cfg[k+'_input']).read_text()) for k in ['group','fine']}
 assert all(d['fitOptions']=='SLQNR0' for d in data.values())
 norm=cfg['normalization'];fine={g['id']:g for g in data['fine']['groups']};metrics={};plot=[];raw={}
 def add(kind,gid,label,series,points):
  for x in points:
   value=x.get('cf');err=x.get('cfError')
   if isinstance(value,dict):value,err=value['value'],value['error']
   plot.append([kind,gid,label,series,x['low'],x['high'],value,err])
 for kind,d in data.items():
  raw[kind]={};metrics[kind]=[]
  for g in d['groups']:
   label=cfg['labels'][g['id']];points=raw_ratio(g,d['settings'],norm);raw[kind][g['id']]=points
   checks=[plateau(points,r) for r in cfg['slope_ranges']]
   edges=g['kstarEdges'];cells=raw_cells(g,d['settings']);ni=aligned_indices(edges,norm)
   total={s:sum(cells[s][i]['value'] for i in ni) for s in ['SE','ME']}
   windows=[]
   for r in cfg['window_ranges']:
    inds=aligned_indices(edges,r);n={s:sum(cells[s][i]['value'] for i in inds) for s in ['SE','ME']}
    value=n['SE']/n['ME']*total['ME']/total['SE'] if n['ME']>0 and total['SE']>0 else None
    windows.append(dict(range=r,counts=n,normalizedRawRatio=value))
   metrics[kind].append(dict(id=g['id'],label=label,slopes=checks,windows=windows))
   add(kind,g['id'],label,'raw',points)
 # The original grouped L-derived results retain the original covariance calculations.
 grouped=json.loads(Path(cfg['group_results']).read_text())
 for g in grouped['groups']:
  for method in ['method3','method5','method6']:add('group',g['id'],cfg['labels'][g['id']],method,g['methods'][method]['points'])
 # Fine M3 makes the centrality diagnostic independent of summed-bin composition.
 fine_results=[]
 for g in data['fine']['groups']:
  r=analyze_group(g,data['fine']['settings'],norm,cfg['plot_range'],cfg['slope_ranges'][0]);m=r['methods']['method3']
  add('fine',g['id'],cfg['labels'][g['id']],'method3',m['points'])
  fine_results.append(dict(id=g['id'],normalization=m['normalization'],validPlotBins=m['valid_plot_bins'],points=m['points']))
 # Per-cent exposure matching is a diagnostic only; count/norm derivatives include reuse.
 composition=[]
 for g in data['group']['groups']:
  selected=[fine['cent'+str(i)] for i in range(g['centMin'],g['centMax']+1)]
  cc=[raw_cells(c,data['fine']['settings']) for c in selected];ni=aligned_indices(g['kstarEdges'],norm)
  totals=[{s:sum(c[s][i]['value'] for i in ni) for s in ['SE','ME']} for c in cc]
  if any(x['SE']<=0 or x['ME']<=0 for x in totals):raise ValueError('Empty fine normalization cell')
  points=[]
  for i,(lo,hi) in enumerate(zip(g['kstarEdges'],g['kstarEdges'][1:])):
   for side in ['SE','ME']:
    assert sum(c[side][i]['value'] for c in cc)==raw['group'][g['id']][i][side]
   numerator=sum(c['SE'][i]['value'] for c in cc)
   denominator=sum(c['ME'][i]['value']*n['SE']/n['ME'] for c,n in zip(cc,totals))
   if denominator<=0:continue
   value=numerator/denominator;variance=0
   for c,n in zip(cc,totals):
    mi=c['ME'][i]['value']
    for j in range(len(g['kstarEdges'])-1):
     indicator=1. if j in ni else 0.
     ds=(1. if i==j else 0.)/denominator-numerator/denominator**2*mi/n['ME']*indicator
     dm=-numerator/denominator**2*(n['SE']/n['ME']*(1. if i==j else 0.)-mi*n['SE']/n['ME']**2*indicator)
     variance+=ds*ds*c['SE'][j]['variance']+dm*dm*c['ME'][j]['variance']
   points.append(dict(low=lo,high=hi,cf=value,cfError=math.sqrt(variance)))
  add('group',g['id'],cfg['labels'][g['id']],'centExposureMatched',points)
  composition.append(dict(id=g['id'],fineCentWeights=[dict(id=c['id'],SEnorm=n['SE'],MEnorm=n['ME'],MEweight=n['SE']/n['ME']) for c,n in zip(selected,totals)],points=points))
 def cell(v):return 'nan' if v is None else str(v).replace('\t',' ')
 (out/'baseline_points.tsv').write_text('kind\tgroup\tlabel\tseries\tlow\thigh\tvalue\terror\n'+''.join('\t'.join(cell(v) for v in row)+'\n' for row in plot))
 dump(out/'baseline_audit.json',dict(phase='analyzed',settings=cfg,sourceHashes={str(Path(cfg[k]).resolve()):sha(cfg[k]) for k in ['group_input','fine_input','group_results']},raw=metrics,composition=composition,fineMethod3=fine_results,caveats=['Normalization 0.4-0.5 is a reference scale, not a validated plateau.','Slope significances use diagonal count-only errors and are descriptive, not calibrated tests.','Per-cent exposure matching is a diagnostic, not an adopted physics correction.','Flow, eventwise conservation, run/vertex dependence and pair/event covariance remain to be checked.']))
 for kind,rows in metrics.items():
  for r in rows:
   x=next(v for v in r['windows'] if v['range']==cfg['window_ranges'][-1]);d=r['slopes'][0]
   print(kind,r['label'],'rawHigh=',x['normalizedRawRatio'],'slope=',d.get('slope'),'error=',d.get('slopeError'))
 print('Saved',len(plot),'plot points')

if __name__=='__main__':main()
