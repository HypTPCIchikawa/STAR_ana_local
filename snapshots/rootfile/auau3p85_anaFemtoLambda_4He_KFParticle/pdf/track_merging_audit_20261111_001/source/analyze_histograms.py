#!/usr/bin/env python3
"""Audit the stored close-pair QA without interpreting raw SE as an efficiency."""
import argparse, csv, hashlib, json, math
from pathlib import Path


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def at(h,x,y):return h['content'][x+h['shape'][1]*y]
def edge_index(edges,value):
 rows=[i for i,x in enumerate(edges) if abs(x-value)<1e-9]
 if len(rows)!=1:raise ValueError('Boundary is not a histogram edge: '+str(value))
 return rows[0]
def square(h,w):
 xa,xb=edge_index(h['xEdges'],-w)+1,edge_index(h['xEdges'],w)+1
 ya,yb=edge_index(h['yEdges'],-w)+1,edge_index(h['yEdges'],w)+1
 return sum(at(h,x,y) for y in range(ya,yb) for x in range(xa,xb))
def audit(h,reference,cfg):
 nx,ny=h['shape'][1]-2,h['shape'][0]-2
 reg=sum(at(h,x,y) for y in range(1,ny+1) for x in range(1,nx+1))
 total=sum(h['content']);maximum=max(at(h,x,y) for y in range(1,ny+1) for x in range(1,nx+1))
 maxima=[{'x':(h['xEdges'][x-1]+h['xEdges'][x])/2,'y':(h['yEdges'][y-1]+h['yEdges'][y])/2,'content':maximum} for y in range(1,ny+1) for x in range(1,nx+1) if at(h,x,y)==maximum]
 boxes=[]
 for w in cfg['origin_halfwidths']:
  c=square(h,w);outer=square(h,2*w);ring=outer-c;area=(2*w)**2;rarea=(4*w)**2-area
  ratio=c/area/(ring/rarea) if ring>0 else None
  ratioerr=ratio*math.sqrt(1/c+1/ring) if c>0 and ring>0 else None
  boxes.append({'halfWidth':w,'counts':c,'area':area,'density':c/area,'adjacentSquareRingCounts':ring,'adjacentSquareRingArea':rarea,'adjacentSquareRingDensity':ring/rarea,'centerToAdjacentDensityRatio':ratio,'countOnlyRatioError':ratioerr})
 return {'name':h['name'],'path':h.get('path'), 'entries':h['entries'],'regularCounts':reg,'flowCounts':total-reg,'sumAllBins':total,'entriesMinusSumBins':h['entries']-total,'regularFraction':reg/total if total>0 else None,'signalPairEntries':reference['entries'],'fillFractionRelativeToSignalPairs':h['entries']/reference['entries'] if reference['entries']>0 else None,'xFlowCounts':sum(at(h,x,y) for y in range(ny+2) for x in [0,nx+1]),'yFlowCounts':sum(at(h,x,y) for y in [0,ny+1] for x in range(nx+2)),'originBoxes':boxes,'maxima':maxima,'xEdges':h['xEdges'],'yEdges':h['yEdges']}

def main():
 p=argparse.ArgumentParser();p.add_argument('--settings',type=Path,required=True);a=p.parse_args();out=a.settings.parent.resolve();cfg=json.loads(a.settings.read_text());doc=json.loads((out/cfg['input']).read_text());hs={h['name']:h for h in doc['histograms'].values()}
 reports=[];histlist=[]
 for sp,suffix in cfg['species'].items():
  ref=hs['hKstarSE_lambda_'+suffix+'_signal']
  for variant,ending in cfg['variants'].items():
   h=hs['hDphiDeta_proton_'+sp+ending]
   if h['dimension']!=2 or h['shape']!=[len(h['yEdges'])+1,len(h['xEdges'])+1]:raise ValueError('Unexpected histogram shape')
   if not all(math.isfinite(x) and x>=0 for x in h['content']+h['error2']):raise ValueError('Invalid content or variance')
   if len(h['content'])!=h['shape'][0]*h['shape'][1] or len(h['error2'])!=len(h['content']):raise ValueError('Wrong payload length')
   rep=audit(h,ref,cfg)
   if abs(rep['entriesMinusSumBins'])>1e-6:raise ValueError('Raw-count Entries closure failed')
   rep['meanRegularBinCount']=rep['regularCounts']/((len(h['xEdges'])-1)*(len(h['yEdges'])-1))
   rep.update(species=sp,variant=variant,id=len(reports));reports.append(rep);histlist.append(h)
 same=[]
 for variant,ending in cfg['variants'].items():
  x,y=hs['hDphiDeta_proton_3He'+ending],hs['hDphiDeta_proton_3HeAs4He'+ending]
  same.append({'variant':variant,'sameAxes':x['xEdges']==y['xEdges'] and x['yEdges']==y['yEdges'],'exactSameContent':x['content']==y['content'],'exactSameError2':x['error2']==y['error2'],'sameEntries':x['entries']==y['entries'],'maxContentDifference':max(abs(a-b) for a,b in zip(x['content'],y['content']))})
 if not all(all(x[k] for k in ['sameAxes','exactSameContent','exactSameError2','sameEntries']) for x in same):raise ValueError('Mass mapping changed geometry QA')
 with (out/'hist_meta.tsv').open('w') as f:
  fields=['id','species','variant','name','nx','xmin','xmax','ny','ymin','ymax','entries','regular','regular_fraction','signal_entries','fill_fraction','origin005','origin010','origin020']
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
  for r,h in zip(reports,histlist):w.writerow([r['id'],r['species'],r['variant'],r['name'],len(h['xEdges'])-1,h['xEdges'][0],h['xEdges'][-1],len(h['yEdges'])-1,h['yEdges'][0],h['yEdges'][-1],r['entries'],r['regularCounts'],r['regularFraction'],r['signalPairEntries'],r['fillFractionRelativeToSignalPairs']]+[b['counts'] for b in r['originBoxes']])
 with (out/'hist_cells.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['id','ix','iy','content','error2'])
  for r,h in zip(reports,histlist):
   for y in range(h['shape'][0]):
    for x in range(h['shape'][1]):
     k=x+h['shape'][1]*y;w.writerow([r['id'],x,y,h['content'][k],h['error2'][k]])
 result={'status':'audited','settings':cfg,'sourceRoot':doc['input'],'sourceUuid':doc['uuid'],'rootVersion':doc['rootVersion'],'histograms':reports,'sameHe3PairsMassMappingChecks':same,'caveats':['Raw same-event close-pair densities are not a two-track efficiency or SE/ME correlation.','No mixed-event counterpart, kstar/centrality/rigidity slicing, pion-daughter counterpart, or skipped-primary-momentum counters are stored in these six histograms.','xFlowCounts and yFlowCounts overlap in corner bins; do not sum them.','Center-to-neighbor ratios are local descriptive shape measures, not formal merging/splitting tests.','Entries include angular flows and are compared to all-kstar signal Entries.']}
 write(out/'hist_audit.json',result)
 write(out/'hist_analysis_validation.json',{'passed':True,'histogramCount':6,'sameMappedHe3Checks':same,'inputSha256':hashlib.sha256((out/cfg['input']).read_bytes()).hexdigest(),'settingsSha256':hashlib.sha256(a.settings.read_bytes()).hexdigest(),'sourceSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
 for h in reports:print(h['name'],'entries',h['entries'],'regular',h['regularCounts'],'fraction',h['regularFraction'],'fill',h['fillFractionRelativeToSignalPairs'],'origin',[(b['halfWidth'],b['counts'],b['centerToAdjacentDensityRatio']) for b in h['originBoxes']])
 print('MERGING_NUMERICS_PASS')
if __name__=='__main__':main()
