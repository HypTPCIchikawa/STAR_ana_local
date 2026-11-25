#!/usr/bin/env python3
"""Build EP8 CF tables from a frozen histogram export; draw with ROOT only."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]


def need(ok,message):
    if not ok:raise RuntimeError(message)


def load(path):
    with Path(path).open() as stream:return json.load(stream)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def resolve(value,base):
    p=Path(value)
    return (p if p.is_absolute() else base/p).resolve()


def stat(path):
    s=path.stat();need(path.is_file() and s.st_size>0,'Missing/empty file: '+str(path))
    return dict(path=str(path),bytes=s.st_size,device=s.st_dev,inode=s.st_ino,mtimeNs=s.st_mtime_ns)


def check_stat(record):
    current=stat(Path(record['path']))
    need(all(current[k]==record[k] for k in ['bytes','device','inode','mtimeNs']),'Input stat changed: '+record['path'])


def pin(path):
    r=stat(path);r['sha256']=sha(path);return r


def save(path,value):
    with path.open('x') as stream:json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')


def f32(x):return struct.unpack('!f',struct.pack('!f',x))[0]
def fkey(x):
    b=struct.unpack('!I',struct.pack('!f',x))[0]
    return ~b&0xffffffff if b&0x80000000 else b^0x80000000
def kval(k):
    b=~k&0xffffffff if k<0x80000000 else k^0x80000000
    return struct.unpack('!f',struct.pack('!I',b))[0]


def epbin(angle):
    half=f32(f32(angle)/2);low,high=-math.pi/2,math.pi/2;width=f32((high-low)/f32(8))
    if not low<half<high:return -1
    b=int((half-low)/width)
    return b if 0<=b<8 else -1


def angles(source):
    body=source.split('inline Int_t XialeiBin(Double_t physicalPsi1) {',1)[1].split('\n}',1)[0]
    need(all(t in body for t in ['const Float_t psi = physicalPsi1;','const Float_t half = psi / 2;',
         'const Float_t width = (high-low)/Float_t(8);','if (!(half > low && half < high)) return -1;',
         'const Int_t bin = Int_t((half-low)/width);']),'Canonical source contract differs')
    lo,hi=fkey(f32(-math.pi)),fkey(f32(math.pi))
    while epbin(kval(lo))<0:lo+=1
    while epbin(kval(hi))<0:hi-=1
    starts=[lo]
    for target in range(1,8):
        a,b=lo,hi
        while a<b:
            m=(a+b)//2
            if epbin(kval(m))>=target:b=m
            else:a=m+1
        starts.append(a)
    result=[]
    for p in range(8):
        last=starts[p+1]-1 if p<7 else hi
        need(epbin(kval(starts[p]))==epbin(kval(last))==p and epbin(kval(starts[p]-1))==(p-1 if p else -1)
             and epbin(kval(last+1))==(p+1 if p<7 else -1),'Canonical float boundary failed')
        result.append(dict(plane=p,nominalLowDeg=-180+45*p,nominalHighDeg=-135+45*p,
                           storedFloatPsi1LowInclusive=kval(starts[p]),storedFloatPsi1HighInclusive=kval(last),
                           nextStoredFloatPsi1UpperExclusive=kval(last+1),
                           label='EP%d: #Psi_{1} %d to %d deg (nominal)'%(p,-180+45*p,-135+45*p)))
    return dict(physicalPlane='corrected first-order EPD Psi1, physical period 2*pi, not Psi2',
                algorithm='Float_t psi; Float_t half=psi/2; Float_t width=pi/8; strict -pi/2<half<pi/2; Int_t((half+pi/2)/width)',
                floatWidthInHalfAngle=f32(math.pi/8),boundaryConvention='Nominal 45-degree labels; exact stored Float_t inclusive endpoints below. Global guard is strict; internal boundaries are lower-inclusive in rounded half-angle space.',ranges=result)


def bounds(edges,window):
    result=[]
    for edge in window:
        found=[i for i,x in enumerate(edges) if math.isclose(edge,x,rel_tol=0,abs_tol=1.e-10)]
        need(len(found)==1,'Window not aligned or outside native axis: '+repr(window));result.append(found[0])
    need(result[0]<result[1],'Reversed interval');return result


def curve(rows,edges,cent,plane,cfg):
    n=len(edges)-1;total={k:[sum(r[k][i] for r in rows) for i in range(n+2)] for k in ['SE','ME','SEVariance','MEVariance']}
    values={k:v[1:-1] for k,v in total.items()};a,b=bounds(edges,cfg['normalizationGeVc'])
    sh,mh=sum(values['SE'][a:b]),sum(values['ME'][a:b]);alpha=sh/mh if sh>0 and mh>0 else None
    out=dict(centralityId=cent['id'],centralityLabel=cent['label'],cent9=cent['cent9'],plane=plane,
             normSE=sh,normME=mh,normAlpha=alpha,available=alpha is not None,highK={},bins=[],
             countsIncludingFlow={k:dict(underflow=total[k][0],regular=sum(values[k]),overflow=total[k][-1],allIncludingFlow=sum(total[k])) for k in ['SE','ME']})
    def ratio(s,m,sv,mv):
        if alpha is None or m<=0:return None,None
        c=s/(alpha*m);v=sv/(alpha*alpha*m*m)+s*s*mv/(alpha*alpha*m**4)
        return c,math.sqrt(max(0,v))
    for name,w in cfg['highKWindowsGeVc'].items():
        a,b=bounds(edges,w);s,m,sv,mv=(sum(values[k][a:b]) for k in ['SE','ME','SEVariance','MEVariance']);c,e=ratio(s,m,sv,mv)
        out['highK'][name]=dict(windowGeVc=w,SE=s,ME=m,SEVariance=sv,MEVariance=mv,integralCF=c,fixedAlphaCountError=e,available=c is not None)
    step=cfg['displayRebin'];need(n%step==0,'Nondivisible display rebin')
    for a in range(0,n,step):
        s,m,sv,mv=(sum(values[k][a:a+step]) for k in ['SE','ME','SEVariance','MEVariance']);c,e=ratio(s,m,sv,mv)
        out['bins'].append(dict(low=edges[a],high=edges[a+step],SE=s,ME=m,SEVariance=sv,MEVariance=mv,CF=c,fixedAlphaCountError=e))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--settings',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--mainconf',type=Path,required=True);ap.add_argument('--plot',action='store_true');args=ap.parse_args();sp=args.settings.resolve();cfg=load(sp)
    need(cfg['schemaVersion']==1,'Unsupported schema');paths={k:resolve(cfg[k],sp.parent) for k in ['rawHistograms','inputProvenance','baselineCompleted','protectedPinProof','independentValidation','canonicalSource']}
    inputs=[pin(p) for p in list(paths.values())+[sp,Path(__file__).resolve(),HERE/'DrawEP8CF.C',resolve(args.mainconf,PROJECT)]]
    raw,prov=load(paths['rawHistograms']),load(paths['inputProvenance']);need(raw['passed'] and raw['treeEntriesRead']==0 and load(paths['baselineCompleted'])['passed'],'Unverified baseline')
    need(len(raw['cases'])==1 and raw['cases'][0]['caseId']==cfg['caseId'],'Wrong case');case=raw['cases'][0]
    need(case['epBins']==8 and case['sources']==len(raw['sources'])==len(prov['sourceStats'])==20,'Baseline coverage differs')
    for r in prov['sourceStats']:check_stat(r)
    protected=load(paths['protectedPinProof'])['pins']
    for r in protected:check_stat(r);need(sha(r['path'])==r['sha256'],'Protected production pin differs: '+r['path'])
    angle=angles(paths['canonicalSource'].read_text());edges=raw['nativeEdges'];n=len(edges)-1;width=edges[1]-edges[0]
    need(n==200 and all(math.isclose(edges[i+1]-edges[i],width,rel_tol=0,abs_tol=1.e-12) for i in range(n)),'Native axis differs')
    lookup={(r['plane'],r['cent9']):r for r in case['slices']};need(len(lookup)==81,'Missing/duplicate native slices');maxdiff=0
    for cent in range(9):
        for key in ['SE','ME']:
            for i in range(n+2):
                inc=lookup[(-1,cent)][key][i];summed=sum(lookup[(p,cent)][key][i] for p in range(8));d=abs(inc-summed);maxdiff=max(maxdiff,d)
                need(d<=cfg['projectionAbsoluteTolerance']+cfg['projectionRelativeTolerance']*max(inc,summed),'Native plane sum closure failed')
    groups=cfg['mainCentralities']+cfg['independentCentralities'];curves=[]
    for cent in groups:
        for p in range(-1,8):curves.append(curve([lookup[(p,c)] for c in range(cent['cent9'][0],cent['cent9'][1]+1)],edges,cent,p,cfg))
    bykey={(c['centralityId'],c['plane']):c for c in curves};independent=load(paths['independentValidation']);crosschecks=0
    need(independent['passed'],'Independent validation not passed')
    for other in independent['curves']:
        p=other['plane'];ours=bykey[('c00_60',p)];need(math.isclose(ours['normAlpha'],other['normAlpha'],rel_tol=1.e-12),'Independent alpha differs')
        for name,w in other['highK'].items():
            match=next((v for v in ours['highK'].values() if v['windowGeVc']==w['windowGeVc']),None);need(match is not None,'Unknown independent window')
            need(match['SE']==w['SE'] and match['ME']==w['ME'] and math.isclose(match['integralCF'],w['integralCF'],rel_tol=1.e-12),'Independent count/CF differs');crosschecks+=1
    summaries=[]
    for cent in groups:
        cc=[bykey[(cent['id'],p)] for p in range(8)];high={}
        for name in cfg['highKWindowsGeVc']:
            vv=[c['highK'][name]['integralCF'] for c in cc if c['highK'][name]['available']]
            high[name]=dict(validPlanes=len(vv),minIntegralCF=min(vv) if vv else None,maxIntegralCF=max(vv) if vv else None,inclusive=bykey[(cent['id'],-1)]['highK'][name])
        summaries.append(dict(centralityId=cent['id'],allEightNormsPositive=all(c['available'] for c in cc),emptyPlaneNorms=[c['plane'] for c in cc if not c['available']],minNormSE=min(c['normSE'] for c in cc),maxNormSE=max(c['normSE'] for c in cc),minNormME=min(c['normME'] for c in cc),maxNormME=max(c['normME'] for c in cc),highK=high))
    totals={k:sum(r[k] for r in raw['sources']) for k in ['inputEvents','completedEvents','epdEventsSeen','epdEventsValid','epdEventsRejected','selectedLambdaCandidates']}
    metrics=dict(schemaVersion=1,passed=True,createdUtc=datetime.now(timezone.utc).isoformat(),caseId=cfg['caseId'],sourceCount=20,nativeBins=n,nativeWidthGeVc=width,displayRebin=cfg['displayRebin'],displayWidthGeVc=width*cfg['displayRebin'],nativeRangeGeVc=[edges[0],edges[-1]],normalizationGeVc=cfg['normalizationGeVc'],eventTotals=totals,canonicalPsi1=angle,planeSumMaxAbsoluteCountDifference=maxdiff,independentWindowCrossChecks=crosschecks,settings=cfg,curves=curves,summaries=summaries,
      definition='Raw signal before background/purity correction; sum raw centrality counts before one norm; alpha=SE(H)/ME(H); CF=SE/(alpha ME). Window values are native raw-count ratios, not arithmetic means of CF points.',
      statisticalLimit='Display/window errors hold alpha fixed. Normalization covariance, shared-event/mixed-pair correlations are omitted. No significance or true asymptotic-convergence claim.',
      eventPopulation='Processed/stored events include EP-invalid events; pair population uses valid corrected EPD plane. See separate seen/valid/rejected totals.',treeEntriesRead=0,ROOTPayloadRead=False)
    out=args.output_dir.resolve();need(HERE in out.parents and not out.exists(),'Output must be a new subdirectory of this study');out.mkdir()
    save(out/'ep8_cf_metrics.json',metrics);save(out/'ep8_source_provenance.json',dict(inputs=inputs,sourceStats=prov['sourceStats'],protectedProductionPins=protected))
    with (out/'ep8_highk_windows.csv').open('x',newline='') as f:
        writer=csv.writer(f);writer.writerow(['centralityId','plane','window','low','high','normSE','normME','alpha','SE','ME','integralCF','fixedAlphaCountError'])
        for c in curves:
            for name,w in c['highK'].items():writer.writerow([c['centralityId'],c['plane'],name]+w['windowGeVc']+[c['normSE'],c['normME'],c['normAlpha'],w['SE'],w['ME'],w['integralCF'],w['fixedAlphaCountError']])
    with (out/'ep8_curve_metadata.tsv').open('x') as f:
        for c in curves:f.write('\t'.join(map(str,[c['centralityId'],c['centralityLabel'],c['plane'],'inclusive raw sum' if c['plane']<0 else angle['ranges'][c['plane']]['label'],c['normSE'],c['normME'],'NA' if c['normAlpha'] is None else c['normAlpha']]))+'\n')
    with (out/'ep8_curve_bins.tsv').open('x') as f:
        for c in curves:
            for i,b in enumerate(c['bins']):f.write('\t'.join(map(str,[c['centralityId'],c['plane'],i,b['low'],b['high'],b['SE'],b['SEVariance'],b['ME'],b['MEVariance'],'NA' if b['CF'] is None else b['CF'],'NA' if b['fixedAlphaCountError'] is None else b['fixedAlphaCountError']]))+'\n')
    with (out/'ep8_plot_windows.tsv').open('x') as f:
        for c in curves:
            for name,w in c['highK'].items():f.write('\t'.join(map(str,[c['centralityId'],c['plane'],name]+w['windowGeVc']+[w['SE'],w['ME'],'NA' if w['integralCF'] is None else w['integralCF'],'NA' if w['fixedAlphaCountError'] is None else w['fixedAlphaCountError']]))+'\n')
    with (out/'ep8_plot_settings.tsv').open('x') as f:
        for c in cfg['mainCentralities']:f.write('mainCentrality\t'+c['id']+'\n')
        for c in cfg['mainCentralities'][1:]:f.write('sixPanel\t'+c['id']+'\n')
        for c in cfg['independentCentralities']:f.write('independent\t'+c['id']+'\n')
        for name in cfg['highKWindowsGeVc']:f.write('window\t'+name+'\n')
    figures=resolve(cfg['figureDirectory'],PROJECT)
    if args.plot:
        figures.mkdir(parents=True,exist_ok=True);inv=str(HERE/'DrawEP8CF.C')+'+('+','.join(json.dumps(str(v)) for v in [out,figures])+')'
        with (out/'ep8_root_plot.log').open('x') as log:subprocess.run(['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(resolve(args.mainconf,PROJECT)),'root4star','-l','-b','-q',inv],cwd=str(PROJECT),stdout=log,stderr=subprocess.STDOUT,check=True)
    for r in inputs:check_stat(r);need(sha(r['path'])==r['sha256'],'Input SHA changed: '+r['path'])
    for r in prov['sourceStats']:check_stat(r)
    for r in protected:check_stat(r);need(sha(r['path'])==r['sha256'],'Protected production pin changed: '+r['path'])
    artifacts=[pin(p) for p in sorted(out.iterdir()) if p.is_file()]
    if args.plot:artifacts += [pin(p) for p in sorted(figures.glob('ep8_*')) if p.suffix in ['.pdf','.png']]
    save(out/'ep8_completed.json',dict(passed=True,plotted=args.plot,inputPinsUnchanged=True,sourceStatsUnchanged=True,protectedProductionPinsUnchanged=True,sourceRoots=20,protectedPinCount=len(protected),treeEntriesRead=0,artifacts=artifacts))
    print(json.dumps(dict(passed=True,output=str(out),eventTotals=totals,summaries=summaries),indent=2))


if __name__=='__main__':
    try:main()
    except Exception as error:print('EP8_CF_FAILED: '+str(error),file=sys.stderr);sys.exit(1)
