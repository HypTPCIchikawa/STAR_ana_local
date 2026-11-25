#!/usr/bin/env python3
"""Recompute completed 0-60% CF comparisons from saved native raw bins only."""
from pathlib import Path
import json, hashlib, math
from datetime import datetime, timezone
HERE=Path(__file__).resolve().parent
CASE=HERE.parent/'aggregate_watch_001/cf_comparison'
CENT9=range(2,9)
NORM=(.4,.5)
WINDOWS={'k060_100':(.6,1.),'k080_100':(.8,1.)}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def same(a,b):return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
def main():
    files=[CASE/n for n in ['raw_histograms.json','cf_comparison.json','completed.json']]
    pins={str(p):sha(p) for p in files}
    raw,cf,done=[json.loads(p.read_text()) for p in files]
    assert raw['passed'] and cf['passed'] and done['passed'] and done['sources']==160
    assert cf['settings']['normalizationGeVc']==list(NORM)
    assert raw['treeEntriesRead']==0 and cf['treeEntriesRead']==0
    edges=raw['nativeEdges'];n=len(edges)-1
    assert n==200 and edges[0]==0 and edges[-1]==1
    checks=0
    def integral(a,lo,hi):
        selected=[i+1 for i in range(n) if edges[i]>=lo-1e-12 and edges[i+1]<=hi+1e-12]
        assert selected and abs(edges[selected[0]-1]-lo)<1e-12 and abs(edges[selected[-1]]-hi)<1e-12
        return sum(a[i] for i in selected)
    result=[];banks={}
    for case in raw['cases']:
        selected=[s for s in case['slices'] if s['plane']==-1 and s['cent9'] in CENT9]
        assert len(selected)==7
        arrays={field:[sum(s[field][i] for s in selected) for i in range(n+2)] for field in ['SE','ME','SEVariance','MEVariance']}
        banks[case['caseId']]=arrays
        ref=next(c for c in cf['curves'] if c['caseId']==case['caseId'] and c['centralityId']=='c00_60' and c['plane']==-1)
        seH=integral(arrays['SE'],*NORM);meH=integral(arrays['ME'],*NORM);alpha=seH/meH
        assert seH==ref['normSE'] and meH==ref['normME'] and same(alpha,ref['normAlpha']);checks+=3
        counts={field:dict(underflow=arrays[field][0],regular=sum(arrays[field][1:-1]),overflow=arrays[field][-1],allIncludingFlow=sum(arrays[field])) for field in ['SE','ME']}
        for field in ['SE','ME']:
            assert counts[field]['regular']==ref[field+'Regular']
            assert counts[field]['underflow']==ref[field+'Underflow']
            assert counts[field]['overflow']==ref[field+'Overflow'];checks+=3
        windows={}
        for key,w in WINDOWS.items():
            se=integral(arrays['SE'],*w);me=integral(arrays['ME'],*w);ratio=se/(alpha*me)
            windows[key]=dict(windowGeVc=list(w),SE=se,ME=me,alphaScaledME=alpha*me,integralCF=ratio)
            if key in ref['highK']:
                assert ref['highK'][key]['SE']==se and ref['highK'][key]['ME']==me and same(ref['highK'][key]['integralCF'],ratio);checks+=3
        for b in ref['bins']:
            se=integral(arrays['SE'],b['low'],b['high']);me=integral(arrays['ME'],b['low'],b['high'])
            assert se==b['SE'] and me==b['ME'];checks+=2
            if me>0:assert same(b['CF'],se/(alpha*me));checks+=1
        result.append(dict(caseId=case['caseId'],epBins=case['epBins'],cutEnabled=case['cutEnabled'],normSE=seH,normME=meH,normAlpha=alpha,counts=counts,highK=windows))
    for cut in [False,True]:
        rows=[row for row in result if row['cutEnabled']==cut]
        first=banks[rows[0]['caseId']]
        for row in rows[1:]:
            for field in ['SE','SEVariance']:
                assert banks[row['caseId']][field]==first[field];checks+=n+2
    by={x['caseId']:x for x in result};cutrows=[]
    for ep in [1,8,16,32]:
        a=by['ep'+str(ep)+'_mergingOff'];b=by['ep'+str(ep)+'_mergingOn'];ret={}
        for field in ['SE','ME']:
            ret[field]={}
            for population in ['regular','allIncludingFlow']:
                off=a['counts'][field][population];on=b['counts'][field][population]
                ret[field][population]=dict(off=off,on=on,rejected=off-on,rejectedFraction=(off-on)/off,retainedFraction=on/off)
        high={}
        for k in WINDOWS:
            off=a['highK'][k]['integralCF'];on=b['highK'][k]['integralCF']
            high[k]=dict(off=off,on=on,onMinusOff=on-off,onOverOff=on/off,relativeChangePercent=100*(on/off-1))
        cutrows.append(dict(epBins=ep,counts=ret,highK=high))
    assert pins=={str(p):sha(p) for p in files}
    result.sort(key=lambda x:(x['epBins'],x['cutEnabled']))
    output=dict(passed=True,createdUTC=datetime.now(timezone.utc).isoformat(),sourceSHA256=pins,checks=checks,centrality='0-60%',cent9=[2,8],normalizationGeVc=list(NORM),nativeAxisGeVc=[0,1],signal='raw Lambda-He4 signal before background subtraction/purity correction',eventPopulation=cf['eventPopulation'],eventsPerCase=cf['inputEventsPerCase'],sameSEAcrossEPAtFixedCut=True,treeEntriesRead=0,curves=result,cutComparisons=cutrows,limitations=cf['statistics'],summary='EP8 moves the normalized high-k count ratio closer to 1 than EP1. EP16/32 make only small further changes; all remain below 1 in the evaluated windows. Proton-He4 origin cut changes the high-k CF by about -0.07%, descriptively; no significance test.')
    with (HERE/'cf_results_summary.json').open('x') as out:json.dump(output,out,indent=2,ensure_ascii=False);out.write('\n')
    with (HERE/'cf_results_summary.tsv').open('x') as out:
        out.write('case\tEPbins\tcut\tCF_060_100\tCF_080_100\tSEall\tMEall\tSEH\tMEH\n')
        for row in result:out.write('\t'.join(map(str,[row['caseId'],row['epBins'],row['cutEnabled'],row['highK']['k060_100']['integralCF'],row['highK']['k080_100']['integralCF'],row['counts']['SE']['allIncludingFlow'],row['counts']['ME']['allIncludingFlow'],row['normSE'],row['normME']]))+'\n')
    print(json.dumps(dict(passed=True,checks=checks,curves=[dict(caseId=x['caseId'],C060_100=x['highK']['k060_100']['integralCF'],C080_100=x['highK']['k080_100']['integralCF']) for x in result],cuts=cutrows),ensure_ascii=False))
if __name__=='__main__':main()
