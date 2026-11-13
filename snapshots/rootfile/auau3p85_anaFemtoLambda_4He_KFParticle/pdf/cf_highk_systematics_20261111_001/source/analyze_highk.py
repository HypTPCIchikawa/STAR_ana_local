#!/usr/bin/env python3
"""Project recorded 1D PID selections and propagate normalization covariance."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def read_histograms(path):
    doc = json.loads(path.read_text())
    hs = {}
    for key, h in doc['histograms'].items():
        name = h['name']
        if name in hs:
            raise ValueError('Ambiguous histogram name: ' + name)
        hs[name] = h
    return doc, hs


def matrix(h, key='content'):
    cache = '_array_' + key
    if cache not in h:
        h[cache] = np.array(h[key], dtype=float).reshape(h['shape'])
    return h[cache]


def aligned(edges, value):
    found = np.flatnonzero(np.isclose(edges, value, rtol=0, atol=1.e-9))
    if len(found) != 1:
        raise ValueError('Selection is not on a stored axis edge: ' + str(value))
    return int(found[0])


def selected_y(h, window):
    n = len(h['yEdges']) - 1
    if window is None:
        return slice(0, n+2)
    lo, hi = window
    return slice(aligned(h['yEdges'], lo)+1, aligned(h['yEdges'], hi)+1)


def project(h, window):
    s = selected_y(h, window)
    return np.array([matrix(h) [s, :].sum(axis=0), matrix(h, 'error2')[s, :].sum(axis=0)])


def region(edges, bounds):
    return slice(aligned(edges, bounds[0]), aligned(edges, bounds[1]))


def ratio(s, m, vs, vm, norm, target):
    sn, mn = s[norm].sum(), m[norm].sum()
    sw, mw = s[target].sum(), m[target].sum()
    if min(sn, mn, sw, mw) <= 0:
        return None, None, None
    c = sw / mw * mn / sn
    ns, ts = np.zeros(len(s)), np.zeros(len(s))
    ns[norm] = 1
    ts[target] = 1
    ds = c * (ts/sw - ns/sn)
    dm = c * (ns/mn - ts/mw)
    variance = float((ds*ds*vs + dm*dm*vm).sum())
    return float(c), math.sqrt(max(0., variance)), (ds, dm)


def analyze(sfull, mfull, edges, cfg, norm_bounds):
    factor = cfg['rebin_factor']
    n = len(edges)-1
    if n % factor:
        raise ValueError('Incomplete final kstar bin')
    rb = lambda a: a[1:-1].reshape(-1, factor).sum(axis=1)
    s, vs = rb(sfull[0]), rb(sfull[1])
    m, vm = rb(mfull[0]), rb(mfull[1])
    e = np.asarray(edges)[::factor]
    norm = region(e, norm_bounds)
    points, jacobians = [], []
    for i in range(len(s)):
        c, err, grad = ratio(s,m,vs,vm,norm,slice(i,i+1))
        points.append({'low':float(e[i]),'high':float(e[i+1]),'value':c,'error':err})
        jacobians.append(grad)
    metrics = {'SEnorm':float(s[norm].sum()),'MEnorm':float(m[norm].sum()),
               'SEtotal':float(s.sum()),'MEtotal':float(m.sum()),
               'SEflow':[float(sfull[0,0]),float(sfull[0,-1])],
               'MEflow':[float(mfull[0,0]),float(mfull[0,-1])],
               'normalization':norm_bounds,'windows':{}}
    for name, bounds in cfg['high_windows'].items():
        c, err, _ = ratio(s,m,vs,vm,norm,region(e,bounds))
        metrics['windows'][name] = {'range':bounds,'value':c,'error':err}
    fit = [i for i,p in enumerate(points) if p['value'] is not None and
           p['low'] >= cfg['slope_range'][0]-1.e-9 and p['high'] <= cfg['slope_range'][1]+1.e-9]
    if len(fit) >= 3:
        js = np.array([jacobians[i][0] for i in fit])
        jm = np.array([jacobians[i][1] for i in fit])
        cov = (js*vs).dot(js.T) + (jm*vm).dot(jm.T)
        inv = np.linalg.pinv(cov, rcond=1.e-12)
        x = np.array([(points[i]['low']+points[i]['high'])/2 for i in fit])
        y = np.array([points[i]['value'] for i in fit])
        design = np.array([np.ones(len(x)),x-float(x.mean())]).T
        pcov = np.linalg.inv(design.T.dot(inv).dot(design))
        pars = pcov.dot(design.T).dot(inv).dot(y)
        residual = y-design.dot(pars)
        metrics.update(slope=float(pars[1]),slope_error=math.sqrt(float(pcov[1,1])),
                       slope_chi2=float(residual.dot(inv).dot(residual)),slope_ndf=len(fit)-2,
                       slope_range=cfg['slope_range'],slope_bin_count=len(fit),slope_bin_edges=[[points[i]['low'],points[i]['high']] for i in fit],intercept_at=float(x.mean()),intercept=float(pars[0]))
    else:
        metrics.update(slope=None,slope_error=None,slope_chi2=None,slope_ndf=None)
    return points, metrics


def main():
    p=argparse.ArgumentParser();p.add_argument('--settings',type=Path,required=True);a=p.parse_args()
    cfg=json.loads(a.settings.read_text());out=a.settings.parent.resolve()
    doc,hs=read_histograms(out/cfg['input_export'])
    olddoc,oldhs=read_histograms(out/cfg['old_export'])
    closure=[];series=[];spectra={}
    suffixes=cfg['mass_suffixes'];baseline_species=cfg['baseline_species']
    def check(label,a,b,legacy_float=False):
        if a.shape != b.shape:raise ValueError('Shape mismatch: '+label)
        diff=np.abs(a-b);rtol=cfg['closure_rtol']
        tol=1.e-9+rtol*np.maximum(np.abs(a),np.abs(b))
        if legacy_float:
            for flow in [0,-1]:
                ulp=np.abs(np.spacing(b[:,flow].astype(np.float32)).astype(float))
                tol[:,flow]=np.maximum(tol[:,flow],cfg['legacy_float_flow_max_ulps']*ulp)
        ok=bool(np.all(diff<=tol))
        closure.append({'check':label,'passed':ok,'maxAbsDiff':float(diff.max()),
                        'maxScaledDiff':float((diff/np.maximum(1.,np.maximum(np.abs(a),np.abs(b)))).max()),
                        'regularMaxAbsDiff':float(diff[:,1:-1].max()),'flowMaxAbsDiff':float(diff[:,[0,-1]].max()),
                        'legacyFloatComparison':legacy_float,'rtol':rtol,'exactContentEqual':bool(np.array_equal(a[0],b[0]))})
        if not ok:
            print('CLOSURE_FAIL',closure[-1],flush=True)
            raise ValueError('Projection closure failed: '+label)
    def qa(species,stage,mass,var,group,window,source=None):
        source=hs if source is None else source
        result=[];axis=None
        for side,prefix in [('SE',''),('ME','Mixed_')]:
            summed=None
            for c in group['bins']:
                name='hKstar'+var+('PrePid' if stage=='prePID' else '')+'_'+prefix+species+suffixes[mass]+'_CentBin'+str(c)
                h=source[name];value=project(h,window)
                if axis is None:axis=h['xEdges']
                if h['xEdges']!=axis:raise ValueError('Axis mismatch')
                summed=value if summed is None else summed+value
            result.append(summed)
        return result[0],result[1],axis
    def baseline(species,mass,group,source=hs):
        result=[];axis=None
        for side in ['SE','ME']:
            name='hKstar'+side+'VsCent_lambda_'+baseline_species[species]+'_'+cfg['baseline_mass_suffixes'][mass]
            h=source[name];s=np.zeros((2,h['shape'][1]))
            for c in group['bins']:
                y=aligned(h['yEdges'],float(c)-0.5)+1
                s[0]+=matrix(h)[y];s[1]+=matrix(h,'error2')[y]
            axis=h['xEdges'];result.append(s)
        return result[0],result[1],axis
    def add(family,species,stage,mass,selection,group,s,m,edges,norm=None):
        meta=dict(family=family,species=species,stage=stage,mass=mass,selection=selection,
                  group=group['id'],label=group['label'])
        sid='|'.join(meta[k] for k in ['family','species','stage','mass','selection','group'])
        pts,met=analyze(s,m,edges,cfg,norm or cfg['normalization'])
        series.append(dict(meta,series_id=sid,points=pts,metrics=met))
        spectra[sid]=(s,m,edges)
    for species in cfg['species']:
        for group in cfg['groups']:
            for mass in suffixes:
                s,m,edges=qa(species,'selected',mass,'NSigma',group,None)
                rs,rm,re=qa(species,'selected',mass,'POverQ',group,None)
                check('NSigma_vs_R_'+species+'_'+mass+'_'+group['id']+'_SE',s,rs)
                check('NSigma_vs_R_'+species+'_'+mass+'_'+group['id']+'_ME',m,rm)
                bs,bm,be=baseline(species,mass,group)
                check('QA_vs_physical_'+species+'_'+mass+'_'+group['id']+'_SE',s,bs,legacy_float=True)
                check('QA_vs_physical_'+species+'_'+mass+'_'+group['id']+'_ME',m,bm,legacy_float=True)
                add('baseline',species,'selected',mass,'all',group,s,m,edges)
            for selection,window in cfg['mt_windows'][species].items():
                result=[]
                for prefix in ['', 'Mixed_']:
                    summed=None
                    for c in group['bins']:
                        h=hs['hKstarMt_'+prefix+species+'_CentBin'+str(c)]
                        value=project(h,window)
                        summed=value if summed is None else summed+value
                    result.append(summed)
                add('mt',species,'selected','signal',selection,group,result[0],result[1],h['xEdges'])
            for stage in cfg['stages']:
                for selection,window in cfg['nsigma_windows'].items():
                    s,m,edges=qa(species,stage,'signal','NSigma',group,window)
                    add('nsigma' if stage=='selected' else 'prepid',species,stage,'signal',selection,group,s,m,edges)
                for selection,window in cfg['rigidity_windows'].items():
                    s,m,edges=qa(species,stage,'signal','POverQ',group,window)
                    add('rigidity' if stage=='selected' else 'prepid_rigidity',species,stage,'signal',selection,group,s,m,edges)
            s,m,edges=qa(species,'selected','signal','NSigma',group,None)
            add('prepid',species,'selected','signal','all',group,s,m,edges)
            for selection,norm in cfg['normalization_windows'].items():
                add('normalization',species,'selected','signal',selection,group,s,m,edges,norm)
            add('production',species,'selected','signal','new1CE',group,s,m,edges)
            os,om,oe=qa(species,'selected','signal','NSigma',group,None,source=oldhs)
            add('production',species,'selected','signal','old04BF',group,os,om,oe)
            for stage in cfg['stages']:
                for var in ['NSigma','POverQ']:
                    all_s,all_m,edges=qa(species,stage,'signal',var,group,None)
                    ls,lm,_=qa(species,stage,'signalLow',var,group,None)
                    us,um,_=qa(species,stage,'signalHigh',var,group,None)
                    check('mass_halves_'+species+'_'+stage+'_'+var+'_'+group['id']+'_SE',all_s,ls+us)
                    check('mass_halves_'+species+'_'+stage+'_'+var+'_'+group['id']+'_ME',all_m,lm+um)
    # Save cut-to-inclusive differences with the overlap covariance of nested counts.
    comparisons=[];byid={s['series_id']:s for s in series}
    for row in series:
        if row['family'] not in ['nsigma','rigidity']:continue
        baseid='|'.join(['baseline',row['species'],'selected','signal','all',row['group']])
        s,m,edges=spectra[row['series_id']];bs,bm,_=spectra[baseid]
        norm=region(edges,cfg['normalization'])
        for name,bounds in cfg['high_windows'].items():
            target=region(edges,bounds)
            cv,ce,cg=ratio(s[0,1:-1],m[0,1:-1],s[1,1:-1],m[1,1:-1],norm,target)
            bv,be,bg=ratio(bs[0,1:-1],bm[0,1:-1],bs[1,1:-1],bm[1,1:-1],norm,target)
            if cv is None or bv is None:continue
            covariance=float((cg[0]*bg[0]*s[1,1:-1]+cg[1]*bg[1]*m[1,1:-1]).sum())
            ve=max(0.,ce*ce+be*be-2*covariance)
            comparisons.append({'series_id':row['series_id'],'window':name,'difference':cv-bv,
                                'error':math.sqrt(ve),'covariance':covariance})
    fields=['series_id','group','label','family','species','stage','mass','selection','low','high','value','error']
    def tab(v):return 'nan' if v is None else v
    with (out/'points.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
        for row in series:
            for pt in row['points']:w.writerow([tab(dict(row,**pt).get(k)) for k in fields])
    mfields=fields[:8]+['high07_08','high07_08_error','high08_10','high08_10_error','slope','slope_error','SEnorm','MEnorm']
    with (out/'metrics.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(mfields)
        for row in series:
            met=row['metrics'];vals=dict(row,**met)
            for name,x in met['windows'].items():vals[name]=x['value'];vals[name+'_error']=x['error']
            w.writerow([tab(vals.get(k)) for k in mfields])
    validation={'passed':True,'seriesCount':len(series),'pointCount':sum(len(x['points']) for x in series),
                'closureChecks':closure,'newRootUUID':doc['uuid'],'oldRootUUID':olddoc['uuid'],
                'sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.settings,Path(__file__),out/cfg['input_export'],out/cfg['old_export']]},
                'statisticalModel':'Stored sumw2 or count variances, independent native SE/ME cells; full derivative normalization covariance; nested cuts use subset overlap covariance. No event/shared-ME-pair covariance is stored.'}
    save(out/'results.json',{'settings':cfg,'series':series,'nestedCutComparisons':comparisons})
    save(out/'validation.json',validation)
    for group in ['cent0_20','cent20_80','cent10_20']:
        for row in series:
            if row['group']==group and row['species']=='4He' and row['family'] in ['baseline','nsigma','rigidity'] and row['mass']=='signal':
                print(row['family'],row['selection'],group,row['metrics']['windows']['high07_08'], 'slope',row['metrics']['slope'])
    print('SERIES',len(series),'CLOSURE',len(closure),'PASS')

if __name__=='__main__':main()
