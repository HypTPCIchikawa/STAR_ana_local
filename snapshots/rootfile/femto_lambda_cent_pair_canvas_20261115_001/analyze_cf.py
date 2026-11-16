#!/usr/bin/env python3
"""Analyze projected raw pairs with the existing normalization Jacobian."""
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RATIO_SOURCE = ROOT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/source/analyze_highk.py'
spec = importlib.util.spec_from_file_location('existing_highk', RATIO_SOURCE)
existing = importlib.util.module_from_spec(spec)
spec.loader.exec_module(existing)


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def window(edges, bounds):
    return existing.region(edges, bounds)


def load_raw():
    groups = {}
    with (HERE / 'raw_projections.tsv').open() as f:
        for row in csv.reader(f, delimiter='\t'):
            if not row or row[0] != 'BIN':
                continue
            if len(row) != 10:
                raise ValueError('Wrong BIN format')
            ds, selection, cent, side = row[1:5]
            key = (ds, cent, selection)
            by_side = groups.setdefault(key, {}).setdefault(side, {})
            b = int(row[5])
            if b in by_side:
                raise ValueError('Duplicate projected bin')
            by_side[b] = tuple(map(float, row[6:]))
    raw = {}
    for key, sides in groups.items():
        if set(sides) != {'SE', 'ME'} or set(sides['SE']) != set(sides['ME']):
            raise ValueError('Incomplete SE/ME')
        n = len(sides['SE']) - 2
        if sorted(sides['SE']) != list(range(n+2)):
            raise ValueError('Missing raw/flow bin')
        edges = [sides['SE'][i][0] for i in range(1,n+1)] + [sides['SE'][n][1]]
        arrays = {}
        for side in ('SE', 'ME'):
            rows = [sides[side][i] for i in range(n+2)]
            if any(rows[i][:2] != sides['SE'][i][:2] for i in range(n+2)):
                raise ValueError('Projected axes differ')
            arrays[side] = np.array([[x[2] for x in rows], [x[3] for x in rows]])
            if np.any(~np.isfinite(arrays[side])) or np.any(arrays[side] < 0):
                raise ValueError('Invalid raw count/variance')
        raw[key] = dict(edges=np.array(edges), **arrays)
    return raw


def closure(raw, cfg):
    checks = []
    def compare(ds, cent, side, label, a, b, partition=False):
        delta = a-b
        tol = cfg['closure']['regular_absolute_tolerance'] + cfg['closure']['regular_relative_tolerance'] * np.maximum(abs(a),abs(b))
        ok = bool(np.all(delta[:,1:-1] >= -tol[:,1:-1])) if partition else bool(np.all(abs(delta[:,1:-1]) <= tol[:,1:-1]))
        checks.append(dict(dataset=ds, centrality=cent, side=side, check=label,
                           passed_regular=ok, relation='nonnegative remainder' if partition else 'equal',
                           max_regular_count_difference=float(abs(delta[0,1:-1]).max()),
                           regular_count_remainder=float(delta[0,1:-1].sum()),
                           kstar_flow_count_difference=[float(delta[0,0]),float(delta[0,-1])]))
        if not ok:
            raise ValueError('Regular-bin projection closure failed: ' + str(checks[-1]))
    for ds in cfg['inputs']:
        for c in cfg['centralities']:
            cent = c['id']
            for side in ('SE','ME'):
                get = lambda name: raw[(ds,cent,name)][side]
                bins = [get('pair_bin%02d' % i) for i in range(1, cfg['pair_rapidity']['n_bins']+1)]
                compare(ds,cent,side,'pairY all including UF/OF vs inclusive',get('pair_all'),get('inclusive'))
                compare(ds,cent,side,'ten regular pairY bins vs regular projection',sum(bins),get('pair_regular'))
                compare(ds,cent,side,'pairY regular + UF + OF vs all',get('pair_regular')+get('pair_under')+get('pair_over'),get('pair_all'))
                compare(ds,cent,side,'inclusive minus strict He4 Below/Above',get('inclusive'),get('he4_below')+get('he4_above'),True)
                for i,(lo,hi) in enumerate(cfg['pair_rapidity']['adjacent_merge_groups'],1):
                    compare(ds,cent,side,'adjacent merge %d' % i,get('pair_merge%02d'%i),sum(bins[lo-1:hi]))
    return checks


def gls(points, gradients, s, m, vs, vm, edges, bounds):
    ids = [i for i,p in enumerate(points) if p['value'] is not None and
           p['low'] >= bounds[0]-1e-9 and p['high'] <= bounds[1]+1e-9]
    if len(ids) < 3:
        return dict(range=bounds, bins=len(ids), slope=None, slope_error=None)
    js = np.array([gradients[i][0] for i in ids]);jm = np.array([gradients[i][1] for i in ids])
    cov = (js*vs).dot(js.T)+(jm*vm).dot(jm.T)
    inv = np.linalg.pinv(cov,rcond=1e-12)
    rank = int(np.linalg.matrix_rank(cov))
    x = np.array([(points[i]['low']+points[i]['high'])/2 for i in ids])
    y = np.array([points[i]['value'] for i in ids]);at=float(x.mean())
    design = np.column_stack((np.ones(len(x)),x-at))
    normal = design.T.dot(inv).dot(design)
    if np.linalg.matrix_rank(normal) < 2:
        return dict(range=bounds,bins=len(ids),slope=None,slope_error=None)
    pcov=np.linalg.inv(normal);pars=pcov.dot(design.T).dot(inv).dot(y)
    residual=y-design.dot(pars);one=np.ones(len(x));den=float(one.dot(inv).dot(one))
    level=float(one.dot(inv).dot(y)/den);rlevel=y-level
    return dict(range=bounds,bins=len(ids),rank=rank,slope=float(pars[1]),
                slope_error=math.sqrt(max(0.,float(pcov[1,1]))),intercept=float(pars[0]),intercept_at=at,
                linear_chi2=float(residual.dot(inv).dot(residual)),linear_ndf=rank-2,
                constant=level,constant_error=math.sqrt(1/den),
                constant_chi2=float(rlevel.dot(inv).dot(rlevel)),constant_ndf=rank-1,
                unity_chi2=float((y-1).dot(inv).dot(y-1)),unity_ndf=rank)


def analyze_one(h,cfg):
    factor=cfg['rebin_factor'];n=len(h['edges'])-1
    if n%factor:
        raise ValueError('Rebin would discard part of the stored range')
    rb=lambda x:x[1:-1].reshape(-1,factor).sum(axis=1)
    s,vs=map(rb,h['SE']);m,vm=map(rb,h['ME']);edges=h['edges'][::factor]
    norm=window(edges,cfg['normalization']);sn=float(s[norm].sum());mn=float(m[norm].sum())
    points=[];grads=[]
    for i in range(len(s)):
        c,error,grad=existing.ratio(s,m,vs,vm,norm,slice(i,i+1))
        points.append(dict(low=float(edges[i]),high=float(edges[i+1]),value=c,error=error))
        grads.append(grad)
    windows={}
    for key,bounds in cfg['high_windows'].items():
        target=window(edges,bounds);c,error,_=existing.ratio(s,m,vs,vm,norm,target)
        windows[key]=dict(range=bounds,value=c,error=error,SE=float(s[target].sum()),ME=float(m[target].sum()))
    middle=window(edges,cfg['high_windows']['k050_070']);tail=window(edges,cfg['high_windows']['k080_100'])
    tilt,tilt_err,_=existing.ratio(s,m,vs,vm,middle,tail)
    valid=[i for i,g in enumerate(grads) if g is not None]
    if valid:
        js=np.array([grads[i][0] for i in valid]);jm=np.array([grads[i][1] for i in valid])
        cov=(js*vs).dot(js.T)+(jm*vm).dot(jm.T)
    else:cov=np.empty((0,0))
    result=dict(edges=edges.tolist(),SE=s.tolist(),ME=m.tolist(),SEvariance=vs.tolist(),MEvariance=vm.tolist(),
                points=points,valid_covariance_bins=valid,normSE=sn,normME=mn,
                MEscaleToSE=(sn/mn if sn>0 and mn>0 else None),
                totalSE=float(s.sum()),totalME=float(m.sum()),
                kstarFlowSE=[float(h['SE'][0,0]),float(h['SE'][0,-1])],
                kstarFlowME=[float(h['ME'][0,0]),float(h['ME'][0,-1])],
                windows=windows,shape_ratio=dict(value=tilt,error=tilt_err,
                  definition='(S/M)[0.8,1.0) divided by (S/M)[0.5,0.7); normalization cancels'),
                fits={name:gls(points,grads,s,m,vs,vm,edges,bounds) for name,bounds in cfg['fit_ranges'].items()})
    return result,cov


def self_check():
    # An independent four-window count identity tests normalization uncertainty.
    s=np.array([40.,60.,80.,120.]);m=np.array([80.,120.,160.,240.])
    c,err,g=existing.ratio(s,m,s,m,slice(0,2),slice(2,4))
    expected=math.sqrt(1/200.+1/400.+1/100.+1/200.)
    assert abs(c-1)<1e-12 and abs(err-expected)<1e-12
    c2,e2,_=existing.ratio(10*s,10*m,10*s,10*m,slice(0,2),slice(2,4))
    assert c2==c and abs(e2-err/math.sqrt(10))<1e-12
    # Shared normalization creates a positive off-diagonal contribution.
    grads=[existing.ratio(s,m,s,m,slice(0,2),slice(i,i+1))[2] for i in [2,3]]
    covariance=float((grads[0][0]*grads[1][0]*s+grads[0][1]*grads[1][1]*m).sum())
    assert abs(covariance-(1/100.+1/200.))<1e-12


def main():
    self_check();cfg=yaml.safe_load((HERE/'study.yaml').read_text());manifest=json.loads((HERE/'input_manifest.json').read_text())
    raw=load_raw();assert len(raw)==manifest['planRows']
    checks=closure(raw,cfg);save(HERE/'closure_checks.json',checks)
    curves={};covariances={}
    with (HERE/'curve_bins.tsv').open('w') as bins_out, (HERE/'curve_metadata.tsv').open('w') as meta_out:
        for (ds,cent,selection),h in sorted(raw.items()):
            result,cov=analyze_one(h,cfg);key='__'.join((ds,cent,selection))
            result.update(dataset=ds,centrality=cent,selection=selection)
            curves[key]=result;covariances[key]=cov
            meta_out.write('\t'.join([ds,cent,selection,str(result['MEscaleToSE']) if result['MEscaleToSE'] is not None else 'NA'])+'\n')
            for i,p in enumerate(result['points']):
                f=[ds,cent,selection,str(i),repr(p['low']),repr(p['high']),repr(result['SE'][i]),repr(result['SEvariance'][i]),
                   repr(result['ME'][i]),repr(result['MEvariance'][i]),'NA' if p['value'] is None else repr(p['value']),
                   'NA' if p['error'] is None else repr(p['error'])]
                bins_out.write('\t'.join(f)+'\n')
    np.savez_compressed(str(HERE/'cf_covariances.npz'),**covariances)
    save(HERE/'analysis_results.json',dict(config=cfg,curves=curves,statistical_assumption='Independent histogram counts with common normalization covariance; event/pair reuse correlations unavailable',
        ratio_source=str(RATIO_SOURCE),ratio_source_sha256=hashlib.sha256(RATIO_SOURCE.read_bytes()).hexdigest()))
    fields=['dataset','centrality','selection','SEtotal','MEtotal','SEnorm','MEnorm','CFtail','CFtailError','shapeRatio','shapeRatioError','slope','slopeError','unityChi2','unityNdf']
    with (HERE/'metrics.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in curves.values():
            t=r['windows']['k080_100'];fit=r['fits']['high'];shape=r['shape_ratio']
            writer.writerow(dict(dataset=r['dataset'],centrality=r['centrality'],selection=r['selection'],SEtotal=r['totalSE'],MEtotal=r['totalME'],SEnorm=r['normSE'],MEnorm=r['normME'],
                CFtail=t['value'],CFtailError=t['error'],shapeRatio=shape['value'],shapeRatioError=shape['error'],slope=fit.get('slope'),slopeError=fit.get('slope_error'),unityChi2=fit.get('unity_chi2'),unityNdf=fit.get('unity_ndf')))
    for ds in cfg['inputs']:
        r=curves[ds+'__c0_20__inclusive'];print(ds,'inclusive tail',r['windows']['k080_100'],'shape',r['shape_ratio'],'slope',r['fits']['high'].get('slope'),flush=True)
    print('PASS numerical self-check,',len(checks),'closure checks,',len(curves),'CF series')


if __name__=='__main__':main()
