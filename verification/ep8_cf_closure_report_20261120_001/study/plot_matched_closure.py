#!/usr/bin/env python3
"""Plot saved matched histograms only; no event processing."""
import pathlib,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=pathlib.Path(__file__).resolve().parent
F=R.parents[1]/'share/figure/20261120_lambda_ep8_cf_closure_report'
F.mkdir(parents=True,exist_ok=True)
v=json.loads((R/'closure_root_audit.json').read_text());assert v['passed']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
artifacts=[];cfdata={}
def get(sp,mode,method='old'):
 return v['histograms']['hKstar'+mode+'_lambda_'+sp+'_signal'][method]
def reb(h):
 y=np.array(h['contentIncludingFlow'][1:-1]).reshape(-1,10).sum(1)
 e=np.sqrt((np.array(h['errorIncludingFlow'][1:-1])**2).reshape(-1,10).sum(1))
 edges=np.array(h['edges'])[::10]
 return edges,(edges[1:]+edges[:-1])/2,y,e
def norm(h):
 e=np.array(h['edges']);y=np.array(h['contentIncludingFlow'][1:-1]);return float(y[(e[:-1]>=.4-1e-12)&(e[1:]<=.5+1e-12)].sum())
def save(fig,name):
 for ext in ['pdf','png']:
  p=F/(name+'.'+ext);assert not p.exists();fig.savefig(p,dpi=180,bbox_inches='tight');artifacts.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 plt.close(fig)
for sp,name in [('he4','closure_matched_seme'),('he3','closure_matched_seme_he3_control')]:
 label=r'$\Lambda$–$^4$He' if sp=='he4' else r'$\Lambda$–$^3$He'
 fig,ax=plt.subplots(2,2,figsize=(12.5,7),sharex='col',gridspec_kw={'height_ratios':[3,1]})
 for col,mode in enumerate(['SE','ME']):
  h=get(sp,mode);n=get(sp,mode,'treeReplay');e,x,y,err=reb(h);_,_,ny,nerr=reb(n);assert np.array_equal(y,ny) and np.array_equal(err,nerr)
  ax[0,col].stairs(y,e,color='black',label='Original direct-Pico histogram analysis');ax[0,col].errorbar(x,y,yerr=err,fmt='none',color='black',lw=.8)
  ax[0,col].plot(x,ny,'o',mfc='none',mec='crimson',ms=5,label='Compact tree writer → replay');ax[0,col].axvspan(.4,.5,color='gold',alpha=.15)
  ax[0,col].set_title(f'{label} signal {mode}: {int(h["entries"])} entries including flow\nVisible: {int(y.sum())}; UF+OF: {int(h["entries"]-y.sum())}; H integral: {int(norm(h))}');ax[0,col].set_ylabel('Raw pairs / 0.05 GeV/c');ax[0,col].grid(alpha=.18)
  ax[1,col].axhline(0,color='black',lw=.8);ax[1,col].plot(x,ny-y,'o',color='crimson',ms=3);ax[1,col].set_ylim(-.5,.5);ax[1,col].set_yticks([-.5,0,.5]);ax[1,col].set_ylabel('Tree − direct\ncounts');ax[1,col].set_xlabel(r'$k^*$ [GeV/c]');ax[1,col].set_xlim(0,1);ax[1,col].grid(alpha=.18)
 ax[0,0].legend(loc='upper left');fig.suptitle('Matched-input closure: original direct PicoDst vs compact tree replay',fontsize=15)
 fig.text(.5,-.02,'Same PicoDst prefix: 10,000 input → 9,761 reconstructed events; replay reads all 9,761 saved events.\nSame cent9 × vertex-z10 × EP8 pools, depth20, bidirectional bufferAll, segment reset; no close-pair veto.\n108 producer QA + 1,186 pair histograms: all 46,336,786 cells including flow match in content/error/entries/axes/raw Sumw2.\nH = 0.4–0.5 GeV/c (yellow). 200 raw bins are summed to 0.05 GeV/c for this display.',ha='center',va='top',fontsize=9)
 fig.tight_layout(rect=[0,.06,1,.96]);save(fig,name)
fig,ax=plt.subplots(2,2,figsize=(12.5,7),sharex='col',gridspec_kw={'height_ratios':[3,1]})
for col,sp in enumerate(['he4','he3']):
 se=get(sp,'SE');me=get(sp,'ME');sH=norm(se);mH=norm(me);assert sH>0 and mH>0
 e,x,s,es=reb(se);_,_,m,em=reb(me);mask=m>0;scale=mH/sH;c=scale*s[mask]/m[mask];err=scale*np.sqrt((es[mask]/m[mask])**2+(s[mask]*em[mask]/m[mask]**2)**2)
 ax[0,col].axhline(1,color='steelblue',ls=':');ax[0,col].errorbar(x[mask],c,yerr=err,fmt='.',color='black',capsize=2,label='Original direct Pico');ax[0,col].plot(x[mask],c,'o',mfc='none',mec='crimson',ms=6,label='Compact tree replay');ax[0,col].axvspan(.4,.5,color='gold',alpha=.15)
 ax[0,col].set_title(f'{sp} signal: C = (ME_H / SE_H) × SE / ME\nH: SE={int(sH)}, ME={int(mH)}, normalization factor={scale:.4g}');ax[0,col].set_ylabel('Raw signal CF');ax[0,col].grid(alpha=.18);ax[0,col].legend()
 ax[1,col].axhline(0,color='black',lw=.8);ax[1,col].plot(x[mask],np.zeros(mask.sum()),'o',color='crimson',ms=3);ax[1,col].set_ylim(-.5,.5);ax[1,col].set_ylabel('Tree − direct\nCF');ax[1,col].set_xlabel(r'$k^*$ [GeV/c]');ax[1,col].set_xlim(0,1);ax[1,col].grid(alpha=.18)
 cfdata[sp]=dict(normSE=sH,normME=mH,normFactor=scale,edges=e.tolist(),validMixedBins=mask.tolist(),cf=c.tolist(),fixedNormalizationError=err.tolist(),oldAndReplayIdentical=True)
fig.suptitle('Small matched sample: exact CF reproduction; no high-k* physics conclusion',fontsize=14)
fig.text(.5,-.02,'Signal window, no background subtraction or purity correction. H = 0.4–0.5 GeV/c; 0.05 GeV/c display bins.\nErrors propagate bin SE/ME errors with normalization fixed. Normalization uncertainty/covariance omitted; He4 H has one SE pair.\nME=0 bins are omitted. Old/replay use the same events; zero difference tests equality of dependent results.',ha='center',va='top',fontsize=9)
fig.tight_layout(rect=[0,.06,1,.96]);save(fig,'closure_matched_cf')
with (R/'closure_plot_artifacts.json').open('x') as out:json.dump(dict(passed=True,sourceAuditSha256=hashlib.sha256((R/'closure_root_audit.json').read_bytes()).hexdigest(),artifacts=artifacts,cf=cfdata,rebinFactor=10,normalizationGeVc=[.4,.5],note='Matched depth20/reset sample only; full-stat periodMix100 depth100/continuation is intentionally different.'),out,indent=2);out.write('\n')
print(json.dumps(dict(passed=True,artifacts=[x['path'] for x in artifacts]),indent=2))
