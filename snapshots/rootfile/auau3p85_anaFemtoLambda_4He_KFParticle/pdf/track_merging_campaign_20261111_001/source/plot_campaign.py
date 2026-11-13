#!/usr/bin/env python3
"""Render one verified raw-SE angular-QA snapshot; no job dispatch or efficiency inference."""
import argparse
import fcntl
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
from datetime import datetime, timezone, timedelta
from pathlib import Path
import yaml


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()


def atomic(path, obj):
    temp=path.with_name(path.name+'.partial')
    temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    os.replace(str(temp),str(path))


def cp(text):return json.dumps(str(text),ensure_ascii=True)


def aligned(edges, value):
    found=[i for i,x in enumerate(edges) if abs(x-value)<1.e-9]
    if len(found)!=1:raise ValueError('Display/origin boundary is not a stored edge: '+str(value))
    return found[0]


def describe(h, widths):
    nx,ny=len(h['xEdges'])-1,len(h['yEdges'])-1
    if h['dimension']!=2 or h['shape']!=[ny+2,nx+2]:raise ValueError('Invalid histogram shape')
    a=h['content'];e=h['error2']
    if len(a)!=(nx+2)*(ny+2) or len(e)!=len(a):raise ValueError('Invalid histogram storage')
    if any(not math.isfinite(x) or x<0 for x in a+e):raise ValueError('Invalid count/variance')
    if abs(sum(a)-h['entries'])>1.e-6:raise ValueError('Histogram Entries/bin sum differs')
    regular=sum(a[y*(nx+2)+x] for y in range(1,ny+1) for x in range(1,nx+1));boxes={}
    for w in widths:
        xl,xh=aligned(h['xEdges'],-w),aligned(h['xEdges'],w)
        yl,yh=aligned(h['yEdges'],-w),aligned(h['yEdges'],w)
        boxes[str(w)]=sum(a[y*(nx+2)+x] for y in range(yl+1,yh+1) for x in range(xl+1,xh+1))
    return {'entries':h['entries'],'regular':regular,'flow':h['entries']-regular,'originBoxes':boxes,'nx':nx,'ny':ny}


def display_max(h, bins, bounds):
    nx,ny=len(h['xEdges'])-1,len(h['yEdges'])-1
    if nx%bins or ny%bins:raise ValueError('Display bins must divide both native axes')
    fx,fy=nx//bins,ny//bins;a=h['content'];maximum=0
    for j in range(bins):
        if h['yEdges'][j*fy]<bounds[0]-1.e-9 or h['yEdges'][(j+1)*fy]>bounds[1]+1.e-9:continue
        for i in range(bins):
            if h['xEdges'][i*fx]<bounds[0]-1.e-9 or h['xEdges'][(i+1)*fx]>bounds[1]+1.e-9:continue
            n=sum(a[y*(nx+2)+x] for y in range(j*fy+1,(j+1)*fy+1) for x in range(i*fx+1,(i+1)*fx+1))
            maximum=max(maximum,n)
    return maximum


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--campaign',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--settings',type=Path,default=Path(__file__).with_name('plot_settings.yaml'))
    ap.add_argument('--once',action='store_true',required=True,help='Render one snapshot; no daemon is started')
    a=ap.parse_args();campaign=a.campaign.resolve(strict=True);out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).resolve();owner=source.parents[1];project=source.parents[5]
    cfg=yaml.safe_load(a.settings.read_text())
    required={'fullRange','zoomRange','fullDisplayBins','zoomDisplayBins','originHalfWidths','canvasWidth','canvasHeight'}
    if set(cfg)!=required:raise ValueError('Plot settings keys differ from the required schema')
    latest=campaign/'aggregate/latest.json'
    if not latest.is_file():print('CAMPAIGN_PLOT_WAIT no verified aggregate yet');return 0
    with (out/'plot.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('CAMPAIGN_PLOT_WAIT another renderer is active');return 0
        raw=latest.read_bytes();meta=json.loads(raw)
        if meta.get('passed') is not True or meta.get('validation',{}).get('passed') is not True:raise ValueError('Aggregate has not passed validation')
        jobs=meta['sourceJobs']
        if meta['sourceCount']!=len(jobs) or not jobs or len(set(jobs))!=len(jobs) or any(not re.fullmatch(r'[0-9]+',str(x)) for x in jobs):raise ValueError('Invalid source-job identity/count')
        if len(meta['inputQaRoots'])!=len(jobs):raise ValueError('Source-QA root list length differs')
        if [Path(x).parent.name for x in meta['inputQaRoots']]!=[str(x) for x in jobs]:raise ValueError('Source-QA root order differs from sourceJobs')
        root=Path(meta['root']).resolve(strict=True);export=Path(meta['histogramExport']).resolve(strict=True)
        if root.parent!=campaign/'aggregate' or export.parent!=root.parent:raise ValueError('Aggregate paths are outside selected campaign')
        if sha(root)!=meta['rootSha256']:raise ValueError('Aggregate ROOT hash differs')
        root_stat=root.stat();finger=(root_stat.st_size,root_stat.st_mtime_ns,root_stat.st_ino)
        data=json.loads(export.read_text())
        if Path(data['input']).resolve()!=root:raise ValueError('Histogram export belongs to another ROOT')
        hist={}
        for path,h in data['histograms'].items():
            if h['name'] in hist:raise ValueError('Ambiguous histogram name')
            hist[h['name']]=h
        if data['histogramCount']!=len(hist) or len(hist)!=meta['validation']['histogramCount']:raise ValueError('Histogram export count mismatch')
        counters=meta['counters'];selected=int(counters.get('signal_pairs_selected_4He',0))
        if counters.get('cache_events')!=meta['sourceAcceptedEventsRepresented']:raise ValueError('Accepted-event counter differs')
        if not 0<=meta['sourceAcceptedEventsRepresented']<=meta['sourceInputEventsRepresented']:raise ValueError('Invalid source coverage')
        cmd_meta=meta['validation']['command'];wrapper=project/'script/singularity_root_batch.sh'
        if len(cmd_meta)<4 or cmd_meta[0]!='bash' or Path(cmd_meta[1]).resolve()!=wrapper or cmd_meta[3]!='root4star':raise ValueError('Unexpected validated ROOT runtime command')
        mainconf=Path(cmd_meta[2]).resolve(strict=True);runtime=mainconf.parents[2]
        if runtime.parent!=owner or not re.fullmatch(r'runtime(?:_v[0-9]+)?',runtime.name):raise ValueError('Mainconf is not a frozen campaign runtime')
        settings_hash=sha(a.settings);source_hash=sha(source);export_hash=sha(export);latest_hash=hashlib.sha256(raw).hexdigest()
        version=hashlib.sha256((source_hash+settings_hash).encode()).hexdigest()[:10]
        target=out/('jobs_%05d_%s_%s'%(len(jobs),meta['rootSha256'][:12],version))
        if target.exists():
            old=json.loads((target/'manifest.json').read_text())
            if old.get('status')!='passed' or old['sourceRootSha256']!=meta['rootSha256'] or old['sourceJobs']!=jobs:raise ValueError('Existing plot snapshot is inconsistent')
            for name,item in old['artifacts'].items():
                if sha(target/name)!=item['sha256']:raise ValueError('Existing plot artifact hash differs')
            if latest.read_bytes()==raw:atomic(out/'latest.json',dict(old,manifest=str(target/'manifest.json')))
            print('CAMPAIGN_PLOT_ALREADY_RENDERED',target);return 0
        staging=Path(tempfile.mkdtemp(prefix='.pending_',dir=str(out)))
        metric={};counter_templates={'GlobalMinusAsin':'global_formula_valid_','HelixFixedGlobalEta':'helix_fixed_valid_','HelixFixedSpatialEta':'helix_fixed_valid_','HelixClosestSpatialEta':'helix_common_valid_'}
        for daughter in ['proton','pion']:
            for method,prefix in counter_templates.items():
                key='hDphiDeta_selected_4He_'+daughter+'_'+method;expected=int(counters.get(prefix+'hDphiDeta_selected_4He_'+daughter,0));h=hist.get(key)
                if h is None:
                    if expected:raise ValueError('Missing nonempty method '+key)
                    continue
                metric[key]=describe(h,cfg['originHalfWidths'])
                if metric[key]['entries']!=expected:raise ValueError('Method count differs from validated counters')
        oldkey='hDphiDeta_proton_4He_EventField'
        if oldkey in hist:
            metric[oldkey]=describe(hist[oldkey],cfg['originHalfWidths'])
            if metric[oldkey]['entries']!=selected-int(counters.get('old_primary_gate_failed_4He',0)):raise ValueError('Old primary-gate count differs')
        method_labels={'GlobalMinusAsin':'Global momentum, corrected -asin','HelixFixedGlobalEta':'Helix at fixed QA radius; global eta','HelixFixedSpatialEta':'Helix at fixed QA radius; spatial eta','HelixClosestSpatialEta':'Closest common radius; spatial eta'}
        corrected=lambda d:['hDphiDeta_selected_4He_'+d+'_'+m for m in method_labels]
        pages=[('proton_global','proton',cfg['fullRange'],cfg['fullDisplayBins'],[oldkey]+corrected('proton')[:2]+[corrected('proton')[-1]]),('proton_zoom','proton',cfg['zoomRange'],cfg['zoomDisplayBins'],corrected('proton')),('pion_global','pion',cfg['fullRange'],cfg['fullDisplayBins'],corrected('pion')),('pion_zoom','pion',cfg['zoomRange'],cfg['zoomDisplayBins'],corrected('pion'))]
        pdf=staging/'merging_maps.pdf';commands=['gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPalette(1);', 'TFile* pf=TFile::Open('+cp(root)+',"READ");','if(!pf || pf->IsZombie() || pf->TestBit(TFile::kRecovered))gSystem->Exit(2);','std::cout << "CAMPAIGN_PLOT_UUID " << pf->GetUUID().AsString() << std::endl;', 'TCanvas cover("plot_cover","Raw SE QA",800,600);','cover.Print('+cp(str(pdf)+'[')+');']
        for page,(name,daughter,bounds,bins,keys) in enumerate(pages):
            lo,hi=bounds;maximum=max([1]+[display_max(hist[k],bins,bounds) for k in keys if k in hist])
            canvas='p%d'%page;commands+=['TCanvas '+canvas+'('+cp(canvas)+',"Raw SE angular QA",%d,%d);'%(cfg['canvasWidth'],cfg['canvasHeight']),canvas+'.Divide(2,2,.003,.003);']
            for panel,key in enumerate(keys):
                pad='%s.cd(%d);'%(canvas,panel+1)
                commands+=[pad,'gPad->SetPad(gPad->GetXlowNDC(),.12+.78*gPad->GetYlowNDC(),gPad->GetXlowNDC()+gPad->GetWNDC(),.12+.78*(gPad->GetYlowNDC()+gPad->GetHNDC()));','gPad->SetLeftMargin(.14);gPad->SetRightMargin(.17);gPad->SetTopMargin(.23);gPad->SetBottomMargin(.15);']
                t='t%d_%d'%(page,panel);commands+=['TLatex* '+t+'=new TLatex();'+t+'->SetNDC();'+t+'->SetTextSize(.032);']
                if key not in hist:commands+=[t+'->DrawLatex(.14,.60,"No allocated map; zero admissible pairs");'];continue
                h=hist[key];m=metric[key];var='h%d_%d'%(page,panel);fx=(len(h['xEdges'])-1)//bins;fy=(len(h['yEdges'])-1)//bins
                method='OLD EventField: +asin / primary gate' if key==oldkey else method_labels[key.rsplit('_',1)[-1]]
                eta='#Delta#eta_{pMom}' if key==oldkey else ('#Delta#eta_{spatial}' if 'SpatialEta' in key else '#Delta#eta_{gMom}')
                commands+=['TH2* '+var+'=(TH2*)((TH2*)pf->Get('+cp(key)+'))->Clone('+cp(var)+');',var+'->SetDirectory(0);',var+'->Rebin2D(%d,%d);'%(fx,fy),var+'->GetXaxis()->SetRangeUser(%g,%g);'%tuple(bounds),var+'->GetYaxis()->SetRangeUser(%g,%g);'%tuple(bounds),var+'->GetXaxis()->SetTitle("#Delta#varphi* (rad)");',var+'->GetYaxis()->SetTitle('+cp(eta)+');',var+'->GetZaxis()->SetTitle("Raw SE / display bin");',var+'->GetYaxis()->SetTitleOffset(1.3);',var+'->GetXaxis()->SetTitleSize(.047);'+var+'->GetYaxis()->SetTitleSize(.047);',var+'->GetXaxis()->SetLabelSize(.04);'+var+'->GetYaxis()->SetLabelSize(.04);',var+'->SetMinimum(0);'+var+'->SetMaximum(%g);'%maximum,var+'->Draw("COLZ");',t+'->DrawLatex(.14,.95,'+cp(method)+');',t+'->SetTextSize(.029);'+t+'->DrawLatex(.14,.88,'+cp('Entries %d; native range %d; origin +/-%.4g: %d'%(m['entries'],m['regular'],cfg['originHalfWidths'][0],m['originBoxes'][str(cfg['originHalfWidths'][0])]))+');']
                z='z%d_%d'%(page,panel);commands+=['TLine* '+z+'a=new TLine(%g,0,%g,0);'%(lo,hi),'TLine* '+z+'b=new TLine(0,%g,0,%g);'%(lo,hi),z+'a->SetLineStyle(3);'+z+'b->SetLineStyle(3);'+z+'a->SetLineColor(kGray+2);'+z+'b->SetLineColor(kGray+2);'+z+'a->Draw();'+z+'b->Draw();']
            commands+=[canvas+'.cd();','TLatex* top%d=new TLatex();top%d->SetNDC();top%d->SetTextSize(.026);'%(page,page,page),'top%d->DrawLatex(.025,.98,'%page+cp('Lambda daughter '+daughter+' vs selected 4He; '+('origin zoom' if 'zoom' in name else 'angular QA comparison'))+');','top%d->SetTextSize(.019);top%d->DrawLatex(.025,.945,'%(page,page)+cp('Source jobs %d; source Pico entries %d; selected He4 signal pairs %d; raw SE, all stored k*'%(len(jobs),meta['sourceInputEventsRepresented'],selected))+');','top%d->SetTextSize(.018);top%d->DrawLatex(.025,.075,'%(page,page)+cp('Origin peak not assessed. No ME reference or pair-efficiency inference. Display binning only; no smoothing.')+');','top%d->DrawLatex(.025,.043,'%page+cp('Admission differs by method. Closest radius is selected per pair by minimum 3D separation; spatial eta is a different ordinate.')+');','top%d->SetTextSize(.017);top%d->DrawLatex(.025,.015,'%(page,page)+cp('Accepted cache frames %d; QA Pico events read %d. Source coverage differs from the QA read count.'%(meta['sourceAcceptedEventsRepresented'],counters['pico_events_read']))+');',canvas+'.Print('+cp(pdf)+');',canvas+'.Print('+cp(staging/(name+'.png'))+');']
        commands+=['cover.Print('+cp(str(pdf)+']')+');','std::cout << "CAMPAIGN_PLOT_ROOT_PASS pages=4" << std::endl;','.q']
        text='\n'.join(commands)+'\n';(staging/'root_commands.txt').write_text(text)
        command=['bash',str(wrapper),str(mainconf),'root4star','-l','-b']
        done=subprocess.run(command,input=text,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (staging/'root.log').write_text(done.stdout)
        if done.returncode or 'CAMPAIGN_PLOT_ROOT_PASS pages=4' not in done.stdout or 'Error:' in done.stdout:raise RuntimeError('ROOT drawing failed; preserved '+str(staging/'root.log'))
        match=re.search(r'CAMPAIGN_PLOT_UUID\s+([A-Za-z0-9-]+)',done.stdout)
        if not match or match.group(1)!=data['uuid']:raise ValueError('ROOT/export UUID mismatch')
        if not pdf.read_bytes().rstrip().endswith(b'%%EOF') or len(re.findall(rb'/Type\s*/Page\b',pdf.read_bytes()))!=4:raise ValueError('Invalid four-page PDF')
        artifacts={}
        for path in [pdf]+[staging/(x[0]+'.png') for x in pages]:
            if path.suffix=='.png' and not path.read_bytes().startswith(b'\x89PNG\r\n\x1a\n'):raise ValueError('Invalid PNG')
            artifacts[path.name]={'bytes':path.stat().st_size,'sha256':sha(path)}
        if sha(root)!=meta['rootSha256'] or (root.stat().st_size,root.stat().st_mtime_ns,root.stat().st_ino)!=finger or sha(export)!=export_hash:raise ValueError('Read-only input changed during plotting')
        result={'status':'passed','createdAtUtc':datetime.now(timezone.utc).isoformat(),'sourceManifest':str(latest),'sourceManifestSha256':latest_hash,'sourceRoot':str(root),'sourceRootSha256':meta['rootSha256'],'histogramExport':str(export),'histogramExportSha256':export_hash,'sourceJobs':jobs,'sourceCount':len(jobs),'selectedHe4Pairs':selected,'sourceInputEventsRepresented':meta['sourceInputEventsRepresented'],'sourceAcceptedEventsRepresented':meta['sourceAcceptedEventsRepresented'],'pages':4,'artifacts':artifacts,'histogramMetrics':metric,'inputUnchanged':True,'sourceSha256':source_hash,'plotSettingsSha256':settings_hash,'rootCommand':command,'rootStdinSha256':hashlib.sha256(text.encode()).hexdigest(),'originPeakAssessment':'not assessed','efficiencyAssessment':'not assessed; raw SE without ME','methodCaveat':'Closest radius is pair-conditioned; spatial eta differs from track eta; old QA has a primary-fit gate.'}
        atomic(staging/'manifest.json',result);os.rename(str(staging),str(target))
        if latest.read_bytes()==raw:atomic(out/'latest.json',dict(result,manifest=str(target/'manifest.json')))
        print('CAMPAIGN_PLOT_PASS sources',len(jobs),'He4',selected,target)
        return 0


if __name__=='__main__':raise SystemExit(main())
