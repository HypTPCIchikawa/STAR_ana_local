#!/usr/bin/env python3
"""Project saved single-He4 lab-rapidity histograms into four/eight raw CFs."""
import argparse,hashlib,json,math,subprocess,sys,types
from pathlib import Path
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]

def need(value,message):
    if not value:raise RuntimeError(message)
def load(path):
    with Path(path).open() as f:return json.load(f)
def save(path,data):
    with Path(path).open('x') as f:json.dump(data,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
def resolve(value,base):
    p=Path(value);return (p if p.is_absolute() else base/p).resolve()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def stat(path):
    p=Path(path);s=p.stat();need(p.is_file() and s.st_size>0,'missing/empty '+str(p));return dict(path=str(p),bytes=s.st_size,device=s.st_dev,inode=s.st_ino,mtimeNs=s.st_mtime_ns)
def same(record):
    actual=stat(record['path']);need(all(actual[k]==record[k] for k in ['bytes','device','inode','mtimeNs']),'input stat changed '+record['path'])
def rows(path):
    out=[]
    for raw in Path(path).read_text().splitlines():
        if not raw.strip() or raw.lstrip().startswith('#'):continue
        need(raw==raw.strip() and len(raw.split())==1,'one path per list row required');out.append(resolve(raw,Path(path).parent))
    need(out and len(out)==len(set(out)),'empty/duplicate list');return out
def module(path,name):
    m=types.ModuleType(name);m.__file__=str(path);exec(compile(path.read_text(),str(path),'exec'),m.__dict__);return m

def validate_settings(path):
    settings=load(path);need(settings['schemaVersion']==1,'settings schema')
    n=settings['rapidityBinning']['nBins'];low=settings['rapidityBinning']['min'];high=settings['rapidityBinning']['max']
    need(n==8 and low<high and settings['coarseBins']==4,'explicit 8/4 binning required')
    need(settings['rapidityBinning']['frame']=='lab','lab frame required')
    need(settings['normalizationInterval']=='half-open [low, high)' and settings['displayRebin']>0,'invalid normalization/display contract')
    need(settings['requiredEPBins']==8 and settings['requiredPairCut'] is False,'reference EP8/cut-Off required')
    for interval in [settings['normalizationGeVc'],settings['displayKGeVc'],settings['displayCFRange'],settings['highKPlotGeVc'],settings['highKCFRange']]+list(settings['highKWindowsGeVc'].values()):
        need(len(interval)==2 and all(math.isfinite(v) for v in interval) and interval[0]<interval[1],'invalid interval')
    need(settings['requiredValidCentralities'] and set(settings['requiredValidCentralities'])<={c['id'] for c in settings['centralities']},'missing main centrality')
    for c in settings['centralities']:need(0<=c['cent9'][0]<=c['cent9'][1]<=8,'invalid cent9 range')
    return settings

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--settings',type=Path,required=True);ap.add_argument('--outputs',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--pilot',action='store_true');ap.add_argument('--pilot-events',type=int);ap.add_argument('--pilot-source',type=Path);ap.add_argument('--pilot-reference',type=Path);args=ap.parse_args()
    sp=args.settings.resolve();s=validate_settings(sp);manifest_path=resolve(s['inputManifest'],sp.parent);manifest=load(manifest_path);sources=manifest['sourceFiles'];need(len(sources)==20,'twenty source periods required')
    paths=rows(args.outputs.resolve());references=rows(resolve(s['referenceOutputsList'],sp.parent));need(len(references)==20,'twenty baseline reference outputs required')
    source_list=resolve(s['sourcePeriodList'],sp.parent);need([str(p) for p in rows(source_list)]==[str(Path(x['path']).resolve()) for x in sources],'source period order differs')
    plan=[]
    if args.pilot:
        need(len(paths)==1 and args.pilot_events and args.pilot_events>0 and args.pilot_source and args.pilot_reference,'pilot requires one output and explicit events/source/reference')
        found=[(i,x) for i,x in enumerate(sources) if Path(x['path']).resolve()==args.pilot_source.resolve()];need(len(found)==1,'pilot source not in period manifest');i,src=found[0];need(args.pilot_events<=src['acceptedEvents'],'pilot events exceed period');plan=[(i,args.pilot_events,src,paths[0],args.pilot_reference.resolve())]
    else:
        need(len(paths)==20,'production requires exactly twenty outputs');plan=[(i,x['acceptedEvents'],x,p,references[i]) for i,(x,p) in enumerate(zip(sources,paths))]
    out=args.output_dir.resolve();need(not out.exists(),'refusing existing output directory');out.mkdir(parents=True)
    pure_path=resolve(s['sharedCFMath'],sp.parent);pure=module(pure_path,'he4_y_cf_math')
    pin_paths=[sp,manifest_path,source_list,args.outputs.resolve(),resolve(s['referenceOutputsList'],sp.parent),pure_path,Path(__file__).resolve(),HERE/'ExtractHe4Y.C',HERE/'DrawHe4Y.C',resolve(s['sharedExtractorUtilities'],sp.parent),resolve(s['mainconfForRoot'],PROJECT)]
    pins=[dict(path=str(p),sha256=sha(p)) for p in pin_paths];source_stats=[stat(x['path']) for x in sources];input_stats=[stat(p) for _,_,_,p,ref in plan]+[stat(ref) for _,_,_,p,ref in plan]
    for original in sources:same(original)
    save(out/'input_provenance.json',dict(pilot=args.pilot,pins=pins,sourceStats=source_stats,inputStats=input_stats,treeEntriesRead=0))
    try:
        with (out/'extraction_plan.tsv').open('x') as f:
            for i,events,src,p,ref in plan:f.write('\t'.join(map(str,[i,events,Path(src['path']).name,p,ref]))+'\n')
        b=s['rapidityBinning'];req=dict(nRapidityBins=b['nBins'],rapidityMin=b['min'],rapidityMax=b['max'],absoluteTolerance=s['projectionAbsoluteTolerance'],relativeTolerance=s['projectionRelativeTolerance'],metadataBins=s['metadata']['nBins'],metadataMin=s['metadata']['min'],metadataMax=s['metadata']['max'],metadataDefinition=s['metadata']['definition'])
        req.update({'mix.'+k:v for k,v in s['requiredMixing'].items()});req.update({'reduced.'+k:v for k,v in s['requiredReducedSettings'].items()})
        req.update({'definition.'+str(i):value for i,value in enumerate(s['requiredDefinitionTokens'])})
        with (out/'requirements.tsv').open('x') as f:
            for k,v in sorted(req.items()):f.write(str(k)+'\t'+str(v)+'\n')
        mainconf=resolve(s['mainconfForRoot'],PROJECT);pure.root_run(mainconf,HERE/'ExtractHe4Y.C',[out/'extraction_plan.tsv',out/'requirements.tsv',out/'raw_histograms.json'],out/'extract.log')
        raw=load(out/'raw_histograms.json');need(raw['passed'] and raw['treeEntriesRead']==0,'extractor failed');need(len(raw['sources'])==len(plan) and raw['events']==sum(x[1] for x in plan),'source/event total differs');edges=raw['nativeEdges'];yedges=raw['rapidityEdges'];need(len(yedges)==b['nBins']+1 and yedges[0]==b['min'] and yedges[-1]==b['max'],'rapidity output edges differ')
        lookup={(x['rootYBin'],x['cent9']):x for x in raw['slices']};need(len(lookup)==(b['nBins']+3)*9,'missing/duplicate native slices')
        curves=[]
        def add(g,cent,index,ys,label,yrange,kind):
            selected=[lookup[(y,c)] for y in ys for c in range(cent['cent9'][0],cent['cent9'][1]+1)];case=dict(caseId='he4Y'+str(g),epBins=g,cutEnabled=False)
            row=pure.curve(case,cent,index,selected,edges,s);row['rapidityBins']=row.pop('epBins');row['rapidityBin']=row.pop('plane');row['mixingEPBins']=s['requiredEPBins'];row['rootYBins']=ys;row['rapidityLabel']=label;row['rapidityRange']=yrange;row['kind']=kind
            for name,w in row['highK'].items():
                a,z=w['windowGeVc'];vs=pure.integral([sum(q['SEVariance'][x] for q in selected) for x in range(1,len(edges))],edges,[a,z]);vm=pure.integral([sum(q['MEVariance'][x] for q in selected) for x in range(1,len(edges))],edges,[a,z]);w['SEVariance']=vs;w['MEVariance']=vm
                alpha=row['normAlpha'];w['fixedAlphaCountError']=(math.sqrt(vs/(alpha*alpha*w['ME']**2)+w['SE']**2*vm/(alpha*alpha*w['ME']**4)) if w['available'] else None)
            curves.append(row)
        for g in [s['coarseBins'],b['nBins']]:
            factor=b['nBins']//g
            for cent in s['centralities']:
                add(g,cent,-1,[-1],'all He4 Y (including Y flow)',None,'inclusive_all_y')
                add(g,cent,-2,list(range(1,b['nBins']+1)),'sum of in-range He4 Y bins',[yedges[0],yedges[-1]],'inclusive_in_range')
                add(g,cent,-3,[0],'He4 Y < '+format(yedges[0],'.6g'),[None,yedges[0]],'rapidity_underflow')
                add(g,cent,-4,[b['nBins']+1],'He4 Y >= '+format(yedges[-1],'.6g'),[yedges[-1],None],'rapidity_overflow')
                for k in range(g):
                    ys=list(range(k*factor+1,(k+1)*factor+1));a,z=yedges[k*factor],yedges[(k+1)*factor];add(g,cent,k,ys,format(a,'.6g')+' <= He4 Y < '+format(z,'.6g'),[a,z],'rapidity_bin')
        by={(x['rapidityBins'],x['centralityId'],x['rapidityBin']):x for x in curves}
        for cent in s['requiredValidCentralities']:
            c=by[(b['nBins'],cent,-1)];need(c['available'] and c['normSE']>0 and c['normME']>0,'main inclusive CF normalization unavailable')
        for cent in s['centralities']:
            a,z=by[(b['nBins'],cent['id'],-1)],by[(s['coarseBins'],cent['id'],-1)];need(a['bins']==z['bins'] and a['normAlpha']==z['normAlpha'],'grouping changes all-Y inclusive result')
        populations=[]
        for cent in s['centralities']:
            control=[by[(b['nBins'],cent['id'],i)] for i in [-1,-2,-3,-4]]
            population=dict(centralityId=cent['id'],centralityLabel=cent['label'],kstarInterval='all native bins including kstar underflow/overflow',pairCut='Off',EPBins=s['requiredEPBins'])
            for mode in ['SE','ME']:
                totals=[x[mode+'Regular']+x[mode+'Underflow']+x[mode+'Overflow'] for x in control]
                need(math.isclose(totals[0],sum(totals[1:]),rel_tol=s['projectionRelativeTolerance'],abs_tol=s['projectionAbsoluteTolerance']),'rapidity population closure failed')
                population[mode]=dict(allY=totals[0],inRange=totals[1],rapidityUnderflow=totals[2],rapidityOverflow=totals[3],rapidityUnderflowFraction=totals[2]/totals[0] if totals[0]>0 else None,rapidityOverflowFraction=totals[3]/totals[0] if totals[0]>0 else None)
            populations.append(population)
        save(out/'he4_y_cf_metrics.json',dict(schemaVersion=1,passed=True,createdUTC=datetime.now(timezone.utc).isoformat(),pilot=args.pilot,events=raw['events'],EPValidEvents=raw['epdEventsValid'],treeEntriesRead=0,normalization='alpha=SE(H)/ME(H); CF=SE/(alpha ME)',rapidityDefinition=raw['rapidityDefinition'],rapidityEdges=yedges,rapidityPairPopulation=populations,settings=s,curves=curves,statisticalLimit='Count errors hold alpha fixed. Normalization covariance and shared-event/shared-ME correlations are not included; no significance test. Y under/overflow are separate controls and remain in all-Y inclusive.'))
        with (out/'curve_metadata.tsv').open('x') as f:
            for c in curves:f.write('\t'.join(map(str,[c['caseId'],c['centralityId'],c['centralityLabel']+' | '+c['rapidityLabel'],c['rapidityBin'],c['rapidityBins'],0,'NA' if c['normAlpha'] is None else repr(c['normAlpha']),c['normSE'],c['normME']]))+'\n')
        with (out/'curve_bins.tsv').open('x') as f:
            for c in curves:
                for i,x in enumerate(c['bins']):f.write('\t'.join(map(str,[c['caseId'],c['centralityId'],c['rapidityBin'],i,x['low'],x['high'],x['SE'],x['SEVariance'],x['ME'],x['MEVariance'],'NA' if x['CF'] is None else repr(x['CF']),'NA' if x['fixedAlphaCountError'] is None else repr(x['fixedAlphaCountError'])]))+'\n')
        with (out/'plot_settings.tsv').open('x') as f:
            for k,z in [('kRange',s['displayKGeVc']),('cfRange',s['displayCFRange']),('normRange',s['normalizationGeVc']),('highKRange',s['highKPlotGeVc']),('highKCFRange',s['highKCFRange'])]:f.write(k+'\t'+'\t'.join(map(str,z))+'\n')
            f.write('mode\the4Y\n')
            for c in s['centralities']:f.write('centrality\t'+c['id']+'\n')
            for k in s['diagnosticCentralities']:f.write('diagnostic\t'+k+'\n')
        pure.root_run(mainconf,HERE/'DrawHe4Y.C',[out/'curve_metadata.tsv',out/'curve_bins.tsv',out/'plot_settings.tsv',out],out/'draw.log')
        for x in source_stats+input_stats:same(x)
        for x in pins:need(sha(x['path'])==x['sha256'],'input/helper pin changed '+x['path'])
        names=['he4_y_cf_metrics.json','raw_histograms.json','he4_y4_CF.pdf','he4_y8_CF.pdf','he4_y4_CF_highk.pdf','he4_y8_CF_highk.pdf','he4_y4_SE_ME_scaled.pdf','he4_y8_SE_ME_scaled.pdf','he4_y_flow_controls.pdf','cf_histograms.root','curve_metadata.tsv','curve_bins.tsv']
        artifacts=[]
        for name in names:rec=stat(out/name);rec['sha256']=sha(out/name);artifacts.append(rec)
        save(out/'completed.json',dict(passed=True,pilot=args.pilot,sources=len(plan),events=raw['events'],curves=len(curves),inputStatsUnchanged=True,pinsUnchanged=True,treeEntriesRead=0,artifacts=artifacts))
        print('HE4_Y_POSTPROCESS_PASS sources=%d curves=%d output=%s'%(len(plan),len(curves),out))
    except Exception as ex:save(out/'failed.json',dict(passed=False,error=str(ex)));raise
if __name__=='__main__':
    try:main()
    except Exception as ex:print('HE4_Y_POSTPROCESS_FAILED: '+str(ex),file=sys.stderr);sys.exit(1)
