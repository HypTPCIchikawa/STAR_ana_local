#!/usr/bin/env python3
"""Create an additive single-He4 lab-rapidity replay profile."""
import argparse,copy,hashlib,json
from pathlib import Path
import yaml
D=Path(__file__).resolve().parent
ROOT=D.parents[1]
P=ROOT/'rootfile/femto_lambda_period_mix100_20261119_001'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return yaml.safe_load(Path(p).read_text())
def write_new(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as out:out.write(data)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bins',type=int,required=True);ap.add_argument('--min',type=float,required=True);ap.add_argument('--max',type=float,required=True);a=ap.parse_args()
    assert a.bins==8 and a.min<a.max
    main_path=ROOT/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml';main=read(main_path)
    refs={k:ROOT/'config'/main[k] for k in ['analysis','maker','femtoHist','mixing','reducedTree','eventPlane','nuclearHist']}
    before={str(p):sha(p) for p in [main_path]+list(refs.values())}
    info=copy.deepcopy(read(refs['analysis']));maker=copy.deepcopy(read(refs['maker']));hist=copy.deepcopy(read(refs['femtoHist']))
    old_hist=copy.deepcopy(hist['histograms']);mix=read(refs['mixing']);reduced=read(refs['reducedTree'])
    assert mix['nEventPlaneBins']==8 and mix['bufferSize']==mix['maxMixEvents']==100
    assert mix['mixingMode']=='bufferAll' and mix['mixBothDirections'] is True
    assert reduced['mode']=='replay' and reduced['resetPoolAtSegmentBoundary'] is False and reduced['requireOriginalMixing'] is False
    assert 'mergingPairCut' not in maker and 'he4RapidityBinning' not in maker
    ana='auau3p85_anaFemtoLambda_4He_tree_he4YLab8_mix100'
    analysis=info['analysis']
    for k in ['anaName','jobName','name','outputFileStem','scratchSubdir']:analysis[k]=ana
    analysis.update(filesPerJob=1,maxEvents=-1,logDir=str(ROOT/'log'/ana),errDir=str(ROOT/'err'/ana),workDir=str(ROOT))
    maker['he4RapidityBinning']=dict(enabled=True,frame='lab',nBins=a.bins,min=a.min,max=a.max)
    axis=dict(nBins=a.bins,min=a.min,max=a.max);hist['axes']['He4SingleYLab']=axis
    names=[]
    for mode in ['SE','ME']:
        for suffix in ['','_signal','_leftSB','_rightSB']:
            ref=hist['histograms']['hKstar'+mode+'VsCent_lambda_he4'+suffix]
            name='hKstar'+mode+'VsHe4YLabVsCent_lambda_he4'+suffix
            assert name not in hist['histograms']
            hist['histograms'][name]=dict(type='TH3D',title=mode+' physical Lambda-He4 '+(suffix or 'full mass')+' accepted pairs;k* [GeV/c];He4 y_lab;cent9',xAxis=copy.deepcopy(ref['xAxis']),yAxis=axis,zAxis=copy.deepcopy(ref['yAxis']))
            names.append(name)
    assert all(hist['histograms'][k]==v for k,v in old_hist.items()) and len(names)==8
    files={
        'analysis':Path('analysis')/('analysis_info_'+ana+'.yaml'),
        'maker':Path('maker')/('maker_'+ana+'.yaml'),
        'femtoHist':Path('hist')/('hist_'+ana+'_kf.yaml')}
    values=dict(analysis=info,maker=maker,femtoHist=hist)
    for k,p in files.items():
        main[k]=p.as_posix();write_new(ROOT/'config'/p,yaml.safe_dump(values[k],sort_keys=False,allow_unicode=True))
    conf=ROOT/'config/mainconf'/('main_'+ana+'.yaml');write_new(conf,yaml.safe_dump(main,sort_keys=False))
    periods=(P/'periods.list').read_bytes();write_new(D/'periods.list',periods.decode())
    manifest=json.loads((P/'input_manifest.json').read_text());sources=manifest['sourceFiles'];assert len(sources)==20
    assert periods.decode().splitlines()==[x['path'] for x in sources]
    for row in sources:
        s=Path(row['path']).stat();assert (s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev)==(row['bytes'],row['mtimeNs'],row['inode'],row['device'])
    assert all(sha(p)==expected for p,expected in before.items())
    settings=dict(schemaVersion=1,frame='lab',nBins=a.bins,min=a.min,max=a.max,coarseBins=4,coarseGroupSize=2,intervalConvention='ROOT half-open [low,high); exact maximum is overflow',retainYUnderflowOverflow=True,referenceEPBins=8,referencePairCutEnabled=False)
    write_new(D/'rapidity_binning.json',json.dumps(settings,indent=2)+'\n')
    case=dict(schemaVersion=1,anaName=ana,outputFileStem=ana,mainconf=str(conf.relative_to(ROOT)),analysisFiles={k:str(ROOT/'config'/v) for k,v in files.items()},histogramNames=names,rapidity=settings,sourceManifest=str(P/'input_manifest.json'),sourceManifestSha256=sha(P/'input_manifest.json'),inputList=str(D/'periods.list'),inputListSha256=sha(D/'periods.list'),sources=sources,storedEvents=sum(x['acceptedEvents'] for x in sources),expectedPeriods=20,outputDirectory=str(ROOT/'rootfile'/ana/'production20261121'),originalReferenceConfigSha256=before,generatedFiles={str(p):sha(p) for p in [conf]+[ROOT/'config'/v for v in files.values()]},baselinePilotReference=str(P/'pilot_B_period100.root'),baselineProductionOutputList=str(P/'output_roots_20261120.list'))
    write_new(D/'analysis_campaign.json',json.dumps(case,indent=2)+'\n')
    scope=json.loads((D/'scope.json').read_text());scope.update(rapidityRangePendingQAAudit=False,rapidity=settings,baseline='EP8/cutOff/depth100/periodContinuous',plannedProductionJobs=20,QAInputCount=127785819,QAInRangeFraction=127671367/127785819,QAPopulation='producer selected single-particle, all cent and EP validity; not pair population');(D/'scope.json').write_text(json.dumps(scope,indent=2)+'\n')
    print(json.dumps(dict(passed=True,mainconf=case['mainconf'],newTH3=8,binning=settings,sourcePeriods=20,storedEvents=case['storedEvents'],baselineConfigsUnchanged=True)))
if __name__=='__main__':main()
