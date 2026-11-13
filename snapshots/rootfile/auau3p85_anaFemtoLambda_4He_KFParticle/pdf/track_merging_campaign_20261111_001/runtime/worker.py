#!/usr/bin/env python3
"""Stage one original candidate job and its PicoDst inputs, then rebuild QA only."""
import argparse, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time, traceback, uuid
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def write(path,obj):
    p=Path(path);tmp=p.with_name(p.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
    os.replace(str(tmp),str(p))


def quote_root(text):
    return '"'+str(text).replace('\\','\\\\').replace('"','\\"')+'"'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['runtime','cache','pico-list','output-dir','policy','scratch-base']:
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args();started=time.time();runtime=args.runtime.resolve(strict=True)
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
    owner=uuid.uuid4().hex;lock=output/'worker.lock'
    with lock.open('x') as f:f.write(owner+'\n')
    work=None;receipt={'status':'starting','startedUnix':started,'pid':os.getpid(),'host':os.uname().nodename,'owner':owner,'runtime':str(runtime),'argv':sys.argv,'sourceCache':str(args.cache.resolve()),'sourcePicoList':str(args.pico_list.resolve()),'sourceList':str(args.pico_list.resolve()),'policy':str(args.policy.resolve())}
    try:
        manifest=json.loads((runtime/'manifest.json').read_text())
        for rel,item in manifest['files'].items():
            p=runtime/rel
            if not p.is_file() or sha(p)!=item['sha256']:raise RuntimeError('Frozen runtime mismatch: '+rel)
        env=json.loads((runtime/'environment.json').read_text());policy=json.loads(args.policy.read_text())
        receipt['runtimeManifestSha256']=sha(runtime/'manifest.json');receipt['policySha256']=sha(args.policy);receipt['inputListSha256']=sha(args.pico_list)
        source=args.cache.resolve(strict=True)
        if not source.is_file() or source.stat().st_size<=0:raise RuntimeError('Missing/empty original candidate ROOT')
        source_stat=source.stat()
        actual_fp={'size':source_stat.st_size,'mtimeNs':source_stat.st_mtime_ns,'ino':source_stat.st_ino,'device':source_stat.st_dev}
        if 'expectedCacheFingerprint' in policy and any(k not in policy['expectedCacheFingerprint'] or actual_fp[k]!=policy['expectedCacheFingerprint'][k] for k in ['size','mtimeNs','ino']):raise RuntimeError('Candidate inventory fingerprint changed')
        if 'expectedInputListSha256' in policy and receipt['inputListSha256']!=policy['expectedInputListSha256']:raise RuntimeError('Original input list changed after inventory')
        receipt['sourceCacheFingerprint']=actual_fp
        receipt['sourceCacheSha256']=sha(source)
        local=bool(policy.get('localValidation',False));allow=json.loads((runtime/'validation/local_allowlist.json').read_text())
        if local and (str(source) not in allow or receipt['sourceCacheSha256']!=allow[str(source)]):raise RuntimeError('Local validation cache is not the frozen smoke input')
        rows=[]
        for text in args.pico_list.read_text().splitlines():
            text=text.split('#',1)[0].strip()
            if not text:continue
            fields=text.split()
            if len(fields) not in [1,2] or (len(fields)==2 and (not fields[1].isdigit() or int(fields[1])<1)):raise RuntimeError('Malformed original Pico input list')
            url=fields[0]
            if local:
                p=Path(url).resolve(strict=True)
                if str(p) not in allow or sha(p)!=allow[str(p)]:raise RuntimeError('Local Pico validation input is not allowlisted')
            elif not (url.startswith('root://') and '/'+env['dataset']+'/' in url and '/'+env['production']+'/' in url and url.endswith('.picoDst.root')):
                raise RuntimeError('Input is not an explicit 2021 P24iy PicoDst URI: '+url)
            rows.append({'source':url,'declaredEntries':int(fields[1]) if len(fields)==2 else None})
        if not rows or len(set(x['source'] for x in rows))!=len(rows):raise RuntimeError('Empty or duplicate Pico input list')
        scratch=args.scratch_base.resolve(strict=True)
        if not scratch.is_dir() or str(scratch).startswith('/gpfs/') or str(scratch).startswith('/star/u/'):
            raise RuntimeError('Scratch base must be a node-local directory')
        work=Path(tempfile.mkdtemp(prefix='mergingqa_',dir=str(scratch)))
        (work/'owner.json').write_text(json.dumps({'owner':owner,'output':str(output)}))
        receipt['scratch']=str(work);write(output/'running.json',receipt)
        stage_started=time.time();staged_cache=work/'candidate.root';shutil.copyfile(str(source),str(staged_cache))
        if sha(staged_cache)!=receipt['sourceCacheSha256']:raise RuntimeError('Candidate byte-copy checksum mismatch')
        staged=[]
        for i,row in enumerate(rows):
            dst=work/('pico_%04d.picoDst.root'%i);one=time.time()
            if local:shutil.copyfile(row['source'],str(dst))
            else:
                cmd=[env['xrdcp'],'--nopbar','--retry','0','--cksum','adler32',row['source'],str(dst)]
                with (output/('stage_%04d.log'%i)).open('xb') as f:
                    done=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=False,timeout=policy['stageTimeoutSeconds'])
                if done.returncode:raise RuntimeError('xrdcp failed for input %d, rc=%d'%(i,done.returncode))
            if not dst.is_file() or dst.stat().st_size<=0:raise RuntimeError('Staged PicoDst is empty')
            staged.append(dict(row,staged=str(dst),bytes=dst.stat().st_size,stageSeconds=time.time()-one,checksumVerification='frozen SHA256 allowlist' if local else 'xrdcp --cksum adler32'))
        receipt['inputs']=staged;receipt['stageSeconds']=time.time()-stage_started;receipt['picoStageBytes']=sum(x['bytes'] for x in staged);receipt['candidateStageBytes']=staged_cache.stat().st_size
        pico_list=work/'pico.list';pico_list.write_text(''.join(x['staged']+'\n' for x in staged))
        cfg=dict(policy)
        for k in ['localValidation']:cfg.pop(k,None)
        main=Path(cfg['sourceMainconf']);main=main if main.is_absolute() else runtime/main
        if not main.is_file():raise RuntimeError('Frozen source mainconf is missing')
        cfg.update(sourceMainconf=str(main),candidateRoot=str(staged_cache),picoList=str(pico_list),outputRoot=str(work/'qa.root'),outputJson=str(work/'summary.json'))
        settings=work/'settings.json';write(settings,cfg);write(output/'effective_settings.json',cfg)
        star_lib=env['starRoot']+'/.'+env['starHostSys']+'/lib'
        ld=':'.join([str(runtime/'lib'),star_lib,env['rootSys']+'/lib']+env['externalLibraryPaths'])
        call=str(runtime/'source/InvokeFrozenQa.C')+'('+quote_root(runtime)+','+quote_root(settings)+')'
        cmd=[env['singularity'],'exec','--pwd',str(work),'-B','/gpfs:/gpfs','-B','/cvmfs:/cvmfs','-B','/star/nfs4/AFS:/star/nfs4/AFS','-B',str(work)+':'+str(work),'--env','ROOTSYS='+env['rootSys'],'--env','STAR='+env['starRoot'],'--env','STAR_HOST_SYS='+env['starHostSys'],'--env','LD_LIBRARY_PATH='+ld,env['image'],env['root4star'],'-l','-b','-q',call]
        receipt['qaCommand']=cmd;qa_start=time.time()
        with (output/'qa.log').open('xb') as f:done=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=False)
        receipt['qaProcessSeconds']=time.time()-qa_start
        if done.returncode:raise RuntimeError('QA process failed, rc='+str(done.returncode))
        if not (work/'summary.json').is_file() or not (work/'qa.root').is_file():raise RuntimeError('QA process did not produce validated outputs')
        summary=json.loads((work/'summary.json').read_text())
        if not summary.get('passed') or summary.get('closureChecks')!=244827:raise RuntimeError('QA closure did not pass')
        if (work/'qa.root').read_bytes()[:4]!=b'root':raise RuntimeError('Invalid output ROOT header')
        receipt['qaSummary']=summary;receipt['outputRootSha256']=sha(work/'qa.root');receipt['outputRootBytes']=(work/'qa.root').stat().st_size
        for name in ['qa.root','summary.json']:
            target=output/name
            if target.exists():raise RuntimeError('Refuse existing published result '+name)
            temp=output/(name+'.'+owner+'.partial');shutil.copyfile(str(work/name),str(temp));os.rename(str(temp),str(target))
        receipt.update(status='completed',wallSeconds=time.time()-started,finishedUnix=time.time(),picoFileCount=len(staged),physicsSelectionChanged=False,kfReconstructionRepeated=False)
        write(output/'receipt.json',receipt);write(output/'worker_result.json',receipt)
        print('CAMPAIGN_WORKER_PASS',output,'Picos',len(staged),'bytes',receipt['picoStageBytes'],'elapsed',receipt['wallSeconds'],flush=True)
        return 0
    except Exception as e:
        receipt.update(status='failed',error=str(e),wallSeconds=time.time()-started,finishedUnix=time.time(),traceback=traceback.format_exc())
        write(output/'failed.json',receipt);write(output/'worker_result.json',receipt);print('CAMPAIGN_WORKER_FAILED',e,file=sys.stderr,flush=True);return 1
    finally:
        if work is not None and work.is_dir():
            marker=json.loads((work/'owner.json').read_text())
            if marker.get('owner')==owner and work.parent==args.scratch_base.resolve() and work.name.startswith('mergingqa_'):
                shutil.rmtree(str(work))
            else:raise RuntimeError('Refusing cleanup of unowned scratch')

if __name__=='__main__':sys.exit(main())
