#!/usr/bin/env python3
"""Validate Poisson mass fits, constant-background fallback and CF publication."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from build_method6_report import pdf_pages


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(a,b):return math.isclose(a,b,rel_tol=2.e-7,abs_tol=2.e-7)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve(strict=True)
    context=json.loads((out/'poisson_context.json').read_text())
    root=Path(context['input']['path']);stat=root.stat()
    assert stat.st_size==context['input']['size'] and stat.st_mtime_ns==context['input']['mtime_ns']
    assert context['negativeBackgroundConstant'] is True
    for record in context['references'].values():assert sha(record['path'])==record['sha256']
    for path,digest in context['sourceCopies'].items():assert sha(path)==digest
    if 'metadataCorrection' in context:
        for record in json.loads(Path(context['metadataCorrection']).read_text())['files']:
            assert sha(record['path'])==record['afterSha256']
    if (out/'publication.json').exists():
        published=json.loads((out/'publication.json').read_text())
        assert sha(published['backup'])==context['targetReport']['sha256']
        assert sha(published['target'])==published['sha256']
    else:
        assert sha(context['targetReport']['path'])==context['targetReport']['sha256']
    checked_backgrounds=0;fallbacks=0;at_zero=0;matched=0
    for name in ('group','fine'):
        data=json.loads((out/(name+'_extract/purity_fit_review.json')).read_text())
        previous=json.loads(Path(context['references'][name]['path']).read_text())
        assert data['settings']==previous['settings']
        assert data['fitOptions']=='SLQNR0'
        partition=data['massPartition']['binned'];width=(data['massEdges'][-1]-data['massEdges'][0])/(len(data['massEdges'])-1)
        for group,old in zip(data['groups'],previous['groups']):
            for method in ('method3','method5'):
                for side in ('SE','ME'):
                    for row,before in zip(group['methods'][method][side],old['methods'][method][side]):
                        for key in ('massCounts','massError2','centralMassCounts','centralMassError2','centralNtot','windowNtot'):
                            assert row[key]==before[key]
                        matched+=1;fit=row['fit']
                        if fit['status']==-1:continue
                        if fit['constantBackgroundFallback']:
                            fallbacks+=1
                            assert fit['primaryBackgroundMinimum']<0
                            assert fit['freeParameters']==4 and fit['ndf']==fit['usedFitBins']-4
                            assert fit['basisParameters'][1:]==[0,0]
                            assert fit['parameters'][4:]==[0,0]
                            assert fit['basisParameters'][0]>=0
                            assert fit['basisParameters'][0]<fit['constantUpperBound']
                            if fit['basisParameters'][0]<1.e-8:at_zero+=1
                        else:
                            assert fit['primaryBackgroundMinimum']>=0
                            assert fit['freeParameters']==6 and fit['ndf']==fit['usedFitBins']-6
                        if row['background'][0] is None:continue
                        b=fit['basisParameters'];cov=fit['basisCovariance'];q=[]
                        for region in ('low','high'):
                            limits=partition[region];lo,hi=limits['lower'],limits['upper']
                            u=(lo-fit['basisMidpoint'])/fit['basisScale'];v=(hi-fit['basisMidpoint'])/fit['basisScale']
                            g=(hi-lo)/width;q.append([g,g*(u+v)/2,g*(u*u+u*v+v*v)/3])
                        for i in range(2):
                            assert close(sum(q[i][j]*b[j] for j in range(3)),row['background'][i])
                            for j in range(2):
                                expected=sum(q[i][a]*cov[a][b]*q[j][b] for a in range(3) for b in range(3))
                                assert close(expected,row['backgroundCovariance'][i][j])
                        checked_backgrounds+=1
                        if row['full']['valid']:
                            assert fit['expectationValid'] and fit['backgroundMinimum']>=0
    assert matched==1960
    analysis=json.loads((out/'analysis_validation.json').read_text());assert analysis['matchedCountRows']==1960
    results=json.loads((out/'methods_comparison_selected.json').read_text());cf_points=0
    for group in results['groups']:
        for method in ('method3','method5','method6'):
            r=group['methods'][method];assert r['normalization']['valid']
            for i,p in enumerate(r['points']):
                if p['cf'] is None:continue
                cf_points+=1;assert close(p['cf']['error']**2,p['cf']['variance'])
                assert close(r['covariance'][i][i],p['cf']['variance'])
        shared=group['methods']['method6']['purity']
        assert shared['SE']==shared['ME']
        for p,q in zip(shared['ME'],group['methods']['method5']['purity']['ME']):
            assert (p is None)==(q is None)
            if p is not None:assert close(p['value'],q['value']) and close(p['error'],q['error'])
    toy=json.loads((out/'empty_bin_toy.json').read_text())
    assert abs(toy['chiSquareEstimate']-1)<.001 and abs(toy['poissonEstimate']-.1)<.001
    metadata=json.loads((out/'poisson_report_metadata.json').read_text())
    if (out/'publication.json').exists():
        assert metadata['report']['sha256']==json.loads((out/'publication.json').read_text())['sha256']
    for product in [metadata['report']]+metadata['figures']:
        assert sha(product['path'])==product['sha256']
        assert pdf_pages(Path(product['path']))==product['pages']
    crops=0
    for path in out.glob('*/crop_manifest.json'):
        for name,record in json.loads(path.read_text())['files'].items():
            assert sha(record['source'])==record['sourceSha256']
            assert sha(path.parent/name)==record['outputSha256'];crops+=1
    record=dict(phase='validated',reportSha256=metadata['report']['sha256'],reportPages=metadata['report']['pages'],
        inputRootUnchanged=True,matchedMassFitRows=matched,backgroundIntegralsAndCovariancesChecked=checked_backgrounds,
        constantFallbackFits=fallbacks,constantAtZeroFits=at_zero,validCfPoints=cf_points,
        cropProducts=crops,independentDevianceChecks=analysis['independentDevianceChecks'],
        sharedMethod6PurityAndErrorsVerified=True,all15SelectedNormalizationsValid=True)
    (out/'validation.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))


if __name__=='__main__':main()
