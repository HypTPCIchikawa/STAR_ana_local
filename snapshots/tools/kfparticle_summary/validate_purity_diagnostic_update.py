#!/usr/bin/env python3
"""Verify staged diagnostic figures, source continuity and report assembly."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from build_method6_report import pdf_pages


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--fine-review', type=Path, required=True)
    parser.add_argument('--reference-review', type=Path, required=True)
    parser.add_argument('--reference-context', type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve(strict=True)
    metadata = json.loads((out/'diagnostic_report_metadata.json').read_text())
    cited = json.loads((out/'cited_sources.json').read_text())
    for source in cited['files'].values():
        assert sha(source['path']) == source['sha256'], source['path']
    context = json.loads(args.reference_context.read_text())
    old = json.loads(args.reference_review.read_text())
    fine = json.loads(args.fine_review.read_text())
    a,b = dict(old['settings']),dict(fine['settings'])
    a.pop('groups');b.pop('groups')
    assert a==b, 'Mass fit or native-bin settings changed'
    assert len(fine['groups'])==9, 'Expected nine fine centrality groups'
    assert len({g['id'] for g in fine['groups']})==9
    fit_count=sum(len(g['methods'][m][s]) for g in fine['groups']
                  for m in ('method3','method5') for s in ('SE','ME'))
    assert fit_count==1260, fit_count
    for key in ('massEdges','originalKstarEdges'):
        assert fine[key]==old[key], key
    # A shared centrality interval provides an exact input-projection cross-check.
    fine6=next(g for g in fine['groups'] if g['centMin']==g['centMax']==6)
    old6=next(g for g in old['groups'] if g['centMin']==g['centMax']==6)
    mass_rows=0
    for method in ('method3','method5'):
        for side in ('SE','ME'):
            first,second=fine6['methods'][method][side],old6['methods'][method][side]
            assert len(first)==len(second)
            for x,y in zip(first,second):
                for key in ('massCounts','centralMassCounts','centralNtot','windowNtot'):
                    assert x[key]==y[key], (method,side,x['bin'],key)
                mass_rows+=1
    root=Path(context['sample']['inputRoot'])
    fingerprint=context['input']
    # The context's input record contains the ROOT fingerprint directly.
    if 'fingerprint' in fingerprint:
        fingerprint=fingerprint['fingerprint']
    stat=root.stat()
    for key,value in [('size',stat.st_size),('mtime_ns',stat.st_mtime_ns),('inode',stat.st_ino)]:
        if key in fingerprint:
            assert fingerprint[key]==value, 'Input ROOT changed: '+key
    assert context['sample']['year']==2021
    assert context['sample']['submittedAnalysis']['starTag']['productionTag']=='P24iy'
    assert context['sample']['jobid']=='04BF2DAE47164B6F040E3A460E539289'
    audit=json.loads((out/'normalization_audit.json').read_text())
    selected=json.loads((out/'methods_comparison_selected.json').read_text())
    chosen=next(c for c in audit['candidates'] if c['id']==audit['settings']['selected'])
    assert chosen['normalization']==[0.4,0.5]
    assert len(selected['groups'])==5
    checked_cf=0
    for group in selected['groups']:
        methods=group['methods']
        for method in ('method3','method5','method6'):
            result=methods[method]
            assert result['normalization']['valid'], (group['id'],method)
            for i,point in enumerate(result['points']):
                if point['cf'] is None:
                    continue
                cf=point['cf'];checked_cf+=1
                assert cf['variance']>=0 and math.isfinite(cf['error'])
                assert math.isclose(cf['error']**2,cf['variance'],rel_tol=1.e-10,abs_tol=1.e-14)
                assert math.isclose(result['covariance'][i][i],cf['variance'],rel_tol=1.e-10,abs_tol=1.e-14)
        shared=methods['method6']['purity']
        assert shared['SE']==shared['ME']
        for x,y in zip(shared['ME'],methods['method5']['purity']['ME']):
            assert (x is None)==(y is None)
            if x is not None:
                assert math.isclose(x['value'],y['value'],rel_tol=1.e-10,abs_tol=1.e-12)
                assert math.isclose(x['error'],y['error'],rel_tol=1.e-10,abs_tol=1.e-12)
        if group['id']=='cent20_80':
            for method in ('method5','method6'):
                for index in (15,16,17):
                    assert methods[method]['points'][index]['cf'] is None, 'Failed purity cell was filled'
    baseline=next(c for c in audit['candidates'] if c['normalization']==[0.3,0.4])
    peripheral=next(g for g in baseline['groups'] if g['id']=='cent20_80')
    assert all(not peripheral['methods'][m]['normalization']['valid'] for m in ('method5','method6'))
    cases=json.loads((out/'selected_purity_cases.json').read_text())
    assert len(cases['cases'])==9
    assert cases['cases'][0]['category']=='normalization_failure'
    assert cases['cases'][0]['ME']['fitRangeObservations']==audit['failedFits'][0]['observationsInFitRange']
    for name,record in json.loads((out/'diagnostic_source_code.json').read_text()).items():
        assert sha(record['path'])==record['sha256'], name
    products=[metadata['report'],metadata['addendum']]+metadata['figures']
    for product in products:
        assert sha(product['path'])==product['sha256'],product['path']
        assert pdf_pages(Path(product['path']))==product['pages']
    for key in ('sourceManifest','content','sourceCode'):
        record=metadata[key];assert sha(record['path'])==record['sha256']
    assert metadata['report']['pages']==1+metadata['originalPages']+metadata['addendum']['pages']
    crop_count=0
    for path in out.glob('*/crop_manifest.json'):
        crop=json.loads(path.read_text())
        for name,record in crop['files'].items():
            assert sha(record['source'])==record['sourceSha256']
            assert sha(path.parent/name)==record['outputSha256']
            assert pdf_pages(path.parent/name)==record['pages']
            crop_count+=1
    report=dict(phase='validated',reportSha256=metadata['report']['sha256'],
        reportPages=metadata['report']['pages'],addendumPages=metadata['addendum']['pages'],
        citedPdfHashes=len(cited['files']),verifiedProducts=len(products),
        cropProducts=crop_count,fineCentralityGroups=9,extractedFits=fit_count,
        matchedFine10_20MassRows=mass_rows,inputRootUnchanged=True,
        validSelectedCfPoints=checked_cf,selectedNormalizationsValid=15,
        sharedMethod6PurityMatchesMethod5ME=True,failedPurityBinsRemainMissing=True,
        massFitSettingsUnchanged=True,fitGateRelaxation=False,
        fineReviewSha256=sha(args.fine_review),referenceReviewSha256=sha(args.reference_review))
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':
    main()
