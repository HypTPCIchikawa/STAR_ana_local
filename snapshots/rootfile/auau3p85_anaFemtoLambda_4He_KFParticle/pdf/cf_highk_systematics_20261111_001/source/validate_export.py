import argparse
import hashlib
import json
import math
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--input',required=True)
p.add_argument('--provenance',required=True)
p.add_argument('--validation',required=True)
p.add_argument('--catalog',required=True)
a=p.parse_args()
f=Path(a.input);d=json.load(f.open())
checks={'histogramCount':len(d['histograms']),'schemaCountMatches':len(d['histograms'])==d['histogramCount'],'allShapesMatch':True,'allValuesFinite':True,'allVariancesNonnegative':True,'allEdgesIncreasing':True,'inputUnchanged':True}
cat={}
for n,h in d['histograms'].items():
 k=math.prod(h['shape'])
 checks['allShapesMatch'] &= len(h['content'])==k and len(h['error2'])==k and h['shape'][-1]==len(h['xEdges'])+1
 if h['dimension']==2:checks['allShapesMatch'] &= h['shape'][0]==len(h['yEdges'])+1
 checks['allValuesFinite'] &= all(math.isfinite(x) for x in h['content']) and all(math.isfinite(x) for x in h['error2'])
 checks['allVariancesNonnegative'] &= all(x>=0 for x in h['error2'])
 for key in ['xEdges','yEdges']:
  if key in h:checks['allEdgesIncreasing'] &= all(y>x for x,y in zip(h[key],h[key][1:]))
 cat[n]={k:v for k,v in h.items() if k not in ['content','error2']}
prov=json.load(open(a.provenance));s=Path(prov['input']).stat()
checks['inputUnchanged']=s.st_size==prov['size'] and s.st_mtime_ns==prov['mtimeNs'] and s.st_ino==prov['ino']
checks['rootVersion']=d['rootVersion'];checks['bytes']=f.stat().st_size
sha=hashlib.sha256()
with f.open('rb') as fd:
 for block in iter(lambda:fd.read(1024*1024),b''):sha.update(block)
checks['sha256']=sha.hexdigest()
Path(a.validation).write_text(json.dumps(checks,indent=2)+'\n')
Path(a.catalog).write_text(json.dumps(cat,indent=2)+'\n')
print(json.dumps(checks,indent=2))
if not all(checks[k] for k in ['schemaCountMatches','allShapesMatch','allValuesFinite','allVariancesNonnegative','allEdgesIncreasing','inputUnchanged']):raise SystemExit(2)
