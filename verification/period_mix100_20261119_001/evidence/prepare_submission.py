#!/usr/bin/env python3
"""Regenerate and pin the one-period-per-worker joblist after pilot validation."""
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import sys
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
STUDY = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    comparison = STUDY / 'pilot_comparison.json'
    assert json.loads(comparison.read_text())['passed'] is True
    assert not (STUDY / 'launch_config.json').exists()
    ana = 'auau3p85_anaFemtoLambda_4He_tree_periodMix100'
    label = 'periodMix10020261119'
    conf = 'config/mainconf/main_' + ana + '.yaml'
    output = ROOT / 'rootfile' / ana / 'production20261119'
    output.mkdir(parents=True, exist_ok=True)
    assert not list(output.glob('*.root'))
    xml_path = ROOT / 'job/joblist' / ('joblist_' + ana + '_' + label + '.xml')
    if xml_path.exists():
        with (STUDY / 'preliminary_joblist.xml').open('xb') as saved:
            saved.write(xml_path.read_bytes())
    subprocess.run([sys.executable, str(ROOT / 'script/generate_joblist_from_list.py'),
                    '--mainconf', conf, '--input', str(STUDY / 'periods.list'),
                    '--files-per-subjob', '1', '--label', label, '--output-dir', str(output)],
                   cwd=str(ROOT), check=True)
    xml = ET.parse(xml_path)
    job = xml.getroot()
    assert job.attrib['maxFilesPerProcess'] == '1'
    job.set('minFilesPerProcess', '1')
    job.set('filesPerHour', '0.025')
    resources = ET.fromstring('<ResourceUsage><Memory><MinMemory>4096</MinMemory><MaxMemory>4096</MaxMemory></Memory>'
                             '<StorageSpace><MinStorage>65536</MinStorage><MaxStorage>65536</MaxStorage></StorageSpace></ResourceUsage>')
    assert job.find('ResourceUsage') is None
    job.insert(list(job).index(job.find('command')) + 1, resources)
    xml.write(str(xml_path), encoding='utf-8', xml_declaration=True)
    manifest = json.loads((STUDY / 'input_manifest.json').read_text())
    for relative, expected in manifest['originalTrackedChangesSha256'].items():
        assert sha(ROOT / relative) == expected, 'Pre-existing user file changed: ' + relative
    fixed = [Path(__file__), STUDY / 'launch_once.py', xml_path,
             ROOT / 'job/run/submit.sh', ROOT / 'script/watch_job_notify.py',
             ROOT / 'analysis/FemtoLambdaRunSupport.h',
             ROOT / 'analysis/anaFemtoLambda_4He_tree.C', ROOT / 'analysis/run_anaFemtoLambda_4He_tree.C',
             ROOT / 'include/FemtoLambdaReducedTree.h', ROOT / 'include/FemtoLambdaReducedTreeConfig.h',
             ROOT / 'include/FemtoLambdaTreeTypes.h', ROOT / 'include/FemtoLambdaTreeLinkDef.h',
             ROOT / 'include/FemtoCandidateEventStore.h', ROOT / 'include/FemtoCandidateEventStoreConfig.h',
             ROOT / 'StMaker/StFemtoMaker/StFemtoMaker.cxx', ROOT / 'StMaker/StFemtoMaker/StFemtoMaker.h',
             ROOT / 'StMaker/common/FemtoLambdaLegacy.cxx', ROOT / 'StMaker/common/FemtoLambdaLegacy.h',
             ROOT / conf]
    for line in (ROOT / conf).read_text().splitlines():
        if ':' in line and not line.lstrip().startswith('#'):
            value = line.split(':', 1)[1].strip()
            if value.endswith('.yaml'):
                fixed.append(ROOT / 'config' / value)
    fixed += [ROOT / 'lib' / name for name in ['libStarAnaConfig.so', 'libStRefMultCorr.so',
              'libKFParticle.so', 'libStKfParticleCommon.so', 'libStCommon.so', 'libStFemtoMaker.so']]
    fixed += [ROOT / 'config/event_plane/xialei_epd_2021P24iy' / name
              for name in ['recenterpar.root', 'shiftpar.root', 'runs_2021.txt']]
    fixed.append(ROOT / 'config/calibration/dedx_he4_auau3p85_2021P24iy_accepted_20261022.yaml')
    config = dict(schemaVersion=1, projectRoot=str(ROOT), directory=str(STUDY), submitHost=socket.getfqdn(),
                  python=sys.executable, anaName=ana, label=label, submitName=ana + '_' + label,
                  outputFileStem=ana, mainconf=conf, joblist=str(xml_path), outputDirectory=str(output),
                  inputList=str(STUDY / 'periods.list'), inputListSha256=sha(STUDY / 'periods.list'),
                  sources=manifest['sourceFiles'], fingerprints={str(p.resolve()): sha(p) for p in fixed},
                  pilotComparison=dict(path=str(comparison), sha256=sha(comparison)),
                  quotaCommand=['/usr/lpp/mmfs/bin/mmlsquota', '-u', 'yichikawa', '-Y', 'rcf.bnl.gov:gpfs:star-pwg'],
                  minimumSoftHeadroomBytes=500000000000,
                  resourceUsageXml=ET.tostring(resources, encoding='unicode'),
                  nativeResources=dict(request_memory=4096, request_disk=65536000),
                  email='yudai.ichikawa.d3@tohoku.ac.jp', emailContent='minimal', pollSeconds=300, settleSeconds=300)
    with (STUDY / 'launch_config.json').open('x') as saved:
        json.dump(config, saved, indent=2)
        saved.write('\n')
    print(json.dumps(dict(joblist=str(xml_path), sourceCount=len(config['sources']),
                          pinnedFiles=len(config['fingerprints']), resourceUsage=config['resourceUsageXml'])))


if __name__ == '__main__':
    main()
