#!/usr/bin/env python3
"""Prepare seven gated cases from the verified one-period launch workflow.

This writes configurations/joblists only. Submissions are separate, one-shot
operations. The aggregate watcher owns completion mail and CF processing.
"""
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import sys
import xml.etree.ElementTree as ET
E = Path(__file__).resolve().parent
ROOT = E.parents[1]
P = ROOT / 'rootfile/femto_lambda_period_mix100_20261119_001'
LABEL = 'epCut20261120'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)

def save(p, obj):
    with Path(p).open('x') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write('\n')

def main():
    assert not (E / 'production_matrix.json').exists()
    for name in ['pilot_comparison.json', 'test_actual_lambda_merging_paircut_result.json']:
        assert json.loads((E / name).read_text())['passed'] is True, name
    cases = json.loads((E / 'cases.json').read_text())
    assert cases['expectedPeriods'] == 20 and cases['plannedProductionJobs'] == 140
    assert sha(cases['inputManifest']) == cases['inputManifestSha256']
    sources = json.loads(Path(cases['inputManifest']).read_text())['sourceFiles']
    template = (P / 'launch_once.py').read_text()
    template = replace_once(template, 'ROOT = STUDY.parent.parent', 'ROOT = Path(' + repr(str(ROOT)) + ')')
    template = replace_once(template, 'STUDY.name == "femto_lambda_period_mix100_20261119_001"',
        'STUDY.parent == Path(' + repr(str(E / 'production')) + ') and STUDY.name == config["caseId"]')
    template = replace_once(template, 'config["anaName"].endswith("_tree_periodMix100")',
        'config["anaName"] == "auau3p85_anaFemtoLambda_4He_tree_" + config["caseId"] + "_mix100"\n'
        '            and re.fullmatch(r"ep(1|8|16|32)_merging(On|Off)", config["caseId"])\n'
        '            and config["caseId"] != "ep8_mergingOff"')
    template = replace_once(template, 'config["label"] == "periodMix10020261119"', 'config["label"] == "' + LABEL + '"')
    template = replace_once(template, 'float(xml.get("filesPerHour", "0")) == 0.025',
        'float(xml.get("filesPerHour", "0")) == 1.0 / 48.0')
    template = replace_once(template, 'reviewed 40-hour-per-period scheduler estimate',
        'reviewed 48-hour-per-period scheduler estimate')
    begin = template.index('        monitor = start_monitor(config, meta, runmeta, event, state)')
    end = template.index('    except Exception as error:', begin)
    template = template[:begin] + ('        print(json.dumps({"phase":"submitted_verified", "jobid":jobid,\n'
        '            "actualJobs":EXPECTED_PERIODS,"runmeta":str(runmeta),"eventLog":event,\n'
        '            "automaticResubmit":False,"automaticMerge":False}), flush=True)\n'
        '        return 0\n') + template[end:]
    template = replace_once(template, 'The existing completion-only notifier sends one minimal email;\nthis launcher does not implement mail delivery or edit shared job tooling.',
        'The aggregate matrix watcher owns completion mail and CF postprocessing;\nthis launcher performs submission verification and returns.')
    prod = E / 'production'
    prod.mkdir(exist_ok=False)
    matrix = dict(schemaVersion=1, projectRoot=str(ROOT), studyDirectory=str(E), expectedPeriods=20,
        plannedProductionJobs=140, casesManifest=str(E/'cases.json'), casesManifestSha256=sha(E/'cases.json'),
        baselineReuseCaseId='ep8_mergingOff', cases=[])
    resource_text = ('<ResourceUsage><Memory><MinMemory>4096</MinMemory><MaxMemory>4096</MaxMemory></Memory>'
        '<StorageSpace><MinStorage>65536</MinStorage><MaxStorage>65536</MaxStorage></StorageSpace></ResourceUsage>')
    base_conf = json.loads((P / 'launch_config.json').read_text())
    for case in cases['cases']:
        if not case['productionRequired']:
            continue
        cid, ana = case['caseId'], case['anaName']
        d = prod / cid
        d.mkdir()
        (d/'launch_once.py').write_text(template)
        (d/'periods.list').write_bytes((E/'periods.list').read_bytes())
        (d/'pilot_comparison.json').write_bytes((E/'pilot_comparison.json').read_bytes())
        output = Path(case['outputDirectory'])
        output.mkdir(parents=True, exist_ok=True)
        assert not list(output.glob('*.root'))
        xml_path = ROOT/'job/joblist'/('joblist_' + ana + '_' + LABEL + '.xml')
        assert not xml_path.exists(), str(xml_path)
        subprocess.run([sys.executable, str(ROOT/'script/generate_joblist_from_list.py'),
            '--mainconf',case['mainconf'],'--input',str(d/'periods.list'),'--files-per-subjob','1',
            '--label',LABEL,'--output-dir',str(output)], cwd=str(ROOT), check=True)
        xml = ET.parse(xml_path)
        job = xml.getroot()
        assert job.get('maxFilesPerProcess') == '1' and job.find('ResourceUsage') is None
        job.set('minFilesPerProcess','1')
        job.set('filesPerHour',str(1.0/48.0))
        resource = ET.fromstring(resource_text)
        job.insert(list(job).index(job.find('command'))+1,resource)
        xml.write(str(xml_path),encoding='utf-8',xml_declaration=True)
        fixed = [Path(p) for p in base_conf['fingerprints']
                 if not str(p).startswith(str(P)) and '/job/joblist/' not in p and '/config/' not in p]
        fixed += [Path(__file__),d/'launch_once.py',xml_path,ROOT/case['mainconf'],E/'cases.json',
            E/'test_actual_lambda_merging_paircut_result.json',ROOT/'include/FemtoLambdaPairGeometry.h',
            ROOT/'include/FemtoEpdEventPlane.h',ROOT/'Makefile']
        for line in (ROOT/case['mainconf']).read_text().splitlines():
            if ':' in line and not line.lstrip().startswith('#'):
                value = line.split(':',1)[1].strip().strip("\"'")
                if value.endswith('.yaml'):
                    fixed.append(ROOT/'config'/value)
        fixed += [Path(p) for p in base_conf['fingerprints']
                  if '/config/event_plane/xialei_epd_2021P24iy/' in p or '/config/calibration/' in p]
        conf = dict(base_conf)
        conf.update(caseId=cid, directory=str(d),python=sys.executable,submitHost=socket.getfqdn(),
            anaName=ana,label=LABEL,submitName=ana+'_'+LABEL,outputFileStem=case['outputFileStem'],
            mainconf=case['mainconf'],joblist=str(xml_path),outputDirectory=str(output),
            inputList=str(d/'periods.list'),inputListSha256=sha(d/'periods.list'),sources=sources,
            fingerprints={str(p.resolve()):sha(p) for p in fixed},
            pilotComparison=dict(path=str(d/'pilot_comparison.json'),sha256=sha(d/'pilot_comparison.json')),
            resourceUsageXml=ET.tostring(resource,encoding='unicode'),
            aggregateMonitor=str(E/'watch_production_matrix.py'))
        save(d/'config.json',conf)
        matrix['cases'].append(dict(caseId=cid,launchDirectory=str(d),config=str(d/'config.json')))
    assert len(matrix['cases']) == 7
    matrix['validation'] = {name:sha(E/name) for name in
        ['pilot_comparison.json','test_actual_lambda_merging_paircut_result.json']}
    save(E/'production_matrix.json',matrix)
    print(json.dumps(dict(prepared=True,productionCases=7,productionJobs=140,submitted=False)))

if __name__ == '__main__':
    main()
