#!/usr/bin/env python3
"""Freeze the input identity and write the explicit raw-projection plan."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cfg = yaml.safe_load((HERE / 'study.yaml').read_text())
    state = json.loads((ROOT / cfg['recovery_state']).read_text())
    assert state['phase'] == 'completed' and state['commonJobs'] > 0
    recovery = json.loads((ROOT / cfg['recovery_config']).read_text())
    known = {e['key'].upper().replace('EPD', 'EP'): e for e in recovery['entries']}
    conditions = [dict(id='inclusive', dimension=2, first=0, last=0, suffix='', label='Inclusive'),
                  dict(id='he4_below', dimension=2, first=0, last=0, suffix='_YLabBelow', label='Y(He4)<-0.5'),
                  dict(id='he4_above', dimension=2, first=0, last=0, suffix='_YLabAbove', label='Y(He4)>-0.5')]
    rapidity = cfg['pair_rapidity']
    for i in range(1, rapidity['n_bins'] + 1):
        conditions.append(dict(id='pair_bin%02d' % i, dimension=3, first=i, last=i, suffix='', label='Pair Y bin %d' % i))
    for i, (first, last) in enumerate(rapidity['adjacent_merge_groups'], 1):
        conditions.append(dict(id='pair_merge%02d' % i, dimension=3, first=first, last=last, suffix='', label='Merged pair Y %d-%d' % (first,last)))
    for name, first, last in [('pair_under',0,0),('pair_over',rapidity['n_bins']+1,rapidity['n_bins']+1),
                              ('pair_regular',1,rapidity['n_bins']),('pair_all',0,rapidity['n_bins']+1)]:
        conditions.append(dict(id=name, dimension=3, first=first, last=last, suffix='', label=name))
    rows, inputs = [], {}
    for dataset, value in cfg['inputs'].items():
        path = ROOT / value['root']; k = known[dataset]
        assert path == Path(k['destination'])
        st = path.stat()
        inputs[dataset] = dict(path=str(path), size=st.st_size, mtimeNs=st.st_mtime_ns,
                               inode=st.st_ino, device=st.st_dev, mainconf=str(ROOT / value['mainconf']))
        for cent in cfg['centralities']:
            for condition in conditions:
                pattern = ('hKstar{side}VsPairYVsCent_' if condition['dimension']==3 else 'hKstar{side}VsCent_') + cfg['channel'] + condition['suffix']
                fields = [dataset,condition['id'],cent['id'],str(cent['min']),str(cent['max']),
                          str(condition['dimension']),str(condition['first']),str(condition['last']),
                          pattern.format(side='SE'),pattern.format(side='ME'),str(path)]
                assert all('\t' not in f and '\n' not in f for f in fields)
                rows.append('\t'.join(fields))
    (HERE / 'projection_plan.tsv').write_text('\n'.join(rows)+'\n')
    manifest = dict(preparedAtJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
                    conditions=conditions, inputs=inputs, planRows=len(rows), commonJobs=state['commonJobs'],
                    sourceFiles={cfg['recovery_config']:sha(ROOT/cfg['recovery_config']),
                                 cfg['recovery_state']:sha(ROOT/cfg['recovery_state']),
                                 'study.yaml':sha(HERE/'study.yaml')},
                    interpretation='EP changes the mixing pool. Rapidity is a common SE/ME pair selection; the two rapidity axes are not stored jointly.')
    (HERE / 'input_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('Prepared',len(rows),'projections;',len(conditions),'conditions per EP/centrality')


if __name__ == '__main__':
    main()
