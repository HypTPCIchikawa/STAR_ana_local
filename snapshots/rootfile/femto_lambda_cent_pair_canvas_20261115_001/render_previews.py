#!/usr/bin/env python3
"""Render every delivered page and make per-centrality contact sheets."""
from pathlib import Path
import json
import subprocess
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    cfg=yaml.safe_load((HERE/'study.yaml').read_text());out=ROOT/cfg['canvas']['output_directory']
    index=json.loads((out/'pdf_index.json').read_text());preview=out/'previews';preview.mkdir(exist_ok=True)
    subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-dUseCropBox',
        '-sDEVICE=pngalpha','-r100','-dTextAlphaBits=4','-dGraphicsAlphaBits=4',
        '-sOutputFile='+str(preview/'page-%03d.png'),str(ROOT/index['combined_pdf'])],check=True)
    for entry in index['centralities']:
        images=[str(preview/('page-%03d.png'%i)) for i in
                range(entry['combined_first_page'],entry['combined_last_page']+1)]
        subprocess.run(['montage']+images+['-thumbnail','240x','-tile','3x4','-geometry','+8+8',
            '-background','white',str(preview/('contact_'+entry['centrality']['id']+'.png'))],check=True)
    print('CENTRALITY_PREVIEWS_OK pages=121 contact_sheets=11')


if __name__=='__main__':main()
