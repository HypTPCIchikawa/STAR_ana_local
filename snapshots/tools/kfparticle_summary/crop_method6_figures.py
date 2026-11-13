"""Set a vector PDF CropBox around visible marks without altering plot data."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import subprocess

from build_report import pdf_pages


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--margin-points', type=float, required=True)
    args = parser.parse_args()
    source = args.input_dir.resolve(strict=True)
    out = args.output_dir.resolve(strict=True)
    if source == out or args.margin_points < 0:
        raise ValueError('Distinct directories and a nonnegative margin are required.')
    paths = sorted(source.glob('*.pdf'))
    if not paths or any((out/path.name).exists() for path in paths):
        raise ValueError('PDF inputs and new output filenames are required.')

    def crop(path):
        target = out/path.name
        before = sha(path)
        bbox = subprocess.run(['gs', '-q', '-dSAFER', '-dBATCH', '-dNOPAUSE',
                               '-sDEVICE=bbox', str(path)], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True, check=True)
        boxes = [list(map(float, row.split())) for row in
                 re.findall(r'%%HiResBoundingBox: ([^\r\n]+)', bbox.stderr)]
        if not boxes or any(len(row) != 4 or row[2] <= row[0] or row[3] <= row[1] for row in boxes):
            raise ValueError('Missing or empty visible-page bounds: '+str(path))
        margin = args.margin_points
        bounds = [min(b[0] for b in boxes)-margin, min(b[1] for b in boxes)-margin,
                  max(b[2] for b in boxes)+margin, max(b[3] for b in boxes)+margin]
        mark = '[/CropBox ['+' '.join(format(v, '.6f') for v in bounds)+'] /PAGES pdfmark'
        result = subprocess.run(['gs', '-q', '-dSAFER', '-dBATCH', '-dNOPAUSE',
                                 '-sDEVICE=pdfwrite', '-dAutoRotatePages=/None',
                                 '-sOutputFile='+str(target), '-c', mark, '-f', str(path)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        pages = pdf_pages(target)
        if pages != len(boxes) or sha(path) != before:
            raise ValueError('Page count or raw PDF changed while cropping.')
        return path.name, dict(source=str(path), sourceSha256=before,
                               outputSha256=sha(target), pages=pages,
                               cropBox=bounds, bytes=target.stat().st_size,
                               writerDiagnostics=(result.stdout+result.stderr).strip())

    with ThreadPoolExecutor(max_workers=2) as pool:
        records = dict(pool.map(crop, paths))
    (out/'crop_manifest.json').write_text(json.dumps(dict(
        phase='complete', operation='Vector CropBox; no fit, bin or numerical changes',
        marginPoints=args.margin_points, files=records), indent=2)+'\n')
    print(json.dumps(dict(pdfs=len(records), pages=sum(v['pages'] for v in records.values()))))


if __name__ == '__main__':
    main()
