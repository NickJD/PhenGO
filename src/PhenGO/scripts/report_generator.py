#!/usr/bin/env python3
"""Simple HTML summariser for PhenGO result directories.

Given an output directory with CSVs and PNGs, generate a small self-contained
HTML index that previews tables and images.
"""
import os as _os, sys as _sys
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here   = _os.path.dirname(_os.path.abspath(__file__))
    _src    = _os.path.normpath(_os.path.join(_here, '..', '..'))
    _phengo = _os.path.dirname(_here)
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',         _os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.scripts', _os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in _sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            _sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.scripts'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo

import os
import argparse
import base64
import logging
import pandas as pd
from ..constants import configure_logger


def embed_image(path):
    try:
        with open(path, 'rb') as fh:
            data = fh.read()
        b64 = base64.b64encode(data).decode('ascii')
        ext = os.path.splitext(path)[1].lstrip('.')
        return f'data:image/{ext};base64,{b64}'
    except Exception:
        return None


def generate_index(input_dir, out_path):
    logger = logging.getLogger('PhenGO.report_generator')
    logger.info(f"Generating HTML report for {input_dir}")
    entries = os.listdir(input_dir)
    imgs = [e for e in entries if e.lower().endswith(('.png', '.jpg', '.jpeg', '.svg'))]
    csvs = [e for e in entries if e.lower().endswith('.csv')]

    html_parts = ["<html><head><meta charset='utf-8'><title>PhenGO Report</title></head><body>"]
    html_parts.append(f"<h1>PhenGO report: {os.path.basename(input_dir)}</h1>")

    if imgs:
        html_parts.append('<h2>Plots</h2>')
        for img in imgs:
            full = os.path.join(input_dir, img)
            datauri = embed_image(full)
            if datauri:
                html_parts.append(f"<div><h3>{img}</h3><img src=\"{datauri}\" style=\"max-width:900px;\"></div>")
            else:
                html_parts.append(f"<div><a href=\"{img}\">{img}</a></div>")

    if csvs:
        html_parts.append('<h2>Tables (preview)</h2>')
        for csv in csvs:
            full = os.path.join(input_dir, csv)
            try:
                df = pd.read_csv(full)
                html_parts.append(f"<h3>{csv}</h3>")
                html_parts.append(df.head(20).to_html(index=False, classes='table'))
            except Exception as e:
                logger.warning(f"Could not read CSV {csv}: {e}")
                html_parts.append(f"<div><a href=\"{csv}\">{csv}</a> (could not preview)</div>")

    html_parts.append('</body></html>')
    html = '\n'.join(html_parts)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(html)
    logger.info(f"Wrote HTML report to {out_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate a simple HTML report from a results directory')
    parser.add_argument('-i', '--input-dir', required=True)
    parser.add_argument('-o', '--output', default='report.html')
    args = parser.parse_args()

    configure_logger('PhenGO.report_generator', enable_file=False)
    generate_index(args.input_dir, args.output)


if __name__ == '__main__':
    main()

