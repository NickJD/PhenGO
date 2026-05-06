#!/usr/bin/env python3
"""ARFF validator for PhenGO

Produces a JSON and text summary describing basic issues in ARFF files.
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
import json
from collections import Counter
import logging
from ..constants import configure_logger


def parse_arff(path):
    """Simple ARFF parser that returns (attribute_names, attribute_types, data_lines).
    data_lines is list of lists (strings).
    """
    attr_names = []
    attr_types = []
    data_lines = []
    in_data = False
    with open(path, 'r', encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('%'):
                continue
            low = line.lower()
            if low.startswith('@relation'):
                continue
            if low.startswith('@attribute'):
                parts = line.split()
                if len(parts) >= 3:
                    name = parts[1].strip('"\'')
                    attr_names.append(name)
                    attr_types.append(' '.join(parts[2:]).strip())
                continue
            if low.startswith('@data'):
                in_data = True
                continue
            if in_data:
                # parse CSV-like respecting quoted strings
                vals = []
                cur = ''
                in_q = False
                for ch in line:
                    if ch == '"':
                        in_q = not in_q
                        cur += ch
                    elif ch == ',' and not in_q:
                        vals.append(cur.strip().strip('"'))
                        cur = ''
                    else:
                        cur += ch
                vals.append(cur.strip().strip('"'))
                data_lines.append(vals)
    return attr_names, attr_types, data_lines


def validate(path):
    logger = logging.getLogger('PhenGO.arff_validator')
    logger.info(f"Validating {path}")
    attr_names, attr_types, data_lines = parse_arff(path)
    n_attrs = len(attr_names)
    n_rows = len(data_lines)
    mismatch_rows = [i for i, row in enumerate(data_lines, 1) if len(row) != n_attrs]

    # missingness per column
    missing = [0] * n_attrs
    for row in data_lines:
        for i in range(n_attrs):
            try:
                v = row[i]
            except Exception:
                v = ''
            if v == '' or v.lower() == 'na' or v == '?':
                missing[i] += 1

    missing_pct = {attr_names[i]: (missing[i], round(missing[i] / max(1, n_rows) * 100, 2))
                   for i in range(n_attrs)}

    # detect numeric-like columns
    numeric_counts = [0] * n_attrs
    for row in data_lines:
        for i in range(n_attrs):
            try:
                v = row[i]
            except Exception:
                v = ''
            try:
                if v is None or v == '' or v == '?':
                    continue
                float(v)
                numeric_counts[i] += 1
            except Exception:
                pass

    numeric_pct = {attr_names[i]: round(numeric_counts[i] / max(1, n_rows) * 100, 2)
                   for i in range(n_attrs)}

    # class distribution (assume last column is class-like)
    class_counts = {}
    if n_rows > 0 and n_attrs > 0:
        lastcol = [row[-1] if len(row) >= 1 else '' for row in data_lines]
        class_counts = dict(Counter(lastcol))

    summary = {
        'file': path,
        'n_attributes': n_attrs,
        'n_rows': n_rows,
        'mismatch_rows_count': len(mismatch_rows),
        'mismatch_row_indices': mismatch_rows[:20],
        'missing_per_column': missing_pct,
        'numeric_pct_per_column': numeric_pct,
        'class_counts': class_counts
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description='PhenGO ARFF validator')
    parser.add_argument('-i', '--input', required=True, help='ARFF file')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    parser.add_argument('--fail-on-error', action='store_true')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    configure_logger('PhenGO.arff_validator', enable_file=True, log_dir=args.output_dir,
                     logfile_name='arff_validator.log')
    logger = logging.getLogger('PhenGO.arff_validator')

    summary = validate(args.input)
    out_json = os.path.join(args.output_dir, 'arff_summary.json')
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2)
    logger.info(f"Wrote {out_json}")

    out_txt = os.path.join(args.output_dir, 'arff_summary.txt')
    with open(out_txt, 'w', encoding='utf-8') as fh:
        fh.write(json.dumps(summary, indent=2))
    logger.info(f"Wrote {out_txt}")

    # if fail on error and mismatches exist, exit 2
    if args.fail_on_error and summary['mismatch_rows_count'] > 0:
        logger.error("Validation failed: mismatched rows present")
        raise SystemExit(2)


if __name__ == '__main__':
    main()

