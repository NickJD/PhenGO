#!/usr/bin/env python3
"""GO term mapper: map GO IDs from an OBO file to name/namespace/def/is_a

Simple, dependency-light OBO parser and mapper.
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

import argparse
import os
import logging
from ..constants import configure_logger


def parse_obo(path):
    terms = {}
    cur = None
    with open(path, 'r', encoding='utf-8') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if line == '[Term]':
                if cur and 'id' in cur:
                    terms[cur['id']] = cur
                cur = {}
                continue
            if not cur:
                continue
            if not line.strip():
                continue
            if ':' not in line:
                continue
            k, v = line.split(':', 1)
            k = k.strip()
            v = v.strip()
            if k in cur:
                # store multiples as list
                if isinstance(cur[k], list):
                    cur[k].append(v)
                else:
                    cur[k] = [cur[k], v]
            else:
                cur[k] = v
    if cur and 'id' in cur:
        terms[cur['id']] = cur
    return terms


def map_go_ids(go_ids, obo_terms):
    results = []
    for gid in go_ids:
        info = obo_terms.get(gid)
        if info:
            results.append((gid, info.get('name', ''), info.get('namespace', ''), info.get('def', ''), info.get('is_a', '')))
        else:
            results.append((gid, '', '', '', ''))
    return results


def main():
    parser = argparse.ArgumentParser(description='Map GO IDs using a local OBO file')
    parser.add_argument('-i', '--input', required=True, help='File with one GO ID per line')
    parser.add_argument('-o', '--output', default='go_mapped.tsv')
    parser.add_argument('--obo', required=True, help='Path to go.obo')
    args = parser.parse_args()

    configure_logger('PhenGO.go_mapper', enable_file=False)
    logger = logging.getLogger('PhenGO.go_mapper')

    if not os.path.exists(args.obo):
        logger.error(f"OBO file not found: {args.obo}")
        raise SystemExit(2)

    obo_terms = parse_obo(args.obo)
    logger.info(f"Parsed {len(obo_terms)} terms from OBO")

    go_ids = []
    with open(args.input, 'r', encoding='utf-8') as fh:
        for raw in fh:
            gid = raw.strip()
            if gid:
                go_ids.append(gid)

    mapped = map_go_ids(go_ids, obo_terms)
    with open(args.output, 'w', encoding='utf-8') as fh:
        fh.write('GO_ID\tName\tNamespace\tDef\tIs_A\n')
        for row in mapped:
            fh.write('\t'.join([str(x).replace('\n',' ') for x in row]) + '\n')

    logger.info(f"Wrote mapping to {args.output}")


if __name__ == '__main__':
    main()

