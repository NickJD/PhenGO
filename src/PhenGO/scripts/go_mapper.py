#!/usr/bin/env python3
"""Map GO IDs to release-matched OBO metadata without guessing replacements."""
import os as _os
import sys as _sys
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
import gzip
import os
import logging
from ..constants import configure_logger


class OBOTerms(dict):
    """Canonical OBO records plus a validated alternate-ID index."""

    def __init__(self):
        super().__init__()
        self.alt_id_to_id = {}


def _open_obo(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    return opener(path, 'rt', encoding='utf-8', errors='replace')


def parse_obo(path):
    terms = OBOTerms()
    cur = None
    with _open_obo(path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if line == '[Term]':
                if cur and 'id' in cur:
                    _store_term(terms, cur)
                cur = {}
                continue
            if line.startswith('['):
                if cur and 'id' in cur:
                    _store_term(terms, cur)
                cur = None
                continue
            if cur is None:
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
        _store_term(terms, cur)

    for term_id, record in terms.items():
        for alt_id in _as_list(record.get('alt_id')):
            previous = terms.alt_id_to_id.get(alt_id)
            if previous is not None and previous != term_id:
                raise ValueError(
                    f"Alternate GO ID {alt_id} maps to both {previous} and {term_id}"
                )
            terms.alt_id_to_id[alt_id] = term_id
    return terms


def _store_term(terms, record):
    term_id = record['id']
    if term_id in terms:
        raise ValueError(f"Duplicate GO term stanza: {term_id}")
    terms[term_id] = record


def _as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _term_id(value):
    return str(value).split()[0]


def map_go_ids(go_ids, obo_terms):
    results = []
    alt_map = getattr(obo_terms, 'alt_id_to_id', {})
    for requested_id in go_ids:
        canonical_id = alt_map.get(requested_id, requested_id)
        used_alt = canonical_id != requested_id
        info = obo_terms.get(canonical_id)
        if info is None:
            results.append((requested_id, '', 'missing', '', '', '', ''))
            continue

        is_obsolete = str(info.get('is_obsolete', 'false')).lower() == 'true'
        if is_obsolete:
            replacements = {
                alt_map.get(_term_id(value), _term_id(value))
                for value in _as_list(info.get('replaced_by'))
            }
            replacements = {
                replacement for replacement in replacements
                if replacement in obo_terms and str(
                    obo_terms[replacement].get('is_obsolete', 'false')
                ).lower() != 'true'
            }
            if len(replacements) != 1:
                results.append((
                    requested_id, '', 'obsolete_unresolved', '', '', '', '',
                ))
                continue
            canonical_id = next(iter(replacements))
            info = obo_terms[canonical_id]
            status = 'alt_id_obsolete_replaced' if used_alt else 'obsolete_replaced'
        else:
            status = 'alt_id' if used_alt else 'exact'

        parents = '|'.join(
            _term_id(value) for value in _as_list(info.get('is_a'))
        )
        results.append((
            requested_id,
            canonical_id,
            status,
            info.get('name', ''),
            info.get('namespace', ''),
            info.get('def', ''),
            parents,
        ))
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
    logger.info(
        "Parsed %d canonical terms and %d alternate IDs from OBO",
        len(obo_terms), len(obo_terms.alt_id_to_id),
    )

    go_ids = []
    with open(args.input, 'r', encoding='utf-8') as fh:
        for raw in fh:
            gid = raw.strip()
            if gid:
                go_ids.append(gid)

    mapped = map_go_ids(go_ids, obo_terms)
    with open(args.output, 'w', encoding='utf-8') as fh:
        fh.write(
            'Requested_GO_ID\tResolved_GO_ID\tMapping_Status\tName\t'
            'Namespace\tDef\tIs_A\n'
        )
        for row in mapped:
            fh.write('\t'.join([str(x).replace('\n',' ') for x in row]) + '\n')

    logger.info(f"Wrote mapping to {args.output}")


if __name__ == '__main__':
    main()
