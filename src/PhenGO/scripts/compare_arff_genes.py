"""
PhenGO Compare-ARFF — compare two ARFF files gene-by-gene.

Usage
-----
    compare-arff -arff_a reference.arff -arff_b comparison.arff -o results.csv
"""
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
import csv
import logging
from ..constants import configure_logger, PhenGO_VERSION
from ..predict.data import load_arff_data



def parse_arff_with_terms(file_path):
    frame, _ = load_arff_data(file_path)
    if frame is None or frame.shape[1] < 3:
        raise ValueError(f"Could not load usable ARFF data from {file_path}")
    terms = list(frame.columns[1:-1])
    genes = {}
    for _, row in frame.iterrows():
        gene = str(row.iloc[0])
        record = {
            'label': str(row.iloc[-1]),
            'features': {term: str(row[term]) for term in terms},
        }
        if gene in genes and genes[gene] != record:
            raise ValueError(f"Conflicting duplicate rows for gene {gene} in {file_path}")
        genes[gene] = record
    return genes, terms

def compare_genes(genes_a, genes_b, all_terms):
    rows = []

    for gene in sorted(set(genes_a) | set(genes_b)):
        a_info = genes_a.get(gene)
        b_info = genes_b.get(gene)
        row = {
            'Gene': gene,
            'Label A': a_info['label'] if a_info else '',
            'Label B': b_info['label'] if b_info else '',
            'GO Terms Differ': '',
            'Status': ''
        }

        statuses = []
        if a_info is None:
            statuses.append('MISSING_IN_A')
        elif b_info is None:
            statuses.append('MISSING_IN_B')
        else:
            if a_info['label'] != b_info['label']:
                statuses.append('LABEL_MISMATCH')

            differing_terms = []
            for term in all_terms:
                va = a_info['features'].get(term, 'NA')
                vb = b_info['features'].get(term, 'NA')
                if va != vb:
                    differing_terms.append(term)

            if differing_terms:
                row['GO Terms Differ'] = ';'.join(differing_terms)
                statuses.append('GO_TERM_MISMATCH')

            if not statuses:
                statuses.append('EXACT_MATCH')

        row['Status'] = ';'.join(statuses)
        rows.append(row)
    return rows

def main():
    parser = argparse.ArgumentParser(description=f"PhenoGO {PhenGO_VERSION} - Compare-ARFF: Compare two ARFF files.")
    parser.add_argument("-arff_a", dest="arff_a", required=True, help="Master ARFF file (reference)")
    parser.add_argument("-arff_b", dest="arff_b", required=True, help="Comparison ARFF file")
    parser.add_argument("-o", dest="output", required=True, help="Output CSV file")
    parser.add_argument("-overwrite", action="store_true")

    args = parser.parse_args()

    missing = [path for path in (args.arff_a, args.arff_b) if not _os.path.isfile(path)]
    if missing:
        parser.error("ARFF file(s) not found: " + ", ".join(missing))
    output = _os.path.abspath(args.output)
    if output in {_os.path.abspath(args.arff_a), _os.path.abspath(args.arff_b)}:
        parser.error("Output must not overwrite an input ARFF")
    if _os.path.exists(output) and not args.overwrite:
        parser.error("Output exists; choose another path or pass -overwrite")
    parent = _os.path.dirname(output)
    if parent:
        _os.makedirs(parent, exist_ok=True)
    args.output = output

    try:
        genes_a, terms_a = parse_arff_with_terms(args.arff_a)
        genes_b, terms_b = parse_arff_with_terms(args.arff_b)
    except ValueError as exc:
        parser.error(str(exc))

    all_terms = sorted(set(terms_a).union(set(terms_b)))
    results = compare_genes(genes_a, genes_b, all_terms)

    # Define desired order of groups
    configure_logger('PhenGO.compare_arff_genes', enable_file=False)
    logger = logging.getLogger('PhenGO.compare_arff_genes')
    status_order = ['MISSING_IN_A', 'MISSING_IN_B', 'LABEL_MISMATCH', 'GO_TERM_MISMATCH', 'EXACT_MATCH']

    # Write structured output grouped by Status
    with open(args.output, 'w', newline='') as csvfile:
        fieldnames = ['Gene', 'Label A', 'Label B', 'GO Terms Differ', 'Status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        writer.writerows(results)

    logger.info(f"\n✅ Grouped comparison complete. Output written to: {args.output}")

    # Print summary
    logger.info("\nSummary of differences:")
    for status in status_order:
        count = sum(status in row['Status'].split(';') for row in results)
        logger.info(f"  {status:17}: {count}")

if __name__ == "__main__":
    main()
