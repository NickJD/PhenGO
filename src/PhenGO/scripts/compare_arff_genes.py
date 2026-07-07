"""
PhenGO Compare-ARFF — compare two ARFF files gene-by-gene.

Usage
-----
    compare-arff -arff_a reference.arff -arff_b comparison.arff -o results.csv
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
import csv
from collections import defaultdict
import logging
from ..constants import configure_logger, PhenGO_VERSION



def parse_arff_with_terms(file_path):

    data_started = False
    genes = {}
    attributes = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            if line.lower().startswith('@attribute'):
                parts = line.split()
                attr_name = parts[1].strip("'\"")  # Clean quotes
                attributes.append(attr_name)
            elif line.lower() == '@data':
                data_started = True
            elif data_started:
                parts = [p.strip().strip('"\'') for p in next(csv.reader([line]))]
                gene = parts[0]
                label = parts[-1]
                values = parts[1:-1]
                feature_dict = {term: val for term, val in zip(attributes[1:-1], values)}
                genes[gene] = {'label': label, 'features': feature_dict}
    return genes, attributes[1:-1]

def compare_genes(genes_a, genes_b, all_terms):
    grouped = defaultdict(list)

    for gene, a_info in genes_a.items():
        row = {
            'Gene': gene,
            'Label A': a_info['label'],
            'Label B': '',
            'GO Terms Differ': '',
            'Status': ''
        }

        statuses = []

        if gene not in genes_b:
            statuses.append('MISSING_IN_B')
        else:
            b_info = genes_b[gene]
            row['Label B'] = b_info['label']

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

        # Add row to each status group for grouping
        for status in statuses:
            grouped[status].append(row)

    return grouped

def main():
    parser = argparse.ArgumentParser(description=f"PhenoGO {PhenGO_VERSION} - Compare-ARFF: Compare two ARFF files.")
    parser.add_argument("-arff_a", dest="arff_a", required=True, help="Master ARFF file (reference)")
    parser.add_argument("-arff_b", dest="arff_b", required=True, help="Comparison ARFF file")
    parser.add_argument("-o", dest="output", required=True, help="Output CSV file")

    args = parser.parse_args()

    genes_a, terms_a = parse_arff_with_terms(args.arff_a)
    genes_b, terms_b = parse_arff_with_terms(args.arff_b)

    all_terms = sorted(set(terms_a).union(set(terms_b)))
    grouped_results = compare_genes(genes_a, genes_b, all_terms)

    # Define desired order of groups
    configure_logger('PhenGO.compare_arff_genes', enable_file=False)
    logger = logging.getLogger('PhenGO.compare_arff_genes')
    status_order = ['MISSING_IN_B', 'LABEL_MISMATCH', 'GO_TERM_MISMATCH', 'EXACT_MATCH']

    # Write structured output grouped by Status
    with open(args.output, 'w', newline='') as csvfile:
        fieldnames = ['Gene', 'Label A', 'Label B', 'GO Terms Differ', 'Status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for status in status_order:
            group = sorted(grouped_results.get(status, []), key=lambda x: x['Gene'])
            for row in group:
                writer.writerow(row)

    logger.info(f"\n✅ Grouped comparison complete. Output written to: {args.output}")

    # Print summary
    logger.info("\nSummary of differences:")
    for status in status_order:
        logger.info(f"  {status:17}: {len(grouped_results.get(status, []))}")

if __name__ == "__main__":
    main()
