import argparse
import csv
from collections import defaultdict

try: # While calling this script through pip
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *



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
                parts = [p.strip() for p in line.split(',')]
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

    print(f"\n✅ Grouped comparison complete. Output written to: {args.output}")

    # Print summary
    print("\nSummary of differences:")
    for status in status_order:
        print(f"  {status:17}: {len(grouped_results.get(status, []))}")

if __name__ == "__main__":
    main()
