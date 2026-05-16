#!/usr/bin/env python3
"""
ARFF GO Term Analysis Script

This script analyzes two ARFF files to compare GO term distributions between
lethal and viable genes, providing comprehensive statistics and comparisons.
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
from collections import Counter
import logging
import os

from ..constants import configure_logger, PhenGO_VERSION

logger = logging.getLogger(__name__)


def parse_arff_with_terms(file_path):
    """Parse ARFF file and extract genes with their GO term features."""
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
                if len(parts) < 3:  # Skip malformed lines
                    continue
                gene = parts[0]
                label = parts[-1]
                values = parts[1:-1]

                # Handle case where there are more values than expected attributes
                feature_dict = {}
                for i, term in enumerate(attributes[1:-1]):
                    if i < len(values):
                        feature_dict[term] = values[i]
                    else:
                        feature_dict[term] = 'NA'

                genes[gene] = {'label': label, 'features': feature_dict}

    return genes, attributes[1:-1]  # Exclude gene name and label columns


def calculate_go_term_stats(genes, go_terms):
    """Calculate GO term statistics for a set of genes."""
    stats = {}

    # Separate genes by label
    lethal_genes = {g: info for g, info in genes.items() if info['label'].lower() in ['lethal', 'essential']}
    viable_genes = {g: info for g, info in genes.items() if info['label'].lower() in ['viable', 'non-essential']}

    # Count GO term occurrences for each group
    lethal_go_counts = Counter()
    viable_go_counts = Counter()
    total_go_counts = Counter()

    for gene, info in genes.items():
        for go_term in go_terms:
            value = info['features'].get(go_term, 'NA')
            if value not in ['NA', '0', 0, '']:  # Consider non-zero, non-NA values as present
                total_go_counts[go_term] += 1
                if info['label'].lower() in ['lethal', 'essential']:
                    lethal_go_counts[go_term] += 1
                elif info['label'].lower() in ['viable', 'non-essential']:
                    viable_go_counts[go_term] += 1

    # Calculate enrichment ratios and statistics
    total_lethal = len(lethal_genes)
    total_viable = len(viable_genes)
    total_genes = len(genes)

    for go_term in go_terms:
        lethal_count = lethal_go_counts[go_term]
        viable_count = viable_go_counts[go_term]
        total_count = total_go_counts[go_term]

        # Calculate frequencies
        lethal_freq = lethal_count / total_lethal if total_lethal > 0 else 0
        viable_freq = viable_count / total_viable if total_viable > 0 else 0
        total_freq = total_count / total_genes if total_genes > 0 else 0

        # Calculate enrichment ratio (lethal vs viable)
        enrichment_ratio = lethal_freq / viable_freq if viable_freq > 0 else float('inf') if lethal_freq > 0 else 0

        # Calculate chi-square-like statistic for significance
        # Using simple odds ratio calculation
        a = lethal_count  # lethal with GO term
        b = total_lethal - lethal_count  # lethal without GO term
        c = viable_count  # viable with GO term
        d = total_viable - viable_count  # viable without GO term

        # Fisher's exact test approximation (odds ratio)
        if b > 0 and c > 0 and d > 0:
            odds_ratio = (a * d) / (b * c)
        else:
            odds_ratio = float('inf') if a > 0 and (c == 0 or b == 0) else 0

        stats[go_term] = {
            'total_count': total_count,
            'total_frequency': total_freq,
            'lethal_count': lethal_count,
            'lethal_frequency': lethal_freq,
            'viable_count': viable_count,
            'viable_frequency': viable_freq,
            'enrichment_ratio': enrichment_ratio,
            'odds_ratio': odds_ratio,
            'lethal_vs_total_ratio': lethal_count / total_count if total_count > 0 else 0
        }

    return stats, {
        'total_genes': total_genes,
        'lethal_genes': total_lethal,
        'viable_genes': total_viable,
        'other_genes': total_genes - total_lethal - total_viable
    }


def compare_go_distributions(stats_a, stats_b, go_terms, summary_a, summary_b):
    """Compare GO term distributions between two datasets."""
    comparisons = {}

    for go_term in go_terms:
        stat_a = stats_a.get(go_term, {})
        stat_b = stats_b.get(go_term, {})

        comparisons[go_term] = {
            # Dataset A stats
            'a_total_count': stat_a.get('total_count', 0),
            'a_total_freq': stat_a.get('total_frequency', 0),
            'a_lethal_count': stat_a.get('lethal_count', 0),
            'a_lethal_freq': stat_a.get('lethal_frequency', 0),
            'a_viable_count': stat_a.get('viable_count', 0),
            'a_viable_freq': stat_a.get('viable_frequency', 0),
            'a_enrichment_ratio': stat_a.get('enrichment_ratio', 0),

            # Dataset B stats
            'b_total_count': stat_b.get('total_count', 0),
            'b_total_freq': stat_b.get('total_frequency', 0),
            'b_lethal_count': stat_b.get('lethal_count', 0),
            'b_lethal_freq': stat_b.get('lethal_frequency', 0),
            'b_viable_count': stat_b.get('viable_count', 0),
            'b_viable_freq': stat_b.get('viable_frequency', 0),
            'b_enrichment_ratio': stat_b.get('enrichment_ratio', 0),

            # Comparisons
            'freq_difference': stat_a.get('total_frequency', 0) - stat_b.get('total_frequency', 0),
            'lethal_freq_difference': stat_a.get('lethal_frequency', 0) - stat_b.get('lethal_frequency', 0),
            'viable_freq_difference': stat_a.get('viable_frequency', 0) - stat_b.get('viable_frequency', 0),
            'enrichment_ratio_difference': stat_a.get('enrichment_ratio', 0) - stat_b.get('enrichment_ratio', 0)
        }

    return comparisons


def write_summary_statistics(output_file, summary_a, summary_b, file_a_name, file_b_name):
    """Write summary statistics to CSV."""
    summary_file = output_file.replace('.csv', '_summary.csv')

    with open(summary_file, 'w', newline='') as csvfile:
        fieldnames = ['Metric', 'Dataset_A', 'Dataset_B', 'Difference', 'Percent_Change']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        metrics = [
            ('Total Genes', 'total_genes'),
            ('Lethal Genes', 'lethal_genes'),
            ('Viable Genes', 'viable_genes'),
            ('Other Genes', 'other_genes')
        ]

        for metric_name, metric_key in metrics:
            val_a = summary_a[metric_key]
            val_b = summary_b[metric_key]
            diff = val_a - val_b
            pct_change = ((val_a - val_b) / val_b * 100) if val_b > 0 else 0

            writer.writerow({
                'Metric': metric_name,
                'Dataset_A': val_a,
                'Dataset_B': val_b,
                'Difference': diff,
                'Percent_Change': f"{pct_change:.2f}%"
            })

        # Add ratio statistics
        lethal_ratio_a = summary_a['lethal_genes'] / summary_a['total_genes'] if summary_a['total_genes'] > 0 else 0
        lethal_ratio_b = summary_b['lethal_genes'] / summary_b['total_genes'] if summary_b['total_genes'] > 0 else 0

        writer.writerow({
            'Metric': 'Lethal Gene Ratio',
            'Dataset_A': f"{lethal_ratio_a:.3f}",
            'Dataset_B': f"{lethal_ratio_b:.3f}",
            'Difference': f"{lethal_ratio_a - lethal_ratio_b:.3f}",
            'Percent_Change': f"{((lethal_ratio_a - lethal_ratio_b) / lethal_ratio_b * 100) if lethal_ratio_b > 0 else 0:.2f}%"
        })

    logger.info(f"Summary statistics written to: {summary_file}")



def write_top_go_terms(output_file, stats_a, stats_b, file_a_name, file_b_name, top_n=50):
    """Write top GO terms by different criteria."""
    top_terms_file = output_file.replace('.csv', '_top_terms.csv')

    with open(top_terms_file, 'w', newline='') as csvfile:
        fieldnames = [
            'GO_Term', 'Category', 'Rank',
            'A_Total_Count', 'A_Total_Freq', 'A_Lethal_Freq', 'A_Viable_Freq', 'A_Enrichment',
            'B_Total_Count', 'B_Total_Freq', 'B_Lethal_Freq', 'B_Viable_Freq', 'B_Enrichment'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Categories to analyze
        categories = [
            ('Most_Common_Overall_A', lambda x: stats_a[x]['total_frequency'], stats_a),
            ('Most_Common_Overall_B', lambda x: stats_b[x]['total_frequency'], stats_b),
            ('Most_Lethal_Enriched_A', lambda x: stats_a[x]['enrichment_ratio'], stats_a),
            ('Most_Lethal_Enriched_B', lambda x: stats_b[x]['enrichment_ratio'], stats_b),
            ('Most_Viable_Enriched_A',
             lambda x: 1 / stats_a[x]['enrichment_ratio'] if stats_a[x]['enrichment_ratio'] > 0 else 0, stats_a),
            ('Most_Viable_Enriched_B',
             lambda x: 1 / stats_b[x]['enrichment_ratio'] if stats_b[x]['enrichment_ratio'] > 0 else 0, stats_b)
        ]

        for category_name, sort_key, stats_dict in categories:
            # Get top terms for this category
            sorted_terms = sorted(stats_dict.keys(), key=sort_key, reverse=True)[:top_n]

            for rank, go_term in enumerate(sorted_terms, 1):
                stat_a = stats_a.get(go_term, {})
                stat_b = stats_b.get(go_term, {})

                writer.writerow({
                    'GO_Term': go_term,
                    'Category': category_name,
                    'Rank': rank,
                    'A_Total_Count': stat_a.get('total_count', 0),
                    'A_Total_Freq': f"{stat_a.get('total_frequency', 0):.4f}",
                    'A_Lethal_Freq': f"{stat_a.get('lethal_frequency', 0):.4f}",
                    'A_Viable_Freq': f"{stat_a.get('viable_frequency', 0):.4f}",
                    'A_Enrichment': f"{stat_a.get('enrichment_ratio', 0):.4f}",
                    'B_Total_Count': stat_b.get('total_count', 0),
                    'B_Total_Freq': f"{stat_b.get('total_frequency', 0):.4f}",
                    'B_Lethal_Freq': f"{stat_b.get('lethal_frequency', 0):.4f}",
                    'B_Viable_Freq': f"{stat_b.get('viable_frequency', 0):.4f}",
                    'B_Enrichment': f"{stat_b.get('enrichment_ratio', 0):.4f}"
                })

    logger.info(f"Top GO terms analysis written to: {top_terms_file}")


def main():
    parser = argparse.ArgumentParser(
        description=f"PhenoGO {PhenGO_VERSION} - ARFF GO Term Analysis: Analyze GO term distributions in ARFF files"
    )
    parser.add_argument("-arff_a", dest="arff_a", required=True,
                        help="First ARFF file for comparison")
    parser.add_argument("-arff_b", dest="arff_b", required=True,
                        help="Second ARFF file for comparison")
    parser.add_argument("-o", dest="output", required=True,
                        help="Output CSV file prefix (multiple files will be generated)")
    parser.add_argument("-output_dir", dest="output_dir", required=False, default='.',
                        help="Directory to write output files (default: current directory)")
    parser.add_argument("--top-n", dest="top_n", type=int, default=50,
                        help="Number of top GO terms to report for each category (default: 50)")

    args = parser.parse_args()

    # Prepare output directory and logging
    args.output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    logger = configure_logger('PhenGO.GO_Compare', enable_file=True, log_dir=args.output_dir, logfile_name='GO_Compare.log', level=logging.INFO)

    logger.info("Loading ARFF files...")
    genes_a, terms_a = parse_arff_with_terms(args.arff_a)
    genes_b, terms_b = parse_arff_with_terms(args.arff_b)

    logger.info(f"Dataset A: {len(genes_a)} genes, {len(terms_a)} GO terms")
    logger.info(f"Dataset B: {len(genes_b)} genes, {len(terms_b)} GO terms")

    # Get union of all GO terms
    all_terms = sorted(set(terms_a).union(set(terms_b)))
    logger.info(f"Total unique GO terms: {len(all_terms)}")

    # Calculate statistics for both datasets
    logger.info("Calculating GO term statistics...")
    stats_a, summary_a = calculate_go_term_stats(genes_a, all_terms)
    stats_b, summary_b = calculate_go_term_stats(genes_b, all_terms)

    # Compare distributions
    logger.info("Comparing GO term distributions...")
    comparisons = compare_go_distributions(stats_a, stats_b, all_terms, summary_a, summary_b)

    # Write detailed comparison
    logger.info("Writing detailed comparison...")
    output_path = os.path.join(args.output_dir, args.output)
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = [
            'GO_Term',
            'A_Total_Count', 'A_Total_Freq', 'A_Lethal_Count', 'A_Lethal_Freq',
            'A_Viable_Count', 'A_Viable_Freq', 'A_Enrichment_Ratio',
            'B_Total_Count', 'B_Total_Freq', 'B_Lethal_Count', 'B_Lethal_Freq',
            'B_Viable_Count', 'B_Viable_Freq', 'B_Enrichment_Ratio',
            'Total_Freq_Diff', 'Lethal_Freq_Diff', 'Viable_Freq_Diff', 'Enrichment_Diff'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Sort by total frequency difference (most changed terms first)
        sorted_terms = sorted(all_terms, key=lambda x: abs(comparisons[x]['freq_difference']), reverse=True)

        for go_term in sorted_terms:
            comp = comparisons[go_term]
            writer.writerow({
                'GO_Term': go_term,
                'A_Total_Count': comp['a_total_count'],
                'A_Total_Freq': f"{comp['a_total_freq']:.4f}",
                'A_Lethal_Count': comp['a_lethal_count'],
                'A_Lethal_Freq': f"{comp['a_lethal_freq']:.4f}",
                'A_Viable_Count': comp['a_viable_count'],
                'A_Viable_Freq': f"{comp['a_viable_freq']:.4f}",
                'A_Enrichment_Ratio': f"{comp['a_enrichment_ratio']:.4f}",
                'B_Total_Count': comp['b_total_count'],
                'B_Total_Freq': f"{comp['b_total_freq']:.4f}",
                'B_Lethal_Count': comp['b_lethal_count'],
                'B_Lethal_Freq': f"{comp['b_lethal_freq']:.4f}",
                'B_Viable_Count': comp['b_viable_count'],
                'B_Viable_Freq': f"{comp['b_viable_freq']:.4f}",
                'B_Enrichment_Ratio': f"{comp['b_enrichment_ratio']:.4f}",
                'Total_Freq_Diff': f"{comp['freq_difference']:.4f}",
                'Lethal_Freq_Diff': f"{comp['lethal_freq_difference']:.4f}",
                'Viable_Freq_Diff': f"{comp['viable_freq_difference']:.4f}",
                'Enrichment_Diff': f"{comp['enrichment_ratio_difference']:.4f}"
            })

    # Write additional analysis files
    write_summary_statistics(output_path, summary_a, summary_b, args.arff_a, args.arff_b)
    write_top_go_terms(output_path, stats_a, stats_b, args.arff_a, args.arff_b, args.top_n)

    logger.info("\n✅ GO term analysis complete!")
    logger.info(f"Main comparison written to: {output_path}")
    logger.info(f"Additional files generated with prefix: {output_path.replace('.csv', '')}")

    # Print quick summary
    logger.info(f"\nQuick Summary:")
    logger.info(f"Dataset A - Total: {summary_a['total_genes']}, Lethal: {summary_a['lethal_genes']}, Viable: {summary_a['viable_genes']}")
    logger.info(f"Dataset B - Total: {summary_b['total_genes']}, Lethal: {summary_b['lethal_genes']}, Viable: {summary_b['viable_genes']}")


if __name__ == "__main__":
    main()