#!/usr/bin/env python3
"""
ARFF GO Term Temporal Analysis Script

This script analyzes multiple ARFF files across time points to track:
1. Changes in top overrepresented GO terms for lethal/viable genes
2. Changes in gene lethal/viable classifications over time
3. GO term ranking stability and trends
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
from collections import defaultdict, Counter
import statistics
import logging
import os
import sys

try:
    from ..constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    try:
        from constants import *
    except:
        PhenGO_VERSION = "1.0"


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
                attr_name = parts[1].strip("'\"")
                attributes.append(attr_name)
            elif line.lower() == '@data':
                data_started = True
            elif data_started:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 3:
                    continue
                gene = parts[0]
                label = parts[-1]
                values = parts[1:-1]

                feature_dict = {}
                for i, term in enumerate(attributes[1:-1]):
                    if i < len(values):
                        feature_dict[term] = values[i]
                    else:
                        feature_dict[term] = 'NA'

                genes[gene] = {'label': label, 'features': feature_dict}

    return genes, attributes[1:-1]


def calculate_enrichment_stats(genes, go_terms):
    """Calculate GO term enrichment statistics."""
    lethal_genes = {g: info for g, info in genes.items()
                    if info['label'].lower() in ['lethal', 'essential']}
    viable_genes = {g: info for g, info in genes.items()
                    if info['label'].lower() in ['viable', 'non-essential']}

    lethal_go_counts = Counter()
    viable_go_counts = Counter()

    for gene, info in genes.items():
        for go_term in go_terms:
            value = info['features'].get(go_term, 'NA')
            if value not in ['NA', '0', 0, '']:
                if info['label'].lower() in ['lethal', 'essential']:
                    lethal_go_counts[go_term] += 1
                elif info['label'].lower() in ['viable', 'non-essential']:
                    viable_go_counts[go_term] += 1

    total_lethal = len(lethal_genes)
    total_viable = len(viable_genes)

    enrichment_stats = {}
    for go_term in go_terms:
        lethal_count = lethal_go_counts[go_term]
        viable_count = viable_go_counts[go_term]

        lethal_freq = lethal_count / total_lethal if total_lethal > 0 else 0
        viable_freq = viable_count / total_viable if total_viable > 0 else 0

        # Enrichment ratio (lethal vs viable)
        enrichment_ratio = lethal_freq / viable_freq if viable_freq > 0 else float('inf') if lethal_freq > 0 else 0

        # Depletion ratio (viable vs lethal) - inverse
        depletion_ratio = viable_freq / lethal_freq if lethal_freq > 0 else float('inf') if viable_freq > 0 else 0

        # Fisher's exact test approximation
        a = lethal_count
        b = total_lethal - lethal_count
        c = viable_count
        d = total_viable - viable_count

        if b > 0 and c > 0 and d > 0:
            odds_ratio = (a * d) / (b * c)
        else:
            odds_ratio = float('inf') if a > 0 and (c == 0 or b == 0) else 0

        enrichment_stats[go_term] = {
            'lethal_count': lethal_count,
            'lethal_freq': lethal_freq,
            'viable_count': viable_count,
            'viable_freq': viable_freq,
            'enrichment_ratio': enrichment_ratio,
            'depletion_ratio': depletion_ratio,
            'odds_ratio': odds_ratio,
            'total_genes_with_term': lethal_count + viable_count
        }

    return enrichment_stats, {
        'total_genes': len(genes),
        'lethal_genes': total_lethal,
        'viable_genes': total_viable
    }


def track_gene_classifications(datasets_dict):
    """Track how gene classifications change across time points."""
    # Get all unique genes across all time points
    all_genes = set()
    for genes_dict, _ in datasets_dict.values():
        all_genes.update(genes_dict.keys())

    # Track classification for each gene at each time point
    gene_timeline = {}
    for gene in all_genes:
        timeline = {}
        for timepoint, (genes_dict, _) in sorted(datasets_dict.items()):
            if gene in genes_dict:
                timeline[timepoint] = genes_dict[gene]['label'].lower()
            else:
                timeline[timepoint] = 'absent'
        gene_timeline[gene] = timeline

    # Identify genes that changed classification
    changing_genes = {}
    stable_genes = {}

    for gene, timeline in gene_timeline.items():
        classifications = [c for c in timeline.values() if c != 'absent']
        unique_classifications = set(classifications)

        # Only consider genes present in at least 2 timepoints
        if len(classifications) >= 2:
            if len(unique_classifications) > 1:
                changing_genes[gene] = timeline
            else:
                stable_genes[gene] = timeline

    return gene_timeline, changing_genes, stable_genes


def get_top_enriched_terms(enrichment_stats, n=20, metric='enrichment_ratio'):
    """Get top N GO terms by enrichment metric."""
    # Filter out infinite values for ranking
    filtered_terms = {
        term: stats for term, stats in enrichment_stats.items()
        if stats[metric] != float('inf') and stats[metric] > 0
    }

    sorted_terms = sorted(
        filtered_terms.items(),
        key=lambda x: x[1][metric],
        reverse=True
    )

    return sorted_terms[:n]


def write_temporal_top_terms(output_prefix, datasets_dict, top_n=20):
    """Write top GO terms over time for lethal and viable enrichment."""

    # Lethal-enriched terms over time
    lethal_file = f"{output_prefix}_lethal_enriched_timeline.csv"
    with open(lethal_file, 'w', newline='') as f:
        # Get all timepoints
        timepoints = sorted(datasets_dict.keys())

        # Build header
        header = ['Rank', 'GO_Term']
        for tp in timepoints:
            header.extend([f'{tp}_EnrichRatio', f'{tp}_LethalFreq', f'{tp}_ViableFreq'])

        writer = csv.writer(f)
        writer.writerow(header)

        # Get top terms for each timepoint
        all_top_terms = set()
        timepoint_rankings = {}

        for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
            top_terms = get_top_enriched_terms(enrichment_stats, top_n, 'enrichment_ratio')
            timepoint_rankings[tp] = {term: (rank + 1, stats)
                                      for rank, (term, stats) in enumerate(top_terms)}
            all_top_terms.update([term for term, _ in top_terms])

        # Write each GO term's trajectory
        for rank, go_term in enumerate(sorted(all_top_terms), 1):
            row = [rank, go_term]
            for tp in timepoints:
                if go_term in timepoint_rankings[tp]:
                    _, stats = timepoint_rankings[tp][go_term]
                    row.extend([
                        f"{stats['enrichment_ratio']:.3f}",
                        f"{stats['lethal_freq']:.4f}",
                        f"{stats['viable_freq']:.4f}"
                    ])
                else:
                    row.extend(['NA', 'NA', 'NA'])
            writer.writerow(row)

    logger.info(f"Lethal-enriched GO terms timeline written to: {lethal_file}")

    # Viable-enriched terms over time
    viable_file = f"{output_prefix}_viable_enriched_timeline.csv"
    with open(viable_file, 'w', newline='') as f:
        header = ['Rank', 'GO_Term']
        for tp in timepoints:
            header.extend([f'{tp}_DepletionRatio', f'{tp}_ViableFreq', f'{tp}_LethalFreq'])

        writer = csv.writer(f)
        writer.writerow(header)

        all_top_terms = set()
        timepoint_rankings = {}

        for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
            top_terms = get_top_enriched_terms(enrichment_stats, top_n, 'depletion_ratio')
            timepoint_rankings[tp] = {term: (rank + 1, stats)
                                      for rank, (term, stats) in enumerate(top_terms)}
            all_top_terms.update([term for term, _ in top_terms])

        for rank, go_term in enumerate(sorted(all_top_terms), 1):
            row = [rank, go_term]
            for tp in timepoints:
                if go_term in timepoint_rankings[tp]:
                    _, stats = timepoint_rankings[tp][go_term]
                    row.extend([
                        f"{stats['depletion_ratio']:.3f}",
                        f"{stats['viable_freq']:.4f}",
                        f"{stats['lethal_freq']:.4f}"
                    ])
                else:
                    row.extend(['NA', 'NA', 'NA'])
            writer.writerow(row)

    logger.info(f"Viable-enriched GO terms timeline written to: {viable_file}")


def write_rank_stability(output_prefix, datasets_dict, top_n=20):
    """Analyze and write GO term rank stability over time."""
    output_file = f"{output_prefix}_rank_stability.csv"

    timepoints = sorted(datasets_dict.keys())

    # Track rankings for each term over time
    term_rankings = defaultdict(dict)

    for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
        # Lethal enrichment rankings
        top_lethal = get_top_enriched_terms(enrichment_stats, top_n, 'enrichment_ratio')
        for rank, (term, stats) in enumerate(top_lethal, 1):
            if term not in term_rankings:
                term_rankings[term] = {'type': 'lethal_enriched', 'rankings': {}}
            term_rankings[term]['rankings'][tp] = rank

        # Viable enrichment rankings
        top_viable = get_top_enriched_terms(enrichment_stats, top_n, 'depletion_ratio')
        for rank, (term, stats) in enumerate(top_viable, 1):
            if term not in term_rankings:
                term_rankings[term] = {'type': 'viable_enriched', 'rankings': {}}
            elif term_rankings[term]['type'] != 'viable_enriched':
                term_rankings[term]['type'] = 'both'
            term_rankings[term]['rankings'][tp] = rank

    with open(output_file, 'w', newline='') as f:
        fieldnames = ['GO_Term', 'Enrichment_Type', 'First_Appearance', 'Last_Appearance',
                      'Times_In_Top20', 'Avg_Rank', 'Rank_StdDev', 'Rank_Range'] + \
                     [f'{tp}_Rank' for tp in timepoints]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for term, data in sorted(term_rankings.items()):
            rankings = data['rankings']
            rank_values = list(rankings.values())

            row = {
                'GO_Term': term,
                'Enrichment_Type': data['type'],
                'First_Appearance': min(rankings.keys()),
                'Last_Appearance': max(rankings.keys()),
                'Times_In_Top20': len(rankings),
                'Avg_Rank': f"{statistics.mean(rank_values):.1f}",
                'Rank_StdDev': f"{statistics.stdev(rank_values):.2f}" if len(rank_values) > 1 else '0',
                'Rank_Range': f"{min(rank_values)}-{max(rank_values)}"
            }

            # Add individual timepoint ranks
            for tp in timepoints:
                row[f'{tp}_Rank'] = rankings.get(tp, 'NA')

            writer.writerow(row)

    logger.info(f"Rank stability analysis written to: {output_file}")


def write_gene_classification_changes(output_prefix, gene_timeline, changing_genes):
    """Write analysis of genes that changed classification."""
    output_file = f"{output_prefix}_gene_classification_changes.csv"

    timepoints = sorted(list(gene_timeline[list(gene_timeline.keys())[0]].keys()))

    with open(output_file, 'w', newline='') as f:
        fieldnames = ['Gene', 'Change_Pattern', 'Num_Changes', 'First_Class', 'Last_Class'] + \
                     [f'{tp}_Classification' for tp in timepoints]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for gene, timeline in sorted(changing_genes.items()):
            classifications = [c for c in timeline.values() if c != 'absent']

            # Count transitions
            transitions = 0
            for i in range(len(classifications) - 1):
                if classifications[i] != classifications[i + 1]:
                    transitions += 1

            # Determine pattern
            unique_classes = set(classifications)
            if len(unique_classes) == 2:
                if 'lethal' in unique_classes and 'viable' in unique_classes:
                    pattern = 'lethal<->viable'
                else:
                    pattern = 'other_changes'
            else:
                pattern = 'multiple_states'

            row = {
                'Gene': gene,
                'Change_Pattern': pattern,
                'Num_Changes': transitions,
                'First_Class': classifications[0] if classifications else 'NA',
                'Last_Class': classifications[-1] if classifications else 'NA'
            }

            for tp in timepoints:
                row[f'{tp}_Classification'] = timeline.get(tp, 'absent')

            writer.writerow(row)

    logger.info(f"Gene classification changes written to: {output_file}")


def write_summary_statistics_timeline(output_prefix, datasets_dict):
    """Write summary statistics over time."""
    output_file = f"{output_prefix}_summary_timeline.csv"

    with open(output_file, 'w', newline='') as f:
        fieldnames = ['Timepoint', 'Total_Genes', 'Lethal_Genes', 'Viable_Genes',
                      'Lethal_Ratio', 'Viable_Ratio', 'Total_GO_Terms']

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for timepoint, (genes_dict, (_, summary)) in sorted(datasets_dict.items()):
            total_genes = summary['total_genes']
            lethal_genes = summary['lethal_genes']
            viable_genes = summary['viable_genes']

            writer.writerow({
                'Timepoint': timepoint,
                'Total_Genes': total_genes,
                'Lethal_Genes': lethal_genes,
                'Viable_Genes': viable_genes,
                'Lethal_Ratio': f"{lethal_genes / total_genes:.3f}" if total_genes > 0 else '0',
                'Viable_Ratio': f"{viable_genes / total_genes:.3f}" if total_genes > 0 else '0',
                'Total_GO_Terms': len([g for g in genes_dict.values()])
            })

    logger.info(f"Summary timeline written to: {output_file}")


def write_classification_stability_summary(output_prefix, gene_timeline, changing_genes, stable_genes):
    """Write summary of classification stability."""
    output_file = f"{output_prefix}_classification_stability_summary.csv"

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Metric', 'Count', 'Percentage'])

        total_genes = len(gene_timeline)
        stable_count = len(stable_genes)
        changing_count = len(changing_genes)

        writer.writerow(['Total Genes Tracked', total_genes, '100%'])
        writer.writerow(['Stable Classification', stable_count,
                         f"{stable_count / total_genes * 100:.1f}%"])
        writer.writerow(['Changed Classification', changing_count,
                         f"{changing_count / total_genes * 100:.1f}%"])

        # Breakdown of changing genes
        writer.writerow([])
        writer.writerow(['Change Type', 'Count', 'Percentage of Changing'])

        lethal_to_viable = 0
        viable_to_lethal = 0
        multiple_changes = 0

        for gene, timeline in changing_genes.items():
            classifications = [c for c in timeline.values() if c != 'absent']
            unique_classes = set(classifications)

            if len(unique_classes) == 2:
                if classifications[0] == 'lethal':
                    lethal_to_viable += 1
                else:
                    viable_to_lethal += 1
            else:
                multiple_changes += 1

        writer.writerow(['Lethal → Viable', lethal_to_viable,
                         f"{lethal_to_viable / changing_count * 100:.1f}%"])
        writer.writerow(['Viable → Lethal', viable_to_lethal,
                         f"{viable_to_lethal / changing_count * 100:.1f}%"])
        writer.writerow(['Multiple Changes', multiple_changes,
                         f"{multiple_changes / changing_count * 100:.1f}%"])

    logger.info(f"Classification stability summary written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=f"PhenoGO {PhenGO_VERSION} - Temporal ARFF Analysis: Track GO terms and gene classifications over time",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s -arff_files file_2015.arff file_2016.arff ... file_2024.arff \\
           -timepoints 2015 2016 ... 2024 \\
           -o output_prefix \\
           --top-n 20
        """
    )

    parser.add_argument("-arff_files", nargs='+', required=True,
                        help="ARFF files in temporal order")
    parser.add_argument("-timepoints", nargs='+', required=True,
                        help="Timepoint labels (years, dates, etc.) corresponding to each ARFF file")
    parser.add_argument("-o", dest="output", required=True,
                        help="Output file prefix (multiple files will be generated)")
    parser.add_argument("-output_dir", dest="output_dir", required=False, default='.',
                        help="Directory to write output files (default: current directory)")
    parser.add_argument("--top-n", dest="top_n", type=int, default=20,
                        help="Number of top GO terms to track (default: 20)")

    args = parser.parse_args()

    # Prepare output directory and logging
    args.output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    logger = configure_logger('PhenGO.GO_temporal_analysis', enable_file=True, log_dir=args.output_dir, logfile_name='GO_temporal_analysis.log', level=logging.INFO)

    # Use output prefix path in given directory
    output_prefix = os.path.join(args.output_dir, args.output)

    # Validate inputs
    if len(args.arff_files) != len(args.timepoints):
        logger.error("Error: Number of ARFF files must match number of timepoints")
        return

    logger.info(f"Loading {len(args.arff_files)} ARFF files...")

    # Load all datasets
    datasets_dict = {}
    all_go_terms = set()

    for arff_file, timepoint in zip(args.arff_files, args.timepoints):
        logger.info(f"  Loading {timepoint}: {arff_file}")
        genes, terms = parse_arff_with_terms(arff_file)
        enrichment_stats, summary = calculate_enrichment_stats(genes, terms)

        datasets_dict[timepoint] = (genes, (enrichment_stats, summary))
        all_go_terms.update(terms)

        logger.info(f"    {len(genes)} genes, {len(terms)} GO terms")

    logger.info(f"\nTotal unique GO terms across all timepoints: {len(all_go_terms)}")

    # Track gene classifications
    logger.info("\nTracking gene classification changes...")
    gene_timeline, changing_genes, stable_genes = track_gene_classifications(datasets_dict)
    logger.info(f"  Total genes tracked: {len(gene_timeline)}")
    logger.info(f"  Stable classifications: {len(stable_genes)}")
    logger.info(f"  Changed classifications: {len(changing_genes)}")

    # Generate outputs
    logger.info("\nGenerating temporal analysis outputs...")

    write_temporal_top_terms(output_prefix, datasets_dict, args.top_n)
    write_rank_stability(output_prefix, datasets_dict, args.top_n)
    write_gene_classification_changes(output_prefix, gene_timeline, changing_genes)
    write_summary_statistics_timeline(output_prefix, datasets_dict)
    write_classification_stability_summary(output_prefix, gene_timeline, changing_genes, stable_genes)

    logger.info("\n✅ Temporal analysis complete!")
    logger.info(f"\nOutput files generated with prefix: {output_prefix}")
    logger.info(f"  - {output_prefix}_lethal_enriched_timeline.csv")
    logger.info(f"  - {output_prefix}_viable_enriched_timeline.csv")
    logger.info(f"  - {output_prefix}_rank_stability.csv")
    logger.info(f"  - {output_prefix}_gene_classification_changes.csv")
    logger.info(f"  - {output_prefix}_summary_timeline.csv")
    logger.info(f"  - {output_prefix}_classification_stability_summary.csv")


if __name__ == "__main__":
    main()