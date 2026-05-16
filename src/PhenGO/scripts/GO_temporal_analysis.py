#!/usr/bin/env python3
"""
ARFF GO Term Temporal Analysis Script

This script analyzes multiple ARFF files across time points to track:
1. Changes in top overrepresented GO terms for lethal/viable genes
2. Changes in gene lethal/viable classifications over time
3. GO term ranking stability and trends
4. Feature-set evolution (GO term churn, Jaccard similarity, cumulative growth)
5. Gene lifecycle (first/last appearance, classification consistency)
6. GO term lifecycle (entry/exit year, enrichment trajectory)
7. Pairwise dataset similarity matrix (gene sets and GO feature sets)
8. Data sparsity over time (GO terms per gene distribution)
9. GO term granularity / tree depth over time

Input modes:
  -input_dir  : Directory containing PhenGO output subdirectories (auto-discovers
                ARFF files; subdirectory name is used as the timepoint label).
  -arff_files : Explicit list of ARFF files (requires matching -timepoints list).
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
import json
from collections import defaultdict, Counter, deque
import statistics
import logging
import os

from ..constants import configure_logger, PhenGO_VERSION

logger = logging.getLogger(__name__)


# ── Input discovery ───────────────────────────────────────────────────────────

def discover_datasets_from_dir(parent_dir):
    """Scan *parent_dir* for PhenGO output subdirectories.

    Each subdirectory that contains a ``*_PhenGO.arff`` file is treated as one
    timepoint.  The subdirectory name is used as the timepoint label (e.g.
    ``Yeast_2025`` → timepoint ``Yeast_2025``).

    Also looks for ``GO_Children&Parents.json`` (present when the run used the
    ``obo_to_graph`` step) to enable GO term depth calculations.

    Returns
    -------
    list of (timepoint_label: str, arff_path: str, go_json_path: str | None)
        Sorted by timepoint label.
    """
    results = []
    parent_dir = os.path.abspath(parent_dir)
    for entry in sorted(os.listdir(parent_dir)):
        subdir = os.path.join(parent_dir, entry)
        if not os.path.isdir(subdir):
            continue
        arff_candidates = [f for f in os.listdir(subdir) if f.endswith('_PhenGO.arff')]
        if not arff_candidates:
            continue
        arff_path = os.path.join(subdir, arff_candidates[0])
        go_json   = os.path.join(subdir, 'GO_Children&Parents.json')
        go_json_path = go_json if os.path.exists(go_json) else None
        results.append((entry, arff_path, go_json_path))
    return results


def load_go_depths(json_path):
    """Compute the tree depth of every GO term from ``GO_Children&Parents.json``.

    The JSON maps each GO term to ``{'p': [parents], 'c': [children]}``.
    Depth 0 is assigned to root terms (no parents); deeper terms are more
    specific / granular.

    A BFS is run forward through children edges from all roots.  Terms that
    are not reachable (disconnected / obsolete) are omitted from the result.

    Returns
    -------
    dict  {go_term: depth_int}
    """
    with open(json_path) as f:
        go_graph = json.load(f)

    # Roots = terms with no parents
    roots = [term for term, data in go_graph.items() if not data.get('p')]
    depths = {}
    queue  = deque()
    for root in roots:
        depths[root] = 0
        queue.append(root)

    while queue:
        term  = queue.popleft()
        depth = depths[term]
        for child in go_graph.get(term, {}).get('c', []):
            if child not in depths:
                depths[child] = depth + 1
                queue.append(child)

    return depths


# ─────────────────────────────────────────────────────────────────────────────

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


def write_summary_statistics_timeline(output_prefix, datasets_dict, depths_map=None):
    """Write summary statistics over time.

    depths_map : dict | None
        Optional mapping  timepoint -> {go_term: depth}  produced by
        :func:`load_go_depths`.  When provided, GO term depth / granularity
        metrics are included in the output.
    """
    output_file = f"{output_prefix}_summary_timeline.csv"

    with open(output_file, 'w', newline='') as f:
        fieldnames = ['Timepoint', 'Total_Genes', 'Lethal_Genes', 'Viable_Genes',
                      'Lethal_Ratio', 'Viable_Ratio', 'Total_GO_Terms',
                      'Avg_GO_Per_Gene', 'Median_GO_Per_Gene',
                      'Mean_GO_Term_Depth', 'Mean_GO_Term_Depth_Weighted', 'Pct_Terms_With_Depth']

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for timepoint, (genes_dict, (enrichment_stats, summary)) in sorted(datasets_dict.items()):
            total_genes = summary['total_genes']
            lethal_genes = summary['lethal_genes']
            viable_genes = summary['viable_genes']
            # Count number of GO terms that appear at least once in this dataset
            total_go_terms = sum(1 for stats in enrichment_stats.values()
                                 if stats['total_genes_with_term'] > 0)

            # GOs per gene statistics (reuse same counting as sparsity analysis)
            go_counts = []
            for gene, info in genes_dict.items():
                n_go = sum(1 for v in info['features'].values() if v not in ('0', 'NA', '', 0))
                go_counts.append(n_go)

            if go_counts:
                avg_go_per_gene = statistics.mean(go_counts)
                median_go_per_gene = statistics.median(go_counts)
            else:
                avg_go_per_gene = 0
                median_go_per_gene = 0

            # GO term depth / granularity metrics (if available)
            depth_map = depths_map.get(timepoint) if depths_map else None
            mean_depth = 'NA'
            mean_depth_weighted = 'NA'
            pct_with_depth = '0'

            if depth_map:
                # Consider only active GO terms in this dataset
                active_terms = [t for t, s in enrichment_stats.items() if s['total_genes_with_term'] > 0]
                depths = [depth_map.get(t) for t in active_terms if depth_map.get(t) is not None]
                if depths:
                    mean_depth = statistics.mean(depths)
                    # Weighted mean: weight by number of genes annotated with term
                    weights = [enrichment_stats[t]['total_genes_with_term'] for t in active_terms if depth_map.get(t) is not None]
                    if sum(weights) > 0:
                        mean_depth_weighted = sum(d * w for d, w in zip(depths, weights)) / sum(weights)
                    else:
                        mean_depth_weighted = 'NA'
                    pct_with_depth = f"{len(depths) / len(active_terms) * 100:.1f}" if active_terms else '0'
                else:
                    mean_depth = 'NA'
                    mean_depth_weighted = 'NA'
                    pct_with_depth = '0'

            writer.writerow({
                'Timepoint': timepoint,
                'Total_Genes': total_genes,
                'Lethal_Genes': lethal_genes,
                'Viable_Genes': viable_genes,
                'Lethal_Ratio': f"{lethal_genes / total_genes:.3f}" if total_genes > 0 else '0',
                'Viable_Ratio': f"{viable_genes / total_genes:.3f}" if total_genes > 0 else '0',
                'Total_GO_Terms': total_go_terms,
                'Avg_GO_Per_Gene': f"{avg_go_per_gene:.2f}",
                'Median_GO_Per_Gene': f"{median_go_per_gene:.1f}",
                'Mean_GO_Term_Depth': f"{mean_depth:.3f}" if isinstance(mean_depth, float) else mean_depth,
                'Mean_GO_Term_Depth_Weighted': f"{mean_depth_weighted:.3f}" if isinstance(mean_depth_weighted, float) else mean_depth_weighted,
                'Pct_Terms_With_Depth': pct_with_depth,
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

        if changing_count > 0:
            writer.writerow(['Lethal → Viable', lethal_to_viable,
                             f"{lethal_to_viable / changing_count * 100:.1f}%"])
            writer.writerow(['Viable → Lethal', viable_to_lethal,
                             f"{viable_to_lethal / changing_count * 100:.1f}%"])
            writer.writerow(['Multiple Changes', multiple_changes,
                             f"{multiple_changes / changing_count * 100:.1f}%"])

    logger.info(f"Classification stability summary written to: {output_file}")


# ── New analyses ──────────────────────────────────────────────────────────────

def analyse_feature_set_evolution(output_prefix, datasets_dict):
    """Track GO feature-set growth, churn, and similarity across timepoints.

    For each consecutive pair of timepoints reports:
      - Total GO terms present (at least once in any gene)
      - New terms (in B not A), Lost terms (in A not B), Retained
      - Jaccard similarity with the previous year
      - Cumulative unique GO terms seen up to this point

    Output: {prefix}_feature_set_evolution.csv
    """
    timepoints = sorted(datasets_dict.keys())

    # For each timepoint collect the set of *active* GO terms
    # (terms that are actually present in at least one gene)
    def active_go_terms(enrichment_stats):
        return {t for t, s in enrichment_stats.items() if s['total_genes_with_term'] > 0}

    tp_go_sets = {}
    for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
        tp_go_sets[tp] = active_go_terms(enrichment_stats)

    cumulative = set()
    rows = []

    for i, tp in enumerate(timepoints):
        current = tp_go_sets[tp]
        cumulative |= current
        prev = tp_go_sets[timepoints[i - 1]] if i > 0 else set()

        new_terms      = current - prev  if i > 0 else current
        lost_terms     = prev - current  if i > 0 else set()
        retained_terms = current & prev  if i > 0 else set()

        union = current | prev
        jaccard = len(retained_terms) / len(union) if union and i > 0 else float('nan')

        rows.append({
            'Timepoint':               tp,
            'Total_GO_Terms':          len(current),
            'New_Terms':               len(new_terms),
            'Lost_Terms':              len(lost_terms),
            'Retained_Terms':          len(retained_terms),
            'Jaccard_Previous':        f"{jaccard:.4f}" if i > 0 else 'NA',
            'Cumulative_Unique_Terms': len(cumulative),
        })

    output_file = f"{output_prefix}_feature_set_evolution.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Feature-set evolution written to: {output_file}")


def analyse_gene_lifecycle(output_prefix, datasets_dict):
    """Per-gene appearance, disappearance, and classification consistency.

    For every gene seen across all timepoints reports:
      - First and last year present
      - Number of years present / absent
      - Classification at each timepoint ('lethal', 'viable', 'absent')
      - Whether the gene ever changed classification

    Output: {prefix}_gene_lifecycle.csv
    """
    timepoints = sorted(datasets_dict.keys())

    all_genes = set()
    for genes_dict, _ in datasets_dict.values():
        all_genes.update(genes_dict.keys())

    rows = []
    for gene in sorted(all_genes):
        timeline = {}
        for tp, (genes_dict, _) in sorted(datasets_dict.items()):
            timeline[tp] = genes_dict[gene]['label'].lower() if gene in genes_dict else 'absent'

        present_tps = [tp for tp, cls in timeline.items() if cls != 'absent']
        classes_seen = [timeline[tp] for tp in present_tps]
        unique_classes = set(classes_seen)

        row = {
            'Gene':                gene,
            'First_Year':          min(present_tps) if present_tps else 'NA',
            'Last_Year':           max(present_tps) if present_tps else 'NA',
            'Years_Present':       len(present_tps),
            'Years_Absent':        len(timepoints) - len(present_tps),
            'Classification_Summary': (
                'consistent_lethal'  if unique_classes == {'lethal'}  else
                'consistent_viable'  if unique_classes == {'viable'}  else
                'mixed'              if len(unique_classes) > 1       else
                'absent_all'
            ),
            'Num_Classification_Changes': sum(
                1 for i in range(len(classes_seen) - 1)
                if classes_seen[i] != classes_seen[i + 1]
            ),
        }
        for tp in timepoints:
            row[f'{tp}_Class'] = timeline[tp]

        rows.append(row)

    if not rows:
        return

    output_file = f"{output_prefix}_gene_lifecycle.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Quick summary log
    n_consistent_lethal = sum(1 for r in rows if r['Classification_Summary'] == 'consistent_lethal')
    n_consistent_viable = sum(1 for r in rows if r['Classification_Summary'] == 'consistent_viable')
    n_mixed             = sum(1 for r in rows if r['Classification_Summary'] == 'mixed')
    n_all_timepoints    = sum(1 for r in rows if r['Years_Present'] == len(timepoints))

    logger.info(f"Gene lifecycle written to: {output_file}")
    logger.info(f"  Consistently lethal:  {n_consistent_lethal}")
    logger.info(f"  Consistently viable:  {n_consistent_viable}")
    logger.info(f"  Mixed classification: {n_mixed}")
    logger.info(f"  Present in all {len(timepoints)} timepoints: {n_all_timepoints}")


def analyse_go_term_lifecycle(output_prefix, datasets_dict):
    """Per-GO-term appearance trajectory and enrichment trend.

    For every GO term seen in any timepoint reports:
      - First and last year present (active in ≥1 gene)
      - Number of years present
      - Mean and trend of lethal enrichment ratio across years
      - Mean lethal frequency and viable frequency

    Output: {prefix}_go_term_lifecycle.csv
    """
    timepoints = sorted(datasets_dict.keys())

    all_terms = set()
    tp_stats = {}
    for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
        tp_stats[tp] = enrichment_stats
        all_terms.update(enrichment_stats.keys())

    rows = []
    for term in sorted(all_terms):
        present_tps = [tp for tp in timepoints
                       if tp_stats[tp].get(term, {}).get('total_genes_with_term', 0) > 0]

        enrichment_vals = [
            tp_stats[tp][term]['enrichment_ratio']
            for tp in present_tps
            if tp_stats[tp].get(term) and tp_stats[tp][term]['enrichment_ratio'] != float('inf')
        ]
        lethal_freqs = [tp_stats[tp][term]['lethal_freq'] for tp in present_tps if tp_stats[tp].get(term)]
        viable_freqs = [tp_stats[tp][term]['viable_freq'] for tp in present_tps if tp_stats[tp].get(term)]

        def _mean(lst): return statistics.mean(lst) if lst else float('nan')

        # Simple linear trend: positive = enrichment increasing over time
        trend = float('nan')
        if len(enrichment_vals) >= 2:
            n = len(enrichment_vals)
            x_mean = (n - 1) / 2.0
            num = sum((i - x_mean) * (v - _mean(enrichment_vals)) for i, v in enumerate(enrichment_vals))
            den = sum((i - x_mean) ** 2 for i in range(n))
            trend = num / den if den else float('nan')

        row = {
            'GO_Term':             term,
            'First_Year':          min(present_tps) if present_tps else 'NA',
            'Last_Year':           max(present_tps) if present_tps else 'NA',
            'Years_Present':       len(present_tps),
            'Mean_Enrichment':     f"{_mean(enrichment_vals):.4f}" if enrichment_vals else 'NA',
            'Enrichment_Trend':    f"{trend:.6f}" if trend == trend else 'NA',  # nan check
            'Mean_Lethal_Freq':    f"{_mean(lethal_freqs):.4f}",
            'Mean_Viable_Freq':    f"{_mean(viable_freqs):.4f}",
        }
        for tp in timepoints:
            s = tp_stats[tp].get(term)
            row[f'{tp}_EnrichRatio'] = (
                f"{s['enrichment_ratio']:.4f}"
                if s and s['total_genes_with_term'] > 0 and s['enrichment_ratio'] != float('inf')
                else 'NA'
            )

        rows.append(row)

    if not rows:
        return

    output_file = f"{output_prefix}_go_term_lifecycle.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"GO term lifecycle written to: {output_file}")
    logger.info(f"  Total unique GO terms tracked: {len(rows)}")
    logger.info(f"  Present in all {len(timepoints)} timepoints: "
                f"{sum(1 for r in rows if r['Years_Present'] == len(timepoints))}")


def compute_pairwise_similarity(output_prefix, datasets_dict):
    """Compute pairwise Jaccard similarity between all dataset pairs.

    Two matrices are written:
      - Gene-set Jaccard:     |genes_A ∩ genes_B| / |genes_A ∪ genes_B|
      - GO feature-set Jaccard: |GO_A ∩ GO_B| / |GO_A ∪ GO_B|

    Output: {prefix}_pairwise_gene_jaccard.csv
            {prefix}_pairwise_go_jaccard.csv
    """
    timepoints = sorted(datasets_dict.keys())

    # Build gene sets and active-GO-term sets per timepoint
    gene_sets = {}
    go_sets   = {}
    for tp, (genes_dict, (enrichment_stats, _)) in sorted(datasets_dict.items()):
        gene_sets[tp] = set(genes_dict.keys())
        go_sets[tp]   = {t for t, s in enrichment_stats.items() if s['total_genes_with_term'] > 0}

    def jaccard(a, b):
        union = a | b
        return len(a & b) / len(union) if union else float('nan')

    def write_matrix(path, sets):
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([''] + timepoints)
            for tp_a in timepoints:
                row = [tp_a] + [
                    f"{jaccard(sets[tp_a], sets[tp_b]):.4f}" for tp_b in timepoints
                ]
                writer.writerow(row)
        logger.info(f"Pairwise similarity matrix written to: {path}")

    write_matrix(f"{output_prefix}_pairwise_gene_jaccard.csv",   gene_sets)
    write_matrix(f"{output_prefix}_pairwise_go_jaccard.csv",     go_sets)

    # Log the consecutive-year Jaccard scores as a quick diagnostic
    logger.info("Consecutive-year Jaccard (gene sets):")
    for i in range(1, len(timepoints)):
        j = jaccard(gene_sets[timepoints[i-1]], gene_sets[timepoints[i]])
        logger.info(f"  {timepoints[i-1]} → {timepoints[i]}: {j:.4f}")
    logger.info("Consecutive-year Jaccard (GO feature sets):")
    for i in range(1, len(timepoints)):
        j = jaccard(go_sets[timepoints[i-1]], go_sets[timepoints[i]])
        logger.info(f"  {timepoints[i-1]} → {timepoints[i]}: {j:.4f}")


def analyse_go_term_depth_distribution(output_prefix, datasets_dict, depths_map):
    """Report how GO term annotation depth (specificity) changes over time.

    GO terms closer to the root (depth 0-2) are very general; deeper terms are
    more specific.  Tracking the depth distribution of *active* terms (those
    present in ≥1 gene in the ARFF) shows whether annotations have become more
    or less specific over time.

    Two files are written:

    1. ``{prefix}_go_depth_distribution.csv``
       Wide table: rows = timepoints, columns = depth levels (0 … max observed).
       Cell value = number of active GO terms at that depth for that timepoint.
       Includes total active terms and % with depth data.

    2. ``{prefix}_go_depth_stats.csv``
       Per-timepoint depth statistics for three subsets:
         - All active terms
         - Lethal-enriched terms (enrichment_ratio > 1)
         - Viable-enriched terms (depletion_ratio > 1 and enrichment_ratio < 1)
       Stats: mean, median, P25, P75, min, max depth.

    Parameters
    ----------
    depths_map : dict
        ``{timepoint: {go_term: depth_int}}`` as built in ``main()``.
        If empty or None, this function logs a warning and returns.
    """
    if not depths_map:
        logger.warning("No GO depth data available — skipping depth distribution analysis.")
        return

    timepoints = sorted(datasets_dict.keys())

    # ── Collect per-timepoint depth data ─────────────────────────────────────
    tp_data = {}   # timepoint -> {'all': [depths], 'lethal': [depths], 'viable': [depths]}

    for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
        dm = depths_map.get(tp)
        if not dm:
            tp_data[tp] = None
            continue

        all_depths     = []
        lethal_depths  = []
        viable_depths  = []

        for term, stats in enrichment_stats.items():
            if stats['total_genes_with_term'] == 0:
                continue
            depth = dm.get(term)
            if depth is None:
                continue
            all_depths.append(depth)
            ratio = stats['enrichment_ratio']
            if ratio != float('inf') and ratio > 1:
                lethal_depths.append(depth)
            elif ratio != float('inf') and ratio < 1:
                viable_depths.append(depth)

        tp_data[tp] = {
            'all':    all_depths,
            'lethal': lethal_depths,
            'viable': viable_depths,
        }

    # ── File 1: depth count matrix ────────────────────────────────────────────
    # Determine all depth levels seen across all timepoints
    all_depth_levels = sorted({
        d for data in tp_data.values()
        if data for d in data['all']
    })

    dist_file = f"{output_prefix}_go_depth_distribution.csv"
    with open(dist_file, 'w', newline='') as f:
        header = ['Timepoint', 'Total_Active_GO_Terms', 'Terms_With_Depth', 'Pct_With_Depth'] + \
                 [f'Depth_{d}' for d in all_depth_levels]
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()

        for tp, (_, (enrichment_stats, _)) in sorted(datasets_dict.items()):
            data = tp_data.get(tp)
            total_active = sum(1 for s in enrichment_stats.values()
                               if s['total_genes_with_term'] > 0)

            if not data:
                row = {'Timepoint': tp, 'Total_Active_GO_Terms': total_active,
                       'Terms_With_Depth': 0, 'Pct_With_Depth': '0.0'}
                for d in all_depth_levels:
                    row[f'Depth_{d}'] = 0
            else:
                depth_counts = Counter(data['all'])
                n_with_depth = len(data['all'])
                row = {
                    'Timepoint':           tp,
                    'Total_Active_GO_Terms': total_active,
                    'Terms_With_Depth':    n_with_depth,
                    'Pct_With_Depth':      f"{n_with_depth / total_active * 100:.1f}" if total_active else '0.0',
                }
                for d in all_depth_levels:
                    row[f'Depth_{d}'] = depth_counts.get(d, 0)

            writer.writerow(row)

    logger.info(f"GO depth distribution matrix written to: {dist_file}")

    # ── File 2: depth statistics split by enrichment group ───────────────────
    def _depth_stats(depths):
        if not depths:
            return {k: 'NA' for k in ('Mean', 'Median', 'P25', 'P75', 'Min', 'Max', 'N')}
        s = sorted(depths)
        n = len(s)
        def pct(p):
            k = (n - 1) * p / 100
            lo, hi = int(k), min(int(k) + 1, n - 1)
            return s[lo] + (s[hi] - s[lo]) * (k - lo)
        return {
            'Mean':   f"{statistics.mean(depths):.3f}",
            'Median': f"{statistics.median(depths):.1f}",
            'P25':    f"{pct(25):.1f}",
            'P75':    f"{pct(75):.1f}",
            'Min':    min(depths),
            'Max':    max(depths),
            'N':      n,
        }

    stats_file = f"{output_prefix}_go_depth_stats.csv"
    fieldnames = [
        'Timepoint',
        'All_Mean_Depth', 'All_Median_Depth', 'All_P25', 'All_P75',
        'All_Min', 'All_Max', 'All_N',
        'Lethal_Enriched_Mean', 'Lethal_Enriched_Median', 'Lethal_Enriched_P25',
        'Lethal_Enriched_P75', 'Lethal_Enriched_Min', 'Lethal_Enriched_Max', 'Lethal_Enriched_N',
        'Viable_Enriched_Mean', 'Viable_Enriched_Median', 'Viable_Enriched_P25',
        'Viable_Enriched_P75', 'Viable_Enriched_Min', 'Viable_Enriched_Max', 'Viable_Enriched_N',
    ]

    with open(stats_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for tp in timepoints:
            data = tp_data.get(tp)
            if not data:
                row = {'Timepoint': tp}
                for field in fieldnames[1:]:
                    row[field] = 'NA'
            else:
                a = _depth_stats(data['all'])
                l = _depth_stats(data['lethal'])
                v = _depth_stats(data['viable'])
                row = {
                    'Timepoint':                tp,
                    'All_Mean_Depth':           a['Mean'],
                    'All_Median_Depth':         a['Median'],
                    'All_P25':                  a['P25'],
                    'All_P75':                  a['P75'],
                    'All_Min':                  a['Min'],
                    'All_Max':                  a['Max'],
                    'All_N':                    a['N'],
                    'Lethal_Enriched_Mean':     l['Mean'],
                    'Lethal_Enriched_Median':   l['Median'],
                    'Lethal_Enriched_P25':      l['P25'],
                    'Lethal_Enriched_P75':      l['P75'],
                    'Lethal_Enriched_Min':      l['Min'],
                    'Lethal_Enriched_Max':      l['Max'],
                    'Lethal_Enriched_N':        l['N'],
                    'Viable_Enriched_Mean':     v['Mean'],
                    'Viable_Enriched_Median':   v['Median'],
                    'Viable_Enriched_P25':      v['P25'],
                    'Viable_Enriched_P75':      v['P75'],
                    'Viable_Enriched_Min':      v['Min'],
                    'Viable_Enriched_Max':      v['Max'],
                    'Viable_Enriched_N':        v['N'],
                }
            writer.writerow(row)

    logger.info(f"GO depth statistics written to: {stats_file}")


def analyse_data_sparsity(output_prefix, datasets_dict):
    """Describe the distribution of GO terms per gene over time.

    A sparse dataset (few GO terms per gene) may produce an unreliable ARFF.
    Tracking sparsity over time reveals whether GO annotation coverage is
    improving or degrading.

    Output: {prefix}_sparsity_over_time.csv
    """
    timepoints = sorted(datasets_dict.keys())
    rows = []

    for tp, (genes_dict, _) in sorted(datasets_dict.items()):
        # Count GO features (value == '1') per gene
        go_counts = []
        for gene, info in genes_dict.items():
            n_go = sum(1 for v in info['features'].values() if v not in ('0', 'NA', '', 0))
            go_counts.append(n_go)

        if not go_counts:
            continue

        go_counts_sorted = sorted(go_counts)
        n = len(go_counts_sorted)

        def percentile(lst, p):
            k = (len(lst) - 1) * p / 100
            lo, hi = int(k), min(int(k) + 1, len(lst) - 1)
            return lst[lo] + (lst[hi] - lst[lo]) * (k - lo)

        zero_count = go_counts.count(0)

        rows.append({
            'Timepoint':           tp,
            'Total_Genes':         n,
            'Genes_With_No_GO':    zero_count,
            'Pct_No_GO':           f"{zero_count / n * 100:.1f}",
            'Mean_GO_Per_Gene':    f"{statistics.mean(go_counts):.2f}",
            'Median_GO_Per_Gene':  f"{statistics.median(go_counts):.1f}",
            'P25_GO_Per_Gene':     f"{percentile(go_counts_sorted, 25):.1f}",
            'P75_GO_Per_Gene':     f"{percentile(go_counts_sorted, 75):.1f}",
            'Max_GO_Per_Gene':     max(go_counts),
            'Sparsity_Pct':        f"{zero_count / n * 100:.1f}",
        })

    if not rows:
        return

    output_file = f"{output_prefix}_sparsity_over_time.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Data sparsity over time written to: {output_file}")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=f"PhenGO {PhenGO_VERSION} - Temporal ARFF Analysis: Track GO terms and gene classifications over time",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-discover from a directory of PhenGO output subdirectories:
  %(prog)s -input_dir /data/PhenGO/Yeast \\
           -output_dir ./yeast_temporal \\
           -o yeast

  # Explicit ARFF files:
  %(prog)s -arff_files 2015.arff 2016.arff 2017.arff \\
           -timepoints 2015 2016 2017 \\
           -output_dir ./results \\
           -o yeast
        """
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-input_dir", dest="input_dir",
        help="Directory containing PhenGO output subdirectories.  Each subdirectory "
             "must contain a *_PhenGO.arff file.  The subdirectory name is used as "
             "the timepoint label (e.g. 'Yeast_2025').  GO_Children&Parents.json is "
             "loaded automatically for depth analysis when present.")
    input_group.add_argument(
        "-arff_files", nargs='+',
        help="Explicit ARFF files in temporal order (requires -timepoints).")

    parser.add_argument(
        "-timepoints", nargs='+',
        help="Timepoint labels corresponding to each -arff_files entry.")
    parser.add_argument(
        "-output_dir", dest="output_dir", required=True,
        help="Directory to write all output files (created if absent).")
    parser.add_argument(
        "-o", dest="output", default="temporal_analysis",
        help="Filename prefix for output files (default: temporal_analysis).")
    parser.add_argument(
        "--top-n", dest="top_n", type=int, default=20,
        help="Number of top GO terms to track (default: 20).")

    args = parser.parse_args()

    # ── Validate arguments ────────────────────────────────────────────────────
    if args.arff_files and not args.timepoints:
        parser.error("-timepoints is required when using -arff_files")
    if args.arff_files and len(args.arff_files) != len(args.timepoints):
        parser.error("Number of -arff_files entries must match number of -timepoints entries")

    # ── Set up output directory and logging ───────────────────────────────────
    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    logger = configure_logger('PhenGO.GO_temporal_analysis', enable_file=True,
                              log_dir=args.output_dir,
                              logfile_name='GO_temporal_analysis.log',
                              level=logging.INFO)

    output_prefix = os.path.join(args.output_dir, args.output)

    # ── Collect ARFF files and GO depth data ──────────────────────────────────
    arff_entries = []   # list of (timepoint_label, arff_path, go_json_path or None)

    if args.input_dir:
        logger.info(f"Discovering PhenGO output directories in: {args.input_dir}")
        arff_entries = discover_datasets_from_dir(args.input_dir)
        if not arff_entries:
            logger.error(f"No PhenGO output subdirectories found in {args.input_dir}")
            return
        logger.info(f"Found {len(arff_entries)} timepoints:")
        for tp, arff, go_json in arff_entries:
            depth_note = " (depths available)" if go_json else " (no depth data)"
            logger.info(f"  {tp}: {arff}{depth_note}")
    else:
        for arff_path, tp in zip(args.arff_files, args.timepoints):
            arff_entries.append((tp, arff_path, None))

    # ── Build depths_map from JSON files ─────────────────────────────────────
    depths_map = {}
    for tp, _, go_json in arff_entries:
        if go_json:
            try:
                logger.info(f"  Loading GO depths for {tp} ...")
                depths_map[tp] = load_go_depths(go_json)
                logger.info(f"    {len(depths_map[tp])} terms with depth data")
            except Exception as e:
                logger.warning(f"  Could not load GO depths for {tp}: {e}")

    # ── Load ARFF files ───────────────────────────────────────────────────────
    logger.info(f"\nLoading {len(arff_entries)} ARFF files...")
    datasets_dict = {}
    all_go_terms  = set()

    for timepoint, arff_file, _ in arff_entries:
        if not os.path.exists(arff_file):
            logger.error(f"  ARFF file not found: {arff_file}  (skipping {timepoint})")
            continue
        logger.info(f"  Loading {timepoint}: {arff_file}")
        genes, terms = parse_arff_with_terms(arff_file)
        enrichment_stats, summary = calculate_enrichment_stats(genes, terms)
        datasets_dict[timepoint] = (genes, (enrichment_stats, summary))
        all_go_terms.update(terms)

        # Per-dataset quick stats
        go_counts = [
            sum(1 for v in info['features'].values() if v not in ('0', 'NA', '', 0))
            for info in genes.values()
        ]
        avg_go = statistics.mean(go_counts) if go_counts else 0
        active_go = sum(1 for s in enrichment_stats.values() if s['total_genes_with_term'] > 0)
        logger.info(f"    {len(genes)} genes | {active_go} active GO terms | "
                    f"avg {avg_go:.1f} GO terms/gene")

        if timepoint in depths_map:
            dm = depths_map[timepoint]
            active_depths = [dm[t] for t in terms if t in dm]
            if active_depths:
                logger.info(f"    GO depth: mean={statistics.mean(active_depths):.2f}  "
                            f"median={statistics.median(active_depths):.1f}  "
                            f"max={max(active_depths)}")

    if len(datasets_dict) < 2:
        logger.error("Need at least 2 valid timepoints to run temporal analysis. Exiting.")
        return

    logger.info(f"\nTotal unique GO terms across all timepoints: {len(all_go_terms)}")

    # ── Gene classification tracking ──────────────────────────────────────────
    logger.info("\nTracking gene classification changes...")
    gene_timeline, changing_genes, stable_genes = track_gene_classifications(datasets_dict)
    logger.info(f"  Total genes tracked: {len(gene_timeline)}")
    logger.info(f"  Stable classifications: {len(stable_genes)}")
    logger.info(f"  Changed classifications: {len(changing_genes)}")

    # ── Standard outputs ──────────────────────────────────────────────────────
    logger.info("\nGenerating temporal analysis outputs...")
    write_temporal_top_terms(output_prefix, datasets_dict, args.top_n)
    write_rank_stability(output_prefix, datasets_dict, args.top_n)
    write_gene_classification_changes(output_prefix, gene_timeline, changing_genes)
    write_summary_statistics_timeline(output_prefix, datasets_dict, depths_map=depths_map)
    write_classification_stability_summary(output_prefix, gene_timeline, changing_genes, stable_genes)

    # ── Extended analyses ──────────────────────────────────────────────────────
    logger.info("\nGenerating extended analyses...")
    analyse_feature_set_evolution(output_prefix, datasets_dict)
    analyse_gene_lifecycle(output_prefix, datasets_dict)
    analyse_go_term_lifecycle(output_prefix, datasets_dict)
    compute_pairwise_similarity(output_prefix, datasets_dict)
    analyse_data_sparsity(output_prefix, datasets_dict)
    analyse_go_term_depth_distribution(output_prefix, datasets_dict, depths_map)

    logger.info("\nTemporal analysis complete!")
    logger.info(f"All outputs written to: {args.output_dir}/")
    logger.info(f"  {output_prefix}_summary_timeline.csv")
    logger.info(f"  {output_prefix}_lethal_enriched_timeline.csv")
    logger.info(f"  {output_prefix}_viable_enriched_timeline.csv")
    logger.info(f"  {output_prefix}_rank_stability.csv")
    logger.info(f"  {output_prefix}_gene_classification_changes.csv")
    logger.info(f"  {output_prefix}_classification_stability_summary.csv")
    logger.info(f"  {output_prefix}_feature_set_evolution.csv")
    logger.info(f"  {output_prefix}_gene_lifecycle.csv")
    logger.info(f"  {output_prefix}_go_term_lifecycle.csv")
    logger.info(f"  {output_prefix}_pairwise_gene_jaccard.csv")
    logger.info(f"  {output_prefix}_pairwise_go_jaccard.csv")
    logger.info(f"  {output_prefix}_sparsity_over_time.csv")
    logger.info(f"  {output_prefix}_go_depth_distribution.csv")
    logger.info(f"  {output_prefix}_go_depth_stats.csv")
    logger.info(f"  GO_temporal_analysis.log")


if __name__ == "__main__":
    main()






