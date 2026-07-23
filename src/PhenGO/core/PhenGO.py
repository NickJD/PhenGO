import argparse
import csv
import sys
import os
import logging
from collections import Counter, defaultdict
import networkx as nx
from datetime import datetime

# ── Allow direct execution: python .../core/PhenGO.py  ─────────────────────
# Manually register parent packages in sys.modules so that relative imports
# (from .X / from ..Y) resolve without error when run as a plain script.
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here      = os.path.dirname(os.path.abspath(__file__))
    _src       = os.path.normpath(os.path.join(_here, '..', '..'))
    _phengo    = os.path.dirname(_here)          # .../src/PhenGO/
    if _src not in sys.path:
        sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',      os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.core', os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.core'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo
# ───────────────────────────────────────────────────────────────────────────

from .obo_to_graph import GO_ROOT_TERMS, obo_to_graph
from ..constants import (
    BUNDLED_FLY_LETHAL_GENES_FILE,
    DEFAULT_FLY_HELPER_LINES_FILE,
    DEFAULT_FLY_SPECIES_FIELDS_FILE,
    DEFAULT_MOUSE_PHENOTYPES_FILE,
    DEFAULT_WORM_LETHAL_GENES_FILE,
    DEFAULT_WORM_LETHAL_PHENOTYPES_FILE,
    PhenGO_VERSION,
    configure_logger,
)
from .phenotype_handling import (
    get_viable_inviable_fish,
    get_viable_inviable_fly,
    get_viable_inviable_mouse,
    get_viable_inviable_worm,
    get_viable_inviable_yeast,
    load_gene_set_labels,
)
from .go_handling import (
    get_viability_go_data_fish,
    get_viability_go_data_fly,
    get_viability_go_data_mouse,
    get_viability_go_data_worm,
    get_viability_go_data_yeast,
)
from ..provenance import (
    assert_safe_output_dir,
    build_run_manifest,
    describe_file,
    prepare_output_dir,
    write_manifest,
)

logger = logging.getLogger(__name__)

def removed_unused_gos(vi_inviable_genes, unique_go_terms):
    # This optional function removes unused GO terms from the ARFF and go_enrichment outputs
    # It rebuilds the GO term list and each gene's bin_vec to only include used terms
    # collect *all* expanded GO features assigned to at least 1 gene
    used_gos = set()
    go_term_index = {go: i for i, go in enumerate(unique_go_terms)}
    for gene_vals in vi_inviable_genes.values():
        bin_vec = gene_vals.get("bin_vec", [])
        if bin_vec:
            used_gos.update(
                go for go, idx in go_term_index.items()
                if idx < len(bin_vec) and bin_vec[idx] == 1
            )
        else:
            used_gos.update(gene_vals.get("expanded_go_list", gene_vals.get("go_list", [])))
    # restrict the master list to only the used terms, in order
    filtered_go_terms = [go for go in unique_go_terms if go in used_gos]
    # now rebuild each gene's expanded_go_list AND its bin_vec *only* over filtered_go_terms
    filtered_genes = {}
    for gene_id, gene_vals in vi_inviable_genes.items():
        old_bin_vec = gene_vals.get("bin_vec", [])
        expanded_go_list = []
        for go in filtered_go_terms:
            old_idx = go_term_index[go]
            if old_idx < len(old_bin_vec) and old_bin_vec[old_idx] == 1:
                expanded_go_list.append(go)

        if not expanded_go_list:
            continue

        direct_go_list = list(dict.fromkeys(gene_vals.get("go_list", [])))
        new_bin_vec = [1 if go in expanded_go_list else 0 for go in filtered_go_terms]
        filtered_genes[gene_id] = dict(gene_vals)
        filtered_genes[gene_id].update({
            "go_list": direct_go_list,
            "expanded_go_list": expanded_go_list,
            "bin_vec": new_bin_vec,
        })
    return filtered_genes, filtered_go_terms


def validate_go_terms_against_graph(vi_inviable_genes, gr, obsolete_go_terms, output_dir):
    """
    Validate all GO terms in gene data against the graph.
    Collect and report all issues.
    """
    logger.info("Validating GO terms against OBO graph...")

    # Collection for reporting
    missing_terms = defaultdict(list)  # term -> [genes using it]
    obsolete_terms_found = defaultdict(list)  # term -> [genes using it]
    valid_terms_count = 0
    total_terms_checked = 0
    alternate_ids_mapped = defaultdict(list)
    obsolete_terms_replaced = defaultdict(list)
    alt_id_to_id = gr.graph.get("alt_id_to_id", {})
    obsolete_replacements = gr.graph.get("obsolete_replacements", {})

    # Process each gene
    for gene, values in vi_inviable_genes.items():
        if 'go_list' not in values:
            continue

        cleaned_go_list = []

        for go_term in values['go_list']:
            total_terms_checked += 1

            canonical = alt_id_to_id.get(go_term, go_term)
            if canonical != go_term:
                alternate_ids_mapped[go_term].append(gene)
                go_term = canonical

            # Check if term is obsolete
            if go_term in obsolete_go_terms:
                replacement = obsolete_replacements.get(go_term)
                if replacement and replacement in gr:
                    obsolete_terms_replaced[go_term].append(gene)
                    go_term = replacement
                else:
                    obsolete_terms_found[go_term].append(gene)
                    continue

            # Check if term exists in graph
            if go_term not in gr:
                missing_terms[go_term].append(gene)
                continue

            # Valid term
            if go_term not in cleaned_go_list:
                cleaned_go_list.append(go_term)
            valid_terms_count += 1

        # Update gene with cleaned GO list
        vi_inviable_genes[gene]['go_list'] = cleaned_go_list

    # Generate comprehensive report
    report = {
        'total_terms_checked': total_terms_checked,
        'valid_terms': valid_terms_count,
        'missing_terms': len(missing_terms),
        'obsolete_terms': len(obsolete_terms_found),
        'obsolete_terms_replaced': len(obsolete_terms_replaced),
        'alternate_ids_mapped': len(alternate_ids_mapped),
        'genes_affected': len(set(
            [g for genes in missing_terms.values() for g in genes] +
            [g for genes in obsolete_terms_found.values() for g in genes]
        ))
    }

    # Write detailed report to file
    report_file = os.path.join(output_dir, 'GO_term_validation_report.txt')
    with open(report_file, 'w') as f:
        f.write("GO TERM VALIDATION REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total GO terms checked: {total_terms_checked}\n")
        f.write(f"Valid terms: {valid_terms_count}\n")
        f.write(f"Missing terms: {len(missing_terms)}\n")
        f.write(f"Obsolete terms found: {len(obsolete_terms_found)}\n")
        f.write(f"Obsolete terms replaced: {len(obsolete_terms_replaced)}\n")
        f.write(f"Alternate IDs canonicalised: {len(alternate_ids_mapped)}\n")
        f.write(f"Genes affected: {report['genes_affected']}\n\n")

        if missing_terms:
            f.write("MISSING TERMS (not in OBO file):\n")
            f.write("-" * 80 + "\n")
            for term, genes in sorted(missing_terms.items()):
                f.write(f"{term} - used by {len(genes)} gene(s)\n")
                if len(genes) <= 10:
                    f.write(f"  Genes: {', '.join(genes)}\n")
                else:
                    f.write(f"  Genes: {', '.join(genes[:10])}... ({len(genes) - 10} more)\n")
            f.write("\n")

        if obsolete_terms_found:
            f.write("OBSOLETE TERMS (deprecated in OBO):\n")
            f.write("-" * 80 + "\n")
            for term, genes in sorted(obsolete_terms_found.items()):
                f.write(f"{term} - used by {len(genes)} gene(s)\n")
                if len(genes) <= 10:
                    f.write(f"  Genes: {', '.join(genes)}\n")
                else:
                    f.write(f"  Genes: {', '.join(genes[:10])}... ({len(genes) - 10} more)\n")
            f.write("\n")

        if obsolete_terms_replaced:
            f.write("OBSOLETE TERMS MAPPED VIA replaced_by:\n")
            f.write("-" * 80 + "\n")
            for term, genes in sorted(obsolete_terms_replaced.items()):
                f.write(f"{term} - mapped for {len(genes)} gene(s)\n")
            f.write("\n")

        if alternate_ids_mapped:
            f.write("ALTERNATE IDS MAPPED TO CANONICAL GO IDS:\n")
            f.write("-" * 80 + "\n")
            for term, genes in sorted(alternate_ids_mapped.items()):
                f.write(f"{term} - mapped for {len(genes)} gene(s)\n")
            f.write("\n")

    # Console output
    logger.info("GO term validation complete:")
    logger.info(f"  Valid: {valid_terms_count}/{total_terms_checked}")
    logger.info(f"  Missing: {len(missing_terms)} unique terms")
    logger.info(f"  Obsolete: {len(obsolete_terms_found)} unique terms")

    if missing_terms or obsolete_terms_found:
        logger.warning(f"  {report['genes_affected']} genes affected by invalid GO terms")
        logger.info(f"  Detailed report written to: {report_file}")

    return vi_inviable_genes, report


def Duplicates(Up):
    seen = set()
    result = []
    for node in Up:
        if node not in seen:
            seen.add(node)
            result.append(node)
    return result


def assign_go_to_vector(vi_inviable_genes, gr, unique_go_terms, propagation="ancestors"):
    # Assumes GO terms have already been validated.
    go_term_index = {go: i for i, go in enumerate(unique_go_terms)}

    bin_vec = [0] * len(unique_go_terms)
    o_rep = []
    genes_processed = 0
    total_genes = len(vi_inviable_genes)

    logger.info(f"Assigning GO terms to binary vectors for {total_genes} genes...")

    for gene, values in vi_inviable_genes.items():
        genes_processed += 1
        if genes_processed % 1000 == 0:
            logger.info(f"  Progress: {genes_processed}/{total_genes} genes")

        if 'go_list' not in values or not values['go_list']:
            continue

        # Collect all ancestor GO terms (including the direct terms themselves)
        all_ancestors = set()

        for go_term in values['go_list']:
            if go_term not in gr:
                # Already filtered by validation, but double-check
                continue

            # Add the term itself
            all_ancestors.add(go_term)

            if propagation == "ancestors":
                try:
                    all_ancestors.update(nx.ancestors(gr, go_term))
                except nx.NetworkXError as e:
                    logger.warning(f"Could not get ancestors for {go_term} in gene {gene}: {e}")
                    continue

        expanded_go_list = [go for go in unique_go_terms if go in all_ancestors]

        # Set binary vector and record overrepresentation output entries (one write per unique term)
        for ancestor in expanded_go_list:
            idx = go_term_index.get(ancestor)
            if idx is not None:
                bin_vec[idx] = 1
                o_rep.append(f"{gene}\t{ancestor}\n")
            else:
                logger.debug(f"Ancestor {ancestor} not in unique_go_terms list")

        vi_inviable_genes[gene]["bin_vec"] = bin_vec.copy()
        vi_inviable_genes[gene]["expanded_go_list"] = expanded_go_list
        o_rep.append('\n')
        bin_vec = [0] * len(unique_go_terms)

    logger.info(f"Vector assignment complete for {genes_processed} genes")

    return vi_inviable_genes, unique_go_terms, o_rep


def filter_go_terms_by_prevalence(genes, go_terms, min_count=1, max_fraction=1.0):
    """Apply feature-frequency filters without changing the gene cohort."""
    if not genes or not go_terms:
        return genes, go_terms
    counts = [0] * len(go_terms)
    for values in genes.values():
        for index, value in enumerate(values.get("bin_vec", [])):
            if index < len(counts) and value == 1:
                counts[index] += 1
    n_genes = len(genes)
    keep_indices = [
        index for index, count in enumerate(counts)
        if count >= min_count and count / n_genes <= max_fraction
    ]
    kept_terms = [go_terms[index] for index in keep_indices]
    filtered = {}
    for gene, values in genes.items():
        new_values = dict(values)
        old = values.get("bin_vec", [])
        new_values["bin_vec"] = [old[index] for index in keep_indices]
        new_values["expanded_go_list"] = [
            term for term, value in zip(kept_terms, new_values["bin_vec"]) if value == 1
        ]
        filtered[gene] = new_values
    return filtered, kept_terms


def get_o_rep_output(vi_inviable_genes, o_rep, output_file):
    # Generate overrepresentation analysis output file
    new_o_rep = []

    for gene, values in vi_inviable_genes.items():
        status = values.get("status")
        if status in {"inviable", "lethal"}:
            label = "1"
        elif status == "viable":
            label = "0"
        else:
            raise ValueError(
                f"Cannot write unsupported phenotype label {status!r} for gene {gene}"
            )
        go_terms = values.get("expanded_go_list") or values.get("go_list", [])
        for go_term in dict.fromkeys(go_terms):
            new_o_rep.append([f"{gene}\t{go_term}\t{label}"])

    with open(output_file, mode='w') as FUNCoutputfile:
        for element in new_o_rep:
            FUNCoutputfile.write(" ".join(element) + "\n")


def write_func_output(vi_inviable_genes, output_dir, species):
    """
    Write GO enrichment / over-representation analysis ready files into a
    dedicated 'go_enrichment/' subdirectory of output_dir.

    Four files are produced under {output_dir}/go_enrichment/:

    1. {species}_binary.tab      — GENE <TAB> GO:XXXX <TAB> 0|1
    2. {species}_annotations.tab — GENE <TAB> GO:XXXX  (all genes)
    3. {species}_study_set.txt   — lethal/inviable genes, one per line
    4. {species}_population.txt  — all genes with GO terms, one per line

    These files are tool-agnostic and compatible with FUNC, TopGO, clusterProfiler,
    g:Profiler, GOrilla, PANTHER, and similar overrepresentation analysis frameworks.
    GO terms are the direct annotations after OBO validation and optional
    unused-GO filtering, consistent with the ARFF feature set.
    Genes without GO terms are excluded.
    """
    func_dir = os.path.join(output_dir, "go_enrichment")
    os.makedirs(func_dir, exist_ok=True)

    binary_path = os.path.join(func_dir, f"{species}_binary.tab")
    annot_path  = os.path.join(func_dir, f"{species}_annotations.tab")
    study_path  = os.path.join(func_dir, f"{species}_study_set.txt")
    pop_path    = os.path.join(func_dir, f"{species}_population.txt")

    binary_lines = []
    annot_lines  = []
    lethal_genes = set()
    all_genes    = set()
    seen_annot   = set()

    for gene, values in vi_inviable_genes.items():
        go_list = values.get("expanded_go_list") or values.get("go_list", [])
        if not go_list:
            continue

        status = values.get("status")
        if status in ("lethal", "inviable"):
            label = "1"
        elif status == "viable":
            label = "0"
        else:
            raise ValueError(
                f"Cannot write unsupported phenotype label {status!r} for gene {gene}"
            )

        all_genes.add(gene)
        if status in ("lethal", "inviable"):
            lethal_genes.add(gene)

        for go_term in go_list:
            annot_key = f"{gene}\t{go_term}"
            if annot_key not in seen_annot:
                seen_annot.add(annot_key)
                annot_lines.append(annot_key)
                binary_lines.append(f"{annot_key}\t{label}")

    with open(binary_path, "w") as f:
        for line in binary_lines:
            f.write(line + "\n")

    with open(annot_path, "w") as f:
        for line in annot_lines:
            f.write(line + "\n")

    with open(study_path, "w") as f:
        for gene in sorted(lethal_genes):
            f.write(gene + "\n")

    with open(pop_path, "w") as f:
        for gene in sorted(all_genes):
            f.write(gene + "\n")

    logger.info(f"GO enrichment files written to {func_dir}/:")
    logger.info(f"  {os.path.basename(binary_path)} — binary annotation ({len(binary_lines)} gene-GO entries)")
    logger.info(f"  {os.path.basename(annot_path)} — GO annotations ({len(annot_lines)} gene-GO entries)")
    logger.info(f"  {os.path.basename(study_path)} — study set / lethal genes ({len(lethal_genes)} genes)")
    logger.info(f"  {os.path.basename(pop_path)} — population set ({len(all_genes)} genes)")



def write_arff_output(vi_inviable_genes, filtered_go_terms, output_file):
    if len(filtered_go_terms) != len(set(filtered_go_terms)):
        raise ValueError("Cannot write ARFF with duplicate GO feature names")

    output_rows = []
    for gene, values in vi_inviable_genes.items():
        status = values["status"]
        if status == 'inviable':
            status = 'lethal'
        if status not in {"viable", "lethal"}:
            raise ValueError(
                f"Cannot write unsupported phenotype label {status!r} for gene {gene}"
            )

        bin_vec = values["bin_vec"]
        if len(bin_vec) != len(filtered_go_terms):
            raise ValueError(
                f"Gene {gene} has {len(bin_vec)} GO values but "
                f"{len(filtered_go_terms)} features were declared"
            )
        try:
            invalid_value = any(value not in {0, 1, "0", "1"} for value in bin_vec)
        except TypeError as exc:
            raise ValueError(f"Gene {gene} has a non-scalar GO value") from exc
        if invalid_value:
            raise ValueError(f"Gene {gene} has a GO value outside the binary 0/1 domain")
        output_rows.append([gene, *(int(value) for value in bin_vec), status])

    with open(output_file, 'w') as f:
        f.write("@RELATION gene_lethality\n\n")
        f.write("@ATTRIBUTE gene string\n")
        for go_term in filtered_go_terms:
            f.write(f"@ATTRIBUTE {go_term} {{0,1}}\n")
        f.write("@ATTRIBUTE class {viable,lethal}\n\n")
        f.write("@DATA\n")
        writer = csv.writer(f, lineterminator="\n")
        writer.writerows(output_rows)


def write_gene_identifier_output(genes, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stable_gene_id", "gene_symbol"])
        for gene, values in genes.items():
            writer.writerow([gene, values.get("gene_symbol", "")])


def join_labels_to_go(species, options, gene_association_file, labels,
                      apply_legacy_fly_override=False):
    """Join one label source to GAF records using canonical stable IDs."""
    if species == "yeast":
        return get_viability_go_data_yeast(gene_association_file, labels, options)
    if species == "fly":
        return get_viability_go_data_fly(
            options,
            gene_association_file,
            labels,
            apply_legacy_lethal_override=apply_legacy_fly_override,
        )
    if species == "fish":
        return get_viability_go_data_fish(gene_association_file, labels, options)
    if species == "worm":
        return get_viability_go_data_worm(gene_association_file, labels, options)
    if species == "mouse":
        return get_viability_go_data_mouse(gene_association_file, labels, options)
    raise ValueError(f"Unsupported species: {species}")


def reconcile_label_records(label_source, release_records=None, gene_set_records=None):
    """Select records for one explicit label-source design and build an audit."""
    release_records = release_records or {}
    gene_set_records = gene_set_records or {}
    selected = {}
    audit_rows = []
    counts = Counter()

    for gene in sorted(set(release_records) | set(gene_set_records)):
        release_record = release_records.get(gene)
        gene_set_record = gene_set_records.get(gene)
        release_label = release_record["status"] if release_record else ""
        gene_set_label = gene_set_record["status"] if gene_set_record else ""

        if label_source == "release_records":
            if not release_record:
                continue
            selected[gene] = release_record
            outcome = "retained_release_record"
        elif label_source == "gene_sets":
            if not gene_set_record:
                continue
            selected[gene] = gene_set_record
            outcome = "retained_gene_set"
        elif release_record and gene_set_record and release_label == gene_set_label:
            selected[gene] = release_record
            outcome = "retained_agreement"
        elif release_record and gene_set_record:
            outcome = "excluded_disagreement"
        elif release_record:
            outcome = "excluded_release_only"
        else:
            outcome = "excluded_gene_set_only"

        counts[outcome] += 1
        record = release_record or gene_set_record
        audit_rows.append({
            "stable_gene_id": gene,
            "gene_symbol": record.get("gene_symbol", ""),
            "release_label": release_label,
            "gene_set_label": gene_set_label,
            "outcome": outcome,
        })

    return selected, audit_rows, dict(sorted(counts.items()))


def write_label_source_audit(rows, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "stable_gene_id", "gene_symbol", "release_label",
            "gene_set_label", "outcome",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_identifier_join_audit(rows, output_file):
    """Write identifier exclusions and precedence decisions for review."""
    fieldnames = [
        "database", "input_identifier", "input_status", "outcome",
        "identifier_type", "stable_gene_ids", "details",
    ]
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fly_parse_audit(rows, output_file):
    """Write fail-closed FlyBase compound decisions for independent review."""
    fieldnames = [
        "line_number", "schema_width", "primary", "phenotype", "primary_gene",
        "components", "component_roles", "outcome",
    ]
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_release_labels(species, options, phenotype_file):
    if species == "yeast":
        return get_viable_inviable_yeast(options, phenotype_file)
    if species == "fly":
        return get_viable_inviable_fly(options, phenotype_file)
    if species == "fish":
        return get_viable_inviable_fish(options, phenotype_file)
    if species == "worm":
        return get_viable_inviable_worm(options, phenotype_file)
    if species == "mouse":
        return get_viable_inviable_mouse(options, phenotype_file)
    raise ValueError(f"Unsupported species: {species}")


def main():
    if sys.argv[1:] == ["--print-defaults"]:
        print(f"Fly assignments: {DEFAULT_FLY_SPECIES_FIELDS_FILE}")
        print(f"Bundled fly lethal collection: {BUNDLED_FLY_LETHAL_GENES_FILE}")
        print(f"Fly helper lines: {DEFAULT_FLY_HELPER_LINES_FILE}")
        print(f"Worm lethal phenotypes: {DEFAULT_WORM_LETHAL_PHENOTYPES_FILE}")
        print(f"Worm lethal genes: {DEFAULT_WORM_LETHAL_GENES_FILE}")
        print(f"Mouse lethal phenotypes: {DEFAULT_MOUSE_PHENOTYPES_FILE}")
        return

    parser = argparse.ArgumentParser(description=f"PhenGO {PhenGO_VERSION} - Convert phenotype and GO data to ARFF format")
    parser.add_argument("--print-defaults", action="store_true",
                        help="Print default files and methods used for each species and exit.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required Options')
    required.add_argument('-species', dest="species", required=True, help='Species tag (e.g., fly, yeast)')
    required.add_argument('-phenotype_file', dest="phenotype_file",
                          help='Phenotype release file; required for release_records and agreement label modes')
    required.add_argument('-gene_association_file', dest="gene_association_file", required=True, help='Path to the gene association file (.gz)')
    required.add_argument(
        '-go_obo_file', dest="go_obo_file", required=True,
        help='Path to the release-matched go-basic.obo file',
    )
    required.add_argument('-output_dir', dest="output_dir", required=True, help='Output directory')

    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-no_filter_unused_gos', dest='filter_unused_gos', action='store_false', required=False,
        help='Disable filtering of unused GO terms from the go_enrichment and ARFF output (filtering is ON by default)')
    optional.add_argument('-filter_mixed_terms', dest='filter_mixed_terms', action='store_true', required=False,
                        help='Deprecated alias for -mixed_label_policy exclude.')
    optional.add_argument('-gene_go_pheno', dest='gene_go_pheno', action='store_true', required=False,
                        help='(Deprecated — go_enrichment files are now always written. Kept for backward compatibility.)'
                             ' Previously output a single Gene-GO-Phenotype file; that file is now'
                             ' superseded by the four {species}_*.tab/txt files in go_enrichment/ generated every run.')
    optional.add_argument('-overwrite', dest='overwrite', action='store_true', required=False,
                        help='Overwrite existing non-log contents in output_dir (default: refuse to delete files).')
    optional.add_argument('-mixed_label_policy', choices=['lethal_wins', 'exclude', 'error'],
                        default='exclude', help='Resolve genes with both lethal and viable observations (default: exclude).')
    optional.add_argument('-nonlethal_policy', choices=['observed_viable', 'explicit_only'],
                        default='explicit_only',
                        help='Require explicit viability evidence (default) or treat other observed phenotypes as viable.')
    optional.add_argument('-ambiguous_label_policy', choices=['exclude', 'viable'], default='exclude',
                        help='Handle partial/semi-lethal observations (default: exclude).')

    label_args = parser.add_argument_group('Phenotype label source')
    label_args.add_argument(
        '-label_source',
        choices=['release_records', 'gene_sets', 'agreement'],
        default='release_records',
        help=(
            'Target definition: parse the supplied release (default), use paired '
            'lethal/viable gene sets, or retain only labels on which both sources agree.'
        ),
    )
    label_args.add_argument(
        '-lethal_gene_set',
        help=(
            'Complete lethal collection for gene_sets/agreement mode. The first two TSV '
            'columns may contain stable IDs and symbols; fly also accepts "bundled".'
        ),
    )
    label_args.add_argument(
        '-viable_gene_set',
        help=(
            'Complete viable collection for gene_sets/agreement mode. Genes absent from '
            'both paired collections are unknown and excluded.'
        ),
    )

    provenance_args = parser.add_argument_group('Snapshot provenance')
    provenance_args.add_argument('-snapshot_id', help='Unique dated release/snapshot identifier')
    provenance_args.add_argument('-phenotype_release', help='Phenotype resource release identifier')
    provenance_args.add_argument('-go_annotation_release', help='GO annotation/GAF release identifier')
    provenance_args.add_argument('-go_ontology_release', help='GO ontology release identifier')
    provenance_args.add_argument('-phenotype_ontology_release', help='Phenotype ontology/helper release identifier')
    provenance_args.add_argument('-retrieval_date', help='Input retrieval date')
    provenance_args.add_argument(
        '-snapshot_semantics', default='component_alignment_unspecified',
        choices=[
            'component_alignment_unspecified',
            'declared_composite_cross_section',
            'synchronized_release',
        ],
        help='Whether component files form a declared cross-section or synchronized release.',
    )
    provenance_args.add_argument(
        '-phenotype_availability', default='available_local',
        help='Controlled provenance status for the phenotype input.',
    )
    provenance_args.add_argument(
        '-go_annotation_availability', default='available_local',
        help='Controlled provenance status for the GAF input.',
    )
    provenance_args.add_argument(
        '-go_ontology_availability', default='available_local',
        help='Controlled provenance status for the GO ontology input.',
    )
    provenance_args.add_argument(
        '-retrieval_route', default='local_archived_or_provider_file',
        help='How the selected resource files were obtained.',
    )
    provenance_args.add_argument('-strict_snapshot', action='store_true',
                        help='Require complete release metadata and explicitly supplied dated auxiliary resources.')

    go_args = parser.add_argument_group('GO feature parameters')
    go_args.add_argument('-go_relations', nargs='+', default=['is_a', 'part_of'],
                        choices=['is_a', 'part_of', 'regulates', 'positively_regulates',
                                 'negatively_regulates', 'occurs_in', 'capable_of',
                                 'capable_of_part_of'],
                        help='Relations used for GO ancestor propagation (default: is_a part_of).')
    go_args.add_argument('-go_namespaces', nargs='+',
                        choices=['biological_process', 'molecular_function', 'cellular_component'],
                        help='Restrict features to selected GO namespaces (default: all).')
    go_args.add_argument('-go_propagation', choices=['ancestors', 'direct'], default='ancestors')
    go_args.add_argument(
        '-allow_cross_namespace_go_edges', action='store_true',
        help=(
            'Explicitly permit selected GO relations to propagate between GO aspects. '
            'Off by default; publication runs should normally use go-basic.obo.'
        ),
    )
    go_args.add_argument('-exclude_go_roots', dest='exclude_go_roots', action='store_true', default=True,
                        help='Exclude the three canonical GO namespace roots (default: on).')
    go_args.add_argument('-include_go_roots', dest='exclude_go_roots', action='store_false',
                        help='Retain ontology roots as features; usually non-informative.')
    go_args.add_argument('-go_evidence_codes', nargs='+',
                        help='Include only GAF annotations with these exact evidence codes.')
    go_args.add_argument('-exclude_go_evidence_codes', nargs='+', default=[],
                        help='Exclude exact GAF evidence codes (for example ND IEA).')
    go_args.add_argument('-min_go_gene_count', type=int, default=2,
                        help='Minimum genes carrying a GO feature (default: 2).')
    go_args.add_argument('-max_go_gene_fraction', type=float, default=0.99,
                        help='Maximum fraction of genes carrying a GO feature (default: 0.99).')

    fly_args = parser.add_argument_group('Fly specific parameters')
    fly_args.add_argument('-fly_lethal_genes', dest='fly_lethal_genes', required=False,
                        help='Deprecated lethal-only override for release_records mode; use paired gene-set mode instead. "default" selects the bundled 2025 collection.')
    fly_args.add_argument('-fly_assignments', dest='fly_assignments', required=False,
                        help='Provide TSV file of fly assignments (file confirming genes are assigned to drosophila melanogaster (default: "data/fly/FlyBase_Fields_2025_07_29.tsv.gz")')
    fly_args.add_argument('-fly_helper_lines', dest='fly_helper_lines', required=False,
                        help='Dated driver/helper-line TSV (bundled 2025 file outside strict mode).')
    fly_args.add_argument('-filter_multigenes', dest='filter_multigenes', action='store_true', default=True,
                        help='Filter out phenotypes involving multiple genes: rows with a "with" tag or "/" in '
                             'the gene field are discarded. In simple mode (default) this is a broad string match; '
                             'in driver-aware mode, second-gene and unresolved compound rows are filtered (default: on).')
    fly_args.add_argument('-allow_multigenes', dest='filter_multigenes', action='store_false',
                        help='Include fly multi-gene phenotypes; not recommended for causal gene-level analysis.')
    fly_args.add_argument('-fly_driver_filtering', dest='fly_driver_filtering', action='store_true', required=False,
                        help='Enable fail-closed driver-aware fly phenotype parsing. Recognized accessory '
                             'components and exact same-gene combinations may be retained; second-gene and '
                             'unresolved components are excluded and audited. Non-descript phenotypes follow -nonlethal_policy. '
                             'Default: OFF (uses simpler broad-match filtering).')

    worm_args = parser.add_argument_group('Worm specific parameters')
    worm_args.add_argument('-worm_phenotypes', dest='worm_phenotypes', required=False,
                        help='Provide TSV file of worm phenotypes (default: "data/worm/lethal_terms_traversed_2025-08-12.tsv.gz")')
    worm_args.add_argument('-worm_lethal_genes', dest='worm_lethal_genes', required=False,
                        help='Provide TSV file of specified lethal worm genes (provide "default" for default lethal genes: "data/worm/genes_direct_and_inferred_for_WBPhenotype_0000062_11-08-2025.txt.gz")')
    worm_args.add_argument('-worm_viable_phenotypes',
                        help='Release-matched TSV of explicitly viable WormBase phenotype IDs.')
    worm_args.add_argument('-worm_viable_genes',
                        help='Release-matched TSV of explicitly viable worm genes for gene-list mode.')
    worm_args.add_argument('-worm_evidence_codes', dest='worm_evidence_codes', nargs='+', required=False,
                        help='Only retain worm phenotype annotations with these evidence codes '
                             '(e.g. IMP RNAi). By default all evidence codes are accepted. '
                             'Example: -worm_evidence_codes IMP IBA')

    mouse_args = parser.add_argument_group('Mouse specific parameters')
    mouse_args.add_argument('-mouse_phenotypes', dest='mouse_phenotypes', required=False,
                        help='Provide TSV file of mouse phenotypes (default: "data/mouse/mouse_lethal_terms.txt.gz")')
    mouse_args.add_argument('-mouse_viable_phenotypes',
                        help='Release-matched TSV of explicitly viable MGI phenotype IDs.')

    misc = parser.add_argument_group('Misc')
    misc.add_argument("-v", "--version", action="version",
                 version=f"PhenGO {PhenGO_VERSION}: Exiting.")


    options = parser.parse_args()

    if options.print_defaults:
        print(f"Fly assignments: {DEFAULT_FLY_SPECIES_FIELDS_FILE}")
        print(f"Bundled fly lethal collection: {BUNDLED_FLY_LETHAL_GENES_FILE}")
        print(f"Fly helper lines: {DEFAULT_FLY_HELPER_LINES_FILE}")
        print(f"Worm lethal phenotypes: {DEFAULT_WORM_LETHAL_PHENOTYPES_FILE}")
        print(f"Worm lethal genes: {DEFAULT_WORM_LETHAL_GENES_FILE}")
        print(f"Mouse lethal phenotypes: {DEFAULT_MOUSE_PHENOTYPES_FILE}")
        return

    if options.filter_mixed_terms:
        options.mixed_label_policy = 'exclude'

    species = options.species.lower()
    if species not in {'fly', 'yeast', 'fish', 'worm', 'mouse'}:
        parser.error(
            f"Unknown species '{options.species}'. Supported: fly, yeast, fish, worm, mouse"
        )
    needs_release_labels = options.label_source in {'release_records', 'agreement'}
    needs_gene_sets = options.label_source in {'gene_sets', 'agreement'}
    if needs_release_labels and not options.phenotype_file:
        parser.error(f"-phenotype_file is required for -label_source {options.label_source}")
    if needs_gene_sets and not (options.lethal_gene_set and options.viable_gene_set):
        parser.error(
            f"-label_source {options.label_source} requires both -lethal_gene_set "
            "and -viable_gene_set"
        )
    if not needs_gene_sets and (options.lethal_gene_set or options.viable_gene_set):
        parser.error(
            "-lethal_gene_set and -viable_gene_set require -label_source gene_sets "
            "or agreement"
        )
    if options.fly_lethal_genes:
        if species != 'fly':
            parser.error('-fly_lethal_genes is only valid for fly datasets')
        if options.label_source != 'release_records':
            parser.error(
                '-fly_lethal_genes is a deprecated release_records-only override; '
                'do not combine it with paired gene-set modes'
            )

    explicit_label_sets = {
        'lethal': bool(
            options.lethal_gene_set and
            str(options.lethal_gene_set).lower() not in {'bundled', 'default'}
        ),
        'viable': bool(options.viable_gene_set),
    }
    if options.lethal_gene_set and str(options.lethal_gene_set).lower() in {'bundled', 'default'}:
        if species != 'fly':
            parser.error('The bundled lethal gene collection is available only for fly')
        options.lethal_gene_set = BUNDLED_FLY_LETHAL_GENES_FILE

    required_inputs = [options.gene_association_file, options.go_obo_file]
    if options.phenotype_file:
        required_inputs.append(options.phenotype_file)
    if needs_gene_sets:
        required_inputs.extend([options.lethal_gene_set, options.viable_gene_set])
    missing_inputs = [path for path in required_inputs if not os.path.isfile(path)]
    if missing_inputs:
        parser.error("Input file(s) not found: " + ", ".join(missing_inputs))
    try:
        gene_set_labels = (
            load_gene_set_labels(options.lethal_gene_set, options.viable_gene_set)
            if needs_gene_sets else None
        )
    except ValueError as exc:
        parser.error(str(exc))
    if options.min_go_gene_count < 1:
        parser.error('-min_go_gene_count must be at least 1')
    if not 0 < options.max_go_gene_fraction <= 1:
        parser.error('-max_go_gene_fraction must be in (0, 1]')
    options.go_evidence_codes = [
        code.upper() for code in (options.go_evidence_codes or [])
    ]
    options.exclude_go_evidence_codes = [
        code.upper() for code in options.exclude_go_evidence_codes
    ]
    options.worm_evidence_codes = [
        code.upper() for code in (options.worm_evidence_codes or [])
    ] or None
    evidence_overlap = set(options.go_evidence_codes) & set(options.exclude_go_evidence_codes)
    if evidence_overlap:
        parser.error(
            'GO evidence codes cannot be both included and excluded: ' +
            ', '.join(sorted(evidence_overlap))
        )

    explicit_aux = {
        'fly_assignments': options.fly_assignments is not None,
        'fly_helper_lines': options.fly_helper_lines is not None,
        'fly_lethal_genes': (
            options.fly_lethal_genes is not None and
            str(options.fly_lethal_genes).lower() != 'default'
        ),
        'worm_phenotypes': options.worm_phenotypes is not None,
        'worm_lethal_genes': (
            options.worm_lethal_genes is not None and
            str(options.worm_lethal_genes).lower() != 'default'
        ),
        'worm_viable_phenotypes': options.worm_viable_phenotypes is not None,
        'worm_viable_genes': options.worm_viable_genes is not None,
        'mouse_phenotypes': options.mouse_phenotypes is not None,
        'mouse_viable_phenotypes': options.mouse_viable_phenotypes is not None,
        'lethal_gene_set': explicit_label_sets['lethal'],
        'viable_gene_set': explicit_label_sets['viable'],
    }

    # Normalise output directory path
    try:
        options.output_dir = assert_safe_output_dir(options.output_dir, required_inputs)
    except ValueError as exc:
        parser.error(str(exc))

    # Ensure output directory exists (create if missing)
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    # Configure logging using centralised helper
    logger = configure_logger('PhenGO', enable_file=True, log_dir=options.output_dir, logfile_name='PhenGO.log', level=logging.INFO)

    logger.info(f"Processing phenotype data for species: {options.species}")
    logger.info(f"Phenotype file: {options.phenotype_file}")
    logger.info(f"Gene association file: {options.gene_association_file}")
    logger.info(f"GO OBO file: {options.go_obo_file}")

### FLY
    if species == "fly":
        ## By default we do not use the fly lethal genes file
        if options.fly_lethal_genes is None:
            pass
        elif options.fly_lethal_genes.lower() == "default":
            logger.warning(
                "Using deprecated bundled fly lethal-only override: %s",
                BUNDLED_FLY_LETHAL_GENES_FILE,
            )
            options.fly_lethal_genes = BUNDLED_FLY_LETHAL_GENES_FILE
        else:
            if not os.path.exists(options.fly_lethal_genes):
                logger.error(f"Error: Fly lethal genes file {options.fly_lethal_genes} does not exist.")
                sys.exit(1)
            else:
                logger.warning(
                    "Using deprecated lethal-only fly override: %s. Prefer paired "
                    "-label_source gene_sets or agreement inputs.",
                    options.fly_lethal_genes,
                )
        ## If the user does not provide a fly species assignments file, we use the default
        if options.fly_assignments is None:
            options.fly_assignments = DEFAULT_FLY_SPECIES_FIELDS_FILE
            logger.warning(
                "Using bundled current FlyBase assignment data. Historical analyses "
                "must supply a release-matched file and use -strict_snapshot."
            )
        else:
            if not os.path.exists(options.fly_assignments):
                logger.error(f"Error: Fly assignments file {options.fly_assignments} does not exist.")
                sys.exit(1)
        logger.info(f"Fly assignments file: {options.fly_assignments}")
        if not needs_release_labels:
            options.fly_helper_lines = None
        elif not options.fly_driver_filtering and options.fly_helper_lines is None:
            options.fly_helper_lines = None
        elif options.fly_helper_lines is None:
            options.fly_helper_lines = DEFAULT_FLY_HELPER_LINES_FILE
            logger.warning(
                "Using bundled current FlyBase helper-line data. Historical analyses "
                "must supply a release-matched file and use -strict_snapshot."
            )
        elif not os.path.exists(options.fly_helper_lines):
            logger.error(f"Error: Fly helper-lines file {options.fly_helper_lines} does not exist.")
            sys.exit(1)

### WORM
    elif species == "worm":
        if needs_release_labels:
            ## By default we do not use the worm lethal genes file
            if options.worm_lethal_genes is None:
                pass
            elif options.worm_lethal_genes.lower() == "default":
                logger.info(f"Using default worm lethal genes file: {DEFAULT_WORM_LETHAL_GENES_FILE}")
                options.worm_lethal_genes = DEFAULT_WORM_LETHAL_GENES_FILE
            else:
                if not os.path.exists(options.worm_lethal_genes):
                    logger.error(f"Error: Worm lethal genes file {options.worm_lethal_genes} does not exist.")
                    sys.exit(1)
                logger.info(f"Worm lethal genes file: {options.worm_lethal_genes}")
            ## Phenotype-term mode needs a release-matched lethal-term set. Gene-list
            ## mode does not use that file and must not silently record a current one.
            if options.worm_lethal_genes is not None:
                if options.worm_phenotypes is not None:
                    logger.warning("-worm_phenotypes is ignored when -worm_lethal_genes is used")
                options.worm_phenotypes = None
            elif options.worm_phenotypes is None:
                options.worm_phenotypes = DEFAULT_WORM_LETHAL_PHENOTYPES_FILE
                logger.info(f"Using default Worm phenotypes file: {options.worm_phenotypes}")
                logger.warning(
                    "The bundled worm phenotype-term list is current, not historical. "
                    "Supply a release-matched list and use -strict_snapshot for publication runs."
                )
            elif not os.path.exists(options.worm_phenotypes):
                logger.error(f"Error: Worm phenotypes file {options.worm_phenotypes} does not exist.")
                sys.exit(1)
            else:
                logger.info(f"Worm phenotypes file: {options.worm_phenotypes}")
            if options.worm_evidence_codes:
                logger.info(f"Worm evidence code filter: {' '.join(options.worm_evidence_codes)}")
            else:
                logger.info("Worm evidence code filter: none (all evidence codes accepted)")
            for label, path in (
                ('Worm viable phenotype', options.worm_viable_phenotypes),
                ('Worm viable gene', options.worm_viable_genes),
            ):
                if path and not os.path.isfile(path):
                    logger.error("%s file does not exist: %s", label, path)
                    sys.exit(1)

### MOUSE
    elif species == "mouse":
        if needs_release_labels:
            if options.mouse_phenotypes is None:
                options.mouse_phenotypes = DEFAULT_MOUSE_PHENOTYPES_FILE
                logger.warning(
                    "The bundled mouse phenotype-term list is current, not historical. "
                    "Supply a release-matched list and use -strict_snapshot for publication runs."
                )
            elif not os.path.exists(options.mouse_phenotypes):
                logger.error(f"Error: Mouse phenotypes file {options.mouse_phenotypes} does not exist.")
                sys.exit(1)
            logger.info(f"Mouse phenotypes file: {options.mouse_phenotypes}")
            if options.mouse_viable_phenotypes and not os.path.isfile(options.mouse_viable_phenotypes):
                logger.error(
                    "Mouse viable phenotype file does not exist: %s",
                    options.mouse_viable_phenotypes,
                )
                sys.exit(1)


    logger.info(f"Output directory: {options.output_dir}")

    if options.strict_snapshot:
        if options.fly_lethal_genes:
            logger.error(
                "Strict publication snapshots do not permit the deprecated "
                "-fly_lethal_genes override. Use paired -label_source gene_sets "
                "or agreement inputs."
            )
            sys.exit(1)
        required_metadata = [
            'snapshot_id', 'phenotype_release', 'go_annotation_release',
            'go_ontology_release', 'retrieval_date',
        ]
        if species == 'fly' or (needs_release_labels and species in {'worm', 'mouse'}):
            required_metadata.append('phenotype_ontology_release')
        missing_metadata = [name for name in required_metadata if not getattr(options, name)]
        if missing_metadata:
            logger.error("Strict snapshot mode requires: %s", ', '.join(missing_metadata))
            sys.exit(1)
        try:
            datetime.strptime(options.retrieval_date, "%Y-%m-%d")
        except ValueError:
            logger.error("Strict snapshot retrieval_date must use YYYY-MM-DD")
            sys.exit(1)
        if needs_release_labels:
            release_aux_ok = (
                species in {'yeast', 'fish'} or
                (species == 'fly' and explicit_aux['fly_assignments'] and
                 (not options.fly_driver_filtering or explicit_aux['fly_helper_lines']) and
                 (options.fly_lethal_genes is None or explicit_aux['fly_lethal_genes'])) or
                (species == 'worm' and (
                    explicit_aux['worm_lethal_genes']
                    if options.worm_lethal_genes is not None
                    else explicit_aux['worm_phenotypes']
                )) or
                (species == 'mouse' and explicit_aux['mouse_phenotypes'])
            )
        else:
            release_aux_ok = species != 'fly' or explicit_aux['fly_assignments']
        gene_set_aux_ok = (
            not needs_gene_sets or
            (explicit_aux['lethal_gene_set'] and explicit_aux['viable_gene_set'])
        )
        strict_aux_ok = release_aux_ok and gene_set_aux_ok
        if (needs_release_labels and options.nonlethal_policy == 'explicit_only'
                and species == 'worm'):
            required_viable_key = (
                'worm_viable_genes' if options.worm_lethal_genes is not None
                else 'worm_viable_phenotypes'
            )
            strict_aux_ok = strict_aux_ok and explicit_aux[required_viable_key]
        if (needs_release_labels and options.nonlethal_policy == 'explicit_only'
                and species == 'mouse'):
            strict_aux_ok = strict_aux_ok and explicit_aux['mouse_viable_phenotypes']
        if not strict_aux_ok:
            logger.error(
                "Strict snapshot mode requires explicitly supplied, release-matched auxiliary files for %s",
                species,
            )
            sys.exit(1)

    auxiliary_inputs = [
        getattr(options, name, None) for name in (
            'fly_assignments', 'fly_helper_lines', 'fly_lethal_genes',
            'worm_phenotypes', 'worm_lethal_genes', 'worm_viable_phenotypes',
            'worm_viable_genes', 'mouse_phenotypes', 'mouse_viable_phenotypes',
            'lethal_gene_set', 'viable_gene_set',
        )
    ]
    try:
        options.output_dir = prepare_output_dir(
            options.output_dir,
            options.overwrite,
            protected_paths=required_inputs + auxiliary_inputs,
            preserve_suffixes=('.log',),
        )
    except ValueError as exc:
        logger.error(str(exc))
        sys.exit(1)

    phenotype_file        = options.phenotype_file
    gene_association_file = options.gene_association_file

    logger.info("Phenotype label source: %s", options.label_source)
    release_labels = {}
    release_records = {}
    if needs_release_labels:
        release_labels = parse_release_labels(species, options, phenotype_file)
        release_records = join_labels_to_go(
            species,
            options,
            gene_association_file,
            release_labels,
            apply_legacy_fly_override=bool(options.fly_lethal_genes),
        )
    elif phenotype_file:
        logger.info(
            "Phenotype file is recorded for provenance but is not used by gene_sets mode"
        )

    fly_parse_audit_rows = list(getattr(options, "_fly_parse_audit_rows", []))
    fly_parse_audit_path = None
    fly_parse_audit_counts = {}
    if species == "fly" and options.fly_driver_filtering:
        fly_parse_audit_path = os.path.join(
            options.output_dir, "fly_compound_row_audit.tsv"
        )
        write_fly_parse_audit(fly_parse_audit_rows, fly_parse_audit_path)
        fly_parse_audit_counts = dict(sorted(Counter(
            row["outcome"] for row in fly_parse_audit_rows
        ).items()))
        logger.info("Fly compound parse audit: %s", fly_parse_audit_counts)

    gene_set_records = {}
    if needs_gene_sets:
        gene_set_records = join_labels_to_go(
            species,
            options,
            gene_association_file,
            gene_set_labels,
            apply_legacy_fly_override=False,
        )

    identifier_audit_rows = list(getattr(options, "_identifier_join_audit_rows", []))
    identifier_audit_path = os.path.join(options.output_dir, "identifier_join_audit.tsv")
    write_identifier_join_audit(identifier_audit_rows, identifier_audit_path)
    identifier_audit_counts = dict(sorted(Counter(
        row["outcome"] for row in identifier_audit_rows
    ).items()))
    logger.info("Identifier-join audit: %s", identifier_audit_counts)

    vi_inviable_genes, label_audit_rows, label_audit_counts = reconcile_label_records(
        options.label_source,
        release_records=release_records,
        gene_set_records=gene_set_records,
    )
    label_audit_path = os.path.join(options.output_dir, "label_source_audit.tsv")
    write_label_source_audit(label_audit_rows, label_audit_path)
    logger.info("Label-source audit: %s", label_audit_counts)
    if not vi_inviable_genes:
        logger.error("No genes remain after applying label source %s", options.label_source)
        sys.exit(1)

    try:
        gr, unique_go_terms, obsolete_go_terms = obo_to_graph(
            options.output_dir,
            options.go_obo_file,
            relations=options.go_relations,
            namespaces=options.go_namespaces,
            allow_cross_namespace=options.allow_cross_namespace_go_edges,
        )
    except ValueError as exc:
        logger.error("GO ontology validation failed: %s", exc)
        sys.exit(1)
    if options.exclude_go_roots:
        unique_go_terms = [term for term in unique_go_terms if term not in GO_ROOT_TERMS]

    # First validate GO terms
    vi_inviable_genes, validation_report = validate_go_terms_against_graph(vi_inviable_genes, gr, obsolete_go_terms, options.output_dir)

    # Then assign to vectors
    vi_inviable_genes, go_terms, Func = assign_go_to_vector(
        vi_inviable_genes,
        gr,
        unique_go_terms,
        propagation=options.go_propagation,
    )

    # Validation may remove every direct annotation for a gene. Do not leave
    # rows without a feature vector when unused-feature filtering is disabled.
    vi_inviable_genes = {
        gene: values for gene, values in vi_inviable_genes.items()
        if "bin_vec" in values
    }

    # Report validation results
    if validation_report['missing_terms'] > 0:
        logger.warning(
            f"WARNING: {validation_report['missing_terms']} GO terms from your gene "
            f"association file are not present in the OBO file. These terms have been "
            f"removed. See GO_term_validation_report.txt for details."
        )

    # Filter out unused GO terms (on by default; disable with -no_filter_unused_gos)
    if options.filter_unused_gos:
        vi_inviable_genes, go_terms = removed_unused_gos(vi_inviable_genes, unique_go_terms)

    vi_inviable_genes, go_terms = filter_go_terms_by_prevalence(
        vi_inviable_genes,
        go_terms,
        min_count=options.min_go_gene_count,
        max_fraction=options.max_go_gene_fraction,
    )
    if not vi_inviable_genes:
        logger.error("No genes remain after phenotype/GO validation and filtering")
        sys.exit(1)
    observed_classes = {values['status'] for values in vi_inviable_genes.values()}
    if not {'viable', 'lethal'} <= observed_classes:
        logger.error(
            "Both viable and lethal classes are required; observed %s. Under "
            "-nonlethal_policy explicit_only, supply the organism-specific viable "
            "term/gene file or choose a documented alternative policy.",
            sorted(observed_classes),
        )
        sys.exit(1)
    if not go_terms:
        logger.error("No GO features remain after relation, namespace, and prevalence filtering")
        sys.exit(1)

    if options.gene_go_pheno:
        get_o_rep_output(vi_inviable_genes, Func, f"{options.output_dir}/{options.species}_OverRepresentation.tab")

    # Always write the four go_enrichment files (binary, annotations, study set, population)
    write_func_output(vi_inviable_genes, options.output_dir, options.species)

    arff_path = f"{options.output_dir}/{options.species}_PhenGO.arff"
    write_arff_output(vi_inviable_genes, go_terms, arff_path)
    identifiers_path = os.path.join(options.output_dir, "gene_identifiers.tsv")
    write_gene_identifier_output(vi_inviable_genes, identifiers_path)

    # Save run parameters
    with open(os.path.join(options.output_dir, "PhenGO_params.txt"), "w") as outfile:
        outfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            if arg.startswith("_"):
                continue
            outfile.write(f"{arg}: {value}\n")

    input_specs = {
        "phenotype": describe_file(options.phenotype_file, options.phenotype_release),
        "go_annotations": describe_file(options.gene_association_file, options.go_annotation_release),
        "go_ontology": describe_file(options.go_obo_file, options.go_ontology_release),
        "fly_assignments": describe_file(getattr(options, 'fly_assignments', None), options.phenotype_ontology_release),
        "fly_helper_lines": describe_file(getattr(options, 'fly_helper_lines', None), options.phenotype_ontology_release),
        "fly_lethal_genes": describe_file(getattr(options, 'fly_lethal_genes', None), options.phenotype_ontology_release),
        "worm_phenotypes": describe_file(getattr(options, 'worm_phenotypes', None), options.phenotype_ontology_release),
        "worm_lethal_genes": describe_file(getattr(options, 'worm_lethal_genes', None), options.phenotype_ontology_release),
        "worm_viable_phenotypes": describe_file(getattr(options, 'worm_viable_phenotypes', None), options.phenotype_ontology_release),
        "worm_viable_genes": describe_file(getattr(options, 'worm_viable_genes', None), options.phenotype_ontology_release),
        "mouse_phenotypes": describe_file(getattr(options, 'mouse_phenotypes', None), options.phenotype_ontology_release),
        "mouse_viable_phenotypes": describe_file(getattr(options, 'mouse_viable_phenotypes', None), options.phenotype_ontology_release),
        "lethal_gene_set": describe_file(getattr(options, 'lethal_gene_set', None), options.phenotype_release),
        "viable_gene_set": describe_file(getattr(options, 'viable_gene_set', None), options.phenotype_release),
    }
    manifest = build_run_manifest(
        repo_dir=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")),
        species=species,
        snapshot_id=options.snapshot_id,
        strict_snapshot=options.strict_snapshot,
        inputs={key: value for key, value in input_specs.items() if value},
        releases={
            "phenotype": options.phenotype_release,
            "go_annotations": options.go_annotation_release,
            "go_ontology": options.go_ontology_release,
            "phenotype_ontology": options.phenotype_ontology_release,
            "retrieval_date": options.retrieval_date,
        },
        resource_context={
            "snapshot_semantics": options.snapshot_semantics,
            "phenotype_availability": options.phenotype_availability,
            "go_annotation_availability": options.go_annotation_availability,
            "go_ontology_availability": options.go_ontology_availability,
            "retrieval_route": options.retrieval_route,
        },
        policies={
            "mixed_label": options.mixed_label_policy,
            "nonlethal": options.nonlethal_policy,
            "ambiguous": options.ambiguous_label_policy,
            "label_source": options.label_source,
            "legacy_fly_lethal_override": bool(options.fly_lethal_genes),
            "filter_multigenes": options.filter_multigenes,
            "fly_driver_filtering": options.fly_driver_filtering,
            "fly_compound_policy": (
                "not_applicable"
                if species != "fly" else
                "external_paired_gene_sets"
                if options.label_source == "gene_sets" else
                "fail_closed_resolved_components_only"
                if options.fly_driver_filtering else
                "exclude_all_compound_rows" if options.filter_multigenes else
                "include_compound_rows"
            ),
            "worm_evidence_codes": options.worm_evidence_codes or "all",
            "go_relations": options.go_relations,
            "go_namespaces": options.go_namespaces or "all",
            "go_propagation": options.go_propagation,
            "allow_cross_namespace_go_edges": options.allow_cross_namespace_go_edges,
            "go_obo_header": gr.graph.get("obo_header", {}),
            "go_cross_namespace_edge_count": gr.graph.get(
                "cross_namespace_edge_count", 0
            ),
            "go_evidence_include": options.go_evidence_codes or "all",
            "go_evidence_exclude": options.exclude_go_evidence_codes,
            "identifier_join_precedence": "stable_id > unique canonical_symbol > unique synonym",
            "strict_identifier_matching": (
                "stable_id > unique canonical_symbol; synonym disabled"
                if options.strict_snapshot else
                "stable_id > unique canonical_symbol > unique synonym"
            ),
            "exclude_go_roots": options.exclude_go_roots,
            "filter_unused_gos": options.filter_unused_gos,
            "min_go_gene_count": options.min_go_gene_count,
            "max_go_gene_fraction": options.max_go_gene_fraction,
        },
        options={
            key: str(value) for key, value in vars(options).items()
            if not key.startswith("_")
        },
        outputs={
            "arff": describe_file(arff_path),
            "gene_identifiers": describe_file(identifiers_path),
            "identifier_join_audit": describe_file(identifier_audit_path),
            "label_source_audit": describe_file(label_audit_path),
            "fly_compound_row_audit": describe_file(fly_parse_audit_path),
            "go_graph": describe_file(os.path.join(options.output_dir, "GO_Children&Parents.json")),
            "go_nodes": describe_file(os.path.join(options.output_dir, "Unique_GO_Nodes.txt")),
            "go_validation": describe_file(os.path.join(options.output_dir, "GO_term_validation_report.txt")),
            "parameters": describe_file(os.path.join(options.output_dir, "PhenGO_params.txt")),
        },
        counts={
            "genes": len(vi_inviable_genes),
            "lethal": sum(v["status"] == "lethal" for v in vi_inviable_genes.values()),
            "viable": sum(v["status"] == "viable" for v in vi_inviable_genes.values()),
            "go_features": len(go_terms),
            "release_label_identifiers": len(release_labels),
            "gene_set_label_identifiers": len(gene_set_labels or {}),
            "release_labels_joined_to_go": len(release_records),
            "gene_set_labels_joined_to_go": len(gene_set_records),
            "identifier_join_lower_priority_conflict_ids": sum(
                run.get("lower_priority_conflict_ids", 0)
                for run in getattr(options, "_identifier_join_runs", [])
            ),
            "identifier_join_ambiguous_identifiers": sum(
                run.get("ambiguous_identifiers", 0)
                for run in getattr(options, "_identifier_join_runs", [])
            ),
            "identifier_join_contradictory_stable_ids": sum(
                run.get("contradictory_stable_ids", 0)
                for run in getattr(options, "_identifier_join_runs", [])
            ),
            "identifier_join_unmatched_identifiers": sum(
                run.get("unmatched_identifiers", 0)
                for run in getattr(options, "_identifier_join_runs", [])
            ),
            **{
                f"identifier_audit_{key}": value
                for key, value in identifier_audit_counts.items()
            },
            **{f"label_audit_{key}": value for key, value in label_audit_counts.items()},
            **{f"fly_parse_audit_{key}": value for key, value in fly_parse_audit_counts.items()},
        },
    )
    write_manifest(os.path.join(options.output_dir, "PhenGO_manifest.json"), manifest)

    logger.info("Thank you for using PhenGO -- A detailed user manual can be found at https://github.com/NickJD/PhenGO\n"
          "Please report any issues to: https://github.com/NickJD/PhenGO\n#####")

if __name__ == "__main__":
    main()
