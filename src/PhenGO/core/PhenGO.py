import argparse
import sys
import shutil
import os
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

from .obo_to_graph import obo_to_graph
from ..constants import *
from .phenotype_handling import *
from .go_handling import *

def removed_unused_gos(vi_inviable_genes, unique_go_terms):
    # This optional function removes unused GO terms from the ARFF and go_enrichment outputs
    # It rebuilds the GO term list and each gene's bin_vec to only include used terms
    # collect *all* GOs assigned to at least 1 gene
    used_gos = set()
    for gene_vals in vi_inviable_genes.values():
        used_gos.update(gene_vals.get("go_list", []))
    # restrict the master list to only the used terms, in order
    filtered_go_terms = [go for go in unique_go_terms if go in used_gos]
    # now rebuild each gene's go_list AND its bin_vec *only* over filtered_go_terms
    filtered_genes = {}
    for gene_id, gene_vals in vi_inviable_genes.items():
        old_go_list = gene_vals.get("go_list", [])
        if not old_go_list:
            continue
        old_go_list = list(dict.fromkeys(old_go_list))
        new_go_list = [go for go in filtered_go_terms if go in old_go_list]
        new_bin_vec = [1 if go in new_go_list else 0 for go in filtered_go_terms]
        filtered_genes[gene_id] = {
            "status" : gene_vals["status"],
            "go_list": new_go_list,
            "bin_vec" : new_bin_vec
        }
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

    # Process each gene
    for gene, values in vi_inviable_genes.items():
        if 'go_list' not in values:
            continue

        cleaned_go_list = []

        for go_term in values['go_list']:
            total_terms_checked += 1

            # Check if term is obsolete
            if go_term in obsolete_go_terms:
                obsolete_terms_found[go_term].append(gene)
                continue

            # Check if term exists in graph
            if go_term not in gr:
                missing_terms[go_term].append(gene)
                continue

            # Valid term
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

    # Console output
    logger.info(f"GO term validation complete:")
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


def assign_go_to_vector(vi_inviable_genes, gr, unique_go_terms):
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

            # Get all ancestors (parents, grandparents, etc.)
            try:
                ancestors = nx.ancestors(gr, go_term)
                all_ancestors.update(ancestors)
            except nx.NetworkXError as e:
                logger.warning(f"Could not get ancestors for {go_term} in gene {gene}: {e}")
                continue

        # Set binary vector and record overrepresentation output entries (one write per unique term)
        for ancestor in all_ancestors:
            idx = go_term_index.get(ancestor)
            if idx is not None:
                bin_vec[idx] = 1
                o_rep.append(f"{gene}\t{ancestor}\n")
            else:
                logger.debug(f"Ancestor {ancestor} not in unique_go_terms list")

        vi_inviable_genes[gene]["bin_vec"] = bin_vec.copy()
        o_rep.append('\n')
        bin_vec = [0] * len(unique_go_terms)

    logger.info(f"Vector assignment complete for {genes_processed} genes")

    return vi_inviable_genes, unique_go_terms, o_rep


def get_o_rep_output(vi_inviable_genes, o_rep, output_file):
    # Generate overrepresentation analysis output file
    new_o_rep = []
    geneSeen = set()
    tempySeen = set()
    Counter = 0

    for gene, values in vi_inviable_genes.items():
        if values['status'] in ['inviable', 'lethal']:  # Accept both for cross-species compatibility
            for line in o_rep:
                if line == "\n":
                    continue
                line_gene = line.split('\t')[0]
                if line_gene == gene and line not in geneSeen:
                    geneSeen.add(line)
                    entry = line.strip() + "\t1"
                    new_o_rep.append([entry])
                if line_gene == gene and gene not in tempySeen:
                    tempySeen.add(gene)
                    Counter += 1

        elif values['status'] == 'viable':
            for line in o_rep:
                if line == "\n":
                    continue
                line_gene = line.split('\t')[0]
                if line_gene == gene and line not in geneSeen:
                    geneSeen.add(line)
                    entry = line.strip() + "\t0"
                    new_o_rep.append([entry])
                if line_gene == gene and gene not in tempySeen:
                    tempySeen.add(gene)
                    Counter += 1

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
        go_list = values.get("go_list", [])
        if not go_list:
            continue

        status = values.get("status", "viable")
        label  = "1" if status in ("lethal", "inviable") else "0"

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
    with open(output_file, 'w') as f:
        f.write("@RELATION gene_lethality\n\n")
        f.write(f"@ATTRIBUTE gene string\n")
        for go_term in filtered_go_terms:
            f.write(f"@ATTRIBUTE {go_term} {{0,1}}\n")
        f.write("@ATTRIBUTE class {viable,lethal}\n\n")
        f.write("@DATA\n")
        for gene, values in vi_inviable_genes.items():
            gene = gene.replace("'", "-")  # Replace single quotes to keep ARFF valid
            bin_vec = ",".join(map(str, values["bin_vec"]))
            status = values["status"]
            if status == 'inviable':
                status = 'lethal'
            f.write(f"{gene},{bin_vec},{status}\n")


def main():
    parser = argparse.ArgumentParser(description=f"PhenGO {PhenGO_VERSION} - Convert phenotype and GO data to ARFF format")
    parser.add_argument("--print-defaults", action="store_true",
                        help="Print default files and methods used for each species and exit.")

    # Parse known args first to check for --print-defaults
    args, remaining_argv = parser.parse_known_args()
    if args.print_defaults:
        configure_logger('PhenGO', enable_file=False)
        logger = logging.getLogger('PhenGO')
        logger.info("\nDefault options and methods for each species:\n")
        logger.info("Fly:")
        logger.info(f"  Assignments file: {DEFAULT_FLY_SPECIES_FIELDS_FILE}")
        logger.info(f"  Lethal genes file: {DEFAULT_FLY_LETHAL_GENES_FILE}")
        logger.info("  Method: get_viable_inviable_fly, get_viability_go_data_fly")
        logger.info("Worm:")
        logger.info(f"  Lethal phenotypes file: {DEFAULT_WORM_LETHAL_PHENOTYPES_FILE}")
        logger.info(f"  Lethal genes file: {DEFAULT_WORM_LETHAL_GENES_FILE}")
        logger.info("  Method: By default we do not use the lethal genes file, we use the default lethal phenotypes file"
                    " to identify the genes in the user-provided phenotype_data file to identify lethal genes. If lethal "
                    "genes file is provided (or 'default' is given, it will be used to identify lethal genes instead.")
        logger.info("Mouse:")
        logger.info(f"  Phenotypes file: {DEFAULT_MOUSE_PHENOTYPES_FILE}")
        logger.info("  Method: get_viable_inviable_mouse, get_viability_go_data_mouse")
        logger.info("Yeast:")
        logger.info("  Method: get_viable_inviable_yeast, get_viability_go_data_yeast")
        logger.info("Fish:")
        logger.info("  Method: get_viable_inviable_fish, get_viability_go_data_fish")
        sys.exit(0)


    parser._action_groups.pop()
    required = parser.add_argument_group('Required Options')
    required.add_argument('-species', dest="species", required=True, help='Species tag (e.g., fly, yeast)')
    required.add_argument('-phenotype_file', dest="phenotype_file", required=True, help='Path to the phenotype data file (.gz)')
    required.add_argument('-gene_association_file', dest="gene_association_file", required=True, help='Path to the gene association file (.gz)')
    required.add_argument('-go_obo_file', dest="go_obo_file", required=True, help='Path to the go.obo file')
    required.add_argument('-output_dir', dest="output_dir", required=True, help='Output directory')

    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-no_filter_unused_gos', dest='filter_unused_gos', action='store_false', required=False,
        help='Disable filtering of unused GO terms from the go_enrichment and ARFF output (filtering is ON by default)')
    optional.add_argument('-filter_mixed_terms', dest='filter_mixed_terms', action='store_true', required=False,
                        help='Filter out genes which have both lethal and viable phenotypes - '
                             'Terms not specifically lethal/viable are not counted in this (default: False)')
    optional.add_argument('-gene_go_pheno', dest='gene_go_pheno', action='store_true', required=False,
                        help='(Deprecated — go_enrichment files are now always written. Kept for backward compatibility.)'
                             ' Previously output a single Gene-GO-Phenotype file; that file is now'
                             ' superseded by the four {species}_*.tab/txt files in go_enrichment/ generated every run.')

    fly_args = parser.add_argument_group('Fly specific parameters')
    fly_args.add_argument('-fly_lethal_genes', dest='fly_lethal_genes', required=False,
                        help='Provide TSV file of specified lethal fly genes (provide "default" for default lethal genes: "data/fly/FlyBase_Lethal_Gene_IDs_2025-08-15.txt.gz")')
    fly_args.add_argument('-fly_assignments', dest='fly_assignments', required=False,
                        help='Provide TSV file of fly assignments (file confirming genes are assigned to drosophila melanogaster (default: "data/fly/FlyBase_Fields_2025_07_29.tsv.gz")')
    fly_args.add_argument('-filter_multigenes', dest='filter_multigenes', action='store_true', required=False,
                        help='Filter out phenotypes involving multiple genes: rows with a "with" tag or "/" in '
                             'the gene field are discarded. In simple mode (default) this is a broad string match; '
                             'in smart mode (-fly_driver_filtering) only genuine second-gene disruptions are filtered.')
    fly_args.add_argument('-fly_driver_filtering', dest='fly_driver_filtering', action='store_true', required=False,
                        help='Enable smart driver-line filtering for fly phenotype parsing. '
                             'Driver lines (GAL4, UAS, lexA, QUAS, etc.) are detected via regex and excluded '
                             'from multi-gene counting, so "viable, with Scer\\GAL4[x]" is correctly kept as '
                             'viable rather than discarded. Non-descript phenotypes are treated as viable. '
                             'Default: OFF (uses simpler broad-match filtering).')

    worm_args = parser.add_argument_group('Worm specific parameters')
    worm_args.add_argument('-worm_phenotypes', dest='worm_phenotypes', required=False,
                        help='Provide TSV file of worm phenotypes (default: "data/worm/lethal_terms_traversed_2025-08-12.tsv.gz")')
    worm_args.add_argument('-worm_lethal_genes', dest='worm_lethal_genes', required=False,
                        help='Provide TSV file of specified lethal worm genes (provide "default" for default lethal genes: "data/worm/genes_direct_and_inferred_for_WBPhenotype_0000062_11-08-2025.txt.gz")')
    worm_args.add_argument('-worm_evidence_codes', dest='worm_evidence_codes', nargs='+', required=False,
                        help='Only retain worm phenotype annotations with these evidence codes '
                             '(e.g. IMP RNAi). By default all evidence codes are accepted. '
                             'Example: -worm_evidence_codes IMP IBA')

    mouse_args = parser.add_argument_group('Mouse specific parameters')
    mouse_args.add_argument('-mouse_phenotypes', dest='mouse_phenotypes', required=False,
                        help='Provide TSV file of mouse phenotypes (default: "data/mouse/mouse_lethal_terms.txt.gz")')

    misc = parser.add_argument_group('Misc')
    misc.add_argument("-v", "--version", action="version",
                 version=f"PhenGO {PhenGO_VERSION}: Exiting.")


    options = parser.parse_args()

    # Normalise output directory path
    options.output_dir = os.path.abspath(options.output_dir)

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
    if options.species.lower() == "fly":
        ## By default we do not use the fly lethal genes file
        if options.fly_lethal_genes is None:
            pass
        elif options.fly_lethal_genes.lower() == "default":
            logger.info(f"Using default fly lethal genes file: {DEFAULT_FLY_LETHAL_GENES_FILE}")
            options.fly_lethal_genes = DEFAULT_FLY_LETHAL_GENES_FILE
        else:
            if not os.path.exists(options.fly_lethal_genes):
                logger.error(f"Error: Fly lethal genes file {options.fly_lethal_genes} does not exist.")
                sys.exit(1)
            else:
                logger.info(f"Fly lethal genes file: {options.fly_lethal_genes}")
        ## If the user does not provide a fly species assignments file, we use the default
        if options.fly_assignments is None:
            options.fly_assignments = DEFAULT_FLY_SPECIES_FIELDS_FILE
        else:
            if not os.path.exists(options.fly_assignments):
                logger.error(f"Error: Fly assignments file {options.fly_assignments} does not exist.")
                sys.exit(1)
        logger.info(f"Fly assignments file: {options.fly_assignments}")

### WORM
    elif options.species.lower() == "worm":
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
            else:
                logger.info(f"Worm lethal genes file: {options.worm_lethal_genes}")
        ## If the user does not provide a worm phenotypes file, we use the default
        if options.worm_phenotypes is None:
            options.worm_phenotypes = DEFAULT_WORM_LETHAL_PHENOTYPES_FILE
            logger.info(f"Using default Worm phenotypes file: {options.worm_phenotypes}")
        else:
            if not os.path.exists(options.worm_phenotypes):
                logger.error(f"Error: Worm phenotypes file {options.worm_phenotypes} does not exist.")
                sys.exit(1)
            else:
                logger.info(f"Worm phenotypes file: {options.worm_phenotypes}")
        if options.worm_evidence_codes:
            logger.info(f"Worm evidence code filter: {' '.join(options.worm_evidence_codes)}")
        else:
            logger.info("Worm evidence code filter: none (all evidence codes accepted)")

### MOUSE
    elif options.species.lower() == "mouse":
        if options.mouse_phenotypes is None:
            options.mouse_phenotypes = DEFAULT_MOUSE_PHENOTYPES_FILE
        else:
            if not os.path.exists(options.mouse_phenotypes):
                logger.error(f"Error: Mouse phenotypes file {options.mouse_phenotypes} does not exist.")
                sys.exit(1)
        logger.info(f"Mouse phenotypes file: {options.mouse_phenotypes}")


    logger.info(f"Output directory: {options.output_dir}")

    # Ensure output directory is empty (preserve existing log files)
    if os.path.exists(options.output_dir):
        non_log_entries = [
            e for e in os.listdir(options.output_dir)
            if not (os.path.isfile(os.path.join(options.output_dir, e)) and e.endswith('.log'))
        ]
        if non_log_entries:
            logger.warning(
                f"Output directory '{options.output_dir}' already contains {len(non_log_entries)} "
                f"file(s)/folder(s). They will be deleted before this run proceeds. "
                f"Pass a new -output_dir to preserve them."
            )
        for entry in os.listdir(options.output_dir):
            path = os.path.join(options.output_dir, entry)
            # Preserve .log files (including PhenGO.log)
            if os.path.isfile(path) and entry.endswith('.log'):
                continue
            try:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
            except Exception:
                logger.warning(f"Could not remove {path} when clearing output directory")
    else:
        os.makedirs(options.output_dir)

    phenotype_file        = options.phenotype_file
    gene_association_file = options.gene_association_file

    species = options.species.lower()
    if species == "yeast":
        vi_inviable_genes = get_viable_inviable_yeast(options, phenotype_file)
    elif species == "fly":
        vi_inviable_genes = get_viable_inviable_fly(options, phenotype_file)
    elif species == "fish":
        vi_inviable_genes = get_viable_inviable_fish(options, phenotype_file)
    elif species == "worm":
        vi_inviable_genes = get_viable_inviable_worm(options, phenotype_file)
    elif species == "mouse":
        vi_inviable_genes = get_viable_inviable_mouse(options, phenotype_file)
    else:
        logger.error(f"Unknown species '{options.species}'. Supported: fly, yeast, fish, worm, mouse")
        sys.exit(1)

    if species == "yeast":
        vi_inviable_genes = get_viability_go_data_yeast(gene_association_file, vi_inviable_genes)
    elif species == "fly":
        vi_inviable_genes = get_viability_go_data_fly(options, gene_association_file, vi_inviable_genes)
    elif species == "fish":
        vi_inviable_genes = get_viability_go_data_fish(gene_association_file, vi_inviable_genes)
    elif species == "worm":
        vi_inviable_genes = get_viability_go_data_worm(gene_association_file, vi_inviable_genes)
    elif species == "mouse":
        vi_inviable_genes = get_viability_go_data_mouse(gene_association_file, vi_inviable_genes)

    gr, unique_go_terms, obsolete_go_terms = obo_to_graph(options.output_dir, options.go_obo_file)

    # First validate GO terms
    vi_inviable_genes, validation_report = validate_go_terms_against_graph(vi_inviable_genes, gr, obsolete_go_terms, options.output_dir)

    # Then assign to vectors
    vi_inviable_genes, go_terms, Func = assign_go_to_vector(vi_inviable_genes, gr, unique_go_terms)

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

    if options.gene_go_pheno:
        get_o_rep_output(vi_inviable_genes, Func, f"{options.output_dir}/{options.species}_OverRepresentation.tab")

    # Always write the four go_enrichment files (binary, annotations, study set, population)
    write_func_output(vi_inviable_genes, options.output_dir, options.species)

    write_arff_output(vi_inviable_genes, go_terms, f"{options.output_dir}/{options.species}_PhenGO.arff")

    # Save run parameters
    with open(os.path.join(options.output_dir, "PhenGO_params.txt"), "w") as outfile:
        outfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")

    logger.info("Thank you for using PhenGO -- A detailed user manual can be found at https://github.com/NickJD/PhenGO\n"
          "Please report any issues to: https://github.com/NickJD/PhenGO\n#####")

if __name__ == "__main__":
    main()







