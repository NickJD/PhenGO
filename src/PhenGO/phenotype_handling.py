import gzip
import csv
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try: # While calling this script through pip
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *


def get_viable_inviable_yeast(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from SGD (yeast) phenotype data.

    Args:
        options: Configuration object with filtering options
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped SGD phenotype TSV file

    Returns:
        dict: {gene_name: "lethal" or "viable"}

    Notes:
        - Column 0 contains gene name, column 1 is gene type, column 6 is mutation type and column 9 contains phenotype description
        - Uses "inviable" as equivalent to lethal in yeast nomenclature
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'mixed_phenotype': 0, 'malformed_rows': 0}

    try:
        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) <= 9:  # Ensure row has enough columns
                    continue

                gene_name = row[0]
                #gene_type = row[1]
                #mutation_type = row[6]
                phenotype_desc = row[9].lower()
                #if gene_type == "ORF" and mutation_type == "null":  - could add checks for gene_type and mutation_type
                # Check both terms in the same column for consistency
                if "inviable" in phenotype_desc or "viable" in phenotype_desc:
                    # Store the actual phenotype description
                    vi_inviable_genes[gene_name].append(row[9])

    except Exception as e:
        logger.error(f"Error reading yeast phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        # Convert to lowercase for comparison
        statuses_lower = [s.lower() for s in statuses]
        has_viable = any("viable" in s and "inviable" not in s for s in statuses_lower)
        has_inviable = any("inviable" in s for s in statuses_lower)

        if options.filter_mixed_terms and has_viable and has_inviable:
            filtered_count['mixed_phenotype'] += 1
            continue

        final_genes[gene] = "lethal" if has_inviable else "viable"

    logger.info(f"Species: Saccharomyces cerevisiae (yeast)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes


def get_viable_inviable_fly(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from FlyBase phenotype data.

    Args:
        options: Configuration object with filtering options
            - filter_multigenes: If True, exclude phenotypes with multiple genes
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped FlyBase phenotype TSV file

    Returns:
        dict: {gene_name: "lethal" or "viable"}

    Notes:
        - Handles two FlyBase formats: concatenated (Format 1) and 'with' clause (Format 2)
        - Excludes "partially lethal" phenotypes as ambiguous
        - If filter_multigenes=True, excludes any phenotype with 'with' clause or '/' in gene field
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'partially_lethal': 0, 'multi_gene': 0, 'mixed_phenotype': 0}

    try:
        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) <= 3 or row[0].startswith('#'):
                    continue

                # Filter partially lethal (ambiguous phenotypes)
                if any("partially lethal" in x.lower() for x in row):
                    filtered_count['partially_lethal'] += 1
                    continue

                # Extract base gene name (before bracket)
                base_gene = row[0].split('[')[0]

                if options.filter_multigenes:
                    # Strict filtering: exclude any complex genotypes
                    has_with_clause = any("with" in x.lower() for x in row)
                    has_compound_het = '/' in row[0]

                    if has_with_clause or has_compound_het:
                        filtered_count['multi_gene'] += 1
                        continue

                    # Classify remaining simple genotypes
                    if any("lethal" in x.lower() for x in row):
                        vi_inviable_genes[base_gene].append("lethal")
                    elif any("viable" in x.lower() for x in row):
                        vi_inviable_genes[base_gene].append("viable")
                    else:
                        vi_inviable_genes[base_gene].append("other")

                else:
                    # Permissive: include all phenotypes regardless of complexity
                    if any("lethal" in x.lower() for x in row):
                        vi_inviable_genes[base_gene].append("lethal")
                    elif any("viable" in x.lower() for x in row):
                        vi_inviable_genes[base_gene].append("viable")
                    else:
                        vi_inviable_genes[base_gene].append("other")

    except Exception as e:
        logger.error(f"Error reading fly phenotype file: {e}")
        raise

    # Resolve genes with multiple phenotype observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        has_viable = "viable" in statuses
        has_lethal = "lethal" in statuses

        if options.filter_mixed_terms and has_viable and has_lethal:
            # Exclude genes with conflicting phenotypes
            filtered_count['mixed_phenotype'] += 1
            continue

        # Priority: lethal > viable > other
        final_genes[gene] = "lethal" if has_lethal else "viable"

    # Report statistics
    logger.info(f"Species: Drosophila melanogaster (fly)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Partially lethal: {filtered_count['partially_lethal']}")
    logger.info(f"Filtered - Multi-gene: {filtered_count['multi_gene']}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes


def get_viable_inviable_fish(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from ZFIN (zebrafish) phenotype data.

    Args:
        options: Configuration object with filtering options
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped ZFIN phenotype TSV file

    Returns:
        dict: {gene_name: "lethal" or "viable"}

    Notes:
        - Assumes column 10 contains phenotype description, column 1 contains gene ID
        - Excludes "semi-lethal" and "semi-viable" phenotypes as ambiguous
        - Includes "dead" as equivalent to lethal, "alive" as equivalent to viable
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'mixed_phenotype': 0}

    try:
        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) <= 10:  # Ensure row has enough columns
                    continue

                phenotype_desc = row[10].lower()
                gene_id = row[1]

                # Classify based on phenotype keywords
                # Exclude semi-lethal and semi-viable (ambiguous)
                if ("lethal" in phenotype_desc and "semi-lethal" not in phenotype_desc) or \
                        "dead" in phenotype_desc:
                    vi_inviable_genes[gene_id].append("lethal")
                elif ("viable" in phenotype_desc and "semi-viable" not in phenotype_desc) or \
                        "alive" in phenotype_desc:
                    vi_inviable_genes[gene_id].append("viable")
                else:
                    # ISSUE: Inconsistent with fly - uses gene_id here but gene_id.split('[')[0] below
                    vi_inviable_genes[gene_id].append("other")

    except Exception as e:
        logger.error(f"Error reading fish phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        has_viable = "viable" in statuses
        has_lethal = "lethal" in statuses

        if options.filter_mixed_terms and has_viable and has_lethal:
            filtered_count['mixed_phenotype'] += 1
            continue

        final_genes[gene] = "lethal" if has_lethal else "viable"

    logger.info(f"Species: Danio rerio (zebrafish)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes


def get_viable_inviable_worm(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from WormBase phenotype data.

    Args:
        options: Configuration object with filtering options
            - worm_lethal_genes: Path to file with pre-defined lethal gene list (optional)
            - worm_phenotypes: Path to file with lethal phenotype terms (used if worm_lethal_genes not provided)
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped WormBase phenotype TSV file

    Returns:
        dict: {gene_name: "lethal" or "viable"}

    Notes:
        - Two modes: gene list mode (uses pre-curated lethal genes) or phenotype term mode
        - In phenotype term mode, matches phenotype terms against lethal term list
        - Filters out annotations with 'NOT' qualifier in phenotype term mode
        - Column 2 assumed to contain gene name, column 4 contains phenotype terms
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'mixed_phenotype': 0}

    try:
        if options.worm_lethal_genes is not None:
            # Mode 1: Use pre-defined lethal gene list
            lethal_genes = set()
            with gzip.open(options.worm_lethal_genes, 'rt', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) > 1:
                        lethal_genes.add(row[1])

            logger.info(f"Loaded {len(lethal_genes)} pre-defined lethal genes")

            with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
                reader = csv.reader(input_file, delimiter='\t')
                for row in reader:
                    if len(row) > 3:
                        gene_name = row[2]
                        # ISSUE: This assumes any gene in the lethal list is always lethal
                        # and everything else is viable, which may not be accurate
                        if gene_name in lethal_genes:
                            vi_inviable_genes[gene_name].append("lethal")
                        else:
                            vi_inviable_genes[gene_name].append("viable")
        else:
            # Mode 2: Use phenotype term matching
            lethal_terms = set()
            with gzip.open(options.worm_phenotypes, 'rt', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if row and not row[0].startswith('#'):
                        lethal_terms.add(row[0])

            logger.info(f"Loaded {len(lethal_terms)} lethal phenotype terms")

            with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
                reader = csv.reader(input_file, delimiter='\t')
                for row in reader:
                    if len(row) > 4:
                        gene_name = row[2]
                        phenotype_terms = row[4]
                        qualifier = row[3] if len(row) > 3 else ""

                        # Match against lethal terms, exclude NOT qualifiers
                        if any(term in phenotype_terms for term in lethal_terms) and \
                                'NOT' not in qualifier:
                            vi_inviable_genes[gene_name].append("lethal")
                        else:
                            vi_inviable_genes[gene_name].append("viable")

    except Exception as e:
        logger.error(f"Error reading worm phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        has_viable = "viable" in statuses
        has_lethal = "lethal" in statuses

        if options.filter_mixed_terms and has_viable and has_lethal:
            filtered_count['mixed_phenotype'] += 1
            continue

        final_genes[gene] = "lethal" if has_lethal else "viable"

    logger.info(f"Species: Caenorhabditis elegans (worm)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes


def get_viable_inviable_mouse(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from MGI (mouse) phenotype data.

    Args:
        options: Configuration object with filtering options
            - mouse_phenotypes: Path to file with lethal phenotype terms
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped MGI phenotype TSV file

    Returns:
        dict: {gene_name: "lethal" or "viable"}

    Notes:
        - Column 3 contains phenotype terms, column 5 contains gene name
        - Filters transgenes by excluding entries with commas in gene name
        - Matches phenotype terms against lethal term list
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'transgenes': 0, 'mixed_phenotype': 0}

    try:
        # Load lethal phenotype terms
        lethal_terms = set()
        with gzip.open(options.mouse_phenotypes, 'rt', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row and row[0] != 'ID':
                    lethal_terms.add(row[0])

        logger.info(f"Loaded {len(lethal_terms)} lethal phenotype terms")

        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            for row in reader:
                if len(row) > 5:
                    gene_name = row[5]
                    phenotype_terms = row[3]

                    # Filter transgenes (simplistic check for multiple genes)
                    # ISSUE: Could miss other multi-gene cases or incorrectly filter valid genes
                    if ',' in gene_name:
                        filtered_count['transgenes'] += 1
                        continue

                    # Classify based on phenotype term matching
                    if any(term in phenotype_terms for term in lethal_terms):
                        vi_inviable_genes[gene_name].append("lethal")
                    else:
                        vi_inviable_genes[gene_name].append("viable")

    except Exception as e:
        logger.error(f"Error reading mouse phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        has_viable = "viable" in statuses
        has_lethal = "lethal" in statuses

        if options.filter_mixed_terms and has_viable and has_lethal:
            filtered_count['mixed_phenotype'] += 1
            continue

        final_genes[gene] = "lethal" if has_lethal else "viable"

    logger.info(f"Species: Mus musculus (mouse)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Transgenes: {filtered_count['transgenes']}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes