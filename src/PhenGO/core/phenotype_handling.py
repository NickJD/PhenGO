import gzip
import csv
import re
from collections import defaultdict
import logging
from .. import constants as _const

logger = logging.getLogger(__name__)

DEFAULT_FLY_HELPER_LINES_FILE = getattr(_const, 'DEFAULT_FLY_HELPER_LINES_FILE', None)



def get_viable_inviable_yeast(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from SGD (yeast) phenotype data.
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
                gene_type = row[1]
                mutation_type = row[6]
                phenotype_desc = row[9].lower()
                if gene_type == "ORF" and mutation_type == "null":
                    # Any observed phenotype for a null ORF is recorded.
                    # Non-inviable observations (including non-descript ones) are
                    # treated as viable — the organism survived the experiment.
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


def get_viable_inviable_worm(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from WormBase phenotype data.
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
                        # This assumes any gene in the lethal list is always lethal
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
                        qualifier = row[3] if len(row) > 3 else ""
                        phenotype_terms = row[4]
                        evidence_codes = row[6] if len(row) > 6 else ""

                        # Match against lethal terms, exclude NOT qualifiers.
                        # Evidence code filtering is only applied when the user
                        # explicitly provides codes via -worm_evidence_codes.
                        evidence_codes_filter = getattr(options, 'worm_evidence_codes', None)
                        evidence_ok = (
                            any(ec in evidence_codes for ec in evidence_codes_filter)
                            if evidence_codes_filter
                            else True
                        )
                        if any(term in phenotype_terms for term in lethal_terms) and \
                                'NOT' not in qualifier and evidence_ok:
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


def get_viable_inviable_fly(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from FlyBase phenotype data.

    Two modes are available, controlled by options.fly_driver_filtering:

    Simple mode (default, fly_driver_filtering=False):
        - Scans all columns for the word "lethal" → lethal; anything else → viable.
        - If filter_multigenes=True, rows with any "with" in any column or "/" in the
          gene field are discarded entirely (old conservative behaviour).
        - Any observed row that passes filters and is not lethal is treated as viable,
          including non-descript phenotypes (neuroanatomy, behaviour, etc.).

    Smart / driver-filtering mode (fly_driver_filtering=True):
        - Handles two FlyBase column formats (4-col Format 1, 7-col Format 2).
        - Driver lines (GAL4, UAS, lexA, QUAS, etc.) are detected via regex patterns
          and an optional helper-lines file.  "viable, with Scer\\GAL4[x]" is stripped
          to just "viable" — the driver is not counted as a second gene.
        - If filter_multigenes=True, only rows where a non-driver, non-self gene
          disruption is present in the 'with' clause are excluded.
        - Phenotype classification checks for the words ``lethal'' and ``viable''
          in the phenotype description; anything else is treated as viable —
          if the fly was observed with a morphological or behavioural phenotype,
          it survived, so the gene is de-facto viable.
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'partially_lethal': 0, 'multi_gene': 0, 'mixed_phenotype': 0}

    use_smart = getattr(options, 'fly_driver_filtering', False)

    # ── Driver-line helpers (smart mode only) ─────────────────────────────────
    if use_smart:
        helper_lines_set = set()
        if DEFAULT_FLY_HELPER_LINES_FILE:
            try:
                with gzip.open(DEFAULT_FLY_HELPER_LINES_FILE, 'rt', encoding='utf-8') as f:
                    reader_h = csv.reader(f, delimiter='\t')
                    helper_lines_set = {row[2] for row in reader_h if len(row) > 2}
            except Exception as e:
                logger.warning(f"Could not load helper lines file: {e}")

        driver_patterns = [
            r'Scer\\GAL4\[', r'Scer\\UAS\.', r'Ecol\\lexA\[', r'Ecol\\lexAop\.',
            r'Ncra\\QF\[', r'Ncra\\QUAS\.', r'Scer\\FLP\d*\[', r'Scer\\GAL80\[',
            r'tetO_', r'rtTA', r'tTA',
        ]
        disruption_patterns = [
            r'\[GD\d+\]', r'\[KK\d+\]', r'\[HMS\d+\]', r'\[JF\d+\]', r'\[GL\d+\]',
            r'\[ex\d+\]', r'\[Δ\d+\]', r'\[del\w*\]', r'\[df\w*\]', r'\[null\]',
            r'\[\w*-\]', r'\[\w*\^\w*\]', r'\[MB\d+\]', r'\[PB\w*\]', r'\[P\w*\]',
            r'\[CR\d+\]', r'\[cas\d+\]',
        ]

        def is_driver_line(component):
            if component in helper_lines_set:
                return True
            return any(re.search(p, component, re.IGNORECASE) for p in driver_patterns)

        def extract_gene_name(allele_string):
            if not allele_string:
                return None
            cleaned = re.sub(r'^(Scer\\|Ecol\\|Ncra\\)', '', str(allele_string))
            match = re.match(r'^([^[\\\s]+)', cleaned)
            return match.group(1) if match else cleaned.split('[')[0].split('\\')[0].strip()

        def parse_with_clause(phenotype_field):
            if not phenotype_field or 'with ' not in phenotype_field:
                return phenotype_field, []
            if ', with ' in phenotype_field:
                clean, with_part = phenotype_field.split(', with ', 1)
            elif ' with ' in phenotype_field:
                clean, with_part = phenotype_field.split(' with ', 1)
            else:
                return phenotype_field, []
            return clean.strip(), [c.strip() for c in with_part.split(',') if c.strip()]

        def classify_phenotype(phenotype_desc):
            if not phenotype_desc:
                return 'other'
            lower = str(phenotype_desc).lower()
            if 'partially' in lower:
                return 'other'
            if 'lethal' in lower:
                return 'lethal'
            if 'viable' in lower:
                return 'viable'
            return 'other'

        def has_real_second_gene(primary_gene, with_components):
            for component in with_components:
                if is_driver_line(component):
                    continue
                if any(re.search(p, component) for p in disruption_patterns):
                    gene = extract_gene_name(component)
                    if gene and gene != primary_gene:
                        return True
            return False

    # ── Parse file ─────────────────────────────────────────────────────────────
    try:
        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue

                if use_smart:
                    # ── Smart mode ─────────────────────────────────────────────
                    is_format1 = len(row) == 4
                    is_format2 = len(row) == 7
                    if not (is_format1 or is_format2):
                        continue

                    if is_format1:
                        primary_gene_allele = row[0]
                        clean_phenotype, with_components = parse_with_clause(row[2])
                    else:
                        primary_gene_allele = row[0].split('/')[0] if '/' in row[0] else row[0]
                        secondary = row[0].split('/')[1:] if '/' in row[0] else []
                        clean_phenotype = row[2]
                        with_components = [s.strip() for s in secondary if s.strip()]

                    primary_gene = extract_gene_name(primary_gene_allele)
                    if not primary_gene or is_driver_line(primary_gene_allele):
                        continue

                    if 'partially' in str(clean_phenotype).lower():
                        filtered_count['partially_lethal'] += 1
                        continue

                    if options.filter_multigenes and has_real_second_gene(primary_gene, with_components):
                        filtered_count['multi_gene'] += 1
                        continue

                    # "other" (non-descript) phenotype → viable: the fly survived
                    phenotype_class = classify_phenotype(clean_phenotype)
                    vi_inviable_genes[primary_gene].append(
                        phenotype_class if phenotype_class != 'other' else 'viable'
                    )

                else:
                    # ── Simple mode (default) ──────────────────────────────────
                    if len(row) <= 3:
                        continue

                    if any("partially lethal" in x.lower() for x in row):
                        filtered_count['partially_lethal'] += 1
                        continue

                    base_gene = row[0].split('[')[0]

                    if options.filter_multigenes:
                        has_with_clause = any("with" in x.lower() for x in row)
                        has_compound_het = '/' in row[0]
                        if has_with_clause or has_compound_het:
                            filtered_count['multi_gene'] += 1
                            continue

                    # Any observed phenotype that is not lethal → viable
                    if any("lethal" in x.lower() for x in row):
                        vi_inviable_genes[base_gene].append("lethal")
                    else:
                        vi_inviable_genes[base_gene].append("viable")

    except Exception as e:
        logger.error(f"Error reading fly phenotype file: {e}")
        raise

    # ── Resolve multiple observations ──────────────────────────────────────────
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        has_viable = "viable" in statuses
        has_lethal = "lethal" in statuses

        if options.filter_mixed_terms and has_viable and has_lethal:
            filtered_count['mixed_phenotype'] += 1
            continue

        final_genes[gene] = "lethal" if has_lethal else "viable"

    mode_label = "smart (driver-filtering)" if use_smart else "simple"
    logger.info(f"Species: Drosophila melanogaster (fly)  [mode: {mode_label}]")
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
                else:
                    # Any non-lethal phenotype (explicit "viable"/"alive" OR morphological
                    # phenotypes like "abnormal fin morphology") is treated as viable.
                    # If the gene was studied in ZFIN and the organism survived with only
                    # non-lethal phenotypes, the gene is de-facto viable.
                    # Semi-lethal/semi-viable are included here as they do not indicate full lethality.
                    vi_inviable_genes[gene_id].append("viable")

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
        - Supports MGI_PhenoGenoMP and MGI_GenePheno-style tabular reports
        - Uses MGI marker accession IDs where available to match GAF column 2
        - Filters multi-marker/transgene rows by excluding entries with commas in gene accession
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

        def extract_mouse_gene_and_mp(row):
            # MGI_GenePheno: allele composition, allele symbol, allele id,
            # background, MP id, PubMed id, marker accession, marker symbol, ...
            if len(row) > 6 and row[4].startswith("MP:"):
                return row[6], row[4]
            # MGI_PhenoGenoMP: allele composition, allele symbol, background,
            # MP id, PubMed id, marker accession, genotype id
            if len(row) > 5 and row[3].startswith("MP:"):
                return row[5], row[3]
            # Legacy/custom fallback used by earlier PhenGO versions.
            if len(row) > 4:
                return row[2], row[4]
            return None, None

        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            for row in reader:
                gene_name, phenotype_term = extract_mouse_gene_and_mp(row)
                if not gene_name or not phenotype_term:
                    continue

                if ',' in gene_name:
                    filtered_count['transgenes'] += 1
                    continue

                if phenotype_term in lethal_terms:
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
