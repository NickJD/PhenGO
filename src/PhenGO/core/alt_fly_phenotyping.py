import gzip
import csv
import re
from collections import defaultdict
import logging
from .. import constants as _const
# Get defaults from constants if available
DEFAULT_FLY_HELPER_LINES_FILE = getattr(_const, 'DEFAULT_FLY_HELPER_LINES_FILE', None)
configure_logger = _const.configure_logger


def get_viable_inviable_fly(options, phenotype_file):
    """
    Simplified function to classify fly genes as viable/inviable using column count detection.
    """
    vi_inviable_genes = defaultdict(list)

    # Load helper lines (driver genes) if a helper-lines file is configured
    helper_lines_set = set()
    if DEFAULT_FLY_HELPER_LINES_FILE:
        try:
            with gzip.open(DEFAULT_FLY_HELPER_LINES_FILE, 'rt', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                helper_lines_set = {row[2] for row in reader if len(row) > 2}
        except Exception as e:
            # Configure logger for this module (console only)
            configure_logger('PhenGO.alt_fly_phenotyping', enable_file=False)
            logger = logging.getLogger('PhenGO.alt_fly_phenotyping')
            logger.warning(f"Could not load helper lines file: {e}")

    # Driver line patterns for better detection
    driver_patterns = [
        r'Scer\\GAL4\[',
        r'Scer\\UAS\.',
        r'Ecol\\lexA\[',
        r'Ecol\\lexAop\.',
        r'Ncra\\QF\[',
        r'Ncra\\QUAS\.',
        r'Scer\\FLP\d*\[',
        r'Scer\\GAL80\[',
        r'tetO_',
        r'rtTA',
        r'tTA'
    ]

    def is_driver_line(component):
        """Check if component is a driver line."""
        if component in helper_lines_set:
            return True
        for pattern in driver_patterns:
            if re.search(pattern, component, re.IGNORECASE):
                return True
        return False

    def extract_gene_name(allele_string):
        """Extract gene name from allele string."""
        if not allele_string:
            return None
        # Remove common prefixes and extract gene name before bracket
        cleaned = re.sub(r'^(Scer\\|Ecol\\|Ncra\\)', '', str(allele_string))
        match = re.match(r'^([^[\\\s]+)', cleaned)
        if match:
            return match.group(1)
        return cleaned.split('[')[0].split('\\')[0].strip()

    def parse_with_clause_from_phenotype(phenotype_field):
        """Extract 'with' components from phenotype field (Format 1)."""
        if not phenotype_field or 'with ' not in phenotype_field:
            return phenotype_field, []

        # Split on ', with ' or ' with '
        if ', with ' in phenotype_field:
            clean_phenotype, with_part = phenotype_field.split(', with ', 1)
        elif ' with ' in phenotype_field:
            clean_phenotype, with_part = phenotype_field.split(' with ', 1)
        else:
            return phenotype_field, []

        # Split with components on commas
        with_components = [c.strip() for c in with_part.split(',') if c.strip()]

        return clean_phenotype.strip(), with_components

    def classify_phenotype(phenotype_desc):
        """Classify phenotype as lethal, viable, or other."""
        if not phenotype_desc:
            return 'other'

        phenotype_lower = str(phenotype_desc).lower()

        # Skip partially lethal/viable
        if 'partially' in phenotype_lower:
            return 'other'

        # Lethal keywords (higher priority)
        lethal_keywords = [
            'lethal', 'dies', 'death', 'dead', 'embryonic lethal', 'larval lethal',
            'pupal lethal', 'adult lethal', 'essential',
            'no eclosion', 'eclosion fail', 'fail to eclose', 'cannot survive',
            'inviable', 'non-viable'
        ]

        for keyword in lethal_keywords:
            if keyword in phenotype_lower:
                return 'lethal'

        # Viable keywords
        viable_keywords = [
            'viable', 'survives', 'normal viability', 'fertile', 'eclosion normal',
            'develops normally', 'no lethality'
        ]

        for keyword in viable_keywords:
            if keyword in phenotype_lower:
                return 'viable'

        return 'other'

    def has_multiple_gene_disruptions(primary_gene, with_components):
        """Check if this involves multiple gene disruptions (multi-gene effect)."""
        # Gene disruption patterns
        disruption_patterns = [
            r'\[GD\d+\]', r'\[KK\d+\]', r'\[HMS\d+\]', r'\[JF\d+\]', r'\[GL\d+\]',
            r'\[ex\d+\]', r'\[Δ\d+\]', r'\[del\w*\]', r'\[df\w*\]', r'\[null\]',
            r'\[\w*-\]', r'\[\w*\^\w*\]', r'\[MB\d+\]', r'\[PB\w*\]', r'\[P\w*\]',
            r'\[CR\d+\]', r'\[cas\d+\]'
        ]

        def is_gene_disruption(component):
            for pattern in disruption_patterns:
                if re.search(pattern, component):
                    return True
            return False

        # Count gene disruptions in with components
        other_gene_disruptions = 0
        for component in with_components:
            if is_gene_disruption(component) and not is_driver_line(component):
                gene_name = extract_gene_name(component)
                if gene_name and gene_name != primary_gene:
                    other_gene_disruptions += 1

        return other_gene_disruptions > 0

    # Process the phenotype file
    with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
        reader = csv.reader(input_file, delimiter='\t')

        for row in reader:
            if not row or row[0].startswith('#'):
                continue

            if 'HT1B[KK112342]' in str(row):
                logger = logging.getLogger('PhenGO.alt_fly_phenotyping')
                logger.debug(f"Debugging row with HT1B: {row}")
            try:
                # Determine format by column count
                is_format1 = len(row) == 4
                is_format2 = len(row) == 7

                if not (is_format1 or is_format2):
                    continue

                # Extract primary gene and phenotype based on format
                if is_format1:
                    # Format 1: 4 columns
                    # Column 0: primary gene allele
                    # Column 2: phenotype (may contain 'with' clause)
                    primary_gene_allele = row[0]
                    phenotype_field = row[2]

                    # Parse 'with' clause from phenotype
                    clean_phenotype, with_components = parse_with_clause_from_phenotype(phenotype_field)
                    if with_components:
                        # minor debug: with components present
                        pass

                else:
                    # Format 2: 7 columns
                    # Column 0: primary gene allele and secondary gene (with clause equivalent) separated by '/'
                    # Column 1:
                    # Column 2: phenotype
                    secondary_genes = row[0].split('/')[1:] if row[0] and '/' in row[0] else []
                    primary_gene_allele = row[0].split('/')[0] if '/' in row[0] else row[0]
                    clean_phenotype = row[2]

                    # Secondary gene becomes the 'with' component
                    #with_components = [secondary_genes] if secondary_genes else []
                    with_components = [s.strip() for s in secondary_genes if s and s.strip()]
                    if with_components:
                        # minor debug: with components present
                        pass

                # Extract primary gene name
                primary_gene = extract_gene_name(primary_gene_allele)
                if not primary_gene:
                    continue

                # Skip if primary gene is a helper/driver line
                if is_driver_line(primary_gene_allele):
                    continue

                # Skip multi-gene interactions for single-gene analysis
                if has_multiple_gene_disruptions(primary_gene, with_components):
                    continue

                # Classify phenotype
                phenotype_class = classify_phenotype(clean_phenotype)

                # Debug for specific genes
                if '5-HT1B' in str(row):
                    logger = logging.getLogger('PhenGO.alt_fly_phenotyping')
                    logger.debug("5-HT1B Debug:")
                    logger.debug(f"  Row: {row}")
                    logger.debug(f"  Format: {'1' if is_format1 else '2'}")
                    logger.debug(f"  Primary gene: {primary_gene}")
                    logger.debug(f"  Clean phenotype: '{clean_phenotype}'")
                    logger.debug(f"  Phenotype class: '{phenotype_class}'")
                    logger.debug(f"  With components: {with_components}")

                #if phenotype_class in ['lethal', 'viable']:
                vi_inviable_genes[primary_gene].append(phenotype_class)

            except Exception as e:
                logger = logging.getLogger('PhenGO.alt_fly_phenotyping')
                logger.error(f"Error processing row {row}: {e}")
                continue

    # Final processing
    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
                continue

        # Set value to a single string: either "viable" or "lethal"
        vi_inviable_genes[gene] = "lethal" if "lethal" in statuses else "viable"

    logger = logging.getLogger('PhenGO.alt_fly_phenotyping')
    logger.info(f"Species: fly")
    logger.info(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")

    return dict(vi_inviable_genes)

