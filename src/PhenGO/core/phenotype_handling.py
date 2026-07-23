import gzip
import csv
import re
from collections import defaultdict
import logging
from .. import constants as _const

logger = logging.getLogger(__name__)

DEFAULT_FLY_HELPER_LINES_FILE = getattr(_const, 'DEFAULT_FLY_HELPER_LINES_FILE', None)


def _open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def _load_values(path, columns=(0,)):
    values = set()
    if not path:
        return values
    with _open_text(path) as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row or row[0].startswith('#'):
                continue
            for column in columns:
                if column < len(row) and row[column]:
                    values.add(row[column])
    return values


def _load_authoritative_gene_values(path):
    """Load one identifier per row, preferring the first (stable-ID) column."""
    values = set()
    if not path:
        return values
    stable_headers = {
        "gene_id", "gene id", "stable_gene_id", "stable gene id", "id",
        "accession", "gene_accession", "gene accession", "marker accession id",
    }
    symbol_headers = {
        "gene", "symbol", "gene_symbol", "gene symbol",
    }
    authority_column = 0
    with _open_text(path) as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row or row[0].startswith('#'):
                continue
            cells = [cell.strip() for cell in row]
            lower = [cell.lower() for cell in cells]
            stable_column = next((
                index for index, value in enumerate(lower[:2])
                if value in stable_headers
            ), None)
            if stable_column is not None or any(
                    value in symbol_headers for value in lower[:2]):
                if stable_column is not None:
                    authority_column = stable_column
                continue
            identifier = (
                cells[authority_column]
                if authority_column < len(cells) else ""
            )
            if not identifier:
                identifier = next((cell for cell in cells[:2] if cell), "")
            if identifier:
                values.add(identifier)
    return values


def load_gene_set_labels(lethal_path, viable_path):
    """Load a complete, explicitly labelled gene-set target definition.

    The first column is authoritative and should contain a stable identifier.
    The second column is used only when the first is blank, allowing explicitly
    documented symbol-only legacy rows without loading both identifiers as if
    they were independent genes. A gene absent from both files remains unknown.
    """
    lethal = _load_authoritative_gene_values(lethal_path)
    viable = _load_authoritative_gene_values(viable_path)
    overlap = lethal & viable
    if overlap:
        examples = ", ".join(sorted(overlap)[:10])
        raise ValueError(
            f"Gene-set labels conflict for {len(overlap)} identifier(s): {examples}"
        )
    if not lethal or not viable:
        raise ValueError(
            "Paired gene sets must each contain at least one identifier: "
            f"lethal={len(lethal)}, viable={len(viable)}"
        )
    labels = {gene: "lethal" for gene in lethal}
    labels.update({gene: "viable" for gene in viable})
    logger.info(
        "Loaded paired gene sets: %d lethal identifiers and %d viable identifiers",
        len(lethal), len(viable),
    )
    return labels


def _resolve_gene_statuses(options, observations, filtered_count):
    """Resolve repeated observations under an explicit, recorded policy."""
    has_viable = "viable" in observations
    has_lethal = "lethal" in observations
    if has_viable and has_lethal:
        policy = getattr(options, "mixed_label_policy", "exclude" if getattr(
            options, "filter_mixed_terms", False) else "lethal_wins")
        if policy == "exclude":
            filtered_count['mixed_phenotype'] += 1
            return None
        if policy == "error":
            raise ValueError("Mixed lethal and viable observations encountered")
    if has_lethal:
        return "lethal"
    if has_viable:
        return "viable"
    return None


def _record_nonlethal(options, observations, explicit_viable=False):
    if explicit_viable or getattr(options, "nonlethal_policy", "observed_viable") == "observed_viable":
        observations.append("viable")


def _field_tokens(value):
    return {token for token in re.split(r"[|,;\s]+", str(value)) if token}



def get_viable_inviable_yeast(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from SGD (yeast) phenotype data.
    Notes:
        - Column 3 contains the stable SGD accession, column 0 the systematic
          name, column 1 the gene type, column 6 the mutation type, and column 9
          the phenotype description
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

                raw_stable_id = row[3].strip() if len(row) > 3 else ""
                stable_id = raw_stable_id if re.fullmatch(r"S\d+", raw_stable_id) else ""
                gene_name = stable_id or row[0].strip()
                if getattr(options, 'strict_snapshot', False) and not stable_id:
                    filtered_count['malformed_rows'] += 1
                    continue
                if not gene_name:
                    filtered_count['malformed_rows'] += 1
                    continue
                gene_type = row[1]
                mutation_type = row[6]
                phenotype_desc = row[9].lower()
                if gene_type == "ORF" and mutation_type == "null":
                    # Non-inviable observations are only labelled viable when
                    # explicit or when the user selects the legacy policy.
                    if "inviable" in phenotype_desc:
                        vi_inviable_genes[gene_name].append("lethal")
                    else:
                        _record_nonlethal(
                            options, vi_inviable_genes[gene_name],
                            explicit_viable=("viable" in phenotype_desc),
                        )

    except Exception as e:
        logger.error(f"Error reading yeast phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        status = _resolve_gene_statuses(options, statuses, filtered_count)
        if status is not None:
            final_genes[gene] = status

    logger.info("Species: Saccharomyces cerevisiae (yeast)")
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
        - Column 1 contains the stable WormBase gene accession, column 2 the
          symbol, and column 4 the phenotype terms
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {'mixed_phenotype': 0}

    try:
        if options.worm_lethal_genes is not None:
            # Mode 1: Use pre-defined lethal gene list
            lethal_genes = _load_authoritative_gene_values(options.worm_lethal_genes)
            viable_genes = _load_authoritative_gene_values(
                getattr(options, 'worm_viable_genes', None),
            )

            logger.info(f"Loaded {len(lethal_genes)} pre-defined lethal genes")

            with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
                reader = csv.reader(input_file, delimiter='\t')
                for row in reader:
                    if len(row) > 3:
                        raw_stable_id = row[1].strip()
                        stable_id = (
                            raw_stable_id
                            if raw_stable_id.startswith("WBGene") else ""
                        )
                        symbol = row[2].strip()
                        gene_name = stable_id or symbol
                        if getattr(options, 'strict_snapshot', False) and not stable_id:
                            continue
                        # Prefer an exact stable-ID label. Symbol-only legacy
                        # lists are supported, but the output key remains the
                        # stable accession supplied by this phenotype record.
                        if stable_id in lethal_genes:
                            vi_inviable_genes[gene_name].append("lethal")
                        elif stable_id in viable_genes:
                            vi_inviable_genes[gene_name].append("viable")
                        elif (not getattr(options, 'strict_snapshot', False)
                              and symbol in lethal_genes):
                            vi_inviable_genes[gene_name].append("lethal")
                        elif (not getattr(options, 'strict_snapshot', False)
                              and symbol in viable_genes):
                            vi_inviable_genes[gene_name].append("viable")
                        else:
                            _record_nonlethal(options, vi_inviable_genes[gene_name])
        else:
            # Mode 2: Use phenotype term matching
            lethal_terms = _load_values(options.worm_phenotypes)
            viable_terms = _load_values(
                getattr(options, 'worm_viable_phenotypes', None),
            )

            logger.info(f"Loaded {len(lethal_terms)} lethal phenotype terms")

            with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
                reader = csv.reader(input_file, delimiter='\t')
                for row in reader:
                    if len(row) > 4:
                        raw_stable_id = row[1].strip()
                        stable_id = (
                            raw_stable_id
                            if raw_stable_id.startswith("WBGene") else ""
                        )
                        gene_name = stable_id or row[2].strip()
                        if getattr(options, 'strict_snapshot', False) and not stable_id:
                            continue
                        if not gene_name:
                            continue
                        qualifier = row[3] if len(row) > 3 else ""
                        phenotype_terms = row[4]
                        evidence_codes = row[6] if len(row) > 6 else ""

                        # Match against lethal terms, exclude NOT qualifiers.
                        # Evidence code filtering is only applied when the user
                        # explicitly provides codes via -worm_evidence_codes.
                        evidence_codes_filter = getattr(options, 'worm_evidence_codes', None)
                        evidence_tokens = _field_tokens(evidence_codes)
                        evidence_ok = (
                            any(ec in evidence_tokens for ec in evidence_codes_filter)
                            if evidence_codes_filter
                            else True
                        )
                        if 'NOT' in _field_tokens(qualifier) or not evidence_ok:
                            continue
                        if lethal_terms & _field_tokens(phenotype_terms):
                            vi_inviable_genes[gene_name].append("lethal")
                        else:
                            _record_nonlethal(
                                options,
                                vi_inviable_genes[gene_name],
                                explicit_viable=bool(
                                    viable_terms & _field_tokens(phenotype_terms)
                                ),
                            )

    except Exception as e:
        logger.error(f"Error reading worm phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        status = _resolve_gene_statuses(options, statuses, filtered_count)
        if status is not None:
            final_genes[gene] = status

    logger.info("Species: Caenorhabditis elegans (worm)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes


def get_viable_inviable_fly(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from FlyBase phenotype data.

    Two modes are available, controlled by options.fly_driver_filtering:

    Simple mode (default, fly_driver_filtering=False):
        - Scans all columns for "lethal"; other observations follow nonlethal_policy.
        - If filter_multigenes=True, rows with any "with" in any column or "/" in the
          gene field are discarded entirely (old conservative behaviour).
        - Explicit viability is retained; non-descript phenotypes are excluded by
          default and can be treated as viable with the legacy policy.

    Fail-closed driver-aware mode (fly_driver_filtering=True):
        - Handles two FlyBase column formats (4-col Format 1, 7-col Format 2).
        - Driver lines (GAL4, UAS, lexA, QUAS, etc.) are detected via regex patterns
          and an optional helper-lines file.  "viable, with Scer\\GAL4[x]" is stripped
          to just "viable" — the driver is not counted as a second gene.
        - If filter_multigenes=True, rows with a second gene or any unresolved
          compound component are excluded. Native gene alleles containing a UAS
          construct are not mistaken for standalone accessory components.
        - Phenotype classification checks for the words ``lethal'' and ``viable'';
          anything else follows nonlethal_policy.
    """
    vi_inviable_genes = defaultdict(list)
    filtered_count = {
        'partially_lethal': 0, 'multi_gene': 0, 'unresolved_compound': 0,
        'accessory_primary': 0, 'unsupported_schema': 0, 'mixed_phenotype': 0,
    }
    parse_audit = []

    use_smart = getattr(options, 'fly_driver_filtering', False)

    # ── Driver-line helpers (smart mode only) ─────────────────────────────────
    if use_smart:
        helper_lines_set = set()
        helper_lines_file = getattr(options, 'fly_helper_lines', None) or DEFAULT_FLY_HELPER_LINES_FILE
        if helper_lines_file:
            try:
                with gzip.open(helper_lines_file, 'rt', encoding='utf-8') as f:
                    reader_h = csv.reader(f, delimiter='\t')
                    helper_lines_set = {row[2] for row in reader_h if len(row) > 2}
            except Exception as e:
                logger.warning(f"Could not load helper lines file: {e}")

        driver_patterns = [
            r'^Scer\\GAL4(?:\[|$)', r'^Scer\\UAS(?:\.|\[|$)',
            r'^Ecol\\lexA(?:\[|$)', r'^Ecol\\lexAop(?:\.|\[|$)',
            r'^Ncra\\QF(?:\[|$)', r'^Ncra\\QUAS(?:\.|\[|$)',
            r'^Scer\\FLP\d*(?:\[|$)', r'^Scer\\GAL80(?:\[|$)',
            r'^tetO(?:_|\[|$)', r'^rtTA(?:\[|$)', r'^tTA(?:\[|$)',
        ]

        def is_driver_line(component):
            if component in helper_lines_set:
                return True
            cleaned = str(component or '').strip().strip('()')
            return any(re.search(p, cleaned, re.IGNORECASE) for p in driver_patterns)

        def extract_gene_name(allele_string):
            if not allele_string:
                return None
            cleaned = str(allele_string).strip().strip('()')
            prefix = cleaned.split('[', 1)[0]
            if '\\' in prefix:
                return None
            match = re.fullmatch(r'([^\[\]/,\s]+)\[[^\]]+\]', cleaned)
            if match:
                return match.group(1).strip()
            return cleaned if re.fullmatch(r'[A-Za-z0-9_.:+-]+', cleaned) else None

        def parse_with_clause(phenotype_field):
            text = str(phenotype_field or '')
            bracket_depth = 0
            depth_at = []
            for character in text:
                depth_at.append(bracket_depth)
                if character == '[':
                    bracket_depth += 1
                elif character == ']' and bracket_depth:
                    bracket_depth -= 1
            match = next(
                (
                    candidate for candidate in re.finditer(
                        r'\bwith\s+', text, flags=re.IGNORECASE,
                    )
                    if depth_at[candidate.start()] == 0
                ),
                None,
            )
            if match is None:
                return phenotype_field, [], False
            clean = text[:match.start()].strip().rstrip(',(').strip()
            with_part = text[match.end():].strip().rstrip(')').strip()
            components, _ = split_components(with_part, delimiters=',')
            return clean, components, True

        def split_components(value, delimiters, split_whitespace=False):
            """Split genotype syntax only at separators outside allele brackets."""
            components = []
            current = []
            bracket_depth = 0
            saw_separator = False
            for character in str(value or ''):
                if character == '[':
                    bracket_depth += 1
                elif character == ']' and bracket_depth:
                    bracket_depth -= 1
                is_separator = bracket_depth == 0 and (
                    character in delimiters
                    or (split_whitespace and character.isspace())
                )
                if is_separator:
                    saw_separator = True
                    component = ''.join(current).strip().strip('()')
                    if component:
                        components.append(component)
                    current = []
                else:
                    current.append(character)
            component = ''.join(current).strip().strip('()')
            if component:
                components.append(component)
            return components, saw_separator

        def classify_phenotype(phenotype_desc):
            if not phenotype_desc:
                return 'other'
            lower = str(phenotype_desc).lower()
            if (
                re.search(r'\b(?:partially|semi[- ]?)\s*(?:lethal|viable)\b', lower)
                or re.search(r'\b(?:non[- ]?lethal|not\s+lethal)\b', lower)
            ):
                return 'ambiguous'
            has_lethal = bool(re.search(r'\blethal\b', lower))
            has_viable = bool(re.search(r'\bviable\b', lower))
            if has_lethal and has_viable:
                return 'ambiguous'
            if has_lethal:
                return 'lethal'
            if has_viable:
                return 'viable'
            return 'other'

        def classify_compound(primary_gene, with_components):
            roles = []
            has_second_gene = False
            has_unresolved = False
            for component in with_components:
                if is_driver_line(component):
                    roles.append(f'accessory:{component}')
                    continue
                gene = extract_gene_name(component)
                if not gene:
                    has_unresolved = True
                    roles.append(f'unresolved:{component}')
                elif gene == primary_gene:
                    roles.append(f'same_gene:{component}')
                else:
                    has_second_gene = True
                    roles.append(f'second_gene:{gene}:{component}')
            if has_second_gene:
                return 'multi_gene', roles
            if has_unresolved:
                return 'unresolved', roles
            return 'resolved', roles

    # ── Parse file ─────────────────────────────────────────────────────────────
    try:
        with gzip.open(phenotype_file, 'rt', encoding='utf-8') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            for line_number, row in enumerate(reader, 1):
                if not row or row[0].startswith('#'):
                    continue

                if use_smart:
                    # ── Smart mode ─────────────────────────────────────────────
                    is_format1 = len(row) == 4
                    is_format2 = len(row) == 7
                    if not (is_format1 or is_format2):
                        filtered_count['unsupported_schema'] += 1
                        parse_audit.append({
                            'line_number': line_number, 'schema_width': len(row),
                            'primary': row[0] if row else '', 'phenotype': '',
                            'primary_gene': '', 'components': '', 'component_roles': '',
                            'outcome': 'excluded_unsupported_schema',
                        })
                        continue

                    if is_format1:
                        primary_gene_allele = row[0]
                        clean_phenotype, with_components, had_compound = parse_with_clause(row[2])
                    else:
                        genotype, had_compound = split_components(
                            row[0], delimiters='/', split_whitespace=True,
                        )
                        primary_gene_allele = genotype[0] if genotype else ''
                        secondary = genotype[1:]
                        clean_phenotype = row[2]
                        with_components = [s.strip() for s in secondary if s.strip()]

                    primary_gene = extract_gene_name(primary_gene_allele)
                    audit = {
                        'line_number': line_number, 'schema_width': len(row),
                        'primary': row[0], 'phenotype': row[2],
                        'primary_gene': primary_gene or '',
                        'components': ' | '.join(with_components),
                        'component_roles': '', 'outcome': '',
                    }
                    if is_driver_line(primary_gene_allele):
                        filtered_count['accessory_primary'] += 1
                        audit['outcome'] = 'excluded_accessory_primary'
                        parse_audit.append(audit)
                        continue
                    if not primary_gene:
                        filtered_count['unresolved_compound'] += 1
                        audit['outcome'] = 'excluded_unresolved_primary'
                        parse_audit.append(audit)
                        continue
                    if had_compound and not with_components:
                        filtered_count['unresolved_compound'] += 1
                        audit['outcome'] = 'excluded_unresolved_compound'
                        parse_audit.append(audit)
                        continue

                    phenotype_class = classify_phenotype(clean_phenotype)
                    if phenotype_class == 'ambiguous':
                        filtered_count['partially_lethal'] += 1
                        audit['outcome'] = 'excluded_ambiguous_phenotype'
                        parse_audit.append(audit)
                        if getattr(options, 'ambiguous_label_policy', 'exclude') == 'viable':
                            vi_inviable_genes[primary_gene].append('viable')
                        continue

                    compound_status, roles = classify_compound(primary_gene, with_components)
                    audit['component_roles'] = ' | '.join(roles)
                    if options.filter_multigenes:
                        if compound_status == 'multi_gene':
                            filtered_count['multi_gene'] += 1
                            audit['outcome'] = 'excluded_multi_gene'
                            parse_audit.append(audit)
                            continue
                        if compound_status == 'unresolved':
                            filtered_count['unresolved_compound'] += 1
                            audit['outcome'] = 'excluded_unresolved_compound'
                            parse_audit.append(audit)
                            continue

                    # Non-descript phenotypes follow the selected label policy.
                    audit['outcome'] = f'retained_{compound_status}_{phenotype_class}'
                    parse_audit.append(audit)
                    if phenotype_class == 'lethal':
                        vi_inviable_genes[primary_gene].append('lethal')
                    else:
                        _record_nonlethal(
                            options, vi_inviable_genes[primary_gene],
                            explicit_viable=(phenotype_class == 'viable'),
                        )

                else:
                    # ── Simple mode (default) ──────────────────────────────────
                    if len(row) <= 3:
                        continue

                    row_text = " ".join(str(value) for value in row).lower()
                    is_ambiguous = bool(
                        re.search(
                            r'\b(?:partially|semi[- ]?)\s*(?:lethal|viable)\b',
                            row_text,
                        )
                        or re.search(r'\b(?:non[- ]?lethal|not\s+lethal)\b', row_text)
                    )
                    if is_ambiguous:
                        filtered_count['partially_lethal'] += 1
                        if getattr(options, 'ambiguous_label_policy', 'exclude') == 'viable':
                            base_gene = row[0].split('[')[0]
                            vi_inviable_genes[base_gene].append('viable')
                        continue

                    base_gene = row[0].split('[')[0]

                    if options.filter_multigenes:
                        has_with_clause = any("with" in x.lower() for x in row)
                        has_compound_het = '/' in row[0]
                        if has_with_clause or has_compound_het:
                            filtered_count['multi_gene'] += 1
                            continue

                    # Nonlethal observations follow the selected label policy.
                    has_lethal = bool(re.search(r'\blethal\b', row_text))
                    has_viable = bool(re.search(r'\bviable\b', row_text))
                    if has_lethal and has_viable:
                        vi_inviable_genes[base_gene].extend(("lethal", "viable"))
                    elif has_lethal:
                        vi_inviable_genes[base_gene].append("lethal")
                    else:
                        _record_nonlethal(
                            options, vi_inviable_genes[base_gene], has_viable,
                        )

    except Exception as e:
        logger.error(f"Error reading fly phenotype file: {e}")
        raise

    # ── Resolve multiple observations ──────────────────────────────────────────
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        status = _resolve_gene_statuses(options, statuses, filtered_count)
        if status is not None:
            final_genes[gene] = status

    mode_label = "smart (driver-filtering)" if use_smart else "simple"
    logger.info(f"Species: Drosophila melanogaster (fly)  [mode: {mode_label}]")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Partially lethal: {filtered_count['partially_lethal']}")
    logger.info(f"Filtered - Multi-gene: {filtered_count['multi_gene']}")
    logger.info(f"Filtered - Unresolved compound: {filtered_count['unresolved_compound']}")
    logger.info(f"Filtered - Accessory primary: {filtered_count['accessory_primary']}")
    logger.info(f"Filtered - Unsupported schema: {filtered_count['unsupported_schema']}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")
    if use_smart:
        options._fly_parse_audit_rows = parse_audit

    return final_genes


def get_viable_inviable_fish(options, phenotype_file):
    """
    Extract gene viability/lethality classifications from ZFIN (zebrafish) phenotype data.

    Args:
        options: Configuration object with filtering options
            - filter_mixed_terms: If True, exclude genes with both viable and lethal phenotypes
        phenotype_file: Path to gzipped ZFIN phenotype TSV file

    Returns:
        dict: {stable_zfin_id: "lethal" or "viable"}

    Notes:
        - Column 10 contains the phenotype description, column 1 the symbol,
          and column 2 the stable ZFIN accession
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
                raw_stable_id = row[2].strip() if len(row) > 2 else ""
                stable_id = raw_stable_id if raw_stable_id.startswith("ZDB-") else ""
                gene_id = stable_id or row[1].strip()
                if getattr(options, 'strict_snapshot', False) and not stable_id:
                    continue
                if not gene_id:
                    continue

                # Classify based on phenotype keywords
                # Exclude semi-lethal and semi-viable (ambiguous)
                is_ambiguous = "semi-lethal" in phenotype_desc or "semi-viable" in phenotype_desc
                if is_ambiguous:
                    if getattr(options, 'ambiguous_label_policy', 'exclude') == 'viable':
                        vi_inviable_genes[gene_id].append("viable")
                    continue
                if re.search(r"\b(?:lethal|dead)\b", phenotype_desc):
                    vi_inviable_genes[gene_id].append("lethal")
                else:
                    # Explicit viability is retained. Other phenotypes only become
                    # viable under the legacy observed_viable policy.
                    _record_nonlethal(
                        options, vi_inviable_genes[gene_id],
                        explicit_viable=bool(re.search(r"\b(?:viable|alive)\b", phenotype_desc)),
                    )

    except Exception as e:
        logger.error(f"Error reading fish phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        status = _resolve_gene_statuses(options, statuses, filtered_count)
        if status is not None:
            final_genes[gene] = status

    logger.info("Species: Danio rerio (zebrafish)")
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
        lethal_terms = _load_values(options.mouse_phenotypes)
        lethal_terms.discard('ID')
        viable_terms = _load_values(
            getattr(options, 'mouse_viable_phenotypes', None),
        )
        viable_terms.discard('ID')

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
                gene_name = gene_name.strip()
                if (getattr(options, 'strict_snapshot', False)
                        and not gene_name.startswith("MGI:")):
                    continue

                if ',' in gene_name:
                    filtered_count['transgenes'] += 1
                    continue

                if phenotype_term in lethal_terms:
                    vi_inviable_genes[gene_name].append("lethal")
                else:
                    _record_nonlethal(
                        options,
                        vi_inviable_genes[gene_name],
                        explicit_viable=phenotype_term in viable_terms,
                    )



    except Exception as e:
        logger.error(f"Error reading mouse phenotype file: {e}")
        raise

    # Resolve multiple observations
    final_genes = {}
    for gene, statuses in vi_inviable_genes.items():
        status = _resolve_gene_statuses(options, statuses, filtered_count)
        if status is not None:
            final_genes[gene] = status

    logger.info("Species: Mus musculus (mouse)")
    logger.info(f"Lethal genes: {sum(1 for v in final_genes.values() if v == 'lethal')}")
    logger.info(f"Viable genes: {sum(1 for v in final_genes.values() if v == 'viable')}")
    logger.info(f"Filtered - Transgenes: {filtered_count['transgenes']}")
    logger.info(f"Filtered - Mixed phenotype: {filtered_count['mixed_phenotype']}")

    return final_genes
