import gzip
import csv
import logging

logger = logging.getLogger(__name__)

def get_viability_go_data_yeast(gene_association_file, vi_inviable_genes):
    # Extract GO terms for yeast genes from gene association file.
    logger.info("Processing yeast GO associations...")

    processed_count = 0
    go_assigned_count = 0
    not_filtered_count = 0

    try:
        with (gzip.open(gene_association_file, 'rt', encoding='utf-8') if gene_association_file.endswith('.gz') else open(gene_association_file,'r', encoding='utf-8')) as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) < 11:  # Validate row has enough columns
                    continue

                if row[0] == "SGD":  # Saccharomyces Genome Database = SGD
                    gene = row[10].partition('|')[0]
                    go = row[4]

                    # Filter NOT qualifiers
                    if 'NOT' in row[3]:
                        not_filtered_count += 1
                        continue

                    if gene in vi_inviable_genes:
                        if isinstance(vi_inviable_genes[gene], tuple):
                            vi_inviable_genes[gene][1].append(go)
                        else:
                            vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
                        go_assigned_count += 1
                    processed_count += 1

    except Exception as e:
        logger.error(f"Error processing yeast GO associations: {e}")
        raise

    # Filter and convert
    genes_with_go = {gene: value for gene, value in vi_inviable_genes.items()
                     if isinstance(value, tuple)}

    for gene, value in genes_with_go.items():
        genes_with_go[gene] = {
            "status": str(value[0]),
            "go_list": list(map(str, value[1]))
        }

    logger.info(f"Processed {processed_count} SGD entries")
    logger.info(f"Filtered {not_filtered_count} NOT annotations")
    logger.info(f"Assigned GO terms to {go_assigned_count} gene instances")
    logger.info(f"Species: yeast")
    logger.info(f"Lethal genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'lethal')}")
    logger.info(f"Viable genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'viable')}")
    logger.info(f"Genes without GO terms: {len(vi_inviable_genes) - len(genes_with_go)}")

    return genes_with_go

def get_viability_go_data_fly(options, gene_association_file, vi_inviable_genes):
    # Extract GO terms for fly genes from gene association file.
    logger.info("Processing fly GO associations...")
    processed_count = 0
    go_assigned_count = 0
    not_filtered_count = 0

    try:
        with (gzip.open(gene_association_file, 'rt', encoding='utf-8') if gene_association_file.endswith('.gz') else open(gene_association_file,'r', encoding='utf-8')) as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            for row in reader:
                if len(row) < 5:
                    continue

                if row[0] == "FB":  # FlyBase
                    gene = row[2].partition('|')[0]
                    go = row[4]

                    # Filter NOT qualifiers
                    if 'NOT' in row[3]:
                        not_filtered_count += 1
                        continue

                    if gene in vi_inviable_genes:
                        if isinstance(vi_inviable_genes[gene], tuple):
                            vi_inviable_genes[gene][1].append(go)
                        else:
                            vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
                        go_assigned_count += 1
                    processed_count += 1

    except Exception as e:
        logger.error(f"Error processing fly GO associations: {e}")
        raise

    # Filter to genes with GO terms FIRST
    genes_with_go = {gene: value for gene, value in vi_inviable_genes.items()
                     if isinstance(value, tuple)}
    genes_without_go = [gene for gene, value in vi_inviable_genes.items()
                        if not isinstance(value, tuple)]

    # Convert ONLY the filtered dict
    for gene, value in genes_with_go.items():
        genes_with_go[gene] = {
            "status": str(value[0]),
            "go_list": list(map(str, value[1]))
        }

    # Filter to melanogaster genes
    if hasattr(options, 'fly_assignments') and options.fly_assignments:
        try:
            with (gzip.open(options.fly_assignments, 'rt', encoding='utf-8') if options.fly_assignments.endswith('.gz') else open(options.fly_assignments, 'r', encoding='utf-8')) as input_file:
                next(input_file)  # Skip header
                valid_genes = set()
                for row in csv.reader(input_file, delimiter='\t'):
                    if (len(row) > 3 and row[0] in genes_with_go and
                            row[1] == 'melanogaster' and row[2] != 'Withdrawn'):
                        valid_genes.add(row[0])

            filtered_out = len(genes_with_go) - len(valid_genes)
            genes_with_go = {gene: value for gene, value in genes_with_go.items()
                             if gene in valid_genes}
            logger.info(f"Filtered to melanogaster genes: {filtered_out} genes removed")
        except Exception as e:
            logger.warning(f"Could not apply melanogaster filter: {e}")

    logger.info(f"Processed {processed_count} FB entries")
    logger.info(f"Filtered {not_filtered_count} NOT annotations")
    logger.info(f"Assigned GO terms to {go_assigned_count} gene instances")
    logger.info(f"Species: fly")
    logger.info(f"Lethal genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'lethal')}")
    logger.info(f"Viable genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'viable')}")
    logger.info(f"Genes without GO terms: {len(genes_without_go)}")

    return genes_with_go

def get_viability_go_data_fish(gene_association_file, vi_inviable_genes):
    # Extract GO terms for fish genes
    logger.info("Processing fish GO associations...")

    processed_count = 0
    go_assigned_count = 0
    not_filtered_count = 0

    try:
        with (
        gzip.open(gene_association_file, 'rt', encoding='utf-8') if gene_association_file.endswith('.gz') else open(gene_association_file, 'r', encoding='utf-8')) as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) < 5:
                    continue

                if row[0] == "ZFIN":
                    gene = row[2]
                    go = row[4]

                    #  Filter NOT qualifiers
                    if len(row) > 3 and 'NOT' in row[3]:
                        not_filtered_count += 1
                        continue

                    if gene in vi_inviable_genes:
                        if isinstance(vi_inviable_genes[gene], tuple):
                            vi_inviable_genes[gene][1].append(go)
                        else:
                            vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
                        go_assigned_count += 1
                    processed_count += 1

    except Exception as e:
        logger.error(f"Error processing fish GO associations: {e}")
        raise

    genes_with_go = {gene: value for gene, value in vi_inviable_genes.items()
                     if isinstance(value, tuple)}

    for gene, value in genes_with_go.items():
        genes_with_go[gene] = {
            "status": str(value[0]),
            "go_list": list(map(str, value[1]))
        }

    logger.info(f"Processed {processed_count} ZFIN entries")
    logger.info(f"Filtered {not_filtered_count} NOT annotations")
    logger.info(f"Species: fish")
    logger.info(f"Lethal genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'lethal')}")
    logger.info(f"Viable genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'viable')}")
    logger.info(f"Genes without GO terms: {len(vi_inviable_genes) - len(genes_with_go)}")

    return genes_with_go

def get_viability_go_data_worm(gene_association_file, vi_inviable_genes):
    # Extract GO terms for worm genes
    logger.info("Processing worm GO associations...")

    processed_count = 0
    go_assigned_count = 0
    not_filtered_count = 0
    try:
        with (gzip.open(gene_association_file, 'rt', encoding='utf-8') if gene_association_file.endswith(
                    '.gz') else open(gene_association_file, 'r', encoding='utf-8')) as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) < 5:
                    continue

                if row[0] == "WB":  # WormBase
                    gene = row[2]
                    go = row[4]

                    #  Filter NOT qualifiers
                    if len(row) > 3 and 'NOT' in row[3]:
                        not_filtered_count += 1
                        continue

                    if gene in vi_inviable_genes:
                        if isinstance(vi_inviable_genes[gene], tuple):
                            vi_inviable_genes[gene][1].append(go)
                        else:
                            vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
                        go_assigned_count += 1
                    processed_count += 1

    except Exception as e:
        logger.error(f"Error processing mouse GO associations: {e}")
        raise

    genes_with_go = {gene: value for gene, value in vi_inviable_genes.items()
                     if isinstance(value, tuple)}

    for gene, value in genes_with_go.items():
        genes_with_go[gene] = {
            "status": str(value[0]),
            "go_list": list(map(str, value[1]))
        }

    logger.info(f"Processed {processed_count} WB entries")
    logger.info(f"Filtered {not_filtered_count} NOT annotations")
    logger.info(f"Species: worm")
    logger.info(f"Lethal genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'lethal')}")
    logger.info(f"Viable genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'viable')}")
    logger.info(f"Genes without GO terms: {len(vi_inviable_genes) - len(genes_with_go)}")

    return genes_with_go

def get_viability_go_data_mouse(gene_association_file, vi_inviable_genes):
    # Extract GO terms for mouse genes
    logger.info("Processing mouse GO associations...")

    processed_count = 0
    go_assigned_count = 0
    not_filtered_count = 0
    try:
        with (
        gzip.open(gene_association_file, 'rt', encoding='utf-8') if gene_association_file.endswith('.gz') else open(gene_association_file, 'r', encoding='utf-8')) as input_file:
            reader = csv.reader(input_file, delimiter='\t')

            for row in reader:
                if len(row) < 5:
                    continue

                if row[0] == "MGI":
                    # FIX: was row[1] which is the MGI accession ID (e.g. "MGI:97490").
                    # The phenotype handler (get_viable_inviable_mouse) uses row[5] which is
                    # the gene *symbol* (e.g. "Trp53"). GAF column 2 (row[2]) is the DB Object
                    # Symbol — the gene name — and this is what matches the phenotype dict keys.
                    gene = row[2]
                    go = row[4]

                    #  Filter NOT qualifiers
                    if len(row) > 3 and 'NOT' in row[3]:
                        not_filtered_count += 1
                        continue

                    if gene in vi_inviable_genes:
                        if isinstance(vi_inviable_genes[gene], tuple):
                            vi_inviable_genes[gene][1].append(go)
                        else:
                            vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
                        go_assigned_count += 1
                    processed_count += 1

    except Exception as e:
        logger.error(f"Error processing mouse GO associations: {e}")
        raise

    genes_with_go = {gene: value for gene, value in vi_inviable_genes.items()
                     if isinstance(value, tuple)}

    for gene, value in genes_with_go.items():
        genes_with_go[gene] = {
            "status": str(value[0]),
            "go_list": list(map(str, value[1]))
        }

    logger.info(f"Processed {processed_count} MGI entries")
    logger.info(f"Filtered {not_filtered_count} NOT annotations")
    logger.info(f"Species: mouse")
    logger.info(f"Lethal genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'lethal')}")
    logger.info(f"Viable genes with GO terms: {sum(1 for v in genes_with_go.values() if v['status'] == 'viable')}")
    logger.info(f"Genes without GO terms: {len(vi_inviable_genes) - len(genes_with_go)}")

    return genes_with_go