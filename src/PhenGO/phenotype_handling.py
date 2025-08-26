import gzip
import csv

try: # While calling this script through pip
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *

def get_viable_inviable_yeast(options, phenotype_file):
    vi_inviable_genes = {}
    input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if "inviable" in row or "viable" in row[9]:
            vi_inviable_genes.setdefault(row[0], []).append(row[9])

    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms == True:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
        # Set value to a single string: either "viable" or "inviable"
        vi_inviable_genes[gene] = "lethal" if "inviable" in statuses else "viable"
    print(f"Species: yeast")
    print(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    print(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")
    return vi_inviable_genes

def get_viable_inviable_fly(options, phenotype_file):
    vi_inviable_genes = {}
    input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    if options.filter_multigenes == True:
        # 'with' if used for Format 1 and '/' in row[0] is used for Format 2
        for row in input:
            if len(row) > 3 and not row[0].startswith('#'):
                if any("partially lethal" in x for x in row):
                    continue
                elif any("lethal" in x and "with" in x for x in row) or '/' in row[0]:
                    continue
                elif any("lethal" in x and "with" not in x for x in row) and  '/' not in row[0]:
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("lethal")
#                    print(str(row ))
                elif not any("lethal" in x or "viable" in x or "with" in x for x in row) and '/' not in row[0]:
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("other")
                elif any("viable" in x and "with" not in x for x in row) and  '/' not in row[0]:
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("viable")
                else:
                    print(f"Unrecognized row format: {row}")
                #     vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("other")
    else:
        for row in input:
            if len(row) > 3 and not row[0].startswith('#'):
                if any("partially lethal" in x for x in row):
                    continue
                elif any("lethal" in x for x in row):
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("lethal")
                elif not any("lethal" in x for x in row) and not any("viable" in x for x in row):
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("other")
                elif any("viable" in x for x in row):
                    vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("viable")
                else:
                    print(f"Unrecognized row format: {row}")
                    #vi_inviable_genes.setdefault(row[0].split('[')[0], []).append("other")

    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms == True:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
        # Set value to a single string: either "viable" or "inviable"
        vi_inviable_genes[gene] = "lethal" if "lethal" in statuses else "viable"
    print(f"Species: fly")
    print(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    print(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")
    return vi_inviable_genes


def get_viable_inviable_fish(options, phenotype_file):
    vi_inviable_genes = {}
    input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    # l = ",lethal" in line
    # d = ",dead" in line
    # v = ",viable" in line
    # a = ",alive" in line
    for row in input:
        if ("lethal" in row[10] and "semi-lethal" not in row[10]) or "dead" in row[10]:
            vi_inviable_genes.setdefault(row[1], []).append("lethal")
        elif ("viable" in row[10] and "semi-viable" not in row[10]) or "alive" in row[10]:
            vi_inviable_genes.setdefault(row[1], []).append("viable")
        else:
            vi_inviable_genes.setdefault(row[1].split('[')[0], []).append("other")

    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms == True:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
        # Set value to a single string: either "viable" or "inviable"
        vi_inviable_genes[gene] = "lethal" if "lethal" in statuses else "viable"
    print(f"Species: fish")
    print(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    print(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")
    return vi_inviable_genes

def get_viable_inviable_worm(options, phenotype_file):
    if options.worm_lethal_genes is not None:
        # Load the worm phenotype terms from the provided file
        lethal_genes = []
        with gzip.open(options.worm_lethal_genes, 'rt', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                lethal_genes.append(row[1])

        vi_inviable_genes = {}
        input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
        input = csv.reader(input, delimiter='\t')
        for row in input:
            if len(row) > 3:
                if any(gene in row[2] for gene in lethal_genes):
                    vi_inviable_genes.setdefault(row[2], []).append("lethal")
                else:
                    vi_inviable_genes.setdefault(row[2], []).append("viable")

    else:
        # Load the worm phenotype terms from the provided file
        terms = []
        with gzip.open(options.worm_phenotypes, 'rt', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row[0].startswith('#'):
                    terms.append(row[0])

        vi_inviable_genes = {}
        input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
        input = csv.reader(input, delimiter='\t')
        for row in input:
            if len(row) > 3:
                if any(term in row[4] for term in terms) and 'NOT' not in row[3]:
                    vi_inviable_genes.setdefault(row[2], []).append("lethal")
                else:
                    vi_inviable_genes.setdefault(row[2], []).append("viable")

    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms == True:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
        # Set value to a single string: either "viable" or "inviable"
        vi_inviable_genes[gene] = "lethal" if "lethal" in statuses else "viable"

    print(f"Species: worm")
    print(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    print(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")
    return vi_inviable_genes

def get_viable_inviable_mouse(options, phenotype_file):
    # Load the worm phenotype terms from the provided file
    terms = []
    with gzip.open(options.mouse_phenotypes, 'rt', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row and row[0] != 'ID':
                terms.append(row[0])

    vi_inviable_genes = {}
    input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if len(row) > 3:
            gene = row[5]
            if ',' not in gene: # Skip rows with multiple genes - These are transgenes which are not considered
                if any(term in row[3] for term in terms):
                    vi_inviable_genes.setdefault(gene, []).append("lethal")
                else:
                    vi_inviable_genes.setdefault(gene, []).append("viable")

    for gene, statuses in list(vi_inviable_genes.items()):
        if options.filter_mixed_terms == True:
            if "viable" in statuses and "lethal" in statuses:
                del vi_inviable_genes[gene]
        # Set value to a single string: either "viable" or "inviable"
        vi_inviable_genes[gene] = "lethal" if "lethal" in statuses else "viable"

    print(f"Species: mouse")
    print(f"Lethal genes: {sum(1 for v in vi_inviable_genes.values() if v == 'lethal')}")
    print(f"Viable genes: {sum(1 for v in vi_inviable_genes.values() if v == 'viable')}")
    return vi_inviable_genes