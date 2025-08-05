import gzip
import csv

def get_viability_go_data_yeast(gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "SGD":  # FlyBase = FB
            gene= row[10].partition('|')[0]
            go = row[4]
            if gene in vi_inviable_genes:
                if isinstance(vi_inviable_genes[gene], tuple):
                    vi_inviable_genes[gene][1].append(go)
                else:
                    # Convert string value to tuple: (original string, [go])
                    vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
    # Filter genes to only those with a GO list (tuple value)
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if isinstance(value, tuple)}
    # Convert tuple values to (string, list of strings)
    for gene, value in vi_inviable_genes.items():
        vi_inviable_genes[gene] = {"status": str(value[0]), "go_list": list(map(str, value[1]))}

    print(f"Species: yeast")
    print(f"Lethal genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'lethal')}")
    print(f"Viable genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'viable')}")
    return vi_inviable_genes

def get_viability_go_data_fly(options, gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "FB":  # FlyBase = FB
            gene= row[2].partition('|')[0]
            go = row[4]
            #dataMarker = row[6] # Current we are not recording the data marker, but it can be used if needed
            if gene in vi_inviable_genes:
                if isinstance(vi_inviable_genes[gene], tuple):
                    vi_inviable_genes[gene][1].append(go)
                else:
                    # Convert string value to tuple: (original string, [go])
                    vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
    # Filter genes to only those with a GO list (tuple value)
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if isinstance(value, tuple)}
    # Convert tuple values to (string, list of strings)
    for gene, value in vi_inviable_genes.items():
        vi_inviable_genes[gene] = {"status": str(value[0]), "go_list": list(map(str, value[1]))}

    # Filter vi_inviable_genes to only those present in the file with 'melanogaster' in the same row
    #Needs to change
    with gzip.open(options.fly_assignments, mode='rt', encoding='utf-8') as f:
        next(f)  # Skip header
        valid_genes = set()
        for row in csv.reader(f, delimiter='\t'):
            if len(row) > 4 and row[0] in vi_inviable_genes and row[3] != 'Withdrawn' and row[4] == 'melanogaster':
                valid_genes.add(row[0])
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if gene in valid_genes}

    print(f"Species: fly")
    print(f"Lethal genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'lethal')}")
    print(f"Viable genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'viable')}")
    return vi_inviable_genes

def get_viability_go_data_fish(gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "ZFIN":
            gene= row[2]
            go = row[4]
            if gene in vi_inviable_genes:
                if isinstance(vi_inviable_genes[gene], tuple):
                    vi_inviable_genes[gene][1].append(go)
                else:
                    # Convert string value to tuple: (original string, [go])
                    vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
    # Filter genes to only those with a GO list (tuple value)
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if isinstance(value, tuple)}
    # Convert tuple values to (string, list of strings)
    for gene, value in vi_inviable_genes.items():
        vi_inviable_genes[gene] = {"status": str(value[0]), "go_list": list(map(str, value[1]))}

    print(f"Species: fish")
    print(f"Lethal genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'lethal')}")
    print(f"Viable genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'viable')}")
    return vi_inviable_genes

def get_viability_go_data_worm(gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "WB":  # WormBase = FB
            gene = row[2]
            go = row[4]
            if gene in vi_inviable_genes:
                if isinstance(vi_inviable_genes[gene], tuple):
                    vi_inviable_genes[gene][1].append(go)
                else:
                    # Convert string value to tuple: (original string, [go])
                    vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
    # Filter genes to only those with a GO list (tuple value)
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if isinstance(value, tuple)}
    # Convert tuple values to (string, list of strings)
    for gene, value in vi_inviable_genes.items():
        vi_inviable_genes[gene] = {"status": str(value[0]), "go_list": list(map(str, value[1]))}

    print(f"Species: worm")
    print(f"Lethal genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'lethal')}")
    print(f"Viable genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'viable')}")
    return vi_inviable_genes

def get_viability_go_data_mouse(gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "MGI":  # FlyBase = FB
            gene = row[1]
            go = row[4]
            if gene in vi_inviable_genes:
                if isinstance(vi_inviable_genes[gene], tuple):
                    vi_inviable_genes[gene][1].append(go)
                else:
                    # Convert string value to tuple: (original string, [go])
                    vi_inviable_genes[gene] = (vi_inviable_genes[gene], [go])
    # Filter genes to only those with a GO list (tuple value)
    vi_inviable_genes = {gene: value for gene, value in vi_inviable_genes.items() if isinstance(value, tuple)}
    # Convert tuple values to (string, list of strings)
    for gene, value in vi_inviable_genes.items():
        vi_inviable_genes[gene] = {"status": str(value[0]), "go_list": list(map(str, value[1]))}

    print(f"Species: mouse")
    print(f"Lethal genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'lethal')}")
    print(f"Viable genes with GO terms: {sum(1 for v in vi_inviable_genes.values() if v['status'] == 'viable')}")

    return vi_inviable_genes