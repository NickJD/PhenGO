import codecs
from itertools import repeat
import csv
import argparse
import sys
import gzip
import glob
import networkx as nx
import json


# def removed_unused_gos(vi_inviable_genes, Refined_GO_Nodes):
#     with open(Refined_GO_Nodes) as input:
#         vec = [line.strip().replace(":", "") for line in input]
#     # Filter vi_inviable_genes to only those with GO terms in vec
#     filtered_genes = {}
#     for gene, values in vi_inviable_genes.items():
#         if "go_list" in values:
#             filtered_go_list = [go for go in values["go_list"] if go.replace(":", "") in vec]
#             if filtered_go_list:
#                 # Trim binVec to only keep indexes corresponding to filtered_go_list
#                 trimmed_binvec = [
#                     values["binVec"][vec.index(go.replace(":", ""))]
#                     for go in filtered_go_list if go.replace(":", "") in vec
#                 ]
#                 filtered_genes[gene] = {
#                     "status": values["status"],
#                     "go_list": filtered_go_list,
#                     "binVec": trimmed_binvec
#                 }
#     return filtered_genes
def removed_unused_gos(vi_inviable_genes, Refined_GO_Nodes):
    # 1) read the master GO‐term list
    with open(Refined_GO_Nodes) as f:
        full_vec = [line.strip() for line in f]

    # 2) collect *all* gos actually used in your genes
    used_gos = set()
    for gene_vals in vi_inviable_genes.values():
        used_gos.update(gene_vals.get("go_list", []))

    # 3) restrict the master list to only the used terms, in order
    filtered_go_terms = [go for go in full_vec if go in used_gos]

    # 4) now rebuild each gene’s go_list AND its binVec *only* over filtered_go_terms
    filtered_genes = {}
    for gene_id, gene_vals in vi_inviable_genes.items():
        old_go_list = gene_vals.get("go_list", [])
        if not old_go_list:
            continue

        # (optional) remove duplicates but keep order
        seen = set()
        old_go_list = [go for go in old_go_list if not (go in seen or seen.add(go))]

        # new_go_list is just the intersection *in the order of* filtered_go_terms
        new_go_list = [go for go in filtered_go_terms if go in old_go_list]

        # new_binvec is a simple 0/1 for each go in filtered_go_terms
        new_binvec = [1 if go in new_go_list else 0 for go in filtered_go_terms]

        filtered_genes[gene_id] = {
            "status" : gene_vals["status"],
            "go_list": new_go_list,
            "binVec" : new_binvec
        }

    # filtered_genes = {}
    # for gene_id, gene_vals in vi_inviable_genes.items():
    #     if "go_list" not in gene_vals:
    #         # If a gene has no go_list, skip or carry it forward with empty structures:
    #         continue
    #
    #     # Remove duplicates from go_list, preserve order
    #     seen = set()
    #     old_go_list = [go for go in gene_vals["go_list"] if not (go in seen or seen.add(go))]
    #
    #     # Restrict the gene's go_list to filtered_vec, preserve order
    #     new_go_list = [go for go in filtered_vec if go in old_go_list]
    #
    #     # Build new binVec: 1 if GO term in old_go_list, else 0
    #     new_binvec = [1 if go in new_go_list else 0 for go in filtered_vec]
    #
    #     filtered_genes[gene_id] = {
    #         "status": gene_vals["status"],
    #         "go_list": new_go_list,
    #         "binVec"  : new_binvec
    #     }

    return filtered_genes, filtered_go_terms


#     with open(Refined_GO_Nodes) as input:
#         vec = [line.strip() for line in input]
#     # Find all GO terms present in any gene
#     used_gos = set()
#     for values in vi_inviable_genes.values():
#         if "go_list" in values:
#             used_gos.update(go for go in values["go_list"])
#     # Filter vec to only GO terms present in any gene
#     filtered_vec = [go for go in vec if go in used_gos]
#     # Update binVecs to only keep indexes corresponding to filtered_vec
#     filtered_genes = {}
#     for gene, values in vi_inviable_genes.items():
#         if "go_list" in values:
#             new_binvec = [values["binVec"][vec.index(go)] for go in filtered_vec]
#             filtered_genes[gene] = {
#                 "status": values["status"],
#                 "go_list": [go for go in values["go_list"] if go.replace(":", "") in filtered_vec],
#                 "binVec": new_binvec
#             }
#     return filtered_genes


def Incidents(Up, Seen, gr):
    i = 0
    for key in Up:
        if key not in Seen:
            parents = []
            try:
                parents.extend(gr.predecessors(key))
                Seen.append(key)
            except (KeyError, ValueError):
                print("Parent Missing")
            if len(parents) > 0:
                i = i + 1
                if i > 0:
                    return False
                else:
                    return True


def Duplicates(Up):
    NewUp = []
    NodesSeen = []
    for node in Up:
        if node not in NodesSeen:  # not a duplicate
            NewUp.append(node)
            NodesSeen.append(node)
    return NewUp


def define_graph_from_file(EdgesInput):
    # Graph creation
    gr = nx.DiGraph()
    count = 0
    # Input
    #GraphInput = codecs.open(GraphInput, encoding='utf-8', mode='rb')
    #EdgesInput = codecs.open(EdgesInput, encoding='utf-8', mode='rb')

    with open(EdgesInput) as json_file:
        json_data = json.load(json_file)
        for key in json_data.keys():
            node = key.replace(":", "")
            gr.add_node(node)
            #print("Added ", node)
        for key in json_data.keys():
            node = key.replace(":", "")
            for parent in json_data[key]['p']:
                #print("PARENT")
                pt = parent.replace(":", "")
                if gr.has_node(node) and gr.has_node(pt):
                    if not gr.has_edge(pt, node):
                        gr.add_edge(pt, node)
            for child in json_data[key]['c']:
                #print("CHILD")
                cd = child.replace(":", "")
                if gr.has_node(node) and gr.has_node(cd):
                    if not gr.has_edge(node, cd):
                        gr.add_edge(node, cd)
            count += 1

    return gr

def assign_go_to_vector(vi_inviable_genes, gr, Refined_GO_Nodes):
    binVec = []
    with open(Refined_GO_Nodes) as input:
        vec = [line.strip().replace(":", "") for line in input]
        binVec = [0] * len(vec)

    counter = 0
    count = 0
    Missing = []
    debug = 0
    #data = open('./Gene&GO_F_With_Lethality.txt', mode="rb")
    outputfile = open('./BinVec.txt', mode='w')
    OutMissing = open('./Missing.txt', mode='w')
    OutParents = open('./Parents.txt', mode='w')
    Testing = open('./Testings.txt', mode='w')
    Func = []
    TempFunc = []
    for gene, values in vi_inviable_genes.items():
        debug = debug + 1
        #csv = line.split(",")
        #Gene = csv[0]
        Continue = True
        goCount = 0
        Ancestors = []
        # Ancestors.append(temp)
        Seen = []
        Up = []
        for t in range(1, len(values['go_list'])):
            goCount = goCount + 1
            if "GO" in values['go_list'][goCount]:
                #print(values[1][goCount])
                temp = values['go_list'][goCount]
                temp = temp.replace(":", "")
                Func.append(gene + "\t" + temp + "\n")
                Up.append(temp)
                # Ancestors.extend(Up)
                #print("Up")
                Continue = True
                Nodes = []
                while Continue == True:
                    l = 0
                    if Incidents(Up,Seen,gr) == False:
                        for node in Up:
                            if node not in Nodes:
                                Nodes.append(node)
                                try:
                                    Up.extend(gr.predecessors(node))
                                    l = l + 1
                                    if l == 1000000:
                                        Up = Duplicates(Up)
                                        # print("Parents Added")
                                except (KeyError, ValueError):
                                    print("Error")
                                    # print("Node Size")
                                    # print(len(Nodes))
                    else:
                        Continue = False
                        #print("Root")

        Ancestors.extend(Up)
        #print("Ancestors")
        # print(Ancestors)
        # print(vec)
        ModifiedAncestors = []
        NodesSeen = []
        for node in Ancestors:
            if node not in NodesSeen:  # not a duplicate
                ModifiedAncestors.append(node)
                NodesSeen.append(node)

        del Ancestors[:]
        for Node in ModifiedAncestors:

            # OutParents.write(Node)
            try:

                #print(vec.index(Node))
                binVec[vec.index(Node)] = 1
            # Parents = gr.incidents(temp)

            except (KeyError, ValueError):
                print("Missing")
                try:
                    Missing.index(Node)
                    print("Already Missing")
                except (KeyError, ValueError):
                    Missing.append(Node)

        for x in ModifiedAncestors:
            Func.append(gene + "\t" + x + "\n")

        vi_inviable_genes[gene]["binVec"] = binVec.copy()

        Func.append('\n')

        #print(gene, binVec)

        binVec = [0] * len(vec)
        try:
            del Seen[:]
            del Up[:]
            del Ancestors[:]
            del ModifiedAncestors[:]
            del Nodes[:]
        except NameError:
            print("NameError")

    return vi_inviable_genes, vec, Func


def get_FUNC_output(vi_inviable_genes, Func, output_file):
    newFUNC = []
    geneSeen = []
    FUNCoutputfile = open(output_file, mode='w')
    tempySeen = []
    Counter = 0

    for gene, values in vi_inviable_genes.items():
        if values['status'] == 'inviable': # May need to change this to 'lethal' for other species
            for line in Func:
                if line == "\n":
                    continue

                tempFUNC = []
                if gene in line and line not in geneSeen:
                    geneSeen.append(line)
                    line = line.strip()
                    line = line.replace("GO", "GO:")
                    tempFUNC.append(str(line) + "\t1")
                    # print tempFUNC
                    newFUNC.append(tempFUNC)
                if gene in line and gene not in tempySeen:
                    tempySeen.append(gene)
                    #tempy.write(gene + ",lethal\n")
                    Counter = Counter + 1

        if values['status'] == 'viable': # May need to change this to 'xxx' for other species
            for line in Func:
                if line == "\n":
                    continue

                tempFUNC = []
                if gene in line and line not in geneSeen:
                    geneSeen.append(line)
                    line = line.strip()
                    line = line.replace("GO", "GO:")
                    tempFUNC.append(str(line) + "\t0")

                    # print tempFUNC
                    newFUNC.append(tempFUNC)
                if gene in line and gene not in tempySeen:
                    tempySeen.append(gene)
                    Counter = Counter + 1

    for element in newFUNC:
        FUNCoutputfile.write(" ".join(element) + "\n")



def get_viable_inviable_yeast(phenotype_file):
    vi_inviable_genes = {}
    input = gzip.open(phenotype_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if "inviable" in row or "viable" in row:
            vi_inviable_genes.setdefault(row[0], []).append(row[9])

    for gene, statuses in list(vi_inviable_genes.items()):
        if "viable" in statuses and "inviable" in statuses:
            del vi_inviable_genes[gene]
        else:
            # Set value to a single string: either "viable" or "inviable"
            vi_inviable_genes[gene] = "viable" if "viable" in statuses else "inviable"

    return vi_inviable_genes

def get_viability_go_data_yeast(gene_association_file, vi_inviable_genes):

    input = gzip.open(gene_association_file, 'rt', encoding='utf-8')
    input = csv.reader(input, delimiter='\t')
    for row in input:
        if row[0] == "SGD":  # FlyBase = FB
            gene= row[10].partition('|')[0]
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

    return vi_inviable_genes

def write_arff_output(vi_inviable_genes, filtered_go_terms, output_file):
    """
    Writes the vi_inviable_genes dictionary to an ARFF formatted file.
    Assumes each gene has a 'binVec' key (list of binary features) and a 'status' key ('viable' or 'inviable').
    """
    with open(output_file, 'w') as f:
        f.write("@RELATION gene_lethality\n\n")
        f.write("@ATTRIBUTE gene {" + ','.join(vi_inviable_genes.keys()) + "}\n")
        for go_term in filtered_go_terms:
            f.write(f"@ATTRIBUTE {go_term} {{0,1}}\n")
        f.write("@ATTRIBUTE class {viable,inviable}\n\n")
        f.write("@DATA\n")
        for gene, values in vi_inviable_genes.items():
            bin_vec = ",".join(map(str, values["binVec"]))
            status = values["status"]
            f.write(f"{gene},{bin_vec},{status}\n")

def main():
    parser = argparse.ArgumentParser(description="Process gene lethality data.")
    parser.add_argument('--species', required=True, help='Species name')
    parser.add_argument('--year', required=True, help='Year')
    args = parser.parse_args()

    print(f"Processing phenotype data for {args.species} in {args.year}...")


    phentoype_data = f"{args.species}/phenotype_data/{args.year}/phenotype_data*.gz"
    gene_association_data = f"{args.species}/gene_association/{args.year}/gene_association*.gz"
    phenotype_files = glob.glob(phentoype_data)
    gene_association_files = glob.glob(gene_association_data)
    if len(phenotype_files) != 1 or len(gene_association_files) != 1:
        print(f"Error: Expected one phenotype/gene_association file.")
        sys.exit(1)
    for file in phenotype_files:
        if args.species.lower() == "yeast":
            vi_inviable_genes = get_viable_inviable_yeast(file)

    for file in gene_association_files:
        if args.species.lower() == "yeast":

            vi_inviable_genes = get_viability_go_data_yeast(file, vi_inviable_genes)

    # Initialise the graph
    gr = define_graph_from_file( f"go/{args.year}/GO_Children&Parents.txt")


    vi_inviable_genes, go_terms, Func = assign_go_to_vector(vi_inviable_genes, gr, f"go/{args.year}/Refined_GO_Nodes.txt")


    vi_inviable_genes, go_terms = removed_unused_gos(vi_inviable_genes, f"go/{args.year}/Refined_GO_Nodes.txt")

    #get_FUNC_output(vi_inviable_genes, Func, f"{args.species}/{args.year}_FUNC.tab")



    write_arff_output(vi_inviable_genes, go_terms, f"{args.species}/{args.species}_{args.year}_GO.arff")




if __name__ == "__main__":
    main()
    print("Complete")






















