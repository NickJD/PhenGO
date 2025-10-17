import argparse
import sys
import shutil
import os
import networkx as nx

try: # Try to import from the package if available
    from .obo_to_graph import obo_to_graph
    from .constants import *
    from .phenotype_handling import *
    from .go_handling import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from obo_to_graph import obo_to_graph
    from constants import *
    from phenotype_handling import *
    #from alt_fly_phenotyping import *
    from go_handling import *



def removed_unused_gos(vi_inviable_genes, unique_go_terms):
    # This optional function removes unused GO terms from the ARFF and FUNC outputs
    # It rebuilds the GO term list and each gene’s bin_vec to only include used terms
    # collect *all* GOs assigned to at least 1 gene
    used_gos = set()
    for gene_vals in vi_inviable_genes.values():
        used_gos.update(gene_vals.get("go_list", []))
    # restrict the master list to only the used terms, in order
    filtered_go_terms = [go for go in unique_go_terms if go in used_gos]
    # now rebuild each gene’s go_list AND its bin_vec *only* over filtered_go_terms
    filtered_genes = {}
    for gene_id, gene_vals in vi_inviable_genes.items():
        old_go_list = gene_vals.get("go_list", [])
        if not old_go_list:
            continue
        seen = set()
        old_go_list = [go for go in old_go_list if not (go in seen or seen.add(go))]
        new_go_list = [go for go in filtered_go_terms if go in old_go_list]
        new_bin_vec = [1 if go in new_go_list else 0 for go in filtered_go_terms]
        filtered_genes[gene_id] = {
            "status" : gene_vals["status"],
            "go_list": new_go_list,
            "bin_vec" : new_bin_vec
        }
    return filtered_genes, filtered_go_terms

def Incidents(options, Up, Seen, gr, obsolete_go_terms):
    # Check if any node in Up has no parents
    i = 0
    for key in Up:
        if key not in Seen:
            parents = []
            try:
                try:
                    Up.extend(gr.predecessors(key))
                except (KeyError, ValueError, nx.exception.NetworkXError):
                    if key in obsolete_go_terms:
                        print(f"Warning: GO term '{key}' is obsolete but present in the gene association file. Skipping.")
                        Seen.append(key)
                        continue
                    print(f"Error: The gene association file '{options.gene_association_file}' contains GO term '{key}' not present in the GO OBO file '{options.go_obo_file}'.")
                    print("Please ensure both files are consistent. Exiting.")
                    sys.exit(1)
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
    # Remove duplicates from Up
    NewUp = []
    NodesSeen = []
    for node in Up:
        if node not in NodesSeen:  # not a duplicate
            NewUp.append(node)
            NodesSeen.append(node)
    return NewUp

def assign_go_to_vector(options, vi_inviable_genes, gr, unique_go_terms, obsolete_go_terms):
    # Assign GO terms to binary vectors for each gene
    bin_vec = [0] * len(unique_go_terms)
    missing = []
    debug = 0
    o_rep = []
    for gene, values in vi_inviable_genes.items():
        debug = debug + 1
        go_count = 0
        Ancestors = []
        Seen = []
        Up = []
        for t in range(1, len(values['go_list'])):
            go_count = go_count + 1
            if "GO" in values['go_list'][go_count]:
                temp = values['go_list'][go_count]
                temp = temp
                o_rep.append(gene + "\t" + temp + "\n")
                Up.append(temp)
                Continue = True
                Nodes = []
                while Continue == True:

                    if Incidents(options, Up, Seen, gr, obsolete_go_terms) == False:
                        for node in Up:
                            if node not in Nodes:
                                Nodes.append(node)
                                try:
                                    Up.extend(gr.predecessors(node))
                                except (KeyError, ValueError):
                                    print(f"Error: Node {node} not found in graph or invalid value")
                        Up = Duplicates(Up)
                    else:
                        Continue = False
        ###########################
        Ancestors.extend(Up)
        ModifiedAncestors = []
        NodesSeen = []
        for node in Ancestors:
            if node not in NodesSeen:  # not a duplicate
                ModifiedAncestors.append(node)
                NodesSeen.append(node)

        del Ancestors[:]
        for Node in ModifiedAncestors:
            try:
                bin_vec[unique_go_terms.index(Node)] = 1
            except (KeyError, ValueError):
                print(f"Missing node: {Node} ")
                try:
                    missing.index(Node)
                    print(f"Node {Node} already marked as missing")
                except (KeyError, ValueError):
                    missing.append(Node)
        for x in ModifiedAncestors:
            o_rep.append(gene + "\t" + x + "\n")
        vi_inviable_genes[gene]["bin_vec"] = bin_vec.copy()
        o_rep.append('\n')
        bin_vec = [0] * len(unique_go_terms)
        try:
            del Seen[:]
            del Up[:]
            del Ancestors[:]
            del ModifiedAncestors[:]
            del Nodes[:]
        except NameError:
            continue

    return vi_inviable_genes, unique_go_terms, o_rep


def get_o_rep_output(vi_inviable_genes, o_rep, output_file):
    # Generate FUNC/Overrepresentation output file
    new_o_rep = []
    geneSeen = []
    FUNCoutputfile = open(output_file, mode='w')
    tempySeen = []
    Counter = 0

    for gene, values in vi_inviable_genes.items():
        if values['status'] in ['inviable','lethal']:  # Accept both 'inviable' and 'lethal' for other species for line in Func:
                if line == "\n":
                    continue
                tmp_o_rep = []
                if gene in line and line not in geneSeen:
                    geneSeen.append(line)
                    line = line.strip()
                    tmp_o_rep.append(str(line) + "\t1")
                    new_o_rep.append(tmp_o_rep)
                if gene in line and gene not in tempySeen:
                    tempySeen.append(gene)
                    Counter = Counter + 1

        if values['status'] == 'viable':
            for line in o_rep:
                if line == "\n":
                    continue
                tempFUNC = []
                if gene in line and line not in geneSeen:
                    geneSeen.append(line)
                    line = line.strip()
                    tempFUNC.append(str(line) + "\t0")
                    new_o_rep.append(tmp_o_rep)
                if gene in line and gene not in tempySeen:
                    tempySeen.append(gene)
                    Counter = Counter + 1

    for element in new_o_rep:
        FUNCoutputfile.write(" ".join(element) + "\n")


def write_arff_output(vi_inviable_genes, filtered_go_terms, output_file):
    with open(output_file, 'w') as f:
        f.write("@RELATION gene_lethality\n\n")
        f.write(f"@ATTRIBUTE gene string\n")
        for go_term in filtered_go_terms:
            f.write(f"@ATTRIBUTE {go_term} {{0,1}}\n")
        f.write("@ATTRIBUTE class {viable,lethal}\n\n")
        f.write("@DATA\n")
        for gene, values in vi_inviable_genes.items():
            gene = gene.replace("'", "-")  # Replace single quotes with underscores
            bin_vec = ",".join(map(str, values["bin_vec"]))
            status = values["status"]
            f.write(f"{gene},{bin_vec},{status}\n")

def main():
    parser = argparse.ArgumentParser(description=f"PhenGO {PhenGO_VERSION} - Convert phenotype and GO data to ARFF format")
    parser.add_argument("--print-defaults", action="store_true",
                        help="Print default files and methods used for each species and exit.")

    # Parse known args first to check for --print-defaults
    args, remaining_argv = parser.parse_known_args()
    if args.print_defaults:
        print("\nDefault options and methods for each species:\n")
        print("Fly:")
        print(f"  Assignments file: {DEFAULT_FLY_SPECIES_FIELDS_FILE}")
        print(f"  Lethal genes file: {DEFAULT_FLY_LETHAL_GENES_FILE}")
        #print(f"  Driver lines file: {DEFAULT_FLY_HELPER_LINES_FILE}")
        print("  Method: get_viable_inviable_fly, get_viability_go_data_fly")
        print("Worm:")
        print(f"  Lethal phenotypes file: {DEFAULT_WORM_LETHAL_PHENOTYPES_FILE}")
        print(f"  Lethal genes file: {DEFAULT_WORM_LETHAL_GENES_FILE}")
        print("  Method: By default we do not use the lethal genes file, we use the default lethal phenotypes file"
              " to identify the genes in the user-provided phenotype_data file to identify lethal genes. If lethal "
              "genes file is provided (or 'default' is given, it will be used to identify lethal genes instead.")
        print("Mouse:")
        print(f"  Phenotypes file: {DEFAULT_MOUSE_PHENOTYPES_FILE}")
        print("  Method: get_viable_inviable_mouse, get_viability_go_data_mouse")
        print("Yeast:")
        print("  Method: get_viable_inviable_yeast, get_viability_go_data_yeast")
        print("Fish:")
        print("  Method: get_viable_inviable_fish, get_viability_go_data_fish")
        sys.exit(0)


    parser._action_groups.pop()
    required = parser.add_argument_group('Required Options')
    required.add_argument('-species', dest="species", required=True, help='Species tag (e.g., fly, yeast)')
    required.add_argument('-phenotype_file', dest="phenotype_file", required=True, help='Path to the phenotype data file (.gz)')
    required.add_argument('-gene_association_file', dest="gene_association_file", required=True, help='Path to the gene association file (.gz)')
    required.add_argument('-go_obo_file', dest="go_obo_file", required=True, help='Path to the go.obo file')
    required.add_argument('-output_dir', dest="output_dir", required=True, help='Output directory')

    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-filter_unused_gos', dest='filter_unused_gos', action='store_false', required=False,
                        help='Filter out unused GO terms from the FUNC and ARFF output (default: True)')
    optional.add_argument('-filter_mixed_terms', dest='filter_mixed_terms', action='store_true', required=False,
                        help='Filter out genes which have both lethal and viable phenotypes - '
                             'Terms not specifically lethal/viable are not counted in this (default: False)')
    optional.add_argument('-gene_go_pheno', dest='gene_go_pheno', action='store_true', required=False,
                        help='Output "Gene-GO-Phenotype" (Rbbp5	GO:0003674	0) file for overrepresentation analysis with tools such as FUNC (default: False)')

    fly_args = parser.add_argument_group('Fly specific parameters')
    fly_args.add_argument('-fly_lethal_genes', dest='fly_lethal_genes', required=False,
                        help='Provide TSV file of specified lethal fly genes (provide "default" for default lethal genes: "data/fly/FlyBase_Lethal_Gene_IDs_2025-08-15.txt.gz")')
    fly_args.add_argument('-fly_assignments', dest='fly_assignments', required=False,
                        help='Provide TSV file of fly assignments (file confirming genes are assignment to drosophila melanogaster (default: "data/fly/FlyBase_Fields_2017.txt.gz")')
    # fly_args.add_argument('-fly_helper_lines', dest='fly_helper_lines', required=False,
    #                     help='Provide TSV file of fly driver lines (file containing the name of driver lines (RNAi) to ignore when present with the "with" tag (default: "data/fly/FlyBase_DriverLine_Fields_2025_08_05.txt.gz")')
    fly_args.add_argument('-filter_multigenes', dest='filter_multigenes', action='store_true', required=False,
                        help='Filter out phenotype with "with" tag (default: DO NOT FILTER)')

    worm_args = parser.add_argument_group('Worm specific parameters')
    worm_args.add_argument('-worm_phenotypes', dest='worm_phenotypes', required=False,
                        help='Provide TSV file of worm phenotypes (default: "data/worm/lethal_terms_traversed_2025-08-12.tsv.gz")')
    worm_args.add_argument('-worm_lethal_genes', dest='worm_lethal_genes', required=False,
                        help='Provide TSV file of specified lethal worm genes (provide "default" for default lethal genes: "data/worm/genes_direct_and_inferred_for_WBPhenotype_0000062_11-08-2025.txt.gz")')

    mouse_args = parser.add_argument_group('Mouse specific parameters')
    mouse_args.add_argument('-mouse_phenotypes', dest='mouse_phenotypes', required=False,
                        help='Provide TSV file of mouse phenotypes (default: "data/mouse/mouse_lethal_terms.txt.gz")')

    misc = parser.add_argument_group('Misc')
    misc.add_argument("-v", "--version", action="version",
                 version=f"PyamilySeq {PhenGO_VERSION}: Exiting.")


    options = parser.parse_args()

    print(f"Processing phenotype data for species: {options.species}")
    print(f"Phenotype file: {options.phenotype_file}")
    print(f"Gene association file: {options.gene_association_file}")
    print(f"GO OBO file: {options.go_obo_file}")

### FLY
    if options.species.lower() == "fly":
        ## By default we do not use the fly lethal genes file
        if options.fly_lethal_genes == None:
            pass
        elif options.fly_lethal_genes.lower() == "default":
            print(f"Using default fly lethal genes file: : {DEFAULT_FLY_LETHAL_GENES_FILE}")
            options.fly_lethal_genes = DEFAULT_FLY_LETHAL_GENES_FILE
        else:
            if not os.path.exists(options.fly_lethal_genes):
                print(f"Error: Fly lethal genes file {options.fly_lethal_genes} does not exist.")
                sys.exit(1)
            else:
                print(f"Fly lethal genes file: {options.fly_lethal_genes}")
        ## If the user does not provide a FlyBase helper lines file, we use the default
        # if options.fly_helper_lines is None:
        #     options.fly_helper_lines = DEFAULT_FLY_HELPER_LINES_FILE
        # else:
        #     if not os.path.exists(options.fly_helper_lines):
        #         print(f"Error: Fly helper lines file {options.fly_helper_lines} does not exist.")
        #         sys.exit(1)
        # print(f"Fly helper lines file: {options.fly_helper_lines}")
        ## If the user does not provide a fly species assignments file, we use the default
        if options.fly_assignments is None:
            options.fly_assignments = DEFAULT_FLY_SPECIES_FIELDS_FILE
        else:
            if not os.path.exists(options.fly_assignments):
                print(f"Error: Fly assignments file {options.fly_assignments} does not exist.")
                sys.exit(1)
        print(f"Fly assignments file: {options.fly_assignments}")


        # if options.driver_lines is None:
        #     options.driver_lines = DEFAULT_FLY_DRIVER_LINES_FILE
        # else:
        #     if not os.path.exists(options.driver_lines):
        #         print(f"Error: Fly driver lines file {options.driver_lines} does not exist.")
        #         sys.exit(1)
        # print(f"Fly driver lines file: {options.driver_lines}")

### WORM
    elif options.species.lower() == "worm":
        ## By default we do not use the worm lethal genes file
        if options.worm_lethal_genes == None:
            pass
        elif options.worm_lethal_genes.lower() == "default":
            print(f"Using default worm lethal genes file: : {DEFAULT_WORM_LETHAL_GENES_FILE}")
            options.worm_lethal_genes = DEFAULT_WORM_LETHAL_GENES_FILE
        else:
            if not os.path.exists(options.worm_lethal_genes):
                print(f"Error: Worm lethal genes file {options.worm_lethal_genes} does not exist.")
                sys.exit(1)
            else:
                print(f"Worm lethal genes file: {options.worm_lethal_genes}")
        ## If the user does not provide a worm phenotypes file, we use the default
        if options.worm_phenotypes is None:
            options.worm_phenotypes = DEFAULT_WORM_LETHAL_PHENOTYPES_FILE
            print(f"Using default Worm phenotypes file: {options.worm_phenotypes}")
        else:
            if not os.path.exists(options.worm_phenotypes):
                print(f"Error: Worm phenotypes file {options.worm_phenotypes} does not exist.")
                sys.exit(1)
            else:
                print(f"Worm phenotypes file: {options.worm_phenotypes}")

### MOUSE
    elif options.species.lower() == "mouse":
        if options.mouse_phenotypes is None:
            options.mouse_phenotypes = DEFAULT_MOUSE_PHENOTYPES_FILE
        else:
            if not os.path.exists(options.mouse_phenotypes):
                print(f"Error: Mouse phenotypes file {options.mouse_phenotypes} does not exist.")
                sys.exit(1)
        print(f"Mouse phenotypes file: {options.mouse_phenotypes}")

###
    print(f"Output directory: {options.output_dir}")

    # Ensure output directory exists and is empty
    if os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir)

    phenotype_files = [options.phenotype_file]
    gene_association_files = [options.gene_association_file]
    if len(phenotype_files) != 1 or len(gene_association_files) != 1:
        print(f"Error: Expected one phenotype/gene_association file.")
        sys.exit(1)
    for file in phenotype_files:
        if __name__ == '__main__':
            if options.species.lower() == "yeast":
                vi_inviable_genes = get_viable_inviable_yeast(options, file)
            elif options.species.lower() == "fly":
                vi_inviable_genes = get_viable_inviable_fly(options, file)
            elif options.species.lower() == "fish":
                vi_inviable_genes = get_viable_inviable_fish(options, file)
            elif options.species.lower() == "worm":
                vi_inviable_genes = get_viable_inviable_worm(options, file)
            elif options.species.lower() == "mouse":
                vi_inviable_genes = get_viable_inviable_mouse(options, file)


    for file in gene_association_files:
        if options.species.lower() == "yeast":
            vi_inviable_genes = get_viability_go_data_yeast(file, vi_inviable_genes)
        elif options.species.lower() == "fly":
            vi_inviable_genes = get_viability_go_data_fly(options, file, vi_inviable_genes)
        elif options.species.lower() == "fish":
            vi_inviable_genes = get_viability_go_data_fish(file, vi_inviable_genes)
        elif options.species.lower() == "worm":
            vi_inviable_genes = get_viability_go_data_worm(file, vi_inviable_genes)
        elif options.species.lower() == "mouse":
            vi_inviable_genes = get_viability_go_data_mouse(file, vi_inviable_genes)

    gr, unique_go_terms, obsolete_go_terms = obo_to_graph(options.output_dir, options.go_obo_file)


    vi_inviable_genes, go_terms, Func = assign_go_to_vector(options, vi_inviable_genes, gr, unique_go_terms, obsolete_go_terms)

    # Filter out unused GO terms
    if options.filter_unused_gos == True:
        vi_inviable_genes, go_terms = removed_unused_gos(vi_inviable_genes, unique_go_terms)

    # Write the filtered GO terms to a file - SLOW
    if options.gene_go_pheno == True:
        get_o_rep_output(vi_inviable_genes, Func, f"{options.output_dir}/{options.species}_OverRepresentation.tab")

    write_arff_output(vi_inviable_genes, go_terms, f"{options.output_dir}/{options.species}_PhenGO.arff")




if __name__ == "__main__":
    main()
    print("Complete")






















