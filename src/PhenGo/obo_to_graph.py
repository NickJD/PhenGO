import codecs
import re
import gzip
import networkx as nx
import json

#from biom.assets.exercise_api import obs_without_top


def define_graph_from_file(EdgesInput):
    # Graph creation
    gr = nx.DiGraph()
    count = 0
    with open(EdgesInput) as json_file:
        json_data = json.load(json_file)
        for key in json_data.keys():
            node = key
            gr.add_node(node)
        for key in json_data.keys():
            node = key
            for parent in json_data[key]['p']:
                pt = parent
                if gr.has_node(node) and gr.has_node(pt):
                    if not gr.has_edge(pt, node):
                        gr.add_edge(pt, node)
            for child in json_data[key]['c']:
                cd = child
                if gr.has_node(node) and gr.has_node(cd):
                    if not gr.has_edge(node, cd):
                        gr.add_edge(node, cd)
            count += 1
    return gr

def getDescendents(goid, terms):
    recursiveArray = [goid]
    if goid in terms:
        children = terms[goid]['c']
        if len(children) > 0:
            for child in children:
                recursiveArray.extend(getDescendents(child,terms))
    return set(recursiveArray)

def getAncestors(goid, terms):
    recursiveArray = [goid]
    if goid in terms:
        parents = terms[goid]['p']
        if len(parents) > 0:
            for parent in parents:
                recursiveArray.extend(getAncestors(parent,terms))
    return set(recursiveArray)

def getTerm(stream):
    block = []
    for line in stream:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        if line.strip() == "[Term]" or line.strip() == "[Typedef]":
            break
        else:
            if line.strip() != "":
                block.append(line.strip())
    return block

def parseTagValue(term, obsolete_go_terms):
    data = {}
    for line in term:
        tag = line.split(': ', 1)[0]
        value = line.split(': ', 1)[1]
        if tag == 'is_obsolete' and value == 'true':
            obsolete_go_terms.append(term[0].split(': ')[1])  # Store the ID of the obsolete term
            return None  # Skip obsolete terms
        if not tag in data:
            data[tag] = []
        data[tag].append(value)
    return data

def obo_to_graph(output_dir, go_obo_file):
    obo_file = gzip.open(go_obo_file, mode='rb')
    go_child_parent_file = output_dir + '/GO_Children&Parents.txt'
    grap_input_file = output_dir + '/GO_Children&Parents.txt'
    refined_nodes_output = output_dir + '/Unique_GO_Nodes.txt'
    terms = {}
    getTerm(obo_file)
    obsolete_go_terms = []
    # infinite loop to go through the obo file.
    # Breaks when the term returned is empty, indicating end of file
    while 1:
        term = parseTagValue(getTerm(obo_file),obsolete_go_terms)
        if term is None:
            continue  # Skip obsolete terms
        elif len(term) != 0:
            termID = term['id'][0]
            alt_ids = term.get('alt_id', [])

            if 'is_a' in term:
                termParents = [p.split()[0] for p in term['is_a']]

                # Create main term entry
                if termID not in terms:
                    terms[termID] = {'p': [], 'c': []}
                terms[termID]['p'] = termParents
                for parent in termParents:
                    if parent not in terms:
                        terms[parent] = {'p': [], 'c': []}
                    terms[parent]['c'].append(termID)

                # Handle alt_ids
                for alt_id in alt_ids:
                    if alt_id not in terms:
                        terms[alt_id] = {'p': [], 'c': []}
                    # Make alt_id point to the same parents
                    terms[alt_id]['p'] = termParents
                    for parent in termParents:
                        if parent not in terms:
                            terms[parent] = {'p': [], 'c': []}
                        terms[parent]['c'].append(alt_id)

        else:
            break
    # while 1:
    #     # get the term using the two parsing functions
    #     term = parseTagValue(getTerm(obo_file))
    #     if len(term) != 0:
    #         termID = term['id'][0]
    #         # only add to the structure if the term has an is_a tag
    #         # the is_a value contain GO ID and term definition
    #         # we only want the GO ID
    #         if 'is_a' in term:
    #             termParents = [p.split()[0] for p in term['is_a']]
    #             if termID not in terms:
    #                 # each GO ID will have two arrays of parents and children
    #                 terms[termID] = {'p': [], 'c': []}
    #             # append parents of the current term
    #             terms[termID]['p'] = termParents
    #             # for every parent term, add this current term as children
    #             for termParent in termParents:
    #                 if termParent not in terms:
    #                     terms[termParent] = {'p': [], 'c': []}
    #                 terms[termParent]['c'].append(termID)
    #     else:
    #         break

    go_child_parent_file = codecs.open(go_child_parent_file, encoding='utf-8', mode='w')
    go_child_parent_file.write(json.dumps(terms, indent=4))
    ############################################

    graph_input = open(grap_input_file, mode='r')
    unique_nodes_out = open(refined_nodes_output, mode='w')
    unique_nodes = []
    go_seen = set()
    attribute_list = []
    for line in graph_input:
        if "GO" in line:
            if line not in go_seen:
                node = line.split('""')
                matches = re.findall(r'\"(.+?)\"', node[0])
                go_t = '\n'.join(matches)
                if go_t not in go_seen:
                    go_seen.add(go_t)
                    unique_nodes_out.write(go_t + '\n')
                    unique_nodes.append(go_t)
                    attribute_list.append("@attribute " + go_t + "{0,1}")

    gr = define_graph_from_file( f"{output_dir}/GO_Children&Parents.txt")
    return gr, unique_nodes, obsolete_go_terms