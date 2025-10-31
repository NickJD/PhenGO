import codecs
import os
import gzip
import networkx as nx
import json
import logging

logger = logging.getLogger(__name__)

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
    # Extract next term block from OBO file.
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
    # Parse tag-value pairs from a term block
    data = {}
    for line in term:
        if ': ' not in line:
            continue
        tag, value = line.split(': ', 1)
        if tag == 'is_obsolete' and value == 'true':
            # Extract ID before marking as obsolete
            for l in term:
                if l.startswith('id: '):
                    obsolete_id = l.split(': ', 1)[1]
                    obsolete_go_terms.add(obsolete_id)
                    break
            return None  # Skip obsolete terms
        if not tag in data:
            data[tag] = []
        data[tag].append(value)
    return data

def obo_to_graph(output_dir, go_obo_file):
    # Parse OBO file and create GO term graph.
    if go_obo_file.endswith('.gz'):
        obo_file = gzip.open(go_obo_file, mode='rb')
    else:
        obo_file = open(go_obo_file, mode='rb')

    logger.info("Parsing GO OBO file...")

    terms = {}
    obsolete_go_terms = set()  # Use set for efficiency

    # Skip to first term
    getTerm(obo_file)

    # Parse all terms
    term_count = 0
    obsolete_count = 0

    while True:
        term = parseTagValue(getTerm(obo_file), obsolete_go_terms)

        if term is None:
            obsolete_count += 1
            continue
        elif len(term) == 0:
            break

        term_count += 1
        termID = term['id'][0]
        alt_ids = term.get('alt_id', [])

        if 'is_a' in term:
            termParents = [p.split()[0] for p in term['is_a']]

            # Create main term entry
            if termID not in terms:
                terms[termID] = {'p': [], 'c': []}
            terms[termID]['p'] = termParents

            # Add this term as child to each parent
            for parent in termParents:
                if parent not in terms:
                    terms[parent] = {'p': [], 'c': []}
                # Check for duplicates before adding
                if termID not in terms[parent]['c']:
                    terms[parent]['c'].append(termID)

            # Handle alt_ids and prevent duplicate children
            for alt_id in alt_ids:
                if alt_id not in terms:
                    terms[alt_id] = {'p': [], 'c': []}
                terms[alt_id]['p'] = termParents

                for parent in termParents:
                    if parent not in terms:
                        terms[parent] = {'p': [], 'c': []}
                    # Only add if not already present (prevent duplicates)
                    if alt_id not in terms[parent]['c']:
                        terms[parent]['c'].append(alt_id)

    logger.info(f"Parsed {term_count} GO terms")
    logger.info(f"Found {obsolete_count} obsolete terms")
    logger.info(f"Total unique GO IDs: {len(terms)}")

    # Write terms to file
    go_child_parent_file = os.path.join(output_dir, 'GO_Children&Parents.txt')
    with codecs.open(go_child_parent_file, encoding='utf-8', mode='w') as f:
        f.write(json.dumps(terms, indent=4))

    # Create NetworkX graph directly (no need to re-read file)
    gr = nx.DiGraph()

    for node_id in terms.keys():
        gr.add_node(node_id)

    for node_id, node_data in terms.items():
        for parent in node_data['p']:
            if parent in gr and node_id in gr:
                if not gr.has_edge(parent, node_id):
                    gr.add_edge(parent, node_id)
        for child in node_data['c']:
            if node_id in gr and child in gr:
                if not gr.has_edge(node_id, child):
                    gr.add_edge(node_id, child)

    logger.info(f"Created graph with {gr.number_of_nodes()} nodes and {gr.number_of_edges()} edges")

    # Extract unique GO terms for attributes
    unique_nodes = sorted(terms.keys())

    # Write unique nodes
    refined_nodes_output = os.path.join(output_dir, 'Unique_GO_Nodes.txt')
    with open(refined_nodes_output, mode='w') as f:
        for node in unique_nodes:
            f.write(node + '\n')

    return gr, unique_nodes, obsolete_go_terms
