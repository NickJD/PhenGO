import codecs
import re
import gzip

Year = './2016'

GO_OBO_File = Year + '/go_2016-01-01.obo.gz'
GO_Child_Parent_File = Year + '/GO_Children&Parents.txt'
GraphInput = Year + '/GO_Children&Parents.txt'
GraphOutput = Year + '/GO_Nodes.txt'
NodesInput = Year + '/GO_Nodes.txt'
RefinedNodesOutput = Year + '/Refined_GO_Nodes.txt'
RefinedNodesInput = Year + '/Refined_GO_Nodes.txt'
AttributeOutput = Year + '/Ontology_Attributes.txt'

def getDescendents(goid):
    recursiveArray = [goid]
    if goid in terms:
        children = terms[goid]['c']
        if len(children) > 0:
            for child in children:
                recursiveArray.extend(getDescendents(child))

    return set(recursiveArray)


def getAncestors(goid):
    recursiveArray = [goid]
    if goid in terms:
        parents = terms[goid]['p']
        if len(parents) > 0:
            for parent in parents:
                recursiveArray.extend(getAncestors(parent))

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


def parseTagValue(term):
    data = {}
    for line in term:
        tag = line.split(': ', 1)[0]
        value = line.split(': ', 1)[1]
        if not tag in data:
            data[tag] = []

        data[tag].append(value)

    return data


oboFile = gzip.open(GO_OBO_File, mode='rb')

# declare a blank dictionary
# keys are the goids
terms = {}

# skip the file header lines
getTerm(oboFile)

# infinite loop to go through the obo file.
# Breaks when the term returned is empty, indicating end of file
while 1:
    # get the term using the two parsing functions
    term = parseTagValue(getTerm(oboFile))
    if len(term) != 0:
        termID = term['id'][0]

        # only add to the structure if the term has a is_a tag
        # the is_a value contain GOID and term definition
        # we only want the GOID
        if 'is_a' in term:
            termParents = [p.split()[0] for p in term['is_a']]

            if termID not in terms:
                # each goid will have two arrays of parents and children
                terms[termID] = {'p': [], 'c': []}

            # append parents of the current term
            terms[termID]['p'] = termParents

            # for every parent term, add this current term as children
            for termParent in termParents:
                if termParent not in terms:
                    terms[termParent] = {'p': [], 'c': []}
                terms[termParent]['c'].append(termID)
    else:
        break

import json


GO_Child_Parent_File = codecs.open(GO_Child_Parent_File, encoding='utf-8', mode='w')

GO_Child_Parent_File.write(json.dumps(terms, indent=4))
############################################



GraphInput = open(GraphInput, mode='r')
GraphOut = open(GraphOutput, mode='w')

GO_Seen = set()

GraphOut.truncate ()
for line in GraphInput:
    if "GO" in line:
        if line not in GO_Seen:
            GO = line.split('""')
            matches = re.findall(r'\"(.+?)\"',GO[0])
            join = '\n'.join(matches)

            print(join)
            GraphOut.write(join+'\n')

            GO_Seen.add(join+'\n')
GraphOut.close()


NodesInput = open(NodesInput, mode='r')

NodesOutput = open(RefinedNodesOutput, mode='w')

Nodes_Seen = set()

NodesOutput.truncate()

for line in NodesInput:
    if line not in Nodes_Seen:
        NodesOutput.write(line)

        Nodes_Seen.add(line)

NodesOutput.close()


###########

AttributeOutput = open(AttributeOutput, mode='w')

line = 0
for Line in open(RefinedNodesInput, mode='r'):
    line = line + 1
    print(line)
    AttributeOutput.write("@attribute ")
    Line = Line.replace("\n", "")
    AttributeOutput.write(Line)
    AttributeOutput.write(" {1,0}")
    AttributeOutput.write("\n")

AttributeOutput.close()