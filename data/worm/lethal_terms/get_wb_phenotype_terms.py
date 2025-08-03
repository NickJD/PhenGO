
import sys
import re

def parse_obo_build_children(obo_path):

    id2name = {}
    children = {}

    current_id = None

    with open(obo_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current_id = None  # start of new term
            elif line.startswith("id: WBPhenotype:"):
                current_id = line.split()[1]
            elif line.startswith("name: ") and current_id:
                id2name[current_id] = line.split('name: ',1)[1].strip()
            elif line.startswith("is_a: ") and current_id:
                parent_id = line.split()[1]
                children.setdefault(parent_id, []).append(current_id)

    return id2name, children

def get_all_descendants(root_id, children_dict):
    """
    Recursively / iteratively gets all descendants of root_id.
    Returns a set including direct children, grandchildren, etc.
    """
    all_descendants = set()
    queue = [root_id]
    while queue:
        current = queue.pop()
        for child in children_dict.get(current, []):
            if child not in all_descendants:
                all_descendants.add(child)
                queue.append(child)
    return all_descendants

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_lethal_descendants.py <root_ids.txt> <wbphenotype.obo> <output.tsv>")
        sys.exit(1)

    root_ids_file, obo_file, output_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

    # Step 1: Read root IDs
    txt = open(root_ids_file).read()
    root_ids = sorted(set(re.findall(r'WBPhenotype:\d{7}', txt)))

    print(f"Loaded {len(root_ids)} root IDs.")

    # Step 2: Parse OBO
    id2name, children = parse_obo_build_children(obo_file)
    print(f"Parsed {len(id2name)} terms from ontology.")

    # Step 3: Collect all descendants
    all_terms = set(root_ids)  # include roots themselves
    for root in root_ids:
        descendants = get_all_descendants(root, children)
        all_terms.update(descendants)
        print(f"Root {root}: found {len(descendants)} descendants.")

    print(f"Total unique lethal terms (including roots): {len(all_terms)}")

    # Step 4: Write TSV
    with open(output_tsv, 'w', newline='') as out:
        out.write("ID\tName\n")
        missing = []
        for wid in sorted(all_terms):
            name = id2name.get(wid)
            if name:
                out.write(f"{wid}\t{name}\n")
            else:
                out.write(f"{wid}\t**NOT-FOUND**\n")
                missing.append(wid)

    if missing:
        print(f"Warning: {len(missing)} IDs were missing from the ontology. They were marked as **NOT-FOUND** in the TSV.")

    print(f"Done. Output written to: {output_tsv}")

if __name__ == "__main__":
    main()
