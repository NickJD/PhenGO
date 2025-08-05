import argparse
import os

def parse_arff(file_path):
    """
    Parse ARFF file into dict: gene → label
    """
    data_started = False
    genes = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            if data_started:
                parts = [p.strip() for p in line.split(',')]
                gene = parts[0].strip()
                label = parts[-1].strip()
                genes[gene] = label
            elif line.lower() == '@data':
                data_started = True
    return genes

def main():
    parser = argparse.ArgumentParser(
        description='Compare genes and labels in two ARFF files')
    parser.add_argument('-arff_a', help='First ARFF file')
    parser.add_argument('-arff_b', help='Second ARFF file')
    args = parser.parse_args()

    genes_a = parse_arff(args.arff_a)
    genes_b = parse_arff(args.arff_b)

    set_a = set(genes_a.keys())
    set_b = set(genes_b.keys())

    shared = set_a & set_b
    only_in_a = set_a - set_b
    only_in_b = set_b - set_a

    print(f"\nGenes only in {args.arff_a}: {len(only_in_a)}")
    for g in sorted(only_in_a):
        print(f"  {g}")

    print(f"\nGenes only in {args.arff_b}: {len(only_in_b)}")
    for g in sorted(only_in_b):
        print(f"  {g}")

    same_label = []
    diff_label = []

    for g in sorted(shared):
        label_a = genes_a[g]
        label_b = genes_b[g]
        if label_a == label_b:
            same_label.append((g, label_a))
        else:
            diff_label.append((g, label_a, label_b))

    # print(f"\nShared genes with SAME label: {len(same_label)}")
    # for g, label in same_label:
    #     print(f"  {g}: {label}")

    print(f"\nShared genes with DIFFERENT label: {len(diff_label)}")
    for g, la, lb in diff_label:
        print(f"  {g}: {os.path.basename(args.arff_a)}={la} vs {os.path.basename(args.arff_b)}={lb}")

if __name__ == '__main__':
    main()
