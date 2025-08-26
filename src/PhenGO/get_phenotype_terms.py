
"""
Phenotype OBO File Parser and Term Traversal Tool

This script loads a phenotype OBO file, reads a list of phenotype terms from a text file,
finds all direct children of those terms, and outputs them with descriptions to a text file.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


class OBOParser:
    """Parser for OBO (Open Biomedical Ontologies) files."""

    def __init__(self):
        self.terms = {}  # term_id -> term_data
        self.children = defaultdict(list)  # parent_id -> [child_ids]

    def parse_obo_file(self, obo_file_path):
        """
        Parse an OBO file and extract term information.

        Args:
            obo_file_path (str): Path to the OBO file
        """
        print(f"Parsing OBO file: {obo_file_path}")

        try:
            with open(obo_file_path, 'r', encoding='utf-8') as f:
                current_term = {}
                in_term_section = False

                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines and comments
                    if not line or line.startswith('!'):
                        continue

                    # Check for section headers
                    if line == '[Term]':
                        # Save previous term if exists
                        if in_term_section and current_term.get('id'):
                            self._save_term(current_term)

                        # Start new term
                        current_term = {}
                        in_term_section = True
                        continue
                    elif line.startswith('[') and line.endswith(']'):
                        # End of term section
                        if in_term_section and current_term.get('id'):
                            self._save_term(current_term)
                        in_term_section = False
                        current_term = {}
                        continue

                    # Parse term fields
                    if in_term_section and ':' in line:
                        try:
                            key, value = line.split(':', 1)
                            key = key.strip()
                            value = value.strip()

                            if key == 'id':
                                current_term['id'] = value
                            elif key == 'name':
                                current_term['name'] = value
                            elif key == 'def':
                                # Extract definition (remove quotes and attribution)
                                if value.startswith('"') and '"' in value[1:]:
                                    def_text = value[1:value.index('"', 1)]
                                    current_term['def'] = def_text
                            elif key == 'is_a':
                                # Extract parent term ID
                                parent_id = value.split('!')[0].strip()
                                if 'is_a' not in current_term:
                                    current_term['is_a'] = []
                                current_term['is_a'].append(parent_id)
                        except ValueError:
                            continue

                # Save last term
                if in_term_section and current_term.get('id'):
                    self._save_term(current_term)

        except FileNotFoundError:
            print(f"Error: OBO file '{obo_file_path}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing OBO file: {e}")
            sys.exit(1)

        print(f"Parsed {len(self.terms)} terms from OBO file")

    def _save_term(self, term_data):
        """Save term data and build parent-child relationships."""
        term_id = term_data['id']
        self.terms[term_id] = term_data

        # Build parent-child relationships
        if 'is_a' in term_data:
            for parent_id in term_data['is_a']:
                self.children[parent_id].append(term_id)

    def get_all_descendants(self, term_id, visited=None):
        """
        Get all descendants (children, grandchildren, etc.) of a given term ID.

        Args:
            term_id (str): The term ID to find descendants for
            visited (set): Set of already visited terms to avoid cycles

        Returns:
            set: Set of all descendant term IDs
        """
        if visited is None:
            visited = set()

        # Avoid infinite loops in case of cycles
        if term_id in visited:
            return set()

        visited.add(term_id)
        descendants = set()

        # Get direct children
        direct_children = self.children.get(term_id, [])

        for child_id in direct_children:
            descendants.add(child_id)
            # Recursively get descendants of each child
            child_descendants = self.get_all_descendants(child_id, visited.copy())
            descendants.update(child_descendants)

        return descendants

    def get_term_info(self, term_id):
        """
        Get term information including name and definition.

        Args:
            term_id (str): The term ID

        Returns:
            dict: Term information or None if not found
        """
        return self.terms.get(term_id)


def load_term_list(term_list_file):
    """
    Load list of phenotype terms from a text file.

    Args:
        term_list_file (str): Path to file containing term IDs (one per line)

    Returns:
        list: List of term IDs
    """
    print(f"Loading term list from: {term_list_file}")

    try:
        with open(term_list_file, 'r', encoding='utf-8') as f:
            #terms = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            terms = [line.split('\t', 1)[0] for line in f if line.strip() and not line.startswith('#')]
        print(f"Loaded {len(terms)} terms from input file")
        return terms

    except FileNotFoundError:
        print(f"Error: Term list file '{term_list_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading term list file: {e}")
        sys.exit(1)


def write_results(output_file, results):
    """
    Write results to tab-delimited output file.

    Args:
        output_file (str): Path to output file
        results (list): List of (term_id, term_name, term_definition) tuples
    """
    print(f"Writing results to: {output_file}")

    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # Write header
            f.write("# Phenotype ID\tDescription\n")

            # Write all descendants
            for term_id, term_name, term_definition in results:
                # Use name if available, otherwise use definition, otherwise indicate missing
                description = term_name if term_name else (
                    term_definition if term_definition else "No description available")
                f.write(f"{term_id}\t{description}\n")

        print(f"Results written successfully to {output_file}")

    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)


def main():
    """Main function to orchestrate the phenotype term traversal."""

    parser = argparse.ArgumentParser(
        description="Traverse phenotype OBO file and find all descendants of specified terms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -o phenotypes.obo -t terms.txt -r results.tsv
  %(prog)s --obo-file hp.obo --term-list my_terms.txt --results descendants.tsv

Input file format:
  The term list file should contain one term ID per line, e.g.:
    HP:0000707
    HP:0001626
    HP:0000118

Output format:
  Tab-delimited file with columns: Phenotype_ID, Description
        """
    )

    parser.add_argument(
        '-o', '--obo-file',
        required=True,
        help='Path to the phenotype OBO file'
    )

    parser.add_argument(
        '-t', '--term-list',
        required=True,
        help='Path to text file containing phenotype term IDs (one per line)'
    )

    parser.add_argument(
        '-r', '--results',
        required=True,
        help='Path to output tab-delimited file for results'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    args = parser.parse_args()

    # Validate input files exist
    if not Path(args.obo_file).exists():
        print(f"Error: OBO file '{args.obo_file}' does not exist.")
        sys.exit(1)

    if not Path(args.term_list).exists():
        print(f"Error: Term list file '{args.term_list}' does not exist.")
        sys.exit(1)

    # Initialize parser and load OBO file
    obo_parser = OBOParser()
    obo_parser.parse_obo_file(args.obo_file)

    # Load term list
    input_terms = load_term_list(args.term_list)

    # Process each term and find all descendants
    all_descendants = set()
    results = []

    print("\nProcessing terms and finding all descendants...")

    for term_id in input_terms:
        if args.verbose:
            print(f"Processing term: {term_id}")

        # Check if term exists in OBO file
        if term_id not in obo_parser.terms:
            print(f"Warning: Term '{term_id}' not found in OBO file")
            continue

        # Get all descendants
        descendant_ids = obo_parser.get_all_descendants(term_id)
        all_descendants.update(descendant_ids)

        if args.verbose:
            print(f"  Found {len(descendant_ids)} total descendants")

    # Collect detailed info for all unique descendants
    print(f"\nCollecting information for {len(all_descendants)} unique descendants...")

    for term_id in sorted(all_descendants):
        term_info = obo_parser.get_term_info(term_id)
        if term_info:
            term_name = term_info.get('name', '')
            term_definition = term_info.get('def', '')
            results.append((term_id, term_name, term_definition))
        else:
            results.append((term_id, '', ''))

    print(f"Found {len(results)} total descendants across {len(input_terms)} input terms")

    # Write results to output file
    write_results(args.results, results)

    print("Processing completed successfully!")


if __name__ == "__main__":
    main()