"""Build auditable, release-specific inputs for the publication workflow."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
from collections import defaultdict


def _open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _write_gene_set(path, records):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stable_gene_id", "gene_symbol"])
        for stable_id, symbol in sorted(records, key=lambda row: (row[0], row[1])):
            writer.writerow([stable_id, symbol])


def _normalise_viability(value):
    value = re.sub(r"\s+", " ", str(value).strip().lower())
    if not value:
        return None
    if "subviable" in value or "incomplete" in value or "conflict" in value:
        return "ambiguous"
    if re.search(r"\b(?:inviable|lethal)\b", value) and not re.search(
            r"\bnon[- ]lethal\b", value):
        return "lethal"
    if re.search(r"\b(?:viable|non[- ]lethal)\b", value):
        return "viable"
    return None


def build_mouse_impc_sets(input_path, lethal_output, viable_output, excluded_output=None):
    """Convert an annual IMPC viability report into paired mouse gene sets.

    The annual exports changed schema several times. This parser starts only at
    the detailed row-level header, accepts symbols with or without MGI IDs, and
    excludes subviable or conflicting genes instead of forcing a binary label.
    """
    observations = defaultdict(set)
    identifiers = defaultdict(lambda: {"ids": set(), "symbols": set()})
    header = None
    column_map = None

    symbol_names = {"gene", "gene symbol"}
    id_names = {"mgi gene id", "gene accession id"}
    label_names = {
        "category", "phenotype", "viability phenotype homs/hemis",
        "viability phenotype homs/hemi", "viability call",
    }

    with _open_text(input_path) as handle:
        for row in csv.reader(handle):
            cleaned = [cell.strip() for cell in row]
            lower = [cell.lower() for cell in cleaned]
            symbol_index = next(
                (index for index, value in enumerate(lower) if value in symbol_names),
                None,
            )
            label_index = next(
                (index for index, value in enumerate(lower) if value in label_names),
                None,
            )
            if symbol_index is not None and label_index is not None:
                id_index = next(
                    (index for index, value in enumerate(lower) if value in id_names),
                    None,
                )
                header = cleaned
                column_map = (symbol_index, id_index, label_index)
                continue
            if column_map is None or not any(cleaned):
                continue

            symbol_index, id_index, label_index = column_map
            if max(symbol_index, label_index) >= len(cleaned):
                continue
            symbol = cleaned[symbol_index]
            stable_id = cleaned[id_index] if id_index is not None and id_index < len(cleaned) else ""
            label = _normalise_viability(cleaned[label_index])
            if not symbol or not label:
                continue
            key = stable_id or symbol
            observations[key].add(label)
            identifiers[key]["symbols"].add(symbol)
            if stable_id:
                identifiers[key]["ids"].add(stable_id)

    if header is None:
        raise ValueError(f"No detailed IMPC viability table header found in {input_path}")

    lethal = []
    viable = []
    excluded = []
    for key, labels in observations.items():
        ids = sorted(identifiers[key]["ids"])
        symbols = sorted(identifiers[key]["symbols"])
        stable_id = ids[0] if ids else ""
        symbol = symbols[0] if symbols else key
        if labels == {"lethal"}:
            lethal.append((stable_id, symbol))
        elif labels == {"viable"}:
            viable.append((stable_id, symbol))
        else:
            excluded.append((stable_id, symbol, ";".join(sorted(labels))))

    if not lethal or not viable:
        raise ValueError(
            f"IMPC conversion did not produce both classes: lethal={len(lethal)}, "
            f"viable={len(viable)}"
        )
    _write_gene_set(lethal_output, lethal)
    _write_gene_set(viable_output, viable)
    if excluded_output:
        os.makedirs(os.path.dirname(os.path.abspath(excluded_output)), exist_ok=True)
        with open(excluded_output, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["stable_gene_id", "gene_symbol", "observed_calls"])
            writer.writerows(sorted(excluded))
    return {
        "source": os.path.abspath(input_path),
        "lethal": len(lethal),
        "viable": len(viable),
        "excluded_ambiguous": len(excluded),
    }


def build_fly_assignments(gaf_path, output_path, taxon="7227", excluded_output=None):
    """Derive unique release-matched FlyBase symbol-to-accession assignments."""
    expected = f"taxon:{taxon}"
    symbol_ids = defaultdict(set)
    missing_taxon = 0
    with _open_text(gaf_path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3 or row[0] != "FB":
                continue
            if len(row) <= 12 or not row[12]:
                missing_taxon += 1
                continue
            taxa = {value.strip() for value in row[12].split("|") if value.strip()}
            if expected not in taxa:
                continue
            symbol = row[2].strip()
            stable_id = row[1].strip()
            if symbol and stable_id:
                symbol_ids[symbol].add(stable_id)
    ambiguous = {
        symbol: stable_ids for symbol, stable_ids in symbol_ids.items()
        if len(stable_ids) > 1
    }
    assignments = {
        symbol: next(iter(stable_ids)) for symbol, stable_ids in symbol_ids.items()
        if len(stable_ids) == 1
    }
    if not assignments:
        raise ValueError(
            f"No unique FB assignments for {expected} found in {gaf_path}; "
            f"ambiguous_symbols={len(ambiguous)}; "
            f"records_without_taxon={missing_taxon}"
        )
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "#SUBMITTED ID", "SPECIES", "GENE_MODEL_STATUS", "FBID_KEY",
            "SOURCE",
        ])
        for symbol, stable_id in sorted(assignments.items()):
            writer.writerow([
                symbol, "melanogaster", "Current", stable_id, f"GAF {expected}",
            ])
    if excluded_output:
        os.makedirs(os.path.dirname(os.path.abspath(excluded_output)), exist_ok=True)
        with open(excluded_output, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_symbol", "stable_gene_ids", "reason"])
            for symbol, stable_ids in sorted(ambiguous.items()):
                writer.writerow([
                    symbol, "|".join(sorted(stable_ids)),
                    "symbol_maps_to_multiple_stable_ids",
                ])
    return {
        "source": os.path.abspath(gaf_path),
        "taxon": expected,
        "assignments": len(assignments),
        "excluded_ambiguous_symbols": len(ambiguous),
        "records_without_taxon": missing_taxon,
    }


def extract_arff_label_sets(arff_path, lethal_output, viable_output, identifiers_path=None):
    """Export paired stable-ID labels from a completed PhenGO ARFF."""
    from PhenGO.predict.data import load_arff_data

    frame, _ = load_arff_data(arff_path)
    if frame is None or frame.empty:
        raise ValueError(f"Could not load a non-empty ARFF: {arff_path}")
    symbol_by_id = {}
    if identifiers_path:
        with open(identifiers_path, encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                symbol_by_id[row.get("stable_gene_id", "")] = row.get("gene_symbol", "")

    gene_column = frame.columns[0]
    label_column = frame.columns[-1]
    lethal = []
    viable = []
    for _, row in frame.iterrows():
        gene = str(row[gene_column])
        label = str(row[label_column]).strip().lower()
        record = (gene, symbol_by_id.get(gene, ""))
        if label in {"lethal", "inviable", "essential"}:
            lethal.append(record)
        elif label in {"viable", "non-essential"}:
            viable.append(record)
        else:
            raise ValueError(f"Unsupported ARFF class {label!r} for {gene}")
    if not lethal or not viable:
        raise ValueError(
            f"ARFF does not contain both classes: lethal={len(lethal)}, viable={len(viable)}"
        )
    _write_gene_set(lethal_output, lethal)
    _write_gene_set(viable_output, viable)
    return {
        "source": os.path.abspath(arff_path),
        "lethal": len(lethal),
        "viable": len(viable),
    }


def _print_summary(summary):
    print(json.dumps(summary, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(
        prog="phengo-publication-inputs",
        description="Prepare release-specific label and assignment inputs for PhenGO.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    mouse = subparsers.add_parser("mouse-impc", help="Convert an IMPC viability CSV")
    mouse.add_argument("--input", required=True)
    mouse.add_argument("--lethal-output", required=True)
    mouse.add_argument("--viable-output", required=True)
    mouse.add_argument("--excluded-output")

    fly = subparsers.add_parser("fly-assignments", help="Build taxon-filtered FlyBase assignments")
    fly.add_argument("--gaf", required=True)
    fly.add_argument("--output", required=True)
    fly.add_argument("--excluded-output")
    fly.add_argument("--taxon", default="7227")

    arff = subparsers.add_parser("arff-labels", help="Export paired labels from an ARFF")
    arff.add_argument("--arff", required=True)
    arff.add_argument("--identifiers")
    arff.add_argument("--lethal-output", required=True)
    arff.add_argument("--viable-output", required=True)

    args = parser.parse_args()
    if args.command == "mouse-impc":
        summary = build_mouse_impc_sets(
            args.input, args.lethal_output, args.viable_output, args.excluded_output,
        )
    elif args.command == "fly-assignments":
        summary = build_fly_assignments(
            args.gaf, args.output, args.taxon, args.excluded_output,
        )
    else:
        summary = extract_arff_label_sets(
            args.arff, args.lethal_output, args.viable_output, args.identifiers,
        )
    _print_summary(summary)


if __name__ == "__main__":
    main()
