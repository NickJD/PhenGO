"""Aggregate a completed PhenGO publication workflow into paper-facing tables."""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from collections import Counter
from datetime import datetime, timezone


def _write_tsv(path, rows):
    if not rows:
        return False
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    return True


def _read_csv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def collect_snapshot_rows(root):
    rows = []
    pattern = os.path.join(root, "02_arff", "*", "*", "*", "PhenGO_manifest.json")
    for path in sorted(glob.glob(pattern)):
        relative = os.path.relpath(path, os.path.join(root, "02_arff"))
        track, organism, timepoint, _ = relative.split(os.sep, 3)
        with open(path, encoding="utf-8") as handle:
            manifest = json.load(handle)
        counts = manifest.get("counts", {})
        releases = manifest.get("releases", {})
        policies = manifest.get("policies", {})
        inputs = manifest.get("inputs", {})
        rows.append({
            "track": track,
            "organism": organism,
            "timepoint": timepoint,
            "snapshot_id": manifest.get("snapshot_id", ""),
            "strict_snapshot": manifest.get("strict_snapshot", ""),
            "genes": counts.get("genes", ""),
            "lethal": counts.get("lethal", ""),
            "viable": counts.get("viable", ""),
            "go_features": counts.get("go_features", ""),
            "phenotype_release": releases.get("phenotype", ""),
            "go_annotation_release": releases.get("go_annotations", ""),
            "go_ontology_release": releases.get("go_ontology", ""),
            "phenotype_ontology_release": releases.get("phenotype_ontology", ""),
            "label_source": policies.get("label_source", ""),
            "nonlethal_policy": policies.get("nonlethal", ""),
            "mixed_label_policy": policies.get("mixed_label", ""),
            "go_relations": ";".join(policies.get("go_relations", [])),
            "go_evidence_exclude": ";".join(policies.get("go_evidence_exclude", [])),
            "arff_sha256": manifest.get("outputs", {}).get("arff", {}).get("sha256", ""),
            "phenotype_sha256": inputs.get("phenotype", {}).get("sha256", ""),
            "gaf_sha256": inputs.get("go_annotations", {}).get("sha256", ""),
            "go_obo_sha256": inputs.get("go_ontology", {}).get("sha256", ""),
            "manifest_path": os.path.abspath(path),
        })
    return rows


def collect_analysis_rows(root, filename):
    rows = []
    pattern = os.path.join(root, "04_ml", "*", "*", filename)
    for path in sorted(glob.glob(pattern)):
        relative = os.path.relpath(path, os.path.join(root, "04_ml"))
        track, organism, _ = relative.split(os.sep, 2)
        for row in _read_csv(path):
            rows.append({"track": track, "organism": organism, **row})
    return rows


def collect_temporal_rows(root, suffix="*_summary_timeline.csv"):
    rows = []
    pattern = os.path.join(root, "05_temporal", "*", "*", suffix)
    for path in sorted(glob.glob(pattern)):
        relative = os.path.relpath(path, os.path.join(root, "05_temporal"))
        track, organism, _ = relative.split(os.sep, 2)
        for row in _read_csv(path):
            rows.append({"track": track, "organism": organism, **row})
    return rows


def collect_single_snapshot_rows(root):
    rows = []
    pattern = os.path.join(
        root, "04_single_snapshot_ml", "*", "*", "Predict", "model_comparison.csv",
    )
    for path in sorted(glob.glob(pattern)):
        relative = os.path.relpath(path, os.path.join(root, "04_single_snapshot_ml"))
        organism, timepoint, _ = relative.split(os.sep, 2)
        for row in _read_csv(path):
            rows.append({"organism": organism, "timepoint": timepoint, **row})
    return rows


def write_index(output_dir, snapshots, table_names):
    track_counts = Counter(row["track"] for row in snapshots)
    organism_counts = Counter(row["organism"] for row in snapshots)
    lines = [
        "# PhenGO publication result index",
        "",
        "This directory was assembled automatically from strict PhenGO manifests and analysis outputs.",
        "Values in the combined tables retain their organism and experimental-track identifiers.",
        "",
        "## Dataset coverage",
        "",
        f"- ARFF snapshots: {len(snapshots)}",
        "- Tracks: " + ", ".join(f"{key} ({value})" for key, value in sorted(track_counts.items())),
        "- Organisms: " + ", ".join(
            f"{key} ({value})" for key, value in sorted(organism_counts.items())
        ),
        "",
        "## Combined tables",
        "",
    ]
    lines.extend(f"- `{name}`" for name in table_names)
    lines.extend([
        "",
        "`snapshot_inventory.tsv` is the provenance spine: it records class counts, feature counts,",
        "release identifiers, policies, and SHA-256 digests for every ARFF. Statistical interpretation",
        "should use the model and temporal tables together with the pre-registered track definitions.",
        "",
    ])
    with open(os.path.join(output_dir, "README.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        prog="phengo-publication-summary",
        description="Aggregate a completed PhenGO publication run.",
    )
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    root = os.path.abspath(args.run_root)
    output = os.path.abspath(args.output_dir)
    os.makedirs(output, exist_ok=True)
    snapshots = collect_snapshot_rows(root)
    if not snapshots:
        raise SystemExit(f"No PhenGO manifests found below {root}/02_arff")

    tables = {
        "snapshot_inventory.tsv": snapshots,
        "within_year_model_performance.tsv": collect_analysis_rows(
            root, "within_year_cv_summary.csv"
        ),
        "cross_year_transfer_performance.tsv": collect_analysis_rows(
            root, "cross_year_transfer_summary.csv"
        ),
        "dataset_drift.tsv": collect_analysis_rows(root, "dataset_drift.csv"),
        "pairwise_dataset_drift.tsv": collect_analysis_rows(root, "pairwise_drift.csv"),
        "evaluation_preflight.tsv": collect_analysis_rows(root, "evaluation_preflight.csv"),
        "prediction_instability_summary.tsv": collect_analysis_rows(
            root, "prediction_instability_summary.csv"
        ),
        "prediction_instability_by_gene.tsv": collect_analysis_rows(
            root, "prediction_instability.csv"
        ),
        "previous_available_snapshot_label_baseline.tsv": collect_analysis_rows(
            root, "previous_available_snapshot_label_baseline.csv"
        ),
        "feature_rank_overlap.tsv": collect_analysis_rows(root, "feature_rank_overlap.csv"),
        "feature_importance_stability.tsv": collect_analysis_rows(
            root, "feature_importance_stability.csv"
        ),
        "single_snapshot_model_performance.tsv": collect_single_snapshot_rows(root),
        "temporal_summary.tsv": collect_temporal_rows(root),
        "temporal_enrichment_statistics.tsv": collect_temporal_rows(
            root, "*_enrichment_statistics.csv"
        ),
    }
    written = [name for name, rows in tables.items() if _write_tsv(os.path.join(output, name), rows)]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": root,
        "snapshots": len(snapshots),
        "tracks": dict(sorted(Counter(row["track"] for row in snapshots).items())),
        "organisms": dict(sorted(Counter(row["organism"] for row in snapshots).items())),
        "tables": written,
    }
    with open(os.path.join(output, "publication_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_index(output, snapshots, written)


if __name__ == "__main__":
    main()
