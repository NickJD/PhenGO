import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).parents[1] / "scripts/v2/publication_multiverse.py"
SPEC = importlib.util.spec_from_file_location("phengo_v2_publication_multiverse", MODULE_PATH)
multiverse = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(multiverse)


def write_arff(path, labels, features=("GO:1", "GO:2")):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["@RELATION test", "@ATTRIBUTE gene STRING"]
    lines.extend(f"@ATTRIBUTE '{feature}' {{0,1}}" for feature in features)
    lines.extend(["@ATTRIBUTE phenotype {viable,lethal}", "@DATA"])
    for index, (gene, label) in enumerate(labels.items()):
        lines.append(f"{gene},{index % 2},{(index + 1) % 2},{label}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_snapshot(root, source_track, organism, timepoint, labels, features=("GO:1", "GO:2")):
    directory = root / "02_arff" / source_track / organism / timepoint
    arff = directory / f"{organism}_PhenGO.arff"
    write_arff(arff, labels, features)
    digest = hashlib.sha256(arff.read_bytes()).hexdigest()
    counts = {
        "genes": len(labels),
        "lethal": sum(value == "lethal" for value in labels.values()),
        "viable": sum(value == "viable" for value in labels.values()),
        "go_features": len(features),
    }
    manifest = {
        "schema_version": 3,
        "snapshot_id": f"{organism}-{timepoint}-{source_track}",
        "counts": counts,
        "policies": {
            "label_source": "release_records", "nonlethal": "explicit_only",
            "mixed_label": "exclude", "ambiguous": "exclude",
            "go_relations": ["is_a", "part_of"],
        },
        "releases": {"phenotype": timepoint, "go_annotations": timepoint, "go_ontology": timepoint},
        "inputs": {},
        "outputs": {"arff": {"sha256": digest}},
        "source_tree_sha256": "source", "tool_version": "test",
        "resource_context": {
            "snapshot_semantics": "declared_composite_cross_section",
            "phenotype_availability": "available_test_fixture",
            "go_annotation_availability": "available_test_fixture",
            "go_ontology_availability": "available_test_fixture",
            "retrieval_route": "generated_test_fixture",
        },
    }
    (directory / "PhenGO_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def write_performance(root, track, organism, dataset, score):
    path = root / "04_ml" / track / organism / "within_year_cv_summary.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset", "model", "panel", "average_precision_mean", "balanced_accuracy_mean"],
        )
        writer.writeheader()
        writer.writerow({
            "dataset": dataset, "model": "lr", "panel": "full",
            "average_precision_mean": score, "balanced_accuracy_mean": score,
        })


def test_arff_contract_and_primary_comparison(tmp_path):
    base = tmp_path / "base"
    extension = tmp_path / "extension"
    write_snapshot(base, "primary", "fly", "2020", {"A": "lethal", "B": "viable"})
    write_snapshot(extension, "alternative", "fly", "2020", {"A": "viable", "C": "viable"}, ("GO:1", "GO:3"))
    base_rows = multiverse.discover_snapshots(base, "base")
    extension_rows = multiverse.discover_snapshots(extension, "extension")
    comparison = multiverse.compare_to_primary(base_rows, extension_rows)[0]
    assert comparison["genes_shared"] == 1
    assert comparison["shared_label_agreement"] == 0.0
    assert comparison["go_features_jaccard"] == 1 / 3


def test_dataset_intervals_distinguish_adjacent_snapshots_from_adjacent_years():
    rows = [{"train_dataset": "2017", "test_dataset": "2024"}]
    multiverse.add_dataset_interval_fields(rows, "train_dataset", "test_dataset")
    assert rows[0]["calendar_gap_years"] == 7
    assert rows[0]["absolute_calendar_gap_years"] == 7
    assert rows[0]["consecutive_calendar_years"] is False


def test_run_integrity_rejects_mixed_source_revisions(tmp_path):
    root = tmp_path / "run"
    (root / "00_run_metadata").mkdir(parents=True)
    (root / "00_run_metadata/run.complete").write_text("complete\n", encoding="utf-8")
    write_snapshot(root, "primary", "fly", "2020", {"A": "lethal"})
    write_snapshot(root, "primary", "fly", "2021", {"B": "viable"})
    manifest_path = root / "02_arff/primary/fly/2021/PhenGO_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_tree_sha256"] = "different-source"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    report = multiverse.audit_scientific_run(root)
    assert report["valid"] is False
    assert report["counts"]["source_tree_hashes"] == 2
    assert any("exactly one source-tree hash" in error for error in report["errors"])


def test_run_integrity_rejects_unverifiable_arff_and_semantics(tmp_path):
    root = tmp_path / "run"
    (root / "00_run_metadata").mkdir(parents=True)
    (root / "00_run_metadata/run.complete").write_text("complete\n", encoding="utf-8")
    write_snapshot(root, "primary", "mouse", "2020", {"A": "lethal"})
    manifest_path = root / "02_arff/primary/mouse/2020/PhenGO_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["resource_context"]["snapshot_semantics"] = "component_alignment_unspecified"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    arff = root / "02_arff/primary/mouse/2020/mouse_PhenGO.arff"
    arff.write_text(arff.read_text(encoding="utf-8") + "% changed\n", encoding="utf-8")

    report = multiverse.audit_scientific_run(root)
    assert report["valid"] is False
    assert any("snapshot semantics" in error for error in report["errors"])
    assert any("ARFF SHA-256" in error for error in report["errors"])


def test_run_integrity_rejects_a_different_launching_source(tmp_path):
    root = tmp_path / "run"
    (root / "00_run_metadata").mkdir(parents=True)
    (root / "00_run_metadata/run.complete").write_text("complete\n", encoding="utf-8")
    write_snapshot(root, "primary", "fish", "2020", {"A": "lethal"})

    report = multiverse.audit_scientific_run(root, expected_source_hash="new-source")
    assert report["valid"] is False
    assert any("launching this workflow" in error for error in report["errors"])


def test_end_to_end_consolidation(tmp_path):
    repo = tmp_path / "repo"
    base = tmp_path / "base"
    extension = tmp_path / "extension"
    output = tmp_path / "consolidated"
    for root in (base, extension):
        metadata = root / "00_run_metadata"
        metadata.mkdir(parents=True)
        (metadata / "run.complete").write_text("complete\n", encoding="utf-8")
    write_snapshot(base, "primary", "mouse", "2020", {"A": "lethal", "B": "viable"})
    write_snapshot(extension, "mouse_impc_assertions", "mouse", "2020", {"A": "lethal", "C": "viable"})
    write_performance(base, "primary", "mouse", "2020", "0.70")
    write_performance(extension, "mouse_impc_assertions", "mouse", "2020", "0.60")

    vocabulary = repo / "data/mouse/phenotype_data/historical_mgi_VOC_MammalianPhenotype"
    vocabulary.mkdir(parents=True)
    (vocabulary / "VOC_MammalianPhenotype_2020.rpt").write_text(
        "MP:1\tterm one\tdefinition\n", encoding="utf-8"
    )
    audit = tmp_path / "audit.tsv"
    audit.write_text("collection\tstatus\trelative_path\nexample\teligible\tfile\n", encoding="utf-8")

    subprocess.run([
        sys.executable, str(MODULE_PATH), "--base-run", str(base),
        "--extension-run", str(extension), "--repo-root", str(repo),
        "--output-dir", str(output), "--input-audit", str(audit),
    ], check=True)
    assert (output / "CONSOLIDATED_REPORT.md").is_file()
    assert (output / "snapshot_inventory.tsv").is_file()
    assert (output / "source_runs/v1/00_run_metadata/run.complete").is_file()
    assert (output / "source_runs/v2_extension/00_run_metadata/run.complete").is_file()
    assert (output / "source_run_file_inventory.tsv").is_file()
    assert (output / "phase_source_hash_inventory.tsv").is_file()
    assert (output / "base_run_integrity.json").is_file()
    assert (output / "extension_run_integrity.json").is_file()
    assert (output / "dataset_drift_all_tracks.csv").is_file()
    assert (output / "previous_available_snapshot_label_baseline_all_tracks.csv").is_file()
    assert (output / "temporal_summary_all_tracks.tsv").is_file()
    deltas = list(csv.DictReader(
        (output / "performance_deltas_vs_primary.csv").open(encoding="utf-8")
    ))
    assert float(deltas[0]["delta_average_precision"]) == pytest.approx(-0.1)
    manifest = json.loads((output / "consolidated_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 2
    assert manifest["counts"]["snapshots"] == 2
    assert manifest["counts"]["copied_source_files"] > 0
    assert manifest["counts"]["phase_source_hash_rows"] == 2
    assert manifest["copy_strategy"] in {"apfs_clone", "copy2"}
    assert manifest["source_runs"]["base"]["copy_verification"]["file_count"] > 0
    assert manifest["source_runs"]["base"]["copy_verification"]["copy_strategy"] == manifest["copy_strategy"]
    assert manifest["scientific_run_integrity"]["base"]["valid"] is True
    assert manifest["scientific_run_integrity"]["extension"]["valid"] is True
    inventory = list(csv.DictReader(
        (output / "snapshot_inventory.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    assert all(str(output / "source_runs") in row["arff_path"] for row in inventory)
    phase_hashes = list(csv.DictReader(
        (output / "phase_source_hash_inventory.tsv").open(encoding="utf-8"),
        delimiter="\t",
    ))
    assert {row["source_run"] for row in phase_hashes} == {"base", "extension"}
    assert {row["source_tree_sha256"] for row in phase_hashes} == {"source"}
