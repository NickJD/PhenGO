import json
import sys
from pathlib import Path

import pytest

from PhenGO.core.PhenGO import write_arff_output
from PhenGO.predict.version_sensitivity import (
    load_version_dataset,
    main,
    validate_dataset_manifests,
)
from PhenGO.provenance import sha256_file


def _write_snapshot(directory: Path, snapshot_id: str, shift: int):
    directory.mkdir()
    arff = directory / "worm_PhenGO.arff"
    genes = {}
    for index in range(8):
        genes[f"WBGene{index:04d}"] = {
            "status": "lethal" if index % 2 else "viable",
            "bin_vec": [index % 2, (index + shift) % 3 == 0],
        }
    write_arff_output(genes, ["GO:0000001", "GO:0000002"], str(arff))
    manifest = {
        "schema_version": 2,
        "tool_version": "test",
        "source_tree_sha256": "test-source-tree",
        "dependencies": {"test": "1"},
        "species": "worm",
        "snapshot_id": snapshot_id,
        "strict_snapshot": True,
        "releases": {
            "phenotype": snapshot_id,
            "go_annotations": snapshot_id,
            "go_ontology": snapshot_id,
            "phenotype_ontology": snapshot_id,
            "retrieval_date": "2026-01-01",
        },
        "inputs": {
            "phenotype": {"sha256": "phenotype"},
            "go_annotations": {"sha256": "gaf"},
            "go_ontology": {"sha256": "obo"},
        },
        "policies": {
            "mixed_label": "exclude",
            "nonlethal": "explicit_only",
            "ambiguous": "exclude",
            "label_source": "release_records",
            "legacy_fly_lethal_override": False,
            "filter_multigenes": True,
            "fly_driver_filtering": False,
            "worm_evidence_codes": "all",
            "go_relations": ["is_a", "part_of"],
            "go_namespaces": "all",
            "go_propagation": "ancestors",
            "allow_cross_namespace_go_edges": False,
            "go_evidence_include": "all",
            "go_evidence_exclude": [],
            "identifier_join_precedence": (
                "stable_id > unique canonical_symbol > unique synonym"
            ),
            "strict_identifier_matching": (
                "stable_id > unique canonical_symbol; synonym disabled"
            ),
            "exclude_go_roots": True,
            "filter_unused_gos": True,
            "min_go_gene_count": 2,
            "max_go_gene_fraction": 0.99,
        },
        "outputs": {"arff": {"sha256": sha256_file(str(arff))}},
    }
    (directory / "PhenGO_manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8",
    )
    return arff


def test_version_sensitivity_cli_writes_oof_summaries(tmp_path: Path, monkeypatch):
    first = _write_snapshot(tmp_path / "WS100", "WS100", 0)
    second = _write_snapshot(tmp_path / "WS101", "WS101", 1)
    output = tmp_path / "analysis"
    monkeypatch.setattr(sys, "argv", [
        "phengo-version-sensitivity",
        "-arff_files", str(first), str(second),
        "-dataset_names", "WS100", "WS101",
        "-models", "nn", "lr",
        "-panels", "matched_both",
        "-cv_folds", "2",
        "-cv_repeats", "2",
        "-calibration", "none",
        "-importance_repeats", "2",
        "-top_k", "1",
        "-n_estimators", "5",
        "-nn_hidden_units", "4",
        "-nn_max_iter", "5",
        "-no_nn_early_stopping",
        "-output_dir", str(output),
    ])

    main()

    assert (output / "within_year_cv_summary.csv").is_file()
    assert (output / "cross_year_transfer_summary.csv").is_file()
    assert (output / "prediction_instability.csv").is_file()
    assert (output / "feature_importance_stability.csv").is_file()
    assert (output / "version_sensitivity_manifest.json").is_file()
    matrices = list((output / "transfer_matrices").glob("*.csv"))
    assert matrices
    transfer_text = (output / "cross_year_transfer_summary.csv").read_text(
        encoding="utf-8",
    )
    assert ",nn," in transfer_text


def test_version_analysis_rejects_policy_drift_between_snapshots(tmp_path: Path):
    first = _write_snapshot(tmp_path / "WS100", "WS100", 0)
    second = _write_snapshot(tmp_path / "WS101", "WS101", 1)
    manifest_path = second.parent / "PhenGO_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["policies"]["go_relations"] = ["is_a"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    datasets = [
        load_version_dataset(str(first), "WS100"),
        load_version_dataset(str(second), "WS101"),
    ]
    with pytest.raises(ValueError, match="publication policies differ"):
        validate_dataset_manifests(datasets)
