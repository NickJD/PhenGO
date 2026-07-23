import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

from PhenGO.core.PhenGO import main
from PhenGO.core.phenotype_handling import load_gene_set_labels
from PhenGO.predict.data import load_arff_data


def _gzip_rows(path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


def _write_rows(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


def _write_common_inputs(tmp_path: Path, statuses):
    phenotype = tmp_path / "phenotypes.tsv.gz"
    gaf = tmp_path / "sgd.gaf.gz"
    obo = tmp_path / "go-basic.obo"
    phenotype_rows = []
    gaf_rows = []
    for index, status in enumerate(statuses, 1):
        systematic = f"YTEST{index:03d}C"
        phenotype_rows.append([
            systematic, "ORF", "", f"S{index:010d}", "", "", "null", "", "", status,
        ])
        gaf_row = ["SGD", f"S{index:010d}", f"SYM{index}", "", "GO:0000002", "", "IMP"]
        gaf_row.extend(["", "", "", systematic])
        gaf_rows.append(gaf_row)
    _gzip_rows(phenotype, phenotype_rows)
    _gzip_rows(gaf, gaf_rows)
    obo.write_text(
        """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: child
namespace: biological_process
is_a: GO:0000001 ! root
""",
        encoding="utf-8",
    )
    return phenotype, gaf, obo


def _base_args(gaf, obo, output):
    return [
        "PhenGO",
        "-species", "yeast",
        "-gene_association_file", str(gaf),
        "-go_obo_file", str(obo),
        "-output_dir", str(output),
        "-snapshot_id", output.name,
        "-phenotype_release", "test-labels",
        "-go_annotation_release", "test-gaf",
        "-go_ontology_release", "test-go-basic",
        "-retrieval_date", "2026-01-01",
        "-strict_snapshot",
        "-include_go_roots",
        "-max_go_gene_fraction", "1.0",
    ]


def test_gene_set_mode_is_complete_and_does_not_require_phenotype_file(
        tmp_path: Path, monkeypatch):
    _, gaf, obo = _write_common_inputs(
        tmp_path, ["viable", "inviable", "viable", "inviable"],
    )
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    _write_rows(lethal, [["stable_gene_id"], ["S0000000002"], ["S0000000004"]])
    _write_rows(viable, [["stable_gene_id"], ["S0000000001"], ["S0000000003"]])
    output = tmp_path / "gene-set-output"
    monkeypatch.setattr(sys, "argv", _base_args(gaf, obo, output) + [
        "-label_source", "gene_sets",
        "-lethal_gene_set", str(lethal),
        "-viable_gene_set", str(viable),
    ])

    main()

    frame, _ = load_arff_data(str(output / "yeast_PhenGO.arff"))
    assert frame.iloc[:, -1].tolist() == ["viable", "lethal", "viable", "lethal"]
    manifest = json.loads((output / "PhenGO_manifest.json").read_text(encoding="utf-8"))
    assert manifest["policies"]["label_source"] == "gene_sets"
    assert manifest["policies"]["go_relations"] == ["is_a", "part_of"]
    assert manifest["inputs"]["lethal_gene_set"]["sha256"]
    audit = (output / "label_source_audit.tsv").read_text(encoding="utf-8")
    assert audit.count("retained_gene_set") == 4


def test_agreement_mode_excludes_disagreement_and_single_source_labels(
        tmp_path: Path, monkeypatch):
    phenotype, gaf, obo = _write_common_inputs(
        tmp_path,
        ["viable", "inviable", "viable", "inviable", "viable", "inviable"],
    )
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    _write_rows(lethal, [["S0000000002"], ["S0000000003"], ["S0000000004"]])
    _write_rows(viable, [["S0000000001"], ["S0000000005"]])
    output = tmp_path / "agreement-output"
    monkeypatch.setattr(sys, "argv", _base_args(gaf, obo, output) + [
        "-phenotype_file", str(phenotype),
        "-label_source", "agreement",
        "-lethal_gene_set", str(lethal),
        "-viable_gene_set", str(viable),
    ])

    main()

    frame, _ = load_arff_data(str(output / "yeast_PhenGO.arff"))
    assert set(frame.iloc[:, 0]) == {
        "S0000000001", "S0000000002", "S0000000004", "S0000000005",
    }
    audit = (output / "label_source_audit.tsv").read_text(encoding="utf-8")
    assert audit.count("retained_agreement") == 4
    assert audit.count("excluded_disagreement") == 1
    assert audit.count("excluded_release_only") == 1


def test_paired_gene_sets_reject_conflicting_identifiers(tmp_path: Path):
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    _write_rows(lethal, [["FBgn0000001"]])
    _write_rows(viable, [["FBgn0000001"]])

    with pytest.raises(ValueError, match="labels conflict"):
        load_gene_set_labels(str(lethal), str(viable))
