import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

from PhenGO.core.PhenGO import main, write_arff_output
from PhenGO.predict.data import load_arff_data


def _gzip_rows(path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


def test_core_writes_stable_ids_relation_features_and_manifest(tmp_path: Path, monkeypatch):
    phenotype = tmp_path / "phenotypes.tsv.gz"
    gaf = tmp_path / "sgd.gaf.gz"
    obo = tmp_path / "go.obo"
    output = tmp_path / "output"

    phenotype_rows = []
    gaf_rows = []
    for index, status in enumerate(["viable", "inviable", "viable", "inviable"], 1):
        systematic = f"YTEST{index:03d}C"
        phenotype_row = [
            systematic, "ORF", "", f"S{index:010d}", "", "", "null", "", "", status,
        ]
        phenotype_rows.append(phenotype_row)
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

    monkeypatch.setattr(sys, "argv", [
        "PhenGO",
        "-species", "yeast",
        "-phenotype_file", str(phenotype),
        "-gene_association_file", str(gaf),
        "-go_obo_file", str(obo),
        "-output_dir", str(output),
        "-snapshot_id", "yeast-test-1",
        "-phenotype_release", "test-phenotype",
        "-go_annotation_release", "test-gaf",
        "-go_ontology_release", "test-obo",
        "-retrieval_date", "2026-01-01",
        "-strict_snapshot",
        "-include_go_roots",
        "-max_go_gene_fraction", "1.0",
    ])
    main()

    frame, _ = load_arff_data(str(output / "yeast_PhenGO.arff"))
    assert list(frame.iloc[:, 0]) == [f"S{index:010d}" for index in range(1, 5)]
    assert list(frame.columns[1:-1]) == ["GO:0000001", "GO:0000002"]
    assert list(frame.iloc[:, -1]) == ["viable", "lethal", "viable", "lethal"]
    manifest = json.loads((output / "PhenGO_manifest.json").read_text(encoding="utf-8"))
    assert manifest["strict_snapshot"] is True
    assert manifest["policies"]["mixed_label"] == "exclude"
    assert manifest["policies"]["nonlethal"] == "explicit_only"
    assert manifest["outputs"]["arff"]["sha256"]


def test_arff_writer_canonicalises_and_validates_binary_features(tmp_path: Path):
    output = tmp_path / "features.arff"
    write_arff_output(
        {"GENE1": {"status": "viable", "bin_vec": [True, False]}},
        ["GO:0000001", "GO:0000002"],
        str(output),
    )
    frame, _ = load_arff_data(str(output))
    assert frame.iloc[0, 1:-1].tolist() == [1, 0]

    invalid_output = tmp_path / "invalid.arff"
    with pytest.raises(ValueError, match="outside the binary 0/1 domain"):
        write_arff_output(
            {"GENE1": {"status": "lethal", "bin_vec": [0.5]}},
            ["GO:0000001"],
            str(invalid_output),
        )
    assert not invalid_output.exists()
