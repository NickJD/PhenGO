import csv
import gzip
import sys
from pathlib import Path

from PhenGO.core.PhenGO import write_arff_output
from PhenGO.scripts.get_phenotype_terms import main as phenotype_terms_main
from PhenGO.scripts.publication_inputs import (
    build_fly_assignments,
    build_mouse_impc_sets,
    extract_arff_label_sets,
)


def test_mouse_impc_conversion_excludes_subviable_and_conflicting(tmp_path: Path):
    source = tmp_path / "viability.csv"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Phenotype", "# Genes", "Gene Symbols"])
        writer.writerow(["Homozygous - Viable", "2", "A, B"])
        writer.writerow([])
        writer.writerow(["Gene Symbol", "MGI Gene Id", "Phenotype"])
        writer.writerow(["A", "MGI:1", "Homozygous - Viable"])
        writer.writerow(["B", "MGI:2", "Homozygous - Lethal"])
        writer.writerow(["C", "MGI:3", "Homozygous - Subviable"])
        writer.writerow(["D", "MGI:4", "Homozygous - Viable"])
        writer.writerow(["D", "MGI:4", "Homozygous - Lethal"])

    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    excluded = tmp_path / "excluded.tsv"
    summary = build_mouse_impc_sets(source, lethal, viable, excluded)

    assert summary == {
        "source": str(source.resolve()),
        "lethal": 1,
        "viable": 1,
        "excluded_ambiguous": 2,
    }
    assert "MGI:2\tB" in lethal.read_text(encoding="utf-8")
    assert "MGI:1\tA" in viable.read_text(encoding="utf-8")
    excluded_text = excluded.read_text(encoding="utf-8")
    assert "MGI:3\tC\tambiguous" in excluded_text
    assert "MGI:4\tD\tlethal;viable" in excluded_text


def test_fly_assignments_are_limited_to_dmel_taxon(tmp_path: Path):
    gaf = tmp_path / "fb.gaf.gz"
    with gzip.open(gaf, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["!gaf-version: 2.2"])
        writer.writerow([
            "FB", "FBgn1", "dmel", "", "GO:1", "", "IMP", "", "P",
            "", "", "gene", "taxon:7227", "20250101", "FB",
        ])
        writer.writerow([
            "FB", "FBgn2", "dvir", "", "GO:1", "", "IMP", "", "P",
            "", "", "gene", "taxon:7244", "20250101", "FB",
        ])

    output = tmp_path / "assignments.tsv"
    summary = build_fly_assignments(gaf, output)
    assert summary["assignments"] == 1
    text = output.read_text(encoding="utf-8")
    assert "dmel\tmelanogaster\tCurrent\tFBgn1" in text
    assert "dvir" not in text


def test_fly_assignments_exclude_symbols_with_multiple_stable_ids(tmp_path: Path):
    gaf = tmp_path / "fb.gaf.gz"
    with gzip.open(gaf, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for stable_id in ("FBgn1", "FBgn2"):
            writer.writerow([
                "FB", stable_id, "duplicate", "", "GO:1", "", "IMP", "", "P",
                "", "", "gene", "taxon:7227", "20250101", "FB",
            ])
        writer.writerow([
            "FB", "FBgn3", "unique", "", "GO:1", "", "IMP", "", "P",
            "", "", "gene", "taxon:7227", "20250101", "FB",
        ])

    output = tmp_path / "assignments.tsv"
    excluded = tmp_path / "excluded.tsv"
    summary = build_fly_assignments(gaf, output, excluded_output=excluded)

    assert summary["assignments"] == 1
    assert summary["excluded_ambiguous_symbols"] == 1
    assert "duplicate" not in output.read_text(encoding="utf-8")
    assert "duplicate\tFBgn1|FBgn2\tsymbol_maps_to_multiple_stable_ids" in (
        excluded.read_text(encoding="utf-8")
    )


def test_extract_arff_label_sets_uses_stable_ids_and_symbols(tmp_path: Path):
    arff = tmp_path / "mouse.arff"
    write_arff_output(
        {
            "MGI:1": {"status": "lethal", "bin_vec": [1]},
            "MGI:2": {"status": "viable", "bin_vec": [0]},
        },
        ["GO:1"],
        str(arff),
    )
    identifiers = tmp_path / "gene_identifiers.tsv"
    identifiers.write_text(
        "stable_gene_id\tgene_symbol\nMGI:1\tOne\nMGI:2\tTwo\n",
        encoding="utf-8",
    )
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    summary = extract_arff_label_sets(arff, lethal, viable, identifiers)

    assert summary["lethal"] == 1
    assert summary["viable"] == 1
    assert "MGI:1\tOne" in lethal.read_text(encoding="utf-8")
    assert "MGI:2\tTwo" in viable.read_text(encoding="utf-8")


def test_phenotype_term_traversal_reads_gzip_and_keeps_root(tmp_path: Path, monkeypatch):
    obo = tmp_path / "phenotypes.obo.gz"
    with gzip.open(obo, "wt", encoding="utf-8") as handle:
        handle.write(
            "[Term]\n"
            "id: WBPhenotype:1\n"
            "name: lethal\n\n"
            "[Term]\n"
            "id: WBPhenotype:2\n"
            "name: larval lethal\n"
            "is_a: WBPhenotype:1 ! lethal\n"
        )
    roots = tmp_path / "roots.txt"
    roots.write_text("WBPhenotype:1\n", encoding="utf-8")
    output = tmp_path / "terms.tsv"
    monkeypatch.setattr(sys, "argv", [
        "get-phenotype-terms", "--obo-file", str(obo),
        "--term-list", str(roots), "--results", str(output),
    ])

    phenotype_terms_main()

    text = output.read_text(encoding="utf-8")
    assert "WBPhenotype:1\tlethal" in text
    assert "WBPhenotype:2\tlarval lethal" in text


def test_worm_complete_lethality_roots_match_curator_confirmed_policy():
    repo = Path(__file__).resolve().parents[1]
    root_file = repo / "data/worm/lethal_terms/root_lethal_phenotype_terms.txt"
    roots = {
        line.split("\t", 1)[0]
        for line in root_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }
    assert roots == {
        "WBPhenotype:0000062",
        "WBPhenotype:0000675",
        "WBPhenotype:0001633",
    }


def test_publication_ledger_uses_historical_go_format_eras():
    repo = Path(__file__).resolve().parents[1]
    ledger = repo / "data/publication_snapshots.tsv"
    with ledger.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    anchors = [row for row in rows if row["cohort"] == "annotation_anchor"]
    early = [row for row in anchors if int(row["year"]) <= 2014]
    assert len(anchors) == 31
    assert len(early) == 25
    assert all((repo / row["gaf_file"]).is_file() for row in early)
    assert all((repo / row["go_obo_file"]).is_file() for row in early)

    for row in early:
        year = int(row["year"])
        ontology_path = row["go_obo_file"]
        if year <= 2012:
            assert ontology_path.endswith(f"go_{year}-01-01.obo.gz")
            assert row["go_ontology_release"] == f"GO full ontology {year}-01-01"
        elif year == 2013:
            assert ontology_path.endswith("go-simple_2013-01-01.obo.gz")
            assert row["go_ontology_release"] == "GO simple 2013-01-01"
        else:
            assert ontology_path.endswith("go-basic_2014-01-01.obo.gz")
            assert row["go_ontology_release"] == "GO basic 2014-01-01"
