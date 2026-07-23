import csv
import gzip
import importlib.util
import json
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts/v2/alternative_inputs.py"
SPEC = importlib.util.spec_from_file_location("phengo_v2_alternative_inputs", MODULE_PATH)
alternative_inputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(alternative_inputs)


def gaf_row(stable_id, symbol, go_id, taxon, qualifier="", evidence="IMP", database="FB"):
    return "\t".join([
        database, stable_id, symbol, qualifier, go_id, "PMID:1", evidence, "",
        "P", symbol, "", "protein", f"taxon:{taxon}", "20250101", database,
        "", "",
    ])


def write_gzip(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def test_audit_year_prefers_collection_directory_over_embedded_file_date():
    path = Path("data/fish/gene_association/zfin/2018/gene_association_2017_12_15.zfin.gz")
    assert alternative_inputs.infer_year(path) == "2018"


def test_logical_file_stats_parse_csv_widths_with_quoted_commas(tmp_path):
    path = tmp_path / "assertions.csv.gz"
    write_gzip(path, 'gene,term,note\nMGI:1,MP:1,"one, two"\n')
    stats = alternative_inputs.logical_file_stats(path)
    assert stats["detected_delimiter"] == "comma"
    assert stats["column_widths"] == "3:2"


def test_filter_gaf_is_taxon_specific_and_reproducible(tmp_path):
    source = tmp_path / "source.gaf.gz"
    write_gzip(
        source,
        "!gaf-version: 2.2\n"
        + gaf_row("FBgn1", "a", "GO:0000001", "7227") + "\n"
        + gaf_row("FBgn2", "b", "GO:0000002", "6239") + "\n"
        + gaf_row("WB1", "c", "GO:0000003", "7227", database="WB") + "\n",
    )
    first = tmp_path / "first.gaf.gz"
    second = tmp_path / "second.gaf.gz"
    summary = alternative_inputs.filter_gaf(str(source), str(first), "7227", "FB")
    alternative_inputs.filter_gaf(str(source), str(second), "7227", "FB")

    assert summary["kept_rows"] == 1
    assert summary["excluded_taxon"] == 1
    assert summary["excluded_database"] == 1
    assert alternative_inputs.sha256_file(first) == alternative_inputs.sha256_file(second)
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        data = handle.read()
    assert "FBgn1" in data
    assert "FBgn2" not in data


def test_gaf_comparison_ignores_relation_qualifier_format_change(tmp_path):
    left = tmp_path / "left.gaf.gz"
    right = tmp_path / "right.gaf.gz"
    write_gzip(left, gaf_row("FBgn1", "a", "GO:0000001", "7227", qualifier="") + "\n")
    write_gzip(right, gaf_row("FBgn1", "a", "GO:0000001", "7227", qualifier="involved_in") + "\n")
    result = alternative_inputs.compare_gaf(str(left), str(right), "7227")
    assert result["stable_id_go_jaccard"] == 1.0
    assert result["stable_id_go_shared"] == 1


def test_impc_and_mgi_label_derivation_use_stable_ids_and_lethal_wins(tmp_path):
    terms = tmp_path / "terms.tsv"
    terms.write_text("ID\tDef\nMP:0000001\tlethal\n", encoding="utf-8")
    impc = tmp_path / "impc.csv.gz"
    with gzip.open(impc, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["marker_accession_id", "marker_symbol", "mp_term_id"])
        writer.writeheader()
        writer.writerow({"marker_accession_id": "MGI:1", "marker_symbol": "One", "mp_term_id": "MP:0000002"})
        writer.writerow({"marker_accession_id": "MGI:1", "marker_symbol": "One", "mp_term_id": "MP:0000001"})
        writer.writerow({"marker_accession_id": "MGI:2", "marker_symbol": "Two", "mp_term_id": "MP:0000002"})
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "other.tsv"
    excluded = tmp_path / "excluded.tsv"
    result = alternative_inputs.build_impc_labels(
        str(impc), str(terms), str(lethal), str(viable), str(excluded)
    )
    assert result["lethal_genes"] == 1
    assert result["operational_nonlethal_genes"] == 1
    assert "MGI:1\tOne" in lethal.read_text(encoding="utf-8")
    assert "MGI:1" not in viable.read_text(encoding="utf-8")

    report = tmp_path / "MGI_GenePheno_2020.rpt.gz"
    write_gzip(
        report,
        "geno\tallele\tMGI:allele\tbackground\tMP:0000001\tPMID:1\tMGI:3\t\tgeno-id\n"
        "geno\tallele\tMGI:allele\tbackground\tMP:0000002\tPMID:1\tMGI:4\t\tgeno-id\n",
    )
    result = alternative_inputs.build_mgi_labels(
        str(report), "genepheno", str(terms), str(lethal), str(viable), str(excluded)
    )
    assert result["lethal_genes"] == 1
    assert result["operational_nonlethal_genes"] == 1


def test_audit_marks_logical_duplicates_and_incomplete_files(tmp_path):
    directory = tmp_path / alternative_inputs.MGI_COLLECTIONS["mgi_genepheno"]
    directory.mkdir(parents=True)
    content = "a\tb\tc\td\tMP:1\tf\tMGI:1\th\ti\n"
    (directory / "MGI_GenePheno_2017.rpt").write_text(content, encoding="utf-8")
    write_gzip(directory / "MGI_GenePheno_2017.rpt.gz", content)
    (directory / "MGI_GenePheno_2025.rpt").write_text(content.rstrip("\n"), encoding="utf-8")

    rows = alternative_inputs.audit_inputs(str(tmp_path))
    selected = {Path(row["path"]).name: row for row in rows if row["collection"] == "mgi_genepheno"}
    assert selected["MGI_GenePheno_2017.rpt.gz"]["status"] == "eligible"
    assert selected["MGI_GenePheno_2017.rpt"]["status"] == "excluded_duplicate"
    assert selected["MGI_GenePheno_2025.rpt"]["status"] == "excluded_truncated"


def test_empty_helper_file_is_explicit_and_empty(tmp_path):
    output = tmp_path / "helper.tsv.gz"
    summary = alternative_inputs.write_empty_fly_helper(str(output), "2017")
    assert summary["helper_records"] == 0
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows == [["#record_id", "release", "helper_line_symbol"]]


def test_fly_fail_closed_labels_exclude_multi_gene_and_unknown_contexts(tmp_path):
    source = tmp_path / "fly.tsv.gz"
    write_gzip(
        source,
        "simpleA[a1]\tFBal1\tlethal\tFBrf1\n"
        "simpleB[b1]\tFBal2\tviable\tFBrf2\n"
        "driverTarget[t1]\tFBal3\tviable, with Scer\\GAL4[elav.PLu]\tFBrf3\n"
        "sameGene[x1]\tFBal4\tlethal, with sameGene[x2]\tFBrf4\n"
        "Dcr-2[Scer\\UAS.cDa]\tFBal5\tviable, with Scer\\GAL4[elav.PLu], sphinx2[GD10874]\tFBrf5\n"
        "unknownGene[u1]\tFBal6\tlethal, with environmental context\tFBrf6\n"
        "Scer\\GAL4[pnr-MD237]\tFBal7\tviable, with sphinx2[GD10874]\tFBrf7\n"
        "uasTarget[u1]\tFBal8\tviable, with otherGene[Scer\\UAS.cDa]\tFBrf8\n"
        "spaceTarget[s1] Scer\\GAL4[elav.PLu]\tFBal9\tviable\tFBcv:1\t\t\tFBrf9\n"
        "spaceMulti[m1] otherGene[o1]\tFBal10\tlethal\tFBcv:2\t\t\tFBrf10\n"
        "commaGene[c1]\tFBal11\tlethal, with commaGene[c,2]\tFBrf11\n"
        "ambiguous[a1]\tFBal12\tlethal and viable\tFBrf12\n"
        "semi[s1]\tFBal13\tsemi viable\tFBrf13\n"
        "nonlethal[n1]\tFBal14\tnon-lethal\tFBrf14\n"
        "bracketWith[w1]\tFBal15\tviable note[with text]\tFBrf15\n",
    )
    lethal = tmp_path / "lethal.tsv.gz"
    viable = tmp_path / "viable.tsv.gz"
    excluded = tmp_path / "excluded.tsv.gz"
    audit = tmp_path / "audit.tsv.gz"
    summary = alternative_inputs.build_fly_fail_closed_labels(
        str(source), str(lethal), str(viable), str(excluded), str(audit)
    )

    assert summary["lethal_genes"] == 3
    assert summary["viable_genes"] == 4
    with gzip.open(lethal, "rt", encoding="utf-8") as handle:
        lethal_text = handle.read()
    with gzip.open(viable, "rt", encoding="utf-8") as handle:
        viable_text = handle.read()
    assert "\tsimpleA" in lethal_text
    assert "\tsameGene" in lethal_text
    assert "\tcommaGene" in lethal_text
    assert "\tsimpleB" in viable_text
    assert "\tdriverTarget" in viable_text
    assert "\tspaceTarget" in viable_text
    assert "\tbracketWith" in viable_text
    assert "Dcr-2" not in viable_text
    assert "uasTarget" not in viable_text

    with gzip.open(audit, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    outcomes = {row["raw_primary"]: row["row_outcome"] for row in rows}
    assert outcomes["Dcr-2[Scer\\UAS.cDa]"] == "excluded_multi_gene"
    assert outcomes["unknownGene[u1]"] == "excluded_unresolved_compound"
    assert outcomes["Scer\\GAL4[pnr-MD237]"] == "excluded_accessory_primary"
    assert outcomes["uasTarget[u1]"] == "excluded_multi_gene"
    assert outcomes["spaceTarget[s1] Scer\\GAL4[elav.PLu]"] == "candidate_accessory_context"
    assert outcomes["spaceMulti[m1] otherGene[o1]"] == "excluded_multi_gene"


def test_fly_fail_closed_labels_exclude_gene_level_label_conflicts(tmp_path):
    source = tmp_path / "fly.tsv.gz"
    write_gzip(
        source,
        "conflict[c1]\tFBal1\tlethal\tFBrf1\n"
        "conflict[c2]\tFBal2\tviable\tFBrf2\n"
        "lethalOnly[l1]\tFBal3\tlethal\tFBrf3\n"
        "viableOnly[v1]\tFBal4\tviable\tFBrf4\n",
    )
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    excluded = tmp_path / "excluded.tsv"
    audit = tmp_path / "audit.tsv"
    summary = alternative_inputs.build_fly_fail_closed_labels(
        str(source), str(lethal), str(viable), str(excluded), str(audit)
    )
    assert summary["mixed_genes_excluded"] == 1
    assert "conflict" not in lethal.read_text(encoding="utf-8")
    assert "conflict" not in viable.read_text(encoding="utf-8")
    assert "mixed_lethal_and_viable_observations" in excluded.read_text(encoding="utf-8")
