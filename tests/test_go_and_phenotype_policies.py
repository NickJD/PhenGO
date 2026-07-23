import csv
import gzip
from pathlib import Path
from types import SimpleNamespace

from PhenGO.core.go_handling import (
    get_viability_go_data_fish,
    get_viability_go_data_worm,
)
from PhenGO.core.phenotype_handling import (
    get_viable_inviable_fly,
    get_viable_inviable_fish,
    get_viable_inviable_worm,
    get_viable_inviable_yeast,
    load_gene_set_labels,
)


def _write_gzip_tsv(path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


def test_fly_driver_aware_parser_fails_closed_and_records_row_audit(tmp_path: Path):
    phenotype = tmp_path / "fly.tsv.gz"
    _write_gzip_tsv(phenotype, [
        ["simpleA[a1]", "FBal1", "lethal", "FBrf1"],
        ["simpleB[b1]", "FBal2", "viable", "FBrf2"],
        ["target[t1]", "FBal3", "viable, with Scer\\GAL4[elav.PLu]", "FBrf3"],
        ["same[s1]", "FBal4", "lethal, with same[s2]", "FBrf4"],
        [
            "Dcr-2[Scer\\UAS.cDa]", "FBal5",
            "viable, with Scer\\GAL4[elav.PLu], sphinx2[GD10874]", "FBrf5",
        ],
        ["unknown[u1]", "FBal6", "lethal, with environmental context", "FBrf6"],
        ["uas[u1]", "FBal7", "viable, with other[Scer\\UAS.cDa]", "FBrf7"],
        ["space[s1] Scer\\GAL4[elav.PLu]", "FBal8", "viable", "FBcv:1", "", "", "FBrf8"],
        ["spaceMulti[m1] other[o1]", "FBal9", "lethal", "FBcv:2", "", "", "FBrf9"],
        ["comma[c1]", "FBal10", "lethal, with comma[c,2]", "FBrf10"],
        ["ambiguous[a1]", "FBal11", "lethal and viable", "FBrf11"],
        ["semi[s1]", "FBal12", "semi viable", "FBrf12"],
        ["nonlethal[n1]", "FBal13", "non-lethal", "FBrf13"],
        ["bracketWith[w1]", "FBal14", "viable note[with text]", "FBrf14"],
    ])
    options = SimpleNamespace(
        fly_driver_filtering=True,
        fly_helper_lines=None,
        filter_multigenes=True,
        ambiguous_label_policy="exclude",
        mixed_label_policy="exclude",
        nonlethal_policy="explicit_only",
        filter_mixed_terms=False,
    )

    result = get_viable_inviable_fly(options, str(phenotype))

    assert result == {
        "simpleA": "lethal", "simpleB": "viable",
        "target": "viable", "same": "lethal", "space": "viable",
        "comma": "lethal", "bracketWith": "viable",
    }
    outcomes = {row["primary"]: row["outcome"] for row in options._fly_parse_audit_rows}
    assert outcomes["Dcr-2[Scer\\UAS.cDa]"] == "excluded_multi_gene"
    assert outcomes["unknown[u1]"] == "excluded_unresolved_compound"
    assert outcomes["uas[u1]"] == "excluded_multi_gene"
    assert outcomes["space[s1] Scer\\GAL4[elav.PLu]"] == "retained_resolved_viable"
    assert outcomes["spaceMulti[m1] other[o1]"] == "excluded_multi_gene"
    assert outcomes["ambiguous[a1]"] == "excluded_ambiguous_phenotype"


def test_fly_simple_parser_resolves_same_row_label_conflict_by_policy(tmp_path: Path):
    phenotype = tmp_path / "fly.tsv.gz"
    _write_gzip_tsv(phenotype, [
        ["conflict[c1]", "FBal1", "lethal and viable", "FBrf1"],
        ["lethal[l1]", "FBal2", "lethal", "FBrf2"],
        ["viable[v1]", "FBal3", "viable", "FBrf3"],
    ])
    options = SimpleNamespace(
        fly_driver_filtering=False,
        filter_multigenes=True,
        ambiguous_label_policy="exclude",
        mixed_label_policy="exclude",
        nonlethal_policy="explicit_only",
        filter_mixed_terms=False,
    )

    result = get_viable_inviable_fly(options, str(phenotype))

    assert result == {"lethal": "lethal", "viable": "viable"}


def test_worm_disallowed_evidence_is_excluded_not_relabelled(tmp_path: Path):
    lethal_terms = tmp_path / "lethal.tsv.gz"
    viable_terms = tmp_path / "viable.tsv.gz"
    phenotype = tmp_path / "phenotype.tsv.gz"
    _write_gzip_tsv(lethal_terms, [["WBPhenotype:0001"]])
    _write_gzip_tsv(viable_terms, [["WBPhenotype:0002"]])
    _write_gzip_tsv(phenotype, [
        ["", "", "gene-a", "", "WBPhenotype:0001", "", "IEA"],
        ["", "", "gene-b", "NOT", "WBPhenotype:0001", "", "IMP"],
        ["", "", "gene-c", "", "WBPhenotype:00010", "", "IMP"],
        ["", "", "gene-d", "", "WBPhenotype:0001", "", "IMP"],
        ["", "", "gene-e", "", "WBPhenotype:0002", "", "IMP"],
    ])
    options = SimpleNamespace(
        worm_lethal_genes=None,
        worm_phenotypes=str(lethal_terms),
        worm_viable_phenotypes=str(viable_terms),
        worm_evidence_codes=["IMP"],
        mixed_label_policy="exclude",
        nonlethal_policy="explicit_only",
        filter_mixed_terms=False,
    )

    result = get_viable_inviable_worm(options, str(phenotype))

    assert result == {"gene-d": "lethal", "gene-e": "viable"}


def test_gaf_uses_stable_id_and_exact_not_and_evidence_filters(tmp_path: Path):
    gaf = tmp_path / "worm.gaf"
    rows = [
        ["WB", "WBGene0001", "gene-a", "", "GO:0000001", "", "IMP"],
        ["WB", "WBGene0001", "gene-a", "NOT", "GO:0000002", "", "IMP"],
        ["WB", "WBGene0001", "gene-a", "contributes_to", "GO:0000003", "", "IEA"],
    ]
    with open(gaf, "w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)
    options = SimpleNamespace(
        go_evidence_codes=["IMP"], exclude_go_evidence_codes=[],
    )

    result = get_viability_go_data_worm(
        str(gaf), {"gene-a": "lethal"}, options,
    )

    assert set(result) == {"WBGene0001"}
    assert result["WBGene0001"]["go_list"] == ["GO:0000001"]
    assert result["WBGene0001"]["gene_symbol"] == "gene-a"


def test_gaf_join_prefers_stable_id_then_symbol_over_conflicting_synonym(tmp_path: Path):
    gaf = tmp_path / "fish.gaf.gz"
    _write_gzip_tsv(gaf, [[
        "ZFIN", "ZDB-GENE-1", "wnt5b", "", "GO:0000001", "", "IMP",
        "", "P", "gene", "wnt5a|old-wnt5", "gene_product", "taxon:7955",
    ]])
    options = SimpleNamespace(go_evidence_codes=[], exclude_go_evidence_codes=[])

    stable_result = get_viability_go_data_fish(
        str(gaf),
        {"ZDB-GENE-1": "lethal", "wnt5b": "viable", "wnt5a": "viable"},
        options,
    )
    assert stable_result["ZDB-GENE-1"]["status"] == "lethal"
    assert options._identifier_join_runs[0]["lower_priority_conflict_ids"] == 1

    symbol_result = get_viability_go_data_fish(
        str(gaf), {"wnt5b": "lethal", "wnt5a": "viable"}, options,
    )
    assert symbol_result["ZDB-GENE-1"]["status"] == "lethal"


def test_gaf_excludes_a_symbol_that_maps_to_multiple_stable_ids(tmp_path: Path):
    gaf = tmp_path / "worm.gaf"
    rows = [
        ["WB", "WBGene0001", "dup", "", "GO:0000001", "", "IMP"],
        ["WB", "WBGene0002", "dup", "", "GO:0000002", "", "IMP"],
    ]
    with open(gaf, "w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)
    options = SimpleNamespace(
        go_evidence_codes=[], exclude_go_evidence_codes=[], strict_snapshot=True,
    )

    result = get_viability_go_data_worm(str(gaf), {"dup": "lethal"}, options)

    assert result == {}
    assert options._identifier_join_runs[0]["ambiguous_identifiers"] == 1
    assert options._identifier_join_audit_rows[0]["outcome"] == (
        "excluded_ambiguous_identifier"
    )


def test_strict_gaf_join_does_not_use_historical_synonyms(tmp_path: Path):
    gaf = tmp_path / "fish.gaf.gz"
    _write_gzip_tsv(gaf, [[
        "ZFIN", "ZDB-GENE-1", "current", "", "GO:0000001", "", "IMP",
        "", "P", "gene", "old-name", "gene_product", "taxon:7955",
    ]])
    options = SimpleNamespace(
        go_evidence_codes=[], exclude_go_evidence_codes=[], strict_snapshot=True,
    )

    result = get_viability_go_data_fish(
        str(gaf), {"old-name": "lethal"}, options,
    )

    assert result == {}
    assert options._identifier_join_runs[0]["synonym"] == 0
    assert options._identifier_join_runs[0]["unmatched_identifiers"] == 1


def test_fish_phenotype_labels_use_stable_accessions(tmp_path: Path):
    phenotype = tmp_path / "fish.tsv.gz"
    lethal = ["1", "wnt5b", "ZDB-GENE-1", "", "", "", "", "", "", "", "lethal"]
    viable = ["2", "wnt5a", "ZDB-GENE-2", "", "", "", "", "", "", "", "abnormal"]
    fallback = ["3", "legacy", "", "", "", "", "", "", "", "", "abnormal"]
    _write_gzip_tsv(phenotype, [lethal, viable, fallback])
    options = SimpleNamespace(
        ambiguous_label_policy="exclude",
        mixed_label_policy="lethal_wins",
        nonlethal_policy="observed_viable",
    )

    result = get_viable_inviable_fish(options, str(phenotype))

    assert result == {
        "ZDB-GENE-1": "lethal",
        "ZDB-GENE-2": "viable",
        "legacy": "viable",
    }


def test_release_parsers_keep_source_stable_accessions(tmp_path: Path):
    yeast = tmp_path / "yeast.tsv.gz"
    worm = tmp_path / "worm.tsv.gz"
    lethal_terms = tmp_path / "lethal.tsv.gz"
    _write_gzip_tsv(yeast, [[
        "YAL001C", "ORF", "TFC3", "S000000001", "", "", "null", "", "",
        "inviable",
    ]])
    _write_gzip_tsv(worm, [[
        "WB", "WBGene00000001", "renamed-symbol", "", "WBPhenotype:0001", "",
        "IMP",
    ]])
    _write_gzip_tsv(lethal_terms, [["WBPhenotype:0001"]])

    yeast_result = get_viable_inviable_yeast(SimpleNamespace(
        mixed_label_policy="exclude", nonlethal_policy="explicit_only",
        strict_snapshot=True,
    ), str(yeast))
    worm_result = get_viable_inviable_worm(SimpleNamespace(
        worm_lethal_genes=None, worm_phenotypes=str(lethal_terms),
        worm_viable_phenotypes=None, worm_evidence_codes=None,
        mixed_label_policy="exclude", nonlethal_policy="explicit_only",
        strict_snapshot=True,
    ), str(worm))

    assert yeast_result == {"S000000001": "lethal"}
    assert worm_result == {"WBGene00000001": "lethal"}


def test_paired_gene_sets_do_not_load_symbols_when_stable_ids_exist(tmp_path: Path):
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    lethal.write_text(
        "stable_gene_id\tgene_symbol\nMGI:1\tshared\n", encoding="utf-8",
    )
    viable.write_text(
        "stable_gene_id\tgene_symbol\nMGI:2\tshared\n", encoding="utf-8",
    )

    assert load_gene_set_labels(lethal, viable) == {
        "MGI:1": "lethal", "MGI:2": "viable",
    }


def test_paired_gene_sets_follow_a_declared_stable_id_column(tmp_path: Path):
    lethal = tmp_path / "lethal.tsv"
    viable = tmp_path / "viable.tsv"
    lethal.write_text(
        "gene_symbol\tstable_gene_id\nOne\tMGI:1\n", encoding="utf-8",
    )
    viable.write_text(
        "gene_symbol\tstable_gene_id\nTwo\tMGI:2\n", encoding="utf-8",
    )

    assert load_gene_set_labels(lethal, viable) == {
        "MGI:1": "lethal", "MGI:2": "viable",
    }
