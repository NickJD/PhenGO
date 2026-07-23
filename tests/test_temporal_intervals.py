import csv

from PhenGO.scripts.GO_temporal_analysis import (
    analyse_feature_set_evolution,
    analyse_gene_lifecycle,
    analyse_go_term_lifecycle,
)


def enrichment(total, ratio):
    return {
        "total_genes_with_term": total,
        "enrichment_ratio": ratio,
        "lethal_freq": 0.5,
        "viable_freq": 0.25,
    }


def test_temporal_outputs_use_available_snapshots_and_calendar_gaps(tmp_path):
    datasets = {
        "2017": (
            {"geneA": {"label": "lethal"}},
            ({"GO:1": enrichment(1, 1.0)}, {}),
        ),
        "2024": (
            {"geneA": {"label": "viable"}, "geneB": {"label": "lethal"}},
            ({"GO:1": enrichment(2, 1.7), "GO:2": enrichment(1, 1.0)}, {}),
        ),
    }
    prefix = str(tmp_path / "yeast")
    analyse_feature_set_evolution(prefix, datasets)
    analyse_gene_lifecycle(prefix, datasets)
    analyse_go_term_lifecycle(prefix, datasets)

    with open(f"{prefix}_feature_set_evolution.csv", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["Previous_Available_Timepoint"] == "2017"
    assert rows[1]["Calendar_Gap_Years"] == "7"
    assert rows[1]["Consecutive_Calendar_Years"] == "False"
    assert "Jaccard_Previous_Available" in rows[1]

    with open(f"{prefix}_gene_lifecycle.csv", encoding="utf-8", newline="") as handle:
        gene_rows = list(csv.DictReader(handle))
    assert "Snapshots_Present" in gene_rows[0]
    assert "Years_Present" not in gene_rows[0]

    with open(f"{prefix}_go_term_lifecycle.csv", encoding="utf-8", newline="") as handle:
        go_rows = {row["GO_Term"]: row for row in csv.DictReader(handle)}
    assert go_rows["GO:1"]["Enrichment_Trend_Per_Calendar_Year"] == "0.100000"


def test_go_trend_gives_each_calendar_year_equal_weight(tmp_path):
    datasets = {
        "2017_plain": ({}, ({"GO:1": enrichment(1, 0.0)}, {})),
        "2017_gzip": ({}, ({"GO:1": enrichment(1, 2.0)}, {})),
        "2019": ({}, ({"GO:1": enrichment(1, 3.0)}, {})),
    }
    prefix = str(tmp_path / "mouse")
    analyse_go_term_lifecycle(prefix, datasets)

    with open(f"{prefix}_go_term_lifecycle.csv", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    # Mean 2017 enrichment is 1.0 and 2019 enrichment is 3.0: slope = 1/year.
    assert row["Enrichment_Trend_Per_Calendar_Year"] == "1.000000"
