from pathlib import Path

from PhenGO.scripts.publication_summary import (
    collect_analysis_rows,
    collect_single_snapshot_rows,
    collect_temporal_rows,
)


def _write_csv(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_publication_collectors_retain_track_organism_and_timepoint(tmp_path):
    _write_csv(
        tmp_path / "04_ml/primary/yeast/prediction_instability_summary.csv",
        "model,prediction_flip_rate\nnn,0.25\n",
    )
    _write_csv(
        tmp_path / "05_temporal/no_iea/worm/worm_enrichment_statistics.csv",
        "GO_Term,p_adj\nGO:0000001,0.01\n",
    )
    _write_csv(
        tmp_path / "04_single_snapshot_ml/fly/2025/Predict/model_comparison.csv",
        "Model,AUC\nNeural Network,0.8\n",
    )

    analysis = collect_analysis_rows(
        str(tmp_path), "prediction_instability_summary.csv"
    )
    temporal = collect_temporal_rows(
        str(tmp_path), "*_enrichment_statistics.csv"
    )
    single = collect_single_snapshot_rows(str(tmp_path))

    assert analysis == [{
        "track": "primary",
        "organism": "yeast",
        "model": "nn",
        "prediction_flip_rate": "0.25",
    }]
    assert temporal == [{
        "track": "no_iea",
        "organism": "worm",
        "GO_Term": "GO:0000001",
        "p_adj": "0.01",
    }]
    assert single == [{
        "organism": "fly",
        "timepoint": "2025",
        "Model": "Neural Network",
        "AUC": "0.8",
    }]
