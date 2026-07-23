from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import PhenGO.predict.version_sensitivity as vs


def _dataset(name, shift=0):
    genes = [f"gene-{index}" for index in range(8)]
    features = pd.DataFrame({
        "GO:1": [(index + shift) % 2 for index in range(8)],
        "GO:2": [int(index >= 4) for index in range(8)],
    }, index=genes, dtype=float)
    labels = pd.Series([index % 2 for index in range(8)], index=genes)
    return vs.VersionDataset(
        name=name,
        path=f"{name}.arff",
        features=features,
        labels=labels,
        label_text=labels.map({0: "viable", 1: "lethal"}),
        feature_names=list(features.columns),
    )


class _DisjointModel:
    def __init__(self, training_genes):
        self.training_genes = set(training_genes)


def test_transfer_and_instability_predictions_are_gene_disjoint(monkeypatch):
    def fake_fit(model_type, X_train, y_train, args, seed, calibrate=True):
        return _DisjointModel(X_train.index)

    def fake_scores(model, X_test):
        assert model.training_genes.isdisjoint(X_test.index)
        return np.full(len(X_test), 0.5)

    monkeypatch.setattr(vs, "fit_sklearn_model", fake_fit)
    monkeypatch.setattr(vs, "model_scores", fake_scores)
    datasets = [_dataset("2019"), _dataset("2020", shift=1)]
    args = SimpleNamespace(cv_folds=2, cv_repeats=2, seed=11)

    rows, summary = vs.run_cross_year_transfer(
        datasets, ["lr"], ["matched_both"], ["GO:1", "GO:2"],
        datasets[0].genes, args,
    )
    per_gene, instability = vs.run_prediction_instability(
        datasets, ["lr"], ["GO:1", "GO:2"], datasets[0].genes, args,
    )

    assert rows and summary and per_gene and instability
    assert {row["evaluation"] for row in rows} == {"fixed_gene_disjoint_oof"}
    assert {row["evaluation"] for row in per_gene} == {"fixed_gene_disjoint_oof"}


def test_publication_preflight_rejects_adaptive_fold_counts():
    datasets = [_dataset("2019"), _dataset("2020", shift=1)]
    args = SimpleNamespace(
        cv_folds=5,
        cv_repeats=2,
        calibration="none",
        calibration_cv=3,
        seed=11,
        allow_incomplete_evaluation=False,
    )

    with pytest.raises(ValueError, match="smallest class has 4 genes"):
        vs.evaluation_preflight(
            datasets, ["full"], ["GO:1", "GO:2"], datasets[0].genes, args,
        )

    args.allow_incomplete_evaluation = True
    rows, errors = vs.evaluation_preflight(
        datasets, ["full"], ["GO:1", "GO:2"], datasets[0].genes, args,
    )
    assert rows and errors
    assert "fail" in {row["status"] for row in rows}


def test_transfer_ignores_go_terms_absent_from_training_release():
    train = _dataset("2019")
    test = _dataset("2020")
    test.features["GO:test-only"] = 1.0
    test.feature_names.append("GO:test-only")

    X_train, _, X_test, _ = vs.pair_panel_data(
        train,
        test,
        "full",
        ["GO:1", "GO:2"],
        train.genes,
    )

    assert list(X_train.columns) == ["GO:1", "GO:2"]
    assert list(X_test.columns) == ["GO:1", "GO:2"]


def test_preflight_checks_nn_early_stopping_split_feasibility():
    datasets = [_dataset("2019"), _dataset("2020")]
    args = SimpleNamespace(
        cv_folds=2,
        cv_repeats=2,
        calibration="none",
        calibration_cv=2,
        seed=11,
        allow_incomplete_evaluation=False,
        nn_early_stopping=True,
        nn_validation_fraction=0.1,
    )

    with pytest.raises(ValueError, match="NN calibration/early-stopping splits"):
        vs.evaluation_preflight(
            datasets,
            ["matched_both"],
            ["GO:1", "GO:2"],
            datasets[0].genes,
            args,
            ["nn"],
        )
