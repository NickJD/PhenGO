from types import SimpleNamespace

import numpy as np

from PhenGO.predict.data import PhenotypeLabelEncoder
from PhenGO.predict.sklearn_models import build_sklearn_model, train_evaluate_sklearn_model
from PhenGO.predict.version_sensitivity import fit_sklearn_model, model_fit_diagnostics, model_scores


def test_calibrated_model_writes_native_importance_and_schema(tmp_path):
    X = np.asarray([[index % 2, int(index >= 10)] for index in range(20)], dtype=float)
    y = np.asarray([index % 2 for index in range(20)], dtype=int)
    train = np.r_[0:8, 10:18]
    test = np.asarray([8, 9, 18, 19])
    options = SimpleNamespace(
        n_estimators=10,
        max_depth=3,
        n_jobs=1,
        seed=7,
        calibration="sigmoid",
        calibration_cv=2,
        perm_repeats=0,
    )

    result = train_evaluate_sklearn_model(
        "lr",
        X[train],
        y[train],
        X[test],
        y[test],
        [f"gene-{index}" for index in test],
        PhenotypeLabelEncoder(),
        ["GO:1", "GO:2"],
        str(tmp_path),
        options,
    )

    assert result["model_type"] == "lr"
    assert (tmp_path / "lr" / "model.joblib").is_file()
    assert (tmp_path / "lr" / "model_schema.json").is_file()
    assert (tmp_path / "lr" / "feature_importance_native.csv").is_file()


def test_version_sensitivity_nn_uses_same_probability_interface():
    X = np.asarray(
        [[index % 2, (index // 2) % 2, int(index >= 20)] for index in range(40)],
        dtype=float,
    )
    y = np.asarray([index % 2 for index in range(40)], dtype=int)
    options = SimpleNamespace(
        n_estimators=10,
        max_depth=None,
        n_jobs=1,
        seed=13,
        calibration="sigmoid",
        calibration_cv=2,
        nn_hidden_units=[6, 3],
        nn_alpha=0.001,
        nn_batch_size=8,
        nn_learning_rate_init=0.001,
        nn_max_iter=20,
        nn_early_stopping=True,
        nn_validation_fraction=0.2,
        nn_n_iter_no_change=5,
    )

    unfitted = build_sklearn_model("nn", options)
    assert unfitted.hidden_layer_sizes == (6, 3)
    model = fit_sklearn_model("nn", X, y, options, seed=13)
    scores = model_scores(model, X)
    diagnostics = model_fit_diagnostics(model)

    assert scores.shape == (40,)
    assert np.all((scores >= 0) & (scores <= 1))
    assert diagnostics["fit_n_iter_max"] <= options.nn_max_iter
    assert diagnostics["fit_convergence_warnings"] >= 0
