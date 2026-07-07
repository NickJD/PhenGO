"""Version-sensitivity analysis for PhenGO ARFF datasets.

This module is designed for the paper-level question:

    Do machine-learning outputs change when the same organism is represented by
    different yearly snapshots of model-organism database resources?

The command intentionally separates several effects:

* within-year model performance
* train-year -> test-year transfer performance
* gene-set churn
* GO-feature-set churn
* label churn among shared genes
* prediction instability for the same genes across years
* shared-feature and matched-gene sensitivity panels
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from types import SimpleNamespace

from ..constants import configure_logger

logger = logging.getLogger("PhenGO.version_sensitivity")

SKLEARN_MODELS = ["dt", "rf", "gb", "lr", "svm"]
DEFAULT_MODELS = ["lr", "rf", "gb", "dt"]
PANELS = ["full", "matched_features", "matched_genes", "matched_both"]
MATRIX_METRICS = ["roc_auc", "average_precision", "f1_lethal", "balanced_accuracy", "mcc"]


@dataclass
class VersionDataset:
    name: str
    path: str
    features: object
    labels: object
    label_text: object
    feature_names: list[str]

    @property
    def genes(self):
        return list(self.features.index)


def natural_key(value: str):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(value))]


def expand_models(model_args: list[str]) -> list[str]:
    out = []
    for model in model_args:
        targets = SKLEARN_MODELS if model == "all" else [model]
        for target in targets:
            if target not in SKLEARN_MODELS:
                raise ValueError(f"Unknown model '{target}'. Valid choices: {SKLEARN_MODELS + ['all']}")
            if target not in out:
                out.append(target)
    return out


def discover_arff_files(input_dir: str) -> tuple[list[str], list[str]]:
    """Discover ARFF datasets from a parent directory.

    Each direct child directory containing ``*_PhenGO.arff`` is treated as one
    dataset. If no child datasets are found, direct ``*.arff`` files in
    ``input_dir`` are used.
    """
    input_dir = os.path.abspath(input_dir)
    discovered = []
    for entry in sorted(os.listdir(input_dir), key=natural_key):
        subdir = os.path.join(input_dir, entry)
        if not os.path.isdir(subdir):
            continue
        arffs = [f for f in os.listdir(subdir) if f.endswith("_PhenGO.arff")]
        if arffs:
            discovered.append((entry, os.path.join(subdir, sorted(arffs, key=natural_key)[0])))

    if not discovered:
        for entry in sorted(os.listdir(input_dir), key=natural_key):
            if entry.endswith(".arff"):
                name = os.path.splitext(entry)[0]
                discovered.append((name, os.path.join(input_dir, entry)))

    names = [name for name, _ in discovered]
    paths = [path for _, path in discovered]
    return paths, names


def ensure_output_dir(path: str, overwrite: bool):
    if os.path.exists(path) and os.listdir(path):
        if not overwrite:
            raise SystemExit(
                f"Output directory '{path}' is not empty. Choose a new directory or pass -overwrite."
            )
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def load_version_dataset(path: str, name: str) -> VersionDataset:
    from .data import PhenotypeLabelEncoder, load_arff_data
    import pandas as pd

    df, _ = load_arff_data(path)
    if df is None:
        raise ValueError(f"Could not load ARFF dataset: {path}")

    gene_col = df.columns[0]
    label_col = df.columns[-1]
    feature_names = list(df.columns[1:-1])
    encoder = PhenotypeLabelEncoder()

    features = df.loc[:, feature_names].astype(float)
    features.index = df[gene_col].astype(str)
    labels = pd.Series(encoder.fit_transform(df[label_col]), index=features.index, name="label")
    label_text = pd.Series(df[label_col].astype(str).str.lower().values,
                           index=features.index, name="label_text")

    if features.index.has_duplicates:
        n_dupes = int(features.index.duplicated().sum())
        logger.warning(f"{name}: {n_dupes} duplicate gene rows found; keeping first occurrence")
        keep = ~features.index.duplicated(keep="first")
        features = features.loc[keep]
        labels = labels.loc[keep]
        label_text = label_text.loc[keep]

    return VersionDataset(
        name=name,
        path=path,
        features=features,
        labels=labels,
        label_text=label_text,
        feature_names=feature_names,
    )


def shared_features(datasets: list[VersionDataset]) -> list[str]:
    common = set(datasets[0].feature_names)
    for dataset in datasets[1:]:
        common &= set(dataset.feature_names)
    return sorted(common)


def shared_genes(datasets: list[VersionDataset]) -> list[str]:
    common = set(datasets[0].genes)
    for dataset in datasets[1:]:
        common &= set(dataset.genes)
    return sorted(common, key=natural_key)


def panel_data(dataset: VersionDataset, panel: str,
               all_shared_features: list[str], all_shared_genes: list[str]):
    features = dataset.feature_names if panel in ("full", "matched_genes") else all_shared_features
    genes = dataset.genes if panel in ("full", "matched_features") else all_shared_genes
    X = dataset.features.reindex(index=genes, columns=features, fill_value=0).astype(float)
    y = dataset.labels.reindex(index=genes)
    valid = y.notna()
    return X.loc[valid], y.loc[valid].astype(int)


def pair_panel_data(train: VersionDataset, test: VersionDataset, panel: str,
                    all_shared_features: list[str], all_shared_genes: list[str]):
    if panel == "full":
        features = sorted(set(train.feature_names) | set(test.feature_names))
        train_genes = train.genes
        test_genes = test.genes
    elif panel == "matched_features":
        features = all_shared_features
        train_genes = train.genes
        test_genes = test.genes
    elif panel == "matched_genes":
        features = sorted(set(train.feature_names) | set(test.feature_names))
        common = sorted(set(train.genes) & set(test.genes), key=natural_key)
        train_genes = common
        test_genes = common
    elif panel == "matched_both":
        features = all_shared_features
        common = all_shared_genes
        train_genes = common
        test_genes = common
    else:
        raise ValueError(f"Unknown panel: {panel}")

    X_train = train.features.reindex(index=train_genes, columns=features, fill_value=0).astype(float)
    y_train = train.labels.reindex(index=train_genes)
    X_test = test.features.reindex(index=test_genes, columns=features, fill_value=0).astype(float)
    y_test = test.labels.reindex(index=test_genes)

    train_valid = y_train.notna()
    test_valid = y_test.notna()
    return (
        X_train.loc[train_valid],
        y_train.loc[train_valid].astype(int),
        X_test.loc[test_valid],
        y_test.loc[test_valid].astype(int),
    )


def can_fit(y) -> bool:
    return len(set(map(int, y))) >= 2


def metric_value(value):
    if value is None:
        return float("nan")
    try:
        if math.isnan(float(value)):
            return float("nan")
    except Exception:
        pass
    return float(value)


def compute_metrics(y_true, y_score, threshold: float = 0.5) -> dict:
    import numpy as np
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        balanced_accuracy_score,
        brier_score_loss,
        f1_score,
        matthews_corrcoef,
        precision_score,
        recall_score,
        roc_auc_score,
    )

    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    y_pred = (y_score >= threshold).astype(int)
    one_class = len(np.unique(y_true)) < 2

    return {
        "n": int(len(y_true)),
        "n_lethal": int(np.sum(y_true == 1)),
        "n_viable": int(np.sum(y_true == 0)),
        "prevalence_lethal": float(np.mean(y_true == 1)) if len(y_true) else float("nan"),
        "accuracy": float(accuracy_score(y_true, y_pred)) if len(y_true) else float("nan"),
        "balanced_accuracy": (
            float(balanced_accuracy_score(y_true, y_pred)) if not one_class else float("nan")
        ),
        "precision_lethal": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall_lethal": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1_lethal": float(f1_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "f1_macro": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, y_pred)) if not one_class else float("nan"),
        "roc_auc": float(roc_auc_score(y_true, y_score)) if not one_class else float("nan"),
        "average_precision": (
            float(average_precision_score(y_true, y_score)) if not one_class else float("nan")
        ),
        "brier_score": float(brier_score_loss(y_true, y_score)) if len(y_true) else float("nan"),
    }


def make_options(args, seed: int):
    return SimpleNamespace(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        random_state=seed,
    )


def fit_sklearn_model(model_type: str, X_train, y_train, args, seed: int):
    from sklearn.utils.class_weight import compute_sample_weight

    from .sklearn_models import build_sklearn_model

    model = build_sklearn_model(model_type, make_options(args, seed))
    if model_type == "gb":
        sample_weight = compute_sample_weight(class_weight="balanced", y=y_train)
        model.fit(X_train, y_train, sample_weight=sample_weight)
    else:
        model.fit(X_train, y_train)
    return model


def model_scores(model, X):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    if hasattr(model, "decision_function"):
        import numpy as np

        scores = model.decision_function(X)
        return 1 / (1 + np.exp(-scores))
    return model.predict(X)


def baseline_scores(kind: str, y_train, n_test: int, seed: int):
    import numpy as np

    prevalence = float(np.mean(y_train == 1)) if len(y_train) else 0.0
    if kind == "baseline_majority":
        return np.full(n_test, 1.0 if prevalence >= 0.5 else 0.0)
    if kind == "baseline_prior":
        return np.full(n_test, prevalence)
    if kind == "baseline_stratified_random":
        rng = np.random.default_rng(seed)
        return rng.binomial(1, prevalence, size=n_test).astype(float)
    raise ValueError(f"Unknown baseline: {kind}")


def dataset_drift_rows(datasets: list[VersionDataset]) -> list[dict]:
    rows = []
    for dataset in datasets:
        X = dataset.features
        y = dataset.labels
        active_features = [col for col in dataset.feature_names if float(X[col].sum()) > 0]
        go_per_gene = X.sum(axis=1)
        rows.append({
            "dataset": dataset.name,
            "path": dataset.path,
            "n_genes": int(X.shape[0]),
            "n_features": int(X.shape[1]),
            "n_active_features": int(len(active_features)),
            "n_lethal": int((y == 1).sum()),
            "n_viable": int((y == 0).sum()),
            "prevalence_lethal": float((y == 1).mean()) if len(y) else float("nan"),
            "mean_go_terms_per_gene": float(go_per_gene.mean()) if len(go_per_gene) else 0.0,
            "median_go_terms_per_gene": float(go_per_gene.median()) if len(go_per_gene) else 0.0,
            "sparsity": float((X.values == 0).sum() / X.size) if X.size else float("nan"),
        })
    return rows


def pairwise_drift_rows(datasets: list[VersionDataset]) -> list[dict]:
    rows = []
    for a in datasets:
        for b in datasets:
            genes_a, genes_b = set(a.genes), set(b.genes)
            feats_a, feats_b = set(a.feature_names), set(b.feature_names)
            common_genes = sorted(genes_a & genes_b, key=natural_key)
            label_churn = float("nan")
            if common_genes:
                label_churn = float((a.labels.loc[common_genes].values != b.labels.loc[common_genes].values).mean())
            rows.append({
                "dataset_a": a.name,
                "dataset_b": b.name,
                "n_common_genes": len(common_genes),
                "gene_jaccard": len(genes_a & genes_b) / len(genes_a | genes_b) if genes_a | genes_b else float("nan"),
                "n_common_features": len(feats_a & feats_b),
                "feature_jaccard": len(feats_a & feats_b) / len(feats_a | feats_b) if feats_a | feats_b else float("nan"),
                "label_churn_common_genes": label_churn,
            })
    return rows


def run_within_year_cv(datasets: list[VersionDataset], models: list[str], panels: list[str],
                       all_features: list[str], all_genes: list[str], args) -> tuple[list[dict], list[dict]]:
    import numpy as np
    import pandas as pd
    from sklearn.model_selection import RepeatedStratifiedKFold

    fold_rows = []
    for panel in panels:
        for dataset in datasets:
            X, y = panel_data(dataset, panel, all_features, all_genes)
            if not can_fit(y):
                logger.warning(f"Skipping CV for {dataset.name}/{panel}: only one class present")
                continue
            min_class = int(np.bincount(y).min())
            n_splits = min(args.cv_folds, min_class)
            if n_splits < 2:
                logger.warning(f"Skipping CV for {dataset.name}/{panel}: fewer than 2 samples in one class")
                continue
            splitter = RepeatedStratifiedKFold(
                n_splits=n_splits,
                n_repeats=args.cv_repeats,
                random_state=args.seed,
            )
            for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(X, y), 1):
                for model_type in models:
                    seed = args.seed + fold_idx
                    model = fit_sklearn_model(
                        model_type,
                        X.iloc[train_idx],
                        y.iloc[train_idx],
                        args,
                        seed,
                    )
                    scores = model_scores(model, X.iloc[test_idx])
                    metrics = compute_metrics(y.iloc[test_idx], scores)
                    fold_rows.append({
                        "analysis": "within_year_cv",
                        "panel": panel,
                        "dataset": dataset.name,
                        "model": model_type,
                        "fold": fold_idx,
                        **metrics,
                    })

                for baseline in ["baseline_majority", "baseline_prior", "baseline_stratified_random"]:
                    scores = baseline_scores(baseline, y.iloc[train_idx], len(test_idx), args.seed + fold_idx)
                    metrics = compute_metrics(y.iloc[test_idx], scores)
                    fold_rows.append({
                        "analysis": "within_year_cv",
                        "panel": panel,
                        "dataset": dataset.name,
                        "model": baseline,
                        "fold": fold_idx,
                        **metrics,
                    })

    if not fold_rows:
        return [], []

    metric_cols = [c for c in fold_rows[0].keys()
                   if c not in {"analysis", "panel", "dataset", "model", "fold"}]
    df = pd.DataFrame(fold_rows)
    summary = []
    for keys, group in df.groupby(["panel", "dataset", "model"], dropna=False):
        row = {"panel": keys[0], "dataset": keys[1], "model": keys[2], "n_folds": int(len(group))}
        for metric in metric_cols:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if len(group) > 1 else 0.0
        summary.append(row)
    return fold_rows, summary


def run_cross_year_transfer(datasets: list[VersionDataset], models: list[str], panels: list[str],
                            all_features: list[str], all_genes: list[str], args) -> list[dict]:
    rows = []
    for panel in panels:
        for train in datasets:
            for test in datasets:
                X_train, y_train, X_test, y_test = pair_panel_data(
                    train, test, panel, all_features, all_genes,
                )
                if len(X_train) == 0 or len(X_test) == 0:
                    continue

                for baseline in ["baseline_majority", "baseline_prior", "baseline_stratified_random"]:
                    scores = baseline_scores(
                        baseline, y_train, len(y_test), args.seed + len(rows)
                    )
                    metrics = compute_metrics(y_test, scores)
                    rows.append({
                        "analysis": "cross_year_transfer",
                        "panel": panel,
                        "train_dataset": train.name,
                        "test_dataset": test.name,
                        "model": baseline,
                        **metrics,
                    })

                if not can_fit(y_train):
                    logger.warning(f"Skipping transfer from {train.name}/{panel}: only one training class")
                    continue
                for model_type in models:
                    model = fit_sklearn_model(model_type, X_train, y_train, args, args.seed)
                    scores = model_scores(model, X_test)
                    metrics = compute_metrics(y_test, scores)
                    rows.append({
                        "analysis": "cross_year_transfer",
                        "panel": panel,
                        "train_dataset": train.name,
                        "test_dataset": test.name,
                        "model": model_type,
                        **metrics,
                    })
    return rows


def run_previous_year_label_baseline(datasets: list[VersionDataset]) -> list[dict]:
    rows = []
    ordered = sorted(datasets, key=lambda ds: natural_key(ds.name))
    for prev, current in zip(ordered, ordered[1:]):
        common = sorted(set(prev.genes) & set(current.genes), key=natural_key)
        if not common:
            continue
        scores = prev.labels.loc[common].astype(float).values
        metrics = compute_metrics(current.labels.loc[common], scores)
        rows.append({
            "analysis": "previous_year_label_baseline",
            "train_dataset": prev.name,
            "test_dataset": current.name,
            "model": "previous_year_labels",
            **metrics,
        })
    return rows


def native_importance(model, model_type: str, feature_names: list[str]) -> dict[str, float] | None:
    import numpy as np

    if model_type in {"dt", "rf", "gb"} and hasattr(model, "feature_importances_"):
        values = model.feature_importances_
    elif model_type == "lr" and hasattr(model, "coef_"):
        values = np.abs(model.coef_[0])
    else:
        return None
    return {feature: float(value) for feature, value in zip(feature_names, values)}


def run_feature_rank_overlap(datasets: list[VersionDataset], models: list[str],
                             all_features: list[str], all_genes: list[str], args) -> list[dict]:
    rows = []
    supported = [m for m in models if m in {"dt", "rf", "gb", "lr"}]
    if not supported or not all_features:
        return rows

    top_sets = {}
    for dataset in datasets:
        X, y = panel_data(dataset, "matched_features", all_features, dataset.genes)
        if not can_fit(y):
            continue
        for model_type in supported:
            model = fit_sklearn_model(model_type, X, y, args, args.seed)
            imp = native_importance(model, model_type, list(X.columns))
            if not imp:
                continue
            ranked = sorted(imp.items(), key=lambda item: item[1], reverse=True)
            top_sets[(dataset.name, model_type)] = {feature for feature, _ in ranked[:args.top_k]}

    for model_type in supported:
        for a in datasets:
            for b in datasets:
                set_a = top_sets.get((a.name, model_type), set())
                set_b = top_sets.get((b.name, model_type), set())
                rows.append({
                    "model": model_type,
                    "dataset_a": a.name,
                    "dataset_b": b.name,
                    "top_k": args.top_k,
                    "n_overlap": len(set_a & set_b),
                    "jaccard": len(set_a & set_b) / len(set_a | set_b) if set_a | set_b else float("nan"),
                })
    return rows


def run_prediction_instability(datasets: list[VersionDataset], models: list[str],
                               all_features: list[str], all_genes: list[str], args) -> tuple[list[dict], list[dict]]:
    import numpy as np
    import pandas as pd

    if not all_features or not all_genes:
        return [], []

    per_gene_rows = []
    summary_rows = []
    ordered = sorted(datasets, key=lambda ds: natural_key(ds.name))
    for model_type in models:
        score_table = {}
        label_table = {}
        pred_table = {}
        for dataset in ordered:
            X, y = panel_data(dataset, "matched_both", all_features, all_genes)
            if not can_fit(y):
                continue
            model = fit_sklearn_model(model_type, X, y, args, args.seed)
            scores = model_scores(model, X)
            score_table[dataset.name] = pd.Series(scores, index=X.index)
            label_table[dataset.name] = y
            pred_table[dataset.name] = pd.Series((scores >= 0.5).astype(int), index=X.index)

        if len(score_table) < 2:
            continue

        scores_df = pd.DataFrame(score_table).loc[all_genes]
        labels_df = pd.DataFrame(label_table).loc[all_genes]
        preds_df = pd.DataFrame(pred_table).loc[all_genes]
        for gene in all_genes:
            score_values = scores_df.loc[gene].astype(float)
            pred_values = preds_df.loc[gene].astype(int)
            label_values = labels_df.loc[gene].astype(int)
            per_gene_rows.append({
                "model": model_type,
                "gene": gene,
                "mean_probability_lethal": float(score_values.mean()),
                "sd_probability_lethal": float(score_values.std(ddof=1)) if len(score_values) > 1 else 0.0,
                "range_probability_lethal": float(score_values.max() - score_values.min()),
                "prediction_flip": int(len(set(pred_values)) > 1),
                "label_flip": int(len(set(label_values)) > 1),
                "min_probability_lethal": float(score_values.min()),
                "max_probability_lethal": float(score_values.max()),
            })

        summary_rows.append({
            "model": model_type,
            "panel": "matched_both",
            "n_genes": len(all_genes),
            "mean_gene_probability_sd": float(np.nanmean(scores_df.std(axis=1, ddof=1))),
            "median_gene_probability_sd": float(np.nanmedian(scores_df.std(axis=1, ddof=1))),
            "mean_gene_probability_range": float(np.nanmean(scores_df.max(axis=1) - scores_df.min(axis=1))),
            "prediction_flip_rate": float((preds_df.nunique(axis=1) > 1).mean()),
            "label_flip_rate": float((labels_df.nunique(axis=1) > 1).mean()),
        })
    return per_gene_rows, summary_rows


def write_csv(path: str, rows: list[dict]):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    preferred = [
        "analysis", "panel", "dataset", "train_dataset", "test_dataset", "model",
        "fold", "gene", "metric",
    ]
    fieldnames = [f for f in preferred if f in fieldnames] + [
        f for f in fieldnames if f not in preferred
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_transfer_matrices(output_dir: str, transfer_rows: list[dict],
                            dataset_names: list[str], metrics: list[str] = MATRIX_METRICS):
    import pandas as pd

    if not transfer_rows:
        return
    matrix_dir = os.path.join(output_dir, "transfer_matrices")
    os.makedirs(matrix_dir, exist_ok=True)
    df = pd.DataFrame(transfer_rows)
    for panel in sorted(df["panel"].dropna().unique(), key=natural_key):
        for model in sorted(df["model"].dropna().unique(), key=natural_key):
            subset = df[(df["panel"] == panel) & (df["model"] == model)]
            for metric in metrics:
                if metric not in subset.columns:
                    continue
                pivot = subset.pivot(index="train_dataset", columns="test_dataset", values=metric)
                pivot = pivot.reindex(index=dataset_names, columns=dataset_names)
                out = os.path.join(matrix_dir, f"{panel}_{model}_{metric}.csv")
                pivot.to_csv(out)


def write_markdown_report(output_dir: str, args, datasets: list[VersionDataset],
                          models: list[str], panels: list[str]):
    report_path = os.path.join(output_dir, "VERSION_SENSITIVITY_REPORT.md")
    lines = [
        "# PhenGO Version-Sensitivity Analysis Report",
        "",
        f"Generated: {datetime.now().isoformat()}",
        "",
        "## Purpose",
        "",
        "This analysis is designed to test whether gene-essentiality ML outputs are stable "
        "across yearly snapshots of model-organism database resources. It separates "
        "performance changes caused by gene coverage, phenotype-label changes, GO feature "
        "churn, and model behavior.",
        "",
        "## Inputs",
        "",
    ]
    for dataset in datasets:
        lines.append(f"- `{dataset.name}`: `{dataset.path}`")
    lines.extend([
        "",
        "## Configuration",
        "",
        f"- Models: {', '.join(models)}",
        f"- Panels: {', '.join(panels)}",
        f"- CV folds: {args.cv_folds}",
        f"- CV repeats: {args.cv_repeats}",
        f"- Random seed: {args.seed}",
        f"- Top-k feature overlap: {args.top_k}",
        "",
        "## Output Files",
        "",
        "- `dataset_drift.csv`: per-year gene counts, GO feature counts, sparsity, and class prevalence.",
        "- `pairwise_drift.csv`: pairwise gene-set Jaccard, GO-feature Jaccard, and label churn among shared genes.",
        "- `within_year_cv.csv`: repeated stratified CV fold-level metrics for each dataset, model, and panel.",
        "- `within_year_cv_summary.csv`: mean and standard deviation for CV metrics.",
        "- `cross_year_transfer.csv`: train-year -> test-year metrics for every pair of datasets.",
        "- `transfer_matrices/`: matrix-form CSVs for key metrics, ready for heatmaps.",
        "- `previous_year_label_baseline.csv`: adjacent-year baseline using the previous year's labels for shared genes.",
        "- `prediction_instability.csv`: per-gene probability variance and prediction flips across years.",
        "- `prediction_instability_summary.csv`: model-level summary of prediction instability.",
        "- `feature_rank_overlap.csv`: pairwise top-k GO importance overlap for models with native importances.",
        "",
        "## Panels",
        "",
        "- `full`: each transfer uses all features from the train/test pair, filling absent GO features with 0.",
        "- `matched_features`: all datasets are restricted to GO features shared by every dataset.",
        "- `matched_genes`: train/test pairs are restricted to genes shared by that pair.",
        "- `matched_both`: all datasets are restricted to genes and GO features shared by every dataset.",
        "",
        "## Recommended Article Use",
        "",
        "1. Use `dataset_drift.csv` and `pairwise_drift.csv` to establish that the source resources themselves change.",
        "2. Use `within_year_cv_summary.csv` to show that ordinary within-snapshot ML performance varies by year.",
        "3. Use `cross_year_transfer.csv` and `transfer_matrices/` as the central evidence: models trained on one snapshot "
        "do not necessarily transfer cleanly to another snapshot.",
        "4. Compare `full`, `matched_features`, `matched_genes`, and `matched_both` panels to separate feature churn, gene churn, "
        "and label/annotation changes.",
        "5. Use `prediction_instability_summary.csv` and `prediction_instability.csv` to show whether the same genes receive "
        "different ML outputs depending on the database year.",
        "6. Use `feature_rank_overlap.csv` to show whether the biological signals apparently learned by the model are stable.",
        "",
        "## Suggested Manuscript Framing",
        "",
        "The models should be described as probes of database-version sensitivity rather than as final production predictors. "
        "If instability appears across simple and complex model families, the most defensible interpretation is that the "
        "changing resource snapshot materially affects ML conclusions.",
        "",
    ])
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        prog="phengo-version-sensitivity",
        description="Analyse ML sensitivity to model-organism database version/year changes.",
    )
    parser.add_argument("-input_dir", help="Parent directory containing PhenGO output subdirectories or ARFF files")
    parser.add_argument("-arff_files", nargs="+", help="Explicit ARFF files in version/time order")
    parser.add_argument("-dataset_names", nargs="+", help="Names matching -arff_files")
    parser.add_argument("-output_dir", required=True)
    parser.add_argument("-models", nargs="+", default=DEFAULT_MODELS,
                        help=f"Models to run: {', '.join(SKLEARN_MODELS)} or all")
    parser.add_argument("-panels", nargs="+", default=PANELS, choices=PANELS)
    parser.add_argument("-cv_folds", type=int, default=5)
    parser.add_argument("-cv_repeats", type=int, default=3)
    parser.add_argument("-n_estimators", type=int, default=200)
    parser.add_argument("-max_depth", type=int, default=None)
    parser.add_argument("-top_k", type=int, default=20)
    parser.add_argument("-seed", type=int, default=42)
    parser.add_argument("-overwrite", action="store_true")
    args = parser.parse_args()

    if args.input_dir:
        arff_files, dataset_names = discover_arff_files(args.input_dir)
    elif args.arff_files:
        arff_files = args.arff_files
        dataset_names = args.dataset_names or [
            os.path.splitext(os.path.basename(path))[0] for path in arff_files
        ]
    else:
        parser.error("Provide either -input_dir or -arff_files")

    if len(arff_files) < 2:
        parser.error("At least two ARFF datasets are required for version-sensitivity analysis")
    if len(dataset_names) != len(arff_files):
        parser.error("Number of -dataset_names values must match number of -arff_files values")

    models = expand_models(args.models)
    ensure_output_dir(args.output_dir, args.overwrite)
    configure_logger(
        "PhenGO.version_sensitivity",
        enable_file=True,
        log_dir=args.output_dir,
        logfile_name="version_sensitivity.log",
    )

    logger.info("Loading datasets")
    datasets = [
        load_version_dataset(path, name)
        for path, name in zip(arff_files, dataset_names)
    ]
    datasets = sorted(datasets, key=lambda ds: natural_key(ds.name))
    dataset_names = [dataset.name for dataset in datasets]

    all_features = shared_features(datasets)
    all_genes = shared_genes(datasets)
    logger.info(f"Datasets: {', '.join(dataset_names)}")
    logger.info(f"Models: {', '.join(models)}")
    logger.info(f"Panels: {', '.join(args.panels)}")
    logger.info(f"Shared genes across all datasets: {len(all_genes)}")
    logger.info(f"Shared GO features across all datasets: {len(all_features)}")

    write_csv(os.path.join(args.output_dir, "dataset_drift.csv"), dataset_drift_rows(datasets))
    write_csv(os.path.join(args.output_dir, "pairwise_drift.csv"), pairwise_drift_rows(datasets))

    logger.info("Running within-year repeated CV")
    cv_rows, cv_summary = run_within_year_cv(
        datasets, models, args.panels, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "within_year_cv.csv"), cv_rows)
    write_csv(os.path.join(args.output_dir, "within_year_cv_summary.csv"), cv_summary)

    logger.info("Running cross-year transfer analysis")
    transfer_rows = run_cross_year_transfer(
        datasets, models, args.panels, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "cross_year_transfer.csv"), transfer_rows)
    write_transfer_matrices(args.output_dir, transfer_rows, dataset_names)

    logger.info("Running previous-year label baseline")
    write_csv(
        os.path.join(args.output_dir, "previous_year_label_baseline.csv"),
        run_previous_year_label_baseline(datasets),
    )

    logger.info("Running matched-gene prediction instability analysis")
    instability_rows, instability_summary = run_prediction_instability(
        datasets, models, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "prediction_instability.csv"), instability_rows)
    write_csv(os.path.join(args.output_dir, "prediction_instability_summary.csv"), instability_summary)

    logger.info("Running native feature-rank overlap analysis")
    write_csv(
        os.path.join(args.output_dir, "feature_rank_overlap.csv"),
        run_feature_rank_overlap(datasets, models, all_features, all_genes, args),
    )

    write_markdown_report(args.output_dir, args, datasets, models, args.panels)
    logger.info(f"Version-sensitivity analysis complete: {args.output_dir}")


if __name__ == "__main__":
    main()
