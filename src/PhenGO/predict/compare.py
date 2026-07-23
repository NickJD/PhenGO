"""Multi-dataset model comparison for PhenGO-Predict.

This convenience workflow honours every model requested by the user. Paper-level
temporal inference should use ``phengo-version-sensitivity``, which adds repeated
out-of-fold transfer estimates and matched cohorts.
"""
from __future__ import annotations

import json
import logging
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split

from .data import load_arff_data, prepare_data
from .sklearn_models import train_evaluate_sklearn_model

logger = logging.getLogger(__name__)


def _safe_name(value):
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "dataset"


def _load_importance(model_dir):
    path = os.path.join(model_dir, "overall_feature_importance.csv")
    if not os.path.isfile(path):
        return None
    result = pd.read_csv(path)
    required = {"GO_Term", "Overall_Importance"}
    return result if required <= set(result.columns) else None


def _train_nn(options, X_train, y_train, X_test, y_test, genes_test,
              label_encoder, go_terms, output_dir):
    try:
        import tensorflow as tf
        from tensorflow import keras
    except ImportError as exc:
        raise RuntimeError(f"TensorFlow is required for -model nn: {exc}") from exc

    from .evaluate import compute_class_weights, find_optimal_threshold
    from .importance import analyse_feature_importance
    from .model import create_model_sparse_optimised

    tf.keras.utils.set_random_seed(options.seed)
    try:
        tf.config.experimental.enable_op_determinism()
    except Exception:
        logger.warning("TensorFlow deterministic operations are not available")

    model_dir = os.path.join(output_dir, "nn")
    os.makedirs(model_dir, exist_ok=True)
    stratify = y_train if len(np.unique(y_train)) > 1 else None
    X_fit, X_valid, y_fit, y_valid = train_test_split(
        X_train,
        y_train,
        test_size=0.2,
        random_state=options.seed,
        stratify=stratify,
    )
    model = create_model_sparse_optimised(
        input_dim=X_train.shape[1],
        hidden_units=options.hidden_units,
        dropout_rate=options.dropout,
    )
    model.fit(
        X_fit,
        y_fit,
        epochs=options.epochs,
        batch_size=options.batch_size,
        validation_data=(X_valid, y_valid),
        class_weight=compute_class_weights(y_fit),
        callbacks=[keras.callbacks.EarlyStopping(
            monitor="val_loss", patience=15, restore_best_weights=True,
        )],
        verbose=0,
    )
    validation_scores = model.predict(X_valid, verbose=0).reshape(-1)
    threshold, _ = find_optimal_threshold(y_valid, validation_scores)
    scores = model.predict(X_test, verbose=0).reshape(-1)
    predictions = (scores >= threshold).astype(int)

    pd.DataFrame({
        "gene": list(genes_test),
        "true_label": label_encoder.inverse_transform(y_test),
        "predicted_label": label_encoder.inverse_transform(predictions),
        "probability_lethal": scores,
        "correct": predictions == y_test,
    }).to_csv(os.path.join(model_dir, "predictions.csv"), index=False)

    importance_options = type("Options", (), vars(options).copy())()
    importance_options.output_dir = model_dir
    analyse_feature_importance(
        importance_options, model, go_terms, X_test, y_test, options.perm_repeats,
    )
    model.save(os.path.join(model_dir, "gene_essentiality_model.keras"))
    with open(os.path.join(model_dir, "model_schema.json"), "w", encoding="utf-8") as handle:
        json.dump({
            "schema_version": 1,
            "model_type": "nn",
            "feature_names": list(go_terms),
            "class_mapping": {"viable": 0, "lethal": 1},
            "positive_class": "lethal",
            "threshold": float(threshold),
        }, handle, indent=2, sort_keys=True)
        handle.write("\n")

    two_classes = len(np.unique(y_test)) > 1
    return {
        "model_type": "nn",
        "label": "Neural Network",
        "accuracy": accuracy_score(y_test, predictions),
        "precision": precision_score(y_test, predictions, zero_division=0),
        "recall": recall_score(y_test, predictions, zero_division=0),
        "auc": roc_auc_score(y_test, scores) if two_classes else float("nan"),
        "average_precision": (
            average_precision_score(y_test, scores) if two_classes else float("nan")
        ),
        "f1_macro": f1_score(y_test, predictions, average="macro", zero_division=0),
        "f1_lethal": f1_score(y_test, predictions, pos_label=1, zero_division=0),
        "f1_viable": f1_score(y_test, predictions, pos_label=0, zero_division=0),
    }


def compare_datasets(options, arff_files, dataset_names, models_to_run):
    """Train each requested model separately on every supplied dataset."""
    logger.info("Comparing %d datasets with models: %s", len(arff_files), ", ".join(models_to_run))
    metric_rows = []
    importances = {}

    for arff_file, dataset_name in zip(arff_files, dataset_names):
        dataset_dir = os.path.join(options.output_dir, _safe_name(dataset_name))
        os.makedirs(dataset_dir, exist_ok=True)
        df, _ = load_arff_data(arff_file)
        if df is None:
            raise ValueError(f"Could not load {arff_file}")
        X, y, gene_names, label_encoder = prepare_data(df)
        if len(np.unique(y)) < 2 or int(np.bincount(y).min()) < 2:
            raise ValueError(f"{dataset_name}: both classes need at least two genes")
        X_train, X_test, y_train, y_test, _, genes_test = train_test_split(
            X,
            y,
            gene_names,
            test_size=options.test_size,
            random_state=options.seed,
            stratify=y,
        )
        go_terms = list(df.columns[1:-1])

        for model_type in models_to_run:
            logger.info("Training %s on %s", model_type, dataset_name)
            if model_type == "nn":
                summary = _train_nn(
                    options, X_train, y_train, X_test, y_test, genes_test,
                    label_encoder, go_terms, dataset_dir,
                )
            else:
                summary = train_evaluate_sklearn_model(
                    model_type,
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    genes_test,
                    label_encoder,
                    go_terms,
                    dataset_dir,
                    options,
                )
            metric_rows.append({"dataset": dataset_name, **summary})
            importance = _load_importance(os.path.join(dataset_dir, model_type))
            if importance is not None:
                importances[(dataset_name, model_type)] = importance

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(os.path.join(options.output_dir, "comparison_metrics.csv"), index=False)
    _write_metric_plot(options.output_dir, metrics_df)
    _write_importance_comparisons(options.output_dir, dataset_names, models_to_run, importances)
    return metric_rows


def _write_metric_plot(output_dir, metrics_df):
    if metrics_df.empty:
        return
    pivot = metrics_df.pivot(index="dataset", columns="model_type", values="auc")
    ax = pivot.plot(kind="bar", figsize=(max(8, len(pivot) * 1.5), 5), ylim=(0, 1))
    ax.set_ylabel("ROC AUC")
    ax.set_xlabel("Dataset")
    ax.set_title("Held-out model performance by dataset")
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=1)
    ax.legend(title="Model")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "model_performance_comparison.png"), dpi=200)
    plt.close()


def _write_importance_comparisons(output_dir, dataset_names, models, importances, top_k=20):
    overlap_rows = []
    for model in models:
        for name_a in dataset_names:
            for name_b in dataset_names:
                imp_a = importances.get((name_a, model))
                imp_b = importances.get((name_b, model))
                if imp_a is None or imp_b is None:
                    continue
                top_a = set(imp_a.nlargest(top_k, "Overall_Importance")["GO_Term"])
                top_b = set(imp_b.nlargest(top_k, "Overall_Importance")["GO_Term"])
                union = top_a | top_b
                overlap_rows.append({
                    "model": model,
                    "dataset_a": name_a,
                    "dataset_b": name_b,
                    "top_k": top_k,
                    "n_overlap": len(top_a & top_b),
                    "jaccard": len(top_a & top_b) / len(union) if union else float("nan"),
                })

        if len(dataset_names) == 2:
            name_a, name_b = dataset_names
            imp_a = importances.get((name_a, model))
            imp_b = importances.get((name_b, model))
            if imp_a is not None and imp_b is not None:
                merged = imp_a[["GO_Term", "Overall_Importance"]].merge(
                    imp_b[["GO_Term", "Overall_Importance"]],
                    on="GO_Term",
                    how="outer",
                    suffixes=(f"_{_safe_name(name_a)}", f"_{_safe_name(name_b)}"),
                ).fillna(0)
                columns = [column for column in merged if column.startswith("Overall_Importance_")]
                merged["Importance_Difference"] = merged[columns[0]] - merged[columns[1]]
                merged["Absolute_Difference"] = merged["Importance_Difference"].abs()
                merged.sort_values("Absolute_Difference", ascending=False).to_csv(
                    os.path.join(output_dir, f"{model}_differential_importance.csv"),
                    index=False,
                )

    if overlap_rows:
        pd.DataFrame(overlap_rows).to_csv(
            os.path.join(output_dir, "feature_importance_overlap.csv"), index=False,
        )
