"""
predict/compare.py — Multi-dataset comparison and differential GO-term importance
reporting for PhenGO-Predict.
"""
import os
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

from .data       import load_arff_data, prepare_data
from .model      import create_model_sparse_optimised
from .evaluate   import compute_class_weights
from .importance import analyse_feature_importance
from .visualise  import plot_roc_curve

logger = logging.getLogger(__name__)

try:
    from tensorflow import keras
except ImportError as e:
    logger.error(f"TensorFlow import failed: {e}")
    raise


def compare_datasets(options, arff_files, dataset_names):
    """Train separate models on each dataset and compare feature importances.

    Args:
        options      : Namespace with output_dir, epochs, batch_size, test_size,
                       hidden_units, dropout, perm_repeats attributes.
        arff_files   : List of ARFF file paths.
        dataset_names: List of display names (one per ARFF file).

    Returns:
        dict mapping dataset_name -> {model, importance, metrics, predictions, ...}
    """
    logger.info("=" * 80)
    logger.info("MULTI-DATASET COMPARISON MODE")
    logger.info(f"Comparing {len(arff_files)} datasets")
    logger.info("=" * 80)

    all_results = {}

    for arff_file, dataset_name in zip(arff_files, dataset_names):
        logger.info(f"Processing dataset: {dataset_name}")
        dataset_dir = os.path.join(options.output_dir, dataset_name)
        os.makedirs(dataset_dir, exist_ok=True)

        df, meta = load_arff_data(arff_file)
        X, y, gene_names, label_encoder = prepare_data(df)

        X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
            X, y, gene_names, test_size=options.test_size, random_state=42, stratify=y,
        )

        model = create_model_sparse_optimised(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout,
        )
        ds_class_weights = compute_class_weights(y_train)

        history = model.fit(
            X_train, y_train,
            epochs=options.epochs,
            batch_size=options.batch_size,
            validation_split=0.2,
            class_weight=ds_class_weights,
            callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=15,
                                                     restore_best_weights=True)],
            verbose=0,
        )

        y_pred_proba = model.predict(X_test, verbose=0)
        test_results = model.evaluate(X_test, y_test, verbose=0)

        dataset_options = argparse.Namespace(**vars(options))
        dataset_options.output_dir = dataset_dir
        go_term_columns = df.columns[1:-1]
        overall_imp, lethal_imp, viable_imp = analyse_feature_importance(
            dataset_options, model, go_term_columns, X_test, y_test, options.perm_repeats,
        )

        all_results[dataset_name] = {
            "model":      model,
            "importance": {"overall": overall_imp, "lethal": lethal_imp, "viable": viable_imp},
            "metrics": {
                "loss":      test_results[0],
                "accuracy":  test_results[1],
                "precision": test_results[2] if len(test_results) > 2 else 0,
                "recall":    test_results[3] if len(test_results) > 3 else 0,
                "auc":       roc_auc_score(y_test, y_pred_proba),
            },
            "predictions": y_pred_proba,
            "true_labels": y_test,
            "go_terms":    list(go_term_columns),
        }

        plot_roc_curve(options, y_test, y_pred_proba, dataset_name)

    generate_comparison_plots(options, all_results, dataset_names)
    generate_differential_importance_report(options, all_results, dataset_names)
    return all_results


def generate_comparison_plots(options, all_results, dataset_names):
    """Metric bar-charts and GO-term importance heatmap across datasets."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for idx, metric in enumerate(["accuracy", "precision", "recall"]):
        values = [all_results[n]["metrics"][metric] for n in dataset_names]
        axes[idx].bar(range(len(dataset_names)), values, alpha=0.7)
        axes[idx].set_xticks(range(len(dataset_names)))
        axes[idx].set_xticklabels(dataset_names, rotation=45, ha="right")
        axes[idx].set_ylabel(metric.capitalize())
        axes[idx].set_title(f"{metric.capitalize()} Comparison")
        axes[idx].set_ylim([0, 1])
        axes[idx].grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(options.output_dir, "model_performance_comparison.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    all_top_terms = set()
    for name in dataset_names:
        top_30 = all_results[name]["importance"]["overall"].head(30)["GO_Term"]
        all_top_terms.update(top_30)
    term_list = sorted(all_top_terms)

    importance_matrix = []
    for term in term_list:
        row = []
        for name in dataset_names:
            imp_df   = all_results[name]["importance"]["overall"]
            term_imp = imp_df[imp_df["GO_Term"] == term]["Overall_Importance"]
            row.append(term_imp.values[0] if len(term_imp) > 0 else 0)
        importance_matrix.append(row)

    fig, ax = plt.subplots(figsize=(12, max(8, len(term_list) * 0.3)))
    im = ax.imshow(importance_matrix, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(dataset_names)))
    ax.set_xticklabels(dataset_names, rotation=45, ha="right")
    ax.set_yticks(range(len(term_list)))
    ax.set_yticklabels(
        [t[:50] + "..." if len(t) > 50 else t for t in term_list], fontsize=8,
    )
    ax.set_title("GO Term Importance Across Datasets", fontsize=14, pad=20)
    plt.colorbar(im, ax=ax, label="Importance Score")
    plt.tight_layout()
    plt.savefig(os.path.join(options.output_dir, "go_term_importance_heatmap.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Comparison plots saved to {options.output_dir}")


def generate_differential_importance_report(options, all_results, dataset_names):
    """Report GO terms with differential importance across exactly two datasets."""
    if len(dataset_names) != 2:
        logger.error("Differential analysis currently only supports 2 datasets")
        return

    name1, name2 = dataset_names
    imp1 = all_results[name1]["importance"]["overall"]
    imp2 = all_results[name2]["importance"]["overall"]

    merged = pd.merge(
        imp1[["GO_Term","Overall_Importance","Overall_Std"]],
        imp2[["GO_Term","Overall_Importance","Overall_Std"]],
        on="GO_Term", suffixes=(f"_{name1}", f"_{name2}"),
    )
    merged["Importance_Diff"] = (merged[f"Overall_Importance_{name1}"] -
                                 merged[f"Overall_Importance_{name2}"])
    merged["Abs_Diff"] = abs(merged["Importance_Diff"])
    merged_sorted = merged.sort_values("Abs_Diff", ascending=False)

    report_path = os.path.join(options.output_dir, "differential_importance_report.csv")
    merged_sorted.to_csv(report_path, index=False)

    logger.info("=" * 80)
    logger.info("TOP 20 DIFFERENTIALLY IMPORTANT GO TERMS")
    logger.info("=" * 80)
    logger.info(f"Positive: more important in {name1} | Negative: more important in {name2}")
    logger.info("\n" + merged_sorted[[
        "GO_Term", "Importance_Diff",
        f"Overall_Importance_{name1}", f"Overall_Importance_{name2}",
    ]].head(20).to_string(index=False))
    logger.info(f"Full report saved to: {report_path}")
