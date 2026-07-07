"""
predict/visualise.py — Plotting utilities for PhenGO-Predict.
"""
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

logger = logging.getLogger(__name__)


def plot_roc_curve(options, y_test, y_pred_proba, dataset_name=""):
    """Plot ROC curve with AUC score and save to output_dir."""
    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    auc_score   = roc_auc_score(y_test, y_pred_proba)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, linewidth=2, label=f"ROC curve (AUC = {auc_score:.3f})")
    plt.plot([0, 1], [0, 1], "k--", linewidth=1, label="Random classifier")
    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.title(f"ROC Curve{' - ' + dataset_name if dataset_name else ''}", fontsize=14)
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    filename = f"roc_curve{'_' + dataset_name if dataset_name else ''}.png"
    output_path = os.path.join(options.output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"ROC curve saved to {output_path}")


def plot_training_history(options, history):
    """Plot training / validation curves (accuracy, loss, precision, recall)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    axes[0, 0].plot(history.history["accuracy"],     label="Training")
    axes[0, 0].plot(history.history["val_accuracy"], label="Validation")
    axes[0, 0].set_title("Model Accuracy")
    axes[0, 0].set_xlabel("Epoch")
    axes[0, 0].set_ylabel("Accuracy")
    axes[0, 0].legend()

    axes[0, 1].plot(history.history["loss"],     label="Training")
    axes[0, 1].plot(history.history["val_loss"], label="Validation")
    axes[0, 1].set_title("Model Loss")
    axes[0, 1].set_xlabel("Epoch")
    axes[0, 1].set_ylabel("Loss")
    axes[0, 1].legend()

    if "precision" in history.history:
        axes[1, 0].plot(history.history["precision"],     label="Training")
        axes[1, 0].plot(history.history["val_precision"], label="Validation")
        axes[1, 0].set_title("Model Precision")
        axes[1, 0].set_xlabel("Epoch")
        axes[1, 0].set_ylabel("Precision")
        axes[1, 0].legend()

    if "recall" in history.history:
        axes[1, 1].plot(history.history["recall"],     label="Training")
        axes[1, 1].plot(history.history["val_recall"], label="Validation")
        axes[1, 1].set_title("Model Recall")
        axes[1, 1].set_xlabel("Epoch")
        axes[1, 1].set_ylabel("Recall")
        axes[1, 1].legend()

    plt.tight_layout()
    output_path = os.path.join(options.output_dir, "training_history.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Training history saved to {output_path}")


def plot_feature_importance_with_errors(options, importance_df, top_n=30,
                                        title="Feature Importance"):
    """Three-panel bar chart: overall / lethal / viable importance with error bars."""
    top_features = importance_df.nlargest(top_n, "Overall_Importance")

    fig, axes = plt.subplots(1, 3, figsize=(20, 8))

    def _bar_panel(ax, df_sorted, col_imp, col_std, colour, panel_title):
        y_pos = np.arange(len(df_sorted))
        ax.barh(y_pos, df_sorted[col_imp], xerr=df_sorted[col_std],
                capsize=3, alpha=0.7, color=colour)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [t[:40] + "..." if len(t) > 40 else t for t in df_sorted["GO_Term"]],
            fontsize=8,
        )
        ax.set_xlabel("Importance (Accuracy Drop)", fontsize=10)
        ax.set_title(panel_title, fontsize=12)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)

    top_lethal = importance_df.nlargest(top_n, "Lethal_Importance")
    top_viable = importance_df.nlargest(top_n, "Viable_Importance")

    _bar_panel(axes[0], top_features, "Overall_Importance", "Overall_Std",
               "steelblue", f"Top {top_n} - Overall")
    _bar_panel(axes[1], top_lethal,   "Lethal_Importance",  "Lethal_Std",
               "red",       f"Top {top_n} - Lethal")
    _bar_panel(axes[2], top_viable,   "Viable_Importance",  "Viable_Std",
               "green",     f"Top {top_n} - Viable")

    plt.tight_layout()
    output_path = os.path.join(options.output_dir, "feature_importance_with_errors.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Feature importance plot saved to {output_path}")
