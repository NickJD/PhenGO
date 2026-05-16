#!/usr/bin/env python3
"""
PhenGO-Predict Report Generator
================================
Produces paper-ready figures and a manuscript draft from PhenGO-Predict outputs.

Two modes
---------
Single dataset::

    predict_report -predict_dir /path/to/Predict \
                   -output_dir  ./report \
                   -dataset_name "C. elegans 2025"

Multi-dataset (compare across species / years)::

    predict_report -predict_dirs /data/Worm/Predict /data/Fly/Predict /data/Yeast/Predict \
                   -dataset_names Worm Fly Yeast \
                   -output_dir ./cross_species_report

    # OR auto-discover: each sub-dir that contains a Predict folder
    predict_report -input_dir /data/PhenGO/Paper_Data \
                   -output_dir ./cross_all_report
"""
import os as _os, sys as _sys
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here   = _os.path.dirname(_os.path.abspath(__file__))
    _src    = _os.path.normpath(_os.path.join(_here, '..', '..'))
    _phengo = _os.path.dirname(_here)
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',         _os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.scripts', _os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in _sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            _sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.scripts'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo

import re
import argparse
import logging
import warnings
import os
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

from sklearn.metrics import roc_curve, auc as sk_auc, f1_score

from ..constants import configure_logger, PhenGO_VERSION

logger = logging.getLogger(__name__)

# ── Publication-quality matplotlib defaults ───────────────────────────────────
plt.rcParams.update({
    'font.family':        'serif',
    'font.size':          11,
    'axes.titlesize':     12,
    'axes.titleweight':   'bold',
    'axes.labelsize':     11,
    'xtick.labelsize':    9,
    'ytick.labelsize':    9,
    'legend.fontsize':    9,
    'legend.framealpha':  0.85,
    'figure.dpi':         150,   # screen preview; saved at 300
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches': 0.05,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
})

# Colorblind-safe palette (Wong 2011)
PALETTE = ['#0072B2', '#E69F00', '#009E73', '#CC79A7',
           '#56B4E9', '#D55E00', '#F0E442', '#000000']

MODEL_LABELS = {
    'nn':  'Neural\nNetwork',
    'rf':  'Random\nForest',
    'gb':  'Gradient\nBoosting',
    'dt':  'Decision\nTree',
    'lr':  'Logistic\nRegression',
    'svm': 'SVM',
}
MODEL_LABELS_SHORT = {k: v.replace('\n', ' ') for k, v in MODEL_LABELS.items()}

METRIC_LABELS = {
    'AUC':       'AUC-ROC',
    'Accuracy':  'Accuracy',
    'Precision': 'Precision (viable)',
    'Recall':    'Recall (viable)',
    'F1_macro':  'F1 (macro)',
    'F1_lethal': 'F1 (lethal)',
    'F1_viable': 'F1 (viable)',
}


# ── Parsing helpers ───────────────────────────────────────────────────────────

def parse_report(path: str) -> dict:
    """Extract metrics dict from a model report.txt."""
    metrics = {}
    if not os.path.exists(path):
        return metrics
    with open(path) as f:
        text = f.read()

    for key in ('Loss', 'Accuracy', 'Precision', 'Recall', 'AUC'):
        m = re.search(rf'{key}:\s*([\d.]+)', text)
        if m:
            metrics[key] = float(m.group(1))

    # Per-class F1 from classification report
    for cls in ('lethal', 'viable'):
        m = re.search(rf'\s+{cls}\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', text)
        if m:
            metrics[f'F1_{cls}'] = float(m.group(3))

    # Macro F1
    m = re.search(r'macro avg\s+[\d.]+\s+[\d.]+\s+([\d.]+)', text)
    if m:
        metrics['F1_macro'] = float(m.group(1))

    # Confusion matrix
    m = re.search(r'\[\[(\d+)\s+(\d+)\]\s*\[(\d+)\s+(\d+)\]\]', text)
    if m:
        metrics['cm'] = np.array([[int(m.group(1)), int(m.group(2))],
                                   [int(m.group(3)), int(m.group(4))]])
    return metrics


def parse_predictions(path: str) -> pd.DataFrame | None:
    """Load predictions.csv; return None if missing."""
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
        # Normalise column names
        df.columns = [c.lower() for c in df.columns]
        return df
    except Exception:
        return None


def parse_feature_importance(model_dir: str, model: str) -> pd.DataFrame | None:
    """Load feature importance file for a model; return top rows."""
    candidates = ['feature_importance.csv', 'overall_feature_importance.csv']
    for fname in candidates:
        path = os.path.join(model_dir, fname)
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                # Normalise: ensure columns GO_Term and Importance exist
                if 'GO_Term' not in df.columns:
                    # try first two columns
                    df.columns = ['GO_Term', 'Importance'] + list(df.columns[2:])
                if 'Overall_Importance' in df.columns:
                    df = df.rename(columns={'Overall_Importance': 'Importance'})
                df = df[['GO_Term', 'Importance']].dropna()
                df['Importance'] = pd.to_numeric(df['Importance'], errors='coerce')
                df = df.dropna().sort_values('Importance', ascending=False).reset_index(drop=True)
                return df
            except Exception:
                pass
    return None


def collect_model_data(predict_dir: str) -> dict:
    """Collect all model outputs from a Predict directory."""
    data = {}
    if not os.path.isdir(predict_dir):
        return data

    model_dirs = [d for d in os.listdir(predict_dir)
                  if os.path.isdir(os.path.join(predict_dir, d))
                  and d in MODEL_LABELS]

    for model in model_dirs:
        mdir = os.path.join(predict_dir, model)
        entry = {
            'metrics':          parse_report(os.path.join(mdir, 'report.txt')),
            'predictions':      parse_predictions(os.path.join(mdir, 'predictions.csv')),
            'feature_importance': parse_feature_importance(mdir, model),
            'roc_png':          os.path.join(mdir, 'roc_curve.png'),
            'fi_png':           os.path.join(mdir, 'feature_importance.png'),
            'training_png':     os.path.join(mdir, 'training_history.png'),
        }
        if entry['metrics']:
            data[model] = entry

    return data


# ── Figure helpers ────────────────────────────────────────────────────────────

def _savefig(fig, output_dir: str, stem: str):
    """Save figure as both PNG and PDF for paper use."""
    for ext in ('png', 'pdf'):
        path = os.path.join(output_dir, f'{stem}.{ext}')
        fig.savefig(path)
    plt.close(fig)
    logger.info(f"  Saved {stem}.png / .pdf")


def fig_metrics_comparison(model_data: dict, output_dir: str, dataset_name: str):
    """Grouped bar chart of all key metrics across models."""
    metrics_to_plot = ['AUC', 'Accuracy', 'F1_macro', 'F1_lethal', 'F1_viable']
    models = sorted(model_data.keys(), key=lambda m: -model_data[m]['metrics'].get('AUC', 0))

    metric_vals = {met: [] for met in metrics_to_plot}
    model_labels = []
    for m in models:
        mets = model_data[m]['metrics']
        model_labels.append(MODEL_LABELS_SHORT.get(m, m))
        for met in metrics_to_plot:
            metric_vals[met].append(mets.get(met, np.nan))

    n_models  = len(models)
    n_metrics = len(metrics_to_plot)
    x         = np.arange(n_models)
    bar_w     = 0.14
    offsets   = np.linspace(-(n_metrics - 1) / 2, (n_metrics - 1) / 2, n_metrics) * bar_w

    fig, ax = plt.subplots(figsize=(max(7, n_models * 1.5), 5))

    for i, (met, offset) in enumerate(zip(metrics_to_plot, offsets)):
        vals = metric_vals[met]
        bars = ax.bar(x + offset, vals, bar_w * 0.9,
                      label=METRIC_LABELS.get(met, met),
                      color=PALETTE[i % len(PALETTE)],
                      edgecolor='white', linewidth=0.5)
        for bar, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.005,
                        f'{v:.3f}', ha='center', va='bottom',
                        fontsize=6.5, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(model_labels, fontsize=9)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel('Score')
    ax.set_title(f'Model performance comparison — {dataset_name}')
    ax.legend(loc='lower right', ncol=2)
    ax.axhline(0.5, color='grey', linewidth=0.7, linestyle='--', alpha=0.5)

    _savefig(fig, output_dir, 'fig1_metrics_comparison')


def fig_roc_curves(model_data: dict, output_dir: str, dataset_name: str):
    """Combined ROC curves recomputed from predictions.csv."""
    fig, ax = plt.subplots(figsize=(5.5, 5))

    plotted = 0
    for i, (model, entry) in enumerate(sorted(model_data.items())):
        preds = entry['predictions']
        if preds is None:
            continue
        if 'probability_lethal' not in preds.columns:
            continue
        y_true = (preds['true_label'].str.lower() == 'lethal').astype(int)
        y_score = preds['probability_lethal']
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = sk_auc(fpr, tpr)
        ax.plot(fpr, tpr,
                color=PALETTE[i % len(PALETTE)],
                linewidth=1.8,
                label=f'{MODEL_LABELS_SHORT.get(model, model)} (AUC = {roc_auc:.3f})')
        plotted += 1

    ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.5, label='Random')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC Curves — {dataset_name}')
    ax.legend(loc='lower right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)

    if plotted > 0:
        _savefig(fig, output_dir, 'fig2_roc_curves')
    else:
        plt.close(fig)


def fig_confusion_matrices(model_data: dict, output_dir: str, dataset_name: str):
    """Tiled normalised confusion matrices, one per model."""
    models = sorted(model_data.keys(), key=lambda m: -model_data[m]['metrics'].get('AUC', 0))
    n = len(models)
    cols = min(3, n)
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols,
                              figsize=(cols * 3.2, rows * 3.0),
                              squeeze=False)

    classes = ['lethal', 'viable']

    for idx, model in enumerate(models):
        ax = axes[idx // cols][idx % cols]
        cm = model_data[model]['metrics'].get('cm')
        if cm is None:
            ax.axis('off')
            continue

        cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
        im = ax.imshow(cm_norm, interpolation='nearest',
                       cmap='Blues', vmin=0, vmax=1)
        ax.set_title(MODEL_LABELS_SHORT.get(model, model), fontsize=10, fontweight='bold')
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(classes, fontsize=8)
        ax.set_yticklabels(classes, fontsize=8)
        ax.set_xlabel('Predicted', fontsize=8)
        ax.set_ylabel('True', fontsize=8)

        thresh = cm_norm.max() / 2
        for i in range(2):
            for j in range(2):
                raw = cm[i, j]
                norm_v = cm_norm[i, j]
                ax.text(j, i, f'{norm_v:.2f}\n({raw})',
                        ha='center', va='center', fontsize=8,
                        color='white' if norm_v > thresh else 'black')

    # Hide unused axes
    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis('off')

    fig.suptitle(f'Confusion Matrices — {dataset_name}', fontsize=12, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _savefig(fig, output_dir, 'fig3_confusion_matrices')


def fig_feature_importance_heatmap(model_data: dict, output_dir: str,
                                   dataset_name: str, top_n: int = 25):
    """Heatmap of top-N GO terms × models with normalised importance."""
    # Collect top_n terms from each model
    candidate_terms = set()
    fi_dfs = {}
    for model, entry in model_data.items():
        fi = entry['feature_importance']
        if fi is not None and not fi.empty:
            top = fi.head(top_n)['GO_Term'].tolist()
            candidate_terms.update(top)
            fi_dfs[model] = fi.set_index('GO_Term')['Importance']

    if not fi_dfs or not candidate_terms:
        logger.warning("  No feature importance data available for heatmap.")
        return

    terms = sorted(candidate_terms)
    models = sorted(fi_dfs.keys())

    mat = pd.DataFrame(index=terms, columns=models, dtype=float)
    for model in models:
        for term in terms:
            mat.loc[term, model] = fi_dfs[model].get(term, 0.0)

    # Normalise each model column to [0, 1] for heatmap colour
    mat_norm = mat.copy()
    for col in mat_norm.columns:
        col_max = mat_norm[col].max()
        if col_max > 0:
            mat_norm[col] = mat_norm[col] / col_max

    # Sort rows by mean normalised importance
    mat_norm['_mean'] = mat_norm[models].mean(axis=1)
    mat_norm = mat_norm.sort_values('_mean', ascending=False).drop(columns='_mean')
    mat = mat.loc[mat_norm.index]

    n_terms  = len(mat_norm)
    n_models = len(models)
    fig_h = max(6, n_terms * 0.32)
    fig_w = max(4, n_models * 1.2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if HAS_SEABORN:
        sns.heatmap(mat_norm.astype(float), ax=ax,
                    cmap='YlOrRd', linewidths=0.3, linecolor='white',
                    xticklabels=[MODEL_LABELS_SHORT.get(m, m) for m in models],
                    yticklabels=mat_norm.index,
                    cbar_kws={'label': 'Normalised importance', 'shrink': 0.6})
    else:
        im = ax.imshow(mat_norm.astype(float).values, aspect='auto', cmap='YlOrRd',
                       vmin=0, vmax=1)
        ax.set_xticks(range(n_models))
        ax.set_xticklabels([MODEL_LABELS_SHORT.get(m, m) for m in models])
        ax.set_yticks(range(n_terms))
        ax.set_yticklabels(mat_norm.index, fontsize=7)
        plt.colorbar(im, ax=ax, label='Normalised importance', shrink=0.6)

    ax.set_title(f'Feature Importance Heatmap (top {top_n} GO terms) — {dataset_name}',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('Model')
    plt.tight_layout()
    _savefig(fig, output_dir, 'fig4_feature_importance_heatmap')


def fig_prediction_agreement(model_data: dict, output_dir: str, dataset_name: str):
    """Pairwise model prediction agreement heatmap + per-gene cross-model consensus."""
    models = sorted(m for m, e in model_data.items() if e['predictions'] is not None)
    if len(models) < 2:
        return

    # Build gene → model prediction matrix
    gene_preds = {}
    for model in models:
        preds = model_data[model]['predictions']
        for _, row in preds.iterrows():
            gene = row.get('gene', row.iloc[0])
            pred = str(row.get('predicted_label', '')).lower()
            gene_preds.setdefault(gene, {})[model] = pred

    genes = sorted(gene_preds.keys())
    pred_mat = pd.DataFrame({m: [gene_preds[g].get(m, np.nan) for g in genes]
                              for m in models}, index=genes)

    # Binary: 1 = lethal predicted
    bin_mat = (pred_mat == 'lethal').astype(float)
    bin_mat[pred_mat.isna()] = np.nan

    # Pairwise agreement
    agreement = pd.DataFrame(index=models, columns=models, dtype=float)
    for m1 in models:
        for m2 in models:
            if m1 == m2:
                agreement.at[m1, m2] = 1.0
                continue
            s1 = bin_mat[m1].dropna()
            s2 = bin_mat[m2].dropna()
            common_idx = s1.index.intersection(s2.index)
            if len(common_idx) > 0:
                agreement.at[m1, m2] = float((s1.loc[common_idx] == s2.loc[common_idx]).mean())
            else:
                agreement.at[m1, m2] = np.nan

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ── Left: agreement heatmap ──
    ax = axes[0]
    agr_arr = agreement.astype(float).values
    if HAS_SEABORN:
        sns.heatmap(agreement.astype(float), ax=ax, annot=True, fmt='.3f',
                    cmap='RdYlGn', vmin=0.5, vmax=1.0,
                    xticklabels=[MODEL_LABELS_SHORT.get(m, m) for m in models],
                    yticklabels=[MODEL_LABELS_SHORT.get(m, m) for m in models],
                    cbar_kws={'label': 'Agreement', 'shrink': 0.7})
    else:
        im = ax.imshow(agr_arr, cmap='RdYlGn', vmin=0.5, vmax=1.0)
        tick_labels = [MODEL_LABELS_SHORT.get(m, m) for m in models]
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        ax.set_yticks(range(len(models)))
        ax.set_yticklabels(tick_labels)
        for i in range(len(models)):
            for j in range(len(models)):
                v = agr_arr[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f'{v:.3f}', ha='center', va='center', fontsize=8)
        plt.colorbar(im, ax=ax, label='Agreement', shrink=0.7)
    ax.set_title('Pairwise prediction agreement', fontsize=10, fontweight='bold')

    # ── Right: per-gene consensus bar by true label ──
    ax2 = axes[1]
    # Count how many models predicted lethal for each gene
    n_models = len(models)
    bin_mat['n_lethal'] = bin_mat[models].sum(axis=1)
    bin_mat['true_label'] = [gene_preds[g].get(models[0], '') for g in genes]
    # Get true label from first model predictions
    first_model_preds = model_data[models[0]]['predictions']
    true_map = {}
    if first_model_preds is not None:
        for _, row in first_model_preds.iterrows():
            gene = row.get('gene', row.iloc[0])
            true_map[gene] = str(row.get('true_label', '')).lower()

    bins = range(n_models + 2)
    lethal_counts = []
    viable_counts = []
    for k in range(n_models + 1):
        genes_at_k = [g for g in genes if bin_mat.loc[g, 'n_lethal'] == k]
        lethal_counts.append(sum(1 for g in genes_at_k if true_map.get(g) == 'lethal'))
        viable_counts.append(sum(1 for g in genes_at_k if true_map.get(g) == 'viable'))

    x = np.arange(n_models + 1)
    ax2.bar(x, lethal_counts, color=PALETTE[0], alpha=0.8, label='True lethal')
    ax2.bar(x, viable_counts, bottom=lethal_counts, color=PALETTE[1], alpha=0.8, label='True viable')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{k}/{n_models}' for k in range(n_models + 1)], fontsize=8)
    ax2.set_xlabel(f'Number of models predicting lethal (out of {n_models})')
    ax2.set_ylabel('Gene count')
    ax2.set_title('Cross-model consensus', fontsize=10, fontweight='bold')
    ax2.legend()

    fig.suptitle(f'Model Prediction Agreement — {dataset_name}',
                 fontsize=12, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _savefig(fig, output_dir, 'fig5_prediction_agreement')


# ── Cross-dataset figures ─────────────────────────────────────────────────────

def fig_cross_dataset_heatmap(all_data: dict, output_dir: str, metric: str = 'AUC'):
    """Heatmap: datasets (rows) × models (columns) for a given metric."""
    datasets = sorted(all_data.keys())
    models_seen = set()
    for d in datasets:
        models_seen.update(all_data[d].keys())
    models = sorted(models_seen, key=lambda m: MODEL_LABELS_SHORT.get(m, m))

    mat = pd.DataFrame(index=datasets, columns=models, dtype=float)
    for ds in datasets:
        for m in models:
            mat.loc[ds, m] = all_data[ds].get(m, {}).get('metrics', {}).get(metric, np.nan)

    fig, ax = plt.subplots(figsize=(max(4, len(models) * 1.4),
                                    max(3, len(datasets) * 0.8 + 1)))
    if HAS_SEABORN:
        sns.heatmap(mat.astype(float), ax=ax, annot=True, fmt='.3f',
                    cmap='YlGnBu', vmin=0.5, vmax=1.0,
                    xticklabels=[MODEL_LABELS_SHORT.get(m, m) for m in models],
                    yticklabels=datasets,
                    cbar_kws={'label': METRIC_LABELS.get(metric, metric), 'shrink': 0.7})
    else:
        arr = mat.astype(float).values
        im = ax.imshow(arr, cmap='YlGnBu', vmin=0.5, vmax=1.0, aspect='auto')
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels([MODEL_LABELS_SHORT.get(m, m) for m in models], rotation=45, ha='right')
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets)
        for i in range(len(datasets)):
            for j in range(len(models)):
                v = arr[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f'{v:.3f}', ha='center', va='center', fontsize=8,
                            color='white' if v > 0.8 else 'black')
        plt.colorbar(im, ax=ax, label=METRIC_LABELS.get(metric, metric), shrink=0.7)

    ax.set_title(f'{METRIC_LABELS.get(metric, metric)} across datasets and models',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    _savefig(fig, output_dir, f'cross_fig1_{metric.lower()}_heatmap')


def fig_cross_dataset_lines(all_data: dict, output_dir: str):
    """Line plots: one per metric, models as separate lines, x = datasets."""
    datasets = sorted(all_data.keys())
    models_seen = set()
    for d in datasets:
        models_seen.update(all_data[d].keys())
    models = sorted(models_seen)

    metrics_to_plot = ['AUC', 'Accuracy', 'F1_macro']
    fig, axes = plt.subplots(1, len(metrics_to_plot),
                             figsize=(5 * len(metrics_to_plot), 4.5), sharey=False)
    if len(metrics_to_plot) == 1:
        axes = [axes]

    x_pos = np.arange(len(datasets))

    for ax, metric in zip(axes, metrics_to_plot):
        for i, model in enumerate(models):
            vals = [all_data[ds].get(model, {}).get('metrics', {}).get(metric, np.nan)
                    for ds in datasets]
            ax.plot(x_pos, vals,
                    marker='o', linewidth=1.8, markersize=5,
                    color=PALETTE[i % len(PALETTE)],
                    label=MODEL_LABELS_SHORT.get(model, model))

        ax.set_xticks(x_pos)
        ax.set_xticklabels(datasets, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(METRIC_LABELS.get(metric, metric))
        ax.set_title(METRIC_LABELS.get(metric, metric), fontweight='bold')
        ax.set_ylim(0.4, 1.02)
        ax.axhline(0.5, color='grey', linewidth=0.7, linestyle='--', alpha=0.4)

    axes[-1].legend(loc='lower right', fontsize=8)
    fig.suptitle('Model performance across datasets', fontsize=12, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    _savefig(fig, output_dir, 'cross_fig2_metric_trends')


def fig_cross_feature_overlap(all_data: dict, output_dir: str, top_n: int = 20):
    """Bar chart: top features shared by multiple datasets and models."""
    # Count how many (dataset, model) pairs each GO term appears in top_n
    term_counts = defaultdict(int)
    total_pairs = 0
    for ds, model_data_d in all_data.items():
        for model, entry in model_data_d.items():
            fi = entry.get('feature_importance')
            if fi is not None and not fi.empty:
                for term in fi.head(top_n)['GO_Term']:
                    term_counts[term] += 1
                total_pairs += 1

    if not term_counts or total_pairs == 0:
        return

    top_terms = sorted(term_counts.items(), key=lambda x: -x[1])[:top_n]
    terms, counts = zip(*top_terms)
    fracs = [c / total_pairs for c in counts]

    fig, ax = plt.subplots(figsize=(7, max(4, len(terms) * 0.35)))
    y = np.arange(len(terms))
    bars = ax.barh(y, fracs, color=PALETTE[0], alpha=0.8, edgecolor='white')
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f'Fraction of (dataset × model) pairs in top {top_n}')
    ax.set_title(f'Most consistently important GO terms across datasets',
                 fontsize=11, fontweight='bold')
    ax.set_xlim(0, 1)
    for bar, frac in zip(bars, fracs):
        ax.text(frac + 0.01, bar.get_y() + bar.get_height() / 2,
                f'{frac:.2f}', va='center', fontsize=7)
    plt.tight_layout()
    _savefig(fig, output_dir, 'cross_fig3_feature_overlap')


def fig_best_model_summary(all_data: dict, output_dir: str):
    """Stacked bar chart: best model by AUC per dataset, coloured by model."""
    datasets = sorted(all_data.keys())
    model_color = {m: PALETTE[i % len(PALETTE)]
                   for i, m in enumerate(sorted(MODEL_LABELS))}

    fig, ax = plt.subplots(figsize=(max(5, len(datasets) * 1.4), 4.5))

    aucs  = []
    colors = []
    best_models = []
    for ds in datasets:
        best_m, best_a = None, 0.0
        for m, entry in all_data[ds].items():
            a = entry.get('metrics', {}).get('AUC', 0.0)
            if a > best_a:
                best_a, best_m = a, m
        aucs.append(best_a)
        colors.append(model_color.get(best_m, '#888888'))
        best_models.append(best_m)

    bars = ax.bar(range(len(datasets)), aucs, color=colors, edgecolor='white')
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels(datasets, rotation=45, ha='right', fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Best AUC-ROC')
    ax.set_title('Best-performing model per dataset', fontsize=11, fontweight='bold')
    ax.axhline(0.5, color='grey', linewidth=0.7, linestyle='--', alpha=0.4)

    for bar, auc_v, model in zip(bars, aucs, best_models):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f'{auc_v:.3f}\n({MODEL_LABELS_SHORT.get(model, model)})',
                ha='center', va='bottom', fontsize=7.5)

    # Legend
    legend_handles = [Patch(color=model_color[m], label=MODEL_LABELS_SHORT[m])
                      for m in sorted(MODEL_LABELS) if m in model_color]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8, ncol=2)

    plt.tight_layout()
    _savefig(fig, output_dir, 'cross_fig4_best_model_summary')


# ── Text report helpers ───────────────────────────────────────────────────────

def _fmt(v, decimals=3):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return 'N/A'
    return f'{v:.{decimals}f}'


def write_metrics_csv(model_data: dict, output_dir: str, dataset_name: str):
    """Write a flat CSV of all metrics for all models."""
    rows = []
    metrics = ['AUC', 'Accuracy', 'Precision', 'Recall', 'F1_macro', 'F1_lethal', 'F1_viable']
    for model in sorted(model_data.keys()):
        row = {'Dataset': dataset_name, 'Model': MODEL_LABELS_SHORT.get(model, model)}
        for met in metrics:
            row[met] = _fmt(model_data[model]['metrics'].get(met), 4)
        rows.append(row)
    df = pd.DataFrame(rows)
    path = os.path.join(output_dir, 'metrics_summary.csv')
    df.to_csv(path, index=False)
    logger.info(f"  Saved metrics_summary.csv")
    return df


def write_results_draft(model_data: dict, output_dir: str, dataset_name: str):
    """Generate a manuscript-ready Results section paragraph."""
    models_by_auc = sorted(
        [(m, e['metrics'].get('AUC', 0)) for m, e in model_data.items()],
        key=lambda x: -x[1]
    )
    best_m, best_auc = models_by_auc[0]
    lines = []
    lines.append(f"### Results draft for: {dataset_name}\n")
    lines.append(
        f"Five machine learning classifiers were trained and evaluated on the PhenGO "
        f"gene essentiality ARFF dataset for {dataset_name} using an 80:20 stratified "
        f"train/test split.  "
        f"Performance was assessed using the area under the receiver operating "
        f"characteristic curve (AUC-ROC), accuracy, precision, recall, and the "
        f"macro-averaged F1 score.\n"
    )

    model_sentences = []
    for model, auc in models_by_auc:
        mets = model_data[model]['metrics']
        acc  = mets.get('Accuracy', float('nan'))
        f1   = mets.get('F1_macro',  float('nan'))
        fl   = mets.get('F1_lethal', float('nan'))
        mname = MODEL_LABELS_SHORT.get(model, model)
        model_sentences.append(
            f"The {mname} achieved an AUC of {_fmt(auc)}, "
            f"accuracy of {_fmt(acc)}, macro-F1 of {_fmt(f1)}, "
            f"and a lethal-class F1 of {_fmt(fl)}."
        )
    lines.append('  '.join(model_sentences) + '\n')

    bname = MODEL_LABELS_SHORT.get(best_m, best_m)
    lines.append(
        f"The {bname} achieved the highest AUC ({_fmt(best_auc)}), "
        f"demonstrating the strongest overall discrimination between lethal and "
        f"viable gene disruptions.  "
    )

    # Feature importance note
    fi = model_data[best_m].get('feature_importance')
    if fi is not None and not fi.empty:
        top3 = fi.head(3)['GO_Term'].tolist()
        lines.append(
            f"Feature importance analysis of the {bname} identified "
            f"{', '.join(top3)} among the most discriminative GO terms, "
            f"suggesting these biological processes are strongly associated "
            f"with gene essentiality in {dataset_name}.\n"
        )

    path = os.path.join(output_dir, 'results_draft.txt')
    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    logger.info(f"  Saved results_draft.txt")


def write_cross_dataset_csv(all_data: dict, output_dir: str):
    """Write full cross-dataset metrics CSV."""
    rows = []
    metrics = ['AUC', 'Accuracy', 'F1_macro', 'F1_lethal', 'F1_viable']
    for ds in sorted(all_data.keys()):
        for model, entry in sorted(all_data[ds].items()):
            row = {'Dataset': ds, 'Model': MODEL_LABELS_SHORT.get(model, model)}
            for met in metrics:
                row[met] = _fmt(entry.get('metrics', {}).get(met), 4)
            rows.append(row)
    df = pd.DataFrame(rows)
    path = os.path.join(output_dir, 'cross_dataset_metrics.csv')
    df.to_csv(path, index=False)
    logger.info(f"  Saved cross_dataset_metrics.csv")
    return df


def write_cross_results_draft(all_data: dict, output_dir: str):
    """Cross-dataset manuscript Results paragraph."""
    datasets = sorted(all_data.keys())
    lines = []
    lines.append(f"### Cross-dataset Results draft\n")
    lines.append(
        f"Gene essentiality prediction with PhenGO-Predict was evaluated across "
        f"{len(datasets)} datasets ({', '.join(datasets)}).  "
        f"All classifiers were trained and evaluated independently for each dataset "
        f"using identical hyperparameters and a stratified 80:20 train/test split.\n"
    )

    # Best model per dataset
    best_rows = []
    for ds in datasets:
        best_m, best_a = None, 0.0
        for m, entry in all_data[ds].items():
            a = entry.get('metrics', {}).get('AUC', 0.0)
            if a > best_a:
                best_a, best_m = a, m
        if best_m:
            best_rows.append(
                f"{ds}: {MODEL_LABELS_SHORT.get(best_m, best_m)} (AUC = {_fmt(best_a)})"
            )
    lines.append(
        "The best-performing model (by AUC-ROC) for each dataset was: "
        + "; ".join(best_rows) + ".\n"
    )

    path = os.path.join(output_dir, 'cross_dataset_results_draft.txt')
    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    logger.info(f"  Saved cross_dataset_results_draft.txt")


# ── Main ──────────────────────────────────────────────────────────────────────

def run_single(predict_dir: str, output_dir: str, dataset_name: str, top_n: int):
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Single-dataset report: {dataset_name}")
    logger.info(f"  Predict dir:  {predict_dir}")
    logger.info(f"  Output dir:   {output_dir}")

    model_data = collect_model_data(predict_dir)
    if not model_data:
        logger.error(f"No model outputs found in {predict_dir}")
        return

    logger.info(f"  Models found: {', '.join(sorted(model_data))}")

    logger.info("Generating figures...")
    fig_metrics_comparison(model_data, output_dir, dataset_name)
    fig_roc_curves(model_data, output_dir, dataset_name)
    fig_confusion_matrices(model_data, output_dir, dataset_name)
    fig_feature_importance_heatmap(model_data, output_dir, dataset_name, top_n)
    fig_prediction_agreement(model_data, output_dir, dataset_name)

    logger.info("Generating text outputs...")
    write_metrics_csv(model_data, output_dir, dataset_name)
    write_results_draft(model_data, output_dir, dataset_name)

    logger.info("Done.")


def run_multi(all_predict_dirs: dict, output_dir: str, top_n: int):
    """all_predict_dirs: {dataset_name: predict_dir_path}"""
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Multi-dataset report: {len(all_predict_dirs)} datasets")

    all_data = {}
    for ds_name, pdir in sorted(all_predict_dirs.items()):
        model_data = collect_model_data(pdir)
        if model_data:
            all_data[ds_name] = model_data
            logger.info(f"  {ds_name}: {', '.join(sorted(model_data))}")
        else:
            logger.warning(f"  {ds_name}: no model outputs found — skipping")

    if not all_data:
        logger.error("No data found for any dataset.")
        return

    # Per-dataset mini-reports inside sub-dirs
    for ds_name, model_data in all_data.items():
        ds_out = os.path.join(output_dir, ds_name.replace(' ', '_'))
        os.makedirs(ds_out, exist_ok=True)
        fig_metrics_comparison(model_data, ds_out, ds_name)
        fig_roc_curves(model_data, ds_out, ds_name)
        fig_confusion_matrices(model_data, ds_out, ds_name)
        fig_feature_importance_heatmap(model_data, ds_out, ds_name, top_n)
        write_metrics_csv(model_data, ds_out, ds_name)

    # Cross-dataset figures at the root
    logger.info("Generating cross-dataset figures...")
    for metric in ('AUC', 'Accuracy', 'F1_macro'):
        fig_cross_dataset_heatmap(all_data, output_dir, metric)
    fig_cross_dataset_lines(all_data, output_dir)
    fig_cross_feature_overlap(all_data, output_dir, top_n)
    fig_best_model_summary(all_data, output_dir)

    logger.info("Generating cross-dataset text outputs...")
    write_cross_dataset_csv(all_data, output_dir)
    write_cross_results_draft(all_data, output_dir)

    logger.info("Done.")


def main():
    parser = argparse.ArgumentParser(
        description=f"PhenGO {PhenGO_VERSION} — Predict Report Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input modes (mutually exclusive groups)
    input_g = parser.add_mutually_exclusive_group(required=True)
    input_g.add_argument('-predict_dir',
                         help='Path to a single Predict output directory (single-dataset mode)')
    input_g.add_argument('-predict_dirs', nargs='+',
                         help='Multiple Predict output dirs (multi-dataset mode)')
    input_g.add_argument('-input_dir',
                         help='Parent directory; auto-discovers sub-dirs that contain a Predict/ folder')

    parser.add_argument('-dataset_name', default=None,
                        help='Display name for single-dataset mode (default: directory basename)')
    parser.add_argument('-dataset_names', nargs='+',
                        help='Names for each dir in -predict_dirs (default: dir basenames)')
    parser.add_argument('-output_dir', required=True,
                        help='Directory to write all report outputs')
    parser.add_argument('--top-n', dest='top_n', type=int, default=25,
                        help='Top N GO terms to include in feature plots (default: 25)')

    args = parser.parse_args()

    log_dir = os.path.abspath(args.output_dir)
    os.makedirs(log_dir, exist_ok=True)
    global logger
    logger = configure_logger('PhenGO.predict_report', enable_file=True,
                              log_dir=log_dir, logfile_name='predict_report.log',
                              level=logging.INFO)

    # ── Single dataset ──────────────────────────────────────────────────────
    if args.predict_dir:
        pdir = os.path.abspath(args.predict_dir)
        name = args.dataset_name or os.path.basename(os.path.dirname(pdir))
        run_single(pdir, os.path.abspath(args.output_dir), name, args.top_n)

    # ── Explicit multi-dataset ──────────────────────────────────────────────
    elif args.predict_dirs:
        pdirs = [os.path.abspath(p) for p in args.predict_dirs]
        if args.dataset_names:
            if len(args.dataset_names) != len(pdirs):
                logger.error("-dataset_names count must match -predict_dirs count")
                return
            names = args.dataset_names
        else:
            names = [os.path.basename(os.path.dirname(p)) for p in pdirs]
        run_multi(dict(zip(names, pdirs)), os.path.abspath(args.output_dir), args.top_n)

    # ── Auto-discover from parent dir ───────────────────────────────────────
    else:
        parent = os.path.abspath(args.input_dir)
        found = {}
        for entry in sorted(os.listdir(parent)):
            candidate = os.path.join(parent, entry, 'Predict')
            if os.path.isdir(candidate):
                found[entry] = candidate
            # Also accept if the dir itself looks like a Predict dir
            elif os.path.isdir(os.path.join(parent, entry)) and any(
                    sub in os.listdir(os.path.join(parent, entry))
                    for sub in ('nn', 'rf', 'dt', 'gb', 'lr', 'svm')):
                found[entry] = os.path.join(parent, entry)
        if not found:
            logger.error(f"No Predict output directories found under {parent}")
            return
        logger.info(f"Auto-discovered {len(found)} dataset(s): {', '.join(sorted(found))}")
        if len(found) == 1:
            name, pdir = list(found.items())[0]
            run_single(pdir, os.path.abspath(args.output_dir), name, args.top_n)
        else:
            run_multi(found, os.path.abspath(args.output_dir), args.top_n)


if __name__ == '__main__':
    main()



