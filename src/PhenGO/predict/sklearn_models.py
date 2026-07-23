"""
sklearn_models.py — scikit-learn model training, evaluation, and feature importance
for PhenGO-Predict.

Supported model types
---------------------
dt  — Decision Tree  (C4.5-style; sklearn DecisionTreeClassifier)
rf  — Random Forest
gb  — Gradient Boosting
lr  — Logistic Regression
svm — Support Vector Machine (RBF kernel, probability=True)
nn  — Feed-forward Neural Network (sklearn MLPClassifier)
"""

import os
import json
import logging
import numpy as np
import pandas as pd

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import (
    roc_auc_score, classification_report, confusion_matrix,
    precision_score, recall_score, accuracy_score, roc_curve, f1_score,
    average_precision_score, balanced_accuracy_score, brier_score_loss,
    matthews_corrcoef,
)
from sklearn.calibration import CalibratedClassifierCV
from sklearn.utils.class_weight import compute_sample_weight

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Human-readable labels for each model type
MODEL_LABELS = {
    'dt':  'Decision Tree',
    'rf':  'Random Forest',
    'gb':  'Gradient Boosting',
    'lr':  'Logistic Regression',
    'svm': 'Support Vector Machine',
    'nn':  'Neural Network (MLP)',
}


def _positive_scores(model, X):
    """Return scores for class 1 without assuming probability-column order."""
    if hasattr(model, 'predict_proba'):
        classes = list(model.classes_)
        return model.predict_proba(X)[:, classes.index(1)]
    decision = np.asarray(model.decision_function(X), dtype=float)
    return 1 / (1 + np.exp(-decision))


def build_sklearn_model(model_type, options):
    """Return an unfitted sklearn classifier for *model_type*.

    Parameters
    ----------
    model_type : str
        One of ``dt``, ``rf``, ``gb``, ``lr``, ``svm``, ``nn``.
    options : argparse.Namespace
        CLI options; reads ``max_depth``, ``n_estimators`` when present.
    """
    max_depth    = getattr(options, 'max_depth',    None)
    n_estimators = getattr(options, 'n_estimators', 200)
    random_state = getattr(options, 'random_state', getattr(options, 'seed', 42))
    n_jobs = getattr(options, 'n_jobs', 1)
    calibration = getattr(options, 'calibration', 'none')

    if model_type == 'dt':
        return DecisionTreeClassifier(
            max_depth=max_depth,
            class_weight='balanced',
            random_state=random_state,
        )
    elif model_type == 'rf':
        return RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            class_weight='balanced',
            n_jobs=n_jobs,
            random_state=random_state,
        )
    elif model_type == 'gb':
        return GradientBoostingClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth or 5,
            random_state=random_state,
        )
    elif model_type == 'lr':
        return LogisticRegression(
            max_iter=2000,
            class_weight='balanced',
            n_jobs=n_jobs,
            random_state=random_state,
        )
    elif model_type == 'svm':
        return SVC(
            kernel='rbf',
            class_weight='balanced',
            # CalibratedClassifierCV consumes decision_function directly. Do
            # not run SVC's internal Platt fit and then calibrate it again.
            probability=(calibration == 'none'),
            random_state=random_state,
        )
    elif model_type == 'nn':
        return MLPClassifier(
            hidden_layer_sizes=tuple(getattr(options, 'nn_hidden_units', (128, 64))),
            activation='relu',
            solver='adam',
            alpha=getattr(options, 'nn_alpha', 0.0001),
            batch_size=getattr(options, 'nn_batch_size', 32),
            learning_rate_init=getattr(options, 'nn_learning_rate_init', 0.001),
            max_iter=getattr(options, 'nn_max_iter', 300),
            early_stopping=getattr(options, 'nn_early_stopping', True),
            validation_fraction=getattr(options, 'nn_validation_fraction', 0.1),
            n_iter_no_change=getattr(options, 'nn_n_iter_no_change', 15),
            random_state=random_state,
        )
    else:
        raise ValueError(f"Unknown sklearn model type: '{model_type}'")


def train_evaluate_sklearn_model(
    model_type,
    X_train, y_train,
    X_test,  y_test,
    genes_test, label_encoder,
    go_term_columns,
    output_dir,
    options,
):
    """Train, evaluate, and write outputs for one sklearn model.

    A subdirectory ``{output_dir}/{model_type}/`` is created for this model's
    outputs.

    Returns
    -------
    dict
        Summary metrics: model_type, label, accuracy, precision, recall, auc.
    """
    label     = MODEL_LABELS.get(model_type, model_type)
    model_dir = os.path.join(output_dir, model_type)
    os.makedirs(model_dir, exist_ok=True)

    logger.info(f"\n{'='*60}")
    logger.info(f"Training: {label}")
    logger.info(f"  Train samples : {len(X_train)}")
    logger.info(f"  Test  samples : {len(X_test)}")
    logger.info('='*60)

    base_model = build_sklearn_model(model_type, options)
    calibration = getattr(options, 'calibration', 'sigmoid')
    calibration_cv = min(
        getattr(options, 'calibration_cv', 3),
        int(np.bincount(np.asarray(y_train, dtype=int)).min()),
    )
    model = base_model
    if calibration != 'none' and calibration_cv >= 2:
        model = CalibratedClassifierCV(
            estimator=base_model,
            method=calibration,
            cv=calibration_cv,
        )
    sample_weight = (
        compute_sample_weight(class_weight='balanced', y=y_train)
        if model_type in {'gb', 'nn'}
        else None
    )
    if sample_weight is not None:
        model.fit(X_train, y_train, sample_weight=sample_weight)
    else:
        model.fit(X_train, y_train)
    logger.info("  Training complete.")

    # ── Metrics ────────────────────────────────────────────────────────────
    y_pred = model.predict(X_test)
    y_proba = _positive_scores(model, X_test)

    acc  = accuracy_score (y_test, y_pred)
    prec = precision_score(y_test, y_pred, zero_division=0)
    rec  = recall_score   (y_test, y_pred, zero_division=0)
    auc  = roc_auc_score(y_test, y_proba) if len(np.unique(y_test)) > 1 else float('nan')
    ap   = average_precision_score(y_test, y_proba) if len(np.unique(y_test)) > 1 else float('nan')
    bal_acc = balanced_accuracy_score(y_test, y_pred) if len(np.unique(y_test)) > 1 else float('nan')
    mcc = matthews_corrcoef(y_test, y_pred) if len(np.unique(y_test)) > 1 else float('nan')
    brier = brier_score_loss(y_test, y_proba)
    f1_macro = f1_score(y_test, y_pred, average='macro', zero_division=0)

    # Per-class F1
    cls_report = classification_report(
        y_test, y_pred, labels=[0, 1], target_names=label_encoder.classes_,
        output_dict=True, zero_division=0,
    )
    f1_lethal = cls_report.get('lethal', {}).get('f1-score', float('nan'))
    f1_viable = cls_report.get('viable', {}).get('f1-score', float('nan'))

    logger.info(f"  Accuracy:  {acc:.3f}")
    logger.info(f"  Precision: {prec:.3f}")
    logger.info(f"  Recall:    {rec:.3f}")
    logger.info(f"  AUC:       {auc:.3f}")
    logger.info(f"  AP:        {ap:.3f}")
    logger.info(f"  F1 macro:  {f1_macro:.3f}  (lethal {f1_lethal:.3f} / viable {f1_viable:.3f})")

    # ── Predictions CSV ────────────────────────────────────────────────────
    pd.DataFrame({
        'gene':               genes_test,
        'true_label':         label_encoder.inverse_transform(y_test),
        'predicted_label':    label_encoder.inverse_transform(y_pred),
        'probability_lethal': np.round(y_proba, 4),
        'correct':            (y_pred == y_test),
    }).to_csv(os.path.join(model_dir, 'predictions.csv'), index=False)

    # ── Text report ────────────────────────────────────────────────────────
    report_path = os.path.join(model_dir, 'report.txt')
    with open(report_path, 'w') as rf:
        rf.write(f"=== {label.upper()} REPORT ===\n\n")
        rf.write(f"Accuracy:  {acc:.3f}\n")
        rf.write(f"Precision: {prec:.3f}\n")
        rf.write(f"Recall:    {rec:.3f}\n")
        rf.write(f"AUC:       {auc:.3f}\n")
        rf.write(f"Average precision: {ap:.3f}\n")
        rf.write(f"Balanced accuracy: {bal_acc:.3f}\n")
        rf.write(f"MCC:       {mcc:.3f}\n")
        rf.write(f"Brier:     {brier:.3f}\n")
        rf.write(f"F1_macro:  {f1_macro:.3f}\n")
        rf.write(f"F1_lethal: {f1_lethal:.3f}\n")
        rf.write(f"F1_viable: {f1_viable:.3f}\n\n")
        rf.write(f"Train samples: {len(X_train)}\n")
        rf.write(f"Test  samples: {len(X_test)}\n\n")
        rf.write("Confusion Matrix (rows=true, cols=predicted):\n")
        rf.write(str(confusion_matrix(y_test, y_pred)) + "\n\n")
        rf.write("Classification Report:\n")
        rf.write(classification_report(
            y_test, y_pred, labels=[0, 1], target_names=label_encoder.classes_,
            zero_division=0,
        ))

    # ── ROC curve ──────────────────────────────────────────────────────────
    if len(np.unique(y_test)) > 1:
        fpr, tpr, _ = roc_curve(y_test, y_proba)
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.plot(fpr, tpr, lw=2, color='steelblue', label=f'AUC = {auc:.3f}')
        ax.plot([0, 1], [0, 1], 'k--', lw=1)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(f'ROC Curve — {label}')
        ax.legend(loc='lower right')
        fig.tight_layout()
        fig.savefig(os.path.join(model_dir, 'roc_curve.png'), dpi=150)
        plt.close(fig)
        logger.info("  ROC curve saved.")

    # ── Feature importance ─────────────────────────────────────────────────
    _write_feature_importance(
        model, model_type, X_test, y_test,
        go_term_columns, model_dir, options,
    )

    import joblib

    joblib.dump(model, os.path.join(model_dir, 'model.joblib'))
    with open(os.path.join(model_dir, 'model_schema.json'), 'w', encoding='utf-8') as handle:
        json.dump({
            'schema_version': 1,
            'model_type': model_type,
            'feature_names': list(go_term_columns),
            'class_mapping': {'viable': 0, 'lethal': 1},
            'positive_class': 'lethal',
            'threshold': 0.5,
            'calibration': calibration if calibration_cv >= 2 else 'none',
        }, handle, indent=2, sort_keys=True)
        handle.write('\n')

    logger.info(f"  Outputs written to: {model_dir}/")
    return {
        'model_type': model_type,
        'label':      label,
        'accuracy':   acc,
        'precision':  prec,
        'recall':     rec,
        'auc':        auc,
        'average_precision': ap,
        'balanced_accuracy': bal_acc,
        'mcc': mcc,
        'brier_score': brier,
        'f1_macro':   f1_macro,
        'f1_lethal':  f1_lethal,
        'f1_viable':  f1_viable,
    }


def _permutation_importance_with_classes(model, X_test, y_test,
                                          go_term_columns, n_repeats, n_top):
    """Run per-feature, per-class permutation importance using predict_proba.

    Mirrors the logic in ``importance.py`` (NN) so that all models produce
    identical output files.  Returns (overall_df, lethal_df, viable_df).
    """
    import time
    n_features = X_test.shape[1]
    total_iters = n_features * n_repeats

    baseline_proba = _positive_scores(model, X_test)
    baseline_pred  = (baseline_proba > 0.5).astype(int)
    baseline_acc   = np.mean(baseline_pred == y_test)

    lethal_mask = y_test == 1
    viable_mask = y_test == 0
    lethal_base = (np.mean(baseline_pred[lethal_mask] == y_test[lethal_mask])
                   if lethal_mask.sum() > 0 else np.nan)
    viable_base = (np.mean(baseline_pred[viable_mask] == y_test[viable_mask])
                   if viable_mask.sum() > 0 else np.nan)

    overall_mean, overall_std = [], []
    lethal_mean,  lethal_std  = [], []
    viable_mean,  viable_std  = [], []

    start = time.time()
    for i in range(n_features):
        drops_o, drops_l, drops_v = [], [], []
        for rep in range(n_repeats):
            X_perm = X_test.copy()
            rng = np.random.default_rng(seed=rep * n_features + i)
            rng.shuffle(X_perm[:, i])

            perm_proba = _positive_scores(model, X_perm)
            perm_pred  = (perm_proba > 0.5).astype(int)
            drops_o.append(baseline_acc - np.mean(perm_pred == y_test))
            if lethal_mask.sum() > 0:
                drops_l.append(lethal_base - np.mean(perm_pred[lethal_mask] == y_test[lethal_mask]))
            if viable_mask.sum() > 0:
                drops_v.append(viable_base - np.mean(perm_pred[viable_mask] == y_test[viable_mask]))

        overall_mean.append(float(np.mean(drops_o)))
        overall_std.append(float(np.std(drops_o)))
        lethal_mean.append(float(np.mean(drops_l)) if drops_l else np.nan)
        lethal_std.append(float(np.std(drops_l))   if drops_l else np.nan)
        viable_mean.append(float(np.mean(drops_v)) if drops_v else np.nan)
        viable_std.append(float(np.std(drops_v))   if drops_v else np.nan)

        done = (i + 1) * n_repeats
        if done % max(1, total_iters // 10) == 0:
            elapsed = time.time() - start
            eta = (total_iters - done) * elapsed / done
            logger.info(f"    Permutation progress: {done}/{total_iters} "
                        f"({100*done/total_iters:.0f}%)  ETA {eta/60:.1f} min")

    logger.info(f"  Permutation importance done in {(time.time()-start)/60:.1f} min")

    base_df = pd.DataFrame({
        'GO_Term':            list(go_term_columns),
        'Overall_Importance': overall_mean,
        'Overall_Std':        overall_std,
        'Lethal_Importance':  lethal_mean,
        'Lethal_Std':         lethal_std,
        'Viable_Importance':  viable_mean,
        'Viable_Std':         viable_std,
    })
    overall_df = base_df.sort_values('Overall_Importance', ascending=False).head(n_top)
    lethal_df  = base_df.sort_values('Lethal_Importance',  ascending=False).head(n_top)
    viable_df  = base_df.sort_values('Viable_Importance',  ascending=False).head(n_top)
    return overall_df, lethal_df, viable_df


def _write_feature_importance(model, model_type, X_test, y_test,
                               go_term_columns, output_dir, options):
    """Write feature importance CSVs and plots, matching the NN output format.

    For all models the following files are produced in *output_dir*:

    * ``overall_feature_importance.csv``   — ranked by overall permutation drop
    * ``lethal_feature_importance.csv``    — ranked by lethal-class drop
    * ``viable_feature_importance.csv``    — ranked by viable-class drop
    * ``feature_importance_with_errors.png`` — bar chart with std error bars

    Tree models (dt / rf / gb) additionally produce:
    * ``feature_importance_native.csv``    — raw Gini/impurity importances (fast
      reference; unaffected by test-set size and permutation noise)
    """
    n_top        = min(50, len(go_term_columns))
    perm_repeats = getattr(options, 'perm_repeats', 5)
    label        = MODEL_LABELS.get(model_type, model_type)

    # ── 1. Native importances (tree / LR only) ─────────────────────────────
    calibrated_estimators = [
        calibrated.estimator for calibrated in
        getattr(model, 'calibrated_classifiers_', [])
    ]
    native_estimators = calibrated_estimators or [model]

    if model_type in ('dt', 'rf', 'gb'):
        available = [
            estimator.feature_importances_ for estimator in native_estimators
            if hasattr(estimator, 'feature_importances_')
        ]
        importances = np.mean(available, axis=0) if available else None
        if importances is None:
            logger.warning("Native feature importance unavailable for %s", label)
        else:
            indices = np.argsort(importances)[::-1][:n_top]
            pd.DataFrame({
                'GO_Term':    [go_term_columns[i] for i in indices],
                'Importance': importances[indices],
            }).to_csv(os.path.join(output_dir, 'feature_importance_native.csv'), index=False)
            logger.info("  Native (Gini) importances written.")

    elif model_type == 'lr':
        available = [
            estimator.coef_[0] for estimator in native_estimators
            if hasattr(estimator, 'coef_')
        ]
        if available:
            coef = np.mean(available, axis=0)
            abs_c = np.abs(coef)
            indices = np.argsort(abs_c)[::-1][:n_top]
            pd.DataFrame({
                'GO_Term':     [go_term_columns[i] for i in indices],
                'Coefficient': coef[indices],
                'Abs_Coef':    abs_c[indices],
            }).to_csv(os.path.join(output_dir, 'feature_importance_native.csv'), index=False)
            logger.info("  Logistic Regression coefficients written.")

    # ── 2. Per-class permutation importance (all models) ───────────────────
    if perm_repeats > 0:
        logger.info(f"  Computing permutation importance "
                    f"({perm_repeats} repeats × {X_test.shape[1]} features) ...")
        overall_df, lethal_df, viable_df = _permutation_importance_with_classes(
            model, X_test, y_test, go_term_columns, perm_repeats, n_top,
        )
        overall_df.to_csv(os.path.join(output_dir, 'overall_feature_importance.csv'), index=False)
        lethal_df.to_csv(os.path.join(output_dir,  'lethal_feature_importance.csv'),  index=False)
        viable_df.to_csv(os.path.join(output_dir,  'viable_feature_importance.csv'),  index=False)
        logger.info("  Per-class importance CSVs written.")

        # ── Error-bar chart (top terms, sorted descending) ─────────────
        plot_df = overall_df.head(min(30, len(overall_df))).iloc[::-1]   # flip so top is highest
        fig, ax = plt.subplots(figsize=(10, max(5, len(plot_df) * 0.22)))
        y_pos   = np.arange(len(plot_df))
        ax.barh(y_pos, plot_df['Overall_Importance'],
                xerr=plot_df['Overall_Std'],
                align='center', color='steelblue', ecolor='black',
                capsize=3, height=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_df['GO_Term'], fontsize=7)
        ax.axvline(0, color='grey', linewidth=0.8, linestyle='--')
        ax.set_xlabel('Mean accuracy drop (permutation importance)')
        ax.set_title(f'Feature Importance with Errors — {label}\n'
                     f'(n_repeats={perm_repeats})', fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, 'feature_importance_with_errors.png'), dpi=150)
        plt.close(fig)
        logger.info("  Error-bar importance plot written.")

        # Also save a plain feature_importance.csv alias for backward compat
        overall_df[['GO_Term', 'Overall_Importance']].rename(
            columns={'Overall_Importance': 'Importance'}
        ).to_csv(os.path.join(output_dir, 'feature_importance.csv'), index=False)

    else:
        # perm_repeats == 0 → only native importances; generate a minimal alias
        if model_type in ('dt', 'rf', 'gb') and available:
            indices     = np.argsort(importances)[::-1][:n_top]
            plot_vals   = importances[indices][::-1]
            plot_names  = [go_term_columns[i] for i in indices][::-1]
            fig, ax = plt.subplots(figsize=(10, max(5, n_top * 0.18)))
            ax.barh(range(len(plot_vals)), plot_vals, color='steelblue')
            ax.set_yticks(range(len(plot_vals)))
            ax.set_yticklabels(plot_names, fontsize=7)
            ax.set_xlabel('Gini importance')
            ax.set_title(f'Feature Importance — {label}')
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, 'feature_importance.png'), dpi=150)
            plt.close(fig)
        logger.info("  Feature importance written (perm_repeats=0, native only).")


def write_model_comparison(results_list, output_dir):
    """Write a summary CSV comparing all models that were run.

    Parameters
    ----------
    results_list : list of dict
        Each dict is the return value of :func:`train_evaluate_sklearn_model`
        or the equivalent NN summary dict.
    output_dir : str
        Root output directory.
    """
    if len(results_list) < 2:
        return

    keep_cols = [c for c in ['label', 'accuracy', 'precision', 'recall', 'auc',
                              'f1_macro', 'f1_lethal', 'f1_viable']
                 if c in results_list[0]]
    df = pd.DataFrame(results_list)[keep_cols].sort_values('auc', ascending=False)
    rename_map = {
        'label': 'Model', 'accuracy': 'Accuracy', 'precision': 'Precision',
        'recall': 'Recall', 'auc': 'AUC',
        'f1_macro': 'F1_macro', 'f1_lethal': 'F1_lethal', 'f1_viable': 'F1_viable',
    }
    df = df.rename(columns=rename_map)

    out_path = os.path.join(output_dir, 'model_comparison.csv')
    df.to_csv(out_path, index=False)
    logger.info(f"\nModel comparison written to: {out_path}")
    logger.info(df.to_string(index=False))

    # Bar chart of AUC scores
    fig, ax = plt.subplots(figsize=(max(6, len(df) * 1.2), 5))
    colors = plt.cm.tab10.colors[:len(df)]
    bars = ax.bar(df['Model'], df['AUC'], color=colors, edgecolor='white', width=0.6)
    ax.bar_label(bars, fmt='%.3f', padding=3, fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('AUC (ROC)')
    ax.set_title('Model Comparison — AUC')
    ax.axhline(0.5, color='grey', linestyle='--', lw=1, label='Random baseline')
    ax.legend()
    plt.xticks(rotation=20, ha='right')
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'model_comparison_auc.png'), dpi=150)
    plt.close(fig)
    logger.info("  Comparison chart: model_comparison_auc.png")
