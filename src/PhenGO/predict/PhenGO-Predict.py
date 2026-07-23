"""
predict/PhenGO-Predict.py — Command-line entry point for PhenGO-Predict.

Usage
-----
Single dataset (logistic-regression baseline, default)::

    PhenGO-predict -arff_file data.arff -output_dir results/

Single dataset with multiple models::

    PhenGO-predict -arff_file data.arff -output_dir results/ \\
                   -model nn rf dt gb lr

Run every built-in model and compare::

    PhenGO-predict -arff_file data.arff -output_dir results/ -model all

Compare two datasets::

    PhenGO-predict -arff_files wt.arff mutant.arff \\
                   -dataset_names "Wild Type" Mutant \\
                   -output_dir comparison/

Available models
----------------
nn  — Neural Network (TensorFlow/Keras)
dt  — Decision Tree  (C4.5-style; interpretable)
rf  — Random Forest  (ensemble; usually best sklearn performer)
gb  — Gradient Boosting
lr  — Logistic Regression (strong sparse-data baseline)
svm — Support Vector Machine (RBF kernel)
all — Run all of the above
"""
import os
import sys
import argparse
import logging
import json
import re
from datetime import datetime

# ── Allow direct execution: python .../predict/PhenGO-Predict.py  ────────────
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here   = os.path.dirname(os.path.abspath(__file__))
    _src    = os.path.normpath(os.path.join(_here, '..', '..'))
    _phengo = os.path.dirname(_here)
    if _src not in sys.path:
        sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',         os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.predict', os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.predict'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo
# ─────────────────────────────────────────────────────────────────────────────

from ..constants import configure_logger
from ..provenance import prepare_output_dir

configure_logger("PhenGO.predict", enable_file=False)
logger = logging.getLogger("PhenGO.predict")

# All supported model type tokens
_ALL_SKLEARN = ['dt', 'rf', 'gb', 'lr', 'svm']
_ALL_MODELS  = ['nn'] + _ALL_SKLEARN


def _resolve_models(model_args):
    """Expand 'all' and deduplicate the model list while preserving order."""
    seen, out = set(), []
    for m in model_args:
        targets = _ALL_MODELS if m == 'all' else [m]
        for t in targets:
            if t not in seen:
                seen.add(t)
                out.append(t)
    return out


def _run_nn(options, X_train, X_test, y_train, y_test,
            genes_test, label_encoder, go_term_columns, output_dir):
    """Train and evaluate the neural network.  Returns a summary dict."""
    try:
        import tensorflow as tf
        from tensorflow import keras
    except ImportError as e:
        logger.error(f"TensorFlow not available — skipping NN: {e}")
        return None

    from .model import create_model_sparse_optimised
    from .evaluate import compute_class_weights, find_optimal_threshold, evaluate_and_analyse_predictions
    from .visualise import plot_training_history, plot_roc_curve
    from .importance import analyse_feature_importance
    from sklearn.metrics import confusion_matrix, classification_report, f1_score
    from sklearn.model_selection import train_test_split
    import numpy as np

    np.random.seed(options.seed)
    tf.keras.utils.set_random_seed(options.seed)
    try:
        tf.config.experimental.enable_op_determinism()
    except Exception:
        logger.warning("TensorFlow deterministic operations are not available")

    nn_dir = os.path.join(output_dir, 'nn')
    os.makedirs(nn_dir, exist_ok=True)

    logger.info(f"\n{'='*60}")
    logger.info("Training: Neural Network")
    logger.info('='*60)

    model = create_model_sparse_optimised(
        input_dim=X_train.shape[1],
        hidden_units=options.hidden_units,
        dropout_rate=options.dropout,
    )
    model.summary(print_fn=logger.info)

    logger.info(f"  Training for up to {options.epochs} epochs ...")
    stratify = y_train if len(np.unique(y_train)) > 1 else None
    X_fit, X_valid, y_fit, y_valid = train_test_split(
        X_train,
        y_train,
        test_size=0.2,
        random_state=options.seed,
        stratify=stratify,
    )
    fit_class_weights = compute_class_weights(y_fit)
    history = model.fit(
        X_fit, y_fit,
        epochs=options.epochs,
        batch_size=options.batch_size,
        validation_data=(X_valid, y_valid),
        class_weight=fit_class_weights,
        callbacks=[
            keras.callbacks.EarlyStopping(
                monitor="val_loss", patience=15, restore_best_weights=True),
            keras.callbacks.ReduceLROnPlateau(
                monitor="val_loss", factor=0.5, patience=10, min_lr=1e-6),
        ],
        verbose=1,
    )
    logger.info("  Training complete.")

    y_pred_proba_valid = model.predict(X_valid, verbose=0).flatten()
    opt_threshold, opt_f1 = find_optimal_threshold(y_valid, y_pred_proba_valid)
    logger.info(f"  Validation-selected threshold {opt_threshold:.2f} (F1={opt_f1:.3f})")

    # Swap output_dir pointer so helpers write into nn/
    class _NNOpts:
        pass
    nn_opts = _NNOpts()
    nn_opts.__dict__.update(vars(options))
    nn_opts.output_dir = nn_dir

    predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score = \
        evaluate_and_analyse_predictions(
            model, X_test, y_test, genes_test, label_encoder,
            threshold=opt_threshold,
        )

    plot_training_history(nn_opts, history)

    y_pred_proba_test = model.predict(X_test, verbose=0)
    plot_roc_curve(nn_opts, y_test, y_pred_proba_test)

    logger.info("  Starting feature importance analysis ...")
    analyse_feature_importance(
        nn_opts, model, go_term_columns, X_test, y_test, options.perm_repeats,
    )

    model.save(os.path.join(nn_dir, "gene_essentiality_model.keras"))
    with open(os.path.join(nn_dir, "model_schema.json"), "w", encoding="utf-8") as handle:
        json.dump({
            "schema_version": 1,
            "model_type": "nn",
            "feature_names": list(go_term_columns),
            "class_mapping": {"viable": 0, "lethal": 1},
            "positive_class": "lethal",
            "threshold": float(opt_threshold),
        }, handle, indent=2, sort_keys=True)
        handle.write("\n")
    predictions_df.to_csv(os.path.join(nn_dir, "predictions.csv"), index=False)

    y_pred_bin = (y_pred_proba_test >= opt_threshold).astype(int).flatten()
    f1_macro = f1_score(y_test, y_pred_bin, average='macro', zero_division=0)
    cls_rep = classification_report(
        y_test, y_pred_bin, labels=[0, 1],
        target_names=label_encoder.classes_, output_dict=True, zero_division=0,
    )
    f1_lethal = cls_rep.get('lethal', {}).get('f1-score', float('nan'))
    f1_viable = cls_rep.get('viable', {}).get('f1-score', float('nan'))

    report_path = os.path.join(nn_dir, "report.txt")
    with open(report_path, "w") as rf:
        rf.write("=== NEURAL NETWORK REPORT ===\n\n")
        rf.write(f"Loss:      {test_loss:.3f}\n")
        rf.write(f"Accuracy:  {test_acc:.3f}\n")
        rf.write(f"Precision: {test_precision:.3f}\n")
        rf.write(f"Recall:    {test_recall:.3f}\n")
        rf.write(f"AUC:       {auc_score:.3f}\n")
        rf.write(f"F1_macro:  {f1_macro:.3f}\n")
        rf.write(f"F1_lethal: {f1_lethal:.3f}\n")
        rf.write(f"F1_viable: {f1_viable:.3f}\n")
        rf.write(f"Validation-selected threshold: {opt_threshold:.2f}\n")
        rf.write(f"Train samples: {len(X_train)}\n")
        rf.write(f"Test  samples: {len(X_test)}\n\n")
        rf.write("Confusion Matrix:\n")
        rf.write(str(confusion_matrix(y_test, y_pred_bin)) + "\n\n")
        rf.write("Classification Report:\n")
        rf.write(classification_report(
            y_test, y_pred_bin, labels=[0, 1],
            target_names=label_encoder.classes_, zero_division=0,
        ))

    logger.info(f"  Outputs written to: {nn_dir}/")
    return {
        'model_type': 'nn',
        'label':      'Neural Network',
        'accuracy':   test_acc,
        'precision':  test_precision,
        'recall':     test_recall,
        'auc':        auc_score,
        'f1_macro':   f1_macro,
        'f1_lethal':  f1_lethal,
        'f1_viable':  f1_viable,
    }


def main():
    parser = argparse.ArgumentParser(
        prog="PhenGO-predict",
        description="PhenGO-Predict — ML gene essentiality prediction platform",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("-arff_file",  help="Path to a single ARFF file")
    parser.add_argument("-arff_files", nargs="+",
                        help="Multiple ARFF files (comparison mode)")
    parser.add_argument("-dataset_names", nargs="+",
                        help="Display names for each ARFF file in comparison mode")

    # ── Model selection ───────────────────────────────────────────────────
    model_choices = _ALL_MODELS + ['all']
    parser.add_argument(
        "-model", nargs="+", default=["lr"],
        metavar="MODEL",
        help=(
            "Model(s) to train. Space-separated list or 'all'. "
            f"Choices: {', '.join(model_choices)}. "
            "Default: lr.  Example: -model nn rf dt"
        ),
    )

    # ── NN hyper-parameters ───────────────────────────────────────────────
    nn_group = parser.add_argument_group("Neural Network parameters")
    nn_group.add_argument("-epochs",       type=int,   default=100)
    nn_group.add_argument("-batch_size",   type=int,   default=32)
    nn_group.add_argument("-hidden_units", nargs="+",  type=int, default=[128, 64])
    nn_group.add_argument("-dropout",      type=float, default=0.3)

    # ── sklearn / shared hyper-parameters ────────────────────────────────
    sk_group = parser.add_argument_group("sklearn model parameters")
    sk_group.add_argument("-n_estimators", type=int, default=200,
                          help="Number of trees for RF and GB (default: 200)")
    sk_group.add_argument("-max_depth", type=int, default=None,
                          help="Max tree depth for DT, RF, GB (default: unlimited)")
    sk_group.add_argument("-calibration", choices=["none", "sigmoid", "isotonic"],
                          default="sigmoid",
                          help="Training-fold probability calibration (default: sigmoid)")
    sk_group.add_argument("-calibration_cv", type=int, default=3)
    sk_group.add_argument("-n_jobs", type=int, default=1,
                          help="Parallel workers inside supported estimators (default: 1)")

    # ── Shared parameters ─────────────────────────────────────────────────
    parser.add_argument("-test_size",    type=float, default=0.2)
    parser.add_argument("-seed", type=int, default=42)
    parser.add_argument("-perm_repeats", type=int,   default=5,
                        help="Permutation repeats for feature importance")
    parser.add_argument("-output_dir",   required=True,
                        help="Output directory")
    parser.add_argument("-overwrite", action="store_true",
                        help="Overwrite an existing non-empty output directory")

    options = parser.parse_args()

    # ── Validate model choices ─────────────────────────────────────────────
    for m in options.model:
        if m not in model_choices:
            parser.error(
                f"Unknown model '{m}'. Valid choices: {', '.join(model_choices)}")
    models_to_run = _resolve_models(options.model)

    if bool(options.arff_file) == bool(options.arff_files):
        parser.error("Provide exactly one of -arff_file or -arff_files")
    input_files = options.arff_files or [options.arff_file]
    missing = [path for path in input_files if not os.path.isfile(path)]
    if missing:
        parser.error("ARFF file(s) not found: " + ", ".join(missing))
    if not 0 < options.test_size < 1:
        parser.error("-test_size must be between 0 and 1")
    if options.calibration_cv < 2:
        parser.error("-calibration_cv must be at least 2")
    if options.n_jobs == 0:
        parser.error("-n_jobs cannot be 0")
    if options.arff_files:
        if not options.dataset_names:
            options.dataset_names = [
                f"Dataset_{index + 1}" for index in range(len(options.arff_files))
            ]
        elif len(options.dataset_names) != len(options.arff_files):
            parser.error("Number of dataset names must match number of ARFF files")
        if len(set(options.dataset_names)) != len(options.dataset_names):
            parser.error("Dataset names must be unique")
        safe_names = [
            re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._") or "dataset"
            for name in options.dataset_names
        ]
        if len(set(safe_names)) != len(safe_names):
            parser.error("Dataset names collide after filesystem-safe normalization")
    from .data import load_arff_data, prepare_data
    from collections import Counter
    from sklearn.model_selection import train_test_split as validate_split

    for path in input_files:
        frame, _ = load_arff_data(path)
        if frame is None:
            parser.error(f"Could not load ARFF file: {path}")
        try:
            _, labels, _, _ = prepare_data(frame)
        except ValueError as exc:
            parser.error(f"{path}: {exc}")
        label_counts = Counter(map(int, labels))
        if len(label_counts) < 2 or min(label_counts.values()) < 2:
            parser.error(f"{path}: both classes require at least two genes")
        try:
            training_labels, _ = validate_split(
                labels,
                test_size=options.test_size,
                random_state=options.seed,
                stratify=labels,
            )
            if 'nn' in models_to_run:
                validate_split(
                    training_labels,
                    test_size=0.2,
                    random_state=options.seed,
                    stratify=training_labels,
                )
        except ValueError as exc:
            parser.error(f"{path}: invalid stratified test/validation split: {exc}")
    try:
        options.output_dir = prepare_output_dir(
            options.output_dir,
            options.overwrite,
            protected_paths=input_files,
        )
    except ValueError as exc:
        parser.error(str(exc))

    configure_logger("PhenGO.predict",
                     enable_file=True,
                     log_dir=options.output_dir,
                     logfile_name="PhenGO_Predict.log")
    logger = logging.getLogger("PhenGO.predict")
    logger.info(f"Models requested: {', '.join(models_to_run)}")

    # ── Comparison mode ───────────────────────────────────────────────────
    if options.arff_files:
        logger.info("=" * 80)
        logger.info("RUNNING IN COMPARISON MODE")
        logger.info(f"Datasets: {', '.join(options.dataset_names)}")
        logger.info("=" * 80)
        from .compare import compare_datasets
        compare_datasets(
            options, options.arff_files, options.dataset_names, models_to_run,
        )
        logger.info("COMPARISON ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {options.output_dir}")

    # ── Single dataset mode ───────────────────────────────────────────────
    elif options.arff_file:
        from .data import load_arff_data, prepare_data
        from .sklearn_models import train_evaluate_sklearn_model, write_model_comparison
        from sklearn.model_selection import train_test_split
        import numpy as np

        logger.info("=" * 80)
        logger.info("RUNNING IN SINGLE DATASET MODE")
        logger.info(f"Models: {', '.join(models_to_run)}")
        logger.info("=" * 80)

        result = load_arff_data(options.arff_file)
        if result[0] is None:
            logger.error("Failed to load data. Exiting.")
            return

        df, meta = result
        X, y, gene_names, label_encoder = prepare_data(df)

        X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
            X, y, gene_names,
            test_size=options.test_size,
            random_state=options.seed,
            stratify=y,
        )

        logger.info(f"Training set : {X_train.shape[0]} samples  "
                    f"{np.bincount(y_train).tolist()}")
        logger.info(f"Test set     : {X_test.shape[0]} samples  "
                    f"{np.bincount(y_test).tolist()}")

        go_term_columns = list(df.columns[1:-1])
        all_results = []

        for model_type in models_to_run:
            if model_type == 'nn':
                summary = _run_nn(
                    options,
                    X_train, X_test, y_train, y_test,
                    genes_test, label_encoder,
                    go_term_columns, options.output_dir,
                )
            else:
                summary = train_evaluate_sklearn_model(
                    model_type,
                    X_train, y_train,
                    X_test,  y_test,
                    genes_test, label_encoder,
                    go_term_columns,
                    options.output_dir,
                    options,
                )

            if summary:
                all_results.append(summary)

        # ── Comparison summary (when >1 model ran) ─────────────────────────
        if len(all_results) > 1:
            write_model_comparison(all_results, options.output_dir)

        logger.info("=" * 80)
        logger.info("ANALYSIS COMPLETE")
        for r in sorted(all_results, key=lambda x: x['auc'], reverse=True):
            logger.info(
                f"  {r['label']:<35} "
                f"Acc={r['accuracy']:.3f}  "
                f"Prec={r['precision']:.3f}  "
                f"Rec={r['recall']:.3f}  "
                f"AUC={r['auc']:.3f}"
            )
        logger.info(f"Results saved to: {options.output_dir}")
        logger.info("=" * 80)

    else:
        logger.error("Must provide -arff_file or -arff_files")
        parser.print_help()
        return

    with open(os.path.join(options.output_dir, "PhenGO_Predict_params.txt"), "w") as pf:
        pf.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            pf.write(f"{arg}: {value}\n")

    logger.info("Thank you for using PhenGO-Predict!")
    logger.info("Documentation: https://github.com/NickJD/PhenGO")


if __name__ == "__main__":
    main()
