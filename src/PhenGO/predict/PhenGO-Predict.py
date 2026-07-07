"""
predict/PhenGO-Predict.py — Command-line entry point for PhenGO-Predict.

Usage
-----
Single dataset (neural network, default)::

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
import shutil
import argparse
import logging
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

    nn_dir = os.path.join(output_dir, 'nn')
    os.makedirs(nn_dir, exist_ok=True)

    logger.info(f"\n{'='*60}")
    logger.info("Training: Neural Network")
    logger.info('='*60)

    class_weights = compute_class_weights(y_train)

    model = create_model_sparse_optimised(
        input_dim=X_train.shape[1],
        hidden_units=options.hidden_units,
        dropout_rate=options.dropout,
    )
    model.summary(print_fn=logger.info)

    logger.info(f"  Training for up to {options.epochs} epochs ...")
    history = model.fit(
        X_train, y_train,
        epochs=options.epochs,
        batch_size=options.batch_size,
        validation_split=0.2,
        class_weight=class_weights,
        callbacks=[
            keras.callbacks.EarlyStopping(
                monitor="val_loss", patience=15, restore_best_weights=True),
            keras.callbacks.ReduceLROnPlateau(
                monitor="val_loss", factor=0.5, patience=10, min_lr=1e-6),
        ],
        verbose=1,
    )
    logger.info("  Training complete.")

    y_pred_proba_train = model.predict(X_train, verbose=0).flatten()
    opt_threshold, opt_f1 = find_optimal_threshold(y_train, y_pred_proba_train)
    logger.info(f"  Threshold {opt_threshold:.2f} (F1={opt_f1:.3f})")

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
    predictions_df.to_csv(os.path.join(nn_dir, "predictions.csv"), index=False)

    y_pred_bin = (y_pred_proba_test >= opt_threshold).astype(int).flatten()
    f1_macro = f1_score(y_test, y_pred_bin, average='macro', zero_division=0)
    cls_rep  = classification_report(y_test, y_pred_bin,
                                     target_names=label_encoder.classes_,
                                     output_dict=True)
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
        rf.write(f"Threshold: {opt_threshold:.2f}\n")
        rf.write(f"Train samples: {len(X_train)}\n")
        rf.write(f"Test  samples: {len(X_test)}\n\n")
        rf.write("Confusion Matrix:\n")
        rf.write(str(confusion_matrix(y_test, y_pred_bin)) + "\n\n")
        rf.write("Classification Report:\n")
        rf.write(classification_report(y_test, y_pred_bin,
                                       target_names=label_encoder.classes_))

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
        "-model", nargs="+", default=["nn"],
        metavar="MODEL",
        help=(
            "Model(s) to train. Space-separated list or 'all'. "
            f"Choices: {', '.join(model_choices)}. "
            "Default: nn.  Example: -model nn rf dt"
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

    # ── Shared parameters ─────────────────────────────────────────────────
    parser.add_argument("-test_size",    type=float, default=0.2)
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

    if os.path.exists(options.output_dir) and os.listdir(options.output_dir):
        if not options.overwrite:
            parser.error(
                f"Output directory '{options.output_dir}' is not empty. "
                "Choose a new directory or pass -overwrite."
            )
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir, exist_ok=True)

    configure_logger("PhenGO.predict",
                     enable_file=True,
                     log_dir=options.output_dir,
                     logfile_name="PhenGO_Predict.log")
    logger = logging.getLogger("PhenGO.predict")
    logger.info(f"Models requested: {', '.join(models_to_run)}")

    # ── Comparison mode ───────────────────────────────────────────────────
    if options.arff_files:
        if not options.dataset_names:
            options.dataset_names = [f"Dataset_{i+1}"
                                     for i in range(len(options.arff_files))]
        elif len(options.dataset_names) != len(options.arff_files):
            logger.error("Number of dataset names must match number of ARFF files")
            return
        logger.info("=" * 80)
        logger.info("RUNNING IN COMPARISON MODE")
        logger.info(f"Datasets: {', '.join(options.dataset_names)}")
        logger.info("=" * 80)
        from .compare import compare_datasets
        compare_datasets(options, options.arff_files, options.dataset_names)
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
            random_state=42,
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
