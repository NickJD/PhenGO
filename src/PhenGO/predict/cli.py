"""
predict/cli.py — Command-line entry point for PhenGO-Predict.

Usage
-----
Single dataset::

    PhenGO-predict -arff_file data.arff -output_dir results/

Compare two datasets::

    PhenGO-predict -arff_files wt.arff mutant.arff \\
                   -dataset_names "Wild Type" Mutant \\
                   -output_dir comparison/
"""
import os
import sys
import shutil
import argparse
import logging
from datetime import datetime

# ── Allow direct execution: python .../predict/cli.py  ─────────────────────
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
# ───────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────

from ..constants import configure_logger

configure_logger("PhenGO.predict", enable_file=False)
logger = logging.getLogger("PhenGO.predict")

try:
    import tensorflow as tf
    from tensorflow import keras
    logger.info(f"TensorFlow version: {tf.__version__}")
except ImportError as e:
    logger.error(f"TensorFlow import failed: {e}")
    logger.error("Install TensorFlow:  pip install tensorflow")
    raise SystemExit(1)

from .data      import load_arff_data, prepare_data
from .model     import create_model_sparse_optimised
from .evaluate  import (compute_class_weights, find_optimal_threshold,
                        evaluate_and_analyse_predictions)
from .visualise import plot_training_history, plot_roc_curve
from .importance import analyse_feature_importance
from .compare   import compare_datasets

from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report


def main():
    parser = argparse.ArgumentParser(
        description="PhenGO-Predict — Neural Network gene essentiality prediction platform",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("-arff_file",  help="Path to a single ARFF file")
    parser.add_argument("-arff_files", nargs="+",
                        help="Multiple ARFF files (comparison mode)")
    parser.add_argument("-dataset_names", nargs="+",
                        help="Display names for each ARFF file in comparison mode")
    parser.add_argument("-epochs",      type=int,   default=100)
    parser.add_argument("-batch_size",  type=int,   default=32)
    parser.add_argument("-test_size",   type=float, default=0.2)
    parser.add_argument("-hidden_units", nargs="+", type=int, default=[128, 64])
    parser.add_argument("-dropout",     type=float, default=0.3)
    parser.add_argument("-perm_repeats", type=int,  default=5,
                        help="Permutation repeats for feature importance")
    parser.add_argument("-output_dir",  required=True,
                        help="Output directory (existing contents will be deleted)")

    options = parser.parse_args()

    if os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir)

    configure_logger("PhenGO.predict",
                     enable_file=True,
                     log_dir=options.output_dir,
                     logfile_name="PhenGO_Predict.log")
    logger = logging.getLogger("PhenGO.predict")

    if options.arff_files:
        if not options.dataset_names:
            options.dataset_names = [f"Dataset_{i+1}" for i in range(len(options.arff_files))]
        elif len(options.dataset_names) != len(options.arff_files):
            logger.error("Number of dataset names must match number of ARFF files")
            return

        logger.info("=" * 80)
        logger.info("RUNNING IN COMPARISON MODE")
        logger.info(f"Datasets: {', '.join(options.dataset_names)}")
        logger.info("=" * 80)
        results = compare_datasets(options, options.arff_files, options.dataset_names)
        logger.info("COMPARISON ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {options.output_dir}")

    elif options.arff_file:
        logger.info("=" * 80)
        logger.info("RUNNING IN SINGLE DATASET MODE")
        logger.info("=" * 80)

        result = load_arff_data(options.arff_file)
        if result[0] is None:
            logger.error("Failed to load data. Exiting.")
            return

        df, meta = result
        X, y, gene_names, label_encoder = prepare_data(df)

        X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
            X, y, gene_names, test_size=options.test_size, random_state=42, stratify=y,
        )

        logger.info(f"Training set: {X_train.shape[0]} samples")
        logger.info(f"Test set:     {X_test.shape[0]} samples")
        logger.info(f"Train class distribution: {np.bincount(y_train)}")
        logger.info(f"Test  class distribution: {np.bincount(y_test)}")

        class_weights = compute_class_weights(y_train)

        model = create_model_sparse_optimised(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout,
        )
        model.summary(print_fn=logger.info)

        logger.info(f"Training for up to {options.epochs} epochs ...")
        try:
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
            logger.info("Training complete.")

            y_pred_proba_train = model.predict(X_train, verbose=0).flatten()
            opt_threshold, opt_f1 = find_optimal_threshold(y_train, y_pred_proba_train)
            logger.info(f"Using threshold {opt_threshold:.2f} (F1={opt_f1:.3f}) for evaluation")

            predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score = \
                evaluate_and_analyse_predictions(
                    model, X_test, y_test, genes_test, label_encoder,
                    threshold=opt_threshold,
                )

            plot_training_history(options, history)

            y_pred_proba_test = model.predict(X_test, verbose=0)
            plot_roc_curve(options, y_test, y_pred_proba_test)

            logger.info("Starting feature importance analysis ...")
            go_term_columns = df.columns[1:-1]
            analyse_feature_importance(
                options, model, go_term_columns, X_test, y_test, options.perm_repeats,
            )

            model.save(os.path.join(options.output_dir, "gene_essentiality_model.keras"))
            predictions_df.to_csv(
                os.path.join(options.output_dir, "test_predictions_detailed.csv"), index=False,
            )

            report_path = os.path.join(options.output_dir, "final_report.txt")
            with open(report_path, "w") as rf:
                rf.write("=== PHENGO PREDICTION REPORT ===\n\n")
                rf.write(f"Model Performance:\n")
                rf.write(f"  Loss:      {test_loss:.3f}\n")
                rf.write(f"  Accuracy:  {test_acc:.3f}\n")
                rf.write(f"  Precision: {test_precision:.3f}\n")
                rf.write(f"  Recall:    {test_recall:.3f}\n")
                rf.write(f"  AUC:       {auc_score:.3f}\n\n")
                rf.write("Confusion Matrix:\n")
                y_pred_bin = (y_pred_proba_test > 0.5).astype(int).flatten()
                rf.write(str(confusion_matrix(y_test, y_pred_bin)) + "\n\n")
                rf.write("Classification Report:\n")
                rf.write(classification_report(y_test, y_pred_bin,
                                               target_names=label_encoder.classes_))

            logger.info("=" * 80)
            logger.info("ANALYSIS COMPLETE")
            logger.info(f"Loss={test_loss:.3f}  Accuracy={test_acc:.3f}  AUC={auc_score:.3f}")
            logger.info(f"Full report: {report_path}")
            logger.info("=" * 80)

        except Exception as e:
            logger.exception(f"Error during training/analysis: {e}")
            return

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
    logger.info("Issues:        https://github.com/NickJD/PhenGO/issues")


if __name__ == "__main__":
    main()
