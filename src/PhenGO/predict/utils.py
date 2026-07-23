"""
predict/utils.py — Cross-validation, unlabelled-gene prediction and
GO-enrichment analysis utilities for PhenGO-Predict.
"""
import os
import logging
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

from .data     import load_arff_data
from .model    import create_model_sparse_optimised
from .evaluate import compute_class_weights

logger = logging.getLogger(__name__)


def _benjamini_hochberg(p_values):
    """Return Benjamini-Hochberg adjusted p-values in original order."""
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], 1):
        rank = len(values) - reverse_rank + 1
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted

def cross_validate_model(X, y, gene_names, options, n_folds=5):
    """Stratified K-fold cross-validation for robust performance estimates.

    Useful for small datasets or when confidence intervals on metrics are needed.

    Args:
        X         : Feature matrix.
        y         : Label array.
        gene_names: Gene name series (for reference).
        options   : Namespace with output_dir, epochs, batch_size, hidden_units, dropout.
        n_folds   : Number of CV folds (default 5).

    Returns:
        dict with lists of per-fold accuracy, precision, recall, auc.
    """
    logger.info("=" * 80)
    logger.info(f"PERFORMING {n_folds}-FOLD CROSS-VALIDATION")
    logger.info("=" * 80)

    try:
        from tensorflow import keras
    except ImportError as exc:
        raise RuntimeError("TensorFlow is required for neural-network CV") from exc

    seed = getattr(options, "seed", 42)
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    cv_results = {"accuracy": [], "precision": [], "recall": [], "auc": []}

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
        logger.info(f"Fold {fold}/{n_folds}")
        X_train_cv, X_test_cv = X[train_idx], X[test_idx]
        y_train_cv, y_test_cv = y[train_idx], y[test_idx]

        model = create_model_sparse_optimised(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout,
        )
        cv_class_weights = compute_class_weights(y_train_cv)
        model.fit(
            X_train_cv, y_train_cv,
            epochs=options.epochs,
            batch_size=options.batch_size,
            validation_split=0.2,
            class_weight=cv_class_weights,
            callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10)],
            verbose=0,
        )

        y_pred_proba = model.predict(X_test_cv, verbose=0)
        test_results  = model.evaluate(X_test_cv, y_test_cv, verbose=0)

        cv_results["accuracy"].append(test_results[1])
        cv_results["precision"].append(test_results[2] if len(test_results) > 2 else 0)
        cv_results["recall"].append(test_results[3]    if len(test_results) > 3 else 0)
        cv_results["auc"].append(roc_auc_score(y_test_cv, y_pred_proba))

        logger.info(f"  Accuracy: {cv_results['accuracy'][-1]:.3f}  "
                    f"AUC: {cv_results['auc'][-1]:.3f}")

    logger.info("=" * 80)
    logger.info("CROSS-VALIDATION RESULTS")
    logger.info("=" * 80)
    for metric, values in cv_results.items():
        logger.info(f"{metric.capitalize()}: {np.mean(values):.3f} +/- {np.std(values):.3f}")

    cv_df = pd.DataFrame(cv_results)
    cv_df.to_csv(os.path.join(options.output_dir, "cross_validation_results.csv"), index=False)
    return cv_results


def predict_unlabeled_genes(model, unlabeled_arff_file, output_file,
                            feature_names=None, threshold=None, schema_file=None):
    """Apply a trained model to genes without phenotype labels.

    Useful for predicting essentiality of newly discovered genes or
    transferring predictions between species.

    Args:
        model               : Trained Keras model.
        unlabeled_arff_file : ARFF file whose last column is absent/ignored.
        output_file         : CSV path to write predictions.

    Returns:
        DataFrame with columns Gene_Name, Predicted_Label, Lethal_Probability, Confidence.
    """
    logger.info("=" * 80)
    logger.info("PREDICTING PHENOTYPES FOR UNLABELLED GENES")
    logger.info("=" * 80)

    schema = {}
    if schema_file:
        with open(schema_file, encoding="utf-8") as handle:
            schema = json.load(handle)
    expected = list(feature_names or schema.get("feature_names") or
                    getattr(model, "feature_names_in_", []))
    threshold = float(threshold if threshold is not None else schema.get("threshold", 0.5))
    if not 0 < threshold < 1:
        raise ValueError("Prediction threshold must be between 0 and 1")
    if not expected:
        raise ValueError(
            "An exact training feature schema is required. Pass feature_names or schema_file."
        )

    df, _ = load_arff_data(unlabeled_arff_file)
    if df is None or df.empty:
        raise ValueError(f"Could not load unlabelled ARFF: {unlabeled_arff_file}")
    gene_names = df.iloc[:, 0].astype(str)
    candidate = df.iloc[:, 1:].copy()
    if candidate.columns[-1].lower() in {"class", "phenotype", "label"}:
        candidate = candidate.iloc[:, :-1]
    missing = sorted(set(expected) - set(candidate.columns))
    extra = sorted(set(candidate.columns) - set(expected))
    if missing or extra:
        raise ValueError(
            "Unlabelled ARFF schema differs from training schema. "
            f"Missing: {missing[:10]}; extra: {extra[:10]}"
        )
    features = candidate.reindex(columns=expected)
    try:
        X_unlabeled = features.astype(float).values
    except (TypeError, ValueError) as exc:
        raise ValueError("Unlabelled GO features must be numeric") from exc
    if not np.isfinite(X_unlabeled).all() or not set(np.unique(X_unlabeled)) <= {0.0, 1.0}:
        raise ValueError("Unlabelled GO features must contain only finite binary 0/1 values")

    if hasattr(model, "predict_proba"):
        classes = list(model.classes_)
        predictions_proba = model.predict_proba(features)[:, classes.index(1)]
    else:
        try:
            predictions_proba = model.predict(X_unlabeled, verbose=0).reshape(-1)
        except TypeError:
            predictions_proba = np.asarray(model.predict(X_unlabeled)).reshape(-1)
    predictions_binary = (predictions_proba >= threshold).astype(int)
    confidence = np.where(
        predictions_proba >= threshold,
        (predictions_proba - threshold) / max(1 - threshold, np.finfo(float).eps),
        (threshold - predictions_proba) / max(threshold, np.finfo(float).eps),
    )

    results_df = pd.DataFrame({
        "Gene_Name":          gene_names,
        "Predicted_Label":    ["lethal" if p == 1 else "viable" for p in predictions_binary],
        "Lethal_Probability": predictions_proba,
        "Confidence":         confidence,
    }).sort_values("Confidence", ascending=False)

    results_df.to_csv(output_file, index=False)

    n_lethal = sum(predictions_binary)
    n_viable = len(predictions_binary) - n_lethal
    logger.info(f"Predictions for {len(results_df)} genes saved to {output_file}")
    logger.info(f"  Predicted lethal: {n_lethal} ({n_lethal/len(predictions_binary)*100:.1f}%)")
    logger.info(f"  Predicted viable: {n_viable} ({n_viable/len(predictions_binary)*100:.1f}%)")
    logger.info(f"  Mean confidence:  {results_df['Confidence'].mean():.3f}")
    return results_df


def analyze_go_term_enrichment_in_predictions(predictions_df, arff_file, output_dir):
    """Identify GO terms enriched in correct vs. incorrect predictions.

    Helps diagnose systematic biases: which GO terms the model relies on
    correctly, and which terms mislead it.

    Args:
        predictions_df: DataFrame returned by evaluate_and_analyse_predictions().
        arff_file     : Original ARFF file (to retrieve GO annotations).
        output_dir    : Directory to save the enrichment CSV.

    Returns:
        DataFrame with Enrichment_Ratio column, sorted descending.
    """
    logger.info("=" * 80)
    logger.info("ANALYSING GO TERM ENRICHMENT IN PREDICTIONS")
    logger.info("=" * 80)

    df, meta = load_arff_data(arff_file)
    go_terms = df.iloc[:, 1:-1]

    merged = predictions_df.merge(
        pd.concat([df.iloc[:, 0], go_terms], axis=1),
        left_on="Gene_Name", right_on=df.columns[0],
    )

    correct_mask = merged["Correct_Prediction"].astype(bool)
    correct = merged[correct_mask]
    incorrect = merged[~correct_mask]

    from scipy.stats import fisher_exact

    enrichment_results = []
    for go_term in go_terms.columns:
        correct_freq   = correct[go_term].sum()   / len(correct)   if len(correct)   > 0 else 0
        incorrect_freq = incorrect[go_term].sum() / len(incorrect) if len(incorrect) > 0 else 0
        correct_with = int(correct[go_term].sum())
        incorrect_with = int(incorrect[go_term].sum())
        odds_ratio, p_value = fisher_exact([
            [correct_with, len(correct) - correct_with],
            [incorrect_with, len(incorrect) - incorrect_with],
        ])
        enrichment_results.append({
            "GO_Term":            go_term,
            "Correct_Frequency":  correct_freq,
            "Incorrect_Frequency": incorrect_freq,
            "Enrichment_Ratio":   (
                (correct_with + 0.5) / (len(correct) + 1) /
                ((incorrect_with + 0.5) / (len(incorrect) + 1))
            ),
            "Odds_Ratio":         odds_ratio,
            "P_Value":            p_value,
            "Difference":         correct_freq - incorrect_freq,
        })

    enrichment_df = pd.DataFrame(enrichment_results)
    if len(enrichment_df):
        enrichment_df["FDR"] = _benjamini_hochberg(enrichment_df["P_Value"])
    enrichment_df = enrichment_df.sort_values(["FDR", "Enrichment_Ratio"], ascending=[True, False])

    output_path = os.path.join(output_dir, "go_enrichment_in_predictions.csv")
    enrichment_df.to_csv(output_path, index=False)
    logger.info(f"GO enrichment results saved to {output_path}")
    return enrichment_df
