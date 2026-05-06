"""
predict/evaluate.py — Class-weight computation, threshold optimisation
and per-gene prediction evaluation for PhenGO-Predict.
"""
import logging
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, classification_report, confusion_matrix, roc_auc_score

logger = logging.getLogger(__name__)

try:
    import tensorflow as tf
except ImportError as e:
    logger.error(f"TensorFlow import failed: {e}")
    raise


def compute_class_weights(y):
    """Compute balanced class weights to counteract class imbalance.

    Gene-essentiality datasets are almost always heavily skewed: viable genes
    typically outnumber lethal genes 3-10:1.  Without compensation the model
    learns to predict "viable" for everything.

    Returns dict {0: w0, 1: w1} suitable for model.fit(class_weight=...).
    Formula: w_i = n_total / (n_classes * n_i)  (sklearn balanced weighting).
    """
    counts = np.bincount(y)
    n_total = len(y)
    n_classes = len(counts)
    weights = {i: n_total / (n_classes * c) for i, c in enumerate(counts)}
    logger.info(f"Class weights: { {k: f'{v:.3f}' for k, v in weights.items()} }")
    return weights


def find_optimal_threshold(y_true, y_pred_proba):
    """Find the probability threshold that maximises the F1 score for the lethal class.

    The default threshold of 0.5 is only optimal when classes are balanced.
    With imbalanced data the optimal cut-point for lethal genes is typically < 0.5.

    Returns (best_threshold float, best_f1 float).
    """
    best_threshold = 0.5
    best_f1 = 0.0

    for t in np.arange(0.1, 0.91, 0.01):
        y_pred = (y_pred_proba >= t).astype(int)
        f1 = f1_score(y_true, y_pred, pos_label=1, zero_division=0)
        if f1 > best_f1:
            best_f1 = f1
            best_threshold = round(float(t), 2)

    logger.info(f"Optimal threshold: {best_threshold:.2f}  (F1={best_f1:.3f} for lethal class)")
    logger.info("  (Default 0.5 threshold only optimal for balanced classes;")
    logger.info("   lethal genes are typically 3-10x rarer than viable genes)")
    return best_threshold, best_f1


def evaluate_and_analyse_predictions(model, X_test, y_test, gene_names_test,
                                     label_encoder, threshold=0.5):
    """Comprehensive per-gene evaluation of a trained model.

    Args:
        model          : Trained Keras model.
        X_test         : Feature matrix for test genes.
        y_test         : Integer labels for test genes.
        gene_names_test: Series/array of gene name strings.
        label_encoder  : Fitted LabelEncoder used during prepare_data().
        threshold      : Decision threshold (use find_optimal_threshold output).

    Returns:
        (predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score)
    """
    logger.info(f"Evaluating model (threshold={threshold:.2f}) ...")

    y_pred_proba = model.predict(X_test, verbose=0)
    y_pred = (y_pred_proba >= threshold).astype(int).flatten()

    test_results = model.evaluate(X_test, y_test, verbose=0)
    test_loss      = test_results[0]
    test_acc       = test_results[1]
    test_precision = test_results[2] if len(test_results) > 2 else 0.0
    test_recall    = test_results[3] if len(test_results) > 3 else 0.0

    try:
        auc_score = roc_auc_score(y_test, y_pred_proba)
    except Exception as e:
        logger.warning(f"Could not calculate AUC: {e}")
        auc_score = 0.0

    logger.info(f"Loss:      {test_loss:.4f}")
    logger.info(f"Accuracy:  {test_acc:.4f}")
    logger.info(f"Precision: {test_precision:.4f}")
    logger.info(f"Recall:    {test_recall:.4f}")
    logger.info(f"AUC:       {auc_score:.4f}")
    logger.info("Classification Report:\n" +
                classification_report(y_test, y_pred, target_names=label_encoder.classes_))
    cm = confusion_matrix(y_test, y_pred)
    logger.info("Confusion Matrix:\n" + np.array2string(cm))

    predictions_df = pd.DataFrame({
        "Gene_Name":                gene_names_test.reset_index(drop=True),
        "True_Label":               [label_encoder.classes_[i] for i in y_test],
        "True_Label_Numeric":       y_test,
        "Predicted_Label":          [label_encoder.classes_[i] for i in y_pred],
        "Predicted_Label_Numeric":  y_pred,
        "Prediction_Probability":   y_pred_proba.flatten(),
        "Confidence":               np.abs(y_pred_proba.flatten() - 0.5),
        "Correct_Prediction":       y_test == y_pred,
    })

    predictions_df["Prediction_Category"] = predictions_df.apply(
        lambda row: "True Positive"  if (row["True_Label_Numeric"] == 1 and row["Predicted_Label_Numeric"] == 1)
               else "True Negative"  if (row["True_Label_Numeric"] == 0 and row["Predicted_Label_Numeric"] == 0)
               else "False Positive" if (row["True_Label_Numeric"] == 0 and row["Predicted_Label_Numeric"] == 1)
               else "False Negative", axis=1,
    )
    predictions_df = predictions_df.sort_values("Confidence", ascending=False)

    logger.info(f"=== PREDICTION SUMMARY ===")
    logger.info(f"Total test samples:   {len(predictions_df)}")
    logger.info(f"Correct predictions:  {predictions_df['Correct_Prediction'].sum()}")
    logger.info(f"Incorrect predictions:{(~predictions_df['Correct_Prediction']).sum()}")
    for category in predictions_df["Prediction_Category"].unique():
        count = (predictions_df["Prediction_Category"] == category).sum()
        logger.info(f"  {category}: {count}")

    logger.info("=== TOP 10 MOST CONFIDENT CORRECT PREDICTIONS ===")
    correct_preds = predictions_df[predictions_df["Correct_Prediction"]].head(10)
    logger.info("\n" + correct_preds[
        ["Gene_Name","True_Label","Prediction_Probability","Confidence"]
    ].to_string(index=False))

    logger.info("=== TOP 10 MOST CONFIDENT INCORRECT PREDICTIONS ===")
    incorrect_preds = predictions_df[~predictions_df["Correct_Prediction"]].head(10)
    if len(incorrect_preds) > 0:
        logger.info("\n" + incorrect_preds[
            ["Gene_Name","True_Label","Predicted_Label","Prediction_Probability","Confidence"]
        ].to_string(index=False))
    else:
        logger.info("No incorrect predictions!")

    logger.info("=== 10 MOST UNCERTAIN PREDICTIONS ===")
    uncertain_preds = predictions_df.nsmallest(10, "Confidence")
    logger.info("\n" + uncertain_preds[
        ["Gene_Name","True_Label","Predicted_Label","Prediction_Probability","Confidence"]
    ].to_string(index=False))

    return predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score
