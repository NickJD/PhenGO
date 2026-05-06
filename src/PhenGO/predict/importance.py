"""
predict/importance.py — Permutation-based GO-term feature importance for PhenGO-Predict.
"""
import os
import time
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def analyse_feature_importance(options, model, feature_names, X_test, y_test,
                               n_repeats, top_n=20):
    """Compute permutation-based feature importance with class-specific breakdown.

    For each GO-term feature, shuffles that column n_repeats times and measures
    the mean accuracy drop compared to the unshuffled baseline.  Reports both
    overall and per-class (lethal / viable) drops.

    Args:
        options      : Namespace with ``output_dir`` attribute.
        model        : Trained Keras model.
        feature_names: Sequence of GO-term column names.
        X_test       : Feature matrix (numpy array).
        y_test       : Integer label array.
        n_repeats    : Number of permutation repeats per feature.
        top_n        : Number of top features to log to console.

    Returns:
        (overall_df, lethal_df, viable_df) — three DataFrames sorted by importance.
    """
    logger.info(f"Analysing feature importance with {n_repeats} repeats ...")
    logger.info(f"({X_test.shape[1] * n_repeats} total permutations)")

    baseline_preds = model.predict(X_test, verbose=0).flatten()
    baseline_acc   = np.mean((baseline_preds > 0.5) == y_test)
    start_time     = time.time()

    lethal_mask  = y_test == 1
    viable_mask  = y_test == 0

    lethal_baseline = (np.mean((baseline_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
                       if np.sum(lethal_mask) > 0 else np.nan)
    viable_baseline = (np.mean((baseline_preds[viable_mask] > 0.5) == y_test[viable_mask])
                       if np.sum(viable_mask) > 0 else np.nan)

    overall_mean, overall_std = [], []
    lethal_mean,  lethal_std  = [], []
    viable_mean,  viable_std  = [], []

    total_iterations = X_test.shape[1] * n_repeats

    for i in range(X_test.shape[1]):
        drops_overall, drops_lethal, drops_viable = [], [], []
        X_base = X_test.copy()

        for rep in range(n_repeats):
            iteration = i * n_repeats + rep + 1
            if iteration % 10 == 0:
                elapsed      = time.time() - start_time
                per_it       = elapsed / iteration
                remaining    = (total_iterations - iteration) * per_it
                logger.debug(f"Progress: {iteration}/{total_iterations} "
                             f"({100*iteration/total_iterations:.1f}%) - "
                             f"ETA: {remaining/60:.1f} min")

            X_permuted = X_base.copy()
            np.random.shuffle(X_permuted[:, i])

            permuted_preds = model.predict(X_permuted, verbose=0).flatten()
            perm_acc       = np.mean((permuted_preds > 0.5) == y_test)
            drops_overall.append(baseline_acc - perm_acc)

            if np.sum(lethal_mask) > 0:
                perm_lethal = np.mean((permuted_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
                drops_lethal.append(lethal_baseline - perm_lethal)
            if np.sum(viable_mask) > 0:
                perm_viable = np.mean((permuted_preds[viable_mask] > 0.5) == y_test[viable_mask])
                drops_viable.append(viable_baseline - perm_viable)

        overall_mean.append(np.mean(drops_overall))
        overall_std.append(np.std(drops_overall))
        lethal_mean.append(np.mean(drops_lethal) if drops_lethal else np.nan)
        lethal_std.append(np.std(drops_lethal)   if drops_lethal else np.nan)
        viable_mean.append(np.mean(drops_viable) if drops_viable else np.nan)
        viable_std.append(np.std(drops_viable)   if drops_viable else np.nan)

    total_time = time.time() - start_time
    logger.info(f"Feature importance completed in {total_time/60:.1f} minutes")

    try:
        weights       = model.get_layer("hidden1").get_weights()[0]
        weight_contrib = np.mean(np.abs(weights), axis=1)
    except Exception:
        weight_contrib = np.zeros(X_test.shape[1])

    n = len(overall_mean)
    importance_df = pd.DataFrame({
        "GO_Term_Index":     range(n),
        "GO_Term":           list(feature_names) if feature_names is not None
                             else [f"GO_Term_{i}" for i in range(n)],
        "Overall_Importance": overall_mean,
        "Overall_Std":        overall_std,
        "Lethal_Importance":  lethal_mean,
        "Lethal_Std":         lethal_std,
        "Viable_Importance":  viable_mean,
        "Viable_Std":         viable_std,
        "Weight_Contribution": weight_contrib,
    })

    overall_df = importance_df.sort_values("Overall_Importance", ascending=False)
    lethal_df  = importance_df.sort_values("Lethal_Importance",  ascending=False)
    viable_df  = importance_df.sort_values("Viable_Importance",  ascending=False)

    overall_df.to_csv(os.path.join(options.output_dir, "overall_feature_importance.csv"), index=False)
    lethal_df.to_csv(os.path.join(options.output_dir,  "lethal_feature_importance.csv"),  index=False)
    viable_df.to_csv(os.path.join(options.output_dir,  "viable_feature_importance.csv"),  index=False)

    logger.info(f"=== TOP {top_n} GO TERMS (Overall Importance) ===")
    logger.info("\n" + overall_df[
        ["GO_Term","Overall_Importance","Overall_Std"]
    ].head(top_n).to_string(index=False))

    from .visualise import plot_feature_importance_with_errors
    plot_feature_importance_with_errors(options, overall_df, top_n=30)

    return overall_df, lethal_df, viable_df
