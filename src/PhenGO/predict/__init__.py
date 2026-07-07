"""
PhenGO Predict — Neural-network gene-essentiality prediction platform.

Submodules
----------
data        : ARFF loading and feature/label preparation
model       : Keras model architectures
evaluate    : Class-weight computation, threshold optimisation, per-gene evaluation
importance  : Permutation-based GO-term feature importance
visualise   : Training-history, ROC and feature-importance plots
compare     : Multi-dataset comparison and differential-importance reporting
utils       : Cross-validation, unlabelled-gene prediction, GO-enrichment in predictions
cli         : Command-line entry point (main())
"""

# Submodule imports are deferred so that merely importing PhenGO.predict
# (e.g. during the direct-execution bootstrap) does not eagerly load
# TensorFlow or NumPy before the user's run configuration is ready.
def __getattr__(name):
    _owners = {
        "load_arff_data": "data",
        "prepare_data": "data",
        "create_model": "model",
        "create_model_sparse_optimised": "model",
        "compute_class_weights": "evaluate",
        "find_optimal_threshold": "evaluate",
        "evaluate_and_analyse_predictions": "evaluate",
        "analyse_feature_importance": "importance",
        "plot_roc_curve": "visualise",
        "plot_training_history": "visualise",
        "plot_feature_importance_with_errors": "visualise",
        "compare_datasets": "compare",
        "predict_unlabeled_genes": "utils",
        "cross_validate_model": "utils",
        "run_cross_year_transfer": "version_sensitivity",
        "run_within_year_cv": "version_sensitivity",
    }
    if name not in _owners:
        raise AttributeError(f"module 'PhenGO.predict' has no attribute {name!r}")
    from importlib import import_module

    module = import_module(f"{__name__}.{_owners[name]}")
    return getattr(module, name)

__all__ = [
    "load_arff_data", "prepare_data",
    "create_model", "create_model_sparse_optimised",
    "compute_class_weights", "find_optimal_threshold",
    "evaluate_and_analyse_predictions",
    "analyse_feature_importance",
    "plot_roc_curve", "plot_training_history",
    "plot_feature_importance_with_errors",
    "compare_datasets",
    "predict_unlabeled_genes", "cross_validate_model",
    "run_cross_year_transfer", "run_within_year_cv",
]
