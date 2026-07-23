"""Version-sensitivity analysis for PhenGO ARFF datasets.

This module is designed for the paper-level question:

    Do machine-learning outputs change when the same organism is represented by
    different yearly snapshots of model-organism database resources?

The command intentionally separates several effects:

* within-year model performance
* train-year -> test-year transfer performance
* gene-set churn
* GO-feature-set churn
* label churn among shared genes
* prediction instability for the same genes across years
* shared-feature and matched-gene sensitivity panels
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from types import SimpleNamespace

from ..constants import PhenGO_VERSION, configure_logger
from ..provenance import (
    dependency_versions,
    git_commit,
    git_status,
    load_manifest_for_arff,
    prepare_output_dir,
    sha256_file,
    source_tree_sha256,
)

logger = logging.getLogger("PhenGO.version_sensitivity")

SKLEARN_MODELS = ["dt", "rf", "gb", "lr", "svm", "nn"]
DEFAULT_MODELS = ["lr", "rf", "gb", "dt"]
PANELS = ["full", "matched_features", "matched_genes", "matched_both"]
MATRIX_METRICS = ["roc_auc", "average_precision", "f1_lethal", "balanced_accuracy", "mcc"]
PUBLICATION_POLICY_KEYS = (
    "mixed_label",
    "nonlethal",
    "ambiguous",
    "label_source",
    "legacy_fly_lethal_override",
    "filter_multigenes",
    "fly_driver_filtering",
    "worm_evidence_codes",
    "go_relations",
    "go_namespaces",
    "go_propagation",
    "allow_cross_namespace_go_edges",
    "go_evidence_include",
    "go_evidence_exclude",
    "identifier_join_precedence",
    "strict_identifier_matching",
    "exclude_go_roots",
    "filter_unused_gos",
    "min_go_gene_count",
    "max_go_gene_fraction",
)


@dataclass
class VersionDataset:
    name: str
    path: str
    features: object
    labels: object
    label_text: object
    feature_names: list[str]
    manifest: dict | None = None

    @property
    def genes(self):
        return list(self.features.index)


def natural_key(value: str):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(value))]


def calendar_year(value: str) -> int | None:
    match = re.search(r"(?:19|20)\d{2}", str(value or ""))
    return int(match.group(0)) if match else None


def add_calendar_interval_fields(rows: list[dict]) -> None:
    for row in rows:
        train_year = calendar_year(row.get("train_dataset", ""))
        test_year = calendar_year(row.get("test_dataset", ""))
        row["train_calendar_year"] = train_year if train_year is not None else ""
        row["test_calendar_year"] = test_year if test_year is not None else ""
        if train_year is None or test_year is None:
            row["calendar_gap_years"] = ""
            row["absolute_calendar_gap_years"] = ""
            row["consecutive_calendar_years"] = ""
            continue
        gap = test_year - train_year
        row["calendar_gap_years"] = gap
        row["absolute_calendar_gap_years"] = abs(gap)
        row["consecutive_calendar_years"] = abs(gap) == 1


def expand_models(model_args: list[str]) -> list[str]:
    out = []
    for model in model_args:
        targets = SKLEARN_MODELS if model == "all" else [model]
        for target in targets:
            if target not in SKLEARN_MODELS:
                raise ValueError(f"Unknown model '{target}'. Valid choices: {SKLEARN_MODELS + ['all']}")
            if target not in out:
                out.append(target)
    return out


def discover_arff_files(input_dir: str) -> tuple[list[str], list[str]]:
    """Discover ARFF datasets from a parent directory.

    Each direct child directory containing ``*_PhenGO.arff`` is treated as one
    dataset. If no child datasets are found, direct ``*.arff`` files in
    ``input_dir`` are used.
    """
    input_dir = os.path.abspath(input_dir)
    discovered = []
    for entry in sorted(os.listdir(input_dir), key=natural_key):
        subdir = os.path.join(input_dir, entry)
        if not os.path.isdir(subdir):
            continue
        arffs = [f for f in os.listdir(subdir) if f.endswith("_PhenGO.arff")]
        if arffs:
            discovered.append((entry, os.path.join(subdir, sorted(arffs, key=natural_key)[0])))

    if not discovered:
        for entry in sorted(os.listdir(input_dir), key=natural_key):
            if entry.endswith(".arff"):
                name = os.path.splitext(entry)[0]
                discovered.append((name, os.path.join(input_dir, entry)))

    names = [name for name, _ in discovered]
    paths = [path for _, path in discovered]
    return paths, names


def ensure_output_dir(path: str, overwrite: bool, protected_paths=()):
    try:
        return prepare_output_dir(path, overwrite, protected_paths=protected_paths)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc


def load_version_dataset(path: str, name: str) -> VersionDataset:
    from .data import PhenotypeLabelEncoder, load_arff_data
    import pandas as pd

    df, _ = load_arff_data(path)
    if df is None:
        raise ValueError(f"Could not load ARFF dataset: {path}")

    gene_col = df.columns[0]
    label_col = df.columns[-1]
    feature_names = list(df.columns[1:-1])
    encoder = PhenotypeLabelEncoder()

    features = df.loc[:, feature_names].astype(float)
    features.index = df[gene_col].astype(str)
    labels = pd.Series(encoder.fit_transform(df[label_col]), index=features.index, name="label")
    label_text = pd.Series(df[label_col].astype(str).str.lower().values,
                           index=features.index, name="label_text")

    if len(set(feature_names)) != len(feature_names):
        raise ValueError(f"{name}: duplicate GO feature names are not permitted")
    if not features.notna().all().all():
        raise ValueError(f"{name}: missing or non-numeric feature values found")
    observed = set(features.to_numpy().ravel())
    if not observed <= {0.0, 1.0}:
        raise ValueError(f"{name}: non-binary GO values found: {sorted(observed)[:10]}")

    if features.index.has_duplicates:
        duplicate_ids = sorted(set(features.index[features.index.duplicated(False)]), key=natural_key)
        for gene in duplicate_ids:
            rows = features.loc[[gene]]
            gene_labels = labels.loc[[gene]]
            if len(set(gene_labels.astype(int))) > 1 or not rows.eq(rows.iloc[0]).all().all():
                raise ValueError(f"{name}: conflicting duplicate rows for gene {gene}")
        logger.warning("%s: removing %d exact duplicate gene rows", name, len(duplicate_ids))
        keep = ~features.index.duplicated(keep="first")
        features, labels, label_text = features.loc[keep], labels.loc[keep], label_text.loc[keep]

    return VersionDataset(
        name=name,
        path=path,
        features=features,
        labels=labels,
        label_text=label_text,
        feature_names=feature_names,
        manifest=load_manifest_for_arff(path),
    )


def validate_dataset_manifests(datasets: list[VersionDataset], allow_missing: bool = False):
    """Require traceable, release-specific inputs for paper-level comparisons."""
    missing = [dataset.name for dataset in datasets if dataset.manifest is None]
    if missing and not allow_missing:
        raise ValueError(
            "Missing PhenGO_manifest.json for dataset(s): " + ", ".join(missing) +
            ". Regenerate them with current PhenGO or pass -allow_missing_manifests "
            "for an explicitly exploratory legacy analysis."
        )

    species = set()
    snapshot_ids = []
    verified_manifests = []
    for dataset in datasets:
        manifest = dataset.manifest
        if manifest is None:
            logger.warning("%s has no manifest; provenance cannot be verified", dataset.name)
            continue
        verified_manifests.append((dataset.name, manifest))
        if int(manifest.get("schema_version", 0)) < 2 and not allow_missing:
            raise ValueError(
                f"{dataset.name}: manifest schema predates the publication label/GO "
                "contract; regenerate this ARFF with the current PhenGO core"
            )
        if manifest.get("species"):
            species.add(str(manifest["species"]).lower())
        snapshot_id = manifest.get("snapshot_id")
        if snapshot_id:
            snapshot_ids.append(str(snapshot_id))
        elif not allow_missing:
            raise ValueError(f"{dataset.name}: manifest has no snapshot_id")
        if not manifest.get("strict_snapshot") and not allow_missing:
            raise ValueError(
                f"{dataset.name}: manifest was not produced with -strict_snapshot"
            )
        recorded = (manifest.get("outputs") or {}).get("arff") or {}
        recorded_hash = recorded.get("sha256")
        if not recorded_hash and not allow_missing:
            raise ValueError(f"{dataset.name}: manifest has no recorded ARFF SHA-256")
        if recorded_hash and recorded_hash != sha256_file(dataset.path):
            raise ValueError(
                f"{dataset.name}: ARFF hash differs from its manifest; the file may "
                "have been changed after generation"
            )
        required_releases = {
            "phenotype", "go_annotations", "go_ontology", "retrieval_date",
        }
        species_name = str(manifest.get("species", "")).lower()
        label_source = (manifest.get("policies") or {}).get(
            "label_source", "release_records"
        )
        if species_name == "fly" or (
            species_name in {"worm", "mouse"} and label_source != "gene_sets"
        ):
            required_releases.add("phenotype_ontology")
        missing_releases = sorted(
            key for key in required_releases
            if not (manifest.get("releases") or {}).get(key)
        )
        if missing_releases and not allow_missing:
            raise ValueError(
                f"{dataset.name}: manifest lacks release metadata: " +
                ", ".join(missing_releases)
            )
        if not allow_missing:
            required_inputs = {"go_annotations", "go_ontology"}
            if label_source in {"release_records", "agreement"}:
                required_inputs.add("phenotype")
            if label_source in {"gene_sets", "agreement"}:
                required_inputs.update({"lethal_gene_set", "viable_gene_set"})
            missing_input_hashes = sorted(
                key for key in required_inputs
                if not ((manifest.get("inputs") or {}).get(key) or {}).get("sha256")
            )
            if missing_input_hashes:
                raise ValueError(
                    f"{dataset.name}: manifest lacks required input hashes: " +
                    ", ".join(missing_input_hashes)
                )
            policies = manifest.get("policies") or {}
            missing_policies = [key for key in PUBLICATION_POLICY_KEYS if key not in policies]
            if missing_policies:
                raise ValueError(
                    f"{dataset.name}: manifest lacks publication policy fields: " +
                    ", ".join(missing_policies)
                )
            if policies.get("legacy_fly_lethal_override"):
                raise ValueError(
                    f"{dataset.name}: deprecated fly lethal-only overrides are not "
                    "permitted in publication analyses"
                )
            if not manifest.get("source_tree_sha256"):
                raise ValueError(f"{dataset.name}: manifest has no source-tree hash")
            if not manifest.get("tool_version"):
                raise ValueError(f"{dataset.name}: manifest has no tool version")

    if len(species) > 1:
        raise ValueError(
            "Version-sensitivity datasets must belong to one organism; found: " +
            ", ".join(sorted(species))
        )
    if len(snapshot_ids) != len(set(snapshot_ids)):
        raise ValueError("Each dataset must have a unique snapshot_id")
    if not allow_missing and verified_manifests:
        reference_name, reference = verified_manifests[0]
        reference_policies = reference["policies"]
        for name, manifest in verified_manifests[1:]:
            changed = [
                key for key in PUBLICATION_POLICY_KEYS
                if manifest["policies"][key] != reference_policies[key]
            ]
            if changed:
                raise ValueError(
                    f"{name}: publication policies differ from {reference_name}: " +
                    ", ".join(changed)
                )
        for field, label in (
            ("tool_version", "tool version"),
            ("source_tree_sha256", "source-tree hash"),
            ("dependencies", "dependency environment"),
        ):
            values = {
                json.dumps(manifest.get(field), sort_keys=True)
                for _, manifest in verified_manifests
            }
            if len(values) > 1:
                raise ValueError(
                    f"Datasets were generated with different {label} values"
                )


def shared_features(datasets: list[VersionDataset]) -> list[str]:
    common = set(datasets[0].feature_names)
    for dataset in datasets[1:]:
        common &= set(dataset.feature_names)
    return sorted(common)


def shared_genes(datasets: list[VersionDataset]) -> list[str]:
    common = set(datasets[0].genes)
    for dataset in datasets[1:]:
        common &= set(dataset.genes)
    return sorted(common, key=natural_key)


def panel_data(dataset: VersionDataset, panel: str,
               all_shared_features: list[str], all_shared_genes: list[str]):
    features = dataset.feature_names if panel in ("full", "matched_genes") else all_shared_features
    genes = dataset.genes if panel in ("full", "matched_features") else all_shared_genes
    X = dataset.features.reindex(index=genes, columns=features, fill_value=0).astype(float)
    y = dataset.labels.reindex(index=genes)
    valid = y.notna()
    return X.loc[valid], y.loc[valid].astype(int)


def pair_panel_data(train: VersionDataset, test: VersionDataset, panel: str,
                    all_shared_features: list[str], all_shared_genes: list[str]):
    if panel == "full":
        features = train.feature_names
        train_genes = train.genes
        test_genes = test.genes
    elif panel == "matched_features":
        features = all_shared_features
        train_genes = train.genes
        test_genes = test.genes
    elif panel == "matched_genes":
        features = train.feature_names
        train_genes = all_shared_genes
        test_genes = all_shared_genes
    elif panel == "matched_both":
        features = all_shared_features
        common = all_shared_genes
        train_genes = common
        test_genes = common
    else:
        raise ValueError(f"Unknown panel: {panel}")

    X_train = train.features.reindex(index=train_genes, columns=features, fill_value=0).astype(float)
    y_train = train.labels.reindex(index=train_genes)
    X_test = test.features.reindex(index=test_genes, columns=features, fill_value=0).astype(float)
    y_test = test.labels.reindex(index=test_genes)

    train_valid = y_train.notna()
    test_valid = y_test.notna()
    return (
        X_train.loc[train_valid],
        y_train.loc[train_valid].astype(int),
        X_test.loc[test_valid],
        y_test.loc[test_valid].astype(int),
    )


def can_fit(y) -> bool:
    return len(set(map(int, y))) >= 2


def cross_year_gene_folds(datasets, panel, all_shared_genes, args):
    """Return fixed gene-ID folds reused across every transfer direction."""
    from sklearn.model_selection import KFold

    if panel in {"matched_genes", "matched_both"}:
        gene_universe = list(all_shared_genes)
    else:
        gene_universe = sorted(
            {gene for dataset in datasets for gene in dataset.genes},
            key=natural_key,
        )
    if len(gene_universe) < args.cv_folds:
        return gene_universe, []

    records = []
    for repeat in range(1, args.cv_repeats + 1):
        splitter = KFold(
            n_splits=args.cv_folds,
            shuffle=True,
            random_state=args.seed + repeat,
        )
        for split, (_, heldout_idx) in enumerate(splitter.split(gene_universe), 1):
            records.append((
                repeat,
                split,
                {gene_universe[index] for index in heldout_idx},
            ))
    return gene_universe, records


def evaluation_preflight(datasets, panels, all_features, all_genes, args, models=None):
    """Verify that the requested evaluation design is complete before fitting."""
    import numpy as np
    from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold

    rows = []
    errors = []

    models = set(models or [])

    def check_training(y_train, context):
        if not can_fit(y_train):
            errors.append(f"{context}: training data contain only one class")
            return
        min_class = int(np.bincount(np.asarray(y_train, dtype=int)).min())
        if args.calibration != "none" and min_class < args.calibration_cv:
            errors.append(
                f"{context}: smallest training class has {min_class} genes, below "
                f"-calibration_cv {args.calibration_cv}"
            )
        if "nn" in models and getattr(args, "nn_early_stopping", True):
            n_inner = len(y_train)
            min_inner_class = min_class
            if args.calibration != "none":
                n_inner -= math.ceil(len(y_train) / args.calibration_cv)
                min_inner_class -= math.ceil(min_class / args.calibration_cv)
            validation_size = math.ceil(
                n_inner * getattr(args, "nn_validation_fraction", 0.1)
            )
            if min_inner_class < 2 or validation_size < 2 or n_inner - validation_size < 2:
                errors.append(
                    f"{context}: NN calibration/early-stopping splits are infeasible "
                    f"(inner n={n_inner}, smallest class={min_inner_class}, "
                    f"validation n={validation_size})"
                )

    for panel in panels:
        for dataset in datasets:
            X, y = panel_data(dataset, panel, all_features, all_genes)
            row = {
                "analysis": "within_year_cv",
                "panel": panel,
                "train_dataset": dataset.name,
                "test_dataset": dataset.name,
                "n_genes": len(y),
                "n_features": X.shape[1],
                "n_viable": int((y == 0).sum()),
                "n_lethal": int((y == 1).sum()),
                "requested_folds": args.cv_folds,
            }
            issue_count = len(errors)
            if X.shape[1] == 0:
                errors.append(f"{dataset.name}/{panel}: no GO features")
            elif len(y) < args.cv_folds:
                errors.append(
                    f"{dataset.name}/{panel}: {len(y)} genes cannot support "
                    f"{args.cv_folds} folds"
                )
            elif not can_fit(y):
                errors.append(f"{dataset.name}/{panel}: only one class is present")
            else:
                min_class = int(np.bincount(y).min())
                if panel not in {"matched_genes", "matched_both"} and min_class < args.cv_folds:
                    errors.append(
                        f"{dataset.name}/{panel}: smallest class has {min_class} genes, "
                        f"below -cv_folds {args.cv_folds}"
                    )
                else:
                    splitter = (
                        RepeatedKFold(
                            n_splits=args.cv_folds,
                            n_repeats=args.cv_repeats,
                            random_state=args.seed,
                        )
                        if panel in {"matched_genes", "matched_both"}
                        else RepeatedStratifiedKFold(
                            n_splits=args.cv_folds,
                            n_repeats=args.cv_repeats,
                            random_state=args.seed,
                        )
                    )
                    split_iter = splitter.split(X) if panel in {
                        "matched_genes", "matched_both"
                    } else splitter.split(X, y)
                    for fold, (train_idx, _) in enumerate(split_iter, 1):
                        check_training(
                            y.iloc[train_idx], f"{dataset.name}/{panel}/fold {fold}",
                        )
            row["status"] = "pass" if len(errors) == issue_count else "fail"
            rows.append(row)

    for panel in panels:
        gene_universe, split_records = cross_year_gene_folds(
            datasets, panel, all_genes, args,
        )
        for train in datasets:
            for test in datasets:
                X_train, y_train, X_test, y_test = pair_panel_data(
                    train, test, panel, all_features, all_genes,
                )
                row = {
                    "analysis": "cross_year_transfer",
                    "panel": panel,
                    "train_dataset": train.name,
                    "test_dataset": test.name,
                    "n_genes": len(y_test),
                    "n_features": X_test.shape[1],
                    "n_viable": int((y_test == 0).sum()),
                    "n_lethal": int((y_test == 1).sum()),
                    "requested_folds": args.cv_folds,
                    "fold_gene_universe": len(gene_universe),
                }
                issue_count = len(errors)
                context = f"{panel} {train.name}->{test.name}"
                if not split_records:
                    errors.append(
                        f"{context}: fixed gene universe cannot support "
                        f"{args.cv_folds} folds"
                    )
                elif not len(X_train) or not len(X_test) or X_train.shape[1] == 0:
                    errors.append(f"{context}: empty train/test cohort or feature set")
                elif not can_fit(y_test):
                    errors.append(f"{context}: test cohort contains only one class")
                else:
                    for repeat, split, heldout_genes in split_records:
                        train_mask = ~X_train.index.isin(heldout_genes)
                        check_training(
                            y_train.loc[train_mask],
                            f"{context}/repeat {repeat}/fold {split}",
                        )
                row["status"] = "pass" if len(errors) == issue_count else "fail"
                rows.append(row)

    if errors and not args.allow_incomplete_evaluation:
        examples = "; ".join(errors[:10])
        suffix = f"; plus {len(errors) - 10} more" if len(errors) > 10 else ""
        raise ValueError(
            "Requested publication evaluation is incomplete: " + examples + suffix +
            ". Reduce folds under a pre-specified protocol or pass "
            "-allow_incomplete_evaluation for explicitly exploratory output."
        )
    return rows, errors


def metric_value(value):
    if value is None:
        return float("nan")
    try:
        if math.isnan(float(value)):
            return float("nan")
    except Exception:
        pass
    return float(value)


def confidence_half_width(values, confidence=0.95):
    import numpy as np
    from scipy.stats import t

    clean = np.asarray([value for value in values if not math.isnan(float(value))], dtype=float)
    if len(clean) < 2:
        return 0.0
    standard_error = float(clean.std(ddof=1) / math.sqrt(len(clean)))
    return float(t.ppf((1 + confidence) / 2, df=len(clean) - 1) * standard_error)


def compute_metrics(y_true, y_score, threshold: float = 0.5) -> dict:
    import numpy as np
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        balanced_accuracy_score,
        brier_score_loss,
        f1_score,
        matthews_corrcoef,
        precision_score,
        recall_score,
        roc_auc_score,
    )

    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    y_pred = (y_score >= threshold).astype(int)
    one_class = len(np.unique(y_true)) < 2

    return {
        "n": int(len(y_true)),
        "n_lethal": int(np.sum(y_true == 1)),
        "n_viable": int(np.sum(y_true == 0)),
        "prevalence_lethal": float(np.mean(y_true == 1)) if len(y_true) else float("nan"),
        "accuracy": float(accuracy_score(y_true, y_pred)) if len(y_true) else float("nan"),
        "balanced_accuracy": (
            float(balanced_accuracy_score(y_true, y_pred)) if not one_class else float("nan")
        ),
        "precision_lethal": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall_lethal": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1_lethal": float(f1_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "f1_macro": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, y_pred)) if not one_class else float("nan"),
        "roc_auc": float(roc_auc_score(y_true, y_score)) if not one_class else float("nan"),
        "average_precision": (
            float(average_precision_score(y_true, y_score)) if not one_class else float("nan")
        ),
        "brier_score": float(brier_score_loss(y_true, y_score)) if len(y_true) else float("nan"),
    }


def make_options(args, seed: int):
    return SimpleNamespace(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        random_state=seed,
        n_jobs=args.n_jobs,
        calibration=getattr(args, "calibration", "none"),
        nn_hidden_units=tuple(getattr(args, "nn_hidden_units", (128, 64))),
        nn_alpha=getattr(args, "nn_alpha", 0.0001),
        nn_batch_size=getattr(args, "nn_batch_size", 32),
        nn_learning_rate_init=getattr(args, "nn_learning_rate_init", 0.001),
        nn_max_iter=getattr(args, "nn_max_iter", 300),
        nn_early_stopping=getattr(args, "nn_early_stopping", True),
        nn_validation_fraction=getattr(args, "nn_validation_fraction", 0.1),
        nn_n_iter_no_change=getattr(args, "nn_n_iter_no_change", 15),
    )


def fit_sklearn_model(model_type: str, X_train, y_train, args, seed: int,
                      calibrate: bool = True):
    import numpy as np
    import warnings
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.utils.class_weight import compute_sample_weight

    from .sklearn_models import build_sklearn_model

    base_model = build_sklearn_model(model_type, make_options(args, seed))
    model = base_model
    calibration = getattr(args, "calibration", "none")
    min_class = int(np.bincount(np.asarray(y_train, dtype=int)).min())
    if calibrate and calibration != "none" and min_class >= 2:
        calibration_cv = min(getattr(args, "calibration_cv", 3), min_class)
        model = CalibratedClassifierCV(
            estimator=base_model,
            method=calibration,
            cv=calibration_cv,
        )
    sample_weight = None
    if model_type in {"gb", "nn"}:
        sample_weight = compute_sample_weight(class_weight="balanced", y=y_train)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        if sample_weight is not None:
            model.fit(X_train, y_train, sample_weight=sample_weight)
        else:
            model.fit(X_train, y_train)
    convergence_warnings = [
        warning for warning in caught
        if issubclass(warning.category, ConvergenceWarning)
    ]
    estimators = [model]
    if hasattr(model, "calibrated_classifiers_"):
        estimators = [item.estimator for item in model.calibrated_classifiers_]
    iterations = [
        int(np.asarray(estimator.n_iter_).max())
        for estimator in estimators
        if hasattr(estimator, "n_iter_")
    ]
    model._phengo_fit_diagnostics = {
        "fit_convergence_warnings": len(convergence_warnings),
        "fit_n_iter_max": max(iterations) if iterations else None,
    }
    if convergence_warnings:
        logger.warning(
            "%s fit emitted %d convergence warning(s) for seed %d",
            model_type,
            len(convergence_warnings),
            seed,
        )
    return model


def model_fit_diagnostics(model) -> dict:
    return dict(getattr(model, "_phengo_fit_diagnostics", {}))


def model_scores(model, X):
    if hasattr(model, "predict_proba"):
        import numpy as np

        classes = list(model.classes_)
        if 1 not in classes:
            return np.zeros(len(X), dtype=float)
        return model.predict_proba(X)[:, classes.index(1)]
    if hasattr(model, "decision_function"):
        import numpy as np

        scores = model.decision_function(X)
        return 1 / (1 + np.exp(-scores))
    return model.predict(X)


def baseline_scores(kind: str, y_train, n_test: int, seed: int):
    import numpy as np

    prevalence = float(np.mean(y_train == 1)) if len(y_train) else 0.0
    if kind == "baseline_majority":
        return np.full(n_test, 1.0 if prevalence >= 0.5 else 0.0)
    if kind == "baseline_prior":
        return np.full(n_test, prevalence)
    if kind == "baseline_stratified_random":
        rng = np.random.default_rng(seed)
        return rng.binomial(1, prevalence, size=n_test).astype(float)
    raise ValueError(f"Unknown baseline: {kind}")


def dataset_drift_rows(datasets: list[VersionDataset]) -> list[dict]:
    rows = []
    for dataset in datasets:
        X = dataset.features
        y = dataset.labels
        active_features = [col for col in dataset.feature_names if float(X[col].sum()) > 0]
        go_per_gene = X.sum(axis=1)
        rows.append({
            "dataset": dataset.name,
            "path": dataset.path,
            "n_genes": int(X.shape[0]),
            "n_features": int(X.shape[1]),
            "n_active_features": int(len(active_features)),
            "n_lethal": int((y == 1).sum()),
            "n_viable": int((y == 0).sum()),
            "prevalence_lethal": float((y == 1).mean()) if len(y) else float("nan"),
            "mean_go_terms_per_gene": float(go_per_gene.mean()) if len(go_per_gene) else 0.0,
            "median_go_terms_per_gene": float(go_per_gene.median()) if len(go_per_gene) else 0.0,
            "sparsity": float((X.values == 0).sum() / X.size) if X.size else float("nan"),
        })
    return rows


def pairwise_drift_rows(datasets: list[VersionDataset]) -> list[dict]:
    rows = []
    for a in datasets:
        for b in datasets:
            genes_a, genes_b = set(a.genes), set(b.genes)
            feats_a, feats_b = set(a.feature_names), set(b.feature_names)
            common_genes = sorted(genes_a & genes_b, key=natural_key)
            label_churn = float("nan")
            if common_genes:
                label_churn = float((a.labels.loc[common_genes].values != b.labels.loc[common_genes].values).mean())
            rows.append({
                "dataset_a": a.name,
                "dataset_b": b.name,
                "n_common_genes": len(common_genes),
                "gene_jaccard": len(genes_a & genes_b) / len(genes_a | genes_b) if genes_a | genes_b else float("nan"),
                "n_common_features": len(feats_a & feats_b),
                "feature_jaccard": len(feats_a & feats_b) / len(feats_a | feats_b) if feats_a | feats_b else float("nan"),
                "label_churn_common_genes": label_churn,
            })
    return rows


def run_within_year_cv(datasets: list[VersionDataset], models: list[str], panels: list[str],
                       all_features: list[str], all_genes: list[str], args) -> tuple[list[dict], list[dict]]:
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold

    fold_rows = []
    repeat_predictions = defaultdict(
        lambda: {
            "y": [],
            "scores": [],
            "splits": set(),
            "fit_count": 0,
            "fit_convergence_warnings": 0,
            "fit_n_iter": [],
        }
    )
    for panel in panels:
        for dataset in datasets:
            X, y = panel_data(dataset, panel, all_features, all_genes)
            if X.shape[1] == 0:
                logger.warning(f"Skipping CV for {dataset.name}/{panel}: no GO features")
                continue
            if not can_fit(y):
                logger.warning(f"Skipping CV for {dataset.name}/{panel}: only one class present")
                continue
            min_class = int(np.bincount(y).min())
            n_splits = (
                min(args.cv_folds, len(y))
                if panel in {"matched_genes", "matched_both"}
                else min(args.cv_folds, min_class)
            )
            if n_splits < 2:
                logger.warning(f"Skipping CV for {dataset.name}/{panel}: fewer than 2 samples in one class")
                continue
            if panel in {"matched_genes", "matched_both"}:
                splitter = RepeatedKFold(
                    n_splits=n_splits, n_repeats=args.cv_repeats, random_state=args.seed,
                )
                split_iter = splitter.split(X)
            else:
                splitter = RepeatedStratifiedKFold(
                    n_splits=n_splits, n_repeats=args.cv_repeats, random_state=args.seed,
                )
                split_iter = splitter.split(X, y)
            for fold_idx, (train_idx, test_idx) in enumerate(split_iter, 1):
                repeat = (fold_idx - 1) // n_splits + 1
                split = (fold_idx - 1) % n_splits + 1
                if not can_fit(y.iloc[train_idx]):
                    logger.warning("Skipping class-deficient training fold for %s/%s", dataset.name, panel)
                    continue
                for model_type in models:
                    seed = args.seed + fold_idx
                    model = fit_sklearn_model(
                        model_type,
                        X.iloc[train_idx],
                        y.iloc[train_idx],
                        args,
                        seed,
                    )
                    scores = model_scores(model, X.iloc[test_idx])
                    metrics = compute_metrics(y.iloc[test_idx], scores)
                    diagnostics = model_fit_diagnostics(model)
                    fold_rows.append({
                        "analysis": "within_year_cv",
                        "panel": panel,
                        "dataset": dataset.name,
                        "model": model_type,
                        "fold": fold_idx,
                        "repeat": repeat,
                        "split": split,
                        **metrics,
                        **diagnostics,
                    })
                    key = (panel, dataset.name, model_type, repeat)
                    repeat_predictions[key]["y"].extend(y.iloc[test_idx].astype(int).tolist())
                    repeat_predictions[key]["scores"].extend(map(float, scores))
                    repeat_predictions[key]["splits"].add(split)
                    repeat_predictions[key]["fit_count"] += 1
                    repeat_predictions[key]["fit_convergence_warnings"] += int(
                        diagnostics.get("fit_convergence_warnings", 0)
                    )
                    if diagnostics.get("fit_n_iter_max") is not None:
                        repeat_predictions[key]["fit_n_iter"].append(
                            int(diagnostics["fit_n_iter_max"])
                        )

                for baseline in ["baseline_majority", "baseline_prior", "baseline_stratified_random"]:
                    scores = baseline_scores(baseline, y.iloc[train_idx], len(test_idx), args.seed + fold_idx)
                    metrics = compute_metrics(y.iloc[test_idx], scores)
                    fold_rows.append({
                        "analysis": "within_year_cv",
                        "panel": panel,
                        "dataset": dataset.name,
                        "model": baseline,
                        "fold": fold_idx,
                        "repeat": repeat,
                        "split": split,
                        **metrics,
                    })
                    key = (panel, dataset.name, baseline, repeat)
                    repeat_predictions[key]["y"].extend(y.iloc[test_idx].astype(int).tolist())
                    repeat_predictions[key]["scores"].extend(map(float, scores))
                    repeat_predictions[key]["splits"].add(split)

    if not fold_rows:
        return [], []

    repeat_rows = []
    for (panel, dataset, model, repeat), values in repeat_predictions.items():
        repeat_rows.append({
            "panel": panel,
            "dataset": dataset,
            "model": model,
            "repeat": repeat,
            "folds_in_repeat": len(values["splits"]),
            "model_fits": values["fit_count"],
            "fit_convergence_warnings": values["fit_convergence_warnings"],
            "fit_n_iter_max": (
                max(values["fit_n_iter"]) if values["fit_n_iter"] else None
            ),
            **compute_metrics(values["y"], values["scores"]),
        })
    metric_cols = [
        column for column in repeat_rows[0]
        if column not in {"panel", "dataset", "model", "repeat", "folds_in_repeat"}
    ]
    df = pd.DataFrame(repeat_rows)
    summary = []
    for keys, group in df.groupby(["panel", "dataset", "model"], dropna=False):
        row = {
            "panel": keys[0],
            "dataset": keys[1],
            "model": keys[2],
            "n_folds": int(group["folds_in_repeat"].sum()),
            "n_repeats": int(group["repeat"].nunique()),
        }
        for metric in metric_cols:
            repeat_values = group[metric].dropna()
            row[f"{metric}_mean"] = float(repeat_values.mean()) if len(repeat_values) else float("nan")
            row[f"{metric}_sd"] = (
                float(repeat_values.std(ddof=1)) if len(repeat_values) > 1 else 0.0
            )
            row[f"{metric}_ci95"] = confidence_half_width(repeat_values)
        summary.append(row)
    return fold_rows, summary


def _summarise_repeated_rows(rows, group_columns):
    import pandas as pd

    if not rows:
        return []
    df = pd.DataFrame(rows)
    excluded = set(group_columns) | {"analysis", "repeat", "evaluation", "n_overlap_genes"}
    metrics = [column for column in df.columns if column not in excluded]
    summaries = []
    for keys, group in df.groupby(group_columns, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_columns, keys))
        row.update({
            "analysis": "cross_year_transfer_summary",
            "evaluation": "fixed_gene_disjoint_oof",
            "n_repeats": int(group["repeat"].nunique()),
            "n_overlap_genes": int(group["n_overlap_genes"].max()),
        })
        for metric in metrics:
            values = group[metric].dropna()
            row[f"{metric}_mean"] = float(values.mean()) if len(values) else float("nan")
            row[f"{metric}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
            row[f"{metric}_ci95"] = confidence_half_width(values)
        summaries.append(row)
    return summaries


def run_cross_year_transfer(datasets: list[VersionDataset], models: list[str], panels: list[str],
                            all_features: list[str], all_genes: list[str], args):
    """Gene-disjoint transfer using fixed folds shared by all test releases."""
    import numpy as np

    rows = []
    model_names = models + ["baseline_majority", "baseline_prior", "baseline_stratified_random"]
    for panel in panels:
        _, split_records = cross_year_gene_folds(datasets, panel, all_genes, args)
        if not split_records:
            logger.warning("Skipping %s transfer: fixed gene folds are unavailable", panel)
            continue
        for train in datasets:
            X_train, y_train = panel_data(train, panel, all_features, all_genes)
            if not len(X_train) or X_train.shape[1] == 0:
                logger.warning("Skipping %s transfer from %s: empty training data", panel, train.name)
                continue
            test_payloads = {}
            for test in datasets:
                _, _, X_test, y_test = pair_panel_data(
                    train, test, panel, all_features, all_genes,
                )
                if len(X_test):
                    test_payloads[test.name] = (test, X_test, y_test)
            for model_name in model_names:
                for repeat in range(1, args.cv_repeats + 1):
                    score_tables = {
                        name: np.full(len(payload[1]), np.nan, dtype=float)
                        for name, payload in test_payloads.items()
                    }
                    fit_count = 0
                    fit_convergence_warnings = 0
                    fit_iterations = []
                    repeat_records = (
                        record for record in split_records if record[0] == repeat
                    )
                    for _, split, heldout_genes in repeat_records:
                        train_mask = ~X_train.index.isin(heldout_genes)
                        fold_X_train = X_train.loc[train_mask]
                        fold_y_train = y_train.loc[train_mask]
                        if not can_fit(fold_y_train):
                            continue
                        seed = args.seed + repeat * 1000 + split
                        model = None
                        if not model_name.startswith("baseline_"):
                            model = fit_sklearn_model(
                                model_name, fold_X_train, fold_y_train, args, seed,
                            )
                            diagnostics = model_fit_diagnostics(model)
                            fit_count += 1
                            fit_convergence_warnings += int(
                                diagnostics.get("fit_convergence_warnings", 0)
                            )
                            if diagnostics.get("fit_n_iter_max") is not None:
                                fit_iterations.append(
                                    int(diagnostics["fit_n_iter_max"])
                                )
                        for test_name, (_, X_test, _) in test_payloads.items():
                            heldout_mask = X_test.index.isin(heldout_genes)
                            if not heldout_mask.any():
                                continue
                            if model is None:
                                fold_scores = baseline_scores(
                                    model_name, fold_y_train, int(heldout_mask.sum()), seed,
                                )
                            else:
                                fold_scores = model_scores(model, X_test.loc[heldout_mask])
                            score_tables[test_name][heldout_mask] = fold_scores

                    for test_name, (_, X_test, y_test) in test_payloads.items():
                        scores = score_tables[test_name]
                        valid = ~np.isnan(scores)
                        if not valid.any():
                            continue
                        metrics = compute_metrics(y_test.iloc[valid], scores[valid])
                        rows.append({
                            "analysis": "cross_year_transfer",
                            "evaluation": "fixed_gene_disjoint_oof",
                            "panel": panel,
                            "train_dataset": train.name,
                            "test_dataset": test_name,
                            "model": model_name,
                            "repeat": repeat,
                            "n_overlap_genes": len(
                                set(X_train.index) & set(X_test.index)
                            ),
                            "model_fits": fit_count,
                            "fit_convergence_warnings": fit_convergence_warnings,
                            "fit_n_iter_max": max(fit_iterations) if fit_iterations else None,
                            **metrics,
                        })
    summary = _summarise_repeated_rows(
        rows, ["panel", "train_dataset", "test_dataset", "model"],
    )
    add_calendar_interval_fields(rows)
    add_calendar_interval_fields(summary)
    return rows, summary


def run_previous_available_snapshot_label_baseline(
    datasets: list[VersionDataset],
) -> list[dict]:
    rows = []
    ordered = sorted(datasets, key=lambda ds: natural_key(ds.name))
    for prev, current in zip(ordered, ordered[1:]):
        common = sorted(set(prev.genes) & set(current.genes), key=natural_key)
        if not common:
            continue
        scores = prev.labels.loc[common].astype(float).values
        metrics = compute_metrics(current.labels.loc[common], scores)
        rows.append({
            "analysis": "previous_available_snapshot_label_baseline",
            "train_dataset": prev.name,
            "test_dataset": current.name,
            "model": "previous_available_snapshot_labels",
            **metrics,
        })
    add_calendar_interval_fields(rows)
    return rows


def run_previous_year_label_baseline(datasets: list[VersionDataset]) -> list[dict]:
    """Compatibility alias; rows explicitly describe previous available snapshots."""
    return run_previous_available_snapshot_label_baseline(datasets)


def native_importance(model, model_type: str, feature_names: list[str]) -> dict[str, float] | None:
    import numpy as np

    if model_type in {"dt", "rf", "gb"} and hasattr(model, "feature_importances_"):
        values = model.feature_importances_
    elif model_type == "lr" and hasattr(model, "coef_"):
        values = np.abs(model.coef_[0])
    else:
        return None
    return {feature: float(value) for feature, value in zip(feature_names, values)}


def run_feature_rank_overlap(datasets: list[VersionDataset], models: list[str],
                             all_features: list[str], all_genes: list[str], args):
    import numpy as np

    overlap_rows = []
    stability_rows = []
    supported = [m for m in models if m in {"dt", "rf", "gb", "lr"}]
    if not supported or not all_features:
        return overlap_rows, stability_rows

    top_sets = {}
    for dataset in datasets:
        X, y = panel_data(dataset, "matched_features", all_features, dataset.genes)
        if not can_fit(y):
            continue
        for model_type in supported:
            importance_samples = {feature: [] for feature in X.columns}
            top_counts = {feature: 0 for feature in X.columns}
            for repeat in range(args.importance_repeats):
                rng = np.random.default_rng(args.seed + repeat)
                bootstrap_parts = [
                    rng.choice(np.flatnonzero(y.to_numpy() == label), size=int((y == label).sum()), replace=True)
                    for label in (0, 1)
                ]
                bootstrap_indices = np.concatenate(bootstrap_parts)
                rng.shuffle(bootstrap_indices)
                model = fit_sklearn_model(
                    model_type,
                    X.iloc[bootstrap_indices],
                    y.iloc[bootstrap_indices],
                    args,
                    args.seed + repeat,
                    calibrate=False,
                )
                importance = native_importance(model, model_type, list(X.columns))
                if not importance:
                    continue
                ranked = sorted(
                    importance.items(), key=lambda item: item[1], reverse=True,
                )
                for feature, value in importance.items():
                    importance_samples[feature].append(value)
                for feature, _ in ranked[:args.top_k]:
                    top_counts[feature] += 1

            mean_importance = {
                feature: float(np.mean(values)) if values else float("nan")
                for feature, values in importance_samples.items()
            }
            ranked_mean = sorted(
                mean_importance.items(),
                key=lambda item: (-item[1] if not math.isnan(item[1]) else float("inf"), item[0]),
            )
            top_sets[(dataset.name, model_type)] = {
                feature for feature, _ in ranked_mean[:args.top_k]
            }
            for rank, (feature, mean_value) in enumerate(ranked_mean, 1):
                values = importance_samples[feature]
                stability_rows.append({
                    "dataset": dataset.name,
                    "model": model_type,
                    "GO_Term": feature,
                    "mean_importance": mean_value,
                    "sd_importance": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
                    "top_k_selection_frequency": top_counts[feature] / max(len(values), 1),
                    "consensus_rank": rank,
                    "bootstrap_repeats": len(values),
                })

    for model_type in supported:
        for a in datasets:
            for b in datasets:
                set_a = top_sets.get((a.name, model_type), set())
                set_b = top_sets.get((b.name, model_type), set())
                overlap_rows.append({
                    "model": model_type,
                    "dataset_a": a.name,
                    "dataset_b": b.name,
                    "top_k": args.top_k,
                    "n_overlap": len(set_a & set_b),
                    "jaccard": len(set_a & set_b) / len(set_a | set_b) if set_a | set_b else float("nan"),
                })
    return overlap_rows, stability_rows


def run_prediction_instability(datasets: list[VersionDataset], models: list[str],
                               all_features: list[str], all_genes: list[str], args) -> tuple[list[dict], list[dict]]:
    import numpy as np
    import pandas as pd
    from sklearn.model_selection import RepeatedKFold

    if not all_features or not all_genes:
        return [], []

    per_gene_rows = []
    summary_rows = []
    ordered = sorted(datasets, key=lambda ds: natural_key(ds.name))
    n_splits = min(args.cv_folds, len(all_genes))
    if n_splits < 2:
        return [], []
    splitter = RepeatedKFold(
        n_splits=n_splits,
        n_repeats=args.cv_repeats,
        random_state=args.seed,
    )
    # The same gene-index folds are reused for every snapshot and model. A gene
    # is removed from training whenever its score is generated.
    split_records = list(splitter.split(np.arange(len(all_genes))))
    for model_type in models:
        score_table = {}
        label_table = {}
        pred_table = {}
        for dataset in ordered:
            X, y = panel_data(dataset, "matched_both", all_features, all_genes)
            if not can_fit(y):
                continue
            score_sum = np.zeros(len(X), dtype=float)
            score_count = np.zeros(len(X), dtype=int)
            for fold_idx, (train_idx, heldout_idx) in enumerate(split_records, 1):
                fold_y = y.iloc[train_idx]
                if not can_fit(fold_y):
                    continue
                model = fit_sklearn_model(
                    model_type,
                    X.iloc[train_idx],
                    fold_y,
                    args,
                    args.seed + fold_idx,
                )
                heldout_scores = model_scores(model, X.iloc[heldout_idx])
                score_sum[heldout_idx] += heldout_scores
                score_count[heldout_idx] += 1
            if not np.all(score_count):
                logger.warning(
                    "Skipping instability scores for %s/%s: some fixed folds lacked both classes",
                    dataset.name,
                    model_type,
                )
                continue
            scores = score_sum / score_count
            score_table[dataset.name] = pd.Series(scores, index=X.index)
            label_table[dataset.name] = y
            pred_table[dataset.name] = pd.Series((scores >= 0.5).astype(int), index=X.index)

        if len(score_table) < 2:
            continue

        scores_df = pd.DataFrame(score_table).loc[all_genes]
        labels_df = pd.DataFrame(label_table).loc[all_genes]
        preds_df = pd.DataFrame(pred_table).loc[all_genes]
        for gene in all_genes:
            score_values = scores_df.loc[gene].astype(float)
            pred_values = preds_df.loc[gene].astype(int)
            label_values = labels_df.loc[gene].astype(int)
            per_gene_rows.append({
                "analysis": "prediction_instability",
                "evaluation": "fixed_gene_disjoint_oof",
                "model": model_type,
                "gene": gene,
                "n_snapshots": int(len(score_values)),
                "n_cv_repeats": int(args.cv_repeats),
                "mean_probability_lethal": float(score_values.mean()),
                "sd_probability_lethal": float(score_values.std(ddof=1)) if len(score_values) > 1 else 0.0,
                "range_probability_lethal": float(score_values.max() - score_values.min()),
                "prediction_flip": int(len(set(pred_values)) > 1),
                "label_flip": int(len(set(label_values)) > 1),
                "min_probability_lethal": float(score_values.min()),
                "max_probability_lethal": float(score_values.max()),
            })

        summary_rows.append({
            "analysis": "prediction_instability_summary",
            "evaluation": "fixed_gene_disjoint_oof",
            "model": model_type,
            "panel": "matched_both",
            "n_genes": len(all_genes),
            "n_snapshots": len(score_table),
            "n_cv_repeats": int(args.cv_repeats),
            "mean_gene_probability_sd": float(np.nanmean(scores_df.std(axis=1, ddof=1))),
            "median_gene_probability_sd": float(np.nanmedian(scores_df.std(axis=1, ddof=1))),
            "mean_gene_probability_range": float(np.nanmean(scores_df.max(axis=1) - scores_df.min(axis=1))),
            "prediction_flip_rate": float((preds_df.nunique(axis=1) > 1).mean()),
            "label_flip_rate": float((labels_df.nunique(axis=1) > 1).mean()),
        })
    return per_gene_rows, summary_rows


def write_csv(path: str, rows: list[dict]):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    preferred = [
        "analysis", "panel", "dataset", "train_dataset", "test_dataset", "model",
        "fold", "gene", "metric",
    ]
    fieldnames = [f for f in preferred if f in fieldnames] + [
        f for f in fieldnames if f not in preferred
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_transfer_matrices(output_dir: str, transfer_rows: list[dict],
                            dataset_names: list[str], metrics: list[str] = MATRIX_METRICS):
    import pandas as pd

    if not transfer_rows:
        return
    matrix_dir = os.path.join(output_dir, "transfer_matrices")
    os.makedirs(matrix_dir, exist_ok=True)
    df = pd.DataFrame(transfer_rows)
    for panel in sorted(df["panel"].dropna().unique(), key=natural_key):
        for model in sorted(df["model"].dropna().unique(), key=natural_key):
            subset = df[(df["panel"] == panel) & (df["model"] == model)]
            for metric in metrics:
                value_column = metric if metric in subset.columns else f"{metric}_mean"
                if value_column not in subset.columns:
                    continue
                pivot = subset.pivot(
                    index="train_dataset", columns="test_dataset", values=value_column,
                )
                pivot = pivot.reindex(index=dataset_names, columns=dataset_names)
                out = os.path.join(matrix_dir, f"{panel}_{model}_{metric}.csv")
                pivot.to_csv(out)


def write_markdown_report(output_dir: str, args, datasets: list[VersionDataset],
                          models: list[str], panels: list[str]):
    report_path = os.path.join(output_dir, "VERSION_SENSITIVITY_REPORT.md")
    lines = [
        "# PhenGO Version-Sensitivity Analysis Report",
        "",
        f"Generated: {datetime.now().isoformat()}",
        "",
        "## Purpose",
        "",
        "This analysis is designed to test whether gene-essentiality ML outputs are stable "
        "across yearly snapshots of model-organism database resources. It separates "
        "performance changes caused by gene coverage, phenotype-label changes, GO feature "
        "churn, and model behavior.",
        "",
        "## Inputs",
        "",
    ]
    for dataset in datasets:
        lines.append(f"- `{dataset.name}`: `{dataset.path}`")
    lines.extend([
        "",
        "## Configuration",
        "",
        f"- Models: {', '.join(models)}",
        f"- Panels: {', '.join(panels)}",
        f"- CV folds: {args.cv_folds}",
        f"- CV repeats: {args.cv_repeats}",
        f"- Probability calibration: {args.calibration}",
        f"- Calibration folds: {args.calibration_cv}",
        "- NN backend: scikit-learn MLPClassifier",
        f"- NN hidden units: {', '.join(map(str, args.nn_hidden_units))}",
        f"- NN maximum iterations: {args.nn_max_iter}",
        f"- NN early stopping: {args.nn_early_stopping}",
        f"- Random seed: {args.seed}",
        f"- Strict manifests required: {not args.allow_missing_manifests}",
        f"- Top-k feature overlap: {args.top_k}",
        f"- Feature-importance bootstrap repeats: {args.importance_repeats}",
        "",
        "## Output Files",
        "",
        "- `evaluation_preflight.csv`: confirms requested folds and calibration are feasible for every panel and transfer cell.",
        "- `dataset_drift.csv`: per-year gene counts, GO feature counts, sparsity, and class prevalence.",
        "- `pairwise_drift.csv`: pairwise gene-set Jaccard, GO-feature Jaccard, and label churn among shared genes.",
        "- `within_year_cv.csv`: repeated stratified CV fold-level metrics for each dataset, model, and panel.",
        "- `within_year_cv_summary.csv`: repeat-level mean, standard deviation, and 95% interval half-width for CV metrics.",
        "- `cross_year_transfer.csv`: repeat-level, gene-disjoint train-year -> test-year metrics.",
        "- `cross_year_transfer_summary.csv`: transfer means, standard deviations, and 95% interval half-widths across repeats.",
        "- `transfer_matrices/`: matrix-form CSVs for key metrics, ready for heatmaps.",
        "- `previous_available_snapshot_label_baseline.csv`: adjacent available-snapshot baseline using the earlier snapshot's labels for shared genes, with the calendar gap reported.",
        "- `prediction_instability.csv`: per-gene probability variance and prediction flips across years.",
        "- `prediction_instability_summary.csv`: model-level summary of prediction instability.",
        "- `feature_rank_overlap.csv`: pairwise top-k GO importance overlap for models with native importances.",
        "- `feature_importance_stability.csv`: bootstrap mean, variability, and top-k selection frequency for each GO feature.",
        "",
        "## Panels",
        "",
        "- `full`: each transfer uses the training release's GO feature vocabulary; GO terms first seen in the test release are ignored.",
        "- `matched_features`: all datasets are restricted to GO features shared by every dataset.",
        "- `matched_genes`: every dataset is restricted to genes shared by all snapshots and each transfer uses its training release's features.",
        "- `matched_both`: all datasets are restricted to genes and GO features shared by every dataset.",
        "",
        "## Recommended Article Use",
        "",
        "1. Use `dataset_drift.csv` and `pairwise_drift.csv` to establish that the source resources themselves change.",
        "2. Use `within_year_cv_summary.csv` to show that ordinary within-snapshot ML performance varies by year.",
        "3. Use `cross_year_transfer_summary.csv` and `transfer_matrices/` as the central evidence: models trained on one snapshot "
        "do not necessarily transfer cleanly to another snapshot.",
        "4. Compare `full`, `matched_features`, `matched_genes`, and `matched_both` panels to separate feature churn, gene churn, "
        "and label/annotation changes.",
        "5. Use `prediction_instability_summary.csv` and `prediction_instability.csv` to show whether the same genes receive "
        "different ML outputs depending on the database year.",
        "6. Use `feature_rank_overlap.csv` to show whether the biological signals apparently learned by the model are stable.",
        "",
        "## Suggested Manuscript Framing",
        "",
        "The models should be described as probes of database-version sensitivity rather than as final production predictors. "
        "If instability appears across simple and complex model families, the most defensible interpretation is that the "
        "changing resource snapshot materially affects ML conclusions.",
        "",
        "All transfer and per-gene instability predictions use fixed gene-ID folds and are gene-disjoint and out-of-fold. "
        "Probability-producing classifiers are calibrated within each training fold unless calibration is disabled explicitly.",
        "The `nn` model is a scikit-learn MLPClassifier evaluated by the same folds, panels, transfer directions, "
        "class-balancing policy, calibration, threshold, seeds, and metrics as the other model families. "
        "Inspect the fit-convergence columns before interpreting neural-network differences.",
        "",
    ])
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        prog="phengo-version-sensitivity",
        description="Analyse ML sensitivity to model-organism database version/year changes.",
    )
    parser.add_argument("-input_dir", help="Parent directory containing PhenGO output subdirectories or ARFF files")
    parser.add_argument("-arff_files", nargs="+", help="Explicit ARFF files in version/time order")
    parser.add_argument("-dataset_names", nargs="+", help="Names matching -arff_files")
    parser.add_argument("-output_dir", required=True)
    parser.add_argument("-models", nargs="+", default=DEFAULT_MODELS,
                        help=f"Models to run: {', '.join(SKLEARN_MODELS)} or all")
    parser.add_argument("-panels", nargs="+", default=PANELS, choices=PANELS)
    parser.add_argument("-cv_folds", type=int, default=5)
    parser.add_argument("-cv_repeats", type=int, default=5)
    parser.add_argument("-calibration", choices=["none", "sigmoid", "isotonic"], default="sigmoid",
                        help="Nested training-fold probability calibration (default: sigmoid)")
    parser.add_argument("-calibration_cv", type=int, default=3)
    parser.add_argument("-n_estimators", type=int, default=200)
    parser.add_argument("-max_depth", type=int, default=None)
    parser.add_argument("-n_jobs", type=int, default=1,
                        help="Parallel workers inside supported estimators (default: 1)")
    nn_group = parser.add_argument_group(
        "Comparable neural-network parameters",
        "Settings for sklearn MLPClassifier in the repeated version-sensitivity design.",
    )
    nn_group.add_argument("-nn_hidden_units", nargs="+", type=int, default=[128, 64])
    nn_group.add_argument("-nn_alpha", type=float, default=0.0001,
                          help="L2 regularization strength (default: 0.0001)")
    nn_group.add_argument("-nn_batch_size", type=int, default=32)
    nn_group.add_argument("-nn_learning_rate_init", type=float, default=0.001)
    nn_group.add_argument("-nn_max_iter", type=int, default=300)
    early_stopping_group = nn_group.add_mutually_exclusive_group()
    early_stopping_group.add_argument(
        "-nn_early_stopping",
        dest="nn_early_stopping",
        action="store_true",
        default=True,
        help="Use a training-only validation fraction for early stopping (default)",
    )
    early_stopping_group.add_argument(
        "-no_nn_early_stopping",
        dest="nn_early_stopping",
        action="store_false",
        help="Disable MLP early stopping and fit for at most -nn_max_iter iterations",
    )
    nn_group.add_argument("-nn_validation_fraction", type=float, default=0.1)
    nn_group.add_argument("-nn_n_iter_no_change", type=int, default=15)
    parser.add_argument("-top_k", type=int, default=20)
    parser.add_argument("-importance_repeats", type=int, default=20)
    parser.add_argument("-seed", type=int, default=42)
    parser.add_argument("-allow_missing_manifests", action="store_true",
                        help="Allow exploratory use of legacy ARFFs without strict snapshot manifests")
    parser.add_argument(
        "-allow_incomplete_evaluation", action="store_true",
        help=(
            "Permit adaptive/skipped folds or cells when requested CV/calibration is "
            "infeasible; exploratory only"
        ),
    )
    parser.add_argument("-overwrite", action="store_true")
    args = parser.parse_args()

    if args.input_dir:
        arff_files, dataset_names = discover_arff_files(args.input_dir)
    elif args.arff_files:
        arff_files = args.arff_files
        dataset_names = args.dataset_names or [
            os.path.splitext(os.path.basename(path))[0] for path in arff_files
        ]
    else:
        parser.error("Provide either -input_dir or -arff_files")

    if len(arff_files) < 2:
        parser.error("At least two ARFF datasets are required for version-sensitivity analysis")
    if len(dataset_names) != len(arff_files):
        parser.error("Number of -dataset_names values must match number of -arff_files values")
    if len(set(dataset_names)) != len(dataset_names):
        parser.error("Dataset names must be unique")
    missing_files = [path for path in arff_files if not os.path.isfile(path)]
    if missing_files:
        parser.error("ARFF file(s) not found: " + ", ".join(missing_files))
    if args.cv_folds < 2:
        parser.error("-cv_folds must be at least 2")
    if args.cv_repeats < 2:
        parser.error("-cv_repeats must be at least 2 for uncertainty estimates")
    if args.calibration_cv < 2:
        parser.error("-calibration_cv must be at least 2")
    if args.top_k < 1 or args.importance_repeats < 2:
        parser.error("-top_k must be at least 1 and -importance_repeats at least 2")
    if args.n_jobs == 0:
        parser.error("-n_jobs cannot be 0")
    if not args.nn_hidden_units or any(units < 1 for units in args.nn_hidden_units):
        parser.error("-nn_hidden_units values must be positive integers")
    if args.nn_alpha < 0:
        parser.error("-nn_alpha cannot be negative")
    if args.nn_batch_size < 1 or args.nn_max_iter < 1 or args.nn_n_iter_no_change < 1:
        parser.error("NN batch size, max iterations, and patience must be positive")
    if args.nn_learning_rate_init <= 0:
        parser.error("-nn_learning_rate_init must be positive")
    if not 0 < args.nn_validation_fraction < 1:
        parser.error("-nn_validation_fraction must be between 0 and 1")

    models = expand_models(args.models)
    try:
        datasets = [
            load_version_dataset(path, name)
            for path, name in zip(arff_files, dataset_names)
        ]
        datasets = sorted(datasets, key=lambda ds: natural_key(ds.name))
        dataset_names = [dataset.name for dataset in datasets]
        validate_dataset_manifests(datasets, args.allow_missing_manifests)
        all_features = shared_features(datasets)
        all_genes = shared_genes(datasets)
        preflight_rows, preflight_errors = evaluation_preflight(
            datasets, args.panels, all_features, all_genes, args, models,
        )
    except ValueError as exc:
        parser.error(str(exc))
    args.output_dir = ensure_output_dir(
        args.output_dir, args.overwrite, protected_paths=arff_files,
    )
    configure_logger(
        "PhenGO.version_sensitivity",
        enable_file=True,
        log_dir=args.output_dir,
        logfile_name="version_sensitivity.log",
    )
    logger.info("Loaded and provenance-validated %d datasets", len(datasets))

    logger.info(f"Datasets: {', '.join(dataset_names)}")
    logger.info(f"Models: {', '.join(models)}")
    logger.info(f"Panels: {', '.join(args.panels)}")
    logger.info(f"Shared genes across all datasets: {len(all_genes)}")
    logger.info(f"Shared GO features across all datasets: {len(all_features)}")
    if preflight_errors:
        logger.warning(
            "Exploratory evaluation preflight retained %d issue(s)",
            len(preflight_errors),
        )

    write_csv(os.path.join(args.output_dir, "evaluation_preflight.csv"), preflight_rows)
    write_csv(os.path.join(args.output_dir, "dataset_drift.csv"), dataset_drift_rows(datasets))
    write_csv(os.path.join(args.output_dir, "pairwise_drift.csv"), pairwise_drift_rows(datasets))

    logger.info("Running within-year repeated CV")
    cv_rows, cv_summary = run_within_year_cv(
        datasets, models, args.panels, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "within_year_cv.csv"), cv_rows)
    write_csv(os.path.join(args.output_dir, "within_year_cv_summary.csv"), cv_summary)

    logger.info("Running cross-year transfer analysis")
    transfer_rows, transfer_summary = run_cross_year_transfer(
        datasets, models, args.panels, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "cross_year_transfer.csv"), transfer_rows)
    write_csv(
        os.path.join(args.output_dir, "cross_year_transfer_summary.csv"),
        transfer_summary,
    )
    write_transfer_matrices(args.output_dir, transfer_summary, dataset_names)

    logger.info("Running previous-available-snapshot label baseline")
    write_csv(
        os.path.join(args.output_dir, "previous_available_snapshot_label_baseline.csv"),
        run_previous_available_snapshot_label_baseline(datasets),
    )

    logger.info("Running matched-gene prediction instability analysis")
    instability_rows, instability_summary = run_prediction_instability(
        datasets, models, all_features, all_genes, args,
    )
    write_csv(os.path.join(args.output_dir, "prediction_instability.csv"), instability_rows)
    write_csv(os.path.join(args.output_dir, "prediction_instability_summary.csv"), instability_summary)

    logger.info("Running native feature-rank overlap analysis")
    feature_overlap, feature_stability = run_feature_rank_overlap(
        datasets, models, all_features, all_genes, args,
    )
    write_csv(
        os.path.join(args.output_dir, "feature_rank_overlap.csv"),
        feature_overlap,
    )
    write_csv(
        os.path.join(args.output_dir, "feature_importance_stability.csv"),
        feature_stability,
    )

    write_markdown_report(args.output_dir, args, datasets, models, args.panels)
    repo_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
    repository_status = git_status(repo_dir)
    output_hashes = {}
    for root, _, filenames in os.walk(args.output_dir):
        for filename in filenames:
            if filename.endswith(".log") or filename == "version_sensitivity_manifest.json":
                continue
            path = os.path.join(root, filename)
            output_hashes[os.path.relpath(path, args.output_dir)] = sha256_file(path)
    analysis_manifest = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "PhenGO version sensitivity",
        "tool_version": PhenGO_VERSION,
        "git_commit": git_commit(repo_dir),
        "git_dirty": bool(repository_status),
        "git_status": repository_status,
        "source_tree_sha256": source_tree_sha256(repo_dir),
        "evaluation": "gene_disjoint_repeated_oof",
        "command": list(os.sys.argv),
        "dependencies": dependency_versions(),
        "configuration": {
            key: value for key, value in vars(args).items()
            if key not in {"arff_files"}
        },
        "datasets": [
            {
                "name": dataset.name,
                "path": os.path.abspath(dataset.path),
                "sha256": sha256_file(dataset.path),
                "snapshot_id": (dataset.manifest or {}).get("snapshot_id"),
                "source_manifest": bool(dataset.manifest),
            }
            for dataset in datasets
        ],
        "outputs": output_hashes,
    }
    with open(
        os.path.join(args.output_dir, "version_sensitivity_manifest.json"),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(analysis_manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    logger.info(f"Version-sensitivity analysis complete: {args.output_dir}")


if __name__ == "__main__":
    main()
