# PhenGO Changelog Since Last Git Commit

Base commit:

```text
d1707db Added additional ML models to Predict, better IO, user options, += data
```

This document summarises the uncommitted source, documentation, and test changes
made after the commit above. The changes are aimed at two related goals:

1. Correct several toolkit issues that could change biological or ML
   interpretation.
2. Add a new analysis workflow for testing whether ML conclusions vary with the
   year/version of model-organism database resources used to generate PhenGO
   ARFF files.

## Executive Summary

The largest conceptual change is that PhenGO now has a dedicated
version-sensitivity analysis command:

```bash
phengo-version-sensitivity
```

This command accepts multiple ARFF files representing different years/releases
for the same organism and generates within-year CV results, train-year/test-year
transfer matrices, source-dataset drift summaries, controlled matched-gene and
matched-feature panels, prediction-instability tables, feature-rank overlap
tables, and a generated Markdown report.

The most important correctness fix is the phenotype label mapping in
PhenGO-Predict. Previously, sklearn's alphabetical `LabelEncoder` mapped
`lethal` to `0` and `viable` to `1`, while the prediction code treated class `1`
as the positive/essential class. The new `PhenotypeLabelEncoder` makes the
mapping explicit:

```text
viable -> 0
lethal -> 1
```

Several ARFF and GO handling fixes were also added so identifiers containing
commas are parsed safely, GO ancestor terms are preserved consistently, mouse
MGI reports are joined more reliably, and output directories are no longer
silently deleted unless `-overwrite` is supplied.

## Version And Packaging Changes

Files changed:

- `setup.cfg`
- `README.md`
- `src/PhenGO/constants.py`

Changes:

- Package version changed from `0.1.2` to `0.2.0`.
- `python_requires` changed from `>=3.6` to `>=3.10`.
- Optional dependency groups were added:
  - `predict`: NumPy, pandas, scikit-learn, matplotlib.
  - `nn`: TensorFlow.
  - `reports`: NumPy, pandas, scikit-learn, matplotlib, seaborn.
  - `all`: all optional dependencies.
- New console scripts were added:
  - `phengo-version-sensitivity`
  - `version-sensitivity`
- README installation guidance now documents optional extras:

```bash
pip install "phengo[predict,reports]"
```

- README usage text now reflects `PhenGO v0.2.0`.
- README now documents the new `-overwrite` option for the core PhenGO command.
- README now includes a Version-Sensitivity Analysis section with example usage.
- Default bundled data paths in `constants.py` now use `importlib.resources`
  instead of fragile relative paths such as `../data/...`.

Why this matters:

- Installed packages can now find bundled default data files reliably.
- Users can install only the dependencies they need.
- TensorFlow is no longer implicitly required for users who only need sklearn
  models or reporting utilities.

## Core PhenGO Pipeline Changes

Files changed:

- `src/PhenGO/core/PhenGO.py`
- `src/PhenGO/core/go_handling.py`
- `src/PhenGO/core/phenotype_handling.py`

### GO Ancestor And Filtering Fixes

The core pipeline now keeps a clearer distinction between:

- direct GO terms from the source annotation file, stored as `go_list`
- expanded GO terms after ancestor traversal, stored as `expanded_go_list`
- binary ARFF feature vector, stored as `bin_vec`

`assign_go_to_vector()` now stores `expanded_go_list` for each gene. This fixes
downstream ambiguity where some outputs used only direct GO terms even though
the ARFF feature vector contained ancestor-expanded GO features.

`removed_unused_gos()` now derives used GO terms from `bin_vec` or
`expanded_go_list`, rather than only the direct `go_list`. This prevents valid
ancestor-derived features from being removed when `-no_filter_unused_gos` is not
used.

`get_o_rep_output()` and `write_func_output()` now use expanded GO terms when
available, keeping enrichment/FUNC-style outputs consistent with the ARFF
feature space.

Why this matters:

- The ARFF matrix, enrichment outputs, and GO feature filtering now describe the
  same expanded GO feature universe.
- Ancestor GO terms should no longer disappear from filtered outputs simply
  because they were not direct annotations.

### Safer ARFF Writing

`write_arff_output()` now writes data rows using Python's `csv.writer`.

Why this matters:

- Gene identifiers containing commas or quote-sensitive characters are emitted
  safely.
- Downstream CSV-aware parsers can round-trip those identifiers correctly.

### Safer Output Directory Handling

The core `PhenGO` command now has an explicit:

```bash
-overwrite
```

If an output directory already contains non-log files, PhenGO now refuses to
delete them unless `-overwrite` is supplied.

Why this matters:

- Previous behaviour could remove existing output contents unexpectedly.
- Runs are now safer and more reproducible.

### Mouse Parser And Join Fixes

Mouse phenotype parsing now supports both common MGI report layouts:

- `MGI_GenePheno`
- `MGI_PhenoGenoMP`

The mouse parser now prefers MGI marker accession IDs where available. Mouse GO
joining now tries both the MGI accession and gene symbol when matching GAF rows
to phenotype-derived gene records.

The transgene/multi-marker filter is now applied to the accession field used for
matching.

Why this matters:

- MGI phenotype files and MGI GAF files do not always line up cleanly by symbol.
- Using marker accessions improves the join between phenotype labels and GO
  annotations.
- Mouse ARFF files generated before this fix may have missed valid GO joins.

### Worm Parser Robustness

The worm phenotype parser now tolerates rows without the expected evidence-code
column.

Why this matters:

- Historical or partially normalised WormBase files are less likely to fail
  during parsing.

## PhenGO-Predict Changes

Files changed:

- `src/PhenGO/predict/PhenGO-Predict.py`
- `src/PhenGO/predict/cli.py`
- `src/PhenGO/predict/__init__.py`
- `src/PhenGO/predict/data.py`
- `src/PhenGO/predict/evaluate.py`
- `src/PhenGO/predict/sklearn_models.py`
- `src/PhenGO/predict/visualise.py`

### Importable CLI Entry Point

A new importable wrapper was added:

```text
src/PhenGO/predict/cli.py
```

It loads the legacy `PhenGO-Predict.py` script through `importlib` and exposes a
standard `main()` function.

Why this matters:

- Console-script entry points can call the predictor reliably even though the
  legacy file name contains a hyphen.

### Lazy Imports And Optional TensorFlow

`PhenGO-Predict.py`, `evaluate.py`, and `predict/__init__.py` were changed to
avoid importing TensorFlow, NumPy, sklearn, or other heavy ML dependencies until
they are actually needed.

Why this matters:

- sklearn-only workflows do not fail just because TensorFlow is missing.
- CLI help and lightweight imports work in more environments.
- The package API avoids eager loading of heavy dependencies.

### Explicit Lethal-Positive Label Encoder

`data.py` now defines:

```python
PhenotypeLabelEncoder
```

Supported lethal/positive aliases include:

- `lethal`
- `inviable`
- `essential`

Supported viable/negative aliases include:

- `viable`
- `non-essential`
- `nonessential`
- `non_essential`

The explicit class order is:

```text
["viable", "lethal"]
```

Why this matters:

- Metrics, thresholds, ROC/AUPRC calculations, predicted probabilities, and
  positive-class interpretation now align with the intended biological meaning:
  lethal/essential is class `1`.

### CSV-Aware ARFF Loading

`load_arff_data()` now uses `csv.reader` for data rows instead of manual comma
splitting.

Why this matters:

- ARFF rows with quoted gene identifiers containing commas are parsed correctly.

### Safer Prediction Output Directories

`PhenGO-Predict.py` now has:

```bash
-overwrite
```

The predictor refuses to delete a non-empty output directory unless this flag is
supplied.

### Threshold Boundary Fix

Neural-network binary predictions now use:

```python
y_pred_proba_test >= opt_threshold
```

instead of:

```python
y_pred_proba_test > opt_threshold
```

Why this matters:

- Predictions exactly equal to the selected threshold are handled consistently
  as positive.

### Model Reproducibility And Class Imbalance

`sklearn_models.py` now uses `options.random_state` or `options.seed` where
available instead of a hard-coded seed of `42`.

Gradient Boosting now receives balanced `sample_weight` during fitting. Before
this change, `dt`, `rf`, `lr`, and `svm` used class balancing, but `gb` did not.

Why this matters:

- Model comparisons are less confounded by class imbalance.
- Repeated analyses can vary seeds in a controlled way.

### Prediction Output Columns

`evaluate_and_analyse_predictions()` now adds:

```text
probability_lethal
```

as a clear alias for the positive-class prediction probability.

Why this matters:

- Downstream reports and version-sensitivity analyses can refer to the lethal
  probability unambiguously.

### Feature Importance Plot Fix

`visualise.py` now uses `nlargest()` rather than `nsmallest()` for top feature
importance plots.

Why this matters:

- Feature-importance plots now show the most important features, not the least
  important ones.

## New Version-Sensitivity Analysis Workflow

Files added:

- `src/PhenGO/predict/version_sensitivity.py`
- `docs/version_sensitivity_workflow.md`

The new command is:

```bash
phengo-version-sensitivity
```

Alias:

```bash
version-sensitivity
```

The command accepts either:

```bash
-arff_files file1.arff file2.arff ...
```

or:

```bash
-input_dir parent_directory
```

When `-input_dir` is used, the command discovers ARFF files from direct child
directories containing `*_PhenGO.arff`. If no child directories match, it falls
back to direct `*.arff` files in the input directory.

### Supported Models

Main models:

- `lr`
- `rf`
- `gb`
- `dt`

Optional model:

- `svm`

Convenience value:

- `all`

The default model set is:

```text
lr rf gb dt
```

### Analysis Panels

The command runs four panel types:

| Panel | Purpose |
| --- | --- |
| `full` | Uses the union of train/test GO features and each dataset's available genes. This captures the full practical effect of database-version changes. |
| `matched_features` | Restricts every dataset to GO features shared by all datasets. This controls for GO feature-space churn. |
| `matched_genes` | Restricts each train/test pair to shared genes. This controls for gene coverage churn. |
| `matched_both` | Restricts all datasets to genes and GO features shared by all datasets. This is the strongest control panel. |

### Output Files

The command writes:

```text
dataset_drift.csv
pairwise_drift.csv
within_year_cv.csv
within_year_cv_summary.csv
cross_year_transfer.csv
transfer_matrices/
previous_year_label_baseline.csv
prediction_instability.csv
prediction_instability_summary.csv
feature_rank_overlap.csv
VERSION_SENSITIVITY_REPORT.md
version_sensitivity.log
```

### Dataset Drift Outputs

`dataset_drift.csv` reports per-snapshot counts and structure:

- number of genes
- number of GO features
- number of active GO features
- lethal and viable counts
- lethal prevalence
- mean and median GO terms per gene
- sparsity

`pairwise_drift.csv` reports pairwise snapshot differences:

- shared genes
- gene-set Jaccard similarity
- shared GO features
- GO-feature Jaccard similarity
- label churn among common genes

### Within-Year CV Outputs

`within_year_cv.csv` contains fold-level repeated stratified CV results for each
dataset, model, and panel.

`within_year_cv_summary.csv` contains mean and standard deviation summaries for
the same metrics.

### Cross-Year Transfer Outputs

`cross_year_transfer.csv` contains train-year to test-year results for every
dataset pair, model, and panel.

`transfer_matrices/` contains matrix-form CSVs for key metrics, ready for
heatmaps.

Metrics include:

- ROC-AUC
- average precision / AUPRC
- balanced accuracy
- MCC
- F1 lethal
- F1 macro
- lethal precision
- lethal recall
- accuracy
- Brier score
- lethal prevalence
- sample counts

### Baselines

The analysis includes:

- `baseline_majority`
- `baseline_prior`
- `baseline_stratified_random`
- `previous_year_labels`

The first three are included in CV and transfer analyses. The previous-year
label baseline is written to `previous_year_label_baseline.csv`.

### Prediction Instability

`prediction_instability.csv` reports per-gene probability variation across
years in the `matched_both` panel.

Important fields:

- `sd_probability_lethal`
- `range_probability_lethal`
- `prediction_flip`
- `label_flip`
- `min_probability_lethal`
- `max_probability_lethal`

`prediction_instability_summary.csv` reports model-level instability:

- mean and median per-gene probability SD
- mean per-gene probability range
- prediction flip rate
- label flip rate

### Feature-Rank Overlap

`feature_rank_overlap.csv` reports top-k GO feature overlap across years for
models with native importance scores:

- `lr`
- `dt`
- `rf`
- `gb`

Low overlap suggests that the GO terms highlighted by the models are not stable
across database versions.

## Script And Report Parser Fixes

Files changed:

- `src/PhenGO/scripts/GO_Compare.py`
- `src/PhenGO/scripts/GO_temporal_analysis.py`
- `src/PhenGO/scripts/arff_validator.py`
- `src/PhenGO/scripts/compare_arff_genes.py`
- `src/PhenGO/scripts/report_generator.py`

Changes:

- ARFF data-row parsing now uses `csv.reader` in multiple scripts.
- Gene identifiers containing commas are now handled consistently.
- `arff_validator.py` now checks whether the input ARFF exists before running.
- `report_generator.py` now recognises both:
  - `report.txt`
  - `final_report.txt`
- Report ARFF stats parsing is now CSV-aware.

Why this matters:

- The ARFF ecosystem is now consistent across generation, validation,
  comparison, temporal analysis, prediction, and reporting.

## Tests Added

Files added:

- `tests/test_predict_labels.py`
- `tests/test_core_go_filtering.py`
- `tests/test_arff_io.py`
- `tests/test_mouse_parser.py`
- `tests/test_version_sensitivity.py`

Coverage added:

- lethal/viable label polarity:
  - confirms `viable -> 0`
  - confirms `lethal`, `inviable`, and `essential -> 1`
- GO filtering:
  - confirms ancestor-expanded GO terms are preserved when unused GO features
    are removed
- ARFF IO:
  - confirms a gene name containing commas can be written and parsed correctly
- mouse parsing:
  - confirms `MGI_GenePheno` marker accessions are parsed
  - confirms `MGI_PhenoGenoMP` marker accessions are parsed
  - confirms mouse GO joins prefer marker accession IDs
- version-sensitivity helpers:
  - confirms `all` model expansion and deduplication
  - confirms natural sorting of year-like dataset names
  - confirms ARFF auto-discovery from output subdirectories

## Recommended Use For The Research Article

### 1. Regenerate ARFF Files

For publication results, regenerate ARFF files from the raw model-organism
database snapshots using the corrected PhenGO core pipeline.

This is important because old ARFFs may have been affected by:

- lethal/viable positive-class interpretation during prediction
- GO ancestor filtering issues
- mouse MGI identifier matching issues
- comma-sensitive ARFF parsing/writing issues
- older output/report inconsistencies

Example:

```bash
PhenGO \
  -species worm \
  -phenotype_file data/worm/phenotype_data/2017/phenotype_association.WS257.wb.clean.col7.gz \
  -gene_association_file data/worm/gene_association/go_archive/2017/wb.gaf.gz \
  -go_obo_file data/go/2017/go_2017-01-01.obo.gz \
  -output_dir results/worm_2017 \
  -overwrite
```

Repeat this for each organism/year or organism/release.

### 2. Run One Version-Sensitivity Analysis Per Organism

Example:

```bash
phengo-version-sensitivity \
  -arff_files \
    results/worm_2017/worm_PhenGO.arff \
    results/worm_2021/worm_PhenGO.arff \
    results/worm_2025/worm_PhenGO.arff \
  -dataset_names 2017 2021 2025 \
  -models lr rf gb dt \
  -cv_folds 5 \
  -cv_repeats 5 \
  -output_dir paper_outputs/worm_version_sensitivity \
  -overwrite
```

Recommended main-manuscript models:

```text
lr rf gb dt
```

Recommended supplementary model:

```text
svm
```

### 3. Structure Results Around The Version Question

Suggested order:

1. Use `dataset_drift.csv` and `pairwise_drift.csv` to show that the source
   resources changed between years.
2. Use `within_year_cv_summary.csv` to show how ordinary model performance
   changes inside each snapshot.
3. Use `cross_year_transfer.csv` and `transfer_matrices/` as the central
   evidence for year/version sensitivity.
4. Compare `full`, `matched_features`, `matched_genes`, and `matched_both`
   panels to identify whether instability is driven by gene churn, GO-feature
   churn, label churn, annotation churn, or all of these.
5. Use `prediction_instability.csv` to show gene-level changes in model output.
6. Use `feature_rank_overlap.csv` to show whether the apparent GO signal learned
   by the model is stable across years.

### 4. Preferred Metrics

Because lethal genes are often the minority class, prioritise:

- AUPRC / average precision
- MCC
- balanced accuracy
- F1 lethal
- lethal recall
- lethal precision

Report class prevalence alongside performance metrics so the reader can see
whether apparent performance changes are partly caused by label distribution
changes.

## Verification Performed

The following checks were run with:

```text
/Users/nicholas/anaconda3/envs/PhenGO/bin/python
```

Unit tests:

```bash
PYTHONPATH=src /Users/nicholas/anaconda3/envs/PhenGO/bin/python -m unittest discover -s tests
```

Result:

```text
Ran 9 tests
OK
```

Version-sensitivity CLI help:

```bash
PYTHONPATH=src /Users/nicholas/anaconda3/envs/PhenGO/bin/python -m PhenGO.predict.version_sensitivity -h
```

Result:

```text
usage: phengo-version-sensitivity ...
```

Whitespace check:

```bash
git diff --check
```

Result:

```text
no whitespace errors reported
```

## Known Environment Issue

The local conda environment currently has a NumPy binary architecture mismatch:

```text
NumPy extension architecture: x86_64
Python architecture needed: arm64 / arm64e
```

This means full ML runs that import NumPy, pandas, or sklearn will fail in that
environment until NumPy is reinstalled for the correct architecture. The
lightweight unit tests and CLI help checks pass, but end-to-end sklearn analyses
need the environment repaired first.

## Generated Files Present In The Working Tree

Running tests created Python bytecode cache files under directories such as:

```text
tests/__pycache__/
src/PhenGO/predict/__pycache__/
```

These are generated artifacts, not source changes. They should normally be
ignored or removed before committing.

## Full Source File Change List

Modified tracked files:

- `README.md`
- `setup.cfg`
- `src/PhenGO/constants.py`
- `src/PhenGO/core/PhenGO.py`
- `src/PhenGO/core/go_handling.py`
- `src/PhenGO/core/phenotype_handling.py`
- `src/PhenGO/predict/PhenGO-Predict.py`
- `src/PhenGO/predict/__init__.py`
- `src/PhenGO/predict/data.py`
- `src/PhenGO/predict/evaluate.py`
- `src/PhenGO/predict/sklearn_models.py`
- `src/PhenGO/predict/visualise.py`
- `src/PhenGO/scripts/GO_Compare.py`
- `src/PhenGO/scripts/GO_temporal_analysis.py`
- `src/PhenGO/scripts/arff_validator.py`
- `src/PhenGO/scripts/compare_arff_genes.py`
- `src/PhenGO/scripts/report_generator.py`

Added source/documentation/test files:

- `CHANGELOG_SINCE_LAST_COMMIT.md`
- `docs/version_sensitivity_workflow.md`
- `src/PhenGO/predict/cli.py`
- `src/PhenGO/predict/version_sensitivity.py`
- `tests/test_arff_io.py`
- `tests/test_core_go_filtering.py`
- `tests/test_mouse_parser.py`
- `tests/test_predict_labels.py`
- `tests/test_version_sensitivity.py`
