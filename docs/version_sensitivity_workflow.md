# PhenGO Version-Sensitivity ML Workflow

## Why This Was Added

The original PhenGO-Predict workflow asked a conventional machine-learning
question: given one ARFF file, how well can a model predict lethal versus viable
genes inside that snapshot?

The research article needs a sharper question:

> Do model-organism database versions materially change machine-learning
> conclusions, potentially making the output reflect resource-version history
> rather than biology alone?

The new `phengo-version-sensitivity` command is designed around that question.
It treats each yearly ARFF as one snapshot of a changing database ecosystem and
compares model behavior across those snapshots.

## What Changed

### 1. Version-Sensitivity Analysis Command

New entry points:

```bash
phengo-version-sensitivity
version-sensitivity
```

Implementation:

```text
src/PhenGO/predict/version_sensitivity.py
```

The command accepts multiple ARFF files or auto-discovers ARFF files from a
parent directory. It produces CSV tables and a generated Markdown report in the
output directory.

### 2. Within-Year Repeated Cross-Validation

For each dataset/year, model, and sensitivity panel, the command runs repeated
stratified cross-validation. This reports the mean and standard deviation of
model performance within each database snapshot.

This is useful for showing that ordinary predictive performance changes across
years.

### 3. Cross-Year Transfer Matrices

For every pair of datasets, models are trained on one year and tested on another
year. This directly measures whether a model trained on one database snapshot
transfers to another snapshot.

This is the central analysis for the article. A strong version effect appears
when off-diagonal train-year/test-year performance drops or behaves
asymmetrically.

### 4. Matched Panels

The command runs four panels:

| Panel | Purpose |
| --- | --- |
| `full` | Uses all GO features in each train/test pair, filling missing terms with 0. Captures the full practical effect of changing resources. |
| `matched_features` | Restricts all datasets to GO features shared by every dataset. Isolates effects beyond GO feature churn. |
| `matched_genes` | Restricts train/test pairs to shared genes. Isolates effects beyond gene coverage changes. |
| `matched_both` | Restricts all datasets to shared genes and shared GO features. Strongest controlled panel for label/annotation instability. |

These panels let the manuscript separate multiple sources of instability:

- genes entering or leaving the resource
- GO terms entering or leaving the feature space
- phenotype label changes
- GO annotation changes for the same genes

### 5. Dataset Drift Tables

The command writes:

```text
dataset_drift.csv
pairwise_drift.csv
```

These quantify the source-resource changes themselves:

- number of genes
- number of GO features
- active GO features
- lethal prevalence
- sparsity
- GO terms per gene
- gene-set Jaccard similarity
- GO-feature Jaccard similarity
- label churn among shared genes

These tables should appear early in the results section, before model outputs.

### 6. Prediction Instability Tables

The command writes:

```text
prediction_instability.csv
prediction_instability_summary.csv
```

These use the `matched_both` panel to ask: for genes shared across all years and
features shared across all years, do model probabilities still move around?

Useful fields include:

- `sd_probability_lethal`
- `range_probability_lethal`
- `prediction_flip`
- `label_flip`
- `prediction_flip_rate`
- `label_flip_rate`

These are strong evidence when the same gene receives materially different
model outputs depending on the database year.

### 7. Feature-Rank Overlap

The command writes:

```text
feature_rank_overlap.csv
```

For models with native feature importance (`lr`, `dt`, `rf`, `gb`), this reports
top-k GO term overlap between years. Low overlap suggests that the apparent
biological signal learned by the model is not stable across database versions.

### 8. Baselines

The command includes:

- `baseline_majority`
- `baseline_prior`
- `baseline_stratified_random`
- `previous_year_labels`

The first three appear in CV and transfer outputs. The previous-year label
baseline appears in:

```text
previous_year_label_baseline.csv
```

Baselines make it easier to argue that performance differences are meaningful
and not simply class imbalance artifacts.

### 9. Expanded Metrics

The new analysis reports:

- ROC-AUC
- average precision / AUPRC
- balanced accuracy
- MCC
- F1 lethal
- F1 macro
- lethal precision
- lethal recall
- Brier score
- lethal prevalence
- sample counts

AUPRC, F1-lethal, balanced accuracy, and MCC are especially important because
lethal genes are usually the minority class.

### 10. Model Improvements

Gradient boosting now uses balanced sample weights. Previously, `dt`, `rf`,
`lr`, and `svm` used class balancing, while gradient boosting did not. That made
gradient boosting more sensitive to class imbalance for reasons unrelated to
database version.

The sklearn model builder now accepts a configurable random seed through the
options namespace, which allows repeated analyses to vary seeds in a controlled
way.

## Recommended Article Workflow

### Step 1. Generate ARFFs For Every Organism-Year

Use the corrected PhenGO core pipeline to generate one ARFF per organism/year or
organism/release. Keep output directories separate and preserve the
`PhenGO_params.txt` files.

Example:

```bash
PhenGO \
  -species worm \
  -phenotype_file data/worm/phenotype_data/2017/phenotype_association.WS257.wb.clean.col7.gz \
  -gene_association_file data/worm/gene_association/go_archive/2017/wb.gaf.gz \
  -go_obo_file data/go/2017/go_2017-01-01.obo.gz \
  -output_dir results/worm_2017
```

Repeat for all years/releases.

Important: because earlier ARFF generation had issues around GO ancestor
filtering and mouse identifier matching, regenerate ARFF files before running
the article analysis.

### Step 2. Run Version-Sensitivity Analysis Per Organism

Example for one organism:

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

Recommended models for the main manuscript:

```text
lr rf gb dt
```

Use `svm` as a supplementary analysis if runtime is acceptable. The RBF SVM can
be slow on large, sparse GO matrices.

### Step 3. Use Dataset Drift Tables First

Start results with:

```text
dataset_drift.csv
pairwise_drift.csv
```

Suggested claims to support:

- database snapshots differ in gene coverage
- database snapshots differ in GO feature coverage
- lethal prevalence changes over time
- labels for shared genes can change
- annotation sparsity changes over time

These establish the resource-version problem before ML is discussed.

### Step 4. Use Within-Year CV To Show Snapshot-Specific Performance

Use:

```text
within_year_cv_summary.csv
```

Report mean ± standard deviation for AUPRC, ROC-AUC, F1-lethal, balanced
accuracy, and MCC.

Suggested figure:

- line plot of AUPRC by year for each model
- line plot of F1-lethal by year for each model

Interpretation:

If within-year performance varies strongly, a model trained on one database
snapshot may not represent a stable biological signal.

### Step 5. Use Cross-Year Transfer As The Central Evidence

Use:

```text
cross_year_transfer.csv
transfer_matrices/
```

Suggested figures:

- heatmap of train-year/test-year AUPRC for each organism and model
- heatmap of train-year/test-year F1-lethal
- heatmap of train-year/test-year MCC

Interpretation:

- strong diagonal performance with weak off-diagonal performance suggests
  snapshot-specific learning
- asymmetric transfer suggests annotation/label evolution rather than generic
  model weakness
- poor transfer even in `matched_both` suggests that labels or annotations for
  shared genes are changing enough to alter ML conclusions

### Step 6. Compare Panels To Separate Causes

Use panel-specific rows in:

```text
cross_year_transfer.csv
within_year_cv_summary.csv
```

Interpretation guide:

| Pattern | Interpretation |
| --- | --- |
| `full` unstable, `matched_features` stable | GO feature churn is a major driver. |
| `full` unstable, `matched_genes` stable | gene coverage churn is a major driver. |
| `matched_both` still unstable | labels or annotations for the same genes/features changed. |
| all panels unstable | version effects are broad and not reducible to one source. |

### Step 7. Use Prediction Instability For Gene-Level Evidence

Use:

```text
prediction_instability.csv
prediction_instability_summary.csv
```

Suggested figure:

- histogram of per-gene `range_probability_lethal`
- bar plot of `prediction_flip_rate` by model
- table of the top genes by probability range

Interpretation:

This shows that not only metrics but also gene-level conclusions can move across
database versions.

### Step 8. Use Feature-Rank Overlap For Biological Signal Stability

Use:

```text
feature_rank_overlap.csv
```

Suggested figure:

- heatmap of top-20 GO term Jaccard overlap between years

Interpretation:

Low overlap suggests that model explanations and highlighted GO processes are
not stable across database versions.

## Recommended Manuscript Framing

Avoid presenting these models as final predictors of gene essentiality. Present
them as standardized probes of database-version sensitivity.

Suggested wording:

> We used multiple common classifier families as probes of annotation-version
> sensitivity. The aim was not to optimize a production predictor, but to test
> whether ML conclusions remain stable when the underlying model-organism
> database snapshot changes.

The strongest claim is supported when:

1. dataset drift is measurable,
2. within-year performance varies by snapshot,
3. cross-year transfer matrices show poor or asymmetric transfer,
4. matched panels still show instability,
5. gene-level predictions flip for shared genes,
6. top GO features learned by models are not stable.

## Output-To-Figure Mapping

| Paper figure/table | Source file |
| --- | --- |
| Dataset drift table | `dataset_drift.csv` |
| Pairwise gene/GO/label drift heatmap | `pairwise_drift.csv` |
| Within-year model performance line plot | `within_year_cv_summary.csv` |
| Cross-year transfer heatmaps | `transfer_matrices/*.csv` |
| Controlled-panel comparison | `cross_year_transfer.csv` |
| Prediction instability histogram | `prediction_instability.csv` |
| Model-level flip-rate bar chart | `prediction_instability_summary.csv` |
| GO importance overlap heatmap | `feature_rank_overlap.csv` |
| Previous-year label baseline table | `previous_year_label_baseline.csv` |

## Practical Notes

- Use `lr` and `rf` as main models; they are easier to defend and interpret.
- Use `gb` and `dt` to show the effect is not limited to one model family.
- Put `svm` in supplementary material unless datasets are small enough for
  comfortable runtime.
- Keep the same ARFF-generation policy for every year in an organism.
- Always report class prevalence alongside performance metrics.
- Prefer AUPRC, MCC, balanced accuracy, and F1-lethal over raw accuracy.
- Report results from the `matched_both` panel prominently; it is the strongest
  control against “you just changed the genes/features” critiques.
