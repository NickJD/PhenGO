# PhenGO Version-Sensitivity Research Workflow

## Research Question

PhenGO is intended here as a controlled instrument for asking:

> How much do machine-learning conclusions change when the model-organism
> database snapshot changes, even though the organism and nominal prediction
> task remain the same?

The objective is not to build the highest-performing essentiality predictor.
The models are standardized probes of version sensitivity. A defensible paper
must therefore preserve release provenance, apply one pre-specified feature and
label policy to every snapshot, use out-of-fold evaluation, and report resource
drift before interpreting model drift.

## Important Compatibility Break

Publication datasets must be regenerated with the current PhenGO core. Existing
ARFFs do not contain the strict snapshot manifest and may differ because the new
pipeline:

- fixes the positive-class mapping to `viable = 0`, `lethal = 1`;
- writes canonical database accessions rather than mutable symbols;
- excludes mixed labels and requires explicit viability evidence by default;
- excludes partial/semi-lethal observations by default;
- filters fly multi-gene phenotypes by default;
- excludes `NOT` and selected GAF evidence annotations exactly;
- canonicalizes alternate GO IDs and one-to-one obsolete replacements;
- supports relation-aware GO propagation, with `is_a + part_of` as the default;
- separates release-derived, paired gene-set, and agreement label sources;
- removes ontology roots and very rare or near-universal features by default;
- records input SHA-256 hashes, releases, policies, software, and outputs.

Do not combine old and regenerated ARFFs in one version analysis.

## 1. Freeze The Protocol

Before generating results, define the following in a protocol or analysis plan:

- organisms and releases to include;
- the exact phenotype, GAF, GO OBO, phenotype ontology, and helper file for each release;
- the label policies;
- GO relation and namespace policies;
- GO evidence-code inclusion/exclusion policy;
- GO prevalence filters;
- primary models and metrics;
- primary panel and secondary sensitivity panels;
- minimum class size and criteria for excluding a snapshot;
- planned figures and statistical comparisons.

Recommended primary settings:

```text
mixed_label_policy       exclude
nonlethal_policy         explicit_only
ambiguous_label_policy   exclude
label_source              release_records
go_relations             is_a part_of
go_ontology_file          release-matched go-basic.obo
go_propagation           ancestors
exclude_go_roots         true
min_go_gene_count        2
max_go_gene_fraction     0.99
models                    dt rf gb lr svm nn
primary metrics           average precision, MCC, balanced accuracy, F1-lethal
primary panel             matched_both
practical-effect panel    full
```

`is_a` and `part_of` are the primary relations because positive GO annotations
can be grouped safely through both. Regulation, occurrence, and capability
relations change the meaning of the propagated gene-term association and are
reserved for separately declared sensitivity analyses.

```bash
-go_relations is_a part_of
```

The toolkit also accepts `regulates`, `positively_regulates`,
`negatively_regulates`, `occurs_in`, `capable_of`, and
`capable_of_part_of`. These relations are not interchangeable with `is_a`.
Keep one relation set fixed across all years and describe its interpretation.
PhenGO rejects cyclic relation graphs and cross-namespace propagation by
default. Use a release-matched `go-basic.obo`; do not infer suitability merely
from a filename.

### Pre-specify the label design

Use `release_records` for the primary as-released analysis. Use paired gene
sets as a second, fixed-label experiment that isolates GO-resource change:

```text
release_records   yearly phenotype labels + yearly GAF + yearly go-basic
gene_sets         fixed paired labels + yearly GAF + yearly go-basic
agreement         yearly labels retained only when paired sets agree
```

Do not merge a lethal-only collection into release labels for manuscript runs.
Both lethal and viable sets are required, and unlisted genes remain unknown.
The same fixed collections should be used across every snapshot in the
fixed-label experiment; applying a later collection while calling it a
historical label source would leak future curation into earlier years.

## 2. Assemble Release-Matched Inputs

Create an input ledger before running PhenGO. One row should represent one
organism snapshot and include:

```text
organism
snapshot_id
phenotype_file
phenotype_release
go_annotation_file
go_annotation_release
go_basic_obo_file
go_ontology_release
phenotype_ontology_or_helper_file
phenotype_ontology_release
explicit_viable_terms_or_genes
label_source
lethal_gene_set_if_used
viable_gene_set_if_used
retrieval_date
archive_url_or_accession
notes
```

Do not use bundled current helper files for historical snapshots. The normal
CLI permits them for exploratory work but emits a warning. `-strict_snapshot`
rejects a historical run unless the required auxiliary files were supplied
explicitly.

## 3. Generate One Strict ARFF Per Snapshot

Example for WormBase phenotype-term mode:

```bash
PhenGO \
  -species worm \
  -phenotype_file inputs/worm/WS257/phenotype_association.WS257.tsv.gz \
  -gene_association_file inputs/worm/WS257/wb.gaf.gz \
  -go_obo_file inputs/worm/WS257/go-basic.WS257.obo.gz \
  -worm_phenotypes inputs/worm/WS257/lethal_terms.WS257.tsv.gz \
  -worm_viable_phenotypes inputs/worm/WS257/viable_terms.WS257.tsv.gz \
  -output_dir results/worm_WS257 \
  -snapshot_id worm_WS257 \
  -phenotype_release WS257 \
  -go_annotation_release WS257 \
  -go_ontology_release WS257 \
  -phenotype_ontology_release WS257 \
  -retrieval_date 2026-07-17 \
  -strict_snapshot \
  -mixed_label_policy exclude \
  -nonlethal_policy explicit_only \
  -ambiguous_label_policy exclude \
  -go_relations is_a part_of \
  -go_propagation ancestors \
  -exclude_go_evidence_codes ND IEA \
  -min_go_gene_count 2 \
  -max_go_gene_fraction 0.99
```

For yeast and fish, `-phenotype_ontology_release` is not required when no
separate phenotype ontology/helper resource is used. Fly, worm, and mouse
strict runs require explicit release-matched auxiliary resources.
With `explicit_only`, worm phenotype-term mode requires
`-worm_viable_phenotypes`, worm lethal-gene mode requires
`-worm_viable_genes`, and mouse requires `-mouse_viable_phenotypes`. These
files must be constructed under a documented, release-matched definition of
explicit viability.

Every output directory should contain at least:

```text
{species}_PhenGO.arff
gene_identifiers.tsv
identifier_join_audit.tsv
label_source_audit.tsv
PhenGO_manifest.json
PhenGO_params.txt
GO_term_validation_report.txt
GO_Children&Parents.json
go_enrichment/
```

Treat `PhenGO_manifest.json` as part of the dataset. It contains input and
output hashes, release identifiers, label and GO policies, the command,
dependency versions, machine architecture, Python executable, git commit, and
row/feature counts.

For a fixed-label run, replace the phenotype-specific inputs with:

```bash
-label_source gene_sets \
-lethal_gene_set inputs/fixed/lethal_genes.tsv \
-viable_gene_set inputs/fixed/viable_genes.tsv
```

For the conservative agreement panel, retain the yearly phenotype inputs and
use `-label_source agreement` with the same paired files. Strict mode rejects
the bundled current fly collection and the deprecated lethal-only override.

## 4. Quality-Control Every Snapshot

Before ML analysis, verify:

1. The manifest has `strict_snapshot: true` and the intended `snapshot_id`.
2. Input hashes and release labels match the input ledger.
3. The ARFF contains stable accessions and `gene_identifiers.tsv` maps them to symbols.
4. Both classes are present with sufficient counts for the planned folds.
5. No duplicate genes or duplicate GO feature names exist.
6. GO values are binary and contain no missing values.
7. The validation report does not show an unexpected spike in missing or unresolved obsolete terms.
8. Label prevalence and GO feature count are biologically plausible for that release.
9. Every year used the identical policy settings.
10. `label_source_audit.tsv` has the expected retained/disagreement/source-only counts.
11. `identifier_join_audit.tsv` contains no unexplained ambiguous or contradictory mappings; unmatched accessions are retained as explicit exclusions.
12. The manifest reports `go_relations: [is_a, part_of]`, no cross-namespace edges, and the expected OBO `data-version` when present.

The version-sensitivity command verifies manifests and ARFF hashes by default.
`-allow_missing_manifests` is only for exploratory legacy analysis and should
not be used for manuscript results. It also rejects differences in source-tree
hash, tool version, dependency environment, label policy, GO policy, prevalence
filtering, or parser settings across snapshots.

## 5. Run One Analysis Per Organism

```bash
phengo-version-sensitivity \
  -arff_files \
    results/worm_WS247/worm_PhenGO.arff \
    results/worm_WS257/worm_PhenGO.arff \
    results/worm_WS270/worm_PhenGO.arff \
    results/worm_WS297/worm_PhenGO.arff \
  -dataset_names WS247 WS257 WS270 WS297 \
  -models all \
  -panels full matched_features matched_genes matched_both \
  -cv_folds 5 \
  -cv_repeats 5 \
  -calibration sigmoid \
  -calibration_cv 3 \
  -importance_repeats 20 \
  -seed 42 \
  -n_jobs 1 \
  -output_dir paper_outputs/worm_version_sensitivity
```

Before model fitting, the command writes `evaluation_preflight.csv` and requires
the requested fold and calibration design to be feasible for every selected
panel and transfer direction. It does not silently lower fold counts or omit
class-deficient cells in publication mode. Reduce the planned fold count for all
snapshots before analysis if necessary. `-allow_incomplete_evaluation` is an
exploratory override and must not be used for primary manuscript results.

Use a new output directory, or pass `-overwrite` deliberately. The command
will not delete an input-containing or otherwise protected directory.

## 6. Understand The Evaluation Design

### Within-snapshot performance

Each model is evaluated by repeated out-of-fold prediction. Metrics are first
computed from the complete out-of-fold predictions for each repeat; means,
standard deviations, and Student-t 95% interval half-widths are then calculated
across repeats. Individual folds are not treated as independent replicates.

### Cross-snapshot transfer

For each repeat, one fixed partition of stable gene identifiers is reused for
every train-release/test-release direction. Every scored test gene is removed
from that fold's training data if it is present in the train release. This makes
diagonal and off-diagonal results gene-disjoint and out-of-fold while allowing
one fitted train-year/fold model to score every test year. Test-year-only GO
terms are ignored because a deployed model cannot learn weights or splits for a
feature absent from its training release.

### Probability calibration

Sigmoid calibration is fitted only within the outer training fold. `isotonic`
is available but generally needs more data. Threshold-free metrics remain
important because any 0.5 decision threshold can obscure ranking changes.

### Comparable neural network

In this command, `nn` is scikit-learn's `MLPClassifier`, not the optional
TensorFlow single-snapshot implementation. It is fitted within exactly the same
outer folds and transfer directions as the other classifiers, receives the
same class-balanced training weights, uses the same nested calibration policy,
and is scored with the same fixed lethal-positive class mapping, threshold, and
metrics. Its default architecture is `128,64` ReLU units with Adam optimization,
L2 regularization, and training-only early stopping. The architecture,
optimizer limits, seed, and validation fraction are fixed across releases and
recorded in `version_sensitivity_manifest.json`.

The fold and transfer tables include `fit_convergence_warnings` and
`fit_n_iter_max`. Neural-network results with convergence warnings must be
resolved by a pre-specified global increase to `-nn_max_iter` and rerunning all
years; do not increase the limit only for an inconvenient release. The optional
TensorFlow `nn` from `phengo-predict` uses a different training and validation
design and is therefore supplementary, not a row in the primary cross-model
comparison.

### Feature importance

Native importance is estimated over stratified bootstrap samples for decision
trees, random forests, gradient boosting, and logistic regression. SVM and MLP
are excluded from this native-importance comparison because their coefficients
are not directly comparable to those estimators. The output
contains mean importance, standard deviation, top-k selection frequency, and a
consensus rank. This exposes unstable GO rankings rather than presenting one
fit as definitive.

## 7. Use The Four Panels Correctly

| Panel | Cohort | Feature space | Interpretation |
| --- | --- | --- | --- |
| `full` | Each snapshot's genes | Training release | Practical effect of deploying each complete training release. |
| `matched_features` | Each snapshot's genes | Features shared by all snapshots | Controls global GO feature-set churn. |
| `matched_genes` | Genes shared by all snapshots | Training release | Controls global gene-set churn. |
| `matched_both` | Genes shared by all snapshots | Features shared by all snapshots | Strongest controlled panel for annotation and label instability. |

Interpret panel differences rather than selecting whichever panel gives the
largest effect. `full` and `matched_both` answer complementary questions.

## 8. Output Files

### Resource drift

```text
dataset_drift.csv
pairwise_drift.csv
```

These contain gene counts, feature counts, class prevalence, sparsity, GO terms
per gene, gene/feature Jaccard similarity, and label churn among shared genes.

### Within-snapshot models

```text
within_year_cv.csv
within_year_cv_summary.csv
```

The first is fold-level diagnostic output. Use the summary for manuscript
tables and plots.

### Temporal transfer

```text
cross_year_transfer.csv
cross_year_transfer_summary.csv
transfer_matrices/
```

Use the summary and matrices for primary results. The raw file contains one
complete out-of-fold estimate per repeat.

### Gene-level stability

```text
prediction_instability.csv
prediction_instability_summary.csv
```

All probabilities are fixed-fold, gene-disjoint out-of-fold predictions in the
`matched_both` panel.

### GO importance stability

```text
feature_rank_overlap.csv
feature_importance_stability.csv
```

Report bootstrap selection frequencies alongside pairwise top-k overlap.

### Baselines and provenance

```text
previous_available_snapshot_label_baseline.csv
version_sensitivity_manifest.json
VERSION_SENSITIVITY_REPORT.md
version_sensitivity.log
```

Majority, class-prior, and stratified-random baselines also appear in within-
and cross-snapshot outputs.

## 9. Primary Results Sequence

1. **Describe source drift.** Report genes, features, prevalence, sparsity, and label churn.
2. **Report within-snapshot performance.** Prioritize average precision, MCC, balanced accuracy, and F1-lethal.
3. **Present transfer matrices.** Show directionality and compare diagonals with off-diagonals.
4. **Compare panels.** Attribute portions of instability to gene churn, feature churn, or within-cohort changes.
5. **Report gene-level movement.** Show probability ranges and prediction-flip rates.
6. **Report explanation stability.** Show bootstrap GO selection frequencies and top-k overlap.
7. **Report baselines and all planned models.** Do not omit models or years after inspecting results.

Always report lethal prevalence and evaluated gene counts with performance.
Accuracy alone is not suitable for an imbalanced target.

## 10. Suggested Figures And Tables

| Item | Source |
| --- | --- |
| Snapshot composition table | `dataset_drift.csv` |
| Gene, GO, and label drift heatmaps | `pairwise_drift.csv` |
| Within-release average precision plot | `within_year_cv_summary.csv` |
| Train-release/test-release heatmaps | `transfer_matrices/*.csv` |
| Full versus matched panel figure | `cross_year_transfer_summary.csv` |
| Per-gene probability-range distribution | `prediction_instability.csv` |
| Prediction-flip-rate table | `prediction_instability_summary.csv` |
| GO top-k overlap heatmap | `feature_rank_overlap.csv` |
| Bootstrap GO stability table | `feature_importance_stability.csv` |

## 11. Temporal GO Enrichment

`temporal-analysis` now computes exact Fisher tests, Haldane-Anscombe corrected
finite effect ratios, and Benjamini-Hochberg FDR values. Ranked timelines only
include terms meeting `--max-fdr` (default `0.05`), while
`*_enrichment_statistics.csv` preserves every tested term.

```bash
temporal-analysis \
  -input_dir results/worm_snapshots \
  -output_dir paper_outputs/worm_temporal_go \
  -o worm \
  --top-n 20 \
  --max-fdr 0.05
```

GO terms remain dependent through the ontology hierarchy, so enrichment is
descriptive even after FDR correction.
The command requires strict source manifests by default and writes
`temporal_analysis_manifest.json`. The legacy override
`-allow_missing_manifests` is not recommended for manuscript outputs.

## 12. Manuscript Language

Suggested methods framing:

> We used several common classifier families as standardized probes of
> database-version sensitivity rather than as optimized production predictors.
> These comprised decision tree, random forest, gradient boosting, logistic
> regression, support vector machine, and a feed-forward multilayer perceptron.
> All scored genes were excluded from their corresponding training folds,
> probability calibration was fitted within training data, and matched panels
> were used to distinguish gene-set and GO-feature-set churn from changes among
> shared genes and features.

Suggested interpretation:

> Differences between releases represent the combined effect of changing
> database content, ontology structure, evidence availability, identifier
> mappings, and curation practice. They do not establish that any individual
> release is biologically incorrect.

A complete manuscript-ready limitations section is provided in
[`manuscript_limitations.md`](manuscript_limitations.md).
