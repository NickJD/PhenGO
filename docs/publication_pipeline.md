# PhenGO One-Shot Publication Pipeline

## Purpose

[`scripts/run_publication_pipeline.sh`](../scripts/run_publication_pipeline.sh)
turns the archived MOD, phenotype, GAF, GO, and phenotype-ontology resources in
`data/` into a provenance-locked collection of ARFF datasets and then runs the
complete version-sensitivity analysis. The pipeline is designed to answer one
question: how much can machine-learning output change when the nominal organism
and prediction task are held constant but the database snapshot changes?

The workflow does not treat the best-performing model as the main result. Its
models are standardized measurement instruments. The principal evidence is the
change in within-year performance, cross-year transfer, per-gene predictions,
and GO associations across archived resource versions.

## Authoritative Inputs

[`data/publication_snapshots.tsv`](../data/publication_snapshots.tsv) is the
machine-readable input authority. It contains 80 rows:

- 49 required release-matched phenotype and GO snapshots;
- 31 annotation anchors used only in fixed-2025-label experiments.

The required series are fish, fly, mouse, and worm for 2015 through 2025, plus
yeast for 2015 through 2017 and 2024 through 2025. The annotation anchors add
2010 through 2014 for all organisms and 2018 through 2023 for yeast where the
required GO resources are present.

`data/File_Overview.ods` is retained as the historical collection notebook. It
is not read by the pipeline because it contains copied paths and does not encode
the label and analysis policies required for reproducibility. The cleaned human
overview is `data/File_Overview_publication.xlsx`.

## Identifier Safety Policy

PhenGO treats the stable database accession as the biological join key. ZFIN,
WormBase, and SGD accessions are retained directly from phenotype records;
release FlyBase symbols are translated through taxon-filtered, release-matched
FBgn assignments; and IMPC MGI accessions are retained when supplied. The GAF
is indexed globally before labels are assigned. A canonical symbol is accepted
only when it maps to one stable accession, and strict snapshots never use the
GAF synonym field. One-to-many symbols and contradictory strongest labels are
excluded rather than resolved heuristically.

The 2015 IMPC export is the one primary-source exception: it contains symbols
but no MGI accessions. Those symbols are matched only to unique exact canonical
symbols in the 2015 MGI GAF. The current archive resolves 1,208 of 1,212 labels
uniquely and leaves four unmatched; no synonym match is attempted. Later IMPC
exports contain MGI accessions. Fly assignment symbols that map to multiple
FBgn IDs are written to `01_derived_labels/fly/{year}/excluded.tsv` and omitted.
Every ARFF directory contains `identifier_join_audit.tsv`, and its exclusions
and counts are part of the manifest review gate.

## Historical GO Archive Policy

The complete planned data collection is present. It includes 2010-2014 GAFs for
all five organisms, full archived GO ontologies for 2010-2012, the first available
`go-simple` ontology for 2013, and `go-basic` from 2014 onward.

PhenGO validated the early ontologies with its production graph parser. Their
`is_a` graphs are acyclic and namespace-safe. Their `part_of` relations contain
cross-namespace links, however, including in the 2013 `go-simple` file. The
pipeline therefore does not pretend that `go-simple` and `go-basic` have the same
propagation guarantees. It runs a continuous 2010-2025 `fixed_2025_is_a` series
and a directly comparable 2014-2025 `fixed_2025` series using `is_a + part_of`.

All anchors are required by default. `--allow-missing-anchors` is available for
exploratory partial runs; skipped rows are recorded in
`00_run_metadata/optional_anchor_status.tsv`, and current files are never
substituted for historical inputs.

## Run Command

The publication run uses the ARM64 environment requested for this project and
writes to the requested Nextcloud directory by default:

```bash
cd /Users/nicholas/Git/PhenGO
./scripts/run_publication_pipeline.sh
```

Equivalent explicit invocation:

```bash
PHENGO_PYTHON=/Users/nicholas/miniconda3/envs/PhenGO/bin/python \
PHENGO_OUTPUT_DIR=/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs \
./scripts/run_publication_pipeline.sh
```

Useful controls:

```bash
# Validate paths and print every planned command.
./scripts/run_publication_pipeline.sh --dry-run

# Small end-to-end smoke run.
./scripts/run_publication_pipeline.sh --quick --output /tmp/phengo-smoke

# Continue an interrupted run without replacing completed units.
./scripts/run_publication_pipeline.sh --resume

# Run the legacy held-out suite for every year instead of 2025 only.
PHENGO_SINGLE_SNAPSHOT_YEARS=all ./scripts/run_publication_pipeline.sh

# Explicitly require every annotation anchor (this is already the default).
./scripts/run_publication_pipeline.sh --require-anchors

# Exploratory only: record and skip missing annotation anchors.
./scripts/run_publication_pipeline.sh --allow-missing-anchors
```

`--resume` compares a fingerprint of source code, the pipeline and ledger,
dependency versions, all authoritative input hashes, and analysis parameters.
If any component changed, it automatically invalidates resume mode and
regenerates every derived label, ARFF, and analysis. A new output root remains
preferable for a distinct manuscript experiment.

## Experimental Tracks

The default run generates five tracks for every selected organism:

| Track | Labels | GO settings | Scientific role |
| --- | --- | --- | --- |
| `primary` | Release-matched | `is_a + part_of` | Combined phenotype, annotation, and ontology version effect. |
| `is_a_only` | Release-matched | `is_a` | GO relation-policy sensitivity. |
| `no_iea` | Release-matched | `is_a + part_of`, IEA excluded | Electronic-annotation sensitivity. |
| `fixed_2025` | 2025 labels held fixed | `is_a + part_of`; `go-basic`; 2014-2025 | Direct GO-resource effect with targets and relation policy controlled. |
| `fixed_2025_is_a` | 2025 labels held fixed | `is_a`; full/simple/basic archive; 2010-2025 | Long historical GO-resource series without unsafe early `part_of` edges. |

The fixed-label series is the appropriate place for GAF-only years. It does not
claim to reconstruct historical phenotypes. Its purpose is to isolate changing
GO resources while holding the label collection fixed.

## Paper-output Design Preview

[`mock_publication_outputs/README.md`](mock_publication_outputs/README.md) is a
pre-run design review containing eight synthetic tables and eight publication-style
figure panels. Every row and image is marked `SYNTHETIC MOCK - NOT RESULTS`.
The preview proposes the manuscript order and visual grammar; it is not an analysis
of the downloaded files and none of its numerical values may be cited.

## Organism-Specific Labels

- **Yeast:** null-mutant observations explicitly reporting viable/non-essential
  or inviable/essential phenotypes; mixed genes are excluded.
- **Fly:** explicit viable or lethal phenotype text; partial lethality, mixed
  genes, and conservative multi-gene records are excluded. D. melanogaster
  assignments are derived from the release GAF using `taxon:7227`.
- **Mouse:** annual IMPC viability calls are converted into paired lethal and
  viable sets. Subviable and conflicting genes are excluded and reported.
- **Fish:** lethal evidence wins over other observations. Genes with only other
  observed phenotype evidence form an operational nonlethal class, not a proven
  viable class.
- **Worm:** lethal terms and descendants are regenerated from each annual
  WormBase phenotype ontology. Genes observed only outside the lethal term set
  form an operational nonlethal class, not an explicit viability assay.

These differences are recorded in every manifest. Cross-organism values should
not be interpreted as if the phenotype evidence had one universal semantics.
The primary statistical comparisons are within organism across resource years.

## Pipeline Stages

1. Validate the ledger, required input paths, ARM64 interpreter, compressed-file
   integrity, scientific Python imports, and the complete test suite.
2. Record the exact command stream, Git state, package environment, machine
   information, and SHA-256 checksums of every selected input.
3. Derive release-specific mouse paired labels, fly taxon-filtered assignments,
   and worm lethal phenotype descendants.
4. Build strict ARFF datasets with stable identifiers and full manifests for
   `primary`, `is_a_only`, and `no_iea`.
5. Export the 2025 labels from each primary ARFF and rebuild the 2014-2025
   `fixed_2025` and 2010-2025 `fixed_2025_is_a` annotation-only series.
6. Validate every ARFF for schema, binary values, duplicate rows/features, and
   class integrity.
7. Run the legacy single-snapshot model suite on each organism's 2025 dataset,
   including the TensorFlow network, as a supplementary implementation check.
8. Run repeated within-year and cross-year analyses for decision tree, random
   forest, gradient boosting, logistic regression, SVM, and neural network.
9. Run temporal GO enrichment and GO-term timeline analysis for each track and
   organism.
10. Aggregate paper-facing tables, generate HTML reports, hash all outputs, and
    write the completion marker.

## Comparable Model Design

The primary `phengo-version-sensitivity` analysis uses six model families:

```text
dt  decision tree
rf  random forest
gb  gradient boosting
lr  logistic regression
svm support vector machine
nn  scikit-learn multilayer perceptron
```

All six use the same repeated gene folds, four data panels, train-year/test-year
directions, seeds, lethal-positive mapping, 0.5 threshold, metrics, and nested
sigmoid calibration. Gradient boosting and the MLP receive class-balanced
sample weights; the remaining estimators use their class-weight mechanisms.

The MLP defaults to two ReLU hidden layers of 128 and 64 units, Adam with a
0.001 initial learning rate, L2 alpha 0.0001, batch size 32, 300 maximum
iterations, and training-only early stopping. Its iteration and convergence
diagnostics are emitted with the result tables. Any convergence failure should
be handled by changing the global pre-specified limit and rerunning every year,
not by tuning an individual release.

Cross-year models use only terms represented in the training release. A GO term
that first appears in the test year cannot have a learned coefficient, split,
or neural-network weight and is therefore ignored for that transfer direction.
One fixed gene-ID fold map is reused across all test years in a repeat. This
allows each train-year/fold model to score every test year while ensuring that a
scored gene is removed from training whenever it occurs in the training
snapshot.

The TensorFlow network under `04_single_snapshot_ml` uses a separate held-out
test/validation design. It is not directly comparable to the repeated temporal
models and should be reported only as a supplementary implementation analysis.
It runs on 2025 by default because the repeated engine already supplies the
within-year result for every release. Set `PHENGO_SINGLE_SNAPSHOT_YEARS=all` or
a space-separated year list only when that additional implementation experiment
has been pre-specified.

## Evaluation Panels

| Panel | Genes | GO features | Main interpretation |
| --- | --- | --- | --- |
| `full` | Native to each release | Training release | Practical consequence of deploying a model trained on that release. |
| `matched_features` | Native to each release | Shared across all years | Controls GO feature-set churn. |
| `matched_genes` | Shared across all years | Training release | Controls gene cohort churn. |
| `matched_both` | Shared across all years | Shared across all years | Strongest controlled version-sensitivity panel. |

The primary controlled result should be `matched_both`; `full` is the practical
effect-size panel. The difference between panels helps distinguish gene-set and
feature-set churn from changes among the same genes and terms.

## Output Layout

```text
New_Outputs/
  00_run_metadata/
    commands.tsv
    environment.txt
    git_commit.txt
    git_status.txt
    input_checksums.sha256
    optional_anchor_status.tsv
    publication_snapshots.tsv
    python_packages.txt
    run_config.env
    run_fingerprint.sha256
    output_checksums.sha256
    run.complete
  01_derived_labels/
    fixed_2025/
    fly/
    mouse/
    worm/
  02_arff/{track}/{organism}/{year}/
    {organism}_PhenGO.arff
    PhenGO_manifest.json
    PhenGO_params.txt
    gene_identifiers.tsv
    identifier_join_audit.tsv
    label_source_audit.tsv
    go_enrichment/
  03_qc/{track}/{organism}/{year}/
  04_ml/{track}/{organism}/
    evaluation_preflight.csv
    dataset_drift.csv
    pairwise_drift.csv
    within_year_cv.csv
    within_year_cv_summary.csv
    cross_year_transfer.csv
    cross_year_transfer_summary.csv
    prediction_instability.csv
    prediction_instability_summary.csv
    previous_available_snapshot_label_baseline.csv
    feature_rank_overlap.csv
    feature_importance_stability.csv
    transfer_matrices/
    version_sensitivity_manifest.json
  04_single_snapshot_ml/{organism}/{year}/Predict/
  05_temporal/{track}/{organism}/
  06_publication_tables/
  07_reports/
  logs/
```

`PhenGO_manifest.json` is part of each ARFF dataset. It identifies every source
release and input hash, all label and GO policies, software and dependency
versions, architecture, source-tree hash, and output hash. The analysis command
rejects mixed manifests by default.

## Paper-Facing Tables

`06_publication_tables/` contains combined TSVs with explicit track and organism
columns:

- `snapshot_inventory.tsv`: class, feature, release, policy, and checksum spine;
- `within_year_model_performance.tsv`: repeated OOF performance by year;
- `cross_year_transfer_performance.tsv`: all directional transfer estimates;
- `dataset_drift.tsv` and `pairwise_dataset_drift.tsv`: source-resource drift;
- `prediction_instability_summary.tsv` and
  `prediction_instability_by_gene.tsv`: probability movement and prediction
  flips for the same genes;
- `feature_rank_overlap.tsv` and `feature_importance_stability.tsv`: stability
  of learned GO associations for models with native importance;
- `previous_available_snapshot_label_baseline.tsv`: label-only baseline against
  the preceding available snapshot, with the calendar interval reported;
- `single_snapshot_model_performance.tsv`: supplementary held-out model results;
- `temporal_summary.tsv` and `temporal_enrichment_statistics.tsv`: GO timeline
  and complete enrichment statistics;
- `evaluation_preflight.tsv`: proof that the requested design was feasible.

Use the summary tables and transfer matrices for figures. Keep fold-level and
per-gene files as archived supplementary data rather than treating folds as
independent biological replicates.

## Paper-Ready Methods Text

> Archived phenotype, Gene Ontology annotation, Gene Ontology structure, and,
> where required, phenotype-ontology releases were organized in a pre-specified
> organism-year ledger. Release-matched phenotype and GO datasets were generated
> for fish, fly, mouse and worm from 2015 to 2025 and for yeast from 2015 to 2017
> and 2024 to 2025. Each dataset was generated under an identical software and
> policy contract and was accompanied by release identifiers, SHA-256 input and
> output hashes, stable gene identifiers, class audits, software dependencies,
> and source-code provenance. Positive GO annotations were propagated through
> `is_a` and `part_of` relationships using the corresponding archived GO
> ontology; ontology roots, terms annotated to fewer than two genes, and terms
> annotated to at least 99% of genes were excluded.
>
> Database-version sensitivity was evaluated separately within each organism.
> Decision tree, random forest, gradient boosting, logistic regression, support
> vector machine, and feed-forward multilayer-perceptron classifiers were used as
> standardized probes rather than year-specifically optimized predictors. The
> same fixed gene-identifier folds, random seeds, feature/cohort panels,
> class-balancing policy, training-fold probability calibration, lethal-positive
> class mapping, threshold, and metrics were used for all model families. For
> each train-year/test-year direction, genes assigned to the evaluated fold were
> removed from the training release whenever present, producing gene-disjoint
> out-of-fold estimates on matrix diagonals and off-diagonals. Models used only
> the GO vocabulary available in their training release.
>
> Results were summarized across five repeated five-fold partitions using
> average precision, ROC AUC, Matthews correlation coefficient, balanced
> accuracy, lethal-class F1, precision and recall, macro F1, Brier score, and
> class prevalence. Full-release analyses were compared with globally matched
> gene, matched GO-feature, and matched-gene-and-feature panels. A fixed-label
> experiment reused the 2025 phenotype labels with each archived GAF and GO
> ontology to separate GO-resource drift from phenotype-label drift. Majority,
> class-prior, stratified-random, and previous-available-snapshot label
> baselines were retained. Calendar intervals were reported explicitly so that
> adjacent available observations were not assumed to represent consecutive
> years.
> Per-gene probability ranges, prediction flips, resource Jaccard similarity,
> label churn, and bootstrap GO feature-rank overlap were reported as descriptive
> measures of instability.

The organism-specific evidence definitions, especially the operational
nonlethal fish and worm classes, must be stated immediately after this common
methods paragraph. The complete limitations text is in
[`manuscript_limitations.md`](manuscript_limitations.md).

## Result Interpretation

1. Establish that the resources changed using snapshot counts, prevalence,
   sparsity, gene Jaccard, GO Jaccard, and shared-gene label churn.
2. Show ordinary within-year performance and uncertainty without selecting a
   favorite year after seeing the results.
3. Present directional transfer heatmaps. Compare each diagonal with transfers
   into earlier and later releases.
4. Contrast `full` with `matched_both`, then use the intermediate panels to
   attribute instability to gene or GO feature churn.
5. Compare `primary` with `fixed_2025` over 2015-2025 to distinguish changing
   phenotype targets from changing GO annotations and ontology structure.
6. Use `fixed_2025_is_a` for the full 2010-2025 historical annotation series;
   do not pool it with the `is_a + part_of` track as though feature policies match.
7. Use `is_a_only` and `no_iea` as pre-declared sensitivity analyses.
8. Report whether conclusions persist across simple and complex model families,
   including the comparable MLP, and report all planned models even when they do
   not support the narrative.
9. Inspect convergence-warning columns before interpreting the MLP and retain
   the evaluation preflight as supplementary evidence.

Differences should be described as sensitivity to the complete database
snapshot, which includes curation policy, evidence availability, ontology
structure, identifier mappings, and phenotype coverage. They do not establish
that one archived release is biologically correct and another is incorrect.

## Completion Gate

The run is complete only when all of the following hold:

- `00_run_metadata/run.complete` exists;
- required ledger rows all have strict ARFF manifests and passing QC;
- every selected track/organism with at least two snapshots has an analysis
  manifest and evaluation preflight;
- no primary MLP convergence warnings remain unresolved;
- input and output checksum files exist;
- `run_fingerprint.sha256` exists and corresponds to the completed run;
- combined publication tables and HTML reports were generated;
- annotation-anchor omissions are explicitly reported in exploratory partial runs;
- the Git state, source-tree hash, package environment, and exact command log
  are archived with the results.

The full default analysis is intentionally computationally intensive. It can be
interrupted and continued with `--resume`; changed code, dependencies, inputs,
ledger content, or parameters invalidate cached units automatically. Keep the
completed output directory immutable and copy it as the article analysis
snapshot before generating revised experiments.
