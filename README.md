[![DOI](https://zenodo.org/badge/1026270084.svg)](https://doi.org/10.5281/zenodo.16748708)

# PhenGO

PhenGO builds traceable machine-learning datasets from model-organism phenotype
records, Gene Ontology annotations, and a release-matched GO ontology. It
supports *Saccharomyces cerevisiae*, *Mus musculus*, *Drosophila melanogaster*,
*Caenorhabditis elegans*, and *Danio rerio*.

The main research use is to test whether ML conclusions change when a different
year or release of a model-organism database is used. PhenGO therefore treats
dataset provenance, label polarity, stable identifiers, GO semantics, and
gene-disjoint evaluation as part of the analysis rather than incidental IO.

## Installation

```bash
pip install phengo
```

Install the scientific prediction and report dependencies with:

```bash
pip install "phengo[predict,reports]"
```

TensorFlow is optional:

```bash
pip install "phengo[nn]"
```

For a development checkout:

```bash
python -m pip install -e ".[predict,reports,test]"
```

Python 3.10 or newer is required.

## Label Definition

The ML class mapping is explicit and fixed:

```text
viable = 0
lethal = 1
```

`inviable` and `essential` are accepted lethal aliases. `non-essential` is an
accepted viable alias. Mixed lethal/viable records and partial/semi-lethal
records are excluded by default. Other observed phenotypes are not assumed to
be viable unless `-nonlethal_policy observed_viable` is selected deliberately.

### Label Sources

The publication default is `-label_source release_records`: lethal and viable
labels are computed from the phenotype file belonging to each database
snapshot. This measures the practical combined effect of phenotype and GO
resource version.

Two paired-collection modes are available for controlled sensitivity analyses:

```bash
-label_source gene_sets -lethal_gene_set lethal.tsv -viable_gene_set viable.tsv
-label_source agreement -lethal_gene_set lethal.tsv -viable_gene_set viable.tsv
```

`gene_sets` uses the paired collections as the complete target definition;
`-phenotype_file` is optional and is not used for classification. `agreement`
requires a phenotype file and retains only genes with identical release-derived
and collection-derived labels. Both collection files are mandatory. Their first
tab-separated column is authoritative and should contain a stable accession;
the optional second column is a display symbol and is used for matching only
when the stable-accession cell is blank. Genes absent from both files are unknown
and excluded; absence from the lethal collection never implies viability.

For fly exploratory work, `-lethal_gene_set bundled` selects the dated 2025
collection packaged with PhenGO, but a bundled current collection cannot satisfy
`-strict_snapshot`. The old `-fly_lethal_genes` lethal-only override is
deprecated and prohibited in strict publication runs.

## Generate A Dataset

Publication runs should use `-strict_snapshot` and release-matched auxiliary
files. Example for WormBase phenotype-term mode:

```bash
PhenGO \
  -species worm \
  -phenotype_file inputs/WS257/phenotype_association.tsv.gz \
  -gene_association_file inputs/WS257/wb.gaf.gz \
  -go_obo_file inputs/WS257/go-basic.obo.gz \
  -worm_phenotypes inputs/WS257/lethal_terms.tsv.gz \
  -worm_viable_phenotypes inputs/WS257/viable_terms.tsv.gz \
  -output_dir results/worm_WS257 \
  -snapshot_id worm_WS257 \
  -phenotype_release WS257 \
  -go_annotation_release WS257 \
  -go_ontology_release WS257 \
  -phenotype_ontology_release WS257 \
  -retrieval_date 2026-07-17 \
  -strict_snapshot \
  -go_relations is_a part_of
```

Key outputs are:

```text
worm_PhenGO.arff
gene_identifiers.tsv
identifier_join_audit.tsv
label_source_audit.tsv
PhenGO_manifest.json
PhenGO_params.txt
GO_term_validation_report.txt
GO_Children&Parents.json
go_enrichment/
```

ARFF rows use stable database accessions. `gene_identifiers.tsv` preserves the
corresponding display symbols. `identifier_join_audit.tsv` records every
unmatched, ambiguous, contradictory, or lower-priority identifier decision.
Strict snapshots disable GAF synonym fallback; canonical-symbol fallback is
accepted only when it resolves to exactly one stable accession.

## GO Feature Policy

The biological and CLI default uses the two GO relations that safely support
positive annotation grouping:

```bash
-go_relations is_a part_of
```

Use a release-matched `go-basic.obo`. PhenGO rejects cyclic selected-relation
graphs and, by default, cross-namespace propagation edges. The latter normally
indicate that a full ontology rather than `go-basic` was supplied. An explicitly
designed cross-aspect sensitivity analysis can opt in with:

```bash
-allow_cross_namespace_go_edges
```

Available relation options are:

```text
is_a
part_of
regulates
positively_regulates
negatively_regulates
occurs_in
capable_of
capable_of_part_of
```

Other GO controls include:

```text
-go_namespaces biological_process molecular_function cellular_component
-go_propagation ancestors|direct
-go_evidence_codes IMP IDA ...
-exclude_go_evidence_codes ND IEA ...
-min_go_gene_count 2
-max_go_gene_fraction 0.99
-include_go_roots
```

By default, canonical GO roots are excluded, annotations propagate to selected
ancestors, alternate IDs are canonicalized, unambiguous obsolete
`replaced_by` mappings are applied, and terms present in fewer than two genes
or more than 99% of genes are removed. Keep the complete GO policy identical
across snapshots. Regulation, occurrence, and capability relations are not
primary defaults because collapsing them into an ordinary binary GO feature
changes the biological meaning of the gene-term association.

## Snapshot Provenance

`PhenGO_manifest.json` records:

- snapshot and release identifiers;
- retrieval date;
- absolute input paths, sizes, timestamps, and SHA-256 hashes;
- label, GO relation, namespace, propagation, and evidence policies;
- Python executable, machine architecture, dependency versions, tool version,
  git commit, dirty-worktree status, source-tree hash, and command;
- output hashes and dataset counts.

Manifest schema 3 records the label source, paired collection hashes, label and
identifier-join audit counts, identifier precedence, OBO header, GO graph
invariants, all feature/phenotype policy settings, and the availability and
retrieval context of each independently released resource. Schema-2 V1
manifests remain readable. Publication version analysis requires every snapshot
within one analysis to have the same tool version, source-tree hash, dependency
environment, and analysis-relevant policy values.

Bundled helper resources are useful for current exploratory runs but are not
appropriate for historical analyses. Strict mode requires explicitly supplied,
release-matched helper resources for fly, worm, and mouse.
Under the default explicit-viability policy, worm term mode additionally needs
`-worm_viable_phenotypes`, worm gene-list mode needs `-worm_viable_genes`, and
mouse needs `-mouse_viable_phenotypes`.

## Version-Sensitivity Analysis

Regenerate all publication ARFFs with the current core before mixing releases.
Then run one analysis per organism:

```bash
phengo-version-sensitivity \
  -arff_files \
    results/worm_WS247/worm_PhenGO.arff \
    results/worm_WS257/worm_PhenGO.arff \
    results/worm_WS297/worm_PhenGO.arff \
  -dataset_names WS247 WS257 WS297 \
  -models all \
  -panels full matched_features matched_genes matched_both \
  -cv_folds 5 \
  -cv_repeats 5 \
  -calibration sigmoid \
  -importance_repeats 20 \
  -seed 42 \
  -output_dir paper_outputs/worm_version_sensitivity
```

The command requires strict manifests by default and verifies each ARFF hash.
Before fitting, `evaluation_preflight.csv` verifies that every requested panel,
transfer direction, fold, and calibration split is feasible. The command fails
instead of silently reducing folds or omitting cells. The
`-allow_incomplete_evaluation` override is for explicitly exploratory work only.
It reports:

- dataset, gene-set, GO-feature, prevalence, sparsity, and label drift;
- repeated within-snapshot out-of-fold performance;
- gene-disjoint train-release to test-release transfer;
- full, matched-feature, matched-gene, and matched-both panels;
- majority, class-prior, random, and previous-available-snapshot label baselines;
- calibrated ROC-AUC, average precision, balanced accuracy, MCC, F1, precision,
  recall, and Brier score;
- fixed-fold per-gene probability and prediction instability;
- bootstrap GO feature-importance stability and between-release overlap;
- repeat-level standard deviations and 95% resampling interval half-widths.

`all` runs decision tree, random forest, gradient boosting, logistic regression,
SVM, and a scikit-learn multilayer perceptron (`nn`). The MLP uses the identical
outer folds, year-to-year transfers, panels, seeds, class balancing,
training-fold calibration, 0.5 lethal-positive threshold, and metrics as the
other estimators. Its fixed architecture and convergence diagnostics are stored
in the analysis manifest and result tables. This is the neural-network result
to use in primary cross-model publication comparisons.

Use `cross_year_transfer_summary.csv` and `transfer_matrices/` for primary
transfer figures. `cross_year_transfer.csv` is repeat-level diagnostic output.

The complete article workflow is in
[`docs/version_sensitivity_workflow.md`](docs/version_sensitivity_workflow.md).
The standalone, one-command dataset-to-paper workflow is described in
[`docs/publication_pipeline_v2.md`](docs/publication_pipeline_v2.md) and
implemented by
[`scripts/run_publication_pipeline_v2.sh`](scripts/run_publication_pipeline_v2.sh).
It builds, resumes, or cleanly repairs the base publication run, adds all
alternative-resource analyses, audits source/manifests/ARFFs before reuse, and
creates one self-contained consolidated result directory. A future run does
not require a separate V1 command.
[`scripts/run_publication_pipeline.sh`](scripts/run_publication_pipeline.sh) is
the lower-level base-run phase and normally does not need to be called directly.
A manuscript-ready limitations section is in
[`docs/manuscript_limitations.md`](docs/manuscript_limitations.md).

## Prediction Models

Single-dataset prediction defaults to logistic regression:

```bash
phengo-predict \
  -arff_file results/worm_WS257/worm_PhenGO.arff \
  -model lr rf gb dt \
  -output_dir predictions/worm_WS257
```

Available models are `lr`, `rf`, `gb`, `dt`, `svm`, and optional TensorFlow
`nn` in this single-snapshot convenience command. Sklearn probabilities use
nested training-fold calibration by default. Saved models include
`model_schema.json`, which fixes feature order, class mapping, positive class,
threshold, and calibration policy. Neural-network thresholds are selected on a
held-out validation subset rather than the training or test set. This
TensorFlow result is a supplementary implementation check; use the
scikit-learn `nn` in `phengo-version-sensitivity` for directly comparable
repeated temporal evaluation.

Multi-dataset comparison mode now honors every requested model:

```bash
phengo-predict \
  -arff_files release_a.arff release_b.arff \
  -dataset_names A B \
  -model lr rf \
  -output_dir comparison
```

Use `phengo-version-sensitivity`, not this convenience comparison mode, for
paper-level temporal inference.

## Temporal GO Analysis

```bash
temporal-analysis \
  -input_dir results/worm_snapshots \
  -output_dir temporal/worm \
  -o worm \
  --top-n 20 \
  --max-fdr 0.05
```

GO enrichment uses Fisher's exact test, Haldane-Anscombe corrected finite
ratios, and Benjamini-Hochberg FDR. All tested terms are retained in
`*_enrichment_statistics.csv`; ranked timelines apply the requested FDR limit.
Strict source manifests and matching ARFF hashes are required by default, and
the analysis writes `temporal_analysis_manifest.json` with input/output hashes.

## Compare ARFF Files

```bash
compare-arff \
  -arff_a release_a.arff \
  -arff_b release_b.arff \
  -o comparison.csv
```

The comparison is symmetric: it reports genes missing in A or B, label
mismatches, GO mismatches, and exact matches. Each gene appears once even when
it has multiple mismatch types. Conflicting duplicate ARFF rows are rejected.

## Safety And Reproducibility

- Non-empty output directories are refused unless `-overwrite` is supplied.
- Output directories that equal or contain protected inputs are rejected.
- Malformed ARFF rows, duplicate features, conflicting duplicate genes,
  non-binary values, and schema mismatches are errors rather than silently
  dropped or converted to zero.
- Use `-n_jobs 1` for the most portable deterministic execution; it is the
  default in the ML commands.
- Set and report `-seed`; the paper workflow defaults to `42`.

Run `PhenGO -h`, `phengo-predict -h`, or
`phengo-version-sensitivity -h` for the full command reference. Run
`PhenGO --print-defaults` to inspect bundled helper resources.
