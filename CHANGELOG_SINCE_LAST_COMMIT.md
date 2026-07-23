# Changelog Since Last Git Commit

Base commit:

```text
0c06311 Major toolkit-wide bug and logic fixes. Will be documented in future
```

This document describes the current uncommitted changes. They complete a
scientific hardening pass of the dataset builder, GO feature construction,
predictors, temporal analysis, provenance system, and version-sensitivity
workflow.

The package version is now `0.3.0`. V1 datasets use schema-2 manifests; newly
generated V2 datasets use schema 3 to add resource-availability and composite-
snapshot context. Both remain readable, but neither should be mixed with
unmanifested `0.2.0` publication datasets.

## Summary

The main behavioral changes are:

1. Publication ARFF generation can now be made release-strict and hash-traceable.
2. ARFF rows now use stable database accessions rather than mutable symbols.
3. Conservative label policies are explicit and are now the defaults.
4. GO parsing is relation-aware, canonicalizes alternate/obsolete IDs, and supports user-selected relations.
5. Version-sensitivity predictions are gene-disjoint and out-of-fold, including transfer diagonals and per-gene instability.
6. Probabilities are calibrated inside training folds, and uncertainty is calculated across repeated complete OOF predictions.
7. A class-balanced MLP neural network now participates in the same repeated,
   gene-disjoint version-sensitivity design as the other model families.
8. Temporal GO enrichment now has exact tests, finite corrected effect sizes, and FDR correction.
9. Destructive output handling, malformed ARFFs, duplicate records, and schema mismatches now fail explicitly.
10. A manuscript-ready workflow and limitations section document what code cannot resolve biologically.
11. Release-derived, paired gene-set, and agreement label designs are now separate and auditable.
12. Publication analysis now enforces identical source code, environments, policies, folds, and calibration feasibility across snapshots.
13. V2 adds fail-closed Fly compound parsing, corrected Worm lethal roots,
    alternative MOD and historical mouse inputs, explicit archive-availability
    context, interval-aware temporal results, and a verified unified output.

These changes alter row identifiers, labels, GO features, and provenance.
Publication ARFF files must be regenerated; legacy and regenerated ARFFs should
not be mixed in one analysis.

## Core Dataset Builder

### Explicit target policies

File: `src/PhenGO/core/PhenGO.py`

New options:

```text
-mixed_label_policy lethal_wins|exclude|error
-nonlethal_policy observed_viable|explicit_only
-ambiguous_label_policy exclude|viable
```

New defaults:

```text
mixed observations       exclude
nonlethal observations   require explicit viability
partial/semi-lethal      exclude
fly multi-gene records   filter
```

`-filter_mixed_terms` remains as a deprecated alias for
`-mixed_label_policy exclude`. `-allow_multigenes` permits the previous fly
behavior explicitly.

### Explicit label-source designs

New options:

```text
-label_source release_records|gene_sets|agreement
-lethal_gene_set
-viable_gene_set
```

`release_records` is the default and computes labels from each supplied MOD
phenotype release. `gene_sets` uses paired lethal and viable collections as the
complete target definition. `agreement` retains only stable genes for which the
release-derived and collection-derived labels agree. Unlisted genes remain
unknown rather than being inferred viable, conflicting identifiers are errors,
and paired files use one authoritative identifier per row. The first column is
the stable accession; the second is a display-symbol fallback only when the
stable-accession cell is blank.

Every run writes `label_source_audit.tsv` with release and collection labels and
an explicit retained/disagreement/source-only outcome. The manifest records the
label source, both collection hashes, source/join counts, and audit counts. The
bundled dated fly lethal collection is available with
`-lethal_gene_set bundled` outside strict mode. The old `-fly_lethal_genes`
lethal-only override remains for exploratory compatibility but is deprecated
and prohibited in strict publication snapshots.

### Phenotype parser corrections

File: `src/PhenGO/core/phenotype_handling.py`

- Centralised repeated-observation resolution.
- Worm `NOT` and evidence fields are tokenized exactly.
- Worm rows excluded by qualifier/evidence policy are skipped; they are no
  longer converted into viable examples.
- Worm phenotype terms are matched as exact tokens rather than substrings.
- Worm lethal-gene-list mode respects `nonlethal_policy`.
- Worm term mode can consume `-worm_viable_phenotypes`; worm gene-list mode can
  consume `-worm_viable_genes`.
- Fish semi-lethal and semi-viable records are excluded by default.
- Fish lethal/dead and viable/alive terms use word-boundary matching.
- Fly partial lethality follows the ambiguous-label policy.
- Yeast, fly, fish, worm, and mouse nonlethal records all follow the same
  explicit policy interface.
- Mouse can consume `-mouse_viable_phenotypes` for explicit negative labels.
- Core generation stops unless both viable and lethal classes remain.
- Mixed labels can be excluded, resolved with lethal priority for legacy
  reproduction, or treated as an error.

### Stable gene identifiers

File: `src/PhenGO/core/go_handling.py`

- GAF DB Object ID is now the canonical ARFF row key.
- GAF joins index the complete filtered annotation file before assigning labels.
- Stable accessions take precedence over canonical symbols; non-stable
  identifiers are accepted only when they resolve to exactly one stable GAF ID.
- Strict snapshots disable historical-synonym fallback. Exploratory runs may use
  a synonym only when it resolves uniquely and no stronger identifier matches.
- Output records retain the current gene symbol for a separate identifier map.
- Contradictory labels at the strongest identifier priority are excluded rather
  than resolved by lethal priority. Lower-priority contradictions are ignored
  and audited without overriding the stronger assignment.
- Fish, worm, and yeast phenotype parsers now retain their source ZFIN,
  WormBase, and SGD stable accessions instead of discarding them for symbols.
- Release-matched fly symbols are translated to unique FBgn accessions before
  the GAF join; ambiguous assignment symbols are excluded.
- `identifier_join_audit.tsv` records unmatched identifiers, one-to-many symbol
  mappings, contradictory labels, and precedence decisions. Counts and policy
  are included in `PhenGO_manifest.json`.
- Exact `NOT` qualifier filtering replaces substring matching.
- User-selected GAF evidence codes are applied consistently across organisms.
- Duplicate GO annotations per gene are removed.
- Deprecated fly lethal-only overrides are applied when requested using one
  authoritative identifier per input row, but paired label sources are the
  publication interface.
- Fly assignment generation excludes any release symbol associated with more
  than one FBgn accession and writes the rejected symbol-to-ID mappings to an
  `excluded.tsv` sidecar.
- `gene_identifiers.tsv` records stable accession to display-symbol mappings.

An all-year primary-input audit confirmed stable source identifiers on every
eligible fish, worm, and yeast row and all IMPC rows from 2016 onward. The 2015
IMPC source contains symbols only; 1,208 of 1,212 labels resolve through unique
exact 2015 GAF canonical symbols, four remain unmatched, and none use synonyms.
Later FlyBase GAFs contain a small number of one-to-many canonical symbols;
these are now explicitly excluded during release-assignment generation.

### Strict snapshot provenance

Files:

```text
src/PhenGO/provenance.py
src/PhenGO/core/PhenGO.py
```

New core options:

```text
-snapshot_id
-phenotype_release
-go_annotation_release
-go_ontology_release
-phenotype_ontology_release
-retrieval_date
-strict_snapshot
```

`-strict_snapshot` requires complete relevant release metadata and explicit,
release-matched helper resources for fly, worm, and mouse. `retrieval_date`
must use `YYYY-MM-DD`. Bundled current helper files emit warnings outside
strict mode and cannot satisfy a strict historical run merely through a
`default` alias.

New runs write schema-3 `PhenGO_manifest.json` (schema-2 V1 manifests remain
readable), containing:

- tool version and git commit;
- dirty-worktree status and a hash of the complete Python source tree;
- command line;
- Python executable, platform, architecture, and dependency versions;
- release identifiers and retrieval date;
- resource availability, retrieval route, and composite-snapshot semantics;
- absolute input paths, sizes, modification times, and SHA-256 hashes;
- label and GO policies;
- output hashes;
- gene, class, and feature counts.

### Output safety

File: `src/PhenGO/provenance.py`

- Root, home, current working directory, and input-containing output paths are rejected.
- Non-empty output directories require `-overwrite`.
- Core log files can be preserved while other output is replaced.
- Auxiliary inputs are protected as well as the three primary inputs.
- Empty post-filter gene cohorts or GO feature sets now stop the run rather than writing unusable ARFFs.
- `PhenGO --print-defaults` works without supplying the otherwise required run arguments.
- Publication `--resume` now fingerprints source, dependencies, pipeline and
  ledger content, input hashes, and analysis parameters. Any change forces a
  complete regeneration instead of mixing stale and current outputs.

## GO Ontology And Features

### Relation-aware OBO parser

File: `src/PhenGO/core/obo_to_graph.py`

The former line-oriented GO parser was replaced by a stanza-aware parser. It
supports:

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

`is_a + part_of` is the default biological grouping policy:

```bash
-go_relations is_a part_of
```

The graph records relation-specific parents and children in
`GO_Children&Parents.json`. Alternate IDs map to one canonical feature rather
than becoming duplicate features. Obsolete terms are followed through
unambiguous `replaced_by` chains; unresolved or ambiguous obsolete terms are
removed and reported. Duplicate alternate-ID mappings are errors.

### GO controls

New options:

```text
-go_namespaces biological_process molecular_function cellular_component
-go_propagation ancestors|direct
-exclude_go_roots
-include_go_roots
-go_evidence_codes ...
-exclude_go_evidence_codes ...
-min_go_gene_count
-max_go_gene_fraction
```

Defaults are ancestor propagation over `is_a + part_of`, exclusion of the three
canonical GO roots, minimum occurrence in two genes, and maximum prevalence of 99%.
Include/exclude evidence-code contradictions are rejected.

Selected relation graphs must be acyclic. Cross-namespace propagation is
rejected by default, which protects publication runs from accidentally using
unsafe edges from a full ontology in place of release-matched `go-basic.obo`.
`-allow_cross_namespace_go_edges` is an explicit sensitivity-analysis override.
The OBO format, ontology, and data-version header fields are stored in the
manifest when present.

### Feature consistency

File: `src/PhenGO/core/PhenGO.py`

- Direct annotations, expanded annotations, and binary vectors remain distinct.
- GO validation canonicalizes alternate IDs and obsolete replacements before vectorization.
- Genes left without a valid vector are removed even when unused-feature filtering is disabled.
- Unused and prevalence filters rebuild every gene vector without discarding identifier metadata.
- ARFF and enrichment outputs use the same expanded feature universe.
- ARFF writing canonicalizes booleans and other valid binary scalars to the
  declared integer `0/1` domain and rejects invalid vector lengths, values,
  labels, or duplicate feature declarations.
- Root exclusion now targets the three canonical GO roots; it does not remove every node lacking a selected-relation parent.

### GO mapping utility

The public `phengo-tools go-map` parser-state bug that could yield an empty
ontology has been fixed. It now reads plain or gzipped OBO files, rejects
duplicate term stanzas and alternate IDs mapped to multiple terms, canonicalizes
alternate IDs, follows only a single direct active `replaced_by` target, and
reports `exact`, `alt_id`, `obsolete_replaced`, `obsolete_unresolved`, or
`missing` for every requested GO ID. Ambiguous obsolete replacements remain
unmapped rather than being guessed.

## Version-Sensitivity ML

File: `src/PhenGO/predict/version_sensitivity.py`

### Comparable neural-network evaluation

- `nn` now selects scikit-learn `MLPClassifier` in the version-sensitivity
  command and is included by `-models all`.
- The MLP uses the same outer folds, year-to-year transfer cells, feature and
  cohort panels, random-seed policy, balanced training weights, nested
  calibration, lethal-positive mapping, threshold, and metrics as the other
  classifiers.
- Architecture, L2 penalty, batch size, learning rate, maximum iterations,
  early-stopping validation fraction, and patience are CLI parameters and are
  captured in `version_sensitivity_manifest.json`.
- Fold and transfer outputs record convergence-warning counts and maximum
  optimizer iterations. Publication runs should be repeated globally if the
  selected iteration limit does not converge.
- Cross-year models now use only the training release's GO vocabulary. This
  prevents test-only terms from contributing untrained random MLP weights or
  altering random-forest feature sampling.
- One fixed gene-ID fold partition is reused across all transfer directions in
  each repeat. A fitted train-year/fold model can therefore score every test
  year while preserving gene-disjoint evaluation, substantially reducing the
  number of redundant fits.
- `scikit-learn>=1.7` is required so MLP training can consume class-balanced
  sample weights.
- The existing TensorFlow network remains available only in the
  single-snapshot runner and is described as a supplementary implementation
  sensitivity analysis.

### Manifest enforcement

- Strict `PhenGO_manifest.json` files are required by default.
- ARFF hashes are checked against manifests.
- All datasets must represent one organism and have unique snapshot IDs.
- Release metadata requirements account for organism-specific auxiliary resources.
- `-allow_missing_manifests` is available only for explicit exploratory legacy analysis.
- `version_sensitivity_manifest.json` records the analysis command,
  dependencies, input hashes, source snapshot IDs, and evaluation design.
- Schema-2 publication manifests must have identical tool versions,
  source-tree hashes, dependency environments, label/GO policies, parser
  settings, and feature filters across snapshots.

### Cohorts and duplicate validation

- Duplicate GO attributes are errors.
- Missing, non-numeric, non-binary, or conflicting duplicate gene rows are errors.
- Exact duplicate rows are deduplicated with a warning.
- `matched_genes` now means genes shared across every snapshot, not a different
  pairwise cohort for every transfer cell.
- Matched panels reuse fixed gene-index partitions across releases.
- Empty matched feature panels are skipped explicitly rather than crashing.

### Out-of-fold evaluation

- Within-snapshot summaries are computed from complete OOF predictions for each repeat.
- Cross-snapshot transfer removes each scored test gene from that fold's training data.
- Diagonal transfer cells are genuine OOF estimates, not resubstitution scores.
- Per-gene prediction-instability probabilities are fixed-fold, gene-disjoint OOF scores.
- Output rows state the evaluation design explicitly.

### Calibration, metrics, and uncertainty

New options:

```text
-calibration none|sigmoid|isotonic
-calibration_cv
-cv_repeats
-n_jobs
-seed
```

Sigmoid calibration is the default and is fitted only inside outer training
data. Reported metrics include ROC-AUC, average precision, balanced accuracy,
MCC, F1-lethal, macro F1, lethal precision/recall, accuracy, Brier score,
prevalence, and sample counts.

Fold metrics remain available for diagnostics, but summaries first combine all
OOF predictions within a repeat. Standard deviations and Student-t 95%
interval half-widths are calculated across repeats. Gradient boosting receives
balanced sample weights. Internal model parallelism defaults to one worker for
portable deterministic execution.

Before fitting, `evaluation_preflight.csv` verifies every selected panel and
transfer direction can use the requested folds and calibration CV without
class-deficient training data. Publication mode fails rather than silently
reducing fold counts or omitting cells. `-allow_incomplete_evaluation` retains
adaptive behavior for explicitly exploratory analyses.

### Transfer outputs

New or changed outputs:

```text
cross_year_transfer.csv
cross_year_transfer_summary.csv
transfer_matrices/
```

The raw file contains repeat-level OOF estimates. The summary and matrices use
means across repeats and include variability/interval columns.

### Feature-importance stability

Native importance is now estimated across stratified bootstrap samples rather
than one full-data fit. New output:

```text
feature_importance_stability.csv
```

It reports mean importance, importance SD, top-k selection frequency, consensus
rank, and completed bootstrap count. `feature_rank_overlap.csv` uses consensus
top-k sets.

## General Predictor

Files:

```text
src/PhenGO/predict/PhenGO-Predict.py
src/PhenGO/predict/compare.py
src/PhenGO/predict/data.py
src/PhenGO/predict/sklearn_models.py
src/PhenGO/predict/utils.py
```

### Safer model defaults and probability handling

- Logistic regression replaces the neural network as the default model.
- Sklearn probabilities use training-fold calibration by default.
- Gradient boosting remains class-balanced through sample weights.
- Neural-network random seeds are explicit and deterministic operations are requested when available.
- Neural-network early stopping and threshold selection use a held-out validation subset.
- The test set is no longer used for threshold selection.
- Single-model output adds average precision, balanced accuracy, MCC, and Brier score.

### Reusable model artifacts

Sklearn models are saved as `model.joblib`. Every saved sklearn and neural model
has `model_schema.json`, recording:

- exact ordered GO feature names;
- `viable = 0`, `lethal = 1`;
- lethal as the positive class;
- decision threshold;
- calibration policy where applicable.

Unlabelled prediction now requires this schema or an explicit feature list,
aligns columns by name, detects whether a class column is actually present,
and rejects missing/extra/non-binary features. It no longer drops the final GO
column unconditionally.

### ARFF validation

- Quoted ARFF attribute names are parsed correctly.
- Malformed rows are errors rather than silently skipped.
- Duplicate attributes are errors.
- Missing/non-finite values are errors rather than converted to zero.
- GO values must be binary.
- Conflicting duplicate gene rows are errors; exact duplicates are removed.

### Multi-dataset comparison

Comparison mode now trains every requested model instead of always invoking
TensorFlow. It writes `comparison_metrics.csv`, model-specific predictions and
importance files, a performance plot, top-k importance overlap, and an outer
join differential-importance table that retains features present in only one
dataset.

## Temporal GO Analysis And ARFF Comparison

Files:

```text
src/PhenGO/scripts/GO_temporal_analysis.py
src/PhenGO/scripts/GO_Compare.py
src/PhenGO/scripts/compare_arff_genes.py
```

Temporal GO enrichment now uses:

- Fisher's exact test;
- Haldane-Anscombe corrected finite enrichment/depletion ratios;
- Benjamini-Hochberg FDR correction;
- `--max-fdr` for ranked timelines;
- `*_enrichment_statistics.csv` for all tested terms.

Infinite one-class ratios are no longer silently discarded. Temporal ARFF
loading also rejects conflicting duplicate genes. Output replacement requires
`-overwrite` and uses protected-path checks.
Strict source manifests and ARFF hashes are required by default. The command
writes `temporal_analysis_manifest.json` with source snapshot IDs, software and
source-tree state, input hashes, configuration, and output hashes.
Temporal analysis uses the same manifested source, environment, and policy
contract as version-sensitivity analysis, so mixed-generation datasets are
rejected consistently.

The two-dataset `GO_Compare` workflow now uses the same corrected finite effect
ratios, Fisher tests, and Benjamini-Hochberg FDR values. Enrichment-ranked terms
respect `--max-fdr`, malformed or conflicting ARFF records are rejected, and a
protected output directory with explicit overwrite behavior is required.

`compare-arff` is now symmetric, reports `MISSING_IN_A` as well as
`MISSING_IN_B`, writes one row per gene with combined statuses, and rejects
conflicting duplicate records. A gene with both a label and GO mismatch no
longer appears twice.

## Dependencies

File: `setup.cfg`

- `scipy` is now explicit in `predict`, `reports`, and `all` extras because the
  statistical workflows use it directly.
- `scikit-learn>=1.7` is required in those extras because comparable MLP
  training uses class-balanced sample weights.
- A `test` extra provides `pytest`.
- `networkx` remains the core ontology graph dependency.

The requested ARM64 environment used for implementation and verification was:

```text
/Users/nicholas/miniconda3/envs/PhenGO
```

`networkx 3.6.1` and `pytest 9.1.1` were installed there. The environment uses
ARM64 Python 3.12 and its existing NumPy, pandas, SciPy, scikit-learn, and
TensorFlow packages.

The current checkout was installed into that environment as editable
`PhenGO 0.3.0`, and the installed `PhenGO`, `phengo-version-sensitivity`, and
`temporal-analysis` entry points were exercised successfully.

## Documentation

Files:

```text
README.md
docs/version_sensitivity_workflow.md
docs/manuscript_limitations.md
docs/publication_pipeline.md
data/PUBLICATION_DATA_README.md
```

The README now describes the current label, GO, provenance, prediction, and
analysis behavior. The research workflow gives a complete protocol from input
ledger through manuscript figures. The limitations file is written as text
that can be inserted into a manuscript and distinguishes software-controlled
choices from unavoidable biological, archival, and inferential limitations.
The publication-pipeline document explains the one-shot command, experimental
tracks, organism-specific label rules, comparable MLP, output tree, completion
gates, result sequence, and includes paper-ready Methods text.

## Publication Data And Automation

Files:

```text
data/publication_snapshots.tsv
data/File_Overview_publication.xlsx
scripts/run_publication_pipeline.sh
src/PhenGO/scripts/publication_inputs.py
src/PhenGO/scripts/publication_summary.py
```

- The 80-row TSV is the runtime authority: 49 required declared composite
  organism-year cross-sections and 31 optional fixed-label annotation anchors.
- A cleaned, human-facing workbook is generated from the ledger. Its
  formula-driven coverage sheets report 49 required rows, 31 annotation-anchor
  rows, and live file-availability checks.
- The master ARM64 script validates inputs, compression, tests, architecture,
  environment, and provenance; derives mouse, fly, and worm inputs; creates and
  validates all ARFF tracks; runs all six comparable models and four baselines;
  runs temporal GO analysis; aggregates paper tables; generates reports; and
  hashes the results.
- The repeated engine evaluates all model families in every year. The older
  held-out TensorFlow/sklearn suite runs all model types on 2025 by default as a
  supplementary implementation check; `PHENGO_SINGLE_SNAPSHOT_YEARS` can
  request additional years without making that weaker design the primary one.
- Required analyses cover the primary annual cross-sections, `is_a`-only,
  no-IEA, fixed-2025
  `is_a + part_of`, and fixed-2025 `is_a` tracks. The 25 newly supplied
  2010-2014 GAFs were validated and added to the complete annotation series.
- Historical GO handling now records the actual format eras: full `go.obo` for
  2010-2012, first available `go-simple` in 2013, and `go-basic` from 2014.
  The long 2010-2025 track uses only `is_a` because early `part_of` edges cross
  namespaces; the directly comparable `is_a + part_of` series starts in 2014.
- Annotation anchors are required by default. Exploratory partial runs must
  explicitly pass `--allow-missing-anchors`; missing inputs are never replaced.
- Publication aggregation now includes prediction instability, feature-rank
  stability, previous-available-snapshot labels with explicit calendar gaps,
  complete temporal enrichment statistics,
  and supplementary single-snapshot model results.
- Dry-run mode now emits downstream ML commands even though ARFFs are not
  physically created, so every model and MLP parameter can be audited before a
  long run.
- `scripts/generate_publication_mock.py` creates a deterministic pre-run design
  review with eight schema-matched synthetic tables and eight publication-style
  figure panels. Every artifact is visibly marked as synthetic and non-citable.

## Publication Pipeline V2 And Correspondence-Informed Controls

Files:

```text
scripts/run_publication_pipeline_v2.sh
scripts/v2/alternative_inputs.py
scripts/v2/publication_multiverse.py
data/resource_availability.tsv
docs/publication_pipeline_v2.md
tests/test_v2_alternative_inputs.py
tests/test_v2_publication_multiverse.py
tests/test_temporal_intervals.py
```

- Provider correspondence is represented as an availability/context ledger,
  including resources that were never generated, were not preserved after a
  platform migration, were temporarily inaccessible, were recovered through a
  web archive, or remain unconfirmed. Missing years are not imputed.
- Organism-year rows are now described as declared composite cross-sections of
  independently released phenotype, GAF, and ontology resources, rather than
  synchronized database releases.
- The Fly sensitivity tracks derive paired gene sets with a fail-closed parser.
  It recognizes top-level space, slash, and `with` syntax outside allele
  brackets; retains only completely resolved accessory or same-gene contexts;
  and excludes accessory primaries, second genes, unknown components, ambiguous
  phenotypes, and gene-level lethal/viable conflicts. Every source row has an
  audit outcome.
- Annual Worm lethal-term sets are regenerated from only the three complete-
  lethality roots confirmed by the supplied WormBase correspondence. `NOT`
  assertions remain exclusions rather than evidence for viability.
- Previous-label baselines and transfer outputs report signed and absolute
  calendar gaps. Temporal lifecycle fields refer to available snapshots, and
  GO enrichment trends are per calendar year. Multiple captures from one
  nominal year are averaged before fitting a trend so that year is not
  inadvertently upweighted.
- V2 runs every configured model family and panel for every eligible track with
  at least two snapshots. Alternative IMPC and recovered MGI negative classes
  are explicitly labelled operational rather than experimentally complete.
- The provider-GAF sensitivity track now includes the two available ZFIN files.
  The 2017 provider/archive stable-ID-plus-GO Jaccard is 0.8271, whereas the
  provider file dated 2017-12-15 is identical to the selected 2018 archive GAF;
  exact source filenames are retained in release metadata.
- The provider-GAF track also includes SGD for the five years with compatible
  archived yeast phenotype snapshots; all 16 annual provider GAFs remain in the
  input audit without fabricating phenotype coverage for missing years.
- V2 is now the single master entry point. It builds a missing V1, resumes an
  inactive partial V1, or audits a completed V1 before running the extension
  and consolidation. A completed run with old schemas, unspecified component
  alignment, orphan or hash-mismatched ARFFs, duplicate snapshot IDs, or mixed
  source revisions is rebuilt cleanly instead of being silently reused. A clean
  V1 from an older source revision is also rebuilt, and consolidation requires
  V1 and V2 to declare the same source hash.
  Date-derived defaults and `--run-root` support
  future one-command runs; process locks and legacy-process discovery prevent
  running or suspended concurrent writers. Recent-log fallback is used only
  where process-table discovery is unavailable, avoiding false active-run
  detection after an intentional stop. Resume now
  preserves command and skipped-input ledgers, while planning-only markers make
  a dry run harmless to a subsequent real run. A fresh-root master dry run
  expands both V1 and V2, while an incomplete V1 with scientific files is left
  untouched.
- Both Bash entry points now change to their repository root before invoking
  relative tests or helper paths, so absolute-path and `launchd` execution do
  not depend on the caller's working directory.
- Consolidation now checks free space, creates independent copies of both
  completed source runs at `source_runs/v1` and
  `source_runs/v2_extension`, hashes and verifies every copied file, and builds
  all combined tables and figures from those copies. APFS uses copy-on-write
  clones to avoid immediately doubling disk blocks; other filesystems use
  ordinary copies. The extension is frozen before copying so its log cannot
  change during verification. Its transient live-process lock is kept outside
  the copied run tree.
- `consolidated_manifest.json`, `source_run_file_inventory.tsv`, and
  `consolidated.complete` make the final directory independently auditable.
  `phase_source_hash_inventory.tsv` additionally exposes the exact source hash
  declared by each ARFF, ML, and temporal-analysis manifest. Mandatory base and
  extension integrity reports are also included. The original V1 and V2 roots
  remain untouched after successful generation.

## Tests

New tests:

```text
tests/test_core_integration.py
tests/test_go_relations.py
tests/test_go_and_phenotype_policies.py
tests/test_go_statistics.py
tests/test_label_sources.py
tests/test_sklearn_calibration.py
tests/test_version_cli_integration.py
tests/test_version_sensitivity_oof.py
tests/test_publication_inputs.py
tests/test_publication_summary.py
tests/test_temporal_intervals.py
tests/test_v2_alternative_inputs.py
tests/test_v2_publication_multiverse.py
```

They add coverage for:

- strict end-to-end ARFF generation and manifests;
- stable accession output and lethal/viable polarity;
- `is_a` and `part_of` propagation;
- alternate GO IDs and obsolete replacement chains;
- exact worm evidence, qualifier, and phenotype-term filtering;
- stable GAF joins;
- global rejection of one-to-many symbols, strict synonym disabling, stable-ID
  phenotype parsing, authoritative paired gene sets, and fly assignment conflicts;
- finite corrected GO effect sizes, Fisher tests, and FDR correction;
- nested probability calibration and calibrated-model feature importance;
- gene-disjoint transfer and prediction-instability scoring;
- complete command-line version-sensitivity output and analysis manifests;
- canonical ARFF binary serialization and invalid feature rejection;
- paired label-set completeness, disagreement filtering, and audit manifests;
- rejection of cross-namespace GO propagation without an explicit override;
- cross-snapshot policy-drift and incomplete-evaluation rejection.
- comparable MLP probability, calibration, class-balancing, convergence, and
  early-stopping preflight behavior;
- training-vocabulary-only transfer and fixed gene-fold reuse;
- publication input conversion, compressed OBO traversal, and result-table
  aggregation.

The original tests continue to cover label aliases, ancestor filtering, quoted
gene IDs, mouse report layouts, mouse accession joins, ARFF discovery, model
expansion, and natural release sorting.

### Verification

Using `/Users/nicholas/miniconda3/envs/PhenGO`:

```text
pytest                         64 passed
focused V2 tests              15 passed
compileall                     passed
Ruff on changed Python files   passed (system Ruff; not installed in the env)
git diff --check               passed
pip check                      no broken requirements
bash -n both pipeline scripts  passed
V1 all-organism dry run        287 ARFF builds and 25 all-model analyses emitted
V2 full dry run                136 ARFF builds and 15 all-model analyses emitted
fresh-root master dry run      423 ARFF builds and 40 all-model analyses emitted
real-data NN smoke             3 yeast releases, 9 transfer directions complete
real-data identifier smoke     2016 snapshots for all five organisms use stable-ID joins
real fish regression           fish 2016 ARFF rebuilt and validated; original collision absent
real FlyBase parser audit      393,904 rows; 7,586 lethal, 3,286 viable, 1,541 mixed genes excluded
real ZFIN provider smoke       2017 ARFF built and validated; 3,591 final genes, manifest schema 3
```

Core, version-sensitivity, and temporal-analysis command help, along with
`PhenGO --print-defaults`, were also exercised successfully.

## Required Regeneration And Migration

Regenerate ARFF files when any manuscript result depends on:

- stable row identifiers;
- conservative label policies;
- worm evidence filtering;
- fish ambiguous-phenotype handling;
- fly lethal or multi-gene filtering;
- GO relation selection;
- GO alternate/obsolete mappings;
- root or prevalence filtering;
- strict provenance manifests.

In practice, regenerate every publication ARFF. The version-sensitivity command
will otherwise reject the files unless legacy manifest checks are disabled.

After regeneration, rerun all model, temporal enrichment, ARFF comparison, and
report outputs. Old trained models are not schema-compatible with a changed GO
feature matrix and should not be reused.
