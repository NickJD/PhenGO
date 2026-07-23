# PhenGO Publication Pipeline V2

## Purpose

`scripts/run_publication_pipeline_v2.sh` is an additive extension of the first
publication run. It tests whether apparently modest data-source and
dataset-construction choices alter the genes, labels, GO features, model
performance, transferability, and apparent biological feature importance.

The script is the sole user-facing publication entry point. If V1 is absent, it
builds V1; if V1 is incomplete and inactive, it resumes it; and if V1 is
complete, it verifies both its checksum ledger and its scientific integrity.
An old completed V1 is rebuilt cleanly when it contains pre-schema-3 manifests,
unspecified component alignment, an ARFF/manifest mismatch, an orphan ARFF,
duplicate snapshot IDs, or more than one source-tree hash. `--extension-only`
instead fails on such a base. A clean V1 from a different source revision is
also rebuilt, so V1 and V2 cannot accidentally compare code versions. A live
lock, a matching running or suspended legacy process, or a recently changing
V1 log prevents concurrent execution.
V2 writes new experimental datasets to a separate extension root.
A final consolidation step copies both completed runs
into a third directory, verifies every copied file, and computes joined
paper-facing tables and figures from those copies. The consolidated directory
is therefore a self-contained publication artifact.

On APFS, consolidation uses independent copy-on-write clones and records that
strategy in the manifest. These are real files that remain valid if a source
run is later removed, but initially share unchanged disk blocks. On other
filesystems V2 uses ordinary copies and requires enough free space for them.

## Run Command

Before reusing the current run root, ensure its earlier pipeline process has
exited. The single command below then audits that partial run, detects that its
source changed during generation, and regenerates V1 before proceeding:

```bash
cd /Users/nicholas/Git/PhenGO

/Users/nicholas/Git/PhenGO/scripts/run_publication_pipeline_v2.sh \
  --run-root /Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs/publication_run_2026-07-22
```

For a completely new future run, use the same command with a new run-root name;
no separate V1 preparation is needed. Invoking the script without arguments
uses `publication_run_YYYY-MM-DD` for the current date. After interruption, use
the same command with `--resume`. Resume is accepted only if the script, source
tree, base fingerprint, track selection, and principal ML settings have the
same fingerprint. `--extension-only` restores the old behavior of requiring a
completed V1.

A two-year preflight is available with `--quick`. It is a smoke analysis, not a
publication analysis. On a fresh root, `--dry-run` expands both V1 and V2 so the
entire command plan can be audited. It records planning-only markers; a later
real run recognizes them and starts cleanly without requiring `--resume`. An
incomplete V1 that already contains scientific outputs is never modified by a
master dry run.

## Experimental Tracks

| Track | Organism(s) | Controlled change |
|---|---|---|
| `primary_mod_gaf` | fish, fly, worm, yeast | Primary annual phenotype policies and annual `go-basic.obo`, but MOD-provided GAF instead of the GO release archive GAF. Fish is limited to the two locally available ZFIN files (2017 and the 2017-12-15 file stored with the 2018 collection); yeast is limited to years with compatible archived phenotype files. |
| `legacy_like_go_archive` | fish, fly, yeast | Approximation of earlier permissive label construction using the GO archive GAF. Fly permits multigene rows; nonlethal observations are treated as viable; lethal labels win. Fish also assigns ambiguous semi-lethal/semi-viable records to the operational viable class. |
| `legacy_like_mod_gaf` | fly | Legacy-like label construction combined with the FlyBase-provided GAF. |
| `fail_closed_driver_aware_go_archive` | fly | Explicit lethal/viable labels with component-by-component genotype parsing and the GO archive GAF. Recognized accessory systems and exact same-gene combinations may be retained; true multi-gene and every unresolved compound record are excluded. |
| `fail_closed_driver_aware_mod_gaf` | fly | The same fail-closed FlyBase parsing combined with the FlyBase-provided GAF. |
| `worm_fixed_2025_terms_go_archive` | worm | Annual phenotype records with the 2025 WormBase lethal-term closure fixed across years and GO archive GAF. |
| `worm_fixed_2025_terms_mod_gaf` | worm | Fixed 2025 lethal-term closure combined with the taxon-filtered WormBase GAF. |
| `mouse_impc_assertions` | mouse | Raw annual IMPC genotype-phenotype assertions classified with the fixed package MP lethal-term list. |
| `mouse_mgi_genepheno_fixed_terms` | mouse | Complete, nonduplicate historical MGI GenePheno captures classified with fixed MP lethal terms. Distinct same-year captures retain separate IDs. |
| `mouse_mgi_phenogenomp_fixed_terms` | mouse | Complete, nonduplicate historical MGI PhenoGenoMP captures classified with fixed MP lethal terms. |

All GO propagation uses `is_a` plus `part_of`, excludes ontology roots, and uses
the same minimum gene count and maximum gene fraction as V1. Every ARFF is
strict-manifested and validated.

With the current data collection, the default V2 plan contains 136 extension
ARFF builds and 15 eligible track-by-organism ML/temporal series. Each of those
15 ML commands requests every model family and all four cohort/feature panels.
Together with V1's planned 287 ARFFs, the consolidated inventory should contain
423 snapshots, subject to every selected input passing the fresh runtime audit.

## Alternative Input Safety

Before ARFF generation, `scripts/v2/alternative_inputs.py` inventories every
alternative GAF, IMPC assertion file, Wayback MGI report, and available MP
vocabulary report. It records physical and decompressed SHA-256 values, file
sizes, row counts, newline completion, column-width distributions, compression,
collection year, years embedded in the source filename, duplicate groups,
eligibility, and exclusion reason.

The current audit detects 92 files: 71 eligible, 10 exact decompressed
duplicates, 5 captures with truncation signatures, 3 vocabulary-only files, and
3 metadata files. The script refuses to use a selected mouse input unless the
fresh audit marks that exact path eligible. A missing selected IMPC or MGI input
also stops the default run; it is never converted into a silent partial series.

The selected MGI captures are:

- GenePheno: 2016 gzip; both distinct 2017 captures; 2020 gzip; both distinct
  2021 captures; 2022, 2024, and 2025 gzip captures.
- PhenoGenoMP: 2017, 2020, 2021, 2024, and 2025 gzip captures.

Exact duplicate nominal releases and truncated files remain visible in the
audit but are excluded from ML. No 2026 mouse ARFF is made because there is no
eligible complete phenotype/GAF/ontology cross-section.

ZFIN, FlyBase, WormBase, and SGD GAFs are filtered to exact taxa `taxon:7955`,
`taxon:7227`, `taxon:6239`, and `taxon:559292`, respectively. This is essential
for the multi-species WormBase files. GAF source overlap is calculated
principally from stable gene ID plus GO ID. Relation qualifiers are not part of
that key because GAF 2.1 and 2.2 encode relations differently;
evidence-code-aware overlap is also reported. The 2017 ZFIN provider/archive stable-ID-plus-GO Jaccard is
0.8271; the ZFIN file dated 2017-12-15 is identical to the selected 2018 archive
GAF on this key. The source filename is retained in release metadata.
The SGD provider and archive files are nearly identical in 2015-2017, whereas
the provider files contain substantially more rows in 2024-2025; their
stable-ID-plus-GO Jaccards are 0.9724 and 0.9751, respectively. Retaining both
periods tests whether apparently small set-level differences mask larger
annotation-count and evidence-composition changes.

FlyBase compound records are parsed separately from the PhenGO core and written
as paired gene sets. The parser never assumes that an unrecognized component is
harmless. Top-level spaces, slashes, and `with` clauses are parsed only outside
allele brackets, so commas or other syntax inside an allele name do not create
false components. It records the original row, schema width, primary gene,
component roles, provisional label, exclusion reason, and final gene-level outcome in
`01_derived_labels/fly_fail_closed/<year>/row_audit.tsv.gz`. A gene with both
lethal and viable candidate observations is excluded from this track.

## Label Meaning

The IMPC assertion and historical MGI association files do not provide a
complete, experimentally established viable class. Their negative class is
therefore explicitly operational: a gene has at least one other phenotype
assertion or association but no selected lethal MP term. A gene with both a
selected lethal term and other phenotypes remains lethal under the declared
`lethal_wins` derivation rule. Genes with invalid or nonunique stable IDs are
excluded rather than symbol-matched.

These tracks test source-definition sensitivity. They must not be described as
definitive viable-versus-lethal training sets. Higher performance in any track
is not evidence that its construction is biologically more correct; it can be
caused by cohort contraction, class balance, easier ascertainment structure, or
curation leakage.

## ML and Temporal Analysis

Each track with at least two snapshots runs all model families from
`PhenGO.predict.version_sensitivity`, including majority, prior, and stratified
random baselines plus logistic regression, random forest, gradient boosting,
SVM, and the comparable scikit-learn MLP. The MLP uses the same fixed gene
folds, panels, repeats, class policy, nested calibration, thresholds, metrics,
and transfer directions as the other learned models.

The default analysis includes:

1. Repeated, gene-disjoint within-snapshot cross-validation.
2. All train-snapshot to test-snapshot transfer directions.
3. `full`, `matched_features`, `matched_genes`, and `matched_both` panels.
4. Previous-available-snapshot label baselines with the calendar interval
   reported explicitly.
5. Fixed-gene out-of-fold prediction and probability instability.
6. Bootstrap feature importance and top-GO rank overlap.
7. GO enrichment, class-label, gene, GO-feature, sparsity, ontology-depth, and
   pairwise temporal analyses.
8. HTML reports and extension-run publication summary tables.

## Output Layout

The V2 extension root contains:

- `00_run_metadata/`: configuration, base fingerprint, environment, commands,
  checksums, skipped-input ledger, and completion marker.
- `01_derived_inputs/`: deterministic taxon-filtered GAFs.
- `01_derived_labels/`: paired label sets, exclusions, assignment tables, and
  derivation summaries.
- `02_arff/`: strict per-track, organism, and timepoint PhenGO datasets.
- `03_qc/`: ARFF validation, full input audit, and archive-versus-MOD GAF
  comparisons.
- `04_ml/`: repeated CV, transfer, instability, and feature-stability outputs.
- `05_temporal/`: GO and phenotype temporal analyses.
- `06_publication_tables/` and `07_reports/`: extension-level summaries.

The consolidated root contains:

- `source_runs/v1/`: verified physical copy of the complete first run.
- `source_runs/v2_extension/`: verified physical copy of the complete extension.
- `source_run_file_inventory.tsv`: path, size, and SHA-256 for every copied file.
- `phase_source_hash_inventory.tsv`: source-tree hash declared by every ARFF,
  version-sensitivity, and temporal-analysis manifest.
- `base_run_integrity.json` and `extension_run_integrity.json`: the mandatory
  fail-closed scientific audits applied before consolidation.
- `resource_availability.tsv`: provider-correspondence and retrieval-history
  context for non-random archive gaps.
- `snapshot_inventory.tsv`
- `primary_alternative_comparisons.tsv`
- `dataset_drift_all_tracks.csv`
- `pairwise_dataset_drift_all_tracks.csv`
- `within_year_cv_all_tracks.csv`
- `performance_deltas_vs_primary.csv`
- `cross_year_transfer_all_tracks.csv`
- `previous_available_snapshot_label_baseline_all_tracks.csv`
- `prediction_instability_all_tracks.csv`
- `feature_rank_overlap_all_tracks.csv`
- `feature_importance_stability_all_tracks.csv`
- `evaluation_preflight_all_tracks.csv`
- `temporal_summary_all_tracks.tsv`
- `temporal_enrichment_statistics_all_tracks.tsv`
- `gaf_source_comparisons.tsv`
- `alternative_input_audit.tsv`
- `mouse_vocabulary_inventory.tsv`
- `mouse_vocabulary_pairwise.tsv`
- `figures/figure_cohort_sizes.png`
- `figures/figure_policy_sensitivity.png`
- `figures/figure_predictability_vs_coverage.png`, when matched performance rows exist
- `CONSOLIDATED_REPORT.md`
- `consolidated_manifest.json` and `consolidated.complete`

The consolidated manifest records the original and copied source roots, full
tree digests and file counts, completion/fingerprint/checksum ledgers, table
counts, the complete command, and a SHA-256 for every newly generated
consolidated output. The workflow checks free disk space before copying and
fails if any copied path, size, or digest differs from its source. The final
completion marker records the consolidated manifest's SHA-256.

This is one Bash entry point, not one monolithic implementation file. It calls
the versioned PhenGO Python modules and V2 helper programs in the repository so
that each parser and analysis remains independently testable. The consolidated
directory contains all computed V1 and V2 results needed for the article; the
immutable raw source archives remain under `data/` and are referenced by path
and SHA-256 rather than duplicated into every publication result.

Nominal years in the consolidated inventory are analytical labels for declared
composite cross-sections. They do not imply that all contributing resources
were released together. The inventory retains each component release and
availability field, while temporal tables report the interval between adjacent
available observations.

## Recommended Paper Use

Start with input, cohort, label, and GO-feature drift before presenting any ML
metric. Use V1's primary matched panels and transfer matrices to establish
year-to-year instability. Then use the V2 tracks as controlled perturbations:
show cohort overlap, label agreement, feature overlap, and prevalence beside
performance deltas. Cross-year transfer, per-gene prediction instability, and
feature-rank instability should carry more interpretive weight than isolated
within-snapshot accuracy.

The central claim supported by this design is not that one historical source is
wrong and another is right. It is that undocumented or weakly versioned choices
about source release, identifier resolution, phenotype inclusion, ontology-term
closure, GAF provider, and filtering policy can materially change the apparent
performance and biological explanation of an ML model, even for intensively
curated model organisms.
