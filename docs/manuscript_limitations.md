# Manuscript-Ready Limitations

## Limitations

This study evaluates the sensitivity of machine-learning results to changing
model-organism database snapshots. It does not estimate a purely biological
effect of calendar time, nor can it determine whether every change between
snapshots represents error, correction, improved biological knowledge, or a
change in curation policy. Database release is therefore treated as a composite
exposure that includes changes in gene coverage, phenotype assertions, Gene
Ontology (GO) annotations, ontology structure, identifier mappings, evidence
availability, and curation practice. Associations between release and model
output should not be interpreted causally.

Historical resources are not equally complete or equally recoverable. Some
archives may lack files, metadata, evidence fields, or release-matched helper
resources that are available for recent releases. Even when a file is dated,
its phenotype data, GO annotation file, GO ontology, and phenotype ontology may
have been released on different dates. PhenGO records release labels, retrieval
dates, paths, hashes, software versions, and analysis parameters, but provenance
metadata cannot reconstruct information that was never archived. Consequently,
some observed temporal differences may reflect archival survivorship or
imperfect synchronization among source resources.

Direct correspondence with database providers showed that these gaps have
different causes and are therefore not missing at random. Some static files
were never generated because an online query service could produce current
data on demand; later platform migration did not preserve the ability to
reconstruct those historical states. Other archived downloads became publicly
inaccessible because of hosting and automated-traffic costs. One historical
mouse-data request was closed without a substantive archival answer in the
supplied correspondence. PhenGO therefore distinguishes files that were not
generated, not preserved, temporarily inaccessible, recovered through a web
archive, or lacking provider confirmation. Adjacent available snapshots are
reported with their calendar interval and are not automatically called
consecutive years.

The binary target is an operational harmonization of organism-specific
phenotype assertions rather than a universal measurement of essentiality.
Lethality can depend on developmental stage, sex, genetic background,
environment, perturbation type, allelic strength, zygosity, and whether a record
describes an allele, genotype, or gene. These dimensions are represented
differently across FlyBase, SGD, WormBase, ZFIN, and MGI and cannot be fully
standardized from the available tabular exports. The pipeline excludes
ambiguous and mixed observations by default and requires explicit viability
evidence, but this conservative policy reduces sample size and does not make
the organism-specific targets biologically identical. Multi-gene and complex
genotype records may also be incompletely identified by source-specific
filters. FlyBase compound records are especially heterogeneous: an accessory
expression system, same-gene transheterozygote, and simultaneous perturbation
of two genes may use similar textual syntax. The primary analysis excludes all
compound contexts. The fail-closed sensitivity analysis retains a compound row
only when every component has a resolved role and excludes all unknown
contexts; this reduces false assignment at the cost of cohort size. Results
should therefore be interpreted within each organism, and
cross-organism differences should not be presented as direct performance
comparisons on an identical endpoint.

Viable genes are not an unbiased set of true negatives. They are genes for
which explicit viability evidence was available under the selected extraction
policy. Genes that were not tested, were tested only under unrepresented
conditions, or had only non-viability phenotypes are excluded rather than
assumed viable. This avoids one form of label contamination but introduces
ascertainment bias: well-studied genes and phenotypes are more likely to enter
the dataset. Publication and curation biases may therefore affect both class
labels and GO feature density.

For WormBase and MGI, explicit negative labels require release-matched viable
phenotype-term or viable-gene lists supplied by the investigator. The criteria
used to construct those lists are themselves an operational definition and
must be documented. If such lists cannot be reconstructed for a historical
release, the less conservative assumption that any observed nonlethal
phenotype implies viability may be used only as a declared sensitivity
analysis; it cannot be considered equivalent to direct viability evidence.

Stable database accessions reduce symbol drift but do not eliminate identifier
uncertainty. Gene records can be merged, split, withdrawn, or remapped between
releases. Restricting analyses to identifiers shared by all snapshots creates a
survivor cohort enriched for stable, well-curated genes and may understate
version effects affecting newly introduced or retired records. Conversely, the
full panel includes gene-set churn and therefore estimates the practical effect
of using each release, not a controlled annotation-only effect. We report both
full and matched panels because neither answers every question alone.

GO features are binary indicators and do not preserve annotation multiplicity,
experimental context, annotation date, or all qualifier semantics. `NOT`
annotations are excluded and evidence-code filters can be fixed by the user,
but evidence composition and annotation depth may still differ by release.
Alternate GO identifiers are canonicalized, one-to-one `replaced_by` mappings
are applied, and obsolete terms without an unambiguous replacement are removed.
These operations avoid duplicate or invalid features but can remove historical
assertions that cannot be translated uniquely into a later ontology.

The primary analysis propagates positive annotations through `is_a` and
`part_of` relations in a release-matched `go-basic` ontology. These relations
support biologically valid grouping to broader classes and necessary wholes.
PhenGO also supports pre-specified combinations including `regulates`,
`positively_regulates`, `negatively_regulates`, `occurs_in`, `capable_of`, and
`capable_of_part_of`. Those relations change the meaning of the gene-term
association and should not be interpreted as equivalent forms of ancestry.
The relation set, namespace set, propagation mode, evidence policy, and
prevalence filters must therefore be fixed before comparing years. Analyses
using alternative relation sets should be reported as sensitivity analyses
rather than selected after inspecting model performance.

The pre-2014 GO archives require an additional historical qualification.
`go-simple` first appears in the 2013 archive, whereas the 2010-2012 anchors use
the full archived ontology. In PhenGO's selected graph, `is_a` is acyclic and
namespace-safe in all four early releases, but `part_of` contains cross-namespace
edges. The continuous 2010-2025 fixed-label analysis is therefore restricted to
`is_a`. A separate 2014-2025 fixed-label analysis uses `go-basic` and
`is_a + part_of`, matching the primary relation policy. These tracks answer
related but non-identical questions and must not be pooled as one homogeneous
time series.

Fixed paired lethal/viable collections allow GO-resource sensitivity to be
examined while holding labels constant, but they introduce a different form of
selection. A collection assembled from a recent release incorporates knowledge
that was unavailable historically and must be described as a fixed contemporary
reference, not as the label set an investigator would have obtained in an
earlier year. Agreement analyses further restrict the cohort to genes supported
by both release records and the fixed collection; their improved apparent label
quality is therefore accompanied by reduced coverage and enrichment for stable,
well-curated genes. Results from release-derived, fixed-label, and agreement
cohorts answer different questions and should not be pooled.

GO terms are hierarchically dependent and often highly correlated. Native tree
importance and absolute logistic-regression coefficients can distribute signal
arbitrarily among correlated ancestors and descendants. Bootstrap selection
frequencies and between-snapshot top-feature overlap quantify instability but do
not turn feature importance into a causal measure of biological mechanism.
Individual GO terms should not be described as causal determinants of
lethality. Temporal GO enrichment is tested with Fisher's exact tests and
Benjamini-Hochberg correction, but dependence among GO terms and repeated
inspection across releases remains; enrichment results are descriptive and
should be interpreted with effect sizes and adjusted p-values together.

The machine-learning models are standardized probes of dataset sensitivity, not
optimized clinical or production predictors. Hyperparameters are intentionally
held largely constant to avoid allowing year-specific tuning to become another
source of variation. This improves comparability but may leave some models below
their best attainable performance. Probability calibration is fitted only
within each outer training fold, and all reported transfer and gene-level
instability predictions are gene-disjoint and out-of-fold. Calibration can
still be unreliable in small or highly imbalanced snapshots, especially when
the number of lethal genes is low.

The primary neural-network analysis uses scikit-learn's multilayer perceptron so
that it can share exactly the same repeated folds, transfer panels,
class-balancing, calibration, threshold, and output metrics as the other model
families. This is not architecture equivalence: optimization behavior and
capacity remain model-specific, and early stopping introduces a training-only
validation subset. Architecture and optimizer limits are fixed across years,
and convergence warnings and iteration counts are reported. The separate
TensorFlow single-snapshot network is retained only as an implementation
sensitivity analysis and should not be substituted into the primary comparative
table.

Repeated cross-validation estimates share the same underlying genes and are not
independent biological replicates. The reported 95% intervals are resampling
intervals across repeated partitions, not confidence intervals over the full
population of genes or future database releases. Random gene-level splits may
also place paralogous or functionally related genes in both training and test
folds. Thus, absolute predictive performance may be optimistic for transfer to
new gene families, even though each scored gene itself is excluded from its
training fold.

Finally, many models, metrics, organisms, panels, releases, GO terms, and
pairwise transfer directions are examined. The analysis pre-specifies primary
models, metrics, and panels, and GO enrichment uses false-discovery-rate
correction, but the complete output remains a high-dimensional descriptive
analysis. The manuscript should identify a limited primary endpoint set before
examining results, report all planned comparisons, and treat secondary patterns
as exploratory. Independent replication using separately archived releases or
another harmonized phenotype source would provide the strongest test of whether
the observed version sensitivity generalizes.
