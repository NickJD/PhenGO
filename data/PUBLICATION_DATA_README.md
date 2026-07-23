# PhenGO publication data collection

`publication_snapshots.tsv` is the machine-readable authority for the publication
workflow. `File_Overview.ods` remains the historical collection notebook; it is not
used to select files at runtime because it contains stale paths, copy-forward entries,
and no explicit label-policy fields.

`File_Overview_publication.xlsx` is the corrected human-facing workbook generated
from the ledger. It includes formula-driven coverage totals, input-availability
checks, the complete snapshot table, and the publication protocol. The TSV remains
authoritative if the workbook and ledger ever differ.

## Current coverage

The required cohort contains 49 declared organism-year cross-sections:

| Organism | Years | Phenotype target source |
|---|---|---|
| fish | 2015-2025 | Annual ZFIN phenotype records |
| fly | 2015-2025 | Annual FlyBase allele/genotype phenotype records |
| mouse | 2015-2025 | Annual IMPC viability calls |
| worm | 2015-2025 | Annual WormBase phenotype associations and phenotype ontology |
| yeast | 2015-2017, 2024-2025 | Archived SGD null-mutant phenotype records |

All required rows use the annual GO archive GAF and the nearest available annual
`go-basic.obo`. Phenotype, GAF, and ontology components were independently
released and are not described as one synchronized database release. Their
exact paths, release identifiers, retrieval-audit date, and label policies are
in the ledger; `resource_availability.tsv` records provider-confirmed archive
gaps and retrieval history.

The fixed-2025-label annotation cohorts reuse 2025 labels while changing the GAF
and GO ontology year. These cohorts can include years for which no phenotype archive
exists because they do not claim to reconstruct historical phenotype labels.

The V2 source-sensitivity analysis additionally audits and uses provider-side GAF
collections when a compatible cross-section exists. This includes ZFIN files for
the 2017 collection and a file dated 2017-12-15 stored with the 2018 collection,
as well as annual FlyBase and WormBase files. Sparse provider coverage is retained
as sparse coverage; missing provider years are recorded and are not filled from the
GO archive under the provider track. Annual SGD provider GAFs are also audited;
only 2015-2017 and 2024-2025 enter this track because those are the years with
compatible archived yeast phenotype files.

## Complete 2010-2025 annotation series

The GO archive GAFs for fish, fly, mouse, worm, and yeast are now present for every
year from 2010 through 2025. The 25 organism-year files added for 2010-2014 were
checked as gzip streams and parsed as GAF records: all records had at least the 15
required columns and every file contained the expected organism taxon. Their exact
paths are recorded in `publication_snapshots.tsv` and their checksums are captured
at run time.

The historical ontology formats require two pre-declared fixed-label analyses:

- `fixed_2025_is_a` spans 2010-2025. It uses full archived `go.obo` in
  2010-2012, the first available `go-simple` release in 2013, and `go-basic`
  from 2014 onward. Propagation is restricted to `is_a` in every year.
- `fixed_2025` spans 2014-2025 and uses `go-basic` with `is_a + part_of` in
  every year, matching the primary ontology policy.

This separation is necessary. PhenGO graph validation found no cross-namespace
`is_a` edges in the 2010-2013 ontologies, but their selected `part_of` edges include
cross-namespace links (422, 650, 802, and 983 respectively). Allowing those edges
would alter the biological meaning of ancestor propagation; silently dropping
`part_of` only in early years would instead introduce a year-dependent feature
policy. The two-track design avoids both problems. `go-simple` first appears in
the 2013 archive, but it does not provide the later `go-basic` cross-namespace
guarantee for the relations used here.

All annotation anchors are required by the default publication run. Use
`--allow-missing-anchors` only for exploratory partial runs; their absence is then
recorded and the affected rows are skipped rather than replaced with current data.

## Label construction

- **Yeast:** null-mutant records explicitly containing `inviable`/`essential` or
  `viable`/`non-essential`; mixed genes are excluded.
- **Fly:** explicit lethal and viable phenotype text; partial lethality, mixed genes,
  and conservative multi-gene records are excluded. A release-specific D. melanogaster
  assignment file is derived from `taxon:7227` GAF records.
- **Mouse:** annual IMPC calls are converted to paired lethal and viable gene sets.
  Subviable and conflicting genes are written to a separate exclusion file.
- **Fish:** lethal/dead evidence takes precedence; genes with only other observed
  phenotype evidence form an operational nonlethal class. This is not equivalent to
  an explicit viability assay.
- **Worm:** lethal descendants are regenerated from each annual WormBase phenotype
  ontology. Genes observed only outside that lethal set form an operational nonlethal
  class. This is not equivalent to an explicit viability assay.

Fish and worm are deliberately labelled `observed_viable` in the ledger so this
assumption is visible in every manifest and enforced consistently within each
organism's year series.

## Integrity rules

The master workflow fails before computation if a required path is missing, the
Python environment is not ARM64, a compressed input fails `gzip -t`, or the toolkit
tests fail. It records SHA-256 checksums for all inputs and outputs. Optional anchor
files are never silently substituted with current data.
