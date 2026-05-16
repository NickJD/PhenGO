[![DOI](https://zenodo.org/badge/1026270084.svg)](https://doi.org/10.5281/zenodo.16748708)

# PhenGO

## Overview

This project provides a unified Python-based tool to generate ready-to-use WEKA ARFF formatted files, specifically designed for machine learning applications involving gene essentiality prediction. 
The tool integrates phenotype data and Gene Ontology (GO) annotations for genes from selected model organisms, streamlining the data preparation process.

## Purpose

The main goal of this project is to simplify and standardise the creation of ARFF files that combine phenotype information with GO-mapped gene data. 
This enables researchers to efficiently apply machine learning techniques (using WEKA or similar platforms) to analyse gene essentiality and related biological questions across various model organisms.

## Features

- **Unified Workflow:** Handles data collection, integration, and formatting in a single pipeline.
- **Model Organism Support:** Designed for commonly studied organisms (*Saccharomyces cerevisiae*, *Mus musculus*, *Drosophila melanogaster*, *Caenorhabditis elegans*, *Danio rerio*).
- **GO Annotation Integration:** Maps genes to their respective GO terms for comprehensive feature representation and traces OBO files to acquire parent terms.
- **Phenotype Data Inclusion:** Incorporates phenotype labels for supervised learning tasks.
- **WEKA ARFF Output:** Produces files in the ARFF format, ready for immediate use in WEKA.

## Installation
To install the PhenGO package, you can use pip:

```bash
pip install phengo
```

## Usage
### PhenGO Package:
#### PhenGO Example:
```commandline
PhenGO -species fly \
       -phenotype_file data/fly/phenotype_data/2017/allele_phenotypic_data_fb_2017_05.tsv.gz \
       -gene_association_file data/fly/gene_association/2017/gene_association_2017_05.fb.gz \
       -go_obo_file data/go/2017/go_2017-05-01.obo.gz \
       -output_dir Documents/PhenGO/fly_2017
```
The output will be saved in the specified output directory, which will contain the ARFF file and other relevant data files.
#### Menu:
```bash
usage: PhenGO [-h] [--print-defaults]
              -species SPECIES -phenotype_file PHENOTYPE_FILE
              -gene_association_file GENE_ASSOCIATION_FILE
              -go_obo_file GO_OBO_FILE -output_dir OUTPUT_DIR
              [-no_filter_unused_gos] [-filter_mixed_terms] [-gene_go_pheno]
              [-fly_lethal_genes FLY_LETHAL_GENES]
              [-fly_assignments FLY_ASSIGNMENTS]
              [-filter_multigenes]
              [-worm_phenotypes WORM_PHENOTYPES]
              [-worm_lethal_genes WORM_LETHAL_GENES]
              [-mouse_phenotypes MOUSE_PHENOTYPES]
              [-v]

PhenGO v0.1.2 - Convert phenotype and GO data to ARFF format

Required Options:
  -species SPECIES      Species tag (fly, yeast, fish, worm, mouse)
  -phenotype_file PHENOTYPE_FILE
                        Path to the phenotype data file (.gz)
  -gene_association_file GENE_ASSOCIATION_FILE
                        Path to the gene association file (.gz)
  -go_obo_file GO_OBO_FILE
                        Path to the go.obo file
  -output_dir OUTPUT_DIR
                        Output directory

Optional parameters:
  -no_filter_unused_gos Disable filtering of unused GO terms from the FUNC
                        and ARFF output (filtering is ON by default)
  -filter_mixed_terms   Filter out genes which have both lethal and viable
                        phenotypes (default: False)
  -gene_go_pheno        Deprecated — FUNC files are always written regardless.
                        Kept for backward compatibility.

Fly specific parameters:
  -fly_lethal_genes FLY_LETHAL_GENES
                        TSV file of specified lethal fly genes; use "default"
                        for bundled FlyBase_Lethal_Gene_IDs file
  -fly_assignments FLY_ASSIGNMENTS
                        TSV file confirming D. melanogaster gene assignments
                        (default: bundled FlyBase_Fields_2025_07_29.tsv.gz)
  -filter_multigenes    Filter out phenotypes involving multiple genes
                        (i.e. with 'with' clause or '/' in gene field)
                        (default: DO NOT FILTER)

Worm specific parameters:
  -worm_phenotypes WORM_PHENOTYPES
                        TSV file of lethal worm phenotype terms
                        (default: bundled lethal_terms_traversed_2025-08-12.tsv.gz)
  -worm_lethal_genes WORM_LETHAL_GENES
                        TSV file of pre-defined lethal worm genes; use "default"
                        for bundled WBPhenotype lethal genes file

Mouse specific parameters:
  -mouse_phenotypes MOUSE_PHENOTYPES
                        TSV file of lethal mouse phenotype terms
                        (default: bundled mouse_lethal_terms.txt.gz)

Misc:
  -v, --version         show program's version number and exit
  --print-defaults      Print default files and methods used for each species
```

### Compare-ARFF:
```commandline
usage: compare-arff [-h] -arff_a ARFF_A -arff_b ARFF_B -o OUTPUT

PhenGO v0.1.2 - Compare-ARFF: Compare two ARFF files.

options:
  -h, --help      show this help message and exit
  -arff_a ARFF_A  Master ARFF file (reference)
  -arff_b ARFF_B  Comparison ARFF file
  -o OUTPUT       Output CSV file

```

**Output:**
The output of the `compare-arff` function is a CSV file that summarizes the comparison between two ARFF files.
```commandline
Gene,Label A,Label B,GO Terms Differ,Status
GeneA,lethal,,,"MISSING_IN_B"
GeneB,lethal,viable,,"LABEL_MISMATCH"
GeneC,viable,viable,GO:0008150;GO:0003674,"GO_TERM_MISMATCH"
GeneD,viable,viable,,"EXACT_MATCH"
```
