# PhenGO

## Overview

This project provides a unified Python-based tool to generate ready-to-use WEKA ARFF formatted files, specifically designed for machine learning applications involving gene essentiality prediction. 
The tool integrates phenotype data and Gene Ontology (GO) annotations for genes from selected model organisms, streamlining the data preparation process.

## Purpose

The main goal of this project is to simplify and standardise the creation of ARFF files that combine phenotype information with GO-mapped gene data. 
This enables researchers to efficiently apply machine learning techniques (using WEKA or similar platforms) to analyse gene essentiality and related biological questions across various model organisms.

## Features

- **Unified Workflow:** Handles data collection, integration, and formatting in a single pipeline.
- **Model Organism Support:** Designed for commonly studied organisms (e.g., *Saccharomyces cerevisiae*, *Mus musculus*).
- **GO Annotation Integration:** Maps genes to their respective GO terms for comprehensive feature representation and traces obo files to acquire parent terms.
- **Phenotype Data Inclusion:** Incorporates phenotype labels for supervised learning tasks.
- **WEKA ARFF Output:** Produces files in the ARFF format, ready for immediate use in WEKA.

