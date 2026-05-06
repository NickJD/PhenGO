"""
PhenGO — A toolkit for building ML-ready ARFF files from model-organism
phenotype databases and Gene Ontology (GO) annotations.

Sub-packages
------------
core      : ARFF-generation pipeline (PhenGO main platform)
predict   : Neural-network gene-essentiality prediction (PhenGO-Predict)
scripts   : Auxiliary analysis and utility tools
"""
from .constants import configure_logger

# Package-level logger (stream only; individual entry points add file handlers).
configure_logger('PhenGO', enable_file=False)

__all__ = []

