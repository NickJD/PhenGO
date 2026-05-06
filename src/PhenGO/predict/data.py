"""
predict/data.py — ARFF loading and data preparation for PhenGO-Predict.
"""
import logging
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder

logger = logging.getLogger(__name__)


def load_arff_data(filepath):
    """Load and parse an ARFF file manually to handle string attributes.

    Returns (DataFrame, metadata_dict) or (None, None) on failure.
    """
    logger.info(f"Loading ARFF file: {filepath}")

    data_lines = []
    attribute_names = []
    attribute_types = []
    in_data_section = False

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('%'):
                    continue
                if line.lower().startswith('@relation'):
                    continue
                if line.lower().startswith('@attribute'):
                    parts = line.split()
                    if len(parts) >= 3:
                        attr_name = parts[1].strip('"\'')
                        attr_type = ' '.join(parts[2:]).strip()
                        attribute_names.append(attr_name)
                        attribute_types.append(attr_type)
                    continue
                if line.lower().startswith('@data'):
                    in_data_section = True
                    continue
                if in_data_section:
                    values = []
                    current_value = ""
                    in_quotes = False
                    for char in line:
                        if char == '"' and not in_quotes:
                            in_quotes = True
                        elif char == '"' and in_quotes:
                            in_quotes = False
                        elif char == ',' and not in_quotes:
                            values.append(current_value.strip().strip('"\''))
                            current_value = ""
                        else:
                            current_value += char
                    values.append(current_value.strip().strip('"\''))
                    if len(values) == len(attribute_names):
                        data_lines.append(values)
                    else:
                        logger.warning(
                            f"Line has {len(values)} values but expected {len(attribute_names)}"
                        )
    except FileNotFoundError:
        logger.error(f"File {filepath} not found")
        return None, None
    except Exception as e:
        logger.error(f"Error reading file: {e}")
        return None, None

    if not data_lines:
        logger.error("No data found in ARFF file")
        return None, None

    df = pd.DataFrame(data_lines, columns=attribute_names)

    for col_name, col_type in zip(attribute_names, attribute_types):
        if col_type.lower() in ['numeric', 'real', 'integer'] or col_type == '{0,1}':
            df[col_name] = pd.to_numeric(df[col_name], errors='coerce')

    logger.info(f"Data shape: {df.shape}")
    logger.debug(f"Columns: {list(df.columns[:5])}... (showing first 5)")
    logger.debug(f"Attribute types: {set(attribute_types)}")

    return df, {'attribute_names': attribute_names, 'attribute_types': attribute_types}


def prepare_data(df):
    """Prepare a loaded ARFF DataFrame for neural-network training.

    Assumes: first column = gene name, last column = phenotype label,
    middle columns = binary GO-term features.

    Returns (X, y, gene_names, label_encoder).
    """
    logger.info("Preparing data...")

    gene_names = df.iloc[:, 0]
    phenotype  = df.iloc[:, -1]
    go_terms   = df.iloc[:, 1:-1]

    logger.debug(f"Gene names column: {df.columns[0]}")
    logger.debug(f"Phenotype column:  {df.columns[-1]}")
    logger.debug(f"Sample gene names: {gene_names.head().tolist()}")
    logger.debug(f"Sample phenotypes: {phenotype.head().tolist()}")
    logger.debug(f"Unique phenotypes: {phenotype.unique()}")

    X = go_terms.astype(float).values
    X = np.nan_to_num(X, nan=0.0)

    le = LabelEncoder()
    y  = le.fit_transform(phenotype)

    logger.info(f"Features shape:    {X.shape}")
    logger.info(f"Number of GO terms:{X.shape[1]}")
    logger.info(f"Class distribution:{np.bincount(y)}")
    logger.info(f"Class labels:      {le.classes_}")

    sparsity = (X == 0).sum() / X.size
    logger.info(f"Data sparsity:     {sparsity:.2%}")

    unique_vals = np.unique(X)
    logger.debug(f"Unique feature values: {unique_vals}")

    return X, y, gene_names, le

