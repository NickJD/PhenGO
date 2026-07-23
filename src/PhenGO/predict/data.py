"""
predict/data.py — ARFF loading and data preparation for PhenGO-Predict.
"""
import logging
import csv
import re

logger = logging.getLogger(__name__)

PHENOTYPE_CLASSES = ["viable", "lethal"]

ATTRIBUTE_RE = re.compile(
    r"^@attribute\s+(?:'([^']*)'|\"([^\"]*)\"|(\S+))\s+(.+)$",
    re.IGNORECASE,
)


class PhenotypeLabelEncoder:
    """Encode phenotype labels with lethal as the positive class.

    sklearn's LabelEncoder sorts labels alphabetically, which maps
    ``lethal -> 0`` and ``viable -> 1``.  The prediction code treats class 1 as
    the positive/essential class, so we keep the mapping explicit:
    ``viable -> 0`` and ``lethal -> 1``.
    """

    classes_ = PHENOTYPE_CLASSES

    _ALIASES = {
        "viable": "viable",
        "non-essential": "viable",
        "nonessential": "viable",
        "non_essential": "viable",
        "lethal": "lethal",
        "inviable": "lethal",
        "essential": "lethal",
    }

    def fit(self, labels):
        unknown = sorted({self._normalise(v) for v in labels} - set(self.classes_))
        if unknown:
            raise ValueError(
                "Unsupported phenotype label(s): "
                + ", ".join(unknown)
                + ". Expected viable/lethal labels."
            )
        return self

    def transform(self, labels):
        mapping = {"viable": 0, "lethal": 1}
        return [mapping[self._normalise(v)] for v in labels]

    def fit_transform(self, labels):
        self.fit(labels)
        return self.transform(labels)

    def inverse_transform(self, values):
        return [self.classes_[int(v)] for v in values]

    def _normalise(self, value):
        label = str(value).strip().lower()
        return self._ALIASES.get(label, label)


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
                    match = ATTRIBUTE_RE.match(line)
                    if not match:
                        raise ValueError(f"Malformed ARFF attribute declaration: {line}")
                    attr_name = next(value for value in match.groups()[:3] if value is not None)
                    attribute_names.append(attr_name)
                    attribute_types.append(match.group(4).strip())
                    continue
                if line.lower().startswith('@data'):
                    in_data_section = True
                    continue
                if in_data_section:
                    values = [v.strip().strip('"\'') for v in next(csv.reader([line]))]
                    if len(values) == len(attribute_names):
                        data_lines.append(values)
                    else:
                        raise ValueError(
                            f"ARFF row has {len(values)} values but expected "
                            f"{len(attribute_names)}: {line[:120]}"
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
    if len(attribute_names) != len(set(attribute_names)):
        duplicates = sorted({name for name in attribute_names if attribute_names.count(name) > 1})
        raise ValueError("Duplicate ARFF attribute(s): " + ", ".join(duplicates))

    import pandas as pd

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

    if df.shape[1] < 3:
        raise ValueError("ARFF must contain a gene ID, at least one GO feature, and a class")
    if df.columns.duplicated().any():
        raise ValueError("Duplicate ARFF feature names are not permitted")

    gene_column = df.columns[0]
    if df[gene_column].astype(str).duplicated().any():
        duplicate_ids = sorted(set(
            df.loc[df[gene_column].astype(str).duplicated(False), gene_column].astype(str)
        ))
        for gene in duplicate_ids:
            rows = df[df[gene_column].astype(str) == gene]
            if not rows.eq(rows.iloc[0]).all().all():
                raise ValueError(f"Conflicting duplicate ARFF rows for gene {gene}")
        logger.warning("Removing %d exact duplicate gene IDs", len(duplicate_ids))
        df = df.drop_duplicates(subset=[gene_column], keep="first").copy()

    gene_names = df.iloc[:, 0].astype(str)
    phenotype  = df.iloc[:, -1]
    go_terms   = df.iloc[:, 1:-1]

    logger.debug(f"Gene names column: {df.columns[0]}")
    logger.debug(f"Phenotype column:  {df.columns[-1]}")
    logger.debug(f"Sample gene names: {gene_names.head().tolist()}")
    logger.debug(f"Sample phenotypes: {phenotype.head().tolist()}")
    logger.debug(f"Unique phenotypes: {phenotype.unique()}")

    import numpy as np

    try:
        X = go_terms.astype(float).values
    except (TypeError, ValueError) as exc:
        raise ValueError("GO features must be numeric binary values") from exc
    if not np.isfinite(X).all():
        raise ValueError("GO features contain missing or non-finite values")
    unique_values = set(np.unique(X))
    if not unique_values <= {0.0, 1.0}:
        raise ValueError(f"GO features must be binary 0/1; found {sorted(unique_values)[:10]}")

    le = PhenotypeLabelEncoder()
    y  = np.asarray(le.fit_transform(phenotype), dtype=int)

    logger.info(f"Features shape:    {X.shape}")
    logger.info(f"Number of GO terms:{X.shape[1]}")
    logger.info(f"Class distribution:{np.bincount(y)}")
    logger.info(f"Class labels:      {le.classes_}")

    sparsity = (X == 0).sum() / X.size
    logger.info(f"Data sparsity:     {sparsity:.2%}")

    logger.debug(f"Unique feature values: {np.unique(X)}")

    return X, y, gene_names, le
