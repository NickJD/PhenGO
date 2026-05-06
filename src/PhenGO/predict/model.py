"""
predict/model.py -- Keras model architectures for PhenGO-Predict.
"""
import logging

logger = logging.getLogger(__name__)

try:
    from tensorflow import keras
    from tensorflow.keras import layers, regularizers
    from tensorflow.keras.metrics import Precision, Recall
except ImportError as e:
    logger.error(f"TensorFlow import failed: {e}")
    raise


def create_model_sparse_optimised(input_dim, hidden_units=None, dropout_rate=0.3,
                                   l2_reg=0.001):
    """Neural network optimised for sparse binary GO-term features.

    Includes L2 regularisation and batch normalisation for stable training on
    high-dimensional sparse ARFF data.
    """
    if hidden_units is None:
        hidden_units = [128, 64]

    model = keras.Sequential([
        layers.Input(shape=(input_dim,)),
        layers.Dense(hidden_units[0], activation="relu",
                     kernel_regularizer=regularizers.l2(l2_reg), name="hidden1"),
        layers.BatchNormalization(),
        layers.Dropout(dropout_rate),
        layers.Dense(hidden_units[1], activation="relu",
                     kernel_regularizer=regularizers.l2(l2_reg), name="hidden2"),
        layers.BatchNormalization(),
        layers.Dropout(dropout_rate),
        layers.Dense(32, activation="relu",
                     kernel_regularizer=regularizers.l2(l2_reg), name="hidden3"),
        layers.Dropout(dropout_rate / 2),
        layers.Dense(1, activation="sigmoid", name="output"),
    ])
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        loss="binary_crossentropy",
        metrics=["accuracy", Precision(name="precision"), Recall(name="recall")],
    )
    return model


def create_model(input_dim, hidden_units=None, dropout_rate=0.3):
    """Simple feedforward neural network without batch-norm or L2 regularisation."""
    if hidden_units is None:
        hidden_units = [128, 64]

    model = keras.Sequential([
        layers.Input(shape=(input_dim,)),
        layers.Dense(hidden_units[0], activation="relu", name="hidden1"),
        layers.Dropout(dropout_rate),
        layers.Dense(hidden_units[1], activation="relu", name="hidden2"),
        layers.Dropout(dropout_rate),
        layers.Dense(32, activation="relu", name="hidden3"),
        layers.Dropout(dropout_rate / 2),
        layers.Dense(1, activation="sigmoid", name="output"),
    ])
    model.compile(
        optimizer="adam",
        loss="binary_crossentropy",
        metrics=["accuracy", Precision(name="precision"), Recall(name="recall")],
    )
    return model

