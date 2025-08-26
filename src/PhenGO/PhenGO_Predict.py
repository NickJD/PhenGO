import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt


def load_arff_data(filepath):
    """Load and parse ARFF file manually to handle string attributes"""
    print(f"Loading ARFF file: {filepath}")

    data_lines = []
    attribute_names = []
    attribute_types = []
    in_data_section = False

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('%'):
                continue

            # Check for relation
            if line.lower().startswith('@relation'):
                continue

            # Check for attributes
            if line.lower().startswith('@attribute'):
                parts = line.split()
                attr_name = parts[1].strip('"\'')  # Remove quotes if present
                attr_type = ' '.join(parts[2:]).strip()
                attribute_names.append(attr_name)
                attribute_types.append(attr_type)
                continue

            # Check for data section
            if line.lower().startswith('@data'):
                in_data_section = True
                continue

            # Process data lines
            if in_data_section:
                # Split by comma, handling quoted strings
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

                # Add the last value
                values.append(current_value.strip().strip('"\''))

                if len(values) == len(attribute_names):
                    data_lines.append(values)
                else:
                    print(f"Warning: Line has {len(values)} values but expected {len(attribute_names)}")

    # Create DataFrame
    df = pd.DataFrame(data_lines, columns=attribute_names)

    # Convert numeric columns
    for i, (col_name, col_type) in enumerate(zip(attribute_names, attribute_types)):
        if col_type.lower() in ['numeric', 'real', 'integer'] or col_type == '{0,1}':
            df[col_name] = pd.to_numeric(df[col_name], errors='coerce')

    print(f"Data shape: {df.shape}")
    print(f"Columns: {list(df.columns[:5])}... (showing first 5)")
    print(f"Attribute types: {set(attribute_types)}")

    return df, {'attribute_names': attribute_names, 'attribute_types': attribute_types}


def prepare_data(df):
    """Prepare data for neural network training"""
    print("\nPreparing data...")

    # Assume first column is gene name, last column is phenotype
    gene_names = df.iloc[:, 0]
    phenotype = df.iloc[:, -1]
    go_terms = df.iloc[:, 1:-1]  # All columns between first and last

    print(f"Gene names column: {df.columns[0]}")
    print(f"Phenotype column: {df.columns[-1]}")
    print(f"Sample gene names: {gene_names.head().tolist()}")
    print(f"Sample phenotypes: {phenotype.head().tolist()}")
    print(f"Unique phenotypes: {phenotype.unique()}")

    # Convert GO terms to numeric (should already be 0/1 but ensure it)
    X = go_terms.astype(float).values

    # Handle potential missing values
    X = np.nan_to_num(X, nan=0.0)

    # Encode phenotype labels
    le = LabelEncoder()
    y = le.fit_transform(phenotype)

    print(f"Features shape: {X.shape}")
    print(f"Number of GO terms: {X.shape[1]}")
    print(f"Class distribution: {np.bincount(y)}")
    print(f"Class labels: {le.classes_}")

    # Check sparsity
    sparsity = (X == 0).sum() / X.size
    print(f"Data sparsity: {sparsity:.2%}")

    # Verify binary data
    unique_vals = np.unique(X)
    print(f"Unique values in features: {unique_vals}")

    return X, y, gene_names, le


def create_model(input_dim, hidden_units=[128, 64], dropout_rate=0.3):
    """Create a simple feedforward neural network"""
    model = keras.Sequential([
        layers.Input(shape=(input_dim,)),

        # First hidden layer
        layers.Dense(hidden_units[0], activation='relu', name='hidden1'),
        layers.Dropout(dropout_rate),

        # Second hidden layer
        layers.Dense(hidden_units[1], activation='relu', name='hidden2'),
        layers.Dropout(dropout_rate),

        # Optional third layer for complex patterns
        layers.Dense(32, activation='relu', name='hidden3'),
        layers.Dropout(dropout_rate / 2),

        # Output layer
        layers.Dense(1, activation='sigmoid', name='output')
    ])

    model.compile(
        optimizer='adam',
        loss='binary_crossentropy',
        metrics=['accuracy', 'precision', 'recall']
    )

    return model


def plot_training_history(history):
    """Plot training history"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Accuracy
    axes[0, 0].plot(history.history['accuracy'], label='Training')
    axes[0, 0].plot(history.history['val_accuracy'], label='Validation')
    axes[0, 0].set_title('Model Accuracy')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Accuracy')
    axes[0, 0].legend()

    # Loss
    axes[0, 1].plot(history.history['loss'], label='Training')
    axes[0, 1].plot(history.history['val_loss'], label='Validation')
    axes[0, 1].set_title('Model Loss')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Loss')
    axes[0, 1].legend()

    # Precision
    axes[1, 0].plot(history.history['precision'], label='Training')
    axes[1, 0].plot(history.history['val_precision'], label='Validation')
    axes[1, 0].set_title('Model Precision')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Precision')
    axes[1, 0].legend()

    # Recall
    axes[1, 1].plot(history.history['recall'], label='Training')
    axes[1, 1].plot(history.history['val_recall'], label='Validation')
    axes[1, 1].set_title('Model Recall')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('Recall')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig('training_history.png', dpi=300, bbox_inches='tight')
    plt.show()


def analyse_feature_importance(model, feature_names, X_test, y_test, top_n=20):
    """Analyse which GO terms are most important for predicting lethal vs viable"""
    print(f"\nAnalysing feature importance (top {top_n} GO terms for each class)...")

    # Method 1: Simple weight-based analysis from first layer
    weights = model.get_layer('hidden1').get_weights()[0]  # Shape: (n_features, n_neurons)

    # Calculate simple importance as mean absolute weight across all first layer neurons
    feature_contributions = np.mean(np.abs(weights), axis=1)

    # Method 2: Permutation-based importance for class-specific analysis
    baseline_preds = model.predict(X_test, verbose=0).flatten()
    baseline_acc = np.mean((baseline_preds > 0.5) == y_test)

    lethal_importance = []
    viable_importance = []
    overall_importance = []

    print("Computing permutation importance for each GO term...")

    for i in range(X_test.shape[1]):
        # Create copy of test data
        X_permuted = X_test.copy()

        # Permute feature i
        np.random.shuffle(X_permuted[:, i])

        # Get predictions with permuted feature
        permuted_preds = model.predict(X_permuted, verbose=0).flatten()
        permuted_acc = np.mean((permuted_preds > 0.5) == y_test)

        # Overall importance
        importance = baseline_acc - permuted_acc
        overall_importance.append(importance)

        # Class-specific importance
        # For lethal class (assuming 1 = lethal)
        lethal_mask = y_test == 1
        if np.sum(lethal_mask) > 0:
            lethal_baseline_acc = np.mean((baseline_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
            lethal_permuted_acc = np.mean((permuted_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
            lethal_importance.append(lethal_baseline_acc - lethal_permuted_acc)
        else:
            lethal_importance.append(0)

        # For viable class (assuming 0 = viable)
        viable_mask = y_test == 0
        if np.sum(viable_mask) > 0:
            viable_baseline_acc = np.mean((baseline_preds[viable_mask] > 0.5) == y_test[viable_mask])
            viable_permuted_acc = np.mean((permuted_preds[viable_mask] > 0.5) == y_test[viable_mask])
            viable_importance.append(viable_baseline_acc - viable_permuted_acc)
        else:
            viable_importance.append(0)

    # Create comprehensive importance DataFrame
    importance_df = pd.DataFrame({
        'GO_Term_Index': range(len(overall_importance)),
        'Overall_Importance': overall_importance,
        'Lethal_Importance': lethal_importance,
        'Viable_Importance': viable_importance,
        'Weight_Contribution': feature_contributions
    })

    if feature_names is not None and len(feature_names) == len(overall_importance):
        importance_df['GO_Term'] = feature_names

    # Sort by different criteria
    overall_df = importance_df.sort_values('Overall_Importance', ascending=False)
    lethal_df = importance_df.sort_values('Lethal_Importance', ascending=False)
    viable_df = importance_df.sort_values('Viable_Importance', ascending=False)

    print(f"\n=== TOP {top_n} GO TERMS FOR OVERALL PREDICTION ===")
    print(overall_df[['GO_Term', 'Overall_Importance', 'Lethal_Importance', 'Viable_Importance']].head(top_n))

    print(f"\n=== TOP {top_n} GO TERMS FOR PREDICTING LETHAL PHENOTYPE ===")
    print(lethal_df[['GO_Term', 'Lethal_Importance', 'Overall_Importance']].head(top_n))

    print(f"\n=== TOP {top_n} GO TERMS FOR PREDICTING VIABLE PHENOTYPE ===")
    print(viable_df[['GO_Term', 'Viable_Importance', 'Overall_Importance']].head(top_n))

    # Create comprehensive visualisation
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Overall importance
    top_overall = overall_df.head(top_n)
    axes[0, 0].barh(range(top_n), top_overall['Overall_Importance'].values)
    axes[0, 0].set_xlabel('Permutation Importance')
    axes[0, 0].set_title(f'Top {top_n} GO Terms - Overall Importance')
    axes[0, 0].set_yticks(range(min(top_n, len(top_overall))))
    if 'GO_Term' in importance_df.columns:
        axes[0, 0].set_yticklabels([f"{name[:30]}..." if len(str(name)) > 30 else str(name)
                                    for name in top_overall['GO_Term'].head(top_n)], fontsize=8)
    axes[0, 0].invert_yaxis()

    # Lethal importance
    top_lethal = lethal_df.head(top_n)
    axes[0, 1].barh(range(top_n), top_lethal['Lethal_Importance'].values, color='red', alpha=0.7)
    axes[0, 1].set_xlabel('Lethal Prediction Importance')
    axes[0, 1].set_title(f'Top {top_n} GO Terms - Lethal Phenotype Prediction')
    axes[0, 1].set_yticks(range(min(top_n, len(top_lethal))))
    if 'GO_Term' in importance_df.columns:
        axes[0, 1].set_yticklabels([f"{name[:30]}..." if len(str(name)) > 30 else str(name)
                                    for name in top_lethal['GO_Term'].head(top_n)], fontsize=8)
    axes[0, 1].invert_yaxis()

    # Viable importance
    top_viable = viable_df.head(top_n)
    axes[1, 0].barh(range(top_n), top_viable['Viable_Importance'].values, color='green', alpha=0.7)
    axes[1, 0].set_xlabel('Viable Prediction Importance')
    axes[1, 0].set_title(f'Top {top_n} GO Terms - Viable Phenotype Prediction')
    axes[1, 0].set_yticks(range(min(top_n, len(top_viable))))
    if 'GO_Term' in importance_df.columns:
        axes[1, 0].set_yticklabels([f"{name[:30]}..." if len(str(name)) > 30 else str(name)
                                    for name in top_viable['GO_Term'].head(top_n)], fontsize=8)
    axes[1, 0].invert_yaxis()

    # Comparison plot
    # Select top features from each category for comparison
    comparison_features = set(list(top_lethal['GO_Term_Index'].head(10)) +
                              list(top_viable['GO_Term_Index'].head(10)))
    comparison_df = importance_df[importance_df['GO_Term_Index'].isin(comparison_features)]

    x_pos = np.arange(len(comparison_df))
    width = 0.35

    axes[1, 1].bar(x_pos - width / 2, comparison_df['Lethal_Importance'], width,
                   label='Lethal Importance', color='red', alpha=0.7)
    axes[1, 1].bar(x_pos + width / 2, comparison_df['Viable_Importance'], width,
                   label='Viable Importance', color='green', alpha=0.7)

    axes[1, 1].set_xlabel('GO Terms')
    axes[1, 1].set_ylabel('Importance Score')
    axes[1, 1].set_title('Lethal vs Viable Importance Comparison')
    axes[1, 1].set_xticks(x_pos)
    if 'GO_Term' in importance_df.columns:
        axes[1, 1].set_xticklabels([f"GO_{idx}" for idx in comparison_df['GO_Term_Index']],
                                   rotation=45, ha='right')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig('class_specific_feature_importance.png', dpi=300, bbox_inches='tight')
    plt.show()

    return overall_df, lethal_df, viable_df


def main():
    parser = argparse.ArgumentParser(description='Neural Network for Gene Essentiality Prediction')
    parser.add_argument('--arff_file', help='Path to ARFF file containing gene data')
    parser.add_argument('--epochs', type=int, default=100, help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for training')
    parser.add_argument('--test_size', type=float, default=0.2, help='Proportion of data for testing')
    parser.add_argument('--hidden_units', nargs='+', type=int, default=[128, 64],
                        help='Hidden layer sizes')
    parser.add_argument('--dropout', type=float, default=0.3, help='Dropout rate')

    args = parser.parse_args()

    # Load data
    df, meta = load_arff_data(args.arff_file)
    X, y, gene_names, label_encoder = prepare_data(df)

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=42, stratify=y
    )

    print(f"\nTraining set: {X_train.shape[0]} samples")
    print(f"Test set: {X_test.shape[0]} samples")

    # Create model
    model = create_model(
        input_dim=X.shape[1],
        hidden_units=args.hidden_units,
        dropout_rate=args.dropout
    )

    print(f"\nModel architecture:")
    model.summary()

    # Train model
    print(f"\nTraining model for {args.epochs} epochs...")

    # Callbacks
    early_stopping = keras.callbacks.EarlyStopping(
        monitor='val_loss', patience=15, restore_best_weights=True
    )

    reduce_lr = keras.callbacks.ReduceLROnPlateau(
        monitor='val_loss', factor=0.5, patience=10, min_lr=1e-6
    )

    history = model.fit(
        X_train, y_train,
        epochs=args.epochs,
        batch_size=args.batch_size,
        validation_split=0.2,
        callbacks=[early_stopping, reduce_lr],
        verbose=1
    )

    # Evaluate model
    print(f"\nEvaluating model...")
    test_loss, test_acc, test_precision, test_recall = model.evaluate(X_test, y_test, verbose=0)

    # Predictions
    y_pred_proba = model.predict(X_test, verbose=0)
    y_pred = (y_pred_proba > 0.5).astype(int).flatten()

    # Calculate AUC
    auc_score = roc_auc_score(y_test, y_pred_proba)

    print(f"\nTest Results:")
    print(f"Accuracy: {test_acc:.4f}")
    print(f"Precision: {test_precision:.4f}")
    print(f"Recall: {test_recall:.4f}")
    print(f"AUC: {auc_score:.4f}")

    print(f"\nDetailed Classification Report:")
    print(classification_report(y_test, y_pred, target_names=label_encoder.classes_))

    print(f"\nConfusion Matrix:")
    print(confusion_matrix(y_test, y_pred))

    # Plot training history
    plot_training_history(history)

    # Analyse feature importance with class-specific analysis
    go_term_columns = df.columns[1:-1]  # GO term column names
    overall_importance, lethal_importance, viable_importance = analyse_feature_importance(
        model, go_term_columns, X_test, y_test
    )

    # Save model and results
    model.save('gene_essentiality_model.keras')
    overall_importance.to_csv('overall_feature_importance.csv', index=False)
    lethal_importance.to_csv('lethal_feature_importance.csv', index=False)
    viable_importance.to_csv('viable_feature_importance.csv', index=False)

    print(f"\nModel saved as 'gene_essentiality_model.keras'")
    print(f"Feature importance saved as:")
    print(f"  - overall_feature_importance.csv")
    print(f"  - lethal_feature_importance.csv")
    print(f"  - viable_feature_importance.csv")


if __name__ == "__main__":
    main()