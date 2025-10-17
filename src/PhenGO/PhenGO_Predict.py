import shutil
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt

# TensorFlow imports with error handling
try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    from tensorflow.keras.metrics import Precision, Recall

    print(f"TensorFlow version: {tf.__version__}")
except ImportError as e:
    print(f"Error importing TensorFlow: {e}")
    print("Please install TensorFlow: pip install tensorflow")
    exit(1)


def load_arff_data(filepath):
    """Load and parse ARFF file manually to handle string attributes"""
    print(f"Loading ARFF file: {filepath}")

    data_lines = []
    attribute_names = []
    attribute_types = []
    in_data_section = False

    try:
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
                    if len(parts) >= 3:
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

    except FileNotFoundError:
        print(f"Error: File {filepath} not found")
        return None, None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None, None

    if not data_lines:
        print("Error: No data found in ARFF file")
        return None, None

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
        metrics=['accuracy', Precision(name='precision'), Recall(name='recall')]
    )

    return model


def plot_training_history(options, history):
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
    if 'precision' in history.history:
        axes[1, 0].plot(history.history['precision'], label='Training')
        axes[1, 0].plot(history.history['val_precision'], label='Validation')
        axes[1, 0].set_title('Model Precision')
        axes[1, 0].set_xlabel('Epoch')
        axes[1, 0].set_ylabel('Precision')
        axes[1, 0].legend()

    # Recall
    if 'recall' in history.history:
        axes[1, 1].plot(history.history['recall'], label='Training')
        axes[1, 1].plot(history.history['val_recall'], label='Validation')
        axes[1, 1].set_title('Model Recall')
        axes[1, 1].set_xlabel('Epoch')
        axes[1, 1].set_ylabel('Recall')
        axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(options.output_dir,'training_history.png'), dpi=300, bbox_inches='tight')
    print("Training history plot saved as 'training_history.png'")
    plt.close()  # Close to free memory

def analyse_feature_importance(options, model, feature_names, X_test, y_test, n_repeats, top_n=20):
    """Permutation-based feature importance with class-specific breakdown"""

    print(f"\nAnalysing feature importance with {n_repeats} repeats...")

    # Baseline predictions
    baseline_preds = model.predict(X_test, verbose=0).flatten()
    baseline_acc = np.mean((baseline_preds > 0.5) == y_test)

    # Class-specific baselines
    lethal_mask = y_test == 1
    viable_mask = y_test == 0

    if np.sum(lethal_mask) > 0:
        lethal_baseline = np.mean((baseline_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
    else:
        lethal_baseline = np.nan

    if np.sum(viable_mask) > 0:
        viable_baseline = np.mean((baseline_preds[viable_mask] > 0.5) == y_test[viable_mask])
    else:
        viable_baseline = np.nan

    # Storage
    overall_mean, overall_std = [], []
    lethal_mean, lethal_std = [], []
    viable_mean, viable_std = [], []

    print("Computing permutation importance for each GO terme...")

    for i in range(X_test.shape[1]):
        drops_overall, drops_lethal, drops_viable = [], [], []

        for _ in range(n_repeats):
            print(f"  Permuting GO term {i+1}/{X_test.shape[1]}, repeat {_+1}/{n_repeats}", end='\r')
            X_permuted = X_test.copy()
            np.random.shuffle(X_permuted[:, i])

            permuted_preds = model.predict(X_permuted, verbose=0).flatten()
            perm_acc = np.mean((permuted_preds > 0.5) == y_test)

            drops_overall.append(baseline_acc - perm_acc)

            if np.sum(lethal_mask) > 0:
                perm_lethal = np.mean((permuted_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
                drops_lethal.append(lethal_baseline - perm_lethal)

            if np.sum(viable_mask) > 0:
                perm_viable = np.mean((permuted_preds[viable_mask] > 0.5) == y_test[viable_mask])
                drops_viable.append(viable_baseline - perm_viable)

        # Store means and stds
        overall_mean.append(np.mean(drops_overall))
        overall_std.append(np.std(drops_overall))
        lethal_mean.append(np.mean(drops_lethal) if drops_lethal else np.nan)
        lethal_std.append(np.std(drops_lethal) if drops_lethal else np.nan)
        viable_mean.append(np.mean(drops_viable) if drops_viable else np.nan)
        viable_std.append(np.std(drops_viable) if drops_viable else np.nan)

    # Weight heuristic (optional)
    try:
        weights = model.get_layer('hidden1').get_weights()[0]
        weight_contrib = np.mean(np.abs(weights), axis=1)
    except Exception:
        weight_contrib = np.zeros(X_test.shape[1])

    # Build DataFrame
    importance_df = pd.DataFrame({
        'GO_Term_Index': range(len(overall_mean)),
        'GO_Term': list(feature_names) if feature_names is not None else [f"GO_Term_{i}" for i in range(len(overall_mean))],
        'Overall_Importance': overall_mean,
        'Overall_Std': overall_std,
        'Lethal_Importance': lethal_mean,
        'Lethal_Std': lethal_std,
        'Viable_Importance': viable_mean,
        'Viable_Std': viable_std,
        'Weight_Contribution': weight_contrib
    })

    # Sort copies
    overall_df = importance_df.sort_values('Overall_Importance', ascending=False)
    lethal_df = importance_df.sort_values('Lethal_Importance', ascending=False)
    viable_df = importance_df.sort_values('Viable_Importance', ascending=False)

    # Save outputs
    overall_df.to_csv(os.path.join(options.output_dir,'overall_feature_importance.csv'), index=False)
    lethal_df.to_csv(os.path.join(options.output_dir,'lethal_feature_importance.csv'), index=False)
    viable_df.to_csv(os.path.join(options.output_dir,'viable_feature_importance.csv'), index=False)

    print(f"\n=== TOP {top_n} GO TERMS (Overall Importance) ===")
    print(overall_df[['GO_Term', 'Overall_Importance', 'Overall_Std']].head(top_n))

    return overall_df, lethal_df, viable_df


#


def evaluate_and_analyse_predictions(model, X_test, y_test, gene_names_test, label_encoder):
    """Comprehensive evaluation with detailed prediction analysis"""
    print(f"\nEvaluating model and analysing predictions...")

    # Get predictions
    y_pred_proba = model.predict(X_test, verbose=0)
    y_pred = (y_pred_proba > 0.5).astype(int).flatten()

    # Basic metrics
    test_results = model.evaluate(X_test, y_test, verbose=0)
    test_loss = test_results[0]
    test_acc = test_results[1]
    test_precision = test_results[2] if len(test_results) > 2 else 0
    test_recall = test_results[3] if len(test_results) > 3 else 0

    try:
        auc_score = roc_auc_score(y_test, y_pred_proba)
    except Exception as e:
        print(f"Warning: Could not calculate AUC score: {e}")
        auc_score = 0.0

    print(f"\nTest Results:")
    print(f"Loss: {test_loss:.4f}")
    print(f"Accuracy: {test_acc:.4f}")
    print(f"Precision: {test_precision:.4f}")
    print(f"Recall: {test_recall:.4f}")
    print(f"AUC: {auc_score:.4f}")

    print(f"\nDetailed Classification Report:")
    print(classification_report(y_test, y_pred, target_names=label_encoder.classes_))

    print(f"\nConfusion Matrix:")
    cm = confusion_matrix(y_test, y_pred)
    print(cm)

    # Create detailed predictions DataFrame
    predictions_df = pd.DataFrame({
        'Gene_Name': gene_names_test.reset_index(drop=True),
        'True_Label': [label_encoder.classes_[i] for i in y_test],
        'True_Label_Numeric': y_test,
        'Predicted_Label': [label_encoder.classes_[i] for i in y_pred],
        'Predicted_Label_Numeric': y_pred,
        'Prediction_Probability': y_pred_proba.flatten(),
        'Confidence': np.abs(y_pred_proba.flatten() - 0.5),  # Distance from decision boundary
        'Correct_Prediction': y_test == y_pred
    })

    # Add prediction categories
    predictions_df['Prediction_Category'] = predictions_df.apply(
        lambda row: 'True Positive' if (row['True_Label_Numeric'] == 1 and row['Predicted_Label_Numeric'] == 1)
        else 'True Negative' if (row['True_Label_Numeric'] == 0 and row['Predicted_Label_Numeric'] == 0)
        else 'False Positive' if (row['True_Label_Numeric'] == 0 and row['Predicted_Label_Numeric'] == 1)
        else 'False Negative', axis=1
    )

    # Sort by confidence (most confident predictions first)
    predictions_df = predictions_df.sort_values('Confidence', ascending=False)

    print(f"\n=== PREDICTION SUMMARY ===")
    print(f"Total test samples: {len(predictions_df)}")
    print(f"Correct predictions: {predictions_df['Correct_Prediction'].sum()}")
    print(f"Incorrect predictions: {(~predictions_df['Correct_Prediction']).sum()}")

    print(f"\nPrediction Categories:")
    for category in predictions_df['Prediction_Category'].unique():
        count = (predictions_df['Prediction_Category'] == category).sum()
        print(f"  {category}: {count}")

    # Show most confident correct predictions
    print(f"\n=== TOP 10 MOST CONFIDENT CORRECT PREDICTIONS ===")
    correct_preds = predictions_df[predictions_df['Correct_Prediction'] == True].head(10)
    print(correct_preds[['Gene_Name', 'True_Label', 'Prediction_Probability', 'Confidence']].to_string(index=False))

    # Show most confident incorrect predictions (these are interesting!)
    print(f"\n=== TOP 10 MOST CONFIDENT INCORRECT PREDICTIONS ===")
    incorrect_preds = predictions_df[predictions_df['Correct_Prediction'] == False].head(10)
    if len(incorrect_preds) > 0:
        print(incorrect_preds[
            ['Gene_Name', 'True_Label', 'Predicted_Label', 'Prediction_Probability', 'Confidence']].to_string(
            index=False))
    else:
        print("No incorrect predictions found!")

    # Show predictions near decision boundary (uncertain predictions)
    print(f"\n=== 10 MOST UNCERTAIN PREDICTIONS (near 0.5 probability) ===")
    uncertain_preds = predictions_df.nsmallest(10, 'Confidence')
    print(uncertain_preds[
        ['Gene_Name', 'True_Label', 'Predicted_Label', 'Prediction_Probability', 'Confidence']].to_string(
        index=False))

    return predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score


def main():
    parser = argparse.ArgumentParser(description='Neural Network for Gene Essentiality Prediction')
    parser.add_argument('-arff_file', required=True, help='Path to ARFF file containing gene data')
    parser.add_argument('-epochs', type=int, default=100, help='Number of training epochs')
    parser.add_argument('-batch_size', type=int, default=32, help='Batch size for training')
    parser.add_argument('-test_size', type=float, default=0.2, help='Proportion of data for testing')
    parser.add_argument('-hidden_units', nargs='+', type=int, default=[128, 64],
                        help='Hidden layer sizes')
    parser.add_argument('-dropout', type=float, default=0.3, help='Dropout rate')
    parser.add_argument('-perm_repeats', type=int, default=1, help='Number of repeats for permutation importance')
    parser.add_argument('-output_dir', dest="output_dir", required=True, help='Output directory (current contents will be deleted)')


    options = parser.parse_args()

    # Load data
    result = load_arff_data(options.arff_file)
    if result[0] is None:
        print("Failed to load data. Exiting.")
        return

    # Ensure output directory exists and is empty
    if os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir)

    df, meta = result
    X, y, gene_names, label_encoder = prepare_data(df)

    # Split data with more control
    X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
        X, y, gene_names, test_size=options.test_size, random_state=42, stratify=y
    )

    print(f"\nData Split:")
    print(f"Training set: {X_train.shape[0]} samples")
    print(f"Test set: {X_test.shape[0]} samples")
    print(f"Training class distribution: {np.bincount(y_train)}")
    print(f"Test class distribution: {np.bincount(y_test)}")

    # Create model
    model = create_model(
        input_dim=X.shape[1],
        hidden_units=options.hidden_units,
        dropout_rate=options.dropout
    )

    print(f"\nModel architecture:")
    model.summary()

    # Train model
    print(f"\nTraining model for {options.epochs} epochs...")

    # Callbacks
    early_stopping = keras.callbacks.EarlyStopping(
        monitor='val_loss', patience=15, restore_best_weights=True
    )

    reduce_lr = keras.callbacks.ReduceLROnPlateau(
        monitor='val_loss', factor=0.5, patience=10, min_lr=1e-6
    )

    try:
        history = model.fit(
            X_train, y_train,
            epochs=options.epochs,
            batch_size=options.batch_size,
            validation_split=0.2,
            callbacks=[early_stopping, reduce_lr],
            verbose=1
        )

        print("\nTraining completed successfully!")

        # Comprehensive evaluation and prediction analysis
        print("\nStarting evaluation phase...")
        predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score = evaluate_and_analyse_predictions(
            model, X_test, y_test, genes_test, label_encoder
        )

        # Plot training history
        print("Plotting training history...")
        plot_training_history(options, history)

        # Analyse feature importance
        print("Starting feature importance analysis...")
        go_term_columns = df.columns[1:-1]  # GO term column names
        overall_importance, lethal_importance, viable_importance = analyse_feature_importance(
            options, model, go_term_columns, X_test, y_test, options.perm_repeats
        )

        # Save model and results
        print("Saving results...")
        model.save(os.path.join(options.output_dir,'gene_essentiality_model.keras'))
        predictions_df.to_csv(os.path.join(options.output_dir,'test_predictions_detailed.csv'), index=False)
        overall_importance.to_csv(os.path.join(options.output_dir,'overall_feature_importance.csv'), index=False)
        lethal_importance.to_csv(os.path.join(options.output_dir,'lethal_feature_importance.csv'), index=False)
        viable_importance.to_csv(os.path.join(options.output_dir,'viable_feature_importance.csv'), index=False)

        print(f"\n=== FILES SAVED ===")
        print(f"Model: gene_essentiality_model.keras")
        print(f"Detailed predictions: test_predictions_detailed.csv")
        print(f"Feature importance files:")
        print(f"  - overall_feature_importance.csv")
        print(f"  - lethal_feature_importance.csv")
        print(f"  - viable_feature_importance.csv")


        # Write final summary and important metrics to a report file
        report_path = os.path.join(options.output_dir, 'final_report.txt')
        with open(report_path, 'w') as report_file:
            report_file.write("=== FINAL SUMMARY ===\n")
            report_file.write(f"Model Performance:\n")
            report_file.write(f"  Loss: {test_loss:.3f}\n")
            report_file.write(f"  Accuracy: {test_acc:.3f}\n")
            report_file.write(f"  Precision: {test_precision:.3f}\n")
            report_file.write(f"  Recall: {test_recall:.3f}\n")
            report_file.write(f"  AUC: {auc_score:.3f}\n")
            report_file.write("\nConfusion Matrix:\n")
            report_file.write(np.array2string(
                confusion_matrix(y_test, (model.predict(X_test, verbose=0) > 0.5).astype(int).flatten())) + "\n")
            report_file.write("\nClassification Report:\n")
            report_file.write(
                classification_report(y_test, (model.predict(X_test, verbose=0) > 0.5).astype(int).flatten(),
                                      target_names=label_encoder.classes_))
        print(f"\n=== FINAL SUMMARY ===")
        print(f"Model Performance: Loss={test_loss:.3f}, Accuracy={test_acc:.3f}, AUC={auc_score:.3f}")
        print(f"Full report written to {report_path}")


    except Exception as e:
        print(f"\nError during training or analysis: {str(e)}")
        import traceback
        traceback.print_exc()

        # Try basic evaluation at least
        try:
            print("\nAttempting basic evaluation...")
            test_results = model.evaluate(X_test, y_test, verbose=1)
            test_acc = test_results[1]
            print(f"Basic accuracy: {test_acc:.4f}")

            # Save model at least
            model.save(os.path.join(options.output_dir,'gene_essentiality_model_basic.keras'))
            print("Model saved successfully")

        except Exception as e2:
            print(f"Even basic evaluation failed: {str(e2)}")


if __name__ == "__main__":
    main()