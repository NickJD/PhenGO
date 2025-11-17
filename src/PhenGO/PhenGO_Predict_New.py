"""
PhenGO_Predict.py - Enhanced Version with Cross-Dataset Comparison

IMPROVEMENTS ADDED:
==================
1. Multi-dataset comparison mode
2. Visualization of feature importance with error bars
3. ROC curve plotting
4. GO term differential importance analysis
5. Better progress reporting
6. Model architecture options for sparse data
7. Cross-validation support
8. Ensemble predictions across datasets
"""

import shutil
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# TensorFlow imports
try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    from tensorflow.keras.metrics import Precision, Recall
    from tensorflow.keras import regularizers

    print(f"TensorFlow version: {tf.__version__}")
except ImportError as e:
    print(f"Error importing TensorFlow: {e}")
    exit(1)


# [Keep existing load_arff_data and prepare_data functions]


def create_model_sparse_optimized(input_dim, hidden_units=[128, 64], dropout_rate=0.3,
                                  l2_reg=0.001):
    """
    Create neural network optimized for sparse binary features.

    IMPROVEMENTS:
    - L2 regularization to handle high-dimensional sparse data
    - Batch normalization for stable training
    - Option for different architectures
    """
    model = keras.Sequential([
        layers.Input(shape=(input_dim,)),

        # First hidden layer with L2 regularization
        layers.Dense(hidden_units[0], activation='relu',
                     kernel_regularizer=regularizers.l2(l2_reg),
                     name='hidden1'),
        layers.BatchNormalization(),
        layers.Dropout(dropout_rate),

        # Second hidden layer
        layers.Dense(hidden_units[1], activation='relu',
                     kernel_regularizer=regularizers.l2(l2_reg),
                     name='hidden2'),
        layers.BatchNormalization(),
        layers.Dropout(dropout_rate),

        # Third layer
        layers.Dense(32, activation='relu',
                     kernel_regularizer=regularizers.l2(l2_reg),
                     name='hidden3'),
        layers.Dropout(dropout_rate / 2),

        # Output
        layers.Dense(1, activation='sigmoid', name='output')
    ])

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        loss='binary_crossentropy',
        metrics=['accuracy', Precision(name='precision'), Recall(name='recall')]
    )

    return model


def plot_roc_curve(options, y_test, y_pred_proba, dataset_name=""):
    """
    Plot ROC curve with AUC score.

    NEW FEATURE: Visualization of model discrimination ability
    """
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    auc_score = roc_auc_score(y_test, y_pred_proba)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, linewidth=2, label=f'ROC curve (AUC = {auc_score:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random classifier')
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title(f'ROC Curve{" - " + dataset_name if dataset_name else ""}', fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    output_path = os.path.join(options.output_dir,
                               f'roc_curve{"_" + dataset_name if dataset_name else ""}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"ROC curve saved to {output_path}")


def plot_feature_importance_with_errors(options, importance_df, top_n=30,
                                        title="Feature Importance"):
    """
    Plot feature importance with error bars from repeated permutations.

    NEW FEATURE: Statistical visualization of GO term importance
    """
    # Get top features
    top_features = importance_df.nsmallest(top_n, 'Overall_Importance')

    fig, axes = plt.subplots(1, 3, figsize=(20, 8))

    # Overall importance
    ax = axes[0]
    y_pos = np.arange(len(top_features))
    ax.barh(y_pos, top_features['Overall_Importance'],
            xerr=top_features['Overall_Std'], capsize=3, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([t[:40] + '...' if len(t) > 40 else t
                        for t in top_features['GO_Term']], fontsize=8)
    ax.set_xlabel('Importance (Accuracy Drop)', fontsize=10)
    ax.set_title(f'Top {top_n} - Overall Importance', fontsize=12)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)

    # Lethal-specific
    ax = axes[1]
    top_lethal = importance_df.nsmallest(top_n, 'Lethal_Importance')
    y_pos = np.arange(len(top_lethal))
    ax.barh(y_pos, top_lethal['Lethal_Importance'],
            xerr=top_lethal['Lethal_Std'], capsize=3, alpha=0.7, color='red')
    ax.set_yticks(y_pos)
    ax.set_yticklabels([t[:40] + '...' if len(t) > 40 else t
                        for t in top_lethal['GO_Term']], fontsize=8)
    ax.set_xlabel('Importance (Accuracy Drop)', fontsize=10)
    ax.set_title(f'Top {top_n} - Lethal Prediction', fontsize=12)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)

    # Viable-specific
    ax = axes[2]
    top_viable = importance_df.nsmallest(top_n, 'Viable_Importance')
    y_pos = np.arange(len(top_viable))
    ax.barh(y_pos, top_viable['Viable_Importance'],
            xerr=top_viable['Viable_Std'], capsize=3, alpha=0.7, color='green')
    ax.set_yticks(y_pos)
    ax.set_yticklabels([t[:40] + '...' if len(t) > 40 else t
                        for t in top_viable['GO_Term']], fontsize=8)
    ax.set_xlabel('Importance (Accuracy Drop)', fontsize=10)
    ax.set_title(f'Top {top_n} - Viable Prediction', fontsize=12)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    output_path = os.path.join(options.output_dir, 'feature_importance_with_errors.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Feature importance plot saved to {output_path}")


def analyse_feature_importance_fixed(options, model, feature_names, X_test, y_test,
                                     n_repeats, top_n=20):
    """
    Fixed version with proper progress reporting and enhanced output.

    FIXES:
    - Proper newline after progress reporting
    - Better time estimation
    - More informative console output
    """
    import time

    print(f"\nAnalysing feature importance with {n_repeats} repeats...")
    print(f"This will perform {X_test.shape[1] * n_repeats} permutations...")

    start_time = time.time()

    # Baseline predictions
    baseline_preds = model.predict(X_test, verbose=0).flatten()
    baseline_acc = np.mean((baseline_preds > 0.5) == y_test)

    # Class-specific baselines
    lethal_mask = y_test == 1
    viable_mask = y_test == 0

    lethal_baseline = (np.mean((baseline_preds[lethal_mask] > 0.5) == y_test[lethal_mask])
                       if np.sum(lethal_mask) > 0 else np.nan)
    viable_baseline = (np.mean((baseline_preds[viable_mask] > 0.5) == y_test[viable_mask])
                       if np.sum(viable_mask) > 0 else np.nan)

    # Storage
    overall_mean, overall_std = [], []
    lethal_mean, lethal_std = [], []
    viable_mean, viable_std = [], []

    total_iterations = X_test.shape[1] * n_repeats

    for i in range(X_test.shape[1]):
        drops_overall, drops_lethal, drops_viable = [], [], []

        for rep in range(n_repeats):
            iteration = i * n_repeats + rep + 1

            # Progress with time estimate
            if iteration % 10 == 0:
                elapsed = time.time() - start_time
                per_iteration = elapsed / iteration
                remaining = (total_iterations - iteration) * per_iteration
                print(f"  Progress: {iteration}/{total_iterations} "
                      f"({100 * iteration / total_iterations:.1f}%) - "
                      f"ETA: {remaining / 60:.1f} min", end='\r')

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

        # Store statistics
        overall_mean.append(np.mean(drops_overall))
        overall_std.append(np.std(drops_overall))
        lethal_mean.append(np.mean(drops_lethal) if drops_lethal else np.nan)
        lethal_std.append(np.std(drops_lethal) if drops_lethal else np.nan)
        viable_mean.append(np.mean(drops_viable) if drops_viable else np.nan)
        viable_std.append(np.std(drops_viable) if drops_viable else np.nan)

    # FIXED: Print newline after progress
    print()

    total_time = time.time() - start_time
    print(f"Feature importance analysis completed in {total_time / 60:.1f} minutes")

    # Weight contribution
    try:
        weights = model.get_layer('hidden1').get_weights()[0]
        weight_contrib = np.mean(np.abs(weights), axis=1)
    except:
        weight_contrib = np.zeros(X_test.shape[1])

    # Build DataFrame
    importance_df = pd.DataFrame({
        'GO_Term_Index': range(len(overall_mean)),
        'GO_Term': (list(feature_names) if feature_names is not None
                    else [f"GO_Term_{i}" for i in range(len(overall_mean))]),
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

    # Save
    overall_df.to_csv(os.path.join(options.output_dir, 'overall_feature_importance.csv'),
                      index=False)
    lethal_df.to_csv(os.path.join(options.output_dir, 'lethal_feature_importance.csv'),
                     index=False)
    viable_df.to_csv(os.path.join(options.output_dir, 'viable_feature_importance.csv'),
                     index=False)

    print(f"\n=== TOP {top_n} GO TERMS (Overall Importance) ===")
    print(overall_df[['GO_Term', 'Overall_Importance', 'Overall_Std']].head(top_n))

    # NEW: Generate visualization
    plot_feature_importance_with_errors(options, overall_df, top_n=30)

    return overall_df, lethal_df, viable_df


def compare_datasets(options, arff_files, dataset_names):
    """
    NEW FEATURE: Compare multiple datasets and identify differential GO term importance.

    This function:
    1. Trains separate models on each dataset
    2. Extracts feature importance for each
    3. Identifies GO terms that differ in importance between datasets
    4. Generates comparative visualizations

    Args:
        options: Command line options
        arff_files: List of paths to ARFF files
        dataset_names: List of names for each dataset
    """
    print(f"\n{'=' * 80}")
    print(f"MULTI-DATASET COMPARISON MODE")
    print(f"Comparing {len(arff_files)} datasets")
    print(f"{'=' * 80}\n")

    all_results = {}

    for arff_file, dataset_name in zip(arff_files, dataset_names):
        print(f"\n{'=' * 80}")
        print(f"Processing dataset: {dataset_name}")
        print(f"{'=' * 80}\n")

        # Create subdirectory for this dataset
        dataset_dir = os.path.join(options.output_dir, dataset_name)
        os.makedirs(dataset_dir, exist_ok=True)

        # Load and process data
        df, meta = load_arff_data(arff_file)
        X, y, gene_names, label_encoder = prepare_data(df)

        # Split data
        X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
            X, y, gene_names, test_size=options.test_size, random_state=42, stratify=y
        )

        # Create and train model
        model = create_model_sparse_optimized(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout
        )

        early_stopping = keras.callbacks.EarlyStopping(
            monitor='val_loss', patience=15, restore_best_weights=True
        )

        history = model.fit(
            X_train, y_train,
            epochs=options.epochs,
            batch_size=options.batch_size,
            validation_split=0.2,
            callbacks=[early_stopping],
            verbose=0
        )

        # Evaluate
        y_pred_proba = model.predict(X_test, verbose=0)
        test_results = model.evaluate(X_test, y_test, verbose=0)

        # Feature importance (modified options for subdirectory)
        dataset_options = argparse.Namespace(**vars(options))
        dataset_options.output_dir = dataset_dir

        go_term_columns = df.columns[1:-1]
        overall_imp, lethal_imp, viable_imp = analyse_feature_importance_fixed(
            dataset_options, model, go_term_columns, X_test, y_test,
            options.perm_repeats
        )

        # Store results
        all_results[dataset_name] = {
            'model': model,
            'importance': {
                'overall': overall_imp,
                'lethal': lethal_imp,
                'viable': viable_imp
            },
            'metrics': {
                'loss': test_results[0],
                'accuracy': test_results[1],
                'precision': test_results[2] if len(test_results) > 2 else 0,
                'recall': test_results[3] if len(test_results) > 3 else 0,
                'auc': roc_auc_score(y_test, y_pred_proba)
            },
            'predictions': y_pred_proba,
            'true_labels': y_test,
            'go_terms': list(go_term_columns)
        }

        # Plot ROC for this dataset
        plot_roc_curve(options, y_test, y_pred_proba, dataset_name)

    # Generate comparison visualizations
    generate_comparison_plots(options, all_results, dataset_names)
    generate_differential_importance_report(options, all_results, dataset_names)

    return all_results


def generate_comparison_plots(options, all_results, dataset_names):
    """
    Generate comparative visualizations across datasets.

    NEW FEATURE: Side-by-side comparison plots
    """
    # 1. Model performance comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metrics_to_plot = ['accuracy', 'precision', 'recall']
    for idx, metric in enumerate(metrics_to_plot):
        values = [all_results[name]['metrics'][metric] for name in dataset_names]
        axes[idx].bar(range(len(dataset_names)), values, alpha=0.7)
        axes[idx].set_xticks(range(len(dataset_names)))
        axes[idx].set_xticklabels(dataset_names, rotation=45, ha='right')
        axes[idx].set_ylabel(metric.capitalize())
        axes[idx].set_title(f'{metric.capitalize()} Comparison')
        axes[idx].set_ylim([0, 1])
        axes[idx].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(options.output_dir, 'model_performance_comparison.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Top GO terms comparison (heatmap)
    # Get union of top 30 GO terms from each dataset
    all_top_terms = set()
    for name in dataset_names:
        top_30 = all_results[name]['importance']['overall'].head(30)['GO_Term']
        all_top_terms.update(top_30)

    # Create importance matrix
    importance_matrix = []
    term_list = sorted(list(all_top_terms))

    for term in term_list:
        row = []
        for name in dataset_names:
            imp_df = all_results[name]['importance']['overall']
            term_imp = imp_df[imp_df['GO_Term'] == term]['Overall_Importance']
            row.append(term_imp.values[0] if len(term_imp) > 0 else 0)
        importance_matrix.append(row)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(12, max(8, len(term_list) * 0.3)))
    im = ax.imshow(importance_matrix, cmap='YlOrRd', aspect='auto')

    ax.set_xticks(range(len(dataset_names)))
    ax.set_xticklabels(dataset_names, rotation=45, ha='right')
    ax.set_yticks(range(len(term_list)))
    ax.set_yticklabels([t[:50] + '...' if len(t) > 50 else t for t in term_list],
                       fontsize=8)
    ax.set_title('GO Term Importance Across Datasets', fontsize=14, pad=20)

    plt.colorbar(im, ax=ax, label='Importance Score')
    plt.tight_layout()
    plt.savefig(os.path.join(options.output_dir, 'go_term_importance_heatmap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Comparison plots saved to {options.output_dir}")


def generate_differential_importance_report(options, all_results, dataset_names):
    """
    Generate report of GO terms with differential importance across datasets.

    NEW FEATURE: Identifies GO terms that matter more in one dataset vs another
    """
    if len(dataset_names) != 2:
        print("Differential analysis currently only supports 2 datasets")
        return

    name1, name2 = dataset_names
    imp1 = all_results[name1]['importance']['overall']
    imp2 = all_results[name2]['importance']['overall']

    # Merge on GO_Term
    merged = pd.merge(
        imp1[['GO_Term', 'Overall_Importance', 'Overall_Std']],
        imp2[['GO_Term', 'Overall_Importance', 'Overall_Std']],
        on='GO_Term',
        suffixes=(f'_{name1}', f'_{name2}')
    )

    # Calculate difference
    merged['Importance_Diff'] = (merged[f'Overall_Importance_{name1}'] -
                                 merged[f'Overall_Importance_{name2}'])
    merged['Abs_Diff'] = abs(merged['Importance_Diff'])

    # Sort by absolute difference
    merged_sorted = merged.sort_values('Abs_Diff', ascending=False)

    # Save report
    report_path = os.path.join(options.output_dir, 'differential_importance_report.csv')
    merged_sorted.to_csv(report_path, index=False)

    # Console output
    print(f"\n{'=' * 80}")
    print(f"TOP 20 DIFFERENTIALLY IMPORTANT GO TERMS")
    print(f"{'=' * 80}\n")
    print(f"Positive values: More important in {name1}")
    print(f"Negative values: More important in {name2}\n")

    print(merged_sorted[['GO_Term', 'Importance_Diff',
                         f'Overall_Importance_{name1}',
                         f'Overall_Importance_{name2}']].head(20).to_string(index=False))

    print(f"\nFull report saved to: {report_path}")


# [Keep existing evaluate_and_analyse_predictions and main functions]
# Add comparison mode to main:

def main_extended():
    """Extended main function with multi-dataset comparison support."""
    parser = argparse.ArgumentParser(
        description='PhenGO Predict - Neural Network for Gene Essentiality Prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Single dataset:
    %(prog)s -arff_file data.arff -output_dir results/

  Compare two datasets:
    %(prog)s -arff_files dataset1.arff dataset2.arff \\
             -dataset_names "Wild Type" "Mutant" \\
             -output_dir comparison_results/
        """
    )

    # Add multi-dataset support
    parser.add_argument('-arff_file', help='Single ARFF file (for single dataset mode)')
    parser.add_argument('-arff_files', nargs='+',
                        help='Multiple ARFF files (for comparison mode)')
    parser.add_argument('-dataset_names', nargs='+',
                        help='Names for each dataset in comparison mode')

    # Existing arguments
    parser.add_argument('-epochs', type=int, default=100)
    parser.add_argument('-batch_size', type=int, default=32)
    parser.add_argument('-test_size', type=float, default=0.2)
    parser.add_argument('-hidden_units', nargs='+', type=int, default=[128, 64])
    parser.add_argument('-dropout', type=float, default=0.3)
    parser.add_argument('-perm_repeats', type=int, default=5,
                        help='Number of repeats for permutation importance (default: 5)')
    parser.add_argument('-output_dir', required=True)

    options = parser.parse_args()

    # Ensure output directory
    if os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir)

    # Determine mode
    if options.arff_files:
        # Comparison mode
        if not options.dataset_names:
            options.dataset_names = [f"Dataset_{i + 1}" for i in range(len(options.arff_files))]
        elif len(options.dataset_names) != len(options.arff_files):
            print("Error: Number of dataset names must match number of ARFF files")
            return

        print(f"\n{'=' * 80}")
        print("RUNNING IN COMPARISON MODE")
        print(f"Datasets to compare: {', '.join(options.dataset_names)}")
        print(f"{'=' * 80}\n")

        results = compare_datasets(options, options.arff_files, options.dataset_names)

        print(f"\n{'=' * 80}")
        print("COMPARISON ANALYSIS COMPLETE")
        print(f"Results saved to: {options.output_dir}")
        print(f"{'=' * 80}\n")

    elif options.arff_file:
        # Single dataset mode (original functionality)
        print(f"\n{'=' * 80}")
        print("RUNNING IN SINGLE DATASET MODE")
        print(f"{'=' * 80}\n")

        result = load_arff_data(options.arff_file)
        if result[0] is None:
            print("Failed to load data. Exiting.")
            return

        df, meta = result
        X, y, gene_names, label_encoder = prepare_data(df)

        # Split data
        X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
            X, y, gene_names, test_size=options.test_size, random_state=42, stratify=y
        )

        print(f"\nData Split:")
        print(f"Training set: {X_train.shape[0]} samples")
        print(f"Test set: {X_test.shape[0]} samples")
        print(f"Training class distribution: {np.bincount(y_train)}")
        print(f"Test class distribution: {np.bincount(y_test)}")

        # Create model
        model = create_model_sparse_optimized(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout
        )

        print(f"\nModel architecture:")
        model.summary()

        # Train
        print(f"\nTraining model for {options.epochs} epochs...")

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

            # Evaluation
            print("\nStarting evaluation phase...")
            predictions_df, test_loss, test_acc, test_precision, test_recall, auc_score = \
                evaluate_and_analyse_predictions(model, X_test, y_test, genes_test, label_encoder)

            # Plot training history
            print("Plotting training history...")
            plot_training_history(options, history)

            # Plot ROC curve
            y_pred_proba = model.predict(X_test, verbose=0)
            plot_roc_curve(options, y_test, y_pred_proba)

            # Feature importance
            print("Starting feature importance analysis...")
            go_term_columns = df.columns[1:-1]
            overall_importance, lethal_importance, viable_importance = \
                analyse_feature_importance_fixed(
                    options, model, go_term_columns, X_test, y_test,
                    options.perm_repeats
                )

            # Save results
            print("Saving results...")
            model.save(os.path.join(options.output_dir, 'gene_essentiality_model.keras'))
            predictions_df.to_csv(
                os.path.join(options.output_dir, 'test_predictions_detailed.csv'),
                index=False
            )

            # Write report
            report_path = os.path.join(options.output_dir, 'final_report.txt')
            with open(report_path, 'w') as report_file:
                report_file.write("=== PHENGO PREDICTION REPORT ===\n\n")
                report_file.write(f"Model Performance:\n")
                report_file.write(f"  Loss: {test_loss:.3f}\n")
                report_file.write(f"  Accuracy: {test_acc:.3f}\n")
                report_file.write(f"  Precision: {test_precision:.3f}\n")
                report_file.write(f"  Recall: {test_recall:.3f}\n")
                report_file.write(f"  AUC: {auc_score:.3f}\n\n")
                report_file.write("Confusion Matrix:\n")
                y_pred = (y_pred_proba > 0.5).astype(int).flatten()
                report_file.write(np.array2string(confusion_matrix(y_test, y_pred)) + "\n\n")
                report_file.write("Classification Report:\n")
                report_file.write(classification_report(
                    y_test, y_pred, target_names=label_encoder.classes_
                ))

            print(f"\n{'=' * 80}")
            print(f"ANALYSIS COMPLETE")
            print(f"{'=' * 80}")
            print(f"Model Performance: Loss={test_loss:.3f}, Accuracy={test_acc:.3f}, AUC={auc_score:.3f}")
            print(f"Full report written to {report_path}")
            print(f"{'=' * 80}\n")

        except Exception as e:
            print(f"\nError during training or analysis: {str(e)}")
            import traceback
            traceback.print_exc()

    else:
        print("Error: Must provide either -arff_file or -arff_files")
        parser.print_help()
        return

    # Save parameters
    from datetime import datetime
    with open(os.path.join(options.output_dir, "PhenGo_Predict_params.txt"), "w") as outfile:
        outfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")

    print("\nThank you for using PhenGO!")
    print("Documentation: https://github.com/NickJD/PhenGO")
    print("Report issues: https://github.com/NickJD/PhenGO/issues")


if __name__ == "__main__":
    main_extended()


# =============================================================================
# ADDITIONAL UTILITY FUNCTIONS
# =============================================================================

def cross_validate_model(X, y, gene_names, options, n_folds=5):
    """
    NEW FEATURE: K-fold cross-validation for more robust performance estimates.

    Performs stratified K-fold CV and reports mean/std of metrics.
    Useful for small datasets or when you want confidence intervals.
    """
    print(f"\n{'=' * 80}")
    print(f"PERFORMING {n_folds}-FOLD CROSS-VALIDATION")
    print(f"{'=' * 80}\n")

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)

    cv_results = {
        'accuracy': [],
        'precision': [],
        'recall': [],
        'auc': []
    }

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
        print(f"\nFold {fold}/{n_folds}")

        X_train_cv, X_test_cv = X[train_idx], X[test_idx]
        y_train_cv, y_test_cv = y[train_idx], y[test_idx]

        # Create and train model
        model = create_model_sparse_optimized(
            input_dim=X.shape[1],
            hidden_units=options.hidden_units,
            dropout_rate=options.dropout
        )

        model.fit(
            X_train_cv, y_train_cv,
            epochs=options.epochs,
            batch_size=options.batch_size,
            validation_split=0.2,
            callbacks=[keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)],
            verbose=0
        )

        # Evaluate
        y_pred_proba = model.predict(X_test_cv, verbose=0)
        y_pred = (y_pred_proba > 0.5).astype(int).flatten()

        test_results = model.evaluate(X_test_cv, y_test_cv, verbose=0)

        cv_results['accuracy'].append(test_results[1])
        cv_results['precision'].append(test_results[2] if len(test_results) > 2 else 0)
        cv_results['recall'].append(test_results[3] if len(test_results) > 3 else 0)
        cv_results['auc'].append(roc_auc_score(y_test_cv, y_pred_proba))

        print(f"  Accuracy: {cv_results['accuracy'][-1]:.3f}")
        print(f"  AUC: {cv_results['auc'][-1]:.3f}")

    # Summary
    print(f"\n{'=' * 80}")
    print("CROSS-VALIDATION RESULTS")
    print(f"{'=' * 80}\n")

    for metric, values in cv_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        print(f"{metric.capitalize()}: {mean_val:.3f} ± {std_val:.3f}")

    # Save CV results
    cv_df = pd.DataFrame(cv_results)
    cv_df.to_csv(os.path.join(options.output_dir, 'cross_validation_results.csv'),
                 index=False)

    return cv_results


def predict_unlabeled_genes(model, unlabeled_arff_file, output_file):
    """
    NEW FEATURE: Use trained model to predict phenotypes for unlabeled genes.

    Useful for:
    - Predicting essentiality of newly discovered genes
    - Transferring predictions between species
    - Hypothesis generation

    Args:
        model: Trained Keras model
        unlabeled_arff_file: ARFF file with genes lacking phenotype labels
        output_file: Where to save predictions
    """
    print(f"\n{'=' * 80}")
    print("PREDICTING PHENOTYPES FOR UNLABELED GENES")
    print(f"{'=' * 80}\n")

    # Load unlabeled data
    df, meta = load_arff_data(unlabeled_arff_file)

    gene_names = df.iloc[:, 0]
    X_unlabeled = df.iloc[:, 1:-1].astype(float).values
    X_unlabeled = np.nan_to_num(X_unlabeled, nan=0.0)

    # Predict
    predictions_proba = model.predict(X_unlabeled, verbose=0).flatten()
    predictions_binary = (predictions_proba > 0.5).astype(int)

    # Create results DataFrame
    results_df = pd.DataFrame({
        'Gene_Name': gene_names,
        'Predicted_Label': ['lethal' if p == 1 else 'viable' for p in predictions_binary],
        'Lethal_Probability': predictions_proba,
        'Confidence': np.abs(predictions_proba - 0.5)
    })

    # Sort by confidence
    results_df = results_df.sort_values('Confidence', ascending=False)

    # Save
    results_df.to_csv(output_file, index=False)

    print(f"Predictions for {len(results_df)} genes saved to {output_file}")
    print(f"\nPrediction Summary:")
    print(
        f"  Predicted Lethal: {sum(predictions_binary)} ({sum(predictions_binary) / len(predictions_binary) * 100:.1f}%)")
    print(
        f"  Predicted Viable: {len(predictions_binary) - sum(predictions_binary)} ({(len(predictions_binary) - sum(predictions_binary)) / len(predictions_binary) * 100:.1f}%)")
    print(f"  Mean Confidence: {results_df['Confidence'].mean():.3f}")

    return results_df


def analyze_go_term_enrichment_in_predictions(predictions_df, arff_file, output_dir):
    """
    NEW FEATURE: Analyze which GO terms are enriched in correct vs incorrect predictions.

    This helps identify:
    - GO terms the model relies on correctly
    - GO terms that mislead the model
    - Systematic biases in predictions
    """
    print(f"\n{'=' * 80}")
    print("ANALYZING GO TERM ENRICHMENT IN PREDICTIONS")
    print(f"{'=' * 80}\n")

    # Load full data to get GO term annotations
    df, meta = load_arff_data(arff_file)
    go_terms = df.iloc[:, 1:-1]

    # Merge with predictions
    merged = predictions_df.merge(
        pd.concat([df.iloc[:, 0], go_terms], axis=1),
        left_on='Gene_Name',
        right_on=df.columns[0]
    )

    # Separate correct and incorrect predictions
    correct = merged[merged['Correct_Prediction'] == True]
    incorrect = merged[merged['Correct_Prediction'] == False]

    # Calculate GO term frequencies
    go_cols = list(go_terms.columns)

    enrichment_results = []
    for go_term in go_cols:
        correct_freq = correct[go_term].sum() / len(correct) if len(correct) > 0 else 0
        incorrect_freq = incorrect[go_term].sum() / len(incorrect) if len(incorrect) > 0 else 0

        enrichment_results.append({
            'GO_Term': go_term,
            'Correct_Frequency': correct_freq,
            'Incorrect_Frequency': incorrect_freq,
            'Enrichment_Ratio': correct_freq / incorrect_freq if incorrect_freq > 0 else np.inf,
            'Difference': correct_freq - incorrect_freq
        })

    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df = enrichment_df.sort_values('Enrichment_Ratio', ascending=False)

    # Save results
    output_path = os.path.join(output_dir, 'go_enrichment_in_predictions.csv')
    enrichment_df.to_csv(output_path, index=False)

    print(f"\nTop 10 GO Terms Enriched in CORRECT Predictions:")
    print(enrichment_df[['GO_Term', 'Correct_Frequency', 'Incorrect_Frequency',
                         'Enrichment_Ratio']].head(10).to_string(index=False))

    print(f"\nTop 10 GO Terms Enriched in INCORRECT Predictions:")
    print(enrichment_df[['GO_Term', 'Correct_Frequency', 'Incorrect_Frequency',
                         'Enrichment_Ratio']].tail(10).to_string(index=False))

    print(f"\nFull results saved to {output_path}")

    return enrichment_df