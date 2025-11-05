# visualization.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class PerformanceVisualizer:
    """
    Visualize model performance metrics and results.
    """
    
    def __init__(self, output_dir: Path):
        """
        Initialize visualizer.
        
        Parameters:
        -----------
        output_dir : Path
            Directory to save visualizations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set style
        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_palette("husl")
    
    def plot_cv_results(self, cv_results: Dict):
        """Plot cross-validation results."""
        
        # Extract metrics across folds
        fold_metrics = pd.DataFrame([
            fold['metrics'] for fold in cv_results['fold_results']
        ])
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Plot 1: F1 Score across folds
        ax = axes[0, 0]
        ax.bar(range(len(fold_metrics)), fold_metrics['f1_score'])
        ax.axhline(cv_results['test_results']['metrics']['f1_score'], 
                  color='red', linestyle='--', label='Test F1')
        ax.set_xlabel('Fold')
        ax.set_ylabel('F1 Score')
        ax.set_title('F1 Score by Fold')
        ax.legend()
        
        # Plot 2: Multiple metrics comparison
        ax = axes[0, 1]
        metrics_to_plot = ['accuracy', 'precision', 'recall', 'f1_score']
        x = np.arange(len(metrics_to_plot))
        
        val_means = [cv_results['summary']['validation'][f'val_{m}_mean'] 
                    for m in metrics_to_plot]
        val_stds = [cv_results['summary']['validation'][f'val_{m}_std'] 
                   for m in metrics_to_plot]
        test_values = [cv_results['test_results']['metrics'][m] 
                      for m in metrics_to_plot]
        
        ax.bar(x - 0.2, val_means, 0.4, yerr=val_stds, 
              label='Validation', alpha=0.8)
        ax.bar(x + 0.2, test_values, 0.4, 
              label='Test', alpha=0.8)
        
        ax.set_xlabel('Metric')
        ax.set_ylabel('Score')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_to_plot, rotation=45)
        ax.set_title('Metrics Comparison')
        ax.legend()
        
        # Plot 3: Confusion Matrix (Test)
        ax = axes[0, 2]
        cm = np.array(cv_results['test_results']['metrics']['confusion_matrix'])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('Actual')
        ax.set_title('Test Confusion Matrix')
        
        # Plot 4: ROC AUC across folds
        ax = axes[1, 0]
        ax.boxplot([fold_metrics['roc_auc']])
        ax.scatter([1], [cv_results['test_results']['metrics']['roc_auc']], 
                  color='red', s=100, label='Test ROC AUC')
        ax.set_ylabel('ROC AUC')
        ax.set_title('ROC AUC Distribution')
        ax.legend()
        
        # Plot 5: Sample distribution
        ax = axes[1, 1]
        n_samples = [fold['metrics']['n_samples'] 
                    for fold in cv_results['fold_results']]
        n_positive = [fold['metrics']['n_positive'] 
                     for fold in cv_results['fold_results']]
        
        ax.bar(range(len(n_samples)), n_samples, label='Total')
        ax.bar(range(len(n_positive)), n_positive, label='Positive')
        ax.set_xlabel('Fold')
        ax.set_ylabel('Number of Samples')
        ax.set_title('Sample Distribution')
        ax.legend()
        
        # Plot 6: Feature Importance (top 10)
        ax = axes[1, 2]
        if hasattr(cv_results.get('final_model'), 'feature_importance'):
            top_features = cv_results['final_model'].feature_importance.head(10)
            ax.barh(range(len(top_features)), top_features['importance'])
            ax.set_yticks(range(len(top_features)))
            ax.set_yticklabels(top_features['feature'], fontsize=8)
            ax.set_xlabel('Importance')
            ax.set_title('Top 10 Feature Importances')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'cv_results.png', dpi=150, bbox_inches='tight')
        plt.show()
    
    def plot_spatial_predictions(self, predictions: np.ndarray, 
                                true_labels: np.ndarray,
                                shape: Tuple[int, int],
                                title: str = "Spatial Predictions"):
        """
        Plot spatial distribution of predictions.
        
        Parameters:
        -----------
        predictions : np.ndarray
            Predicted labels (flattened)
        true_labels : np.ndarray
            True labels (flattened)
        shape : Tuple[int, int]
            Original spatial shape (H, W)
        """
        # Reshape to spatial dimensions
        pred_spatial = predictions.reshape(shape)
        true_spatial = true_labels.reshape(shape)
        
        # Create figure
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # True labels
        im1 = axes[0].imshow(true_spatial, cmap='Blues', vmin=0, vmax=1)
        axes[0].set_title('True Flood Extent')
        plt.colorbar(im1, ax=axes[0])
        
        # Predictions
        im2 = axes[1].imshow(pred_spatial, cmap='Blues', vmin=0, vmax=1)
        axes[1].set_title('Predicted Flood Extent')
        plt.colorbar(im2, ax=axes[1])
        
        # Difference
        diff = pred_spatial - true_spatial
        im3 = axes[2].imshow(diff, cmap='RdBu', vmin=-1, vmax=1)
        axes[2].set_title('Prediction Error')
        plt.colorbar(im3, ax=axes[2])
        
        for ax in axes:
            ax.axis('off')
        
        plt.suptitle(title)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'spatial_predictions.png', 
                   dpi=150, bbox_inches='tight')
        plt.show()
    
    def create_interactive_metrics_dashboard(self, cv_results: Dict):
        """Create interactive dashboard with Plotly."""
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('F1 Score Across Folds', 
                          'Metrics Comparison',
                          'Feature Importance', 
                          'Confusion Matrix'),
            specs=[[{'type': 'bar'}, {'type': 'bar'}],
                  [{'type': 'bar'}, {'type': 'heatmap'}]]
        )
        
        # F1 scores
        fold_f1 = [fold['metrics']['f1_score'] 
                  for fold in cv_results['fold_results']]
        fig.add_trace(
            go.Bar(y=fold_f1, name='Validation F1'),
            row=1, col=1
        )
        
        # Metrics comparison
        metrics = ['accuracy', 'precision', 'recall', 'f1_score']
        val_means = [cv_results['summary']['validation'][f'val_{m}_mean'] 
                    for m in metrics]
        test_values = [cv_results['test_results']['metrics'][m] 
                      for m in metrics]
        
        fig.add_trace(
            go.Bar(x=metrics, y=val_means, name='Validation'),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(x=metrics, y=test_values, name='Test'),
            row=1, col=2
        )
        
        # Feature importance
        if hasattr(cv_results.get('final_model'), 'feature_importance'):
            top_features = cv_results['final_model'].feature_importance.head(10)
            fig.add_trace(
                go.Bar(x=top_features['importance'], 
                      y=top_features['feature'],
                      orientation='h'),
                row=2, col=1
            )
        
        # Confusion matrix
        cm = np.array(cv_results['test_results']['metrics']['confusion_matrix'])
        fig.add_trace(
            go.Heatmap(z=cm, colorscale='Blues'),
            row=2, col=2
        )
        
        fig.update_layout(height=800, showlegend=True, 
                         title_text="Model Performance Dashboard")
        
        fig.write_html(self.output_dir / 'interactive_dashboard.html')
        return fig
    
    def save_metrics_report(self, cv_results: Dict):
        """Save detailed metrics report."""
        
        report = {
            'summary': cv_results['summary'],
            'fold_details': [],
            'test_details': {
                'metrics': cv_results['test_results']['metrics'],
                'watershed': cv_results['test_results'].get('test_watershed')
            }
        }
        
        # Add fold details
        for fold_result in cv_results['fold_results']:
            report['fold_details'].append({
                'fold': fold_result['fold'],
                'watershed': fold_result.get('val_watershed'),
                'metrics': fold_result['metrics']
            })
        
        # Save as JSON
        with open(self.output_dir / 'metrics_report.json', 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        # Save as formatted text
        with open(self.output_dir / 'metrics_report.txt', 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("FLOOD PREDICTION MODEL PERFORMANCE REPORT\n")
            f.write("=" * 60 + "\n\n")
            
            f.write("VALIDATION PERFORMANCE (Cross-Validation)\n")
            f.write("-" * 40 + "\n")
            for metric, value in report['summary']['validation'].items():
                f.write(f"{metric:20s}: {value:.4f}\n")
            
            f.write("\nTEST PERFORMANCE\n")
            f.write("-" * 40 + "\n")
            for metric, value in report['test_details']['metrics'].items():
                if metric != 'confusion_matrix':
                    f.write(f"{metric:20s}: {value}\n")
            
            f.write("\nCONFUSION MATRIX (Test Set)\n")
            f.write("-" * 40 + "\n")
            cm = np.array(report['test_details']['metrics']['confusion_matrix'])
            f.write(f"{'':10s} {'Pred 0':10s} {'Pred 1':10s}\n")
            f.write(f"{'True 0':10s} {cm[0,0]:10d} {cm[0,1]:10d}\n")
            f.write(f"{'True 1':10s} {cm[1,0]:10d} {cm[1,1]:10d}\n")
