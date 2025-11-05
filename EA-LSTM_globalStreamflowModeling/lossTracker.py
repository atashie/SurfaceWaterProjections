################# Downloading and Analysis for SageMaker Training Outputs
import os
import json
import tarfile
import boto3
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Optional: for timestamps etc.
from datetime import datetime

print("=" * 60)
print("DOWNLOADING SAGEMAKER TRAINING OUTPUTS")
print("=" * 60)

# -----------------------------------------------------------------------------
# Configuration (update these three as needed)
# -----------------------------------------------------------------------------
bucket_name = 'sagemaker-us-east-2-672980278503'
training_job_name = 'pytorch-training-2025-08-24-15-04-22-070'
experiment_name = 'experiment_lightning_116_pcaShuffle_yt2lyr5dialNoNoiseComplxNorm3xLoss'

# Exact S3 base prefix for this experiment (without 's3://bucket/')
# Example: s3://sagemaker-us-east-2-672980278503/training-outputs_lightning_114_pcaShuffle_nt
s3_base_prefix = 'training-outputs_lightning_116_pcaShuffle_yt2lyr5dialNoNoiseComplxNorm3xLoss'

# -----------------------------------------------------------------------------
# Local directories
# -----------------------------------------------------------------------------
local_dir = f'./model_analysis_{training_job_name}'
Path(local_dir).mkdir(exist_ok=True)

s3 = boto3.client('s3')

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def head_exists(s3_client, bucket, key) -> bool:
    """Check if s3://bucket/key exists via HEAD."""
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except Exception:
        return False

def extract_and_list_tar(tar_path, extract_to):
    """Extract tar file and list all contents with their paths."""
    print(f"\nExtracting {os.path.basename(tar_path)}...")
    with tarfile.open(tar_path, 'r:gz') as tar:
        # List all members first
        print("  Archive contents:")
        for member in tar.getmembers():
            print(f"    {member.name} ({member.size} bytes)")
        # Extract all
        tar.extractall(extract_to)

    # Walk through extracted files
    print("\n  Extracted file structure:")
    for root, dirs, files in os.walk(extract_to):
        level = root.replace(extract_to, '').count(os.sep)
        indent = '    ' * (level + 1)
        print(f'{indent}{os.path.basename(root)}/')
        sub_indent = '    ' * (level + 2)
        for file in files:
            file_path = os.path.join(root, file)
            size = os.path.getsize(file_path)
            if size < 1024:
                size_str = f"{size} B"
            elif size < 1024 * 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size / 1024 / 1024:.1f} MB"
            print(f'{sub_indent}{file} ({size_str})')

# -----------------------------------------------------------------------------
# 1) Download model.tar.gz
# -----------------------------------------------------------------------------
print("\n1. Downloading model.tar.gz...")
model_tar_path = f'{local_dir}/model.tar.gz'
model_key = f'{s3_base_prefix}/{training_job_name}/output/model.tar.gz'

if head_exists(s3, bucket_name, model_key):
    try:
        s3.download_file(bucket_name, model_key, model_tar_path)
        print(f"✓ Downloaded model.tar.gz from s3://{bucket_name}/{model_key} "
              f"({os.path.getsize(model_tar_path) / 1024 / 1024:.1f} MB)")
        model_extract_dir = f'{local_dir}/model_files'
        Path(model_extract_dir).mkdir(exist_ok=True)
        extract_and_list_tar(model_tar_path, model_extract_dir)
    except Exception as e:
        print(f"✗ Failed to download model.tar.gz: {e}")
else:
    print("✗ Error with model.tar.gz: Not Found at:")
    print(f"  s3://{bucket_name}/{model_key}")

# -----------------------------------------------------------------------------
# 2) Download output.tar.gz
# -----------------------------------------------------------------------------
print("\n2. Downloading output.tar.gz...")
output_tar_path = f'{local_dir}/output.tar.gz'
output_key = f'{s3_base_prefix}/{training_job_name}/output/output.tar.gz'

if head_exists(s3, bucket_name, output_key):
    try:
        s3.download_file(bucket_name, output_key, output_tar_path)
        print(f"✓ Downloaded output.tar.gz from s3://{bucket_name}/{output_key} "
              f"({os.path.getsize(output_tar_path) / 1024:.1f} KB)")
        output_extract_dir = f'{local_dir}/output_files'
        Path(output_extract_dir).mkdir(exist_ok=True)
        extract_and_list_tar(output_tar_path, output_extract_dir)
    except Exception as e:
        print(f"✗ Failed to download output.tar.gz: {e}")
else:
    print("✗ Error with output.tar.gz: Not Found at:")
    print(f"  s3://{bucket_name}/{output_key}")

# -----------------------------------------------------------------------------
# 3) Find and load all JSON files with a comprehensive search
# -----------------------------------------------------------------------------
print("\n" + "=" * 60)
print("LOADING TRAINING RESULTS")
print("=" * 60)

# Define all possible locations where files might be
search_dirs = [
    local_dir,
    f'{local_dir}/model_files',
    f'{local_dir}/output_files',
    f'{local_dir}/output_files/output',
    f'{local_dir}/output_files/output/data',
    f'{local_dir}/output_files/opt/ml/output/data',
    f'{local_dir}/model_files/opt/ml/model',
]

json_files_to_find = {
    'results.json': None,
    'model_info.json': None,
    'normalization_params.json': None
}

for json_file in json_files_to_find.keys():
    print(f"\nSearching for {json_file}...")
    found = False
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            search_path = Path(search_dir)
            for path in search_path.rglob(json_file):
                print(f"  Found at: {path}")
                try:
                    with open(path, 'r') as f:
                        json_files_to_find[json_file] = json.load(f)
                    print(f"  ✓ Successfully loaded {json_file}")
                    found = True
                    break
                except Exception as e:
                    print(f"  ✗ Error loading {path}: {e}")
        if found:
            break
    if not found:
        print(f"  ✗ {json_file} not found in any location")

# Safe defaults if JSONs are missing
results = json_files_to_find.get('results.json') or {}
model_info = json_files_to_find.get('model_info.json') or {}
norm_params = json_files_to_find.get('normalization_params.json') or {}

test_results = results.get('test_results', []) if isinstance(results, dict) else []
hyperparams = results.get('hyperparameters', {}) if isinstance(results, dict) else {}
data_info = results.get('data_info', {}) if isinstance(results, dict) else {}

# -----------------------------------------------------------------------------
# 4) Print training summary (console)
# -----------------------------------------------------------------------------
print("\n" + "=" * 60)
print("TRAINING SUMMARY")
print("=" * 60)

if isinstance(results, dict) and results:
    print(f"Experiment: {results.get('experiment_name', 'N/A')}")
    print(f"Total epochs trained: {results.get('total_epochs_trained', 'N/A')}")
    if test_results and isinstance(test_results, list) and len(test_results) > 0 and isinstance(test_results[0], dict):
        print(f"\nTest Results:")
        for key, value in test_results[0].items():
            if isinstance(value, (int, float)):
                print(f"  {key}: {value:.4f}")
            else:
                print(f"  {key}: {value}")

if isinstance(model_info, dict) and model_info:
    print(f"\nModel Checkpoint Info:")
    print(f"  Best Validation Score: {model_info.get('best_score', 'N/A')}")
    print(f"  Monitor metric: {model_info.get('monitor_metric', 'N/A')}")
    if 'best_epoch' in model_info:
        be = model_info.get('best_epoch', {})
        try:
            print(f"  Best checkpoints saved: {len(be)}")
        except Exception:
            pass

if isinstance(results, dict) and results:
    print(f"\nModel Architecture:")
    for key, value in hyperparams.items():
        print(f"  {key}: {value}")
    print(f"\nDataset Split:")
    for key, value in data_info.items():
        print(f"  {key}: {value}")

# -----------------------------------------------------------------------------
# 5) Locate TensorBoard logs
# -----------------------------------------------------------------------------
tb_search_dirs = [
    f'{local_dir}/output_files/{experiment_name}/lightning',
    f'{local_dir}/output_files/output/data/{experiment_name}/lightning',
    f'{local_dir}/output_files/opt/ml/output/data/{experiment_name}/lightning',
]

tb_log_dir = None
for tb_dir in tb_search_dirs:
    if os.path.exists(tb_dir):
        tb_log_dir = tb_dir
        print(f"\n✓ TensorBoard logs found at: {tb_log_dir}")
        break

# Fallback: discover event files recursively
if not tb_log_dir:
    search_root = Path(f'{local_dir}/output_files')
    if search_root.exists():
        event_files = list(search_root.rglob('events.out.tfevents.*'))
        if event_files:
            tb_log_dir = str(event_files[0].parent)
            print(f"\n✓ TensorBoard logs found at (discovered): {tb_log_dir}")

if not tb_log_dir:
    print(f"\n✗ TensorBoard logs not found. Searched in:")
    for tb_dir in tb_search_dirs:
        print(f"  - {tb_dir}")

################# Plotting
if tb_log_dir:
    from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

    def parse_tensorboard_logs_with_epoch_conversion(log_dir):
        """Extract training metrics from TensorBoard logs with proper epoch conversion."""
        ea = EventAccumulator(log_dir)
        ea.Reload()

        scalar_tags = ea.Tags().get('scalars', [])
        print(f"\nAvailable TensorBoard metrics: {len(scalar_tags)}")
        print(f"Metrics: {scalar_tags}")

        metrics = {}
        steps_per_epoch = None

        # Method 1: Use epoch metric to determine steps per epoch
        if 'epoch' in scalar_tags:
            epoch_events = ea.Scalars('epoch')
            if len(epoch_events) > 10:
                for i in range(1, len(epoch_events)):
                    if epoch_events[i].value > epoch_events[i - 1].value:
                        steps_per_epoch = epoch_events[i].step
                        break

        # Method 2: Use validation metrics gaps
        if steps_per_epoch is None and 'val_loss' in scalar_tags:
            val_events = ea.Scalars('val_loss')
            if len(val_events) > 1:
                steps_per_epoch = val_events[1].step - val_events[0].step

        # Method 3: Estimate from total steps and epochs.json (if present)
        if steps_per_epoch is None:
            max_step = 0
            for tag in scalar_tags:
                events = ea.Scalars(tag)
                if events:
                    max_step = max(max_step, max(e.step for e in events))
            total_epochs = results.get('total_epochs_trained', 0) if isinstance(results, dict) else 0
            if max_step > 0 and total_epochs:
                steps_per_epoch = max_step / total_epochs

        if steps_per_epoch is not None:
            print(f"Estimated steps per epoch: {steps_per_epoch:.1f}")
        else:
            print("Estimated steps per epoch: Unknown")

        # Extract metrics
        for tag in scalar_tags:
            events = ea.Scalars(tag)
            if not events:
                continue

            steps = [e.step for e in events]
            values = [e.value for e in events]

            # Treat validation/test metrics as epoch-level; also *_epoch tags
            if '_epoch' in tag or any(keyword in tag for keyword in ['val_', 'test_']):
                if steps and steps_per_epoch:
                    epochs = [s / steps_per_epoch for s in steps]
                    metrics[tag] = {'steps': epochs, 'values': values}
                else:
                    metrics[tag] = {'steps': steps, 'values': values}
            else:
                metrics[tag] = {'steps': steps, 'values': values}

        return metrics, steps_per_epoch

    # Parse metrics
    metrics, steps_per_epoch = parse_tensorboard_logs_with_epoch_conversion(tb_log_dir)

    # Create comprehensive training visualization
    if metrics:
        # Plot configuration
        plot_configs = [
            # Primary Loss (MSE)
            {
                'metrics': [('train_loss_epoch', 'Training Loss (MSE)'), ('val_loss', 'Validation Loss (MSE)')],
                'title': 'Loss (MSE - Primary Loss Function)',
                'ylabel': 'MSE',
                'invert_y': False
            },
            # KGE Score
            {
                'metrics': [('val_kge', 'Validation KGE')],
                'title': 'KGE Score',
                'ylabel': 'KGE',
                'invert_y': False,
                'add_lines': [(0.5, 'orange', '--', 'KGE=0.5'),
                              (0.7, 'red', '--', 'KGE=0.7')]
            },
            # NSE Score (Top-right)
            {
                'metrics': [('val_nse', 'Validation NSE')],
                'title': 'NSE Score',
                'ylabel': 'NSE',
                'invert_y': False,
                'add_lines': [(0.0, 'gray', '--', 'NSE=0.0'),
                              (0.5, 'orange', '--', 'NSE=0.5'),
                              (0.7, 'red', '--', 'NSE=0.7')]
            },
            # RMSE
            {
                'metrics': [('val_rmse', 'Validation RMSE')],
                'title': 'Root Mean Squared Error',
                'ylabel': 'RMSE (mm/day)',
                'invert_y': False
            },
            # MAE
            {
                'metrics': [('val_mae', 'Validation MAE')],
                'title': 'Mean Absolute Error',
                'ylabel': 'MAE (mm/day)',
                'invert_y': False
            },
            # KGE Component: Correlation (r)
            {
                'metrics': [('val_kge_r', 'Correlation (r)')],
                'title': 'KGE Component: Correlation Coefficient',
                'ylabel': 'r',
                'invert_y': False,
                'add_lines': [(1.0, 'green', '--', 'Perfect Correlation')]
            },
            # KGE Component: Variability (alpha)
            {
                'metrics': [('val_kge_alpha', 'Variability (α)')],
                'title': 'KGE Component: Variability Ratio',
                'ylabel': 'α (σ_sim/σ_obs)',
                'invert_y': False,
                'add_lines': [(1.0, 'green', '--', 'Perfect Variability')]
            },
            # KGE Component: Bias (beta)
            {
                'metrics': [('val_kge_beta', 'Bias (β)')],
                'title': 'KGE Component: Bias Ratio',
                'ylabel': 'β (μ_sim/μ_obs)',
                'invert_y': False,
                'add_lines': [(1.0, 'green', '--', 'No Bias')]
            },
            # Learning Rate
            {
                'metrics': [('lr-Adam', 'Learning Rate')],
                'title': 'Learning Rate Schedule',
                'ylabel': 'Learning Rate',
                'invert_y': False,
                'log_scale': True
            }
        ]

        # Compute layout
        n_plots = len(plot_configs) + 1  # +1 for summary panel
        n_cols = 3
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows))
        axes = axes.flatten()

        # Plot each configured panel
        for idx, config in enumerate(plot_configs):
            ax = axes[idx]
            plotted_something = False

            for metric_key, label in config['metrics']:
                if metric_key in metrics:
                    data = metrics[metric_key]
                    epochs = data['steps']
                    values = data['values']

                    color = 'red' if 'train' in metric_key else 'blue'
                    linestyle = '--' if 'train' in metric_key else '-'
                    ax.plot(epochs, values, color=color, linestyle=linestyle,
                            linewidth=2, label=label, alpha=0.8)
                    plotted_something = True

                    # Annotate best/closest values
                    if len(values) > 0:
                        if metric_key in ['val_kge_r', 'val_kge_alpha', 'val_kge_beta']:
                            distances_from_one = np.abs(np.array(values) - 1.0)
                            best_idx = np.argmin(distances_from_one)
                            ax.scatter(epochs[best_idx], values[best_idx],
                                       color='red', s=100, zorder=5,
                                       label=f'Closest to 1: {values[best_idx]:.3f}')
                        elif 'val_kge' in metric_key:
                            best_idx = np.argmax(values)
                            ax.scatter(epochs[best_idx], values[best_idx],
                                       color='red', s=100, zorder=5,
                                       label=f'Best: {values[best_idx]:.3f}')
                        elif metric_key == 'val_nse':
                            best_idx = np.argmax(values)
                            ax.scatter(epochs[best_idx], values[best_idx],
                                       color='red', s=100, zorder=5,
                                       label=f'Best: {values[best_idx]:.3f}')

            if plotted_something:
                ax.set_xlabel('Epoch')
                ax.set_ylabel(config['ylabel'])
                ax.set_title(config['title'])
                ax.grid(True, alpha=0.3)
                ax.legend()

                # Add reference lines
                if 'add_lines' in config:
                    for y_val, color, style, label in config['add_lines']:
                        ax.axhline(y=y_val, color=color, linestyle=style, alpha=0.5, label=label)

                # Log scale if requested
                if config.get('log_scale', False):
                    ax.set_yscale('log')

                # Invert y if requested
                if config.get('invert_y', False):
                    ax.invert_yaxis()
            else:
                ax.text(0.5, 0.5, f"No data for {config['title']}",
                        ha='center', va='center', transform=ax.transAxes)
                ax.set_title(config['title'])

        # Summary panel
        ax = axes[-1]
        ax.axis('off')

        # Gather final validation metrics (include NSE)
        final_metrics = {}
        for key in ['val_kge', 'val_rmse', 'val_mae', 'val_kge_r', 'val_kge_alpha', 'val_kge_beta', 'val_nse']:
            if key in metrics and len(metrics[key]['values']) > 0:
                final_metrics[key] = metrics[key]['values'][-1]

        # Format validation metrics
        val_kge_str = f"{final_metrics['val_kge']:.4f}" if 'val_kge' in final_metrics else 'N/A'
        val_rmse_str = f"{final_metrics['val_rmse']:.4f}" if 'val_rmse' in final_metrics else 'N/A'
        val_mae_str = f"{final_metrics['val_mae']:.4f}" if 'val_mae' in final_metrics else 'N/A'
        val_kge_r_str = f"{final_metrics['val_kge_r']:.4f}" if 'val_kge_r' in final_metrics else 'N/A'
        val_kge_alpha_str = f"{final_metrics['val_kge_alpha']:.4f}" if 'val_kge_alpha' in final_metrics else 'N/A'
        val_kge_beta_str = f"{final_metrics['val_kge_beta']:.4f}" if 'val_kge_beta' in final_metrics else 'N/A'
        val_nse_str = f"{final_metrics['val_nse']:.4f}" if 'val_nse' in final_metrics else 'N/A'

        # Test metrics
        test_kge_str = 'N/A'
        test_rmse_str = 'N/A'
        test_mae_str = 'N/A'
        test_nse_str = 'N/A'
        if test_results and len(test_results) > 0 and isinstance(test_results[0], dict):
            if 'test_kge' in test_results[0]:
                tv = test_results[0]['test_kge']
                test_kge_str = f"{tv:.4f}" if isinstance(tv, (int, float)) else str(tv)
            if 'test_rmse' in test_results[0]:
                tv = test_results[0]['test_rmse']
                test_rmse_str = f"{tv:.4f}" if isinstance(tv, (int, float)) else str(tv)
            if 'test_mae' in test_results[0]:
                tv = test_results[0]['test_mae']
                test_mae_str = f"{tv:.4f}" if isinstance(tv, (int, float)) else str(tv)
            if 'test_nse' in test_results[0]:
                tv = test_results[0]['test_nse']
                test_nse_str = f"{tv:.4f}" if isinstance(tv, (int, float)) else str(tv)

        summary_text = f"""Training Summary

Experiment: {results.get('experiment_name', 'N/A') if isinstance(results, dict) else 'N/A'}
Total Epochs: {results.get('total_epochs_trained', 'N/A') if isinstance(results, dict) else 'N/A'}

Test Performance:
  Test KGE: {test_kge_str}
  Test RMSE: {test_rmse_str}
  Test MAE: {test_mae_str}
  Test NSE: {test_nse_str}

Final Validation Metrics:
  KGE: {val_kge_str}
  RMSE: {val_rmse_str} mm/day
  MAE: {val_mae_str} mm/day
  NSE: {val_nse_str}

KGE Components:
  r (correlation): {val_kge_r_str}
  α (variability): {val_kge_alpha_str}
  β (bias): {val_kge_beta_str}

Model Configuration:
  Hidden Size: {hyperparams.get('hidden_size', 'N/A')}
  Layers: {hyperparams.get('num_layers', 'N/A')}
  Batch Size: {hyperparams.get('batch_size', 'N/A')}
  Learning Rate: {hyperparams.get('learning_rate', 'N/A')}

Dataset:
  Total: {data_info.get('total_watersheds', 'N/A')} watersheds
  Train/Val/Test: {data_info.get('train_watersheds', '?')}/{data_info.get('val_watersheds', '?')}/{data_info.get('test_watersheds', '?')}"""

        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Hide any unused subplots (if grid bigger than needed)
        used_positions = list(range(len(plot_configs))) + [len(axes) - 1]
        for idx in range(len(axes)):
            if idx not in used_positions:
                axes[idx].set_visible(False)

        plt.tight_layout()
        out_png = f'{local_dir}/training_analysis_comprehensive.png'
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"\n✓ Saved comprehensive analysis plot to {out_png}")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
print(f"All files saved to: {local_dir}/")
