import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import json
from typing import Dict, List, Tuple
from tqdm import tqdm
import warnings
import os
from pathlib import Path
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_model_and_data(model_dir, watershed_df, bucket_name, base_data_dir, base_attr_dir):
    """Load trained model and prepare test data"""
    
    # Updated paths to match your actual file structure
    model_path = os.path.join(model_dir, 'model_files', 'model.ckpt')  # Added 'model_files'
    model_info_path = os.path.join(model_dir, 'model_files', 'model_info.json')  # Added 'model_files'
    norm_params_path = os.path.join(model_dir, 'output_files', 'normalization_params.json')  # Changed from 'output_data' to 'output_files'
    results_path = os.path.join(model_dir, 'output_files', 'results.json')  # Changed from 'output_data' to 'output_files'
    
    # Load JSON files
    with open(norm_params_path, 'r') as f:
        norm_params = json.load(f)
    
    with open(model_info_path, 'r') as f:
        model_info = json.load(f)
        
    with open(results_path, 'r') as f:
        results = json.load(f)
    
    # Import classes from train.py
    import sys
    sys.path.append('.')  # Add current directory to path
    from train import StreamflowLightningModule, CaravanDataModule, EntityAwareLSTM
    
    # Extract hyperparameters from results
    hyperparams = results.get('hyperparameters', {})
    
    # Create model instance
    model = StreamflowLightningModule.load_from_checkpoint(
        model_path,
        dynamic_input_size=len(model_info['feature_info']['dynamic_features']),
        static_input_size=len(model_info['feature_info']['static_features']),
        hidden_size=hyperparams.get('hidden_size', 256),
        num_layers=hyperparams.get('num_layers', 2),
        dropout=hyperparams.get('dropout', 0.3),
        learning_rate=hyperparams.get('learning_rate', 0.0001),
        norm_params=norm_params
    )
    model.eval()
    
    # Create data module
    data_module = CaravanDataModule(
        watersheds_df=watershed_df,
        bucket_name=bucket_name,
        base_data_dir=base_data_dir,
        base_attr_dir=base_attr_dir,
        sequence_length=365,
        batch_size=32,
        num_workers=4,
        chunk_size=50,
        train_split=0.6,
        val_split=0.2,
        random_seed=42,
        norm_params=norm_params
    )
    
    data_module.setup('fit')
    
    return model, data_module, norm_params, model_info, results

def generate_predictions_for_watersheds(model, data_module, n_watersheds=4):
    """Generate predictions for the first n test watersheds"""
    
    from train import prepare_ptf_dataframe, unnormalize_predictions
    
    # Get test watersheds
    test_watersheds = data_module.test_watersheds.iloc[:n_watersheds]
    predictions_dict = {}
    
    device = next(model.parameters()).device
    
    for idx, (_, watershed) in enumerate(test_watersheds.iterrows()):
        watershed_id = f"{watershed['subdirectory_name']}_{watershed['watershedID']}"
        print(f"Generating predictions for {watershed_id} ({idx+1}/{n_watersheds})")
        
        # Load data for this watershed
        watershed_subset = pd.DataFrame([watershed])
        df, static_cols, dynamic_cols, _ = prepare_ptf_dataframe(
            watershed_subset,
            data_module.bucket_name,
            data_module.base_data_dir,
            data_module.base_attr_dir,
            norm_params=data_module.norm_params
        )
        
        if df.empty:
            continue
        
        # Prepare sequences
        sequence_length = 365
        predictions = []
        actuals = []
        dates = []
        
        # Get data for this watershed
        watershed_data = df[df['group_id'] == watershed_id].sort_values('date')
        
        # Generate predictions for each valid sequence
        for i in range(len(watershed_data) - sequence_length):
            # Extract sequence
            seq_data = watershed_data.iloc[i:i+sequence_length]
            target_data = watershed_data.iloc[i+sequence_length]
            
            # Prepare inputs
            dynamic_seq = torch.FloatTensor(seq_data[data_module.dynamic_cols_no_target].values).unsqueeze(0).to(device)
            static_feat = torch.FloatTensor(seq_data[static_cols].iloc[0].values).unsqueeze(0).to(device)
            
            # Get prediction
            with torch.no_grad():
                pred = model(dynamic_seq, static_feat)
            
            predictions.append(pred.cpu().item())
            actuals.append(target_data['streamflow'])
            dates.append(target_data['date'])
        
        # Unnormalize predictions and actuals
        predictions_array = np.array(predictions)
        actuals_array = np.array(actuals)
        
        predictions_original = unnormalize_predictions(predictions_array, data_module.norm_params)
        actuals_original = unnormalize_predictions(actuals_array, data_module.norm_params)
        
        predictions_dict[watershed_id] = {
            'dates': dates,
            'predictions': predictions_original,
            'actuals': actuals_original,
            'watershed_info': watershed.to_dict()
        }
    
    return predictions_dict

def plot_hydrographs(predictions_dict, save_dir='./hydrograph_plots'):
    """Plot hydrographs for test watersheds in various formats"""
    
    os.makedirs(save_dir, exist_ok=True)
    
    from train import calculate_kge
    
    for watershed_id, data in predictions_dict.items():
        dates = pd.to_datetime(data['dates'])
        predictions = data['predictions']
        actuals = data['actuals']
        
        # Create figure with 4 subplots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'Hydrograph Analysis - {watershed_id}', fontsize=16, fontweight='bold')
        
        # a) Linear space - full period
        ax = axes[0, 0]
        ax.plot(dates, actuals, 'b-', label='Observed', alpha=0.7, linewidth=1)
        ax.plot(dates, predictions, 'r-', label='Predicted', alpha=0.7, linewidth=1)
        ax.set_xlabel('Date')
        ax.set_ylabel('Streamflow (mm/day)')
        ax.set_title('Full Period - Linear Scale')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Calculate and display KGE for full period
        kge_full = calculate_kge(actuals, predictions)
        ax.text(0.02, 0.95, f'KGE: {kge_full:.3f}', transform=ax.transAxes, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                fontsize=12, fontweight='bold')
        
        # b) Log space - full period
        ax = axes[0, 1]
        # Add small constant to avoid log(0)
        actuals_log = actuals + 1e-6
        predictions_log = predictions + 1e-6
        ax.semilogy(dates, actuals_log, 'b-', label='Observed', alpha=0.7, linewidth=1)
        ax.semilogy(dates, predictions_log, 'r-', label='Predicted', alpha=0.7, linewidth=1)
        ax.set_xlabel('Date')
        ax.set_ylabel('Streamflow (mm/day) - Log Scale')
        ax.set_title('Full Period - Log Scale')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        # c) Linear space - last 2 years
        ax = axes[1, 0]
        last_2_years_mask = dates >= (dates.max() - pd.Timedelta(days=730))
        dates_2y = dates[last_2_years_mask]
        actuals_2y = actuals[last_2_years_mask]
        predictions_2y = predictions[last_2_years_mask]
        
        if len(dates_2y) > 0:
            ax.plot(dates_2y, actuals_2y, 'b-', label='Observed', alpha=0.7, linewidth=1.5)
            ax.plot(dates_2y, predictions_2y, 'r-', label='Predicted', alpha=0.7, linewidth=1.5)
            
            # Calculate and display KGE for last 2 years
            kge_2y = calculate_kge(actuals_2y, predictions_2y)
            ax.text(0.02, 0.95, f'KGE: {kge_2y:.3f}', transform=ax.transAxes,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                    fontsize=12, fontweight='bold')
        
        ax.set_xlabel('Date')
        ax.set_ylabel('Streamflow (mm/day)')
        ax.set_title('Last 2 Years - Linear Scale')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # d) Log space - last 2 years
        ax = axes[1, 1]
        if len(dates_2y) > 0:
            actuals_2y_log = actuals_2y + 1e-6
            predictions_2y_log = predictions_2y + 1e-6
            ax.semilogy(dates_2y, actuals_2y_log, 'b-', label='Observed', alpha=0.7, linewidth=1.5)
            ax.semilogy(dates_2y, predictions_2y_log, 'r-', label='Predicted', alpha=0.7, linewidth=1.5)
        
        ax.set_xlabel('Date')
        ax.set_ylabel('Streamflow (mm/day) - Log Scale')
        ax.set_title('Last 2 Years - Log Scale')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        plt.savefig(f'{save_dir}/hydrograph_{watershed_id}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Saved hydrograph plot for {watershed_id}")

def calculate_feature_importance_perturbation(model, data_module, n_samples=1000, debug=True):
    """Calculate feature importance using perturbation method with debugging"""
    
    from train import calculate_kge, unnormalize_predictions
    
    device = next(model.parameters()).device
    
    # Get validation dataloader
    val_loader = data_module.train_dataloader()
    
    # Collect baseline predictions
    baseline_preds = []
    all_targets = []
    all_dynamic_seqs = []
    all_static_feats = []
    
    print("Collecting baseline predictions...")
    sample_count = 0
    
    with torch.no_grad():
        for batch in val_loader:
            if sample_count >= n_samples:
                break
                
            dynamic_seq, static_feat, target = batch
            dynamic_seq = dynamic_seq.to(device)
            static_feat = static_feat.to(device)
            
            pred = model(dynamic_seq, static_feat)
            
            baseline_preds.append(pred.cpu())
            all_targets.append(target)
            all_dynamic_seqs.append(dynamic_seq.cpu())
            all_static_feats.append(static_feat.cpu())
            
            sample_count += len(target)
    
    # Concatenate all batches
    baseline_preds = torch.cat(baseline_preds)[:n_samples]
    all_targets = torch.cat(all_targets)[:n_samples]
    all_dynamic_seqs = torch.cat(all_dynamic_seqs)[:n_samples]
    all_static_feats = torch.cat(all_static_feats)[:n_samples]
    
    # Debug: Check static feature statistics
    if debug:
        print(f"\nDEBUG: Static features shape: {all_static_feats.shape}")
        print(f"Number of unique watersheds in validation: {len(data_module.val_watersheds)}")
        
        # Check variance of static features
        static_features = data_module.static_cols
        print("\nStatic feature statistics (normalized values):")
        for idx, feat_name in enumerate(static_features[:10]):  # Show first 10
            feat_values = all_static_feats[:, idx]
            print(f"  {feat_name}:")
            print(f"    Mean: {feat_values.mean():.4f}, Std: {feat_values.std():.4f}")
            print(f"    Min: {feat_values.min():.4f}, Max: {feat_values.max():.4f}")
            print(f"    Unique values: {len(torch.unique(feat_values))}")
    
    # Calculate baseline KGE
    baseline_preds_np = baseline_preds.numpy().flatten()
    targets_np = all_targets.numpy().flatten()
    
    # Unnormalize for KGE calculation
    baseline_preds_orig = unnormalize_predictions(baseline_preds_np, data_module.norm_params)
    targets_orig = unnormalize_predictions(targets_np, data_module.norm_params)
    baseline_kge = calculate_kge(targets_orig, baseline_preds_orig)
    
    print(f"\nBaseline KGE: {baseline_kge:.4f}")
    
    # Calculate importance for each feature
    feature_importance = {}
    
    # Dynamic features (weather variables)
    dynamic_features = data_module.dynamic_cols_no_target
    print("\nCalculating dynamic feature importance...")
    
    for idx, feat_name in enumerate(tqdm(dynamic_features)):
        # Perturb this feature
        perturbed_dynamic = all_dynamic_seqs.clone()
        
        # Shuffle values across samples to break relationship
        perm_indices = torch.randperm(n_samples)
        perturbed_dynamic[:, :, idx] = perturbed_dynamic[perm_indices, :, idx]
        
        # Get predictions with perturbed feature
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_samples, batch_size):
                end_idx = min(i + batch_size, n_samples)
                dynamic_batch = perturbed_dynamic[i:end_idx].to(device)
                static_batch = all_static_feats[i:end_idx].to(device)
                pred = model(dynamic_batch, static_batch)
                perturbed_preds.append(pred.cpu())
            perturbed_preds = torch.cat(perturbed_preds)
        
        # Calculate KGE with perturbation
        perturbed_preds_np = perturbed_preds.numpy().flatten()
        perturbed_preds_orig = unnormalize_predictions(perturbed_preds_np, data_module.norm_params)
        perturbed_kge = calculate_kge(targets_orig, perturbed_preds_orig)
        
        # Importance is the drop in performance
        importance = baseline_kge - perturbed_kge
        feature_importance[feat_name] = importance
    
    # Static features - Try multiple perturbation methods
    static_features = data_module.static_cols
    print("\nCalculating static feature importance with multiple methods...")
    
    # Method 1: Shuffle (original method)
    print("  Method 1: Shuffling across samples")
    static_importance_shuffle = {}
    
    for idx, feat_name in enumerate(tqdm(static_features[:30])):
        # Check if feature has any variance
        feat_values = all_static_feats[:, idx]
        if feat_values.std() < 1e-6:
            print(f"    Warning: {feat_name} has near-zero variance (std={feat_values.std():.6f})")
            static_importance_shuffle[feat_name] = 0.0
            continue
        
        # Perturb this feature
        perturbed_static = all_static_feats.clone()
        
        # Shuffle values across samples
        perm_indices = torch.randperm(n_samples)
        perturbed_static[:, idx] = perturbed_static[perm_indices, idx]
        
        # Get predictions with perturbed feature
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_samples, batch_size):
                end_idx = min(i + batch_size, n_samples)
                dynamic_batch = all_dynamic_seqs[i:end_idx].to(device)
                static_batch = perturbed_static[i:end_idx].to(device)
                pred = model(dynamic_batch, static_batch)
                perturbed_preds.append(pred.cpu())
            perturbed_preds = torch.cat(perturbed_preds)
        
        # Calculate KGE with perturbation
        perturbed_preds_np = perturbed_preds.numpy().flatten()
        perturbed_preds_orig = unnormalize_predictions(perturbed_preds_np, data_module.norm_params)
        perturbed_kge = calculate_kge(targets_orig, perturbed_preds_orig)
        
        importance = baseline_kge - perturbed_kge
        static_importance_shuffle[feat_name] = importance
    
    # Method 2: Set to mean value
    print("\n  Method 2: Setting to mean value")
    static_importance_mean = {}
    
    for idx, feat_name in enumerate(tqdm(static_features[:30])):
        # Set feature to its mean value
        perturbed_static = all_static_feats.clone()
        mean_val = all_static_feats[:, idx].mean()
        perturbed_static[:, idx] = mean_val
        
        # Get predictions
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_samples, batch_size):
                end_idx = min(i + batch_size, n_samples)
                dynamic_batch = all_dynamic_seqs[i:end_idx].to(device)
                static_batch = perturbed_static[i:end_idx].to(device)
                pred = model(dynamic_batch, static_batch)
                perturbed_preds.append(pred.cpu())
            perturbed_preds = torch.cat(perturbed_preds)
        
        perturbed_preds_np = perturbed_preds.numpy().flatten()
        perturbed_preds_orig = unnormalize_predictions(perturbed_preds_np, data_module.norm_params)
        perturbed_kge = calculate_kge(targets_orig, perturbed_preds_orig)
        
        importance = baseline_kge - perturbed_kge
        static_importance_mean[feat_name] = importance
    
    # Method 3: Add noise
    print("\n  Method 3: Adding Gaussian noise")
    static_importance_noise = {}
    
    for idx, feat_name in enumerate(tqdm(static_features[:30])):
        # Add Gaussian noise
        perturbed_static = all_static_feats.clone()
        feat_std = all_static_feats[:, idx].std()
        if feat_std > 1e-6:
            noise = torch.randn(n_samples) * feat_std
            perturbed_static[:, idx] = all_static_feats[:, idx] + noise
        else:
            # If no variance, add small noise
            noise = torch.randn(n_samples) * 0.1
            perturbed_static[:, idx] = all_static_feats[:, idx] + noise
        
        # Get predictions
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_samples, batch_size):
                end_idx = min(i + batch_size, n_samples)
                dynamic_batch = all_dynamic_seqs[i:end_idx].to(device)
                static_batch = perturbed_static[i:end_idx].to(device)
                pred = model(dynamic_batch, static_batch)
                perturbed_preds.append(pred.cpu())
            perturbed_preds = torch.cat(perturbed_preds)
        
        perturbed_preds_np = perturbed_preds.numpy().flatten()
        perturbed_preds_orig = unnormalize_predictions(perturbed_preds_np, data_module.norm_params)
        perturbed_kge = calculate_kge(targets_orig, perturbed_preds_orig)
        
        importance = baseline_kge - perturbed_kge
        static_importance_noise[feat_name] = importance
    
    # Combine results - use maximum importance from all methods
    print("\n  Combining static feature importance from all methods...")
    for feat_name in static_features[:30]:
        importance_values = [
            static_importance_shuffle.get(feat_name, 0),
            static_importance_mean.get(feat_name, 0),
            static_importance_noise.get(feat_name, 0)
        ]
        # Use maximum absolute importance
        max_importance = max(importance_values, key=abs)
        feature_importance[feat_name] = max_importance
        
        if debug and abs(max_importance) > 1e-6:
            print(f"    {feat_name}: shuffle={importance_values[0]:.4f}, "
                  f"mean={importance_values[1]:.4f}, noise={importance_values[2]:.4f}")
    
    # Debug: Check if model uses static features at all
    if debug:
        print("\n  Testing if model responds to static features...")
        # Set all static features to zero
        zero_static = torch.zeros_like(all_static_feats)
        
        with torch.no_grad():
            zero_preds = []
            batch_size = 64
            for i in range(0, min(200, n_samples), batch_size):  # Test on subset
                end_idx = min(i + batch_size, n_samples)
                dynamic_batch = all_dynamic_seqs[i:end_idx].to(device)
                static_batch = zero_static[i:end_idx].to(device)
                pred = model(dynamic_batch, static_batch)
                zero_preds.append(pred.cpu())
            zero_preds = torch.cat(zero_preds)
        
        # Compare with baseline
        baseline_subset = baseline_preds[:len(zero_preds)]
        diff = (zero_preds - baseline_subset).abs().mean()
        print(f"    Mean absolute difference when static features = 0: {diff:.6f}")
        
        if diff < 1e-4:
            print("    WARNING: Model appears to be ignoring static features!")
    
    return feature_importance, baseline_kge


def diagnose_static_features(model, data_module, watershed_df):
    """Diagnose why static features might have zero importance"""
    
    print("="*60)
    print("STATIC FEATURE DIAGNOSTIC")
    print("="*60)
    
    # 1. Check variance in raw static features
    print("\n1. Checking variance in raw (unnormalized) static features:")
    
    # Load raw static features for all watersheds
    static_features_raw = []
    for _, row in watershed_df.iterrows():
        attrs = get_watershed_attributes(
            data_module.bucket_name,
            data_module.base_attr_dir,
            row['subdirectory_name'],
            row['watershedID']
        )
        if not attrs.empty:
            static_features_raw.append(attrs)
    
    if static_features_raw:
        static_df = pd.DataFrame(static_features_raw)
        
        # Show statistics for first 10 static features
        for col in list(static_df.columns)[:10]:
            if col in data_module.static_cols:
                values = static_df[col].dropna()
                if len(values) > 0:
                    print(f"\n  {col}:")
                    print(f"    Count: {len(values)}, Missing: {static_df[col].isna().sum()}")
                    print(f"    Mean: {values.mean():.4f}, Std: {values.std():.4f}")
                    print(f"    Min: {values.min():.4f}, Max: {values.max():.4f}")
                    print(f"    CV (std/mean): {values.std()/values.mean():.4f}" if values.mean() != 0 else "    CV: N/A")
    
    # 2. Check if static features are connected in the model
    print("\n2. Checking model's use of static features:")
    
    # Get a sample batch
    val_loader = data_module.val_dataloader()
    for batch in val_loader:
        dynamic_seq, static_feat, target = batch
        break
    
    device = next(model.parameters()).device
    dynamic_seq = dynamic_seq.to(device)
    static_feat = static_feat.to(device)
    
    # Test 1: Normal prediction
    with torch.no_grad():
        normal_pred = model(dynamic_seq, static_feat)
    
    # Test 2: Zero static features
    zero_static = torch.zeros_like(static_feat)
    with torch.no_grad():
        zero_pred = model(dynamic_seq, zero_static)
    
    # Test 3: Random static features
    random_static = torch.randn_like(static_feat)
    with torch.no_grad():
        random_pred = model(dynamic_seq, random_static)
    
    print(f"\n  Predictions with different static inputs:")
    print(f"    Normal: mean={normal_pred.mean():.4f}, std={normal_pred.std():.4f}")
    print(f"    Zero static: mean={zero_pred.mean():.4f}, std={zero_pred.std():.4f}")
    print(f"    Random static: mean={random_pred.mean():.4f}, std={random_pred.std():.4f}")
    print(f"    Diff (normal-zero): {(normal_pred - zero_pred).abs().mean():.6f}")
    print(f"    Diff (normal-random): {(normal_pred - random_pred).abs().mean():.6f}")
    
    # 3. Check model weights related to static features
    print("\n3. Checking model weights for static features:")
    
    # Check static encoder weights
    static_encoder_weights = []
    for name, param in model.model.static_encoder.named_parameters():
        if 'weight' in name:
            weight_stats = {
                'layer': name,
                'mean': param.abs().mean().item(),
                'std': param.std().item(),
                'max': param.abs().max().item()
            }
            static_encoder_weights.append(weight_stats)
            print(f"    {name}: mean_abs={weight_stats['mean']:.6f}, "
                  f"std={weight_stats['std']:.6f}, max_abs={weight_stats['max']:.6f}")
    
    # Check output layer weights for static features
    output_weight = model.model.output_layer[0].weight  # First linear layer
    hidden_size = model.model.hidden_size
    static_size = len(data_module.static_cols)
    
    # The output layer concatenates [lstm_output, static_features]
    # So weights[:, :hidden_size] are for LSTM output
    # and weights[:, hidden_size:] are for static features
    static_weights = output_weight[:, hidden_size:hidden_size+static_size]
    lstm_weights = output_weight[:, :hidden_size]
    
    print(f"\n  Output layer weight magnitudes:")
    print(f"    LSTM weights: mean_abs={lstm_weights.abs().mean():.6f}, "
          f"std={lstm_weights.std():.6f}")
    print(f"    Static weights: mean_abs={static_weights.abs().mean():.6f}, "
          f"std={static_weights.std():.6f}")
    print(f"    Ratio (static/lstm): {static_weights.abs().mean() / lstm_weights.abs().mean():.4f}")
    
    return static_df if 'static_df' in locals() else None



def diagnose_static_features(model, data_module, watershed_df):
    """Diagnose why static features might have zero importance"""
    
    print("="*60)
    print("STATIC FEATURE DIAGNOSTIC")
    print("="*60)
    
    # 1. Check variance in raw static features
    print("\n1. Checking variance in raw (unnormalized) static features:")
    
    # Load raw static features for all watersheds
    static_features_raw = []
    for _, row in watershed_df.iterrows():
        attrs = get_watershed_attributes(
            data_module.bucket_name,
            data_module.base_attr_dir,
            row['subdirectory_name'],
            row['watershedID']
        )
        if not attrs.empty:
            static_features_raw.append(attrs)
    
    if static_features_raw:
        static_df = pd.DataFrame(static_features_raw)
        
        # Show statistics for first 10 static features
        for col in list(static_df.columns)[:10]:
            if col in data_module.static_cols:
                values = static_df[col].dropna()
                if len(values) > 0:
                    print(f"\n  {col}:")
                    print(f"    Count: {len(values)}, Missing: {static_df[col].isna().sum()}")
                    print(f"    Mean: {values.mean():.4f}, Std: {values.std():.4f}")
                    print(f"    Min: {values.min():.4f}, Max: {values.max():.4f}")
                    print(f"    CV (std/mean): {values.std()/values.mean():.4f}" if values.mean() != 0 else "    CV: N/A")
    
    # 2. Check if static features are connected in the model
    print("\n2. Checking model's use of static features:")
    
    # Get a sample batch
    val_loader = data_module.val_dataloader()
    for batch in val_loader:
        dynamic_seq, static_feat, target = batch
        break
    
    device = next(model.parameters()).device
    dynamic_seq = dynamic_seq.to(device)
    static_feat = static_feat.to(device)
    
    # Test 1: Normal prediction
    with torch.no_grad():
        normal_pred = model(dynamic_seq, static_feat)
    
    # Test 2: Zero static features
    zero_static = torch.zeros_like(static_feat)
    with torch.no_grad():
        zero_pred = model(dynamic_seq, zero_static)
    
    # Test 3: Random static features
    random_static = torch.randn_like(static_feat)
    with torch.no_grad():
        random_pred = model(dynamic_seq, random_static)
    
    print(f"\n  Predictions with different static inputs:")
    print(f"    Normal: mean={normal_pred.mean():.4f}, std={normal_pred.std():.4f}")
    print(f"    Zero static: mean={zero_pred.mean():.4f}, std={zero_pred.std():.4f}")
    print(f"    Random static: mean={random_pred.mean():.4f}, std={random_pred.std():.4f}")
    print(f"    Diff (normal-zero): {(normal_pred - zero_pred).abs().mean():.6f}")
    print(f"    Diff (normal-random): {(normal_pred - random_pred).abs().mean():.6f}")
    
    # 3. Check model weights related to static features
    print("\n3. Checking model weights for static features:")
    
    # Check static encoder weights
    static_encoder_weights = []
    for name, param in model.model.static_encoder.named_parameters():
        if 'weight' in name:
            weight_stats = {
                'layer': name,
                'mean': param.abs().mean().item(),
                'std': param.std().item(),
                'max': param.abs().max().item()
            }
            static_encoder_weights.append(weight_stats)
            print(f"    {name}: mean_abs={weight_stats['mean']:.6f}, "
                  f"std={weight_stats['std']:.6f}, max_abs={weight_stats['max']:.6f}")
    
    # Check output layer weights for static features
    output_weight = model.model.output_layer[0].weight  # First linear layer
    hidden_size = model.model.hidden_size
    static_size = len(data_module.static_cols)
    
    # The output layer concatenates [lstm_output, static_features]
    # So weights[:, :hidden_size] are for LSTM output
    # and weights[:, hidden_size:] are for static features
    static_weights = output_weight[:, hidden_size:hidden_size+static_size]
    lstm_weights = output_weight[:, :hidden_size]
    
    print(f"\n  Output layer weight magnitudes:")
    print(f"    LSTM weights: mean_abs={lstm_weights.abs().mean():.6f}, "
          f"std={lstm_weights.std():.6f}")
    print(f"    Static weights: mean_abs={static_weights.abs().mean():.6f}, "
          f"std={static_weights.std():.6f}")
    print(f"    Ratio (static/lstm): {static_weights.abs().mean() / lstm_weights.abs().mean():.4f}")
    
    return static_df if 'static_df' in locals() else None

def plot_feature_importance_bars(feature_importance, data_module, save_dir='./feature_importance_plots'):
    """Create barplots for feature importance - both raw and normalized versions"""
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Define weather features (including aridity)
    weather_features = [
#        'potential_evaporation_sum_ERA5_LAND', 
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 
        'temperature_2m_mean', 
        'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean',
        'Precip_smoothed_5day',
        'Precip_lagged_90day'
    ]
    
    static_features = data_module.static_cols
    
    # Filter importance scores
    weather_importance = {k: v for k, v in feature_importance.items() if k in weather_features}
    static_importance = {k: v for k, v in feature_importance.items() if k in static_features}
    
    # Function to normalize importance scores
    def normalize_importance(importance_dict, method='minmax'):
        if not importance_dict:
            return {}
        
        values = np.array(list(importance_dict.values()))
        
        if method == 'minmax':
            # Min-max normalization to [0, 1]
            min_val = values.min()
            max_val = values.max()
            if max_val > min_val:
                normalized = (values - min_val) / (max_val - min_val)
            else:
                normalized = np.ones_like(values)
        elif method == 'zscore':
            # Z-score normalization
            mean_val = values.mean()
            std_val = values.std()
            if std_val > 0:
                normalized = (values - mean_val) / std_val
            else:
                normalized = np.zeros_like(values)
        elif method == 'relative':
            # Relative to sum (percentage)
            total = values.sum()
            if total > 0:
                normalized = values / total * 100
            else:
                normalized = np.zeros_like(values)
        
        return {k: norm_val for k, norm_val in zip(importance_dict.keys(), normalized)}
    
    # Create figure with RAW importance values (original)
    fig, axes = plt.subplots(3, 1, figsize=(12, 16))
    fig.suptitle('Feature Importance - Raw KGE Drop', fontsize=16, fontweight='bold')
    
    # 1. Weather variables importance (RAW)
    ax = axes[0]
    if weather_importance:
        weather_df = pd.DataFrame(list(weather_importance.items()), columns=['Feature', 'Importance'])
        weather_df = weather_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(weather_df['Feature'], weather_df['Importance'], color='skyblue', edgecolor='navy')
        ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
        ax.set_title('Weather Variables - Raw Values', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                    ha='left', va='center', fontsize=9)
    
    # 2. Static watershed features importance (RAW)
    ax = axes[1]
    if static_importance:
        static_df = pd.DataFrame(list(static_importance.items()), columns=['Feature', 'Importance'])
        static_df = static_df.sort_values('Importance', ascending=False).head(20)
        static_df = static_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(static_df['Feature'], static_df['Importance'], color='lightcoral', edgecolor='darkred')
        ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
        ax.set_title('Top 20 Static Watershed Features - Raw Values', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                    ha='left', va='center', fontsize=9)
    
    # 3. All features combined (RAW)
    ax = axes[2]
    all_df = pd.DataFrame(list(feature_importance.items()), columns=['Feature', 'Importance'])
    all_df = all_df.sort_values('Importance', ascending=False).head(25)
    all_df = all_df.sort_values('Importance', ascending=True)
    
    # Color code by type
    colors = ['skyblue' if feat in weather_features else 'lightcoral' for feat in all_df['Feature']]
    
    bars = ax.barh(all_df['Feature'], all_df['Importance'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
    ax.set_title('Top 25 All Features - Raw Values', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                ha='left', va='center', fontsize=9)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='skyblue', label='Weather Variables'),
                      Patch(facecolor='lightcoral', label='Static Features')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_barplots_raw.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure with NORMALIZED importance values
    fig, axes = plt.subplots(3, 1, figsize=(12, 16))
    fig.suptitle('Feature Importance - Normalized Values', fontsize=16, fontweight='bold')
    
    # Normalize within each category
    weather_importance_norm = normalize_importance(weather_importance, method='minmax')
    static_importance_norm = normalize_importance(static_importance, method='minmax')
    all_importance_norm = normalize_importance(feature_importance, method='minmax')
    
    # 1. Weather variables importance (NORMALIZED)
    ax = axes[0]
    if weather_importance_norm:
        weather_df = pd.DataFrame(list(weather_importance_norm.items()), columns=['Feature', 'Importance'])
        weather_df = weather_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(weather_df['Feature'], weather_df['Importance'], color='skyblue', edgecolor='navy')
        ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax.set_title('Weather Variables - Normalized within Category', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, 1.1)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                    ha='left', va='center', fontsize=9)
    
    # 2. Static watershed features importance (NORMALIZED)
    ax = axes[1]
    if static_importance_norm:
        static_df = pd.DataFrame(list(static_importance_norm.items()), columns=['Feature', 'Importance'])
        static_df = static_df.sort_values('Importance', ascending=False).head(20)
        static_df = static_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(static_df['Feature'], static_df['Importance'], color='lightcoral', edgecolor='darkred')
        ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax.set_title('Top 20 Static Watershed Features - Normalized within Category', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, 1.1)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                    ha='left', va='center', fontsize=9)
    
    # 3. All features combined (NORMALIZED)
    ax = axes[2]
    all_df = pd.DataFrame(list(all_importance_norm.items()), columns=['Feature', 'Importance'])
    all_df = all_df.sort_values('Importance', ascending=False).head(25)
    all_df = all_df.sort_values('Importance', ascending=True)
    
    # Color code by type
    colors = ['skyblue' if feat in weather_features else 'lightcoral' for feat in all_df['Feature']]
    
    bars = ax.barh(all_df['Feature'], all_df['Importance'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
    ax.set_title('Top 25 All Features - Normalized across All Features', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, 1.1)
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                ha='left', va='center', fontsize=9)
    
    # Add legend
    legend_elements = [Patch(facecolor='skyblue', label='Weather Variables'),
                      Patch(facecolor='lightcoral', label='Static Features')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_barplots_normalized.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create comparison plot showing both raw and normalized for static features
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))
    fig.suptitle('Static Feature Importance Comparison: Raw vs Normalized', fontsize=16, fontweight='bold')
    
    if static_importance:
        # Get top 15 static features by raw importance
        static_df = pd.DataFrame(list(static_importance.items()), columns=['Feature', 'RawImportance'])
        static_df = static_df.sort_values('RawImportance', ascending=False).head(15)
        
        # Add normalized values
        static_df['NormImportance'] = static_df['Feature'].map(static_importance_norm)
        
        # Sort by normalized importance for plotting
        static_df = static_df.sort_values('NormImportance', ascending=True)
        
        # Raw values plot
        bars1 = ax1.barh(static_df['Feature'], static_df['RawImportance'], 
                         color='lightcoral', edgecolor='darkred')
        ax1.set_xlabel('Raw Feature Importance (KGE Drop)', fontsize=12)
        ax1.set_title('Raw Values', fontsize=14)
        ax1.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar, val in zip(bars1, static_df['RawImportance']):
            ax1.text(bar.get_width() + 0.0001, bar.get_y() + bar.get_height()/2, 
                    f'{val:.5f}', ha='left', va='center', fontsize=9)
        
        # Normalized values plot
        bars2 = ax2.barh(static_df['Feature'], static_df['NormImportance'], 
                         color='darkseagreen', edgecolor='darkgreen')
        ax2.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax2.set_title('Normalized Values', fontsize=14)
        ax2.grid(axis='x', alpha=0.3)
        ax2.set_xlim(0, 1.1)
        
        # Add value labels
        for bar, val in zip(bars2, static_df['NormImportance']):
            ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                    f'{val:.3f}', ha='left', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_static_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save feature importance to CSV for further analysis
    all_features_df = pd.DataFrame([
        {
            'Feature': feat, 
            'RawImportance': imp,
            'NormalizedImportance': all_importance_norm.get(feat, 0),
            'Type': 'Weather' if feat in weather_features else 'Static'
        }
        for feat, imp in feature_importance.items()
    ])
    all_features_df = all_features_df.sort_values('RawImportance', ascending=False)
    all_features_df.to_csv(f'{save_dir}/feature_importance_values.csv', index=False)
    
    print(f"Saved feature importance barplots to {save_dir}")
    print(f"  - Raw values: feature_importance_barplots_raw.png")
    print(f"  - Normalized values: feature_importance_barplots_normalized.png")
    print(f"  - Static comparison: feature_importance_static_comparison.png")
    print(f"  - CSV data: feature_importance_values.csv")

def calculate_temporal_feature_importance(model, data_module, n_samples=500, lag_days=None):
    """Calculate feature importance at different time lags"""
    
    from train import calculate_kge, unnormalize_predictions
    
    device = next(model.parameters()).device
    
    if lag_days is None:
        # Analyze every 10 days for computational efficiency
        lag_days = list(range(0, 365, 10)) + [364]  # Include last day
    
    # Get validation dataloader
    val_loader = data_module.val_dataloader()
    
    # Collect data
    all_dynamic_seqs = []
    all_static_feats = []
    all_targets = []
    
    print("Collecting validation data...")
    sample_count = 0
    
    with torch.no_grad():
        for batch in val_loader:
            if sample_count >= n_samples:
                break
                
            dynamic_seq, static_feat, target = batch
            
            all_dynamic_seqs.append(dynamic_seq)
            all_static_feats.append(static_feat)
            all_targets.append(target)
            
            sample_count += len(target)
    
    # Concatenate
    all_dynamic_seqs = torch.cat(all_dynamic_seqs)[:n_samples]
    all_static_feats = torch.cat(all_static_feats)[:n_samples]
    all_targets = torch.cat(all_targets)[:n_samples]
    
    # Get baseline predictions
    with torch.no_grad():
        baseline_preds = []
        batch_size = 64
        for i in range(0, n_samples, batch_size):
            end_idx = min(i + batch_size, n_samples)
            dynamic_batch = all_dynamic_seqs[i:end_idx].to(device)
            static_batch = all_static_feats[i:end_idx].to(device)
            pred = model(dynamic_batch, static_batch)
            baseline_preds.append(pred.cpu())
        baseline_preds = torch.cat(baseline_preds)
    
    # Calculate baseline KGE
    baseline_preds_np = baseline_preds.numpy().flatten()
    targets_np = all_targets.numpy().flatten()
    baseline_preds_orig = unnormalize_predictions(baseline_preds_np, data_module.norm_params)
    targets_orig = unnormalize_predictions(targets_np, data_module.norm_params)
    baseline_kge = calculate_kge(targets_orig, baseline_preds_orig)
    
    # Weather features (excluding aridity for temporal analysis)
    weather_features = [
 #       'potential_evaporation_sum_ERA5_LAND', 
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 
        'temperature_2m_mean', 
        'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean',
        'Precip_smoothed_5day',
        'Precip_lagged_90day'
    ]
    
    # Find indices of weather features
    dynamic_features = data_module.dynamic_cols_no_target
    weather_indices = [i for i, feat in enumerate(dynamic_features) if feat in weather_features]
    weather_names = [feat for feat in dynamic_features if feat in weather_features]
    
    # Calculate importance at each lag
    temporal_importance = {feat: [] for feat in weather_names}
    
    print(f"\nCalculating temporal feature importance for {len(lag_days)} lags...")
    
    for lag in tqdm(lag_days):
        time_idx = 364 - lag  # 364 is the last timestep (0-indexed for 365 days)
        
        for feat_idx, feat_name in zip(weather_indices, weather_names):
            # Perturb feature at specific time lag
            perturbed_dynamic = all_dynamic_seqs.clone()
            
            # Shuffle values at this specific timestep
            perm_indices = torch.randperm(n_samples)
            perturbed_dynamic[:, time_idx, feat_idx] = perturbed_dynamic[perm_indices, time_idx, feat_idx]
            
            # Get predictions
            with torch.no_grad():
                perturbed_preds = []
                batch_size = 64
                for i in range(0, n_samples, batch_size):
                    end_idx = min(i + batch_size, n_samples)
                    dynamic_batch = perturbed_dynamic[i:end_idx].to(device)
                    static_batch = all_static_feats[i:end_idx].to(device)
                    pred = model(dynamic_batch, static_batch)
                    perturbed_preds.append(pred.cpu())
                perturbed_preds = torch.cat(perturbed_preds)
            
            # Calculate importance
            perturbed_preds_np = perturbed_preds.numpy().flatten()
            perturbed_preds_orig = unnormalize_predictions(perturbed_preds_np, data_module.norm_params)
            perturbed_kge = calculate_kge(targets_orig, perturbed_preds_orig)
            
            importance = baseline_kge - perturbed_kge
            temporal_importance[feat_name].append(importance)
    
    return temporal_importance, lag_days, baseline_kge

def plot_feature_importance_bars(feature_importance, data_module, save_dir='./feature_importance_plots'):
    """Create barplots for feature importance - both raw and normalized versions"""
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Define weather features (including aridity)
    weather_features = [
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 
        'temperature_2m_mean', 
        'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean',
        'Precip_smoothed_5day',
        'Precip_lagged_90day'
    ]
    
    static_features = data_module.static_cols
    
    # Filter importance scores
    weather_importance = {k: v for k, v in feature_importance.items() if k in weather_features}
    static_importance = {k: v for k, v in feature_importance.items() if k in static_features}
    
    # Function to normalize importance scores
    def normalize_importance(importance_dict, method='minmax'):
        if not importance_dict:
            return {}
        
        values = np.array(list(importance_dict.values()))
        
        if method == 'minmax':
            # Min-max normalization to [0, 1]
            min_val = values.min()
            max_val = values.max()
            if max_val > min_val:
                normalized = (values - min_val) / (max_val - min_val)
            else:
                normalized = np.ones_like(values)
        elif method == 'zscore':
            # Z-score normalization
            mean_val = values.mean()
            std_val = values.std()
            if std_val > 0:
                normalized = (values - mean_val) / std_val
            else:
                normalized = np.zeros_like(values)
        elif method == 'relative':
            # Relative to sum (percentage)
            total = values.sum()
            if total > 0:
                normalized = values / total * 100
            else:
                normalized = np.zeros_like(values)
        
        return {k: norm_val for k, norm_val in zip(importance_dict.keys(), normalized)}
    
    # Create figure with RAW importance values (original)
    fig, axes = plt.subplots(3, 1, figsize=(12, 16))
    fig.suptitle('Feature Importance - Raw KGE Drop', fontsize=16, fontweight='bold')
    
    # 1. Weather variables importance (RAW)
    ax = axes[0]
    if weather_importance:
        weather_df = pd.DataFrame(list(weather_importance.items()), columns=['Feature', 'Importance'])
        weather_df = weather_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(weather_df['Feature'], weather_df['Importance'], color='skyblue', edgecolor='navy')
        ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
        ax.set_title('Weather Variables - Raw Values', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                    ha='left', va='center', fontsize=9)
    
    # 2. Static watershed features importance (RAW)
    ax = axes[1]
    if static_importance:
        static_df = pd.DataFrame(list(static_importance.items()), columns=['Feature', 'Importance'])
        static_df = static_df.sort_values('Importance', ascending=False).head(20)
        static_df = static_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(static_df['Feature'], static_df['Importance'], color='lightcoral', edgecolor='darkred')
        ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
        ax.set_title('Top 20 Static Watershed Features - Raw Values', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                    ha='left', va='center', fontsize=9)
    
    # 3. All features combined (RAW)
    ax = axes[2]
    all_df = pd.DataFrame(list(feature_importance.items()), columns=['Feature', 'Importance'])
    all_df = all_df.sort_values('Importance', ascending=False).head(25)
    all_df = all_df.sort_values('Importance', ascending=True)
    
    # Color code by type
    colors = ['skyblue' if feat in weather_features else 'lightcoral' for feat in all_df['Feature']]
    
    bars = ax.barh(all_df['Feature'], all_df['Importance'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Feature Importance (KGE Drop)', fontsize=12)
    ax.set_title('Top 25 All Features - Raw Values', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.0005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                ha='left', va='center', fontsize=9)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='skyblue', label='Weather Variables'),
                      Patch(facecolor='lightcoral', label='Static Features')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_barplots_raw.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure with NORMALIZED importance values
    fig, axes = plt.subplots(3, 1, figsize=(12, 16))
    fig.suptitle('Feature Importance - Normalized Values', fontsize=16, fontweight='bold')
    
    # Normalize within each category
    weather_importance_norm = normalize_importance(weather_importance, method='minmax')
    static_importance_norm = normalize_importance(static_importance, method='minmax')
    all_importance_norm = normalize_importance(feature_importance, method='minmax')
    
    # 1. Weather variables importance (NORMALIZED)
    ax = axes[0]
    if weather_importance_norm:
        weather_df = pd.DataFrame(list(weather_importance_norm.items()), columns=['Feature', 'Importance'])
        weather_df = weather_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(weather_df['Feature'], weather_df['Importance'], color='skyblue', edgecolor='navy')
        ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax.set_title('Weather Variables - Normalized within Category', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, 1.1)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                    ha='left', va='center', fontsize=9)
    
    # 2. Static watershed features importance (NORMALIZED)
    ax = axes[1]
    if static_importance_norm:
        static_df = pd.DataFrame(list(static_importance_norm.items()), columns=['Feature', 'Importance'])
        static_df = static_df.sort_values('Importance', ascending=False).head(20)
        static_df = static_df.sort_values('Importance', ascending=True)
        
        bars = ax.barh(static_df['Feature'], static_df['Importance'], color='lightcoral', edgecolor='darkred')
        ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax.set_title('Top 20 Static Watershed Features - Normalized within Category', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, 1.1)
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                    ha='left', va='center', fontsize=9)
    
    # 3. All features combined (NORMALIZED)
    ax = axes[2]
    all_df = pd.DataFrame(list(all_importance_norm.items()), columns=['Feature', 'Importance'])
    all_df = all_df.sort_values('Importance', ascending=False).head(25)
    all_df = all_df.sort_values('Importance', ascending=True)
    
    # Color code by type
    colors = ['skyblue' if feat in weather_features else 'lightcoral' for feat in all_df['Feature']]
    
    bars = ax.barh(all_df['Feature'], all_df['Importance'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
    ax.set_title('Top 25 All Features - Normalized across All Features', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, 1.1)
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
                ha='left', va='center', fontsize=9)
    
    # Add legend
    legend_elements = [Patch(facecolor='skyblue', label='Weather Variables'),
                      Patch(facecolor='lightcoral', label='Static Features')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_barplots_normalized.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create comparison plot showing both raw and normalized for static features
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))
    fig.suptitle('Static Feature Importance Comparison: Raw vs Normalized', fontsize=16, fontweight='bold')
    
    if static_importance:
        # Get top 15 static features by raw importance
        static_df = pd.DataFrame(list(static_importance.items()), columns=['Feature', 'RawImportance'])
        static_df = static_df.sort_values('RawImportance', ascending=False).head(15)
        
        # Add normalized values
        static_df['NormImportance'] = static_df['Feature'].map(static_importance_norm)
        
        # Sort by normalized importance for plotting
        static_df = static_df.sort_values('NormImportance', ascending=True)
        
        # Raw values plot
        bars1 = ax1.barh(static_df['Feature'], static_df['RawImportance'], 
                         color='lightcoral', edgecolor='darkred')
        ax1.set_xlabel('Raw Feature Importance (KGE Drop)', fontsize=12)
        ax1.set_title('Raw Values', fontsize=14)
        ax1.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar, val in zip(bars1, static_df['RawImportance']):
            ax1.text(bar.get_width() + 0.0001, bar.get_y() + bar.get_height()/2, 
                    f'{val:.5f}', ha='left', va='center', fontsize=9)
        
        # Normalized values plot
        bars2 = ax2.barh(static_df['Feature'], static_df['NormImportance'], 
                         color='darkseagreen', edgecolor='darkgreen')
        ax2.set_xlabel('Normalized Feature Importance (0-1 scale)', fontsize=12)
        ax2.set_title('Normalized Values', fontsize=14)
        ax2.grid(axis='x', alpha=0.3)
        ax2.set_xlim(0, 1.1)
        
        # Add value labels
        for bar, val in zip(bars2, static_df['NormImportance']):
            ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                    f'{val:.3f}', ha='left', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/feature_importance_static_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save feature importance to CSV for further analysis
    all_features_df = pd.DataFrame([
        {
            'Feature': feat, 
            'RawImportance': imp,
            'NormalizedImportance': all_importance_norm.get(feat, 0),
            'Type': 'Weather' if feat in weather_features else 'Static'
        }
        for feat, imp in feature_importance.items()
    ])
    all_features_df = all_features_df.sort_values('RawImportance', ascending=False)
    all_features_df.to_csv(f'{save_dir}/feature_importance_values.csv', index=False)
    
    print(f"Saved feature importance barplots to {save_dir}")
    print(f"  - Raw values: feature_importance_barplots_raw.png")
    print(f"  - Normalized values: feature_importance_barplots_normalized.png")
    print(f"  - Static comparison: feature_importance_static_comparison.png")
    print(f"  - CSV data: feature_importance_values.csv")


def plot_temporal_feature_importance(temporal_importance, lag_days, save_dir='./feature_importance_plots'):
    """Plot temporal feature importance in both raw and normalized scales"""
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Define a minimum value for log transformation
    # Use 1% of the maximum importance as the floor, or 1e-4 if all values are very small
    all_values = []
    for values in temporal_importance.values():
        all_values.extend(values)
    max_importance = max(all_values) if all_values else 1.0
    log_floor = max(max_importance * 0.01, 1e-4)
    
    print(f"Using log floor value: {log_floor:.6f} (1% of max importance: {max_importance:.4f})")
    
    # Create figure with RAW values (original)
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    fig.suptitle('Temporal Feature Importance - Raw KGE Drop', fontsize=16, fontweight='bold')
    
    # Color palette
    colors = plt.cm.tab10(np.linspace(0, 1, len(temporal_importance)))
    
    # Linear scale (RAW)
    ax = axes[0]
    for idx, (feature, importance_values) in enumerate(temporal_importance.items()):
        # Reverse lag days for intuitive display (1 day ago to 365 days ago)
        display_lags = [365 - lag for lag in lag_days]
        ax.plot(display_lags, importance_values, 'o-', label=feature.replace('_', ' ').title(), 
                color=colors[idx], linewidth=2, markersize=6, markeredgecolor='white')
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Feature Importance (KGE Drop)', fontsize=12)
    ax.set_title('Linear Scale - Raw Values', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 365)
    
    # Add vertical lines for key periods
    ax.axvline(x=7, color='gray', linestyle='--', alpha=0.5, label='1 week')
    ax.axvline(x=30, color='gray', linestyle='--', alpha=0.5, label='1 month')
    ax.axvline(x=90, color='gray', linestyle='--', alpha=0.5, label='3 months')
    
    # Log scale (RAW) - with improved handling of zeros
    ax = axes[1]
    for idx, (feature, importance_values) in enumerate(temporal_importance.items()):
        display_lags = [365 - lag for lag in lag_days]
        # Add log floor to avoid log(0) and make small values visible
        importance_log = np.array(importance_values) + log_floor
        
        # Mark points that were originally zero
        zero_mask = np.array(importance_values) == 0
        
        # Plot the line
        ax.semilogy(display_lags, importance_log, 'o-', label=feature.replace('_', ' ').title(), 
                    color=colors[idx], linewidth=2, markersize=6, markeredgecolor='white')
        
        # Highlight zero values with different marker
        if np.any(zero_mask):
            zero_lags = [dl for dl, zm in zip(display_lags, zero_mask) if zm]
            zero_vals = [log_floor] * len(zero_lags)
            ax.semilogy(zero_lags, zero_vals, 'x', color=colors[idx], markersize=8, 
                       markeredgecolor='black', markeredgewidth=1)
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel(f'Feature Importance (KGE Drop) + {log_floor:.4f} - Log Scale', fontsize=12)
    ax.set_title('Log Scale - Raw Values (x marks indicate originally zero values)', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0, 365)
    
    # Add horizontal line at log floor
    ax.axhline(y=log_floor, color='red', linestyle=':', alpha=0.5, label=f'Log floor ({log_floor:.4f})')
    
    # Add vertical lines
    ax.axvline(x=7, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=30, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=90, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/temporal_feature_importance_raw.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure with NORMALIZED values
    fig, axes = plt.subplots(3, 1, figsize=(14, 16))
    fig.suptitle('Temporal Feature Importance - Normalized Values', fontsize=16, fontweight='bold')
    
    # Normalize within each feature (across time lags)
    normalized_importance_within = {}
    for feature, values in temporal_importance.items():
        values_array = np.array(values)
        min_val = values_array.min()
        max_val = values_array.max()
        if max_val > min_val:
            normalized = (values_array - min_val) / (max_val - min_val)
        else:
            normalized = np.ones_like(values_array) * 0.5  # Set to 0.5 if all values are the same
        normalized_importance_within[feature] = normalized.tolist()
    
    # Normalize across features at each lag
    normalized_importance_across = {feat: [] for feat in temporal_importance.keys()}
    for i, lag in enumerate(lag_days):
        values_at_lag = [temporal_importance[feat][i] for feat in temporal_importance.keys()]
        values_array = np.array(values_at_lag)
        min_val = values_array.min()
        max_val = values_array.max()
        if max_val > min_val:
            normalized = (values_array - min_val) / (max_val - min_val)
        else:
            normalized = np.ones_like(values_array) * (1.0 / len(values_array))  # Equal importance
        
        for j, feat in enumerate(temporal_importance.keys()):
            normalized_importance_across[feat].append(normalized[j])
    
    # Plot 1: Normalized within each feature
    ax = axes[0]
    for idx, (feature, importance_values) in enumerate(normalized_importance_within.items()):
        display_lags = [365 - lag for lag in lag_days]
        ax.plot(display_lags, importance_values, 'o-', label=feature.replace('_', ' ').title(), 
                color=colors[idx], linewidth=2, markersize=6, markeredgecolor='white')
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Normalized Importance (0-1)', fontsize=12)
    ax.set_title('Normalized Within Each Feature (Shows temporal pattern per feature)', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 365)
    ax.set_ylim(-0.05, 1.05)
    
    # Add vertical lines
    ax.axvline(x=7, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=30, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=90, color='gray', linestyle='--', alpha=0.5)
    
    # Plot 2: Normalized across features at each lag
    ax = axes[1]
    for idx, (feature, importance_values) in enumerate(normalized_importance_across.items()):
        display_lags = [365 - lag for lag in lag_days]
        ax.plot(display_lags, importance_values, 'o-', label=feature.replace('_', ' ').title(), 
                color=colors[idx], linewidth=2, markersize=6, markeredgecolor='white')
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Normalized Importance (0-1)', fontsize=12)
    ax.set_title('Normalized Across Features at Each Lag (Shows relative importance between features)', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 365)
    ax.set_ylim(-0.05, 1.05)
    
    # Add vertical lines
    ax.axvline(x=7, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=30, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=90, color='gray', linestyle='--', alpha=0.5)
    
    # Plot 3: Heatmap of normalized importance
    ax = axes[2]
    
    # Create matrix for heatmap
    feature_names = list(temporal_importance.keys())
    importance_matrix = np.array([temporal_importance[feat] for feat in feature_names])
    
    # Add log floor before normalization for better visualization
    importance_matrix_log = importance_matrix + log_floor
    
    # Normalize each row (feature) to [0, 1]
    importance_matrix_norm = np.zeros_like(importance_matrix)
    for i in range(importance_matrix.shape[0]):
        row = importance_matrix[i]
        row_min, row_max = row.min(), row.max()
        if row_max > row_min:
            importance_matrix_norm[i] = (row - row_min) / (row_max - row_min)
        else:
            importance_matrix_norm[i] = np.ones_like(row) * 0.5
    
    # Create heatmap
    display_lags = [365 - lag for lag in lag_days]
    im = ax.imshow(importance_matrix_norm, aspect='auto', cmap='YlOrRd', interpolation='nearest')
    
    # Set ticks and labels
    ax.set_yticks(range(len(feature_names)))
    ax.set_yticklabels([feat.replace('_', ' ').title() for feat in feature_names])
    
    # Set x-axis ticks at specific intervals
    tick_indices = [i for i, lag in enumerate(display_lags) if lag in [1, 7, 30, 60, 90, 180, 365]]
    tick_labels = [display_lags[i] for i in tick_indices]
    ax.set_xticks(tick_indices)
    ax.set_xticklabels(tick_labels)
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Weather Variable', fontsize=12)
    ax.set_title('Temporal Importance Heatmap (Normalized by Feature)', fontsize=14)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Normalized Importance', rotation=270, labelpad=20)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/temporal_feature_importance_normalized.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create focused plots for recent time periods with log scale
    fig, axes = plt.subplots(3, 1, figsize=(14, 14))
    fig.suptitle('Temporal Feature Importance - Focus on Recent Time Periods', fontsize=16, fontweight='bold')
    
    # Filter for recent 90 days
    recent_mask = np.array([365 - lag <= 90 for lag in lag_days])
    recent_lags = [lag for lag, mask in zip(lag_days, recent_mask) if mask]
    
    # Plot 1: Raw values for recent 90 days
    ax = axes[0]
    for idx, (feature, importance_values) in enumerate(temporal_importance.items()):
        recent_values = [val for val, mask in zip(importance_values, recent_mask) if mask]
        recent_display_lags = [365 - lag for lag in recent_lags]
        ax.plot(recent_display_lags, recent_values, 'o-', label=feature.replace('_', ' ').title(), 
                color=colors[idx], linewidth=2.5, markersize=8, markeredgecolor='white')
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Feature Importance (KGE Drop)', fontsize=12)
    ax.set_title('Last 90 Days - Raw Values', fontsize=14)
    ax.legend(frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 90)
    
    # Add vertical lines for key periods
    ax.axvline(x=1, color='red', linestyle=':', alpha=0.5, label='1 day')
    ax.axvline(x=7, color='orange', linestyle='--', alpha=0.5, label='1 week')
    ax.axvline(x=30, color='green', linestyle='--', alpha=0.5, label='1 month')
    
    # Plot 2: Log scale for recent 90 days
    ax = axes[1]
    for idx, (feature, importance_values) in enumerate(temporal_importance.items()):
        recent_values = [val for val, mask in zip(importance_values, recent_mask) if mask]
        recent_display_lags = [365 - lag for lag in recent_lags]
        
        # Add log floor
        recent_values_log = np.array(recent_values) + log_floor
        
        # Mark zeros
        zero_mask_recent = np.array(recent_values) == 0
        
        ax.semilogy(recent_display_lags, recent_values_log, 'o-', label=feature.replace('_', ' ').title(), 
                    color=colors[idx], linewidth=2.5, markersize=8, markeredgecolor='white')
        
        # Mark zero values
        if np.any(zero_mask_recent):
            zero_lags = [dl for dl, zm in zip(recent_display_lags, zero_mask_recent) if zm]
            zero_vals = [log_floor] * len(zero_lags)
            ax.semilogy(zero_lags, zero_vals, 'x', color=colors[idx], markersize=10, 
                       markeredgecolor='black', markeredgewidth=1.5)
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel(f'Feature Importance + {log_floor:.4f} (Log Scale)', fontsize=12)
    ax.set_title('Last 90 Days - Log Scale (x marks indicate zero importance)', fontsize=14)
    ax.legend(frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0, 90)
    
    # Add reference line at log floor
    ax.axhline(y=log_floor, color='red', linestyle=':', alpha=0.5)
    
    # Add vertical lines
    ax.axvline(x=1, color='red', linestyle=':', alpha=0.5)
    ax.axvline(x=7, color='orange', linestyle='--', alpha=0.5)
    ax.axvline(x=30, color='green', linestyle='--', alpha=0.5)
    
    # Plot 3: Normalized values for recent 90 days
    ax = axes[2]
    for idx, (feature, importance_values) in enumerate(normalized_importance_within.items()):
        recent_values = [val for val, mask in zip(importance_values, recent_mask) if mask]
        recent_display_lags = [365 - lag for lag in recent_lags]
        ax.plot(recent_display_lags, recent_values, 'o-', label=feature.replace('_', ' ').title(), 
                color=colors[idx], linewidth=2.5, markersize=8, markeredgecolor='white')
    
    ax.set_xlabel('Days Before Prediction', fontsize=12)
    ax.set_ylabel('Normalized Importance (0-1)', fontsize=12)
    ax.set_title('Last 90 Days - Normalized Within Each Feature', fontsize=14)
    ax.legend(frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 90)
    ax.set_ylim(-0.05, 1.05)
    
    # Add vertical lines
    ax.axvline(x=1, color='red', linestyle=':', alpha=0.5)
    ax.axvline(x=7, color='orange', linestyle='--', alpha=0.5)
    ax.axvline(x=30, color='green', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/temporal_feature_importance_recent.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save temporal importance data to CSV with additional info
    temporal_df_data = []
    for lag_idx, lag in enumerate(lag_days):
        days_before = 365 - lag
        row = {'days_before_prediction': days_before}
        
        for feature in temporal_importance.keys():
            raw_value = temporal_importance[feature][lag_idx]
            row[f'{feature}_raw'] = raw_value
            row[f'{feature}_raw_plus_floor'] = raw_value + log_floor
            row[f'{feature}_is_zero'] = raw_value == 0
            row[f'{feature}_norm_within'] = normalized_importance_within[feature][lag_idx]
            row[f'{feature}_norm_across'] = normalized_importance_across[feature][lag_idx]
        
        temporal_df_data.append(row)
    
    temporal_df = pd.DataFrame(temporal_df_data)
    temporal_df.to_csv(f'{save_dir}/temporal_feature_importance_values.csv', index=False)
    
    # Print summary statistics
    print("\nTemporal Feature Importance Summary:")
    print(f"Log floor used: {log_floor:.6f}")
    for feature in temporal_importance.keys():
        values = temporal_importance[feature]
        n_zeros = sum(1 for v in values if v == 0)
        print(f"\n{feature}:")
        print(f"  Zero values: {n_zeros}/{len(values)} ({100*n_zeros/len(values):.1f}%)")
        print(f"  Range: [{min(values):.6f}, {max(values):.6f}]")
        print(f"  Mean: {np.mean(values):.6f}, Std: {np.std(values):.6f}")
    
    print(f"\nSaved temporal feature importance plots to {save_dir}")
    print(f"  - Raw values: temporal_feature_importance_raw.png")
    print(f"  - Normalized values: temporal_feature_importance_normalized.png")
    print(f"  - Recent period focus: temporal_feature_importance_recent.png")
    print(f"  - CSV data: temporal_feature_importance_values.csv")


# Main execution function
def run_comprehensive_analysis(model_dir, watershed_df, bucket_name, base_data_dir, base_attr_dir, n_test_watersheds=4):
    """Run all analyses"""
    
    # Load model and data
    print("Loading model and data...")
    model, data_module, norm_params, model_info, results = load_model_and_data(
        model_dir, watershed_df, bucket_name, base_data_dir, base_attr_dir
    )
    
    print(f"\nModel loaded successfully!")
    print(f"Test KGE from training: {results['test_results'][0]['test_kge']:.4f}")
    print(f"Best Validation KGE: {model_info['best_score']:.4f}")
    
    # 1. Generate predictions and plot hydrographs
    print(f"\n1. Generating predictions for {n_test_watersheds} test watersheds...")
    predictions_dict = generate_predictions_for_watersheds(model, data_module, n_test_watersheds)
    
    print("\n2. Plotting hydrographs...")
    plot_hydrographs(predictions_dict)
    
    # 2. Calculate feature importance
    print("\n3. Calculating feature importance (this may take a few minutes)...")
    feature_importance, baseline_kge = calculate_feature_importance_perturbation(
        model, data_module, n_samples=1000
    )
    
    print(f"\nBaseline validation KGE: {baseline_kge:.4f}")
    print("\nTop 10 most important features:")
    sorted_importance = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)[:10]
    for feat, imp in sorted_importance:
        print(f"  {feat}: {imp:.4f}")
    
    print("\n4. Plotting feature importance barplots...")
    plot_feature_importance_bars(feature_importance, data_module)
    
    # 3. Calculate temporal feature importance
    print("\n5. Calculating temporal feature importance...")
    temporal_importance, lag_days, baseline_kge_temporal = calculate_temporal_feature_importance(
        model, data_module, n_samples=500
    )
    
    print("\n6. Plotting temporal feature importance...")
    plot_temporal_feature_importance(temporal_importance, lag_days)
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("Generated outputs:")
    print("  - Hydrograph plots in ./hydrograph_plots/")
    print("  - Feature importance plots in ./feature_importance_plots/")
    print("="*60)
    
    return predictions_dict, feature_importance, temporal_importance


import shutil
import os

# Your model directory
model_dir = './model_analysis_pytorch-training-2025-08-06-23-42-19-294'

# Create output_data directory if it doesn't exist
os.makedirs(os.path.join(model_dir, 'output_data'), exist_ok=True)

# Copy files from output_files to output_data
files_to_copy = ['normalization_params.json', 'results.json']

for file in files_to_copy:
    src = os.path.join(model_dir, 'output_files', file)
    dst = os.path.join(model_dir, 'output_data', file)
    
    if os.path.exists(src):
        shutil.copy2(src, dst)
        print(f"✓ Copied {file} to output_data/")
    else:
        print(f"✗ {file} not found in output_files/")


import sys
import os
import pandas as pd
from pathlib import Path

# Add the directory containing train.py to Python path
sys.path.append('.')

# Import the analysis functions
#from analyze_outputs import run_comprehensive_analysis

# Import necessary functions from train.py
from train import identify_all_available_watersheds, get_watershed_attributes

# Configuration - UPDATE THESE WITH YOUR ACTUAL VALUES
BUCKET_NAME = 'climate-ai-data-science-shiny-app-data'#bucket_name#'sagemaker-us-east-2-672980278503'  # Your actual S3 bucket
BASE_DATA_DIR = 'Caravan-Jan25-nc/timeseries/netcdf'
BASE_ATTR_DIR = 'Caravan-Jan25-nc/attributes'
MODEL_DIR = './model_analysis_pytorch-training-2025-08-06-23-42-19-294'  # Your model directory

def main():
    # First, verify the model directory exists and has required files
    print(f"Checking model directory: {MODEL_DIR}")
    if not os.path.exists(MODEL_DIR):
        print(f"ERROR: Model directory {MODEL_DIR} not found!")
        print("Make sure you've run lossTracker.py first to download the model files.")
        return
    
    # Check for required files
    required_files = [
        'model_files/model.ckpt',
        'model_files/model_info.json',
        'output_data/normalization_params.json',
        'output_data/results.json'
    ]
    
    missing_files = []
    for file in required_files:
        file_path = os.path.join(MODEL_DIR, file)
        if not os.path.exists(file_path):
            missing_files.append(file)
        else:
            print(f"✓ Found {file}")
    
    if missing_files:
        print(f"\nERROR: Missing required files:")
        for file in missing_files:
            print(f"  ✗ {file}")
        print("\nMake sure you've run lossTracker.py successfully first.")
        return
    
    # Load watershed information
    print("\nLoading watershed information from S3...")
    print(f"Bucket: {BUCKET_NAME}")
    print(f"Data directory: {BASE_DATA_DIR}")
    
    try:
        watershed_df = identify_all_available_watersheds(BUCKET_NAME, BASE_DATA_DIR)
        print(f"Total watersheds found: {len(watershed_df)}")
        
        # Filter to CAMELS watersheds
        watershed_df = watershed_df[watershed_df['subdirectory_name'] == 'camels']
        print(f"CAMELS watersheds: {len(watershed_df)}")
        
        # Use same subset as training (first 100)
        watershed_df = watershed_df.iloc[0:100]
        print(f"Using first {len(watershed_df)} watersheds for analysis")
        
        if len(watershed_df) == 0:
            print("ERROR: No watersheds found! Check your S3 bucket and path.")
            return
            
    except Exception as e:
        print(f"ERROR loading watershed data: {e}")
        print("\nTroubleshooting:")
        print("1. Check your AWS credentials are configured")
        print("2. Verify the S3 bucket name and paths are correct")
        print("3. Ensure you have read access to the S3 bucket")
        return
    
    # Run comprehensive analysis
    try:
        predictions_dict, feature_importance, temporal_importance = run_comprehensive_analysis(
            model_dir=MODEL_DIR,
            watershed_df=watershed_df,
            bucket_name=BUCKET_NAME,
            base_data_dir=BASE_DATA_DIR,
            base_attr_dir=BASE_ATTR_DIR,
            n_test_watersheds=4  # Analyze 4 test watersheds
        )
        
        print("\nAnalysis complete! Check the output directories for results.")
        
    except Exception as e:
        print(f"\nERROR during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
