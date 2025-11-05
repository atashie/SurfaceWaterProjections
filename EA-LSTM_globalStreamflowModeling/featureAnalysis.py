import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import json
from typing import Dict, List, Tuple, Optional, Any
from tqdm import tqdm
import warnings
import os
from pathlib import Path
import gc
import sys
from train import get_watershed_attributes
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")



def find_split_csv(model_dir: str, experiment_name: str) -> str | None:
    """
    Try to locate the saved PCA split CSV from training outputs. We search:
    - ./output_files/{experiment_name}_pca_split_log.csv
    - Any file matching *_pca_split_log.csv under ./output_files
    Returns the first match or None.
    """
    # Primary expected path
    primary = os.path.join(model_dir, 'output_files', f'{experiment_name}_pca_split_log.csv')
    if os.path.exists(primary):
        return primary

    # Fallback: any split log under output_files
    candidates = list(Path(os.path.join(model_dir, 'output_files')).rglob('*_pca_split_log.csv'))
    if candidates:
        # Prefer the one whose filename contains the experiment_name substring
        for c in candidates:
            if experiment_name in c.name:
                return str(c)
        # Else return the first one
        return str(candidates[0])

    return None


def apply_training_split(watershed_df: pd.DataFrame, split_csv_path: str) -> pd.DataFrame:
    """
    Merge the training split into watershed_df. Align types and join keys robustly.
    Preferred join key is 'gauge_id' if present; else use subdirectory_name + watershedID.
    Returns a new DataFrame with a 'split' column and logs merge stats.
    """
    split_df = pd.read_csv(split_csv_path)
    merged = watershed_df.copy()

    # Build/align gauge_id on both sides if available
    if 'gauge_id' in split_df.columns:
        # Normalize types
        merged['watershedID'] = merged['watershedID'].astype(str)
        merged['subdirectory_name'] = merged['subdirectory_name'].astype(str)
        merged['gauge_id'] = merged['subdirectory_name'] + '_' + merged['watershedID']

        split_df['gauge_id'] = split_df['gauge_id'].astype(str)

        before = len(merged)
        merged = merged.merge(split_df[['gauge_id', 'split']], on='gauge_id', how='left')
        after = merged['split'].notna().sum()
        missing = merged['split'].isna().sum()
        print(f"Applied split via 'gauge_id': matched={after}, missing={missing}, total={before}")
    else:
        # Fallback join on subdirectory_name + watershedID
        # Align dtypes
        merged['watershedID'] = merged['watershedID'].astype(str)
        split_df['watershedID'] = split_df['watershedID'].astype(str)
        merged['subdirectory_name'] = merged['subdirectory_name'].astype(str)
        split_df['subdirectory_name'] = split_df['subdirectory_name'].astype(str)

        before = len(merged)
        merged = merged.merge(split_df[['subdirectory_name', 'watershedID', 'split']],
                              on=['subdirectory_name', 'watershedID'],
                              how='left')
        after = merged['split'].notna().sum()
        missing = merged['split'].isna().sum()
        print(f"Applied split via ['subdirectory_name','watershedID']: matched={after}, missing={missing}, total={before}")

    # Drop rows without a split (they weren’t part of the training run)
    if merged['split'].isna().any():
        missing = merged['split'].isna().sum()
        print(f"Warning: {missing} watersheds missing split after merge. Dropping them for analysis.")
        merged = merged.dropna(subset=['split']).reset_index(drop=True)

    # Show per-split counts
    print("Split counts after merge:", merged['split'].value_counts().to_dict())

    return merged


def apply_experiment_split_to_watersheds(watershed_df, model_dir, results, bucket_name, base_attr_dir):
    """
    Ensure watershed_df carries the 'split' column consistent with the training run:
    - Try to load the saved PCA split log from the extracted output_files dir
    - If not found, recompute the PCA split using the same function as train.py
    Returns a copy of watershed_df with a 'split' column.
    """
    # Determine experiment_name
    exp_name = None
    if isinstance(results, dict):
        exp_name = results.get('experiment_name', None)

    # Candidate path to the split log saved during training
    split_csv = None
    if exp_name:
        candidate = os.path.join(model_dir, 'output_files', f'{exp_name}_pca_split_log.csv')
        if os.path.exists(candidate):
            split_csv = candidate

    if split_csv:
        print(f"Loading training split from: {split_csv}")
        split_df = pd.read_csv(split_csv)
        # Only keep the join keys and split column
        cols_keep = ['subdirectory_name', 'watershedID', 'split']
        # Backward compatibility if IDs are strings/ints
        split_df['watershedID'] = split_df['watershedID'].astype(str)
        merged = watershed_df.copy()
        merged['watershedID'] = merged['watershedID'].astype(str)
        merged = merged.merge(split_df[cols_keep], on=['subdirectory_name', 'watershedID'], how='left')

        if merged['split'].isna().any():
            n_missing = merged['split'].isna().sum()
            print(f"Warning: {n_missing} watersheds missing split after merge. You may have filtered differently.")
            # Drop those or assign 'train' as default; here we drop them:
            merged = merged.dropna(subset=['split']).reset_index(drop=True)
        return merged

    # Fallback: recompute PCA-based split to best approximate training
    print("Training split CSV not found in outputs. Recomputing PCA-based split to approximate training...")
    # Import function from train.py
    import sys
    sys.path.append('.')
    from train import pca_stratified_reorder_watersheds

    recomputed = pca_stratified_reorder_watersheds(
        watershed_df=watershed_df,
        bucket_name=bucket_name,
        base_attr_dir=base_attr_dir,
        grid_size=5,             # default; optionally parse from results if saved
        train_frac=0.6,
        val_frac=0.2,
        test_frac=0.2,
        random_state=42,
        subdirectory_col='subdirectory_name',
        id_col='watershedID',
        log_path=None
    )
    return recomputed


from pathlib import Path

def _find_checkpoint(model_dir: str, preference: str = 'best') -> str:
    """
    Return a checkpoint path under model_dir based on preference.
    - 'best': use ./model_files/model.ckpt (packaged best)
    - 'last': try to locate any '*last*.ckpt' recursively; fallback to 'best'
    - explicit file path: if 'preference' is a file path, return it if exists
    """
    # If an explicit file path is provided
    if os.path.isfile(preference):
        return preference

    # Prefer a "last" checkpoint if requested and present
    if preference.lower() == 'last':
        search_roots = [
            os.path.join(model_dir, 'model_files'),
            os.path.join(model_dir, 'output_files'),
        ]
        for root in search_roots:
            if os.path.isdir(root):
                for ckpt in Path(root).rglob('*.ckpt'):
                    if 'last' in ckpt.name.lower():
                        return str(ckpt)

    # Fallback: packaged best checkpoint
    best_ckpt = os.path.join(model_dir, 'model_files', 'model.ckpt')
    if os.path.isfile(best_ckpt):
        return best_ckpt

    # Last resort: any .ckpt under model_dir
    for ckpt in Path(model_dir).rglob('*.ckpt'):
        return str(ckpt)

    raise FileNotFoundError("No checkpoint file found under model_dir.")


def load_model_and_data(
    model_dir: str,
    watershed_df: pd.DataFrame,
    bucket_name: str,
    base_data_dir: str,
    base_attr_dir: str,
    ckpt_preference: str = 'best',
    use_cache: bool = False,
    cache_dir: str | None = None,
    batch_size: int = 32,
    num_workers: int = 4
):
    """Load trained model and prepare data respecting the training split."""
    
    # Handle the nested model_files issue
    # Check if files exist at the expected locations first
    model_info_candidates = [
        os.path.join(model_dir, 'model_files', 'model_files', 'model_info.json'),  # Double nested
        os.path.join(model_dir, 'model_files', 'model_info.json'),  # Single nested
        os.path.join(model_dir, 'model_info.json'),  # Direct
    ]
    
    model_info_path = None
    for candidate in model_info_candidates:
        if os.path.exists(candidate):
            model_info_path = candidate
            print(f"Found model_info.json at: {model_info_path}")
            break
    
    if not model_info_path:
        raise FileNotFoundError(f"Could not find model_info.json in any expected location under {model_dir}")
    
    # Determine the actual base path based on where we found model_info.json
    if 'model_files/model_files' in model_info_path:
        # Double nested case
        model_files_base = os.path.join(model_dir, 'model_files', 'model_files')
        output_files_base = os.path.join(model_dir, 'output_files')
    elif model_info_path.endswith('model_files/model_info.json'):
        # Single nested case (expected)
        model_files_base = os.path.join(model_dir, 'model_files')
        output_files_base = os.path.join(model_dir, 'output_files')
    else:
        # Direct case
        model_files_base = model_dir
        output_files_base = model_dir
    
    # Now look for the other files relative to the determined base paths
    norm_params_candidates = [
        os.path.join(output_files_base, 'normalization_params.json'),
        os.path.join(model_dir, 'normalization_params.json'),
        os.path.join(model_dir, 'output_files', 'normalization_params.json'),
    ]
    
    norm_params_path = None
    for candidate in norm_params_candidates:
        if os.path.exists(candidate):
            norm_params_path = candidate
            print(f"Found normalization_params.json at: {norm_params_path}")
            break
    
    if not norm_params_path:
        raise FileNotFoundError(f"Could not find normalization_params.json in any expected location")
    
    results_candidates = [
        os.path.join(output_files_base, 'results.json'),
        os.path.join(model_dir, 'results.json'),
        os.path.join(model_dir, 'output_files', 'results.json'),
    ]
    
    results_path = None
    for candidate in results_candidates:
        if os.path.exists(candidate):
            results_path = candidate
            print(f"Found results.json at: {results_path}")
            break
    
    # It's okay if results.json doesn't exist
    if not results_path:
        print("Warning: results.json not found, continuing without it")
        results = {}
    else:
        with open(results_path, 'r') as f:
            results = json.load(f)
    
    # Load the JSON files
    with open(norm_params_path, 'r') as f:
        norm_params = json.load(f)
    with open(model_info_path, 'r') as f:
        model_info = json.load(f)
    
    # Import classes from train.py
    import sys
    sys.path.append('.')
    from train import StreamflowLightningModule, CaravanDataModule
    
    # Find the checkpoint file
    ckpt_candidates = [
        os.path.join(model_files_base, 'model.ckpt'),
        os.path.join(model_dir, 'model.ckpt'),
        os.path.join(model_dir, 'model_files', 'model.ckpt'),
    ]
    
    ckpt_path = None
    for candidate in ckpt_candidates:
        if os.path.exists(candidate):
            ckpt_path = candidate
            print(f"Found model checkpoint at: {ckpt_path}")
            break
    
    if not ckpt_path:
        # Try to find any .ckpt file
        print("Warning: model.ckpt not found at expected locations, searching for any .ckpt file...")
        for root, dirs, files in os.walk(model_dir):
            for file in files:
                if file.endswith('.ckpt'):
                    ckpt_path = os.path.join(root, file)
                    print(f"Found checkpoint: {ckpt_path}")
                    break
            if ckpt_path:
                break
    
    if not ckpt_path:
        # Last resort: try .pt file
        pt_candidates = [
            os.path.join(model_files_base, 'model.pt'),
            os.path.join(model_dir, 'model.pt'),
            os.path.join(model_dir, 'model_files', 'model.pt'),
        ]
        
        pt_path = None
        for candidate in pt_candidates:
            if os.path.exists(candidate):
                pt_path = candidate
                print(f"Found model.pt at: {pt_path}")
                break
        
        if pt_path:
            print(f"Loading from state dict: {pt_path}")
            # Create model and load state dict
            model_config = {
                'dynamic_input_size': len(model_info.get('feature_info', {}).get('dynamic_features', [])) or 8,
                'static_input_size': len(model_info.get('feature_info', {}).get('static_features', [])) or 50,
                'hidden_size': 256,
                'num_layers': 1,
                'dropout': 0.2,
                'learning_rate': 0.0001,
                'norm_params': norm_params
            }
            model = StreamflowLightningModule(**model_config)
            state_dict = torch.load(pt_path, map_location='cpu')
            model.load_state_dict(state_dict)
        else:
            raise FileNotFoundError(f"Could not find any model checkpoint (.ckpt or .pt) in {model_dir}")
    else:
        print(f"Loading checkpoint: {ckpt_path}")
        model = StreamflowLightningModule.load_from_checkpoint(ckpt_path, map_location='cpu')
    
    setattr(model, 'norm_params', norm_params)
    model.eval()
    
    # Apply the exact training split (if split CSV is present and needed)
    exp_name = results.get('experiment_name', None) if isinstance(results, dict) else None
    watersheds_df_with_split = watershed_df.copy()
    
    if 'split' in watersheds_df_with_split.columns:
        print("Info: Input watershed_df already contains a 'split' column. Using it as-is.")
    else:
        split_csv = find_split_csv(model_dir, exp_name) if exp_name else None
        if split_csv:
            print(f"Loading training split from: {split_csv}")
            watersheds_df_with_split = apply_training_split(watersheds_df_with_split, split_csv)
        else:
            print("Warning: Training split CSV not found; DataModule will fallback to its own split.")

    
    # Cache directory
    if cache_dir is not None and use_cache:
        os.makedirs(cache_dir, exist_ok=True)
    
    # Create data module (it honors the 'split' column if present)
    data_module = CaravanDataModule(
        watersheds_df=watersheds_df_with_split,
        bucket_name=bucket_name,
        base_data_dir=base_data_dir,
        base_attr_dir=base_attr_dir,
        sequence_length=365,
        batch_size=batch_size,
        num_workers=num_workers,
        chunk_size=50,
        train_split=0.6,
        val_split=0.2,
        random_seed=42,
        norm_params=norm_params,
        cache_data=use_cache,
        cache_dir=cache_dir
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

def calculate_feature_importance_perturbation(model, data_module, n_eff=1000, debug=True):
    """Calculate feature importance using perturbation method with debugging"""

    from train import calculate_kge, unnormalize_predictions

    device = next(model.parameters()).device

    # Get validation dataloader (FIXED)
    val_loader = data_module.val_dataloader()

    # Collect baseline predictions
    baseline_preds = []
    all_targets = []
    all_dynamic_seqs = []
    all_static_feats = []

    print("Collecting baseline predictions...")
    sample_count = 0

    with torch.no_grad():
        for batch in val_loader:
            if sample_count >= n_eff:
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
    if len(baseline_preds) == 0:
        raise RuntimeError("No validation samples collected for feature importance.")
    
    baseline_preds = torch.cat(baseline_preds)
    all_targets = torch.cat(all_targets)
    all_dynamic_seqs = torch.cat(all_dynamic_seqs)
    all_static_feats = torch.cat(all_static_feats)
    
    n_eff = int(min(n_eff, baseline_preds.shape[0]))
    baseline_preds = baseline_preds[:n_eff]
    all_targets = all_targets[:n_eff]
    all_dynamic_seqs = all_dynamic_seqs[:n_eff]
    all_static_feats = all_static_feats[:n_eff]

    
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
        perm_indices = torch.randperm(n_eff)
        perturbed_dynamic[:, :, idx] = perturbed_dynamic[perm_indices, :, idx]
        
        # Get predictions with perturbed feature
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_eff, batch_size):
                end_idx = min(i + batch_size, n_eff)
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
    
    for idx, feat_name in enumerate(tqdm(static_features)):#[:30])):
        # Check if feature has any variance
        feat_values = all_static_feats[:, idx]
        if feat_values.std() < 1e-6:
            print(f"    Warning: {feat_name} has near-zero variance (std={feat_values.std():.6f})")
            static_importance_shuffle[feat_name] = 0.0
            continue
        
        # Perturb this feature
        perturbed_static = all_static_feats.clone()
        
        # Shuffle values across samples
        perm_indices = torch.randperm(n_eff)
        perturbed_static[:, idx] = perturbed_static[perm_indices, idx]
        
        # Get predictions with perturbed feature
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_eff, batch_size):
                end_idx = min(i + batch_size, n_eff)
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
    
    for idx, feat_name in enumerate(tqdm(static_features)):#[:30])):
        # Set feature to its mean value
        perturbed_static = all_static_feats.clone()
        mean_val = all_static_feats[:, idx].mean()
        perturbed_static[:, idx] = mean_val
        
        # Get predictions
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_eff, batch_size):
                end_idx = min(i + batch_size, n_eff)
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
    
    for idx, feat_name in enumerate(tqdm(static_features)):#[:30])):
        # Add Gaussian noise
        perturbed_static = all_static_feats.clone()
        feat_std = all_static_feats[:, idx].std()
        if feat_std > 1e-6:
            noise = torch.randn(n_eff) * feat_std
            perturbed_static[:, idx] = all_static_feats[:, idx] + noise
        else:
            # If no variance, add small noise
            noise = torch.randn(n_eff) * 0.1
            perturbed_static[:, idx] = all_static_feats[:, idx] + noise
        
        # Get predictions
        with torch.no_grad():
            perturbed_preds = []
            batch_size = 64
            for i in range(0, n_eff, batch_size):
                end_idx = min(i + batch_size, n_eff)
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
    for feat_name in static_features:#[:30]:
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
            for i in range(0, min(200, n_eff), batch_size):  # Test on subset
                end_idx = min(i + batch_size, n_eff)
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
        'total_likely_snow_sum',  
        'total_likely_rain_sum',   
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

def calculate_temporal_feature_importance(model, data_module, n_eff=500, lag_days=None):
    """Calculate feature importance at different time lags"""
    from train import calculate_kge, unnormalize_predictions

    device = next(model.parameters()).device

    if lag_days is None:
        lag_days = list(range(0, 365, 10)) + [364]

    val_loader = data_module.val_dataloader()

    # Collect CPU tensors in lists
    dyn_list, stat_list, targ_list = [], [], []
    print("Collecting validation data...")
    collected = 0
    with torch.no_grad():
        for batch in val_loader:
            if collected >= n_eff:
                break
            dyn, stat, targ = batch
            dyn_list.append(dyn)    # keep on CPU
            stat_list.append(stat)  # keep on CPU
            targ_list.append(targ)  # keep on CPU
            collected += len(targ)

    if not dyn_list:
        raise RuntimeError("No validation samples collected for temporal importance.")

    # Concatenate and trim to effective count
    all_dynamic_seqs = torch.cat(dyn_list, dim=0)
    all_static_feats = torch.cat(stat_list, dim=0)
    all_targets = torch.cat(targ_list, dim=0)

    n_eff = int(min(n_eff, all_dynamic_seqs.shape[0]))
    all_dynamic_seqs = all_dynamic_seqs[:n_eff]
    all_static_feats = all_static_feats[:n_eff]
    all_targets = all_targets[:n_eff]

    # Baseline predictions
    with torch.no_grad():
        baseline_preds = []
        batch_size = 64
        for i in range(0, n_eff, batch_size):
            j = min(i + batch_size, n_eff)
            dyn_b = all_dynamic_seqs[i:j].to(device)
            stat_b = all_static_feats[i:j].to(device)
            pr = model(dyn_b, stat_b)
            baseline_preds.append(pr.cpu())
        baseline_preds = torch.cat(baseline_preds, dim=0)

    # Baseline KGE (original scale)
    baseline_preds_np = baseline_preds.numpy().flatten()
    targets_np = all_targets.numpy().flatten()
    baseline_preds_org = unnormalize_predictions(baseline_preds_np, data_module.norm_params)
    targets_org = unnormalize_predictions(targets_np, data_module.norm_params)
    baseline_kge = calculate_kge(targets_org, baseline_preds_org)

    # Weather features to analyze
    weather_features = [
#        'potential_evaporation_sum_ERA5_LAND', 
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 
        'temperature_2m_mean', 
        'temperature_2m_min',
        'total_precipitation_sum',
        'total_likely_snow_sum',  
        'total_likely_rain_sum',   
        'dewpoint_temperature_2m_mean',
        'Precip_lagged_90day',
        'Precip_smoothed_5day'
    ]
    dynamic_features = data_module.dynamic_cols_no_target
    weather_indices = [i for i, f in enumerate(dynamic_features) if f in weather_features]
    weather_names = [dynamic_features[i] for i in weather_indices]

    temporal_importance = {name: [] for name in weather_names}
    print(f"\nCalculating temporal feature importance for {len(lag_days)} lags...")

    for lag in tqdm(lag_days):
        time_idx = 364 - lag  # last step is 364 in a 365-window
        for feat_idx, feat_name in zip(weather_indices, weather_names):
            perturbed_dynamic = all_dynamic_seqs.clone()
            perm = torch.randperm(n_eff)
            perturbed_dynamic[:, time_idx, feat_idx] = perturbed_dynamic[perm, time_idx, feat_idx]

            with torch.no_grad():
                perturbed_preds = []
                batch_size = 64
                for i in range(0, n_eff, batch_size):
                    j = min(i + batch_size, n_eff)
                    dyn_b = perturbed_dynamic[i:j].to(device)
                    stat_b = all_static_feats[i:j].to(device)
                    pr = model(dyn_b, stat_b)
                    perturbed_preds.append(pr.cpu())
                perturbed_preds = torch.cat(perturbed_preds, dim=0)

            pert_np = perturbed_preds.numpy().flatten()
            pert_org = unnormalize_predictions(pert_np, data_module.norm_params)
            pert_kge = calculate_kge(targets_org, pert_org)
            temporal_importance[feat_name].append(baseline_kge - pert_kge)

    return temporal_importance, lag_days, baseline_kge


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


def calculate_rmce(observed: np.ndarray, simulated: np.ndarray) -> float:
    """
    Calculate Cube-root Mean Cubed Error (RMCE).
    RMCE = (mean(|sim - obs|^3))^(1/3)
    """
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()
    
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated) & \
                   np.isfinite(observed) & np.isfinite(simulated)
    
    if np.sum(valid_indices) < 1:
        return np.nan
    
    obs_valid = observed[valid_indices]
    sim_valid = simulated[valid_indices]
    
    # Calculate cubed errors (preserve sign)
    cubed_errors = (sim_valid - obs_valid) ** 3
    mean_cubed = np.mean(cubed_errors)
    
    # Take cube root preserving sign
    if mean_cubed < 0:
        rmce = -((-mean_cubed) ** (1/3))
    else:
        rmce = mean_cubed ** (1/3)
    
    return abs(rmce)  # Return absolute value as error metric


def compute_train_style_metrics_streaming(model, dataloader, norm_params):
    """
    Match train.py epoch-end behavior but stream over batches:
      - loss_norm = MAE + RMSE + RMCE in normalized space (val_loss/test_loss)
      - KGE, r, alpha, beta, NSE, RMSE, MAE, RMCE in original scale
    No giant tensors are stored; everything is accumulated online.
    """
    device = next(model.parameters()).device

    # Normalized-space accumulators
    n_norm = 0
    sum_abs_err = 0.0
    sum_sq_err = 0.0
    sum_abs_err_cubed = 0.0  # for RMCE

    # Original-scale accumulators (for KGE/NSE/RMSE/MAE/RMCE)
    n_org = 0
    sum_obs = 0.0
    sum_sim = 0.0
    sum_obs2 = 0.0
    sum_sim2 = 0.0
    sum_obs_sim = 0.0
    sum_abs_err_org = 0.0
    sum_sq_err_org = 0.0
    sum_abs_err_cubed_org = 0.0

    from train import unnormalize_predictions

    model.eval()
    with torch.inference_mode():
        for batch in tqdm(dataloader, desc="Eval (train-style)", unit="batch"):
            dyn, stat, targ = batch
            dyn = dyn.to(device, non_blocking=True)
            stat = stat.to(device, non_blocking=True)

            pred = model(dyn, stat).detach().squeeze(-1)  # [B]
            targ = targ.squeeze(-1)  # [B]

            # Normalized space errors
            e = (pred - targ).float()
            bsz = e.numel()
            n_norm += bsz
            sum_abs_err += torch.sum(torch.abs(e)).item()
            sum_sq_err += torch.sum(e * e).item()
            sum_abs_err_cubed += torch.sum(torch.abs(e) ** 3).item()

            # Original scale (unnormalize batch-wise)
            pred_np = pred.cpu().numpy().astype(np.float64)
            targ_np = targ.cpu().numpy().astype(np.float64)
            sim_org = unnormalize_predictions(pred_np, norm_params).astype(np.float64)
            obs_org = unnormalize_predictions(targ_np, norm_params).astype(np.float64)

            # Pairwise valid mask
            valid = np.isfinite(obs_org) & np.isfinite(sim_org)
            if not np.any(valid):
                continue

            obs = obs_org[valid]
            sim = sim_org[valid]
            m = obs.size

            n_org += m
            sum_obs += float(np.sum(obs))
            sum_sim += float(np.sum(sim))
            sum_obs2 += float(np.sum(obs * obs))
            sum_sim2 += float(np.sum(sim * sim))
            sum_obs_sim += float(np.sum(obs * sim))
            abs_err_org = np.abs(sim - obs)
            sum_abs_err_org += float(np.sum(abs_err_org))
            sq_err_org = (sim - obs) ** 2
            sum_sq_err_org += float(np.sum(sq_err_org))
            sum_abs_err_cubed_org += float(np.sum(abs_err_org ** 3))

    # Normalized-space metrics
    if n_norm > 0:
        mae_norm = sum_abs_err / n_norm
        rmse_norm = np.sqrt(sum_sq_err / n_norm)
        rmce_norm = (sum_abs_err_cubed / n_norm) ** (1.0 / 3.0)
        loss_norm = mae_norm + rmse_norm + rmce_norm
    else:
        loss_norm = mae_norm = rmse_norm = rmce_norm = np.nan

    # Original-scale metrics (KGE components and others)
    if n_org > 1:
        mean_obs = sum_obs / n_org
        mean_sim = sum_sim / n_org
        var_obs = max(sum_obs2 / n_org - mean_obs ** 2, 0.0)
        var_sim = max(sum_sim2 / n_org - mean_sim ** 2, 0.0)
        std_obs = np.sqrt(var_obs)
        std_sim = np.sqrt(var_sim)

        if std_obs > 1e-12 and std_sim > 1e-12:
            cov = sum_obs_sim / n_org - mean_obs * mean_sim
            r = cov / (std_obs * std_sim)
        else:
            r = 0.0

        alpha = (std_sim / std_obs) if std_obs > 1e-12 else np.nan
        beta = (mean_sim / mean_obs) if abs(mean_obs) > 1e-12 else np.nan
        kge = 1.0 - np.sqrt((r - 1.0) ** 2 + (alpha - 1.0) ** 2 + (beta - 1.0) ** 2)

        rmse = np.sqrt(sum_sq_err_org / n_org)
        mae = sum_abs_err_org / n_org

        # NSE
        denom = sum_obs2 - n_org * (mean_obs ** 2)
        if denom > 1e-12:
            nse = 1.0 - (sum_sq_err_org / denom)
        else:
            nse = np.nan

        # RMCE (original scale)
        rmce = (sum_abs_err_cubed_org / n_org) ** (1.0 / 3.0)
    else:
        kge = r = alpha = beta = rmse = mae = nse = rmce = np.nan

    return {
        'loss_norm': float(loss_norm),
        'mae_norm': float(mae_norm),
        'rmse_norm': float(rmse_norm),
        'rmce_norm': float(rmce_norm),
        'kge': float(kge),
        'r': float(r),
        'alpha': float(alpha),
        'beta': float(beta),
        'rmse': float(rmse),
        'mae': float(mae),
        'nse': float(nse),
        'rmce': float(rmce),
        'n_eff': int(n_org)
    }



def calculate_comprehensive_metrics_optimized(model, data_module, output_dir='./metrics_analysis'):
    """
    Efficient metrics: median and percentile metrics from per-watershed evaluation only.
    Removed the continuous evaluation that was causing memory issues.
    """
    os.makedirs(output_dir, exist_ok=True)
    comprehensive_rows = []
    last_per_ws_df = None

    for split in ['val', 'test']:
        print(f"\n{'='*60}")
        print(f"Processing {split.upper()} split")
        print('='*60)

        # Per-watershed metrics only (no more continuous evaluation)
        print(f"\n1. Calculating per-watershed {split} metrics (optimized)...")
        per_ws_df = calculate_watershed_metrics_optimized(model, data_module, split=split)
        last_per_ws_df = per_ws_df

        if not per_ws_df.empty:
            ws_csv = os.path.join(output_dir, f'{split}_per_watershed_metrics.csv')
            per_ws_df.to_csv(ws_csv, index=False)
            print(f"Saved per-watershed metrics to {ws_csv}")

            metric_columns = ['kge', 'r', 'alpha', 'beta', 'rmse', 'mae', 'nse', 'rmce']
            median_metrics = per_ws_df[metric_columns].median()
            comprehensive_rows.append({
                'split': split,
                'aggregation': 'median_watersheds',
                **{m: float(median_metrics[m]) for m in metric_columns}
            })
            for percentile in [10, 25, 75, 90]:
                q = per_ws_df[metric_columns].quantile(percentile / 100.0)
                comprehensive_rows.append({
                    'split': split,
                    'aggregation': f'p{percentile}_watersheds',
                    **{m: float(q[m]) for m in metric_columns}
                })

    comprehensive_df = pd.DataFrame(comprehensive_rows)
    out_csv = os.path.join(output_dir, 'comprehensive_metrics.csv')
    comprehensive_df.to_csv(out_csv, index=False)
    print(f"\nSaved comprehensive metrics to {out_csv}")

    return comprehensive_df, last_per_ws_df


def calculate_watershed_metrics_optimized(model, data_module, split='val'):
    """
    Compute per-watershed metrics efficiently by loading each watershed once
    and batching sequences. Uses same normalization and dataset logic as training.
    """
    from train import (
        calculate_kge_components, calculate_rmse, calculate_mae,
        calculate_nse, unnormalize_predictions, WatershedDataset, prepare_ptf_dataframe
    )

    device = next(model.parameters()).device
    if split == 'val':
        watersheds = data_module.val_watersheds
    else:
        watersheds = data_module.test_watersheds

    # If cached dataframe exists, use it; else load once for the split
    if data_module.cache_data and split in data_module._cached_data and data_module._cached_data[split] is not None:
        all_data = data_module._cached_data[split]
    else:
        print(f"Loading all {split} watersheds data (single pass)...")
        all_data, _, _, _ = prepare_ptf_dataframe(
            watersheds,
            data_module.bucket_name,
            data_module.base_data_dir,
            data_module.base_attr_dir,
            norm_params=data_module.norm_params
        )

    if all_data is None or all_data.empty:
        return pd.DataFrame([])

    groups = all_data.groupby('group_id')
    results = []

    for _, ws_row in tqdm(watersheds.iterrows(), total=len(watersheds), desc=f"Per-watershed {split}"):
        ws_id = f"{ws_row['subdirectory_name']}_{ws_row['watershedID']}"
        if ws_id not in groups.groups:
            continue

        ws_df = groups.get_group(ws_id).sort_values('date')
        if len(ws_df) <= data_module.sequence_length:
            continue

        # Dataset and loader for this one watershed
        ws_dataset = WatershedDataset(
            ws_df,
            data_module.static_cols,
            data_module.dynamic_cols_no_target,
            'streamflow',
            data_module.sequence_length
        )
        if len(ws_dataset) == 0:
            continue

        ws_loader = torch.utils.data.DataLoader(
            ws_dataset, batch_size=256, shuffle=False, num_workers=0
        )

        preds_norm, targs_norm = [], []
        with torch.no_grad():
            for batch in ws_loader:
                dyn, stat, targ = batch
                dyn = dyn.to(device, non_blocking=True)
                stat = stat.to(device, non_blocking=True)
                pr = model(dyn, stat)
                preds_norm.append(pr.detach().cpu().numpy().flatten())
                targs_norm.append(targ.detach().cpu().numpy().flatten())

        preds_norm = np.concatenate(preds_norm, axis=0)
        targs_norm = np.concatenate(targs_norm, axis=0)

        if preds_norm.size < 10:
            continue

        # Unnormalize
        preds_orig = unnormalize_predictions(preds_norm, data_module.norm_params)
        targs_orig = unnormalize_predictions(targs_norm, data_module.norm_params)

        # Metrics
        kge, r, alpha, beta = calculate_kge_components(targs_orig, preds_orig)
        rmse = calculate_rmse(targs_orig, preds_orig)
        mae = calculate_mae(targs_orig, preds_orig)
        nse = calculate_nse(targs_orig, preds_orig)
        rmce = calculate_rmce(targs_orig, preds_orig)

        results.append({
            'watershed_id': ws_id,
            'subdirectory': ws_row['subdirectory_name'],
            'watershedID': ws_row['watershedID'],
            'split': split,
            'n_eff': int(preds_norm.size),
            'kge': float(kge),
            'r': float(r),
            'alpha': float(alpha),
            'beta': float(beta),
            'rmse': float(rmse),
            'mae': float(mae),
            'nse': float(nse),
            'rmce': float(rmce)
        })

        # Periodic cleanup
        if len(results) % 10 == 0:
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    return pd.DataFrame(results)




# Main execution function
def run_comprehensive_analysis(
    model_dir: str,
    watershed_df: pd.DataFrame,
    bucket_name: str,
    base_data_dir: str,
    base_attr_dir: str,
    n_test_watersheds: int = 4,
    ckpt_preference: str = 'best',  # <- NEW
    use_cache: bool = False,        # <- optional
    cache_dir: str | None = None,   # <- optional
    batch_size: int = 32,           # <- optional
    num_workers: int = 4            # <- optional
):
    """Run all analyses including comprehensive metrics"""

    print("Loading model and data...")
    model, data_module, norm_params, model_info, results = load_model_and_data(
        model_dir=model_dir,
        watershed_df=watershed_df,
        bucket_name=bucket_name,
        base_data_dir=base_data_dir,
        base_attr_dir=base_attr_dir,
        ckpt_preference=ckpt_preference,  # <- pass through
        use_cache=use_cache,              # <- pass through
        cache_dir=cache_dir,              # <- pass through
        batch_size=batch_size,            # <- pass through
        num_workers=num_workers           # <- pass through
    )

    # Quick consistency diagnostics
    def _count_samples(dataloader):
        total = 0
        for batch in dataloader:
            total += len(batch[2])  # target length
        return total

    val_loader = data_module.val_dataloader()
    test_loader = data_module.test_dataloader()
    val_samples = _count_samples(val_loader)
    test_samples = _count_samples(test_loader)
    print(f"Val samples (sequences): {val_samples}")
    print(f"Test samples (sequences): {test_samples}")

    print("\nModel loaded successfully!")
    if isinstance(results, dict) and results.get('test_results'):
        try:
            print(f"Test KGE from training: {results['test_results'][0].get('test_kge', float('nan')):.4f}")
        except Exception:
            pass
    try:
        print(f"Best Validation Score (monitor={model_info.get('monitor_metric', 'val_loss')}): {float(model_info.get('best_score', float('nan'))):.6f}")
    except Exception:
        pass

    # 1) Comprehensive metrics (aggregate + per-watershed)
    print("\n1. Calculating comprehensive metrics (optimized)...")
    comprehensive_metrics_df, watershed_metrics_df = calculate_comprehensive_metrics_optimized(
        model, data_module, output_dir='./metrics_analysis'
    )

    # 2) (Optional) Generate a few hydrographs
    print(f"\n2. Generating predictions for {n_test_watersheds} test watersheds...")
    predictions_dict = generate_predictions_for_watersheds(model, data_module, n_test_watersheds)

    print("\n3. Plotting hydrographs...")
    plot_hydrographs(predictions_dict)

    # 3) Feature importance
    print("\n4. Calculating feature importance (perturbation)...")
    feature_importance, baseline_kge = calculate_feature_importance_perturbation(
        model, data_module, n_eff=500
    )
    print(f"\nBaseline validation KGE: {baseline_kge:.4f}")

    print("\n5. Plotting feature importance barplots...")
    plot_feature_importance_bars(feature_importance, data_module)

    # 4) Temporal feature importance
    print("\n6. Calculating temporal feature importance...")
    temporal_importance, lag_days, baseline_kge_temporal = calculate_temporal_feature_importance(
        model, data_module, n_eff=500
    )

    print("\n7. Plotting temporal feature importance...")
    plot_temporal_feature_importance(temporal_importance, lag_days)

    print("\n" + "="*60)
    print("Analysis complete!")
    print("Generated outputs:")
    print("  - Comprehensive metrics in ./metrics_analysis/ (includes mse_norm for 1:1 with val_loss)")
    print("  - Hydrograph plots in ./hydrograph_plots/")
    print("  - Feature importance plots in ./feature_importance_plots/")
    print("="*60)

    # Add smoothing filter analysis if model has it
    if hasattr(model.model, 'feature_smoother'):
        print("\n8. Analyzing learned smoothing filters...")
        filter_analysis = analyze_learned_smoothing_filters(
            model, 
            data_module,
            save_dir='./filter_analysis'
        )
        print("Smoothing filter analysis complete!")
    else:
        print("\n8. Model does not use feature smoothing, skipping filter analysis")
        filter_analysis = None

    return predictions_dict, feature_importance, temporal_importance, filter_analysis



def analyze_learned_smoothing_filters(model, data_module, save_dir='./filter_analysis'):
    """
    Analyze the learned smoothing filters and compare with uniform filters.
    Includes frequency domain analysis.
    """
    import scipy.signal as signal
    import scipy.fft as fft
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Extract learned filters
    if not hasattr(model.model, 'feature_smoother'):
        print("Model does not have learnable smoothing. Skipping filter analysis.")
        return None
    
    filters = model.model.feature_smoother.get_filter_weights()
    n_features = model.model.feature_smoother.n_features
    kernel_sizes = model.model.feature_smoother.kernel_sizes
    feature_names = data_module.dynamic_cols_no_target
    
    # Analyze each feature's filters
    filter_analysis = {}
    
    for feat_idx in range(min(n_features, len(feature_names))):
        feat_name = feature_names[feat_idx] if feat_idx < len(feature_names) else f'feature_{feat_idx}'
        feat_filters = filters[f'feature_{feat_idx}']
        
        fig, axes = plt.subplots(len(kernel_sizes), 3, figsize=(15, 4*len(kernel_sizes)))
        if len(kernel_sizes) == 1:
            axes = axes.reshape(1, -1)
        
        fig.suptitle(f'Learned Smoothing Filters - {feat_name}', fontsize=14, fontweight='bold')
        
        analysis_data = {}
        
        for scale_idx, kernel_size in enumerate(kernel_sizes):
            learned_filter = feat_filters[f'kernel_{kernel_size}']
            uniform_filter = np.ones(kernel_size) / kernel_size
            
            # 1. Filter coefficients comparison
            ax = axes[scale_idx, 0]
            x = np.arange(kernel_size)
            ax.bar(x - 0.2, uniform_filter, 0.4, label='Uniform', alpha=0.7, color='blue')
            ax.bar(x + 0.2, learned_filter, 0.4, label='Learned', alpha=0.7, color='red')
            ax.set_xlabel('Coefficient Index')
            ax.set_ylabel('Weight')
            ax.set_title(f'Kernel Size {kernel_size} - Coefficients')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # 2. Frequency response
            ax = axes[scale_idx, 1]
            
            # Compute frequency response
            w_uniform, h_uniform = signal.freqz(uniform_filter, 1, worN=512)
            w_learned, h_learned = signal.freqz(learned_filter, 1, worN=512)
            
            # Convert to Hz (assuming daily sampling)
            freq_uniform = w_uniform / (2 * np.pi)  # Normalized frequency
            freq_learned = w_learned / (2 * np.pi)
            
            # Plot magnitude response
            ax.plot(freq_uniform, np.abs(h_uniform), 'b-', label='Uniform', linewidth=2)
            ax.plot(freq_learned, np.abs(h_learned), 'r-', label='Learned', linewidth=2)
            ax.set_xlabel('Normalized Frequency (cycles/day)')
            ax.set_ylabel('Magnitude')
            ax.set_title(f'Frequency Response')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_xlim([0, 0.5])  # Nyquist frequency
            
            # 3. Phase response
            ax = axes[scale_idx, 2]
            ax.plot(freq_uniform, np.angle(h_uniform), 'b-', label='Uniform', linewidth=2)
            ax.plot(freq_learned, np.angle(h_learned), 'r-', label='Learned', linewidth=2)
            ax.set_xlabel('Normalized Frequency (cycles/day)')
            ax.set_ylabel('Phase (radians)')
            ax.set_title(f'Phase Response')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_xlim([0, 0.5])
            
            # Store analysis metrics
            analysis_data[f'kernel_{kernel_size}'] = {
                'learned_filter': learned_filter.tolist(),
                'deviation_from_uniform': float(np.mean(np.abs(learned_filter - uniform_filter))),
                'filter_variance': float(np.var(learned_filter)),
                'effective_bandwidth': float(np.sum(np.abs(h_learned) > 0.5) / len(h_learned)),
                'is_causal': bool(np.all(learned_filter >= -1e-6)),  # Check if non-negative
            }
        
        plt.tight_layout()
        plt.savefig(f'{save_dir}/filters_{feat_name}.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        filter_analysis[feat_name] = analysis_data
    
    # Plot scale mixing weights heatmap
    scale_weights = filters['scale_weights']
    
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(scale_weights, aspect='auto', cmap='YlOrRd')
    
    # Labels
    ax.set_xticks(range(len(kernel_sizes)))
    ax.set_xticklabels([f'K={k}' for k in kernel_sizes])
    ax.set_yticks(range(len(feature_names)))
    ax.set_yticklabels(feature_names)
    
    ax.set_xlabel('Kernel Size')
    ax.set_ylabel('Feature')
    ax.set_title('Multi-Scale Mixing Weights (Softmax)', fontsize=14, fontweight='bold')
    
    # Add values
    for i in range(len(feature_names)):
        for j in range(len(kernel_sizes)):
            text = ax.text(j, i, f'{scale_weights[i, j]:.2f}',
                         ha='center', va='center', color='white' if scale_weights[i, j] > 0.5 else 'black')
    
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/scale_mixing_weights.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save analysis to JSON
    with open(f'{save_dir}/filter_analysis.json', 'w') as f:
        json.dump(filter_analysis, f, indent=2)
    
    # Summary statistics
    print("\n" + "="*60)
    print("LEARNED SMOOTHING FILTER ANALYSIS")
    print("="*60)
    
    for feat_name, feat_data in filter_analysis.items():
        print(f"\n{feat_name}:")
        for kernel_name, metrics in feat_data.items():
            deviation = metrics['deviation_from_uniform']
            print(f"  {kernel_name}: deviation from uniform = {deviation:.4f}")
    
    # Identify which features prefer which kernel sizes
    print("\n\nPreferred kernel sizes (by scale weight):")
    for feat_idx, feat_name in enumerate(feature_names):
        weights = scale_weights[feat_idx]
        preferred_idx = np.argmax(weights)
        preferred_size = kernel_sizes[preferred_idx]
        print(f"  {feat_name}: {preferred_size} (weight={weights[preferred_idx]:.3f})")
    
    return filter_analysis


def compare_smoothing_effectiveness(
    model_with_smoothing, 
    model_without_smoothing,
    data_module,
    n_samples=1000
):
    """
    Compare model performance with and without learned smoothing.
    """
    from train import calculate_kge, unnormalize_predictions
    
    device = next(model_with_smoothing.parameters()).device
    val_loader = data_module.val_dataloader()
    
    results = {'with_smoothing': [], 'without_smoothing': []}
    
    print("Comparing models with and without learned smoothing...")
    
    sample_count = 0
    with torch.no_grad():
        for batch in tqdm(val_loader, desc="Evaluating"):
            if sample_count >= n_samples:
                break
            
            dynamic_seq, static_feat, target = batch
            dynamic_seq = dynamic_seq.to(device)
            static_feat = static_feat.to(device)
            
            # Predictions with smoothing
            pred_with = model_with_smoothing(dynamic_seq, static_feat)
            pred_with_np = pred_with.cpu().numpy().flatten()
            
            # Predictions without smoothing
            pred_without = model_without_smoothing(dynamic_seq, static_feat)
            pred_without_np = pred_without.cpu().numpy().flatten()
            
            # Unnormalize
            target_np = target.numpy().flatten()
            pred_with_orig = unnormalize_predictions(pred_with_np, data_module.norm_params)
            pred_without_orig = unnormalize_predictions(pred_without_np, data_module.norm_params)
            target_orig = unnormalize_predictions(target_np, data_module.norm_params)
            
            # Calculate KGE for each
            kge_with = calculate_kge(target_orig, pred_with_orig)
            kge_without = calculate_kge(target_orig, pred_without_orig)
            
            results['with_smoothing'].append(kge_with)
            results['without_smoothing'].append(kge_without)
            
            sample_count += len(target)
    
    # Summary statistics
    mean_kge_with = np.mean(results['with_smoothing'])
    mean_kge_without = np.mean(results['without_smoothing'])
    improvement = mean_kge_with - mean_kge_without
    
    print(f"\nResults:")
    print(f"  Mean KGE with smoothing: {mean_kge_with:.4f}")
    print(f"  Mean KGE without smoothing: {mean_kge_without:.4f}")
    print(f"  Improvement: {improvement:.4f} ({100*improvement/abs(mean_kge_without):.1f}%)")
    
    return results

