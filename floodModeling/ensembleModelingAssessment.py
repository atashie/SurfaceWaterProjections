###############################################################
###############################################################
###############################################################
### Ensemble Modeling Assessment - Enhanced Aesthetics
###############################################################
###############################################################
###############################################################

"""
ensembleModelingAssessment.py

Compares model ensemble approaches:
1. XGBoost Baseline
2. Architecture Ensemble (6 models)
3. Feature Ensemble (5 configurations)
4. Transfer Ensemble (7 source basins)
5. Kitchen Sink (all unique models)

Uses majority voting from binarized predictions.
"""

import os
import pickle
import tempfile
import re
from typing import List, Dict, Tuple
import numpy as np
import pandas as pd
from osgeo import gdal
import boto3

import pyarrow as pa
import pyarrow.dataset as ds

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.patches import PathPatch

from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    matthews_corrcoef, roc_auc_score, average_precision_score
)

gdal.UseExceptions()

# ============ Configuration ============

BASE_S3_PREFIX = "climate-ai-data-science-datasets/arrakis-data/floodOutputs"

basin_locations = [
    {"latitude": "10.7875",  "longitude": "0.7479",  "name": "EMSR470"},
    {"latitude": "48.1747",  "longitude": "23.3978", "name": "EMSR444"},
    {"latitude": "53.6583",  "longitude": "-1.1398", "name": "EMSR407"},
    {"latitude": "-16.5646", "longitude": "46.7857", "name": "EMSR424"},
    {"latitude": "-1.8038",  "longitude": "33.3857", "name": "EMSR438"},
    {"latitude": "40.9437",  "longitude": "24.7669", "name": "EMSR292"},
    {"latitude": "27.2981",  "longitude": "68.1875", "name": "EMSR629"},
    {"latitude": "41.9275",  "longitude": "-1.3792", "name": "EMSR279"},
]

ARCHITECTURE_SUBFOLDERS = {
    "GLM":        "assessGLMModel",
    "RandomForest": "assessRandomForestModel",
    "XGBoost":    "assessXGBoostModel",
    "MLP-FiLM":   "assessMLPModel",
    "U-Net-FiLM": "assessUNetModel",
    "TabNet":     "assessTabNetModel",
}

FEATURE_SUBFOLDERS = {
    "XGB - All Features":    "assessXGBoostModel",
    "XGB - Basic Features":  "assessXGBoostTraditionalModel",
    "XGB - No Streamflow":   "assessXGBoostWeatherModel",
    "XGB - Multi Objective": "assessXGBoostMultiLearnModel",
    "XGB - 25% Noise":       "assessXGBoostNoisyModel"
}

ENSEMBLE_ORDER = ["XGBoost Baseline", "Architecture Ensemble", "Feature Ensemble", 
                  "Transfer Ensemble", "Kitchen Sink"]

METRICS = ["MCC", "F1", "Accuracy", "Precision", "Recall", "ROC_AUC", "PR_AUC"]
OUTPUT_DIR = "ensembleComparisons"
ALPHA = 0.05


# ============ Parquet Utilities ============

def get_s3fs():
    try:
        import pyarrow.fs as pafs
        return pafs.S3FileSystem()
    except:
        from pyarrow import filesystem as pafs
        return pafs.S3FileSystem()


def open_parquet_dataset(parquet_uri: str):
    if parquet_uri.startswith("s3://"):
        bkt, key = parquet_uri[5:].split("/", 1)
        s3fs = get_s3fs()
        return ds.dataset(f"{bkt}/{key}", format="parquet", filesystem=s3fs)
    else:
        return ds.dataset(parquet_uri, format="parquet")


def load_test_data_for_location(latitude: str, longitude: str,
                                bucket: str = "climate-ai-data-science-datasets",
                                base_prefix: str = "arrakis-data/floodOutputs") -> pd.DataFrame:
    """Load test set data from parquet for a specific target watershed."""
    
    parquet_uri = (
        f"s3://{bucket}/{base_prefix}/"
        f"trainingFolderFor_Lat_{latitude}_Lon_{longitude}/"
        f"postProcessedTable/flood_training.parquet"
    )
    
    dataset = open_parquet_dataset(parquet_uri)
    target_watershed = f"Lat_{latitude}_Lon_{longitude}"
    
    # Load data for target watershed
    filt = ds.field("Watershed_Loc") == target_watershed
    scanner = dataset.scanner(filter=filt)
    
    pixel_data = []
    for batch in scanner.to_batches():
        if batch.num_rows > 0:
            df = pa.Table.from_batches([batch]).to_pandas()
            pixel_data.append(df)
    
    if not pixel_data:
        return pd.DataFrame()
    
    df_all = pd.concat(pixel_data, ignore_index=True)
    
    # Normalize dates
    def normalize_date(d):
        if pd.isna(d):
            return None
        if isinstance(d, str):
            return d.replace('-', '').replace('/', '')[:8]
        if isinstance(d, (pd.Timestamp, np.datetime64)):
            return pd.Timestamp(d).strftime('%Y%m%d')
        return str(d)
    
    df_all['Date_Normalized'] = df_all['Date'].apply(normalize_date)
    
    return df_all


# ============ TIFF Loading ============

def load_prediction_tiff(tiff_path: str) -> Tuple[Dict[str, np.ndarray], Dict]:
    """Load prediction TIFF and return dict of date -> 2D prediction array."""
    
    vsi_path = tiff_path.replace('s3://', '/vsis3/')
    
    try:
        ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
        if ds is None:
            print(f"      ✗ GDAL cannot open file")
            return {}, None
        
        predictions_by_date = {}
        
        # Try multiple date patterns
        date_patterns = [
            r'(\d{4}-\d{2}-\d{2})',  # YYYY-MM-DD
            r'(\d{4}_\d{2}_\d{2})',  # YYYY_MM_DD  
            r'(\d{8})',              # YYYYMMDD
        ]
        
        for b in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(b)
            desc = band.GetDescription()
            
            date_found = None
            
            # Try each pattern
            for pattern in date_patterns:
                match = re.search(pattern, desc)
                if match:
                    date_str = match.group(1)
                    # Normalize to YYYYMMDD format
                    date_found = date_str.replace('-', '').replace('_', '')
                    break
            
            if date_found and len(date_found) == 8:
                pred_array = band.ReadAsArray()
                predictions_by_date[date_found] = pred_array
        
        # Store metadata
        metadata = {
            'geotransform': ds.GetGeoTransform(),
            'projection': ds.GetProjection(),
            'width': ds.RasterXSize,
            'height': ds.RasterYSize,
            'n_bands': ds.RasterCount
        }
        
        ds = None
        
        if predictions_by_date:
            return predictions_by_date, metadata
        else:
            print(f"      ✗ No dates extracted from {ds.RasterCount} bands")
            return {}, None
        
    except Exception as e:
        print(f"      ✗ Error loading: {e}")
        return {}, None

# ============ Transfer Model Predictions ============

def load_model_from_location(latitude: str, longitude: str,
                             bucket: str = "climate-ai-data-science-datasets",
                             base_prefix: str = "arrakis-data/floodOutputs") -> Dict:
    """Load pre-trained XGBoost model."""
    
    model_path = (
        f"s3://{bucket}/{base_prefix}/"
        f"trainingFolderFor_Lat_{latitude}_Lon_{longitude}/"
        f"modelXGBoost/xgb_model.pkl"
    )
    
    s3_client = boto3.client('s3')
    model_key = model_path.replace(f"s3://{bucket}/", "")
    
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        s3_client.download_file(bucket, model_key, temp_path)
        
        with open(temp_path, 'rb') as f:
            results = pickle.load(f)
        
        os.remove(temp_path)
        return results
        
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return None


def predict_with_transfer_model(source_model: Dict, target_data: pd.DataFrame,
                                feature_cols: List[str], means: np.ndarray,
                                stds: np.ndarray) -> np.ndarray:
    """Generate predictions using transfer model on target data."""
    
    # Align features
    missing_features = [f for f in feature_cols if f not in target_data.columns]
    for feat in missing_features:
        target_data[feat] = np.nan
    
    X = target_data[feature_cols].to_numpy(dtype=np.float32, copy=False)
    
    # Filter rows with sufficient features
    non_nan_counts = np.sum(~np.isnan(X), axis=1)
    feature_threshold = int(0.70 * len(feature_cols))
    valid_mask = non_nan_counts >= feature_threshold
    
    predictions = np.full(len(target_data), np.nan, dtype=np.float32)
    
    print(f"      Transfer predict: {valid_mask.sum()}/{len(valid_mask)} rows have ≥{feature_threshold} features")
    
    if valid_mask.sum() == 0:
        print(f"      ⚠️ No valid rows - all predictions will be NaN")
        return predictions
    
    X_valid = X[valid_mask]
    X_valid = np.where(np.isnan(X_valid), means, X_valid)
    X_valid = (X_valid - means) / stds
    
    try:
        y_proba = source_model['model'].predict_proba(X_valid)[:, 1]
        predictions[valid_mask] = y_proba
        print(f"      ✓ Generated {valid_mask.sum()} predictions")
    except Exception as e:
        print(f"      ✗ Prediction failed: {type(e).__name__}: {str(e)[:100]}")
    
    return predictions


# ============ Ensemble Logic ============

def majority_vote_ensemble(predictions_list: List[np.ndarray], 
                           threshold: float = 0.25) -> np.ndarray:
    """
    Majority vote ensemble from list of prediction arrays.
    
    - Binarize each at threshold
    - Take majority vote
    - Ties default to 1 (water)
    - Returns binary predictions
    """
    
    if not predictions_list:
        return np.array([])
    
    # Stack predictions
    stacked = np.stack(predictions_list, axis=-1)  # shape: (pixels, n_models)
    
    # Binarize
    binary = (stacked >= threshold).astype(np.float32)
    
    # Count votes
    votes = np.sum(binary, axis=-1)
    n_models = stacked.shape[-1]
    
    # Majority vote (ties go to 1)
    ensemble_pred = (votes >= (n_models / 2.0)).astype(np.uint8)
    
    return ensemble_pred


# ============ Metrics Calculation ============

def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Compute all 7 metrics."""
    
    metrics = {
        'Accuracy': float(accuracy_score(y_true, y_pred)),
        'Precision': float(precision_score(y_true, y_pred, zero_division=0)),
        'Recall': float(recall_score(y_true, y_pred, zero_division=0)),
        'F1': float(f1_score(y_true, y_pred, zero_division=0)),
        'MCC': float(matthews_corrcoef(y_true, y_pred)) if len(np.unique(y_true)) > 1 else 0.0,
    }
    
    try:
        metrics['ROC_AUC'] = float(roc_auc_score(y_true, y_pred))
    except:
        metrics['ROC_AUC'] = np.nan
    
    try:
        metrics['PR_AUC'] = float(average_precision_score(y_true, y_pred))
    except:
        metrics['PR_AUC'] = np.nan
    
    return metrics


# ============ Main Ensemble Generation ============

def generate_ensemble_predictions_for_location(
    location: Dict,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs"
) -> pd.DataFrame:
    """Generate ensemble predictions and metrics for one location."""
    
    lat = location['latitude']
    lon = location['longitude']
    name = location['name']
    
    print(f"\n{'='*60}")
    print(f"Processing {name} (Lat_{lat}_Lon_{lon})")
    print(f"{'='*60}")
    
    # Load test data
    print("  Loading test data...")
    df_test = load_test_data_for_location(lat, lon, bucket, base_prefix)
    
    if df_test.empty:
        print("  ✗ No test data found")
        return pd.DataFrame()
    
    dates = sorted(df_test['Date_Normalized'].dropna().unique())
    print(f"  Found {len(dates)} unique dates")
    
    # Load architecture predictions
    print("  Loading architecture predictions...")
    
    arch_preds = {}
    for arch_name, subfolder in ARCHITECTURE_SUBFOLDERS.items():
        tiff_path = (
            f"s3://{bucket}/{base_prefix}/"
            f"trainingFolderFor_Lat_{lat}_Lon_{lon}/"
            f"{subfolder}/"
            f"{subfolder.replace('assess', '').replace('Model', '').lower()}_predictions_Lat_{lat}_Lon_{lon}.tif"
        )
        
        print(f"    Trying {arch_name}:")
        print(f"      Path: .../{subfolder}/{os.path.basename(tiff_path)}")
        
        preds, _ = load_prediction_tiff(tiff_path)
        if preds:
            arch_preds[arch_name] = preds
            print(f"      ✓ Loaded {len(preds)} dates")
        else:
            print(f"      ✗ Failed")
    
    print(f"    Total loaded: {len(arch_preds)}/{len(ARCHITECTURE_SUBFOLDERS)} architectures")
    
    # Load feature predictions
    print("  Loading feature predictions...")
    
    FEATURE_FILENAME_PREFIXES = {
        "assessXGBoostModel": "xgboost_predictions",
        "assessXGBoostTraditionalModel": "xgboost_predictions",
        "assessXGBoostWeatherModel": "xgboost_predictions",
        "assessXGBoostMultiLearnModel": "xgboostmultilearn_predictions",
        "assessXGBoostNoisyModel": "xgboost_predictions"
    }
    
    feat_preds = {}
    for feat_name, subfolder in FEATURE_SUBFOLDERS.items():
        if feat_name == "XGB - All Features":
            continue
        
        filename_prefix = FEATURE_FILENAME_PREFIXES.get(subfolder, "xgboost_predictions")
        
        tiff_path = (
            f"s3://{bucket}/{base_prefix}/"
            f"trainingFolderFor_Lat_{lat}_Lon_{lon}/"
            f"{subfolder}/"
            f"{filename_prefix}_Lat_{lat}_Lon_{lon}.tif"
        )
        preds, _ = load_prediction_tiff(tiff_path)
        if preds:
            feat_preds[feat_name] = preds
    
    print(f"    Loaded {len(feat_preds)} feature configs")
    
    # Load transfer models and generate predictions
    print("  Loading transfer models...")
    transfer_preds = {}
    source_locs = [loc for loc in basin_locations if loc != location]
    
    for source_loc in source_locs:
        source_model = load_model_from_location(
            source_loc['latitude'], source_loc['longitude'], bucket, base_prefix
        )
        
        if source_model:
            transfer_preds[source_loc['name']] = {
                'model': source_model['model'],
                'feature_cols': source_model['feature_cols'],
                'means': source_model['means'],
                'stds': source_model['stds']
            }
    
    print(f"    Loaded {len(transfer_preds)} transfer models")
    
    # Process each date
    all_metrics = []
    
    for date_str in dates:
        date_data = df_test[df_test['Date_Normalized'] == date_str].copy()
        
        if len(date_data) < 100:
            continue
        
        y_true = date_data['Target'].to_numpy(dtype=np.uint8)
        pixel_rows = date_data['PixelRow'].to_numpy()
        pixel_cols = date_data['PixelCol'].to_numpy()
        
        # 1. XGBoost Baseline
        xgb_preds_2d = arch_preds.get('XGBoost', {}).get(date_str)
        
        if xgb_preds_2d is not None:
            xgb_1d = xgb_preds_2d[pixel_rows, pixel_cols]
            xgb_binary = (xgb_1d >= 0.5).astype(np.uint8)
            
            valid = ~np.isnan(xgb_1d)
            if valid.sum() > 100:
                metrics = compute_metrics(y_true[valid], xgb_binary[valid])
                metrics['date'] = date_str
                metrics['ensemble'] = 'XGBoost Baseline'
                all_metrics.append(metrics)
        
        # 2. Architecture Ensemble
        arch_preds_1d = []
        for arch_name, preds_dict in arch_preds.items():
            if date_str in preds_dict:
                pred_2d = preds_dict[date_str]
                pred_1d = pred_2d[pixel_rows, pixel_cols]
                valid = ~np.isnan(pred_1d)
                if valid.sum() > 0:
                    arch_preds_1d.append(pred_1d)
        
        if len(arch_preds_1d) > 0:
            arch_ensemble = majority_vote_ensemble(arch_preds_1d)
            valid = np.all([~np.isnan(p) for p in arch_preds_1d], axis=0)
            
            total_pixels = len(valid)
            valid_pixels = valid.sum()
            retention_pct = (valid_pixels / total_pixels * 100) if total_pixels > 0 else 0
            
            print(f"    [Architecture] Date {date_str}: {len(arch_preds_1d)} models, "
                  f"{valid_pixels}/{total_pixels} pixels retained ({retention_pct:.1f}%)")
            
            if retention_pct < 50:
                print(f"      ⚠️  Warning: Strict mask removed >{100-retention_pct:.0f}% of pixels")
            
            if valid.sum() > 100:
                metrics = compute_metrics(y_true[valid], arch_ensemble[valid])
                metrics['date'] = date_str
                metrics['ensemble'] = 'Architecture Ensemble'
                all_metrics.append(metrics)
        
        # 3. Feature Ensemble
        feat_preds_1d = []
        if xgb_preds_2d is not None:
            feat_preds_1d.append(xgb_preds_2d[pixel_rows, pixel_cols])
        
        for feat_name, preds_dict in feat_preds.items():
            if date_str in preds_dict:
                pred_2d = preds_dict[date_str]
                pred_1d = pred_2d[pixel_rows, pixel_cols]
                valid = ~np.isnan(pred_1d)
                if valid.sum() > 0:
                    feat_preds_1d.append(pred_1d)
        
        if len(feat_preds_1d) > 0:
            feat_ensemble = majority_vote_ensemble(feat_preds_1d)
            valid = np.all([~np.isnan(p) for p in feat_preds_1d], axis=0)
            
            total_pixels = len(valid)
            valid_pixels = valid.sum()
            retention_pct = (valid_pixels / total_pixels * 100) if total_pixels > 0 else 0
            
            print(f"    [Feature] Date {date_str}: {len(feat_preds_1d)} models, "
                  f"{valid_pixels}/{total_pixels} pixels retained ({retention_pct:.1f}%)")
            
            if retention_pct < 50:
                print(f"      ⚠️  Warning: Strict mask removed >{100-retention_pct:.0f}% of pixels")
            
            if valid.sum() > 100:
                metrics = compute_metrics(y_true[valid], feat_ensemble[valid])
                metrics['date'] = date_str
                metrics['ensemble'] = 'Feature Ensemble'
                all_metrics.append(metrics)
       
        # 4. Transfer Ensemble
        transfer_preds_1d = []
        for source_name, source_info in transfer_preds.items():
            pred_1d = predict_with_transfer_model(
                source_info, date_data,
                source_info['feature_cols'],
                source_info['means'],
                source_info['stds']
            )
            valid = ~np.isnan(pred_1d)
            if valid.sum() > 0:
                transfer_preds_1d.append(pred_1d)
        
        if len(transfer_preds_1d) > 0:
            transfer_ensemble = majority_vote_ensemble(transfer_preds_1d)
            valid = np.all([~np.isnan(p) for p in transfer_preds_1d], axis=0)
            
            total_pixels = len(valid)
            valid_pixels = valid.sum()
            retention_pct = (valid_pixels / total_pixels * 100) if total_pixels > 0 else 0
            
            print(f"    [Transfer] Date {date_str}: {len(transfer_preds_1d)} models, "
                  f"{valid_pixels}/{total_pixels} pixels retained ({retention_pct:.1f}%)")
            
            if retention_pct < 50:
                print(f"      ⚠️  Warning: Strict mask removed >{100-retention_pct:.0f}% of pixels")
            
            if valid.sum() > 100:
                metrics = compute_metrics(y_true[valid], transfer_ensemble[valid])
                metrics['date'] = date_str
                metrics['ensemble'] = 'Transfer Ensemble'
                all_metrics.append(metrics)

        if len(transfer_preds_1d) > 1:
            transfer_binary = [(p >= 0.5).astype(int) for p in transfer_preds_1d]
            agreements = []
            for i in range(len(transfer_binary)):
                for j in range(i+1, len(transfer_binary)):
                    agree = (transfer_binary[i] == transfer_binary[j]).mean()
                    agreements.append(agree)
            print(f"    Transfer model avg agreement: {np.mean(agreements):.3f}")
                
        # 5. Kitchen Sink
        kitchen_sink_preds = []
        kitchen_sink_preds.extend(arch_preds_1d)
        
        for feat_name, preds_dict in feat_preds.items():
            if date_str in preds_dict:
                pred_2d = preds_dict[date_str]
                pred_1d = pred_2d[pixel_rows, pixel_cols]
                valid = ~np.isnan(pred_1d)
                if valid.sum() > 0:
                    kitchen_sink_preds.append(pred_1d)
        
        kitchen_sink_preds.extend(transfer_preds_1d)
        
        if len(kitchen_sink_preds) > 0:
            kitchen_ensemble = majority_vote_ensemble(kitchen_sink_preds)
            valid = np.all([~np.isnan(p) for p in kitchen_sink_preds], axis=0)
            
            total_pixels = len(valid)
            valid_pixels = valid.sum()
            retention_pct = (valid_pixels / total_pixels * 100) if total_pixels > 0 else 0
            
            print(f"    [Kitchen Sink] Date {date_str}: {len(kitchen_sink_preds)} models, "
                  f"{valid_pixels}/{total_pixels} pixels retained ({retention_pct:.1f}%)")
            
            if retention_pct < 50:
                print(f"      ⚠️  Warning: Strict mask removed >{100-retention_pct:.0f}% of pixels")
            
            if valid.sum() > 100:
                metrics = compute_metrics(y_true[valid], kitchen_ensemble[valid])
                metrics['date'] = date_str
                metrics['ensemble'] = 'Kitchen Sink'
                all_metrics.append(metrics)
    
    if all_metrics:
        df_metrics = pd.DataFrame(all_metrics)
        df_metrics['location'] = name
        df_metrics['latitude'] = lat
        df_metrics['longitude'] = lon
        print(f"  ✓ Generated {len(df_metrics)} ensemble comparisons")
        return df_metrics
    else:
        print(f"  ✗ No metrics generated")
        return pd.DataFrame()


# ============ Visualization ============

def _get_box_patches(ax):
    boxes = [a for a in ax.artists if isinstance(a, PathPatch)]
    if not boxes:
        boxes = [p for p in ax.patches if isinstance(p, PathPatch)]
    n_cats = len(ax.get_xticklabels())
    return boxes[:n_cats]


def _highlight_boxes_vs_baseline(ax, df_plot: pd.DataFrame, metric: str, alpha: float = ALPHA):
    df_plot = df_plot[df_plot["ensemble"].isin(ENSEMBLE_ORDER)].dropna(subset=[metric])
    
    baseline_vals = df_plot.loc[df_plot["ensemble"] == "XGBoost Baseline", metric].dropna()
    if baseline_vals.empty:
        return
    
    xticklabels = [t.get_text() for t in ax.get_xticklabels()]
    box_patches = _get_box_patches(ax)
    
    if len(box_patches) != len(xticklabels):
        return
    
    model_to_box = {label: box for label, box in zip(xticklabels, box_patches)}
    
    for ensemble in ENSEMBLE_ORDER:
        if ensemble == "XGBoost Baseline":
            continue
        if ensemble not in df_plot["ensemble"].unique():
            continue
        if ensemble not in model_to_box:
            continue
        
        ensemble_vals = df_plot.loc[df_plot["ensemble"] == ensemble, metric].dropna()
        if ensemble_vals.empty:
            continue
        
        box = model_to_box[ensemble]
        
        try:
            _, p_greater = mannwhitneyu(ensemble_vals, baseline_vals, alternative="greater")
        except ValueError:
            p_greater = 1.0
        
        try:
            _, p_less = mannwhitneyu(ensemble_vals, baseline_vals, alternative="less")
        except ValueError:
            p_less = 1.0
        
        if p_greater < alpha:
            box.set_edgecolor("#FFD700")
            box.set_linewidth(3.0)
        elif p_less < alpha:
            box.set_edgecolor("#DC143C")
            box.set_linewidth(3.0)


def plot_all_metrics_summary(metrics_df: pd.DataFrame, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for idx, metric in enumerate(METRICS):
        ax = axes[idx]
        
        df_plot = metrics_df[["ensemble", metric]].dropna()
        df_plot = df_plot[df_plot["ensemble"].isin(ENSEMBLE_ORDER)]
        
        if df_plot.empty:
            ax.set_title(f"{metric} (no data)", fontsize=34)  # ✅ 70% bigger
            ax.axis('off')
            continue
        
        sns.boxplot(
            data=df_plot,
            x="ensemble",
            y=metric,
            order=ENSEMBLE_ORDER,
            showfliers=False,
            ax=ax,
        )
        
        ax.set_ylim(0.0, 1.0)
        # ✅ Remove "Ensemble Method" label
        ax.set_xlabel("", fontsize=29)  # ✅ 70% bigger
        # ✅ Remove metric name from y-axis
        ax.set_ylabel("", fontsize=34, fontweight='bold')  # ✅ 70% bigger
        ax.set_title(f"{metric}", fontsize=34, fontweight='bold')  # ✅ 70% bigger
        ax.tick_params(axis='x', rotation=45, labelsize=26)  # ✅ 70% bigger
        ax.tick_params(axis='y', labelsize=26)  # ✅ 70% bigger
        ax.grid(True, alpha=0.3, axis='y')
        
        _highlight_boxes_vs_baseline(ax, df_plot, metric)
    
    axes[7].axis('off')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', linewidth=1, label='No significant difference'),
        Patch(facecolor='lightblue', edgecolor='#FFD700', linewidth=3, label='Significantly better (p<0.05)'),
        Patch(facecolor='lightblue', edgecolor='#DC143C', linewidth=3, label='Significantly worse (p<0.05)'),
    ]
    axes[7].legend(handles=legend_elements, loc='center', 
                   fontsize=34,  # ✅ 70% bigger
                   frameon=True,
                   title='vs. XGBoost Baseline', 
                   title_fontsize=34)  # ✅ 70% bigger
    
    fig.suptitle('Model Ensemble Comparison\nAggregated Performance Across All Locations',
                fontsize=46, fontweight='bold', y=0.98)  # ✅ 70% bigger
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    out_path = os.path.join(output_dir, "all_metrics_ensemble_summary.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Saved summary figure to: {out_path}")


def plot_per_metric_all_locations_ensemble(
    metrics_df: pd.DataFrame,
    output_dir: str,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs"
):
    """
    For each metric, create a PNG with 2x4 subplots (one per location).
    Each subplot shows 5 boxplots: XGBoost Baseline + 4 ensemble methods.
    Uses gold/red highlighting for significance vs baseline.
    """
    
    n_locs = len(basin_locations)
    ncols = 4
    nrows = 2  # Fixed 2x4 grid for 8 locations
    
    for metric in METRICS:
        print(f"  Generating per-location figure for {metric}...")
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(20, 10))
        axes = axes.flatten()
        
        for idx, location in enumerate(basin_locations):
            ax = axes[idx]
            loc_name = location['name']
            
            # Filter data for this location and metric
            df_loc = metrics_df[metrics_df['location'] == loc_name].copy()
            
            if df_loc.empty or metric not in df_loc.columns:
                ax.set_title(f"{loc_name} (no data)", fontsize=32)  # ✅ 70% bigger
                ax.axis('off')
                continue
            
            df_plot = df_loc[["ensemble", metric]].dropna()
            df_plot = df_plot[df_plot["ensemble"].isin(ENSEMBLE_ORDER)]
            
            if df_plot.empty:
                ax.set_title(f"{loc_name} (no data)", fontsize=32)  # ✅ 70% bigger
                ax.axis('off')
                continue
            
            # Create boxplot
            sns.boxplot(
                data=df_plot,
                x="ensemble",
                y=metric,
                order=ENSEMBLE_ORDER,
                showfliers=False,
                ax=ax,
            )
            
            ax.set_ylim(0.0, 1.0)
            # ✅ Remove "Ensemble Method" label
            ax.set_xlabel("", fontsize=29)  # ✅ 70% bigger
            # ✅ Remove metric name from y-axis
            ax.set_ylabel("", fontsize=34)  # ✅ 70% bigger
            ax.set_title(loc_name, fontsize=32, fontweight='bold')  # ✅ 70% bigger
            ax.tick_params(axis='x', rotation=45, labelsize=24)  # ✅ 70% bigger
            ax.tick_params(axis='y', labelsize=26)  # ✅ 70% bigger
            ax.grid(True, alpha=0.3, axis='y')
            
            # Apply gold/red highlighting vs baseline
            _highlight_boxes_vs_baseline(ax, df_plot, metric, alpha=ALPHA)
        
        # Overall title
        fig.suptitle(
            f'Ensemble Performance by Location: {metric}\n'
            f'Gold outline = significantly better than XGBoost Baseline (p<0.05) | '
            f'Red outline = significantly worse (p<0.05)',
            fontsize=41, fontweight='bold', y=0.98  # ✅ 70% bigger
        )
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Save figure
        output_path = f"s3://{bucket}/{base_prefix}/{OUTPUT_DIR}/per_location_{metric}_ensemble_comparison.png"
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            temp_png = tmp.name
        
        plt.savefig(temp_png, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Upload to S3
        s3_client = boto3.client('s3')
        output_key = output_path.replace(f"s3://{bucket}/", "")
        s3_client.upload_file(temp_png, bucket, output_key)
        os.remove(temp_png)
        
        print(f"    ✓ Saved: {output_path}")

# ============ Main Execution ============

def run_ensemble_comparison(
    basin_locations: List[Dict],
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs"
):
    """Run complete ensemble comparison analysis."""
    
    print("=" * 70)
    print("ENSEMBLE COMPARISON ANALYSIS")
    print("=" * 70)
    
    all_metrics = []
    
    for location in basin_locations:
        metrics_df = generate_ensemble_predictions_for_location(
            location, bucket, base_prefix
        )
        
        if not metrics_df.empty:
            all_metrics.append(metrics_df)
    
    if not all_metrics:
        print("\n❌ No metrics generated")
        return
    
    combined_df = pd.concat(all_metrics, ignore_index=True)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total comparisons: {len(combined_df)}")
    print(f"Locations: {combined_df['location'].nunique()}")
    print(f"\nBy ensemble:")
    print(combined_df.groupby('ensemble').size())
    
    # Save CSV
    output_dir = f"s3://{bucket}/{base_prefix}/{OUTPUT_DIR}"
    csv_path = f"{output_dir}/ensemble_per_date_metrics.csv"
    
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
        temp_csv = tmp.name
    combined_df.to_csv(temp_csv, index=False)
    
    s3_client = boto3.client('s3')
    csv_key = csv_path.replace(f"s3://{bucket}/", "")
    s3_client.upload_file(temp_csv, bucket, csv_key)
    os.remove(temp_csv)
    
    print(f"\n✓ Saved CSV: {csv_path}")
    
    # Generate summary figure (all metrics, aggregated)
    plot_all_metrics_summary(combined_df, OUTPUT_DIR)
    
    # Generate per-location figures (one per metric)
    print(f"\n{'='*70}")
    print("GENERATING PER-LOCATION FIGURES")
    print(f"{'='*70}")
    plot_per_metric_all_locations_ensemble(combined_df, OUTPUT_DIR, bucket, base_prefix)
    
    print(f"\n✅ Ensemble comparison complete!")



if __name__ == "__main__":
    run_ensemble_comparison(basin_locations)