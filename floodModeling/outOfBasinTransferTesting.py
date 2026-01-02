"""
outOfBasinTransferTesting.py - H2 Transfer Testing

Tests model transferability across basins without importing xgboostFloods.py
(to avoid statsmodels compatibility issues).

Basin locations are passed as parameters rather than hardcoded.
"""

import os
import pickle
import tempfile
from typing import Dict, List
import numpy as np
import pandas as pd
import boto3

import pyarrow as pa
import pyarrow.dataset as ds

from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    cohen_kappa_score, matthews_corrcoef, roc_auc_score,
    average_precision_score, confusion_matrix
)


# ============ Helper Functions (copied to avoid statsmodels import) ============

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


def make_scanner(dataset, columns=None, filter=None, batch_size=None):
    kwargs = {}
    if columns is not None:
        kwargs["columns"] = columns
    if filter is not None:
        kwargs["filter"] = filter
    if batch_size is not None:
        try:
            kwargs["batch_size"] = batch_size
        except:
            pass
    return dataset.scanner(**kwargs)


def scanner_batches(scanner):
    if hasattr(scanner, "to_batches"):
        return scanner.to_batches()
    elif hasattr(scanner, "scan_batches"):
        return scanner.scan_batches()
    else:
        tbl = scanner.to_table()
        return tbl.to_batches(max_chunksize=65536)


def arrow_batch_to_numpy(batch, columns: List[str]) -> np.ndarray:
    """Efficiently convert Arrow batch to NumPy array."""
    if batch.num_rows == 0:
        return np.empty((0, len(columns)), dtype=np.float32)
    
    arrays = []
    for i, col_name in enumerate(columns):
        arr = batch.column(i).to_numpy(zero_copy_only=False)
        if arr.dtype != np.float32:
            arr = arr.astype(np.float32, copy=False)
        arrays.append(arr)
    
    return np.column_stack(arrays) if arrays else np.empty((batch.num_rows, 0), dtype=np.float32)


def compute_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics."""
    y_pred = (y_proba >= threshold).astype(np.uint8)
    
    prevalence = float(np.mean(y_true))
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    kappa = cohen_kappa_score(y_true, y_pred)
    
    if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1:
        mcc = matthews_corrcoef(y_true, y_pred)
    else:
        mcc = 0.0
    
    try:
        roc_auc = roc_auc_score(y_true, y_proba)
    except:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_proba)
    except:
        pr_auc = np.nan
    
    if len(np.unique(y_true)) > 1:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    else:
        if y_true[0] == 0:
            tn, fp, fn, tp = len(y_true), 0, 0, 0
        else:
            tn, fp, fn, tp = 0, 0, 0, len(y_true)
    
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    tnr = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    
    return {
        "threshold": float(threshold),
        "prevalence": prevalence,
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "kappa": kappa,
        "mcc": mcc,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "tpr": tpr,
        "tnr": tnr,
    }


def get_unique_watersheds(dataset) -> List[str]:
    """Extract unique watershed IDs from dataset."""
    uniques = set()
    scanner = make_scanner(dataset, columns=["Watershed_Loc"])
    for batch in scanner_batches(scanner):
        if batch.num_columns == 0 or batch.num_rows == 0:
            continue
        s = pa.Table.from_batches([batch]).column(0).to_pandas()
        uniques.update(s.dropna().unique().tolist())
    return sorted(uniques)


# ============ Main Transfer Testing Functions ============

def load_model_from_location(latitude: str, longitude: str, 
                             bucket: str = "climate-ai-data-science-datasets",
                             base_prefix: str = "arrakis-data/floodOutputs") -> Dict:
    """Load pre-trained XGBoost model from a specific basin location."""
    model_path = (
        f"s3://{bucket}/{base_prefix}/"
        f"trainingFolderFor_Lat_{latitude}_Lon_{longitude}/"
        f"modelXGBoost/xgb_model.pkl"
    )
    
    print(f"  Loading model from Lat_{latitude}_Lon_{longitude}...")
    
    s3_client = boto3.client('s3')
    model_key = model_path.replace(f"s3://{bucket}/", "")
    
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        s3_client.download_file(bucket, model_key, temp_path)
        
        with open(temp_path, 'rb') as f:
            results = pickle.load(f)
        
        os.remove(temp_path)
        print(f"    ✓ Loaded successfully ({len(results['feature_cols'])} features)")
        return results
        
    except Exception as e:
        print(f"    ✗ Failed to load: {e}")
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return None


def validate_and_align_features(df: pd.DataFrame, feature_cols: List[str]) -> np.ndarray:
    """Align features from target data to match source model's feature order."""
    missing_features = [f for f in feature_cols if f not in df.columns]
    
    if missing_features:
        print(f"    Warning: {len(missing_features)} features missing, filling with NaN")
        for feat in missing_features:
            df[feat] = np.nan
    
    X = df[feature_cols].to_numpy(dtype=np.float32, copy=False)
    return X


def evaluate_model_on_target_location(
    source_model_results: Dict,
    target_latitude: str,
    target_longitude: str,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs"
) -> pd.DataFrame:
    """Evaluate a source basin's model on target basin's data."""
    
    model = source_model_results['model']
    threshold = source_model_results['threshold']
    feature_cols = source_model_results['feature_cols']
    means = source_model_results['means']
    stds = source_model_results['stds']
    
    target_parquet = (
        f"s3://{bucket}/{base_prefix}/"
        f"trainingFolderFor_Lat_{target_latitude}_Lon_{target_longitude}/"
        f"postProcessedTable/flood_training.parquet"
    )
    
    print(f"    Loading target data from Lat_{target_latitude}_Lon_{target_longitude}...")
    
    try:
        dataset = open_parquet_dataset(target_parquet)
    except Exception as e:
        print(f"    ✗ Cannot open target parquet: {e}")
        return pd.DataFrame()
    
    target_watershed = f"Lat_{target_latitude}_Lon_{target_longitude}"
    
    filt = ds.field("Watershed_Loc") == target_watershed
    scanner = make_scanner(dataset, filter=filt)
    
    pixel_data = []
    for batch in scanner_batches(scanner):
        if batch.num_rows > 0:
            df = pa.Table.from_batches([batch]).to_pandas()
            pixel_data.append(df)
    
    if not pixel_data:
        print(f"    ✗ No data found for target watershed")
        return pd.DataFrame()
    
    df_all = pd.concat(pixel_data, ignore_index=True)
    
    def normalize_date(d):
        if pd.isna(d):
            return None
        if isinstance(d, str):
            return d.replace('-', '').replace('/', '')[:8]
        if isinstance(d, (pd.Timestamp, np.datetime64)):
            return pd.Timestamp(d).strftime('%Y%m%d')
        return str(d)
    
    df_all['Date_Normalized'] = df_all['Date'].apply(normalize_date)
    dates = sorted(df_all['Date_Normalized'].dropna().unique())
    
    print(f"    Processing {len(dates)} dates...")
    
    per_date_metrics = []
    
    for date_str in dates:
        date_data = df_all[df_all['Date_Normalized'] == date_str].copy()
        
        X = validate_and_align_features(date_data, feature_cols)
        y_true = date_data['Target'].to_numpy(dtype=np.uint8, copy=False)
        
        non_nan_counts = np.sum(~np.isnan(X), axis=1)
        feature_threshold = int(0.70 * len(feature_cols))
        valid_mask = non_nan_counts >= feature_threshold
        
        if valid_mask.sum() < 100:
            continue
        
        X_valid = X[valid_mask]
        y_true_valid = y_true[valid_mask]
        
        X_valid = np.where(np.isnan(X_valid), means, X_valid)
        X_valid = (X_valid - means) / stds
        
        try:
            y_proba = model.predict_proba(X_valid)[:, 1]
            
            metrics = compute_metrics(y_true_valid, y_proba, threshold)
            metrics['date'] = date_str
            metrics['valid_pixel_count'] = int(valid_mask.sum())
            metrics['total_pixel_count'] = len(date_data)
            
            per_date_metrics.append(metrics)
            
        except Exception as e:
            print(f"      Warning: Failed for {date_str}: {e}")
            continue
    
    if per_date_metrics:
        print(f"    ✓ Generated metrics for {len(per_date_metrics)} dates")
        return pd.DataFrame(per_date_metrics)
    else:
        print(f"    ✗ No valid predictions generated")
        return pd.DataFrame()


def run_out_of_basin_transfer_test(
    target_latitude: str,
    target_longitude: str,
    basin_locations: List[Dict[str, str]],
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs",
    model_type: str = "xgboost"
):
    """Run complete out-of-basin transfer test for one target location.
    
    Args:
        target_latitude: Target basin latitude
        target_longitude: Target basin longitude
        basin_locations: List of dicts with 'latitude', 'longitude', 'name' for each basin
        bucket: S3 bucket
        base_prefix: S3 prefix
        model_type: Model type (currently only 'xgboost')
    """
    
    print("=" * 70)
    print("OUT-OF-BASIN TRANSFER TEST FOR H2")
    print(f"Target Location: Lat_{target_latitude}_Lon_{target_longitude}")
    print("=" * 70)
    
    target_loc = None
    for loc in basin_locations:
        if loc['latitude'] == target_latitude and loc['longitude'] == target_longitude:
            target_loc = loc
            break
    
    if target_loc is None:
        raise ValueError(f"Target location Lat_{target_latitude}_Lon_{target_longitude} not found in basin_locations list")
    
    source_locations = [loc for loc in basin_locations if loc != target_loc]
    
    print(f"\nTesting {len(source_locations)} source models on target location")
    print(f"Models will be saved to: outOfBasinsTests/\n")
    
    output_folder = (
        f"s3://{bucket}/{base_prefix}/"
        f"trainingFolderFor_Lat_{target_latitude}_Lon_{target_longitude}/"
        f"outOfBasinsTests"
    )
    
    s3_client = boto3.client('s3')
    results_summary = []
    
    for idx, source_loc in enumerate(source_locations, 1):
        print(f"\n{'='*60}")
        print(f"Source Model {idx}/{len(source_locations)}: {source_loc['name']}")
        print(f"  Lat_{source_loc['latitude']}_Lon_{source_loc['longitude']}")
        print(f"{'='*60}")
        
        source_model = load_model_from_location(
            source_loc['latitude'], 
            source_loc['longitude'],
            bucket,
            base_prefix
        )
        
        if source_model is None:
            print(f"  ✗ Skipping due to load failure\n")
            continue
        
        metrics_df = evaluate_model_on_target_location(
            source_model,
            target_latitude,
            target_longitude,
            bucket,
            base_prefix
        )
        
        if metrics_df.empty:
            print(f"  ✗ No metrics generated\n")
            continue
        
        source_name = f"Lat_{source_loc['latitude']}_Lon_{source_loc['longitude']}"
        target_name = f"Lat_{target_latitude}_Lon_{target_longitude}"
        
        csv_filename = (
            f"{model_type}_{source_name}_per_day_metrics_{target_name}.csv"
        )
        csv_path = f"{output_folder}/{csv_filename}"
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
            temp_csv = tmp.name
        
        metrics_df.to_csv(temp_csv, index=False)
        
        csv_key = csv_path.replace(f"s3://{bucket}/", "")
        s3_client.upload_file(temp_csv, bucket, csv_key)
        os.remove(temp_csv)
        
        print(f"  ✓ Saved: {csv_filename}")
        
        summary = {
            'source_location': source_loc['name'],
            'source_lat': source_loc['latitude'],
            'source_lon': source_loc['longitude'],
            'target_lat': target_latitude,
            'target_lon': target_longitude,
            'n_dates': len(metrics_df),
            'mean_mcc': metrics_df['mcc'].mean(),
            'mean_f1': metrics_df['f1'].mean(),
            'mean_accuracy': metrics_df['accuracy'].mean(),
            'mean_precision': metrics_df['precision'].mean(),
            'mean_recall': metrics_df['recall'].mean(),
            'csv_path': csv_path
        }
        results_summary.append(summary)
    
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        
        summary_csv = f"{output_folder}/transfer_test_summary_{target_name}.csv"
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
            temp_csv = tmp.name
        
        summary_df.to_csv(temp_csv, index=False)
        
        summary_key = summary_csv.replace(f"s3://{bucket}/", "")
        s3_client.upload_file(temp_csv, bucket, summary_key)
        os.remove(temp_csv)
        
        print(f"\n{'='*60}")
        print("SUMMARY OF TRANSFER TEST")
        print(f"{'='*60}")
        print(summary_df[['source_location', 'n_dates', 'mean_mcc', 'mean_f1']].to_string(index=False))
        print(f"\nFull summary saved to: transfer_test_summary_{target_name}.csv")
        
        print(f"\n{'='*60}")
        print("HYPOTHESIS H2 ANALYSIS")
        print(f"{'='*60}")
        print(f"Within-basin models typically achieve MCC ≈ 0.62-0.72")
        print(f"Transfer models achieved:")
        print(f"  Mean MCC: {summary_df['mean_mcc'].mean():.3f} ± {summary_df['mean_mcc'].std():.3f}")
        print(f"  Mean F1:  {summary_df['mean_f1'].mean():.3f} ± {summary_df['mean_f1'].std():.3f}")
        print(f"\nH2 is SUPPORTED if transfer performance is substantially worse,")
        print(f"demonstrating that basin-specific training is necessary.")
    
    print(f"\n✅ Out-of-basin transfer test complete!")
    print(f"All outputs saved to: {output_folder}")


if __name__ == "__main__":
    # Example usage
    basin_locations = [
        {"latitude": "10.7875",  "longitude": "0.7479", "name": "EMSR470"},
        {"latitude": "48.1747",  "longitude": "23.3978", "name": "EMSR444"},
        {"latitude": "41.9275",  "longitude": "-1.3792", "name": "EMSR279"},
    ]
    
    target_lat = "10.7875"
    target_lon = "0.7479"
    
    run_out_of_basin_transfer_test(
        target_latitude=target_lat,
        target_longitude=target_lon,
        basin_locations=basin_locations
    )