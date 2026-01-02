###############################################################
###############################################################
###############################################################
### STURM-Floods validation
###############################################################
###############################################################
###############################################################


"""
sturmPerformanceAggregator.py - Aggregate STURM Performance Analysis

Analyzes XGBoost model performance against STURM ground truth across multiple locations.
Creates:
1. Per-location boxplots (8 locations x 7 metrics)
2. Aggregate boxplots (7 metrics across all locations)
3. CSV exports of all metrics
"""

import os
import re
import tempfile
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from osgeo import gdal
import boto3
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    matthews_corrcoef, average_precision_score, roc_auc_score
)

gdal.UseExceptions()


# ======================== GDAL Configuration ========================

def configure_gdal():
    """Configure GDAL environment."""
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        os.environ["PROJ_LIB"] = proj_dir
        gdal.SetConfigOption("PROJ_LIB", proj_dir)
    except:
        pass

    gdal.SetConfigOption("PROJ_NETWORK", "ON")
    gdal.SetConfigOption("GDAL_CACHEMAX", "512")
    gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,tiff")
    gdal.SetConfigOption("AWS_NO_SIGN_REQUEST", "NO")
    gdal.PushErrorHandler('CPLQuietErrorHandler')

configure_gdal()


# ======================== Discovery Functions ========================

def discover_confusion_matrix_tiffs(location: Dict) -> List[str]:
    """
    Discover stitched confusion matrix TIFFs for a location.
    Only finds files of the form: confusion_matrix_S{1-or-2}_{date}.tif
    Excludes per-watershed files: confusion_matrix_S{1-or-2}_{date}_Lat_X_Lon_Y.tif
    """
    
    lat = location['latitude']
    lon = location['longitude']
    emsr = location['name']
    
    sturm_folder = f"s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/trainingFolderFor_Lat_{lat}_Lon_{lon}/assessXGBoostModel/sturmComparison"
    
    s3_client = boto3.client('s3')
    bucket, prefix = sturm_folder.replace('s3://', '').split('/', 1)
    
    cm_tiffs = []
    
    try:
        response = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix + '/')
        
        if 'Contents' not in response:
            return []
        
        for obj in response['Contents']:
            key = obj['Key']
            filename = os.path.basename(key)
            
            # Match: confusion_matrix_S1_2020-05-15.tif or confusion_matrix_S2_2020-05-15.tif
            # Exclude: confusion_matrix_S1_2020-05-15_Lat_X_Lon_Y.tif
            if filename.startswith('confusion_matrix_S') and filename.endswith('.tif'):
                # Check if it contains "_Lat_" (which would indicate per-watershed file)
                if '_Lat_' not in filename:
                    cm_tiffs.append(f"s3://{bucket}/{key}")
        
    except Exception as e:
        print(f"  ⚠️ Error discovering TIFFs for {emsr}: {e}")
    
    return cm_tiffs


# ======================== Metrics Calculation ========================

def calculate_metrics_from_cm_tiff(cm_tiff_uri: str) -> Dict[str, float]:
    """
    Calculate all 7 metrics from confusion matrix TIFF.
    CM values: 0=TN, 1=TP, 2=FN, 3=FP, 255=NoData
    """
    
    vsi_path = cm_tiff_uri.replace('s3://', '/vsis3/')
    
    try:
        ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
        if ds is None:
            return None
        
        cm_array = ds.GetRasterBand(1).ReadAsArray()
        ds = None
        
        # Extract counts
        valid_mask = cm_array != 255
        
        if valid_mask.sum() == 0:
            return None
        
        tn = int((cm_array == 0).sum())
        tp = int((cm_array == 1).sum())
        fn = int((cm_array == 2).sum())
        fp = int((cm_array == 3).sum())
        
        total = tn + tp + fn + fp
        
        if total == 0:
            return None
        
        # Reconstruct binary arrays for metric calculation
        y_true = []
        y_pred = []
        
        # TN: true=0, pred=0
        y_true.extend([0] * tn)
        y_pred.extend([0] * tn)
        
        # TP: true=1, pred=1
        y_true.extend([1] * tp)
        y_pred.extend([1] * tp)
        
        # FN: true=1, pred=0
        y_true.extend([1] * fn)
        y_pred.extend([0] * fn)
        
        # FP: true=0, pred=1
        y_true.extend([0] * fp)
        y_pred.extend([1] * fp)
        
        y_true = np.array(y_true)
        y_pred = np.array(y_pred)
        
        # Calculate metrics
        metrics = {
            'accuracy': float(accuracy_score(y_true, y_pred)),
            'precision': float(precision_score(y_true, y_pred, zero_division=0)),
            'recall': float(recall_score(y_true, y_pred, zero_division=0)),
            'f1': float(f1_score(y_true, y_pred, zero_division=0)),
            'mcc': float(matthews_corrcoef(y_true, y_pred)) if len(np.unique(y_true)) > 1 else 0.0,
        }
        
        try:
            metrics['roc_auc'] = float(roc_auc_score(y_true, y_pred))
        except:
            metrics['roc_auc'] = np.nan
        
        try:
            metrics['pr_auc'] = float(average_precision_score(y_true, y_pred))
        except:
            metrics['pr_auc'] = np.nan
        
        metrics['valid_pixels'] = total
        metrics['tp'] = tp
        metrics['fp'] = fp
        metrics['tn'] = tn
        metrics['fn'] = fn
        
        return metrics
        
    except Exception as e:
        print(f"  ⚠️ Error calculating metrics: {e}")
        return None


# ======================== Visualization Functions ========================

# Location labels mapping
LOCATION_LABELS = {
    "EMSR470": "1 (Togo)",
    "EMSR444": "2 (Ukraine)",
    "EMSR407": "3 (UK)",
    "EMSR424": "4 (Madagascar)",
    "EMSR438": "5 (Uganda)",
    "EMSR292": "6 (Greece)",
    "EMSR629": "7 (Pakistan)",
    "EMSR279": "8 (Spain)",
}




LOCATION_ORDER = ["EMSR470", "EMSR444", "EMSR407", "EMSR424", "EMSR438", "EMSR292", "EMSR629", "EMSR279"]

def create_per_location_boxplots(df_all: pd.DataFrame, output_png: str):
    """
    Create 2x4 grid of boxplots showing 7 metrics across 8 locations.
    One subplot per metric, with 8 boxes per subplot (one per location).
    """
    
    metrics = ['accuracy', 'precision', 'recall', 'f1', 'mcc', 'roc_auc', 'pr_auc']
    metric_labels = {
        'accuracy': 'Accuracy',
        'precision': 'Precision',
        'recall': 'Recall',
        'f1': 'F1 Score',
        'mcc': 'MCC',
        'roc_auc': 'ROC-AUC',
        'pr_auc': 'PR-AUC'
    }
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    # Get unique locations in consistent order
    locations = [loc for loc in LOCATION_ORDER if loc in df_all['location'].unique()]
    
    # Create display labels for locations
    location_display_labels = [LOCATION_LABELS.get(loc, loc) for loc in locations]
    
    # Map locations to display labels in the dataframe
    df_plot_all = df_all.copy()
    df_plot_all['location_label'] = df_plot_all['location'].map(LOCATION_LABELS)
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        # Filter data for this metric
        df_plot = df_plot_all[['location_label', metric]].dropna()
        
        if df_plot.empty:
            ax.set_title(f"{metric_labels[metric]} (no data)", fontsize=20, fontweight='bold')
            ax.axis('off')
            continue
        
        # Create boxplot using seaborn
        sns.boxplot(
            data=df_plot,
            x='location_label',
            y=metric,
            order=location_display_labels,
            showfliers=False,
            ax=ax,
        )
        
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel("", fontsize=17)
        ax.set_ylabel("", fontsize=20)
        ax.set_title(metric_labels[metric], fontsize=20, fontweight='bold')
        ax.tick_params(axis='x', labelsize=15, rotation=45) 
        ax.tick_params(axis='y', labelsize=15)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Hide the 8th subplot (bottom right) and add legend
    axes[7].axis('off')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', linewidth=1, label='IQR'),
    ]
#    axes[7].legend(handles=legend_elements, loc='center', fontsize=20, frameon=True)
    
    # Overall title
    fig.suptitle('Model Performance Across Locations\nXGBoost vs STURM Ground Truth',
                fontsize=27, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    
    plt.savefig(temp_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Upload to S3
    if output_png.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_png.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_png, bucket, key)
    
    os.remove(temp_png)
    print(f"  ✓ Created: {output_png}")


def create_aggregate_boxplots(df_all: pd.DataFrame, output_png: str):
    """
    Create single figure with 7 boxplots showing aggregate performance across all locations.
    X-axis: Metrics, Y-axis: Performance (0-1)
    """
    
    metrics = ['accuracy', 'precision', 'recall', 'f1', 'mcc', 'roc_auc', 'pr_auc']
    metric_labels = ['Accuracy', 'Precision', 'Recall', 'F1', 'MCC', 'ROC-AUC', 'PR-AUC']
    
    # Reshape data for seaborn
    df_melted = df_all[metrics].melt(var_name='metric', value_name='value')
    
    # Map metric names to labels
    metric_map = dict(zip(metrics, metric_labels))
    df_melted['metric'] = df_melted['metric'].map(metric_map)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create boxplot using seaborn
    sns.boxplot(
        data=df_melted,
        x='metric',
        y='value',
        order=metric_labels,
        showfliers=False,
        ax=ax,
    )
    
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("", fontsize=17)
    ax.set_ylabel("", fontsize=20)
    ax.set_title('Aggregate Model Performance Across All Locations\nXGBoost vs STURM Ground Truth',
                fontsize=20, fontweight='bold')
    ax.tick_params(axis='x',labelsize=15, rotation=45)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    
    plt.savefig(temp_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Upload to S3
    if output_png.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_png.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_png, bucket, key)
    
    os.remove(temp_png)
    print(f"  ✓ Created: {output_png}")

# ======================== Main Aggregation Function ========================

def aggregate_sturm_performance(basin_locations: List[Dict],
                                output_folder: str = "s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/model_performance_against_sturm") -> Dict:
    """
    Aggregate STURM performance metrics across multiple locations.
    
    Args:
        basin_locations: List of dicts with 'latitude', 'longitude', 'name' (EMSR code)
        output_folder: S3 folder for outputs
    
    Returns:
        Dictionary with results and output paths
    """
    
    print(f"\n{'='*70}")
    print("STURM PERFORMANCE AGGREGATION")
    print(f"{'='*70}")
    print(f"Locations: {len(basin_locations)}")
    print(f"Output: {output_folder}")
    
    all_metrics = []
    
    # Process each location
    for loc_idx, location in enumerate(basin_locations, 1):
        emsr = location['name']
        lat = location['latitude']
        lon = location['longitude']
        
        print(f"\n[{loc_idx}/{len(basin_locations)}] Processing {emsr} (Lat {lat}, Lon {lon})...")
        
        # Discover confusion matrix TIFFs
        cm_tiffs = discover_confusion_matrix_tiffs(location)
        
        if not cm_tiffs:
            print(f"  ⚠️ No confusion matrix TIFFs found")
            continue
        
        print(f"  Found {len(cm_tiffs)} confusion matrix TIFFs")
        
        # Calculate metrics from each TIFF
        for cm_tiff in cm_tiffs:
            filename = os.path.basename(cm_tiff)
            
            # Extract date and sensor from filename
            # Match YYYYMMDD format (e.g., confusion_matrix_S1_20201014.tif)
            match = re.search(r'confusion_matrix_(S\d)_(\d{8})\.tif', filename)
            
            if not match:
                print(f"  ⚠️ Could not parse filename: {filename}")
                continue
            
            sensor = match.group(1)
            date_yyyymmdd = match.group(2)
            
            # Convert YYYYMMDD to YYYY-MM-DD for consistency
            date = f"{date_yyyymmdd[:4]}-{date_yyyymmdd[4:6]}-{date_yyyymmdd[6:]}"
            
            metrics = calculate_metrics_from_cm_tiff(cm_tiff)
            
            if metrics:
                metrics['location'] = emsr
                metrics['sensor'] = sensor
                metrics['date'] = date
                metrics['latitude'] = lat
                metrics['longitude'] = lon
                
                all_metrics.append(metrics)
                print(f"    ✓ {sensor} {date}: F1={metrics['f1']:.3f}, MCC={metrics['mcc']:.3f}")
    
    if not all_metrics:
        print("\n❌ No metrics calculated - no confusion matrix TIFFs found")
        return {}
    
    # Create DataFrame
    df_all = pd.DataFrame(all_metrics)
    
    print(f"\n{'='*70}")
    print("SUMMARY STATISTICS")
    print(f"{'='*70}")
    print(f"Total comparisons: {len(df_all)}")
    print(f"Locations: {df_all['location'].nunique()}")
    print(f"Unique dates: {df_all['date'].nunique()}")
    print(f"\nBy location:")
    print(df_all.groupby('location').size())
    
    # Save aggregated CSV
    print(f"\n{'='*70}")
    print("SAVING OUTPUTS")
    print(f"{'='*70}")
    
    agg_csv = f"{output_folder}/per_location_metrics.csv"
    
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
        temp_csv = tmp.name
    df_all.to_csv(temp_csv, index=False)
    
    if agg_csv.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = agg_csv.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_csv, bucket, key)
    
    os.remove(temp_csv)
    print(f"  ✓ Created: {agg_csv}")
    
    # Create per-location boxplots
    per_loc_png = f"{output_folder}/per_location_boxplots.png"
    create_per_location_boxplots(df_all, per_loc_png)
    
    # Create aggregate boxplots
    agg_png = f"{output_folder}/aggregate_boxplots.png"
    create_aggregate_boxplots(df_all, agg_png)
    
    # Calculate summary statistics
    print(f"\n{'='*70}")
    print("AGGREGATE STATISTICS")
    print(f"{'='*70}")
    
    summary_stats = []
    for metric in ['accuracy', 'precision', 'recall', 'f1', 'mcc', 'roc_auc', 'pr_auc']:
        stats = {
            'metric': metric,
            'mean': df_all[metric].mean(),
            'std': df_all[metric].std(),
            'median': df_all[metric].median(),
            'min': df_all[metric].min(),
            'max': df_all[metric].max(),
            'q25': df_all[metric].quantile(0.25),
            'q75': df_all[metric].quantile(0.75)
        }
        summary_stats.append(stats)
        
        print(f"{metric.upper():12s}: μ={stats['mean']:.3f} ± {stats['std']:.3f}, "
              f"M={stats['median']:.3f}, range=[{stats['min']:.3f}, {stats['max']:.3f}]")
    
    df_summary = pd.DataFrame(summary_stats)
    
    summary_csv = f"{output_folder}/aggregate_metrics.csv"
    
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
        temp_csv = tmp.name
    df_summary.to_csv(temp_csv, index=False)
    
    if summary_csv.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = summary_csv.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_csv, bucket, key)
    
    os.remove(temp_csv)
    print(f"\n  ✓ Created: {summary_csv}")
    
    print(f"\n{'='*70}")
    print("✅ AGGREGATION COMPLETE")
    print(f"{'='*70}")
    
    return {
        'per_location_csv': agg_csv,
        'per_location_png': per_loc_png,
        'aggregate_csv': summary_csv,
        'aggregate_png': agg_png,
        'dataframe': df_all,
        'summary': df_summary
    }
