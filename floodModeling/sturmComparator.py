"""
sturmComparator.py - Complete Automated STURM Comparison

Automatically constructs paths based on sturmFloodStitcher.py naming convention
and compares XGBoost predictions to STURM ground truth.
"""

import os
import gc
import re
import tempfile
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from osgeo import gdal, osr
import boto3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    matthews_corrcoef, average_precision_score, roc_auc_score,
    confusion_matrix
)

# ======================== GDAL/PROJ Configuration ========================

def configure_gdal_proj():
    """Configure GDAL and PROJ library paths."""
    
    # Suppress specific warnings
    gdal.SetConfigOption('CPL_DEBUG', 'OFF')
    
    # Try to find PROJ data directory
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        if proj_dir and os.path.isdir(proj_dir):
            os.environ["PROJ_LIB"] = proj_dir
            gdal.SetConfigOption("PROJ_LIB", proj_dir)
            print(f"  Using PROJ from pyproj: {proj_dir}")
            return
    except:
        pass
    
    # Check common locations
    for cand in ["/opt/conda/share/proj", "/usr/share/proj", "/usr/local/share/proj"]:
        if os.path.isdir(cand):
            os.environ["PROJ_LIB"] = cand
            gdal.SetConfigOption("PROJ_LIB", cand)
            print(f"  Using PROJ from: {cand}")
            return
    
    print("  ⚠️  Could not locate PROJ data directory")
    print("  Continuing with default PROJ configuration...")

configure_gdal_proj()

# Enable PROJ network for datum transformations
gdal.SetConfigOption("PROJ_NETWORK", "ON")
gdal.SetConfigOption("GDAL_CACHEMAX", "256")
gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,tiff")
gdal.SetConfigOption("AWS_NO_SIGN_REQUEST", "NO")

# Suppress warnings about PROJ database
gdal.PushErrorHandler('CPLQuietErrorHandler')



# ======================== GDAL Utilities ========================

def to_vsis3(s3_uri: str) -> str:
    """Convert s3:// URI to /vsis3/ path."""
    return s3_uri.replace('s3://', '/vsis3/')


def extract_band_by_description(tiff_uri: str, description_pattern: str) -> Optional[Tuple[np.ndarray, object]]:
    """Extract band matching description pattern."""
    vsi_path = to_vsis3(tiff_uri)
    ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
    
    if ds is None:
        return None, None
    
    for b in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(b)
        desc = band.GetDescription()
        if description_pattern in desc:
            arr = band.ReadAsArray()
            return arr, ds
    
    # If no match, try first band
    if ds.RasterCount > 0:
        arr = ds.GetRasterBand(1).ReadAsArray()
        return arr, ds
    
    ds = None
    return None, None


def warp_to_match(src_uri: str, ref_ds: object,
                 resample_method: str = 'nearest',
                 src_nodata: int = 255,
                 dst_nodata: int = 255) -> Optional[np.ndarray]:
    """Warp source raster to match reference dataset with aggressive cleanup."""
    
    vsi_src = to_vsis3(src_uri)
    
    # Get reference parameters
    ref_gt = ref_ds.GetGeoTransform()
    ref_proj = ref_ds.GetProjection()
    ref_width = ref_ds.RasterXSize
    ref_height = ref_ds.RasterYSize
    
    # Calculate bounds
    minx = ref_gt[0]
    maxy = ref_gt[3]
    maxx = minx + ref_gt[1] * ref_width
    miny = maxy + ref_gt[5] * ref_height
    
    output_bounds = [minx, miny, maxx, maxy]
    
    # Create temporary output file
    temp_output = tempfile.mktemp(suffix='.tif')
    
    result_ds = None
    arr = None
    
    try:
        warp_options = gdal.WarpOptions(
            format='GTiff',
            creationOptions=['COMPRESS=DEFLATE', 'TILED=YES'],
            dstSRS=ref_proj,
            outputBounds=output_bounds,
            width=ref_width,
            height=ref_height,
            resampleAlg=resample_method,
            srcNodata=src_nodata,
            dstNodata=dst_nodata,
            multithread=False  # Disable to reduce memory
        )
        
        # Warp to temporary file
        result_ds = gdal.Warp(temp_output, vsi_src, options=warp_options)
        
        # Check if valid
        if result_ds is None or not isinstance(result_ds, gdal.Dataset):
            return None
        
        # Read array immediately
        band = result_ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        
        # Close dataset IMMEDIATELY
        result_ds = None
        band = None
        
        # Delete temp file IMMEDIATELY
        if os.path.exists(temp_output):
            try:
                os.remove(temp_output)
            except:
                pass
        
        return arr
        
    except Exception as e:
        print(f"    ⚠️  Warp error: {e}")
        return None
        
    finally:
        # Ensure cleanup in all cases
        if result_ds is not None:
            result_ds = None
        
        if os.path.exists(temp_output):
            try:
                os.remove(temp_output)
            except:
                pass
        
        # Force garbage collection
        gc.collect()


# ======================== Metrics ========================

def calculate_binary_metrics(y_pred: np.ndarray, y_true: np.ndarray,
                            pred_nodata: int = 255,
                            true_nodata: int = 255) -> Dict[str, float]:
    """Calculate comprehensive binary classification metrics."""
    
    valid_mask = (y_pred != pred_nodata) & (y_true != true_nodata)
    
    if valid_mask.sum() == 0:
        return {'valid_pixels': 0}
    
    y_pred_valid = y_pred[valid_mask]
    y_true_valid = y_true[valid_mask]
    
    metrics = {
        'valid_pixels': int(valid_mask.sum()),
        'accuracy': float(accuracy_score(y_true_valid, y_pred_valid)),
        'precision': float(precision_score(y_true_valid, y_pred_valid, zero_division=0)),
        'recall': float(recall_score(y_true_valid, y_pred_valid, zero_division=0)),
        'f1': float(f1_score(y_true_valid, y_pred_valid, zero_division=0)),
    }
    
    if len(np.unique(y_true_valid)) > 1 and len(np.unique(y_pred_valid)) > 1:
        metrics['mcc'] = float(matthews_corrcoef(y_true_valid, y_pred_valid))
    else:
        metrics['mcc'] = 0.0
    
    try:
        metrics['roc_auc'] = float(roc_auc_score(y_true_valid, y_pred_valid))
    except:
        metrics['roc_auc'] = np.nan
    
    try:
        metrics['pr_auc'] = float(average_precision_score(y_true_valid, y_pred_valid))
    except:
        metrics['pr_auc'] = np.nan
    
    cm = confusion_matrix(y_true_valid, y_pred_valid, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    
    metrics.update({
        'tp': int(tp), 'fp': int(fp), 'tn': int(tn), 'fn': int(fn),
        'pred_positive_frac': float(y_pred_valid.mean()),
        'true_positive_frac': float(y_true_valid.mean())
    })
    
    return metrics


def create_confusion_matrix_array(y_pred: np.ndarray, y_true: np.ndarray,
                                 pred_nodata: int = 255,
                                 true_nodata: int = 255,
                                 output_nodata: int = 255) -> np.ndarray:
    """Create confusion matrix array: 0=TN, 1=TP, 2=FN, 3=FP."""
    
    cm = np.full(y_pred.shape, output_nodata, dtype=np.uint8)
    valid = (y_pred != pred_nodata) & (y_true != true_nodata)
    
    if valid.sum() == 0:
        return cm
    
    cm[(y_pred == 0) & (y_true == 0) & valid] = 0  # TN
    cm[(y_pred == 1) & (y_true == 1) & valid] = 1  # TP
    cm[(y_pred == 0) & (y_true == 1) & valid] = 2  # FN
    cm[(y_pred == 1) & (y_true == 0) & valid] = 3  # FP
    
    return cm


# ======================== Visualization ========================

def create_confusion_matrix_png(cm_array: np.ndarray, metrics: Dict,
                               output_uri: str, sensor: str, watershed: str, date_str: str):
    """Create confusion matrix PNG visualization."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    cm_disp = cm_array.astype(float)
    cm_disp[cm_array == 255] = np.nan
    
    cmap = matplotlib.colors.ListedColormap(['green', 'blue', 'orange', 'red'])
    ax.imshow(cm_disp, cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
    
    legend_elements = [
        mpatches.Patch(color='green', label='True Negative (Land correct)'),
        mpatches.Patch(color='blue', label='True Positive (Water correct)'),
        mpatches.Patch(color='orange', label='False Negative (Missed water)'),
        mpatches.Patch(color='red', label='False Positive (False alarm)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    ax.set_title(f'{sensor} vs XGBoost Confusion Matrix\n{watershed} - {date_str}',
                fontsize=14, fontweight='bold')
    
    text = (
        f"Watershed: {watershed}\n"
        f"Sensor: {sensor}\n"
        f"Date: {date_str}\n"
        f"Valid pixels: {metrics.get('valid_pixels', 0):,}\n"
        f"MCC: {metrics.get('mcc', 0):.3f}\n"
        f"F1 Score: {metrics.get('f1', 0):.3f}\n"
        f"Accuracy: {metrics.get('accuracy', 0):.3f}\n"
        f"Precision: {metrics.get('precision', 0):.3f}\n"
        f"Recall: {metrics.get('recall', 0):.3f}\n"
        f"PR-AUC: {metrics.get('pr_auc', 0):.3f}\n"
        f"ROC-AUC: {metrics.get('roc_auc', 0):.3f}"
    )
    props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11,
           va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    plt.savefig(temp_png, dpi=150, bbox_inches='tight')
    plt.close()
    
    if output_uri.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_uri.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_png, bucket, key)
    
    os.remove(temp_png)


def save_tiff(array: np.ndarray, geotransform: tuple, projection: str,
             output_uri: str, nodata_value: int = 255):
    """Save numpy array as GeoTIFF."""
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_tiff = tmp.name
    
    driver = gdal.GetDriverByName('GTiff')
    H, W = array.shape
    out_ds = driver.Create(temp_tiff, W, H, 1, gdal.GDT_Byte,
                          options=['COMPRESS=DEFLATE', 'TILED=YES'])
    
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection(projection)
    
    band = out_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata_value)
    band.WriteArray(array)
    
    out_ds = None
    
    if output_uri.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_uri.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_tiff, bucket, key)
    
    os.remove(temp_tiff)


# ======================== Path Discovery ========================

def discover_watershed_predictions(assess_folder: str) -> Dict[str, str]:
    """Discover all watershed prediction TIFFs."""
    
    s3_client = boto3.client('s3')
    bucket, prefix = assess_folder.replace('s3://', '').split('/', 1)
    
    response = s3_client.list_objects_v2(
        Bucket=bucket,
        Prefix=prefix + '/',
        Delimiter='/'
    )
    
    watershed_predictions = {}
    
    if 'CommonPrefixes' in response:
        for prefix_obj in response['CommonPrefixes']:
            ws_prefix = prefix_obj['Prefix']
            ws_name = ws_prefix.rstrip('/').split('/')[-1]
            
            if ws_name.startswith('Lat_'):
                pred_tiff = f"s3://{bucket}/{ws_prefix}xgboost_predictions_{ws_name}.tif"
                
                try:
                    tiff_bucket, tiff_key = pred_tiff.replace('s3://', '').split('/', 1)
                    s3_client.head_object(Bucket=tiff_bucket, Key=tiff_key)
                    watershed_predictions[ws_name] = pred_tiff
                except:
                    pass
    
    return watershed_predictions

def discover_available_sturm_dates(sturm_folder: str, sturm_code: str) -> Dict[str, List[str]]:
    """
    Auto-discover what dates actually have STURM data.
    
    Returns:
        {'S1': ['2019-11-12', ...], 'S2': ['2019-11-13', ...]}
    """
    s3_client = boto3.client('s3')
    bucket, prefix = sturm_folder.replace('s3://', '').split('/', 1)
    
    available = {'S1': [], 'S2': []}
    
    for sensor in ['S1', 'S2']:
        sensor_prefix = f"{prefix}/{sensor}/"
        
        try:
            response = s3_client.list_objects_v2(Bucket=bucket, Prefix=sensor_prefix)
            
            if 'Contents' in response:
                for obj in response['Contents']:
                    filename = os.path.basename(obj['Key'])
                    
                    # Extract date from filename: EMSR###_YYYY-MM-DD_mask.tif
                    if sturm_code in filename and '_mask.tif' in filename:
                        match = re.search(r'(\d{4}-\d{2}-\d{2})', filename)
                        if match:
                            date_str = match.group(1)
                            if date_str not in available[sensor]:
                                available[sensor].append(date_str)
        except:
            pass
    
    return available

# Add after the existing visualization functions, before compare_all_watersheds_to_sturm

def stitch_confusion_matrices(watershed_cm_tiffs: List[str], 
                              output_tiff_uri: str,
                              target_crs: str = "EPSG:4326",
                              resolution: Optional[float] = None) -> bool:
    """Stitch individual watershed confusion matrix TIFFs into one mosaic."""
    
    if not watershed_cm_tiffs:
        return False
    
    print(f"  Stitching {len(watershed_cm_tiffs)} confusion matrix TIFFs...")
    
    vsis3_paths = [uri.replace('s3://', '/vsis3/') for uri in watershed_cm_tiffs]
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        warp_options = gdal.WarpOptions(
            format='GTiff',
            creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
            dstSRS=target_crs,
            xRes=resolution,
            yRes=resolution,
            resampleAlg='nearest',
            dstNodata=255,
            srcNodata=255,
            multithread=False
        )
        
        result = gdal.Warp(temp_path, vsis3_paths, options=warp_options)
        
        if result is None:
            print(f"    ⚠️  Failed to create confusion matrix mosaic")
            os.remove(temp_path)
            return False
        
        print(f"    ✓ Created mosaic: {result.RasterXSize}x{result.RasterYSize}")
        
        result = None
        
        # Upload to S3
        if output_tiff_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_tiff_uri.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_path, bucket, key)
            
            # Verify with S3 API
            import time
            max_retries = 5
            for attempt in range(max_retries):
                try:
                    s3_client.head_object(Bucket=bucket, Key=key)
                    print(f"    ✓ Verified upload to S3")
                    break
                except:
                    if attempt < max_retries - 1:
                        time.sleep(2)
                    else:
                        os.remove(temp_path)
                        return False
            
            # ✅ NEW: Verify GDAL can open it via /vsis3/
            vsi_output_path = output_tiff_uri.replace('s3://', '/vsis3/')
            gdal_verified = False
            
            for attempt in range(5):
                try:
                    test_ds = gdal.Open(vsi_output_path, gdal.GA_ReadOnly)
                    if test_ds is not None:
                        # Successfully opened
                        test_ds = None
                        gdal_verified = True
                        print(f"    ✓ Verified GDAL can access via /vsis3/")
                        break
                    else:
                        if attempt < 4:
                            print(f"    ⏳ GDAL cache not ready, waiting... (attempt {attempt + 1}/5)")
                            time.sleep(3)  # Longer wait for GDAL
                except Exception as e:
                    if attempt < 4:
                        print(f"    ⏳ GDAL open failed, waiting... (attempt {attempt + 1}/5)")
                        time.sleep(3)
                    
            if not gdal_verified:
                print(f"    ⚠️  File uploaded to S3 but GDAL cannot access it yet")
                os.remove(temp_path)
                return False
        
        os.remove(temp_path)
        print(f"    ✓ Saved to: {output_tiff_uri}")
        
        gc.collect()
        return True
        
    except Exception as e:
        print(f"    ⚠️  Error: {e}")
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return False


def create_stitched_confusion_matrix_png(stitched_tiff_uri: str,
                                        output_png_uri: str,
                                        date_str: str,
                                        aggregate_metrics: Dict[str, float]):
    """Create PNG visualization of stitched confusion matrix."""
    
    print(f"  Creating stitched confusion matrix PNG...")
    
    # Load stitched TIFF
    vsi_path = stitched_tiff_uri.replace('s3://', '/vsis3/')
    ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
    
    if ds is None:
        print(f"    ⚠️  Could not open stitched TIFF")
        return
    
    cm_array = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    
    # ✅ FIX: Find extent of valid (non-NoData) pixels
    valid_mask = cm_array != 255
    
    if not valid_mask.any():
        print(f"    ⚠️  No valid data in confusion matrix")
        return
    
    # Get bounding box of valid data
    rows, cols = np.where(valid_mask)
    row_min, row_max = rows.min(), rows.max()
    col_min, col_max = cols.min(), cols.max()
    
    # Add 5% padding on each side
    padding_rows = max(1, int(0.05 * (row_max - row_min)))
    padding_cols = max(1, int(0.05 * (col_max - col_min)))
    
    row_min = max(0, row_min - padding_rows)
    row_max = min(cm_array.shape[0] - 1, row_max + padding_rows)
    col_min = max(0, col_min - padding_cols)
    col_max = min(cm_array.shape[1] - 1, col_max + padding_cols)
    
    print(f"    Data extent: rows [{row_min}:{row_max}], cols [{col_min}:{col_max}]")
    
    # Create visualization
    fig, ax = plt.subplots(figsize=(16, 12))
    
    cm_disp = cm_array.astype(float)
    cm_disp[cm_array == 255] = np.nan
    
    cmap = matplotlib.colors.ListedColormap(['green', 'blue', 'orange', 'red'])
    
    # Display full array with extent
    H, W = cm_array.shape
    ax.imshow(cm_disp, cmap=cmap, vmin=0, vmax=3, interpolation='nearest',
             extent=[0, W, H, 0], origin='upper')
    
    # ✅ Zoom to data extent
    ax.set_xlim(col_min, col_max)
    ax.set_ylim(row_max, row_min)  # Reversed because origin='upper'
    
    # Legend
    legend_elements = [
        mpatches.Patch(color='green', label='True Negative (Land correct)'),
        mpatches.Patch(color='blue', label='True Positive (Water correct)'),
        mpatches.Patch(color='orange', label='False Negative (Missed water)'),
        mpatches.Patch(color='red', label='False Positive (False alarm)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    ax.set_title(f'All Watersheds Confusion Matrix - {date_str}\nXGBoost vs STURM Ground Truth',
                fontsize=16, fontweight='bold', pad=20)
    
    # Aggregate metrics text box
    if aggregate_metrics:
        text = (
            f"Date: {date_str}\n"
            f"Aggregated Metrics:\n"
            f"F1 Score: {aggregate_metrics.get('f1_mean', 0):.3f} ± {aggregate_metrics.get('f1_std', 0):.3f}\n"
            f"MCC: {aggregate_metrics.get('mcc_mean', 0):.3f} ± {aggregate_metrics.get('mcc_std', 0):.3f}\n"
            f"Accuracy: {aggregate_metrics.get('accuracy_mean', 0):.3f} ± {aggregate_metrics.get('accuracy_std', 0):.3f}\n"
            f"Precision: {aggregate_metrics.get('precision_mean', 0):.3f} ± {aggregate_metrics.get('precision_std', 0):.3f}\n"
            f"Recall: {aggregate_metrics.get('recall_mean', 0):.3f} ± {aggregate_metrics.get('recall_std', 0):.3f}\n"
            f"Watersheds: {aggregate_metrics.get('n_watersheds', 0)}"
        )
        props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='black', linewidth=2)
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=12,
               va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save (rest of code unchanged)
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    
    try:
        plt.savefig(temp_png, dpi=200, bbox_inches='tight')
        plt.close()
        
        if not os.path.exists(temp_png):
            print(f"    ❌ ERROR: Failed to create temp PNG file")
            return
        
        file_size = os.path.getsize(temp_png)
        print(f"    ✓ Temp PNG created: {file_size:,} bytes")
        
        if output_png_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_png_uri.replace('s3://', '').split('/', 1)
            
            print(f"    Uploading to s3://{bucket}/{key}...")
            s3_client.upload_file(temp_png, bucket, key)
            
            response = s3_client.head_object(Bucket=bucket, Key=key)
            uploaded_size = response['ContentLength']
            
            if uploaded_size != file_size:
                print(f"    ⚠️  WARNING: Size mismatch! Local={file_size}, S3={uploaded_size}")
            else:
                print(f"    ✓ Upload verified: {uploaded_size:,} bytes")
            
            print(f"    ✓ Created: {output_png_uri}")
        
        os.remove(temp_png)
        
    except Exception as e:
        print(f"    ❌ ERROR creating PNG: {e}")
        import traceback
        traceback.print_exc()
        
        if os.path.exists(temp_png):
            try:
                os.remove(temp_png)
            except:
                pass

# Add after create_stitched_confusion_matrix_png(), before compare_all_watersheds_to_sturm


def stitch_prediction_tiffs(watershed_pred_tiffs: List[str],
                           output_tiff_uri: str,
                           target_crs: str = "EPSG:4326",
                           resolution: Optional[float] = None) -> bool:
    """Stitch individual watershed prediction TIFFs into one mosaic with non-blocking verification.
    
    Args:
        watershed_pred_tiffs: List of S3 URIs to individual watershed prediction TIFFs
        output_tiff_uri: S3 URI for output stitched prediction TIFF
        target_crs: Target coordinate reference system
        resolution: Target resolution (auto if None)
        
    Returns:
        True if successful, False otherwise
    """
    
    if not watershed_pred_tiffs:
        return False
    
    print(f"  Stitching {len(watershed_pred_tiffs)} prediction TIFFs...")
    
    vsis3_paths = [uri.replace('s3://', '/vsis3/') for uri in watershed_pred_tiffs]
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        warp_options = gdal.WarpOptions(
            format='GTiff',
            creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
            dstSRS=target_crs,
            xRes=resolution,
            yRes=resolution,
            resampleAlg='nearest',  # Nearest neighbor for discrete classifications
            dstNodata=255,
            srcNodata=255,
            multithread=False
        )
        
        result = gdal.Warp(temp_path, vsis3_paths, options=warp_options)
        
        if result is None:
            print(f"    ⚠️  Failed to create prediction mosaic")
            if os.path.exists(temp_path):
                os.remove(temp_path)
            return False
        
        print(f"    ✓ Created prediction mosaic: {result.RasterXSize}x{result.RasterYSize}")
        result = None
        
        # Upload to S3
        if output_tiff_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_tiff_uri.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_path, bucket, key)
            
            # Verify with S3 API (REQUIRED - must pass)
            import time
            s3_verified = False
            for attempt in range(5):
                try:
                    s3_client.head_object(Bucket=bucket, Key=key)
                    print(f"    ✓ Verified upload to S3")
                    s3_verified = True
                    break
                except:
                    if attempt < 4:
                        print(f"    ⏳ Waiting for S3 propagation (attempt {attempt + 1}/5)...")
                        time.sleep(2)
            
            if not s3_verified:
                print(f"    ❌ S3 upload verification failed")
                if os.path.exists(temp_path):
                    os.remove(temp_path)
                return False
            
            # GDAL /vsis3/ verification (OPTIONAL - warn but don't fail)
            vsi_output_path = output_tiff_uri.replace('s3://', '/vsis3/')
            gdal_verified = False
            
            for attempt in range(5):
                try:
                    test_ds = gdal.Open(vsi_output_path, gdal.GA_ReadOnly)
                    if test_ds is not None:
                        test_ds = None
                        gdal_verified = True
                        print(f"    ✓ Verified GDAL can access via /vsis3/")
                        break
                except:
                    pass
                
                if attempt < 4:
                    print(f"    ⏳ GDAL cache not ready, waiting... (attempt {attempt + 1}/5)")
                    time.sleep(3)
            
            if not gdal_verified:
                # ✅ CRITICAL CHANGE: Don't fail, just warn
                print(f"    ⚠️  GDAL verification timed out (cache lag)")
                print(f"    ➜  Continuing anyway - downstream operations may succeed")
                print(f"    ➜  File is confirmed in S3 and will be accessible shortly")
        
        # Clean up temp file
        if os.path.exists(temp_path):
            os.remove(temp_path)
        
        print(f"    ✓ Saved to: {output_tiff_uri}")
        gc.collect()
        
        # ✅ Return True even if GDAL verification failed (S3 verification is enough)
        return True
        
    except Exception as e:
        print(f"    ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return False


def create_stitched_prediction_png(stitched_tiff_uri: str,
                                   output_png_uri: str,
                                   date_str: str,
                                   n_watersheds: int):
    """Create PNG visualization of stitched binary predictions."""
    
    print(f"  Creating stitched prediction PNG...")
    
    # Load stitched TIFF
    vsi_path = stitched_tiff_uri.replace('s3://', '/vsis3/')
    ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
    
    if ds is None:
        print(f"    ⚠️  Could not open stitched prediction TIFF")
        return
    
    pred_array = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    
    # ✅ FIX: Find extent of valid data
    valid_mask = pred_array != 255
    
    if not valid_mask.any():
        print(f"    ⚠️  No valid data in predictions")
        return
    
    rows, cols = np.where(valid_mask)
    row_min, row_max = rows.min(), rows.max()
    col_min, col_max = cols.min(), cols.max()
    
    # Add 5% padding
    padding_rows = max(1, int(0.05 * (row_max - row_min)))
    padding_cols = max(1, int(0.05 * (col_max - col_min)))
    
    row_min = max(0, row_min - padding_rows)
    row_max = min(pred_array.shape[0] - 1, row_max + padding_rows)
    col_min = max(0, col_min - padding_cols)
    col_max = min(pred_array.shape[1] - 1, col_max + padding_cols)
    
    print(f"    Data extent: rows [{row_min}:{row_max}], cols [{col_min}:{col_max}]")
    
    # Create visualization
    fig, ax = plt.subplots(figsize=(16, 12))
    
    pred_disp = pred_array.astype(float)
    pred_disp[pred_array == 255] = np.nan
    
    cmap = matplotlib.colors.ListedColormap(['tan', 'blue'])
    
    H, W = pred_array.shape
    ax.imshow(pred_disp, cmap=cmap, vmin=0, vmax=1, interpolation='nearest',
             extent=[0, W, H, 0], origin='upper')
    
    # ✅ Zoom to data extent
    ax.set_xlim(col_min, col_max)
    ax.set_ylim(row_max, row_min)
    
    # Legend
    legend_elements = [
        mpatches.Patch(color='tan', label='Predicted Land'),
        mpatches.Patch(color='blue', label='Predicted Water')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    ax.set_title(f'All Watersheds XGBoost Predictions - {date_str}',
                fontsize=16, fontweight='bold', pad=20)
    
    # Calculate statistics
    if valid_mask.sum() > 0:
        water_fraction = pred_array[valid_mask].mean()
        
        text = (
            f"Date: {date_str}\n"
            f"Model: XGBoost\n"
            f"Watersheds: {n_watersheds}\n"
            f"Valid pixels: {valid_mask.sum():,}\n"
            f"Predicted water: {water_fraction:.3f}"
        )
        props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='black', linewidth=2)
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=12,
               va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save (rest unchanged)
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    
    try:
        plt.savefig(temp_png, dpi=200, bbox_inches='tight')
        plt.close()
        
        if not os.path.exists(temp_png):
            print(f"    ❌ ERROR: Failed to create temp PNG file")
            return
        
        file_size = os.path.getsize(temp_png)
        print(f"    ✓ Temp PNG created: {file_size:,} bytes")
        
        if output_png_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_png_uri.replace('s3://', '').split('/', 1)
            
            print(f"    Uploading to s3://{bucket}/{key}...")
            s3_client.upload_file(temp_png, bucket, key)
            
            response = s3_client.head_object(Bucket=bucket, Key=key)
            uploaded_size = response['ContentLength']
            
            if uploaded_size != file_size:
                print(f"    ⚠️  WARNING: Size mismatch! Local={file_size}, S3={uploaded_size}")
            else:
                print(f"    ✓ Upload verified: {uploaded_size:,} bytes")
            
            print(f"    ✓ Created: {output_png_uri}")
        
        os.remove(temp_png)
        
    except Exception as e:
        print(f"    ❌ ERROR creating prediction PNG: {e}")
        import traceback
        traceback.print_exc()
        
        if os.path.exists(temp_png):
            try:
                os.remove(temp_png)
            except:
                pass

def create_floodmask_png(floodmask_array: np.ndarray, output_uri: str, 
                        sensor: str, watershed: str, date_str: str):
    """Create PNG visualization of aligned STURM floodmask."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    fm_disp = floodmask_array.astype(float)
    fm_disp[floodmask_array == 255] = np.nan
    
    cmap = matplotlib.colors.ListedColormap(['tan', 'blue'])
    ax.imshow(fm_disp, cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
    
    legend_elements = [
        mpatches.Patch(color='tan', label='Land (Not Flooded)'),
        mpatches.Patch(color='blue', label='Water (Flooded)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    ax.set_title(f'STURM {sensor} Flood Mask\n{watershed} - {date_str}',
                fontsize=14, fontweight='bold')
    
    water_fraction = floodmask_array[floodmask_array != 255].mean() if (floodmask_array != 255).any() else 0.0
    
    text = (
        f"Watershed: {watershed}\n"
        f"Sensor: {sensor}\n"
        f"Date: {date_str}\n"
        f"Source: STURM Ground Truth\n"
        f"Water Fraction: {water_fraction:.3f}"
    )
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.9)
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11,
           va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    plt.savefig(temp_png, dpi=150, bbox_inches='tight')
    plt.close()
    
    if output_uri.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_uri.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_png, bucket, key)
    
    os.remove(temp_png)


def create_imagery_png(imagery_array: np.ndarray, output_uri: str,
                      sensor: str, watershed: str, date_str: str):
    """Create PNG visualization of aligned STURM imagery."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    img_disp = imagery_array.astype(float)
    img_disp[imagery_array == 255] = np.nan
    
    # Use grayscale for imagery
    vmin, vmax = np.nanpercentile(img_disp, [2, 98])
    
    ax.imshow(img_disp, cmap='gray', vmin=vmin, vmax=vmax, interpolation='nearest')
    
    ax.set_title(f'STURM {sensor} Imagery\n{watershed} - {date_str}',
                fontsize=14, fontweight='bold')
    
    text = (
        f"Watershed: {watershed}\n"
        f"Sensor: {sensor}\n"
        f"Date: {date_str}\n"
        f"Source: STURM\n"
        f"Display: 2-98% stretch"
    )
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.9)
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11,
           va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    plt.savefig(temp_png, dpi=150, bbox_inches='tight')
    plt.close()
    
    if output_uri.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = output_uri.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_png, bucket, key)
    
    os.remove(temp_png)

def stitch_imagery_tiffs(watershed_img_tiffs: List[str],
                        output_tiff_uri: str,
                        target_crs: str = "EPSG:4326",
                        resolution: Optional[float] = None) -> bool:
    """Stitch individual watershed imagery TIFFs into one mosaic with S3 verification.
    
    Args:
        watershed_img_tiffs: List of S3 URIs to individual watershed imagery TIFFs
        output_tiff_uri: S3 URI for output stitched TIFF
        target_crs: Target coordinate reference system
        resolution: Target resolution (auto if None)
        
    Returns:
        True if successful, False otherwise
    """
    
    if not watershed_img_tiffs:
        return False
    
    print(f"  Stitching {len(watershed_img_tiffs)} imagery TIFFs...")
    
    vsis3_paths = [uri.replace('s3://', '/vsis3/') for uri in watershed_img_tiffs]
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        warp_options = gdal.WarpOptions(
            format='GTiff',
            creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
            dstSRS=target_crs,
            xRes=resolution,
            yRes=resolution,
            resampleAlg='bilinear',  # Better for continuous imagery
            dstNodata=255,
            srcNodata=255,
            multithread=False
        )
        
        result = gdal.Warp(temp_path, vsis3_paths, options=warp_options)
        
        if result is None:
            print(f"    ⚠️  Failed to create imagery mosaic")
            if os.path.exists(temp_path):
                os.remove(temp_path)
            return False
        
        print(f"    ✓ Created imagery mosaic: {result.RasterXSize}x{result.RasterYSize}")
        
        result = None
        
        # Upload to S3
        if output_tiff_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_tiff_uri.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_path, bucket, key)
            
            # ✅ FIX: Verify file exists in S3 before returning (S3 eventual consistency)
            import time
            max_retries = 5
            for attempt in range(max_retries):
                try:
                    s3_client.head_object(Bucket=bucket, Key=key)
                    print(f"    ✓ Verified upload to S3")
                    break
                except:
                    if attempt < max_retries - 1:
                        print(f"    ⏳ Waiting for S3 propagation (attempt {attempt + 1}/{max_retries})...")
                        time.sleep(2)  # Wait 2 seconds
                    else:
                        print(f"    ⚠️  File uploaded but not yet visible in S3")
                        if os.path.exists(temp_path):
                            os.remove(temp_path)
                        return False
        
        if os.path.exists(temp_path):
            os.remove(temp_path)
        
        print(f"    ✓ Saved to: {output_tiff_uri}")
        
        gc.collect()
        
        return True
        
    except Exception as e:
        print(f"    ⚠️  Error: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return False



def create_sturm_imagery_png_direct(sturm_stitched_tiff: str,
                                   output_png_uri: str,
                                   sensor: str,
                                   date_str: str):
    """Create PNG directly from stitched STURM imagery file.
    For S2, reads bands 2,3,4 as Blue,Green,Red for true color.
    For S1, uses grayscale.
    """
    
    print(f"  Creating {sensor} imagery PNG from stitched STURM file...")
    
    vsi_path = sturm_stitched_tiff.replace('s3://', '/vsis3/')
    ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
    
    if ds is None:
        print(f"    ⚠️  Could not open {sturm_stitched_tiff}")
        return
    
    n_bands = ds.RasterCount
    print(f"    TIFF has {n_bands} band(s)")
    
    # S2 TRUE COLOR (RGB) - Requires at least 4 bands to access band 4
    if sensor == 'S2' and n_bands >= 4:  # ✅ FIXED: Changed from >= 3 to >= 4
        print(f"    Processing as true color RGB (bands 2,3,4)...")
        
        # Read Sentinel-2 bands for natural color
        # Band 2 = Blue (490nm), Band 3 = Green (560nm), Band 4 = Red (665nm)
        blue = ds.GetRasterBand(1).ReadAsArray().astype(float)   # S2 Band 2 (Blue, 490nm)
        green = ds.GetRasterBand(2).ReadAsArray().astype(float)  # S2 Band 3 (Green, 560nm)
        red = ds.GetRasterBand(3).ReadAsArray().astype(float)    # S2 Band 4 (Red, 665nm)
        ds = None
        
        # Mask NoData
        red[red == 255] = np.nan
        green[green == 255] = np.nan
        blue[blue == 255] = np.nan
        
        valid_mask = ~(np.isnan(red) | np.isnan(green) | np.isnan(blue))
        
        if not valid_mask.any():
            print(f"    ⚠️  No valid data")
            return
        
        # Find data extent
        rows, cols = np.where(valid_mask)
        row_min, row_max = rows.min(), rows.max()
        col_min, col_max = cols.min(), cols.max()
        
        padding = max(1, int(0.05 * max(row_max - row_min, col_max - col_min)))
        row_min = max(0, row_min - padding)
        row_max = min(red.shape[0] - 1, row_max + padding)
        col_min = max(0, col_min - padding)
        col_max = min(red.shape[1] - 1, col_max + padding)
        
        # Stack as RGB
        rgb_array = np.dstack([red, green, blue])
        
        # Per-channel 2-98% stretch
        rgb_stretched = np.zeros_like(rgb_array)
        
        for i, color_name in enumerate(['Red', 'Green', 'Blue']):
            channel = rgb_array[:, :, i]
            vmin, vmax = np.nanpercentile(channel, [2, 98])
            print(f"      {color_name} channel: range=[{vmin:.3f}, {vmax:.3f}]")
            
            # Normalize to 0-1
            channel_norm = (channel - vmin) / (vmax - vmin + 1e-8)
            channel_norm = np.clip(channel_norm, 0, 1)
            rgb_stretched[:, :, i] = channel_norm
        
        # Create figure
        fig, ax = plt.subplots(figsize=(16, 12))
        H, W = red.shape
        
        ax.imshow(rgb_stretched, interpolation='nearest',
                 extent=[0, W, H, 0], origin='upper')
        
        ax.set_xlim(col_min, col_max)
        ax.set_ylim(row_max, row_min)
        
        ax.set_title(f'STURM S2 True Color Imagery - {date_str}',
                    fontsize=16, fontweight='bold', pad=20)
        
        text = (
            f"Date: {date_str}\n"
            f"Sensor: S2 (Optical)\n"
            f"Source: STURM Ground Truth\n"
            f"Display: Natural Color RGB\n"
            f"Bands: 4 (Red), 3 (Green), 2 (Blue)\n"
            f"Stretch: 2-98% per channel"
        )
        props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.9,
                    edgecolor='darkgreen', linewidth=2)
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=12,
               va='top', bbox=props, family='monospace')
    
    # S1 GRAYSCALE (or S2 fallback if insufficient bands)
    else:
        print(f"    Processing as grayscale...")
        
        img = ds.GetRasterBand(1).ReadAsArray().astype(float)
        ds = None
        
        img[img == 255] = np.nan
        valid_mask = ~np.isnan(img)
        
        if not valid_mask.any():
            return
        
        rows, cols = np.where(valid_mask)
        row_min, row_max = rows.min(), rows.max()
        col_min, col_max = cols.min(), cols.max()
        
        padding = max(1, int(0.05 * max(row_max - row_min, col_max - col_min)))
        row_min = max(0, row_min - padding)
        row_max = min(img.shape[0] - 1, row_max + padding)
        col_min = max(0, col_min - padding)
        col_max = min(img.shape[1] - 1, col_max + padding)
        
        vmin, vmax = np.nanpercentile(img, [2, 98])
        print(f"      Grayscale range: [{vmin:.3f}, {vmax:.3f}]")
        
        fig, ax = plt.subplots(figsize=(16, 12))
        H, W = img.shape
        
        ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, interpolation='nearest',
                 extent=[0, W, H, 0], origin='upper')
        
        ax.set_xlim(col_min, col_max)
        ax.set_ylim(row_max, row_min)
        
        sensor_type = "SAR" if sensor == "S1" else "Optical"
        ax.set_title(f'STURM {sensor} Imagery - {date_str}',
                    fontsize=16, fontweight='bold', pad=20)
        
        text = (
            f"Date: {date_str}\n"
            f"Sensor: {sensor} ({sensor_type})\n"
            f"Source: STURM Ground Truth\n"
            f"Display: Grayscale\n"
            f"Range: [{vmin:.2f}, {vmax:.2f}]"
        )
        props = dict(boxstyle='round', facecolor='lightgray', alpha=0.9)
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=12,
               va='top', bbox=props, family='monospace')
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_png = tmp.name
    
    try:
        plt.savefig(temp_png, dpi=200, bbox_inches='tight')
        plt.close()
        
        if output_png_uri.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = output_png_uri.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_png, bucket, key)
            print(f"    ✓ Created: {output_png_uri}")
        
        os.remove(temp_png)
    except Exception as e:
        print(f"    ❌ ERROR: {e}")
        if os.path.exists(temp_png):
            os.remove(temp_png)


# ======================== Main Comparison Function ========================

def compare_all_watersheds_to_sturm(
    latitude: str,
    longitude: str,
    sturm_code: str,
    date_list: List[str],
    bucket: str = "climate-ai-data-science-datasets",
    flood_base: str = "arrakis-data/floodOutputs",
    sturm_nodata: int = 255,
    pred_nodata: int = 255
) -> Dict:
    """Complete automated comparison workflow with stitched outputs.
    
    Args:
        latitude: Base latitude (e.g., "10.7875")
        longitude: Base longitude (e.g., "0.7479")
        sturm_code: STURM event code (e.g., "EMSR470")
        date_list: List of dates in YYYYMMDD format
        bucket: S3 bucket name
        flood_base: Base prefix for flood outputs
        sturm_nodata: NoData value in STURM data
        pred_nodata: NoData value in predictions
        
    Returns:
        Dictionary with comparison results including stitched outputs
    """
    
    print(f"\n{'='*70}")
    print(f"AUTOMATED STURM COMPARISON WORKFLOW")
    print(f"{'='*70}")
    print(f"Base Location: Lat {latitude}, Lon {longitude}")
    print(f"STURM Code: {sturm_code}")
    print(f"Dates: {date_list}")
    
    # Construct paths
    training_folder = f"s3://{bucket}/{flood_base}/trainingFolderFor_Lat_{latitude}_Lon_{longitude}"
    assess_folder = f"{training_folder}/assessXGBoostModel"
    sturm_folder = f"{training_folder}/sturm_stitched"
    
    print(f"\nTraining folder: {training_folder}")
    print(f"Assessment folder: {assess_folder}")
    print(f"STURM folder: {sturm_folder}")
    
    # Discover watershed predictions
    print(f"\n{'='*70}")
    print("DISCOVERING WATERSHED PREDICTIONS")
    print(f"{'='*70}")
    
    watershed_predictions = discover_watershed_predictions(assess_folder)
    
    if not watershed_predictions:
        print("❌ No watershed predictions found")
        return {}
    
    print(f"Found {len(watershed_predictions)} watersheds:")
    for ws_id in watershed_predictions:
        print(f"  - {ws_id}")
    
    # Auto-discover available STURM dates
    print(f"\n{'='*70}")
    print("DISCOVERING AVAILABLE STURM DATES")
    print(f"{'='*70}")
    
    available_sturm = discover_available_sturm_dates(sturm_folder, sturm_code)
    
    print(f"Available STURM data:")
    for sensor, dates in available_sturm.items():
        if dates:
            print(f"  {sensor}: {', '.join(sorted(dates))}")
        else:
            print(f"  {sensor}: None")
    
    # Find intersection with requested dates
    date_list_formatted = [f"{d[:4]}-{d[4:6]}-{d[6:]}" for d in date_list]
    
    all_available_dates = set(available_sturm['S1'] + available_sturm['S2'])
    requested_dates_set = set(date_list_formatted)
    
    matching_dates = all_available_dates & requested_dates_set
    
    if not matching_dates:
        print(f"\n❌ No overlap between requested dates and available STURM data!")
        print(f"Requested: {sorted(requested_dates_set)}")
        print(f"Available: {sorted(all_available_dates)}")
        print(f"\nSuggestion: Update your date_list to match available STURM dates")
        return {}
    
    print(f"\n✓ Processing {len(matching_dates)} matching dates: {sorted(matching_dates)}")
    
    # Build sturm_files dict only for available dates
    print(f"\n{'='*70}")
    print("LOCATING STURM FLOOD MAPS")
    print(f"{'='*70}")

    s3_client = boto3.client('s3') 
    
    sturm_files = {}
    
    for date_fmt in sorted(matching_dates):
        # Convert back to YYYYMMDD for consistency
        date_orig = date_fmt.replace('-', '')
        sturm_files[date_orig] = {}
        
        # Check S1
        if date_fmt in available_sturm['S1']:
            s1_fm_path = f"{sturm_folder}/S1/{sturm_code}_{date_fmt}_mask.tif"
            s1_img_path = f"{sturm_folder}/S1/{sturm_code}_{date_fmt}_imagery.tif"
            
            try:
                p_bucket, p_key = s1_fm_path.replace('s3://', '').split('/', 1)
                s3_client.head_object(Bucket=p_bucket, Key=p_key)
                sturm_files[date_orig]['s1_floodmask'] = s1_fm_path
                print(f"  ✓ Found S1 mask: {date_fmt}")
            except:
                print(f"  ⚠️  Expected S1 mask missing: {date_fmt}")
            
            try:
                p_bucket, p_key = s1_img_path.replace('s3://', '').split('/', 1)
                s3_client.head_object(Bucket=p_bucket, Key=p_key)
                sturm_files[date_orig]['s1_imagery'] = s1_img_path
                print(f"  ✓ Found S1 imagery: {date_fmt}")
            except:
                pass
        
        # Check S2
        if date_fmt in available_sturm['S2']:
            s2_fm_path = f"{sturm_folder}/S2/{sturm_code}_{date_fmt}_mask.tif"
            s2_img_path = f"{sturm_folder}/S2/{sturm_code}_{date_fmt}_imagery.tif"
            
            try:
                p_bucket, p_key = s2_fm_path.replace('s3://', '').split('/', 1)
                s3_client.head_object(Bucket=p_bucket, Key=p_key)
                sturm_files[date_orig]['s2_floodmask'] = s2_fm_path
                print(f"  ✓ Found S2 mask: {date_fmt}")
            except:
                print(f"  ⚠️  Expected S2 mask missing: {date_fmt}")
            
            try:
                p_bucket, p_key = s2_img_path.replace('s3://', '').split('/', 1)
                s3_client.head_object(Bucket=p_bucket, Key=p_key)
                sturm_files[date_orig]['s2_imagery'] = s2_img_path
                print(f"  ✓ Found S2 imagery: {date_fmt}")
            except:
                pass
    
    # ================================================================
    # MAIN PROCESSING: Compare predictions to STURM
    # ================================================================
    
    all_results = {}
    all_metrics = []
    
    for date_str in date_list:
        if not sturm_files.get(date_str):
            print(f"\n⚠️  No STURM files for date {date_str}, skipping")
            continue
        
        print(f"\n{'='*70}")
        print(f"PROCESSING DATE: {date_str}")
        print(f"{'='*70}")
        
        date_output_folder = f"{assess_folder}/sturmComparison"
        
        for ws_idx, (ws_id, pred_tiff) in enumerate(watershed_predictions.items(), 1):
            print(f"\n[{ws_idx}/{len(watershed_predictions)}] Watershed: {ws_id}")
            
            # Extract prediction band for this date
            pred_array, ref_ds = extract_band_by_description(pred_tiff, date_str)
            
            if pred_array is None or ref_ds is None:
                print(f"  ⚠️  Could not load prediction for {date_str}")
                continue
            
            ref_geotransform = ref_ds.GetGeoTransform()
            ref_projection = ref_ds.GetProjection()
            H, W = pred_array.shape
            
            print(f"  Prediction shape: {H}x{W}")
            print(f"  Valid pixels: {(pred_array != pred_nodata).sum():,}")
            
            ws_output_folder = date_output_folder
            
            ws_results = {
                'watershed_id': ws_id,
                'date': date_str,
                'prediction_tiff': pred_tiff,
                'metrics': {},
                'confusion_matrices': {}
            }
            
            # Process each STURM product
            for sturm_key, sturm_path in sturm_files[date_str].items():
                sensor = 'S1' if 's1' in sturm_key else 'S2'
                is_floodmask = 'floodmask' in sturm_key
                is_imagery = 'imagery' in sturm_key
                
                if is_floodmask:
                    # Process floodmask: generate confusion matrices
                    print(f"  Processing STURM {sensor} floodmask...")
                    
                    aligned = warp_to_match(sturm_path, ref_ds,
                                           resample_method='nearest',
                                           src_nodata=255,
                                           dst_nodata=pred_nodata)
                    
                    if aligned is not None:
                        print(f"    Reprojected: {aligned.shape}")
                        print(f"    STURM valid pixels: {(aligned != pred_nodata).sum():,}")
                        
                        # Create validity mask
                        pred_valid = pred_array != pred_nodata
                        sturm_valid = (aligned != pred_nodata)
                        valid_overlap = pred_valid & sturm_valid
                        
                        print(f"    Prediction valid pixels: {pred_valid.sum():,}")
                        print(f"    STURM valid pixels: {sturm_valid.sum():,}")
                        print(f"    Overlapping valid pixels: {valid_overlap.sum():,}")
                        
                        if valid_overlap.sum() < 100:
                            print(f"    ⚠️  Insufficient overlap, skipping")
                            del aligned
                            gc.collect()
                            continue
                        
                        # Binarize: 0 stays 0 (land), 1+ becomes 1 (water)
                        sturm_binary = np.full_like(aligned, pred_nodata, dtype=np.uint8)
                        sturm_binary[sturm_valid] = np.where(aligned[sturm_valid] > 0, 1, 0)
                        
                        # Save aligned floodmask TIFF
                        fm_tiff = f"{ws_output_folder}/sturm_{sensor.lower()}_floodmask_{date_str}.tif"
                        save_tiff(sturm_binary, ref_geotransform, ref_projection, fm_tiff, pred_nodata)
                        
                        # Create floodmask PNG
                        fm_png = f"{ws_output_folder}/sturm_{sensor.lower()}_floodmask_{date_str}.png"
                        create_floodmask_png(sturm_binary, fm_png, sensor, ws_id, date_str)
                        
                        # Calculate metrics
                        metrics = calculate_binary_metrics(pred_array, sturm_binary,
                                                           pred_nodata, pred_nodata)
                        
                        if metrics.get('valid_pixels', 0) > 0:
                            print(f"    F1={metrics['f1']:.3f}, MCC={metrics['mcc']:.3f}")
                            
                            ws_results['metrics'][sensor] = metrics
                            
                            # Create confusion matrix
                            cm = create_confusion_matrix_array(pred_array, sturm_binary,
                                                               pred_nodata, pred_nodata)
                            
                            # Save confusion matrix TIFF
                            cm_tiff = f"{ws_output_folder}/confusion_matrix_{sensor}_{date_str}_{ws_id}.tif"
                            save_tiff(cm, ref_geotransform, ref_projection, cm_tiff, 255)
                            ws_results['confusion_matrices'][sensor] = cm_tiff
                            
                            # Create confusion matrix PNG
                            cm_png = f"{ws_output_folder}/confusion_matrix_{sensor}_{date_str}_{ws_id}.png"
                            create_confusion_matrix_png(cm, metrics, cm_png, sensor, ws_id, date_str)
                            
                            # Collect for aggregation
                            metrics.update({
                                'watershed': ws_id,
                                'date': date_str,
                                'sensor': sensor
                            })
                            all_metrics.append(metrics)
                            
                            del sturm_binary, cm
                            gc.collect()
                        
                        del aligned
                    else:
                        print(f"    ⚠️  Failed to reproject STURM {sensor} floodmask")
                
                elif is_imagery:
                    # Process imagery: just save aligned imagery
                    print(f"  Processing STURM {sensor} imagery...")
                    
                    aligned = warp_to_match(sturm_path, ref_ds,
                                           resample_method='bilinear',
                                           src_nodata=255,
                                           dst_nodata=pred_nodata)
                    
                    if aligned is not None:
                        print(f"    Reprojected: {aligned.shape}")
                        
                        # Save aligned imagery TIFF
                        img_tiff = f"{ws_output_folder}/sturm_{sensor.lower()}_imagery_{date_str}.tif"
                        save_tiff(aligned, ref_geotransform, ref_projection, img_tiff, pred_nodata)
                        
                        # Create imagery PNG
                        img_png = f"{ws_output_folder}/sturm_{sensor.lower()}_imagery_{date_str}.png"
                        create_imagery_png(aligned, img_png, sensor, ws_id, date_str)
                        
                        del aligned
                        gc.collect()
                    else:
                        print(f"    ⚠️  Failed to reproject STURM {sensor} imagery")
                
                gc.collect()
            
            # Cleanup after each watershed
            ref_ds = None
            del pred_array
            gc.collect()
            
            all_results[f"{ws_id}_{date_str}"] = ws_results
            print(f"  ✅ Completed {ws_id}")
            
            # Memory check
            try:
                import psutil
                process = psutil.Process()
                mem_mb = process.memory_info().rss / 1024 / 1024
                print(f"  Memory usage: {mem_mb:.0f} MB")
                
                if mem_mb > 8000:
                    print(f"  ⚠️  High memory usage, forcing cleanup...")
                    gc.collect()
                    gc.collect()
            except:
                pass
    
    # ================================================================
    # AGGREGATED SUMMARY (only if STURM comparisons were made)
    # ================================================================
    
    if all_metrics:
        print(f"\n{'='*70}")
        print("CREATING AGGREGATED SUMMARY")
        print(f"{'='*70}")
        
        df_metrics = pd.DataFrame(all_metrics)
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
            temp_csv = tmp.name
        df_metrics.to_csv(temp_csv, index=False)
        
        df_s1 = df_metrics[df_metrics['sensor'] == 'S1']
        df_s2 = df_metrics[df_metrics['sensor'] == 'S2']
        
        if not df_s1.empty:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
                temp_csv_s1 = tmp.name
            df_s1.to_csv(temp_csv_s1, index=False)
            
            summary_csv_s1 = f"{assess_folder}/sturmComparison/S1_per_date_metrics.csv"
            if summary_csv_s1.startswith('s3://'):
                bucket, key = summary_csv_s1.replace('s3://', '').split('/', 1)
                s3_client.upload_file(temp_csv_s1, bucket, key)
            os.remove(temp_csv_s1)
            print(f"  Created S1 metrics: {summary_csv_s1}")
        
        if not df_s2.empty:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
                temp_csv_s2 = tmp.name
            df_s2.to_csv(temp_csv_s2, index=False)
            
            summary_csv_s2 = f"{assess_folder}/sturmComparison/S2_per_date_metrics.csv"
            if summary_csv_s2.startswith('s3://'):
                bucket, key = summary_csv_s2.replace('s3://', '').split('/', 1)
                s3_client.upload_file(temp_csv_s2, bucket, key)
            os.remove(temp_csv_s2)
            print(f"  Created S2 metrics: {summary_csv_s2}")
        
        # Keep old combined CSV for backwards compatibility
        summary_csv = f"{assess_folder}/sturmComparison/all_metrics.csv"
        if summary_csv.startswith('s3://'):
            bucket, key = summary_csv.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_csv, bucket, key)
        
        os.remove(temp_csv)
        print(f"  Created: {summary_csv}")
        
        all_results['summary_csv'] = summary_csv
        all_results['summary_dataframe'] = df_metrics
        
        # ================================================================
        # STITCHED CONFUSION MATRICES (by sensor)
        # ================================================================
        print(f"\n{'='*70}")
        print("CREATING STITCHED CONFUSION MATRICES")
        print(f"{'='*70}")
        
        for date_str in date_list:
            print(f"\nProcessing date: {date_str}")
            
            for sensor_name in ['S1', 'S2']:
                print(f"\n  Processing {sensor_name} confusion matrices...")
                
                cm_tiffs_for_sensor = []
                sensor_metrics_list = []
                
                for ws_id, ws_results in all_results.items():
                    # Skip non-dictionary entries
                    if not isinstance(ws_results, dict):
                        continue
                    if ws_id in ['summary_csv', 'summary_dataframe']:
                        continue
                    if 'date' not in ws_results:
                        continue
                    
                    if ws_results['date'] == date_str and 'confusion_matrices' in ws_results:
                        if sensor_name in ws_results['confusion_matrices']:
                            cm_tiff = ws_results['confusion_matrices'][sensor_name]
                            cm_tiffs_for_sensor.append(cm_tiff)
                            
                            if 'metrics' in ws_results and sensor_name in ws_results['metrics']:
                                sensor_metrics_list.append(ws_results['metrics'][sensor_name])
                
                if not cm_tiffs_for_sensor:
                    print(f"    ⚠️  No {sensor_name} confusion matrices found for {date_str}")
                    continue
                
                print(f"    Found {len(cm_tiffs_for_sensor)} {sensor_name} confusion matrices")
                
                stitched_cm_tiff = f"{assess_folder}/sturmComparison/confusion_matrix_{sensor_name}_{date_str}.tif"
                
                if stitch_confusion_matrices(cm_tiffs_for_sensor, stitched_cm_tiff):
                    all_results[f'stitched_cm_tiff_{sensor_name}_{date_str}'] = stitched_cm_tiff
                    
                    aggregate_metrics = {}
                    
                    if sensor_metrics_list:
                        df_sensor = pd.DataFrame(sensor_metrics_list)
                        
                        for metric in ['f1', 'mcc', 'accuracy', 'precision', 'recall']:
                            if metric in df_sensor.columns:
                                aggregate_metrics[f'{metric}_mean'] = float(df_sensor[metric].mean())
                                aggregate_metrics[f'{metric}_std'] = float(df_sensor[metric].std())
                        
                        aggregate_metrics['n_watersheds'] = len(cm_tiffs_for_sensor)
                        aggregate_metrics['sensor'] = sensor_name
                    
                    stitched_cm_png = f"{assess_folder}/sturmComparison/confusion_matrix_{sensor_name}_{date_str}.png"
                    create_stitched_confusion_matrix_png(
                        stitched_cm_tiff, 
                        stitched_cm_png, 
                        date_str,
                        aggregate_metrics
                    )
                    
                    all_results[f'stitched_cm_png_{sensor_name}_{date_str}'] = stitched_cm_png
                    
                    gc.collect()
        
    
    # ================================================================
    # STITCHED PREDICTIONS (always generated, regardless of STURM)
    # ================================================================
    print(f"\n{'='*70}")
    print("CREATING STITCHED PREDICTIONS")
    print(f"{'='*70}")
    
    for date_str in date_list:
        print(f"\nProcessing predictions for date: {date_str}")
        
        # Extract date-specific bands from each watershed's multi-band TIFF
        temp_tiffs_for_date = []
        
        for ws_id, pred_tiff in watershed_predictions.items():
            # Open the multi-band prediction TIFF
            vsi_path = pred_tiff.replace('s3://', '/vsis3/')
            ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
            
            if ds is None:
                print(f"  ⚠️  Could not open prediction TIFF for {ws_id}")
                continue
            
            # Find the band for this date
            date_band_idx = None
            for b in range(1, ds.RasterCount + 1):
                band = ds.GetRasterBand(b)
                desc = band.GetDescription()
                if date_str in desc:
                    date_band_idx = b
                    break
            
            if date_band_idx is None:
                print(f"  ⚠️  Date {date_str} not found in {ws_id} prediction TIFF")
                ds = None
                continue
            
            # Extract this band to a temporary single-band TIFF
            with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
                temp_tiff = tmp.name
            
            # Create single-band TIFF
            driver = gdal.GetDriverByName('GTiff')
            out_ds = driver.Create(
                temp_tiff, 
                ds.RasterXSize, 
                ds.RasterYSize, 
                1,
                gdal.GDT_Byte,
                options=['COMPRESS=DEFLATE', 'TILED=YES']
            )
            
            out_ds.SetGeoTransform(ds.GetGeoTransform())
            out_ds.SetProjection(ds.GetProjection())
            
            # Copy the specific date band
            in_band = ds.GetRasterBand(date_band_idx)
            out_band = out_ds.GetRasterBand(1)
            out_band.SetNoDataValue(255)
            out_band.WriteArray(in_band.ReadAsArray())
            
            out_ds = None
            ds = None
            
            temp_tiffs_for_date.append(temp_tiff)
        
        if not temp_tiffs_for_date:
            print(f"  ⚠️  No prediction TIFFs found for {date_str}")
            continue
        
        # Stitch the single-band TIFFs
        stitched_pred_tiff = f"{assess_folder}/xgboost_predictions_stitched_{date_str}.tif"
        
        if stitch_prediction_tiffs(temp_tiffs_for_date, stitched_pred_tiff):
            all_results[f'stitched_pred_tiff_{date_str}'] = stitched_pred_tiff
            
            # Create PNG visualization
            stitched_pred_png = f"{assess_folder}/xgboost_predictions_stitched_{date_str}.png"
            create_stitched_prediction_png(
                stitched_pred_tiff,
                stitched_pred_png,
                date_str,
                len(watershed_predictions)
            )
            
            all_results[f'stitched_pred_png_{date_str}'] = stitched_pred_png
        
        # Clean up temporary files
        for temp_tiff in temp_tiffs_for_date:
            try:
                os.remove(temp_tiff)
            except:
                pass
        
        gc.collect()

    # ================================================================
    # CREATE STURM IMAGERY PNGs DIRECTLY FROM STITCHED FILES
    # ================================================================
    print(f"\n{'='*70}")
    print("CREATING STURM IMAGERY PNGs FROM STITCHED FILES")
    print(f"{'='*70}")
    
    # Get available stitched STURM files
    available_sturm = {}
    
    for sensor in ['S1', 'S2']:
        img_path = f"{sturm_folder}/sturm_{sturm_code}_{sensor.lower()}_imagery.tif"
        try:
            img_bucket, img_key = img_path.replace('s3://', '').split('/', 1)
            s3_client.head_object(Bucket=img_bucket, Key=img_key)
            available_sturm[sensor] = img_path
            print(f"  Found {sensor} imagery: {img_path}")
        except:
            print(f"  ⚠️  Missing {sensor} imagery")
    
    # Create PNGs for each date
    for date_str in date_list:
        print(f"\nCreating imagery PNGs for {date_str}...")
        
        for sensor, stitched_tiff in available_sturm.items():
            output_png = f"{assess_folder}/sturm_imagery_{sensor.lower()}_{date_str}.png"
            create_sturm_imagery_png_direct(stitched_tiff, output_png, sensor, date_str)
    
    # ================================================================
    # FINAL SUMMARY
    # ================================================================
    
    print(f"\n{'='*70}")
    print(f"✅ COMPARISON COMPLETE")
    print(f"{'='*70}")
    print(f"Watersheds processed: {len(watershed_predictions)}")
    print(f"Total comparisons: {len(all_metrics)}")
    
    if all_metrics:
        print(f"\nOutputs created:")
        print(f"  - Summary CSV: {all_results.get('summary_csv')}")
        for date_str in date_list:
            # Confusion matrices by sensor
            for sensor_name in ['S1', 'S2']:
                if f'stitched_cm_tiff_{sensor_name}_{date_str}' in all_results:
                    print(f"  - Stitched {sensor_name} CM TIFF ({date_str}): {all_results[f'stitched_cm_tiff_{sensor_name}_{date_str}']}")
                    print(f"  - Stitched {sensor_name} CM PNG ({date_str}): {all_results[f'stitched_cm_png_{sensor_name}_{date_str}']}")
            
            # Imagery by sensor
            for sensor_name in ['S1', 'S2']:
                if f'stitched_imagery_tiff_{sensor_name}_{date_str}' in all_results:
                    print(f"  - Stitched {sensor_name} Imagery TIFF ({date_str}): {all_results[f'stitched_imagery_tiff_{sensor_name}_{date_str}']}")
                    print(f"  - Stitched {sensor_name} Imagery PNG ({date_str}): {all_results[f'stitched_imagery_png_{sensor_name}_{date_str}']}")
            
            # Predictions
            if f'stitched_pred_tiff_{date_str}' in all_results:
                print(f"  - Stitched Predictions TIFF ({date_str}): {all_results[f'stitched_pred_tiff_{date_str}']}")
                print(f"  - Stitched Predictions PNG ({date_str}): {all_results[f'stitched_pred_png_{date_str}']}")
    
    return all_results
