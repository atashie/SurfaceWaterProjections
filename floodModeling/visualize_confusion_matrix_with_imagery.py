###############################################################
###############################################################
###############################################################
### Confusion Matrix Visualization - Metrics & Layout Update
###############################################################
###############################################################
###############################################################

"""
visualize_confusion_matrix_with_imagery.py

Generates a 3x2 diagnostic figure:
Row 1: S2 True Color | S2 NDWI | S2 Scene Classification
Row 2: XGB-OPERA CM | OPERA-S2 CM | XGB-S2 CM

Features:
- Metrics (MCC, F1, ROC, PR) calculated and displayed
- Legends moved to bottom-left inside axes
- Font sizes increased (Titles +20%, Text +50%)
- Robust RGB/Raw detection
- Sentinel Hub standard colormaps
"""

import os
import tempfile
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from matplotlib.patches import Patch
from osgeo import gdal, ogr, osr
import boto3
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score, average_precision_score

gdal.UseExceptions()


# ======================== Configuration ========================

LAT = "40.9437"
LON = "24.7669"
LOCATION_NAME = "EMSR292"
DATE_STR = "2024-11-20"
DATE_BAND = "20241120"

BASE_S3_PREFIX = "climate-ai-data-science-datasets/arrakis-data/floodOutputs"
LOCATION_FOLDER = f"trainingFolderFor_Lat_{LAT}_Lon_{LON}"

CONFUSION_MATRIX_TIFF = (
    f"s3://{BASE_S3_PREFIX}/{LOCATION_FOLDER}/"
    f"assessXGBoostModel/xgboost_confusion_matrix_Lat_{LAT}_Lon_{LON}.tif"
)

S2_FOLDER = f"s3://{BASE_S3_PREFIX}/{LOCATION_FOLDER}/s2_{DATE_STR}/"

S2_TRUE_COLOR = f"{S2_FOLDER}{DATE_STR}-00_00_{DATE_STR}-23_59_Sentinel-2_L2A_True_color.tiff"
S2_NDWI = f"{S2_FOLDER}{DATE_STR}-00_00_{DATE_STR}-23_59_Sentinel-2_L2A_NDWI.tiff"
S2_SCL = f"{S2_FOLDER}{DATE_STR}-00_00_{DATE_STR}-23_59_Sentinel-2_L2A_Scene_classification_map_.tiff"

WATERSHED_GPKG = (
    f"s3://{BASE_S3_PREFIX}/{LOCATION_FOLDER}/"
    f"watershedBoundaries/processedWatershedAndBasins_Lon{LON}_Lat{LAT}.gpkg"
)

OUTPUT_PATH = f"confusion_matrix_visualization_{LOCATION_NAME}_{DATE_STR}.png"


# ======================== Data Loading & Diagnostics ========================

def analyze_band_data(data: np.ndarray, name: str):
    """Log detailed statistics about the data array."""
    print(f"   [{name}] Analysis:")
    print(f"     - Shape: {data.shape}")
    print(f"     - Dtype: {data.dtype}")
    
    if data.ndim == 2:
        mn, mx = np.nanmin(data), np.nanmax(data)
        mean = np.nanmean(data)
        print(f"     - Range: [{mn:.4f}, {mx:.4f}]")
        print(f"     - Mean: {mean:.4f}")
        
        uniques = np.unique(data[~np.isnan(data)])
        if len(uniques) < 25:
            print(f"     - Unique Values: {uniques}")
            
    elif data.ndim == 3:
        print(f"     - Bands: {data.shape[0]}")
        if data.shape[0] >= 3:
            b1, b2, b3 = data[0], data[1], data[2]
            if np.array_equal(b1, b2) and np.array_equal(b2, b3):
                print("     - ⚠️  WARNING: Bands 1, 2, and 3 are IDENTICAL (Grayscale as RGB).")
            else:
                print("     - Bands are distinct (True Color/Multispectral).")


def load_tiff_from_s3_via_download(s3_path: str, target_extent: list = None, target_shape: tuple = None) -> tuple:
    """Download and load TIFF. Resamples if target_shape provided."""
    s3_client = boto3.client('s3')
    bucket, key = s3_path.replace('s3://', '').split('/', 1)
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        print(f"      Downloading {os.path.basename(s3_path)}...")
        s3_client.download_file(bucket, key, temp_path)
        
        ds = gdal.Open(temp_path, gdal.GA_ReadOnly)
        if ds is None:
            print(f"      ✗ Failed to open downloaded file")
            os.remove(temp_path)
            return None, None, None
        
        geotransform = ds.GetGeoTransform()
        width = ds.RasterXSize
        height = ds.RasterYSize
        n_bands = ds.RasterCount
        
        left = geotransform[0]
        top = geotransform[3]
        right = left + width * geotransform[1]
        bottom = top + height * geotransform[5]
        full_extent = [left, right, bottom, top]
        
        if target_extent is not None:
            target_left, target_right, target_bottom, target_top = target_extent
            x_off = int((target_left - left) / geotransform[1])
            y_off = int((target_top - top) / geotransform[5])
            x_size = int((target_right - target_left) / geotransform[1])
            y_size = int((target_bottom - target_top) / geotransform[5])
            
            x_off = max(0, min(x_off, width - 1))
            y_off = max(0, min(y_off, height - 1))
            x_size = max(1, min(x_size, width - x_off))
            y_size = max(1, min(y_size, height - y_off))
            
            actual_extent = [
                left + x_off * geotransform[1],
                left + (x_off + x_size) * geotransform[1],
                top + (y_off + y_size) * geotransform[5],
                top + y_off * geotransform[5]
            ]
        else:
            x_off, y_off, x_size, y_size = 0, 0, width, height
            actual_extent = full_extent

        if target_shape is not None:
            buf_ysize, buf_xsize = target_shape
        else:
            buf_xsize, buf_ysize = x_size, y_size

        data = ds.ReadAsArray(x_off, y_off, x_size, y_size, buf_xsize=buf_xsize, buf_ysize=buf_ysize)
        
        if n_bands >= 3:
            if np.array_equal(data[0], data[1]) and np.array_equal(data[1], data[2]):
                print("      -> Detected identical bands, collapsing to single band.")
                data = data[0]
        
        ds = None
        os.remove(temp_path)
        
        return data, actual_extent, geotransform
        
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        print(f"      ✗ Error: {e}")
        return None, None, None


def load_tiff_band_by_name(s3_path: str, band_name_pattern: str) -> tuple:
    """Load specific band from TIFF by matching band description."""
    vsi_path = s3_path.replace('s3://', '/vsis3/')
    
    try:
        ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
        if ds is None:
            print(f"✗ Failed to open: {s3_path}")
            return None, None, None, None
        
        for b in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(b)
            desc = band.GetDescription()
            
            if band_name_pattern in desc:
                print(f"✓ Found band {b}: {desc}")
                data = band.ReadAsArray()
                geotransform = ds.GetGeoTransform()
                projection = ds.GetProjection()
                
                width = ds.RasterXSize
                height = ds.RasterYSize
                left = geotransform[0]
                top = geotransform[3]
                right = left + width * geotransform[1]
                bottom = top + height * geotransform[5]
                extent = [left, right, bottom, top]
                
                ds = None
                return data, geotransform, projection, extent
        
        print(f"✗ Band matching '{band_name_pattern}' not found")
        ds = None
        return None, None, None, None
        
    except Exception as e:
        print(f"✗ Error loading {s3_path}: {e}")
        return None, None, None, None


def load_watershed_boundary(s3_path: str) -> gpd.GeoDataFrame:
    """Load watershed boundary from GPKG (first row = subwatershed)."""
    s3_client = boto3.client('s3')
    bucket, key = s3_path.replace('s3://', '').split('/', 1)
    
    with tempfile.NamedTemporaryFile(suffix='.gpkg', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        s3_client.download_file(bucket, key, temp_path)
        gdf = gpd.read_file(temp_path)
        
        if len(gdf) >= 1:
            watershed = gdf.iloc[[0]].copy()
            os.remove(temp_path)
            print(f"✓ Loaded watershed boundary from row 1 (subwatershed)")
            return watershed
        else:
            print(f"✗ GPKG is empty")
            os.remove(temp_path)
            return None
        
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        print(f"✗ Error loading watershed: {e}")
        return None


def rasterize_watershed_mask(watershed_gdf: gpd.GeoDataFrame, 
                             geotransform: tuple, 
                             width: int, height: int) -> np.ndarray:
    """Rasterize watershed polygon to create boolean mask."""
    
    mem_driver = gdal.GetDriverByName('MEM')
    mem_ds = mem_driver.Create('', width, height, 1, gdal.GDT_Byte)
    mem_ds.SetGeoTransform(geotransform)
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    mem_ds.SetProjection(srs.ExportToWkt())
    
    mem_layer_ds = ogr.GetDriverByName('Memory').CreateDataSource('')
    mem_layer = mem_layer_ds.CreateLayer('watershed', srs, ogr.wkbPolygon)
    
    for geom in watershed_gdf.geometry:
        feature = ogr.Feature(mem_layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkb(geom.wkb))
        mem_layer.CreateFeature(feature)
    
    gdal.RasterizeLayer(mem_ds, [1], mem_layer, burn_values=[1])
    
    mask = mem_ds.GetRasterBand(1).ReadAsArray()
    
    mem_ds = None
    mem_layer_ds = None
    
    return mask.astype(bool)


# ======================== Colormap Definitions ========================

def get_scl_colormap():
    """Sentinel Hub Scene Classification colormap."""
    colors_rgb = [
        [0, 0, 0],           # 0: No Data
        [255, 0, 0],         # 1: Saturated
        [47, 47, 47],        # 2: Dark Features
        [100, 50, 0],        # 3: Cloud Shadows
        [0, 160, 0],         # 4: Vegetation
        [255, 230, 90],      # 5: Not Vegetated
        [0, 0, 255],         # 6: Water
        [128, 128, 128],     # 7: Unclassified
        [192, 192, 192],     # 8: Cloud Med
        [255, 255, 255],     # 9: Cloud High
        [100, 200, 255],     # 10: Cirrus
        [255, 150, 255],     # 11: Snow
    ]
    colors_norm = [[c[0]/255, c[1]/255, c[2]/255] for c in colors_rgb]
    cmap = ListedColormap(colors_norm)
    norm = BoundaryNorm(range(13), cmap.N)
    labels = [
        'No Data', 'Saturated', 'Dark Features', 'Cloud Shadows',
        'Vegetation', 'Not Vegetated', 'Water', 'Unclassified',
        'Cloud Med', 'Cloud High', 'Cirrus', 'Snow/Ice'
    ]
    return cmap, norm, labels


def get_ndwi_colormap():
    """Sentinel Hub NDWI colormap."""
    colors = [
        (0.0, '#FF0000'),   # -1.0: Red (Dry)
        (0.4, '#FFFF00'),   # -0.2: Yellow
        (0.5, '#FFFFFF'),   #  0.0: White
        (0.6, '#00FFFF'),   #  0.2: Cyan
        (1.0, '#0000FF'),   #  1.0: Blue (Water)
    ]
    cmap = LinearSegmentedColormap.from_list('ndwi_custom', colors)
    return cmap


def get_confusion_matrix_colormap():
    """RandomForest confusion matrix encoding."""
    colors = ['green', 'blue', 'orange', 'red']
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(range(5), cmap.N)
    return cmap, norm


# ======================== Visualization ========================

def prepare_rgb_image(data: np.ndarray) -> np.ndarray:
    """Prepare 3-band data for RGB display."""
    if data.ndim != 3:
        return data
    if data.shape[0] == 3:
        data = np.transpose(data, (1, 2, 0))
    if data.max() > 1.0:
        data = data.astype(float) / 255.0
    return np.clip(data, 0, 1)


def calculate_metrics_from_cm(cm_array: np.ndarray) -> str:
    """
    Calculate metrics from confusion matrix array.
    Codes: 0=TN, 1=TP, 2=FN, 3=FP, 255=Mask
    """
    valid_mask = cm_array != 255
    vals = cm_array[valid_mask]
    
    if len(vals) == 0:
        return "No Data"
    
    # Reconstruct y_true and y_pred
    # TP(1) or FN(2) -> Actual Water (1)
    y_true = np.isin(vals, [1, 2]).astype(int)
    # TP(1) or FP(3) -> Predicted Water (1)
    y_pred = np.isin(vals, [1, 3]).astype(int)
    
    if len(np.unique(y_true)) < 2:
        return "Single Class"
        
    mcc = matthews_corrcoef(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    roc = roc_auc_score(y_true, y_pred)
    pr = average_precision_score(y_true, y_pred)
    
    return f"MCC: {mcc:.2f}\nF1: {f1:.2f}\nROC: {roc:.2f}\nPR: {pr:.2f}"


def calculate_derived_confusion_matrix(
    base_cm: np.ndarray, 
    scl: np.ndarray, 
    pred_source: str
) -> np.ndarray:
    """
    Calculate new confusion matrix using SCL as Ground Truth.
    """
    
    # --- 1. Define Ground Truth (S2) ---
    if scl.ndim == 3:
        # RGB Logic (Sentinel Hub Colors)
        R = scl[0]
        G = scl[1]
        B = scl[2]
        
        # Water: Blue channel is dominant
        s2_water = (B > R) & (B > G)
        
        # Land: Green or Red channel is dominant (covers Green, Yellow, Brown)
        # Also exclude No Data (Black) and Clouds (White/Grey)
        is_not_blue_dominant = ~s2_water
        is_not_black = (R > 0) | (G > 0) | (B > 0)
        
        # Simple Cloud Mask: High brightness, low saturation
        is_bright = (R > 200) & (G > 200) & (B > 200)
        is_grey = (np.abs(R - G) < 20) & (np.abs(G - B) < 20) & (R > 100)
        is_cloud = is_bright | is_grey
        
        s2_land = is_not_blue_dominant & is_not_black & (~is_cloud)
        s2_valid = s2_water | s2_land
        
    else:
        # 2D Class ID Logic
        s2_water = (scl == 6)
        s2_land = np.isin(scl, [2, 3, 4, 5, 7, 11])
        s2_valid = s2_water | s2_land
    
    # --- 2. Define Prediction from Base CM ---
    if pred_source == 'OPERA':
        pred_water = (base_cm == 1) | (base_cm == 2)
        pred_land = (base_cm == 0) | (base_cm == 3)
    elif pred_source == 'XGBoost':
        pred_water = (base_cm == 1) | (base_cm == 3)
        pred_land = (base_cm == 0) | (base_cm == 2)
    else:
        return None

    # --- 3. Compute New CM ---
    new_cm = np.full_like(base_cm, 255)
    valid_mask = s2_valid & (base_cm != 255) & ((pred_water | pred_land))
    
    new_cm[valid_mask & s2_land & pred_land] = 0  # TN
    new_cm[valid_mask & s2_water & pred_water] = 1  # TP
    new_cm[valid_mask & s2_water & pred_land] = 2  # FN
    new_cm[valid_mask & s2_land & pred_water] = 3  # FP
    
    return new_cm


def create_confusion_matrix_figure(
    confusion_matrix: np.ndarray,
    cm_extent: list,
    cm_geotransform: tuple,
    true_color: np.ndarray,
    tc_extent: list,
    ndwi: np.ndarray,
    ndwi_extent: list,
    scl: np.ndarray,
    scl_extent: list,
    watershed: gpd.GeoDataFrame,
    output_path: str
):
    """Create 3x2 figure with satellite imagery and 3 confusion matrices."""
    
    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    
    # Create watershed mask for original CM
    height, width = confusion_matrix.shape
    ws_mask = rasterize_watershed_mask(watershed, cm_geotransform, width, height)
    
    # Apply mask to original CM
    cm_masked = confusion_matrix.copy()
    cm_masked[~ws_mask] = 255
    
    # --- Calculate Derived CMs ---
    scl_calc = scl if scl.ndim == 2 else scl[0]
    cm_opera_s2 = calculate_derived_confusion_matrix(cm_masked, scl, 'OPERA')
    cm_xgb_s2 = calculate_derived_confusion_matrix(cm_masked, scl, 'XGBoost')
    
    # --- Plotting Helper ---
    def plot_cm(ax, data, title):
        cm_disp = data.astype(float)
        cm_disp[data == 255] = np.nan
        cmap, norm = get_confusion_matrix_colormap()
        ax.imshow(cm_disp, extent=cm_extent, origin='upper', cmap=cmap, norm=norm, interpolation='nearest')
        
        # Legend (Bottom Left)
        legend_elements = [
            Patch(facecolor='green', edgecolor='black', label='TN'),
            Patch(facecolor='blue', edgecolor='black', label='TP'),
            Patch(facecolor='orange', edgecolor='black', label='FN'),
            Patch(facecolor='red', edgecolor='black', label='FP'),
        ]
        ax.legend(handles=legend_elements, loc='lower left', fontsize=16, frameon=True, framealpha=0.8)
        
        # Metrics (Top Right)
        metrics_text = calculate_metrics_from_cm(data)
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.96, 0.96, metrics_text, transform=ax.transAxes, fontsize=18,
                verticalalignment='top', horizontalalignment='right', bbox=props)
        
        if watershed is not None:
            watershed.boundary.plot(ax=ax, color='#2a2a2a', linewidth=2, zorder=2)
        
        ax.set_title(title, fontsize=22, fontweight='bold')
        ax.set_xlabel('Longitude', fontsize=21)
        ax.set_ylabel('Latitude', fontsize=21)
        ax.tick_params(labelsize=22)
        ax.grid(alpha=0.3)

    # ========== Row 1: Imagery ==========
    
    # 1. True Color
    ax = axes[0, 0]
    if true_color is not None:
        analyze_band_data(true_color, "True Color")
        tc_rgb = prepare_rgb_image(true_color)
        ax.imshow(tc_rgb, extent=tc_extent, origin='upper')
    if watershed is not None:
        watershed.boundary.plot(ax=ax, color='#2a2a2a', linewidth=2, zorder=2)
    ax.set_title('S2 True Color', fontsize=22, fontweight='bold')
    ax.set_xlabel('Longitude', fontsize=21)
    ax.set_ylabel('Latitude', fontsize=21)
    ax.tick_params(labelsize=22)
    ax.grid(alpha=0.3)
    
    # 2. NDWI
    ax = axes[0, 1]
    if ndwi is not None:
        analyze_band_data(ndwi, "NDWI")
        if ndwi.ndim == 3:
            ax.imshow(prepare_rgb_image(ndwi), extent=ndwi_extent, origin='upper')
        else:
            cmap = get_ndwi_colormap()
            im = ax.imshow(ndwi, extent=ndwi_extent, origin='upper', cmap=cmap, vmin=-1, vmax=1)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, label='NDWI')
            cbar.ax.tick_params(labelsize=16)
            cbar.set_label('NDWI', fontsize=18)
    if watershed is not None:
        watershed.boundary.plot(ax=ax, color='#2a2a2a', linewidth=2, zorder=2)
    ax.set_title('S2 NDWI', fontsize=22, fontweight='bold')
    ax.set_xlabel('Longitude', fontsize=21)
    ax.set_ylabel('Latitude', fontsize=21)
    ax.tick_params(labelsize=22)
    ax.grid(alpha=0.3)
    
    # 3. Scene Classification
    ax = axes[0, 2]
    if scl is not None:
        analyze_band_data(scl, "SCL")
        if scl.ndim == 3:
            ax.imshow(prepare_rgb_image(scl), extent=scl_extent, origin='upper')
        else:
            cmap, norm, labels = get_scl_colormap()
            im = ax.imshow(scl, extent=scl_extent, origin='upper', cmap=cmap, norm=norm, interpolation='nearest')
            legend_elements = [Patch(facecolor=cmap(i), label=f'{i}: {labels[i]}') for i in range(len(labels))]
            ax.legend(handles=legend_elements, loc='lower left', fontsize=12, frameon=True, framealpha=0.8)
    if watershed is not None:
        watershed.boundary.plot(ax=ax, color='#2a2a2a', linewidth=2, zorder=2)
    ax.set_title('S2 Classification', fontsize=22, fontweight='bold')
    ax.set_xlabel('Longitude', fontsize=21)
    ax.set_ylabel('Latitude', fontsize=21)
    ax.tick_params(labelsize=22)
    ax.grid(alpha=0.3)
    
    # ========== Row 2: Confusion Matrices ==========
    
    # 4. XGBoost - OPERA (Original)
    plot_cm(axes[1, 0], cm_masked, 'XGBoost -> OPERA')
    
    # 5. OPERA - S2
    plot_cm(axes[1, 1], cm_opera_s2, 'OPERA -> S2')
    
    # 6. XGBoost - S2
    plot_cm(axes[1, 2], cm_xgb_s2, 'XGBoost -> S2')
    
    # ========== Overall Title ==========
    fig.suptitle(
        f'Flood Detection Analysis: {LOCATION_NAME} - {DATE_STR}',
        fontsize=29, fontweight='bold', y=0.98
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved figure: {output_path}")


# ======================== Main Workflow ========================

def generate_confusion_matrix_visualization():
    """Main workflow."""
    print("=" * 70)
    print("CONFUSION MATRIX VISUALIZATION (3x2 Grid)")
    print(f"Location: {LOCATION_NAME} ({LAT}, {LON})")
    print(f"Date: {DATE_STR}")
    print("=" * 70)
    
    # 1. Load confusion matrix (Reference)
    print("\n1. Loading confusion matrix...")
    cm_data, cm_geotransform, cm_projection, cm_extent = load_tiff_band_by_name(
        CONFUSION_MATRIX_TIFF, f"XGBoost_CM_{DATE_BAND}"
    )
    
    if cm_data is None:
        print("✗ Failed to load confusion matrix")
        return
    
    target_shape = cm_data.shape  # (height, width)
    print(f"   Target Shape: {target_shape}")
    
    # 2. Load Sentinel-2 imagery (Resampled to match CM)
    print("\n2. Loading Sentinel-2 imagery...")
    
    print("   - True Color...")
    tc_data, tc_extent, _ = load_tiff_from_s3_via_download(
        S2_TRUE_COLOR, target_extent=cm_extent, target_shape=target_shape
    )
    
    print("   - NDWI...")
    ndwi_data, ndwi_extent, _ = load_tiff_from_s3_via_download(
        S2_NDWI, target_extent=cm_extent, target_shape=target_shape
    )
    
    print("   - Scene Classification...")
    scl_data, scl_extent, _ = load_tiff_from_s3_via_download(
        S2_SCL, target_extent=cm_extent, target_shape=target_shape
    )
    
    # 3. Load watershed boundary
    print("\n3. Loading watershed boundary...")
    watershed = load_watershed_boundary(WATERSHED_GPKG)
    
    # 4. Create visualization
    print("\n4. Creating visualization...")
    create_confusion_matrix_figure(
        confusion_matrix=cm_data,
        cm_extent=cm_extent,
        cm_geotransform=cm_geotransform,
        true_color=tc_data,
        tc_extent=tc_extent,
        ndwi=ndwi_data,
        ndwi_extent=ndwi_extent,
        scl=scl_data,
        scl_extent=scl_extent,
        watershed=watershed,
        output_path=OUTPUT_PATH
    )
    
    print("\n" + "=" * 70)
    print("✅ VISUALIZATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    generate_confusion_matrix_visualization()