"""
sturmFloodStitcher.py - STURM Flood Dataset Stitcher (With Bbox Filtering & NoData Validation)

Stitches STURM flood data with:
- Bbox filtering to clip to predictions area
- Validation to skip mosaics containing only NoData
"""

import os
import gc
import tempfile
import io
import numpy as np
from typing import Dict, List, Optional, Tuple
import pandas as pd
from osgeo import gdal, osr
import boto3

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
    gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,tiff,csv")
    gdal.SetConfigOption("AWS_NO_SIGN_REQUEST", "NO")
    gdal.SetConfigOption("VSI_CACHE", "TRUE")

configure_gdal()


# ======================== Validation Utilities ========================

def raster_contains_valid_data(file_path: str) -> bool:
    """
    Check if a raster contains any valid (non-NoData) data.
    Returns True if valid data exists, False if all NoData/null.
    """
    try:
        ds = gdal.Open(file_path, gdal.GA_ReadOnly)
        if not ds:
            return False
        
        # Check first band (sufficient for our purposes)
        band = ds.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        
        # Read entire band to check for valid data
        data = band.ReadAsArray()
        
        if data is None:
            ds = None
            return False
        
        # Count valid pixels
        if nodata is not None:
            valid_pixels = np.count_nonzero(data != nodata)
        else:
            # If no nodata value defined, check for non-zero values
            valid_pixels = np.count_nonzero(data != 0)
        
        ds = None
        
        return valid_pixels > 0
        
    except Exception as e:
        print(f"        ⚠️ Error checking raster data: {e}")
        return False


# ======================== Metadata Loading ========================

def load_metadata_from_s3(bucket: str, key: str) -> pd.DataFrame:
    """Load a CSV from S3 into a Pandas DataFrame."""
    s3_client = boto3.client('s3')
    try:
        print(f"  Loading metadata: s3://{bucket}/{key}")
        obj = s3_client.get_object(Bucket=bucket, Key=key)
        df = pd.read_csv(io.BytesIO(obj['Body'].read()))
        df.columns = [c.lower() for c in df.columns]
        return df
    except Exception as e:
        print(f"  ⚠️ Error loading metadata from s3://{bucket}/{key}: {e}")
        return pd.DataFrame()


def get_tiles_from_metadata(ems_code: str, metadata_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Filter metadata for a specific EMSR code and group tile_ids by date."""
    if metadata_df.empty:
        return {}

    code_col = 'ems_code' if 'ems_code' in metadata_df.columns else 'emsr_code'
    if code_col not in metadata_df.columns:
        print(f"  ⚠️ Metadata missing '{code_col}' column.")
        return {}

    event_df = metadata_df[metadata_df[code_col] == ems_code].copy()
    if event_df.empty:
        print(f"  ⚠️ No entries found for event code: {ems_code}")
        return {}

    date_col = None
    if 'floodmap_date' in event_df.columns:
        date_col = 'floodmap_date'
    elif 'sentinel_timestamp' in event_df.columns:
        date_col = 'sentinel_timestamp'
    elif 'sentinel_date' in event_df.columns:
        date_col = 'sentinel_date'
    
    if not date_col:
        print(f"  ⚠️ Metadata missing date column.")
        return {}

    try:
        event_df['date_str'] = pd.to_datetime(event_df[date_col]).dt.strftime('%Y-%m-%d')
    except Exception as e:
        print(f"  ⚠️ Error parsing dates: {e}")
        return {}
    
    tiles_by_date = {}
    if 'tile_id' not in event_df.columns:
        print("  ⚠️ Metadata missing 'tile_id' column.")
        return {}

    for date_str, group in event_df.groupby('date_str'):
        tiles = group['tile_id'].unique().tolist()
        tiles_by_date[date_str] = tiles
        
    return tiles_by_date


# ======================== Spatial Filtering ========================

def get_bbox_from_tiff(tiff_uri: str, target_crs: str = "EPSG:4326") -> Optional[Tuple[float, float, float, float]]:
    """Extract bounding box from TIFF."""
    vsis3_path = tiff_uri.replace('s3://', '/vsis3/')
    try:
        ds = gdal.Open(vsis3_path, gdal.GA_ReadOnly)
        if not ds: return None
        
        gt = ds.GetGeoTransform()
        width = ds.RasterXSize
        height = ds.RasterYSize
        
        minx = gt[0]
        maxy = gt[3]
        maxx = minx + (gt[1] * width)
        miny = maxy + (gt[5] * height)
        
        return (minx, miny, maxx, maxy)
    except:
        return None


# ======================== Stitching with Validation ========================

def stitch_date_mosaic_with_validation(
    tile_paths: List[str], 
    output_uri: str,
    is_mask: bool,
    target_crs: str = "EPSG:4326",
    predictions_bbox: Optional[Tuple[float, float, float, float]] = None,
    bbox_buffer: float = 0.0
) -> bool:
    """
    Stitch tiles with bbox filtering and validate output contains data.
    Only saves file if it contains valid (non-NoData) pixels.
    """
    if not tile_paths:
        return False
    
    print(f"      Mosaicking {len(tile_paths)} tiles with bbox filtering...")
    
    vsis3_paths = [p.replace('s3://', '/vsis3/') for p in tile_paths]
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_output = tmp.name
    
    try:
        # Setup output bounds for clipping
        output_bounds = None
        if predictions_bbox:
            minx, miny, maxx, maxy = predictions_bbox
            output_bounds = [minx - bbox_buffer, miny - bbox_buffer, maxx + bbox_buffer, maxy + bbox_buffer]
            print(f"        Clipping to bbox: {output_bounds}")
        
        # Warp with bbox clipping
        warp_options = gdal.WarpOptions(
            format='GTiff',
            creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
            dstSRS=target_crs,
            resampleAlg='near' if is_mask else 'bilinear',
            outputType=gdal.GDT_Byte if is_mask else gdal.GDT_Float32,
            outputBounds=output_bounds,
            srcNodata=None,
            dstNodata=255 if is_mask else None,
            multithread=True,
            warpMemoryLimit=2048
        )
        
        result = gdal.Warp(temp_output, vsis3_paths, options=warp_options)
        
        if result is not None:
            result = None
            
            # Validate output contains data
            print(f"        Validating mosaic contains valid data...")
            has_data = raster_contains_valid_data(temp_output)
            
            if not has_data:
                print(f"        ⚠️ Mosaic contains only NoData/null values - skipping save")
                os.remove(temp_output)
                return False
            
            # Get mosaic info
            ds = gdal.Open(temp_output, gdal.GA_ReadOnly)
            if ds:
                width = ds.RasterXSize
                height = ds.RasterYSize
                bands = ds.RasterCount
                print(f"        ✓ Valid mosaic created: {width}x{height}, {bands} bands")
                ds = None
            
            # Upload to S3
            s3_client = boto3.client('s3')
            bucket, key = output_uri.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_output, bucket, key)
            
            os.remove(temp_output)
            return True
        else:
            print(f"        ⚠️ Warp failed")
            if os.path.exists(temp_output):
                os.remove(temp_output)
            return False
            
    except Exception as e:
        print(f"        ⚠️ Error: {e}")
        if os.path.exists(temp_output):
            os.remove(temp_output)
        return False


# ======================== Main Function ========================

def stitch_sturm_floods(
    flood_event_code: str,
    predictions_tiff_uri: str,
    output_folder: str,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs/curated_flood_maps",
    target_crs: str = "EPSG:4326",
    bbox_buffer: float = 0.1
) -> Dict:
    """
    Main function to stitch STURM flood data with bbox filtering.
    Only saves mosaics that contain valid data within the bbox.
    """
    
    print(f"\n{'='*70}")
    print(f"STITCHING STURM FLOODS: {flood_event_code}")
    print("With bbox filtering and NoData validation")
    print(f"{'='*70}")
    
    # 1. Get Reference BBox
    print(f"\n[1] Extracting reference bounding box...")
    ref_bbox = get_bbox_from_tiff(predictions_tiff_uri, target_crs)
    if ref_bbox:
        print(f"  ✓ Reference BBox: {ref_bbox}")
    else:
        print("  ⚠️ Reference BBox not found, stitching full extent.")

    # 2. Load Metadata
    print(f"\n[2] Loading Metadata...")
    s1_meta_key = f"{base_prefix}/Sentinel1_metadata.csv"
    s2_meta_key = f"{base_prefix}/Sentinel2_metadata.csv"
    
    df_s1 = load_metadata_from_s3(bucket, s1_meta_key)
    df_s2 = load_metadata_from_s3(bucket, s2_meta_key)
    
    # 3. Group Tiles by Date
    print(f"\n[3] Grouping Tiles by Date...")
    s1_tiles_by_date = get_tiles_from_metadata(flood_event_code, df_s1)
    s2_tiles_by_date = get_tiles_from_metadata(flood_event_code, df_s2)
    
    results = {'S1': {}, 'S2': {}}
    
    # 4. Process Sensors
    sturm_root = f"s3://{bucket}/{base_prefix}/STURM_Floods"
    
    sensor_configs = [
        {
            'name': 'S1',
            'tiles_map': s1_tiles_by_date,
            'mask_dir': f"{sturm_root}/Sentinel1/Floodmaps",
            'img_dir': f"{sturm_root}/Sentinel1/S1"
        },
        {
            'name': 'S2',
            'tiles_map': s2_tiles_by_date,
            'mask_dir': f"{sturm_root}/Sentinel2/Floodmaps",
            'img_dir': f"{sturm_root}/Sentinel2/S2"
        }
    ]
    
    for config in sensor_configs:
        s_name = config['name']
        tiles_map = config['tiles_map']
        
        print(f"\n[Processing {s_name}]")
        if not tiles_map:
            print(f"  ⚠️ No tiles found in metadata for {s_name} / {flood_event_code}")
            continue
            
        print(f"  Found {len(tiles_map)} unique dates.")
        
        for date_str, tile_ids in tiles_map.items():
            print(f"  Date: {date_str} ({len(tile_ids)} tiles)")
            
            # --- Stitch Masks ---
            mask_paths = [f"{config['mask_dir']}/{tid}" for tid in tile_ids]
            mask_out = f"{output_folder}/{s_name}/{flood_event_code}_{date_str}_mask.tif"
            
            print(f"    > Stitching Masks...")
            if stitch_date_mosaic_with_validation(mask_paths, mask_out, is_mask=True, 
                                                   target_crs=target_crs, predictions_bbox=ref_bbox, bbox_buffer=bbox_buffer):
                if date_str not in results[s_name]: results[s_name][date_str] = {}
                results[s_name][date_str]['mask'] = mask_out
            else:
                print(f"    ⚠️ Mask for {date_str} not saved (no data in bbox)")
            
            # --- Stitch Imagery ---
            img_paths = [f"{config['img_dir']}/{tid}" for tid in tile_ids]
            img_out = f"{output_folder}/{s_name}/{flood_event_code}_{date_str}_imagery.tif"
            
            print(f"    > Stitching Imagery...")
            if stitch_date_mosaic_with_validation(img_paths, img_out, is_mask=False, 
                                                   target_crs=target_crs, predictions_bbox=ref_bbox, bbox_buffer=bbox_buffer):
                if date_str not in results[s_name]: results[s_name][date_str] = {}
                results[s_name][date_str]['imagery'] = img_out
            else:
                print(f"    ⚠️ Imagery for {date_str} not saved (no data in bbox)")

    print(f"\n{'='*70}")
    print("STITCHING COMPLETE")
    print(f"{'='*70}")
    return results