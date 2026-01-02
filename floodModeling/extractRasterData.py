# extractRasterData.py

import json
import numpy as np
from osgeo import gdal
import logging
import tempfile
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from scipy import ndimage
import boto3
from typing import Dict, Optional, List, Tuple
import warnings
from scipy.ndimage import generic_filter
from pyproj import CRS
DSWX_DEBUG = False  # set True for detailed per-tile diagnostics

# Import ClimateAI utilities
from climateai.geolately.db_utils import Geospatial_DB
from climateai.geolately.raster_utils import (
    clip_raster_to_vector, 
    mosaic_and_fn, 
    resample_raster_to_raster
)
from coordinate_utils import format_coordinate
from delineateWatersheds import read_s3_gpkg

# Configure logging
logging.basicConfig(level=logging.INFO)
warnings.filterwarnings('ignore')

# Define soil feature products
SOIL_FEATURES = [
    'CLAY_0_TO_5CM',
#    'CLAY_5_TO_15CM',
#    'CLAY_15_TO_30CM', 
    'CLAY_30_TO_60CM',
#    'CLAY_60_TO_100CM',
#    'CLAY_100_TO_200CM',
#    'SAND_0_TO_5CM', 
    'SAND_5_TO_15CM', 
#    'SAND_15_TO_30CM',
#    'SAND_30_TO_60CM', 
    'SAND_60_TO_100CM', 
#    'SAND_100_TO_200CM',
#    'SILT_0_TO_5CM', 
#    'SILT_5_TO_15CM', 
    'SILT_15_TO_30CM',
#    'SILT_30_TO_60CM', 
#    'SILT_60_TO_100CM', 
    'SILT_100_TO_200CM'
]


def prepare_watershed_polygon(longitude: float, latitude: float, 
                            polygon_type: str = 'subwatershed',
                            training_folder_path: Optional[str] = None) -> Dict:
    """
    Read and process watershed polygons from pre-processed HydroBASINs data.
    """
    # Format coordinates immediately
    lat_str = format_coordinate(latitude)
    lon_str = format_coordinate(longitude)
    
    print(f"\n1. Loading HydroBASINs polygons for location ({lat_str}, {lon_str})...")
    
    gpkg_filename = f"processedWatershedAndBasins_Lon{lon_str}_Lat{lat_str}.gpkg"
    
    # Determine where to look for the file
    if training_folder_path:
        gpkg_path = f"{training_folder_path}/watershedBoundaries/{gpkg_filename}"
    else:
        gpkg_path = f"s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/watershedBoundaries/{gpkg_filename}"
    
    print(f"   Reading from: {gpkg_path}")
    
    try:
        watershed_gdf = read_s3_gpkg(gpkg_path)
        print(f"   ✅ Successfully loaded {len(watershed_gdf)} polygons")
    except Exception as e:
        error_msg = f"""
        ❌ Failed to load watershed file.
        Expected file: {gpkg_filename}
        Full path: {gpkg_path}
        
        Please ensure you're using the exact coordinates that were used when creating the file.
        Error: {str(e)}
        """
        raise FileNotFoundError(error_msg)
    
    print(f"\n2. Extracting {polygon_type} polygon...")
    polygon_gdf = watershed_gdf[watershed_gdf['feature_type'] == polygon_type].copy()
    
    if len(polygon_gdf) == 0:
        raise ValueError(f"No {polygon_type} found in the GeoPackage")
    
    polygon_row = polygon_gdf.iloc[0]
    metadata = {
        'HYBAS_ID': polygon_row.get('HYBAS_ID', 'N/A'),
        'SUB_AREA': polygon_row.get('SUB_AREA', None),
        'feature_type': polygon_type,
        'crs': str(polygon_gdf.crs),
        'longitude': lon_str,
        'latitude': lat_str
    }
    
    print(f"   ✅ {polygon_type.capitalize()} extracted")
    print(f"      HYBAS_ID: {metadata['HYBAS_ID']}")
    if metadata['SUB_AREA']:
        print(f"      Area: {metadata['SUB_AREA']:.2f} km²")
    
    print(f"\n3. Converting to GeoJSON...")
    
    # Robust CRS comparison and conversion
    try:
        crs_obj = polygon_gdf.crs
        if crs_obj is not None and crs_obj != CRS.from_epsg(4326):
            polygon_gdf = polygon_gdf.to_crs('EPSG:4326')
    except Exception:
        # Fallback to previous behavior if CRS objects are not comparable
        if str(polygon_gdf.crs) != 'EPSG:4326':
            polygon_gdf = polygon_gdf.to_crs('EPSG:4326')
    
    geojson_path = f"/tmp/{polygon_type}_{lon_str}_{lat_str}.geojson" 
    polygon_gdf.to_file(geojson_path, driver='GeoJSON')
    print(f"   ✅ GeoJSON saved to: {geojson_path}")
    
    return {
        'watershed_gdf': watershed_gdf,
        'polygon_gdf': polygon_gdf,
        'geojson_path': geojson_path,
        'metadata': metadata
    }

def align_array_to_reference(array: np.ndarray, ref_shape: Tuple[int, int], 
                            method: str = 'crop_or_pad') -> np.ndarray:
    """
    Align an array to match reference dimensions.
    
    Parameters:
    -----------
    array : np.ndarray
        Input array to align
    ref_shape : tuple
        Reference shape (height, width)
    method : str
        Method for alignment: 'crop_or_pad', 'resample'
        
    Returns:
    --------
    np.ndarray
        Aligned array matching ref_shape
    """
    ref_height, ref_width = ref_shape
    arr_height, arr_width = array.shape
    
    if method == 'crop_or_pad':
        # Create output array filled with nodata
        aligned = np.full(ref_shape, -9999, dtype=array.dtype)
        
        # Calculate overlap region
        min_height = min(arr_height, ref_height)
        min_width = min(arr_width, ref_width)
        
        # Copy valid data
        aligned[:min_height, :min_width] = array[:min_height, :min_width]
        
        return aligned
    
    elif method == 'resample':
        # Use scipy zoom for resampling
        from scipy.ndimage import zoom
        zoom_factors = (ref_height / arr_height, ref_width / arr_width)
        return zoom(array, zoom_factors, order=1)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def query_and_process_hand(geojson_path: str, polygon_gdf: gpd.GeoDataFrame, 
                          reference_dem: gdal.Dataset) -> Dict:
    """
    Query and process HAND (Height Above Nearest Drainage) data.
    HAND represents the vertical distance to the nearest drainage in meters.
    """
    print(f"\n3d. Querying and processing HAND data...")
    
    geo_db = Geospatial_DB(enable_s3=True)
    
    hand_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product='HAND',
        start_date='2000-01-01',
        end_date='2025-01-01'
    )
    
    if hand_results is None or len(hand_results) == 0:
        print(f"   ⚠️ No HAND data found for this region")
        return None
    
    print(f"   ✅ Found {len(hand_results)} HAND tiles")
    
    # Handle multiple tiles with mosaicing
    hand_uris = hand_results['URI'].tolist()
    
    if len(hand_uris) > 1:
        print(f"   Mosaicking {len(hand_uris)} HAND tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in hand_uris]
        temp_mosaic = f"/tmp/HAND_mosaic_temp.tif"
        
        # Use max for HAND data (taking maximum height above drainage)
        hand_ds = mosaic_and_fn(
            input_filelist=vsi_paths,
            output_pth=temp_mosaic,
            fn="max",
            dtype="float"
        )
    else:
        hand_uri = hand_uris[0]
        vsi_path = hand_uri.replace('s3://', '/vsis3/')
        
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        hand_ds = gdal.Open(vsi_path)
    
    if hand_ds is None:
        print(f"   ❌ Failed to open HAND data")
        return None
    
    print(f"   Opened HAND: {hand_ds.RasterXSize} x {hand_ds.RasterYSize} pixels")
    
    # First clip to polygon
    clipped_hand_path = f"/tmp/HAND_clipped_temp.tif"
    clipped_hand = clip_raster_to_vector(
        gd_raster=hand_ds,
        gpd_shp=polygon_gdf,
        out_pth=clipped_hand_path,
        crop=True,
        alltouched=False,
        crs=4326
    )
    
    # Then resample to match reference DEM
    aligned_hand_path = f"/tmp/HAND_aligned.tif"
    aligned_hand = resample_raster_to_raster(
        gd_raster=clipped_hand,
        ref_raster=reference_dem,
        interp='bilinear',  # Use bilinear for continuous HAND data
        out_pth=aligned_hand_path
    )
    
    # Compute and report statistics
    if aligned_hand:
        band = aligned_hand.GetRasterBand(1)
        stats = band.ComputeStatistics(False)
        print(f"   HAND stats: Min={stats[0]:.1f}m, Max={stats[1]:.1f}m, Mean={stats[2]:.1f}m")
        
        # Check for reasonable HAND values (0-100m typical range)
        if stats[1] > 500:
            print(f"   ⚠️ Warning: Unusually high HAND values detected (max={stats[1]:.1f}m)")
    
    print(f"   ✅ HAND processed and aligned to reference DEM")
    
    return {
        'raster': aligned_hand,
        'path': aligned_hand_path,
        'name': 'HAND',
        'description': 'Height Above Nearest Drainage (m)',
        'units': 'meters'
    }



def query_and_process_soil_features(geojson_path: str, polygon_gdf: gpd.GeoDataFrame, 
                                   reference_dem: gdal.Dataset) -> Dict:
    """
    Query and process soil features using max aggregation to exclude NoData.
    SoilGrids provides values in g/kg (0-1000 range), with NoData = -32768.
    """
    print(f"\n3c. Querying and processing soil features...")
    
    geo_db = Geospatial_DB(enable_s3=True)
    soil_features_dict = {}
    
    for soil_feature in SOIL_FEATURES:
        print(f"   Processing {soil_feature}...")
        
        soil_results = geo_db.query_asset_catalog(
            geojson_pth=geojson_path,
            product=soil_feature,
            start_date='2000-01-01',
            end_date='2025-01-01'
        )
        
        if soil_results is None or len(soil_results) == 0:
            print(f"      ⚠️ No {soil_feature} data found")
            continue
        
        soil_uris = soil_results['URI'].tolist()
        
        if len(soil_uris) > 1:
            print(f"      Mosaicking {len(soil_uris)} tiles for {soil_feature}...")
            vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in soil_uris]
            temp_mosaic = f"/tmp/{soil_feature}_mosaic_temp.tif"
            
            # Use 'max' to effectively exclude NoData (-32768) values
            soil_ds = mosaic_and_fn(
                input_filelist=vsi_paths,
                output_pth=temp_mosaic,
                fn="max",  # This works because valid values (0-1000) > NoData (-32768)
                dtype="float"
            )
        else:
            soil_uri = soil_uris[0]
            vsi_path = soil_uri.replace('s3://', '/vsis3/')
            soil_ds = gdal.Open(vsi_path)
        
        if soil_ds is None:
            print(f"      ❌ Failed to open {soil_feature}")
            continue
        
        # Set NoData value explicitly
        band = soil_ds.GetRasterBand(1)
        band.SetNoDataValue(-32768)
        
        # Clip to polygon
        clipped_soil_path = f"/tmp/{soil_feature}_clipped_temp.tif"
        clipped_soil = clip_raster_to_vector(
            gd_raster=soil_ds,
            gpd_shp=polygon_gdf,
            out_pth=clipped_soil_path,
            crop=True,
            alltouched=False,
            crs=4326
        )
        
        # Resample to match reference DEM
        aligned_path = f"/tmp/{soil_feature}_aligned.tif"
        aligned_soil = resample_raster_to_raster(
            gd_raster=clipped_soil,
            ref_raster=reference_dem,
            interp='bilinear',
            out_pth=aligned_path
        )
        
        # Verify stats are in correct range (g/kg not %)
        if aligned_soil:
            band = aligned_soil.GetRasterBand(1)
            band.SetNoDataValue(-32768)
            stats = band.ComputeStatistics(False)
            print(f"      Stats: Min={stats[0]:.1f}, Max={stats[1]:.1f}, Mean={stats[2]:.1f} g/kg")
            
            # Check for CORRECT SoilGrids range (0-1000 g/kg)
            if stats[0] < -100 or stats[1] > 1000:
                print(f"      ⚠️ Warning: Values outside expected range (0-1000 g/kg)")
            else:
                # Convert to percentage for display only
                print(f"      ({stats[2]/10:.1f}% average content)")
        
        soil_features_dict[soil_feature] = {
            'raster': aligned_soil,
            'path': aligned_path,
            'name': soil_feature,
            'description': f'Soil {soil_feature.split("_")[0].lower()} content (g/kg)',
            'units': 'g/kg'
        }
        
        print(f"      ✅ {soil_feature} processed and aligned")
    
    print(f"   ✅ Processed {len(soil_features_dict)} soil features")
    return soil_features_dict


def derive_dem_features(dem_ds: gdal.Dataset) -> Dict:
    """
    Calculate multiple types of local elevation anomalies - OPTIMIZED VERSION.
    """
    print(f"\n   Deriving DEM features with multiple anomaly types...")
    
    dem_band = dem_ds.GetRasterBand(1)
    dem_array = dem_band.ReadAsArray().astype(np.float32)
    nodata = dem_band.GetNoDataValue()
    
    if nodata is not None:
        valid_mask = dem_array != nodata
        dem_working = np.where(valid_mask, dem_array, np.nan)
    else:
        valid_mask = np.ones_like(dem_array, dtype=bool)
        dem_working = dem_array.copy()
    
    derived_features = {}
    
    def out_with_nodata(arr: np.ndarray) -> np.ndarray:
        if nodata is not None:
            return np.where(valid_mask & np.isfinite(arr), arr, nodata).astype(np.float32)
        else:
            return np.where(valid_mask & np.isfinite(arr), arr, -9999).astype(np.float32)
    
    # Helper functions remain the same...
    def fast_nanmin_filter(arr: np.ndarray, size: int) -> np.ndarray:
        """Fast minimum filter that handles NaN values."""
        arr_filled = np.where(np.isnan(arr), np.inf, arr)
        result = ndimage.minimum_filter(arr_filled, size=size, mode='constant', cval=np.inf)
        return np.where(np.isinf(result), np.nan, result)
    
    def fast_nanmax_filter(arr: np.ndarray, size: int) -> np.ndarray:
        """Fast maximum filter that handles NaN values."""
        arr_filled = np.where(np.isnan(arr), -np.inf, arr)
        result = ndimage.maximum_filter(arr_filled, size=size, mode='constant', cval=-np.inf)
        return np.where(np.isinf(result), np.nan, result)
    
    def fast_nanmedian_filter(arr: np.ndarray, size: int) -> np.ndarray:
        """Optimized median filter for large windows."""
        if size <= 15:
            return ndimage.median_filter(arr, size=size, mode='constant', cval=np.nan)
        else:
            downsample_factor = max(1, size // 15)
            if downsample_factor > 1:
                small_arr = arr[::downsample_factor, ::downsample_factor]
                small_size = max(3, size // downsample_factor)
                small_median = ndimage.median_filter(small_arr, size=small_size, 
                                                    mode='constant', cval=np.nan)
                from scipy.ndimage import zoom
                zoom_factor = arr.shape[0] / small_median.shape[0], arr.shape[1] / small_median.shape[1]
                result = zoom(small_median, zoom_factor, order=1)
                result = result[:arr.shape[0], :arr.shape[1]]
                return result
            else:
                return ndimage.median_filter(arr, size=size, mode='constant', cval=np.nan)
    
    # 1. Standard elevation anomaly (elevation - local minimum)
    anomaly_windows = [3, 7, 15, 31]
    print(f"      Computing elevation_anomaly for windows: {anomaly_windows}")
    
    for window in anomaly_windows:
        print(f"         Processing {window}x{window} anomaly...")
        local_min = fast_nanmin_filter(dem_working, window)
        anomaly = dem_array - local_min
        anomaly = out_with_nodata(anomaly)
        
        feature_name = f"elevation_anomaly_{window}x{window}"
        derived_features[feature_name] = {
            'array': anomaly,
            'window_size': window,
            'description': f'Elevation above lowest pixel in {window}x{window} window',
            'type': 'anomaly'
        }
        
        vm = (valid_mask & np.isfinite(anomaly))
        if vm.any():
            print(f"            Mean: {np.nanmean(anomaly[vm]):.2f}m, Max: {np.nanmax(anomaly[vm]):.2f}m")
    
    # 2. Inverse elevation anomaly (local maximum - elevation)
    inverse_windows = [5, 9, 17, 33]
    print(f"      Computing elevation_anomaly_inverse for windows: {inverse_windows}")
    
    for window in inverse_windows:
        print(f"         Processing {window}x{window} inverse anomaly...")
        local_max = fast_nanmax_filter(dem_working, window)
        inverse_anomaly = local_max - dem_array
        inverse_anomaly = out_with_nodata(inverse_anomaly)
        
        feature_name = f"elevation_anomaly_inverse_{window}x{window}"
        derived_features[feature_name] = {
            'array': inverse_anomaly,
            'window_size': window,
            'description': f'Highest pixel in {window}x{window} window above elevation',
            'type': 'anomaly_inverse'
        }
        
        vm = (valid_mask & np.isfinite(inverse_anomaly))
        if vm.any():
            print(f"            Mean: {np.nanmean(inverse_anomaly[vm]):.2f}m, Max: {np.nanmax(inverse_anomaly[vm]):.2f}m")
    
    # 3. Relative elevation anomaly (elevation - local median)
    relative_windows = [7, 19, 35]
    print(f"      Computing elevation_anomaly_relative for windows: {relative_windows}")
    
    for window in relative_windows:
        print(f"         Processing {window}x{window} relative anomaly...")
        if window > 20:
            print(f"            Using optimized method for large window...")
        
        local_median = fast_nanmedian_filter(dem_working, window)
        relative_anomaly = dem_array - local_median
        relative_anomaly = out_with_nodata(relative_anomaly)
        
        feature_name = f"elevation_anomaly_relative_{window}x{window}"
        derived_features[feature_name] = {
            'array': relative_anomaly,
            'window_size': window,
            'description': f'Elevation relative to median in {window}x{window} window',
            'type': 'anomaly_relative'
        }
        
        vm = (valid_mask & np.isfinite(relative_anomaly))
        if vm.any():
            print(f"            Mean: {np.nanmean(relative_anomaly[vm]):.2f}m, Std: {np.nanstd(relative_anomaly[vm]):.2f}m")
    
    # 4. NEW: Ratio of upslope area
    upslope_windows = [3, 9, 21]
    print(f"      Computing ratio_upslope_area for windows: {upslope_windows}")
    
    def count_upslope_ratio(arr: np.ndarray, window_size: int) -> np.ndarray:
        """
        OPTIMIZED: Vectorized upslope ratio calculation using broadcasting.
        Computes for each pixel: (# pixels higher than center) / (# valid pixels in window)
        """
        from scipy.ndimage import uniform_filter, generic_filter
        
        result = np.full(arr.shape, np.nan, dtype=np.float32)
        half_window = window_size // 2
        
        print(f"         Computing upslope ratios (vectorized method)...")
        
        # Create a mask of valid pixels
        valid_mask = np.isfinite(arr)
        
        # For small windows, use generic_filter (still faster than nested loops)
        if window_size <= 9:
            def upslope_func(values):
                # Center value is the middle element
                center_idx = len(values) // 2
                center = values[center_idx]
                if not np.isfinite(center):
                    return np.nan
                
                valid_vals = values[np.isfinite(values)]
                if len(valid_vals) == 0:
                    return np.nan
                
                upslope_count = np.sum(valid_vals > center)
                return upslope_count / len(valid_vals)
            
            result = generic_filter(arr, upslope_func, size=window_size, 
                                   mode='constant', cval=np.nan)
        else:
            # For large windows, use strided approach with chunking
            print(f"            Using chunked processing for large window...")
            padded = np.pad(arr, half_window, mode='constant', constant_values=np.nan)
            
            # Process in chunks for memory efficiency
            chunk_size = 50
            n_rows = arr.shape[0]
            n_chunks = (n_rows + chunk_size - 1) // chunk_size
            
            for chunk_idx in range(n_chunks):
                start_row = chunk_idx * chunk_size
                end_row = min(start_row + chunk_size, n_rows)
                
                if chunk_idx % max(1, n_chunks // 10) == 0:
                    progress = 100 * chunk_idx // n_chunks
                    print(f"            Progress: {progress}%")
                
                # Process this chunk of rows
                for i in range(start_row, end_row):
                    center_vals = arr[i, :]
                    valid_centers = np.isfinite(center_vals)
                    
                    if not valid_centers.any():
                        continue
                    
                    # Vectorized: extract all windows for this row at once
                    for j in np.where(valid_centers)[0]:
                        window = padded[i:i+window_size, j:j+window_size]
                        valid_window = window[np.isfinite(window)]
                        
                        if len(valid_window) == 0:
                            continue
                        
                        upslope_count = np.sum(valid_window > center_vals[j])
                        result[i, j] = upslope_count / len(valid_window)
        
        print(f"            Progress: 100%")
        return result.astype(np.float32)
    
    for window in upslope_windows:
        print(f"         Processing {window}x{window} upslope ratio...")
        
        upslope_ratio = count_upslope_ratio(dem_working, window)
        upslope_ratio = out_with_nodata(upslope_ratio)
        
        feature_name = f"ratio_upslope_area_{window}x{window}"
        derived_features[feature_name] = {
            'array': upslope_ratio,
            'window_size': window,
            'description': f'Ratio of upslope pixels in {window}x{window} window',
            'type': 'upslope_ratio'
        }
        
        vm = (valid_mask & np.isfinite(upslope_ratio))
        if vm.any():
            mean_ratio = np.nanmean(upslope_ratio[vm])
            print(f"            Mean ratio: {mean_ratio:.3f} ({mean_ratio*100:.1f}%)")
    
    # Original elevation
    derived_features['elevation'] = {
        'array': out_with_nodata(dem_array),
        'window_size': 1,
        'description': 'Original elevation',
        'type': 'base'
    }
    
    derived_features['metadata'] = {
        'geotransform': dem_ds.GetGeoTransform(),
        'projection': dem_ds.GetProjection(),
        'width': dem_ds.RasterXSize,
        'height': dem_ds.RasterYSize,
        'nodata': nodata if nodata is not None else -9999
    }
    
    print(f"      ✅ Generated {len(derived_features)-1} DEM-derived features")
    return derived_features

def compute_built_area_ratio(arr: np.ndarray, window_size: int, nodata_val) -> np.ndarray:
    """
    Compute ratio of Built Area (class 7) pixels in sliding window.
    Returns values from 0.0 to 1.0 representing the fraction of built area.
    """
    from scipy import ndimage
    
    result = np.full(arr.shape, nodata_val if nodata_val is not None else -9999, dtype=np.float32)
    
    print(f"         Computing built area ratios for {arr.shape[0]}x{arr.shape[1]} pixels...")
    
    # For small windows, use optimized generic_filter
    if window_size <= 9:
        def built_ratio_func(values):
            # Filter nodata
            if nodata_val is not None:
                valid = values[values != nodata_val]
            else:
                valid = values[values >= 0]
            
            if len(valid) == 0:
                return nodata_val if nodata_val is not None else -9999
            
            # Count Built Area pixels (class 7)
            built_count = np.sum(valid == 7)
            return float(built_count) / len(valid)
        
        result = ndimage.generic_filter(
            arr,
            built_ratio_func,
            size=window_size,
            mode='constant',
            cval=nodata_val if nodata_val is not None else -9999
        )
    else:
        # For large windows, use batch processing
        print(f"            Using batch processing for large window...")
        
        half_window = window_size // 2
        pad_value = nodata_val if nodata_val is not None else -9999
        padded = np.pad(arr, half_window, mode='constant', constant_values=pad_value)
        
        chunk_size = 100
        n_rows = arr.shape[0]
        n_chunks = (n_rows + chunk_size - 1) // chunk_size
        
        for chunk_idx in range(n_chunks):
            start_row = chunk_idx * chunk_size
            end_row = min(start_row + chunk_size, n_rows)
            
            if chunk_idx % max(1, n_chunks // 10) == 0:
                progress = 100 * chunk_idx // n_chunks
                print(f"            Progress: {progress}%")
            
            # Process chunk rows
            for i in range(start_row, end_row):
                for j in range(arr.shape[1]):
                    # Extract window
                    window = padded[i:i+window_size, j:j+window_size]
                    
                    # Filter nodata
                    if nodata_val is not None:
                        valid_values = window[window != nodata_val]
                    else:
                        valid_values = window[window >= 0]
                    
                    if len(valid_values) == 0:
                        continue
                    
                    # Calculate ratio of Built Area (class 7)
                    built_count = np.sum(valid_values == 7)
                    result[i, j] = float(built_count) / len(valid_values)
        
        print(f"            Progress: 100%")
    
    return result.astype(np.float32)


def derive_lulc_features(lulc_ds: gdal.Dataset) -> Dict:
    """
    Calculate neighborhood features from LULC data.
    Computes the most common (modal) LULC class in various window sizes.
    """
    print(f"\n   Deriving LULC neighborhood features...")
    
    lulc_band = lulc_ds.GetRasterBand(1)
    lulc_array = lulc_band.ReadAsArray().astype(np.float32)
    nodata = lulc_band.GetNoDataValue()
    
    if nodata is not None:
        valid_mask = lulc_array != nodata
    else:
        valid_mask = np.ones_like(lulc_array, dtype=bool)
    
    derived_features = {}
    
    def compute_modal_filter(arr: np.ndarray, window_size: int, nodata_val) -> np.ndarray:
        """
        OPTIMIZED: Compute modal (most common) value using vectorized operations.
        ~100x faster than nested loops for large arrays.
        """
        from scipy import ndimage
        from scipy.stats import mode as scipy_mode
        
        result = np.full(arr.shape, nodata_val if nodata_val is not None else -9999, dtype=np.float32)
        
        print(f"         Computing modal values for {arr.shape[0]}x{arr.shape[1]} pixels...")
        
        # For small windows (<=9), use optimized generic_filter
        if window_size <= 9:
            # Define mode function for generic_filter
            def mode_func(values):
                # Filter nodata
                if nodata_val is not None:
                    valid = values[values != nodata_val]
                else:
                    valid = values[values >= 0]
                
                if len(valid) == 0:
                    return nodata_val if nodata_val is not None else -9999
                
                # Fast mode using bincount
                valid_int = valid.astype(np.int32)
                # Handle potential negative or very large values
                if valid_int.min() < 0 or valid_int.max() > 100:
                    # Fallback for unusual ranges
                    from collections import Counter
                    counts = Counter(valid_int)
                    return float(counts.most_common(1)[0][0])
                
                counts = np.bincount(valid_int)
                return float(np.argmax(counts))
            
            # Use generic_filter with mode function
            result = ndimage.generic_filter(
                arr, 
                mode_func, 
                size=window_size, 
                mode='constant',
                cval=nodata_val if nodata_val is not None else -9999
            )
            
        else:
            # For large windows, use strided approach with batch processing
            print(f"            Using optimized batch processing for large window...")
            
            half_window = window_size // 2
            pad_value = nodata_val if nodata_val is not None else -9999
            padded = np.pad(arr, half_window, mode='constant', constant_values=pad_value)
            
            # Process in chunks for better cache efficiency
            chunk_size = 100
            n_rows = arr.shape[0]
            n_chunks = (n_rows + chunk_size - 1) // chunk_size
            
            for chunk_idx in range(n_chunks):
                start_row = chunk_idx * chunk_size
                end_row = min(start_row + chunk_size, n_rows)
                
                if chunk_idx % max(1, n_chunks // 10) == 0:
                    progress = 100 * chunk_idx // n_chunks
                    print(f"            Progress: {progress}%")
                
                # Process chunk rows
                for i in range(start_row, end_row):
                    for j in range(arr.shape[1]):
                        # Extract window
                        window = padded[i:i+window_size, j:j+window_size]
                        
                        # Filter nodata
                        if nodata_val is not None:
                            valid_values = window[window != nodata_val]
                        else:
                            valid_values = window[window >= 0]
                        
                        if len(valid_values) == 0:
                            continue
                        
                        # Compute mode
                        valid_int = valid_values.astype(np.int32)
                        
                        # Optimize bincount for typical LULC range (0-255)
                        if valid_int.min() >= 0 and valid_int.max() <= 255:
                            counts = np.bincount(valid_int, minlength=256)
                            mode_value = np.argmax(counts[:256])
                        else:
                            # Fallback for unusual ranges
                            counts = np.bincount(valid_int)
                            mode_value = np.argmax(counts)
                        
                        result[i, j] = float(mode_value)
        
        print(f"            Progress: 100%")
        return result.astype(np.float32)    
        
    # Modal LULC in various window sizes
    lulc_windows = [3, 9, 27, 81]
    print(f"      Computing most_common_lulc_neighbor for windows: {lulc_windows}")
    
    for window in lulc_windows:
        print(f"         Processing {window}x{window} modal LULC...")
        
        modal_lulc = compute_modal_filter(lulc_array, window, nodata)
        
        feature_name = f"most_common_lulc_neighbor_{window}x{window}"
        derived_features[feature_name] = {
            'array': modal_lulc,
            'window_size': window,
            'description': f'Most common LULC class in {window}x{window} window',
            'type': 'lulc_modal'
        }
        
        # Statistics
        if nodata is not None:
            valid_data = modal_lulc[modal_lulc != nodata]
        else:
            valid_data = modal_lulc[modal_lulc >= 0]
        
        if len(valid_data) > 0:
            unique_classes = np.unique(valid_data)
            print(f"            Found {len(unique_classes)} unique LULC classes in output")

    # Built Area ratio in various window sizes
    built_area_windows = [3, 9, 27]
    print(f"      Computing ratio_built_area for windows: {built_area_windows}")
    
    for window in built_area_windows:
        print(f"         Processing {window}x{window} built area ratio...")
        
        built_ratio = compute_built_area_ratio(lulc_array, window, nodata)
        
        feature_name = f"ratio_built_area_{window}x{window}"
        derived_features[feature_name] = {
            'array': built_ratio,
            'window_size': window,
            'description': f'Ratio of Built Area (class 7) in {window}x{window} window',
            'type': 'lulc_built_ratio'
        }
        
        # Statistics
        valid_data = built_ratio[(built_ratio != nodata) & np.isfinite(built_ratio)] if nodata is not None else built_ratio[np.isfinite(built_ratio)]
        
        if len(valid_data) > 0:
            mean_ratio = np.mean(valid_data)
            max_ratio = np.max(valid_data)
            print(f"            Mean ratio: {mean_ratio:.3f} ({mean_ratio*100:.1f}%), Max: {max_ratio:.3f}")
    
    print(f"      ✅ Generated {len(derived_features)} LULC features (modal + built area ratio)")
    return derived_features
    
def create_multiband_raster(
    all_features: Dict, 
    output_path: str, 
    metadata: Dict,
    use_bigtiff: bool = False  # NEW parameter
) -> Tuple[gdal.Dataset, List[str]]:
    """
    Stack all features into a multi-band GeoTIFF with proper NoData handling per band.
    Includes all three types of elevation anomalies plus other features.
    
    Parameters:
    -----------
    all_features : Dict
        Dictionary of all features to include
    output_path : str
        Path for output file
    metadata : Dict
        Raster metadata (geotransform, projection, dimensions, nodata)
    use_bigtiff : bool
        If True, create BigTIFF format (for files > 4GB)
    """
    print(f"\n   Creating multi-band raster with enhanced DEM features...")
    if use_bigtiff:
        print(f"   Using BigTIFF format for large file support")
    
    # Reference dims
    ref_height = metadata['height']
    ref_width = metadata['width']
    ref_shape = (ref_height, ref_width)
    
    band_order = []
    band_arrays = []
    band_nodata_values = []
    
    default_nodata = metadata.get('nodata', -9999)
     
    # Elevation first
    if 'elevation' in all_features:
        band_order.append('elevation')
        array = all_features['elevation']['array']
        if array.shape != ref_shape:
            print(f"      Aligning elevation: {array.shape} -> {ref_shape}")
            array = align_array_to_reference(array, ref_shape)
        band_arrays.append(array.astype(np.float32))
        band_nodata_values.append(default_nodata)
    
    # Standard anomalies
    for window in [3, 7, 15, 31]:
        key = f'elevation_anomaly_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            band_nodata_values.append(default_nodata)
    
    # Inverse anomalies
    for window in [5, 9, 17, 33]:
        key = f'elevation_anomaly_inverse_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            band_nodata_values.append(default_nodata)
    
    # Relative anomalies
    for window in [7, 19, 35]:
        key = f'elevation_anomaly_relative_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            band_nodata_values.append(default_nodata)
    
    # Upslope ratio features
    for window in [3, 9, 21]:
        key = f'ratio_upslope_area_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            band_nodata_values.append(default_nodata)
    
    # HAND
    if 'HAND' in all_features:
        band_order.append('HAND')
        hand_ds = all_features['HAND']['raster']
        band = hand_ds.GetRasterBand(1)
        array = band.ReadAsArray().astype(np.float32)
        nodata = band.GetNoDataValue()
        if array.shape != ref_shape:
            print(f"      Aligning HAND: {array.shape} -> {ref_shape}")
            array = align_array_to_reference(array, ref_shape)
        band_arrays.append(array)
        band_nodata_values.append(nodata if nodata is not None else default_nodata)
    
    # LULC
    if 'IO_LULC' in all_features:
        band_order.append('IO_LULC')
        lulc_ds = all_features['IO_LULC']['raster']
        band = lulc_ds.GetRasterBand(1)
        array = band.ReadAsArray().astype(np.float32)
        nodata = band.GetNoDataValue()
        if array.shape != ref_shape:
            print(f"      Aligning IO_LULC: {array.shape} -> {ref_shape}")
            array = align_array_to_reference(array, ref_shape)
        band_arrays.append(array)
        band_nodata_values.append(nodata if nodata is not None else default_nodata)
 
    # LULC neighborhood features
    for window in [3, 9, 27, 81]:
        key = f'most_common_lulc_neighbor_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            nodata = all_features['IO_LULC']['raster'].GetRasterBand(1).GetNoDataValue()
            band_nodata_values.append(nodata if nodata is not None else default_nodata)

    # LULC Built Area ratios
    for window in [3, 9, 27]:
        key = f'ratio_built_area_{window}x{window}'
        if key in all_features:
            band_order.append(key)
            array = all_features[key]['array']
            if array.shape != ref_shape:
                print(f"      Aligning {key}: {array.shape} -> {ref_shape}")
                array = align_array_to_reference(array, ref_shape)
            band_arrays.append(array.astype(np.float32))
            nodata = all_features['IO_LULC']['raster'].GetRasterBand(1).GetNoDataValue()
            band_nodata_values.append(nodata if nodata is not None else default_nodata)

    
    # Soil features
    for soil_feature in SOIL_FEATURES:
        if soil_feature in all_features:
            band_order.append(soil_feature)
            soil_ds = all_features[soil_feature]['raster']
            band = soil_ds.GetRasterBand(1)
            array = band.ReadAsArray().astype(np.float32)
            nodata = -32768.0
            if array.shape != ref_shape:
                print(f"      Aligning {soil_feature}: {array.shape} -> {ref_shape}")
                aligned = np.full(ref_shape, nodata, dtype=np.float32)
                min_h = min(array.shape[0], ref_shape[0])
                min_w = min(array.shape[1], ref_shape[1])
                aligned[:min_h, :min_w] = array[:min_h, :min_w]
                array = aligned
            band_arrays.append(array)
            band_nodata_values.append(nodata)
    
    n_bands = len(band_order)
    print(f"      Stacking {n_bands} bands with dimensions {ref_shape}")
    
    # Create TIFF with BigTIFF support if requested
    driver = gdal.GetDriverByName('GTiff')
    
    # Build creation options
    creation_options = ['COMPRESS=DEFLATE', 'TILED=YES']
    if use_bigtiff:
        creation_options.append('BIGTIFF=YES')
    
    out_ds = driver.Create(
        output_path,
        ref_width,
        ref_height,
        n_bands,
        gdal.GDT_Float32,
        options=creation_options
    )
    out_ds.SetGeoTransform(metadata['geotransform'])
    out_ds.SetProjection(metadata['projection'])
    
    # Write each band
    for band_idx, (band_name, array, nodata_val) in enumerate(zip(band_order, band_arrays, band_nodata_values), 1):
        out_band = out_ds.GetRasterBand(band_idx)
        arr = array.astype(np.float32, copy=True)
        if nodata_val is not None:
            nan_mask = np.isnan(arr)
            if nan_mask.any():
                arr[nan_mask] = float(nodata_val)
            out_band.SetNoDataValue(float(nodata_val))
        out_band.WriteArray(arr)
        out_band.SetDescription(band_name)
        out_band.FlushCache()
    
    print(f"   ✅ Created multi-band raster with {n_bands} bands")
    
    # Optional band summary (unchanged)
    dem_bands = [b for b in band_order if 'elevation' in b.lower() or 'anomaly' in b.lower()]
    other_bands = [b for b in band_order if b not in dem_bands]
    
    print(f"      DEM-derived bands ({len(dem_bands)}):")
    for i, band in enumerate(dem_bands[:5], 1):
        print(f"        Band {i}: {band}")
    if len(dem_bands) > 5:
        print(f"        ... and {len(dem_bands)-5} more")
    
    if other_bands:
        print(f"      Other bands ({len(other_bands)}):")
        for band in other_bands[:3]:
            idx = band_order.index(band) + 1
            print(f"        Band {idx}: {band}")
        if len(other_bands) > 3:
            print(f"        ... and {len(other_bands)-3} more")
    
    return out_ds, band_order

def _finalize_output(
    multiband_ds: gdal.Dataset,
    multiband_path: str,
    band_names: List[str],
    base_filename: str,
    output_dir: str,
    geojson_path: str,
    all_features: Dict,
    dem_results
) -> Dict:
    """
    Finalize outputs: plot, save, and return results.
    """
    # Plot all bands
    print(f"\nPlotting all bands...")
    plot_path = f"/tmp/multiband_{base_filename}_plot.png"
    figure = plot_all_bands(multiband_ds, band_names, save_path=plot_path)
    plt.show()

    # CRITICAL: Close the dataset to flush all data to disk
    multiband_ds.FlushCache()
    multiband_ds = None
    print("   ✅ Dataset closed and data written to disk")

    # Save to final location
    print(f"\nSaving multi-band raster...")
    
    if output_dir.startswith('s3://'):
        s3_client = boto3.client('s3')
        s3_parts = output_dir.replace('s3://', '').rstrip('/').split('/', 1)
        bucket = s3_parts[0]
        key_prefix = s3_parts[1] if len(s3_parts) > 1 else ''
        
        multiband_filename = f"multiband_{base_filename}.tif"
        multiband_key = f"{key_prefix}/{multiband_filename}" if key_prefix else multiband_filename
        s3_client.upload_file(multiband_path, bucket, multiband_key)
        final_multiband_path = f"s3://{bucket}/{multiband_key}"
        
        plot_filename = f"multiband_{base_filename}_plot.png"
        plot_key = f"{key_prefix}/{plot_filename}" if key_prefix else plot_filename
        s3_client.upload_file(plot_path, bucket, plot_key)
        final_plot_path = f"s3://{bucket}/{plot_key}"
        
        print(f"   ✅ Uploaded to S3")
    else:
        import shutil
        final_multiband_path = os.path.join(output_dir, f"multiband_{base_filename}.tif")
        final_plot_path = os.path.join(output_dir, f"multiband_{base_filename}_plot.png")
        shutil.move(multiband_path, final_multiband_path)
        shutil.move(plot_path, final_plot_path)
        print(f"   ✅ Saved locally")
    
    # Cleanup
    if os.path.exists(geojson_path):
        os.remove(geojson_path)
    
    # Count feature types
    n_anomaly = len([k for k in all_features if 'anomaly' in k and 'inverse' not in k and 'relative' not in k])
    n_inverse = len([k for k in all_features if 'inverse' in k])
    n_relative = len([k for k in all_features if 'relative' in k])
    
    return {
        'multiband_path': final_multiband_path,
        'plot_path': final_plot_path,
        'band_names': band_names,
        'n_bands': len(band_names),
        'base_filename': base_filename,
        'dem_results': dem_results,
        'has_hand': 'HAND' in all_features,
        'has_lulc': 'IO_LULC' in all_features,
        'n_soil_features': len([k for k in all_features.keys() if k in SOIL_FEATURES]),
        'n_elevation_anomaly': n_anomaly,
        'n_elevation_inverse': n_inverse,
        'n_elevation_relative': n_relative
    }

def plot_all_bands(multiband_raster: gdal.Dataset, band_names: List[str], 
                   save_path: Optional[str] = None, figsize: Tuple = (24, 20)) -> plt.Figure:
    """
    Plot all bands from the multi-band raster including the three types of elevation anomalies.
    """
    n_bands = len(band_names)
    n_cols = 5
    n_rows = int(np.ceil(n_bands / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes
    
    for band_idx, band_name in enumerate(band_names):
        ax = axes[band_idx]
        band = multiband_raster.GetRasterBand(band_idx + 1)
        array = band.ReadAsArray()
        nodata = band.GetNoDataValue()
        
        if nodata is not None:
            array = np.ma.masked_equal(array, nodata)
            if any(soil_type in band_name for soil_type in ['CLAY', 'SAND', 'SILT']):
                array = np.ma.masked_where((array < 0) | (array > 1000), array)
        
        # Colormap selection
        if band_name == 'HAND':
            cmap = 'Blues_r'
        elif 'elevation' in band_name.lower() and 'anomaly' not in band_name.lower():
            cmap = 'terrain'
        elif 'anomaly_inverse' in band_name.lower():
            cmap = 'Reds'
        elif 'anomaly_relative' in band_name.lower():
            cmap = 'RdBu_r'
        elif 'anomaly' in band_name.lower():
            cmap = 'YlOrRd'
        elif 'ratio_upslope_area' in band_name.lower():  # NEW
            cmap = 'viridis'
        elif 'most_common_lulc' in band_name.lower():  # NEW
            cmap = 'tab20'
        elif 'ratio_built_area' in band_name.lower():
            cmap = 'Reds'  # Red intensity for urbanization
        elif 'lulc' in band_name.lower():
            cmap = 'tab20'
        elif 'clay' in band_name.lower():
            cmap = 'YlOrBr'
        elif 'sand' in band_name.lower():
            cmap = 'YlGn'
        elif 'silt' in band_name.lower():
            cmap = 'PuBu'
        else:
            cmap = 'viridis'
        
        # Center relative anomaly around 0
        if 'anomaly_relative' in band_name.lower() and hasattr(array, 'compressed'):
            valid_data = array.compressed()
            if len(valid_data) > 0:
                vmax = np.percentile(np.abs(valid_data), 95)
                im = ax.imshow(array, cmap=cmap, aspect='auto', vmin=-vmax, vmax=vmax)
            else:
                im = ax.imshow(array, cmap=cmap, aspect='auto')
        else:
            im = ax.imshow(array, cmap=cmap, aspect='auto')
        
        # Titles
        if band_name == 'HAND':
            title = 'Height Above\nNearest Drainage'
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.1f}m mean)"
        elif 'elevation_anomaly_inverse' in band_name:
            window = band_name.split('_')[-1]
            title = f"Inverse Anomaly {window}\n(max - elevation)"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.1f}m mean)"
        elif 'elevation_anomaly_relative' in band_name:
            window = band_name.split('_')[-1]
            title = f"Relative Anomaly {window}\n(elevation - median)"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    std_val = valid_data.std()
                    title += f"\n(μ={mean_val:.1f}m, σ={std_val:.1f}m)"
        elif 'elevation_anomaly' in band_name:
            window = band_name.split('_')[-1]
            title = f"Elevation Anomaly {window}\n(elevation - min)"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.1f}m mean)"
        elif band_name in SOIL_FEATURES:
            # Improve soil label formatting: e.g., CLAY_0_TO_5CM -> "Clay 0 to 5 cm"
            parts = band_name.split('_')  # e.g., ['CLAY','0','TO','5CM']
            soil_type = parts[0].title()
            depth = ' '.join(parts[1:]).lower().replace('_', ' ')  # "0 to 5cm"
            depth = depth.replace('cm', ' cm')  # "0 to 5 cm"
            title = f"{soil_type} {depth}"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.0f} g/kg mean)"
        # Titles - add new cases
        elif 'ratio_upslope_area' in band_name:
            window = band_name.split('_')[-1]
            title = f"Upslope Ratio {window}"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.2%} mean)"
        elif 'most_common_lulc' in band_name:
            window = band_name.split('_')[-1]
            title = f"Modal LULC {window}"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    n_classes = len(np.unique(valid_data))
                    title += f"\n({n_classes} classes)"
        elif 'ratio_built_area' in band_name:
            window = band_name.split('_')[-1]
            title = f"Built Area Ratio {window}"
            if hasattr(array, 'compressed'):
                valid_data = array.compressed()
                if len(valid_data) > 0:
                    mean_val = valid_data.mean()
                    title += f"\n({mean_val:.2%} mean)"
                    
        else:
            title = band_name.replace('_', ' ').title()
        
        ax.set_title(title, fontsize=9)
        ax.axis('off')
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046)
        cbar.ax.tick_params(labelsize=7)
        if band_name == 'HAND':
            cbar.set_label('Height (m)', rotation=270, labelpad=12, fontsize=8)
        elif 'anomaly_inverse' in band_name:
            cbar.set_label('Height diff (m)', rotation=270, labelpad=12, fontsize=8)
        elif 'anomaly_relative' in band_name:
            cbar.set_label('Rel. height (m)', rotation=270, labelpad=12, fontsize=8)
        elif 'anomaly' in band_name and 'inverse' not in band_name and 'relative' not in band_name:
            cbar.set_label('Height above min (m)', rotation=270, labelpad=12, fontsize=8)
        elif any(soil_type in band_name for soil_type in ['CLAY', 'SAND', 'SILT']):
            cbar.set_label('g/kg', rotation=270, labelpad=12, fontsize=8)
        elif 'ratio_upslope_area' in band_name:
            cbar.set_label('Ratio (0-1)', rotation=270, labelpad=12, fontsize=8)
        elif 'most_common_lulc' in band_name:
            cbar.set_label('LULC class', rotation=270, labelpad=12, fontsize=8)            
        elif 'ratio_built_area' in band_name:
            cbar.set_label('Built Ratio (0-1)', rotation=270, labelpad=12, fontsize=8)
    
    # Hide extras
    for idx in range(n_bands, len(axes)):
        axes[idx].set_visible(False)
    
    n_anomaly = len([b for b in band_names if 'anomaly' in b and 'inverse' not in b and 'relative' not in b])
    n_inverse = len([b for b in band_names if 'inverse' in b])
    n_relative = len([b for b in band_names if 'relative' in b])
    title_text = f'Multi-Band Raster Features ({n_bands} bands total)\n'
    title_text += f'DEM Features: {n_anomaly} anomaly, {n_inverse} inverse, {n_relative} relative'
    plt.suptitle(title_text, fontsize=14, y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"   ✅ Plot saved to: {save_path}")
    
    return fig

# Update the main workflow function to use the new features
def extract_raster_data_for_polygon(
    polygon_dict: Dict, 
    output_dir: str = "./", 
    product: str = 'NASADEM',
    reference_raster_type: str = 'DEM'  # NEW: 'DEM' or 'LULC'
) -> Dict:
    """
    Complete workflow to extract, process, clip, plot and save all raster data for a polygon.
    Now includes three types of elevation anomalies.
    
    Parameters:
    -----------
    polygon_dict : Dict
        Dictionary containing watershed polygon data
    output_dir : str
        Output directory path (local or S3)
    product : str
        DEM product to use (default: 'NASADEM')
    reference_raster_type : str
        Which raster to use as reference for reprojection/resampling.
        Options: 'DEM' (30m resolution) or 'LULC' (10m resolution)
        Default: 'DEM'
    """
    # Validate reference_raster_type
    if reference_raster_type not in ['DEM', 'LULC']:
        raise ValueError(f"reference_raster_type must be 'DEM' or 'LULC', got: {reference_raster_type}")
    
    print(f"\n{'='*60}")
    print(f"Reference Raster: {reference_raster_type}")
    if reference_raster_type == 'LULC':
        print("  Resolution: 10m (High resolution - outputs will be saved as BigTIFF)")
    else:
        print("  Resolution: 30m")
    print(f"{'='*60}\n")
    
    polygon_gdf = polygon_dict['polygon_gdf']
    geojson_path = polygon_dict['geojson_path']
    metadata = polygon_dict['metadata']
    
    all_features = {}
    
    # Branch workflow based on reference raster choice
    if reference_raster_type == 'LULC':
        # LULC-first workflow
        result = _extract_with_lulc_reference(
            polygon_gdf, geojson_path, metadata, 
            all_features, product, output_dir
        )
    else:
        # DEM-first workflow (original)
        result = _extract_with_dem_reference(
            polygon_gdf, geojson_path, metadata, 
            all_features, product, output_dir
        )
    
    return result

def query_and_clip_lulc_native(
    geojson_path: str, 
    polygon_gdf: gpd.GeoDataFrame,
    max_retries: int = 2
) -> Optional[Dict]:
    """
    Query and clip LULC at native 10m resolution.
    Validates ALL source tiles (both same-proj and cross-proj) before processing.
    """
    print(f"\n3b. Querying and clipping IO_LULC at native resolution...")
    
    geo_db = Geospatial_DB(enable_s3=True)
    
    lulc_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product='IO_LULC',
        start_date='2015-01-01',
        end_date='2025-01-01'
    )
    
    if lulc_results is None or len(lulc_results) == 0:
        print(f"   ⚠️ No IO_LULC data found for this region")
        return None
    
    print(f"   ✅ Found {len(lulc_results)} IO_LULC tiles")
    
    lulc_uris = lulc_results['URI'].tolist()
    
    if len(lulc_uris) > 1:
        print(f"   Mosaicking {len(lulc_uris)} LULC tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in lulc_uris]
        
        # Read projection of first VALID tile (not just first in list)
        reference_proj = None
        first_valid_idx = -1
        
        for idx, vsi_path in enumerate(vsi_paths):
            try:
                test_ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
                if test_ds is not None:
                    reference_proj = test_ds.GetProjection()
                    first_valid_idx = idx
                    test_ds = None
                    break
            except:
                continue
        
        if reference_proj is None:
            print(f"   ❌ Cannot open any LULC tiles")
            return None
        
        print(f"   Using tile {first_valid_idx+1} as reference projection")
        
        # Process ALL tiles (validate originals, reproject if needed)
        unified_paths = []
        temp_files = []
        
        for idx, vsi_path in enumerate(vsi_paths):
            # CRITICAL: Validate the original tile FIRST
            try:
                print(f"      Validating tile {idx+1}...")
                tile_ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
                if tile_ds is None:
                    print(f"         ⚠️ Cannot open tile {idx+1}, skipping")
                    continue
                
                # Test read a small sample
                try:
                    test_band = tile_ds.GetRasterBand(1)
                    test_sample = test_band.ReadAsArray(0, 0, 
                        min(100, tile_ds.RasterXSize),
                        min(100, tile_ds.RasterYSize))
                    
                    if test_sample is None:
                        print(f"         ⚠️ Cannot read tile {idx+1}, skipping")
                        tile_ds = None
                        continue
                except Exception as e:
                    print(f"         ⚠️ Read error on tile {idx+1}: {str(e)[:50]}, skipping")
                    tile_ds = None
                    continue
                
                # Get projection
                tile_proj = tile_ds.GetProjection()
                tile_ds = None
                
                # If projection matches, use original path
                if tile_proj == reference_proj:
                    print(f"         ✓ Tile {idx+1} validated (same projection)")
                    unified_paths.append(vsi_path)
                else:
                    # Reproject to match reference projection
                    print(f"         Reprojecting tile {idx+1} to match reference...")
                    temp_reproj = f"/tmp/lulc_reproj_{idx}_{os.getpid()}.tif"
                    
                    warp_options = gdal.WarpOptions(
                        dstSRS=reference_proj,
                        resampleAlg='near',
                        format='GTiff',
                        creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']
                    )
                    
                    result_ds = gdal.Warp(temp_reproj, vsi_path, options=warp_options)
                    
                    if result_ds is not None:
                        result_ds.FlushCache()
                        result_ds = None
                        
                        # Validate reprojected file
                        try:
                            validate_ds = gdal.Open(temp_reproj, gdal.GA_ReadOnly)
                            if validate_ds is None:
                                print(f"            ⚠️ Tile {idx+1}: cannot reopen reprojected file, skipping")
                                try:
                                    os.remove(temp_reproj)
                                except:
                                    pass
                                continue
                            
                            test_band = validate_ds.GetRasterBand(1)
                            test_read = test_band.ReadAsArray(0, 0,
                                min(100, validate_ds.RasterXSize),
                                min(100, validate_ds.RasterYSize))
                            
                            validate_ds = None
                            
                            if test_read is None:
                                print(f"            ⚠️ Tile {idx+1}: cannot read reprojected file, skipping")
                                try:
                                    os.remove(temp_reproj)
                                except:
                                    pass
                                continue
                            
                            # Success
                            unified_paths.append(temp_reproj)
                            temp_files.append(temp_reproj)
                            print(f"            ✓ Tile {idx+1} reprojected successfully")
                            
                        except Exception as e:
                            print(f"            ⚠️ Tile {idx+1}: validation failed ({str(e)[:50]}), skipping")
                            try:
                                os.remove(temp_reproj)
                            except:
                                pass
                            continue
                    else:
                        print(f"         ⚠️ Failed to reproject tile {idx+1}, skipping")
            
            except Exception as e:
                print(f"      ⚠️ Error processing tile {idx+1}: {str(e)[:50]}")
                continue
        
        if len(unified_paths) == 0:
            print(f"   ❌ No tiles could be processed")
            return None
        
        print(f"   Using {len(unified_paths)} tiles ({len(temp_files)} reprojected to common projection)")
        
        # Now mosaic with unified projection
        temp_mosaic = f"/tmp/IO_LULC_mosaic_native_{os.getpid()}.tif"
        
        lulc_ds = None
        for attempt in range(max_retries):
            try:
                lulc_ds = mosaic_and_fn(
                    input_filelist=unified_paths,
                    output_pth=temp_mosaic,
                    fn="max", 
                    dtype="float"
                )
                
                if lulc_ds is not None:
                    band = lulc_ds.GetRasterBand(1)
                    if band is not None:
                        print(f"   ✅ Mosaic successful (attempt {attempt + 1}/{max_retries})")
                        break
                    else:
                        print(f"   ⚠️ Mosaic dataset invalid (attempt {attempt + 1}/{max_retries})")
                        lulc_ds = None
                else:
                    print(f"   ⚠️ Mosaic returned None (attempt {attempt + 1}/{max_retries})")
                
            except Exception as e:
                print(f"   ⚠️ Mosaic attempt {attempt + 1}/{max_retries} failed: {e}")
                lulc_ds = None
            
            if lulc_ds is None and attempt < max_retries - 1:
                print(f"   Retrying in 2 seconds...")
                import time
                time.sleep(2)
                try:
                    if os.path.exists(temp_mosaic):
                        os.remove(temp_mosaic)
                except:
                    pass
        
        # Clean up reprojected temp files
        for temp_file in temp_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            except:
                pass
        
        if lulc_ds is None:
            print(f"   ❌ Failed to mosaic LULC after {max_retries} attempts")
            return None
            
    else:
        # Single tile
        lulc_uri = lulc_uris[0]
        vsi_path = lulc_uri.replace('s3://', '/vsis3/')
        
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        lulc_ds = gdal.Open(vsi_path)
    
    if lulc_ds is None:
        print(f"   ❌ Failed to open IO_LULC")
        return None
    
    print(f"   Opened LULC: {lulc_ds.RasterXSize} x {lulc_ds.RasterYSize} pixels")
    
    # Clip to polygon
    clipped_lulc_path = f"/tmp/IO_LULC_native_clipped_{os.getpid()}.tif"
    
    try:
        clipped_lulc = clip_raster_to_vector(
            gd_raster=lulc_ds,
            gpd_shp=polygon_gdf,
            out_pth=clipped_lulc_path,
            crop=True,
            alltouched=False,
            crs=4326
        )
    except Exception as e:
        print(f"   ❌ Failed to clip LULC: {e}")
        lulc_ds = None
        return None
    
    lulc_ds = None
    
    if clipped_lulc is None:
        print(f"   ❌ Failed to clip LULC")
        return None
    
    # Validate clipped output
    band = clipped_lulc.GetRasterBand(1)
    array = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    
    if nodata is not None:
        valid_data = array[array != nodata]
    else:
        valid_data = array[array >= 0]
    
    if len(valid_data) > 0:
        unique_classes = np.unique(valid_data)
        print(f"   ✅ LULC clipped: found {len(unique_classes)} unique classes: {unique_classes[:10]}")
    else:
        print(f"   ⚠️ Warning: No valid LULC data after clipping")
    
    return {
        'raster': clipped_lulc,
        'path': clipped_lulc_path,
        'name': 'IO_LULC',
        'description': 'Land Use Land Cover (native 10m)',
        'resolution': '10m'
    }


def query_and_clip_lulc_native_og(
    geojson_path: str, 
    polygon_gdf: gpd.GeoDataFrame,
    max_retries: int = 2
) -> Optional[Dict]:
    """
    Query and clip LULC at native 10m resolution WITHOUT resampling.
    Includes retry logic for VRT mosaicking failures.
    """
    print(f"\n3b. Querying and clipping IO_LULC at native resolution...")
    
    geo_db = Geospatial_DB(enable_s3=True)
    
    lulc_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product='IO_LULC',
        start_date='2015-01-01',
        end_date='2025-01-01'
    )
    
    if lulc_results is None or len(lulc_results) == 0:
        print(f"   ⚠️ No IO_LULC data found for this region")
        return None
    
    print(f"   ✅ Found {len(lulc_results)} IO_LULC tiles")
    
    lulc_uris = lulc_results['URI'].tolist()
    
    if len(lulc_uris) > 1:
        print(f"   Mosaicking {len(lulc_uris)} LULC tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in lulc_uris]
        temp_mosaic = f"/tmp/IO_LULC_mosaic_native_{os.getpid()}.tif"
        
        lulc_ds = None
        for attempt in range(max_retries):
            try:
                lulc_ds = mosaic_and_fn(
                    input_filelist=vsi_paths,
                    output_pth=temp_mosaic,
                    fn="max", 
                    dtype="float"
                )
                
                # Validate the dataset was created successfully
                if lulc_ds is not None:
                    band = lulc_ds.GetRasterBand(1)
                    if band is not None:
                        print(f"   ✅ Mosaic successful (attempt {attempt + 1}/{max_retries})")
                        break
                    else:
                        print(f"   ⚠️ Mosaic dataset invalid (attempt {attempt + 1}/{max_retries})")
                        lulc_ds = None
                else:
                    print(f"   ⚠️ Mosaic returned None (attempt {attempt + 1}/{max_retries})")
                
            except Exception as e:
                print(f"   ⚠️ Mosaic attempt {attempt + 1}/{max_retries} failed: {e}")
                lulc_ds = None
            
            if lulc_ds is None and attempt < max_retries - 1:
                print(f"   Retrying in 2 seconds...")
                import time
                time.sleep(2)
                # Clean up temp file before retry
                try:
                    if os.path.exists(temp_mosaic):
                        os.remove(temp_mosaic)
                except:
                    pass
        
        if lulc_ds is None:
            print(f"   ❌ Failed to mosaic LULC after {max_retries} attempts")
            return None
            
    else:
        # Single tile - direct open
        lulc_uri = lulc_uris[0]
        vsi_path = lulc_uri.replace('s3://', '/vsis3/')
        
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        lulc_ds = gdal.Open(vsi_path)
    
    if lulc_ds is None:
        print(f"   ❌ Failed to open IO_LULC")
        return None
    
    print(f"   Opened LULC: {lulc_ds.RasterXSize} x {lulc_ds.RasterYSize} pixels at native 10m resolution")
    
    # Clip to polygon
    clipped_lulc_path = f"/tmp/IO_LULC_native_clipped_{os.getpid()}.tif"
    
    try:
        clipped_lulc = clip_raster_to_vector(
            gd_raster=lulc_ds,
            gpd_shp=polygon_gdf,
            out_pth=clipped_lulc_path,
            crop=True,
            alltouched=False,
            crs=4326
        )
    except Exception as e:
        print(f"   ❌ Failed to clip LULC: {e}")
        lulc_ds = None
        return None
    
    lulc_ds = None  # Close source
    
    if clipped_lulc is None:
        print(f"   ❌ Failed to clip LULC")
        return None
    
    print(f"   ✅ LULC clipped at native 10m resolution")
    
    return {
        'raster': clipped_lulc,
        'path': clipped_lulc_path,
        'name': 'IO_LULC',
        'description': 'Land Use Land Cover (native 10m)',
        'resolution': '10m'
    }


def resample_lulc_feature_to_reference(
    feature_array: np.ndarray,
    native_lulc_ds: gdal.Dataset,
    reference_ds: gdal.Dataset,
    feature_name: str
) -> np.ndarray:
    """
    Resample a LULC-derived feature array from native 10m to reference resolution.
    Uses nearest neighbor to preserve categorical values.
    """
    # Create temporary raster with the feature array
    temp_path = f"/tmp/temp_lulc_feature_{feature_name}.tif"
    
    driver = gdal.GetDriverByName('GTiff')
    temp_ds = driver.Create(
        temp_path,
        native_lulc_ds.RasterXSize,
        native_lulc_ds.RasterYSize,
        1,
        gdal.GDT_Float32
    )
    
    temp_ds.SetGeoTransform(native_lulc_ds.GetGeoTransform())
    temp_ds.SetProjection(native_lulc_ds.GetProjection())
    
    band = temp_ds.GetRasterBand(1)
    band.WriteArray(feature_array)
    
    # Set nodata if present
    native_band = native_lulc_ds.GetRasterBand(1)
    nodata = native_band.GetNoDataValue()
    if nodata is not None:
        band.SetNoDataValue(nodata)
    
    temp_ds.FlushCache()
    temp_ds = None
    
    # Reopen for reading
    temp_ds = gdal.Open(temp_path)
    
    # Resample to reference using nearest neighbor (categorical data)
    resampled_path = f"/tmp/resampled_lulc_feature_{feature_name}.tif"
    resampled_ds = resample_raster_to_raster(
        gd_raster=temp_ds,
        ref_raster=reference_ds,
        interp='near',  # Preserve categorical values
        out_pth=resampled_path
    )
    
    if resampled_ds is None:
        print(f"      ⚠️ Failed to resample {feature_name}")
        temp_ds = None
        try:
            os.remove(temp_path)
        except:
            pass
        return None
    
    # Read resampled array
    resampled_array = resampled_ds.GetRasterBand(1).ReadAsArray()
    
    # Cleanup
    temp_ds = None
    resampled_ds = None
    try:
        os.remove(temp_path)
        os.remove(resampled_path)
    except:
        pass
    
    return resampled_array


def _extract_with_dem_reference(
    polygon_gdf: gpd.GeoDataFrame,
    geojson_path: str,
    metadata: Dict,
    all_features: Dict,
    product: str,
    output_dir: str
) -> Dict:
    """
    Extract raster data using DEM as reference (30m resolution).
    Now processes LULC features at native 10m resolution for consistency.
    """
    # Step 4: Query NASADEM
    print(f"\n4. Querying {product} from Geospatial Asset Catalog...")
    geo_db = Geospatial_DB(enable_s3=True)
    
    dem_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product=product,
        start_date='2000-01-01',
        end_date='2025-01-01'
    )
    
    if dem_results is None or len(dem_results) == 0:
        raise ValueError(f"No {product} data found for this region")
    
    print(f"   ✅ Found {len(dem_results)} {product} tiles")
    
    dem_uris = dem_results['URI'].tolist()
    
    if len(dem_uris) > 1:
        print(f"   Mosaicking {len(dem_uris)} tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in dem_uris]
        temp_mosaic = f"/tmp/{product}_mosaic_temp.tif"
        dem_ds = mosaic_and_fn(
            input_filelist=vsi_paths,
            output_pth=temp_mosaic,
            fn="max",
            dtype="float"
        )
    else:
        dem_uri = dem_uris[0]
        print(f"   Opening single tile: {dem_uri}")
        vsi_path = dem_uri.replace('s3://', '/vsis3/')
        
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        gdal.SetConfigOption('CPL_VSIL_CURL_CHUNK_SIZE', '524288')
        
        dem_ds = gdal.Open(vsi_path)
        
        if dem_ds is None:
            raise RuntimeError(f"Failed to open {product}")
    
    print(f"   ✅ DEM opened: {dem_ds.RasterXSize} x {dem_ds.RasterYSize} pixels")
    
    # Step 5: Derive features
    print(f"\n5. Deriving enhanced DEM features from full tile...")
    derived_features = derive_dem_features(dem_ds)
    print(f"   ✅ Derived {len(derived_features)-1} feature layers")
    
    # Step 6: Clip DEM to polygon
    print(f"\n6. Clipping DEM to polygon boundary...")
    
    base_filename = f"{metadata['feature_type']}_Lon{metadata['longitude']}_Lat{metadata['latitude']}"
    clipped_dem_path = f"/tmp/{product}_{base_filename}.tif"
    
    clipped_dem_ds = clip_raster_to_vector(
        gd_raster=dem_ds,
        gpd_shp=polygon_gdf,
        out_pth=clipped_dem_path,
        crop=True,
        alltouched=False,
        crs=4326
    )
    
    print(f"   ✅ DEM clipped")
    
    # NEW: Query and clip LULC at native 10m resolution (don't resample yet)
    lulc_native = query_and_clip_lulc_native(geojson_path, polygon_gdf)
    
    if lulc_native:
        # Derive LULC features at native 10m resolution
        print(f"\n   Deriving LULC neighborhood features at native 10m resolution...")
        print(f"   (This ensures consistent spatial scales regardless of reference choice)")
        lulc_features = derive_lulc_features(lulc_native['raster'])
        
        # Now resample the base LULC to DEM reference
        print(f"\n   Resampling LULC from 10m to 30m to match DEM reference...")
        aligned_lulc_path = f"/tmp/IO_LULC_aligned_to_dem.tif"
        aligned_lulc = resample_raster_to_raster(
            gd_raster=lulc_native['raster'],
            ref_raster=clipped_dem_ds,
            interp='near',  # Preserve categorical values
            out_pth=aligned_lulc_path
        )
        
        if aligned_lulc:
            all_features['IO_LULC'] = {
                'raster': aligned_lulc,
                'path': aligned_lulc_path,
                'name': 'IO_LULC',
                'description': 'Land Use Land Cover (resampled to 30m)'
            }
            print(f"   ✅ Base LULC resampled to 30m")
        
        # Resample each derived LULC feature to DEM reference
        print(f"\n   Resampling LULC-derived features from 10m to 30m...")
        for feature_name, feature_data in lulc_features.items():
            resampled_array = resample_lulc_feature_to_reference(
                feature_data['array'],
                lulc_native['raster'],
                clipped_dem_ds,
                feature_name
            )
            
            if resampled_array is not None:
                all_features[feature_name] = {
                    'array': resampled_array,
                    'window_size': feature_data['window_size'],
                    'description': feature_data['description'] + ' (computed at 10m, resampled to 30m)',
                    'type': feature_data['type']
                }
        
        print(f"   ✅ All LULC features resampled to match DEM reference")
    
    # Process soil and HAND features (unchanged)
    soil_features = query_and_process_soil_features(geojson_path, polygon_gdf, clipped_dem_ds)
    all_features.update(soil_features)
    
    hand_data = query_and_process_hand(geojson_path, polygon_gdf, clipped_dem_ds)
    if hand_data:
        all_features['HAND'] = hand_data
    
    # Process derived DEM features
    print(f"\n   Clipping DEM-derived features to polygon...")
    clipped_metadata = {
        'geotransform': clipped_dem_ds.GetGeoTransform(),
        'projection': clipped_dem_ds.GetProjection(),
        'width': clipped_dem_ds.RasterXSize,
        'height': clipped_dem_ds.RasterYSize,
        'nodata': derived_features['metadata'].get('nodata')
    }
    
    _process_derived_features(
        derived_features, all_features, polygon_gdf, 
        clipped_dem_ds, clipped_metadata
    )
    
    # Create multi-band raster (standard TIFF)
    print(f"\n7. Creating multi-band raster...")
    multiband_path = f"/tmp/multiband_{base_filename}.tif"
    multiband_ds, band_names = create_multiband_raster(
        all_features, multiband_path, clipped_metadata, use_bigtiff=False
    )
    
    return _finalize_output(
        multiband_ds, multiband_path, band_names, base_filename, 
        output_dir, geojson_path, all_features, dem_results
    )

def _extract_with_lulc_reference(
    polygon_gdf: gpd.GeoDataFrame,
    geojson_path: str,
    metadata: Dict,
    all_features: Dict,
    product: str,
    output_dir: str
) -> Dict:
    """
    Extract raster data using LULC as reference (10m resolution).
    This will produce larger files requiring BigTIFF format.
    """
    # Step 3b: Query LULC FIRST (becomes reference)
    print(f"\n3. Querying IO_LULC as reference raster...")
    geo_db = Geospatial_DB(enable_s3=True)
    
    lulc_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product='IO_LULC',
        start_date='2015-01-01',
        end_date='2025-01-01'
    )
    
    if lulc_results is None or len(lulc_results) == 0:
        raise ValueError("No IO_LULC data found for this region - cannot use as reference")
    
    print(f"   ✅ Found {len(lulc_results)} IO_LULC tiles")
    
    # Process LULC
    lulc_uris = lulc_results['URI'].tolist()
    
    if len(lulc_uris) > 1:
        print(f"   Mosaicking {len(lulc_uris)} LULC tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in lulc_uris]
        temp_mosaic = f"/tmp/IO_LULC_mosaic_temp.tif"
        
        lulc_ds = mosaic_and_fn(
            input_filelist=vsi_paths,
            output_pth=temp_mosaic,
            fn="max",
            dtype="float"
        )
    else:
        lulc_uri = lulc_uris[0]
        vsi_path = lulc_uri.replace('s3://', '/vsis3/')
        
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        lulc_ds = gdal.Open(vsi_path)
    
    if lulc_ds is None:
        raise RuntimeError("Failed to open IO_LULC")
    
    print(f"   ✅ LULC opened: {lulc_ds.RasterXSize} x {lulc_ds.RasterYSize} pixels")
    
    # Clip LULC to polygon - THIS BECOMES THE REFERENCE
    print(f"\n4. Clipping LULC to polygon (establishing reference grid)...")
    base_filename = f"{metadata['feature_type']}_Lon{metadata['longitude']}_Lat{metadata['latitude']}"
    clipped_lulc_path = f"/tmp/IO_LULC_reference_{base_filename}.tif"
    
    clipped_lulc_ds = clip_raster_to_vector(
        gd_raster=lulc_ds,
        gpd_shp=polygon_gdf,
        out_pth=clipped_lulc_path,
        crop=True,
        alltouched=False,
        crs=4326
    )
    
    if clipped_lulc_ds is None:
        raise RuntimeError("Failed to clip LULC reference raster")
    
    print(f"   ✅ LULC reference clipped: {clipped_lulc_ds.RasterXSize} x {clipped_lulc_ds.RasterYSize} pixels")
    print(f"   ⚠️  High resolution detected - will use BigTIFF format for outputs")
    
    # Store LULC as feature
    lulc_band = clipped_lulc_ds.GetRasterBand(1)
    lulc_array = lulc_band.ReadAsArray()
    all_features['IO_LULC'] = {
        'raster': clipped_lulc_ds,
        'path': clipped_lulc_path,
        'name': 'IO_LULC',
        'description': 'Land Use Land Cover (Reference)'
    }
    
    # NEW: Derive LULC neighborhood features from the clipped LULC reference
    print(f"\n   Deriving LULC neighborhood features...")
    lulc_features = derive_lulc_features(clipped_lulc_ds)
    
    # Add LULC features to all_features
    for feature_name, feature_data in lulc_features.items():
        all_features[feature_name] = feature_data
    
    # Step 5: Query DEM and resample to LULC reference
    print(f"\n5. Querying {product} and resampling to LULC resolution...")
    
    dem_results = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product=product,
        start_date='2000-01-01',
        end_date='2025-01-01'
    )
    
    if dem_results is None or len(dem_results) == 0:
        raise ValueError(f"No {product} data found for this region")
    
    print(f"   ✅ Found {len(dem_results)} {product} tiles")
    
    dem_uris = dem_results['URI'].tolist()
    
    if len(dem_uris) > 1:
        print(f"   Mosaicking {len(dem_uris)} DEM tiles...")
        vsi_paths = [uri.replace('s3://', '/vsis3/') for uri in dem_uris]
        temp_mosaic = f"/tmp/{product}_mosaic_temp.tif"
        dem_ds = mosaic_and_fn(
            input_filelist=vsi_paths,
            output_pth=temp_mosaic,
            fn="max",
            dtype="float"
        )
    else:
        dem_uri = dem_uris[0]
        vsi_path = dem_uri.replace('s3://', '/vsis3/')
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        dem_ds = gdal.Open(vsi_path)
    
    if dem_ds is None:
        raise RuntimeError(f"Failed to open {product}")
    
    # Derive features from full-resolution DEM FIRST
    print(f"\n6. Deriving DEM features at original resolution...")
    derived_features = derive_dem_features(dem_ds)
    print(f"   ✅ Derived {len(derived_features)-1} feature layers")
    
    # Clip and resample DEM-derived features to LULC reference
    print(f"\n7. Resampling DEM features to LULC resolution...")
    
    reference_metadata = {
        'geotransform': clipped_lulc_ds.GetGeoTransform(),
        'projection': clipped_lulc_ds.GetProjection(),
        'width': clipped_lulc_ds.RasterXSize,
        'height': clipped_lulc_ds.RasterYSize,
        'nodata': derived_features['metadata'].get('nodata')
    }
    
    _process_derived_features(
        derived_features, all_features, polygon_gdf,
        clipped_lulc_ds, reference_metadata
    )
    
    # Query other datasets and align to LULC reference
    print(f"\n8. Processing additional datasets aligned to LULC...")
    
    soil_features = query_and_process_soil_features(geojson_path, polygon_gdf, clipped_lulc_ds)
    all_features.update(soil_features)
    
    hand_data = query_and_process_hand(geojson_path, polygon_gdf, clipped_lulc_ds)
    if hand_data:
        all_features['HAND'] = hand_data
    
    # Create multi-band raster with BigTIFF
    print(f"\n9. Creating multi-band raster (BigTIFF format)...")
    multiband_path = f"/tmp/multiband_{base_filename}.tif"
    multiband_ds, band_names = create_multiband_raster(
        all_features, multiband_path, reference_metadata, use_bigtiff=True
    )
    
    return _finalize_output(
        multiband_ds, multiband_path, band_names, base_filename,
        output_dir, geojson_path, all_features, dem_results
    )

def _process_derived_features(
    derived_features: Dict,
    all_features: Dict,
    polygon_gdf: gpd.GeoDataFrame,
    reference_ds: gdal.Dataset,
    reference_metadata: Dict
) -> None:
    """
    Clip and align derived DEM features to reference raster.
    """
    for feature_name, feature_data in derived_features.items():
        if feature_name == 'metadata':
            continue
        
        temp_path = f"/tmp/temp_{feature_name}.tif"
        driver = gdal.GetDriverByName('GTiff')
        temp_ds = driver.Create(
            temp_path,
            derived_features['metadata']['width'],
            derived_features['metadata']['height'],
            1,
            gdal.GDT_Float32
        )
        
        temp_ds.SetGeoTransform(derived_features['metadata']['geotransform'])
        temp_ds.SetProjection(derived_features['metadata']['projection'])
        temp_ds.GetRasterBand(1).WriteArray(feature_data['array'])
        
        if derived_features['metadata']['nodata'] is not None:
            temp_ds.GetRasterBand(1).SetNoDataValue(derived_features['metadata']['nodata'])
        
        temp_ds.FlushCache()
        temp_ds = None
        temp_ds = gdal.Open(temp_path)
        
        # Clip to polygon
        clipped_path = f"/tmp/clipped_{feature_name}.tif"
        clipped_feature = clip_raster_to_vector(
            gd_raster=temp_ds,
            gpd_shp=polygon_gdf,
            out_pth=clipped_path,
            crop=True,
            alltouched=False,
            crs=4326
        )
        
        if clipped_feature:
            # Resample to reference resolution
            aligned_path = f"/tmp/aligned_{feature_name}.tif"
            aligned_feature = resample_raster_to_raster(
                gd_raster=clipped_feature,
                ref_raster=reference_ds,
                interp='bilinear',
                out_pth=aligned_path
            )
            
            if aligned_feature:
                clipped_array = aligned_feature.GetRasterBand(1).ReadAsArray()
                all_features[feature_name] = {
                    'array': clipped_array,
                    'window_size': feature_data.get('window_size', 1),
                    'description': feature_data['description'],
                    'type': feature_data.get('type', 'unknown')
                }
            
            aligned_feature = None
        
        clipped_feature = None
        temp_ds = None
    
    print(f"   ✅ All DEM features aligned to reference")


def run_complete_workflow(
    longitude: float, 
    latitude: float, 
    output_dir: str = "s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/rasters",
    reference_raster_type: str = 'DEM'  # NEW parameter
):
    """
    Run the complete raster extraction and processing workflow.
    Now includes HAND data processing and choice of reference raster.
    
    Parameters:
    -----------
    longitude : float
        Longitude coordinate
    latitude : float
        Latitude coordinate
    output_dir : str
        Output directory (local or S3)
    reference_raster_type : str
        Which raster to use as reference: 'DEM' (30m) or 'LULC' (10m)
        Default: 'DEM'
    """
    print("="*60)
    print("Starting Multi-Raster Extraction and Processing Workflow")
    print("Including: DEM, LULC, Soils, and HAND data")
    print(f"Reference Raster: {reference_raster_type}")
    print("="*60)
    
    watershed_data = prepare_watershed_polygon(longitude, latitude, polygon_type='subwatershed')
    results = extract_raster_data_for_polygon(
        watershed_data, 
        output_dir=output_dir,
        reference_raster_type=reference_raster_type
    )
    
    print("\n" + "="*60)
    print("Workflow Complete!")
    print("="*60)
    print(f"\nOutputs saved to: {output_dir}")
    print(f"Base filename: {results['base_filename']}")
    print(f"Reference raster: {reference_raster_type}")
    print(f"\nFiles created:")
    print(f"  - Multi-band raster: {results['multiband_path']}")
    print(f"    Contains {results['n_bands']} bands:")
    
    # Group bands by type for cleaner display
    elevation_bands = [b for b in results['band_names'] if 'elevation' in b.lower()]
    hand_bands = [b for b in results['band_names'] if b == 'HAND']
    lulc_bands = [b for b in results['band_names'] if 'lulc' in b.lower()]
    soil_bands = [b for b in results['band_names'] if b in SOIL_FEATURES]
    
    band_counter = 1
    if elevation_bands:
        print(f"    Elevation features ({len(elevation_bands)} bands):")
        for band in elevation_bands:
            print(f"      Band {band_counter}: {band}")
            band_counter += 1
    
    if hand_bands:
        print(f"    HAND features ({len(hand_bands)} band):")
        for band in hand_bands:
            print(f"      Band {band_counter}: {band}")
            band_counter += 1
    
    if lulc_bands:
        print(f"    Land use/Land cover ({len(lulc_bands)} band):")
        for band in lulc_bands:
            print(f"      Band {band_counter}: {band}")
            band_counter += 1
    
    if soil_bands:
        print(f"    Soil features ({len(soil_bands)} bands):")
        for band in soil_bands:
            print(f"      Band {band_counter}: {band}")
            band_counter += 1
    
    print(f"  - Plot: {results['plot_path']}")
    
    # Data availability summary
    print(f"\nData availability:")
    print(f"  ✅ DEM data: Found")
    print(f"  {'✅' if results['has_hand'] else '⚠️'} HAND data: {'Found' if results['has_hand'] else 'Not available'}")
    print(f"  {'✅' if results['has_lulc'] else '⚠️'} LULC data: {'Found' if results['has_lulc'] else 'Not available'}")
    print(f"  {'✅' if results['n_soil_features'] > 0 else '⚠️'} Soil data: {results['n_soil_features']} features found")
    
    return results




############################################################################################################################
#### STEP 2: process opera data day-by-day

def _normalize_dswx_classes(raw: np.ndarray) -> np.ndarray:
    """
    Normalize DSWx-S1 class raster to canonical classes:
      0 = Not water
      1 = Open water
      2 = Partial water
      3 = High-confidence water
      255 = NoData

    Handles common coding variants robustly (0..3 or 1..4).
    """
    a = raw.astype(np.int32, copy=True)
    out = np.full(a.shape, 255, dtype=np.uint8)

    valid = (a >= 0) & (a <= 250)
    if not valid.any():
        return out

    vals = np.unique(a[valid])

    # Common encodings
    if set(vals).issubset({0, 1, 2, 3}):
        mapping = {0: 0, 1: 1, 2: 2, 3: 3}
    elif set(vals).issubset({1, 2, 3, 4}):
        mapping = {1: 0, 2: 1, 3: 2, 4: 3}
    else:
        # Fallback: map min to land (0), subsequent unique values to 1..3
        mapping = {}
        sorted_vals = sorted([int(v) for v in vals if v <= 4])
        if len(sorted_vals) == 0:
            return out
        mapping[sorted_vals[0]] = 0
        next_class = 1
        for v in sorted_vals[1:]:
            mapping[v] = min(3, next_class)
            next_class += 1

    for v, m in mapping.items():
        out[a == v] = np.uint8(m)

    return out


def _resample_dswx_preserving_partial(src_ds: gdal.Dataset, ref_ds: gdal.Dataset, out_path: str) -> Optional[gdal.Dataset]:
    """
    Resample src_ds onto the exact grid of ref_ds (same width/height/gt/proj).
    - If source is finer (downsampling): try 'mode' (category-preserving), else fallback to 'near'
    - Otherwise (same/coarser): use 'near'
    """
    try:
        src_gt = src_ds.GetGeoTransform()
        ref_gt = ref_ds.GetGeoTransform()
        src_px = abs(src_gt[1])
        ref_px = abs(ref_gt[1])
    except Exception:
        src_px = 0.0
        ref_px = 1.0

    prefer_mode = (src_px > 0) and (ref_px > 0) and (src_px < ref_px * 0.9)

    # Try mode where available
    if prefer_mode:
        try:
            ds = resample_raster_to_raster(
                gd_raster=src_ds,
                ref_raster=ref_ds,
                interp='mode',
                out_pth=out_path
            )
            if ds is not None:
                return ds
        except Exception:
            pass

    # Fallback to nearest neighbor
    try:
        ds = resample_raster_to_raster(
            gd_raster=src_ds,
            ref_raster=ref_ds,
            interp='near',
            out_pth=out_path
        )
        return ds
    except Exception:
        return None


def _process_tile_with_reference(uri: str, polygon_gdf: gpd.GeoDataFrame,
                                 ref_ds: gdal.Dataset, date_str: str) -> Optional[np.ndarray]:
    """
    Clip a single DSWx tile to the polygon and align it onto ref_ds grid.
    Adds diagnostics when DSWX_DEBUG=True.
    """
    try:
        vsi_path = uri.replace('s3://', '/vsis3/')
        gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'NO')
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')

        tile_ds = gdal.Open(vsi_path, gdal.GA_ReadOnly)
        if tile_ds is None:
            if DSWX_DEBUG:
                print(f"      [tile] Could not open: {uri}")
            return None

        clipped_path = f"/tmp/dswx_{date_str}_clipped.tif"
        clipped_ds = clip_raster_to_vector(
            gd_raster=tile_ds,
            gpd_shp=polygon_gdf,
            out_pth=clipped_path,
            crop=True,
            alltouched=True,  # include narrow edges along polygon
            crs=4326
        )
        tile_ds = None
        if clipped_ds is None:
            if DSWX_DEBUG:
                print(f"      [tile] Clip returned None for: {uri}")
            try:
                os.remove(clipped_path)
            except Exception:
                pass
            return None

        if DSWX_DEBUG:
            try:
                pre = clipped_ds.GetRasterBand(1).ReadAsArray()
                uq = np.unique(pre)
                print(f"      [tile] Post-clip uniques: {uq[:10]}{'..' if uq.size > 10 else ''}")
            except Exception:
                pass

        aligned_path = f"/tmp/dswx_{date_str}_aligned.tif"
        aligned_ds = _resample_dswx_preserving_partial(clipped_ds, ref_ds, aligned_path)
        clipped_ds = None
        if aligned_ds is None:
            if DSWX_DEBUG:
                print(f"      [tile] Resample failed for: {uri}")
            try:
                os.remove(clipped_path)
            except Exception:
                pass
            try:
                os.remove(aligned_path)
            except Exception:
                pass
            return None

        if (aligned_ds.RasterXSize != ref_ds.RasterXSize) or (aligned_ds.RasterYSize != ref_ds.RasterYSize):
            if DSWX_DEBUG:
                print(f"      [tile] Shape mismatch: {aligned_ds.RasterXSize}x{aligned_ds.RasterYSize} vs ref {ref_ds.RasterXSize}x{ref_ds.RasterYSize}")
            aligned_ds = None
            try:
                os.remove(clipped_path)
            except Exception:
                pass
            try:
                os.remove(aligned_path)
            except Exception:
                pass
            return None

        arr = aligned_ds.GetRasterBand(1).ReadAsArray()
        aligned_ds = None

        # Cleanup
        for tmp in (clipped_path, aligned_path):
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass

        arr_norm = _normalize_dswx_classes(arr)
        if DSWX_DEBUG:
            uq2 = np.unique(arr_norm)
            print(f"      [tile] Post-norm uniques: {uq2[:10]}{'..' if uq2.size > 10 else ''}")

        # If everything is nodata, return None to allow other tiles to fill
        if (arr_norm == 255).all():
            if DSWX_DEBUG:
                print(f"      [tile] All nodata after normalization: {uri}")
            return None

        return arr_norm
    except Exception as e:
        if DSWX_DEBUG:
            print(f"      [tile] Exception for {uri}: {e}")
        return None


def _mosaic_dswx_tiles_for_date(uris: List[str], polygon_gdf: gpd.GeoDataFrame,
                                ref_ds: gdal.Dataset, date_str: str,
                                max_tiles: int = 32) -> Optional[np.ndarray]:
    """
    Process and mosaic up to max_tiles DSWx tiles to cover the watershed for a single date.
    Combine by class priority: 3 > 2 > 1 > 0; keep 255 as NoData where all tiles are nodata.
    """
    h, w = ref_ds.RasterYSize, ref_ds.RasterXSize
    mosaic = np.full((h, w), 255, dtype=np.uint8)
    have_any = False

    for idx, uri in enumerate(uris[:max_tiles]):
        arr = _process_tile_with_reference(uri, polygon_gdf, ref_ds, date_str)
        if arr is None:
            continue

        both_valid = (mosaic != 255) & (arr != 255)
        mosaic[both_valid] = np.maximum(mosaic[both_valid], arr[both_valid])
        only_arr = (mosaic == 255) & (arr != 255)
        mosaic[only_arr] = arr[only_arr]
        have_any = have_any or (arr != 255).any()

    if not have_any:
        return None

    # Require minimal valid pixels (very permissive)
    if (mosaic != 255).sum() < 50:
        return None

    print(f"      Mosaicked {min(len(uris), max_tiles)} tile(s)")
    return mosaic


def _plot_dswx_bands(tiff_path: str, dates: List[str], save_path: str):
    """
    Plot DSWx bands (canonical-coded: 0 Not Water, 1 Open, 2 Partial, 3 High Conf, 255 NoData).
    """
    import matplotlib.patches as mpatches

    ds = gdal.Open(tiff_path)
    n_bands = min(len(dates), ds.RasterCount)
    n_cols = min(5, n_bands)
    n_rows = int(np.ceil(n_bands / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3))

    if n_bands == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = list(axes)
    else:
        axes = axes.flatten()

    colors = ['#8B4513', '#0000FF', '#00BFFF', '#000080', '#FFFFFF']  # land/open/partial/high/nodata
    boundaries = [-0.5, 0.5, 1.5, 2.5, 3.5, 255.5]
    cmap = plt.matplotlib.colors.ListedColormap(colors)
    norm = plt.matplotlib.colors.BoundaryNorm(boundaries, cmap.N)

    for band_idx in range(n_bands):
        ax = axes[band_idx]
        band = ds.GetRasterBand(band_idx + 1)
        arr = band.ReadAsArray()
        arr = _normalize_dswx_classes(arr)
        masked = np.ma.masked_equal(arr, 255)
        ax.imshow(masked, cmap=cmap, norm=norm, aspect='auto')
        ax.set_title(dates[band_idx], fontsize=9)
        ax.axis('off')

    for idx in range(n_bands, len(axes)):
        axes[idx].set_visible(False)

    legend_labels = ['Not Water', 'Open Water', 'Partial', 'High Conf', 'No Data']
    patches = [mpatches.Patch(color=colors[i], label=legend_labels[i]) for i in range(len(colors))]
    fig.legend(handles=patches, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.05))
    plt.suptitle('OPERA DSWx-S1 (SAR) Water Time Series', fontsize=12)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()


def process_opera_dswx_to_tiff(longitude: float, latitude: float,
                               output_dir: str = "./",
                               start_date: str = "2023-01-01",
                               end_date: str = None,
                               training_folder_path: Optional[str] = None) -> Dict:
    """
    Process OPERA DSWx-S1 (SAR) water extent into a multi-band GeoTIFF aligned to the static raster.
    Output bands: one per date, NoData=255, band description 'DSWx_S1_YYYY-MM-DD'.
    Saves TIFF to output_dir (e.g., .../floodMaps).
    """
    import pandas as pd
    from datetime import datetime

    if end_date is None:
        end_date = datetime.utcnow().strftime("%Y-%m-%d")

    print("\n" + "=" * 60)
    print("Processing OPERA DSWx-S1 (SAR) Water Extent Time Series")
    print(f"Date range: {start_date} to {end_date}")
    print("DSWx Classification (canonical):")
    print("  0 = Not Water (land)")
    print("  1 = Open Water")
    print("  2 = Partial Surface Water")
    print("  3 = High Confidence Water")
    print("  255 = NoData")
    print("=" * 60 + "\n")

    # 1) Reference static raster (target grid)
    lat_str = format_coordinate(latitude)
    lon_str = format_coordinate(longitude)
    print("1. Loading reference static raster...")
    if training_folder_path:
        static_raster_path = f"{training_folder_path}/rasters/multiband_subwatershed_Lon{lon_str}_Lat{lat_str}.tif"
    else:
        static_raster_path = f"s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/rasters/multiband_subwatershed_Lon{lon_str}_Lat{lat_str}.tif"
    print(f"   Reference raster: {static_raster_path}")

    if static_raster_path.startswith("s3://"):
        vsi = static_raster_path.replace("s3://", "/vsis3/")
        gdal.SetConfigOption("AWS_NO_SIGN_REQUEST", "NO")
        gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "YES")
        reference_ds = gdal.Open(vsi, gdal.GA_ReadOnly)
    else:
        reference_ds = gdal.Open(static_raster_path, gdal.GA_ReadOnly)

    if reference_ds is None:
        alt = static_raster_path.replace("multiband_subwatershed", "subwatershed")
        print(f"   Trying alternative: {alt}")
        if alt.startswith("s3://"):
            reference_ds = gdal.Open(alt.replace("s3://", "/vsis3/"), gdal.GA_ReadOnly)
        else:
            reference_ds = gdal.Open(alt, gdal.GA_ReadOnly)
        if reference_ds is None:
            raise FileNotFoundError("Reference static raster not found. Please run static raster extraction first.")

    ref_w = reference_ds.RasterXSize
    ref_h = reference_ds.RasterYSize
    ref_gt = reference_ds.GetGeoTransform()
    ref_proj = reference_ds.GetProjection()
    print(f"   ✅ Reference loaded: {ref_w} x {ref_h}")

    # Create on-disk empty reference (for warp alignment)
    ref_tmp = f"/tmp/dswx_reference_{lon_str}_{lat_str}.tif"
    drv = gdal.GetDriverByName("GTiff")
    ref_copy = drv.Create(ref_tmp, ref_w, ref_h, 1, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=DEFLATE"])
    ref_copy.SetGeoTransform(ref_gt)
    ref_copy.SetProjection(ref_proj)
    ref_copy.GetRasterBand(1).WriteArray(np.zeros((ref_h, ref_w), dtype=np.uint8))
    ref_copy.FlushCache()
    ref_copy = None
    ref_ds = gdal.Open(ref_tmp, gdal.GA_ReadOnly)

    # 2) Watershed polygon
    print("\n2. Preparing watershed polygon...")
    ws = prepare_watershed_polygon(longitude=longitude, latitude=latitude,
                                   polygon_type="subwatershed",
                                   training_folder_path=training_folder_path)
    polygon_gdf = ws["polygon_gdf"]
    geojson_path = ws["geojson_path"]

    # 3) Query DSWx-S1
    print("\n3. Querying OPERA DSWx-S1 (SAR) data...")
    geo_db = Geospatial_DB(enable_s3=True)
    dswx = geo_db.query_asset_catalog(
        geojson_pth=geojson_path,
        product="OPERA_DSWX_S1",
        start_date=start_date,
        end_date=end_date
    )
    if dswx is None or len(dswx) == 0:
        print("   ⚠️ No OPERA DSWx-S1 data found")
        reference_ds = None
        ref_ds = None
        return {"status": "no_data"}

    print(f"   ✅ Found {len(dswx)} DSWx tiles")

    # 4) Prepare dates and iterate
    print("\n4. Filtering and grouping tiles by date...")
    dswx["date"] = pd.to_datetime(dswx["START_DATE"]).dt.date
    dswx = dswx[dswx["date"].notna()].copy()
    dates = sorted(dswx["date"].unique())
    print(f"   Dates to process: {len(dates)}")

    # 5) Process dates with MOSAIC of all intersecting tiles
    print("\n5. Processing dates with improved tile mosaicking...")
    processed_bands: List[np.ndarray] = []
    band_dates: List[str] = []
    skipped = 0

    # Permissive coverage thresholds (tuneable)
    MIN_VALID_PIXELS = 200       # allow sparse passes
    MIN_VALID_RATIO = 0.002      # 0.2% of ref grid

    for i, d in enumerate(dates, start=1):
        if i % 10 == 0 or i == 1:
            print(f"   Date {i}/{len(dates)}: {d}")

        tiles = dswx[dswx["date"] == d]
        all_uris = tiles["URI"].tolist()
        use_uris = all_uris  # Use all accessible URIs; do not restrict to a single bucket

        arr = _mosaic_dswx_tiles_for_date(use_uris, polygon_gdf, ref_ds, str(d), max_tiles=32)
        if arr is None or arr.shape != (ref_h, ref_w):
            skipped += 1
            continue

        valid_mask = (arr >= 0) & (arr <= 3)
        valid_count = int(valid_mask.sum())
        min_required = max(MIN_VALID_PIXELS, int(ref_h * ref_w * MIN_VALID_RATIO))
        if valid_count < min_required:
            print(f"      ⚠️ Insufficient coverage: {valid_count}/{ref_h * ref_w}")
            skipped += 1
            continue

        # Quick diagnostics
        water_mask = (arr >= 1) & (arr <= 3)
        partial_count = int((arr == 2).sum())
        if valid_count > 0 and (i % 10 == 0 or partial_count > 0):
            water_pct = 100.0 * water_mask.sum() / valid_count
            print(f"      Valid: {valid_count:,} ({water_pct:.1f}% water) | Partial: {partial_count:,}")

        processed_bands.append(arr.astype(np.uint8))
        band_dates.append(str(d))

    reference_ds = None
    ref_ds = None
    if len(processed_bands) == 0:
        print("   ❌ No valid date produced a usable array")
        return {"status": "processing_failed", "message": "No valid DSWx pixels inside watershed"}

    print(f"   ✅ Processed {len(processed_bands)} dates; skipped {skipped}")

    # 6) Write multi-band DSWx TIFF (Byte), nodata=255, band desc = 'DSWx_S1_YYYY-MM-DD'
    print("\n6. Writing DSWx multi-band raster...")
    base_filename = f"dswx_s1_timeseries_subwatershed_Lon{lon_str}_Lat{lat_str}"
    local_tif = f"/tmp/{base_filename}.tif"
    drv = gdal.GetDriverByName("GTiff")
    out = drv.Create(local_tif, ref_w, ref_h, len(processed_bands), gdal.GDT_Byte,
                     options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"])
    out.SetGeoTransform(ref_gt)
    out.SetProjection(ref_proj)

    for b_idx, (arr, d) in enumerate(zip(processed_bands, band_dates), start=1):
        band = out.GetRasterBand(b_idx)
        band.WriteArray(arr)
        band.SetNoDataValue(255)
        band.SetDescription(f"DSWx_S1_{d}")
        band.FlushCache()

    out.FlushCache()
    out = None  # close to persist

    # 7) Plot PNG
    print("\n7. Plotting bands...")
    local_png = f"/tmp/{base_filename}_plot.png"
    _plot_dswx_bands(local_tif, band_dates, local_png)

    # 8) Upload/save outputs
    print("\n8. Saving outputs...")
    if output_dir.startswith("s3://"):
        s3 = boto3.client("s3")
        path = output_dir.replace("s3://", "").rstrip("/")
        parts = path.split("/", 1)
        bucket = parts[0]
        prefix = parts[1] if len(parts) > 1 else ""
        tif_key = f"{prefix}/{base_filename}.tif" if prefix else f"{base_filename}.tif"
        png_key = f"{prefix}/{base_filename}_plot.png" if prefix else f"{base_filename}_plot.png"
        s3.upload_file(local_tif, bucket, tif_key)
        s3.upload_file(local_png, bucket, png_key)
        final_tif = f"s3://{bucket}/{tif_key}"
        final_png = f"s3://{bucket}/{png_key}"
        # Clean local temp files
        try:
            if os.path.exists(local_tif):
                os.remove(local_tif)
            if os.path.exists(local_png):
                os.remove(local_png)
        except Exception:
            pass
        print("   ✅ Uploaded to S3 and removed local temp files")
    else:
        import shutil
        os.makedirs(output_dir, exist_ok=True)
        final_tif = os.path.join(output_dir, f"{base_filename}.tif")
        final_png = os.path.join(output_dir, f"{base_filename}_plot.png")
        shutil.move(local_tif, final_tif)
        shutil.move(local_png, final_png)
        print("   ✅ Saved locally")

    # Cleanup temp
    try:
        if os.path.exists(geojson_path):
            os.remove(geojson_path)
    except Exception:
        pass
    try:
        if os.path.exists(ref_tmp):
            os.remove(ref_tmp)
    except Exception:
        pass

    print("\n" + "=" * 60)
    print("✅ DSWx-S1 Processing Complete")
    print("=" * 60)
    print(f"Summary:")
    print(f"  - Output size: {ref_w} x {ref_h}")
    print(f"  - Dates: {band_dates[0]} .. {band_dates[-1]} ({len(band_dates)} bands)")
    print(f"  - Output: {final_tif}")

    # Close dataset handle
    try:
        ds = None
    except Exception:
        pass

    return {
        "status": "success",
        "tiff_path": final_tif,
        "plot_path": final_png,
        "band_dates": band_dates,
        "n_dates_processed": len(band_dates)
    } 