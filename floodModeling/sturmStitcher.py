"""
sturmFloodStitcher.py - STURM Flood Dataset Stitcher (Memory-Efficient)

Stitches STURM flood data with spatial filtering and memory-efficient processing.
"""

import os
import re
import gc
import tempfile
from typing import Dict, List, Optional, Tuple
import numpy as np
from osgeo import gdal, osr
import boto3

gdal.UseExceptions()


# ======================== GDAL/PROJ Configuration ========================

def configure_gdal_proj():
    """Configure GDAL and PROJ library paths."""
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
    
    for cand in ["/opt/conda/share/proj", "/usr/share/proj", "/usr/local/share/proj"]:
        if os.path.isdir(cand):
            os.environ["PROJ_LIB"] = cand
            gdal.SetConfigOption("PROJ_LIB", cand)
            print(f"  Using PROJ from: {cand}")
            return
    
    print("  ⚠️  Could not locate PROJ data directory")

configure_gdal_proj()
gdal.SetConfigOption("PROJ_NETWORK", "ON")
gdal.SetConfigOption("GDAL_CACHEMAX", "256")  # Increased from 64 for better performance
gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,tiff")
gdal.SetConfigOption("AWS_NO_SIGN_REQUEST", "NO")


# ======================== Spatial Filtering Functions ========================

def get_bbox_from_tiff(tiff_uri: str, target_crs: str = "EPSG:4326") -> Optional[Tuple[float, float, float, float]]:
    """Extract bounding box from TIFF and transform to target CRS."""
    vsis3_path = tiff_uri.replace('s3://', '/vsis3/')
    
    try:
        ds = gdal.Open(vsis3_path, gdal.GA_ReadOnly)
        if ds is None:
            return None
        
        gt = ds.GetGeoTransform()
        width = ds.RasterXSize
        height = ds.RasterYSize
        src_proj_wkt = ds.GetProjection()
        
        minx_src = gt[0]
        maxy_src = gt[3]
        maxx_src = minx_src + gt[1] * width
        miny_src = maxy_src + gt[5] * height
        
        ds = None
        
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt(src_proj_wkt)
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        
        tgt_srs = osr.SpatialReference()
        if target_crs.upper().startswith('EPSG:'):
            epsg_code = int(target_crs.split(':')[1])
            tgt_srs.ImportFromEPSG(epsg_code)
        else:
            tgt_srs.SetFromUserInput(target_crs)
        tgt_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        
        if src_srs.IsSame(tgt_srs):
            return (minx_src, miny_src, maxx_src, maxy_src)
        
        transform = osr.CoordinateTransformation(src_srs, tgt_srs)
        
        corners = [
            (minx_src, miny_src),
            (maxx_src, miny_src),
            (maxx_src, maxy_src),
            (minx_src, maxy_src)
        ]
        
        transformed_x = []
        transformed_y = []
        
        for x, y in corners:
            try:
                tx, ty, _ = transform.TransformPoint(x, y)
                transformed_x.append(tx)
                transformed_y.append(ty)
            except Exception as e:
                return None
        
        return (min(transformed_x), min(transformed_y), max(transformed_x), max(transformed_y))
        
    except Exception as e:
        return None


def bboxes_intersect(bbox1: Tuple[float, float, float, float],
                    bbox2: Tuple[float, float, float, float],
                    buffer: float = 0.0) -> bool:
    """Check if two bounding boxes intersect."""
    minx1, miny1, maxx1, maxy1 = bbox1
    minx2, miny2, maxx2, maxy2 = bbox2
    
    minx1 -= buffer
    miny1 -= buffer
    maxx1 += buffer
    maxy1 += buffer
    
    no_overlap = (maxx1 < minx2 or minx1 > maxx2 or maxy1 < miny2 or miny1 > maxy2)
    return not no_overlap


def filter_tiles_by_bbox(tile_uris: List[str], 
                         reference_bbox: Tuple[float, float, float, float],
                         target_crs: str = "EPSG:4326",
                         buffer: float = 0.1,
                         debug: bool = False) -> List[str]:
    """Filter tiles to only those intersecting reference bbox."""
    if debug:
        ref_minx, ref_miny, ref_maxx, ref_maxy = reference_bbox
        print(f"  Reference bbox: [{ref_minx:.4f}, {ref_miny:.4f}, {ref_maxx:.4f}, {ref_maxy:.4f}]")
        print(f"  With {buffer}° buffer: [{ref_minx-buffer:.4f}, {ref_miny-buffer:.4f}, {ref_maxx+buffer:.4f}, {ref_maxy+buffer:.4f}]")
        print(f"  Checking {len(tile_uris)} tiles for intersection...")
    
    filtered_tiles = []
    failed_tiles = 0
    
    for i, tile_uri in enumerate(tile_uris):
        tile_bbox = get_bbox_from_tiff(tile_uri, target_crs)
        
        if tile_bbox is None:
            filtered_tiles.append(tile_uri)
            failed_tiles += 1
            continue
        
        if debug and i < 3:
            tile_minx, tile_miny, tile_maxx, tile_maxy = tile_bbox
            print(f"    Tile {i+1}: [{tile_minx:.4f}, {tile_miny:.4f}, {tile_maxx:.4f}, {tile_maxy:.4f}]")
            intersects = bboxes_intersect(reference_bbox, tile_bbox, buffer=buffer)
            print(f"      Intersects: {intersects}")
        
        if bboxes_intersect(reference_bbox, tile_bbox, buffer=buffer):
            filtered_tiles.append(tile_uri)
    
    if debug:
        print(f"  Result: {len(filtered_tiles)} tiles kept, {len(tile_uris) - len(filtered_tiles)} excluded")
        if failed_tiles > 0:
            print(f"  Note: {failed_tiles} tiles kept due to bbox read failure")
    
    return filtered_tiles


# ======================== S3 Utilities ========================

def list_s3_files(bucket: str, prefix: str, pattern: str = None) -> List[str]:
    """List files in S3 matching optional pattern."""
    s3_client = boto3.client('s3')
    files = []
    
    paginator = s3_client.get_paginator('list_objects_v2')
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        if 'Contents' not in page:
            continue
        
        for obj in page['Contents']:
            key = obj['Key']
            if pattern:
                filename = os.path.basename(key)
                if not re.search(pattern, filename):
                    continue
            files.append(f"s3://{bucket}/{key}")
    
    return files


# ======================== Memory-Efficient Raster Processing ========================

def stitch_and_save_to_disk(file_list: List[str], 
                            output_uri: str,
                            data_type: str = "floodmap",  # ✅ NEW: "floodmap" or "imagery"
                            target_crs: str = "EPSG:4326",
                            resolution: Optional[float] = None,
                            nodata_value: int = 255) -> bool:
    """Stitch rasters and save directly to disk (memory efficient).
    
    Args:
        file_list: List of S3 URIs to rasters
        output_uri: S3 URI for output
        data_type: "floodmap" (integer, nearest) or "imagery" (float, bilinear)
        target_crs: Target CRS
        resolution: Target resolution
        nodata_value: NoData value
    
    Returns:
        True if successful, False otherwise
    """
    if not file_list:
        return False
    
    print(f"    Stitching {len(file_list)} {data_type} tiles...")
    
    vsis3_paths = [f.replace('s3://', '/vsis3/') for f in file_list]
    
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        # Configure based on data type
        if data_type == "imagery":
            # ✅ IMAGERY: Preserve float reflectance, use bilinear
            warp_options = gdal.WarpOptions(
                format='GTiff',
                creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
                dstSRS=target_crs,
                xRes=resolution,
                yRes=resolution,
                resampleAlg='bilinear',  # Better for continuous data
                outputType=gdal.GDT_Float32,  # ✅ Preserve float!
                dstNodata=-9999,  # Use float-compatible NoData
                multithread=True,
                warpMemoryLimit=512
            )
        else:
            # ✅ FLOODMAP: Integer classification, use nearest neighbor
            warp_options = gdal.WarpOptions(
                format='GTiff',
                creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
                dstSRS=target_crs,
                xRes=resolution,
                yRes=resolution,
                resampleAlg='nearest',  # Preserve discrete classes
                outputType=gdal.GDT_Byte,  # Integer 0-255
                dstNodata=nodata_value,
                srcNodata=99,  # STURM uses 99 for NoData
                multithread=True,
                warpMemoryLimit=512
            )
        
        result = gdal.Warp(temp_path, vsis3_paths, options=warp_options)
        
        if result is None:
            print(f"    ⚠️  Failed to create mosaic")
            os.remove(temp_path)
            return False
        
        print(f"    ✓ Created mosaic: {result.RasterXSize}x{result.RasterYSize}, {result.RasterCount} bands")
        result = None
        
        # Upload to S3
        s3_client = boto3.client('s3')
        bucket, key = output_uri.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        
        print(f"    ✓ Saved to: {output_uri}")
        gc.collect()
        
        return True
        
    except Exception as e:
        print(f"    ⚠️  Error: {e}")
        if os.path.exists(temp_path):
            os.remove(temp_path)
        return False


# ======================== Main Stitching Function ========================

def stitch_sturm_floods(
    flood_event_code: str,
    predictions_tiff_uri: Optional[str] = None,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs/curated_flood_maps/STURM_Floods",
    output_folder: Optional[str] = None,
    target_crs: str = "EPSG:4326",
    resolution: Optional[float] = None,
    bbox_buffer: float = 0.1,
    save_to_disk: bool = True
) -> Dict:
    """Stitch STURM flood maps with memory-efficient processing and proper data type handling.
    
    Args:
        flood_event_code: Flood event code (e.g., "EMSR279")
        predictions_tiff_uri: Optional predictions TIFF for spatial filtering
        bucket: S3 bucket name
        base_prefix: Base S3 prefix for STURM dataset
        output_folder: S3 folder to save stitched outputs (required if save_to_disk=True)
        target_crs: Target CRS (default EPSG:4326)
        resolution: Target resolution (auto if None)
        bbox_buffer: Buffer in degrees for spatial filtering
        save_to_disk: If True, save to disk instead of loading to memory
        
    Returns:
        Dictionary with paths to saved files or None for each dataset type
    """
    
    print(f"\n{'='*70}")
    print(f"STITCHING STURM FLOODS: {flood_event_code}")
    print(f"{'='*70}")
    print(f"Memory mode: {'DISK' if save_to_disk else 'RAM'}")
    
    if save_to_disk and output_folder is None:
        raise ValueError("output_folder must be specified when save_to_disk=True")
    
    # Extract reference bbox
    reference_bbox = None
    if predictions_tiff_uri:
        print(f"\n[0] Extracting reference bounding box...")
        reference_bbox = get_bbox_from_tiff(predictions_tiff_uri, target_crs)
        
        if reference_bbox:
            minx, miny, maxx, maxy = reference_bbox
            print(f"    Reference bbox: [{minx:.4f}, {miny:.4f}, {maxx:.4f}, {maxy:.4f}]")
            print(f"    Buffer: {bbox_buffer}°")
        else:
            print(f"    ⚠️  Could not extract bbox - no spatial filtering")
    
    results = {
        's1_floodmap': None,
        's1_imagery': None,
        's2_floodmap': None,
        's2_imagery': None,
        'metadata': {
            'event_code': flood_event_code,
            'target_crs': target_crs,
            'reference_bbox': reference_bbox,
            'bbox_buffer': bbox_buffer,
            's1_floodmap_tiles': [],
            's1_imagery_tiles': [],
            's2_floodmap_tiles': [],
            's2_imagery_tiles': [],
            'save_to_disk': save_to_disk,
            'output_folder': output_folder
        }
    }
    
    # ============================================================
    # [1/4] S1 FLOOD MAPS (Integer Classification)
    # ============================================================
    
    print("\n[1/4] Processing Sentinel-1 flood maps...")
    s1_fm_files = list_s3_files(bucket, f"{base_prefix}/Sentinel1/Floodmaps/", f"^{flood_event_code}_.*\\.tif$")
    
    if s1_fm_files:
        print(f"  Found {len(s1_fm_files)} tiles")
        
        if reference_bbox:
            print(f"  Filtering by spatial intersection...")
            s1_fm_files = filter_tiles_by_bbox(s1_fm_files, reference_bbox, target_crs, bbox_buffer, debug=True)
        
        if s1_fm_files:
            results['metadata']['s1_floodmap_tiles'] = s1_fm_files
            
            if save_to_disk:
                output_uri = f"{output_folder}/sturm_{flood_event_code}_s1_floodmap.tif"
                # ✅ Use data_type="floodmap" for integer classification
                if stitch_and_save_to_disk(s1_fm_files, output_uri, data_type="floodmap",
                                          target_crs=target_crs, resolution=resolution, nodata_value=255):
                    results['s1_floodmap'] = output_uri
            else:
                print("    ⚠️  Memory mode not recommended - use save_to_disk=True")
        else:
            print(f"  ⚠️  No tiles intersect reference area")
    else:
        print(f"  ⚠️  No S1 flood maps found")
    
    gc.collect()
    
    # ============================================================
    # [2/4] S1 IMAGERY (Float Reflectance)
    # ============================================================
    
    print("\n[2/4] Processing Sentinel-1 imagery...")
    s1_img_files = list_s3_files(bucket, f"{base_prefix}/Sentinel1/S1/", f"^{flood_event_code}_.*\\.tif$")
    
    if s1_img_files:
        print(f"  Found {len(s1_img_files)} tiles")
        
        if reference_bbox:
            print(f"  Filtering by spatial intersection...")
            s1_img_files = filter_tiles_by_bbox(s1_img_files, reference_bbox, target_crs, bbox_buffer, debug=False)
        
        if s1_img_files:
            results['metadata']['s1_imagery_tiles'] = s1_img_files
            
            if save_to_disk:
                output_uri = f"{output_folder}/sturm_{flood_event_code}_s1_imagery.tif"
                # ✅ Use data_type="imagery" for float reflectance
                if stitch_and_save_to_disk(s1_img_files, output_uri, data_type="imagery",
                                          target_crs=target_crs, resolution=resolution, nodata_value=-9999):
                    results['s1_imagery'] = output_uri
        else:
            print(f"  ⚠️  No tiles intersect reference area")
    else:
        print(f"  ⚠️  No S1 imagery found")
    
    gc.collect()
    
    # ============================================================
    # [3/4] S2 FLOOD MAPS (Integer Classification)
    # ============================================================
    
    print("\n[3/4] Processing Sentinel-2 flood maps...")
    s2_fm_files = list_s3_files(bucket, f"{base_prefix}/Sentinel2/Floodmaps/", f"^{flood_event_code}_.*\\.tif$")
    
    if s2_fm_files:
        print(f"  Found {len(s2_fm_files)} tiles")
        
        if reference_bbox:
            print(f"  Filtering by spatial intersection...")
            s2_fm_files = filter_tiles_by_bbox(s2_fm_files, reference_bbox, target_crs, bbox_buffer, debug=False)
        
        if s2_fm_files:
            results['metadata']['s2_floodmap_tiles'] = s2_fm_files
            
            if save_to_disk:
                output_uri = f"{output_folder}/sturm_{flood_event_code}_s2_floodmap.tif"
                # ✅ Use data_type="floodmap" for integer classification
                if stitch_and_save_to_disk(s2_fm_files, output_uri, data_type="floodmap",
                                          target_crs=target_crs, resolution=resolution, nodata_value=255):
                    results['s2_floodmap'] = output_uri
        else:
            print(f"  ⚠️  No tiles intersect reference area")
    else:
        print(f"  ⚠️  No S2 flood maps found")
    
    gc.collect()
    
    # ============================================================
    # [4/4] S2 IMAGERY (Float Reflectance - Multi-band RGB)
    # ============================================================
    
    print("\n[4/4] Processing Sentinel-2 imagery...")
    s2_img_files = list_s3_files(bucket, f"{base_prefix}/Sentinel2/S2/", f"^{flood_event_code}_.*\\.tif$")
    
    if s2_img_files:
        print(f"  Found {len(s2_img_files)} tiles")
        
        if reference_bbox:
            print(f"  Filtering by spatial intersection...")
            s2_img_files = filter_tiles_by_bbox(s2_img_files, reference_bbox, target_crs, bbox_buffer, debug=False)
        
        if s2_img_files:
            results['metadata']['s2_imagery_tiles'] = s2_img_files
            
            if save_to_disk:
                output_uri = f"{output_folder}/sturm_{flood_event_code}_s2_imagery.tif"
                # ✅ Use data_type="imagery" for float reflectance (multi-band)
                if stitch_and_save_to_disk(s2_img_files, output_uri, data_type="imagery",
                                          target_crs=target_crs, resolution=resolution, nodata_value=-9999):
                    results['s2_imagery'] = output_uri
        else:
            print(f"  ⚠️  No tiles intersect reference area")
    else:
        print(f"  ⚠️  No S2 imagery found")
    
    gc.collect()
    
    # ============================================================
    # SUMMARY
    # ============================================================
    
    print(f"\n{'='*70}")
    print("STITCHING COMPLETE")
    print(f"{'='*70}")
    
    for key in ['s1_floodmap', 's1_imagery', 's2_floodmap', 's2_imagery']:
        val = results[key]
        if val:
            print(f"  ✓ {key}: {val}")
        else:
            print(f"  ✗ {key}: None")
    
    return results

