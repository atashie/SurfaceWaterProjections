# processAllTrainingData.py - UPDATED WITH OPERA DSWX INTEGRATION

"""
Optimized basin-wide training data preparation pipeline with OPERA DSWx integration.
Training folders now include floodMaps/ with water extent time series for each watershed.
"""

import os
import logging
import pandas as pd
import numpy as np
from shapely.geometry import Point
from typing import List, Dict, Tuple, Optional
import boto3
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import warnings
import json
import geopandas as gpd
from coordinate_utils import format_coordinate

# Import optimized modules
from delineateWatersheds import (
    find_containing_subwatershed_optimized, 
    delineateHydroBasins, 
    read_s3_gpkg
)
from extractRasterData import (
    run_complete_workflow as extract_rasters,
    process_opera_dswx_to_tiff
)
from extractClimateData import ClimateDataExtractor

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OptimizedBasinTrainingPipeline:
    """
    Optimized training data preparation for multiple watersheds in a basin.
    Now includes OPERA DSWx water extent time series processing.
    """
    
    def __init__(self, base_output_path: str = "s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs"):
        self.base_output_path = base_output_path.rstrip('/')
        self.s3_client = boto3.client('s3')
        self._basin_data_cache = {}
        
    def define_all_watersheds_in_basin(
        self,
        latitude: float,
        longitude: float,
        n_neighbors: int = 5,
        level: int = 12,
        parallel_processing: bool = True,
        skip_existing: bool = True,
        include_dswx: bool = True,
        dswx_start_date: str = "2023-01-01",
        dswx_end_date: str = None,
        reference_raster_type: str = 'DEM'  # NEW: 'DEM' (30m) or 'LULC' (10m)
    ) -> Dict:
        """
        Main workflow with watershed-based folder naming and DSWx integration.

        Parameters:
        -----------
        latitude : float
            Latitude of the point of interest
        longitude : float
            Longitude of the point of interest
        n_neighbors : int
            Number of nearest neighbor watersheds to process
        level : int
            BasinATLAS level (default 12)
        parallel_processing : bool
            Whether to use parallel processing
        skip_existing : bool
            Whether to skip existing training folders
        include_dswx : bool
            Whether to process OPERA DSWx water extent time series
        dswx_start_date : str
            Start date for DSWx data query (YYYY-MM-DD)
        dswx_end_date : str
            End date for DSWx data query (None = today)
        reference_raster_type : str
            Reference raster for spatial alignment: 'DEM' (30m) or 'LULC' (10m)
            Default: 'DEM'
        """        
        # Format user input coordinates for display
        user_lat_str = format_coordinate(latitude)
        user_lon_str = format_coordinate(longitude)
        
        print(f"\n{'='*60}")
        print(f"OPTIMIZED BASIN TRAINING DATA PREPARATION")
        print(f"User Input Location: ({user_lat_str}, {user_lon_str})")
        print(f"Neighbors to process: {n_neighbors}")
        print(f"Reference Raster: {reference_raster_type} "
              f"({'30m resolution' if reference_raster_type == 'DEM' else '10m resolution'})")
        if include_dswx:
            print(f"OPERA DSWx Processing: ENABLED")
            print(f"  Date range: {dswx_start_date} to {dswx_end_date or 'today'}")
        print(f"{'='*60}\n")
        
        # Step 1: Find containing subwatershed
        print("\n🔍 Step 1: Finding containing subwatershed (optimized)...")
        gdf = self._load_basin_data(level)
        
        point = Point(longitude, latitude)
        initial_watershed_gdf = find_containing_subwatershed_optimized(gdf, point, use_spatial_index=True)
        
        if initial_watershed_gdf is None:
            error_msg = f"No subwatershed found containing point ({user_lat_str}, {user_lon_str})"
            logger.error(error_msg)
            return {'error': error_msg, 'status': 'error'}
        
        initial_watershed = initial_watershed_gdf.iloc[0]
        print(f"   ✅ Found subwatershed with HYBAS_ID: {initial_watershed['HYBAS_ID']}")
        
        # Get watershed representative coordinates for folder naming
        watershed_coords = self._get_watershed_coordinates(initial_watershed)
        watershed_lat = watershed_coords['latitude']
        watershed_lon = watershed_coords['longitude']
        watershed_lat_str = format_coordinate(watershed_lat)
        watershed_lon_str = format_coordinate(watershed_lon)
        
        print(f"   📍 Watershed representative location: ({watershed_lat_str}, {watershed_lon_str})")
        
        # Step 2: Check for existing training folder
        print("\n📁 Step 2: Checking for existing training folder...")
        print(f"   Note: Folder named by watershed location, not user input")
        
        training_folder = f"trainingFolderFor_Lat_{watershed_lat_str}_Lon_{watershed_lon_str}"
        training_folder_path = f"{self.base_output_path}/{training_folder}"
        
        # Store metadata about the query
        query_metadata = {
            'user_input_lat': float(user_lat_str),
            'user_input_lon': float(user_lon_str),
            'watershed_lat': float(watershed_lat_str),
            'watershed_lon': float(watershed_lon_str),
            'hybas_id': initial_watershed['HYBAS_ID'],
            'query_time': pd.Timestamp.now().isoformat(),
            'include_dswx': include_dswx,
            'dswx_start_date': dswx_start_date if include_dswx else None,
            'dswx_end_date': dswx_end_date if include_dswx else None,
            'reference_raster_type': reference_raster_type
        }

        
        if self._check_s3_folder_exists(training_folder_path):
            if skip_existing:
                print(f"   ⚠️ Training folder already exists for this watershed: {training_folder_path}")
                print(f"      User queried: ({user_lat_str}, {user_lon_str})")
                print(f"      Watershed center: ({watershed_lat_str}, {watershed_lon_str})")
                return {
                    'status': 'exists',
                    'folder': training_folder_path,
                    'training_folder': training_folder_path,
                    'query_metadata': query_metadata,
                    'message': 'Training data already exists for this watershed. Set skip_existing=False to reprocess.'
                }
            else:
                print(f"   ⚠️ Training folder exists for watershed but will be reprocessed")
        else:
            print(f"   ✅ No existing training folder found for this watershed")
        
        print(f"   Training folder path: {training_folder_path}")
        
        # Save query metadata
        self._save_query_metadata(training_folder_path, query_metadata)
        
        # Step 3: Find nearest subwatersheds
        neighbors_to_fetch = max(0, n_neighbors)
        print(f"\n🗺️ Step 3: Target total watersheds: {n_neighbors} (including the seed)")
        print(f"   Finding {neighbors_to_fetch} nearest neighbors (excluding the seed)...")

        all_watersheds = self._find_nearest_subwatersheds_optimized(
            initial_watershed, neighbors_to_fetch, gdf
        )
        
        # Extract locations with formatted coordinates
        watershed_locations = self._extract_watershed_locations_fast(all_watersheds)
        
        # Add user query info to the first watershed
        watershed_locations[0]['user_query_lat'] = float(user_lat_str)
        watershed_locations[0]['user_query_lon'] = float(user_lon_str)
        
        # Save locations CSV
        locations_df = pd.DataFrame(watershed_locations)
        locations_path = f"{training_folder_path}/list_of_training_locations.csv"
        self._save_df_to_s3(locations_df, locations_path)
        print(f"   ✅ Saved watershed locations to: {locations_path}")
        
        # Step 4: Process watersheds with DSWx integration
        print(f"\n⚙️ Step 4: Processing {len(all_watersheds)} watersheds...")
        if include_dswx:
            print(f"   Including OPERA DSWx water extent time series")
        
        if parallel_processing and len(all_watersheds) > 1:
            results = self._process_watersheds_parallel_optimized(
                watershed_locations, training_folder_path, include_dswx, 
                dswx_start_date, dswx_end_date, reference_raster_type
            )
        else:
            results = self._process_watersheds_sequential_optimized(
                watershed_locations, training_folder_path, include_dswx, 
                dswx_start_date, dswx_end_date, reference_raster_type
            )
        
        # Step 5: Create summary
        print(f"\n📊 Step 5: Creating processing summary...")
        summary = self._create_summary(
            results, watershed_locations, training_folder_path, query_metadata
        )
        
        print(f"\n{'='*60}")
        print(f"✅ BASIN TRAINING DATA PREPARATION COMPLETE")
        print(f"{'='*60}")
        print(f"\nSummary:")
        print(f"  - User queried location: ({user_lat_str}, {user_lon_str})")
        print(f"  - Watershed center: ({watershed_lat_str}, {watershed_lon_str})")
        print(f"  - Total watersheds processed: {summary['total_watersheds']}")
        print(f"  - Successful: {summary['successful']}")
        print(f"  - Failed: {summary['failed']}")
        if include_dswx:
            dswx_summary = summary.get('dswx_summary', {})
            print(f"  - DSWx processing: {dswx_summary.get('successful', 0)} successful, "
                  f"{dswx_summary.get('failed', 0)} failed, "
                  f"{dswx_summary.get('no_data', 0)} no data")
            if dswx_summary.get('total_dates_processed', 0) > 0:
                print(f"    Total dates processed: {dswx_summary['total_dates_processed']}")
        print(f"  - Output folder: {training_folder_path}")
        
        return summary
    
    def _process_single_watershed_optimized(
        self,
        location: Dict,
        training_folder_path: str,
        include_dswx: bool = True,
        dswx_start_date: str = "2023-01-01",
        dswx_end_date: str = None,
        reference_raster_type: str = 'DEM',
        climate_extractor: Optional[ClimateDataExtractor] = None
    ) -> Dict:
        """
        Process a single watershed with DSWx water extent integration.

        Parameters:
        -----------
        location : Dict
            Watershed location information
        training_folder_path : str
            Path to training folder
        include_dswx : bool
            Whether to process DSWx data
        dswx_start_date : str
            Start date for DSWx query
        dswx_end_date : str
            End date for DSWx query
        reference_raster_type : str
            Reference raster: 'DEM' (30m) or 'LULC' (10m)
        climate_extractor : Optional[ClimateDataExtractor]
            Climate data extractor (will create new if None)
        """
        lat = location['latitude']
        lon = location['longitude']
        lat_str = location.get('lat_str', format_coordinate(lat))
        lon_str = location.get('lon_str', format_coordinate(lon))
        hybas_id = location['hybas_id']
        
        result = {
            'location': location,
            'status': 'processing',
            'errors': []
        }

        max_retries = 3
        for attempt in range(max_retries):
            try:
                print(f"\n   Processing watershed {hybas_id} at ({lat_str}, {lon_str})...")
                
                # Step 1: Delineate watersheds
                print(f"      Delineating watersheds...")
                watershed_path = delineateHydroBasins(
                    latitude=lat,
                    longitude=lon,
                    output_location=f"{training_folder_path}/watershedBoundaries",
                    level=12,
                    use_spatial_index=True
                )
                result['watershed_path'] = watershed_path
                
                # Step 2: Extract climate data
                print(f"      Extracting climate data...")
                local_climate_extractor = ClimateDataExtractor(cache_masks=False)
                
                try:
                    climate_results = local_climate_extractor.process_watersheds(
                        gpkg_path=watershed_path,
                        latitude=lat,
                        longitude=lon,
                        dataset='era5',
                        output_location=f"{training_folder_path}/climateData",
                        parallel_members=False
                    )
                    result['climate_data'] = climate_results
                finally:
                    # Force cleanup of climate extractor
                    import gc
                    del local_climate_extractor
                    gc.collect()
                
                # Step 3: Create model metadata
                print(f"      Creating model metadata...")
                self._create_model_metadata(hybas_id, lat, lon, training_folder_path)
                
                # Step 4: Extract raster data
                print(f"      Extracting raster data...")
                import time
                time.sleep(0.5)  # Small delay to prevent connection overload
                
                try:
                    from extractRasterData import prepare_watershed_polygon, extract_raster_data_for_polygon
                    
                    watershed_data = prepare_watershed_polygon(
                        longitude=lon,
                        latitude=lat,
                        polygon_type='subwatershed',
                        training_folder_path=training_folder_path
                    )
                    
                    raster_results = extract_raster_data_for_polygon(
                        watershed_data,
                        output_dir=f"{training_folder_path}/rasters",
                        reference_raster_type=reference_raster_type
                    )
                    result['raster_data'] = raster_results
                except Exception as raster_error:
                    print(f"      ⚠️ Raster extraction failed: {raster_error}")
                    result['raster_data'] = {'status': 'failed', 'error': str(raster_error)}
                
                # Step 5: Process OPERA DSWx water extent time series
                if include_dswx:
                    print(f"      Processing OPERA DSWx water extent time series...")
                    try:
                        # Process DSWx data for this watershed
                        dswx_results = process_opera_dswx_to_tiff(
                            longitude=lon,
                            latitude=lat,
                            output_dir=f"{training_folder_path}/floodMaps",
                            start_date=dswx_start_date,
                            end_date=dswx_end_date,
                            training_folder_path=training_folder_path
                        )
                        
                        result['dswx_data'] = dswx_results
                        
                        if dswx_results.get('status') == 'success':
                            n_dates = dswx_results.get('n_dates_processed', 0)
                            print(f"      ✅ DSWx processing successful: {n_dates} dates")
                        else:
                            status = dswx_results.get('status', 'unknown')
                            print(f"      ⚠️ DSWx processing status: {status}")
                            
                    except Exception as dswx_error:
                        print(f"      ⚠️ DSWx processing failed: {dswx_error}")
                        result['dswx_data'] = {'status': 'failed', 'error': str(dswx_error)}
                else:
                    result['dswx_data'] = {'status': 'skipped'}
                
                result['status'] = 'success'
                print(f"      ✅ Successfully processed watershed {hybas_id}")
                break

            except Exception as e:
                result['status'] = 'failed'
                result['errors'].append(str(e))
                if "libgdal" in str(e) and attempt < max_retries - 1:
                    print(f"      ⚠️ GDAL error, retrying (attempt {attempt + 2}/{max_retries})...")
                    import time
                    time.sleep(10)  # Wait before retry
                    continue
                else:
                    raise  # Re-raise if not GDAL error or max retries reached
            
            finally:  # ← FIXED: Proper indentation (aligned with try/except)
                # Cleanup any temporary resources
                import gc
                gc.collect()
        
        return result
    
    def _process_watersheds_parallel_optimized(
        self,
        watershed_locations: List[Dict],
        training_folder_path: str,
        include_dswx: bool = True,
        dswx_start_date: str = "2023-01-01",
        dswx_end_date: str = None,
        reference_raster_type: str = 'DEM'
    ) -> List[Dict]:
        """Parallel processing with DSWx integration."""
        results = []
        max_workers = min(2, len(watershed_locations))
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_single_watershed_optimized, 
                    location, 
                    training_folder_path,
                    include_dswx,
                    dswx_start_date,
                    dswx_end_date,
                    reference_raster_type,
                    None
                ): location 
                for location in watershed_locations
            }
            
            for future in tqdm(as_completed(futures), total=len(futures), 
                             desc="Processing watersheds"):
                try:
                    result = future.result(timeout=1200)  # 20 minute timeout
                    results.append(result)
                except Exception as e:
                    location = futures[future]
                    logger.error(f"Failed processing watershed {location['hybas_id']}: {e}")
                    results.append({
                        'location': location,
                        'status': 'failed',
                        'errors': [str(e)]
                    })
        
        return results
    
    def _process_watersheds_sequential_optimized(
        self,
        watershed_locations: List[Dict],
        training_folder_path: str,
        include_dswx: bool = True,
        dswx_start_date: str = "2023-01-01",
        dswx_end_date: str = None,
        reference_raster_type: str = 'DEM'
    ) -> List[Dict]:
        """Sequential processing with DSWx integration."""
        results = []
        
        for location in tqdm(watershed_locations, desc="Processing watersheds"):
            result = self._process_single_watershed_optimized(
                location, 
                training_folder_path, 
                include_dswx, 
                dswx_start_date, 
                dswx_end_date,
                reference_raster_type,
                None
            )
            results.append(result)
            
            # Force garbage collection between watersheds
            import gc
            gc.collect()
            
            # Small delay to prevent connection overload
            import time
            time.sleep(1)
        
        return results
    
    def _create_summary(
        self,
        results: List[Dict],
        watershed_locations: List[Dict],
        training_folder_path: str,
        query_metadata: Dict
    ) -> Dict:
        """Create summary including DSWx processing results."""
        successful = sum(1 for r in results if r['status'] == 'success')
        failed = sum(1 for r in results if r['status'] == 'failed')
        
        # Count DSWx processing results
        dswx_successful = sum(1 for r in results 
                             if r.get('dswx_data', {}).get('status') == 'success')
        dswx_failed = sum(1 for r in results 
                         if r.get('dswx_data', {}).get('status') == 'failed')
        dswx_no_data = sum(1 for r in results 
                          if r.get('dswx_data', {}).get('status') in ['no_data', 'no_accessible_data'])
        dswx_skipped = sum(1 for r in results
                          if r.get('dswx_data', {}).get('status') == 'skipped')
        
        # Calculate total dates processed across all watersheds
        total_dates_processed = sum(
            r.get('dswx_data', {}).get('n_dates_processed', 0) 
            for r in results
        )
        
        summary = {
            'status': 'success' if successful > 0 else 'failed',
            'total_watersheds': len(watershed_locations),
            'successful': successful,
            'failed': failed,
            'training_folder': training_folder_path,
            'query_metadata': query_metadata,
            'watershed_locations': watershed_locations,
            'processing_results': results,
            'dswx_summary': {
                'successful': dswx_successful,
                'failed': dswx_failed,
                'no_data': dswx_no_data,
                'skipped': dswx_skipped,
                'total_dates_processed': total_dates_processed
            },
            'created_at': pd.Timestamp.now().isoformat()
        }
        
        # Save summary
        summary_json = json.dumps(summary, indent=2, default=str)
        summary_path = f"{training_folder_path}/processing_summary.json"
        
        bucket, key = self._parse_s3_path(summary_path)
        self.s3_client.put_object(
            Bucket=bucket,
            Key=key,
            Body=summary_json.encode('utf-8')
        )
        
        print(f"   ✅ Saved processing summary to: {summary_path}")
        
        return summary
    
    def _get_watershed_coordinates(self, watershed: pd.Series) -> Dict:
        """Get representative coordinates for a watershed."""
        geom = watershed.geometry
        bounds = geom.bounds
        lon = (bounds[0] + bounds[2]) / 2
        lat = (bounds[1] + bounds[3]) / 2
        
        return {
            'latitude': lat,
            'longitude': lon
        }
    
    def _save_query_metadata(self, training_folder_path: str, metadata: Dict):
        """Save metadata about queries to this watershed's training folder."""
        queries_path = f"{training_folder_path}/query_history.jsonl"
        query_record = json.dumps(metadata) + "\n"
        
        bucket, key = self._parse_s3_path(queries_path)
        
        try:
            response = self.s3_client.get_object(Bucket=bucket, Key=key)
            existing_content = response['Body'].read().decode('utf-8')
            new_content = existing_content + query_record
        except self.s3_client.exceptions.NoSuchKey:
            new_content = query_record
        except Exception as e:
            logger.warning(f"Error reading query history: {e}")
            new_content = query_record
        
        self.s3_client.put_object(
            Bucket=bucket,
            Key=key,
            Body=new_content.encode('utf-8')
        )
    
    def _extract_watershed_locations_fast(self, all_watersheds: List[pd.Series]) -> List[Dict]:
        """Extract locations with properly formatted coordinates."""
        watershed_locations = []
        
        for idx, watershed in enumerate(all_watersheds):
            geom = watershed.geometry
            bounds = geom.bounds
            lon = (bounds[0] + bounds[2]) / 2
            lat = (bounds[1] + bounds[3]) / 2
            
            lat_str = format_coordinate(lat)
            lon_str = format_coordinate(lon)
            lat_float = float(lat_str)
            lon_float = float(lon_str)
            
            watershed_locations.append({
                'watershed_id': idx,
                'hybas_id': watershed['HYBAS_ID'],
                'latitude': lat_float,
                'longitude': lon_float,
                'lat_str': lat_str,
                'lon_str': lon_str,
                'area_km2': watershed.get('SUB_AREA', 0),
                'is_initial': idx == 0
            })
            
            if idx == 0:
                print(f"   Initial watershed: HYBAS_ID={watershed['HYBAS_ID']}, loc=({lat_str}, {lon_str})")
            else:
                print(f"   Neighbor {idx}: HYBAS_ID={watershed['HYBAS_ID']}, loc=({lat_str}, {lon_str})")
        
        return watershed_locations
    
    def _load_basin_data(self, level: int) -> gpd.GeoDataFrame:
        """Load and cache BasinATLAS data with spatial index."""
        if level not in self._basin_data_cache:
            s3_path = "s3://climate-ai-data-science-datasets/arrakis-data/BasinATLAS/BasinATLAS_v10.gdb"
            layer_name = f"BasinATLAS_v10_lev{level:02d}"
            
            print(f"   Loading BasinATLAS level {level} data...")
            gdf = gpd.read_file(s3_path, layer=layer_name)
            
            print(f"   Building spatial index...")
            _ = gdf.sindex
            
            self._basin_data_cache[level] = gdf
            print(f"   ✅ Loaded and indexed {len(gdf)} subwatersheds")
        
        return self._basin_data_cache[level]
    
    def _find_nearest_subwatersheds_optimized(
        self, 
        initial_watershed: pd.Series,
        n_neighbors: int,
        gdf: gpd.GeoDataFrame
    ) -> List[pd.Series]:
        """
        Optimized nearest neighbor search (index-safe).
        Returns a list: [initial_watershed] + n_neighbors nearest neighbors.
        """
        initial_point = initial_watershed.geometry.representative_point()

        print(f"   Calculating distances using vectorized operations...")
        # distances is a Series aligned to gdf.index (labels preserved)
        rep_points = gdf.geometry.representative_point()
        distances = rep_points.distance(initial_point)  # index-aligned Series

        # Build DataFrame indexed by gdf.index to keep alignment
        distance_df = pd.DataFrame(
            {"distance": distances, "HYBAS_ID": gdf["HYBAS_ID"]},
            index=gdf.index
        )

        # Exclude the seed by HYBAS_ID
        before_filter = len(distance_df)
        distance_df = distance_df.loc[distance_df["HYBAS_ID"] != initial_watershed["HYBAS_ID"]]
        after_filter = len(distance_df)
        print(f"   Excluded {before_filter - after_filter} watershed(s) with matching HYBAS_ID")

        if after_filter < n_neighbors:
            print(f"   WARNING: Only {after_filter} watersheds available after filtering!")
            print(f"      Requested {n_neighbors} neighbors but will return {after_filter}")

        # Select n nearest by distance
        k = min(n_neighbors, len(distance_df))
        nearest = distance_df.nsmallest(k, "distance")
        nearest_indices = nearest.index  # index labels

        # Use .loc (labels), not .iloc (positions)
        neighbor_rows = gdf.loc[nearest_indices]

        # Return [seed] + neighbors
        all_watersheds = [initial_watershed] + [row for _, row in neighbor_rows.iterrows()]

        print(f"   ✓ Returning {len(all_watersheds)} watersheds: 1 initial + {len(nearest_indices)} neighbors")
        for i, (idx, row) in enumerate(neighbor_rows.head(3).iterrows(), start=1):
            print(f"      Neighbor {i}: HYBAS_ID={row['HYBAS_ID']}, distance={nearest.loc[idx, 'distance']:.4f} degrees")

        return all_watersheds

   
    def _check_s3_folder_exists(self, training_folder_path: str) -> bool:
        """Check if the training folder exists."""
        bucket, prefix = self._parse_s3_path(training_folder_path)
        
        if not prefix.endswith('/'):
            prefix = prefix + '/'
        
        try:
            response = self.s3_client.list_objects_v2(
                Bucket=bucket,
                Prefix=prefix,
                MaxKeys=1
            )
            
            exists = 'Contents' in response and len(response['Contents']) > 0
            
            if exists:
                print(f"   Found existing content in training folder")
            
            return exists
            
        except Exception as e:
            logger.debug(f"Error checking folder existence: {e}")
            return False
    
    def _save_df_to_s3(self, df: pd.DataFrame, s3_path: str):
        """Save DataFrame to S3 as CSV."""
        csv_buffer = df.to_csv(index=False)
        bucket, key = self._parse_s3_path(s3_path)
        
        self.s3_client.put_object(
            Bucket=bucket,
            Key=key,
            Body=csv_buffer.encode('utf-8')
        )
    
    def _parse_s3_path(self, s3_path: str) -> Tuple[str, str]:
        """Parse S3 path into bucket and key."""
        path = s3_path.replace('s3://', '').rstrip('/')
        parts = path.split('/', 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ''
        return bucket, key
    
    def _create_model_metadata(
        self,
        hybas_id: float,
        lat: float,
        lon: float,
        training_folder_path: str
    ):
        """Create and save model metadata with properly formatted coordinates."""
        lat_str = format_coordinate(lat)
        lon_str = format_coordinate(lon)
        
        metadata = {
            'hybas_id': float(hybas_id),
            'latitude': float(lat_str),
            'longitude': float(lon_str),
            'lat_str': lat_str,
            'lon_str': lon_str,
            'created_at': pd.Timestamp.now().isoformat(),
            'model_status': 'pending_training',
            'model_version': '1.0.0'
        }
        
        metadata_json = json.dumps(metadata, indent=2)
        path = f"{training_folder_path}/models/metadata_{int(hybas_id)}.json"
        
        bucket, key = self._parse_s3_path(path)
        
        try:
            self.s3_client.put_object(
                Bucket=bucket,
                Key=key,
                Body=metadata_json.encode('utf-8'),
                ContentType='application/json'
            )
            logger.debug(f"Saved model metadata for HYBAS_ID {hybas_id}")
        except Exception as e:
            logger.error(f"Failed to save model metadata: {e}")
            raise


# Main execution function with DSWx parameters
def define_all_watersheds_in_basin(
    latitude: float,
    longitude: float,
    n_neighbors: int = 5,
    include_dswx: bool = True,
    dswx_start_date: str = "2023-01-01",
    dswx_end_date: str = None,
    reference_raster_type: str = 'DEM',
    **kwargs
) -> Dict:
    """
    Wrapper function for basin-wide watershed processing with DSWx integration.

    Parameters:
    -----------
    latitude : float
        Latitude of the point of interest
    longitude : float
        Longitude of the point of interest
    n_neighbors : int
        Number of nearest neighbor watersheds to process
    include_dswx : bool
        Whether to process OPERA DSWx water extent time series
    dswx_start_date : str
        Start date for DSWx data query (YYYY-MM-DD)
    dswx_end_date : str
        End date for DSWx data query (None = today)
    reference_raster_type : str
        Reference raster for spatial alignment: 'DEM' (30m) or 'LULC' (10m)
        Default: 'DEM'
    **kwargs
        Additional parameters passed to define_all_watersheds_in_basin
        
    Returns:
    --------
    Dict
        Processing summary with results
    """
    pipeline = OptimizedBasinTrainingPipeline()
    return pipeline.define_all_watersheds_in_basin(
        latitude=latitude,
        longitude=longitude,
        n_neighbors=n_neighbors,
        include_dswx=include_dswx,
        dswx_start_date=dswx_start_date,
        dswx_end_date=dswx_end_date,
        reference_raster_type=reference_raster_type,
        **kwargs
    )

