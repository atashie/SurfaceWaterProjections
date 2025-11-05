# processAllStreamflowDataForTraining.py
"""
Wrapper for running streamflow predictions using existing watershed boundaries and climate data.
Reuses data from define_all_watersheds_in_basin() instead of re-delineating.
"""

import os
import json
import logging
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import pandas as pd
import numpy as np
import boto3
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings
import geopandas as gpd

from runStreamflowPrediction import OperationalStreamflowPredictor
from coordinate_utils import format_coordinate
from delineateWatersheds import read_s3_gpkg

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BasinStreamflowPredictor:
    """
    Run streamflow predictions for all watersheds in a basin training folder.
    REUSES existing watershed boundaries and climate data instead of re-delineating.
    """
    
    def __init__(self):
        self.s3_client = boto3.client('s3')
        
    def run_predictions_for_basin(
        self,
        training_folder_path: str,
        training_job_name: str,
        output_base_uri: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        source_tar_path: Optional[str] = None,
        checkpoint_base_uri: Optional[str] = None,
        parallel_processing: bool = False,
        max_workers: int = 2
    ) -> Dict:
        """
        Run streamflow predictions for all watersheds in a training folder.
        REUSES existing watershed boundaries and climate data from the training folder.
        
        Parameters:
        -----------
        training_folder_path : str
            S3 path to training folder (e.g., s3://bucket/trainingFolderFor_Lat_35.9167_Lon_-78.9937)
        training_job_name : str
            Name of the SageMaker training job
        output_base_uri : str
            Base S3 URI for model outputs
        start_date : str, optional
            Start date for predictions (YYYY-MM-DD)
        end_date : str, optional
            End date for predictions (YYYY-MM-DD)
        source_tar_path : str, optional
            Path to source code tar file
        checkpoint_base_uri : str, optional
            Base URI for model checkpoints
        parallel_processing : bool
            Whether to process watersheds in parallel
        max_workers : int
            Maximum number of parallel workers
            
        Returns:
        --------
        Dict
            Summary of prediction results
        """
        
        print(f"\n{'='*60}")
        print(f"BASIN-WIDE STREAMFLOW PREDICTIONS")
        print(f"Training folder: {training_folder_path}")
        print(f"Model: {training_job_name}")
        print(f"{'='*60}\n")
        
        # Step 1: Load watershed locations
        print("📍 Step 1: Loading watershed locations...")
        watershed_locations = self._load_watershed_locations(training_folder_path)
        
        if not watershed_locations:
            logger.error("No watershed locations found")
            return {'error': 'No watershed locations found'}
        
        print(f"   Found {len(watershed_locations)} watersheds to process")
        
        # Step 2: Load model once for all watersheds
        print("\n🔧 Step 2: Loading model...")
        predictor = OperationalStreamflowPredictor()
        
        try:
            predictor.load_model_from_training_job(
                training_job_name=training_job_name,
                output_base_uri=output_base_uri,
                checkpoint_base_uri=checkpoint_base_uri,
                source_tar_path=source_tar_path
            )
            print("   ✅ Model loaded successfully")
            
            # Step 3: Create output folder for predictions
            predictions_folder = f"{training_folder_path}/streamflowPredictions"
            print(f"\n📁 Step 3: Predictions will be saved to: {predictions_folder}")
            
            # Step 4: Process each watershed using existing data
            print(f"\n🌊 Step 4: Running streamflow predictions using existing data...")
            
            if parallel_processing and len(watershed_locations) > 1:
                results = self._process_watersheds_parallel(
                    watershed_locations,
                    training_folder_path,
                    predictions_folder,
                    predictor,
                    start_date,
                    end_date,
                    max_workers
                )
            else:
                results = self._process_watersheds_sequential(
                    watershed_locations,
                    training_folder_path,
                    predictions_folder,
                    predictor,
                    start_date,
                    end_date
                )
            
            # Step 5: Create summary
            print(f"\n📊 Step 5: Creating prediction summary...")
            summary = self._create_prediction_summary(
                results, 
                training_folder_path, 
                predictions_folder,
                training_job_name
            )
            
            print(f"\n{'='*60}")
            print(f"✅ STREAMFLOW PREDICTIONS COMPLETE")
            print(f"{'='*60}")
            print(f"\nSummary:")
            print(f"  - Total watersheds: {summary['total_watersheds']}")
            print(f"  - Successful predictions: {summary['successful_predictions']}")
            print(f"  - Failed predictions: {summary['failed_predictions']}")
            print(f"  - Output folder: {predictions_folder}")
            
            return summary
            
        finally:
            # Cleanup
            predictor.cleanup()
    
    def _process_single_watershed_optimized(
        self,
        location: Dict,
        training_folder_path: str,
        predictions_folder: str,
        predictor: OperationalStreamflowPredictor,
        start_date: Optional[str],
        end_date: Optional[str]
    ) -> Dict:
        """
        Process streamflow predictions using EXISTING watershed boundaries and climate data.
        """
        
        hybas_id = location['hybas_id']
        lat = location['latitude']
        lon = location['longitude']
        lat_str = format_coordinate(lat)
        lon_str = format_coordinate(lon)
        
        result = {
            'hybas_id': hybas_id,
            'location': location,
            'status': 'processing',
            'predictions': {},
            'errors': []
        }
        
        print(f"\n   Processing watershed {hybas_id} at ({lat}, {lon})...")
        
        try:
            # Step 1: Load EXISTING watershed boundary file
            watershed_file = f"{training_folder_path}/watershedBoundaries/processedWatershedAndBasins_Lon{lon_str}_Lat{lat_str}.gpkg"
            
            # Check if watershed file exists
            if not self._check_s3_file_exists(watershed_file):
                print(f"      ⚠️ Watershed boundary not found: {watershed_file}")
                result['status'] = 'no_watershed'
                return result
            
            print(f"      ✓ Using existing watershed boundary")
            # Compute and cache areas (km^2) for both subwatershed and river basin
            areas_km2 = self._get_areas_from_gpkg(watershed_file)
            
            # Step 2: Load EXISTING climate data
            print(f"      Loading existing climate data...")
            climate_data = {}
            
            climate_folder = f"{training_folder_path}/climateData"
            
            # Load subwatershed climate
            subwatershed_climate = f"{climate_folder}/era5_subwatershed_Lat{lat_str}_Lon{lon_str}.csv"
            if self._check_s3_file_exists(subwatershed_climate):
                climate_data['subwatershed'] = self._load_csv_from_s3(subwatershed_climate)
                print(f"      ✓ Loaded subwatershed climate data ({len(climate_data['subwatershed'])} days)")
            
            # Load river basin climate
            basin_climate = f"{climate_folder}/era5_river_basin_Lat{lat_str}_Lon{lon_str}.csv"
            if self._check_s3_file_exists(basin_climate):
                climate_data['river_basin'] = self._load_csv_from_s3(basin_climate)
                print(f"      ✓ Loaded river basin climate data ({len(climate_data['river_basin'])} days)")
            
            if not climate_data:
                print(f"      ⚠️ No climate data found for watershed {hybas_id}")
                result['status'] = 'no_climate_data'
                return result
            
            # Step 3: Extract static features from existing watershed
            print(f"      Extracting watershed attributes...")
            static_features = {}
            
            # Extract attributes for both watershed types
            static_features['subwatershed'] = predictor._extract_watershed_attributes(
                watershed_file, 'subwatershed'
            )
            static_features['river_basin'] = predictor._extract_watershed_attributes(
                watershed_file, 'riverBasin'
            )
            
            # Step 4: Process predictions for each watershed type
            print(f"      Running predictions...")
            watershed_results = {}
            
            for watershed_type in ['subwatershed', 'river_basin']:
                if watershed_type not in climate_data:
                    continue
                
                print(f"        Processing {watershed_type}...")
                
                climate_df = climate_data[watershed_type]
                
                # Store original climate data for merging later
                original_climate_df = climate_df.copy()
                
                # Ensure Date column is datetime
                if 'Date' in climate_df.columns:
                    climate_df['Date'] = pd.to_datetime(climate_df['Date'])
                    original_climate_df['Date'] = pd.to_datetime(original_climate_df['Date'])
                
                # Filter dates if specified
                if start_date:
                    climate_df = climate_df[climate_df['Date'] >= pd.to_datetime(start_date)]
                    original_climate_df = original_climate_df[original_climate_df['Date'] >= pd.to_datetime(start_date)]
                if end_date:
                    climate_df = climate_df[climate_df['Date'] <= pd.to_datetime(end_date)]
                    original_climate_df = original_climate_df[original_climate_df['Date'] <= pd.to_datetime(end_date)]
                
                if climate_df.empty:
                    print(f"        ⚠️ No data in date range for {watershed_type}")
                    continue
                
                # Process climate data (unit conversions, etc.)
                processed_climate = predictor._process_climate_dataframe(
                    climate_df, start_date, end_date
                )
                
                if processed_climate.empty:
                    print(f"        ⚠️ No data after processing for {watershed_type}")
                    continue
                
                # Prepare input sequences
                input_sequences = predictor._prepare_input_sequences(
                    processed_climate, 
                    static_features[watershed_type]
                )
                
                if not input_sequences:
                    print(f"        ⚠️ Insufficient data for {watershed_type} (need 365+ days)")
                    continue
                
                print(f"        Created {len(input_sequences)} sequences")
                
                # Run inference
                predictions = predictor._run_inference(input_sequences)

                # Convert LSTM output from mm to km^3 using area (km^2)
                area_km2 = areas_km2.get(watershed_type)
                if area_km2 is None or not np.isfinite(area_km2):
                    print(f"        ⚠️ Missing SUB_AREA for {watershed_type}; leaving units as mm/day")
                else:
                    # 1 mm = 1e-6 km; volume (km^3) = depth_km * area_km2
                    # => streamflow_km3 = streamflow_mm * 1e-6 * area_km2
                    predictions["streamflow_pred"] = (
                        predictions["streamflow_pred"].astype(np.float32) * 1e-6 * float(area_km2)
                    )
                
                # Recompute statistics after potential unit conversion
                stats = {
                    'mean': float(predictions['streamflow_pred'].mean()),
                    'median': float(predictions['streamflow_pred'].median()),
                    'max': float(predictions['streamflow_pred'].max()),
                    'min': float(predictions['streamflow_pred'].min()),
                    'std': float(predictions['streamflow_pred'].std()),
                    'zero_flow_days': int((predictions['streamflow_pred'] == 0).sum()),
                    'n_predictions': len(predictions)
                }
                
                # Save predictions (unchanged)
                output_path = f"{predictions_folder}/{watershed_type}_predictions_{int(hybas_id)}.csv"
                self._save_df_to_s3(predictions, output_path)
                
                watershed_results[watershed_type] = {
                    'predictions': predictions,
                    'output_path': output_path,
                    'statistics': stats
                }
                
                print(f"        ✅ Saved {watershed_type} predictions ({len(predictions)} days)")
                # Update units message based on conversion
                if area_km2 is None or not np.isfinite(area_km2):
                    print(f"           Mean flow: {stats['mean']:.6f} mm/day (no SUB_AREA found)")
                else:
                    print(f"           Mean flow: {stats['mean']:.8f} km^3/day (area={area_km2:.2f} km^2)")
                
                # Fix column name mismatch: predictions has 'date', climate has 'Date'
                if 'date' in predictions.columns:
                    predictions = predictions.rename(columns={'date': 'Date'})
                
                # Align predictions with original climate data by Date
                if 'Date' in predictions.columns and 'Date' in original_climate_df.columns:
                    # Merge predictions with ALL original climate columns
                    merged_df = pd.merge(
                        original_climate_df,
                        predictions[['Date', 'streamflow_pred']],
                        on='Date',
                        how='inner'
                    )
                    
                    # Reorder columns to have Date first, then weather features, then streamflow prediction
                    date_col = ['Date']
                    weather_cols = [col for col in original_climate_df.columns if col != 'Date']
                    streamflow_cols = ['streamflow_pred']
                    
                    merged_df = merged_df[date_col + weather_cols + streamflow_cols]
                    
                    print(f"        ✅ Merged {len(weather_cols)} weather features with predictions")
                else:
                    # Fallback if Date column is missing
                    merged_df = predictions
                    print(f"        ⚠️ Could not merge weather data (Date column missing)")
                
                # Calculate statistics on streamflow predictions only
                stats = {
                    'mean': float(merged_df['streamflow_pred'].mean()),
                    'median': float(merged_df['streamflow_pred'].median()),
                    'max': float(merged_df['streamflow_pred'].max()),
                    'min': float(merged_df['streamflow_pred'].min()),
                    'std': float(merged_df['streamflow_pred'].std()),
                    'zero_flow_days': int((merged_df['streamflow_pred'] == 0).sum()),
                    'n_predictions': len(merged_df)
                }
                
                # CHANGE NAMING CONVENTION: Use Lat_Lon format instead of hybas_id
                output_path = f"{predictions_folder}/{watershed_type}_predictions_Lat{lat_str}_Lon{lon_str}.csv"
                self._save_df_to_s3(merged_df, output_path)
                
                watershed_results[watershed_type] = {
                    'predictions': merged_df,  # Now includes weather data
                    'output_path': output_path,
                    'statistics': stats,
                    'n_weather_features': len(weather_cols) if 'weather_cols' in locals() else 0
                }
                
                print(f"        ✅ Saved {watershed_type} predictions with weather data ({len(merged_df)} days)")
                print(f"           Mean flow: {stats['mean']:.8f} km3/day")
                print(f"           Weather features included: {watershed_results[watershed_type]['n_weather_features']}")
            
            if watershed_results:
                result['predictions'] = watershed_results
                result['status'] = 'success'
                print(f"      ✅ Successfully processed watershed {hybas_id}")
            else:
                result['status'] = 'no_valid_predictions'
                print(f"      ⚠️ No valid predictions for watershed {hybas_id}")
            
        except Exception as e:
            result['status'] = 'failed'
            result['errors'].append(str(e))
            logger.error(f"Failed to process watershed {hybas_id}: {e}")
            print(f"      ❌ Failed to process watershed {hybas_id}: {e}")
        
        return result    
    def _process_watersheds_sequential(
        self,
        watershed_locations: List[Dict],
        training_folder_path: str,
        predictions_folder: str,
        predictor: OperationalStreamflowPredictor,
        start_date: Optional[str],
        end_date: Optional[str]
    ) -> List[Dict]:
        """Process watersheds sequentially using existing data."""
        results = []
        
        for location in tqdm(watershed_locations, desc="Processing watersheds"):
            result = self._process_single_watershed_optimized(
                location,
                training_folder_path,
                predictions_folder,
                predictor,
                start_date,
                end_date
            )
            results.append(result)
        
        return results
    
    def _process_watersheds_parallel(
        self,
        watershed_locations: List[Dict],
        training_folder_path: str,
        predictions_folder: str,
        predictor: OperationalStreamflowPredictor,
        start_date: Optional[str],
        end_date: Optional[str],
        max_workers: int
    ) -> List[Dict]:
        """Process watersheds in parallel using existing data."""
        results = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_single_watershed_optimized,
                    location,
                    training_folder_path,
                    predictions_folder,
                    predictor,
                    start_date,
                    end_date
                ): location
                for location in watershed_locations
            }
            
            for future in tqdm(as_completed(futures), total=len(futures),
                             desc="Processing watersheds"):
                try:
                    result = future.result(timeout=600)  # 10 minute timeout
                    results.append(result)
                except Exception as e:
                    location = futures[future]
                    logger.error(f"Failed processing watershed {location['hybas_id']}: {e}")
                    results.append({
                        'hybas_id': location['hybas_id'],
                        'location': location,
                        'status': 'failed',
                        'errors': [str(e)]
                    })
        
        return results
    
    def _load_watershed_locations(self, training_folder_path: str) -> List[Dict]:
        """Load watershed locations from the training folder."""
        locations_path = f"{training_folder_path}/list_of_training_locations.csv"
        
        try:
            df = self._load_csv_from_s3(locations_path)
            locations = df.to_dict('records')
            return locations
            
        except Exception as e:
            logger.error(f"Failed to load watershed locations: {e}")
            return []
    
    def _load_csv_from_s3(self, s3_path: str) -> pd.DataFrame:
        """Load CSV from S3."""
        bucket, key = self._parse_s3_path(s3_path)
        response = self.s3_client.get_object(Bucket=bucket, Key=key)
        return pd.read_csv(response['Body'])
    
    def _check_s3_file_exists(self, s3_path: str) -> bool:
        """Check if an S3 file exists."""
        bucket, key = self._parse_s3_path(s3_path)
        
        try:
            self.s3_client.head_object(Bucket=bucket, Key=key)
            return True
        except:
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
    
    def _save_json_to_s3(self, data: Dict, s3_path: str):
        """Save dictionary to S3 as JSON."""
        json_str = json.dumps(data, indent=2, default=str)
        bucket, key = self._parse_s3_path(s3_path)
        
        self.s3_client.put_object(
            Bucket=bucket,
            Key=key,
            Body=json_str.encode('utf-8'),
            ContentType='application/json'
        )
    
    def _parse_s3_path(self, s3_path: str) -> Tuple[str, str]:
        """Parse S3 path into bucket and key."""
        path = s3_path.replace('s3://', '').rstrip('/')
        parts = path.split('/', 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ''
        return bucket, key

    def _get_areas_from_gpkg(self, watershed_file: str) -> Dict[str, Optional[float]]:
        """
        Read the GeoPackage and return areas (km^2) for subwatershed and river basin.
        - SUB_AREA is expected to be in km^2 (HydroBASINS / HydroATLAS).
        - For river basin, if SUB_AREA is missing, falls back to summing subwatershed SUB_AREA.
        """
        areas = {"subwatershed": None, "river_basin": None}
        try:
            gdf = read_s3_gpkg(watershed_file)
        except Exception as e:
            print(f"      ⚠️ Failed to read GPkg for areas: {e}")
            return areas
    
        # Subwatershed area
        try:
            if "feature_type" in gdf.columns and "SUB_AREA" in gdf.columns:
                sw = gdf[gdf["feature_type"] == "subwatershed"]
                if not sw.empty and pd.notna(sw.iloc[0]["SUB_AREA"]):
                    areas["subwatershed"] = float(sw.iloc[0]["SUB_AREA"])
        except Exception:
            pass
    
        # River basin area
        try:
            if "feature_type" in gdf.columns and "SUB_AREA" in gdf.columns:
                rb = gdf[gdf["feature_type"] == "riverBasin"]
                if not rb.empty and pd.notna(rb.iloc[0]["SUB_AREA"]):
                    areas["river_basin"] = float(rb.iloc[0]["SUB_AREA"])
                else:
                    # Fallback: sum all subwatershed SUB_AREA values
                    sw = gdf[gdf["feature_type"] == "subwatershed"]
                    if not sw.empty and "SUB_AREA" in sw.columns:
                        areas["river_basin"] = float(sw["SUB_AREA"].sum())
        except Exception:
            pass
    
        return areas
    
    def _create_prediction_summary(
        self, 
        results: List[Dict],
        training_folder_path: str,
        predictions_folder: str,
        training_job_name: str
    ) -> Dict:
        """Create and save a summary of all predictions."""
        
        successful = sum(1 for r in results if r['status'] == 'success')
        failed = sum(1 for r in results if r['status'] == 'failed')
        no_data = sum(1 for r in results if r['status'] in ['no_watershed', 'no_climate_data', 'no_valid_predictions'])
        
        # Aggregate statistics across all successful predictions
        all_predictions = []
        for result in results:
            if result['status'] == 'success' and 'predictions' in result:
                for watershed_type, pred_data in result['predictions'].items():
                    all_predictions.extend(pred_data['predictions']['streamflow_pred'].values)
        
        if all_predictions:
            overall_stats = {
                'mean_flow': float(np.mean(all_predictions)),
                'median_flow': float(np.median(all_predictions)),
                'max_flow': float(np.max(all_predictions)),
                'min_flow': float(np.min(all_predictions)),
                'std_flow': float(np.std(all_predictions)),
                'total_predictions': len(all_predictions)
            }
        else:
            overall_stats = None

        if overall_stats:
            print(f"   ℹ️ Overall mean flow (km^3/day): {overall_stats['mean_flow']:.8f}")
            
        summary = {
            'total_watersheds': len(results),
            'successful_predictions': successful,
            'failed_predictions': failed,
            'no_data_watersheds': no_data,
            'training_folder': training_folder_path,
            'predictions_folder': predictions_folder,
            'training_job_name': training_job_name,
            'created_at': datetime.now().isoformat(),
            'overall_statistics': overall_stats,
            'watershed_results': []
        }
        
        # Add individual watershed summaries
        for result in results:
            watershed_summary = {
                'hybas_id': result['hybas_id'],
                'location': result['location'],
                'status': result['status'],
                'latitude': result['location']['latitude'],
                'longitude': result['location']['longitude'],
                'lat_str': format_coordinate(result['location']['latitude']),
                'lon_str': format_coordinate(result['location']['longitude'])
            }
            
            if result['status'] == 'success' and 'predictions' in result:
                watershed_summary['predictions'] = {}
                for watershed_type, pred_data in result['predictions'].items():
                    watershed_summary['predictions'][watershed_type] = {
                        'output_path': pred_data['output_path'],
                        'statistics': pred_data['statistics'],
                        'n_weather_features': pred_data.get('n_weather_features', 0)
                    }
            elif 'errors' in result:
                watershed_summary['errors'] = result['errors']
            
            summary['watershed_results'].append(watershed_summary)
        
        # Save summary
        summary_path = f"{predictions_folder}/prediction_summary.json"
        self._save_json_to_s3(summary, summary_path)
        print(f"   ✅ Saved prediction summary to: {summary_path}")
        
        return summary


# Convenience function
def run_basin_streamflow_predictions(
    training_folder_path: str,
    training_job_name: str,
    output_base_uri: str,
    **kwargs
) -> Dict:
    """
    Convenience function to run streamflow predictions for all watersheds in a basin.
    
    Example usage:
    --------------
    results = run_basin_streamflow_predictions(
        training_folder_path="s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/trainingFolderFor_Lat_35.9167_Lon_-78.9937",
        training_job_name="pytorch-training-2025-09-05-16-19-21-357",
        output_base_uri="s3://sagemaker-us-east-2-672980278503/training-outputs_lightning_164_pcaShuffle_yt2lyr5dialNoNoiseComplxNorm3xLossVal",
        start_date="2010-01-01",
        end_date="2025-09-17",
        source_tar_path="s3://sagemaker-us-east-2-672980278503/pytorch-training-2025-09-05-16-19-21-357/source/sourcedir.tar.gz",
        checkpoint_base_uri="s3://sagemaker-us-east-2-672980278503/checkpoints_lightning_164_pcaShuffle_yt2lyr5dialNoNoiseComplxNorm3xLossVal",
        parallel_processing=False
    )
    """
    predictor = BasinStreamflowPredictor()
    return predictor.run_predictions_for_basin(
        training_folder_path=training_folder_path,
        training_job_name=training_job_name,
        output_base_uri=output_base_uri,
        **kwargs
    )
