"""
Optimized operational streamflow prediction pipeline.
Key improvements:
- No redundant ERA5 file saving
- Watershed caching to avoid re-delineation
- Consolidated methods and reduced redundancy
- Better memory management
"""

import os
import sys
import json
import warnings
import tempfile
import tarfile
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import shutil
from pathlib import Path

import boto3
import numpy as np
import pandas as pd
import geopandas as gpd
import torch
import xarray as xr
from shapely.geometry import Point
import s3fs

# Import your existing modules
from delineateWatersheds import delineateHydroBasins, read_s3_gpkg
from extractClimateData import ClimateDataExtractor

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


class OperationalStreamflowPredictor:
    """
    Optimized operational streamflow predictor with caching and efficiency improvements.
    """
    
    def __init__(
        self,
        metadata_output_path: str = "s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/models/",
        historical_output_path: str = "s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/historicalOutputs/",
        watershed_cache_path: str = "s3://climate-ai-data-science-datasets/arrakis-data/temp_watersheds"
    ):
        """Initialize with optimizations."""
        self.metadata_output_path = metadata_output_path.rstrip('/')
        self.historical_output_path = historical_output_path.rstrip('/')
        self.watershed_cache_path = watershed_cache_path.rstrip('/')
        
        self.s3_client = boto3.client('s3')
        self.s3_fs = s3fs.S3FileSystem()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        # Model components
        self.model = None
        self.norm_params = None
        self.static_feature_names = None
        self.dynamic_feature_names = None
        self.lstm_static_features = None
        self.training_job_name = None
        self.model_tar_path = None
        self.source_tar_path = None
        self.checkpoint_base_uri = None
        self.work_dir = None
        self.model_dir = None
        self.checkpoint_path = None
        
        # Climate extractor for efficient processing
        self.climate_extractor = ClimateDataExtractor(cache_masks=True)
    
    def load_model_from_training_job(
        self,
        training_job_name: str,
        output_base_uri: str,
        checkpoint_base_uri: str = None,
        source_tar_path: Optional[str] = None,
        **kwargs  # Ignore unused parameters for backward compatibility
    ):
        """Load model components from a specific training job."""
        self.training_job_name = training_job_name
        
        # Construct paths
        self.model_tar_path = f"{output_base_uri}/{training_job_name}/output/output.tar.gz"
        
        if source_tar_path is None:
            bucket = output_base_uri.split('/')[2]
            self.source_tar_path = f"s3://{bucket}/{training_job_name}/source/sourcedir.tar.gz"
        else:
            self.source_tar_path = source_tar_path
        
        # Infer checkpoint URI if not provided
        if checkpoint_base_uri is None:
            output_parts = output_base_uri.split('/')[-1]
            if 'training-outputs_' in output_parts:
                experiment_suffix = output_parts.replace('training-outputs_', '')
                bucket = output_base_uri.split('/')[2]
                checkpoint_base_uri = f"s3://{bucket}/checkpoints_{experiment_suffix}"
        
        self.checkpoint_base_uri = checkpoint_base_uri
        
        print("Loading model components from S3...")
        print(f"  Training job: {training_job_name}")
        print(f"  Output tar path: {self.model_tar_path}")
        if checkpoint_base_uri:
            print(f"  Checkpoint base URI: {checkpoint_base_uri}")
        print(f"  Source tar path: {self.source_tar_path}")
        
        # Create working directory
        self.work_dir = tempfile.mkdtemp(prefix='streamflow_model_')
        print(f"  Working directory: {self.work_dir}")
        
        try:
            # Download components
            print("\n1. Downloading and extracting model outputs...")
            self._download_and_extract_model()
            
            print("\n2. Downloading and extracting source code...")
            self._download_and_extract_source()
            
            if self.checkpoint_base_uri:
                print("\n3. Downloading checkpoint from checkpoint storage...")
                self._download_checkpoint_from_s3()
            
            print("\n4. Loading model components...")
            self._load_model_components()
            
            print("\n✅ Model loaded successfully")
            
        except Exception as e:
            print(f"❌ Error loading model: {e}")
            self.cleanup()
            raise
    
    def _check_watershed_cache(self, latitude: float, longitude: float) -> Optional[str]:
        """Check if watershed has already been delineated for this location."""
        # Round coordinates for cache matching
        lat_rounded = round(latitude, 4)
        lon_rounded = round(longitude, 4)
        
        # Check common patterns
        patterns = [
            f"processedWatershedAndBasins_Lon{longitude}_Lat{latitude}.gpkg",
            f"processedWatershedAndBasins_Lon{lon_rounded}_Lat{lat_rounded}.gpkg",
        ]
        
        for pattern in patterns:
            test_path = f"{self.watershed_cache_path}/{pattern}"
            try:
                bucket, key = self._parse_s3_path(test_path)
                self.s3_client.head_object(Bucket=bucket, Key=key)
                print(f"  ✓ Found cached watershed: {pattern}")
                return test_path
            except:
                continue
        
        return None
    
    def _delineate_watersheds(self, latitude: float, longitude: float, level: int) -> Optional[str]:
        """Delineate or retrieve cached watersheds."""
        # First check cache
        cached_path = self._check_watershed_cache(latitude, longitude)
        if cached_path:
            return cached_path
        
        # If not cached, delineate new
        print("  Delineating new watersheds...")
        try:
            gpkg_path = delineateHydroBasins(
                latitude=latitude,
                longitude=longitude,
                output_location=self.watershed_cache_path,
                level=level
            )
            return gpkg_path
        except Exception as e:
            print(f"  ⚠️ Watershed delineation failed: {e}")
            return None
    
    def _extract_climate_without_saving(
        self,
        gpkg_path: str,
        start_date: Optional[str],
        end_date: Optional[str]
    ) -> Dict[str, pd.DataFrame]:
        """Extract climate data WITHOUT saving separate ERA5 files."""
        # Load watersheds
        gdf = read_s3_gpkg(gpkg_path)
        
        # Open ERA5 dataset efficiently
        ds = xr.open_zarr(
            's3://climate-ai-data-eng-prod-processed/era5/v1/daily_time_chunk_optimized.zarr',
            chunks={'time': 365},
            consolidated=True
        )
        
        # Filter to needed variables
        variables = ['d2m', 'ssrd', 't2m', 't2m_max', 't2m_min', 'tp']
        ds = ds[variables]
        
        climate_data = {}
        
        # Process each watershed
        for feature_type in ['subwatershed', 'riverBasin']:
            features = gdf[gdf['feature_type'] == feature_type]
            if features.empty:
                continue
            
            polygon = features.iloc[0].geometry
            
            # Extract using the efficient extractor
            df = self.climate_extractor._extract_optimized(
                ds, polygon, gdf.crs, 'time'
            )
            
            # Clean data
            df = self.climate_extractor._clean_climate_data(df)
            
            # Filter dates if specified
            if start_date:
                df = df[df['Date'] >= pd.to_datetime(start_date)]
            if end_date:
                df = df[df['Date'] <= pd.to_datetime(end_date)]
            
            # Map to expected key names
            climate_key = 'river_basin' if feature_type == 'riverBasin' else feature_type
            climate_data[climate_key] = df
        
        return climate_data
    
    def predict_for_location(
        self,
        latitude: float,
        longitude: float,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        delineate_level: int = 12
    ) -> Dict[str, Any]:
        """Optimized prediction pipeline."""
        print(f"\n{'='*60}")
        print(f"OPERATIONAL STREAMFLOW PREDICTION")
        print(f"Location: ({latitude}, {longitude})")
        print(f"{'='*60}\n")
        
        if self.model is None:
            raise RuntimeError("Model not loaded. Call load_model_from_training_job() first.")
        
        # Step 1: Get watersheds (cached or new)
        print("\n📍 Step 1: Delineating watersheds...")
        watershed_gpkg = self._delineate_watersheds(latitude, longitude, delineate_level)
        
        if watershed_gpkg is None:
            return {'error': 'Failed to delineate watershed'}
        
        # Step 2: Extract watershed attributes for both types
        print("\n📊 Step 2: Extracting watershed attributes...")
        static_features = {}
        static_features['subwatershed'] = self._extract_watershed_attributes(watershed_gpkg, 'subwatershed')
        static_features['river_basin'] = self._extract_watershed_attributes(watershed_gpkg, 'riverBasin')
        
        # Step 3: Get climate data WITHOUT saving separate files
        print("\n🌡️ Step 3: Extracting climate data...")
        climate_data = self._extract_climate_without_saving(watershed_gpkg, start_date, end_date)
        
        # Step 4: Process predictions for each watershed
        watershed_results = {}
        
        for watershed_type in ['subwatershed', 'river_basin']:
            print(f"\n🔧 Processing predictions for {watershed_type}...")
            
            # Get the correct climate data key
            climate_key = 'river_basin' if watershed_type == 'river_basin' else watershed_type
            
            if climate_key not in climate_data:
                print(f"  ⚠️ No climate data for {watershed_type}")
                continue
            
            climate_df = climate_data[climate_key]
            
            # Process climate data
            processed_climate = self._process_climate_dataframe(climate_df, start_date, end_date)
            
            if processed_climate.empty:
                print(f"  ⚠️ No data for {watershed_type} after filtering")
                continue
            
            # Prepare sequences with correct static features
            print(f"  Preparing input sequences for {watershed_type}...")
            input_sequences = self._prepare_input_sequences(
                processed_climate, 
                static_features[watershed_type]
            )
            
            if not input_sequences:
                print(f"  ⚠️ Insufficient data for {watershed_type}")
                continue
            
            print(f"  Created {len(input_sequences)} sequences for {watershed_type}")
            
            # Run inference
            print(f"  Running model inference for {watershed_type}...")
            predictions = self._run_inference(input_sequences)
            
            # Print statistics
            print(f"\n  Prediction Statistics for {watershed_type}:")
            print(f"    Mean: {predictions['streamflow_pred'].mean():.2f} mm/day")
            print(f"    Median: {predictions['streamflow_pred'].median():.2f} mm/day")
            print(f"    Max: {predictions['streamflow_pred'].max():.2f} mm/day")
            print(f"    % zeros: {100*(predictions['streamflow_pred']==0).mean():.1f}%")
            
            # Store results
            watershed_results[watershed_type] = {
                'predictions': predictions,
                'climate_data': processed_climate,
                'static_features': static_features[watershed_type]
            }
        
        # Step 5: Save consolidated results
        print("\n💾 Step 5: Saving results...")
        results = self._save_all_results(
            watershed_results, latitude, longitude, watershed_gpkg
        )
        
        print(f"\n{'='*60}")
        print("✅ PREDICTION COMPLETE")
        print(f"{'='*60}\n")
        
        return results
    
    def _download_and_extract_model(self):
        """Download and extract the model output tar.gz file."""
        local_tar = os.path.join(self.work_dir, 'output.tar.gz')
        
        bucket, key = self._parse_s3_path(self.model_tar_path)
        print(f"  Downloading from s3://{bucket}/{key}")
        self.s3_client.download_file(bucket, key, local_tar)
        
        self.model_dir = os.path.join(self.work_dir, 'model')
        os.makedirs(self.model_dir, exist_ok=True)
        
        print("  Extracting model outputs...")
        with tarfile.open(local_tar, 'r:gz') as tar:
            tar.extractall(self.model_dir)
        
        print(f"  Extracted to {self.model_dir}")
        os.remove(local_tar)
    
    def _download_and_extract_source(self):
        """Download and extract the source code tar.gz file."""
        source_local = os.path.join(self.work_dir, 'source.tar.gz')
        
        bucket, key = self._parse_s3_path(self.source_tar_path)
        print(f"  Downloading source from s3://{bucket}/{key}")
        self.s3_client.download_file(bucket, key, source_local)
        
        source_dir = os.path.join(self.work_dir, 'source')
        os.makedirs(source_dir, exist_ok=True)
        
        print("  Extracting source code...")
        with tarfile.open(source_local, 'r:gz') as tar:
            tar.extractall(source_dir)
        
        sys.path.insert(0, source_dir)
        print(f"  Added {source_dir} to Python path")
        os.remove(source_local)
    
    def _download_checkpoint_from_s3(self):
        """Download the best checkpoint based on loss."""
        bucket, prefix = self._parse_s3_path(self.checkpoint_base_uri)
        
        print(f"  Listing checkpoints in s3://{bucket}/{prefix}")
        response = self.s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix)
        
        checkpoint_files = []
        if 'Contents' in response:
            for obj in response['Contents']:
                key = obj['Key']
                if key.endswith('.ckpt'):
                    file_name = key.split('/')[-1]
                    size_mb = obj['Size'] / (1024 * 1024)
                    
                    # Try to extract loss from filename
                    loss = float('inf')
                    try:
                        if 'loss_' in file_name:
                            loss_str = file_name.split('loss_')[-1].replace('.ckpt', '').split('_')[0]
                            loss = float(loss_str)
                    except:
                        pass
                    
                    checkpoint_files.append({
                        'key': key,
                        'name': file_name,
                        'size_mb': size_mb,
                        'loss': loss
                    })
                    print(f"    Found: {file_name} ({size_mb:.2f} MB)")
        
        if not checkpoint_files:
            raise FileNotFoundError("No checkpoints found")
        
        # Sort by loss and take best
        checkpoint_files.sort(key=lambda x: x['loss'])
        best_checkpoint = checkpoint_files[0]
        
        checkpoint_dir = os.path.join(self.work_dir, 'checkpoints')
        os.makedirs(checkpoint_dir, exist_ok=True)
        
        print(f"\n  Downloading checkpoint: {best_checkpoint['name']}")
        local_checkpoint_path = os.path.join(checkpoint_dir, best_checkpoint['name'])
        
        self.s3_client.download_file(
            Bucket=bucket,
            Key=best_checkpoint['key'],
            Filename=local_checkpoint_path
        )
        
        self.checkpoint_path = local_checkpoint_path
        print(f"  ✓ Downloaded checkpoint to: {local_checkpoint_path}")
    
    def _load_model_components(self):
        """Load model components with proper feature name inference."""
        import train
        
        # Load normalization parameters
        norm_params_candidates = [
            os.path.join(self.model_dir, 'output_files', 'normalization_params.json'),
            os.path.join(self.model_dir, 'normalization_params.json'),
        ]
        
        norm_params_path = None
        for candidate in norm_params_candidates:
            if os.path.exists(candidate):
                norm_params_path = candidate
                print(f"  Found normalization_params.json at: {norm_params_path}")
                break
        
        if not norm_params_path:
            raise FileNotFoundError("Could not find normalization_params.json")
        
        with open(norm_params_path, 'r') as f:
            self.norm_params = json.load(f)
        
        # Infer feature names from normalization parameters
        if 'static_means' in self.norm_params:
            self.static_feature_names = list(self.norm_params['static_means'].keys())
            print(f"  Inferred {len(self.static_feature_names)} static features from norm params")
        
        if 'dynamic_means' in self.norm_params:
            self.dynamic_feature_names = list(self.norm_params['dynamic_means'].keys())
            print(f"  Inferred {len(self.dynamic_feature_names)} dynamic features from norm params")
            print(f"  Dynamic features: {self.dynamic_feature_names}")
        
        # Check for checkpoint
        if not (hasattr(self, 'checkpoint_path') and self.checkpoint_path and os.path.exists(self.checkpoint_path)):
            raise FileNotFoundError("Could not find model checkpoint")
        
        ckpt_path = self.checkpoint_path
        print(f"  Using checkpoint from checkpoint storage: {ckpt_path}")
        
        # Define lstm_static_features
        lstm_static_features = ['pre_mm_syr', 'pet_mm_syr', 'tmp_dc_syr']
        
        print(f"  Using LSTM static features from training: {lstm_static_features}")
        print(f"  Static feature names: {self.static_feature_names[:5] if self.static_feature_names else None}...")
        
        # Load model
        print(f"  Loading model from checkpoint...")
        self.model = train.StreamflowLightningModule.load_from_checkpoint(
            ckpt_path,
            map_location=self.device,
            norm_params=self.norm_params,
            lstm_static_features=lstm_static_features,
            static_feature_names=self.static_feature_names
        )
        self.model.eval()
        self.model.to(self.device)
        
        self.lstm_static_features = lstm_static_features
        
        print(f"  ✓ Model loaded on {self.device}")
        print(f"  ✓ LSTM expects {len(self.dynamic_feature_names) + len(lstm_static_features)} total input features")
    
    def _process_climate_dataframe(
        self,
        climate_df: pd.DataFrame,
        start_date: Optional[str],
        end_date: Optional[str]
    ) -> pd.DataFrame:
        """Process climate dataframe with unit conversions and date filtering."""
        df = climate_df.copy()
        
        # Debug original values
        print("\n  DEBUG - Original ERA5 values:")
        for col in ['d2m', 't2m', 'tp', 'ssrd']:
            if col in df.columns:
                print(f"    {col}: mean={df[col].mean():.3f}, min={df[col].min():.3f}, max={df[col].max():.3f}")
        
        # Filter dates if not already done
        if start_date and 'Date' in df.columns:
            df = df[df['Date'] >= pd.to_datetime(start_date)]
        if end_date and 'Date' in df.columns:
            df = df[df['Date'] <= pd.to_datetime(end_date)]
        
        # Rename columns to match Caravan format
        rename_map = {
            'd2m': 'dewpoint_temperature_2m_mean',
            'ssrd': 'surface_net_solar_radiation_mean',
            't2m_max': 'temperature_2m_max',
            't2m': 'temperature_2m_mean',
            't2m_min': 'temperature_2m_min',
            'tp': 'total_precipitation_sum'
        }
        df = df.rename(columns=rename_map)
        
        # Unit conversions from ERA5 to Caravan format
        temp_cols = [
            'temperature_2m_mean', 
            'temperature_2m_max', 
            'temperature_2m_min', 
            'dewpoint_temperature_2m_mean'
        ]
        for col in temp_cols:
            if col in df.columns and df[col].mean() > 200:
                print(f"    Converting {col} from Kelvin to Celsius")
                df[col] = df[col] - 273.15
        
        if 'total_precipitation_sum' in df.columns and df['total_precipitation_sum'].mean() < 0.1:
            print(f"    Converting precipitation from meters to mm")
            df['total_precipitation_sum'] = df['total_precipitation_sum'] * 1000
        
        if 'surface_net_solar_radiation_mean' in df.columns and df['surface_net_solar_radiation_mean'].mean() > 1000000:
            print(f"    Converting solar radiation from J/m² to W/m²")
            df['surface_net_solar_radiation_mean'] = df['surface_net_solar_radiation_mean'] / 86400
        
        # Debug after conversion
        print("\n  DEBUG - After unit conversion:")
        for col in ['temperature_2m_mean', 'total_precipitation_sum', 'surface_net_solar_radiation_mean']:
            if col in df.columns:
                print(f"    {col}: mean={df[col].mean():.3f}, min={df[col].min():.3f}, max={df[col].max():.3f}")
        
        return df
    
    def _add_precipitation_phase_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add engineered precipitation features matching the training pipeline."""
        df = df.copy()
        
        precip_col = 'total_precipitation_sum'
        temp_col = 'temperature_2m_mean'
        
        # Safety check
        if precip_col not in df.columns:
            df[precip_col] = 0.0
        if temp_col not in df.columns:
            df[temp_col] = 0.0
        
        print("  Computing rolling precipitation features...")
        
        # Rolling features
        df['Precip_smoothed_5day'] = df[precip_col].rolling(window=5, min_periods=1, center=False).mean()
        df['Precip_lagged_90day'] = df[precip_col].rolling(window=90, min_periods=1, center=False).mean()
        
        # Fill NaN values
        df['Precip_smoothed_5day'] = df['Precip_smoothed_5day'].fillna(df[precip_col])
        df['Precip_lagged_90day'] = df['Precip_lagged_90day'].fillna(df[precip_col])
        df['Precip_smoothed_5day'] = df['Precip_smoothed_5day'].ffill().bfill()
        df['Precip_lagged_90day'] = df['Precip_lagged_90day'].ffill().bfill()
        
        # Snow/rain partitioning
        temp = df[temp_col].values
        precip = df[precip_col].values
        
        df['total_likely_snow_sum'] = np.where(temp < 0, precip, 0.0)
        df['total_likely_rain_sum'] = np.where(temp > 0, precip, 0.0)
        
        # Debug
        print("  DEBUG - Engineered precipitation features:")
        engineered_cols = ['Precip_smoothed_5day', 'Precip_lagged_90day', 
                          'total_likely_snow_sum', 'total_likely_rain_sum']
        for col in engineered_cols:
            if col in df.columns:
                values = df[col]
                print(f"    {col}: mean={values.mean():.3f}, min={values.min():.3f}, "
                      f"max={values.max():.3f}, %zero={100*(values==0).mean():.1f}%")
        
        return df
    
    def _prepare_input_sequences(
        self,
        climate_df: pd.DataFrame,
        static_features: pd.DataFrame
    ) -> List[Dict]:
        """Prepare input sequences for the model."""
        sequence_length = 365
        
        # Add engineered features BEFORE normalization
        climate_df = self._add_precipitation_phase_features(climate_df)
        
        # Normalize features
        normalized_climate = self._normalize_climate_features(climate_df)
        normalized_static = self._normalize_static_features(static_features)
        
        # Use the expected dynamic features
        if self.dynamic_feature_names:
            dynamic_cols = self.dynamic_feature_names
        else:
            dynamic_cols = [
                'surface_net_solar_radiation_mean',
                'temperature_2m_max',
                'temperature_2m_min', 
                'dewpoint_temperature_2m_mean',
                'Precip_smoothed_5day',
                'Precip_lagged_90day',
                'total_likely_snow_sum',
                'total_likely_rain_sum'
            ]
        
        # Ensure all expected dynamic features exist
        for feat in dynamic_cols:
            if feat not in normalized_climate.columns:
                normalized_climate[feat] = 0.0
        
        # Ensure static features are complete
        if self.static_feature_names:
            for feat in self.static_feature_names:
                if feat not in normalized_static.columns:
                    normalized_static[feat] = 0.0
            static_cols = self.static_feature_names
        else:
            static_cols = list(normalized_static.columns)
        
        # Create sequences
        sequences = []
        
        if len(normalized_climate) < sequence_length:
            print(f"  ⚠️ Warning: Only {len(normalized_climate)} days of data")
            return sequences
        
        for i in range(len(normalized_climate) - sequence_length + 1):
            dynamic_seq = normalized_climate.iloc[i:i+sequence_length][dynamic_cols].values
            static_seq = normalized_static[static_cols].values[0]
            
            sequences.append({
                'dynamic': torch.FloatTensor(dynamic_seq),
                'static': torch.FloatTensor(static_seq),
                'date': climate_df.iloc[i+sequence_length-1]['Date']
            })
        
        return sequences
    
    def _normalize_climate_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize climate features using training statistics."""
        normalized = df.copy()
        
        dynamic_means = self.norm_params.get('dynamic_means', {})
        dynamic_stds = self.norm_params.get('dynamic_stds', {})
        dynamic_mins = self.norm_params.get('dynamic_mins', {})
        dynamic_maxs = self.norm_params.get('dynamic_maxs', {})
        
        minmax_cols = self.norm_params.get('minmax_cols', [])
        std_only_cols = self.norm_params.get('std_only_cols', [])
        
        print("\n  DEBUG - Normalizing climate features:")
        
        for col in df.columns:
            if col == 'Date':
                continue
            
            if col in minmax_cols and col in dynamic_mins:
                min_val = dynamic_mins[col]
                max_val = dynamic_maxs[col]
                if max_val > min_val:
                    normalized[col] = (df[col] - min_val) / (max_val - min_val)
                    print(f"    {col}: min-max normalized [{min_val:.3f}, {max_val:.3f}]")
            elif col in std_only_cols and col in dynamic_stds:
                std_val = dynamic_stds[col]
                normalized[col] = df[col] / (std_val + 1e-8)
                print(f"    {col}: std-only normalized (std={std_val:.3f})")
            elif col in dynamic_means and col in dynamic_stds:
                mean_val = dynamic_means[col]
                std_val = dynamic_stds[col]
                normalized[col] = (df[col] - mean_val) / (std_val + 1e-8)
                print(f"    {col}: z-score normalized (mean={mean_val:.3f}, std={std_val:.3f})")
        
        return normalized
    
    def _normalize_static_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize static features using training statistics."""
        normalized = df.copy()
        
        static_means = self.norm_params.get('static_means', {})
        static_stds = self.norm_params.get('static_stds', {})
        
        for col in df.columns:
            if col in static_means and col in static_stds:
                normalized[col] = (df[col] - static_means[col]) / (static_stds[col] + 1e-8)
        
        print("\n  DEBUG - Static feature normalization:")
        for feat in self.static_feature_names[:10] if self.static_feature_names else []:
            if feat in df.columns:
                raw_val = df[feat].values[0]
                norm_val = normalized[feat].values[0]
                mean = static_means.get(feat, 0)
                std = static_stds.get(feat, 1)
                print(f"    {feat}: raw={raw_val:.3f}, norm={norm_val:.3f} (mean={mean:.3f}, std={std:.3f})")
        
        return normalized
    
    def _extract_watershed_attributes(self, gpkg_path: str, feature_type: str) -> pd.DataFrame:
        """Extract HydroATLAS attributes for a specific watershed type."""
        gdf = read_s3_gpkg(gpkg_path)
        
        # Get the appropriate feature
        if feature_type == 'subwatershed':
            features = gdf[gdf['feature_type'] == 'subwatershed']
            print(f"  Extracting attributes for single subwatershed")
        else:
            features = gdf[gdf['feature_type'] == 'riverBasin']
            print(f"  Extracting attributes for aggregated river basin")
        
        if features.empty:
            print(f"  ⚠️ No {feature_type} found")
            watershed_data = gdf.iloc[0]
        else:
            watershed_data = features.iloc[0]
        
        # Consolidated column mapping
        column_mapping = {
            'SUB_AREA': 'area', 'sub_area': 'area',
            'PET_MM_SYR': 'pet_mm_syr', 'pet_mm_syr': 'pet_mm_syr',
            'AET_MM_SYR': 'aet_mm_syr', 'aet_mm_syr': 'aet_mm_syr',
            'PRE_MM_SYR': 'pre_mm_syr', 'pre_mm_syr': 'pre_mm_syr',
            'TMP_DC_SYR': 'tmp_dc_syr', 'tmp_dc_syr': 'tmp_dc_syr',
            'TMP_DC_SMN': 'tmp_dc_smn', 'tmp_dc_smn': 'tmp_dc_smn',
            'TMP_DC_SMX': 'tmp_dc_smx', 'tmp_dc_smx': 'tmp_dc_smx',
            'SNW_PC_SYR': 'snw_pc_syr', 'snw_pc_syr': 'snw_pc_syr',
            'SNW_PC_SMX': 'snw_pc_smx', 'snw_pc_smx': 'snw_pc_smx',
            'SWC_PC_SYR': 'swc_pc_syr', 'swc_pc_syr': 'swc_pc_syr',
            'INU_PC_SLT': 'inu_pc_slt', 'inu_pc_slt': 'inu_pc_slt',
            'INU_PC_SMN': 'inu_pc_smn', 'inu_pc_smn': 'inu_pc_smn',
            'INU_PC_SMX': 'inu_pc_smx', 'inu_pc_smx': 'inu_pc_smx',
            'RUN_MM_SYR': 'run_mm_syr', 'run_mm_syr': 'run_mm_syr',
            'DIS_M3_PMN': 'dis_m3_pmn', 'dis_m3_pmn': 'dis_m3_pmn',
            'DIS_M3_PMX': 'dis_m3_pmx', 'dis_m3_pmx': 'dis_m3_pmx',
            'DIS_M3_PYR': 'dis_m3_pyr', 'dis_m3_pyr': 'dis_m3_pyr',
            'LKA_PC_SSE': 'lka_pc_sse', 'lka_pc_sse': 'lka_pc_sse',
            'LKV_MC_USU': 'lkv_mc_usu', 'lkv_mc_usu': 'lkv_mc_usu',
            'REV_MC_USU': 'rev_mc_usu', 'rev_mc_usu': 'rev_mc_usu',
            'RIV_TC_USU': 'riv_tc_usu', 'riv_tc_usu': 'riv_tc_usu',
            'DOR_PC_PVA': 'dor_pc_pva', 'dor_pc_pva': 'dor_pc_pva',
            'GWT_CM_SAV': 'gwt_cm_sav', 'gwt_cm_sav': 'gwt_cm_sav',
            'ELE_MT_SAV': 'ele_mt_sav', 'ele_mt_sav': 'ele_mt_sav',
            'ELE_MT_SMN': 'ele_mt_smn', 'ele_mt_smn': 'ele_mt_smn',
            'SGR_DK_SAV': 'sgr_dk_sav', 'sgr_dk_sav': 'sgr_dk_sav',
            'ARI_IX_SAV': 'ari_ix_sav', 'ari_ix_sav': 'ari_ix_sav',
            'GLC_CL_SMJ': 'glc_cl_smj', 'glc_cl_smj': 'glc_cl_smj',
            'FOR_PC_SSE': 'for_pc_sse', 'for_pc_sse': 'for_pc_sse',
            'PST_PC_SSE': 'pst_pc_sse', 'pst_pc_sse': 'pst_pc_sse',
            'CRP_PC_SSE': 'crp_pc_sse', 'crp_pc_sse': 'crp_pc_sse',
            'IRE_PC_SSE': 'ire_pc_sse', 'ire_pc_sse': 'ire_pc_sse',
            'GLA_PC_SSE': 'gla_pc_sse', 'gla_pc_sse': 'gla_pc_sse',
            'PRM_PC_SSE': 'prm_pc_sse', 'prm_pc_sse': 'prm_pc_sse',
            'CLY_PC_SAV': 'cly_pc_sav', 'cly_pc_sav': 'cly_pc_sav',
            'SLT_PC_SAV': 'slt_pc_sav', 'slt_pc_sav': 'slt_pc_sav',
            'SND_PC_SAV': 'snd_pc_sav', 'snd_pc_sav': 'snd_pc_sav',
            'SOC_TH_SAV': 'soc_th_sav', 'soc_th_sav': 'soc_th_sav',
            'LIT_CL_SMJ': 'lit_cl_smj', 'lit_cl_smj': 'lit_cl_smj',
            'KAR_PC_SSE': 'kar_pc_sse', 'kar_pc_sse': 'kar_pc_sse',
            'URB_PC_SSE': 'urb_pc_sse', 'urb_pc_sse': 'urb_pc_sse',
            'HDI_IX_SAV': 'hdi_ix_sav', 'hdi_ix_sav': 'hdi_ix_sav',
            'NLI_IX_SAV': 'nli_ix_sav', 'nli_ix_sav': 'nli_ix_sav',
            'PPD_PK_SAV': 'ppd_pk_sav', 'ppd_pk_sav': 'ppd_pk_sav',
            'CLS_CL_SMJ': 'cls_cl_smj', 'cls_cl_smj': 'cls_cl_smj',
            'CLZ_CL_SMJ': 'clz_cl_smj', 'clz_cl_smj': 'clz_cl_smj'
        }
        
        # Extract features
        static_features = {}
        for source_col, model_col in column_mapping.items():
            if source_col in watershed_data:
                value = watershed_data[source_col]
                if pd.notna(value) and value != 0:
                    static_features[model_col] = value
        
        # Debug output
        print(f"\n  DEBUG - {feature_type} static features extracted:")
        for feat in ['area', 'pre_mm_syr', 'pet_mm_syr', 'tmp_dc_syr']:
            if feat in static_features:
                print(f"    {feat}: {static_features[feat]:.2f}")
        
        # Ensure all expected features exist
        if self.static_feature_names:
            for feature in self.static_feature_names:
                if feature not in static_features:
                    static_features[feature] = 0.0
        
        return pd.DataFrame([static_features])
    
    def _run_inference(self, sequences: List[Dict]) -> pd.DataFrame:
        """Run model inference on prepared sequences."""
        import train
        
        predictions = []
        
        # Debug first few sequences
        if sequences:
            first_seq = sequences[0]
            print(f"\n  DEBUG - First sequence shape: dynamic {first_seq['dynamic'].shape}, static {first_seq['static'].shape}")
            
            for i, feat_name in enumerate(self.dynamic_feature_names[:5]):
                feat_values = first_seq['dynamic'][:, i].numpy()
                print(f"    {feat_name}: min={feat_values.min():.3f}, max={feat_values.max():.3f}, mean={feat_values.mean():.3f}")
        
        with torch.no_grad():
            for i, seq_data in enumerate(sequences):
                if i < 5:  # Debug first 5
                    dynamic_input = seq_data['dynamic'].unsqueeze(0).to(self.device)
                    static_input = seq_data['static'].unsqueeze(0).to(self.device)
                    
                    if hasattr(self.model, 'model'):
                        pred = self.model.model(dynamic_input, static_input, self.static_feature_names)
                    else:
                        pred = self.model(dynamic_input, static_input)
                    
                    pred_norm = pred.cpu().numpy()[0, 0]
                    
                    # Debug output
                    print(f"\n  DEBUG - Sequence {i+1}:")
                    print(f"    Raw model output (normalized): {pred_norm:.6f}")
                    
                    target_mean = self.norm_params.get('target_mean', 0.0)
                    target_std = self.norm_params.get('target_std', 1.0)
                    
                    pred_original = train.unnormalize_predictions(
                        np.array([pred_norm]), 
                        self.norm_params
                    )[0]
                    
                    print(f"    After denormalization: {pred_original:.6f}")
                    print(f"    Final: {max(0, pred_original):.6f}")
                    
                    predictions.append({
                        'date': seq_data['date'],
                        'streamflow_pred': max(0, pred_original)
                    })
                else:
                    # Process rest without debug
                    if i % 100 == 0:
                        print(f"  Processing sequence {i+1}/{len(sequences)}...")
                    
                    dynamic_input = seq_data['dynamic'].unsqueeze(0).to(self.device)
                    static_input = seq_data['static'].unsqueeze(0).to(self.device)
                    
                    if hasattr(self.model, 'model'):
                        pred = self.model.model(dynamic_input, static_input, self.static_feature_names)
                    else:
                        pred = self.model(dynamic_input, static_input)
                    
                    pred_norm = pred.cpu().numpy()[0, 0]
                    pred_original = train.unnormalize_predictions(
                        np.array([pred_norm]), 
                        self.norm_params
                    )[0]
                    
                    predictions.append({
                        'date': seq_data['date'],
                        'streamflow_pred': max(0, pred_original)
                    })
        
        return pd.DataFrame(predictions)
    
    def _save_all_results(self, watershed_results: Dict, latitude: float, longitude: float, watershed_gpkg: str) -> Dict[str, Any]:
        """Save predictions and metadata for all processed watersheds."""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        historical_path = self.historical_output_path.rstrip('/')
        metadata_path = self.metadata_output_path.rstrip('/')
        
        all_results = {}
        
        for watershed_type, data in watershed_results.items():
            location_tag = f"Lat{latitude}_Lon{longitude}_{watershed_type}"
            
            # Merge predictions with climate data
            results_df = pd.merge(
                data['predictions'],
                data['climate_data'],
                left_on='date',
                right_on='Date',
                how='left'
            )
            
            # Save predictions (with climate data included)
            predictions_filename = f"streamflow_predictions_{location_tag}_{timestamp}.csv"
            predictions_path = f"{historical_path}/{predictions_filename}"
            
            self._save_df_to_s3(results_df, predictions_path)
            
            # Calculate statistics
            stats = {
                'mean_flow': float(data['predictions']['streamflow_pred'].mean()),
                'median_flow': float(data['predictions']['streamflow_pred'].median()),
                'max_flow': float(data['predictions']['streamflow_pred'].max()),
                'min_flow': float(data['predictions']['streamflow_pred'].min()),
                'std_flow': float(data['predictions']['streamflow_pred'].std())
            }
            
            # Get area
            static_df = data.get('static_features')
            area = float(static_df['area'].values[0]) if static_df is not None and 'area' in static_df.columns else None
            
            all_results[watershed_type] = {
                'predictions': results_df,
                'output_path': predictions_path,
                'statistics': stats,
                'area_km2': area
            }
            
            print(f"  - Saved {watershed_type}: {predictions_path}")
        
        # Save metadata
        metadata = {
            'timestamp': timestamp,
            'location': {'latitude': latitude, 'longitude': longitude},
            'watershed': {'gpkg_path': watershed_gpkg},
            'model_info': {
                'training_job': self.training_job_name,
                'model_tar_path': self.model_tar_path,
                'source_tar_path': self.source_tar_path
            },
            'predictions': {}
        }
        
        for watershed_type, data in all_results.items():
            metadata['predictions'][watershed_type] = {
                'output_path': data['output_path'],
                'n_predictions': len(data['predictions']),
                'area_km2': data['area_km2'],
                'date_range': {
                    'start': data['predictions']['date'].min().strftime('%Y-%m-%d'),
                    'end': data['predictions']['date'].max().strftime('%Y-%m-%d')
                },
                'statistics': data['statistics']
            }
        
        metadata_filename = f"metadata_Lat{latitude}_Lon{longitude}_{timestamp}.json"
        metadata_path_full = f"{metadata_path}/{metadata_filename}"
        
        self._save_json_to_s3(metadata, metadata_path_full)
        
        print(f"\n📁 Results saved:")
        print(f"  - Metadata: {metadata_path_full}")
        
        return {
            'watersheds': all_results,
            'metadata': metadata,
            'output_paths': {
                'metadata': metadata_path_full,
                'watershed': watershed_gpkg
            }
        }
    
    def cleanup(self):
        """Clean up temporary files."""
        if self.work_dir and os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
            print(f"Cleaned up working directory: {self.work_dir}")
            self.work_dir = None
    
    def _parse_s3_path(self, s3_path: str) -> Tuple[str, str]:
        """Parse S3 path into bucket and key."""
        s3_path = s3_path.replace('s3://', '')
        parts = s3_path.split('/', 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ''
        return bucket, key
    
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
            Body=json_str.encode('utf-8')
        )


def run_streamflow_prediction(
    training_job_name: str,
    output_base_uri: str,
    latitude: float,
    longitude: float,
    start_date: Optional[str] = None,
    end_date: Optional[str] = None,
    source_tar_path: Optional[str] = None,
    checkpoint_base_uri: Optional[str] = None,
    **kwargs  # Accept and ignore extra parameters for backward compatibility
) -> Dict[str, Any]:
    """Convenience function to run streamflow prediction for a location."""
    predictor = OperationalStreamflowPredictor()
    
    try:
        predictor.load_model_from_training_job(
            training_job_name=training_job_name,
            output_base_uri=output_base_uri,
            checkpoint_base_uri=checkpoint_base_uri,
            source_tar_path=source_tar_path
        )
        
        results = predictor.predict_for_location(
            latitude=latitude,
            longitude=longitude,
            start_date=start_date,
            end_date=end_date
        )
        
        return results
        
    finally:
        predictor.cleanup()
