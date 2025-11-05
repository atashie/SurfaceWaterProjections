# flood_data_container.py

import numpy as np
import pandas as pd
import torch
import rasterio
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from datetime import datetime
import logging
from dataclasses import dataclass
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class WatershedMetadata:
    """Metadata for a single watershed."""
    hybas_id: float
    latitude: float
    longitude: float
    lat_str: str
    lon_str: str
    is_initial: bool
    area_km2: float


class FloodDataContainer:
    """
    Unified container for all flood modeling data components.
    Supports multiple output formats for different ML algorithms.
    """
    
    def __init__(self, watershed_metadata: WatershedMetadata, base_path: str):
        """
        Initialize container for a single watershed.
        
        Parameters:
        -----------
        watershed_metadata : WatershedMetadata
            Metadata for the watershed
        base_path : str
            Base S3 or local path to training folder
        """
        self.metadata = watershed_metadata
        self.base_path = Path(base_path)
        
        # Spatial dimensions (will be set when loading raster)
        self.height: Optional[int] = None
        self.width: Optional[int] = None
        self.n_pixels: Optional[int] = None
        self.crs: Optional[str] = None
        self.transform: Optional[object] = None
        
        # Static features
        self.static_raster: Optional[np.ndarray] = None  # (H, W, C_static)
        self.static_bands: List[str] = []
        self.static_tabular: Optional[np.ndarray] = None  # (N_pixels, C_static)
        
        # Dynamic features
        self.dynamic_subwatershed: Optional[pd.DataFrame] = None
        self.dynamic_basin: Optional[pd.DataFrame] = None
        
        # Targets
        self.flood_masks_spatial: Optional[np.ndarray] = None  # (T, H, W)
        self.flood_masks_tabular: Optional[np.ndarray] = None  # (T, N_pixels)
        self.flood_dates: List[datetime] = []
        
        # Temporal alignment
        self.timestamps: List[datetime] = []
        self.valid_timesteps: Optional[np.ndarray] = None
        
        # Spatial metadata
        self.pixel_coords: Optional[np.ndarray] = None  # (N_pixels, 2) lat/lon
        self.pixel_valid_mask: Optional[np.ndarray] = None  # (N_pixels,)
        
        # Placeholder for graph structure (GNN)
        self.adjacency_matrix: Optional[torch.sparse.Tensor] = None
        self.edge_features: Optional[torch.Tensor] = None
        
        logger.info(f"Initialized FloodDataContainer for watershed {watershed_metadata.hybas_id}")
    
    def load_static_raster(self, raster_path: Optional[str] = None):
        """
        Load static raster features from multiband TIFF.
        Properly handles per-band NoData values as saved by extractRasterData.py
        """
        if raster_path is None:
            raster_path = (self.base_path / 'rasters' / 
                          f'multiband_subwatershed_Lon{self.metadata.lon_str}_'
                          f'Lat{self.metadata.lat_str}.tif')
        
        raster_path = str(raster_path)
        logger.info(f"Loading static raster from {raster_path}")
        
        # Handle S3 paths
        if raster_path.startswith('s3://'):
            vsi_path = '/vsis3/' + raster_path.replace('s3://', '')
            vsi_path = vsi_path.replace('//', '/')
            
            from rasterio.session import AWSSession
            with rasterio.Env(AWSSession(requester_pays=False)):
                src_path = vsi_path
        else:
            src_path = raster_path
        
        # Open and read with rasterio
        with rasterio.open(src_path) as src:
            self.height = src.height
            self.width = src.width
            self.n_pixels = self.height * self.width
            self.crs = src.crs
            self.transform = src.transform
            
            bands_data = []
            band_names = []
            
            # Known NoData values from extractRasterData.py
            KNOWN_NODATA_VALUES = {
                -9999: 'default/elevation',
                -32768: 'soil',
                -3.4028235e+38: 'float32 nodata',
                255: 'byte nodata',
                0: 'some LULC/mask nodata'
            }
            
            for i in range(1, src.count + 1):
                # Get band data and metadata
                band_data = src.read(i).astype(np.float32)
                band_name = src.descriptions[i-1] if src.descriptions[i-1] else f'band_{i}'
                band_nodata = src.nodatavals[i-1]  # Get per-band NoData value
                
                # Convert NoData to NaN based on band type
                if band_nodata is not None:
                    band_data = np.where(band_data == band_nodata, np.nan, band_data)
                
                # Handle soil bands specifically (they use -32768)
                if any(soil_type in band_name.upper() for soil_type in ['CLAY', 'SAND', 'SILT']):
                    band_data = np.where(band_data == -32768, np.nan, band_data)
                    # Also mask impossible soil values
                    band_data = np.where((band_data < 0) | (band_data > 1000), np.nan, band_data)
                
                # Handle LULC (categorical data, but 0 might be valid)
                elif 'LULC' in band_name.upper():
                    # Only mask explicit NoData, not 0 (which might be a valid class)
                    if band_nodata is not None and band_nodata != 0:
                        pass  # Already handled above
                
                # Handle elevation anomalies (can be negative)
                elif 'anomaly' in band_name.lower():
                    # These can have negative values, only mask explicit NoData
                    pass  # Already handled above
                
                # Additional check for unrealistic values
                if band_name == 'elevation' or 'elevation' in band_name.lower():
                    # Mask unrealistic elevation values
                    band_data = np.where(band_data < -500, np.nan, band_data)  # Below Dead Sea
                    band_data = np.where(band_data > 9000, np.nan, band_data)  # Above Everest
                
                bands_data.append(band_data)
                band_names.append(band_name)
            
            # Stack all bands
            self.static_raster = np.stack(bands_data, axis=-1)
            self.static_bands = band_names
            
            # Create tabular version
            self.static_tabular = self.static_raster.reshape(-1, len(band_names))
            
            # Generate pixel coordinates
            self._generate_pixel_coordinates()
            
            # CRITICAL FIX: A pixel is valid if it has valid data in ANY band
            # Not all bands need to be valid (e.g., soil might be missing in water areas)
            self.pixel_valid_mask = ~np.all(np.isnan(self.static_tabular), axis=1)
            
            # Alternative: Require at least the elevation band to be valid
            # elevation_idx = band_names.index('elevation') if 'elevation' in band_names else 0
            # self.pixel_valid_mask = ~np.isnan(self.static_tabular[:, elevation_idx])
        
        logger.info(f"Loaded {len(self.static_bands)} static bands: {self.static_bands}")
        logger.info(f"Raster shape: {self.static_raster.shape}")
        logger.info(f"Valid pixels: {np.sum(self.pixel_valid_mask)}/{self.n_pixels}")
        
        # Debug information if very few valid pixels
        if np.sum(self.pixel_valid_mask) < self.n_pixels * 0.1:  # Less than 10% valid
            logger.warning(f"Only {np.sum(self.pixel_valid_mask)/self.n_pixels*100:.1f}% pixels are valid")
            for i, band_name in enumerate(self.static_bands):
                band_valid = ~np.isnan(self.static_tabular[:, i])
                logger.info(f"  Band {band_name}: {np.sum(band_valid)}/{self.n_pixels} valid pixels")
   
    def _load_from_local_file(self, filepath):
        """Helper to load from local file."""
        with rasterio.open(filepath) as src:
            self.height = src.height
            self.width = src.width
            self.n_pixels = self.height * self.width
            self.crs = src.crs
            self.transform = src.transform
            
            bands_data = []
            band_names = []
            
            for i in range(1, src.count + 1):
                band_data = src.read(i)
                band_name = src.descriptions[i-1] if src.descriptions[i-1] else f'band_{i}'
                
                if src.nodata is not None:
                    band_data = np.where(band_data == src.nodata, np.nan, band_data)
                
                bands_data.append(band_data)
                band_names.append(band_name)
            
            self.static_raster = np.stack(bands_data, axis=-1)
            self.static_bands = band_names
            self.static_tabular = self.static_raster.reshape(-1, len(band_names))
            self._generate_pixel_coordinates()
            self.pixel_valid_mask = ~np.any(np.isnan(self.static_tabular), axis=1)
            
    def load_dynamic_features(self, streamflow_path: Optional[str] = None):
        """
        Load dynamic climate, weather, and streamflow features from CSV files.
        Now uses Lat{latitude}_Lon{longitude} naming convention and loads all weather features.
        """
        if streamflow_path is None:
            streamflow_path = self.base_path / 'streamflowPredictions'
        else:
            # Ensure streamflow_path is a Path object
            streamflow_path = Path(streamflow_path)
        
        # UPDATED NAMING CONVENTION: Use Lat_Lon format instead of hybas_id
        # Load subwatershed predictions
        subwatershed_file = streamflow_path / f'subwatershed_predictions_Lat{self.metadata.lat_str}_Lon{self.metadata.lon_str}.csv'
        
        if subwatershed_file.exists():
            self.dynamic_subwatershed = pd.read_csv(subwatershed_file)
            
            # Ensure Date column exists and is datetime
            if 'Date' in self.dynamic_subwatershed.columns:
                self.dynamic_subwatershed['Date'] = pd.to_datetime(self.dynamic_subwatershed['Date'])
                self.dynamic_subwatershed.set_index('Date', inplace=True)
            elif 'date' in self.dynamic_subwatershed.columns:
                # Handle lowercase date column
                self.dynamic_subwatershed.rename(columns={'date': 'Date'}, inplace=True)
                self.dynamic_subwatershed['Date'] = pd.to_datetime(self.dynamic_subwatershed['Date'])
                self.dynamic_subwatershed.set_index('Date', inplace=True)
            
            # Log what columns were loaded (weather features + streamflow)
            weather_cols = [col for col in self.dynamic_subwatershed.columns if col != 'streamflow_pred']
            logger.info(f"Loaded subwatershed data: {self.dynamic_subwatershed.shape}")
            logger.info(f"  - Weather features: {weather_cols}")
            logger.info(f"  - Contains streamflow_pred: {'streamflow_pred' in self.dynamic_subwatershed.columns}")
        else:
            logger.warning(f"Subwatershed file not found: {subwatershed_file}")
        
        # UPDATED NAMING CONVENTION: Load river basin predictions
        basin_file = streamflow_path / f'river_basin_predictions_Lat{self.metadata.lat_str}_Lon{self.metadata.lon_str}.csv'
        
        if basin_file.exists():
            self.dynamic_basin = pd.read_csv(basin_file)
            
            # Ensure Date column exists and is datetime
            if 'Date' in self.dynamic_basin.columns:
                self.dynamic_basin['Date'] = pd.to_datetime(self.dynamic_basin['Date'])
                self.dynamic_basin.set_index('Date', inplace=True)
            elif 'date' in self.dynamic_basin.columns:
                # Handle lowercase date column
                self.dynamic_basin.rename(columns={'date': 'Date'}, inplace=True)
                self.dynamic_basin['Date'] = pd.to_datetime(self.dynamic_basin['Date'])
                self.dynamic_basin.set_index('Date', inplace=True)
            
            # Log what columns were loaded
            weather_cols = [col for col in self.dynamic_basin.columns if col != 'streamflow_pred']
            logger.info(f"Loaded river basin data: {self.dynamic_basin.shape}")
            logger.info(f"  - Weather features: {weather_cols}")
            logger.info(f"  - Contains streamflow_pred: {'streamflow_pred' in self.dynamic_basin.columns}")
        else:
            logger.warning(f"River basin file not found: {basin_file}")
    
    def get_dynamic_features_summary(self) -> Dict:
        """
        Get summary of all available dynamic features.
        
        Returns:
        --------
        Dict with feature information for both subwatershed and basin levels
        """
        summary = {}
        
        if self.dynamic_subwatershed is not None:
            summary['subwatershed'] = {
                'n_timesteps': len(self.dynamic_subwatershed),
                'weather_features': [col for col in self.dynamic_subwatershed.columns 
                                    if col not in ['streamflow_pred']],
                'has_streamflow': 'streamflow_pred' in self.dynamic_subwatershed.columns,
                'date_range': (self.dynamic_subwatershed.index.min(), 
                              self.dynamic_subwatershed.index.max()) if len(self.dynamic_subwatershed) > 0 else None
            }
        
        if self.dynamic_basin is not None:
            summary['river_basin'] = {
                'n_timesteps': len(self.dynamic_basin),
                'weather_features': [col for col in self.dynamic_basin.columns 
                                    if col not in ['streamflow_pred']],
                'has_streamflow': 'streamflow_pred' in self.dynamic_basin.columns,
                'date_range': (self.dynamic_basin.index.min(), 
                              self.dynamic_basin.index.max()) if len(self.dynamic_basin) > 0 else None
            }
        
        return summary

   
    def load_flood_targets(self, flood_path: Optional[str] = None):
        """
        Load and process OPERA DSWx flood maps.
        Convert to binary classification: water (1,2,3) vs not_water (0).
        """
        if flood_path is None:
            flood_file = (self.base_path / 'floodMaps' / 
                         f'dswx_s1_timeseries_subwatershed_Lon{self.metadata.lon_str}_'
                         f'Lat{self.metadata.lat_str}.tif')
        else:
            flood_file = Path(flood_path)
        
        if not flood_file.exists():
            logger.warning(f"Flood map file not found: {flood_file}")
            # Initialize with empty arrays so feature processor doesn't fail
            self.flood_masks_spatial = np.array([])
            self.flood_masks_tabular = np.array([])
            self.flood_dates = []
            self.valid_timesteps = np.array([])
            return
        
        logger.info(f"Loading flood targets from {flood_file}")
        
        with rasterio.open(str(flood_file)) as src:
            flood_masks = []
            flood_dates = []
            
            for i in range(1, src.count + 1):
                band_data = src.read(i)
                band_desc = src.descriptions[i-1] if src.descriptions[i-1] else f'band_{i}'
                
                # Extract date from band description (format: DSWx_S1_YYYY-MM-DD)
                try:
                    date_str = band_desc.split('_')[-1]
                    flood_date = pd.to_datetime(date_str)
                    flood_dates.append(flood_date)
                except:
                    flood_dates.append(None)
                
                # Convert to binary classification
                # Water classes: 1, 2, 3
                # Not water: 0
                # NoData: 255 -> -1 (will be masked)
                
                binary_mask = np.full_like(band_data, -1, dtype=np.float32)
                binary_mask[np.isin(band_data, [1, 2, 3])] = 1.0  # Water
                binary_mask[band_data == 0] = 0.0  # Not water
                
                flood_masks.append(binary_mask)
            
            # Stack temporal dimension (T, H, W)
            self.flood_masks_spatial = np.stack(flood_masks, axis=0)
            self.flood_dates = flood_dates
            
            # Create flattened version (T, N_pixels)
            self.flood_masks_tabular = self.flood_masks_spatial.reshape(
                len(flood_masks), -1
            )
            
            # Create valid timesteps mask (timesteps with at least some valid data)
            self.valid_timesteps = np.any(
                self.flood_masks_tabular != -1, axis=1
            )
            
            logger.info(f"Loaded {len(flood_masks)} flood masks")
            logger.info(f"Valid timesteps: {np.sum(self.valid_timesteps)}/{len(flood_masks)}")
    
    def _generate_pixel_coordinates(self):
        """Generate lat/lon coordinates for each pixel."""
        rows, cols = np.meshgrid(range(self.height), range(self.width), indexing='ij')
        
        # Convert pixel indices to coordinates using transform
        xs, ys = rasterio.transform.xy(
            self.transform, rows.flatten(), cols.flatten()
        )
        
        self.pixel_coords = np.column_stack([xs, ys])
    
    def get_aligned_features_for_date(self, target_date: datetime, 
                                     lag_days: int = 3) -> Dict:
        """
        Get aligned features for a specific flood observation date.
        Includes all weather features and streamflow predictions.
        
        Parameters:
        -----------
        target_date : datetime
            Date of flood observation
        lag_days : int
            Number of days to look back for dynamic features
            
        Returns:
        --------
        Dict with aligned features including weather and streamflow data
        """
        aligned_data = {
            'static': self.static_tabular,
            'target_date': target_date
        }
        
        # Get dynamic features with lag
        if self.dynamic_subwatershed is not None:
            date_range = pd.date_range(
                end=target_date, periods=lag_days + 1, freq='D'
            )
            
            # Extract features for date range
            try:
                # Get all columns (weather + streamflow) for the date range
                subwatershed_features = self.dynamic_subwatershed.loc[date_range]
                
                # Also get basin features if available
                basin_features = None
                if self.dynamic_basin is not None:
                    basin_features = self.dynamic_basin.loc[date_range]
                
                aligned_data['dynamic_subwatershed'] = subwatershed_features
                aligned_data['dynamic_basin'] = basin_features
                
                # Add summary of what features are included
                aligned_data['feature_summary'] = {
                    'subwatershed_features': list(subwatershed_features.columns),
                    'basin_features': list(basin_features.columns) if basin_features is not None else [],
                    'n_days': len(date_range)
                }
                
            except KeyError as e:
                logger.warning(f"Missing dynamic data for {target_date}: {e}")
                # Try to get partial data if some dates are missing
                available_dates = []
                for date in date_range:
                    if date in self.dynamic_subwatershed.index:
                        available_dates.append(date)
                
                if available_dates:
                    logger.info(f"Found {len(available_dates)}/{len(date_range)} dates")
                    aligned_data['dynamic_subwatershed'] = self.dynamic_subwatershed.loc[available_dates]
                    if self.dynamic_basin is not None:
                        aligned_data['dynamic_basin'] = self.dynamic_basin.loc[available_dates]
                else:
                    aligned_data['dynamic_subwatershed'] = None
                    aligned_data['dynamic_basin'] = None
        
        return aligned_data
