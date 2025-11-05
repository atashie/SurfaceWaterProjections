# feature_processor.py

import numpy as np
import pandas as pd
from scipy import ndimage
from typing import Dict, List, Optional, Tuple
import logging
from datetime import datetime

# Import FloodDataContainer from flood_data_container module
from flood_data_container import FloodDataContainer

logger = logging.getLogger(__name__)


class TabularFeatureProcessor:
    """
    Feature engineering for tree-based and traditional ML models.
    """
    
    def __init__(self, 
                 lag_days: int = 3,
                 window_sizes: List[int] = [3, 5, 7],
                 include_temporal_features: bool = True,
                 include_spatial_context: bool = True):
        """
        Initialize feature processor.
        
        Parameters:
        -----------
        lag_days : int
            Number of days to look back for dynamic features
        window_sizes : List[int]
            Window sizes for spatial context features
        include_temporal_features : bool
            Whether to include temporal encoding features
        include_spatial_context : bool
            Whether to include neighborhood statistics
        """
        self.lag_days = lag_days
        self.window_sizes = window_sizes
        self.include_temporal_features = include_temporal_features
        self.include_spatial_context = include_spatial_context
        
    def process_watershed(self, container: FloodDataContainer) -> pd.DataFrame:
        """
        Process all features for a watershed into tabular format.
        
        Returns:
        --------
        DataFrame with features and targets for all valid samples
        """
        all_samples = []
        
        # Iterate through valid timesteps
        for t_idx, (flood_date, is_valid) in enumerate(
            zip(container.flood_dates, container.valid_timesteps)
        ):
            if not is_valid or flood_date is None:
                continue
            
            # Get aligned features for this date
            aligned_features = container.get_aligned_features_for_date(
                flood_date, self.lag_days
            )
            
            if aligned_features['dynamic_subwatershed'] is None:
                continue
            
            # Process features for this timestep
            timestep_features = self._process_timestep(
                container, aligned_features, t_idx, flood_date
            )
            
            if timestep_features is not None:
                all_samples.append(timestep_features)
        
        if not all_samples:
            logger.warning(f"No valid samples for watershed {container.metadata.hybas_id}")
            return pd.DataFrame()
        
        # Combine all samples
        df = pd.concat(all_samples, ignore_index=True)
        
        logger.info(f"Processed {len(df)} samples for watershed {container.metadata.hybas_id}")
        logger.info(f"Feature dimensions: {df.shape}")
        
        return df
    
    def _process_timestep(self, container: FloodDataContainer,
                         aligned_features: Dict,
                         t_idx: int,
                         flood_date: datetime) -> Optional[pd.DataFrame]:
        """Process features for a single timestep."""
        
        # Get target for this timestep
        target_flat = container.flood_masks_tabular[t_idx]
        
        # Get valid pixels (not NoData in static features or target)
        valid_mask = container.pixel_valid_mask & (target_flat != -1)
        
        if np.sum(valid_mask) == 0:
            return None
        
        # Extract valid pixels
        valid_indices = np.where(valid_mask)[0]
        n_valid = len(valid_indices)
        
        # Initialize feature dictionary
        features = {
            'pixel_idx': valid_indices,
            'timestep_idx': [t_idx] * n_valid,
            'date': [flood_date] * n_valid,
            'watershed_id': [container.metadata.hybas_id] * n_valid,
            'target': target_flat[valid_mask]
        }
        
        # Add static features
        static_features = aligned_features['static'][valid_mask]
        for i, band_name in enumerate(container.static_bands):
            features[f'static_{band_name}'] = static_features[:, i]
        
        # Add spatial context if requested
        if self.include_spatial_context:
            spatial_features = self._compute_spatial_context(
                container, valid_indices
            )
            features.update(spatial_features)
        
        # Add dynamic features
        dynamic_features = self._process_dynamic_features(
            aligned_features, n_valid
        )
        features.update(dynamic_features)
        
        # Add temporal features if requested
        if self.include_temporal_features:
            temporal_features = self._compute_temporal_features(
                flood_date, n_valid
            )
            features.update(temporal_features)
        
        return pd.DataFrame(features)
    
    def _compute_spatial_context(self, container: FloodDataContainer,
                                valid_indices: np.ndarray) -> Dict:
        """Compute neighborhood statistics for each pixel."""
        
        spatial_features = {}
        
        # Convert flat indices back to 2D coordinates
        rows = valid_indices // container.width
        cols = valid_indices % container.width
        
        # Process key static bands for spatial context
        key_bands = ['elevation', 'slope', 'HAND']
        
        for band_name in key_bands:
            if band_name not in container.static_bands:
                continue
            
            band_idx = container.static_bands.index(band_name)
            band_data = container.static_raster[:, :, band_idx]
            
            for window_size in self.window_sizes:
                # Compute neighborhood statistics
                mean_filter = ndimage.uniform_filter(
                    band_data, size=window_size, mode='constant', cval=np.nan
                )
                std_filter = ndimage.generic_filter(
                    band_data, np.nanstd, size=window_size, mode='constant', cval=np.nan
                )
                
                # Extract values at valid pixels
                spatial_features[f'{band_name}_mean_{window_size}x{window_size}'] = mean_filter[rows, cols]
                spatial_features[f'{band_name}_std_{window_size}x{window_size}'] = std_filter[rows, cols]
        
        return spatial_features
    
    def _process_dynamic_features(self, aligned_features: Dict, 
                                 n_valid: int) -> Dict:
        """Process dynamic climate and streamflow features."""
        
        dynamic_features = {}
        
        # Subwatershed-level features
        if aligned_features['dynamic_subwatershed'] is not None:
            df_sub = aligned_features['dynamic_subwatershed']
            
            for col in df_sub.columns:
                if col == 'Date':
                    continue
                
                # Current value (at target date)
                dynamic_features[f'sub_{col}_t0'] = [df_sub[col].iloc[-1]] * n_valid
                
                # Lagged values
                for lag in range(1, min(self.lag_days + 1, len(df_sub))):
                    dynamic_features[f'sub_{col}_t-{lag}'] = [df_sub[col].iloc[-1-lag]] * n_valid
                
                # Aggregations over window
                dynamic_features[f'sub_{col}_mean_{self.lag_days}d'] = [df_sub[col].mean()] * n_valid
                dynamic_features[f'sub_{col}_max_{self.lag_days}d'] = [df_sub[col].max()] * n_valid
                dynamic_features[f'sub_{col}_sum_{self.lag_days}d'] = [df_sub[col].sum()] * n_valid
        
        # Basin-level features
        if aligned_features['dynamic_basin'] is not None:
            df_basin = aligned_features['dynamic_basin']
            
            for col in df_basin.columns:
                if col == 'Date':
                    continue
                
                # Current value
                dynamic_features[f'basin_{col}_t0'] = [df_basin[col].iloc[-1]] * n_valid
                
                # Aggregations
                dynamic_features[f'basin_{col}_mean_{self.lag_days}d'] = [df_basin[col].mean()] * n_valid
        
        return dynamic_features
    
    def _compute_temporal_features(self, flood_date: datetime, 
                                  n_valid: int) -> Dict:
        """Compute temporal encoding features."""
        
        temporal_features = {
            'year': [flood_date.year] * n_valid,
            'month': [flood_date.month] * n_valid,
            'day_of_year': [flood_date.timetuple().tm_yday] * n_valid,
            'season': [self._get_season(flood_date.month)] * n_valid,
            'is_wet_season': [self._is_wet_season(flood_date.month)] * n_valid
        }
        
        # Add cyclical encoding for month and day of year
        temporal_features['month_sin'] = [np.sin(2 * np.pi * flood_date.month / 12)] * n_valid
        temporal_features['month_cos'] = [np.cos(2 * np.pi * flood_date.month / 12)] * n_valid
        
        day_of_year = flood_date.timetuple().tm_yday
        temporal_features['doy_sin'] = [np.sin(2 * np.pi * day_of_year / 365)] * n_valid
        temporal_features['doy_cos'] = [np.cos(2 * np.pi * day_of_year / 365)] * n_valid
        
        return temporal_features
    
    @staticmethod
    def _get_season(month: int) -> int:
        """Get season (0-3) from month."""
        if month in [12, 1, 2]:
            return 0  # Winter
        elif month in [3, 4, 5]:
            return 1  # Spring
        elif month in [6, 7, 8]:
            return 2  # Summer
        else:
            return 3  # Fall
    
    @staticmethod
    def _is_wet_season(month: int) -> int:
        """Determine if month is in wet season (varies by region)."""
        # This is a simplified example - adjust based on your region
        return 1 if month in [3, 4, 5, 6, 7, 8] else 0


# Placeholder classes for future implementation
class SpatialFeatureProcessor:
    """Placeholder for CNN/U-Net feature processing."""
    def __init__(self):
        raise NotImplementedError("CNN/U-Net processing to be implemented")


class GraphFeatureProcessor:
    """Placeholder for GNN feature processing."""
    def __init__(self):
        raise NotImplementedError("GNN processing to be implemented")
