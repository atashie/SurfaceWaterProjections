# train.py
from torch.cuda.amp import autocast, GradScaler
import argparse
import json
import os
import sys
import psutil
import torch
import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple, Union
import pytorch_lightning as pl
from torch.optim.lr_scheduler import ReduceLROnPlateau
#torch.set_num_threads(16)  # Limit PyTorch threads
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from datetime import datetime
#os.environ["OMP_NUM_THREADS"] = "4"  # Limit OpenMP threads
#os.environ["MKL_NUM_THREADS"] = "4"   # Limit MKL threads
#os.environ["NUMEXPR_NUM_THREADS"] = "4"  # Limit NumExpr threads
import matplotlib.pyplot as plt
import time
import seaborn as sns
from tqdm import tqdm  
import pickle
import signal  # signal handling for spot instances

#########################################################################################################################
# Caravan data processing
#########################################################################################################################

import boto3
import re
import xarray as xr
import s3fs  # Required by xarray to open S3 files directly
import warnings
import gc


def identify_all_available_watersheds(bucket_name, base_directory):
    """
    Identify all available watersheds in the specified S3 bucket and directory,
    returning their subdirectory and ID in a DataFrame.

    Args:
        bucket_name (str): Name of the S3 bucket.
        base_directory (str): Base directory in the S3 bucket where watershed files are stored.
                               Expected structure: .../netcdf/{subdirectory_name}/{subdirectory_name}_{watershedID}.nc

    Returns:
        pandas.DataFrame: A DataFrame with columns 'subdirectory_name' and 'watershedID'.
                          Returns an empty DataFrame if no matching files are found.
    """
    # Initialize S3 client
    s3 = boto3.client('s3')

    # List to store watershed data dictionaries
    watershed_data = []

    # Regex to capture subdirectory name and watershed ID
    # Pattern breakdown:
    # netcdf/       : Matches the literal 'netcdf/'
    # ([^/]+)       : Captures the subdirectory name (group 1) - one or more characters not being '/'
    # /             : Matches the literal '/'
    # [^/]+         : Matches the filename prefix (e.g., subdirectory_name again)
    # _             : Matches the literal '_'
    # (\d+)         : Captures the watershed ID (group 2) - one or more digits
    # \.nc$         : Matches the literal '.nc' at the end of the string
    pattern = re.compile(rf'{base_directory.strip("/")}/([^/]+)/[^/]+_(\d+)\.nc$')

    # Paginate through the objects in the bucket
    paginator = s3.get_paginator('list_objects_v2')
    operation_parameters = {'Bucket': bucket_name, 'Prefix': base_directory}

    print(f"Searching for watershed files in s3://{bucket_name}/{base_directory}...")

    for page in paginator.paginate(**operation_parameters):
        if 'Contents' in page:
            for obj in page['Contents']:
                file_path = obj['Key']
                # Match the pattern
                match = pattern.search(file_path)
                if match:
                    subdirectory_name = match.group(1)
                    watershed_id = match.group(2)
                    watershed_data.append({
                        'subdirectory_name': subdirectory_name,
                        'watershedID': watershed_id
                    })

    # Convert the list of dictionaries to a Pandas DataFrame
    if watershed_data:
        df_watersheds = pd.DataFrame(watershed_data)
        # Optional: Remove potential duplicate entries if the same file path was listed multiple times
        df_watersheds = df_watersheds.drop_duplicates().reset_index(drop=True)
        print(f"Found {len(df_watersheds)} unique watershed files.")
    else:
        df_watersheds = pd.DataFrame(columns=['subdirectory_name', 'watershedID'])
        print("No watershed files found matching the pattern.")

    return df_watersheds



import scipy # Required for the engine
import traceback # For detailed error reporting

def load_and_process_watershed_data(bucket_name, base_directory, subdirectory_name, watershedID, variables_to_extract=None):
    """
    Loads and processes timeseries data for a specific watershed from an S3 NetCDF file.
    Allows specifying variables to extract, with a default list if none provided.
    (Uses explicit s3fs object and scipy engine for NetCDF Classic files)

    Args:
        bucket_name (str): Name of the S3 bucket.
        base_directory (str): Base directory in the S3 bucket where watershed files are stored.
        subdirectory_name (str): The specific subdirectory identifier (e.g., 'hysets').
        watershedID (str or int): The ID of the watershed.
        variables_to_extract (list, optional): A list of variable names (strings) to extract
                                               from the NetCDF file. Defaults to a predefined list
                                               if None is provided.

    Returns:
        pandas.DataFrame: A DataFrame containing the processed timeseries data for the
                          specified variables. Returns an empty DataFrame on error or if
                          no valid data is found.
    """
    s3_path = f"s3://{bucket_name}/{base_directory.strip('/')}/{subdirectory_name}/{subdirectory_name}_{watershedID}.nc"
    print(f"Attempting to open: {s3_path}")
    print(f"Using s3fs and engine='scipy'") # Specify engine being used

    # --- Define default variables and handle user input ---
    default_variables = [
        'streamflow', 
#        'potential_evaporation_sum_ERA5_LAND', 
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min', 
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean'
    ]
    if variables_to_extract is None:
        # Use the default list if the user didn't provide one
        variables_to_extract = default_variables
        print(f"Using default variables: {variables_to_extract}")
    else:
        # Use the list provided by the user
        print(f"Using user-specified variables: {variables_to_extract}")
    # --- End variable handling ---

    processed_df = pd.DataFrame()

    try:
        # --- Explicitly create an s3fs filesystem object ---
        s3 = s3fs.S3FileSystem(anon=False) # Assumes credentials from environment/role

        # --- Open the dataset using the s3fs object and scipy engine ---
        with s3.open(s3_path, mode='rb') as s3file: # Open in binary read mode
            # Pass the file-like object to xarray with the scipy engine
            with xr.open_dataset(s3file, engine='scipy', cache=False) as ds:
                print(f"Successfully opened via s3fs and scipy: {s3_path}")

                # --- Check if variables exist before trying to subset ---
                available_vars = list(ds.data_vars)
                # Use the (potentially user-defined) variables_to_extract list here
                missing_vars = [v for v in variables_to_extract if v not in available_vars]
                if missing_vars:
                    print(f"Warning: The following requested variables are missing from the file: {missing_vars}")
                    # Adjust the list of variables to extract to only include available ones for this operation
                    current_variables_to_extract = [v for v in variables_to_extract if v in available_vars]
                    if not current_variables_to_extract:
                         print("Error: None of the requested variables are available. Cannot proceed.")
                         return pd.DataFrame()
                else:
                    # All requested variables are available
                    current_variables_to_extract = variables_to_extract

                # --- Subset using the validated list ---
                ds_subset = ds[current_variables_to_extract]
                df = ds_subset.to_dataframe().reset_index()

                # --- Date Handling ---
                if 'date' not in df.columns:
                     # Scipy engine might name the time dimension differently, check coordinates
                     time_coord_name = None
                     if 'time' in ds.coords: # Common name
                         time_coord_name = 'time'
                     elif ds.coords: # If not 'time', maybe there's only one coord?
                         potential_time_coords = [c for c in ds.coords if ds[c].ndim == 1] # Look for 1D coords
                         if len(potential_time_coords) == 1:
                             time_coord_name = potential_time_coords[0]
                             print(f"Using coordinate '{time_coord_name}' as date.")
                             df = df.rename(columns={time_coord_name: 'date'}) # Rename it

                     if 'date' not in df.columns:
                         raise ValueError("Date coordinate not found or could not be identified.")

                # Convert date if it's not already datetime (scipy might not auto-decode)
                if not pd.api.types.is_datetime64_any_dtype(df['date']):
                    print("Attempting to decode date coordinate...")
                    try:
                        # Get units and calendar from dataset attributes if possible
                        # Use the identified time coordinate name or default to 'date'
                        time_var_name = time_coord_name if time_coord_name else 'date'
                        time_var = ds[time_var_name]
                        units = time_var.attrs.get('units', None)
                        calendar = time_var.attrs.get('calendar', 'standard') # Default calendar

                        if units:
                            print(f"Using units: '{units}', calendar: '{calendar}'")
                            # Use xarray's utility for robust decoding
                            # Pass df['date'].values for potentially better performance
                            decoded_dates = xr.coding.times.decode_cf_datetime(df['date'].values, units, calendar)
                            df['date'] = decoded_dates
                        else:
                            raise ValueError(f"Could not find 'units' attribute for date coordinate '{time_var_name}'.")

                        if not pd.api.types.is_datetime64_any_dtype(df['date']):
                             raise TypeError("Date conversion failed.")
                        print("Date decoding successful.")

                    except Exception as date_err:
                         print(f"Error decoding date: {date_err}")
                         warnings.warn("Date column may not be in datetime format.")

                print(f"Initial data loaded. Shape: {df.shape}")

                # --- Streamflow Trimming (only if streamflow was requested and exists) ---
                # Check if 'streamflow' is in the list of variables we actually extracted
                if 'streamflow' in current_variables_to_extract and 'streamflow' in df.columns:
                    first_valid_idx = df['streamflow'].first_valid_index()
                    last_valid_idx = df['streamflow'].last_valid_index()

                    if first_valid_idx is None or last_valid_idx is None:
                        print("Warning: No valid streamflow data found. Returning empty DataFrame.")
                        return pd.DataFrame()

                    print(f"Trimming data based on streamflow: First valid at index {first_valid_idx}, Last valid at index {last_valid_idx}")
                    df_trimmed = df.loc[first_valid_idx:last_valid_idx].reset_index(drop=True)
                    print(f"Data shape after trimming: {df_trimmed.shape}")
                else:
                    # Streamflow wasn't requested or wasn't found in the loaded data, skip trimming
                    df_trimmed = df.copy()
                    if 'streamflow' not in current_variables_to_extract:
                        print("Skipping streamflow NaN trimming as 'streamflow' was not requested.")
                    else: # Streamflow was requested but not found in df.columns (should be caught earlier, but defensive check)
                        print("Warning: 'streamflow' requested but not found in loaded columns. Skipping trimming.")


                # --- Check remaining NaNs (on df_trimmed) ---
                remaining_nans = df_trimmed.isnull().sum()
                remaining_nans = remaining_nans[remaining_nans > 0]

                if not remaining_nans.empty:
                    print("\nWarning: Remaining NaNs found after trimming:")
                else:
                    print("\nNo remaining NaNs found in the trimmed data.")

                processed_df = df_trimmed

                    # --- Fill any remaining NaNs with column means ---
                if processed_df.isnull().values.any():
                    print("Filling remaining NaNs with column means...")
                    # Only fill numeric columns
                    numeric_cols = processed_df.select_dtypes(include=[np.number]).columns
                    for col in numeric_cols:
                        if processed_df[col].isnull().any():
                            mean_val = processed_df[col].mean()
                            processed_df[col] = processed_df[col].fillna(mean_val)
                    # Optionally, warn if any NaNs remain (shouldn't happen unless all values were NaN)
                    if processed_df.isnull().values.any():
                        print("Warning: Some NaNs remain even after filling with means.")
                    else:
                        print("All NaNs filled with column means.")
                

    except FileNotFoundError:
        print(f"Error: File not found at {s3_path} (reported by s3fs)")
        return pd.DataFrame()
    except KeyError as e:
        # This might occur if a variable is missing AND the check above fails somehow,
        # or if date handling fails to find the coordinate.
        print(f"Error: Key {e} not found or issue during processing.")
        traceback.print_exc()
        return pd.DataFrame()
    except ValueError as e:
        # Catch specific ValueErrors like the date coordinate issue
        print(f"Error: {e}")
        traceback.print_exc()
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred while processing {s3_path}: {e}")
        traceback.print_exc() # Print full traceback for debugging
        return pd.DataFrame()

    return processed_df



def get_watershed_attributes(bucket_name, base_attributes_directory, subdirectory_name, watershedID, hydroatlas_cols_to_extract=None, caravan_cols_to_extract=None):
    """
    Extracts static attributes for a specific watershed from up to three CSV files on S3.
    Allows specifying HydroATLAS and Caravan columns to extract, with defaults.

    Args:
        bucket_name (str): Name of the S3 bucket.
        base_attributes_directory (str): Base directory in S3 for attribute files
                                         (e.g., "Caravan-Jan25-nc/attributes/").
        subdirectory_name (str): The specific subdirectory identifier (e.g., 'hysets', 'camels_br').
        watershedID (str or int): The ID of the watershed.
        hydroatlas_cols_to_extract (list, optional): A list of HydroATLAS column names (strings)
                                                     to extract from the hydroatlas attributes file.
                                                     Defaults to a predefined list if None is provided.
        caravan_cols_to_extract (list, optional): A list of Caravan column names (strings)
                                                  to extract from the caravan attributes file.
                                                  Defaults to None (file is not read).

    Returns:
        pandas.Series: A Series containing the combined attributes for the watershed from
                       all requested and successfully processed files.
                       Returns an empty Series if the watershed is not found in any requested file
                       or if required files/columns are missing.
    """
    # --- Construct Paths and Target ID ---
    # Check if base_attributes_directory is a local path or S3 path
    if base_attributes_directory.startswith('/'):  # Local path
        # For local paths, construct direct file paths
        base_path = f"{base_attributes_directory}/{subdirectory_name}"
        file1_path = f"{base_path}/attributes_hydroatlas_{subdirectory_name}.csv"
        file2_path = f"{base_path}/attributes_other_{subdirectory_name}.csv"
        file3_path = f"{base_path}/attributes_caravan_{subdirectory_name}.csv"
    else:  # S3 path
        # For S3 paths, construct S3 URLs
        base_path = f"s3://{bucket_name}/{base_attributes_directory.strip('/')}/{subdirectory_name}"
        file1_path = f"{base_path}/attributes_hydroatlas_{subdirectory_name}.csv"
        file2_path = f"{base_path}/attributes_other_{subdirectory_name}.csv"
        file3_path = f"{base_path}/attributes_caravan_{subdirectory_name}.csv"
    
    target_gauge_id = f"{subdirectory_name}_{watershedID}"

    print(f"Searching for attributes for gauge_id: {target_gauge_id}")

    # --- Define Default HydroATLAS Columns and Handle User Input ---
    default_hydroatlas_cols = [
        'pet_mm_syr', 'aet_mm_syr', 
        'pre_mm_syr', 
        'tmp_dc_syr', 'tmp_dc_smn', 'tmp_dc_smx',
        'snw_pc_syr', 'snw_pc_smx', 
        'swc_pc_syr',
        'gdp_ud_sav', 'gdp_ud_ssu',
        'inu_pc_slt', 'inu_pc_smn', 'inu_pc_smx',
        'run_mm_syr',
        'dis_m3_pmn', 'dis_m3_pmx', 'dis_m3_pyr',
        'lka_pc_sse',
        'lkv_mc_usu',
        'rev_mc_usu',
        'ria_ha_usu', 'riv_tc_usu', #river area and volume
        'dor_pc_pva', # degree of regulation
        'gwt_cm_sav',
        'ele_mt_sav', 'slp_dg_sav', 
        'sgr_dk_sav', #stream gradient
        'ari_ix_sav', 'glc_cl_smj', 'for_pc_sse',
        'crp_pc_sse', 'ire_pc_sse', 'gla_pc_sse', 'prm_pc_sse',
        'cly_pc_sav', 'slt_pc_sav', 'snd_pc_sav', 
        'soc_th_sav', 'lit_cl_smj', 'kar_pc_sse', 'pop_ct_usu',
        'rdd_mk_sav', 'urb_pc_sse',
    ]
    if hydroatlas_cols_to_extract is None:
        hydroatlas_cols_to_use = default_hydroatlas_cols
        print(f"Using default HydroATLAS columns.")
    else:
        hydroatlas_cols_to_use = hydroatlas_cols_to_extract
        print(f"Using user-specified HydroATLAS columns.")
    # --- End HydroATLAS column handling ---

    # --- Define Other Columns (remains fixed) ---
    other_cols = ['gauge_lat', 'gauge_lon', 'area']

    # --- Initialize empty Series for results ---
    attributes_series1 = pd.Series(dtype=object) # HydroATLAS
    attributes_series2 = pd.Series(dtype=object) # Other
    attributes_series3 = pd.Series(dtype=object) # Caravan (Optional)

    # --- Process File 1 (HydroATLAS Attributes) ---
    print(f"Reading HydroATLAS attributes from: {file1_path}")
    try:
        df1 = pd.read_csv(file1_path, low_memory=False)
        if 'gauge_id' not in df1.columns: raise KeyError("'gauge_id' column not found")

        missing_cols1 = [col for col in hydroatlas_cols_to_use if col not in df1.columns]
        current_hydroatlas_cols = [col for col in hydroatlas_cols_to_use if col in df1.columns]

        if missing_cols1:
            warnings.warn(f"Missing requested HydroATLAS columns in {file1_path}: {missing_cols1}")
        if not current_hydroatlas_cols:
            print("Error: None of the requested HydroATLAS columns are available. Skipping HydroATLAS file.")
        else:
            df1_filtered = df1[df1['gauge_id'] == target_gauge_id]
            if len(df1_filtered) == 1:
                attributes_series1 = df1_filtered[current_hydroatlas_cols].iloc[0]
                print("Found matching watershed in HydroATLAS attributes file.")
            elif len(df1_filtered) > 1:
                print(f"Warning: Multiple entries for gauge_id '{target_gauge_id}' in {file1_path}. Using first.")
                attributes_series1 = df1_filtered[current_hydroatlas_cols].iloc[0]
            else: # len == 0
                 print(f"Warning: gauge_id '{target_gauge_id}' not found in {file1_path}")

    except FileNotFoundError: print(f"Error: File not found: {file1_path}")
    except KeyError as e: print(f"Error: {e} in {file1_path}.")
    except Exception as e: print(f"An error occurred while processing {file1_path}: {e}")

    # --- Process File 2 (Other Attributes) ---
    print(f"Reading Other attributes from: {file2_path}")
    try:
        df2 = pd.read_csv(file2_path, low_memory=False)
        if 'gauge_id' not in df2.columns: raise KeyError("'gauge_id' column not found")

        missing_cols2 = [col for col in other_cols if col not in df2.columns]
        current_other_cols = [col for col in other_cols if col in df2.columns]

        if missing_cols2:
            warnings.warn(f"Missing required 'other' columns in {file2_path}: {missing_cols2}")
        if not current_other_cols:
            print("Error: None of the required 'other' columns are available. Skipping Other attributes file.")
        else:
            df2_filtered = df2[df2['gauge_id'] == target_gauge_id]
            if len(df2_filtered) == 1:
                attributes_series2 = df2_filtered[current_other_cols].iloc[0]
                print("Found matching watershed in Other attributes file.")
            elif len(df2_filtered) > 1:
                print(f"Warning: Multiple entries for gauge_id '{target_gauge_id}' in {file2_path}. Using first.")
                attributes_series2 = df2_filtered[current_other_cols].iloc[0]
            else: # len == 0
                print(f"Warning: gauge_id '{target_gauge_id}' not found in {file2_path}")

    except FileNotFoundError: print(f"Error: File not found: {file2_path}")
    except KeyError as e: print(f"Error: {e} in {file2_path}.")
    except Exception as e: print(f"An error occurred while processing {file2_path}: {e}")

    # --- Process File 3 (Caravan Attributes - Optional) ---
    if caravan_cols_to_extract is not None:
        print(f"Reading Caravan attributes from: {file3_path}")
        if not isinstance(caravan_cols_to_extract, list) or not all(isinstance(item, str) for item in caravan_cols_to_extract):
             print("Warning: 'caravan_cols_to_extract' provided but is not a list of strings. Skipping Caravan file.")
        else:
            try:
                df3 = pd.read_csv(file3_path, low_memory=False)
                if 'gauge_id' not in df3.columns: raise KeyError("'gauge_id' column not found")

                missing_cols3 = [col for col in caravan_cols_to_extract if col not in df3.columns]
                current_caravan_cols = [col for col in caravan_cols_to_extract if col in df3.columns]

                if missing_cols3:
                    warnings.warn(f"Missing requested Caravan columns in {file3_path}: {missing_cols3}")
                if not current_caravan_cols:
                    print("Error: None of the requested Caravan columns are available. Skipping Caravan file.")
                else:
                    df3_filtered = df3[df3['gauge_id'] == target_gauge_id]
                    if len(df3_filtered) == 1:
                        attributes_series3 = df3_filtered[current_caravan_cols].iloc[0]
                        print("Found matching watershed in Caravan attributes file.")
                    elif len(df3_filtered) > 1:
                        print(f"Warning: Multiple entries for gauge_id '{target_gauge_id}' in {file3_path}. Using first.")
                        attributes_series3 = df3_filtered[current_caravan_cols].iloc[0]
                    else: # len == 0
                        print(f"Warning: gauge_id '{target_gauge_id}' not found in {file3_path}")

            except FileNotFoundError: print(f"Error: File not found: {file3_path}")
            except KeyError as e: print(f"Error: {e} in {file3_path}.")
            except Exception as e: print(f"An error occurred while processing {file3_path}: {e}")
    else:
        print("Skipping Caravan attributes file as 'caravan_cols_to_extract' was not provided.")


    # --- Combine Results ---
    # Collect all non-empty series into a list
    all_series = []
    if not attributes_series1.empty:
        all_series.append(attributes_series1)
    if not attributes_series2.empty:
        all_series.append(attributes_series2)
    if not attributes_series3.empty: # Will only be non-empty if requested and successful
        all_series.append(attributes_series3)

    # Concatenate the collected series
    if all_series:
        # Check for duplicate index entries (column names) before concatenating
        combined_index = pd.concat(all_series).index
        if combined_index.duplicated().any():
            duplicates = combined_index[combined_index.duplicated()].unique().tolist()
            warnings.warn(f"Duplicate attribute names found across files: {duplicates}. Values from later files in the sequence (HydroATLAS -> Other -> Caravan) will overwrite earlier ones.")
            # Keep the last occurrence in case of duplicates
            combined_attributes = pd.concat(all_series)
            combined_attributes = combined_attributes[~combined_attributes.index.duplicated(keep='last')]
        else:
            combined_attributes = pd.concat(all_series)

        print(f"Successfully combined attributes from {len(all_series)} file(s).")
        return combined_attributes
    else:
        print("Error: Could not retrieve attributes from any requested file for the specified watershed.")
        return pd.Series(dtype=object) # Return empty Series on complete failure






#########################################################################################################################
# Modeling - Helper Functions
#########################################################################################################################

def calculate_kge(observed: np.ndarray, simulated: np.ndarray) -> float:
    """
    Calculate Kling-Gupta Efficiency (KGE) with robust handling of extreme values.
    """
    # Ensure inputs are numpy arrays
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()
    
    # Remove NaN and inf values pairwise
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated) & \
                   np.isfinite(observed) & np.isfinite(simulated)
    
    if np.sum(valid_indices) < 2:
        warnings.warn("Not enough valid pairs to calculate KGE after removing NaNs/Infs.")
        return np.nan
    
    obs_valid = observed[valid_indices]
    sim_valid = simulated[valid_indices]
    
    # Additional check for extreme values
    if np.any(np.abs(obs_valid) > 1e10) or np.any(np.abs(sim_valid) > 1e10):
        warnings.warn("Extreme values detected in KGE calculation. Clipping to reasonable range.")
        obs_valid = np.clip(obs_valid, -1e10, 1e10)
        sim_valid = np.clip(sim_valid, -1e10, 1e10)
    
    # Calculate components with float64 for better precision
    mean_obs = np.mean(obs_valid.astype(np.float64))
    mean_sim = np.mean(sim_valid.astype(np.float64))
    
    if np.isnan(mean_obs) or np.isnan(mean_sim) or np.isinf(mean_obs) or np.isinf(mean_sim):
        warnings.warn("Mean calculation resulted in NaN or Inf in KGE.")
        return np.nan
    
    std_obs = np.std(obs_valid.astype(np.float64))
    std_sim = np.std(sim_valid.astype(np.float64))
    
    # Avoid division by zero
    if std_obs < 1e-10 or mean_obs == 0:
        warnings.warn("Zero or near-zero standard deviation or mean observed value in KGE.")
        return np.nan
    
    # Pearson correlation coefficient
    if std_sim < 1e-10:
        r = 0.0  # No variation in simulated values
    else:
        correlation_matrix = np.corrcoef(obs_valid, sim_valid)
        r = correlation_matrix[0, 1]
        if np.isnan(r):
            r = 0.0
    
    # Alpha (ratio of standard deviations)
    alpha = std_sim / std_obs
    
    # Beta (ratio of means)
    beta = mean_sim / mean_obs
    
    # Calculate KGE
    kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    
    return kge


class KGELoss(nn.Module):
    """Differentiable KGE loss function optimized for normalized space training"""
    def __init__(self, epsilon=1e-6, beta_scale=2.0):
        super().__init__()
        self.epsilon = epsilon  
        self.beta_scale = beta_scale
        
    def soft_sign_transform(self, x):
        """Maps values to positive domain [0, 2] while preserving relative differences"""
        return 1.0 + torch.tanh(x / self.beta_scale)
        
    def forward(self, predictions, targets):
        # Flatten tensors
        predictions = predictions.view(-1)
        targets = targets.view(-1)
        
        # Calculate means
        mean_pred = torch.mean(predictions)
        mean_target = torch.mean(targets)
        
        # Calculate standard deviations
        std_pred = torch.std(predictions, unbiased=False)
        std_target = torch.std(targets, unbiased=False)
        
        # Correlation calculation
        if std_pred < self.epsilon or std_target < self.epsilon:
            r = torch.tensor(0.0, device=predictions.device)
        else:
            pred_standardized = (predictions - mean_pred) / (std_pred + self.epsilon)
            target_standardized = (targets - mean_target) / (std_target + self.epsilon)
            r = torch.mean(pred_standardized * target_standardized)
            r = torch.clamp(r, min=-0.999, max=0.999)
        
        # Beta calculation with soft sign transformation
        transformed_mean_pred = self.soft_sign_transform(mean_pred)
        transformed_mean_target = self.soft_sign_transform(mean_target)
        
        # Calculate beta as ratio of transformed values
        beta = transformed_mean_pred / (transformed_mean_target + self.epsilon)
        
        # Alpha calculation - ENSURE IT'S ALWAYS A TENSOR
        if std_target < self.epsilon:
            # Create tensor instead of using float
            if std_pred < self.epsilon:
                alpha = torch.tensor(1.0, device=predictions.device)
            else:
                alpha = torch.tensor(10.0, device=predictions.device)
        else:
            alpha = (std_pred + self.epsilon) / (std_target + self.epsilon)
        
        # Apply bounds - now alpha and beta are guaranteed to be tensors
        alpha = torch.clamp(alpha, min=0.01, max=100.0)
        beta = torch.clamp(beta, min=0.01, max=100.0)
        
        # Calculate KGE
        kge = 1 - torch.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
        kge = torch.clamp(kge, min=-10.0)
        
        return -kge


def calculate_kge_components(observed: np.ndarray, simulated: np.ndarray) -> tuple:
    """
    Calculate KGE and its three components.
    
    Args:
        observed: Array of observed values
        simulated: Array of simulated values
        
    Returns:
        tuple: (kge, r, alpha, beta)
    """
    # Ensure inputs are numpy arrays
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()
    
    # Remove NaN values pairwise
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated) & np.isfinite(observed) & np.isfinite(simulated)
    if np.sum(valid_indices) < 2:
        return np.nan, np.nan, np.nan, np.nan
    
    obs_valid = observed[valid_indices]
    sim_valid = simulated[valid_indices]
    
    # Calculate components
    mean_obs = np.mean(obs_valid)
    mean_sim = np.mean(sim_valid)
    std_obs = np.std(obs_valid)
    std_sim = np.std(sim_valid)
    
    # Avoid division by zero
    if std_obs < 1e-10 or mean_obs == 0:
        return np.nan, np.nan, np.nan, np.nan
    
    # Pearson correlation coefficient
    if std_sim < 1e-10:
        r = 0.0
    else:
        correlation_matrix = np.corrcoef(obs_valid, sim_valid)
        r = correlation_matrix[0, 1]
        if np.isnan(r):
            r = 0.0
    
    # Alpha (ratio of standard deviations)
    alpha = std_sim / std_obs
    
    # Beta (ratio of means)
    beta = mean_sim / mean_obs
    
    # Calculate KGE
    kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    
    return kge, r, alpha, beta


def calculate_rmse(observed: np.ndarray, simulated: np.ndarray) -> float:
    """Calculate Root Mean Squared Error"""
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()
    
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated)
    if np.sum(valid_indices) < 1:
        return np.nan
        
    return np.sqrt(np.mean((observed[valid_indices] - simulated[valid_indices])**2))


def calculate_mae(observed: np.ndarray, simulated: np.ndarray) -> float:
    """Calculate Mean Absolute Error"""
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()
    
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated)
    if np.sum(valid_indices) < 1:
        return np.nan
        
    return np.mean(np.abs(observed[valid_indices] - simulated[valid_indices]))




import traceback

def unnormalize_predictions(predictions, norm_params):
    """
    Reverse normalization for streamflow predictions with safety bounds.
    """
    # Clip predictions to reasonable range before unnormalization
    predictions = np.clip(predictions, -10, 10)
    
    # Reverse standardization
    log_predictions = predictions * norm_params['target_std'] + norm_params['target_mean']
    
    # Clip log predictions to prevent extreme exponentials
    # exp(20) ≈ 485 million, which should be more than enough for any streamflow
    log_predictions = np.clip(log_predictions, -10, 20)
    
    # Reverse log transform
    original_predictions = np.expm1(log_predictions)  # exp(x) - 1
    
    # Final safety clip for streamflow (mm/day)
    # Maximum reasonable streamflow is ~10,000 mm/day
    original_predictions = np.clip(original_predictions, 0, 10000)
    
    return original_predictions



def unnormalize_features(features_df, norm_params, feature_type='dynamic'):
    """
    Reverse normalization for features (useful for interpretation).
    """
    unnormalized_df = features_df.copy()
    
    if feature_type == 'dynamic':
        means = norm_params['dynamic_means']
        stds = norm_params['dynamic_stds']
    else:  # static
        means = norm_params['static_means']
        stds = norm_params['static_stds']
    
    for col in features_df.columns:
        if col in means:
            if stds[col] > 0:
                unnormalized_df[col] = features_df[col] * stds[col] + means[col]
    
    return unnormalized_df



#################################################################################################






def prepare_ptf_dataframe(
    watershed_subset_df,
    bucket_name,
    base_data_dir,
    base_attr_dir,
    include_static_attributes=True,
    norm_params=None,  
    compute_norm_params=False 
    ):
    """
    Loads data for multiple watersheds and prepares a single DataFrame suitable
    for pytorch-forecasting's TimeSeriesDataSet.

    Args:
        watershed_subset_df (pd.DataFrame): DataFrame with 'subdirectory_name' and 'watershedID'.
        bucket_name (str): S3 bucket name.
        base_data_dir (str): S3 base directory for timeseries NetCDF files.
        base_attr_dir (str): S3 base directory for attribute CSV files.
        include_static_attributes (bool): Whether to load and include static attributes.

    Returns:
        pd.DataFrame: A single DataFrame containing timeseries data, static attributes (if included),
                      a 'group_id' column, and an integer 'time_idx' column.
                      Returns an empty DataFrame if loading fails for all watersheds.
        list: List of column names identified as static real attributes.
        list: List of column names identified as time-varying real attributes (covariates).
    """
    # Define these at the start to avoid undefined variable errors if the loop is skipped
    weather_cols = [
#        'potential_evaporation_sum_ERA5_LAND', 
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean'
    ]
    # ADD: Include new smoothed precip features in weather columns
    precip_cols = ['Precip_smoothed_5day', 'Precip_lagged_90day']
    time_varying_cols = weather_cols + precip_cols  # Modified to include smoothed precip
    
    all_data_list = []
    static_attribute_names = []

    for _, row in watershed_subset_df.iterrows():
        subdir = row['subdirectory_name']
        ws_id = row['watershedID']
        group_id = f"{subdir}_{ws_id}"
        print(f"  Loading data for {group_id}...")

        try:
            # 1. Load Timeseries Data
            timeseries_df = load_and_process_watershed_data(
                bucket_name, base_data_dir, subdir, ws_id
            )
            if timeseries_df.empty:
                print(f"  Warning: No timeseries data loaded for {group_id}. Skipping.")
                continue
                
            weather_cols = [
#                'potential_evaporation_sum_ERA5_LAND',
                'surface_net_solar_radiation_mean',
                'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min',
                'total_precipitation_sum',
                'dewpoint_temperature_2m_mean'
            ]
            target_col = 'streamflow'
            
            # streamflow validataion (check for unrealistic values, or values of 0 that will be problematic for log transformation)
            invalid_mask = timeseries_df[target_col] <= 0
            if invalid_mask.any():
                n_invalid = invalid_mask.sum()
                print(f"WARNING: Found {n_invalid} non-positive streamflow values in {group_id}")
                timeseries_df.loc[invalid_mask, target_col] = 1e-6
            
            
            

            # Calculate precipitation features BEFORE subsetting columns
            if 'total_precipitation_sum' in timeseries_df.columns:
                # Calculate Precip_smoothed_5day (5-day centered smoothing)
                # Use center=True for centered window, min_periods=1 to handle edges
                timeseries_df['Precip_smoothed_5day'] = (
                    timeseries_df['total_precipitation_sum']
                    .rolling(window=5, min_periods=1, center=True)
                    .mean()
                )
                
                # Calculate Precip_lagged_90day (90-day lagged average)
                # This is the average precipitation over the preceding 90 days
                timeseries_df['Precip_lagged_90day'] = (
                    timeseries_df['total_precipitation_sum']
                    .rolling(window=90, min_periods=1, center=False)
                    .mean()
                )
                
                # Handle any potential NaN values at the beginning of the series
                # For smoothed: should have minimal NaNs due to min_periods=1
                timeseries_df['Precip_smoothed_5day'] = timeseries_df['Precip_smoothed_5day'].ffill().bfill()
                # For lagged: forward fill then backward fill
                timeseries_df['Precip_lagged_90day'] = timeseries_df['Precip_lagged_90day'].ffill().bfill()
                
                # Log extreme values if any
                for col in precip_cols:
                    col_max = timeseries_df[col].max()
                    col_min = timeseries_df[col].min()
                    col_mean = timeseries_df[col].mean()
                    if col_max > 1000 or col_min < 0:  # Precipitation shouldn't be negative or extremely high
                        print(f"  Warning: Extreme {col} values in {group_id}: min={col_min:.3f}, max={col_max:.3f}, mean={col_mean:.3f}")
            else:
                print(f"  Warning: Missing total_precipitation_sum column for precipitation features in {group_id}")
                # Create dummy columns with appropriate default values (using typical precipitation values)
                timeseries_df['Precip_smoothed_5day'] = 0  
                timeseries_df['Precip_lagged_90day'] = 0

            
            # Now update the columns to keep
            date_col = 'date'
            target_col = 'streamflow'
            keep_cols = [date_col, target_col] + weather_cols + precip_cols
            timeseries_df = timeseries_df[keep_cols]

            
            

            timeseries_df[date_col] = pd.to_datetime(timeseries_df[date_col])
            timeseries_df = timeseries_df.sort_values(by=date_col)

            # Create full date range and reindex BEFORE filling NaNs
            if not timeseries_df.empty:
                date_range = pd.date_range(start=timeseries_df[date_col].min(),
                                           end=timeseries_df[date_col].max(), freq='D')
                timeseries_df = timeseries_df.set_index(date_col).reindex(date_range)

                # --- Fill NaNs ---
                # Forward fill is often suitable for time series
                # Fill target first
                timeseries_df[target_col] = timeseries_df[target_col].ffill()
                # Fill covariates
                for col in weather_cols:
                    timeseries_df[col] = timeseries_df[col].ffill()
                # Optional: Backward fill any remaining NaNs at the beginning
                timeseries_df = timeseries_df.bfill()
                # Check for extreme values that could cause instability
                for col in weather_cols + [target_col]:
                    if timeseries_df[col].max() > 1e6 or timeseries_df[col].abs().mean() > 1e4:
                        warnings.warn(f"Potential outlier detected in {group_id} for column {col}. Max: {timeseries_df[col].max()}, Mean: {timeseries_df[col].mean()}")

                # Check if any NaNs remain after filling (shouldn't happen with ffill+bfill unless all were NaN)
                if timeseries_df[[target_col] + weather_cols].isnull().any().any():
                     print(f"  Warning: NaNs still present in timeseries for {group_id} after filling. Skipping.")
                     continue

                timeseries_df = timeseries_df.reset_index().rename(columns={'index': date_col})
            else:
                # Handle case where timeseries_df was empty after loading/trimming
                 print(f"  Warning: Empty timeseries after loading/trimming for {group_id}. Skipping.")
                 continue


            timeseries_df['time_idx'] = (timeseries_df[date_col] - timeseries_df[date_col].min()).dt.days
            timeseries_df['group_id'] = group_id

            # 2. Load and Merge Static Attributes
            if include_static_attributes:
                # (Static attribute loading logic remains the same)
                print(f"  Searching for attributes for gage_id: {group_id}")
                static_attrs = get_watershed_attributes(
                    bucket_name, base_attr_dir, subdir, ws_id
                )
                if not static_attrs.empty:
                    current_static_names = []
                    for col_name, value in static_attrs.items():
                        col_name_str = str(col_name)
                        try:
                            timeseries_df[col_name_str] = pd.to_numeric(value)
                            current_static_names.append(col_name_str)
                        except (ValueError, TypeError):
                             print(f"  Warning: Skipping non-numeric static attribute '{col_name_str}'.")
                    if not static_attribute_names:
                        static_attribute_names = current_static_names
                    elif set(static_attribute_names) != set(current_static_names):
                        warnings.warn(f"Inconsistent static attributes found for {group_id}. Using attributes from the first watershed.")
                        timeseries_df = timeseries_df.drop(columns=[c for c in current_static_names if c not in static_attribute_names])
                else:
                    print(f"  Warning: Could not load static attributes for {group_id}. Filling with NaN.")
                    for col_name_str in static_attribute_names:
                         timeseries_df[col_name_str] = np.nan
                 # Fill NaNs potentially introduced in static columns if some groups lacked them
#                for col in static_attribute_names:
#                     timeseries_df[col] = timeseries_df[col].ffill().bfill()
                for col in static_attribute_names:
                    if timeseries_df[col].isnull().any():
                        # Use a domain-appropriate default or median from all watersheds
                        default_val = 0.0
                        timeseries_df[col] = timeseries_df[col].fillna(default_val)

            all_data_list.append(timeseries_df)

        except Exception as e:
            print(f"  Error processing watershed {group_id}: {e}")
            traceback.print_exc()
            continue

    if not all_data_list:
        print("Error: No data loaded successfully.")
        return pd.DataFrame(), [], []

    combined_df = pd.concat(all_data_list, ignore_index=True)
    time_varying_cols = weather_cols


    # In prepare_ptf_dataframe, after loading data:
    # Check for extreme values
    for col in weather_cols + [target_col]:
        q99 = timeseries_df[col].quantile(0.99)
        q01 = timeseries_df[col].quantile(0.01)
        if q99 / (q01 + 1e-6) > 1000:  # Huge range
            print(f"WARNING: Extreme values in {col} for {group_id}")
            # Consider winsorizing or robust scaling


    
    # NORMALIZE FEATURES
     # Compute or apply normalization parameters
    if compute_norm_params:
        # Only compute params, don't normalize yet
        # Update time_varying_cols to include smoothed precip
        all_time_varying_cols = weather_cols + precip_cols
        
        norm_params = {
            'static_means': {col: combined_df[col].mean() for col in static_attribute_names},
            'static_stds': {col: combined_df[col].std() for col in static_attribute_names},
            'dynamic_means': {col: combined_df[col].mean() for col in all_time_varying_cols},
            'dynamic_stds': {col: combined_df[col].std() for col in all_time_varying_cols},
            'target_mean': np.mean(np.log1p(combined_df[target_col])),
            'target_std': np.std(np.log1p(combined_df[target_col]))
        }
    
    if norm_params is not None:
        # Apply normalization using provided parameters
        # Static features
        for col in static_attribute_names:
            if col in norm_params['static_means']:
                mean_val = norm_params['static_means'][col]
                std_val = norm_params['static_stds'][col]
                if std_val > 0:
                    combined_df[col] = (combined_df[col] - mean_val) / std_val
                else:
                    combined_df[col] = 0.0
        
        # Dynamic features (including smoothed precip)
        all_time_varying_cols = weather_cols + precip_cols
        for col in all_time_varying_cols:
            if col in norm_params['dynamic_means']:
                mean_val = norm_params['dynamic_means'][col]
                std_val = norm_params['dynamic_stds'][col]
                if std_val > 0:
                    combined_df[col] = (combined_df[col] - mean_val) / std_val
                else:
                    combined_df[col] = 0.0
        
        # Target (log-transform first, then normalize)
        combined_df[target_col] = np.log1p(combined_df[target_col])
        mean_target = norm_params['target_mean']
        std_target = norm_params['target_std']
        combined_df[target_col] = (combined_df[target_col] - mean_target) / std_target
        
        # ========== ADD VALIDATION DATA CHECKS HERE ==========
        # Check for extreme normalized values
        print("\nChecking for extreme normalized values...")
        extreme_found = False
        
        # Check target values
        target_normalized = combined_df[target_col]
        if target_normalized.max() > 20 or target_normalized.min() < -20:
            extreme_found = True
            print(f"WARNING: Extreme normalized target values detected!")
            print(f"  Min: {target_normalized.min():.2f}, Max: {target_normalized.max():.2f}")
            print(f"  Mean: {target_normalized.mean():.2f}, Std: {target_normalized.std():.2f}")
            # Clip extreme values
            combined_df[target_col] = combined_df[target_col].clip(-20, 20)
            print(f"  Clipped target values to [-20, 20]")
        
        # Check dynamic features
        for col in all_time_varying_cols:
            if col in combined_df.columns:
                col_values = combined_df[col]
                if col_values.max() > 20 or col_values.min() < -20:
                    extreme_found = True
                    print(f"WARNING: Extreme normalized values in {col}!")
                    print(f"  Min: {col_values.min():.2f}, Max: {col_values.max():.2f}")
                    combined_df[col] = combined_df[col].clip(-20, 20)
        
        # Check static features
        for col in static_attribute_names:
            if col in combined_df.columns:
                col_values = combined_df[col]
                if col_values.max() > 20 or col_values.min() < -20:
                    extreme_found = True
                    print(f"WARNING: Extreme normalized values in static feature {col}!")
                    print(f"  Min: {col_values.min():.2f}, Max: {col_values.max():.2f}")
                    combined_df[col] = combined_df[col].clip(-20, 20)
        
        if not extreme_found:
            print("  No extreme values found - normalization looks good!")
        
        # Additional check for NaN or Inf values - FIXED VERSION
        # Only check numeric columns
        numeric_cols = combined_df.select_dtypes(include=[np.number]).columns
        
        # Check for NaN values in numeric columns
        if combined_df[numeric_cols].isna().any().any():
            print(f"WARNING: NaN values detected after normalization!")
            # Fill NaN values in numeric columns only
            for col in numeric_cols:
                if combined_df[col].isna().any():
                    combined_df[col] = combined_df[col].fillna(0)
            print("  Filled NaN values with 0")
        
        # Check for Inf values in numeric columns
        inf_mask = np.isinf(combined_df[numeric_cols].values)
        if inf_mask.any():
            print(f"WARNING: Inf values detected after normalization!")
            # Replace inf values with large but finite values
            for col in numeric_cols:
                combined_df[col] = combined_df[col].replace([np.inf, -np.inf], [10, -10])
            print("  Replaced Inf values with ±10")
        # ========== END OF VALIDATION DATA CHECKS ==========
            
    # Ensure correct types
    combined_df[target_col] = combined_df[target_col].astype(float)
    for col in time_varying_cols:
        combined_df[col] = combined_df[col].astype(float)
    for col in static_attribute_names:
        combined_df[col] = combined_df[col].astype(float)


    time_varying_cols = weather_cols + precip_cols
    
    print(f"Successfully prepared combined DataFrame with shape: {combined_df.shape}")
    print(f"Identified static reals: {static_attribute_names}")
    print(f"Identified time-varying reals: {time_varying_cols}")

    return combined_df, static_attribute_names, time_varying_cols, norm_params




class WatershedDataset(Dataset):
    """PyTorch Dataset for watershed timeseries data"""
    def __init__(self, data_df, static_cols, dynamic_cols, target_col, sequence_length):
        self.data_df = data_df.copy()
        self.static_cols = static_cols
        self.dynamic_cols = dynamic_cols
        self.target_col = target_col
        self.sequence_length = sequence_length
        
        # Group data by watershed
        self.groups = self.data_df.groupby('group_id')
        self.group_keys = list(self.groups.groups.keys())
        
        # Calculate valid indices for each group
        self.valid_indices = []
        for group_id in self.group_keys:
            group_data = self.groups.get_group(group_id)
            n_samples = len(group_data) - self.sequence_length
            if n_samples > 0:
                self.valid_indices.extend([(group_id, i) for i in range(n_samples)])
    
    def __len__(self):
        return len(self.valid_indices)
    
    def __getitem__(self, idx):
        group_id, start_idx = self.valid_indices[idx]
        group_data = self.groups.get_group(group_id)
        
        # Extract sequence of dynamic features
        end_idx = start_idx + self.sequence_length
        dynamic_seq = group_data[self.dynamic_cols].iloc[start_idx:end_idx].values
        
        # Extract static features (same for all timesteps in watershed)
        static_features = group_data[self.static_cols].iloc[0].values
        
        # Extract target value (next timestep after sequence)
        target = group_data[self.target_col].iloc[end_idx:end_idx+1].values[0]
        
        return (torch.FloatTensor(dynamic_seq), 
                torch.FloatTensor(static_features), 
                torch.FloatTensor([target]))



class EntityAwareLSTM(nn.Module):
    """Entity-Aware LSTM with optional Time Modulation"""
    def __init__(self, dynamic_input_size, static_input_size, hidden_size=256, 
                 num_layers=1, dropout=0.2, use_time_modulation=True):
        super().__init__()
        
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.use_time_modulation = use_time_modulation
        
        # Standard LSTM with dynamic inputs
        self.lstm = nn.LSTM(
            input_size=dynamic_input_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            dropout=dropout if num_layers > 1 else 0
        )
        
        # Static feature processing for input gate modulation
        self.static_encoder = nn.Sequential(
            nn.Linear(static_input_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
            nn.Sigmoid()  # Output between 0 and 1 for gate modulation
        )
        
        # Optional: Additional static feature integration
        self.static_hidden_init = nn.Sequential(
            nn.Linear(static_input_size, hidden_size * num_layers),
            nn.GELU()
        )
        
        # Time Modulation Network (NEW)
        if self.use_time_modulation:
            self.time_modulator = nn.Sequential(
                nn.Conv1d(hidden_size, hidden_size, kernel_size=1),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Conv1d(hidden_size, hidden_size, kernel_size=1),
                nn.Sigmoid()  # Output between 0 and 1 for temporal gating
            )
        
        # Output layers
        self.output_layer = nn.Sequential(
            nn.Linear(hidden_size + static_input_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, 1)
        )
        
        self._initialize_weights()
        
    def forward(self, dynamic_inputs, static_inputs):
        batch_size = dynamic_inputs.size(0)
        
        # Initialize hidden state with static features
        h_0 = self.static_hidden_init(static_inputs)
        h_0 = h_0.view(self.num_layers, batch_size, self.hidden_size)
        c_0 = torch.zeros_like(h_0)
        
        # Generate input gate modulation from static features
        input_gate_modulation = self.static_encoder(static_inputs)
        
        # Run LSTM
        lstm_out, (h_n, c_n) = self.lstm(dynamic_inputs, (h_0, c_0))
        
        # Apply time modulation if enabled
        if self.use_time_modulation:
            # Reshape for Conv1d: (batch, channels, length)
            lstm_out_reshaped = lstm_out.transpose(1, 2)  # (batch, hidden_size, seq_len)
            
            # Apply temporal modulation
            time_modulation = self.time_modulator(lstm_out_reshaped)
            
            # Reshape back: (batch, seq_len, hidden_size)
            time_modulation = time_modulation.transpose(1, 2)
            
            # Combine static and temporal modulation
            # This allows the model to learn both entity-specific and time-varying patterns
            combined_modulation = input_gate_modulation.unsqueeze(1) * time_modulation
            lstm_out = lstm_out * combined_modulation
        else:
            # Original entity-aware modulation only
            lstm_out = lstm_out * input_gate_modulation.unsqueeze(1)
        
        # Take last timestep output
        last_output = lstm_out[:, -1, :]
        
        # Concatenate with static features for final prediction
        combined = torch.cat([last_output, static_inputs], dim=1)
        prediction = self.output_layer(combined)
        
        return prediction
    
    def _initialize_weights(self):  
        """Initialize weights with smaller values to prevent gradient explosion"""
        for name, param in self.named_parameters():
            if 'weight' in name:
                if 'lstm' in name:
                    # Use Xavier initialization for LSTM weights
                    nn.init.xavier_uniform_(param, gain=0.5)  # Small gain
                elif 'conv' in name:
                    # Use He initialization for convolutional layers
                    nn.init.kaiming_uniform_(param, a=0, mode='fan_in', nonlinearity='leaky_relu')
                    param.data *= 0.5  # Scale down for stability
                else:
                    # Use He initialization for other layers
                    nn.init.kaiming_uniform_(param, a=0, mode='fan_in', nonlinearity='leaky_relu')
                    param.data *= 0.1  # Scale down
            elif 'bias' in name:
                nn.init.constant_(param, 0.0)






class StreamflowLightningModule(pl.LightningModule):
    """Lightning wrapper for EA-LSTM streamflow model with KGE loss"""
    
    def __init__(
        self,
        dynamic_input_size: int,
        static_input_size: int,
        hidden_size: int = 256,
        num_layers: int = 1,
        dropout: float = 0.2,
        learning_rate: float = 1e-4,
        use_time_modulation: bool = True,
        # Training parameters
        target_noise_std: float = 0.1,
        target_noise_type: str = 'log-normal',
        gradient_clip_val: float = 5.0,
        # LR scheduler parameters
        lr_scheduler_params: Optional[Dict[str, Any]] = None,
        # Normalization parameters for metrics
        norm_params: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        super().__init__()
        
        # Save all hyperparameters for checkpointing and logging
        self.save_hyperparameters()
        
        # Initialize the EA-LSTM model
        self.model = EntityAwareLSTM(
            dynamic_input_size=dynamic_input_size,
            static_input_size=static_input_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            dropout=dropout,
            use_time_modulation=use_time_modulation
        )
        
        # Loss functions
        self.criterion = KGELoss()  # Base KGE loss
        self.mse_criterion = nn.MSELoss()  # MSE loss for additional term
        
        # Target noise wrapper
        self.target_noise = TargetNoiseWrapper(
            noise_std=target_noise_std,
            noise_type=target_noise_type
        )
        
        # Store parameters
        self.learning_rate = learning_rate
        self.gradient_clip_val = gradient_clip_val
        self.lr_scheduler_params = lr_scheduler_params or {
            'factor': 0.5,
            'patience': 2,
            'min_lr': 1e-7,
            'threshold': 0.001,
            'cooldown': 1,
            'verbose': True
        }
        self.norm_params = norm_params
        
        # Track validation outputs for epoch-level metrics
        self.validation_step_outputs = []
        self.test_step_outputs = []

    def forward(self, dynamic_inputs, static_inputs):
        """Forward pass through the model"""
        predictions = self.model(dynamic_inputs, static_inputs)
        # Clip predictions to prevent extreme values during any forward pass
        predictions = torch.clamp(predictions, min=-10, max=15)
        return predictions
        
    def training_step_og(self, batch, batch_idx):
        """Training step with progressive KGE + MSE loss"""
        dynamic_seq, static_feat, target = batch
        
        # Forward pass
        predictions = self(dynamic_seq, static_feat)
        
        # Add noise to target during training
        noisy_target = self.target_noise.add_noise(target, training=True)
        
        # Calculate base losses
        kge_loss = self.criterion(predictions, noisy_target)
        mse_loss = self.mse_criterion(predictions, noisy_target)
        
        # Calculate MSE weight: ramp from 0 to 0.5 over first 50 epochs
        mse_weight = min(0.5, self.current_epoch / 50.0)
        
        # Combined loss: (1-weight)*KGE + weight*MSE
        loss = (1 - mse_weight) * kge_loss + mse_weight * mse_loss
        
        # Check for NaN/Inf
        if torch.isnan(loss) or torch.isinf(loss):
            self.log('train_nan_count', 1.0, on_step=True, prog_bar=False)
            return None  # Skip this batch
        
        # Log metrics
        self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_kge_loss', kge_loss, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('train_mse_loss', mse_loss, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('mse_weight', mse_weight, on_step=False, on_epoch=True, prog_bar=False, logger=True)
        
        # Log gradient norm if available
        if self.global_step > 0 and self.global_step % 50 == 0:
            total_norm = 0
            for p in self.model.parameters():
                if p.grad is not None:
                    param_norm = p.grad.data.norm(2)
                    total_norm += param_norm.item() ** 2
            total_norm = total_norm ** 0.5
            self.log('grad_norm', total_norm, on_step=True, prog_bar=False)
        
        return loss
        
    def training_step(self, batch, batch_idx):
        """Training step with scale-balanced KGE + MSE loss"""
        dynamic_seq, static_feat, target = batch
        
        # Forward pass
        predictions = self(dynamic_seq, static_feat)
        
        # Add noise to target during training
        noisy_target = self.target_noise.add_noise(target, training=True)
        
        # Calculate base losses
        kge_loss = self.criterion(predictions, noisy_target)
        mse_loss = self.mse_criterion(predictions, noisy_target)
        
        # Normalize MSE by target variance with safety checks
        target_var = torch.var(noisy_target)
        
        # Add safety check for very small variance
        if target_var < 1e-6:
            # If variance is too small, use standard deviation instead
            target_std = torch.std(noisy_target) + 1e-6
            normalized_mse = mse_loss / (target_std ** 2)
        else:
            normalized_mse = mse_loss / target_var
        
        # Additional safety: clamp normalized MSE to reasonable range
        normalized_mse = torch.clamp(normalized_mse, min=0.0, max=10.0)
        
        # Calculate MSE weight: ramp from 0 to 0.5 over first 1000 epochs
        mse_weight = min(0.5, self.current_epoch / 1000.0)
        
        # Combined loss
        loss = (1 - mse_weight) * kge_loss + mse_weight * normalized_mse
        
        # Check for NaN/Inf
        if torch.isnan(loss) or torch.isinf(loss):
            self.log('train_nan_count', 1.0, on_step=True, prog_bar=False)
            # Log debug info
            print(f"NaN/Inf in loss. KGE: {kge_loss.item()}, MSE: {mse_loss.item()}, "
                  f"Var: {target_var.item()}, Normalized MSE: {normalized_mse.item()}")
            return None  # Skip this batch
        
        # Log metrics
        self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_kge_loss', kge_loss, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('train_mse_loss', mse_loss, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('train_normalized_mse', normalized_mse, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('train_target_var', target_var, on_step=True, on_epoch=True, prog_bar=False, logger=True)
        self.log('mse_weight', mse_weight, on_step=False, on_epoch=True, prog_bar=False, logger=True)
        
        # Log gradient norm if available
        if self.global_step > 0 and self.global_step % 50 == 0:
            total_norm = 0
            for p in self.model.parameters():
                if p.grad is not None:
                    param_norm = p.grad.data.norm(2)
                    total_norm += param_norm.item() ** 2
            total_norm = total_norm ** 0.5
            self.log('grad_norm', total_norm, on_step=True, prog_bar=False)
        
        return loss
    
    def validation_step(self, batch, batch_idx):
        """Validation step with normalized MSE"""
        dynamic_seq, static_feat, target = batch
        
        # Forward pass
        predictions = self(dynamic_seq, static_feat)
        
        # Calculate losses
        kge_loss = self.criterion(predictions, target)
        mse_loss = self.mse_criterion(predictions, target)
        
        # Normalize MSE by target variance (same as training)
        target_var = torch.var(target)
        if target_var < 1e-6:
            target_std = torch.std(target) + 1e-6
            normalized_mse = mse_loss / (target_std ** 2)
        else:
            normalized_mse = mse_loss / target_var
        
        normalized_mse = torch.clamp(normalized_mse, min=0.0, max=10.0)
        
        # Calculate MSE weight: ramp from 0 to 0.5 over first 1000 epochs
        mse_weight = min(0.5, self.current_epoch / 1000.0)
        loss = (1 - mse_weight) * kge_loss + mse_weight * normalized_mse
        
        # Store outputs for epoch-level metrics calculation
        self.validation_step_outputs.append({
            'loss': loss,
            'kge_loss': kge_loss,
            'mse_loss': mse_loss,
            'normalized_mse': normalized_mse,
            'predictions': predictions.detach().cpu(),
            'targets': target.detach().cpu()
        })
        
        return loss


    def on_validation_epoch_end(self):
        """Calculate comprehensive epoch-level validation metrics"""
        # Gather all outputs
        all_losses = [x['loss'] for x in self.validation_step_outputs]
        all_kge_losses = [x['kge_loss'] for x in self.validation_step_outputs]
        all_mse_losses = [x['mse_loss'] for x in self.validation_step_outputs]
        all_preds = torch.cat([x['predictions'] for x in self.validation_step_outputs])
        all_targets = torch.cat([x['targets'] for x in self.validation_step_outputs])
        
        # Calculate mean losses
        avg_loss = torch.stack(all_losses).mean()
        avg_kge_loss = torch.stack(all_kge_losses).mean()
        avg_mse_loss = torch.stack(all_mse_losses).mean()
        
        # Convert to numpy for metric calculations
        preds_np = all_preds.numpy().flatten()
        targets_np = all_targets.numpy().flatten()
        
        # Calculate metrics on original scale
        if self.norm_params is not None:
            # Unnormalize predictions and targets
            preds_original = unnormalize_predictions(preds_np, self.norm_params)
            targets_original = unnormalize_predictions(targets_np, self.norm_params)
        else:
            preds_original = preds_np
            targets_original = targets_np
        
        # Calculate all metrics
        kge, r, alpha, beta = calculate_kge_components(targets_original, preds_original)
        rmse = calculate_rmse(targets_original, preds_original)
        mae = calculate_mae(targets_original, preds_original)
        
        # Log all metrics
        self.log('val_loss', avg_loss, on_epoch=True, prog_bar=True, logger=True)
        self.log('val_kge_loss', avg_kge_loss, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_mse_loss', avg_mse_loss, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_kge', kge, on_epoch=True, prog_bar=True, logger=True)
        self.log('val_kge_r', r, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_kge_alpha', alpha, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_kge_beta', beta, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_rmse', rmse, on_epoch=True, prog_bar=False, logger=True)
        self.log('val_mae', mae, on_epoch=True, prog_bar=False, logger=True)
        
        # Clear stored outputs
        self.validation_step_outputs.clear()
        
    def test_step(self, batch, batch_idx):
        """Test step with comprehensive metrics"""
        dynamic_seq, static_feat, target = batch
        
        # Forward pass
        predictions = self(dynamic_seq, static_feat)
        
        # Calculate losses
        kge_loss = self.criterion(predictions, target)
        mse_loss = self.mse_criterion(predictions, target)
        
        # Use same weighting as training
        mse_weight = min(0.5, self.current_epoch / 50.0)
        loss = (1 - mse_weight) * kge_loss + mse_weight * mse_loss
        
        # Store outputs for metrics calculation
        self.test_step_outputs.append({
            'loss': loss,
            'kge_loss': kge_loss,
            'mse_loss': mse_loss,
            'predictions': predictions.detach().cpu(),
            'targets': target.detach().cpu()
        })
        
        return loss
    
    def on_test_epoch_end(self):
        """Calculate comprehensive test metrics"""
        # Gather all outputs
        all_losses = [x['loss'] for x in self.test_step_outputs]
        all_kge_losses = [x['kge_loss'] for x in self.test_step_outputs]
        all_mse_losses = [x['mse_loss'] for x in self.test_step_outputs]
        all_preds = torch.cat([x['predictions'] for x in self.test_step_outputs])
        all_targets = torch.cat([x['targets'] for x in self.test_step_outputs])
        
        # Calculate mean losses
        avg_loss = torch.stack(all_losses).mean()
        avg_kge_loss = torch.stack(all_kge_losses).mean()
        avg_mse_loss = torch.stack(all_mse_losses).mean()
        
        # Convert to numpy
        preds_np = all_preds.numpy().flatten()
        targets_np = all_targets.numpy().flatten()
        
        # Calculate metrics on original scale
        if self.norm_params is not None:
            preds_original = unnormalize_predictions(preds_np, self.norm_params)
            targets_original = unnormalize_predictions(targets_np, self.norm_params)
        else:
            preds_original = preds_np
            targets_original = targets_np
        
        # Calculate all metrics
        kge, r, alpha, beta = calculate_kge_components(targets_original, preds_original)
        rmse = calculate_rmse(targets_original, preds_original)
        mae = calculate_mae(targets_original, preds_original)
        
        # Log all metrics
        self.log('test_loss', avg_loss, on_epoch=True, logger=True)
        self.log('test_kge_loss', avg_kge_loss, on_epoch=True, logger=True)
        self.log('test_mse_loss', avg_mse_loss, on_epoch=True, logger=True)
        self.log('test_kge', kge, on_epoch=True, logger=True)
        self.log('test_kge_r', r, on_epoch=True, logger=True)
        self.log('test_kge_alpha', alpha, on_epoch=True, logger=True)
        self.log('test_kge_beta', beta, on_epoch=True, logger=True)
        self.log('test_rmse', rmse, on_epoch=True, logger=True)
        self.log('test_mae', mae, on_epoch=True, logger=True)
        
        # Clear stored outputs
        self.test_step_outputs.clear()
    
    def configure_optimizers(self):
        """Configure optimizer and learning rate scheduler"""
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        
        # Configure scheduler - now monitoring combined loss
        scheduler = ReduceLROnPlateau(
            optimizer,
            mode='min',
            factor=self.lr_scheduler_params['factor'],
            patience=self.lr_scheduler_params['patience'],
            min_lr=self.lr_scheduler_params['min_lr'],
            threshold=self.lr_scheduler_params['threshold'],
            cooldown=self.lr_scheduler_params['cooldown'],
            verbose=self.lr_scheduler_params['verbose']
        )
        
        return {
            'optimizer': optimizer,
            'lr_scheduler': {
                'scheduler': scheduler,
                'monitor': 'val_loss',
                'interval': 'epoch',
                'frequency': 1
            }
        }
    
    def on_after_backward(self):
        """Custom gradient clipping and monitoring"""
        if self.global_step % 100 == 0:
            total_norm = 0
            max_grad = 0
            for p in self.model.parameters():
                if p.grad is not None:
                    param_norm = p.grad.data.norm(2)
                    total_norm += param_norm.item() ** 2
                    max_grad = max(max_grad, p.grad.data.abs().max().item())
            total_norm = total_norm ** 0.5
            
            # Log large gradient warnings
            if total_norm > 100 or max_grad > 50:
                self.log('large_gradient_warning', 1.0, on_step=True)
                print(f"Warning: Large gradients detected - norm: {total_norm:.2f}, max: {max_grad:.2f}")





class CaravanDataModule(pl.LightningDataModule):
    """Lightning DataModule for Caravan watershed data with data caching for improved GPU utilization"""
    
    def __init__(
        self,
        watersheds_df: pd.DataFrame,
        bucket_name: str,
        base_data_dir: str,
        base_attr_dir: str,
        sequence_length: int = 365,
        batch_size: int = 256,
        num_workers: int = 16,  # Increased from 4
        chunk_size: int = 50,
        train_split: float = 0.6,
        val_split: float = 0.2,
        random_seed: int = 42,
        pin_memory: bool = True,
        persistent_workers: bool = True,
        prefetch_factor: int = 4,  # Increased from 2
        norm_params: Optional[Dict[str, Any]] = None,
        cache_data: bool = True,  # NEW: Enable data caching
        cache_dir: str = '/opt/ml/input/data/cache'  # NEW: Cache directory
    ):
        super().__init__()
        
        # Save parameters
        self.watersheds_df = watersheds_df
        self.bucket_name = bucket_name
        self.base_data_dir = base_data_dir
        self.base_attr_dir = base_attr_dir
        self.sequence_length = sequence_length
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.chunk_size = chunk_size
        self.train_split = train_split
        self.val_split = val_split
        self.random_seed = random_seed
        self.pin_memory = pin_memory
        self.persistent_workers = persistent_workers
        self.prefetch_factor = prefetch_factor
        self.cache_data = cache_data
        self.cache_dir = cache_dir
        
        # Data attributes (will be set during setup)
        self.train_watersheds = None
        self.val_watersheds = None
        self.test_watersheds = None
        self.norm_params = norm_params
        self.static_cols = None
        self.dynamic_cols_no_target = None
        
        # NEW: Cached data storage
        self._cached_data = {}
        self._cached_datasets = {}
        
        # Current chunk tracking for memory-efficient loading (kept for fallback)
        self.current_train_chunk = 0
        self.train_chunks = []
        
        # Create cache directory if caching is enabled
        if self.cache_data and self.cache_dir:
            os.makedirs(self.cache_dir, exist_ok=True)
    
    def setup(self, stage: Optional[str] = None):
        """Prepare data splits and compute normalization parameters"""
        
        if stage == 'fit' or stage is None:
            # Split watersheds
            np.random.seed(self.random_seed)
            train_watersheds, temp_watersheds = train_test_split(
                self.watersheds_df, test_size=(1 - self.train_split), random_state=self.random_seed
            )
            val_watersheds, test_watersheds = train_test_split(
                temp_watersheds, test_size=0.5, random_state=self.random_seed
            )
            
            self.train_watersheds = train_watersheds
            self.val_watersheds = val_watersheds
            self.test_watersheds = test_watersheds
            
            print(f"Train watersheds: {len(self.train_watersheds)}")
            print(f"Val watersheds: {len(self.val_watersheds)}")
            print(f"Test watersheds: {len(self.test_watersheds)}")
            
            # Compute normalization parameters if not provided
            if self.norm_params is None:
                print("Computing normalization parameters from training data...")
                _, _, _, self.norm_params = prepare_ptf_dataframe(
                    self.train_watersheds,
                    self.bucket_name,
                    self.base_data_dir,
                    self.base_attr_dir,
                    compute_norm_params=True
                )
            
            # Get feature columns from a sample
            sample_df, static_cols, dynamic_cols, _ = prepare_ptf_dataframe(
                self.train_watersheds.iloc[:5],
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params
            )
            
            self.static_cols = static_cols
            self.dynamic_cols_no_target = [col for col in dynamic_cols if col != 'streamflow']
            
            # NEW: Pre-cache all data if enabled
            if self.cache_data and stage == 'fit':
                self._cache_all_data()
            else:
                # Original chunking approach as fallback
                self.train_chunks = [
                    self.train_watersheds[i:i+self.chunk_size]
                    for i in range(0, len(self.train_watersheds), self.chunk_size)
                ]
            
            del sample_df
            gc.collect()
            
        if stage == 'test':
            # For testing, we already have test_watersheds from fit stage
            pass
    
    def _cache_all_data(self):
        """Pre-cache all training and validation data for faster loading"""
        print("\n" + "="*60)
        print("PRE-CACHING ALL DATA FOR FASTER TRAINING")
        print("="*60)
        
        # Try to load from disk cache first
        cache_file_train = os.path.join(self.cache_dir, 'train_data.pkl') if self.cache_dir else None
        cache_file_val = os.path.join(self.cache_dir, 'val_data.pkl') if self.cache_dir else None
        
        # Cache training data
        if cache_file_train and os.path.exists(cache_file_train):
            print("Loading cached training data from disk...")
            try:
                with open(cache_file_train, 'rb') as f:
                    self._cached_data['train'] = pickle.load(f)
                print(f"✓ Loaded training data from cache: {cache_file_train}")
            except Exception as e:
                print(f"✗ Failed to load cache: {e}")
                self._cached_data['train'] = None
        else:
            print(f"Loading {len(self.train_watersheds)} training watersheds from S3...")
            start_time = time.time()
            
            self._cached_data['train'], _, _, _ = prepare_ptf_dataframe(
                self.train_watersheds,
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params
            )
            
            load_time = time.time() - start_time
            print(f"✓ Training data loaded in {load_time:.1f} seconds")
            
            # Save to disk cache if path provided
            if cache_file_train:
                try:
                    print(f"Saving training data to cache: {cache_file_train}")
                    with open(cache_file_train, 'wb') as f:
                        pickle.dump(self._cached_data['train'], f, protocol=pickle.HIGHEST_PROTOCOL)
                    print("✓ Training data cached to disk")
                except Exception as e:
                    print(f"✗ Failed to save cache: {e}")
        
        # Cache validation data
        if cache_file_val and os.path.exists(cache_file_val):
            print("\nLoading cached validation data from disk...")
            try:
                with open(cache_file_val, 'rb') as f:
                    self._cached_data['val'] = pickle.load(f)
                print(f"✓ Loaded validation data from cache: {cache_file_val}")
            except Exception as e:
                print(f"✗ Failed to load cache: {e}")
                self._cached_data['val'] = None
        else:
            print(f"\nLoading {len(self.val_watersheds)} validation watersheds from S3...")
            start_time = time.time()
            
            self._cached_data['val'], _, _, _ = prepare_ptf_dataframe(
                self.val_watersheds,
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params
            )
            
            load_time = time.time() - start_time
            print(f"✓ Validation data loaded in {load_time:.1f} seconds")
            
            # Save to disk cache if path provided
            if cache_file_val:
                try:
                    print(f"Saving validation data to cache: {cache_file_val}")
                    with open(cache_file_val, 'wb') as f:
                        pickle.dump(self._cached_data['val'], f, protocol=pickle.HIGHEST_PROTOCOL)
                    print("✓ Validation data cached to disk")
                except Exception as e:
                    print(f"✗ Failed to save cache: {e}")
        
        # Pre-create datasets for even faster access
        print("\nPre-creating PyTorch datasets...")
        
        if 'train' in self._cached_data and self._cached_data['train'] is not None:
            self._cached_datasets['train'] = OptimizedWatershedDataset(
                self._cached_data['train'],
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length
            )
            print(f"✓ Training dataset created with {len(self._cached_datasets['train'])} samples")
        
        if 'val' in self._cached_data and self._cached_data['val'] is not None:
            self._cached_datasets['val'] = OptimizedWatershedDataset(
                self._cached_data['val'],
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length
            )
            print(f"✓ Validation dataset created with {len(self._cached_datasets['val'])} samples")
        
        # Memory report
        if 'train' in self._cached_data and self._cached_data['train'] is not None:
            train_memory = self._cached_data['train'].memory_usage(deep=True).sum() / 1e9
            print(f"\nMemory usage - Training data: {train_memory:.2f} GB")
        if 'val' in self._cached_data and self._cached_data['val'] is not None:
            val_memory = self._cached_data['val'].memory_usage(deep=True).sum() / 1e9
            print(f"Memory usage - Validation data: {val_memory:.2f} GB")
        
        print("="*60 + "\n")
    
    def train_dataloader(self):
        """Create training dataloader using cached data when available"""
        if self.cache_data and 'train' in self._cached_datasets:
            # Use pre-created cached dataset
            print("Using cached training dataset")
            dataset = self._cached_datasets['train']
        elif self.cache_data and 'train' in self._cached_data and self._cached_data['train'] is not None:
            # Create dataset from cached data
            dataset = OptimizedWatershedDataset(
                self._cached_data['train'],
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length
            )
        else:
            # Fallback to original chunking approach
            print("Warning: Using chunked loading (slower). Consider enabling caching.")
            chunk_watersheds = self.train_chunks[self.current_train_chunk]
            self.current_train_chunk = (self.current_train_chunk + 1) % len(self.train_chunks)
            
            chunk_df, _, _, _ = prepare_ptf_dataframe(
                chunk_watersheds,
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params
            )
            
            dataset = WatershedDataset(
                chunk_df,
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length
            )
        
        # Create dataloader with optimized settings
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory and torch.cuda.is_available(),
            persistent_workers=self.persistent_workers and self.num_workers > 0,
            prefetch_factor=self.prefetch_factor if self.num_workers > 0 else 2,
            drop_last=True
        )
    
    def val_dataloader(self):
        """Create validation dataloader using cached data when available"""
        if self.cache_data and 'val' in self._cached_datasets:
            # Use pre-created cached dataset
            print("Using cached validation dataset")
            dataset = self._cached_datasets['val']
            
            return DataLoader(
                dataset,
                batch_size=self.batch_size,
                shuffle=False,
                num_workers=min(self.num_workers, 4),  # Fewer workers for validation
                pin_memory=self.pin_memory and torch.cuda.is_available(),
                persistent_workers=self.persistent_workers and self.num_workers > 0,
                prefetch_factor=self.prefetch_factor if self.num_workers > 0 else 2
            )
        else:
            # Fallback to original approach
            return self._create_eval_dataloader(self.val_watersheds)
    
    def test_dataloader(self):
        """Create test dataloader"""
        # Test data is not cached by default to save memory
        return self._create_eval_dataloader(self.test_watersheds)
    
    def _create_eval_dataloader(self, watersheds_df):
        """Helper to create evaluation dataloaders with chunking (fallback method)"""
        chunks = [
            watersheds_df[i:i+self.chunk_size]
            for i in range(0, len(watersheds_df), self.chunk_size)
        ]
        
        dataloaders = []
        for chunk in chunks:
            chunk_df, _, _, _ = prepare_ptf_dataframe(
                chunk,
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params
            )
            
            if not chunk_df.empty:
                dataset = WatershedDataset(
                    chunk_df,
                    self.static_cols,
                    self.dynamic_cols_no_target,
                    'streamflow',
                    self.sequence_length
                )
                
                if len(dataset) > 0:
                    dataloader = DataLoader(
                        dataset,
                        batch_size=self.batch_size,
                        shuffle=False,
                        num_workers=min(self.num_workers, 4),
                        pin_memory=self.pin_memory and torch.cuda.is_available(),
                        persistent_workers=self.persistent_workers and self.num_workers > 0,
                        prefetch_factor=self.prefetch_factor if self.num_workers > 0 else 2
                    )
                    dataloaders.append(dataloader)
        
        # Return combined dataloader
        return CombinedDataLoader(dataloaders)
    
    def on_after_batch_transfer(self, batch, dataloader_idx):
        """Clean up memory after batch transfer"""
        if hasattr(self, '_cleanup_counter'):
            self._cleanup_counter += 1
        else:
            self._cleanup_counter = 0
            
        if self._cleanup_counter % 100 == 0:
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        
        return batch


class OptimizedWatershedDataset(Dataset):
    """Optimized PyTorch Dataset with pre-computed arrays for faster data loading"""
    
    def __init__(self, data_df, static_cols, dynamic_cols, target_col, sequence_length):
        self.data_df = data_df.copy()
        self.static_cols = static_cols
        self.dynamic_cols = dynamic_cols
        self.target_col = target_col
        self.sequence_length = sequence_length
        
        # Group data by watershed
        self.groups = self.data_df.groupby('group_id')
        self.group_keys = list(self.groups.groups.keys())
        
        # Pre-compute all valid samples for faster access
        print("Pre-computing dataset arrays for faster access...")
        self._precompute_all_samples()
        print(f"Dataset ready with {len(self)} samples")
    
    def _precompute_all_samples(self):
        """Pre-compute and cache all samples as numpy arrays"""
        self.all_samples = []
        self.valid_indices = []
        
        # Pre-convert all groups to numpy arrays
        self.group_arrays = {}
        
        for group_id in tqdm(self.group_keys, desc="Pre-processing watersheds"):
            group_data = self.groups.get_group(group_id)
            
            # Convert to numpy arrays once (much faster than repeated pandas indexing)
            dynamic_array = group_data[self.dynamic_cols].values.astype(np.float32)
            static_array = group_data[self.static_cols].iloc[0].values.astype(np.float32)
            target_array = group_data[self.target_col].values.astype(np.float32)
            
            self.group_arrays[group_id] = {
                'dynamic': dynamic_array,
                'static': static_array,
                'target': target_array
            }
            
            # Calculate valid indices for this group
            n_samples = len(group_data) - self.sequence_length
            if n_samples > 0:
                for i in range(n_samples):
                    self.valid_indices.append((group_id, i))
                    
                    # Pre-extract sequences
                    end_idx = i + self.sequence_length
                    dynamic_seq = dynamic_array[i:end_idx]
                    target = target_array[end_idx]
                    
                    # Store pre-extracted samples
                    self.all_samples.append((
                        dynamic_seq,
                        static_array,
                        np.array([target], dtype=np.float32)
                    ))
    
    def __len__(self):
        return len(self.all_samples)
    
    def __getitem__(self, idx):
        # Direct access to pre-computed samples (very fast)
        dynamic_seq, static_features, target = self.all_samples[idx]
        
        # Convert to tensors
        return (
            torch.from_numpy(dynamic_seq),
            torch.from_numpy(static_features),
            torch.from_numpy(target)
        )


class CombinedDataLoader:
    """Helper class to combine multiple dataloaders for chunked evaluation"""
    
    def __init__(self, dataloaders):
        self.dataloaders = dataloaders
        
    def __iter__(self):
        for dataloader in self.dataloaders:
            for batch in dataloader:
                yield batch
                
    def __len__(self):
        return sum(len(dl) for dl in self.dataloaders)






class ValidationFrequencyCallback(pl.Callback):
    def __init__(self, initial_frequency=1, later_frequency=10, switch_epoch=10):
        self.initial_frequency = initial_frequency
        self.later_frequency = later_frequency
        self.switch_epoch = switch_epoch
        
    def on_train_epoch_start(self, trainer, pl_module):
        if trainer.current_epoch == self.switch_epoch:
            trainer.check_val_every_n_epoch = self.later_frequency
            print(f"Switching validation frequency to every {self.later_frequency} epochs")




class TargetNoiseWrapper:
    """Wrapper to add noise to targets during training"""
    def __init__(self, noise_std=0.1, noise_type='multiplicative'):
        """
        Args:
            noise_std: Standard deviation of noise
            noise_type: 'additive' or 'multiplicative'
        """
        self.noise_std = noise_std
        self.noise_type = noise_type
        
    def add_noise(self, target, training=True):
        """Add noise to target values during training"""
        if not training or self.noise_std == 0:
            return target
            
        if self.noise_type == 'multiplicative':
            # Multiplicative noise: noise scales with target magnitude
            # Good for streamflow where variance often scales with mean
            noise = torch.randn_like(target) * self.noise_std
            noisy_target = target * (1 + noise)
            # Ensure non-negative streamflow
            noisy_target = torch.clamp(noisy_target, min=0)
            
        elif self.noise_type == 'additive':
            # Additive noise: constant variance
            noise = torch.randn_like(target) * self.noise_std * target.mean()
            noisy_target = target + noise
            noisy_target = torch.clamp(noisy_target, min=0)
            
        elif self.noise_type == 'log-normal':
            # Log-normal noise: better for skewed distributions like streamflow
            log_target = torch.log(target + 1e-6)  # Avoid log(0)
            noise = torch.randn_like(log_target) * self.noise_std
            noisy_target = torch.exp(log_target + noise) - 1e-6
            noisy_target = torch.clamp(noisy_target, min=0)  # 
            
        return noisy_target









#####################################################################################################################
#####################################################################################################################



def plot_training_history_with_lr(history, save_dir):
    """Plot training history including learning rate changes"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    epochs = history['epoch']
    
    # Loss plot
    ax1.plot(epochs, history['val_loss'], 'b-', label='Validation Loss')
    if history['train_loss'] and not all(np.isnan(history['train_loss'])):
        ax1.plot(epochs, history['train_loss'], 'r--', label='Train Loss')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')
    ax1.set_title('Training History - Loss')
    ax1.legend()
    ax1.grid(True)
    
    # KGE plot
    ax2.plot(epochs, history['val_kge'], 'g-', label='Validation KGE')
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('KGE')
    ax2.set_title('Training History - KGE')
    ax2.legend()
    ax2.grid(True)
    
    # Learning rate plot
    ax3.plot(epochs, history['learning_rate'], 'orange', marker='o')
    ax3.set_xlabel('Epoch')
    ax3.set_ylabel('Learning Rate')
    ax3.set_title('Learning Rate Schedule')
    ax3.set_yscale('log')
    ax3.grid(True)
    
    # Combined normalized plot
    if len(epochs) > 0:
        val_loss_norm = (history['val_loss'] - np.min(history['val_loss'])) / (np.max(history['val_loss']) - np.min(history['val_loss']))
        val_kge_norm = (history['val_kge'] - np.min(history['val_kge'])) / (np.max(history['val_kge']) - np.min(history['val_kge']))
        lr_norm = (history['learning_rate'] - np.min(history['learning_rate'])) / (np.max(history['learning_rate']) - np.min(history['learning_rate']))
        
        ax4.plot(epochs, val_loss_norm, 'b-', label='Val Loss (norm)')
        ax4.plot(epochs, val_kge_norm, 'g-', label='Val KGE (norm)')
        ax4.plot(epochs, lr_norm, 'orange', label='LR (norm)')
        ax4.set_xlabel('Epoch')
        ax4.set_ylabel('Normalized Value')
        ax4.set_title('All Metrics (Normalized)')
        ax4.legend()
        ax4.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'training_history_with_lr.png'), dpi=150)
    plt.close()












################################################################


    
def tensor_to_python(obj):
    """Recursively convert PyTorch tensors to Python native types"""
    if torch.is_tensor(obj):
        if obj.dim() == 0:  # scalar tensor
            return obj.item()
        else:
            return obj.tolist()
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: tensor_to_python(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [tensor_to_python(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(tensor_to_python(v) for v in obj)
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    else:
        return obj







def main():
    """Main entry point for SageMaker training with PyTorch Lightning"""
    
    # CUDA Memory Management
    if torch.cuda.is_available():
        # Clear any existing allocations
        torch.cuda.empty_cache()
        
        # Set memory fraction to prevent OOM
        torch.cuda.set_per_process_memory_fraction(0.95)  # Use 95% of GPU memory
        
        # Reset peak memory stats
        torch.cuda.reset_peak_memory_stats()
        
        # Set CUDA allocator settings
        os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'
    
    # Enable garbage collection
    import gc
    gc.collect()    
    # Enable CUDNN optimizations
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True
    
    # Enable TF32 on Ampere GPUs
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    
    # Parse arguments
    parser = argparse.ArgumentParser()
    
    # SageMaker specific arguments
    parser.add_argument('--model-dir', type=str, default=os.environ.get('SM_MODEL_DIR'))
    parser.add_argument('--output-data-dir', type=str, default=os.environ.get('SM_OUTPUT_DATA_DIR'))
    
    # Training arguments
    parser.add_argument('--bucket-name', type=str, required=True)
    parser.add_argument('--base-data-dir', type=str, required=True)
    parser.add_argument('--base-attr-dir', type=str, required=True)
    parser.add_argument('--experiment-name', type=str, default='experiment_1')
    parser.add_argument('--chunk-size', type=int, default=200)
    parser.add_argument('--total-epochs', type=int, default=200)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--hidden-size', type=int, default=256)
    parser.add_argument('--num-layers', type=int, default=2)
    parser.add_argument('--learning-rate', type=float, default=0.0001)
    parser.add_argument('--dropout', type=float, default=0.3)
    
    args = parser.parse_args()
    
    # Define checkpoint directory EARLY, before any usage
    checkpoint_dir = os.environ.get('SM_CHECKPOINT_DIR', '/opt/ml/checkpoints')
    
    # Log configuration
    print(f"Starting PyTorch Lightning training with experiment: {args.experiment_name}")
    print(f"Model directory: {args.model_dir}")
    print(f"Output directory: {args.output_data_dir}")
    print(f"Checkpoint directory: {checkpoint_dir}")
    
    # Detect device
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")
    if device == 'cuda':
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_data_dir, exist_ok=True)
    if args.model_dir:
        os.makedirs(args.model_dir, exist_ok=True)
    
    # Load watershed data
    print("Loading watershed data...")
    watershed_df = identify_all_available_watersheds(args.bucket_name, args.base_data_dir)
    watershed_df = watershed_df[watershed_df['subdirectory_name'] == 'camels']
    watershed_df = watershed_df#.iloc[0:100]  # Adjust as needed
    print(f"Found {len(watershed_df)} watersheds")
    
    # Create data module
    print("Creating Lightning DataModule...")
    data_module = CaravanDataModule(
        watersheds_df=watershed_df,
        bucket_name=args.bucket_name,
        base_data_dir=args.base_data_dir,
        base_attr_dir=args.base_attr_dir,
        sequence_length=365,
        batch_size=args.batch_size,
        num_workers=8,
        chunk_size=args.chunk_size,
        train_split=0.6,
        val_split=0.2,
        random_seed=42
    )
    
    # Setup data (computes normalization params)
    print("Setting up data splits and computing normalization parameters...")
    data_module.setup('fit')
    
    # Save normalization parameters
    norm_params_path = os.path.join(args.output_data_dir, 'normalization_params.json')
    with open(norm_params_path, 'w') as f:
        json_params = {
            'static_means': {k: float(v) for k, v in data_module.norm_params['static_means'].items()},
            'static_stds': {k: float(v) for k, v in data_module.norm_params['static_stds'].items()},
            'dynamic_means': {k: float(v) for k, v in data_module.norm_params['dynamic_means'].items()},
            'dynamic_stds': {k: float(v) for k, v in data_module.norm_params['dynamic_stds'].items()},
            'target_mean': float(data_module.norm_params['target_mean']),
            'target_std': float(data_module.norm_params['target_std'])
        }
        json.dump(json_params, f, indent=2)
    print(f"Saved normalization parameters to {norm_params_path}")
    
    # Create model
    print("Creating Lightning model...")
    model = StreamflowLightningModule(
        dynamic_input_size=len(data_module.dynamic_cols_no_target),
        static_input_size=len(data_module.static_cols),
        hidden_size=args.hidden_size,
        num_layers=args.num_layers,
        dropout=args.dropout,
        learning_rate=args.learning_rate,
        norm_params=data_module.norm_params,
        target_noise_std=0.1,
        target_noise_type='multiplicative',
        gradient_clip_val=5.0,
        lr_scheduler_params={
            'factor': 0.5,
            'patience': 2,
            'min_lr': 1e-7,
            'threshold': 0.001,
            'cooldown': 1,
            'verbose': True
        }
    )
    
    # Create trainer with SageMaker-compatible settings
    trainer = pl.Trainer(
        max_epochs=args.total_epochs,
        accelerator='gpu' if device == 'cuda' else 'cpu',
        devices=1,
        precision=16 if device == 'cuda' else 32,  # Mixed precision on GPU only
        gradient_clip_val=5.0,
        accumulate_grad_batches=4, # Accumulate 4 batches before update which increases the effective batch size
        num_sanity_val_steps=2,  # Run validation at start to initialize metrics
#        val_check_interval=1 if args.total_epochs <= 10 else None,  # Validate every epoch for short runs
        
        # Callbacks
        callbacks=[
            ValidationFrequencyCallback(initial_frequency=1, later_frequency=10, switch_epoch=10),
            # Checkpoint callback with SageMaker paths
            pl.callbacks.ModelCheckpoint(
                dirpath=checkpoint_dir,
                filename=f'{args.experiment_name}_epoch_{{epoch}}_kge_{{val_kge:.3f}}',
                monitor='val_kge',  # Still monitor KGE for saving best model
                mode='max',  # Maximize KGE
                save_top_k=3,
                save_last=True,
                save_on_train_epoch_end=True,
                every_n_epochs=1
            ),
            # Update Early stopping to use KGE loss
            pl.callbacks.EarlyStopping(
                monitor='val_kge',  # This is now negative KGE; wait now positive kge
                patience=10,
                mode='max',  # Minimize negative KGE; wait now maximize kge
                verbose=True
            ),
            # Learning rate monitor
            pl.callbacks.LearningRateMonitor(logging_interval='epoch'),
            # Progress bar refresh rate
            pl.callbacks.TQDMProgressBar(refresh_rate=50),
        ],
        
        # Logging - TensorBoard with SageMaker output path
        logger=pl.loggers.TensorBoardLogger(
            save_dir=args.output_data_dir,
            name=args.experiment_name,
            version='lightning'
        ),
        
        # Validation frequency
        check_val_every_n_epoch=1, #if args.total_epochs > 10 else None,  # Modified: validate every epoch initially
        limit_val_batches=1.0,  # Use all validation data
        
        # Enable better error messages
        enable_model_summary=True,
        
        # For distributed training readiness
        strategy='auto',
        
        # Other settings
        log_every_n_steps=50,
        enable_progress_bar=True,
        enable_checkpointing=True,
        deterministic=False,  # Faster training
        benchmark=True,  # CUDNN auto-tuner
    )
    
    # Check for existing checkpoints
    ckpt_path = None  # Initialize checkpoint path
    
    if os.path.exists(checkpoint_dir):
        checkpoint_files = sorted([f for f in os.listdir(checkpoint_dir) 
                                 if f.endswith('.ckpt') and 'last' in f])
        if checkpoint_files:
            ckpt_path = os.path.join(checkpoint_dir, checkpoint_files[-1])
            print(f"Found checkpoint: {ckpt_path}")
            print("Will resume training from this checkpoint")
        else:
            print("No checkpoint files found, starting fresh training")
    else:
        print(f"Checkpoint directory {checkpoint_dir} does not exist, starting fresh training")
    
    # Train the model
    print("Starting training...")
    trainer.fit(model, data_module, ckpt_path=ckpt_path)
    
    # Test the model
    print("Running test evaluation...")
    test_results = trainer.test(model, data_module)
    
    # Save the best model to SageMaker model directory
    if hasattr(trainer, 'checkpoint_callback') and trainer.checkpoint_callback:
        best_model_path = trainer.checkpoint_callback.best_model_path
        if best_model_path and os.path.exists(best_model_path):
            print(f"Copying best model from {best_model_path} to {args.model_dir}")
            import shutil
            shutil.copy(best_model_path, os.path.join(args.model_dir, 'model.ckpt'))
            
            # Also save as state dict for compatibility
            checkpoint = torch.load(best_model_path, map_location=device)
            torch.save(checkpoint['state_dict'], os.path.join(args.model_dir, 'model.pt'))
            
            # Save model info - APPLY tensor_to_python HERE
            best_score = trainer.checkpoint_callback.best_model_score
            if best_score is not None:
                best_score = tensor_to_python(best_score)

            # Convert best_k_models as well - it might contain tensors
            best_k_models = trainer.checkpoint_callback.best_k_models
            if best_k_models is not None:
                best_k_models = tensor_to_python(best_k_models)

            model_info = {
                'best_epoch': best_k_models,  # Now converted
                'best_score': best_score,
                'monitor_metric': trainer.checkpoint_callback.monitor,
                'feature_info': {
                    'static_features': data_module.static_cols,
                    'dynamic_features': data_module.dynamic_cols_no_target,
                    'target': 'streamflow'
                }
            }
            # Extra safety: convert the entire dict
            model_info = tensor_to_python(model_info)
            
            with open(os.path.join(args.model_dir, 'model_info.json'), 'w') as f:
                json.dump(model_info, f, indent=2)
        else:
            print("Warning: No best model checkpoint found")
    
    # Save training results
    results = {
        'test_results': tensor_to_python(test_results),  # Apply here
        'total_epochs_trained': trainer.current_epoch,
        'experiment_name': args.experiment_name,
        'hyperparameters': {
            'batch_size': args.batch_size,
            'hidden_size': args.hidden_size,
            'num_layers': args.num_layers,
            'learning_rate': args.learning_rate,
            'dropout': args.dropout,
            'chunk_size': args.chunk_size
        },
        'feature_info': {
            'static_features': data_module.static_cols,
            'dynamic_features': data_module.dynamic_cols_no_target,
            'target': 'streamflow',
            'num_static_features': len(data_module.static_cols),
            'num_dynamic_features': len(data_module.dynamic_cols_no_target)
        },
        'data_info': {
            'total_watersheds': len(watershed_df),
            'train_watersheds': len(data_module.train_watersheds),
            'val_watersheds': len(data_module.val_watersheds),
            'test_watersheds': len(data_module.test_watersheds)
        }
    }

    # Add timing information if available
    if hasattr(trainer, 'logged_metrics'):
        results['final_metrics'] = tensor_to_python(trainer.logged_metrics)  # Use tensor_to_python here too

    results_path = os.path.join(args.output_data_dir, 'results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\n" + "="*50)
    print("Training completed successfully!")
    print(f"Best model saved to: {args.model_dir}")
    print(f"Results saved to: {results_path}")
    print(f"TensorBoard logs available at: {args.output_data_dir}/{args.experiment_name}")
    print(f"Normalization parameters saved to: {norm_params_path}")
    print("="*50)


if __name__ == '__main__':
    main()
