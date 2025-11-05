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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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
import shutil
import subprocess

from sagemaker import Session
from sagemaker.pytorch import PyTorch



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

def load_and_process_watershed_data(
    bucket_name,
    base_directory,
    subdirectory_name,
    watershedID,
    variables_to_extract=None,
    max_missing_pct: float = 1.0
):
    """
    Loads and processes timeseries data for a specific watershed from an S3 NetCDF file.
    Handles case variations in filenames (e.g., GRDC vs grdc).
    Handles both NetCDF3 (scipy) and NetCDF4 (netcdf4/h5netcdf) formats.
    Excludes watersheds with > max_missing_pct percent missing streamflow BEFORE any filling.
    """
    # Build path variants to try
    default_path = f"s3://{bucket_name}/{base_directory.strip('/')}/{subdirectory_name}/{subdirectory_name}_{watershedID}.nc"
    
    # For GRDC, also try uppercase filename variant
    path_variants = [default_path]
    if subdirectory_name.lower() == 'grdc':
        uppercase_path = f"s3://{bucket_name}/{base_directory.strip('/')}/{subdirectory_name}/GRDC_{watershedID}.nc"
        path_variants.append(uppercase_path)
    
    # Define default variables if not specified
    default_variables = [
        'streamflow',
        'surface_net_solar_radiation_mean',
        'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean'
    ]
    
    if variables_to_extract is None:
        variables_to_extract = default_variables
        print(f"Using default variables: {variables_to_extract}")
    else:
        print(f"Using user-specified variables: {variables_to_extract}")
    
    # Initialize S3 filesystem
    s3 = s3fs.S3FileSystem(anon=False)
    
    # Try each path variant
    processed_df = pd.DataFrame()
    successful_path = None
    
    for path_to_try in path_variants:
        print(f"Attempting to open: {path_to_try}")
        
        try:
            # Check if file exists
            if not s3.exists(path_to_try):
                if path_to_try == path_variants[-1]:  # Last variant
                    print(f"Error: None of the path variants exist.")
                    print(f"Tried: {path_variants}")
                continue
            
            # Try different engines to handle both NetCDF3 and NetCDF4
            engines_to_try = ['h5netcdf', 'netcdf4', 'scipy']
            dataset_opened = False
            
            for engine in engines_to_try:
                try:
                    with s3.open(path_to_try, mode='rb') as s3file:
                        # Try opening with current engine
                        with xr.open_dataset(s3file, engine=engine, cache=False) as ds:
                            print(f"Successfully opened via s3fs and {engine}: {path_to_try}")
                            successful_path = path_to_try
                            dataset_opened = True
                            
                            # Check available variables
                            available_vars = list(ds.data_vars)
                            missing_vars = [v for v in variables_to_extract if v not in available_vars]
                            
                            if missing_vars:
                                print(f"Warning: The following requested variables are missing from the file: {missing_vars}")
                                current_variables_to_extract = [v for v in variables_to_extract if v in available_vars]
                                if not current_variables_to_extract:
                                    print("Error: None of the requested variables are available. Cannot proceed.")
                                    return pd.DataFrame()
                            else:
                                current_variables_to_extract = variables_to_extract
                            
                            # Extract subset of variables
                            ds_subset = ds[current_variables_to_extract]
                            df = ds_subset.to_dataframe().reset_index()
                            
                            # Handle date coordinate
                            if 'date' not in df.columns:
                                time_coord_name = None
                                if 'time' in ds.coords:
                                    time_coord_name = 'time'
                                elif ds.coords:
                                    potential_time_coords = [c for c in ds.coords if ds[c].ndim == 1]
                                    if len(potential_time_coords) == 1:
                                        time_coord_name = potential_time_coords[0]
                                        print(f"Using coordinate '{time_coord_name}' as date.")
                                        df = df.rename(columns={time_coord_name: 'date'})
                                
                                if 'date' not in df.columns:
                                    raise ValueError("Date coordinate not found or could not be identified.")
                            
                            # Decode date if necessary
                            if not pd.api.types.is_datetime64_any_dtype(df['date']):
                                print("Attempting to decode date coordinate...")
                                try:
                                    time_var_name = time_coord_name if 'time_coord_name' in locals() and time_coord_name else 'date'
                                    time_var = ds[time_var_name]
                                    units = time_var.attrs.get('units', None)
                                    calendar = time_var.attrs.get('calendar', 'standard')
                                    
                                    if units:
                                        print(f"Using units: '{units}', calendar: '{calendar}'")
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
                            
                            # Process streamflow if present
                            if 'streamflow' in current_variables_to_extract and 'streamflow' in df.columns:
                                # Trim to first/last valid streamflow
                                first_valid_idx = df['streamflow'].first_valid_index()
                                last_valid_idx = df['streamflow'].last_valid_index()
                                
                                if first_valid_idx is None or last_valid_idx is None:
                                    print("Warning: No valid streamflow data found. Returning empty DataFrame.")
                                    return pd.DataFrame()
                                
                                print(f"Trimming data based on streamflow: First valid at index {first_valid_idx}, Last valid at index {last_valid_idx}")
                                df_trimmed = df.loc[first_valid_idx:last_valid_idx].reset_index(drop=True)
                                
                                # Check missing percentage after trimming
                                missing_after_trim = df_trimmed['streamflow'].isna().sum()
                                total_after_trim = len(df_trimmed)
                                missing_pct_after = (missing_after_trim / max(total_after_trim, 1)) * 100.0
                                
                                if missing_pct_after > max_missing_pct:
                                    print(f"WARNING: {missing_pct_after:.2f}% of streamflow still missing after trimming (>{max_missing_pct}% threshold)")
                                    print(f"  Watershed {subdirectory_name}_{watershedID} will be excluded")
                                    return pd.DataFrame()
                                
                                print(f"Data shape after trimming: {df_trimmed.shape}")
                                print(f"Missing streamflow after trimming: {missing_after_trim} days ({missing_pct_after:.3f}%)")
                            else:
                                df_trimmed = df.copy()
                                if 'streamflow' not in current_variables_to_extract:
                                    print("Skipping streamflow NaN trimming as 'streamflow' was not requested.")
                                else:
                                    print("Warning: 'streamflow' requested but not found. Skipping trimming.")
                            
                            # Report remaining NaNs (no filling here)
                            rem = df_trimmed.isnull().sum()
                            rem = rem[rem > 0]
                            if not rem.empty:
                                print("\nWarning: Remaining NaNs found after trimming:")
                                for k, v in rem.items():
                                    print(f"  {k}: {int(v)}")
                            
                            processed_df = df_trimmed
                            break  # Successfully processed with this engine
                            
                except ImportError as ie:
                    if engine == 'h5netcdf':
                        print(f"h5netcdf not available, trying next engine...")
                        continue
                    elif engine == 'netcdf4':
                        print(f"netcdf4 not available, trying next engine...")
                        continue
                    else:
                        print(f"Engine {engine} import error: {ie}")
                        continue
                except Exception as engine_error:
                    if engine == engines_to_try[-1]:  # Last engine
                        print(f"Failed with all engines. Last error with {engine}: {engine_error}")
                    else:
                        print(f"Failed with {engine}, trying next engine...")
                    continue
            
            if dataset_opened:
                break  # Successfully opened with one of the engines
                
        except FileNotFoundError:
            if path_to_try == path_variants[-1]:  # Last variant
                print(f"Error: File not found at any of the attempted paths")
                print(f"Tried: {path_variants}")
            continue
        except Exception as e:
            if path_to_try == path_variants[-1]:  # Last variant
                print(f"An unexpected error occurred while processing {path_to_try}: {e}")
                import traceback
                traceback.print_exc()
            continue
    
    # If no successful path was found and we have an empty dataframe
    if successful_path is None and processed_df.empty:
        print(f"Error: Could not successfully open any variant for {subdirectory_name}_{watershedID}")
        print(f"Attempted paths: {path_variants}")
    
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
        'pet_mm_syr',
        'aet_mm_syr',  # optional aet_mm_s01,...,12
        'pre_mm_syr',
        'tmp_dc_syr', 'tmp_dc_smn', 'tmp_dc_smx',
        'snw_pc_syr', 'snw_pc_smx',
        'swc_pc_syr',
#        'gdp_ud_sav',
        'inu_pc_slt', 'inu_pc_smn', 'inu_pc_smx',  # inundation percent
        'run_mm_syr',
        'dis_m3_pmn', 'dis_m3_pmx', 'dis_m3_pyr',
        'lka_pc_sse',  # limnicity pct area
        'lkv_mc_usu',  # lake vol
        'rev_mc_usu',  # reservoir volume
#        'ria_ha_usu',  # river area
        'riv_tc_usu',  # river volume
        'dor_pc_pva',  # degree of regulation
        'gwt_cm_sav',  # groundwater table depth
        'ele_mt_sav', 'ele_mt_smn', 'ele_mt_mx',
#        'slp_dg_sav',  # terrain slope
        'sgr_dk_sav',  # stream gradient
        'ari_ix_sav',
        'glc_cl_smj',  # land cover classes
        #'glc_cl_s01','glc_cl_s02','glc_cl_s03','glc_cl_s04','glc_cl_s05','glc_cl_s06','glc_cl_s07','glc_cl_s08','glc_cl_s09','glc_cl_s10','glc_cl_s11',
        #'glc_cl_s12','glc_cl_s13','glc_cl_s14','glc_cl_s15','glc_cl_s16','glc_cl_s17','glc_cl_s18','glc_cl_s19','glc_cl_s20','glc_cl_s21','glc_cl_s22',
        'for_pc_sse','pst_pc_sse',
        'crp_pc_sse', 'ire_pc_sse', 'gla_pc_sse', 'prm_pc_sse',
        'cly_pc_sav', 'slt_pc_sav', 'snd_pc_sav',
        'soc_th_sav', 'lit_cl_smj', 'kar_pc_sse',
#        'rdd_mk_sav',
        'urb_pc_sse',
        'hdi_ix_sav', 'nli_ix_sav', 'ppd_pk_sav',  # human development index, nighttime lights, pop density
        'cls_cl_smj', 'clz_cl_smj',  # climate strata and zones
    ]
    
    if hydroatlas_cols_to_extract is None:
        hydroatlas_cols_to_use = default_hydroatlas_cols
        print(f"Using default HydroATLAS columns.")
    else:
        hydroatlas_cols_to_use = hydroatlas_cols_to_extract
        print(f"Using user-specified HydroATLAS columns.")
    # --- End HydroATLAS column handling ---

    # --- Define Other Columns (remains fixed) ---
    other_cols = ['area']#,'gauge_lat', 'gauge_lon']

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

def calculate_nse(observed: np.ndarray, simulated: np.ndarray) -> float:
    """
    Calculate Nash–Sutcliffe Efficiency (NSE).
    Returns NaN if insufficient valid data or near-zero variance in observations.
    """
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()

    valid = ~np.isnan(observed) & ~np.isnan(simulated) & np.isfinite(observed) & np.isfinite(simulated)
    if np.sum(valid) < 2:
        return np.nan

    obs = observed[valid].astype(np.float64)
    sim = simulated[valid].astype(np.float64)

    denom = np.sum((obs - np.mean(obs))**2)
    if denom <= 1e-12:
        return np.nan

    num = np.sum((sim - obs)**2)
    return 1.0 - (num / denom)




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


class ProgressiveExtremeLoss(nn.Module):
    """
    Multi-component loss that progressively shifts focus from general to extreme errors.
    Combines MAE (L1), RMSE (L2), and CMCE (L3^(1/3)) with time-varying weights.
    
    All three components are on the same scale (normalized units), making them directly comparable.
    - MAE: treats all errors equally
    - RMSE: moderately emphasizes larger errors  
    - CMCE: strongly emphasizes extreme errors
    """
    def __init__(self, warmup_epochs=20, transition_epochs=30, 
                 final_weights=(0.2, 0.3, 0.5), outlier_percentile=99):
        super().__init__()
        self.warmup_epochs = warmup_epochs
        self.transition_epochs = transition_epochs
        self.final_weights = final_weights
        self.outlier_percentile = outlier_percentile
        
    def forward(self, pred, target, epoch=0):
        pred = pred.view(-1)
        target = target.view(-1)
        
        residuals = pred - target
        abs_residuals = torch.abs(residuals)
        
        # Optional: Clip extreme outliers for stability (still included but dampened)
        if self.outlier_percentile < 100:
            threshold = torch.quantile(abs_residuals, self.outlier_percentile / 100)
            clipped_residuals = torch.sign(residuals) * torch.minimum(abs_residuals, threshold * 2)
        else:
            clipped_residuals = residuals
        
        # MAE
        mae = torch.mean(abs_residuals)
        
        # RMSE (use clipped for stability in high-power terms)
        mse = torch.mean(clipped_residuals ** 2)
        rmse = torch.sqrt(mse + 1e-8)
        
        # CMCE with numerical stability
        # Use clipped residuals and add small epsilon to prevent overflow
        safe_residuals = torch.clamp(clipped_residuals, -10, 10)  # Prevent extreme overflow
        cubed = safe_residuals ** 3
        mean_cubed = torch.mean(cubed)
        
        # Stable cube root that preserves sign
        if mean_cubed.abs() < 1e-8:
            cmce_loss = torch.zeros_like(mean_cubed)
        else:
            sign = torch.sign(mean_cubed)
            magnitude = torch.abs(mean_cubed) ** (1/3)
            cmce_loss = torch.abs(sign * magnitude)
        
        # Rest of the weighting logic remains the same...
        if epoch < self.warmup_epochs:
            return mse
        else:
            progress = min(1.0, (epoch - self.warmup_epochs) / self.transition_epochs)
            
            w_mae = (1 - progress) * 0.6 + progress * self.final_weights[0]
            w_rmse = (1 - progress) * 0.3 + progress * self.final_weights[1]  
            w_cmce = (1 - progress) * 0.1 + progress * self.final_weights[2]
            
            total = w_mae + w_rmse + w_cmce
            w_mae, w_rmse, w_cmce = w_mae/total, w_rmse/total, w_cmce/total
            
            loss = w_mae * mae + w_rmse * rmse + w_cmce * cmce_loss
            
            return loss, {'mae': mae, 'rmse': rmse, 'cmce': cmce_loss,
                         'w_mae': w_mae, 'w_rmse': w_rmse, 'w_cmce': w_cmce}

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
    Reverse normalization for streamflow predictions when targets are log-transformed then z-scored.
    
    Steps:
    1. Reverse z-score: multiply by std and add mean
    2. Reverse log1p: apply expm1 (exp(x) - 1)
    """
    preds = np.asarray(predictions, dtype=float).flatten()
    
    # Get normalization parameters
    target_mean = norm_params.get('target_mean', 0.0)
    target_std = norm_params.get('target_std', 1.0)
    
    # Step 1: Reverse z-score normalization
    preds_log_space = preds * target_std + target_mean
    
    # Step 2: Reverse log1p transform (inverse is expm1)
    original = np.expm1(preds_log_space)  # exp(x) - 1, inverse of log(1 + x)
    
    # Safety: clip to reasonable streamflow range and handle any numerical issues
    original = np.clip(original, 0.0, 100000.0)
    
    # Handle any potential NaN/Inf from exp
    original = np.nan_to_num(original, nan=0.0, posinf=100000.0, neginf=0.0)
    
    return original




def unnormalize_features(features_df, norm_params, feature_type='dynamic'):
    """
    Reverse normalization for features (useful for interpretation).
    Handles:
      - Dynamic z-score features (mean/std)
      - Dynamic min–max features (to [0,1])
      - Dynamic std-only features (multiply by std)
      - Static z-score features
    """
    unnormalized_df = features_df.copy()

    if feature_type == 'dynamic':
        means = norm_params.get('dynamic_means', {})
        stds = norm_params.get('dynamic_stds', {})
        mins = norm_params.get('dynamic_mins', {})
        maxs = norm_params.get('dynamic_maxs', {})
        std_only_cols = set(norm_params.get('std_only_cols', []))
        minmax_cols = set(norm_params.get('minmax_cols', []))

        for col in features_df.columns:
            if col in minmax_cols and col in mins and col in maxs:
                # Invert min–max: x = x_norm * (max - min) + min
                rng = maxs[col] - mins[col]
                if rng > 0:
                    unnormalized_df[col] = features_df[col] * rng + mins[col]
                else:
                    unnormalized_df[col] = mins[col]
            elif col in std_only_cols and col in stds:
                # Invert std-only: x = x_norm * std
                std_val = stds[col]
                if std_val and std_val > 0:
                    unnormalized_df[col] = features_df[col] * std_val
                else:
                    unnormalized_df[col] = 0.0
            else:
                # Default z-score inversion: x = x_norm * std + mean
                if col in means and col in stds and stds[col] > 0:
                    unnormalized_df[col] = features_df[col] * stds[col] + means[col]
    else:
        # Static: z-score inversion
        means = norm_params.get('static_means', {})
        stds = norm_params.get('static_stds', {})
        for col in features_df.columns:
            if col in means and col in stds and stds[col] > 0:
                unnormalized_df[col] = features_df[col] * stds[col] + means[col]

    return unnormalized_df




#################################################################################################



import importlib.util

def _import_feature_analysis_safe(preferred_paths=None):
    """
    Try `import featureAnalysis` normally; if that fails, search a few
    likely locations and import by file path. Returns the module or raises.
    """
    try:
        import featureAnalysis as fa
        return fa
    except ModuleNotFoundError:
        pass

    # Build candidate paths to featureAnalysis.py
    candidates = []
    # 1) Same directory as this script
    try:
        this_dir = os.path.dirname(os.path.abspath(__file__))
        candidates.append(os.path.join(this_dir, 'featureAnalysis.py'))
    except Exception:
        pass
    # 2) Current working dir
    candidates.append(os.path.join(os.getcwd(), 'featureAnalysis.py'))
    # 3) Any preferred paths provided by caller
    if preferred_paths:
        for p in preferred_paths:
            candidates.append(os.path.join(p, 'featureAnalysis.py'))

    for path in candidates:
        if path and os.path.exists(path):
            spec = importlib.util.spec_from_file_location('featureAnalysis', path)
            if spec and spec.loader:
                fa = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(fa)
                print(f"Imported featureAnalysis from file: {path}")
                return fa

    raise ModuleNotFoundError(
        "featureAnalysis module not found. Ensure featureAnalysis.py is included "
        "in your training source_dir, or place it alongside train.py."
    )



def prepare_ptf_dataframe(
    watershed_subset_df,
    bucket_name,
    base_data_dir,
    base_attr_dir,
    include_static_attributes=True,
    norm_params=None,
    compute_norm_params=False,
    max_missing_pct: float = 1.0  # NEW
):
    """
    Loads data for multiple watersheds and prepares a single DataFrame suitable for modeling.
    Excludes watersheds with > max_missing_pct missing streamflow (handled upstream).
    """
    weather_cols = [
        'surface_net_solar_radiation_mean',
        'temperature_2m_max',
        'temperature_2m_mean',
        'temperature_2m_min',
        'total_precipitation_sum',
        'dewpoint_temperature_2m_mean'
    ]
    precip_cols = [#'Precip_smoothed_5day',
                   #'Precip_lagged_90day',
                   'total_likely_snow_sum',
                   'total_likely_rain_sum']
    time_varying_cols = weather_cols + precip_cols

    all_data_list = []
    static_attribute_names = []
    excluded_watersheds = []

    for _, row in watershed_subset_df.iterrows():
        subdir = row['subdirectory_name']
        ws_id = row['watershedID']
        group_id = f"{subdir}_{ws_id}"
        print(f"  Loading data for {group_id}...")

        try:
            # Load timeseries with missing % restriction
            timeseries_df = load_and_process_watershed_data(
                bucket_name, base_data_dir, subdir, ws_id, max_missing_pct=max_missing_pct
            )
            if timeseries_df.empty:
                excluded_watersheds.append({
                    'watershed_id': group_id,
                    'reason': 'Missing data exceeds threshold or no valid data'
                })
                print(f"  Warning: No valid timeseries data for {group_id}. Skipping.")
                continue

            date_col = 'date'
            target_col = 'streamflow'
            # Fix non-positive target (for log1p later)
            invalid_mask = timeseries_df[target_col] <= 0
            if invalid_mask.any():
                n_invalid = int(invalid_mask.sum())
                print(f"WARNING: Found {n_invalid} non-positive streamflow values in {group_id}")
                timeseries_df.loc[invalid_mask, target_col] = 1e-6

            # Add precip features if available; else create defaults
            if 'total_precipitation_sum' in timeseries_df.columns and 'temperature_2m_mean' in timeseries_df.columns:
                timeseries_df['Precip_smoothed_5day'] = (
                    timeseries_df['total_precipitation_sum']
                    .rolling(window=5, min_periods=1, center=False)
                    .mean()
                )
                timeseries_df['Precip_lagged_90day'] = (
                    timeseries_df['total_precipitation_sum']
                    .rolling(window=90, min_periods=1, center=False)
                    .mean()
                )
                timeseries_df['Precip_smoothed_5day'] = timeseries_df['Precip_smoothed_5day'].ffill().bfill()
                timeseries_df['Precip_lagged_90day'] = timeseries_df['Precip_lagged_90day'].ffill().bfill()

                # NEW: Snow/rain partitioning based on temperature (vectorized)
                temp = timeseries_df['temperature_2m_mean'].values
                precip = timeseries_df['total_precipitation_sum'].values
                
                # Snow: precip when temp < 0°C
                timeseries_df['total_likely_snow_sum'] = np.where(temp < 0, precip, 0.0)
                
                # Rain: precip when temp > 0°C
                timeseries_df['total_likely_rain_sum'] = np.where(temp > 0, precip, 0.0)

                
            else:
                print(f"  Warning: Missing required columns for {group_id} (precip/temp features set to 0)")
                timeseries_df['Precip_smoothed_5day'] = 0.0
                timeseries_df['Precip_lagged_90day'] = 0.0
                timeseries_df['total_likely_snow_sum'] = 0.0
                timeseries_df['total_likely_rain_sum'] = 0.0

            # Subset to required columns
            keep_cols = [date_col, target_col] + weather_cols + precip_cols
            timeseries_df = timeseries_df[keep_cols]

            # Sort and reindex to full daily range
            timeseries_df[date_col] = pd.to_datetime(timeseries_df[date_col])
            timeseries_df = timeseries_df.sort_values(by=date_col)

            if timeseries_df.empty:
                print(f"  Warning: Empty timeseries after preprocessing for {group_id}. Skipping.")
                continue

            date_range = pd.date_range(start=timeseries_df[date_col].min(),
                                       end=timeseries_df[date_col].max(), freq='D')
            timeseries_df = timeseries_df.set_index(date_col).reindex(date_range)

            # Fill NaNs now that watershed passed missing % checks
            # 1) target
            timeseries_df[target_col] = timeseries_df[target_col].ffill()
            # 2) covariates
            for col in weather_cols + precip_cols:
                if col in timeseries_df.columns:
                    timeseries_df[col] = timeseries_df[col].ffill()
            # Back-fill any leading NaNs
            timeseries_df = timeseries_df.bfill()

            # Final NaN check on required columns
            if timeseries_df[[target_col] + weather_cols].isnull().any().any():
                print(f"  Warning: NaNs still present after filling for {group_id}. Skipping.")
                excluded_watersheds.append({
                    'watershed_id': group_id,
                    'reason': 'NaNs remained after filling'
                })
                continue

            # Restore date column
            timeseries_df = timeseries_df.reset_index().rename(columns={'index': date_col})
            timeseries_df['time_idx'] = (timeseries_df[date_col] - timeseries_df[date_col].min()).dt.days
            timeseries_df['group_id'] = group_id

            # Static attributes
            if include_static_attributes:
                print(f"  Searching for attributes for gage_id: {group_id}")
                static_attrs = get_watershed_attributes(bucket_name, base_attr_dir, subdir, ws_id)
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
                        # Drop extras not in the initial set
                        drop_cols = [c for c in current_static_names if c not in static_attribute_names]
                        if drop_cols:
                            timeseries_df = timeseries_df.drop(columns=drop_cols, errors='ignore')
                else:
                    print(f"  Warning: Could not load static attributes for {group_id}. Filling with 0.")
                    for col_name_str in static_attribute_names:
                        timeseries_df[col_name_str] = 0.0

            all_data_list.append(timeseries_df)

        except Exception as e:
            print(f"  Error processing watershed {group_id}: {e}")
            traceback.print_exc()
            excluded_watersheds.append({'watershed_id': group_id, 'reason': str(e)})
            continue

    # Final checks after loop
    if not all_data_list:
        print("Error: No data loaded successfully.")
        if excluded_watersheds:
            print(f"\nExcluded {len(excluded_watersheds)} watersheds due to missing data or errors:")
            for ws in excluded_watersheds[:10]:
                print(f"  - {ws['watershed_id']}: {ws['reason']}")
            if len(excluded_watersheds) > 10:
                print(f"  ... and {len(excluded_watersheds) - 10} more")
        return pd.DataFrame(), [], [], norm_params

    combined_df = pd.concat(all_data_list, ignore_index=True)

    # Summary
    print(f"\nData Loading Summary:")
    print(f"  Successfully loaded: {len(all_data_list)} watersheds")
    print(f"  Excluded: {len(excluded_watersheds)} watersheds")
    print(f"  Total attempted: {len(watershed_subset_df)} watersheds")

    # Optional diagnostic: extreme values on combined_df
    target_col = 'streamflow'
    for col in [c for c in weather_cols + [target_col] if c in combined_df.columns]:
        q99 = combined_df[col].quantile(0.99)
        q01 = combined_df[col].quantile(0.01)
        denom = max(q01 + 1e-6, 1e-6)
        if (q99 / denom) > 1000:
            print(f"WARNING: Extreme values in {col}: q01={q01:.4f}, q99={q99:.4f}")

    # NORMALIZE FEATURES
    if compute_norm_params:
        all_time_varying_cols = weather_cols + precip_cols
        dynamic_means = {col: combined_df[col].mean() for col in all_time_varying_cols}
        dynamic_stds = {col: combined_df[col].std() for col in all_time_varying_cols}
        std_only_cols = ['temperature_2m_min']
        minmax_cols = ['total_precipitation_sum', 'total_likely_snow_sum', 'total_likely_rain_sum']
        dynamic_mins = {col: float(combined_df[col].min()) for col in minmax_cols if col in combined_df.columns}
        dynamic_maxs = {col: float(combined_df[col].max()) for col in minmax_cols if col in combined_df.columns}

        norm_params = {
            'static_means': {col: combined_df[col].mean() for col in static_attribute_names},
            'static_stds': {col: combined_df[col].std() for col in static_attribute_names},
            'dynamic_means': dynamic_means,
            'dynamic_stds': dynamic_stds,
            'target_mean': float(np.mean(np.log1p(combined_df[target_col]))),
            'target_std': float(np.std(np.log1p(combined_df[target_col]))),
            'dynamic_mins': dynamic_mins,
            'dynamic_maxs': dynamic_maxs,
            'std_only_cols': std_only_cols,
            'minmax_cols': minmax_cols
        }

    if norm_params is not None:
        # Static z-score
        for col in static_attribute_names:
            mean_val = norm_params['static_means'].get(col, 0.0)
            std_val = norm_params['static_stds'].get(col, 1.0)
            combined_df[col] = (combined_df[col] - mean_val) / std_val if std_val > 0 else 0.0

        # Dynamic per-column strategy
        std_only_cols = set(norm_params.get('std_only_cols', []))
        minmax_cols = set(norm_params.get('minmax_cols', []))
        dyn_means = norm_params.get('dynamic_means', {})
        dyn_stds = norm_params.get('dynamic_stds', {})
        dyn_mins = norm_params.get('dynamic_mins', {})
        dyn_maxs = norm_params.get('dynamic_maxs', {})

        for col in weather_cols + precip_cols:
            if col not in combined_df.columns:
                continue
            if col in minmax_cols:
                min_v = dyn_mins.get(col, None)
                max_v = dyn_maxs.get(col, None)
                if min_v is not None and max_v is not None and max_v > min_v:
                    combined_df[col] = (combined_df[col] - min_v) / (max_v - min_v)
                    combined_df[col] = combined_df[col].clip(0.0, 1.0)
                else:
                    combined_df[col] = 0.0
            elif col in std_only_cols:
                std_val = dyn_stds.get(col, None)
                combined_df[col] = (combined_df[col] / std_val) if (std_val and std_val > 0) else 0.0
            else:
                mean_val = dyn_means.get(col, None)
                std_val = dyn_stds.get(col, None)
                if std_val and std_val > 0 and mean_val is not None:
                    combined_df[col] = (combined_df[col] - mean_val) / std_val
                else:
                    combined_df[col] = 0.0

        # Target: log1p then z-score
        combined_df[target_col] = np.log1p(combined_df[target_col])
        mean_target = norm_params['target_mean']
        std_target = norm_params['target_std']
        combined_df[target_col] = (combined_df[target_col] - mean_target) / std_target

        # Validation checks
        print("\nChecking for extreme normalized values...")
        extreme_found = False

        tgt_vals = combined_df[target_col]
        if tgt_vals.max() > 20 or tgt_vals.min() < -20:
            extreme_found = True
            print(f"WARNING: Extreme normalized target values detected!")
            print(f"  Min: {tgt_vals.min():.2f}, Max: {tgt_vals.max():.2f}")
            print(f"  Mean: {tgt_vals.mean():.2f}, Std: {tgt_vals.std():.2f}")
            combined_df[target_col] = combined_df[target_col].clip(-20, 20)

        for col in weather_cols + precip_cols:
            if col not in combined_df.columns:
                continue
            col_values = combined_df[col]
            if col in minmax_cols:
                if col_values.min() < -1e-6 or col_values.max() > 1 + 1e-6:
                    extreme_found = True
                    print(f"WARNING: Min-max normalized {col} out of [0,1] bounds! Clipping.")
                    combined_df[col] = col_values.clip(0.0, 1.0)
            else:
                if col_values.max() > 20 or col_values.min() < -20:
                    extreme_found = True
                    print(f"WARNING: Extreme normalized values in {col}!")
                    print(f"  Min: {col_values.min():.2f}, Max: {col_values.max():.2f}")
                    combined_df[col] = col_values.clip(-20, 20)

        for col in static_attribute_names:
            if col in combined_df.columns:
                col_values = combined_df[col]
                if col_values.max() > 20 or col_values.min() < -20:
                    extreme_found = True
                    print(f"WARNING: Extreme normalized values in static feature {col}!")
                    print(f"  Min: {col_values.min():.2f}, Max: {col_values.max():.2f}")
                    combined_df[col] = col_values.clip(-20, 20)

        if not extreme_found:
            print("  No extreme values found - normalization looks good!")

        numeric_cols = combined_df.select_dtypes(include=[np.number]).columns
        if combined_df[numeric_cols].isna().any().any():
            print(f"WARNING: NaN values detected after normalization! Filling with 0.")
            for col in numeric_cols:
                if combined_df[col].isna().any():
                    combined_df[col] = combined_df[col].fillna(0)

        if np.isinf(combined_df[numeric_cols].values).any():
            print(f"WARNING: Inf values detected after normalization! Replacing with ±10.")
            for col in numeric_cols:
                combined_df[col] = combined_df[col].replace([np.inf, -np.inf], [10, -10])

    # Ensure dtypes
    combined_df[target_col] = combined_df[target_col].astype(float)
    for col in time_varying_cols:
        if col in combined_df.columns:
            combined_df[col] = combined_df[col].astype(float)
    for col in static_attribute_names:
        if col in combined_df.columns:
            combined_df[col] = combined_df[col].astype(float)

    time_varying_cols = weather_cols + precip_cols

    print(f"Successfully prepared combined DataFrame with shape: {combined_df.shape}")
    print(f"Identified static reals: {static_attribute_names}")
    print(f"Identified time-varying reals: {time_varying_cols}")

    return combined_df, static_attribute_names, time_varying_cols, norm_params




class WatershedDataset(Dataset):
    """PyTorch Dataset for watershed timeseries data"""
    def __init__(
        self,
        data_df: pd.DataFrame,
        static_cols: List[str],
        dynamic_cols: List[str],
        target_col: str,
        sequence_length: int,
        stride: int = 1,  # NEW
    ):
        self.data_df = data_df.copy()
        if 'group_id' not in self.data_df.columns:
            raise KeyError(
                "Expected 'group_id' column in data_df but it was missing. "
                "Upstream prepare_ptf_dataframe likely returned an empty/malformed DataFrame."
            )
        
        self.static_cols = static_cols
        self.dynamic_cols = dynamic_cols
        self.target_col = target_col
        self.sequence_length = sequence_length
        self.stride = max(int(stride), 1)  # NEW

        # Group data by watershed
        self.groups = self.data_df.groupby('group_id')
        self.group_keys = list(self.groups.groups.keys())

        # Valid start indices for all groups, using configured stride
        self.valid_indices: List[Tuple[str, int]] = []
        for group_id in self.group_keys:
            group_data = self.groups.get_group(group_id)
            n_samples = len(group_data) - self.sequence_length
            if n_samples > 0:
                for i in range(0, n_samples, self.stride):
                    self.valid_indices.append((group_id, i))

    def __len__(self) -> int:
        return len(self.valid_indices)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        group_id, start_idx = self.valid_indices[idx]
        group_data = self.groups.get_group(group_id)

        end_idx = start_idx + self.sequence_length

        dynamic_seq = group_data[self.dynamic_cols].iloc[start_idx:end_idx].values
        static_features = group_data[self.static_cols].iloc[0].values
        target = group_data[self.target_col].iloc[end_idx:end_idx + 1].values[0]

        return (
            torch.as_tensor(dynamic_seq, dtype=torch.float32),
            torch.as_tensor(static_features, dtype=torch.float32),
            torch.as_tensor([target], dtype=torch.float32),
        )


class CausalConv1d(nn.Module):
    """Causal convolution with dilation for capturing long-range dependencies"""
    def __init__(self, in_channels, out_channels, kernel_size, dilation=1):
        super().__init__()
        self.padding = (kernel_size - 1) * dilation
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, 
                             padding=self.padding, dilation=dilation)
    
    def forward(self, x):
        # Fix: Handle edge case when padding = 0
        if self.padding == 0:
            return self.conv(x)
        # Remove future timesteps that were added by padding
        return self.conv(x)[:, :, :-self.padding]
    

class LearnableFeatureSmoothing(nn.Module):
    """
    Learnable smoothing filters for each input feature.
    Applies different 1D convolutions to each feature channel independently.
    """
    def __init__(self, n_features, kernel_sizes=None, initialize_uniform=True):
        super().__init__()
        
        # Default kernel sizes to test
        if kernel_sizes is None:
            kernel_sizes = [3, 5, 7, 15, 31]  # Various time scales
        
        self.n_features = n_features
        self.kernel_sizes = kernel_sizes
        
        # Create a Conv1d for each feature with learnable weights
        self.smoothing_convs = nn.ModuleDict()
        
        for feat_idx in range(n_features):
            # Each feature gets its own set of multi-scale smoothing filters
            feat_convs = nn.ModuleList()
            for kernel_size in kernel_sizes:
                conv = nn.Conv1d(
                    in_channels=1,
                    out_channels=1, 
                    kernel_size=kernel_size,
                    padding=kernel_size-1,  # Padding for same-size output
                    bias=False
                )
                
                if initialize_uniform:
                    # Initialize to uniform (standard moving average)
                    with torch.no_grad():
                        conv.weight.fill_(1.0 / kernel_size)
                else:
                    # Xavier/Kaiming initialization for learning from scratch
                    nn.init.kaiming_normal_(conv.weight)
                
                feat_convs.append(conv)
            
            self.smoothing_convs[f'feature_{feat_idx}'] = feat_convs
        
        # Learnable mixing weights for multi-scale features
        self.scale_mixer = nn.Parameter(
            torch.ones(n_features, len(kernel_sizes)) / len(kernel_sizes)
        )
        
    def forward(self, x):
        """
        Args:
            x: (batch, seq_len, n_features)
        Returns:
            smoothed: (batch, seq_len, n_features)
        """
        batch_size, seq_len, _ = x.shape
        smoothed_features = []
        
        for feat_idx in range(self.n_features):
            # Extract single feature and reshape for Conv1d
            feat = x[:, :, feat_idx:feat_idx+1]  # (batch, seq_len, 1)
            feat = feat.transpose(1, 2)  # (batch, 1, seq_len)
            
            # Apply multi-scale smoothing
            multiscale_outputs = []
            for conv_idx, conv in enumerate(self.smoothing_convs[f'feature_{feat_idx}']):
                smoothed = conv(feat)
                # Trim to original sequence length (remove extra padding)
                kernel_size = self.kernel_sizes[conv_idx]
                smoothed = smoothed[:, :, :seq_len]
                multiscale_outputs.append(smoothed)
            
            # Stack and mix scales
            multiscale_stack = torch.stack(multiscale_outputs, dim=-1)  # (batch, 1, seq_len, n_scales)
            scale_weights = torch.softmax(self.scale_mixer[feat_idx], dim=0)
            mixed = torch.sum(multiscale_stack * scale_weights, dim=-1)  # (batch, 1, seq_len)
            
            mixed = mixed.transpose(1, 2)  # (batch, seq_len, 1)
            smoothed_features.append(mixed)
        
        return torch.cat(smoothed_features, dim=-1)
    
    def get_filter_weights(self):
        """Extract learned filter coefficients for analysis"""
        filters = {}
        for feat_idx in range(self.n_features):
            feat_filters = {}
            for scale_idx, conv in enumerate(self.smoothing_convs[f'feature_{feat_idx}']):
                kernel_size = self.kernel_sizes[scale_idx]
                weight = conv.weight.data.cpu().numpy().squeeze()
                feat_filters[f'kernel_{kernel_size}'] = weight
            filters[f'feature_{feat_idx}'] = feat_filters
        
        # Also return scale mixing weights
        filters['scale_weights'] = self.scale_mixer.data.cpu().numpy()
        
        return filters


class FiLM_LSTM(nn.Module):
    """FiLM LSTM with hybrid static feature processing"""

    def _initialize_weights(self):
        """Initialize weights with sensible defaults."""
        for name, param in self.named_parameters():
            if 'weight' in name:
                if 'lstm' in name:
                    # LSTM weights - use orthogonal initialization
                    if param.ndim >= 2:
                        nn.init.orthogonal_(param)
                elif param.ndim >= 2:
                    # Linear layers - use Xavier/Glorot initialization
                    nn.init.xavier_uniform_(param)
            elif 'bias' in name:
                nn.init.zeros_(param)

    
    def __init__(self, 
                 dynamic_input_size, 
                 static_input_size,
                 lstm_static_features=None,  # NEW: list of static feature names for LSTM
                 hidden_size=256, 
                 num_layers=1, 
                 dropout=0.2, 
                 use_time_modulation=False,
                 use_feature_smoothing=True,  # NEW
                 smoothing_kernel_sizes=None):  # NEW
        super().__init__()
        
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.use_time_modulation = use_time_modulation
        
        # NEW: Store which features go where
        self.lstm_static_features = lstm_static_features or []
        self.n_lstm_static = len(self.lstm_static_features)
        # NEW: Add learnable smoothing module
        self.use_feature_smoothing = use_feature_smoothing
        if use_feature_smoothing:
            self.feature_smoother = LearnableFeatureSmoothing(
                n_features=dynamic_input_size,
                kernel_sizes=smoothing_kernel_sizes or [3, 5, 7, 15, 31],
                initialize_uniform=True
            )
            print(f"Added learnable smoothing with kernels: {smoothing_kernel_sizes or [3, 5, 7, 15, 31]}")
        
        # LSTM now takes dynamic + selected static features
        lstm_input_size = dynamic_input_size + self.n_lstm_static
        
        self.lstm = nn.LSTM(
            input_size=lstm_input_size,  # Expanded input
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            dropout=dropout if num_layers > 1 else 0
        )
        
        # Static encoder for FiLM (uses ALL static features)
        self.static_encoder = nn.Sequential(
            nn.Linear(static_input_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
        )
        
        # Hidden initialization with all static features
        self.static_hidden_init = nn.Sequential(
            nn.Linear(static_input_size, hidden_size * num_layers),
            nn.GELU()
        )
        
        # Time Modulation Network (if enabled)
        if self.use_time_modulation:
            num_branches = 5
            base_channels = hidden_size // num_branches
            remainder = hidden_size % num_branches
            channel_sizes = [base_channels] * num_branches
            for i in range(remainder):
                channel_sizes[i] += 1
            
            self.time_modulator = nn.ModuleList([
                CausalConv1d(hidden_size, channel_sizes[0], kernel_size=3, dilation=1),
                CausalConv1d(hidden_size, channel_sizes[1], kernel_size=3, dilation=4),
                CausalConv1d(hidden_size, channel_sizes[2], kernel_size=3, dilation=16),
                CausalConv1d(hidden_size, channel_sizes[3], kernel_size=3, dilation=64),
                CausalConv1d(hidden_size, channel_sizes[4], kernel_size=3, dilation=182),
            ])
            
            self.temporal_norm = nn.LayerNorm(hidden_size)
            
            self.time_gate = nn.Sequential(
                nn.Linear(hidden_size, hidden_size),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size, hidden_size),
            )
        
        # Output layers
        self.output_layer = nn.Sequential(
            nn.Linear(hidden_size + static_input_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, 1)
        )
        
        self._initialize_weights()
    
    def forward(self, dynamic_inputs, static_inputs, static_feature_names=None):
        batch_size, seq_len, _ = dynamic_inputs.shape
        
        # NEW: Apply learnable smoothing to dynamic features
        if self.use_feature_smoothing:
            # Smooth the dynamic inputs
            dynamic_smoothed = self.feature_smoother(dynamic_inputs)
            
            # Option 1: Replace original with smoothed
            dynamic_inputs = dynamic_smoothed
            
            # Option 2: Concatenate original and smoothed (doubles feature dim)
            #dynamic_inputs = torch.cat([dynamic_inputs, dynamic_smoothed], dim=-1)
            
            # Option 3: Weighted combination (requires additional learnable weight)
            # dynamic_inputs = 0.5 * dynamic_inputs + 0.5 * dynamic_smoothed
        
        # NEW: Extract static features for LSTM if indices provided
        if self.n_lstm_static > 0 and static_feature_names is not None:
            # Find indices of LSTM static features
            lstm_static_indices = []
            for feat_name in self.lstm_static_features:
                if feat_name in static_feature_names:
                    idx = static_feature_names.index(feat_name)
                    lstm_static_indices.append(idx)
            
            if lstm_static_indices:
                # Extract and broadcast static features for LSTM
                lstm_static = static_inputs[:, lstm_static_indices]
                lstm_static_expanded = lstm_static.unsqueeze(1).expand(-1, seq_len, -1)
                
                # Concatenate with dynamic features
                lstm_inputs = torch.cat([dynamic_inputs, lstm_static_expanded], dim=-1)
            else:
                # Fallback if features not found
                lstm_inputs = dynamic_inputs
        else:
            lstm_inputs = dynamic_inputs
        
        # Initialize hidden state with ALL static features
        h_0 = self.static_hidden_init(static_inputs)
        # FIX: Correct reshaping with proper dimension ordering
        h_0 = h_0.view(batch_size, self.num_layers, self.hidden_size).transpose(0, 1)
        c_0 = torch.zeros_like(h_0)
        
        # Run LSTM with expanded inputs
        lstm_out, (h_n, c_n) = self.lstm(lstm_inputs, (h_0, c_0))
        
        # Generate static modulation for FiLM (using ALL static features)
        static_logits = self.static_encoder(static_inputs)
        
        # Apply modulation
        if self.use_time_modulation:
            lstm_out_reshaped = lstm_out.transpose(1, 2)
            
            multi_scale_features = []
            for conv in self.time_modulator:
                multi_scale_features.append(conv(lstm_out_reshaped))
            
            combined_features = torch.cat(multi_scale_features, dim=1)
            combined_features = combined_features.transpose(1, 2)
            combined_features = self.temporal_norm(combined_features)
            
            temporal_logits = self.time_gate(combined_features)
            combined_logits = static_logits.unsqueeze(1) + temporal_logits
            combined_modulation = torch.sigmoid(combined_logits)
            
            lstm_out = lstm_out * combined_modulation
        else:
            static_modulation = torch.sigmoid(static_logits)
            lstm_out = lstm_out * static_modulation.unsqueeze(1)
        
        # Take last timestep output
        last_output = lstm_out[:, -1, :]
        
        # Concatenate with ALL static features for final prediction
        combined = torch.cat([last_output, static_inputs], dim=1)
        prediction = self.output_layer(combined)
        
        return prediction




class StreamflowLightningModule(pl.LightningModule):
    """
    Lightning wrapper for the FiLM-LSTM streamflow model.

    Training/validation/test loss (normalized space):
        loss = MAE + RMSE + RMCE

    Reporting metrics (original streamflow scale, computed at epoch end):
        - KGE and components (r, alpha, beta)
        - NSE
        - RMSE
        - MAE
    """
    def __init__(
        self,
        dynamic_input_size: int,
        static_input_size: int,
        lstm_static_features: List[str] = None,
        static_feature_names: List[str] = None,
        hidden_size: int = 256,
        num_layers: int = 1,
        dropout: float = 0.2,
        learning_rate: float = 1e-4,
        use_time_modulation: bool = False,
        norm_params: Optional[Dict] = None, 
        target_noise_std: float = 0.0,       
        target_noise_type: str = 'additive', 
        gradient_clip_val: float = 5.0,      
        lr_scheduler_params: Optional[Dict] = None,  
        use_feature_smoothing: bool = True,  # NEW
        smoothing_kernel_sizes: List[int] = None,  # NEW
        **kwargs
    ):
        super().__init__()
        
        # Default to the three climate features if not specified
        if lstm_static_features is None:
            lstm_static_features = ['pre_mm_syr', 'pet_mm_syr', 'tmp_dc_syr', 'snw_pc_syr', 'area']
        
        self.lstm_static_features = lstm_static_features
        self.static_feature_names = static_feature_names
        
        # Save all hyperparameters
        self.save_hyperparameters()
        
        # Adjust dynamic input size for LSTM
        n_lstm_static = len(lstm_static_features)
        
       # Adjust input size if concatenating smoothed features
        effective_dynamic_size = dynamic_input_size
        #if use_feature_smoothing:
        #    # If concatenating, double the size
        #    effective_dynamic_size = dynamic_input_size * 2
        
        self.model = FiLM_LSTM(
            dynamic_input_size=dynamic_input_size,  # Original size for smoother
            static_input_size=static_input_size,
            lstm_static_features=lstm_static_features,
            hidden_size=hidden_size,
            num_layers=num_layers,
            dropout=dropout,
            use_time_modulation=use_time_modulation,
            use_feature_smoothing=use_feature_smoothing,  # NEW
            smoothing_kernel_sizes=smoothing_kernel_sizes  # NEW
        )

        # Target noise wrapper
        self.target_noise = TargetNoiseWrapper(
            noise_std=target_noise_std,
            noise_type=target_noise_type
        )

        # Store params
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

        # Epoch accumulators
        self.validation_step_outputs: List[Dict[str, torch.Tensor]] = []
        self.test_step_outputs: List[Dict[str, torch.Tensor]] = []

    # --------------------------
    # Loss helpers (normalized space)
    # --------------------------
    def _loss_components_norm(
        self,
        predictions: torch.Tensor,
        targets: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Compute MAE, RMSE, and RMCE in normalized space.
        Returns: (mae_norm, rmse_norm, rmce_norm)
        """
        err = predictions - targets
        mae_norm = torch.mean(torch.abs(err))
        rmse_norm = torch.sqrt(torch.mean(err ** 2))
        rmce_norm = torch.mean(torch.abs(err) ** 3) ** (1.0 / 3.0)
        return mae_norm, rmse_norm, rmce_norm

    def _combined_loss_norm(
        self,
        predictions: torch.Tensor,
        targets: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Combined loss used for optimization and monitoring:
            loss = mae_norm + rmse_norm + rmce_norm
        Returns: (loss, mae_norm, rmse_norm, rmce_norm)
        """
        mae_norm, rmse_norm, rmce_norm = self._loss_components_norm(predictions, targets)
        loss = mae_norm + rmse_norm + rmce_norm
        return loss, mae_norm, rmse_norm, rmce_norm

    # --------------------------
    # Lightning hooks
    # --------------------------
    def forward(self, dynamic_inputs: torch.Tensor, static_inputs: torch.Tensor) -> torch.Tensor:
        return self.model(dynamic_inputs, static_inputs, self.static_feature_names)

    def training_step(self, batch: Tuple[torch.Tensor, torch.Tensor, torch.Tensor], batch_idx: int) -> torch.Tensor:
        dynamic_seq, static_feat, target = batch
        predictions = self(dynamic_seq, static_feat)
        noisy_target = self.target_noise.add_noise(target, training=True)
        loss, mae_norm, rmse_norm, rmce_norm = self._combined_loss_norm(predictions, noisy_target)

        if torch.isnan(loss) or torch.isinf(loss):
            self.log('train_nan_count', 1.0, on_step=True, prog_bar=False, logger=True)
            return torch.zeros(1, device=loss.device, requires_grad=True)

        self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_mae_norm', mae_norm, on_step=False, on_epoch=True, logger=True)
        self.log('train_rmse_norm', rmse_norm, on_step=False, on_epoch=True, logger=True)
        self.log('train_rmce_norm', rmce_norm, on_step=False, on_epoch=True, logger=True)
        return loss

    def validation_step(self, batch: Tuple[torch.Tensor, torch.Tensor, torch.Tensor], batch_idx: int) -> torch.Tensor:
        """
        Validation step: compute combined loss (normalized space), store preds/targets
        for original-scale metrics at epoch end. Also log val_* on_epoch so EarlyStopping
        sees 'val_loss' during sanity check.
        """
        dynamic_seq, static_feat, target = batch
        predictions = self(dynamic_seq, static_feat)

        loss, mae_norm, rmse_norm, rmce_norm = self._combined_loss_norm(predictions, target)

        # Log on_epoch so Lightning aggregates across val batches (including sanity check)
        self.log('val_loss', loss, on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('val_mae_norm', mae_norm, on_step=False, on_epoch=True, logger=True)
        self.log('val_rmse_norm', rmse_norm, on_step=False, on_epoch=True, logger=True)
        self.log('val_rmce_norm', rmce_norm, on_step=False, on_epoch=True, logger=True)

        # Store for original-scale metrics at epoch end
        self.validation_step_outputs.append({
            'loss': loss.detach(),
            'mae_norm': mae_norm.detach(),
            'rmse_norm': rmse_norm.detach(),
            'rmce_norm': rmce_norm.detach(),
            'predictions': predictions.detach().cpu(),
            'targets': target.detach().cpu()
        })
        return loss

    def on_validation_epoch_end(self) -> None:
        """
        Aggregate validation metrics:
          - val_loss (combined loss in normalized space)
          - Original-scale KGE and components, NSE, RMSE, MAE
        """
        if not self.validation_step_outputs:
            return

        avg_loss = torch.stack([x['loss'] for x in self.validation_step_outputs]).mean()
        avg_mae_norm = torch.stack([x['mae_norm'] for x in self.validation_step_outputs]).mean()
        avg_rmse_norm = torch.stack([x['rmse_norm'] for x in self.validation_step_outputs]).mean()
        avg_rmce_norm = torch.stack([x['rmce_norm'] for x in self.validation_step_outputs]).mean()

        all_preds = torch.cat([x['predictions'] for x in self.validation_step_outputs])
        all_targets = torch.cat([x['targets'] for x in self.validation_step_outputs])

        preds_np = all_preds.numpy().flatten()
        targets_np = all_targets.numpy().flatten()

        if self.norm_params is not None:
            preds_original = unnormalize_predictions(preds_np, self.norm_params)
            targets_original = unnormalize_predictions(targets_np, self.norm_params)
        else:
            preds_original = preds_np
            targets_original = targets_np

        kge, r, alpha, beta = calculate_kge_components(targets_original, preds_original)
        rmse = calculate_rmse(targets_original, preds_original)
        mae = calculate_mae(targets_original, preds_original)
        nse = calculate_nse(targets_original, preds_original)

        # Log epoch-level metrics (val_loss is monitored by callbacks/scheduler)
        self.log('val_loss', avg_loss, on_epoch=True, prog_bar=True, logger=True)
        self.log('val_mae_norm', avg_mae_norm, on_epoch=True, logger=True)
        self.log('val_rmse_norm', avg_rmse_norm, on_epoch=True, logger=True)
        self.log('val_rmce_norm', avg_rmce_norm, on_epoch=True, logger=True)

        self.log('val_kge', kge, on_epoch=True, prog_bar=True, logger=True)
        self.log('val_kge_r', r, on_epoch=True, logger=True)
        self.log('val_kge_alpha', alpha, on_epoch=True, logger=True)
        self.log('val_kge_beta', beta, on_epoch=True, logger=True)
        self.log('val_rmse', rmse, on_epoch=True, logger=True)
        self.log('val_mae', mae, on_epoch=True, logger=True)
        self.log('val_nse', nse, on_epoch=True, logger=True)

        self.validation_step_outputs.clear()

    def test_step(self, batch: Tuple[torch.Tensor, torch.Tensor, torch.Tensor], batch_idx: int) -> torch.Tensor:
        dynamic_seq, static_feat, target = batch
        predictions = self(dynamic_seq, static_feat)
        loss, mae_norm, rmse_norm, rmce_norm = self._combined_loss_norm(predictions, target)

        self.test_step_outputs.append({
            'loss': loss.detach(),
            'mae_norm': mae_norm.detach(),
            'rmse_norm': rmse_norm.detach(),
            'rmce_norm': rmce_norm.detach(),
            'predictions': predictions.detach().cpu(),
            'targets': target.detach().cpu()
        })
        return loss

    def on_test_epoch_end(self) -> None:
        if not self.test_step_outputs:
            return

        avg_loss = torch.stack([x['loss'] for x in self.test_step_outputs]).mean()
        avg_mae_norm = torch.stack([x['mae_norm'] for x in self.test_step_outputs]).mean()
        avg_rmse_norm = torch.stack([x['rmse_norm'] for x in self.test_step_outputs]).mean()
        avg_rmce_norm = torch.stack([x['rmce_norm'] for x in self.test_step_outputs]).mean()

        all_preds = torch.cat([x['predictions'] for x in self.test_step_outputs]).numpy().flatten()
        all_targets = torch.cat([x['targets'] for x in self.test_step_outputs]).numpy().flatten()

        if self.norm_params is not None:
            preds_original = unnormalize_predictions(all_preds, self.norm_params)
            targets_original = unnormalize_predictions(all_targets, self.norm_params)
        else:
            preds_original = all_preds
            targets_original = all_targets

        kge, r, alpha, beta = calculate_kge_components(targets_original, preds_original)
        rmse = calculate_rmse(targets_original, preds_original)
        mae = calculate_mae(targets_original, preds_original)
        nse = calculate_nse(targets_original, preds_original)

        self.log('test_loss', avg_loss, on_epoch=True, logger=True)
        self.log('test_mae_norm', avg_mae_norm, on_epoch=True, logger=True)
        self.log('test_rmse_norm', avg_rmse_norm, on_epoch=True, logger=True)
        self.log('test_rmce_norm', avg_rmce_norm, on_epoch=True, logger=True)
        self.log('test_kge', kge, on_epoch=True, logger=True)
        self.log('test_kge_r', r, on_epoch=True, logger=True)
        self.log('test_kge_alpha', alpha, on_epoch=True, logger=True)
        self.log('test_kge_beta', beta, on_epoch=True, logger=True)
        self.log('test_rmse', rmse, on_epoch=True, logger=True)
        self.log('test_mae', mae, on_epoch=True, logger=True)
        self.log('test_nse', nse, on_epoch=True, logger=True)

        self.test_step_outputs.clear()

    def configure_optimizers(self):
        """Configure optimizer and learning rate scheduler."""
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)

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
                'monitor': 'val_loss',  # Combined loss (MAE+RMSE+RMCE) in normalized space
                'interval': 'epoch',
                'frequency': 1
            }
        }

    def on_after_backward(self) -> None:
        """Custom gradient diagnostics (optional)."""
        if self.global_step % 100 == 0:
            total_norm = 0.0
            max_grad = 0.0
            for p in self.model.parameters():
                if p.grad is not None:
                    param_norm = p.grad.data.norm(2)
                    total_norm += param_norm.item() ** 2
                    max_grad = max(max_grad, p.grad.data.abs().max().item())
            total_norm = total_norm ** 0.5
            if total_norm > 100 or max_grad > 50:
                self.log('large_gradient_warning', 1.0, on_step=True)
                print(f"Warning: Large gradients detected - norm: {total_norm:.2f}, max: {max_grad:.2f}")




class CaravanDataModule(pl.LightningDataModule):
    def __init__(
        self,
        watersheds_df: pd.DataFrame,
        bucket_name: str,
        base_data_dir: str,
        base_attr_dir: str,
        sequence_length: int = 365,
        batch_size: int = 256,
        num_workers: int = 16,
        chunk_size: int = 50,
        train_split: float = 0.6,
        val_split: float = 0.2,
        random_seed: int = 42,
        pin_memory: bool = True,
        persistent_workers: bool = True,
        prefetch_factor: int = 4,
        norm_params: Optional[Dict[str, Any]] = None,
        cache_data: bool = True,
        cache_dir: str = '/opt/ml/input/data/cache',
        max_missing_pct: float = 1.0,
        train_stride: int = 1,   # NEW
        val_stride: int = 1,    # NEW
        test_stride: int = 1,  # NEW
    ):
        super().__init__()
        # existing assignments...
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
        self.norm_params = norm_params
        self.cache_data = cache_data
        self.cache_dir = cache_dir
        self.max_missing_pct = max_missing_pct

        # NEW: strides
        self.train_stride = max(int(train_stride), 1)
        self.val_stride = max(int(val_stride), 1)
        self.test_stride = max(int(test_stride), 1)

        # Data attrs, caches...
        self.train_watersheds = None
        self.val_watersheds = None
        self.test_watersheds = None
        self.static_cols = None
        self.dynamic_cols_no_target = None
        self._cached_data = {}
        self._cached_datasets = {}
        self.current_train_chunk = 0
        self.train_chunks = []

        # Create cache directory if caching is enabled
        if self.cache_data and self.cache_dir:
            os.makedirs(self.cache_dir, exist_ok=True)
        
   
    def setup(self, stage: Optional[str] = None):
        """Prepare data splits and compute normalization parameters"""
        
        if stage == 'fit' or stage is None:
            # First, do a test load to identify which watersheds will be excluded
            print("Performing initial data quality check...")
            test_df, _, _, _ = prepare_ptf_dataframe(
                self.watersheds_df.iloc[:min(10, len(self.watersheds_df))],
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                compute_norm_params=False,
                max_missing_pct=self.max_missing_pct
            )
            
            # If a preassigned split is provided, honor it; else fall back to random splitting
            if 'split' in self.watersheds_df.columns:
                print("Using preassigned PCA-stratified splits from 'split' column.")
                self.train_watersheds = self.watersheds_df[self.watersheds_df['split'] == 'train'].reset_index(drop=True)
                temp_val = self.watersheds_df[self.watersheds_df['split'] == 'val'].reset_index(drop=True)
                temp_test = self.watersheds_df[self.watersheds_df['split'] == 'test'].reset_index(drop=True)
    
                self.val_watersheds = temp_val
                self.test_watersheds = temp_test
    
                print(f"Train watersheds: {len(self.train_watersheds)}")
                print(f"Val watersheds: {len(self.val_watersheds)}")
                print(f"Test watersheds: {len(self.test_watersheds)}")
            else:
                # Original behavior: random split (ignores order)
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
                    compute_norm_params=True,
                    max_missing_pct=self.max_missing_pct
                )
    
            # Get feature columns from a sample
            sample_df, static_cols, dynamic_cols, _ = prepare_ptf_dataframe(
                self.train_watersheds.iloc[:5],
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params,
                max_missing_pct=self.max_missing_pct
            )
            
            self.static_cols = static_cols
            self.dynamic_cols_no_target = [col for col in dynamic_cols if col != 'streamflow']
    
            # Pre-cache if enabled (unchanged)
            if self.cache_data and stage == 'fit':
                self._cache_all_data()
            else:
                self.train_chunks = [
                    self.train_watersheds[i:i+self.chunk_size]
                    for i in range(0, len(self.train_watersheds), self.chunk_size)
                ]
            
            del sample_df
            gc.collect()
    
        if stage == 'test':
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
                norm_params=self.norm_params,
                max_missing_pct=self.max_missing_pct
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
                norm_params=self.norm_params,
                max_missing_pct=self.max_missing_pct
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
        
        def _df_ok(df: pd.DataFrame) -> bool:
            return isinstance(df, pd.DataFrame) and (not df.empty) and ('group_id' in df.columns)
        
        train_df = self._cached_data.get('train')
        val_df = self._cached_data.get('val')
        
        if _df_ok(train_df):
            self._cached_datasets['train'] = OptimizedWatershedDataset(
                train_df,
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length,
                stride=self.train_stride,
            )
            print(f"✓ Training dataset created with {len(self._cached_datasets['train'])} samples")
        else:
            print("WARNING: Cached training DataFrame is empty or missing 'group_id'. Skipping pre-create of training dataset.")
        
        if _df_ok(val_df):
            self._cached_datasets['val'] = OptimizedWatershedDataset(
                val_df,
                self.static_cols,
                self.dynamic_cols_no_target,
                'streamflow',
                self.sequence_length,
                stride=self.val_stride,
            )
            print(f"✓ Validation dataset created with {len(self._cached_datasets['val'])} samples")
        else:
            print("WARNING: Cached validation DataFrame is empty or missing 'group_id'. Skipping pre-create of validation dataset.")
        
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
            # Create dataset from cached data if valid, else fall back
            train_df = self._cached_data['train']
            if isinstance(train_df, pd.DataFrame) and (not train_df.empty) and ('group_id' in train_df.columns):
                dataset = OptimizedWatershedDataset(
                    train_df,
                    self.static_cols,
                    self.dynamic_cols_no_target,
                    'streamflow',
                    self.sequence_length,
                    stride=self.train_stride,  # keep stride consistent
                )
            else:
                print("WARNING: Cached training DataFrame is empty or missing 'group_id'. Falling back to chunked loading.")
        
                # Fallback to original chunking approach
                attempts = 0
                dataset = None
                total_chunks = len(self.train_chunks)
                while attempts < total_chunks:
                    chunk_watersheds = self.train_chunks[self.current_train_chunk]
                    self.current_train_chunk = (self.current_train_chunk + 1) % total_chunks
                    attempts += 1
        
                    chunk_df, _, _, _ = prepare_ptf_dataframe(
                        chunk_watersheds,
                        self.bucket_name,
                        self.base_data_dir,
                        self.base_attr_dir,
                        norm_params=self.norm_params,
                        max_missing_pct=self.max_missing_pct
                    )
        
                    if isinstance(chunk_df, pd.DataFrame) and (not chunk_df.empty) and ('group_id' in chunk_df.columns):
                        candidate = WatershedDataset(
                            chunk_df,
                            self.static_cols,
                            self.dynamic_cols_no_target,
                            'streamflow',
                            self.sequence_length,
                            stride=self.train_stride
                        )
                        if len(candidate) > 0:
                            dataset = candidate
                            break
                    else:
                        print("Info: Skipping empty/malformed training chunk (no data or missing 'group_id').")
        
                if dataset is None:
                    raise RuntimeError(
                        "No non-empty training data available after filtering. "
                        "All chunks are empty or missing 'group_id'. "
                        "Consider relaxing --max-missing-pct or verifying input data."
                    )

        else:
            # Fallback to original chunking approach
            print("Warning: Using chunked loading (slower). Consider enabling caching.")
        
            attempts = 0
            dataset = None
            total_chunks = len(self.train_chunks)
            while attempts < total_chunks:
                chunk_watersheds = self.train_chunks[self.current_train_chunk]
                self.current_train_chunk = (self.current_train_chunk + 1) % total_chunks
                attempts += 1
        
                chunk_df, _, _, _ = prepare_ptf_dataframe(
                    chunk_watersheds,
                    self.bucket_name,
                    self.base_data_dir,
                    self.base_attr_dir,
                    norm_params=self.norm_params,
                    max_missing_pct=self.max_missing_pct
                )
        
                if isinstance(chunk_df, pd.DataFrame) and (not chunk_df.empty) and ('group_id' in chunk_df.columns):
                    candidate = WatershedDataset(
                        chunk_df,
                        self.static_cols,
                        self.dynamic_cols_no_target,
                        'streamflow',
                        self.sequence_length,
                        stride=self.train_stride
                    )
                    if len(candidate) > 0:
                        dataset = candidate
                        break
                else:
                    print("Info: Skipping empty/malformed training chunk (no data or missing 'group_id').")
        
            if dataset is None:
                raise RuntimeError(
                    "No non-empty training data available after filtering. "
                    "All chunks are empty or missing 'group_id'. "
                    "Consider relaxing --max-missing-pct or verifying input data."
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
        if self.cache_data and 'val' in self._cached_datasets:
            dataset = self._cached_datasets['val']
            return DataLoader(
                dataset,
                batch_size=self.batch_size,
                shuffle=False,
                num_workers=min(self.num_workers, 4),
                pin_memory=self.pin_memory and torch.cuda.is_available(),
                persistent_workers=self.persistent_workers and self.num_workers > 0,
                prefetch_factor=self.prefetch_factor if self.num_workers > 0 else 2,
            )
        # Fallback
        return self._create_eval_dataloader(self.val_watersheds, stride=self.val_stride)  # NEW
    
    def test_dataloader(self):
        return self._create_eval_dataloader(self.test_watersheds, stride=self.test_stride)  # NEW

    
    def _create_eval_dataloader(self, watersheds_df: pd.DataFrame, stride: int) -> DataLoader:
        chunks = [watersheds_df[i:i + self.chunk_size] for i in range(0, len(watersheds_df), self.chunk_size)]
        dataloaders = []
        for chunk in chunks:
            chunk_df, _, _, _ = prepare_ptf_dataframe(
                chunk,
                self.bucket_name,
                self.base_data_dir,
                self.base_attr_dir,
                norm_params=self.norm_params,
                max_missing_pct=self.max_missing_pct,
            )
            if isinstance(chunk_df, pd.DataFrame) and (not chunk_df.empty) and ('group_id' in chunk_df.columns):
                dataset = WatershedDataset(
                    chunk_df,
                    self.static_cols,
                    self.dynamic_cols_no_target,
                    'streamflow',
                    self.sequence_length,
                    stride=stride,
                )
                if len(dataset) > 0:
                    dataloaders.append(DataLoader(
                        dataset,
                        batch_size=self.batch_size,
                        shuffle=False,
                        num_workers=min(self.num_workers, 4),
                        pin_memory=self.pin_memory and torch.cuda.is_available(),
                        persistent_workers=self.persistent_workers and self.num_workers > 0,
                        prefetch_factor=self.prefetch_factor if self.num_workers > 0 else 2,
                    ))
            else:
                print("Info: Skipping empty/malformed eval chunk (no data or missing 'group_id').")
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

    def __init__(
        self,
        data_df: pd.DataFrame,
        static_cols: List[str],
        dynamic_cols: List[str],
        target_col: str,
        sequence_length: int,
        stride: int = 1,  # NEW
    ):
        self.data_df = data_df.copy()
        if 'group_id' not in self.data_df.columns:
            raise KeyError(
                "Expected 'group_id' column in data_df but it was missing. "
                "Upstream prepare_ptf_dataframe likely returned an empty/malformed DataFrame."
            )
        
        self.static_cols = static_cols
        self.dynamic_cols = dynamic_cols
        self.target_col = target_col
        self.sequence_length = sequence_length
        self.stride = max(int(stride), 1)  # NEW

        self.groups = self.data_df.groupby('group_id')
        self.group_keys = list(self.groups.groups.keys())

        print("Pre-computing dataset arrays for faster access...")
        self._precompute_all_samples()
        print(f"Dataset ready with {len(self)} samples")

    def _precompute_all_samples(self) -> None:
        self.all_samples: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        self.valid_indices: List[Tuple[str, int]] = []
        self.group_arrays: Dict[str, Dict[str, np.ndarray]] = {}

        for group_id in tqdm(self.group_keys, desc="Pre-processing watersheds"):
            group_data = self.groups.get_group(group_id)

            dynamic_array = group_data[self.dynamic_cols].values.astype(np.float32)
            static_array = group_data[self.static_cols].iloc[0].values.astype(np.float32)
            target_array = group_data[self.target_col].values.astype(np.float32)

            self.group_arrays[group_id] = {
                'dynamic': dynamic_array,
                'static': static_array,
                'target': target_array,
            }

            n_samples = len(group_data) - self.sequence_length
            if n_samples > 0:
                for i in range(0, n_samples, self.stride):  # NEW: use stride
                    self.valid_indices.append((group_id, i))

                    end_idx = i + self.sequence_length
                    dynamic_seq = dynamic_array[i:end_idx]
                    target = target_array[end_idx]

                    self.all_samples.append(
                        (dynamic_seq, static_array, np.array([target], dtype=np.float32))
                    )

    def __len__(self) -> int:
        return len(self.all_samples)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        dynamic_seq, static_features, target = self.all_samples[idx]
        return (
            torch.from_numpy(dynamic_seq),
            torch.from_numpy(static_features),
            torch.from_numpy(target),
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
    def __init__(self, initial_frequency=1, later_frequency=2, switch_epoch=20):
        self.initial_frequency = initial_frequency
        self.later_frequency = later_frequency
        self.switch_epoch = switch_epoch
        
    def on_train_epoch_start(self, trainer, pl_module):
        if trainer.current_epoch == self.switch_epoch:
            trainer.check_val_every_n_epoch = self.later_frequency
            print(f"Switching validation frequency to every {self.later_frequency} epochs")

class LossHistoryCallback(pl.Callback):
    """
    Collect per-epoch history for training/validation loss and LR.
    Uses trainer.callback_metrics to gather values at the end of each validation epoch.
    """
    def __init__(self):
        super().__init__()
        self.history = {
            'epoch': [],
            'train_loss': [],
            'val_loss': [],
            'val_kge': [],
            'learning_rate': []
        }

    def on_validation_epoch_end(self, trainer, pl_module):
        epoch = int(trainer.current_epoch)
        cm = trainer.callback_metrics  # dict-like

        # Extract metrics if present
        val_loss = cm.get('val_loss', None)
        val_kge = cm.get('val_kge', None)
        train_loss = cm.get('train_loss_epoch', cm.get('train_loss', None))

        # Current LR from optimizer
        try:
            lr = float(trainer.optimizers[0].param_groups[0]['lr'])
        except Exception:
            lr = np.nan

        self.history['epoch'].append(epoch)
        self.history['val_loss'].append(float(val_loss.item()) if val_loss is not None else np.nan)
        self.history['val_kge'].append(float(val_kge.item()) if val_kge is not None else np.nan)
        self.history['train_loss'].append(float(train_loss.item()) if train_loss is not None and hasattr(train_loss, 'item') else (float(train_loss) if train_loss is not None else np.nan))
        self.history['learning_rate'].append(lr)



class TargetNoiseWrapper:
    """Wrapper to add noise to targets during training (for log-normalized targets)"""
    def __init__(self, noise_std=0.1, noise_type='additive'):
        self.noise_std = noise_std
        self.noise_type = noise_type
        
    def add_noise(self, target, training=True):
        """Add noise to target values during training"""
        if not training or self.noise_std == 0:
            return target
            
        if self.noise_type == 'additive':
            # Simple additive Gaussian noise in normalized log space
            noise = torch.randn_like(target) * self.noise_std
            noisy_target = target + noise
            
        elif self.noise_type == 'multiplicative':
            # Multiplicative noise (scales with magnitude)
            noise = torch.randn_like(target) * self.noise_std
            noisy_target = target * (1 + noise)
            
        # No clamping needed in log-normalized space
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



def save_loss_history(history: Dict[str, list], out_dir: str):
    """
    Save loss history to CSV and plot using existing plot_training_history_with_lr.
    """
    os.makedirs(out_dir, exist_ok=True)
    # Save CSV
    hist_df = pd.DataFrame(history)
    csv_path = os.path.join(out_dir, 'training_history.csv')
    hist_df.to_csv(csv_path, index=False)
    print(f"Saved training history CSV to {csv_path}")
    # Save plots
    try:
        plot_training_history_with_lr(history, out_dir)
        print(f"Saved training history plots to {out_dir}")
    except Exception as e:
        print(f"Warning: Failed to plot training history: {e}")









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

def _ensure_output_files_dir(model_dir: Optional[str]) -> Optional[str]:
    if not model_dir:
        return None
    out_dir = os.path.join(model_dir, 'output_files')
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

def copy_to_output_files(src_path: str, model_dir: Optional[str], dst_name: Optional[str] = None):
    """
    Copy a file or directory into model_dir/output_files (SageMaker will upload this to S3).
    """
    out_dir = _ensure_output_files_dir(model_dir)
    if out_dir is None:
        return
    if not os.path.exists(src_path):
        return
    dst_path = os.path.join(out_dir, dst_name if dst_name else os.path.basename(src_path))
    if os.path.isdir(src_path):
        # Copy a directory tree (merge if exists)
        if os.path.exists(dst_path):
            shutil.rmtree(dst_path)
        shutil.copytree(src_path, dst_path)
    else:
        shutil.copy2(src_path, dst_path)
    print(f"Copied to output_files: {dst_path}")

class TotalEpochTracker(pl.Callback):
    """Track total epochs across multiple training runs"""
    
    def __init__(self, initial_epoch=0, checkpoint_dir='/opt/ml/checkpoints', s3_uri=None):
        self.total_epochs = initial_epoch
        self.checkpoint_dir = checkpoint_dir
        self.s3_uri = s3_uri
        self.start_time = time.time()
    
    def on_train_epoch_end(self, trainer, pl_module):
        self.total_epochs += 1
        
        # Log with proper epoch number
        pl_module.log('total_epochs_trained', self.total_epochs)
        
        # Save state periodically
        if trainer.current_epoch % 5 == 0:
            self._save_state()
    
    def on_train_end(self, trainer, pl_module):
        self._save_state()
        
        # Final sync to S3
        if self.s3_uri:
            import subprocess
            subprocess.run([
                'aws', 's3', 'cp',
                os.path.join(self.checkpoint_dir, 'training_state.json'),
                f"{self.s3_uri}/training_state.json"
            ])
    
    def _save_state(self):
        state = {
            'completed_epochs': self.total_epochs,
            'timestamp': datetime.now().isoformat(),
            'runtime_hours': (time.time() - self.start_time) / 3600
        }
        
        state_file = os.path.join(self.checkpoint_dir, 'training_state.json')
        with open(state_file, 'w') as f:
            json.dump(state, f, indent=2)

def find_best_checkpoint(checkpoint_dir, epoch_offset=0):
    """Find the best checkpoint considering epoch offset"""
    
    # First try last.ckpt
    last_ckpt = os.path.join(checkpoint_dir, 'last.ckpt')
    if os.path.exists(last_ckpt):
        print(f"Resuming from last checkpoint: {last_ckpt}")
        return last_ckpt
    
    # Look for checkpoint with highest epoch number
    import glob
    pattern = os.path.join(checkpoint_dir, '*.ckpt')
    checkpoints = glob.glob(pattern)
    
    if checkpoints:
        # Sort by modification time and take the most recent
        checkpoints.sort(key=os.path.getmtime, reverse=True)
        print(f"Resuming from checkpoint: {checkpoints[0]}")
        return checkpoints[0]
    
    print("No checkpoint found, starting fresh")
    return None



def monitor_and_continue_training(initial_job_name, max_total_days=30):
    """Monitor training and automatically continue when hitting time limits"""
    
    sagemaker_client = boto3.client('sagemaker')
    total_days_trained = 0
    current_job = initial_job_name
    
    while total_days_trained < max_total_days:
        # Monitor current job
        while True:
            response = sagemaker_client.describe_training_job(
                TrainingJobName=current_job
            )
            status = response['TrainingJobStatus']
            
            if status == 'Completed':
                print(f"Training completed successfully!")
                return
            elif status in ['Failed', 'Stopped']:
                # Check if it's a MaxRuntimeExceeded failure
                if 'MaxRuntimeExceeded' in response.get('FailureReason', ''):
                    print(f"Job {current_job} hit time limit. Starting continuation...")
                    break
                else:
                    print(f"Job failed: {response.get('FailureReason')}")
                    return
            
            time.sleep(300)  # Check every 5 minutes
        
        # Start continuation job
        total_days_trained += 5
        if total_days_trained < max_total_days:
            # Launch new training job with same configuration
            # (implement continuation logic here similar to section 1)
            pass


def pca_stratified_reorder_watersheds(
    watershed_df: pd.DataFrame,
    bucket_name: str,
    base_attr_dir: str,
    grid_size: int = 5,
    train_frac: float = 0.6,
    val_frac: float = 0.2,
    test_frac: float = 0.2,
    random_state: int = 42,
    subdirectory_col: str = 'subdirectory_name',
    id_col: str = 'watershedID',
    log_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Perform PCA on HydroATLAS static attributes and create a stratified train/val/test split
    across a grid in the first two PCs. Returns a reshuffled DataFrame with ['pc1','pc2','stratum','split'].

    - Loads attributes_hydroatlas_{subdirectory}.csv once from S3/local.
    - Uses the default HydroATLAS columns (intersection with available).
    - Imputes missing values with per-feature median.
    - Standardizes features, computes PCA (2 components).
    - Bins PC1 and PC2 into quantile bins (grid_size each) and defines strata.
    - Within each stratum, shuffles and allocates train/val/test by target fractions.
    - Returns a reshuffled watershed_df and logs the order.

    Note: CaravanDataModule.setup must honor the 'split' column to use this stratification.
    """

    rng = np.random.default_rng(random_state)
    df = watershed_df.copy()

    # Build gauge_id and ensure consistent dtype
    df['gauge_id'] = df[subdirectory_col].astype(str) + '_' + df[id_col].astype(str)
    df['gauge_id'] = df['gauge_id'].str.lower()  # Normalize to lowercase

    # Default HydroATLAS columns (same list used in get_watershed_attributes)
    default_hydroatlas_cols = [
        'pet_mm_syr',
        'aet_mm_syr',  # optional aet_mm_s01,...,12
        'pre_mm_syr',
        'tmp_dc_syr', 'tmp_dc_smn', 'tmp_dc_smx',
        'snw_pc_syr', 'snw_pc_smx',
        'swc_pc_syr',
#        'gdp_ud_sav',
        'inu_pc_slt', 'inu_pc_smn', 'inu_pc_smx',  # inundation percent
        'run_mm_syr',
        'dis_m3_pmn', 'dis_m3_pmx', 'dis_m3_pyr',
        'lka_pc_sse',  # limnicity pct area
        'lkv_mc_usu',  # lake vol
        'rev_mc_usu',  # reservoir volume
#        'ria_ha_usu',  # river area
        'riv_tc_usu',  # river volume
        'dor_pc_pva',  # degree of regulation
        'gwt_cm_sav',  # groundwater table depth
        'ele_mt_sav', 'ele_mt_smn', 'ele_mt_mx',
#        'slp_dg_sav',  # terrain slope
        'sgr_dk_sav',  # stream gradient
        'ari_ix_sav',
        'glc_cl_smj',  # land cover classes
        #'glc_cl_s01','glc_cl_s02','glc_cl_s03','glc_cl_s04','glc_cl_s05','glc_cl_s06','glc_cl_s07','glc_cl_s08','glc_cl_s09','glc_cl_s10','glc_cl_s11',
        #'glc_cl_s12','glc_cl_s13','glc_cl_s14','glc_cl_s15','glc_cl_s16','glc_cl_s17','glc_cl_s18','glc_cl_s19','glc_cl_s20','glc_cl_s21','glc_cl_s22',
        'for_pc_sse','pst_pc_sse',
        'crp_pc_sse', 'ire_pc_sse', 'gla_pc_sse', 'prm_pc_sse',
        'cly_pc_sav', 'slt_pc_sav', 'snd_pc_sav',
        'soc_th_sav', 'lit_cl_smj', 'kar_pc_sse',
#        'rdd_mk_sav',
        'urb_pc_sse',
        'hdi_ix_sav', 'nli_ix_sav', 'ppd_pk_sav',  # human development index, nighttime lights, pop density
        'cls_cl_smj', 'clz_cl_smj',  # climate strata and zones
    ]
    
    # We assume df is already filtered to a single subdirectory (e.g., camels)
    unique_subdirs = df[subdirectory_col].unique().tolist()
    all_parts = []

    for subdir in unique_subdirs:
        sub_df = df[df[subdirectory_col] == subdir].copy()
        if sub_df.empty:
            continue

        # Resolve path to attributes_hydroatlas_{subdir}.csv
        if base_attr_dir.startswith('/'):
            attr_path = f"{base_attr_dir.rstrip('/')}/{subdir}/attributes_hydroatlas_{subdir}.csv"
        else:
            attr_path = f"s3://{bucket_name}/{base_attr_dir.strip('/')}/{subdir}/attributes_hydroatlas_{subdir}.csv"

        print(f"[PCA Split] Loading attributes from: {attr_path}")
        attrs = pd.read_csv(attr_path, low_memory=False)

        if 'gauge_id' not in attrs.columns:
            raise ValueError(f"gauge_id not found in {attr_path}")

       # FIRST: Normalize CSV gauge_ids to lowercase
        attrs['gauge_id'] = attrs['gauge_id'].str.lower()
        
        # THEN: Filter to selected gauges using lowercase comparison
        attrs = attrs[attrs['gauge_id'].isin(sub_df['gauge_id'])]
        
        if attrs.empty:
            # Debug output to help diagnose
            print(f"[DEBUG] No matches found!")
            print(f"[DEBUG] First 5 gauge_ids we're looking for: {sub_df['gauge_id'].head().tolist()}")
            print(f"[DEBUG] First 5 gauge_ids in CSV (after lowercase): {attrs['gauge_id'].head().tolist()}")
            raise ValueError(f"No matching gauge_ids found in attributes for subdirectory: {subdir}")
    
        # Select available HydroATLAS columns (intersection)
        avail_cols = [c for c in default_hydroatlas_cols if c in attrs.columns]
        if not avail_cols:
            raise ValueError(f"No default HydroATLAS columns found in {attr_path}")

        X = attrs[avail_cols].copy()

        # Impute per-feature median
        for c in avail_cols:
            median_val = X[c].median()
            X[c] = X[c].fillna(median_val)

        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X.values)

        # PCA (2 components)
        pca = PCA(n_components=2, random_state=random_state)
        pcs = pca.fit_transform(X_scaled)
        attrs['pc1'] = pcs[:, 0]
        attrs['pc2'] = pcs[:, 1]

        # Merge PCs back to sub_df
        sub_df = sub_df.merge(attrs[['gauge_id', 'pc1', 'pc2']], on='gauge_id', how='left')

        # Bin PC1/PC2 into quantile bins
        def quantile_bins(series: pd.Series, bins: int) -> Tuple[pd.Series, int]:
            # Use qcut with duplicates='drop' to be robust to ties
            try:
                binned = pd.qcut(series, q=bins, labels=False, duplicates='drop')
                n_bins_eff = binned.nunique()
                binned = binned.astype(int)
                return binned, n_bins_eff
            except Exception:
                # Fallback: rank-based uniform bins
                ranks = series.rank(method='average', pct=True)
                binned = np.minimum((ranks * bins).astype(int), bins - 1)
                return pd.Series(binned.values, index=series.index), bins

        sub_df['pc1_bin'], n1 = quantile_bins(sub_df['pc1'], grid_size)
        sub_df['pc2_bin'], n2 = quantile_bins(sub_df['pc2'], grid_size)
        sub_df['stratum'] = sub_df['pc1_bin'] * n2 + sub_df['pc2_bin']

        # Allocate splits per stratum
        split_labels = []
        for stratum, grp in sub_df.groupby('stratum', sort=False):
            idx = grp.index.to_numpy()
            # Shuffle within stratum
            rng.shuffle(idx)
            n = len(idx)

            # Integer allocation by rounding, then adjust to sum n
            n_train = int(round(n * train_frac))
            n_val = int(round(n * val_frac))
            n_test = n - n_train - n_val

            # Fix edge cases to avoid negatives
            if n_test < 0:
                # Reduce validation first
                reduce = -n_test
                take_from_val = min(reduce, n_val)
                n_val -= take_from_val
                reduce -= take_from_val
                n_train = max(0, n_train - reduce)
                n_test = n - n_train - n_val

            # Assign
            split = pd.Series(index=idx, dtype=object)
            split.iloc[:n_train] = 'train'
            split.iloc[n_train:n_train + n_val] = 'val'
            split.iloc[n_train + n_val:] = 'test'
            split_labels.append(split)

        sub_df['split'] = pd.concat(split_labels).sort_index()

        # Reorder rows: train first, then val, then test; shuffle within split for randomness
        def reorder_within_split(frame: pd.DataFrame) -> pd.DataFrame:
            parts = []
            for split_name in ['train', 'val', 'test']:
                part = frame[frame['split'] == split_name].sample(frac=1.0, random_state=random_state)
                parts.append(part)
            return pd.concat(parts, axis=0)

        sub_df = reorder_within_split(sub_df)

        # Keep relevant columns and append
        all_parts.append(sub_df)

    # Combine all subdirectories (if more than one)
    out_df = pd.concat(all_parts, axis=0).reset_index(drop=True)

    # Log reshuffled order and per-split counts
    print("[PCA Split] Final split counts:")
    print(out_df['split'].value_counts(dropna=False).to_dict())

    # Optional: save a log CSV with the new order and PCs
    if log_path is not None:
        try:
            cols_to_save = [
                subdirectory_col, id_col, 'gauge_id', 'pc1', 'pc2', 'pc1_bin', 'pc2_bin', 'stratum', 'split'
            ]
            out_df[cols_to_save].to_csv(log_path, index=False)
            print(f"[PCA Split] Saved reshuffled order log to {log_path}")
        except Exception as e:
            print(f"[PCA Split] Failed to save reshuffle log: {e}")

    # Return only the columns needed downstream plus the split/PCA info
    # Preserve original columns and add new ones
    reorder_cols = watershed_df.columns.tolist() + ['pc1', 'pc2', 'stratum', 'split', 'gauge_id', 'pc1_bin', 'pc2_bin']
    reorder_cols = [c for c in reorder_cols if c in out_df.columns]
    out_df = out_df[reorder_cols]

    return out_df





def main():
    """Main entry point with continuation support for spot instances and 5-day limits"""
    
    # Parse arguments first to get access to them throughout
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
    parser.add_argument('--pca-grid-size', type=int, default=5)
    parser.add_argument('--split-seed', type=int, default=42)
    parser.add_argument('--max-missing-pct', type=float, default=1.0)
    parser.add_argument('--train-stride', type=int, default=1)
    parser.add_argument('--val-stride', type=int, default=30)
    parser.add_argument('--test-stride', type=int, default=365)
    parser.add_argument('--analysis-batch-size', type=int, default=256)
    parser.add_argument('--analysis-use-cache', action='store_true', default=False)
    parser.add_argument('--analysis-cache-dir', type=str, default=None)
    
    # Continuation-specific arguments
    parser.add_argument('--completed-epochs', type=int, default=0,
                        help='Number of epochs completed in previous runs')
    parser.add_argument('--is-continuation', type=str, default='False',
                        help='Whether this is a continuation run')
    parser.add_argument('--checkpoint-s3-uri', type=str, default=None,
                        help='S3 URI for checkpoint storage')
    
    args = parser.parse_args()
    
    # Start time for tracking wall time
    start_time = time.time()
    
    # Checkpoint management with S3 sync
    checkpoint_dir = os.environ.get('SM_CHECKPOINT_DIR', '/opt/ml/checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Download checkpoints from S3 if continuation
    if args.is_continuation == 'True' and args.checkpoint_s3_uri:
        print(f"Downloading checkpoints from {args.checkpoint_s3_uri}")
        import subprocess
        try:
            subprocess.run([
                'aws', 's3', 'sync', 
                args.checkpoint_s3_uri, 
                checkpoint_dir,
                '--exclude', '*',
                '--include', '*.ckpt',
                '--include', 'training_state.json'
            ], check=True)
            print(f"Downloaded checkpoints to {checkpoint_dir}")
        except Exception as e:
            print(f"Warning: Failed to download checkpoints: {e}")
    
    # Load training state
    training_state_file = os.path.join(checkpoint_dir, 'training_state.json')
    if os.path.exists(training_state_file):
        with open(training_state_file, 'r') as f:
            training_state = json.load(f)
        print(f"Loaded training state: {training_state}")
    else:
        training_state = {
            'completed_epochs': args.completed_epochs,
            'best_metric': float('inf'),
            'total_wall_time': 0,
            'last_learning_rate': args.learning_rate
        }
    
    # Calculate remaining epochs
    actual_completed_epochs = training_state['completed_epochs']
    remaining_epochs = args.total_epochs - actual_completed_epochs
    
    if remaining_epochs <= 0:
        print(f"Training already complete! {actual_completed_epochs}/{args.total_epochs} epochs done.")
        # Sync final checkpoints to S3 if needed
        if args.checkpoint_s3_uri:
            subprocess.run([
                'aws', 's3', 'sync',
                checkpoint_dir,
                args.checkpoint_s3_uri
            ])
        sys.exit(0)
    
    print(f"Continuing training: {actual_completed_epochs}/{args.total_epochs} epochs completed")
    print(f"Will train for {remaining_epochs} more epochs")
    
    # Signal handler for spot instance interruption
    def signal_handler(signum, frame):
        """Handle spot instance interruption gracefully"""
        print(f"\nReceived signal {signum}. Attempting to save checkpoint and state...")
        try:
            if 'trainer' in locals() and trainer is not None:
                # Save checkpoint
                checkpoint_path = os.path.join(checkpoint_dir, "interrupted.ckpt")
                trainer.save_checkpoint(checkpoint_path)
                print(f"Checkpoint saved to {checkpoint_path}")
                
                # Update and save training state
                current_training_state = {
                    'completed_epochs': actual_completed_epochs + trainer.current_epoch,
                    'best_metric': float(trainer.checkpoint_callback.best_model_score) 
                                   if trainer.checkpoint_callback else training_state.get('best_metric'),
                    'total_wall_time': training_state.get('total_wall_time', 0) + time.time() - start_time,
                    'last_learning_rate': float(trainer.optimizers[0].param_groups[0]['lr']),
                    'last_updated': datetime.now().isoformat()
                }
                
                with open(training_state_file, 'w') as f:
                    json.dump(current_training_state, f, indent=2)
                
                # Sync to S3
                if args.checkpoint_s3_uri:
                    subprocess.run([
                        'aws', 's3', 'sync',
                        checkpoint_dir,
                        args.checkpoint_s3_uri
                    ])
                    print(f"Checkpoints synced to {args.checkpoint_s3_uri}")
                    
        except Exception as e:
            print(f"Failed to save checkpoint: {e}")
        sys.exit(0)
    
    # Register signal handlers
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    # CUDA Memory Management
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.set_per_process_memory_fraction(0.95)
        torch.cuda.reset_peak_memory_stats()
        os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'
    
    gc.collect()
    
    # Enable CUDNN optimizations
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    
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
    
    # Create output directories
    os.makedirs(args.output_data_dir, exist_ok=True)
    if args.model_dir:
        os.makedirs(args.model_dir, exist_ok=True)
    
    # Load watershed data
    print("Loading watershed data...")
    watershed_df = identify_all_available_watersheds(args.bucket_name, args.base_data_dir)
    # subset the watersheds to camels, camelsbr, grdc, hysets, etc.
    watershed_df = watershed_df[watershed_df['subdirectory_name'] != 'grdc']
    print(f"Found {len(watershed_df)} watersheds")
    
    # Apply PCA-based stratified split
    print("Applying PCA-based stratified split and reshuffle...")
    reshuffle_log = os.path.join(args.output_data_dir, f"{args.experiment_name}_pca_split_log.csv")
    watershed_df = pca_stratified_reorder_watersheds(
        watershed_df=watershed_df,
        bucket_name=args.bucket_name,
        base_attr_dir=args.base_attr_dir,
        grid_size=args.pca_grid_size,
        train_frac=0.6,
        val_frac=0.2,
        test_frac=0.2,
        random_state=args.split_seed,
        subdirectory_col='subdirectory_name',
        id_col='watershedID',
        log_path=reshuffle_log
    )
    print("PCA-based reshuffle complete.")
    
    # Copy split log to model_dir/output_files
    if os.path.exists(reshuffle_log):
        copy_to_output_files(reshuffle_log, args.model_dir, dst_name=os.path.basename(reshuffle_log))
    
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
        random_seed=42,
        max_missing_pct=args.max_missing_pct,
        train_stride=args.train_stride,
        val_stride=args.val_stride,
        test_stride=args.test_stride
    )
    
    # Setup data
    print("Setting up data splits and computing normalization parameters...")
    data_module.setup('fit')
    
    # Validate normalization params
    if data_module.norm_params is None or not isinstance(data_module.norm_params, dict):
        raise RuntimeError(
            "\n" + "="*60 + "\n"
            "CRITICAL ERROR: Failed to compute normalization parameters!\n"
            "This typically means all training watersheds were excluded due to:\n"
            "  1. More than 1% missing streamflow data\n"
            "  2. No valid streamflow data after trimming\n"
            "  3. Data loading errors\n"
            "\n"
            "Debugging steps:\n"
            "  - Try increasing --max-missing-pct (e.g., 2.0 or 5.0) temporarily\n"
            "  - Check S3 paths and permissions\n"
            "  - Verify the grdc subset has valid .nc files\n"
            "  - Look for 'WARNING' and 'ERROR' messages above\n"
            + "="*60
        )
    
    if not hasattr(data_module, 'train_watersheds') or len(data_module.train_watersheds) == 0:
        raise RuntimeError("No training watersheds available after splitting")
    
    print(f"Normalization parameters computed successfully")
    print(f"  Target mean: {data_module.norm_params.get('target_mean', 'N/A')}")
    print(f"  Target std: {data_module.norm_params.get('target_std', 'N/A')}")
    
    # Save normalization parameters
    norm_params_path = os.path.join(args.output_data_dir, 'normalization_params.json')
    with open(norm_params_path, 'w') as f:
        json_params = {
            'static_means': {k: float(v) for k, v in data_module.norm_params.get('static_means', {}).items()},
            'static_stds': {k: float(v) for k, v in data_module.norm_params.get('static_stds', {}).items()},
            'dynamic_means': {k: float(v) for k, v in data_module.norm_params.get('dynamic_means', {}).items()},
            'dynamic_stds': {k: float(v) for k, v in data_module.norm_params.get('dynamic_stds', {}).items()},
            'target_mean': float(data_module.norm_params.get('target_mean', 0.0)),
            'target_std': float(data_module.norm_params.get('target_std', 1.0)),
            'dynamic_mins': {k: float(v) for k, v in data_module.norm_params.get('dynamic_mins', {}).items()},
            'dynamic_maxs': {k: float(v) for k, v in data_module.norm_params.get('dynamic_maxs', {}).items()},
            'std_only_cols': data_module.norm_params.get('std_only_cols', []),
            'minmax_cols': data_module.norm_params.get('minmax_cols', []),
        }
        json.dump(json_params, f, indent=2)
    
    print(f"Saved normalization parameters to {norm_params_path}")
    copy_to_output_files(norm_params_path, args.model_dir, dst_name='normalization_params.json')
    
    # Define static features for LSTM
    lstm_static_features = ['pre_mm_syr', 'pet_mm_syr', 'tmp_dc_syr']
    
    missing_features = [f for f in lstm_static_features if f not in data_module.static_cols]
    if missing_features:
        print(f"WARNING: LSTM static features {missing_features} not found in static columns")
        print(f"Available static features: {data_module.static_cols[:20]}...")
        lstm_static_features = []
    
    # Create or load model
    print("Creating Lightning model...")
    
    # Adjust learning rate if continuation (can decay or use last LR)
    initial_lr = training_state.get('last_learning_rate', args.learning_rate)
    
    model = StreamflowLightningModule(
        dynamic_input_size=len(data_module.dynamic_cols_no_target),
        static_input_size=len(data_module.static_cols),
        lstm_static_features=lstm_static_features,
        static_feature_names=data_module.static_cols,
        hidden_size=args.hidden_size,
        num_layers=args.num_layers,
        dropout=args.dropout,
        learning_rate=initial_lr,  # Use last LR or initial
        norm_params=data_module.norm_params,
        use_time_modulation=False,
        target_noise_std=0.0,
        target_noise_type='additive',
        gradient_clip_val=5.0,
        use_feature_smoothing=True,
        smoothing_kernel_sizes=[1,2,3,5,7,30],
        lr_scheduler_params={
            'factor': 0.5,
            'patience': 2,
            'min_lr': 1e-7,
            'threshold': 0.001,
            'cooldown': 1,
            'verbose': True
        }
    )
    
    # Custom callback for tracking total epochs
    class TotalEpochTracker(pl.Callback):
        """Track total epochs across multiple training runs"""
        
        def __init__(self, initial_epoch=0, checkpoint_dir='/opt/ml/checkpoints', s3_uri=None):
            self.total_epochs = initial_epoch
            self.checkpoint_dir = checkpoint_dir
            self.s3_uri = s3_uri
            self.start_time = time.time()
        
        def on_train_epoch_end(self, trainer, pl_module):
            self.total_epochs += 1
            
            # Log with proper epoch number
            pl_module.log('total_epochs_trained', self.total_epochs)
            
            # Save state periodically
            if trainer.current_epoch % 5 == 0:
                self._save_state()
        
        def on_train_end(self, trainer, pl_module):
            self._save_state()
            
            # Final sync to S3
            if self.s3_uri:
                import subprocess
                subprocess.run([
                    'aws', 's3', 'sync',
                    self.checkpoint_dir,
                    self.s3_uri
                ])
        
        def _save_state(self):
            state = {
                'completed_epochs': self.total_epochs,
                'timestamp': datetime.now().isoformat(),
                'runtime_hours': (time.time() - self.start_time) / 3600
            }
            
            state_file = os.path.join(self.checkpoint_dir, 'training_state.json')
            with open(state_file, 'w') as f:
                json.dump(state, f, indent=2)
    
    # Enhanced checkpoint callback with epoch offset
    class ContinuationModelCheckpoint(pl.callbacks.ModelCheckpoint):
        def __init__(self, epoch_offset=0, **kwargs):
            super().__init__(**kwargs)
            self.epoch_offset = epoch_offset
        
        def on_save_checkpoint(self, trainer, pl_module, checkpoint):
            # Add continuation metadata
            checkpoint['continuation_metadata'] = {
                'total_epochs_trained': self.epoch_offset + trainer.current_epoch,
                'continuation_run': True,
                'timestamp': datetime.now().isoformat()
            }
            return checkpoint
    
    # Create callbacks
    loss_history_cb = LossHistoryCallback()
    
    # Create trainer with continuation support
    trainer = pl.Trainer(
        max_epochs=remaining_epochs,  # Only train remaining epochs
        accelerator='gpu' if device == 'cuda' else 'cpu',
        devices=1,
        precision=32,  # Full precision for stability
        gradient_clip_val=1.0,  # More aggressive clipping
        accumulate_grad_batches=1,  # No accumulation to reduce instability
        num_sanity_val_steps=2,
        
        callbacks=[
            ValidationFrequencyCallback(initial_frequency=1, later_frequency=2, switch_epoch=20),
            ContinuationModelCheckpoint(
                epoch_offset=actual_completed_epochs,
                dirpath=checkpoint_dir,
                filename=f'{args.experiment_name}_total_epoch_{{epoch:03d}}_loss_{{val_loss:.4f}}',
                monitor='val_loss',
                mode='min',
                save_top_k=3,
                save_last=True,
                save_on_train_epoch_end=True,
                every_n_epochs=1,
                auto_insert_metric_name=False
            ),
            pl.callbacks.EarlyStopping(
                monitor='val_loss',
                patience=40,
                mode='min',
                verbose=True
            ),
            pl.callbacks.LearningRateMonitor(logging_interval='epoch'),
            pl.callbacks.TQDMProgressBar(refresh_rate=50),
            loss_history_cb,
            TotalEpochTracker(
                initial_epoch=actual_completed_epochs,
                checkpoint_dir=checkpoint_dir,
                s3_uri=args.checkpoint_s3_uri
            ),
        ],
        
        logger=pl.loggers.TensorBoardLogger(
            save_dir=args.output_data_dir,
            name=args.experiment_name,
            version=f'lightning_run_{actual_completed_epochs}'
        ),
        
        check_val_every_n_epoch=1,
        limit_val_batches=1.0,
        enable_model_summary=True,
        strategy='auto',
        log_every_n_steps=50,
        enable_progress_bar=True,
        enable_checkpointing=True,
        deterministic=False,
        benchmark=True,
    )
    
    # Find the best checkpoint to resume from
    def find_best_checkpoint(checkpoint_dir, epoch_offset=0):
        """Find the best checkpoint considering epoch offset"""
        
        # First try last.ckpt
        last_ckpt = os.path.join(checkpoint_dir, 'last.ckpt')
        if os.path.exists(last_ckpt):
            print(f"Resuming from last checkpoint: {last_ckpt}")
            return last_ckpt
        
        # Look for checkpoint with highest epoch number
        import glob
        pattern = os.path.join(checkpoint_dir, '*.ckpt')
        checkpoints = glob.glob(pattern)
        
        if checkpoints:
            # Sort by modification time and take the most recent
            checkpoints.sort(key=os.path.getmtime, reverse=True)
            print(f"Resuming from checkpoint: {checkpoints[0]}")
            return checkpoints[0]
        
        print("No checkpoint found, starting fresh")
        return None
    
    ckpt_path = find_best_checkpoint(checkpoint_dir, actual_completed_epochs)
    
    # Train the model
    print("Starting training...")
    trainer.fit(model, data_module, ckpt_path=ckpt_path)
    
    # Save updated training state
    final_training_state = {
        'completed_epochs': actual_completed_epochs + trainer.current_epoch,
        'best_metric': float(trainer.checkpoint_callback.best_model_score) 
                       if trainer.checkpoint_callback else training_state.get('best_metric'),
        'total_wall_time': training_state.get('total_wall_time', 0) + time.time() - start_time,
        'last_learning_rate': float(trainer.optimizers[0].param_groups[0]['lr']),
        'last_updated': datetime.now().isoformat()
    }
    
    with open(training_state_file, 'w') as f:
        json.dump(final_training_state, f, indent=2)
    
    # Sync checkpoints back to S3
    if args.checkpoint_s3_uri:
        print(f"Syncing checkpoints to {args.checkpoint_s3_uri}")
        import subprocess
        subprocess.run([
            'aws', 's3', 'sync',
            checkpoint_dir,
            args.checkpoint_s3_uri,
            '--exclude', '*',
            '--include', '*.ckpt',
            '--include', 'training_state.json'
        ], check=True)
    
    # Continue with evaluation and analysis (rest of original main() code)
    best_model_path = None
    best_model_score = None
    best_epoch = None
    
    if hasattr(trainer, 'checkpoint_callback') and trainer.checkpoint_callback:
        best_model_path = trainer.checkpoint_callback.best_model_path
        best_model_score = trainer.checkpoint_callback.best_model_score
        if best_model_score is not None:
            best_model_score = float(best_model_score)
    
        if best_model_path and os.path.exists(best_model_path):
            print(f"\nReloading best checkpoint for evaluation: {best_model_path}")
            checkpoint = torch.load(best_model_path, map_location='cpu')
            best_epoch = checkpoint.get('epoch', None)
            
            # Account for continuation metadata if present
            if 'continuation_metadata' in checkpoint:
                total_epochs = checkpoint['continuation_metadata']['total_epochs_trained']
                print(f"Best checkpoint from total epoch {total_epochs}")
            
            model.load_state_dict(checkpoint['state_dict'])
            if best_epoch is not None:
                print(f"Best checkpoint was from epoch {best_epoch} (this run)")
            if best_model_score is not None:
                print(f"Best validation loss (combined): {best_model_score:.6f}")
        else:
            print("Warning: Best checkpoint not found; evaluating with final model state")
    else:
        print("Warning: No checkpoint callback found; evaluating with final model state")
    
    # Validate and test
    print("\nRunning validation on best model...")
    val_results = trainer.validate(model, datamodule=data_module, verbose=False)
    
    print("\nRunning test evaluation on best model...")
    test_results = trainer.test(model, datamodule=data_module, verbose=False)
    
    # Export best checkpoint and save results
    if best_model_path and os.path.exists(best_model_path):
        import shutil
        
        model_files_dir = os.path.join(args.model_dir, 'model_files')
        os.makedirs(model_files_dir, exist_ok=True)
        
        dst_ckpt = os.path.join(model_files_dir, 'model.ckpt')
        print(f"\nCopying best model to: {dst_ckpt}")
        shutil.copy(best_model_path, dst_ckpt)
        
        checkpoint = torch.load(best_model_path, map_location='cpu')
        torch.save(checkpoint['state_dict'], os.path.join(model_files_dir, 'model.pt'))
        
        # Write model_info.json
        model_info = {
            'best_epoch': best_epoch,
            'total_epochs_trained': final_training_state['completed_epochs'],
            'best_score': best_model_score,
            'best_model_path': 'model_files/model.ckpt',
            'monitor_metric': 'val_loss',
            'feature_info': {
                'static_features': data_module.static_cols,
                'dynamic_features': data_module.dynamic_cols_no_target,
                'target': 'streamflow'
            },
            'continuation_info': {
                'total_runs': training_state.get('total_runs', 1),
                'total_wall_time_hours': final_training_state['total_wall_time'] / 3600
            }
        }
        
        model_info_path = os.path.join(model_files_dir, 'model_info.json')
        with open(model_info_path, 'w') as f:
            json.dump(tensor_to_python(model_info), f, indent=2)
        print(f"Saved model info to {model_info_path}")
        
        # Build and write results.json
        results = {
            'experiment_name': args.experiment_name,
            'total_epochs_trained': final_training_state['completed_epochs'],
            'best_model_path': best_model_path,
            'best_model_score': best_model_score,
            'best_model_epoch': best_epoch,
            'best_val_results': tensor_to_python(val_results),
            'test_results': tensor_to_python(test_results),
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
            },
            'continuation_info': final_training_state
        }
        
        if hasattr(trainer, 'logged_metrics'):
            results['final_metrics'] = tensor_to_python(trainer.logged_metrics)
        
        results_path = os.path.join(args.output_data_dir, 'results.json')
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        copy_to_output_files(results_path, args.model_dir, dst_name='results.json')
        print(f"Saved results to {results_path} and copied to output_files/")
        
        # Run post-training analyses
        analysis_root = os.path.join(args.output_data_dir, 'analysis')
        os.makedirs(analysis_root, exist_ok=True)
        
        # 1) Save loss history
        try:
            loss_out = os.path.join(analysis_root, 'loss_tracker')
            save_loss_history(loss_history_cb.history, loss_out)
            copy_to_output_files(loss_out, args.model_dir, dst_name=os.path.join('analysis', 'loss_tracker'))
        except Exception as e:
            print(f"WARNING: loss-tracker failed: {e}")
        
        # 2) Run feature analysis
        try:
            try:
                import featureAnalysis as fa
                print("Imported featureAnalysis via standard import")
            except ModuleNotFoundError:
                fa = _import_feature_analysis_safe(
                    preferred_paths=[args.model_dir, args.output_data_dir]
                )
            
            cwd = os.getcwd()
            os.chdir(analysis_root)
            try:
                _ = fa.run_comprehensive_analysis(
                    model_dir=os.path.abspath(args.model_dir),
                    watershed_df=watershed_df,
                    bucket_name=args.bucket_name,
                    base_data_dir=args.base_data_dir,
                    base_attr_dir=args.base_attr_dir,
                    ckpt_preference='best',
                    use_cache=args.analysis_use_cache,
                    cache_dir=args.analysis_cache_dir,
                    batch_size=args.analysis_batch_size,
                    num_workers=min(4, os.cpu_count() or 4)
                )
            finally:
                os.chdir(cwd)
            
            # Mirror analysis outputs
            for subdir in ['metrics_analysis', 'hydrograph_plots', 'feature_importance_plots']:
                src = os.path.join(analysis_root, subdir)
                if os.path.exists(src) and os.path.isdir(src):
                    copy_to_output_files(src, args.model_dir, dst_name=os.path.join('analysis', subdir))
                    print(f"Mirrored {subdir} to model_dir/output_files/analysis/")
            
            print("featureAnalysis finished and outputs were saved/mirrored.")
            
        except Exception as e:
            import traceback
            print(f"WARNING: featureAnalysis failed: {e}")
            traceback.print_exc()
    
    else:
        print("Warning: No best model checkpoint found to export")
    
    # Final checkpoint sync to S3
    if args.checkpoint_s3_uri:
        print(f"Final sync of all outputs to {args.checkpoint_s3_uri}")
        subprocess.run([
            'aws', 's3', 'sync',
            checkpoint_dir,
            args.checkpoint_s3_uri
        ])
    
    print("\n" + "="*50)
    print("Training completed successfully!")
    print(f"Total epochs trained: {final_training_state['completed_epochs']}/{args.total_epochs}")
    print(f"Total wall time: {final_training_state['total_wall_time']/3600:.2f} hours")
    print(f"Best model saved to: {os.path.join(args.model_dir, 'model_files')}")
    print(f"Results saved to: {os.path.join(args.output_data_dir, 'results.json')}")
    print(f"TensorBoard logs available at: {args.output_data_dir}/{args.experiment_name}")
    print("="*50)


if __name__ == '__main__':
    main()
