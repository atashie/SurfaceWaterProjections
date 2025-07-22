# train.py
from torch.cuda.amp import autocast, GradScaler
import argparse
import json
import os
import sys

# Add all your imports here
import torch
import numpy as np
import pandas as pd


#########################################################################################################################
# Caravan data processing
#########################################################################################################################

import boto3
import re
import pandas as pd
import xarray as xr
import s3fs  # Required by xarray to open S3 files directly
import numpy as np
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



import pandas as pd
import xarray as xr
import s3fs # Required
import scipy # Required for the engine
import numpy as np
import warnings
import boto3
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
        'streamflow', 'potential_evaporation_sum_ERA5_LAND', 'surface_net_solar_radiation_mean',
        'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min', 'total_precipitation_sum'
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


import pandas as pd
import s3fs # Required for pandas to read directly from S3
import warnings

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
    base_path = f"s3://{bucket_name}/{base_attributes_directory.strip('/')}/{subdirectory_name}"
    file1_path = f"{base_path}/attributes_hydroatlas_{subdirectory_name}.csv" # HydroATLAS
    file2_path = f"{base_path}/attributes_other_{subdirectory_name}.csv"      # Other
    file3_path = f"{base_path}/attributes_caravan_{subdirectory_name}.csv"    # Caravan (Optional)
    target_gauge_id = f"{subdirectory_name}_{watershedID}"

    print(f"Searching for attributes for gauge_id: {target_gauge_id}")

    # --- Define Default HydroATLAS Columns and Handle User Input ---
    default_hydroatlas_cols = [
        'pet_mm_syr', 'aet_mm_syr', 
        'pre_mm_syr', 
        'tmp_dc_syr', 'tmp_dc_smn', 'tmp_dc_smx'
        'snw_pc_syr', 'snw_pc_smx', 
        'swc_pc_syr',
        'gdp_ud_sav', 'gdp_ud_ssu', 'gdp_ud_usu',
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

import numpy as np
import pandas as pd
import warnings

def calculate_kge(observed: np.ndarray, simulated: np.ndarray) -> float:
    """
    Calculate Kling-Gupta Efficiency (KGE).

    Args:
        observed (np.ndarray): Array of observed values.
        simulated (np.ndarray): Array of simulated values.

    Returns:
        float: KGE value. Returns np.nan if calculation fails.
    """
    # Ensure inputs are numpy arrays
    observed = np.asarray(observed).flatten()
    simulated = np.asarray(simulated).flatten()

    # Remove NaN values pairwise
    valid_indices = ~np.isnan(observed) & ~np.isnan(simulated)
    if np.sum(valid_indices) < 2: # Need at least 2 pairs for std dev and correlation
        warnings.warn("Not enough valid pairs to calculate KGE after removing NaNs.")
        return np.nan

    obs_valid = observed[valid_indices]
    sim_valid = simulated[valid_indices]

    # Calculate components
    mean_obs = np.mean(obs_valid)
    mean_sim = np.mean(sim_valid)
    std_obs = np.std(obs_valid)
    std_sim = np.std(sim_valid)

    # Avoid division by zero for standard deviations
    if std_obs == 0 or std_sim == 0 or mean_obs == 0:
         warnings.warn("Zero standard deviation or mean observed value encountered in KGE calculation.")
         # Return NaN or a default low value depending on desired behavior
         return np.nan # Or perhaps -np.inf or a very negative number

    # Pearson correlation coefficient
    correlation_matrix = np.corrcoef(obs_valid, sim_valid)
    r = correlation_matrix[0, 1]

    # Alpha (ratio of standard deviations)
    alpha = std_sim / std_obs

    # Beta (ratio of means)
    beta = mean_sim / mean_obs

    # Calculate KGE
    kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)

    return kge


import traceback

def unnormalize_predictions(predictions, norm_params):
    """
    Reverse normalization for streamflow predictions.
    
    Args:
        predictions: Normalized log-transformed predictions
        norm_params: Dictionary with target_mean and target_std
        
    Returns:
        Original scale predictions
    """
    # Reverse standardization
    log_predictions = predictions * norm_params['target_std'] + norm_params['target_mean']
    
    # Reverse log transform
    original_predictions = np.expm1(log_predictions)  # exp(x) - 1
    
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


import psutil

def get_ram_usage():
    """
    Returns current RAM usage statistics.
    """
    virtual_mem = psutil.virtual_memory()
    total_ram_gb = virtual_mem.total / (1024**3)  # Total RAM in GB
    used_ram_gb = virtual_mem.used / (1024**3)    # Used RAM in GB
    percent_used = virtual_mem.percent             # Percentage of RAM used
    
    return {
        "total_ram_gb": total_ram_gb,
        "used_ram_gb": used_ram_gb,
        "percent_used": percent_used
    }

def print_ram_usage(context_message="Current RAM Usage"):
    """
    Prints the current RAM usage statistics in a formatted way.
    """
    usage = get_ram_usage()
    print(f"--- {context_message} ---")
    print(f"  Total RAM: {usage['total_ram_gb']:.2f} GB")
    print(f"  Used RAM:  {usage['used_ram_gb']:.2f} GB")
    print(f"  Percent Used: {usage['percent_used']:.1f}%")
    print("--------------------------------------")


class EarlyStopping:
    """Early stopping to stop training when validation loss doesn't improve"""
    def __init__(self, patience=7, verbose=True, delta=0, save_path='checkpoint.pt'):
        """
        Args:
            patience (int): How many epochs to wait after last improvement
            verbose (bool): If True, prints messages
            delta (float): Minimum change to qualify as improvement
            save_path (str): Path to save the best model
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.save_path = save_path
        
    def __call__(self, val_loss, model, optimizer, epoch, additional_info=None):
        score = -val_loss  # We want to minimize loss, so negate it
        
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model, optimizer, epoch, additional_info)
        elif score < self.best_score + self.delta:
            self.counter += 1
            if self.verbose:
                print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model, optimizer, epoch, additional_info)
            self.counter = 0
            
    def save_checkpoint(self, val_loss, model, optimizer, epoch, additional_info):
        """Saves model when validation loss decreases"""
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}). Saving model...')
        
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'val_loss': val_loss,
        }
        if additional_info:
            checkpoint.update(additional_info)
            
        torch.save(checkpoint, self.save_path)
        self.val_loss_min = val_loss


def get_learning_rate_scheduler(optimizer, scheduler_type='reduce_on_plateau', **kwargs):
    """
    Get a learning rate scheduler based on type.
    
    Args:
        optimizer: PyTorch optimizer
        scheduler_type: Type of scheduler ('reduce_on_plateau', 'cosine', 'step', 'exponential')
        **kwargs: Additional arguments for the scheduler
    
    Returns:
        scheduler object
    """
    if scheduler_type == 'reduce_on_plateau':
        return torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode='min',
            factor=kwargs.get('factor', 0.5),
            patience=kwargs.get('patience', 5),
            verbose=True,
            min_lr=kwargs.get('min_lr', 1e-6)
        )
    elif scheduler_type == 'cosine':
        return torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer,
            T_max=kwargs.get('T_max', 100),
            eta_min=kwargs.get('min_lr', 1e-6)
        )
    elif scheduler_type == 'step':
        return torch.optim.lr_scheduler.StepLR(
            optimizer,
            step_size=kwargs.get('step_size', 30),
            gamma=kwargs.get('gamma', 0.1)
        )
    elif scheduler_type == 'exponential':
        return torch.optim.lr_scheduler.ExponentialLR(
            optimizer,
            gamma=kwargs.get('gamma', 0.95)
        )
    elif scheduler_type == 'warmup_cosine':
        # Custom warmup + cosine annealing
        return WarmupCosineScheduler(
            optimizer,
            warmup_epochs=kwargs.get('warmup_epochs', 5),
            max_epochs=kwargs.get('max_epochs', 100),
            min_lr=kwargs.get('min_lr', 1e-6)
        )
    else:
        return None

class WarmupCosineScheduler:
    """Learning rate scheduler with linear warmup and cosine annealing"""
    def __init__(self, optimizer, warmup_epochs, max_epochs, min_lr=1e-6):
        self.optimizer = optimizer
        self.warmup_epochs = warmup_epochs
        self.max_epochs = max_epochs
        self.min_lr = min_lr
        self.base_lr = optimizer.param_groups[0]['lr']
        self.current_epoch = 0
        
    def step(self, epoch=None):
        if epoch is not None:
            self.current_epoch = epoch
        else:
            self.current_epoch += 1
            
        if self.current_epoch < self.warmup_epochs:
            # Linear warmup
            lr = self.base_lr * (self.current_epoch + 1) / self.warmup_epochs
        else:
            # Cosine annealing
            progress = (self.current_epoch - self.warmup_epochs) / (self.max_epochs - self.warmup_epochs)
            lr = self.min_lr + (self.base_lr - self.min_lr) * 0.5 * (1 + np.cos(np.pi * progress))
            
        for param_group in self.optimizer.param_groups:
            param_group['lr'] = lr
            
    def get_last_lr(self):
        return [param_group['lr'] for param_group in self.optimizer.param_groups]


# Assume load_and_process_watershed_data and get_watershed_attributes
# functions from previous steps are defined and imported here.
# from your_module import load_and_process_watershed_data, get_watershed_attributes

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
        'potential_evaporation_sum_ERA5_LAND', 'surface_net_solar_radiation_mean',
        'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min',
        'total_precipitation_sum'
    ]
    time_varying_cols = weather_cols  # Define early
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
                'potential_evaporation_sum_ERA5_LAND', 'surface_net_solar_radiation_mean',
                'temperature_2m_max', 'temperature_2m_mean', 'temperature_2m_min',
                'total_precipitation_sum'
            ]
            target_col = 'streamflow'
            
            # streamflow validataion (check for unrealistic values, or values of 0 that will be problematic for log transformation)
            invalid_mask = timeseries_df[target_col] <= 0
            if invalid_mask.any():
                n_invalid = invalid_mask.sum()
                print(f"WARNING: Found {n_invalid} non-positive streamflow values in {group_id}")
                timeseries_df.loc[invalid_mask, target_col] = 1e-6
            
            
            date_col = 'date'
            keep_cols = [date_col, target_col] + weather_cols
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
        norm_params = {
            'static_means': {col: combined_df[col].mean() for col in static_attribute_names},
            'static_stds': {col: combined_df[col].std() for col in static_attribute_names},
            'dynamic_means': {col: combined_df[col].mean() for col in time_varying_cols},
            'dynamic_stds': {col: combined_df[col].std() for col in time_varying_cols},
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
        
        # Dynamic features
        for col in time_varying_cols:
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
            
    # Ensure correct types
    combined_df[target_col] = combined_df[target_col].astype(float)
    for col in time_varying_cols:
        combined_df[col] = combined_df[col].astype(float)
    for col in static_attribute_names:
        combined_df[col] = combined_df[col].astype(float)

    print(f"Successfully prepared combined DataFrame with shape: {combined_df.shape}")
    print(f"Identified static reals: {static_attribute_names}")
    print(f"Identified time-varying reals: {time_varying_cols}")

    return combined_df, static_attribute_names, time_varying_cols, norm_params




import torch
#torch.set_num_threads(16)  # Limit PyTorch threads
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import warnings
from datetime import datetime
import json
import gc
import os
#os.environ["OMP_NUM_THREADS"] = "4"  # Limit OpenMP threads
#os.environ["MKL_NUM_THREADS"] = "4"   # Limit MKL threads
#os.environ["NUMEXPR_NUM_THREADS"] = "4"  # Limit NumExpr threads
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm  
import pickle
import signal  # signal handling for spot instances



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
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
            nn.Sigmoid()  # Output between 0 and 1 for gate modulation
        )
        
        # Optional: Additional static feature integration
        self.static_hidden_init = nn.Sequential(
            nn.Linear(static_input_size, hidden_size * num_layers),
            nn.ReLU()
        )
        
        # Time Modulation Network (NEW)
        if self.use_time_modulation:
            self.time_modulator = nn.Sequential(
                nn.Conv1d(hidden_size, hidden_size, kernel_size=1),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Conv1d(hidden_size, hidden_size, kernel_size=1),
                nn.Sigmoid()  # Output between 0 and 1 for temporal gating
            )
        
        # Output layers
        self.output_layer = nn.Sequential(
            nn.Linear(hidden_size + static_input_size, hidden_size),
            nn.ReLU(),
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
                    nn.init.kaiming_uniform_(param, a=0, mode='fan_in', nonlinearity='relu')
                    param.data *= 0.5  # Scale down for stability
                else:
                    # Use He initialization for other layers
                    nn.init.kaiming_uniform_(param, a=0, mode='fan_in', nonlinearity='relu')
                    param.data *= 0.1  # Scale down
            elif 'bias' in name:
                nn.init.constant_(param, 0.0)




class EntityAwareLSTM_old(nn.Module):
    """Entity-Aware LSTM following Kratzert et al. (2019) approach"""
    def __init__(self, dynamic_input_size, static_input_size, hidden_size=256, 
                 num_layers=1, dropout=0.2):
        super().__init__()
        
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        
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
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
            nn.Sigmoid()  # Output between 0 and 1 for gate modulation
        )
        
        # Optional: Additional static feature integration
        self.static_hidden_init = nn.Sequential(
            nn.Linear(static_input_size, hidden_size * num_layers),
            nn.ReLU()
        )
        
        # Output layers
        self.output_layer = nn.Sequential(
            nn.Linear(hidden_size + static_input_size, hidden_size),
            nn.ReLU(),
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
        
        # Apply entity-aware modulation to LSTM outputs
        # This scales the LSTM output based on watershed characteristics
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
                    nn.init.xavier_uniform_(param, gain=0.5)  # Small gain is 0.1
                else:
                    # Use He initialization for other layers
                    nn.init.kaiming_uniform_(param, a=0, mode='fan_in', nonlinearity='relu')
                    param.data *= 0.1  # Scale down
            elif 'bias' in name:
                nn.init.constant_(param, 0.0)



def print_gpu_memory(context_message="Current GPU Usage"):
    """
    Prints current GPU memory usage statistics.
    """
    if torch.cuda.is_available():
        # Get current GPU memory stats
        allocated = torch.cuda.memory_allocated() / (1024**3)  # GB
        reserved = torch.cuda.memory_reserved() / (1024**3)    # GB
        max_allocated = torch.cuda.max_memory_allocated() / (1024**3)  # GB
        
        print(f"--- {context_message} ---")
        print(f"  GPU Memory Allocated: {allocated:.2f} GB")
        print(f"  GPU Memory Reserved:  {reserved:.2f} GB")
        print(f"  Max GPU Memory Used:  {max_allocated:.2f} GB")
        print("--------------------------------------")
        
        # Optional: Print per-GPU stats if multiple GPUs
        if torch.cuda.device_count() > 1:
            for i in range(torch.cuda.device_count()):
                print(f"  GPU {i}: {torch.cuda.memory_allocated(i)/(1024**3):.2f} GB")
    else:
        print(f"--- {context_message} ---")
        print("  No GPU available")
        print("--------------------------------------")

def print_memory_usage(context_message="Current Memory Usage"):
    """
    Prints current CPU memory usage and process information.
    """
    import psutil
    import os
    
    # System-wide memory
    virtual_mem = psutil.virtual_memory()
    
    # Current process memory
    process = psutil.Process(os.getpid())
    process_info = process.memory_info()
    
    print(f"--- {context_message} ---")
    print(f"System Memory:")
    print(f"  Total: {virtual_mem.total / (1024**3):.2f} GB")
    print(f"  Used:  {virtual_mem.used / (1024**3):.2f} GB ({virtual_mem.percent:.1f}%)")
    print(f"  Available: {virtual_mem.available / (1024**3):.2f} GB")
    print(f"Process Memory:")
    print(f"  RSS (Resident): {process_info.rss / (1024**3):.2f} GB")
    print(f"  VMS (Virtual):  {process_info.vms / (1024**3):.2f} GB")
    print(f"  CPU Count: {psutil.cpu_count(logical=False)} physical, {psutil.cpu_count()} logical")
    print(f"  Process Threads: {process.num_threads()}")
    print("--------------------------------------")


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



def train_ea_lstm_model_with_scheduling(
    watersheds_df,
    bucket_name,
    base_data_dir,
    base_attr_dir,
    hyperparameters,
    save_dir,
    chunk_size=50,
    epochs_per_chunk=10,
    device='cpu',
    random_seed=42,
    gradient_accumulation_steps=8,
    memory_cleanup_frequency=5,
    # New parameters for early stopping and LR scheduling
    early_stopping_patience=15,
    lr_scheduler_type='reduce_on_plateau',
    lr_scheduler_params=None,
    target_noise_std=0.1,  # New parameter
    target_noise_type='multiplicative'  # New parameter
):
    """
    Train EntityAwareLSTM model with early stopping and learning rate scheduling.
    
    Additional Args:
        early_stopping_patience: Number of epochs without improvement before stopping
        lr_scheduler_type: Type of learning rate scheduler to use
        lr_scheduler_params: Dict of parameters for the learning rate scheduler
    """

    # anomaly detection for debugging
    torch.autograd.set_detect_anomaly(True)
    
    # CPU-specific setup
    n_threads = min(4, psutil.cpu_count(logical=False) // 2)  # Half the physical cores
    torch.set_num_threads(n_threads)
    os.environ["OMP_NUM_THREADS"] = str(n_threads)
    os.environ["MKL_NUM_THREADS"] = str(n_threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(n_threads)
    torch.set_flush_denormal(True)

    scaler = GradScaler(
        init_scale=1024,  # Lower initial scale
        growth_factor=2,
        backoff_factor=0.5,
        growth_interval=2000,
        enabled=device == 'cuda'
    ) if device == 'cuda' else None
    
    # Set random seeds
    np.random.seed(random_seed)
    torch.manual_seed(random_seed)
    
    # Create save directory
    os.makedirs(save_dir, exist_ok=True)
    
    # Log setup
    log_file = os.path.join(save_dir, 'training_log.txt')
    log_buffer = []
    
    def log_message(msg, force_flush=False):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        full_msg = f"[{timestamp}] {msg}"
        print(full_msg)
        log_buffer.append(full_msg)
        
        if force_flush or len(log_buffer) >= 10:
            with open(log_file, 'a') as f:
                f.write('\n'.join(log_buffer) + '\n')
            log_buffer.clear()
    
    log_message("Starting EntityAwareLSTM training with early stopping and LR scheduling")
    log_message(f"Early stopping patience: {early_stopping_patience}")
    log_message(f"LR scheduler: {lr_scheduler_type}")
    log_message(f"Total watersheds available: {len(watersheds_df)}")

    # Initialize target noise wrapper
    target_noise = TargetNoiseWrapper(
        noise_std=target_noise_std,
        noise_type=target_noise_type
    )
    
    log_message(f"Using target noise: std={target_noise_std}, type={target_noise_type}")
    
    # Initial memory state
    print_memory_usage("Initial state")

    # Resume from checkpoint if exists
    start_epoch = 0
    resume_checkpoint = None
    checkpoint_dir = '/opt/ml/checkpoints'
    
    if os.path.exists(checkpoint_dir):
        # Look for checkpoint files
        checkpoint_files = sorted([f for f in os.listdir(checkpoint_dir) if f.endswith('.pt')])
        
        if checkpoint_files:
            # Use the most recent checkpoint
            latest_checkpoint = os.path.join(checkpoint_dir, checkpoint_files[-1])
            log_message(f"Found checkpoint: {latest_checkpoint}")
            
            try:
                resume_checkpoint = torch.load(latest_checkpoint, map_location=device)
                if 'norm_params' in resume_checkpoint:
                    train_norm_params = resume_checkpoint['norm_params']
                    log_message("Loaded normalization parameters from checkpoint")
                else:
                    # If not in checkpoint, need to recompute or error
                    log_message("WARNING: No normalization parameters in checkpoint!")

                # Check for a poisoned checkpoint from a previous NaN state
                if 'early_stopping_state' in resume_checkpoint and np.isnan(resume_checkpoint['early_stopping_state'].get('best_score', 0.0)):
                    log_message("!!! WARNING: POISONED CHECKPOINT DETECTED (NaN state) !!!")
                    log_message("Starting with fresh optimizer, scheduler, and history, but keeping model weights.")
                    # Invalidate the parts of the checkpoint that cause cascading failure
                    resume_checkpoint['optimizer_state_dict'] = None
                    resume_checkpoint['scheduler_state_dict'] = None
                    resume_checkpoint['early_stopping_state'] = None
                    training_history = {} # Force re-initialization
                
                # Extract state from checkpoint
                start_epoch = resume_checkpoint['epoch']
                training_history = resume_checkpoint['training_history']
                chunks_processed = resume_checkpoint.get('chunks_processed', 0)
                
                # Feature info
                if 'feature_info' in resume_checkpoint:
                    static_cols = resume_checkpoint['feature_info']['static_cols']
                    dynamic_cols_no_target = resume_checkpoint['feature_info']['dynamic_cols_no_target']
                    dynamic_input_size = resume_checkpoint['feature_info']['dynamic_input_size']
                    static_input_size = resume_checkpoint['feature_info']['static_input_size']
                    
                log_message(f"Resuming from epoch {start_epoch}")
                log_message(f"Training history has {len(training_history['epoch'])} entries")
                
            except Exception as e:
                log_message(f"Error loading checkpoint: {e}. Starting fresh.")
                resume_checkpoint = None
                start_epoch = 0
                training_history = {
                    'train_loss': [],
                    'val_loss': [],
                    'train_kge': [],
                    'val_kge': [],
                    'epoch': [],
                    'learning_rate': []
                }
        else:
            log_message("No checkpoint files found in checkpoint directory")
            training_history = {
                'train_loss': [],
                'val_loss': [],
                'train_kge': [],
                'val_kge': [],
                'epoch': [],
                'learning_rate': []
            }
    else:
        log_message("Checkpoint directory does not exist, starting fresh")
        training_history = {
            'train_loss': [],
            'val_loss': [],
            'train_kge': [],
            'val_kge': [],
            'epoch': [],
            'learning_rate': []
        }



    
    # Split watersheds
    train_watersheds, temp_watersheds = train_test_split(
        watersheds_df, test_size=0.4, random_state=random_seed
    )
    val_watersheds, test_watersheds = train_test_split(
        temp_watersheds, test_size=0.5, random_state=random_seed
    )
    
    log_message(f"Train watersheds: {len(train_watersheds)}")
    log_message(f"Val watersheds: {len(val_watersheds)}")
    log_message(f"Test watersheds: {len(test_watersheds)}")
    
    # Save splits
    train_watersheds.to_csv(os.path.join(save_dir, 'train_watersheds.csv'), index=False)
    val_watersheds.to_csv(os.path.join(save_dir, 'val_watersheds.csv'), index=False)
    test_watersheds.to_csv(os.path.join(save_dir, 'test_watersheds.csv'), index=False)
    
    # Compute normalization parameters from full training data
    log_message("Computing normalization parameters from full training data...")
    _, _, _, train_norm_params = prepare_ptf_dataframe(
        train_watersheds, bucket_name, base_data_dir, base_attr_dir,
        compute_norm_params=True  # Only compute params, don't normalize
    )

    # Save normalization parameters
    norm_params_path = os.path.join(save_dir, 'normalization_params.json')
    with open(norm_params_path, 'w') as f:
        json_params = {
            'static_means': {k: float(v) for k, v in train_norm_params['static_means'].items()},
            'static_stds': {k: float(v) for k, v in train_norm_params['static_stds'].items()},
            'dynamic_means': {k: float(v) for k, v in train_norm_params['dynamic_means'].items()},
            'dynamic_stds': {k: float(v) for k, v in train_norm_params['dynamic_stds'].items()},
            'target_mean': float(train_norm_params['target_mean']),
            'target_std': float(train_norm_params['target_std'])
        }
        json.dump(json_params, f, indent=2)
    log_message(f"Saved normalization parameters to {norm_params_path}")
    
    # Extract hyperparameters
    hp = hyperparameters
    sequence_length = hp.get('sequence_length', 365)
    hidden_size = hp.get('hidden_size', 256)
    learning_rate = hp.get('learning_rate', 0.0001)
    batch_size = hp.get('batch_size', 256)
    dropout = hp.get('dropout', 0.2)
    num_layers = hp.get('num_layers', 1)
    total_epochs = hp.get('total_epochs', 100)
    
    # Setup default LR scheduler parameters if not provided
    if lr_scheduler_params is None:
        lr_scheduler_params = {
            'factor': 0.5,
            'patience': 5,
            'min_lr': 1e-6,
            'T_max': total_epochs,
            'warmup_epochs': 5,
            'max_epochs': total_epochs
        }
    
    # Save all hyperparameters
    hp_to_save = hyperparameters.copy()
    hp_to_save.update({
        'gradient_accumulation_steps': gradient_accumulation_steps,
        'effective_batch_size': batch_size * gradient_accumulation_steps,
        'early_stopping_patience': early_stopping_patience,
        'lr_scheduler_type': lr_scheduler_type,
        'lr_scheduler_params': lr_scheduler_params
    })
    
    with open(os.path.join(save_dir, 'hyperparameters.json'), 'w') as f:
        json.dump(hp_to_save, f, indent=2)
    
    # Prepare data
    log_message("Loading and processing training data...")
    train_chunks = [train_watersheds[i:i+chunk_size] 
                   for i in range(0, len(train_watersheds), chunk_size)]
    
    # Determine feature dimensions
    sample_chunk = train_chunks[0][:5]
    sample_df, static_cols, dynamic_cols, norm_params_sample = prepare_ptf_dataframe(
        sample_chunk, bucket_name, base_data_dir, base_attr_dir, norm_params=train_norm_params
    )
    
    if sample_df.empty:
        raise ValueError("Failed to load sample data")
    
    del sample_df
    gc.collect()
    
    dynamic_cols_no_target = [col for col in dynamic_cols if col != 'streamflow']
    dynamic_input_size = len(dynamic_cols_no_target)
    static_input_size = len(static_cols)
    
    log_message(f"Dynamic input size: {dynamic_input_size}")
    log_message(f"Static input size: {static_input_size}")
    
    gc.collect()

    # Save normalization parameters to a JSON file
    norm_params_path = os.path.join(save_dir, 'normalization_params.json')
    with open(norm_params_path, 'w') as f:
        # Convert numpy types to standard Python types for JSON serialization
        json_params = {
            'static_means': {k: float(v) for k, v in train_norm_params['static_means'].items()},
            'static_stds': {k: float(v) for k, v in train_norm_params['static_stds'].items()},
            'dynamic_means': {k: float(v) for k, v in train_norm_params['dynamic_means'].items()},
            'dynamic_stds': {k: float(v) for k, v in train_norm_params['dynamic_stds'].items()},
            'target_mean': float(train_norm_params['target_mean']),
            'target_std': float(train_norm_params['target_std'])
        }
        json.dump(json_params, f, indent=2)
    log_message(f"Saved normalization parameters to {norm_params_path}")

    
    
    
    # Initialize model
    model = EntityAwareLSTM(
        dynamic_input_size=dynamic_input_size,
        static_input_size=static_input_size,
        hidden_size=hidden_size,
        num_layers=num_layers,
        dropout=dropout,
        use_time_modulation=True 
    ).to(device)
    
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    
    # Initialize learning rate scheduler
    scheduler = get_learning_rate_scheduler(optimizer, lr_scheduler_type, **lr_scheduler_params)
    
    # Initialize early stopping
    early_stopping = EarlyStopping(
        patience=early_stopping_patience,
        verbose=True,
        save_path=os.path.join(save_dir, 'best_model.pt')
    )


    # Restore states from checkpoint if resuming
    if resume_checkpoint is not None:
        log_message("Restoring model and optimizer states from checkpoint...")
        
        # Restore model state
        model.load_state_dict(resume_checkpoint['model_state_dict'])
        
        # Restore optimizer state
        optimizer.load_state_dict(resume_checkpoint['optimizer_state_dict'])
        
        # Restore scheduler state if it exists
        if scheduler and 'scheduler_state_dict' in resume_checkpoint and resume_checkpoint['scheduler_state_dict'] is not None:
            try:
                scheduler.load_state_dict(resume_checkpoint['scheduler_state_dict'])
                log_message("Restored learning rate scheduler state")
            except Exception as e:
                log_message(f"Warning: Could not restore scheduler state: {e}")
        
        # Restore early stopping state
        if 'early_stopping_state' in resume_checkpoint:
            es_state = resume_checkpoint['early_stopping_state']
            early_stopping.counter = es_state['counter']
            early_stopping.best_score = es_state['best_score']
            early_stopping.val_loss_min = es_state['val_loss_min']
            early_stopping.early_stop = es_state.get('early_stop', False)
            log_message(f"Restored early stopping state: counter={es_state['counter']}, best_score={es_state['best_score']}")
        
        log_message("Successfully restored all states from checkpoint")

    
    
    # Training history
    if resume_checkpoint is None:
        training_history = {
            'train_loss': [],
            'val_loss': [],
            'train_kge': [],
            'val_kge': [],
            'epoch': [],
            'learning_rate': []
        }
    else:
        # training_history was already loaded from checkpoint
        log_message(f"Continuing with existing training history ({len(training_history['epoch'])} epochs)")
                    
    # Training loop
        # Initialize tracking variables
    gradient_explosion_count = 0
    nan_batch_count = 0
    global_step = 0
    if resume_checkpoint is None:
        chunks_processed = 0
    
    for epoch in range(0, total_epochs, epochs_per_chunk):
        # Check early stopping
        if early_stopping.early_stop:
            log_message("Early stopping triggered! Stopping training.", force_flush=True)
            break

 
        # Manually reduce LR based on epoch
        current_epoch = epoch + epochs_per_chunk - 1
#        if current_epoch > 100 and optimizer.param_groups[0]['lr'] > 1e-4:
#            for param_group in optimizer.param_groups:
#                param_group['lr'] = 1e-4
#            log_message(f"Manually reduced learning rate to 1e-4 at epoch {current_epoch}")
#        elif current_epoch > 200 and optimizer.param_groups[0]['lr'] > 1e-5:
#            for param_group in optimizer.param_groups:
#                param_group['lr'] = 1e-5
#            log_message(f"Manually reduced learning rate to 1e-5 at epoch {current_epoch}")
        
        log_message(f"\n=== Global epochs {epoch} to {min(epoch + epochs_per_chunk - 1, total_epochs - 1)} ===")
        
        # Log current learning rate
        current_lr = optimizer.param_groups[0]['lr']
        log_message(f"Current learning rate: {current_lr:.6f}")
        
        gc.collect()
        np.random.shuffle(train_chunks)
        
        epoch_train_losses = []
        
        for chunk_idx, chunk_watersheds in enumerate(train_chunks):
            chunks_processed += 1
            log_message(f"\nProcessing chunk {chunk_idx + 1}/{len(train_chunks)}")
            
            try:
                # Check for termination request
                termination_file = '/opt/ml/checkpoints/TERMINATION_REQUESTED'
                if os.path.exists(termination_file):
                    log_message("!!! Spot instance termination detected! Saving emergency checkpoint...", force_flush=True)
                    
                    # Save complete state
                    emergency_checkpoint = os.path.join(checkpoint_dir, 'emergency_checkpoint.pt')
                    torch.save({
                        'epoch': current_epoch,
                        'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': optimizer.state_dict(),
                        'scheduler_state_dict': scheduler.state_dict() if scheduler else None,
                        'training_history': training_history,
                        'norm_params': train_norm_params,
                        'chunks_processed': chunks_processed,
                        'chunk_idx': chunk_idx,
                        'early_stopping_state': {
                            'counter': early_stopping.counter,
                            'best_score': early_stopping.best_score,
                            'val_loss_min': early_stopping.val_loss_min,
                            'early_stop': early_stopping.early_stop
                        },
                        'feature_info': {
                            'static_cols': static_cols,
                            'dynamic_cols_no_target': dynamic_cols_no_target,
                            'dynamic_input_size': dynamic_input_size,
                            'static_input_size': static_input_size
                        }
                    }, emergency_checkpoint)
                    
                    log_message(f"Emergency checkpoint saved to {emergency_checkpoint}", force_flush=True)
                    log_message("Exiting due to spot instance termination...", force_flush=True)
                    
                    # Clean exit
                    sys.exit(0)

                
                # Load chunk data
                chunk_df, _, _, _ = prepare_ptf_dataframe(
                    chunk_watersheds, bucket_name, base_data_dir, base_attr_dir, norm_params=train_norm_params  
                )
                
                if chunk_df.empty:
                    log_message(f"Warning: Empty data for chunk {chunk_idx + 1}, skipping")
                    continue
                
                # Create dataset
                train_dataset = WatershedDataset(
                    chunk_df, static_cols, dynamic_cols_no_target, 
                    'streamflow', sequence_length
                )
                
                del chunk_df
                gc.collect()
                
                if len(train_dataset) == 0:
                    log_message(f"Warning: No valid sequences in chunk {chunk_idx + 1}, skipping")
                    del train_dataset
                    continue
                
                train_loader = DataLoader(
                    train_dataset, 
                    batch_size=batch_size, 
                    shuffle=True, 
                    num_workers=4,  # Reduce from 8 to 4 for GPU
                    pin_memory=True,  # Important for GPU!
                    persistent_workers=True,
                    prefetch_factor=2,  # Added for speed
                    drop_last=True  # Add this to avoid variable batch sizes
                )
                
                # Train for epochs_per_chunk
                for chunk_epoch in range(epochs_per_chunk):
                    current_epoch = epoch + chunk_epoch
                    if current_epoch >= total_epochs:
                        break
                    
                    if early_stopping.early_stop:
                        break
                        
                    model.train()
                    chunk_losses = []
                    
                    show_progress = (current_epoch == 0) or (current_epoch % 10 == 0)
                    
                    if show_progress:
                        progress_bar = tqdm(train_loader, 
                                          desc=f'Epoch {current_epoch}, Chunk {chunk_idx+1}')
                    else:
                        progress_bar = train_loader
                    
                    optimizer.zero_grad()
                    
                    for batch_idx, (dynamic_seq, static_feat, target) in enumerate(progress_bar):
                            # Move to GPU
                        dynamic_seq = dynamic_seq.to(device)
                        static_feat = static_feat.to(device)
                        target = target.to(device)
                        
                        if device == 'cuda':
                            with autocast():
                                predictions = model(dynamic_seq, static_feat)
                                noisy_target = target_noise.add_noise(target, training=model.training)
                                loss = criterion(predictions, noisy_target)

                                # Add this NaN propagation check here:
                                if torch.isnan(loss) or torch.isinf(loss):
                                    nan_batch_count += 1
                                    log_message(f"Warning: Invalid loss detected: {loss.item()} at epoch {current_epoch}, batch {batch_idx}")
                                    optimizer.zero_grad()  # Clear any partial gradients
                                    continue  # Skip this batch
                            
                            scaler.scale(loss / gradient_accumulation_steps).backward()
                            
                            if (batch_idx + 1) % gradient_accumulation_steps == 0:
                                scaler.unscale_(optimizer)

                                total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)#1.0)
                                if total_norm > 10:  # Warning threshold
                                    log_message(f"Warning: Large gradient norm: {total_norm:.2f}")
                                if torch.isnan(total_norm):
                                    log_message("Warning: NaN gradients detected! Skipping batch.")
                                    optimizer.zero_grad()
                                    scaler.update()
                                    continue

                                # Emergency LR reduction on extreme gradients
                                if total_norm > 50:
                                    gradient_explosion_count += 1 
                                    current_lr = optimizer.param_groups[0]['lr']
                                    new_lr = max(current_lr * 0.5, 1e-7)
                                    for param_group in optimizer.param_groups:
                                        param_group['lr'] = new_lr
                                    log_message(f"Emergency LR reduction due to gradient explosion: {current_lr} -> {new_lr}")

#                                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                                scaler.step(optimizer)
                                scaler.update()
                                optimizer.zero_grad()
                        else:
                            # CPU path (no mixed precision)
                            predictions = model(dynamic_seq, static_feat)
                            noisy_target = target_noise.add_noise(target, training=model.training)
                            loss = criterion(predictions, noisy_target)
                            loss = loss / gradient_accumulation_steps
                            loss.backward()
                            
                            if (batch_idx + 1) % gradient_accumulation_steps == 0:
                                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
                                optimizer.step()
                                optimizer.zero_grad()          
                                
                        loss_value = loss.item() * gradient_accumulation_steps
                        chunk_losses.append(loss_value)

                        
                        if batch_idx % 50 == 0:
                            del predictions, loss
                        
                        if show_progress and hasattr(progress_bar, 'set_postfix'):
                            progress_bar.set_postfix({'loss': np.mean(chunk_losses[-100:]) if len(chunk_losses) > 100 else np.mean(chunk_losses)})
                    
                        if batch_idx % 10 == 0 and device == 'cuda':
                            torch.cuda.empty_cache()

                    
                    # Handle remaining gradients
                    if (batch_idx + 1) % gradient_accumulation_steps != 0:
                        total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
                        if total_norm > 100:  # Emergency brake
                            optimizer.zero_grad()  # Skip this update entirely
                    
                    mean_chunk_loss = np.mean(chunk_losses)
                    epoch_train_losses.append(mean_chunk_loss)
                    
                    log_message(f"Chunk {chunk_idx+1} - Epoch {current_epoch}: Loss={mean_chunk_loss:.4f}")
                
                # Cleanup
                del train_dataset, train_loader
                if 'dynamic_seq' in locals():
                    del dynamic_seq, static_feat, target
                gc.collect()
                
                # Periodic deep cleanup
                if chunks_processed % memory_cleanup_frequency == 0 or chunks_processed == len(train_chunks):
                    log_message("Performing deep memory cleanup...")
                    gc.collect()
                    import ctypes
                    if hasattr(ctypes, 'CDLL'):
                        log_message("confirming hasattr(ctypes, 'CDLL')...")
                        try:
                            libc = ctypes.CDLL("libc.so.6")
                            libc.malloc_trim(0)
                            log_message("confirming libc.malloc_trim(0) was performed...")
                        except:
                            log_message("confirming libc.malloc_trim(0) has failed...")
                            pass
                
            except Exception as e:
                log_message(f"Error processing chunk {chunk_idx + 1}: {str(e)}", force_flush=True)
                import traceback
                traceback.print_exc()
                gc.collect()
                continue
        
        # Validation and scheduler step
        if (epoch + epochs_per_chunk) % 5 == 0 or early_stopping.early_stop:
            log_message("\nRunning validation...", force_flush=True)
            gc.collect()
            
            val_loss, val_kge = evaluate_model_gpu_optimized(
                model, val_watersheds, bucket_name, base_data_dir, base_attr_dir,
                static_cols, dynamic_cols_no_target, sequence_length, 
                batch_size, device, chunk_size, norm_params=train_norm_params
            )
            
            log_message(f"Validation - Loss: {val_loss:.4f}, KGE: {val_kge:.4f}")
            
            # Update history
            current_epoch_num = min(epoch + epochs_per_chunk, total_epochs)
            training_history['epoch'].append(current_epoch_num)
            training_history['val_loss'].append(val_loss)
            training_history['val_kge'].append(val_kge)
            training_history['train_loss'].append(np.mean(epoch_train_losses) if epoch_train_losses else np.nan)
            training_history['learning_rate'].append(optimizer.param_groups[0]['lr'])
            
            # Early stopping check
            early_stopping(val_loss, model, optimizer, current_epoch_num, 
                         additional_info={'val_kge': val_kge})
            
            # Learning rate scheduler step
            if scheduler:
                if lr_scheduler_type == 'reduce_on_plateau':
                    scheduler.step(val_loss)
                    print(f"Scheduler step called with val_loss: {val_loss}")  # ADD THIS LINE HERE
                elif hasattr(scheduler, 'step'):
                    scheduler.step()
                    
                # Log new learning rate if changed
                new_lr = optimizer.param_groups[0]['lr']
                if new_lr != current_lr:
                    log_message(f"Learning rate changed: {current_lr:.6f} -> {new_lr:.6f}")
            
            # Save intermediate results
            with open(os.path.join(save_dir, 'training_history.pkl'), 'wb') as f:
                pickle.dump(training_history, f)

            checkpoint_dir = '/opt/ml/checkpoints'
            if os.path.exists(checkpoint_dir):
                checkpoint_path = os.path.join(checkpoint_dir, f'checkpoint_epoch_{current_epoch_num}.pt')
                torch.save({
                    'epoch': current_epoch_num,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'scheduler_state_dict': scheduler.state_dict() if scheduler else None,
                    'training_history': training_history,
                    'norm_params': train_norm_params,
                    'early_stopping_state': {
                        'counter': early_stopping.counter,
                        'best_score': early_stopping.best_score,
                        'val_loss_min': early_stopping.val_loss_min,
                        'early_stop': early_stopping.early_stop
                    },
                    'feature_info': {
                        'static_cols': static_cols,
                        'dynamic_cols_no_target': dynamic_cols_no_target,
                        'dynamic_input_size': dynamic_input_size,
                        'static_input_size': static_input_size
                    }
                }, checkpoint_path)
                print(f"Checkpoint saved to {checkpoint_path}")
    
    # Final evaluation
    log_message("\n=== Final evaluation on test set ===", force_flush=True)
    gc.collect()
    
    # Load best model
    checkpoint = torch.load(os.path.join(save_dir, 'best_model.pt'), map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    
    test_loss, test_kge = evaluate_model_gpu_optimized(
        model, test_watersheds, bucket_name, base_data_dir, base_attr_dir,
        static_cols, dynamic_cols_no_target, sequence_length, 
        batch_size, device, chunk_size, norm_params=train_norm_params
    )
    
    log_message(f"Test set - Loss: {test_loss:.4f}, KGE: {test_kge:.4f}", force_flush=True)
    
    # Save final results
    results = {
        'test_loss': test_loss,
        'test_kge': test_kge,
        'best_val_loss': checkpoint['val_loss'],
        'best_val_kge': checkpoint.get('val_kge', np.nan),
        'best_epoch': checkpoint['epoch'],
        'total_epochs_trained': training_history['epoch'][-1] if training_history['epoch'] else 0,
        'early_stopped': early_stopping.early_stop,
        'training_history': training_history,
        'feature_info': {
            'static_features': static_cols,
            'dynamic_features': dynamic_cols_no_target,
            'target': 'streamflow'
        }
    }
    
    with open(os.path.join(save_dir, 'results.json'), 'w') as f:
        json.dump(results, f, indent=2)
    
    # Enhanced plotting with learning rate
    plot_training_history_with_lr(training_history, save_dir)
    
    print_memory_usage("Final state")
    
    log_message("\n=== Training Summary ===")
    log_message(f"Total gradient explosions: {gradient_explosion_count}")
    log_message(f"Total NaN batches skipped: {nan_batch_count}")
    log_message(f"Final learning rate: {optimizer.param_groups[0]['lr']}")
    log_message(f"Best validation KGE: {early_stopping.best_score}")
    
    log_message("\nTraining completed successfully!", force_flush=True)
    
    return results




def sigterm_handler(signum, frame):
    """
    Handle SIGTERM signal from spot instance termination.
    This gives us ~2 minutes to save our progress before shutdown.
    """
    print("\n!!! SIGTERM received - Spot instance terminating !!!")
    print("Attempting to save emergency checkpoint...")
    
    # Try to access global training state if available
    # In your case, these would need to be made accessible
    checkpoint_path = '/opt/ml/checkpoints/emergency_checkpoint.pt'
    
    # Since we can't easily access the training loop variables here,
    # we'll create a marker file that the training loop can check
    with open('/opt/ml/checkpoints/TERMINATION_REQUESTED', 'w') as f:
        f.write('SIGTERM received')
    
    print(f"Termination marker created. Training loop will save state on next iteration.")
    # Don't exit immediately - let the training loop handle it











def evaluate_model_gpu_optimized(model, watersheds_df, bucket_name, base_data_dir, base_attr_dir,
                   static_cols, dynamic_cols, sequence_length, batch_size, 
                   device, chunk_size=50, norm_params=None):
    """Evaluate model with GPU support"""
    model.eval()
    all_losses = []
    all_predictions = []
    all_targets = []
    
    chunks = [watersheds_df[i:i+chunk_size] 
              for i in range(0, len(watersheds_df), chunk_size)]
    
    with torch.no_grad():
        for chunk_idx, chunk_watersheds in enumerate(chunks):
            try:
                chunk_df, _, _, _ = prepare_ptf_dataframe(
                    chunk_watersheds, bucket_name, base_data_dir, base_attr_dir, norm_params=norm_params
                )
                
                if chunk_df.empty:
                    continue
                
                dataset = WatershedDataset(
                    chunk_df, static_cols, dynamic_cols, 
                    'streamflow', sequence_length
                )
                
                del chunk_df
                gc.collect()
                
                if len(dataset) == 0:
                    del dataset
                    continue
                
                loader = DataLoader(
                    dataset, batch_size=batch_size, 
                    shuffle=False, 
                    num_workers=4,  # Reduced for GPU
                    pin_memory=True  # Important for GPU
                )
                
                chunk_losses = []
                for batch_idx, (dynamic_seq, static_feat, target) in enumerate(loader):
                    # Move to GPU
                    dynamic_seq = dynamic_seq.to(device)
                    static_feat = static_feat.to(device)
                    target = target.to(device)
                    
                    predictions = model(dynamic_seq, static_feat)
                    loss = nn.MSELoss()(predictions, target)
                    
                    chunk_losses.append(loss.item())
                    # Move back to CPU for storage
                    all_predictions.extend(predictions.cpu().numpy())
                    all_targets.extend(target.cpu().numpy())
                    
                    del predictions, loss
                    
                    if batch_idx % 20 == 0:
                        if device == 'cuda':
                            torch.cuda.empty_cache()
                
                all_losses.extend(chunk_losses)
                
                del dataset, loader
                gc.collect()
                if device == 'cuda':
                    torch.cuda.empty_cache()
                
            except Exception as e:
                print(f"Error evaluating chunk {chunk_idx}: {str(e)}")
                gc.collect()
                continue
    
    mean_loss = np.mean(all_losses) if all_losses else float('inf')
    # Before calculating metrics (line ~1863)
    if norm_params is not None:
        # Unnormalize predictions and targets for KGE calculation
        predictions_original = unnormalize_predictions(
            np.array(all_predictions).flatten(), 
            norm_params
        )
        targets_original = unnormalize_predictions(
            np.array(all_targets).flatten(), 
            norm_params
        )
        
        # Calculate KGE on original scale
        kge = calculate_kge(targets_original, predictions_original) if all_targets else float('-inf')
    else:
        # Fallback to normalized KGE
        kge = calculate_kge(
            np.array(all_targets).flatten(), 
            np.array(all_predictions).flatten()
        ) if all_targets else float('-inf')

    mean_loss = np.mean(all_losses) if all_losses else float('inf')
    
    return mean_loss, kge









def save_training_state(save_path, model, optimizer, scheduler, training_history, 
                       chunk_idx, epoch, chunks_processed, early_stopping, 
                       static_cols, dynamic_cols_no_target, additional_info=None):
    """Save complete training state to disk"""
    state = {
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict() if scheduler else None,
        'training_history': training_history,
        'chunk_idx': chunk_idx,
        'current_epoch': epoch,
        'chunks_processed': chunks_processed,
        'early_stopping_state': {
            'counter': early_stopping.counter,
            'best_score': early_stopping.best_score,
            'val_loss_min': early_stopping.val_loss_min,
            'early_stop': early_stopping.early_stop
        },
        'feature_info': {
            'static_cols': static_cols,
            'dynamic_cols_no_target': dynamic_cols_no_target,
            'dynamic_input_size': len(dynamic_cols_no_target),
            'static_input_size': len(static_cols)
        },
        'timestamp': datetime.now().isoformat()
    }
    
    if additional_info:
        state.update(additional_info)
    
    # Save with atomic write (write to temp, then rename)
    temp_path = save_path + '.tmp'
    torch.save(state, temp_path)
    os.rename(temp_path, save_path)
    
    print(f"Saved training state to {save_path}")
    return save_path



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

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    
    # SageMaker specific arguments
    parser.add_argument('--model-dir', type=str, default=os.environ.get('SM_MODEL_DIR'))
    parser.add_argument('--output-data-dir', type=str, default=os.environ.get('SM_OUTPUT_DATA_DIR'))
    
    # Your training arguments
    parser.add_argument('--bucket-name', type=str, required=True)
    parser.add_argument('--base-data-dir', type=str, required=True)
    parser.add_argument('--base-attr-dir', type=str, required=True)
    parser.add_argument('--experiment-name', type=str, default='experiment_1')
    parser.add_argument('--chunk-size', type=int, default=25)
    parser.add_argument('--epochs-per-chunk', type=int, default=100)#10)
    parser.add_argument('--batch-size', type=int, default=256)#32)
    parser.add_argument('--hidden-size', type=int, default=128)
    parser.add_argument('--num-layers', type=int, default=2)
    parser.add_argument('--total-epochs', type=int, default=500)
    parser.add_argument('--learning-rate', type=float, default=0.0001)
    parser.add_argument('--dropout', type=float, default=0.3)
    
    args = parser.parse_args()
    
    # Register SIGTERM handler for spot instance termination
    signal.signal(signal.SIGTERM, sigterm_handler)
    print("SIGTERM handler registered for spot instance termination")
    
    # Create hyperparameters dict
    hyperparameters = {
        'sequence_length': 365,
        'hidden_size': args.hidden_size,
        'learning_rate': args.learning_rate,
        'batch_size': args.batch_size,
        'dropout': args.dropout,
        'num_layers': args.num_layers,
        'total_epochs': args.total_epochs
    }
    
    lr_scheduler_params = {
        'factor': 0.1,  # More aggressive reduction (was 0.5)
        'patience': 3,  # Less patience (was 5)
        'min_lr': 1e-7,  # Lower minimum (was 1e-6)
        'threshold': 0.01,  # Add threshold for significant improvement
        'cooldown': 5  # Add cooldown period after reduction
    }

    
    # Load watersheds data
    # You'll need to modify this to load from S3 or pass as argument
    watershed_df = identify_all_available_watersheds(args.bucket_name, args.base_data_dir)
    watershed_df = watershed_df[watershed_df['subdirectory_name'] == 'camels']
    watershed_df = watershed_df.iloc[0:20]
    
    # Use the output directory for saving (this persists after training)
    save_dir = args.output_data_dir


        # Detect if GPU is available
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")
    if device == 'cuda':
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")


    
    # Run your training function
    results = train_ea_lstm_model_with_scheduling(
        watersheds_df=watershed_df,
        bucket_name=args.bucket_name,
        base_data_dir=args.base_data_dir,
        base_attr_dir=args.base_attr_dir,
        hyperparameters=hyperparameters,
        save_dir=save_dir,
        chunk_size=args.chunk_size,
        epochs_per_chunk=args.epochs_per_chunk,
        gradient_accumulation_steps=8,#2
        early_stopping_patience=10,
        lr_scheduler_type='reduce_on_plateau',
        lr_scheduler_params=lr_scheduler_params,
        target_noise_std=0.1,
        target_noise_type='multiplicative', #'log-normal',
        device=device, 
    )
    
    # Save final model to model directory
    if os.path.exists(os.path.join(save_dir, 'best_model.pt')):
        import shutil
        shutil.copy(
            os.path.join(save_dir, 'best_model.pt'),
            os.path.join(args.model_dir, 'model.pt')
        )
    
    print("Training completed successfully!")

if __name__ == '__main__':
    # Copy all your function definitions here
    # (all the functions from your lstmCode.txt)
    
    main()