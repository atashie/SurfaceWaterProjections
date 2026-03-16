"""
Human interference metadata module.

Loads GAGES-II and HYDAT metadata for watershed characterization.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import numpy as np

from .config import (
    INCLUDE_HUMAN_INTERFERENCE,
    GAGES_II_DIR,
    HYDAT_PATH,
    INTERFERENCE_COLUMNS,
)


# GAGES-II file definitions (matching R's config.R)
GAGES_II_FILES_CONUS = {
    "hydromod_dams": "conterm_hydromod_dams.txt",
    "pop_infrastr": "conterm_pop_infrastr.txt",
    "hydromod_other": "conterm_hydromod_other.txt",
    "bas_classif": "conterm_bas_classif.txt",
    "lc06_basin": "conterm_lc06_basin.txt",
}

GAGES_II_FILES_AKHIPR = {
    "hydromod_dams": "AKHIPR_hydromod_dams.txt",
    "pop_infrastr": "AKHIPR_pop_infrastr.txt",
    "hydromod_other": "AKHIPR_hydromod_other.txt",
    "bas_classif": "AKHIPR_bas_classif.txt",
    # Note: AKHIPR does not have lc06_basin file
}

# Sentinel value for missing data in GAGES-II
GAGES_II_MISSING_VALUE = -999


def load_gages_ii_interference(gages_dir: str = GAGES_II_DIR) -> pd.DataFrame:
    """
    Load GAGES-II human interference metadata.

    Parameters
    ----------
    gages_dir : str
        Path to GAGES-II metadata directory

    Returns
    -------
    pd.DataFrame
        GAGES-II interference columns for USGS gages
    """
    if gages_dir is None or gages_dir == "":
        print("Warning: GAGES-II directory not configured")
        return pd.DataFrame(columns=[
            "STAID", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
            "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
            "HYDRO_DISTURB_INDX", "CLASS"
        ])

    if not os.path.isdir(gages_dir):
        print(f"Warning: GAGES-II directory not found: {gages_dir}")
        return pd.DataFrame()

    # Load CONUS data
    conus_data = _load_gages_ii_region(gages_dir, GAGES_II_FILES_CONUS)

    # Load AKHIPR data
    akhipr_data = _load_gages_ii_region(gages_dir, GAGES_II_FILES_AKHIPR)

    # Combine
    if len(conus_data) > 0 and len(akhipr_data) > 0:
        combined = pd.concat([conus_data, akhipr_data], ignore_index=True)
    elif len(conus_data) > 0:
        combined = conus_data
    elif len(akhipr_data) > 0:
        combined = akhipr_data
    else:
        return pd.DataFrame()

    # Replace -999 (missing value sentinel) with NaN
    numeric_cols = combined.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        combined.loc[combined[col] == GAGES_II_MISSING_VALUE, col] = np.nan

    # Standardize CLASS values
    if "CLASS" in combined.columns:
        combined["CLASS"] = combined["CLASS"].str.strip()
        combined.loc[combined["CLASS"] == "", "CLASS"] = np.nan

    return combined


def _load_gages_ii_region(gages_dir: str, files: Dict[str, str]) -> pd.DataFrame:
    """
    Load GAGES-II data for a specific region (CONUS or AKHIPR).

    Parameters
    ----------
    gages_dir : str
        Path to GAGES-II metadata directory
    files : dict
        Dictionary mapping file types to filenames

    Returns
    -------
    pd.DataFrame
        Combined data from all region files
    """
    result = pd.DataFrame()

    for file_type, filename in files.items():
        filepath = os.path.join(gages_dir, filename)
        if not os.path.isfile(filepath):
            continue

        try:
            df = pd.read_csv(
                filepath,
                na_values=["-999", "NA", ""],
                dtype={"STAID": str},
                encoding="latin-1",
            )

            # Ensure STAID is string
            if "STAID" in df.columns:
                df["STAID"] = df["STAID"].astype(str)

            if len(result) == 0:
                result = df
            else:
                # Merge on STAID
                result = result.merge(df, on="STAID", how="left", suffixes=("", "_dup"))
                # Remove duplicate columns
                dup_cols = [c for c in result.columns if c.endswith("_dup")]
                result = result.drop(columns=dup_cols)

        except Exception as e:
            print(f"Warning: Failed to load {filename}: {e}")

    # Select relevant columns
    relevant_cols = [
        "STAID", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
        "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
        "HYDRO_DISTURB_INDX", "CLASS"
    ]
    available_cols = [c for c in relevant_cols if c in result.columns]

    if len(available_cols) > 0:
        return result[available_cols]
    else:
        return pd.DataFrame()


def load_canadian_interference(hydat_path: Optional[str] = HYDAT_PATH) -> pd.DataFrame:
    """
    Load Canadian HYDAT interference metadata.

    For Canadian gages, returns RHBN and REGULATED status.
    Note: This requires a pre-exported CSV since tidyhydat is an R package.

    Parameters
    ----------
    hydat_path : str or None
        Path to HYDAT metadata CSV (pre-exported from tidyhydat)

    Returns
    -------
    pd.DataFrame
        Canadian interference columns
    """
    if hydat_path is None or hydat_path == "":
        print("Warning: HYDAT path not configured - Canadian interference metadata unavailable")
        return pd.DataFrame(columns=[
            "gage_id", "RHBN", "REGULATED", "human_interference_class"
        ])

    if not os.path.isfile(hydat_path):
        print(f"Warning: HYDAT file not found: {hydat_path}")
        return pd.DataFrame(columns=[
            "gage_id", "RHBN", "REGULATED", "human_interference_class"
        ])

    try:
        df = pd.read_csv(hydat_path)

        # Ensure required columns exist
        if "STATION_NUMBER" not in df.columns:
            print("Warning: HYDAT file missing STATION_NUMBER column")
            return pd.DataFrame()

        # Rename column
        df = df.rename(columns={"STATION_NUMBER": "gage_id"})

        # Calculate human_interference_class
        def classify_interference(row):
            rhbn = row.get("RHBN", np.nan)
            if pd.notna(rhbn) and rhbn == True:
                return "reference"
            elif pd.notna(rhbn) and rhbn == False:
                return "non-reference"
            else:
                return "unknown"

        df["human_interference_class"] = df.apply(classify_interference, axis=1)

        return df

    except Exception as e:
        print(f"Warning: Failed to load HYDAT data: {e}")
        return pd.DataFrame()


def enrich_signatures_with_metadata(
    signatures: pd.DataFrame,
    metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    Add human interference columns to signature output.

    Parameters
    ----------
    signatures : pd.DataFrame
        Signature results with gage_id column
    metadata : pd.DataFrame
        Watershed metadata with interference columns

    Returns
    -------
    pd.DataFrame
        Signatures with interference columns added
    """
    if not INCLUDE_HUMAN_INTERFERENCE:
        return signatures

    # Ensure both have gage_id column
    if "gage_id" not in signatures.columns or "gage_id" not in metadata.columns:
        print("Warning: Missing gage_id column - cannot enrich with metadata")
        return signatures

    # Select interference columns from metadata
    interference_cols = [c for c in INTERFERENCE_COLUMNS if c in metadata.columns]
    cols_to_join = ["gage_id"] + interference_cols

    if len(cols_to_join) <= 1:
        print("Warning: No interference columns found in metadata")
        return signatures

    # Join metadata to signatures
    metadata_subset = metadata[cols_to_join].copy()
    result = signatures.merge(metadata_subset, on="gage_id", how="left")

    return result
