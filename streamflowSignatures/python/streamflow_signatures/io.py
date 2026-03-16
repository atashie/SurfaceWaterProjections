"""
I/O functions for reading and writing streamflow signature data.

Handles parquet input, CSV output, and schema validation.
"""

from typing import Optional, List, Dict, Tuple, Union
import calendar
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

from .config import MIN_NUM_YEARS, MIN_FRAC_GOOD_DATA, MIN_Q_VALUE, MIN_DAYS_ABOVE_THRESHOLD


# Schema constants matching R config.R
STAT_SUFFIXES = [
    "_senn_slp",
    "_linear_slp",
    "_spearman_rho",
    "_spearman_pval",
    "_mk_rho",
    "_mk_pval",
    "_mean",
    "_median",
]

# Expected signature base names (matches R config.R EXPECTED_SIGNATURE_BASES)
EXPECTED_SIGNATURE_BASES = [
    # Flow volumes (22 metrics)
    "Qann", "Qwin", "Qspr", "Qsum", "Qfal",
    "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
    "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99", "Q95_Q10",
    # FDC (3 metrics)
    "FDCall", "FDC90th", "FDCmid",
    # Baseflow (2 metrics)
    "BFI_Eckhardt", "BFI_LyneHollick",
    # Recession (5 metrics + 6 seasonality)
    "log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity",
    # Pulses (10 metrics - year-specific + 4 period-of-record)
    "n_high_pulses_year", "n_low_pulses_year",
    "dur_high_pulses_year", "dur_low_pulses_year",
    "n_high_pulses_all", "n_low_pulses_all",
    "dur_high_pulses_all", "dur_low_pulses_all",
    "TQmean",
    "Flow_Reversals_annual", "Flow_Reversals_winter", "Flow_Reversals_spring",
    "Flow_Reversals_summer", "Flow_Reversals_fall",
    "Flow_Reversals_rising_annual", "Flow_Reversals_falling_annual",
    "Flow_Reversals_rising_winter", "Flow_Reversals_falling_winter",
    "Flow_Reversals_rising_spring",
    # Flashiness (1 metric)
    "flashinessRB",
    # Timing (13 metrics)
    "D5_day", "D10_day", "D20_day", "D30_day", "D40_day", "D50_day",
    "D60_day", "D70_day", "D80_day", "D90_day", "D95_day",
    "D25_to_D75", "Dmax",
    # Climate-dependent signatures
    "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
    "summer_runoff_ratio", "fall_runoff_ratio",
    "elasticity",
    "qp_slope_sd", "qp_bimodality",
    "avg_storage",
]

# Metadata columns expected in output
METADATA_COLUMNS = [
    "gage_id",
    "latitude",
    "longitude",
    "basin_area_km2",
    "gage_type",
    "processing_status",
    "num_water_years",
    "start_year",
    "end_year",
]


def read_parquet(
    file_path: Union[str, Path],
    columns: Optional[List[str]] = None,
    gage_ids: Optional[List[str]] = None,
    normalize_columns: bool = True,
) -> pd.DataFrame:
    """
    Read streamflow data from a parquet file.

    Parameters
    ----------
    file_path : str or Path
        Path to the parquet file
    columns : list of str, optional
        Specific columns to read. If None, reads all columns.
    gage_ids : list of str, optional
        Filter to specific gage IDs. If None, reads all gages.
    normalize_columns : bool, default True
        Auto-rename common column variants to standard names:
        Date->date, site_id->gage_id, prcp->PPT. Also ensures
        date column is pd.Timestamp.

    Returns
    -------
    pd.DataFrame
        Streamflow data with columns including 'gage_id', 'date', 'Q',
        'water_year', 'month', 'dowy' (day of water year)

    Raises
    ------
    FileNotFoundError
        If the parquet file does not exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    # Read parquet
    if columns is not None:
        df = pq.read_table(file_path, columns=columns).to_pandas()
    else:
        df = pq.read_table(file_path).to_pandas()

    # Normalize column names
    if normalize_columns:
        rename_map = {}
        if "Date" in df.columns and "date" not in df.columns:
            rename_map["Date"] = "date"
        if "site_id" in df.columns and "gage_id" not in df.columns:
            rename_map["site_id"] = "gage_id"
        if "prcp" in df.columns and "PPT" not in df.columns:
            rename_map["prcp"] = "PPT"
        if rename_map:
            df = df.rename(columns=rename_map)

        # Ensure date column is datetime
        if "date" in df.columns and not pd.api.types.is_datetime64_any_dtype(df["date"]):
            df["date"] = pd.to_datetime(df["date"])

    # Filter to specific gages if requested
    if gage_ids is not None and "gage_id" in df.columns:
        df = df[df["gage_id"].isin(gage_ids)]

    return df


def write_signatures(
    signatures: Dict[str, Dict[str, float]],
    metadata: pd.DataFrame,
    output_path: Union[str, Path],
) -> None:
    """
    Write signature results to CSV.

    Parameters
    ----------
    signatures : dict
        Dictionary mapping gage_id to signature statistics dict
    metadata : pd.DataFrame
        Metadata DataFrame with gage information
    output_path : str or Path
        Output CSV path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert signatures dict to DataFrame
    sig_df = pd.DataFrame.from_dict(signatures, orient="index")
    sig_df.index.name = "gage_id"
    sig_df = sig_df.reset_index()

    # Merge with metadata
    if metadata is not None and len(metadata) > 0:
        # Ensure gage_id is string in both
        sig_df["gage_id"] = sig_df["gage_id"].astype(str)
        metadata = metadata.copy()
        metadata["gage_id"] = metadata["gage_id"].astype(str)

        # Select metadata columns that exist
        meta_cols = [c for c in METADATA_COLUMNS if c in metadata.columns]
        result = sig_df.merge(metadata[meta_cols], on="gage_id", how="left")
    else:
        result = sig_df

    # Write to CSV
    result.to_csv(output_path, index=False)


def validate_schema(
    df: pd.DataFrame,
    required_columns: Optional[List[str]] = None,
    signature_bases: Optional[List[str]] = None,
) -> Dict[str, Union[bool, List[str]]]:
    """
    Validate that a DataFrame has the expected schema.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate
    required_columns : list of str, optional
        Required column names. Defaults to METADATA_COLUMNS.
    signature_bases : list of str, optional
        Expected signature base names. Defaults to EXPECTED_SIGNATURE_BASES.

    Returns
    -------
    dict
        Validation result with keys:
            - valid: bool, True if all validations pass
            - missing_metadata: list of missing metadata columns
            - missing_signatures: list of missing signature base names
            - extra_columns: list of unexpected columns
    """
    if required_columns is None:
        required_columns = METADATA_COLUMNS
    if signature_bases is None:
        signature_bases = EXPECTED_SIGNATURE_BASES

    # Check metadata columns
    missing_metadata = [c for c in required_columns if c not in df.columns]

    # Check signature columns
    missing_signatures = []
    for base in signature_bases:
        expected_cols = [f"{base}{suffix}" for suffix in STAT_SUFFIXES]
        if not any(c in df.columns for c in expected_cols):
            missing_signatures.append(base)

    # Identify extra columns
    expected_all = set(required_columns)
    for base in signature_bases:
        expected_all.update(f"{base}{suffix}" for suffix in STAT_SUFFIXES)
    extra_columns = [c for c in df.columns if c not in expected_all]

    return {
        "valid": len(missing_metadata) == 0 and len(missing_signatures) == 0,
        "missing_metadata": missing_metadata,
        "missing_signatures": missing_signatures,
        "extra_columns": extra_columns,
    }


def add_water_year_columns(df: pd.DataFrame, date_col: str = "date") -> pd.DataFrame:
    """
    Add water year related columns to a streamflow DataFrame.

    Water year runs from October 1 to September 30. Day 1 is October 1.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a date column
    date_col : str, default "date"
        Name of the date column

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns:
            - water_year: int, the water year (year of the end date)
            - month: int, calendar month (1-12)
            - dowy: int, day of water year (1-366)
    """
    # Auto-detect date column if specified name not found
    if date_col not in df.columns:
        variants = {c for c in df.columns if c.lower() == date_col.lower()}
        if len(variants) == 1:
            date_col = variants.pop()
        else:
            raise KeyError(
                f"Date column '{date_col}' not found. Available columns: {list(df.columns)}"
            )

    df = df.copy()

    # Ensure date column is datetime
    if not pd.api.types.is_datetime64_any_dtype(df[date_col]):
        df[date_col] = pd.to_datetime(df[date_col])

    # Extract month
    df["month"] = df[date_col].dt.month

    # Calculate water year (Oct-Dec = next year, Jan-Sep = current year)
    df["water_year"] = df[date_col].dt.year
    df.loc[df["month"] >= 10, "water_year"] += 1

    # Calculate day of water year (vectorized)
    # Water year starts Oct 1, so subtract Oct 1 of (water_year - 1)
    wy_start = pd.to_datetime(
        (df["water_year"] - 1).astype(str) + "-10-01"
    )
    df["dowy"] = (df[date_col] - wy_start).dt.days + 1

    return df


def filter_qualifying_years(
    gage_data: pd.DataFrame,
    min_q_value: float = MIN_Q_VALUE,
    min_days_above: int = MIN_DAYS_ABOVE_THRESHOLD,
    min_frac_good: float = MIN_FRAC_GOOD_DATA,
    min_num_years: int = MIN_NUM_YEARS,
) -> Tuple[List[int], bool]:
    """Filter water years per-gage matching R's process_signatures_from_parquet().

    Three-stage per-year filtering:
    1. Per water year, check at least min_days_above days with Q > min_q_value
    2. Per water year, check data completeness (>= min_frac_good of expected days)
    3. Gage qualifies if >= min_num_years pass both sub-checks

    Parameters
    ----------
    gage_data : pd.DataFrame
        DataFrame for a single gage with water_year and Q columns.
    min_q_value : float
        Minimum Q threshold for sub-check 1.
    min_days_above : int
        Minimum days above min_q_value per water year.
    min_frac_good : float
        Minimum fraction of expected days with non-NA Q per water year.
    min_num_years : int
        Minimum number of qualifying water years for gage to pass.

    Returns
    -------
    tuple of (list[int], bool)
        (qualifying_years, gage_qualifies)
    """
    qualifying_years = []
    for wy in gage_data["water_year"].unique():
        yr_data = gage_data[gage_data["water_year"] == wy]
        q = yr_data["Q"]
        n_nona = q.notna().sum()

        # Sub-check 1: days above threshold
        n_above = (q > min_q_value).sum()
        if n_above < min_days_above:
            continue

        # Sub-check 2: data completeness (accounting for leap years)
        # Water year Y spans Oct 1 (Y-1) to Sep 30 (Y)
        expected_days = 366 if calendar.isleap(wy) else 365
        min_good_days = int(expected_days * min_frac_good)
        if n_nona < min_good_days:
            continue

        qualifying_years.append(wy)

    return qualifying_years, len(qualifying_years) >= min_num_years
