"""
I/O functions for reading and writing streamflow signature data.

Handles parquet input, CSV output, and schema validation.

Schema model (August 2026 product — 1,653 columns):

    1,653 = 20 metadata columns
          + 100 signature bases x 16 statistic suffixes   (= 1,600)
          + 21 per-gage scalars
          + 12 QA flag columns

The 16 suffixes are the 8 ordinary statistics (``stats.STAT_SUFFIXES``) plus
the 8 Pettitt changepoint fields (``stats.CP_SUFFIXES``). The 20 metadata
columns mirror ``metadata_order`` in docs/benchmarks/run_python_benchmark.py
(the production runner's output assembly); the 12 flags come from
``qa_qc.get_flag_columns()``. Use ``expected_output_columns()`` for the full
list and ``validate_schema(df, strict=True)`` to gate a complete product CSV.
"""

from typing import Optional, List, Dict, Tuple, Union
import calendar
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

import numpy as np

from .config import (
    MIN_NUM_YEARS, MIN_FRAC_GOOD_DATA, MIN_Q_VALUE, MIN_DAYS_ABOVE_THRESHOLD,
    NA_MAX_GAP_DAYS, NA_MAX_RAW_NA_PER_YEAR, NA_REJECT_NEGATIVE_FLOW,
    NA_REJECT_RESIDUAL_NA, NA_CONSTANT_SD_ENABLED, NA_CONSTANT_SD_MIN_DAYS,
    NA_CONSTANT_SD_MAX_UNIQUE, NA_SEASONAL_MIN_FRACTION, NA_SEASONAL_DEFINITIONS,
    NA_MAX_RAW_NA_PPT, NA_MAX_GAP_PPT, NA_REJECT_NEGATIVE_PPT,
    NA_MAX_RAW_NA_SWE, NA_MAX_GAP_SWE, NA_REJECT_NEGATIVE_SWE,
)
from .stats import STAT_SUFFIXES, CP_SUFFIXES
from .qa_qc import get_flag_columns


# The 16 per-base statistic suffixes of the full product: the 8 ordinary
# statistics plus the 8 Pettitt changepoint fields.
PETTITT_SUFFIXES = list(CP_SUFFIXES)
ALL_STAT_SUFFIXES = list(STAT_SUFFIXES) + PETTITT_SUFFIXES

# Expected signature base names — the 100 time-series signature bases of the
# current (Aug 2026) 1,653-column product. Mirrors the Julia annual-collector
# registry: EXPECTED_DENSE_SIGNATURES in julia/test/test_annual_collector.jl
# (78 dense bases) plus the sparse recession / parameterized-BFI bases
# exercised there (8) plus the 14 snow metrics (julia/src/snow.jl SNOW_METRICS).
EXPECTED_SIGNATURE_BASES = [
    # Flow volumes (21 metrics)
    "Qann", "Qwin", "Qspr", "Qsum", "Qfal",
    "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
    "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99", "Q95_Q10",
    # FDC (3 metrics)
    "FDCall", "FDC90th", "FDCmid",
    # Baseflow fixed-parameter (2 metrics)
    "BFI_Eckhardt", "BFI_LyneHollick",
    # Baseflow recession-parameterized (2 metrics)
    "BFI_Eckhardt_param", "BFI_LyneHollick_param",
    # Recession (7 annual metrics; the 6 seasonality scalars are not
    # time-series bases)
    "log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events",
    "concavity", "n_recession_events", "alpha_linear",
    # Pulses (14 metrics)
    "n_high_pulses_year", "n_low_pulses_year",
    "dur_high_pulses_year", "dur_low_pulses_year",
    "n_high_pulses_all", "n_low_pulses_all",
    "dur_high_pulses_all", "dur_low_pulses_all",
    "TQmean",
    "Flow_Reversals_annual", "Flow_Reversals_winter", "Flow_Reversals_spring",
    "Flow_Reversals_summer", "Flow_Reversals_fall",
    # Negative flow days (1 metric)
    "negative_ann",
    # Flashiness (1 metric)
    "flashinessRB",
    # Timing (15 metrics: 13 D-days D1-D99 + D25_to_D75 + Dmax)
    "D1_day", "D5_day", "D10_day", "D20_day", "D30_day", "D40_day", "D50_day",
    "D60_day", "D70_day", "D80_day", "D90_day", "D95_day", "D99_day",
    "D25_to_D75", "Dmax",
    # Runoff ratios (5 metrics)
    "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
    "summer_runoff_ratio", "fall_runoff_ratio",
    # Elasticity (2 time-series metrics; elasticity_static is a scalar)
    "elasticity_rolling", "elasticity_annual",
    # Q-P seasonality (2 metrics)
    "qp_slope_sd", "qp_bimodality",
    # Storage (1 metric)
    "avg_storage",
    # Snow (14 metrics)
    "swe_max", "swe_max_dowy", "snow_cover_days", "snow_on_dowy",
    "snow_off_dowy", "melt_season_days", "melt_rate", "ssm", "swe_apr1",
    "melt_before_peak", "melt_before_peak_pct", "melt_before_peak_to_max_swe",
    "melt_com_dowy", "swe_max_to_ppt",
    # Streamflow drought (10 metrics)
    "drought_duration_fixed_p2", "drought_duration_fixed_p5",
    "drought_duration_fixed_p10", "drought_duration_fixed_p20",
    "drought_duration_fixed_p30",
    "drought_deficit_fixed_p2", "drought_deficit_fixed_p5",
    "drought_deficit_fixed_p10", "drought_deficit_fixed_p20",
    "drought_deficit_fixed_p30",
]

# Metadata columns of the production output (20), in the exact order the
# benchmark runner emits them — mirrors `metadata_order` in
# docs/benchmarks/run_python_benchmark.py: gage id / coords / area / type /
# record info / area_normalized, the 8 GAGES-II interference attributes, the
# 2 Canadian HYDAT interference fields, and the unified classification.
METADATA_COLUMNS = [
    "gage_id",
    "latitude",
    "longitude",
    "basin_area",
    "gage_type",
    "num_water_years",
    "start_water_year",
    "end_water_year",
    "area_normalized",
    "NDAMS_2009",
    "MAJ_DDENS_2009",
    "STOR_NID_2009",
    "IMPNLCD06",
    "DEVNLCD06",
    "FRESHW_WITHDRAWAL",
    "HYDRO_DISTURB_INDX",
    "CLASS",
    "RHBN",
    "REGULATED",
    "human_interference_class",
]

# Per-gage scalar outputs (21) — documented exceptions to the 8-statistic rule.
# Each name is emitted as-is (no suffix). Sources: elasticity.py, recession.py,
# runoff_ratios.py, signatures.py (season exclusion counts), drought.py
# (DROUGHT_THRESHOLD_SCALARS), and the benchmark runner (ice_affected_days_total,
# aggregated from preprocessor diagnostics).
PER_GAGE_SCALARS = [
    "elasticity_static",
    "log_a_seasonality_amplitude_all",
    "log_a_seasonality_amplitude_first_half",
    "log_a_seasonality_amplitude_last_half",
    "log_a_seasonality_minimum_all",
    "log_a_seasonality_minimum_first_half",
    "log_a_seasonality_minimum_last_half",
    "runoff_ratio_high_count",
    "elasticity_years_total",
    "elasticity_years_low_ppt",
    "ice_affected_days_total",
    "recession_alpha_point_cloud_linear_reservoir",
    "season_excluded_years_winter",
    "season_excluded_years_spring",
    "season_excluded_years_summer",
    "season_excluded_years_fall",
    "drought_threshold_fixed_p2",
    "drought_threshold_fixed_p5",
    "drought_threshold_fixed_p10",
    "drought_threshold_fixed_p20",
    "drought_threshold_fixed_p30",
]

# QA/QC flag columns (12) — single source of truth is qa_qc.get_flag_columns().
QA_FLAG_COLUMNS = get_flag_columns()


def expected_output_columns(
    metadata_columns: Optional[List[str]] = None,
    signature_bases: Optional[List[str]] = None,
    scalars: Optional[List[str]] = None,
    flag_columns: Optional[List[str]] = None,
) -> List[str]:
    """
    Return the full expected column list of the production signatures CSV.

    With the defaults this is the 1,653-column August 2026 product:
    20 metadata + 100 bases x 16 suffixes + 21 scalars + 12 QA flags,
    ordered as the benchmark runner writes it (metadata first, then the
    signature + scalar columns sorted, then the flags sorted).
    """
    if metadata_columns is None:
        metadata_columns = METADATA_COLUMNS
    if signature_bases is None:
        signature_bases = EXPECTED_SIGNATURE_BASES
    if scalars is None:
        scalars = PER_GAGE_SCALARS
    if flag_columns is None:
        flag_columns = QA_FLAG_COLUMNS

    signature_cols = [
        f"{base}{suffix}" for base in signature_bases for suffix in ALL_STAT_SUFFIXES
    ]
    return (
        list(metadata_columns)
        + sorted(signature_cols + list(scalars))
        + sorted(flag_columns)
    )


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
        if "swe" in df.columns and "SWE" not in df.columns:
            rename_map["swe"] = "SWE"
        if rename_map:
            df.rename(columns=rename_map, inplace=True)

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
    strict: bool = False,
) -> Dict[str, Union[bool, List[str]]]:
    """
    Validate a signatures DataFrame against the full product schema.

    The expected column set is metadata + every signature base x all 16
    statistic suffixes (8 ordinary stats + 8 Pettitt changepoint fields) +
    the 21 per-gage scalars + the 12 QA flag columns — 1,653 columns with
    the defaults (see the module docstring).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate (only its columns are inspected).
    required_columns : list of str, optional
        Required metadata column names. Defaults to METADATA_COLUMNS (20).
    signature_bases : list of str, optional
        Expected signature base names. Defaults to EXPECTED_SIGNATURE_BASES (100).
    strict : bool, default False
        strict=True: every expected column must be present and no unexpected
        column may exist — the gate for a complete production CSV.
        strict=False: missing metadata / Pettitt fields / scalars / flags are
        reported but do not fail validation, so partial outputs (e.g. a bare
        ``calculate_all_signatures`` result without changepoint or metadata)
        still validate as long as no signature base is entirely absent.

    Returns
    -------
    dict
        Validation result with keys:
            - valid: bool. strict=True: no missing and no unexpected columns.
              strict=False: no signature base entirely absent.
            - missing_metadata: metadata columns absent from df
            - missing_signatures: bases with NONE of their 16 fields present
            - missing_columns: every expected column absent from df (sorted)
            - extra_columns: columns of df outside the expected set
            - n_expected: total size of the expected column set
    """
    if required_columns is None:
        required_columns = METADATA_COLUMNS
    if signature_bases is None:
        signature_bases = EXPECTED_SIGNATURE_BASES

    expected = expected_output_columns(
        metadata_columns=required_columns, signature_bases=signature_bases
    )
    expected_set = set(expected)
    present = set(df.columns)

    # Metadata columns absent from df
    missing_metadata = [c for c in required_columns if c not in present]

    # Bases with none of their 16 fields present (entirely absent signature)
    missing_signatures = []
    for base in signature_bases:
        if not any(f"{base}{suffix}" in present for suffix in ALL_STAT_SUFFIXES):
            missing_signatures.append(base)

    # Full accounting against the expected set
    missing_columns = sorted(expected_set - present)
    extra_columns = [c for c in df.columns if c not in expected_set]

    if strict:
        valid = len(missing_columns) == 0 and len(extra_columns) == 0
    else:
        valid = len(missing_signatures) == 0

    return {
        "valid": valid,
        "missing_metadata": missing_metadata,
        "missing_signatures": missing_signatures,
        "missing_columns": missing_columns,
        "extra_columns": extra_columns,
        "n_expected": len(expected),
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


def preprocess_daily_data(
    gage_flow: pd.DataFrame,
    config: Optional[Dict] = None,
) -> Dict:
    """
    Preprocess daily streamflow data with standardized NA handling.

    For each water year, this function:
      1. Normalizes to a complete daily grid (Oct 1 - Sep 30), filling missing
         dates with NaN.
      2. Computes raw diagnostics: NA count, max consecutive NA run, seasonal
         completeness (from raw data), negative-flow flag, constant-SD flag.
      3. Rejects years that exceed raw-NA limits, have gaps longer than
         max_gap_days, or contain negative Q values.
      4. Interpolates internal gaps <= max_gap_days using linear interpolation
         (internal-only, i.e. leading/trailing NAs are not filled).
      5. After interpolation, checks for residual NAs and optionally rejects
         years that still contain them.
      6. Applies the same interpolation logic to PPT if present.

    Parameters
    ----------
    gage_flow : pd.DataFrame
        Daily data with required columns: date, Q, water_year, month, dowy.
        Optional column: PPT (precipitation).
    config : dict, optional
        Override configuration parameters. Keys (all optional):
            - max_gap_days: int (default from NA_MAX_GAP_DAYS)
            - max_raw_na_per_year: int (default from NA_MAX_RAW_NA_PER_YEAR)
            - reject_negative_flow: bool (default from NA_REJECT_NEGATIVE_FLOW)
            - reject_residual_na: bool (default from NA_REJECT_RESIDUAL_NA)
            - constant_sd_enabled: bool (default from NA_CONSTANT_SD_ENABLED)
            - constant_sd_min_days: int (default from NA_CONSTANT_SD_MIN_DAYS)
            - constant_sd_max_unique: int (default from NA_CONSTANT_SD_MAX_UNIQUE)
            - seasonal_min_fraction: float (default from NA_SEASONAL_MIN_FRACTION)
            - seasonal_definitions: dict (default from NA_SEASONAL_DEFINITIONS)
            - max_raw_na_ppt: int (default from NA_MAX_RAW_NA_PPT)
            - max_gap_ppt: int (default from NA_MAX_GAP_PPT)
            - reject_negative_ppt: bool (default from NA_REJECT_NEGATIVE_PPT)

    Returns
    -------
    dict
        Keys:
            - data: pd.DataFrame with cleaned/interpolated daily data
            - valid_years: list of int, water years that passed all Q checks
            - valid_climate_years: list of int, water years that also passed PPT checks
            - rejected_years: list of dict, each with keys (water_year, reason)
            - seasonal_flags: list of dict, each with keys
              (water_year, winter_complete, spring_complete, summer_complete,
              fall_complete)
            - diagnostics: list of dict, each with keys
              (water_year, raw_na_count, max_na_run, negative_flag,
              constant_sd_flag, post_interp_na_count)
    """
    if config is None:
        config = {}

    max_gap_days = config.get("max_gap_days", NA_MAX_GAP_DAYS)
    max_raw_na = config.get("max_raw_na_per_year", NA_MAX_RAW_NA_PER_YEAR)
    reject_negative = config.get("reject_negative_flow", NA_REJECT_NEGATIVE_FLOW)
    reject_residual = config.get("reject_residual_na", NA_REJECT_RESIDUAL_NA)
    const_sd_enabled = config.get("constant_sd_enabled", NA_CONSTANT_SD_ENABLED)
    const_sd_min_days = config.get("constant_sd_min_days", NA_CONSTANT_SD_MIN_DAYS)
    const_sd_max_unique = config.get("constant_sd_max_unique", NA_CONSTANT_SD_MAX_UNIQUE)
    seasonal_min_frac = config.get("seasonal_min_fraction", NA_SEASONAL_MIN_FRACTION)
    seasonal_defs = config.get("seasonal_definitions", NA_SEASONAL_DEFINITIONS)
    max_raw_na_swe = config.get("max_raw_na_swe", NA_MAX_RAW_NA_SWE)
    max_gap_swe = config.get("max_gap_swe", NA_MAX_GAP_SWE)
    reject_neg_swe = config.get("reject_negative_swe", NA_REJECT_NEGATIVE_SWE)
    max_raw_na_ppt = config.get("max_raw_na_ppt", NA_MAX_RAW_NA_PPT)
    max_gap_ppt = config.get("max_gap_ppt", NA_MAX_GAP_PPT)
    reject_neg_ppt = config.get("reject_negative_ppt", NA_REJECT_NEGATIVE_PPT)

    has_ppt = "PPT" in gage_flow.columns
    has_swe = "SWE" in gage_flow.columns

    valid_years: List[int] = []
    valid_climate_years: List[int] = []
    valid_swe_years: List[int] = []
    rejected_years: List[Dict] = []
    seasonal_flags: List[Dict] = []
    diagnostics: List[Dict] = []
    cleaned_frames: List[pd.DataFrame] = []

    water_years = sorted(gage_flow["water_year"].unique())

    for wy in water_years:
        wy = int(wy)
        # --- (a) Normalize to complete daily grid ---
        # Water year wy spans Oct 1 (wy-1) to Sep 30 (wy)
        wy_start = pd.Timestamp(year=wy - 1, month=10, day=1)
        wy_end = pd.Timestamp(year=wy, month=9, day=30)
        full_dates = pd.date_range(start=wy_start, end=wy_end, freq="D")

        yr_data = gage_flow[gage_flow["water_year"] == wy].copy()

        # Ensure date is datetime for merge
        if "date" in yr_data.columns and not pd.api.types.is_datetime64_any_dtype(yr_data["date"]):
            yr_data["date"] = pd.to_datetime(yr_data["date"])

        # Build complete grid and merge
        grid = pd.DataFrame({"date": full_dates})
        grid["water_year"] = wy
        grid["month"] = grid["date"].dt.month
        grid["dowy"] = (grid["date"] - wy_start).dt.days + 1

        yr_data = grid.merge(yr_data.drop(columns=["water_year", "month", "dowy"], errors="ignore"),
                             on="date", how="left")

        n_days = len(yr_data)

        # --- (b) Raw diagnostics ---
        q_raw = yr_data["Q"].values.copy()
        raw_na_count = int(np.sum(np.isnan(q_raw.astype(float))))

        # Max consecutive NA run
        max_na_run = _max_consecutive_na(yr_data["Q"])

        # Seasonal completeness (from raw observations)
        season_flags_yr = {"water_year": wy}
        for season_name, season_months in seasonal_defs.items():
            season_mask = yr_data["month"].isin(season_months)
            season_q = yr_data.loc[season_mask, "Q"]
            n_season = len(season_q)
            n_valid = int(season_q.notna().sum())
            frac = n_valid / n_season if n_season > 0 else 0.0
            season_flags_yr[f"{season_name}_complete"] = frac >= seasonal_min_frac
        seasonal_flags.append(season_flags_yr)

        # Negative flow flag
        q_finite = q_raw[~np.isnan(q_raw.astype(float))]
        negative_flag = bool(np.any(q_finite < 0)) if len(q_finite) > 0 else False

        # Constant-SD flag: check if any month has suspiciously constant values
        constant_sd_flag = False
        if const_sd_enabled:
            for m in range(1, 13):
                m_mask = yr_data["month"] == m
                m_q = yr_data.loc[m_mask & (yr_data["Q"] > 0), "Q"].dropna()
                if len(m_q) >= const_sd_min_days:
                    n_unique = m_q.nunique()
                    if n_unique <= const_sd_max_unique:
                        constant_sd_flag = True
                        break

        # NA-cause tracking (ice-affected days from USGS qualifier flags)
        na_cause_ice = 0
        na_cause_other = 0
        if raw_na_count > 0 and "flag" in yr_data.columns:
            na_rows = yr_data.loc[yr_data["Q"].isna()]
            if len(na_rows) > 0 and "flag" in na_rows.columns:
                na_flags = na_rows["flag"].fillna("")
                na_cause_ice = int(na_flags.str.contains("Ice|ice", case=False, regex=True).sum())
                na_cause_other = raw_na_count - na_cause_ice

        diag = {
            "water_year": wy,
            "raw_na_count": raw_na_count,
            "max_na_run": max_na_run,
            "negative_flag": negative_flag,
            "constant_sd_flag": constant_sd_flag,
            "post_interp_na_count": raw_na_count,  # updated after interp
            "na_cause_ice": na_cause_ice,
            "na_cause_other": na_cause_other,
        }

        # --- (c) Reject years ---
        rejected = False

        if raw_na_count > max_raw_na:
            rejected_years.append({"water_year": wy, "reason": f"raw_na_count={raw_na_count} > {max_raw_na}"})
            rejected = True

        if not rejected and max_na_run > max_gap_days:
            rejected_years.append({"water_year": wy, "reason": f"max_na_run={max_na_run} > {max_gap_days}"})
            rejected = True

        if not rejected and reject_negative and negative_flag:
            rejected_years.append({"water_year": wy, "reason": "negative_flow"})
            rejected = True

        if rejected:
            diagnostics.append(diag)
            continue

        # --- (d) Interpolate internal gaps <= max_gap_days ---
        yr_data["Q"] = yr_data["Q"].interpolate(
            method="linear", limit=max_gap_days, limit_area="inside"
        )

        # --- (e) Post-interpolation residual NA check ---
        post_na_count = int(yr_data["Q"].isna().sum())
        diag["post_interp_na_count"] = post_na_count

        if reject_residual and post_na_count > 0:
            rejected_years.append({"water_year": wy, "reason": f"residual_na={post_na_count}"})
            diagnostics.append(diag)
            continue

        diagnostics.append(diag)
        valid_years.append(wy)

        # --- (f) PPT handling ---
        ppt_valid = False
        if has_ppt:
            ppt_raw = yr_data["PPT"].values.copy()
            ppt_raw_na = int(np.sum(np.isnan(ppt_raw.astype(float))))
            ppt_max_run = _max_consecutive_na(yr_data["PPT"])
            ppt_finite = ppt_raw[~np.isnan(ppt_raw.astype(float))]
            ppt_negative = bool(np.any(ppt_finite < 0)) if len(ppt_finite) > 0 else False

            ppt_reject = False
            if ppt_raw_na > max_raw_na_ppt:
                ppt_reject = True
            if ppt_max_run > max_gap_ppt:
                ppt_reject = True
            if reject_neg_ppt and ppt_negative:
                ppt_reject = True

            if not ppt_reject:
                yr_data["PPT"] = yr_data["PPT"].interpolate(
                    method="linear", limit=max_gap_ppt, limit_area="inside"
                )
                ppt_post_na = int(yr_data["PPT"].isna().sum())
                if reject_residual and ppt_post_na > 0:
                    ppt_reject = True

            if not ppt_reject:
                ppt_valid = True

        if ppt_valid:
            valid_climate_years.append(wy)

        # --- (f2) SWE handling (same rules as PPT, tracked separately) ---
        swe_valid = False
        if has_swe:
            swe_raw = yr_data["SWE"].values.copy()
            swe_raw_na = int(np.sum(np.isnan(swe_raw.astype(float))))
            swe_max_run = _max_consecutive_na(yr_data["SWE"])
            swe_finite = swe_raw[~np.isnan(swe_raw.astype(float))]
            swe_negative = bool(np.any(swe_finite < 0)) if len(swe_finite) > 0 else False

            swe_reject = False
            if swe_raw_na > max_raw_na_swe:
                swe_reject = True
            if swe_max_run > max_gap_swe:
                swe_reject = True
            if reject_neg_swe and swe_negative:
                swe_reject = True

            if not swe_reject:
                yr_data["SWE"] = yr_data["SWE"].interpolate(
                    method="linear", limit=max_gap_swe, limit_area="inside"
                )
                swe_post_na = int(yr_data["SWE"].isna().sum())
                if reject_residual and swe_post_na > 0:
                    swe_reject = True

            if not swe_reject:
                swe_valid = True

        if swe_valid:
            valid_swe_years.append(wy)

        cleaned_frames.append(yr_data)

    # --- (g) Assemble output ---
    if len(cleaned_frames) > 0:
        cleaned_data = pd.concat(cleaned_frames, ignore_index=True)
    else:
        cleaned_data = pd.DataFrame(columns=gage_flow.columns)

    return {
        "data": cleaned_data,
        "valid_years": valid_years,
        "valid_climate_years": valid_climate_years,
        "valid_swe_years": valid_swe_years,
        "rejected_years": rejected_years,
        "seasonal_flags": seasonal_flags,
        "diagnostics": diagnostics,
    }


def _max_consecutive_na(series: pd.Series) -> int:
    """Compute the length of the longest consecutive NaN run in a Series."""
    is_na = series.isna().astype(int)
    if is_na.sum() == 0:
        return 0
    # Group consecutive NAs: when value is not-NA, cumsum increments,
    # so NAs between non-NAs share the same group label.
    groups = (~series.isna()).astype(int).cumsum()
    # Within each group, cumsum the is_na flag to get run lengths
    run_lengths = is_na.groupby(groups).cumsum()
    return int(run_lengths.max())
