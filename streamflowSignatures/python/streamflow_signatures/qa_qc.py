"""
QA/QC flag computation for streamflow signatures.

Computes quality control flags based on expected ranges and consistency checks.
"""

from typing import Dict
import numpy as np
import pandas as pd

from .config import (
    QAQC_QANN_RANGE,
    QAQC_BFI_RANGE,
    QAQC_FLASHINESS_RANGE,
    QAQC_TQMEAN_RANGE,
    QAQC_D50_RANGE,
    QAQC_ELASTICITY_RANGE,
    QAQC_RUNOFF_RATIO_RANGE,
    QAQC_SEASONAL_SUM_TOLERANCE,
    QAQC_MAX_NA_FRACTION,
    QAQC_HIGH_NA_SUFFIXES,
    QAQC_HIGH_NA_SCALARS,
)


def compute_qa_flags(signatures: pd.DataFrame) -> pd.DataFrame:
    """
    Compute QA/QC flags for a DataFrame of streamflow signatures.

    Parameters
    ----------
    signatures : pd.DataFrame
        DataFrame with signature columns (one row per gage)

    Returns
    -------
    pd.DataFrame
        DataFrame with QA/QC flag columns added

    Notes
    -----
    Flags are computed based on expected ranges from config.
    A flag value of True indicates a potential quality issue.
    """
    df = signatures.copy()

    # Initialize all flag columns to False
    flag_columns = [
        "flagged_for_qann_range",
        "flagged_for_bfi_eckhardt_range",
        "flagged_for_bfi_lynehollick_range",
        "flagged_for_flashiness_range",
        "flagged_for_tqmean_range",
        "flagged_for_d50_range",
        "flagged_for_elasticity_range",
        "flagged_for_runoff_ratio_range",
        "flagged_for_seasonal_sum",
        "flagged_for_percentile_order",
        "flagged_for_timing_order",
        "flagged_for_high_na",
    ]

    for col in flag_columns:
        df[col] = False

    # 1. Qann_mean outside 0-2000 mm
    if "Qann_mean" in df.columns:
        df["flagged_for_qann_range"] = (
            (df["Qann_mean"] < QAQC_QANN_RANGE[0]) |
            (df["Qann_mean"] > QAQC_QANN_RANGE[1])
        ).fillna(False)

    # 2. BFI_Eckhardt outside 0-1
    if "BFI_Eckhardt_mean" in df.columns:
        df["flagged_for_bfi_eckhardt_range"] = (
            (df["BFI_Eckhardt_mean"] < QAQC_BFI_RANGE[0]) |
            (df["BFI_Eckhardt_mean"] > QAQC_BFI_RANGE[1])
        ).fillna(False)

    # 3. BFI_LyneHollick outside 0-1
    if "BFI_LyneHollick_mean" in df.columns:
        df["flagged_for_bfi_lynehollick_range"] = (
            (df["BFI_LyneHollick_mean"] < QAQC_BFI_RANGE[0]) |
            (df["BFI_LyneHollick_mean"] > QAQC_BFI_RANGE[1])
        ).fillna(False)

    # 4. flashinessRB outside 0-2
    if "flashinessRB_mean" in df.columns:
        df["flagged_for_flashiness_range"] = (
            (df["flashinessRB_mean"] < QAQC_FLASHINESS_RANGE[0]) |
            (df["flashinessRB_mean"] > QAQC_FLASHINESS_RANGE[1])
        ).fillna(False)

    # 5. TQmean outside 0-100%
    if "TQmean_mean" in df.columns:
        df["flagged_for_tqmean_range"] = (
            (df["TQmean_mean"] < QAQC_TQMEAN_RANGE[0]) |
            (df["TQmean_mean"] > QAQC_TQMEAN_RANGE[1])
        ).fillna(False)

    # 6. D50_day outside 1-366
    if "D50_day_mean" in df.columns:
        df["flagged_for_d50_range"] = (
            (df["D50_day_mean"] < QAQC_D50_RANGE[0]) |
            (df["D50_day_mean"] > QAQC_D50_RANGE[1])
        ).fillna(False)

    # 7. elasticity_static outside 0.1-5
    if "elasticity_static" in df.columns:
        df["flagged_for_elasticity_range"] = (
            (df["elasticity_static"] < QAQC_ELASTICITY_RANGE[0]) |
            (df["elasticity_static"] > QAQC_ELASTICITY_RANGE[1])
        ).fillna(False)

    # 8. annual_runoff_ratio outside 0.01-1.5
    if "annual_runoff_ratio_mean" in df.columns:
        df["flagged_for_runoff_ratio_range"] = (
            (df["annual_runoff_ratio_mean"] < QAQC_RUNOFF_RATIO_RANGE[0]) |
            (df["annual_runoff_ratio_mean"] > QAQC_RUNOFF_RATIO_RANGE[1])
        ).fillna(False)

    # 9. Seasonal sum deviation > tolerance
    # |sum of seasonal / annual - 1| > tolerance
    seasonal_cols = ["Qwin_mean", "Qspr_mean", "Qsum_mean", "Qfal_mean"]
    if all(col in df.columns for col in seasonal_cols) and "Qann_mean" in df.columns:
        seasonal_sum = (
            df["Qwin_mean"] + df["Qspr_mean"] + df["Qsum_mean"] + df["Qfal_mean"]
        )
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = seasonal_sum / df["Qann_mean"]
            deviation = np.abs(ratio - 1)
        df["flagged_for_seasonal_sum"] = (
            deviation > QAQC_SEASONAL_SUM_TOLERANCE
        ).fillna(False)

    # 10. Percentile order violation (Q5 < Q25 < Q50 < Q75 < Q95)
    pct_cols = ["Q5_mean", "Q25_mean", "Q50_mean", "Q75_mean", "Q95_mean"]
    if all(col in df.columns for col in pct_cols):
        df["flagged_for_percentile_order"] = (
            (df["Q5_mean"] > df["Q25_mean"]) |
            (df["Q25_mean"] > df["Q50_mean"]) |
            (df["Q50_mean"] > df["Q75_mean"]) |
            (df["Q75_mean"] > df["Q95_mean"])
        ).fillna(False)

    # 11. Timing order violation (D5 < D50 < D95)
    timing_cols = ["D5_day_mean", "D50_day_mean", "D95_day_mean"]
    if all(col in df.columns for col in timing_cols):
        df["flagged_for_timing_order"] = (
            (df["D5_day_mean"] > df["D50_day_mean"]) |
            (df["D50_day_mean"] > df["D95_day_mean"])
        ).fillna(False)

    # 12. High NA fraction: > max_na_fraction of the gage's SIGNATURE columns are NA.
    # Denominator = the signature columns PRESENT in the table (16 statistic suffixes
    # + the per-gage scalar registry, config qa_qc.high_na_denominator — one manifest
    # shared with julia/src/qa_qc.jl and rpkg/R/qa_qc.R). Metadata, unregistered
    # diagnostics and flag columns never enter it; families NA-filled by the table
    # union DO count. Before 2026-09-04 this counted every non-flag column (incl.
    # string metadata) and differed from Julia (metadata only) and rpkg (mean/median).
    signature_cols = high_na_denominator_columns(df.columns)
    if signature_cols:
        na_fraction = df[signature_cols].isna().mean(axis=1)
        df["flagged_for_high_na"] = (na_fraction > QAQC_MAX_NA_FRACTION)

    return df


def high_na_denominator_columns(columns) -> list:
    """
    Columns that enter the ``flagged_for_high_na`` denominator: every name ending in
    one of the statistic suffixes (``QAQC_HIGH_NA_SUFFIXES``) plus the per-gage scalar
    registry (``QAQC_HIGH_NA_SCALARS``), in input order; ``flagged_*`` columns excluded.
    Shared definition with ``julia/src/qa_qc.jl`` and ``rpkg/R/qa_qc.R``.
    """
    scalars = set(QAQC_HIGH_NA_SCALARS)
    out = []
    for c in columns:
        name = str(c)
        if name.startswith("flagged_"):
            continue
        if name in scalars or any(name.endswith(sfx) for sfx in QAQC_HIGH_NA_SUFFIXES):
            out.append(name)
    return out


def get_flag_columns() -> list:
    """Return list of QA/QC flag column names."""
    return [
        "flagged_for_qann_range",
        "flagged_for_bfi_eckhardt_range",
        "flagged_for_bfi_lynehollick_range",
        "flagged_for_flashiness_range",
        "flagged_for_tqmean_range",
        "flagged_for_d50_range",
        "flagged_for_elasticity_range",
        "flagged_for_runoff_ratio_range",
        "flagged_for_seasonal_sum",
        "flagged_for_percentile_order",
        "flagged_for_timing_order",
        "flagged_for_high_na",
    ]
