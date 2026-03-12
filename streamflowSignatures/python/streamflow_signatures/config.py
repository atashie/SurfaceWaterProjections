"""
Centralized configuration for streamflow signature calculations.

Loads configuration from the shared JSON config file and exposes
parameters as module-level constants for use across all signature modules.
"""

import json
from pathlib import Path
from typing import Dict, Any, List

# Find config file relative to this module
_CONFIG_DIR = Path(__file__).parent.parent.parent / "config"
_CONFIG_FILE = _CONFIG_DIR / "signatures_config.json"


def load_config() -> Dict[str, Any]:
    """Load configuration from JSON file."""
    if not _CONFIG_FILE.exists():
        raise FileNotFoundError(f"Config file not found: {_CONFIG_FILE}")

    with open(_CONFIG_FILE, "r") as f:
        return json.load(f)


# Load config at module import time
_config = load_config()

# =============================================================================
# Filtering Parameters
# =============================================================================
MIN_NUM_YEARS: int = _config["filtering"]["min_num_years"]
MIN_FRAC_GOOD_DATA: float = _config["filtering"]["min_frac_good_data"]
MIN_Q_VALUE: float = _config["filtering"]["min_q_value"]
MIN_DAYS_ABOVE_THRESHOLD: int = _config["filtering"]["min_days_above_threshold"]
MIN_NONA_DAYS_ANNUAL: int = _config["filtering"]["min_nona_days_annual"]

# =============================================================================
# Water Year
# =============================================================================
WATER_YEAR_START_MONTH: int = _config["water_year"]["start_month"]

# =============================================================================
# Baseflow Parameters
# =============================================================================
ECKHARDT_BFIMAX: float = _config["baseflow"]["eckhardt_bfimax"]
ECKHARDT_ALPHA: float = _config["baseflow"]["eckhardt_alpha"]
LYNE_HOLLICK_ALPHA: float = _config["baseflow"]["lyne_hollick_alpha"]
LYNE_HOLLICK_PASSES: int = _config["baseflow"]["lyne_hollick_passes"]
BASEFLOW_MIN_DAYS: int = _config["baseflow"]["min_days"]
BASEFLOW_MAX_MISSING_FRAC: float = _config["baseflow"]["max_missing_frac"]

# =============================================================================
# Recession Parameters
# =============================================================================
RECESSION_MIN_LENGTH: int = _config["recession"]["min_length"]
RECESSION_MIN_EVENTS: int = _config["recession"]["min_events"]

# =============================================================================
# Pulse Parameters
# =============================================================================
HIGH_PULSE_PERCENTILE: float = _config["pulses"]["high_percentile"]
LOW_PULSE_PERCENTILE: float = _config["pulses"]["low_percentile"]
FLOW_REVERSAL_THRESHOLD: float = _config["pulses"]["flow_reversal_threshold"]
PULSES_MIN_DAYS: int = _config["pulses"]["min_days"]

# =============================================================================
# Timing Parameters
# =============================================================================
TIMING_MIN_DAYS: int = _config["timing"]["min_days"]
D_PERCENTILES: List[int] = _config["timing"]["d_percentiles"]

# =============================================================================
# Elasticity Parameters
# =============================================================================
ELASTICITY_WINDOW_YEARS: int = _config["elasticity"]["window_years"]
ELASTICITY_MIN_YEARS: int = _config["elasticity"]["min_years"]
ELASTICITY_MIN_DATA_COMPLETENESS: float = _config["elasticity"]["min_data_completeness"]
ELASTICITY_MIN_ANNUAL_PPT: float = _config["elasticity"]["min_annual_ppt"]

# =============================================================================
# Q-P Seasonality Parameters
# =============================================================================
QP_SLOPE_WINDOW_DAYS: int = _config["qp_seasonality"]["slope_window_days"]
QP_MIN_YEARS: int = _config["qp_seasonality"]["min_years"]
QP_MIN_DAYS_PER_YEAR: int = _config["qp_seasonality"]["min_days_per_year"]
QP_MAX_NA_FRAC: float = _config["qp_seasonality"]["max_na_frac"]

# =============================================================================
# Runoff Ratio Parameters
# =============================================================================
RUNOFF_MIN_ANNUAL_PPT: float = _config["runoff_ratios"]["min_annual_ppt"]
RUNOFF_MIN_SEASONAL_PPT: float = _config["runoff_ratios"]["min_seasonal_ppt"]

# =============================================================================
# Flow Volumes Parameters
# =============================================================================
FLOW_VOLUMES_MIN_DAYS: int = _config["flow_volumes"]["min_days"]
FLOW_PERCENTILES: List[int] = _config["flow_volumes"]["percentiles"]

# =============================================================================
# FDC Parameters
# =============================================================================
FDC_MIN_DAYS: int = _config["fdc"]["min_days"]

# =============================================================================
# Flashiness Parameters
# =============================================================================
FLASHINESS_MIN_DAYS: int = _config["flashiness"]["min_days"]
FLASHINESS_MAX_MISSING_FRAC: float = _config["flashiness"]["max_missing_frac"]

# =============================================================================
# Storage Parameters
# =============================================================================
STORAGE_MIN_YEARS: int = _config["storage"]["min_years"]
STORAGE_MIN_DAYS_PER_YEAR: int = _config["storage"]["min_days_per_year"]
STORAGE_MAX_NA_FRAC: float = _config["storage"]["max_na_frac"]

# =============================================================================
# QA/QC Parameters
# =============================================================================
QAQC_QANN_RANGE: List[float] = _config["qa_qc"]["qann_range"]
QAQC_BFI_RANGE: List[float] = _config["qa_qc"]["bfi_range"]
QAQC_FLASHINESS_RANGE: List[float] = _config["qa_qc"]["flashiness_range"]
QAQC_TQMEAN_RANGE: List[float] = _config["qa_qc"]["tqmean_range"]
QAQC_D50_RANGE: List[int] = _config["qa_qc"]["d50_range"]
QAQC_ELASTICITY_RANGE: List[float] = _config["qa_qc"]["elasticity_range"]
QAQC_RUNOFF_RATIO_RANGE: List[float] = _config["qa_qc"]["runoff_ratio_range"]
QAQC_SEASONAL_SUM_TOLERANCE: float = _config["qa_qc"]["seasonal_sum_tolerance"]
QAQC_MAX_NA_FRACTION: float = _config["qa_qc"]["max_na_fraction"]


# =============================================================================
# Metadata Parameters
# =============================================================================
INCLUDE_HUMAN_INTERFERENCE: bool = _config["metadata"]["include_human_interference"]
GAGES_II_DIR: str = _config["metadata"]["gages_ii_dir"]
HYDAT_PATH = _config["metadata"]["hydat_path"]
INTERFERENCE_COLUMNS: List[str] = _config["metadata"]["interference_columns"]


def get_config() -> Dict[str, Any]:
    """Return the full configuration dictionary."""
    return _config.copy()
