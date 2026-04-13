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

# =============================================================================
# Timing Parameters
# =============================================================================
D_PERCENTILES: List[int] = _config["timing"]["d_percentiles"]

# =============================================================================
# Elasticity Parameters
# =============================================================================
ELASTICITY_WINDOW_YEARS: int = _config["elasticity"]["window_years"]
ELASTICITY_MIN_YEARS: int = _config["elasticity"]["min_years"]
ELASTICITY_MIN_ANNUAL_PPT: float = _config["elasticity"]["min_annual_ppt"]

# =============================================================================
# Q-P Seasonality Parameters
# =============================================================================
QP_SLOPE_WINDOW_DAYS: int = _config["qp_seasonality"]["slope_window_days"]
QP_MIN_YEARS: int = _config["qp_seasonality"]["min_years"]

# =============================================================================
# Runoff Ratio Parameters
# =============================================================================
RUNOFF_MIN_ANNUAL_PPT: float = _config["runoff_ratios"]["min_annual_ppt"]
RUNOFF_MIN_SEASONAL_PPT: float = _config["runoff_ratios"]["min_seasonal_ppt"]

# =============================================================================
# Flow Volumes Parameters
# =============================================================================
FLOW_PERCENTILES: List[int] = _config["flow_volumes"]["percentiles"]

# =============================================================================
# Storage Parameters
# =============================================================================
STORAGE_MIN_YEARS: int = _config["storage"]["min_years"]

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

# ==============================================================================
# NA HANDLING
# ==============================================================================
_na_handling = _config.get("na_handling", {})
_na_interp = _na_handling.get("interpolation", {})
_na_reject = _na_handling.get("year_rejection", {})
_na_const_sd = _na_handling.get("constant_sd_flag", {})
_na_trend = _na_handling.get("trend_completeness", {})
_na_seasonal = _na_handling.get("seasonal_completeness", {})
_na_climate = _na_handling.get("climate_na_policy", {})

NA_MAX_GAP_DAYS = int(_na_interp.get("max_gap_days", 3))
NA_INTERPOLATION_METHOD = _na_interp.get("method", "linear")
NA_INTERNAL_ONLY = _na_interp.get("internal_only", True)
NA_MAX_RAW_NA_PER_YEAR = int(_na_reject.get("max_raw_na_per_year", 30))
NA_REJECT_NEGATIVE_FLOW = _na_reject.get("reject_negative_flow", False)
NA_REJECT_RESIDUAL_NA = _na_reject.get("reject_residual_na", True)
NA_CONSTANT_SD_ENABLED = _na_const_sd.get("enabled", True)
NA_CONSTANT_SD_MIN_DAYS = int(_na_const_sd.get("min_nonzero_days_per_month", 15))
NA_CONSTANT_SD_MAX_UNIQUE = int(_na_const_sd.get("max_unique_values", 1))
NA_TREND_MIN_FRACTION = float(_na_trend.get("min_fraction", 0.80))
NA_DECADE_MIN_FRACTION = float(_na_trend.get("decade_min_fraction", 0.80))
NA_SEASONAL_MIN_FRACTION = float(_na_seasonal.get("min_fraction", 0.80))
NA_SEASONAL_USE_RAW = _na_seasonal.get("use_raw_observations", True)
NA_SEASONAL_DEFINITIONS = _na_seasonal.get("season_definitions", {
    "winter": [12, 1, 2], "spring": [3, 4, 5],
    "summer": [6, 7, 8], "fall": [9, 10, 11]
})
NA_MAX_RAW_NA_PPT = int(_na_climate.get("max_raw_na_per_year_ppt", 30))
NA_MAX_GAP_PPT = int(_na_climate.get("max_interpolation_gap_ppt", 3))
NA_REJECT_NEGATIVE_PPT = _na_climate.get("reject_negative_ppt", True)

# Legacy filtering flag
USE_LEGACY_FILTERING = _config.get("filtering", {}).get("use_legacy_filtering", True)


def get_config() -> Dict[str, Any]:
    """Return the full configuration dictionary."""
    return _config.copy()
