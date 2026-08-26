"""
Centralized configuration for streamflow signature calculations.

Loads configuration from the shared JSON config file and exposes
parameters as module-level constants for use across all signature modules.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional

# Config resolution order (2026-08-25, Codex Phase-1 review MAJOR):
#   1. STREAMFLOW_CONFIG env var (explicit override — the same variable the
#      Julia package honors and the benchmark runners record in provenance)
#   2. the repo's canonical config/signatures_config.json (source tree /
#      editable install — the source of truth during development)
#   3. the copy bundled inside the package (installed-wheel fallback; kept in
#      sync with the canonical file, guarded by a unit test)
_REPO_CONFIG = Path(__file__).parent.parent.parent / "config" / "signatures_config.json"
_PACKAGED_CONFIG = Path(__file__).parent / "data" / "signatures_config.json"


def _resolve_config_file() -> Path:
    env = os.environ.get("STREAMFLOW_CONFIG", "")
    if env:
        p = Path(env)
        if not p.is_file():
            raise FileNotFoundError(f"STREAMFLOW_CONFIG points at a missing file: {p}")
        return p
    if _REPO_CONFIG.is_file():
        return _REPO_CONFIG
    if _PACKAGED_CONFIG.is_file():
        return _PACKAGED_CONFIG
    raise FileNotFoundError(
        "signatures_config.json not found: set STREAMFLOW_CONFIG, or run from the "
        f"repo (looked at {_REPO_CONFIG}), or reinstall the package "
        f"(bundled copy missing at {_PACKAGED_CONFIG})")


_CONFIG_FILE = _resolve_config_file()


def load_config() -> Dict[str, Any]:
    """Load configuration from the resolved JSON file (see resolution order above)."""
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
_na_snow = _na_handling.get("snow_na_policy", {})
NA_MAX_RAW_NA_SWE = int(_na_snow.get("max_raw_na_per_year_swe", 30))
NA_MAX_GAP_SWE = int(_na_snow.get("max_interpolation_gap_swe", 3))
NA_REJECT_NEGATIVE_SWE = _na_snow.get("reject_negative_swe", True)
NA_MAX_RAW_NA_PPT = int(_na_climate.get("max_raw_na_per_year_ppt", 30))
NA_MAX_GAP_PPT = int(_na_climate.get("max_interpolation_gap_ppt", 3))
NA_REJECT_NEGATIVE_PPT = _na_climate.get("reject_negative_ppt", True)

# Legacy filtering flag
USE_LEGACY_FILTERING = _config.get("filtering", {}).get("use_legacy_filtering", True)

# ==============================================================================
# STATS FLOOR — None when the section is absent (no floor; legacy behavior).
# Exemptions (recession, elasticity) are enforced at the orchestration layer.
# ==============================================================================
_stats_floor = _config.get("stats_floor", {})
_mv = _stats_floor.get("min_values_for_stats", None)
MIN_VALUES_FOR_STATS: Optional[int] = int(_mv) if _mv is not None else None

# ==============================================================================
# DROUGHT (Adelsperger et al.) — fallbacks mirror julia/src/config.jl.
# An ABSENT section disables the family entirely (no drought columns).
# ==============================================================================
_drought_config = _config.get("drought", {})
DROUGHT_ENABLED: bool = bool(_drought_config)
DROUGHT_SMOOTH_WINDOW: int = int(_drought_config.get("smoothing_window_days", 7))
DROUGHT_SMOOTH_ALIGNMENT: str = str(_drought_config.get("smoothing_alignment", "center"))
DROUGHT_SMOOTH_MIN_VALID: int = int(_drought_config.get("smoothing_min_valid_days", 4))
DROUGHT_THRESHOLD_METHOD: str = str(_drought_config.get("threshold_method", "fixed"))
DROUGHT_PERCENTILES: List[int] = [int(p) for p in
                                  _drought_config.get("threshold_percentiles", [2, 5, 10, 20, 30])]
DROUGHT_PLOTTING_POSITION: str = str(_drought_config.get("plotting_position", "weibull"))
DROUGHT_BELOW_RANGE_POLICY: str = str(_drought_config.get("below_plotting_range_policy", "na"))
DROUGHT_MIN_YEARS: int = int(_drought_config.get("min_years_for_threshold", 10))

if DROUGHT_ENABLED:
    # Fail fast on unsupported/invalid configuration (mirrors julia/src/config.jl)
    if DROUGHT_THRESHOLD_METHOD != "fixed":
        raise ValueError(f"drought.threshold_method must be 'fixed', got "
                         f"{DROUGHT_THRESHOLD_METHOD!r} (the variable day-of-year "
                         f"method is intentionally not implemented)")
    if DROUGHT_PLOTTING_POSITION != "weibull":
        raise ValueError(f"drought.plotting_position must be 'weibull', got "
                         f"{DROUGHT_PLOTTING_POSITION!r}")
    if DROUGHT_BELOW_RANGE_POLICY not in ("na", "clamp"):
        raise ValueError(f"drought.below_plotting_range_policy must be 'na' or "
                         f"'clamp', got {DROUGHT_BELOW_RANGE_POLICY!r}")
    if DROUGHT_SMOOTH_WINDOW < 1:
        raise ValueError("drought.smoothing_window_days must be >= 1")
    if DROUGHT_SMOOTH_ALIGNMENT == "center" and DROUGHT_SMOOTH_WINDOW % 2 == 0:
        raise ValueError("drought.smoothing_window_days must be ODD for centered alignment")
    if DROUGHT_SMOOTH_MIN_VALID > DROUGHT_SMOOTH_WINDOW:
        raise ValueError("drought.smoothing_min_valid_days must be <= smoothing_window_days")
    if DROUGHT_SMOOTH_ALIGNMENT not in ("center", "trailing"):
        raise ValueError(f"drought.smoothing_alignment must be 'center' or "
                         f"'trailing', got {DROUGHT_SMOOTH_ALIGNMENT!r}")
    if DROUGHT_SMOOTH_MIN_VALID < 1:
        raise ValueError("drought.smoothing_min_valid_days must be >= 1")
    if not DROUGHT_PERCENTILES or any(not (0 < p < 100) for p in DROUGHT_PERCENTILES):
        raise ValueError("drought.threshold_percentiles must be non-empty and all in (0, 100)")
    if list(DROUGHT_PERCENTILES) != sorted(DROUGHT_PERCENTILES):
        raise ValueError("drought.threshold_percentiles must be sorted ascending")
    if len(set(DROUGHT_PERCENTILES)) != len(DROUGHT_PERCENTILES):
        raise ValueError("drought.threshold_percentiles must be unique")
    if DROUGHT_MIN_YEARS < 1:
        raise ValueError("drought.min_years_for_threshold must be >= 1")

# ==============================================================================
# SNOW — fallbacks mirror julia/src/config.jl (absent section => gate disabled)
# ==============================================================================
_snow_config = _config.get("snow", {})
SNOW_SWE_THRESHOLD_MM: float = float(_snow_config.get("swe_day_threshold_mm", 10.0))
SNOW_SEASONAL_MIN_DAYS: int = int(_snow_config.get("seasonal_spell_min_days", 60))
SNOW_MELT_COM_FRACTION: float = float(_snow_config.get("melt_com_fraction", 0.5))
SNOW_MIN_ANNUAL_PPT_MM: float = float(_snow_config.get("min_annual_ppt_mm", 10))
# Record-anchored decade gate: threshold-dependent snow metrics need
# >= decade_completeness of the SWE-valid years in the record's first AND last
# decade to be computable before trend stats are computed. Threshold is the SAME
# decade_min_fraction knob the streamflow gate uses (linked by design).
SNOW_RECORD_DECADE_GATE: bool = bool(_snow_config.get("record_anchored_decade_gate", False))

# ==============================================================================
# ANNUAL VALUES EXPORT — defaults to False when the section is absent
# (matching julia/src/config.jl CFG_SAVE_ANNUAL_VALUES)
# ==============================================================================
SAVE_ANNUAL_VALUES: bool = bool(_config.get("annual_values", {}).get("save", False))

# ==============================================================================
# CHANGEPOINT (Pettitt) — fallbacks mirror julia/src/config.jl
# ==============================================================================
_cp_config = _config.get("changepoint", {})
CHANGEPOINT_ENABLED = bool(_cp_config.get("enabled", False))
CP_START_WATER_YEAR = int(_cp_config.get("start_water_year", 1980))
CP_END_WATER_YEAR = int(_cp_config.get("end_water_year", 2024))
CP_MIN_TOTAL_OBS = int(_cp_config.get("min_total_obs", 20))
CP_MIN_SEGMENT_OBS = int(_cp_config.get("min_segment_obs", 10))


def get_config() -> Dict[str, Any]:
    """Return the full configuration dictionary."""
    return _config.copy()
