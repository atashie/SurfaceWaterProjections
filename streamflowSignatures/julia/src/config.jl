"""
Centralized configuration for streamflow signature calculations.

Loads configuration from the shared JSON config file and exports
parameters as constants for use across all signature modules.
"""

using JSON3

# Find config file relative to this module
const _CONFIG_DIR = joinpath(@__DIR__, "..", "..", "config")
const _DEFAULT_CONFIG_FILE = joinpath(_CONFIG_DIR, "signatures_config.json")

# ⚠️ CACHING GOTCHA (verified 2026-07-28): every `CFG_*` constant below is evaluated when
# this module is PRECOMPILED, and Julia does NOT invalidate the precompile cache when an
# ENV var changes. So setting `STREAMFLOW_CONFIG` (or editing the JSON) has NO EFFECT on
# an already-precompiled package — a config-variant experiment can silently run with the
# previous config. Symptom: output identical to the prior run despite a changed config.
# To make a config change take effect the cache must be invalidated by CONTENT — `touch`
# alone does NOT work (verified 2026-07-28). Delete the compiled cache:
#     rm -rf ~/.julia/compiled/v1.12/StreamflowSignatures
# then run, and verify with:
#     julia --project=julia -e 'using StreamflowSignatures; println(CFG_DROUGHT_ENABLED)'
# (`julia --compiled-modules=no` also works but is much slower). The cache then holds the
# VARIANT build, so purge it again afterwards. This is also why the benchmark runner
# reads STREAMFLOW_START_WATER_YEAR / _END_WATER_YEAR / _MIN_QUALIFYING_DATA_FRACTION at
# RUNTIME inside main() rather than from constants here.
const _CONFIG_FILE = get(ENV, "STREAMFLOW_CONFIG", _DEFAULT_CONFIG_FILE)

function load_config()
    if !isfile(_CONFIG_FILE)
        error("Config file not found: $_CONFIG_FILE")
    end
    return JSON3.read(read(_CONFIG_FILE, String))
end

# Load config at module load time
const _config = load_config()

# =============================================================================
# Filtering Parameters
# =============================================================================
const CFG_MIN_NUM_YEARS = _config.filtering.min_num_years
const CFG_MIN_FRAC_GOOD_DATA = _config.filtering.min_frac_good_data
const CFG_MIN_Q_VALUE = _config.filtering.min_q_value
const CFG_MIN_DAYS_ABOVE_THRESHOLD = _config.filtering.min_days_above_threshold
const CFG_MIN_NONA_DAYS_ANNUAL = _config.filtering.min_nona_days_annual

# Experiment controls — JSON defaults with ENV overrides
const CFG_START_WATER_YEAR = let
    json_val = get(get(_config, "filtering", Dict()), "start_water_year", nothing)
    env_val = get(ENV, "STREAMFLOW_START_WATER_YEAR", "")
    v = if env_val != ""
        try parse(Int, env_val) catch e
            error("Invalid STREAMFLOW_START_WATER_YEAR='$env_val': must be an integer")
        end
    else
        json_val
    end
    v === nothing ? nothing : Int(v)
end

const CFG_MIN_QUALIFYING_DATA_FRACTION = let
    json_val = get(get(_config, "filtering", Dict()), "min_qualifying_data_fraction", nothing)
    env_val = get(ENV, "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "")
    v = if env_val != ""
        try parse(Float64, env_val) catch e
            error("Invalid STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION='$env_val': must be a number")
        end
    else
        json_val
    end
    v === nothing ? nothing : Float64(v)
end

# =============================================================================
# Water Year
# =============================================================================
const CFG_WATER_YEAR_START_MONTH = _config.water_year.start_month

# =============================================================================
# Baseflow Parameters
# =============================================================================
const CFG_ECKHARDT_BFIMAX = _config.baseflow.eckhardt_bfimax
const CFG_ECKHARDT_ALPHA = _config.baseflow.eckhardt_alpha
const CFG_LYNE_HOLLICK_ALPHA = _config.baseflow.lyne_hollick_alpha
const CFG_LYNE_HOLLICK_PASSES = _config.baseflow.lyne_hollick_passes

# =============================================================================
# Recession Parameters
# =============================================================================
const CFG_RECESSION_MIN_LENGTH = _config.recession.min_length
const CFG_RECESSION_MIN_EVENTS = _config.recession.min_events

# =============================================================================
# Pulse Parameters
# =============================================================================
const CFG_HIGH_PULSE_PERCENTILE = _config.pulses.high_percentile
const CFG_LOW_PULSE_PERCENTILE = _config.pulses.low_percentile
const CFG_FLOW_REVERSAL_THRESHOLD = _config.pulses.flow_reversal_threshold

# =============================================================================
# Timing Parameters
# =============================================================================
const CFG_D_PERCENTILES = collect(_config.timing.d_percentiles)

# =============================================================================
# Elasticity Parameters
# =============================================================================
const CFG_ELASTICITY_WINDOW_YEARS = _config.elasticity.window_years
const CFG_ELASTICITY_MIN_YEARS = _config.elasticity.min_years
const CFG_ELASTICITY_MIN_ANNUAL_PPT = _config.elasticity.min_annual_ppt

# =============================================================================
# Q-P Seasonality Parameters
# =============================================================================
const CFG_QP_SLOPE_WINDOW_DAYS = _config.qp_seasonality.slope_window_days
const CFG_QP_MIN_YEARS = _config.qp_seasonality.min_years

# =============================================================================
# Runoff Ratio Parameters
# =============================================================================
const CFG_RUNOFF_MIN_ANNUAL_PPT = _config.runoff_ratios.min_annual_ppt
const CFG_RUNOFF_MIN_SEASONAL_PPT = _config.runoff_ratios.min_seasonal_ppt

# =============================================================================
# Flow Volumes Parameters
# =============================================================================
const CFG_FLOW_PERCENTILES = collect(_config.flow_volumes.percentiles)

# =============================================================================
# Storage Parameters
# =============================================================================
const CFG_STORAGE_MIN_YEARS = _config.storage.min_years

# =============================================================================
# Snow Parameters (absent-section defaults match the shipped config)
# =============================================================================
const _snow_config = get(_config, "snow", Dict())
const CFG_SNOW_SWE_THRESHOLD_MM = Float64(get(_snow_config, "swe_day_threshold_mm", 10.0))
const CFG_SNOW_SEASONAL_MIN_DAYS = Int(get(_snow_config, "seasonal_spell_min_days", 60))
const CFG_SNOW_MELT_COM_FRACTION = Float64(get(_snow_config, "melt_com_fraction", 0.5))
const CFG_SNOW_MIN_ANNUAL_PPT_MM = Float64(get(_snow_config, "min_annual_ppt_mm", 10))
# Record-anchored decade gate (July 2026): threshold-dependent snow metrics need
# >= decade_completeness of the SWE-valid years in the record's first AND last
# decade to be computable (snowy) before trend stats are computed. NOTE: shipped
# config = true; absent key => disabled (legacy). Threshold is the SAME
# decade_min_fraction knob the streamflow gate uses (linked by design).
const CFG_SNOW_RECORD_DECADE_GATE = Bool(get(_snow_config, "record_anchored_decade_gate", false))

# =============================================================================
# Climate Parameters (guidelines 3.9, 2026-10-12)
#
# 39 metrics from basin-mean Daymet prcp / tmin / tmax: 20 precipitation bases on
# PPT-valid years, 19 temperature bases on temperature-valid years. Absent section
# => family disabled (no climate columns emitted) — the additivity-proof switch.
# =============================================================================
const _climate_config = get(_config, "climate", Dict())
const CFG_CLIMATE_ENABLED = !isempty(_climate_config)

# Strict scalar readers (review 2026-10-12): a numeric setting must be a finite Real that
# is NOT a Bool (Julia would otherwise accept `true` as 1); an integer setting must also
# be whole; a flag must be an actual Bool (not 0 / 1). Shared with the Python and R ports.
function _cfg_num(section, key::String, default)::Float64
    v = get(section, key, default)
    (v isa Real && !(v isa Bool) && isfinite(v)) || error("config $key: expected a finite number, got $(repr(v))")
    return Float64(v)
end
function _cfg_int(section, key::String, default)::Int
    v = get(section, key, default)
    (v isa Real && !(v isa Bool) && isfinite(v) && isinteger(v)) || error("config $key: expected an integer, got $(repr(v))")
    return Int(v)
end
function _cfg_bool(section, key::String, default)::Bool
    v = get(section, key, default)
    (v isa Bool) || error("config $key: expected true/false, got $(repr(v))")
    return v
end
function _cfg_numvec(section, key::String, default)::Vector{Float64}
    v = collect(get(section, key, default))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x), v) || error("config $key: expected a list of finite numbers, got $(repr(v))")
    return Float64.(v)
end

const CFG_CLIMATE_DRY_DAY_MM = _cfg_num(_climate_config, "dry_day_threshold_mm", 0.1)
const CFG_CLIMATE_WET_THRESHOLDS_MM = _cfg_numvec(_climate_config, "wet_day_thresholds_mm", [10, 25, 50])
const CFG_CLIMATE_MULTIDAY_WINDOW = _cfg_int(_climate_config, "multi_day_max_window_days", 5)
const _climate_cum_fraction_raw = _cfg_num(_climate_config, "cumulative_fraction", 0.5)
const CFG_CLIMATE_FROST_TEMP_C = _cfg_num(_climate_config, "frost_temp_c", 0.0)
const CFG_CLIMATE_FROST_LAST_END_MONTH = _cfg_int(_climate_config, "frost_last_end_month", 6)
const CFG_CLIMATE_GDD_BASE_C = _cfg_num(_climate_config, "gdd_base_c", 0.0)
# Record-anchored decade gate for frost_first_dowy / frost_last_dowy (snow precedent;
# user decision 2026-10-12). Shipped config = true; absent key => disabled.
const CFG_CLIMATE_RECORD_DECADE_GATE = _cfg_bool(_climate_config, "record_anchored_decade_gate", false)

if CFG_CLIMATE_ENABLED
    if CFG_CLIMATE_DRY_DAY_MM < 0
        error("climate.dry_day_threshold_mm must be >= 0, got $(CFG_CLIMATE_DRY_DAY_MM).")
    end
    if isempty(CFG_CLIMATE_WET_THRESHOLDS_MM) || any(x -> x <= 0, CFG_CLIMATE_WET_THRESHOLDS_MM) ||
       !issorted(CFG_CLIMATE_WET_THRESHOLDS_MM) ||
       length(unique(CFG_CLIMATE_WET_THRESHOLDS_MM)) != length(CFG_CLIMATE_WET_THRESHOLDS_MM)
        error("climate.wet_day_thresholds_mm must be a non-empty, ascending, unique list of positive values (column names derive from them), got $(CFG_CLIMATE_WET_THRESHOLDS_MM).")
    end
    # Window 1 would duplicate prcp_1d_max under a second name; > 366 can never be
    # computed on a water year (review finding 6, 2026-10-12)
    if !(2 <= CFG_CLIMATE_MULTIDAY_WINDOW <= 366)
        error("climate.multi_day_max_window_days must lie in 2..366, got $(CFG_CLIMATE_MULTIDAY_WINDOW).")
    end
    if !(0 < _climate_cum_fraction_raw < 1)
        error("climate.cumulative_fraction must lie strictly inside (0, 1), got $(_climate_cum_fraction_raw).")
    end
    # The column name carries the percentage (prcp_<pct>_dowy), so only whole
    # percentages are representable (review finding 5)
    if !isinteger(round(100 * _climate_cum_fraction_raw; digits=9))
        error("climate.cumulative_fraction must be a whole percentage (e.g. 0.5 -> prcp_50_dowy), got $(_climate_cum_fraction_raw).")
    end
end
# The fraction USED is the canonical whole percentage the column name shows (review
# round 2: validating a rounded value but computing with the raw one gave two
# definitions under one name)
const CFG_CLIMATE_CUM_FRACTION = CFG_CLIMATE_ENABLED ? round(100 * _climate_cum_fraction_raw) / 100 : _climate_cum_fraction_raw
if CFG_CLIMATE_ENABLED
    # the canonical percentage must stay strictly inside the open interval (review round 3)
    if !(1 <= round(Int, 100 * CFG_CLIMATE_CUM_FRACTION) <= 99)
        error("climate.cumulative_fraction must canonicalise to a whole percentage between 1 and 99, got $(_climate_cum_fraction_raw).")
    end
    if !(1 <= CFG_CLIMATE_FROST_LAST_END_MONTH <= 12)
        error("climate.frost_last_end_month must be a calendar month 1-12, got $(CFG_CLIMATE_FROST_LAST_END_MONTH).")
    end
end

# =============================================================================
# Drought Parameters (Adelsperger et al., in review)
#
# Fixed (whole-record) percentile thresholds on the 7-day-smoothed flow series.
# Percentiles are MAGNITUDE (non-exceedance) levels — 10 == the flow exceeded 90%
# of the time — consistent with the paper and this repo's Q{n} columns.
# Absent section => family disabled (no drought columns emitted).
# =============================================================================
const _drought_config = get(_config, "drought", Dict())
const CFG_DROUGHT_ENABLED = !isempty(_drought_config)
const CFG_DROUGHT_SMOOTH_WINDOW = Int(get(_drought_config, "smoothing_window_days", 7))
const CFG_DROUGHT_SMOOTH_ALIGNMENT = String(get(_drought_config, "smoothing_alignment", "center"))
const CFG_DROUGHT_SMOOTH_MIN_VALID = Int(get(_drought_config, "smoothing_min_valid_days", 4))
const CFG_DROUGHT_THRESHOLD_METHOD = String(get(_drought_config, "threshold_method", "fixed"))
const CFG_DROUGHT_PERCENTILES = Int.(collect(get(_drought_config, "threshold_percentiles", [2, 5, 10, 20, 30])))
const CFG_DROUGHT_PLOTTING_POSITION = String(get(_drought_config, "plotting_position", "weibull"))
const CFG_DROUGHT_BELOW_RANGE_POLICY = String(get(_drought_config, "below_plotting_range_policy", "na"))
const CFG_DROUGHT_MIN_YEARS = Int(get(_drought_config, "min_years_for_threshold", 10))

# Fail fast on unimplemented options rather than silently falling back — a config
# asking for a method we do not implement must not quietly produce columns computed
# by a different one.
if CFG_DROUGHT_ENABLED
    if CFG_DROUGHT_THRESHOLD_METHOD != "fixed"
        error("Unsupported drought.threshold_method='$(CFG_DROUGHT_THRESHOLD_METHOD)': only \"fixed\" is implemented. " *
              "The paper's variable day-of-year method is intentionally out of scope " *
              "(docs/plans/2026-07-27-drought-signatures-plan.md §2.5).")
    end
    if CFG_DROUGHT_PLOTTING_POSITION != "weibull"
        error("Unsupported drought.plotting_position='$(CFG_DROUGHT_PLOTTING_POSITION)': only \"weibull\" (i/(n+1), Hyndman-Fan definition 6) is implemented.")
    end
    if !(CFG_DROUGHT_SMOOTH_ALIGNMENT in ("center", "trailing"))
        error("Unsupported drought.smoothing_alignment='$(CFG_DROUGHT_SMOOTH_ALIGNMENT)': expected \"center\" or \"trailing\".")
    end
    if !(CFG_DROUGHT_BELOW_RANGE_POLICY in ("na", "clamp"))
        error("Unsupported drought.below_plotting_range_policy='$(CFG_DROUGHT_BELOW_RANGE_POLICY)': expected \"na\" or \"clamp\".")
    end
    # Numeric sanity — a nonsensical value here would silently produce a plausible but
    # wrong metric (e.g. an EVEN "centered" window is asymmetric: (w-1)÷2 each side
    # spans w-1 days, not w).
    if CFG_DROUGHT_SMOOTH_WINDOW < 1
        error("drought.smoothing_window_days must be >= 1, got $(CFG_DROUGHT_SMOOTH_WINDOW).")
    end
    if CFG_DROUGHT_SMOOTH_ALIGNMENT == "center" && iseven(CFG_DROUGHT_SMOOTH_WINDOW)
        error("drought.smoothing_window_days must be ODD for alignment=\"center\" (an even window cannot be centered), got $(CFG_DROUGHT_SMOOTH_WINDOW).")
    end
    if !(1 <= CFG_DROUGHT_SMOOTH_MIN_VALID <= CFG_DROUGHT_SMOOTH_WINDOW)
        error("drought.smoothing_min_valid_days must satisfy 1 <= min_valid <= smoothing_window_days, got $(CFG_DROUGHT_SMOOTH_MIN_VALID) with window $(CFG_DROUGHT_SMOOTH_WINDOW).")
    end
    if isempty(CFG_DROUGHT_PERCENTILES)
        error("drought.threshold_percentiles must be non-empty.")
    end
    if any(p -> !(0 < p < 100), CFG_DROUGHT_PERCENTILES)
        error("drought.threshold_percentiles must all lie strictly inside (0, 100), got $(CFG_DROUGHT_PERCENTILES).")
    end
    if length(unique(CFG_DROUGHT_PERCENTILES)) != length(CFG_DROUGHT_PERCENTILES)
        error("drought.threshold_percentiles must be unique, got $(CFG_DROUGHT_PERCENTILES).")
    end
    if !issorted(CFG_DROUGHT_PERCENTILES)
        error("drought.threshold_percentiles must be ascending (column order and the level-monotonicity checks depend on it), got $(CFG_DROUGHT_PERCENTILES).")
    end
    if CFG_DROUGHT_MIN_YEARS < 1
        error("drought.min_years_for_threshold must be >= 1, got $(CFG_DROUGHT_MIN_YEARS).")
    end
end

# =============================================================================
# QA/QC Parameters
# =============================================================================
const CFG_QAQC_QANN_RANGE = Tuple(_config.qa_qc.qann_range)
const CFG_QAQC_BFI_RANGE = Tuple(_config.qa_qc.bfi_range)
const CFG_QAQC_FLASHINESS_RANGE = Tuple(_config.qa_qc.flashiness_range)
const CFG_QAQC_TQMEAN_RANGE = Tuple(_config.qa_qc.tqmean_range)
const CFG_QAQC_D50_RANGE = Tuple(_config.qa_qc.d50_range)
const CFG_QAQC_ELASTICITY_RANGE = Tuple(_config.qa_qc.elasticity_range)
const CFG_QAQC_RUNOFF_RATIO_RANGE = Tuple(_config.qa_qc.runoff_ratio_range)
const CFG_QAQC_SEASONAL_SUM_TOLERANCE = _config.qa_qc.seasonal_sum_tolerance
const CFG_QAQC_MAX_NA_FRACTION = _config.qa_qc.max_na_fraction
# flagged_for_high_na denominator manifest (2026-09-04). Absent section -> the
# same defaults, so old config variants keep the documented semantics.
const _high_na_cfg = get(_config.qa_qc, "high_na_denominator", Dict())
const CFG_QAQC_HIGH_NA_SUFFIXES = String.(collect(get(_high_na_cfg, "statistic_suffixes",
    ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
     "_mean", "_median", "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
     "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
     "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"])))
const CFG_QAQC_HIGH_NA_SCALARS = String.(collect(get(_high_na_cfg, "per_gage_scalars",
    ["elasticity_static", "log_a_seasonality_amplitude_all",
     "log_a_seasonality_amplitude_first_half", "log_a_seasonality_amplitude_last_half",
     "log_a_seasonality_minimum_all", "log_a_seasonality_minimum_first_half",
     "log_a_seasonality_minimum_last_half", "runoff_ratio_high_count",
     "elasticity_years_total", "elasticity_years_low_ppt", "ice_affected_days_total",
     "recession_alpha_point_cloud_linear_reservoir", "season_excluded_years_winter",
     "season_excluded_years_spring", "season_excluded_years_summer",
     "season_excluded_years_fall", "drought_threshold_fixed_p2", "drought_threshold_fixed_p5",
     "drought_threshold_fixed_p10", "drought_threshold_fixed_p20", "drought_threshold_fixed_p30"])))

# =============================================================================
# Metadata Parameters
# =============================================================================
const CFG_INCLUDE_HUMAN_INTERFERENCE = _config.metadata.include_human_interference

"""
    gages_ii_dir() -> String

GAGES-II metadata directory, resolved **at call time** from
`STREAMFLOW_GAGES_II_DIR`, falling back to `metadata.gages_ii_dir` in the config.

⚠️ This is a FUNCTION, not a `const`, deliberately. An ENV-derived constant is baked in
at precompile time (see the caching note at the top of this file), so a wrapper that sets
`ENV["STREAMFLOW_GAGES_II_DIR"]` at runtime would be silently ignored and the 11
human-interference columns would vanish from the output with only a warning. That
happened twice on 2026-07-28 before this was made a runtime lookup.
"""
gages_ii_dir() = let env_val = get(ENV, "STREAMFLOW_GAGES_II_DIR", "")
    env_val != "" ? env_val : String(_config.metadata.gages_ii_dir)
end

# Back-compat constant (JSON value only — does NOT see a runtime ENV override).
# Prefer `gages_ii_dir()`.
const CFG_GAGES_II_DIR = String(_config.metadata.gages_ii_dir)
const CFG_HYDAT_PATH = _config.metadata.hydat_path
const CFG_INTERFERENCE_COLUMNS = collect(_config.metadata.interference_columns)

# ==============================================================================
# NA HANDLING
# ==============================================================================
const _na_handling = get(_config, "na_handling", Dict())
const _na_interp = get(_na_handling, "interpolation", Dict())
const _na_reject = get(_na_handling, "year_rejection", Dict())
const _na_const_sd = get(_na_handling, "constant_sd_flag", Dict())
const _na_trend = get(_na_handling, "trend_completeness", Dict())
const _na_seasonal = get(_na_handling, "seasonal_completeness", Dict())
const _na_climate = get(_na_handling, "climate_na_policy", Dict())

const CFG_NA_MAX_GAP_DAYS = Int(get(_na_interp, "max_gap_days", 3))
const CFG_NA_INTERPOLATION_METHOD = get(_na_interp, "method", "linear")
const CFG_NA_INTERNAL_ONLY = get(_na_interp, "internal_only", true)
const CFG_NA_MAX_RAW_NA_PER_YEAR = Int(get(_na_reject, "max_raw_na_per_year", 30))
const CFG_NA_REJECT_NEGATIVE_FLOW = get(_na_reject, "reject_negative_flow", false)
const CFG_NA_REJECT_RESIDUAL_NA = get(_na_reject, "reject_residual_na", true)
const CFG_NA_CONSTANT_SD_ENABLED = get(_na_const_sd, "enabled", true)
const CFG_NA_CONSTANT_SD_MIN_DAYS = Int(get(_na_const_sd, "min_nonzero_days_per_month", 15))
const CFG_NA_CONSTANT_SD_MAX_UNIQUE = Int(get(_na_const_sd, "max_unique_values", 1))
const CFG_NA_TREND_MIN_FRACTION = Float64(get(_na_trend, "min_fraction", 0.80))
const CFG_NA_DECADE_MIN_FRACTION = Float64(get(_na_trend, "decade_min_fraction", 0.80))
const CFG_NA_SEASONAL_MIN_FRACTION = Float64(get(_na_seasonal, "min_fraction", 0.80))
const CFG_NA_SEASONAL_USE_RAW = get(_na_seasonal, "use_raw_observations", true)
const _na_season_defs = get(_na_seasonal, "season_definitions", Dict())
# Season months must be whole calendar months 1..12 (review round 2: a fractional or
# out-of-range month must fail loading, not be truncated or silently dropped)
function _cfg_months(section, key::String, default)::Vector{Int}
    v = collect(get(section, key, default))
    (!isempty(v) && all(x -> x isa Real && !(x isa Bool) && isfinite(x) && isinteger(x) && 1 <= x <= 12, v)) ||
        error("na_handling.seasonal_completeness.season_definitions.$key: expected calendar months 1..12, got $(repr(v))")
    return Int.(v)
end
const CFG_NA_SEASON_MONTHS = Dict{String, Vector{Int}}(
    "win" => _cfg_months(_na_season_defs, "winter", [12, 1, 2]),
    "spr" => _cfg_months(_na_season_defs, "spring", [3, 4, 5]),
    "sum" => _cfg_months(_na_season_defs, "summer", [6, 7, 8]),
    "fal" => _cfg_months(_na_season_defs, "fall", [9, 10, 11])
)
const CFG_NA_MAX_RAW_NA_PPT = Int(get(_na_climate, "max_raw_na_per_year_ppt", 30))
const CFG_NA_MAX_GAP_PPT = Int(get(_na_climate, "max_interpolation_gap_ppt", 3))
const CFG_NA_REJECT_NEGATIVE_PPT = get(_na_climate, "reject_negative_ppt", true)
const _na_snow = get(_na_handling, "snow_na_policy", Dict())
const CFG_NA_MAX_RAW_NA_SWE = Int(get(_na_snow, "max_raw_na_per_year_swe", 30))
const CFG_NA_MAX_GAP_SWE = Int(get(_na_snow, "max_interpolation_gap_swe", 3))
const CFG_NA_REJECT_NEGATIVE_SWE = get(_na_snow, "reject_negative_swe", true)
# Temperature (tmin/tmax) policy for the climate metrics (2026-10-12): same shape as
# the PPT/SWE policies, no negative-value check (temperatures are signed).
const _na_temp = get(_na_handling, "temperature_na_policy", Dict())
const CFG_NA_MAX_RAW_NA_TEMP = _cfg_int(_na_temp, "max_raw_na_per_year_temp", 30)
const CFG_NA_MAX_GAP_TEMP = _cfg_int(_na_temp, "max_interpolation_gap_temp", 3)
if CFG_NA_MAX_RAW_NA_TEMP < 0 || CFG_NA_MAX_GAP_TEMP < 0
    error("na_handling.temperature_na_policy: max_raw_na_per_year_temp and max_interpolation_gap_temp must be >= 0, got $(CFG_NA_MAX_RAW_NA_TEMP) / $(CFG_NA_MAX_GAP_TEMP).")
end

# Changepoint analysis configuration
const _cp_config = get(_config, "changepoint", Dict())
const CFG_CHANGEPOINT_ENABLED = Bool(get(_cp_config, "enabled", false))
const CFG_CP_START_WATER_YEAR = Int(get(_cp_config, "start_water_year", 1980))
const CFG_CP_END_WATER_YEAR = Int(get(_cp_config, "end_water_year", 2024))
const CFG_CP_MIN_TOTAL_OBS = Int(get(_cp_config, "min_total_obs", 20))
const CFG_CP_MIN_SEGMENT_OBS = Int(get(_cp_config, "min_segment_obs", 10))

# Legacy filtering flag
const CFG_USE_LEGACY_FILTERING = get(get(_config, "filtering", Dict()), "use_legacy_filtering", true)

# Stats floor (July 2026): minimum non-NA annual values a metric needs before ANY
# of its 8 statistics (+ changepoint fields) are computed. `nothing` when the
# section is absent => no floor (legacy behavior). Exemptions (recession,
# elasticity) are enforced at the orchestration layer, like trend_completeness.
const _stats_floor = get(_config, "stats_floor", Dict())
const CFG_MIN_VALUES_FOR_STATS = let
    v = get(_stats_floor, "min_values_for_stats", nothing)
    v === nothing ? nothing : Int(v)
end

# Annual values export (per-year annual signature values, long-format parquet).
# Default FALSE when the section is absent: config is loaded once at module import,
# so a permissive fallback would silently enable a large new artifact for old configs.
const _annual_values = get(_config, "annual_values", Dict())
const CFG_SAVE_ANNUAL_VALUES = Bool(get(_annual_values, "save", false))

"""
    get_config()

Return the full configuration dictionary.
"""
function get_config()
    return _config
end
