"""
Core statistical functions for streamflow signature analysis.

This module provides the statistics engine that produces 8 standard statistics
for each signature metric, matching the R implementation's generate_stats().
"""

from typing import Optional, Dict, List, Union
import numpy as np
import pandas as pd
from scipy import stats

from .changepoint import pettitt_test, segment_differential_metrics

# Stat suffixes matching the canonical Julia implementation (stats.jl)
STAT_SUFFIXES = [
    "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
    "_mk_rho", "_mk_pval", "_mean", "_median",
]

# Changepoint column suffixes (Pettitt test only)
CP_SUFFIXES = [
    "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean",
    "_pettitt_delta_mean", "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval",
]


class AnnualCollector:
    """Opt-in accumulator for per-year signature values (long format).

    Port of the `AnnualCollector` struct in julia/src/stats.jl. When passed to
    ``generate_stats``, the annual series is appended EXACTLY as received —
    before any min_rows, stats-floor or trend-completeness gating — so the
    collected rows are precisely what the summary statistics were computed
    from. Collection is read-only with respect to the statistics.

    One instance per gage; the benchmark runner drains it into a global
    long-format table tagged with gage_id.
    """

    __slots__ = ("signature", "water_year", "value")

    def __init__(self):
        self.signature = []
        self.water_year = []
        self.value = []

    def __len__(self):
        return len(self.value)


def _collect_annual(collector: "AnnualCollector", data: pd.DataFrame,
                    value_cols: List[str], year_col: str) -> None:
    """Append the annual series as passed to generate_stats (long format).

    Rows whose year is missing or non-finite are skipped (they cannot be
    keyed); missing values become NaN. NaN values are KEPT — they record
    "year present, metric not computable".
    """
    years_raw = pd.to_numeric(data[year_col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(years_raw)
    if not keep.any():
        return
    keyed = np.rint(years_raw[keep]).astype(np.int32)
    n_keep = int(keep.sum())

    for col in value_cols:
        if col not in data.columns:
            continue
        vals = pd.to_numeric(data[col], errors="coerce").to_numpy(dtype=float)[keep]
        collector.signature.extend([col] * n_keep)
        collector.water_year.extend(keyed.tolist())
        collector.value.extend(vals.tolist())


def empty_stats(metric: str) -> Dict[str, float]:
    """
    Return a dict with NaN values for all 8 statistics for a metric.

    Mirrors Julia's ``empty_stats`` exactly: 8 stat keys ONLY — no changepoint
    keys (gages routed through empty-stats branches lack the ``_pettitt_*``
    keys, which surface as missing/NaN via the runner's union schema, identical
    to the canonical output).
    """
    return {f"{metric}{suffix}": np.nan for suffix in STAT_SUFFIXES}


def _run_changepoint_block(
    results: Dict[str, float],
    col: str,
    valid_years: np.ndarray,
    valid_values: np.ndarray,
    changepoint: dict,
) -> None:
    """Run the Pettitt test + differential metrics; write 8 keys into results.

    Mirrors julia/src/stats.jl `_run_changepoint_block!`: filter to the
    changepoint analysis window, Pettitt (4 keys), differential metrics (4 keys).
    """
    cp_start = changepoint.get("start_year", 1980)
    cp_end = changepoint.get("end_year", 2024)
    cp_min_obs = changepoint.get("min_total_obs", 20)
    cp_min_seg = changepoint.get("min_segment_obs", 10)

    cp_mask = (valid_years >= cp_start) & (valid_years <= cp_end)
    cp_years = valid_years[cp_mask]
    cp_vals = valid_values[cp_mask]

    pettitt = pettitt_test(cp_years, cp_vals,
                           min_total_obs=cp_min_obs, min_segment_obs=cp_min_seg)
    results[f"{col}_pettitt_cp_year"] = pettitt["cp_year"]
    results[f"{col}_pettitt_pval"] = pettitt["pval"]
    results[f"{col}_pettitt_pre_mean"] = pettitt["pre_mean"]
    results[f"{col}_pettitt_post_mean"] = pettitt["post_mean"]

    pettitt_diff = segment_differential_metrics(cp_years, cp_vals, pettitt["cp_year"])
    results[f"{col}_pettitt_delta_mean"] = pettitt_diff["delta_mean"]
    results[f"{col}_pettitt_pct_change"] = pettitt_diff["pct_change"]
    results[f"{col}_pettitt_pre_mk_pval"] = pettitt_diff["pre_mk_pval"]
    results[f"{col}_pettitt_post_mk_pval"] = pettitt_diff["post_mk_pval"]


def theil_sen_slope(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the Theil-Sen slope estimator.

    This is a robust non-parametric estimator of linear trend, calculated as
    the median of all pairwise slopes between points.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (typically years)
    y : np.ndarray
        Dependent variable (metric values)

    Returns
    -------
    float
        Theil-Sen slope estimate, or np.nan if calculation fails

    Notes
    -----
    Matches R's zyp::zyp.sen implementation.
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 2:
        return np.nan

    try:
        result = stats.theilslopes(y_clean, x_clean)
        return result.slope
    except Exception:
        return np.nan


def mann_kendall_test(y: np.ndarray) -> tuple:
    """
    Perform the Mann-Kendall trend test.

    The Mann-Kendall test is a non-parametric test for monotonic trend.

    Parameters
    ----------
    y : np.ndarray
        Time series values (assumed to be in temporal order)

    Returns
    -------
    tuple
        (tau, p_value) where tau is Kendall's tau and p_value is the
        two-sided p-value for the test. Returns (np.nan, np.nan) if
        calculation fails.

    Notes
    -----
    Computed directly (not via ``scipy.stats.kendalltau``) so that BOTH outputs
    match canonical Julia — which is itself a faithful port of R's
    ``Kendall::MannKendall`` (verified identical to 1.000 at n = 10…33):

    - tau is tau-b: ``S / sqrt((n_pairs - T_y) * n_pairs)``. The time index has
      no ties, so only ties in ``y`` enter. Algebraically identical to scipy's
      tau-b for this input; computing it here removes the last-ulp difference.
    - the p-value uses the tie-corrected variance AND the **continuity
      correction** ``z = (S ∓ 1)/sqrt(var_S)`` — the standard adjustment when
      approximating a discrete statistic with a normal.

    Why not scipy's p-value (changed 2026-08-26, port campaign): scipy selects
    on TIES, not sample size — with no ties it returns the EXACT Kendall
    p-value, and with ties it uses an asymptotic branch with NO continuity
    correction. The second case dominates here, because day counts, pulse
    counts and drought durations are tie-heavy: it shifted ~0.24 % of
    main-path and ~0.45 % of Pettitt-segment significance calls at alpha = 0.05
    relative to Julia/R.
    """
    # Remove NaN values
    y_clean = np.asarray(y, dtype=float)
    y_clean = y_clean[~np.isnan(y_clean)]
    n = len(y_clean)

    if n < 3:
        return np.nan, np.nan

    try:
        # S = sum over i<j of sign(y_j - y_i). n <= ~46 here, so the pairwise
        # matrix is trivial; np.subtract.outer(y, y)[i, j] = y_i - y_j, hence
        # the sign flip on the strict upper triangle.
        S = float(-np.triu(np.sign(np.subtract.outer(y_clean, y_clean)), 1).sum())

        # Tie-corrected variance of S
        _, tie_groups = np.unique(y_clean, return_counts=True)
        ties = tie_groups[tie_groups > 1].astype(float)
        tie_correction = float((ties * (ties - 1) * (2 * ties + 5)).sum())
        var_S = (n * (n - 1) * (2 * n + 5) - tie_correction) / 18.0
        if var_S <= 0:
            return np.nan, np.nan

        # tau-b (x = time indices, so only ties in y contribute)
        n_pairs = n * (n - 1) / 2.0
        T_y = float((tie_groups * (tie_groups - 1) / 2.0).sum())
        denom = np.sqrt((n_pairs - T_y) * n_pairs)
        tau = S / denom if denom > 0.0 else np.nan

        # Continuity-corrected normal approximation
        if S > 0:
            z = (S - 1) / np.sqrt(var_S)
        elif S < 0:
            z = (S + 1) / np.sqrt(var_S)
        else:
            z = 0.0
        p_value = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

        # R's Kendall::MannKendall fails for n<=3 (IFAULT=12, returns p=1.0);
        # Julia mirrors that. Match it.
        if n <= 3:
            p_value = 1.0

        return tau, p_value
    except Exception:
        return np.nan, np.nan


def spearman_correlation(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Calculate Spearman rank correlation and p-value.

    Parameters
    ----------
    x : np.ndarray
        First variable (typically years)
    y : np.ndarray
        Second variable (metric values)

    Returns
    -------
    tuple
        (rho, p_value) where rho is Spearman's correlation coefficient
        and p_value is the two-sided p-value. Returns (np.nan, np.nan)
        if calculation fails.
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 3:
        return np.nan, np.nan

    try:
        rho, p_value = stats.spearmanr(x_clean, y_clean)
        # Handle constant input: scipy returns (NaN, NaN) for constant series.
        # R's cor.test returns rho=NA but we match by returning (NaN, NaN) —
        # the key fix is for MK where R returns (0, 1.0) instead of NaN.
        return rho, p_value
    except Exception:
        return np.nan, np.nan


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate ordinary least squares linear regression slope.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (typically years)
    y : np.ndarray
        Dependent variable (metric values)

    Returns
    -------
    float
        Linear regression slope, or np.nan if calculation fails
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 2:
        return np.nan

    try:
        result = stats.linregress(x_clean, y_clean)
        return result.slope
    except Exception:
        return np.nan


def generate_stats(
    data: pd.DataFrame,
    value_cols: Optional[List[str]] = None,
    year_col: str = "water_year",
    min_rows: int = 3,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    changepoint: Optional[dict] = None,
    min_values_for_stats: Optional[int] = None,
    force_skip_trends: Optional[set] = None,
    collector: Optional["AnnualCollector"] = None,
) -> Dict[str, float]:
    """
    Generate 8 standard statistics for each metric column.

    This is the core statistics engine that produces standardized trend
    and summary statistics for each signature metric.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing annual metric values with a year column
    value_cols : list of str, optional
        Columns to analyze. If None, uses all numeric columns except year_col
    year_col : str, default "water_year"
        Name of the year column
    min_rows : int, default 3
        Minimum number of non-NA rows required for statistics
    trend_completeness : float, optional
        Minimum fraction of non-NaN values (0-1) required across the full
        record for trend statistics to be computed. If the fraction of
        valid values is below this threshold, the 6 trend statistics
        (slopes, correlations, p-values) are set to NaN while mean and
        median are still computed. If None, no completeness check is applied.
    decade_completeness : float, optional
        Minimum fraction of non-NaN values required in both the first and
        last decade of the record for trend statistics. This prevents
        trends being driven by data concentrated in one end of the record.
        Skipped if the year span is < 10 years. If None, no decade check
        is applied.
    changepoint : dict, optional
        Changepoint (Pettitt) analysis configuration with keys ``start_year``,
        ``end_year``, ``min_total_obs``, ``min_segment_obs``. When provided,
        8 additional ``_pettitt_*`` fields are emitted per metric. Changepoint
        analysis runs INDEPENDENTLY of the trend-completeness gates (matching
        Julia): a trend-gated metric still gets its Pettitt fields, while a
        metric below ``min_rows`` gets NaN Pettitt fields.
    min_values_for_stats : int, optional
        Stats floor (July 2026 Julia feature): a metric with fewer non-NaN
        annual values than this emits NaN for ALL 8 statistics AND its
        changepoint fields — a 4-point clustered series must not carry a
        Theil-Sen slope. Recession and elasticity are exempt at the
        orchestration layer (the kwarg is never passed to them).
    force_skip_trends : set of str, optional
        Metric names whose 6 trend statistics are suppressed through the SAME
        path as trend_completeness — mean/median and changepoint unaffected.
        Used by caller-computed gates (e.g. the snow record-anchored decade
        gate).
    collector : AnnualCollector, optional
        When supplied, the annual series is appended to it in long format
        BEFORE any gating (export contract). Read-only with respect to the
        statistics: passing a collector never changes the returned values.

    Returns
    -------
    dict
        Dictionary with keys for each statistic. For each metric column,
        8 statistics are produced with suffixes:
            - _senn_slp: Theil-Sen slope
            - _linear_slp: Linear regression slope
            - _spearman_rho: Spearman correlation coefficient
            - _spearman_pval: Spearman p-value
            - _mk_rho: Mann-Kendall tau
            - _mk_pval: Mann-Kendall p-value
            - _mean: Arithmetic mean
            - _median: Median

    Examples
    --------
    >>> data = pd.DataFrame({
    ...     'water_year': [2000, 2001, 2002, 2003, 2004],
    ...     'Qann': [100, 110, 105, 120, 115]
    ... })
    >>> stats = generate_stats(data, value_cols=['Qann'])
    >>> print(stats['Qann_mean'])
    110.0
    """
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input 'data' must be a DataFrame")

    if year_col not in data.columns:
        raise ValueError(f"Year column '{year_col}' not found in data")

    # Resolve value columns BEFORE sorting/gating (mirrors Julia's hoisted
    # resolution — eltypes are unaffected by sorting)
    if value_cols is None:
        numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
        value_cols = [c for c in numeric_cols if c != year_col]

    results = {}

    # Opt-in annual-values collection: capture the series exactly as passed,
    # BEFORE any gating below (matching Julia's placement).
    if collector is not None:
        _collect_annual(collector, data, value_cols, year_col)

    # Frame-level gate (mirrors Julia's earliest nrow < min_rows return): every
    # REQUESTED metric — including requested-but-absent columns — emits the 8
    # NaN statistics (+ NaN changepoint fields when enabled). Added 2026-08-25
    # (Codex Phase-1 review): the per-column path below silently skips absent
    # columns, which diverges from Julia on short frames.
    if len(data) < min_rows:
        for col in value_cols:
            for suffix in STAT_SUFFIXES:
                results[f"{col}{suffix}"] = np.nan
            if changepoint is not None:
                for cp_suffix in CP_SUFFIXES:
                    results[f"{col}{cp_suffix}"] = np.nan
        return results

    # Sort by year for trend calculations
    data = data.sort_values(year_col).copy()

    for col in value_cols:
        if col not in data.columns:
            # Skip missing columns silently (no stats emitted for them)
            continue

        # Extract non-NA data
        mask = ~data[col].isna()
        years = data.loc[mask, year_col].values.astype(float)
        values = data.loc[mask, col].values.astype(float)

        # Check minimum data requirement. The stats floor (min_values_for_stats)
        # gates through the same branch as min_rows: below the floor, ALL 8
        # statistics AND the changepoint fields are NaN (matching Julia; the
        # annual-values collector — Phase 3 — collects BEFORE this gate).
        if len(values) < min_rows or (
                min_values_for_stats is not None and len(values) < min_values_for_stats):
            results[f"{col}_senn_slp"] = np.nan
            results[f"{col}_linear_slp"] = np.nan
            results[f"{col}_spearman_rho"] = np.nan
            results[f"{col}_spearman_pval"] = np.nan
            results[f"{col}_mk_rho"] = np.nan
            results[f"{col}_mk_pval"] = np.nan
            results[f"{col}_mean"] = np.nan
            results[f"{col}_median"] = np.nan
            if changepoint is not None:
                for cp_suffix in CP_SUFFIXES:
                    results[f"{col}{cp_suffix}"] = np.nan
            continue

        # Trend completeness check: suppress trend stats if data is too sparse.
        # Uses the metric's own non-NaN year range for decade boundaries, so
        # climate signatures use climate years and recession uses recession years.
        suppress_trends = False

        # Externally-forced trend skip (snow record-anchored decade gate etc.):
        # suppresses the 6 trend stats; mean/median + changepoint unaffected.
        if force_skip_trends is not None and col in force_skip_trends:
            suppress_trends = True

        if not suppress_trends and trend_completeness is not None and len(years) > 0:
            metric_yr_min = float(years.min())
            metric_yr_max = float(years.max())
            metric_span = int(metric_yr_max - metric_yr_min) + 1

            # Overall: non-NaN values / span of metric's year range
            overall_frac = len(values) / metric_span if metric_span > 0 else 0.0
            if overall_frac < trend_completeness:
                suppress_trends = True

        if not suppress_trends and decade_completeness is not None and len(years) > 0:
            metric_yr_min = float(years.min())
            metric_yr_max = float(years.max())
            metric_span = int(metric_yr_max - metric_yr_min) + 1

            if metric_span >= 10:
                # First decade: metric_yr_min to metric_yr_min + 9
                first_expected = min(10, metric_span)
                first_decade_mask = years <= metric_yr_min + 9
                first_decade_n = int(np.sum(first_decade_mask))
                if first_expected > 0 and first_decade_n / first_expected < decade_completeness:
                    suppress_trends = True

                # Last decade: metric_yr_max - 9 to metric_yr_max
                if not suppress_trends:
                    last_start = metric_yr_max - 9
                    last_expected = min(10, int(metric_yr_max - max(last_start, metric_yr_min)) + 1)
                    last_decade_mask = years >= last_start
                    last_decade_n = int(np.sum(last_decade_mask))
                    if last_expected > 0 and last_decade_n / last_expected < decade_completeness:
                        suppress_trends = True

        if suppress_trends:
            # Set 6 trend stats to NaN, still compute mean/median.
            # Changepoint runs INDEPENDENTLY of trend gating (matching Julia).
            results[f"{col}_senn_slp"] = np.nan
            results[f"{col}_linear_slp"] = np.nan
            results[f"{col}_spearman_rho"] = np.nan
            results[f"{col}_spearman_pval"] = np.nan
            results[f"{col}_mk_rho"] = np.nan
            results[f"{col}_mk_pval"] = np.nan
            results[f"{col}_mean"] = np.mean(values)
            results[f"{col}_median"] = np.median(values)
            if changepoint is not None:
                _run_changepoint_block(results, col, years, values, changepoint)
            continue

        # Calculate Theil-Sen slope
        results[f"{col}_senn_slp"] = theil_sen_slope(years, values)

        # Calculate linear regression slope
        results[f"{col}_linear_slp"] = linear_slope(years, values)

        # Calculate Spearman correlation
        rho, pval = spearman_correlation(years, values)
        results[f"{col}_spearman_rho"] = rho
        results[f"{col}_spearman_pval"] = pval

        # Calculate Mann-Kendall test
        mk_tau, mk_pval = mann_kendall_test(values)
        results[f"{col}_mk_rho"] = mk_tau
        results[f"{col}_mk_pval"] = mk_pval

        # Calculate mean and median
        results[f"{col}_mean"] = np.mean(values)
        results[f"{col}_median"] = np.median(values)

        # Changepoint analysis (opt-in)
        if changepoint is not None:
            _run_changepoint_block(results, col, years, values, changepoint)

    return results
