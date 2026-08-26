"""
Streamflow drought signatures (Adelsperger et al., in review).

Port of julia/src/drought.jl (canonical). Per-water-year drought duration and
deficit against FIXED whole-record percentile thresholds, at five severity
levels mirroring the U.S. Drought Monitor ladder (2/5/10/20/30 % =
D4/D3/D2/D1/D0).

Conventions pinned to the canonical implementation:
- 7-day CENTERED smoothing applied to the continuous date-indexed series within
  maximal runs of consecutive dates, so the window never averages across a
  temporal gap (rejected years, duplicate or missing dates all break a run);
  the window shrinks at run edges and a day whose in-run window holds fewer
  than `min_valid` values is NaN.
- Thresholds are percentiles of the smoothed values pooled over the whole
  record, computed with the unbiased Weibull plotting position `p_i = i/(n+1)`
  (Hyndman & Fan definition 6) — deliberately NOT the type-7 default used
  elsewhere in this project, which diverges most in exactly this low tail.
  A probability below the smallest plotting position is NOT interpolable and
  yields NaN rather than clamping to the sample minimum.
- Per water year and level: count days with `Q_smooth < thr` (STRICT) and sum
  `thr - Q_smooth` over those days.
- Both metrics are dense series whose ZEROS are meaningful (a wet year has no
  drought days), so the trend gate and stats floor apply normally — this family
  is NOT exempt, and no snow-style record-anchored gate is needed.
- `drought_deficit_*` is unit-carrying (mm only when Q is area-normalized); the
  duration metrics are scale-invariant.

References: Adelsperger et al. (in review); Laaha, G. et al. (2017) — unbiased
plotting positions for low-flow frequency analysis.
"""

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .stats import generate_stats
from .config import (
    DROUGHT_PERCENTILES,
    DROUGHT_SMOOTH_WINDOW,
    DROUGHT_SMOOTH_ALIGNMENT,
    DROUGHT_SMOOTH_MIN_VALID,
    DROUGHT_MIN_YEARS,
    DROUGHT_BELOW_RANGE_POLICY,
)

# Metric names are derived from the configured levels so config and columns
# cannot drift apart. Durations first, then deficits (level order as configured).
DROUGHT_DURATION_METRICS = [f"drought_duration_fixed_p{p}" for p in DROUGHT_PERCENTILES]
DROUGHT_DEFICIT_METRICS = [f"drought_deficit_fixed_p{p}" for p in DROUGHT_PERCENTILES]
DROUGHT_METRICS = DROUGHT_DURATION_METRICS + DROUGHT_DEFICIT_METRICS

# Per-gage threshold values (mm/day), emitted as scalar diagnostics — not 8-stat
# metrics — so the thresholds behind the metrics are auditable in the delivered
# CSV (and gages whose low-level threshold is exactly 0 are identifiable).
DROUGHT_THRESHOLD_SCALARS = [f"drought_threshold_fixed_p{p}" for p in DROUGHT_PERCENTILES]


def weibull_quantile(x, p: float, below_range: str = "na") -> float:
    """
    Quantile at non-exceedance probability `p` using the unbiased Weibull
    plotting position `p_i = i/(n+1)` (Hyndman & Fan definition 6).

    NaNs are dropped. `p` below the smallest plotting position `1/(n+1)` is NOT
    interpolable: with ``below_range="na"`` the result is NaN (honest — the
    sample cannot resolve that probability), with ``"clamp"`` it is the sample
    minimum. Above `n/(n+1)` the sample maximum is returned.
    """
    vals = np.asarray(x, dtype=float)
    vals = np.sort(vals[~np.isnan(vals)])
    n = len(vals)
    if n == 0:
        return np.nan

    h = (n + 1) * float(p)
    if h < 1.0:
        return float(vals[0]) if below_range == "clamp" else np.nan
    if h >= n:
        return float(vals[n - 1])

    fl = int(np.floor(h))
    frac = h - fl
    # 1-based -> 0-based indexing
    return float(vals[fl - 1] + frac * (vals[fl] - vals[fl - 1]))


def _day_numbers(dates):
    """Integer day numbers plus a validity mask. Unparseable/missing dates are
    marked invalid and never treated as contiguous with a neighbour."""
    d = pd.to_datetime(pd.Series(list(dates)), errors="coerce")
    ok = d.notna().to_numpy()
    nums = np.zeros(len(d), dtype=np.int64)
    if ok.any():
        # days since epoch; only the DIFFERENCE between neighbours matters
        nums[ok] = (d[ok].to_numpy().astype("datetime64[D]").astype(np.int64))
    return nums, ok


def smooth_daily_flow(dates, Q, window: int = None, alignment: str = None,
                      min_valid: int = None) -> np.ndarray:
    """
    Moving-average smoothing of a daily series, applied within maximal runs of
    consecutive calendar dates so the window never averages across a temporal
    gap (rejected years, duplicated dates, or missing dates all break a run).

    ``alignment="center"`` uses ±(window-1)//2 days; ``"trailing"`` uses the
    `window` days ending on the day itself. The window shrinks at run edges; a
    day whose in-run window contains fewer than `min_valid` non-NaN values
    yields NaN. `dates` must be sorted ascending (callers sort first).
    """
    window = DROUGHT_SMOOTH_WINDOW if window is None else window
    alignment = DROUGHT_SMOOTH_ALIGNMENT if alignment is None else alignment
    min_valid = DROUGHT_SMOOTH_MIN_VALID if min_valid is None else min_valid

    Q = np.asarray(Q, dtype=float)
    n = len(Q)
    out = np.full(n, np.nan)
    if n == 0:
        return out
    if len(dates) != n:
        raise ValueError(f"smooth_daily_flow: dates and Q length mismatch "
                         f"({len(dates)} vs {n})")

    nums, ok = _day_numbers(dates)
    half = (window - 1) // 2

    run_start = 0
    for i in range(n):
        # A run ends at i when i is the last row, either endpoint is invalid, or
        # the next date is not exactly one day later (gap OR duplicate).
        ends_run = (i == n - 1) or (not ok[i]) or (not ok[i + 1]) \
            or (nums[i + 1] != nums[i] + 1)
        if not ends_run:
            continue
        # A row with a missing/unparseable date belongs to NO run: it must not be
        # averaged with its neighbours, so the run stops before it and its own
        # output stays NaN.
        run_end = i if ok[i] else i - 1
        for j in range(run_start, run_end + 1):
            if alignment == "center":
                lo, hi = max(run_start, j - half), min(run_end, j + half)
            else:  # trailing
                lo, hi = max(run_start, j - window + 1), j
            seg = Q[lo:hi + 1]
            valid = seg[~np.isnan(seg)]
            out[j] = valid.mean() if len(valid) >= min_valid else np.nan
        run_start = i + 1

    return out


def calculate_drought_metrics(
    streamflow_data: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate streamflow drought duration and deficit with trend statistics.

    For each configured percentile level `p` (default 2, 5, 10, 20, 30):
    - ``drought_duration_fixed_p{p}`` — days per water year below the threshold
    - ``drought_deficit_fixed_p{p}``  — summed departures below it (mm)

    Each produces 8 statistics via generate_stats() (+ 8 Pettitt fields when
    enabled), plus per-gage scalar diagnostics ``drought_threshold_fixed_p{p}``.
    """
    result: Dict[str, float] = {}

    def _all_nan_result():
        """Every early-return path emits the IDENTICAL key set (including
        changepoint fields when enabled) via the same explicit-value_cols
        machinery — schema contract."""
        empty_annual = pd.DataFrame({"water_year": pd.Series([], dtype="int64")})
        for m in DROUGHT_METRICS:
            empty_annual[m] = pd.Series([], dtype="float64")
        result.update(generate_stats(
            empty_annual, value_cols=DROUGHT_METRICS, year_col="water_year",
            trend_completeness=trend_completeness,
            decade_completeness=decade_completeness,
            min_values_for_stats=min_values_for_stats,
            changepoint=changepoint, collector=collector))
        for s in DROUGHT_THRESHOLD_SCALARS:
            result[s] = np.nan
        return result

    missing = [c for c in ("Q", "water_year", "date") if c not in streamflow_data.columns]
    if missing:
        import logging
        logging.getLogger(__name__).warning(
            "calculate_drought_metrics: Missing columns: %s", missing)
        return _all_nan_result()
    if len(streamflow_data) == 0:
        return _all_nan_result()

    dfs = streamflow_data.sort_values("date")
    Q_smooth = smooth_daily_flow(dfs["date"].to_numpy(),
                                 pd.to_numeric(dfs["Q"], errors="coerce").to_numpy(dtype=float))

    wy = pd.to_numeric(dfs["water_year"], errors="coerce").to_numpy(dtype=float)
    years = sorted({int(y) for y in wy if not np.isnan(y)})
    if not years:
        return _all_nan_result()

    # Fixed thresholds: whole-record percentiles of the smoothed series. Short
    # records get NaN thresholds (and hence NaN metrics) rather than an unstable
    # low tail.
    n_levels = len(DROUGHT_PERCENTILES)
    thresholds = np.full(n_levels, np.nan)
    if len(years) >= DROUGHT_MIN_YEARS:
        pool = Q_smooth[~np.isnan(Q_smooth)]
        if pool.size:
            for k, p in enumerate(DROUGHT_PERCENTILES):
                thresholds[k] = weibull_quantile(
                    pool, p / 100.0, below_range=DROUGHT_BELOW_RANGE_POLICY)

    annual = pd.DataFrame({"water_year": years})
    for m in DROUGHT_METRICS:
        annual[m] = np.nan

    for yr_idx, yr in enumerate(years):
        year_mask = wy == yr
        if not year_mask.any():
            continue
        qs = Q_smooth[year_mask]
        qs_valid = qs[~np.isnan(qs)]

        for k in range(n_levels):
            thr = thresholds[k]
            if np.isnan(thr):
                continue
            below = qs_valid[qs_valid < thr]       # STRICT comparison
            annual.at[yr_idx, DROUGHT_DURATION_METRICS[k]] = float(below.size)
            annual.at[yr_idx, DROUGHT_DEFICIT_METRICS[k]] = float((thr - below).sum())

    result.update(generate_stats(
        annual, value_cols=DROUGHT_METRICS, year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness,
        min_values_for_stats=min_values_for_stats,
        changepoint=changepoint, collector=collector))

    for k, s in enumerate(DROUGHT_THRESHOLD_SCALARS):
        result[s] = float(thresholds[k])

    return result
