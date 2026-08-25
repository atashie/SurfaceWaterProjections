"""
Changepoint detection for annual signature series (Pettitt test).

Port of the pipeline-reachable subset of julia/src/changepoint.jl (canonical):
``pettitt_test`` and ``segment_differential_metrics``, called from
``generate_stats``'s changepoint block. The Julia module's BIC model-selection
path (``detect_changepoint`` + helpers) is exported Julia-only API that the
production pipeline never calls — intentionally excluded from parity scope
(port plan 2026-08-24 §3; Codex F11).

Semantics pinned to the canonical implementation:
- NaN pairs filtered, then sorted by year (annual series have unique years).
- U_t computed via the recursive relation U_t = U_{t-1} + V_t with
  V_t = Σ_j sgn(x_t − x_j).
- The changepoint index is the FIRST t attaining the strict maximum of |U_t|
  over t ∈ [min_segment_obs, n − min_segment_obs] (1-based; Julia updates only
  on ``>``) — ties select the earliest candidate.
- Asymptotic p-value p = min(1, 2·exp(−6K²/(n³+n²))) (Pettitt 1979).
- Differential metrics split by year value (≤ cp_year vs > cp_year), require
  ≥3 values per segment, and reuse this package's own ``mann_kendall_test``
  (exactly as Julia reuses its own) for the per-segment p-values.
"""

from typing import Dict

import numpy as np


def pettitt_test(
    years,
    values,
    min_total_obs: int = 20,
    min_segment_obs: int = 10,
) -> Dict[str, float]:
    """
    Run the Pettitt non-parametric changepoint test on an annual time series.

    Parameters
    ----------
    years, values : array-like of float
        Annual series (NaN pairs are filtered; sorted by year internally).
    min_total_obs : int, default 20
        Minimum number of valid observations for the test to run.
    min_segment_obs : int, default 10
        Minimum observations on each side of a candidate changepoint.

    Returns
    -------
    dict with keys ``cp_year``, ``pval``, ``pre_mean``, ``post_mean``
    (all NaN when insufficient data).
    """
    na_result = {"cp_year": np.nan, "pval": np.nan,
                 "pre_mean": np.nan, "post_mean": np.nan}

    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)

    valid = ~(np.isnan(years) | np.isnan(values))
    yrs = years[valid]
    vals = values[valid]
    n = len(yrs)

    if n < min_total_obs:
        return na_result

    order = np.argsort(yrs, kind="stable")
    yrs = yrs[order]
    vals = vals[order]

    # U_t = cumsum(V_t), V_t = Σ_j sgn(x_t - x_j). n <= ~46 for annual series,
    # so the n x n sign matrix is trivial.
    V = np.sign(vals[:, None] - vals[None, :]).sum(axis=1)
    U = np.cumsum(V)

    # Candidate range t in [min_segment_obs, n - min_segment_obs] (1-based).
    # np.argmax returns the FIRST maximum — same tie rule as Julia's strict `>`.
    seg = np.abs(U[min_segment_obs - 1: n - min_segment_obs])
    if seg.size == 0:
        return na_result
    offset = int(np.argmax(seg))
    best_t = min_segment_obs + offset  # 1-based changepoint index
    K = float(seg[offset])

    pval = 2.0 * np.exp(-6.0 * K ** 2 / (float(n) ** 3 + float(n) ** 2))
    pval = min(pval, 1.0)

    pre_vals = vals[:best_t]
    post_vals = vals[best_t:]

    return {
        "cp_year": float(yrs[best_t - 1]),
        "pval": float(pval),
        "pre_mean": float(np.mean(pre_vals)),
        "post_mean": float(np.mean(post_vals)),
    }


def segment_differential_metrics(years, values, cp_year: float) -> Dict[str, float]:
    """
    Differential metrics at a changepoint: delta_mean, pct_change, and pre/post
    Mann-Kendall p-values. NaN dict when cp_year is NaN or a segment has <3
    values.
    """
    na_result = {"delta_mean": np.nan, "pct_change": np.nan,
                 "pre_mk_pval": np.nan, "post_mk_pval": np.nan}

    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)

    if np.isnan(cp_year) or len(years) == 0:
        return na_result

    pre_mask = years <= cp_year
    post_mask = years > cp_year
    pre_vals = values[pre_mask]
    post_vals = values[post_mask]

    if len(pre_vals) < 3 or len(post_vals) < 3:
        return na_result

    pre_m = float(np.mean(pre_vals))
    post_m = float(np.mean(post_vals))
    delta_mean = post_m - pre_m
    pct_change = delta_mean / abs(pre_m) * 100.0 if abs(pre_m) > 1e-10 else np.nan

    # Deferred import avoids a stats <-> changepoint import cycle; reusing the
    # package's own MK implementation mirrors Julia's structure exactly.
    from .stats import mann_kendall_test
    _, pre_mk_pval = mann_kendall_test(pre_vals)
    _, post_mk_pval = mann_kendall_test(post_vals)

    return {"delta_mean": delta_mean, "pct_change": pct_change,
            "pre_mk_pval": pre_mk_pval, "post_mk_pval": post_mk_pval}
