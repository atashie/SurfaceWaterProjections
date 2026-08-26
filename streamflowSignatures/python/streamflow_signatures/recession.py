"""
Recession analysis signatures.

Analyzes streamflow recession behavior using:
- Power-law fitting: dQ/dt = a * Q^b
- Point cloud and event-based methods
- Sinusoidal seasonality of recession parameters
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from .stats import generate_stats
from .config import RECESSION_MIN_EVENTS, RECESSION_MIN_LENGTH


def identify_recession_events(
    Q: np.ndarray,
    min_length: int = RECESSION_MIN_LENGTH
) -> List[Dict]:
    """
    Identify recession events using position-level marking.

    Marks each position where both recession criteria hold (Q decreasing AND
    |dQ/dt| decreasing), then finds contiguous runs of valid positions >=
    min_length. This captures the full extent of each recession event, unlike
    the previous look-ahead algorithm which could truncate events.

    Parameters
    ----------
    Q : np.ndarray
        Daily discharge values
    min_length : int, default 5
        Minimum number of valid transitions for a recession event

    Returns
    -------
    list of dict
        Each dict contains 'start', 'end', and 'indices' keys
    """
    n = len(Q)
    if n < min_length + 1:
        return []

    # Calculate dQ/dt (forward difference) with NaN at end
    dQ_dt = np.append(np.diff(Q), np.nan)

    # Step 1: Mark each position where both recession criteria hold
    # Position i is "valid" if Q[i+1] < Q[i] AND |dQdt[i+1]| < |dQdt[i]|
    is_valid = np.zeros(n, dtype=bool)
    for i in range(n - 1):
        if (not np.isnan(Q[i]) and not np.isnan(Q[i + 1]) and
                not np.isnan(dQ_dt[i]) and not np.isnan(dQ_dt[i + 1])):
            if Q[i + 1] < Q[i] and abs(dQ_dt[i + 1]) < abs(dQ_dt[i]):
                is_valid[i] = True

    # Step 2: Find contiguous runs of valid positions >= min_length
    # A run of k valid transitions at positions start to start+k-1 covers
    # days start to start+k (k+1 days).
    recession_events = []
    run_start = -1
    for i in range(n):
        if is_valid[i]:
            if run_start < 0:
                run_start = i
        else:
            if run_start >= 0:
                run_length = i - run_start  # number of valid transitions
                if run_length >= min_length:
                    # Event covers days run_start through i (inclusive)
                    recession_events.append({
                        'start': run_start,
                        'end': i,
                        'indices': list(range(run_start, i + 1))
                    })
                run_start = -1

    # Handle run ending at data boundary
    if run_start >= 0:
        run_length = n - run_start
        if run_length >= min_length:
            recession_events.append({
                'start': run_start,
                'end': n - 1,
                'indices': list(range(run_start, n))
            })

    return recession_events


def fit_recession_event(
    Q_values: np.ndarray,
    remove_first_day: bool = True
) -> Tuple[float, float]:
    """
    Fit recession parameters for a single event.

    Fits the power-law relationship: -dQ/dt = a * Q^b
    in log-log space: log(-dQ/dt) = log(a) + b*log(Q)

    Parameters
    ----------
    Q_values : np.ndarray
        Discharge values for the recession event
    remove_first_day : bool, default True
        Whether to remove the first day (often influenced by peak)

    Returns
    -------
    tuple
        (log_a, b) parameters, or (np.nan, np.nan) if fitting fails
    """
    if remove_first_day and len(Q_values) > 1:
        Q_values = Q_values[1:]

    n = len(Q_values)
    if n < 3:
        return np.nan, np.nan

    # Calculate -dQ/dt
    dQ_dt = -np.diff(Q_values)
    Q_subset = Q_values[:-1]

    # Remove non-positive values for log transformation
    valid_mask = (Q_subset > 0) & (dQ_dt > 0)
    if valid_mask.sum() < 2:
        return np.nan, np.nan

    Q_valid = Q_subset[valid_mask]
    dQ_dt_valid = dQ_dt[valid_mask]

    try:
        # Fit in log-log space
        result = scipy_stats.linregress(np.log(Q_valid), np.log(dQ_dt_valid))
        b = result.slope
        log_a = result.intercept
        return log_a, b
    except Exception:
        return np.nan, np.nan


def _event_log_a_b1(Q_event: np.ndarray) -> float:
    """Per-event log(a) under a FIXED linear reservoir (b = 1): the median of
    log(-dQ/dt) - log(Q) over the event's valid pairs — no regression.

    Port of julia/src/recession.jl `event_log_a_b1` (July 2026 convention,
    Phase 2 of the port campaign): fixing b = 1 removes the slope-intercept
    convolution of the free power-law fit. Conventions match the fitting path:
    first day removed (storm peak), pairs require Q > 0 and -dQ > 0, natural
    log. NaN when there are no valid pairs or fewer than 3 days.
    """
    Q_event = np.asarray(Q_event, dtype=float)
    if len(Q_event) < 3:
        return np.nan
    Q_work = Q_event[1:]            # Remove first day (storm peak)
    dQdt = -np.diff(Q_work)
    Q_mid = Q_work[:-1]             # Q at start of each interval
    valid = (Q_mid > 0) & (dQdt > 0) & ~np.isnan(Q_mid) & ~np.isnan(dQdt)
    if not valid.any():
        return np.nan
    return float(np.median(np.log(dQdt[valid]) - np.log(Q_mid[valid])))


def fit_sinusoidal_model(
    doy_values: np.ndarray,
    log_a_values: np.ndarray
) -> Tuple[float, float]:
    """
    Fit sinusoidal model to log(a) values across day of year.

    Fits: log(a) = A * sin(2*pi/365 * (doy - phi)) + C

    Parameters
    ----------
    doy_values : np.ndarray
        Day of year values
    log_a_values : np.ndarray
        Corresponding log(a) values

    Returns
    -------
    tuple
        (amplitude, minimum_doy) or (np.nan, np.nan) if fitting fails
    """
    # Remove NA values
    valid_mask = ~(np.isnan(log_a_values) | np.isnan(doy_values))
    if valid_mask.sum() < 10:
        return np.nan, np.nan

    doy_clean = doy_values[valid_mask]
    log_a_clean = log_a_values[valid_mask]

    try:
        # Create design matrix for sinusoidal fit
        X = np.column_stack([
            np.sin(2 * np.pi * doy_clean / 365),
            np.cos(2 * np.pi * doy_clean / 365),
            np.ones(len(doy_clean))
        ])

        # Fit linear model
        result = np.linalg.lstsq(X, log_a_clean, rcond=None)
        coeffs = result[0]

        B1 = coeffs[0]  # sin coefficient
        B2 = coeffs[1]  # cos coefficient
        # C = coeffs[2]  # intercept (not used)

        # Calculate amplitude and phase
        amplitude = np.sqrt(B1**2 + B2**2)

        # Calculate phase (in days)
        phase_rad = np.arctan2(-B2, B1)
        phase_days = phase_rad * 365 / (2 * np.pi)

        # Ensure phase is between 0 and 365
        if phase_days < 0:
            phase_days += 365

        # Minimum occurs at phase + 273.75 days (3/4 of a cycle)
        minimum_doy = phase_days + 273.75
        if minimum_doy > 365:
            minimum_doy -= 365

        return amplitude, minimum_doy

    except Exception:
        return np.nan, np.nan


def analyze_recession_parameters(
    streamflow_data: pd.DataFrame,
    min_events: int = RECESSION_MIN_EVENTS,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate recession parameter trends.

    Analyzes recession behavior using power-law fitting and
    calculates trends in recession parameters over time.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - dowy: int, day of water year
    min_events : int, default 25
        Minimum number of recession events required

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - log_a_pointcloud_*: Point cloud log(a) statistics
            - log_a_events_*: Event-based log(a) statistics
            - b_pointcloud_*: Point cloud b statistics
            - b_events_*: Event-based b statistics
            - concavity_*: Concavity statistics
            - log_a_seasonality_*: Seasonality metrics (single values)

    Notes
    -----
    - Power-law: dQ/dt = a * Q^b
    - Point cloud: Fits all recession points together
    - Events: Fits each event separately, takes median
    - Concavity: Difference in b between first and second half of event
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Define signatures
    signatures_with_stats = ["log_a_pointcloud", "log_a_events", "b_pointcloud",
                             "b_events", "concavity", "alpha_linear"]
    seasonality_signatures = [
        "log_a_seasonality_amplitude_all", "log_a_seasonality_minimum_all",
        "log_a_seasonality_amplitude_first_half", "log_a_seasonality_minimum_first_half",
        "log_a_seasonality_amplitude_last_half", "log_a_seasonality_minimum_last_half"
    ]

    # Sort by year and dowy
    df = df.sort_values(["water_year", "dowy"])

    years = df["water_year"].unique()

    # Initialize annual metrics
    annual_metrics = pd.DataFrame({
        "water_year": years,
        "log_a_pointcloud": np.nan,
        "log_a_events": np.nan,
        "b_pointcloud": np.nan,
        "b_events": np.nan,
        "concavity": np.nan,
        "n_recession_events": np.nan,
        "alpha_linear": np.nan,
    })

    # Store all recession events with timing
    all_recession_events = []
    all_alpha_linear = []

    # Process each year (df already sorted by water_year and dowy)
    for yr, year_data in df.groupby("water_year", sort=False):

        Q = year_data["Q"].values
        dowy = year_data["dowy"].values

        # Identify recession events
        recession_events = identify_recession_events(Q)

        # Store event count for this year (INDEPENDENT of min_events gate)
        yr_idx = annual_metrics.index[annual_metrics["water_year"] == yr]
        annual_metrics.loc[yr_idx, "n_recession_events"] = len(recession_events)

        # Collect all recession data for point cloud
        all_Q = []
        all_dQ_dt = []
        year_alpha_linear = []

        # Store individual event parameters. log_a values are the per-event
        # FIXED-b=1 values (linear reservoir; July 2026 canonical convention) —
        # the free-fit intercept and the median-b recalculation are gone: alpha
        # is decoupled from b everywhere. b itself stays a free fit.
        event_log_a_values = []
        event_b_values = []
        event_concavities = []

        # Process each recession event
        for event in recession_events:
            Q_event = Q[event['indices']]

            # Get middle day of recession for timing (feeds the seasonality
            # sinusoid). Canonical is the FLOOR midpoint of the index range
            # (julia/src/recession.jl `div(start_idx + end_idx, 2)`; rpkg's
            # 1-based `ceiling(n/2)` is algebraically the same). The previous
            # `indices[len // 2]` took the UPPER middle, putting even-length
            # events one day late — 8 of 12 events/year on a monthly-recession
            # gage (found 2026-08-25 by the Phase-2 cross-language check, which
            # the b=1 fix finally made visible).
            ev_idx = event['indices']
            mid_idx = (ev_idx[0] + ev_idx[-1]) // 2
            event_dowy = dowy[mid_idx]

            # Per-event FREE power-law fit (remove_first_day=True) — used for b
            # only; fit success also gates event acceptance (unchanged counts)
            log_a, b = fit_recession_event(Q_event, remove_first_day=True)

            # Per-event alpha under fixed b = 1 (same valid-pair definition as
            # the free fit, so non-NaN whenever the fit succeeds)
            log_a_b1 = _event_log_a_b1(Q_event)

            # Compute discrete recession constant alpha = Q_{i+1}/Q_i (b=1 linear reservoir)
            # INDEPENDENT of power-law fit — depends only on raw Q pairs
            # Remove first day (storm peak) consistent with point-cloud fitting
            if len(Q_event) > 2:  # Need at least 3 days (2 after removing first)
                Q_alpha = Q_event[1:]  # Remove first day (0-indexed)
                for j in range(len(Q_alpha) - 1):
                    if Q_alpha[j] > 0 and not np.isnan(Q_alpha[j]) and not np.isnan(Q_alpha[j + 1]):
                        a_i = Q_alpha[j + 1] / Q_alpha[j]
                        if 0 < a_i < 1:  # Must be valid recession (decreasing Q)
                            year_alpha_linear.append(a_i)
                            all_alpha_linear.append(a_i)

            if not np.isnan(log_a) and not np.isnan(b):
                event_log_a_values.append(log_a_b1)
                event_b_values.append(b)

                # Store event with timing (log_a is the b=1 value — the
                # seasonality sinusoid consumes these, per the July 2026
                # convention)
                all_recession_events.append({
                    'water_year': yr,
                    'dowy': event_dowy,
                    'log_a': log_a_b1,
                    'b': b
                })

                # Calculate concavity (overlapping split matching R)
                # R: first=Q[1:mid_point], second=Q[mid_point:n] (1-indexed, overlap at mid_point)
                # Python 0-indexed: first=Q[:mid_point], second=Q[mid_point-1:] (overlap at mid_point-1)
                if len(Q_event) >= 6:
                    mid_point = len(Q_event) // 2
                    first_half = Q_event[:mid_point]        # matches R's Q[1:mid_point]
                    second_half = Q_event[mid_point - 1:]   # matches R's Q[mid_point:n], overlap at mid_point-1

                    _, b_first = fit_recession_event(first_half, remove_first_day=False)
                    _, b_second = fit_recession_event(second_half, remove_first_day=False)

                    if not np.isnan(b_first) and not np.isnan(b_second):
                        concavity = b_second - b_first
                        event_concavities.append(concavity)

                # Add to point cloud data
                if len(Q_event) > 1:
                    Q_subset = Q_event[1:]  # Remove first day
                    dQ_subset = -np.diff(Q_event[1:])

                    valid_mask = (Q_subset[:-1] > 0) & (dQ_subset > 0)
                    if valid_mask.any():
                        all_Q.extend(Q_subset[:-1][valid_mask])
                        all_dQ_dt.extend(dQ_subset[valid_mask])

        # Point cloud analysis
        if len(all_Q) > 10:
            all_Q = np.array(all_Q)
            all_dQ_dt = np.array(all_dQ_dt)
            log_Q = np.log(all_Q)
            log_dQ_dt = np.log(all_dQ_dt)

            # log_a under fixed b = 1 (linear reservoir): log(a) = log(-dQ/dt)
            # - log(Q). No regression is involved, so this needs no singularity
            # gate and is independent of the free b fit below (July 2026
            # canonical convention; Phase 2 of the port campaign).
            annual_metrics.loc[annual_metrics["water_year"] == yr, "log_a_pointcloud"] = \
                np.median(log_dQ_dt - log_Q)

            # b from the FREE point-cloud fit (unchanged).
            # Skip near-singular data (matching R's lm() QR rank check);
            # R's tolerance is .Machine$double.eps^0.5 ≈ 1.49e-8
            if np.var(log_Q) >= 1e-8:
                result = scipy_stats.linregress(log_Q, log_dQ_dt)
                if not np.isnan(result.slope):
                    annual_metrics.loc[annual_metrics["water_year"] == yr, "b_pointcloud"] = result.slope

        # Event-based metrics. log_a_events uses the per-event fixed-b=1 values
        # — the previous median-b recalculation is gone: alpha is decoupled
        # from b by assuming a linear reservoir everywhere.
        if len(event_b_values) > 0:
            valid_b1 = [v for v in event_log_a_values if not np.isnan(v)]
            annual_metrics.loc[annual_metrics["water_year"] == yr, "log_a_events"] = \
                np.median(valid_b1) if valid_b1 else np.nan
            annual_metrics.loc[annual_metrics["water_year"] == yr, "b_events"] = np.median(event_b_values)

        # Concavity
        if len(event_concavities) > 0:
            annual_metrics.loc[annual_metrics["water_year"] == yr, "concavity"] = np.mean(event_concavities)

        # Per-year alpha_linear (linear reservoir assumption)
        if len(year_alpha_linear) > 10:
            annual_metrics.loc[annual_metrics["water_year"] == yr, "alpha_linear"] = np.median(year_alpha_linear)

    # n_recession_events stats computed INDEPENDENTLY of min_events gate
    events_df = annual_metrics[["water_year", "n_recession_events"]].dropna()
    if len(events_df) >= 3:
        n_events_stats = generate_stats(
            events_df, value_cols=["n_recession_events"], year_col="water_year",
            trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint, collector=collector,
        )
    else:
        n_events_stats = {}
        for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                       "_mk_rho", "_mk_pval", "_mean", "_median"]:
            n_events_stats[f"n_recession_events{suffix}"] = np.nan

    # Whole-record recession alpha (scalar — independent of min_events gate)
    recession_alpha_scalar = np.median(all_alpha_linear) if len(all_alpha_linear) > 10 else np.nan

    # Check minimum events requirement
    if len(all_recession_events) < min_events:
        # Return all NAs (but include n_recession_events stats)
        result = {}
        for sig in signatures_with_stats:
            for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                           "_mk_rho", "_mk_pval", "_mean", "_median"]:
                result[f"{sig}{suffix}"] = np.nan
        for sig in seasonality_signatures:
            result[sig] = np.nan
        result.update(n_events_stats)
        result["recession_alpha_point_cloud_linear_reservoir"] = recession_alpha_scalar
        return result

    # Generate statistics for main signatures
    result = generate_stats(
        annual_metrics,
        value_cols=signatures_with_stats,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, changepoint=changepoint, collector=collector,
    )

    # Add n_recession_events stats (computed independently above)
    result.update(n_events_stats)

    # Add seasonality signatures (single values, not trends)
    for sig in seasonality_signatures:
        result[sig] = np.nan

    # Calculate seasonality of log(a)
    if len(all_recession_events) >= 10:
        event_dowys = np.array([e['dowy'] for e in all_recession_events])
        event_log_a_values = np.array([e['log_a'] for e in all_recession_events])
        event_water_years = np.array([e['water_year'] for e in all_recession_events])

        # Fit to all data
        amplitude, minimum_doy = fit_sinusoidal_model(event_dowys, event_log_a_values)
        result["log_a_seasonality_amplitude_all"] = amplitude
        result["log_a_seasonality_minimum_all"] = minimum_doy

        # Split into first and last half of years
        median_water_year = np.median(np.unique(event_water_years))
        first_half_mask = event_water_years <= median_water_year
        last_half_mask = event_water_years > median_water_year

        # First half
        if first_half_mask.sum() >= 10:
            amp, min_doy = fit_sinusoidal_model(
                event_dowys[first_half_mask],
                event_log_a_values[first_half_mask]
            )
            result["log_a_seasonality_amplitude_first_half"] = amp
            result["log_a_seasonality_minimum_first_half"] = min_doy

        # Last half
        if last_half_mask.sum() >= 10:
            amp, min_doy = fit_sinusoidal_model(
                event_dowys[last_half_mask],
                event_log_a_values[last_half_mask]
            )
            result["log_a_seasonality_amplitude_last_half"] = amp
            result["log_a_seasonality_minimum_last_half"] = min_doy

    result["recession_alpha_point_cloud_linear_reservoir"] = recession_alpha_scalar
    return result
