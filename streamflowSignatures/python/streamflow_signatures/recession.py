"""
Recession analysis signatures.

Analyzes streamflow recession behavior using:
- Power-law fitting: dQ/dt = a * Q^b
- Point cloud and event-based methods
- Sinusoidal seasonality of recession parameters
"""

from typing import Dict, List, Tuple, Optional
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
    Identify recession events in a streamflow time series.

    A recession event is a period of monotonically decreasing flow
    with decreasing rate of change (|dQ/dt| also decreasing).

    Parameters
    ----------
    Q : np.ndarray
        Daily discharge values
    min_length : int, default 5
        Minimum number of consecutive days for a recession event

    Returns
    -------
    list of dict
        Each dict contains 'start', 'end', and 'indices' keys
    """
    n = len(Q)
    if n < min_length + 1:
        return []

    # Calculate dQ/dt (forward difference)
    dQ_dt = np.diff(Q)
    dQ_dt = np.append(dQ_dt, np.nan)

    recession_events = []
    in_recession = False
    start_idx = None

    for i in range(n - min_length):
        # Check if we have valid data
        if np.isnan(Q[i]) or np.isnan(dQ_dt[i]):
            if in_recession:
                # End current recession if we hit NA
                if start_idx is not None and (i - start_idx) >= min_length:
                    recession_events.append({
                        'start': start_idx,
                        'end': i - 1,
                        'indices': list(range(start_idx, i))
                    })
                in_recession = False
                start_idx = None
            continue

        # Check for monotonic decrease in both Q and |dQ/dt|
        if i < n - 1:
            is_recession = True

            # Check next min_length consecutive decreases
            for j in range(min_length):
                idx = i + j
                if idx + 1 >= n:
                    is_recession = False
                    break
                if (np.isnan(Q[idx]) or np.isnan(Q[idx + 1]) or
                    np.isnan(dQ_dt[idx]) or np.isnan(dQ_dt[idx + 1])):
                    is_recession = False
                    break

                # Check if Q is decreasing
                if Q[idx + 1] >= Q[idx]:
                    is_recession = False
                    break

                # Check if |dQ/dt| is decreasing (becoming less negative)
                if abs(dQ_dt[idx + 1]) >= abs(dQ_dt[idx]):
                    is_recession = False
                    break

            if is_recession and not in_recession:
                # Start new recession
                in_recession = True
                start_idx = i
            elif not is_recession and in_recession:
                # End current recession
                if start_idx is not None and (i - start_idx) >= min_length:
                    recession_events.append({
                        'start': start_idx,
                        'end': i - 1,
                        'indices': list(range(start_idx, i))
                    })
                in_recession = False
                start_idx = None

    # Check if we ended in a recession
    if in_recession and start_idx is not None and (n - start_idx) >= min_length:
        recession_events.append({
            'start': start_idx,
            'end': n - 1,
            'indices': list(range(start_idx, n))
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
                             "b_events", "concavity"]
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
    })

    # Store all recession events with timing
    all_recession_events = []

    # Process each year (df already sorted by water_year and dowy)
    for yr, year_data in df.groupby("water_year", sort=False):

        Q = year_data["Q"].values
        dowy = year_data["dowy"].values

        # Identify recession events
        recession_events = identify_recession_events(Q)

        # Collect all recession data for point cloud
        all_Q = []
        all_dQ_dt = []

        # Store individual event parameters
        event_log_a_values = []
        event_b_values = []
        event_concavities = []
        successful_events_Q = []  # Store Q data for successful events (for log_a recalc)

        # Process each recession event
        for event in recession_events:
            Q_event = Q[event['indices']]

            # Get middle day of recession for timing
            mid_idx = event['indices'][len(event['indices']) // 2]
            event_dowy = dowy[mid_idx]

            # Fit parameters for this event
            log_a, b = fit_recession_event(Q_event, remove_first_day=True)

            if not np.isnan(log_a) and not np.isnan(b):
                event_log_a_values.append(log_a)
                event_b_values.append(b)
                successful_events_Q.append(Q_event.copy())

                # Store event with timing
                all_recession_events.append({
                    'water_year': yr,
                    'dowy': event_dowy,
                    'log_a': log_a,
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
            try:
                all_Q = np.array(all_Q)
                all_dQ_dt = np.array(all_dQ_dt)
                log_Q = np.log(all_Q)
                log_dQ_dt = np.log(all_dQ_dt)

                # Skip near-singular data (matching R's lm() QR rank check)
                # R's tolerance is .Machine$double.eps^0.5 ≈ 1.49e-8
                if np.var(log_Q) < 1e-8:
                    pass  # Year remains NA
                else:
                    result = scipy_stats.linregress(log_Q, log_dQ_dt)
                    b_pointcloud = result.slope

                    # Calculate log(a) using b_pointcloud
                    log_a_values_pc = log_dQ_dt - b_pointcloud * log_Q
                    log_a_pointcloud = np.median(log_a_values_pc)

                    annual_metrics.loc[annual_metrics["water_year"] == yr, "b_pointcloud"] = b_pointcloud
                    annual_metrics.loc[annual_metrics["water_year"] == yr, "log_a_pointcloud"] = log_a_pointcloud
            except Exception:
                pass

        # Individual events analysis — recalculate log_a using median_b (matching R)
        if len(event_b_values) > 0:
            median_b = np.median(event_b_values)

            # Recalculate log_a for each event using the ensemble median_b
            log_a_events_recalc = []
            for Q_ev in successful_events_Q:
                if len(Q_ev) > 1:
                    Q_subset = Q_ev[1:]  # Remove first day
                    dQ_subset = -np.diff(Q_ev[1:])
                    Q_for_calc = Q_subset[:-1]
                    valid_mask = (Q_for_calc > 0) & (dQ_subset > 0)
                    if valid_mask.any():
                        log_a_vals = np.log(dQ_subset[valid_mask]) - median_b * np.log(Q_for_calc[valid_mask])
                        log_a_events_recalc.append(np.median(log_a_vals))

            if log_a_events_recalc:
                annual_metrics.loc[annual_metrics["water_year"] == yr, "log_a_events"] = np.median(log_a_events_recalc)
            annual_metrics.loc[annual_metrics["water_year"] == yr, "b_events"] = np.median(event_b_values)

        # Concavity
        if len(event_concavities) > 0:
            annual_metrics.loc[annual_metrics["water_year"] == yr, "concavity"] = np.mean(event_concavities)

    # Check minimum events requirement
    if len(all_recession_events) < min_events:
        # Return all NAs
        result = {}
        for sig in signatures_with_stats:
            for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                           "_mk_rho", "_mk_pval", "_mean", "_median"]:
                result[f"{sig}{suffix}"] = np.nan
        for sig in seasonality_signatures:
            result[sig] = np.nan
        return result

    # Generate statistics for main signatures
    result = generate_stats(
        annual_metrics,
        value_cols=signatures_with_stats,
        year_col="water_year"
    )

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

    return result
