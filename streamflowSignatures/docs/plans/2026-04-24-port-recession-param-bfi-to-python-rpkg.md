# Port Recession-Parameterized BFI to Python & rpkg

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port 25 new recession-parameterized BFI columns (alpha_linear×8, BFI_Eckhardt_param×8, BFI_LyneHollick_param×8, recession_alpha_point_cloud_linear_reservoir×1) from Julia canonical to Python and rpkg, then validate against Julia golden output.

**Architecture:** Three code changes per language: (1) add alpha_linear computation + whole-record scalar to recession module, (2) add `analyze_baseflow_indices_with_parameters()` to baseflow module, (3) wire orchestration + exports. Then run benchmarks and compare vs Julia golden (624 cols).

**Tech Stack:** Python (numpy, pandas, scipy), R (rpkg package), Julia golden output for validation.

**Julia reference files (canonical):**
- `julia/src/recession.jl:342-478` — alpha_linear per-year + whole-record scalar
- `julia/src/baseflow.jl:267-345` — `analyze_baseflow_indices_with_parameters()`
- `julia/src/signatures.jl:80-98` — orchestration wiring

---

### Task 1: Add alpha_linear + scalar to Python recession.py

**Files:**
- Modify: `python/streamflow_signatures/recession.py:270-484`

**Context:** The Python recession function (`analyze_recession_parameters`) processes events per year. We need to add Q_{i+1}/Q_i ratio computation inside the event loop (after `fit_recession_event`), store per-year and whole-record collections, and output both per-year stats and the whole-record scalar.

**Step 1: Add alpha_linear to annual_metrics initialization**

In `recession.py`, find line 273 (the `annual_metrics` DataFrame initialization). Add `"alpha_linear": np.nan` to the dict:

```python
    annual_metrics = pd.DataFrame({
        "water_year": years,
        "log_a_pointcloud": np.nan,
        "log_a_events": np.nan,
        "b_pointcloud": np.nan,
        "b_events": np.nan,
        "concavity": np.nan,
        "n_recession_events": np.nan,
        "alpha_linear": np.nan,          # NEW
    })
```

**Step 2: Add whole-record alpha collection vector**

After `all_recession_events = []` (line 284), add:

```python
    all_alpha_linear = []  # Whole-record Q_{i+1}/Q_i ratios for scalar
```

**Step 3: Add alpha computation inside the event loop**

Inside the `for event in recession_events:` loop, right after `log_a, b = fit_recession_event(Q_event, remove_first_day=True)` (line 318), BEFORE the `if not np.isnan(log_a)` check (line 320), add the alpha computation. This is INDEPENDENT of fit success:

```python
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
```

**Step 4: Initialize per-year alpha list at top of year loop**

At the start of the `for yr, year_data in df.groupby(...)` loop (around line 287), add `year_alpha_linear = []` alongside the existing `all_Q = []` etc. Place it right after `all_dQ_dt = []` (line 301):

```python
        year_alpha_linear = []
```

**Step 5: Store per-year alpha_linear median**

After the concavity block (after line 405), add:

```python
        # Per-year alpha_linear (linear reservoir assumption)
        if len(year_alpha_linear) > 10:
            annual_metrics.loc[annual_metrics["water_year"] == yr, "alpha_linear"] = np.median(year_alpha_linear)
```

**Step 6: Add alpha_linear to signatures_with_stats**

Find line ~258 where `signatures_with_stats` is defined. Add `"alpha_linear"`:

```python
    signatures_with_stats = ["log_a_pointcloud", "log_a_events",
                             "b_pointcloud", "b_events", "concavity",
                             "alpha_linear"]
```

**Step 7: Add whole-record scalar output**

Right before the min_events check (line ~420 `if len(all_recession_events) < min_events:`), after `n_events_stats`, add:

```python
    # Whole-record recession alpha (scalar — independent of min_events gate)
    recession_alpha_scalar = np.median(all_alpha_linear) if len(all_alpha_linear) > 10 else np.nan
```

Then in BOTH branches (the `< min_events` NA-return path around line 431, AND the normal return at the end), make sure the scalar is added to the result:

In the `< min_events` branch, add before `return result`:
```python
        result["recession_alpha_point_cloud_linear_reservoir"] = recession_alpha_scalar
```

At the very end of the function (after all seasonality), add before `return result`:
```python
    result["recession_alpha_point_cloud_linear_reservoir"] = recession_alpha_scalar
```

**Step 8: Add alpha_linear to NA template**

In the `< min_events` branch (lines 423-431), the loop `for sig in signatures_with_stats` already covers all stats. Since we added `"alpha_linear"` to `signatures_with_stats`, this is automatic. Verify.

**Step 9: Verify the function returns the new columns**

Run: `python -c "from streamflow_signatures.recession import analyze_recession_parameters; print('import ok')"`
Expected: `import ok`

**Step 10: Commit**

```bash
git add python/streamflow_signatures/recession.py
git commit -m "feat(python): add alpha_linear + recession_alpha scalar to recession analysis"
```

---

### Task 2: Add analyze_baseflow_indices_with_parameters() to Python baseflow.py

**Files:**
- Modify: `python/streamflow_signatures/baseflow.py` (append after line 223)

**Context:** This function mirrors `analyze_baseflow_indices()` (lines 136-223) but takes a recession-derived `alpha` parameter instead of using defaults. The existing `eckhardt_filter()` and `lyne_hollick_filter()` already accept custom parameters.

**Step 1: Add the new function**

Append after the existing `analyze_baseflow_indices()` (after line 223):

```python
def analyze_baseflow_indices_with_parameters(
    streamflow_data: pd.DataFrame,
    alpha: float,
    BFImax: float = ECKHARDT_BFIMAX,
    passes: int = LYNE_HOLLICK_PASSES,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calculate baseflow indices using recession-derived filter parameters.

    Same as analyze_baseflow_indices() but uses recession-derived alpha
    instead of fixed defaults. BFImax remains fixed at 0.8.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns: water_year, Q, dowy.
    alpha : float
        Recession-derived filter parameter (from recession_alpha_point_cloud_linear_reservoir).
        Must be in (0, 1). NaN or out-of-range returns all NaN stats.
    BFImax : float, default 0.8
        Maximum baseflow index for Eckhardt filter.
    passes : int, default 2
        Number of passes for Lyne-Hollick filter.

    Returns
    -------
    dict
        Dictionary with BFI_Eckhardt_param_* and BFI_LyneHollick_param_* statistics.
    """
    metrics = ["BFI_Eckhardt_param", "BFI_LyneHollick_param"]

    # Validate alpha
    if np.isnan(alpha) or alpha <= 0.0 or alpha >= 1.0:
        result = {}
        for m in metrics:
            for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                           "_mk_rho", "_mk_pval", "_mean", "_median"]:
                result[f"{m}{suffix}"] = np.nan
        return result

    # Validate required columns
    required_cols = ["water_year", "Q", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data
    years = df["water_year"].unique()

    bfi_by_year = pd.DataFrame({
        "water_year": years,
        "BFI_Eckhardt_param": np.nan,
        "BFI_LyneHollick_param": np.nan,
    })

    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()
        year_data = year_data.sort_values("dowy")
        Q = year_data["Q"].values.astype(float)

        # Apply filters with recession-derived alpha
        bf_eck = eckhardt_filter(Q, BFImax=BFImax, a=alpha)
        bf_lh = lyne_hollick_filter(Q, alpha=alpha, passes=passes)

        # Eckhardt BFI with paired masking
        valid_eck = ~np.isnan(Q) & ~np.isnan(bf_eck)
        if valid_eck.sum() > 0 and np.sum(Q[valid_eck]) > 0:
            bfi_eck = np.sum(bf_eck[valid_eck]) / np.sum(Q[valid_eck])
            bfi_eck = np.clip(bfi_eck, 0.0, 1.0)
            bfi_by_year.loc[bfi_by_year["water_year"] == yr, "BFI_Eckhardt_param"] = bfi_eck

        # Lyne-Hollick BFI with paired masking
        valid_lh = ~np.isnan(Q) & ~np.isnan(bf_lh)
        if valid_lh.sum() > 0 and np.sum(Q[valid_lh]) > 0:
            bfi_lh = np.sum(bf_lh[valid_lh]) / np.sum(Q[valid_lh])
            bfi_lh = np.clip(bfi_lh, 0.0, 1.0)
            bfi_by_year.loc[bfi_by_year["water_year"] == yr, "BFI_LyneHollick_param"] = bfi_lh

    return generate_stats(
        bfi_by_year,
        value_cols=metrics,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness,
    )
```

**Step 2: Verify import works**

Run: `python -c "from streamflow_signatures.baseflow import analyze_baseflow_indices_with_parameters; print('import ok')"`
Expected: `import ok`

**Step 3: Commit**

```bash
git add python/streamflow_signatures/baseflow.py
git commit -m "feat(python): add analyze_baseflow_indices_with_parameters()"
```

---

### Task 3: Wire Python orchestration + exports

**Files:**
- Modify: `python/streamflow_signatures/signatures.py:13-14,115-119`
- Modify: `python/streamflow_signatures/__init__.py:27,55-58`

**Step 1: Add import in signatures.py**

Find line 13: `from .baseflow import analyze_baseflow_indices`
Change to:
```python
from .baseflow import analyze_baseflow_indices, analyze_baseflow_indices_with_parameters
```

**Step 2: Wire recession alpha extraction + parameterized BFI call**

Replace the recession block (lines 115-119):
```python
    # Recession: inherently sparse (event-based), no trend completeness
    try:
        results.update(analyze_recession_parameters(gage_data))
    except Exception:
        pass
```

With:
```python
    # Recession: inherently sparse (event-based), no trend completeness
    recession_alpha = float('nan')
    try:
        recession_results = analyze_recession_parameters(gage_data)
        # Extract scalar for parameterized BFI before merging
        recession_alpha = recession_results.get(
            "recession_alpha_point_cloud_linear_reservoir", float('nan'))
        results.update(recession_results)
    except Exception:
        pass

    # Parameterized BFI using recession-derived alpha (requires recession to have run first)
    # Uses trend completeness (same as fixed-parameter BFI)
    try:
        results.update(analyze_baseflow_indices_with_parameters(
            gage_data, recession_alpha, **trend_kwargs))
    except Exception:
        pass
```

**Step 3: Add export in __init__.py**

Find line 27: `from .baseflow import analyze_baseflow_indices, eckhardt_filter, lyne_hollick_filter`
Change to:
```python
from .baseflow import analyze_baseflow_indices, analyze_baseflow_indices_with_parameters, eckhardt_filter, lyne_hollick_filter
```

Add `"analyze_baseflow_indices_with_parameters"` to `__all__` after `"analyze_baseflow_indices"` (around line 55):
```python
    "analyze_baseflow_indices",
    "analyze_baseflow_indices_with_parameters",
```

**Step 4: Verify import chain**

Run: `python -c "from streamflow_signatures import calculate_all_signatures; print('ok')"`
Expected: `ok`

**Step 5: Commit**

```bash
git add python/streamflow_signatures/signatures.py python/streamflow_signatures/__init__.py
git commit -m "feat(python): wire recession-parameterized BFI into orchestration + exports"
```

---

### Task 4: Add alpha_linear + scalar to rpkg recession.R

**Files:**
- Modify: `rpkg/R/recession.R:109-302`

**Context:** Mirrors Task 1 (Python) but in R. The rpkg recession function processes events per year in a `for (yr in years)` loop. We add the same Q_{i+1}/Q_i computation.

**Step 1: Add alpha_linear to annual_metrics initialization**

Find line 136 (the `annual_metrics` data.frame). Add `alpha_linear`:

```r
  annual_metrics <- data.frame(water_year = years, log_a_pointcloud = NA_real_,
                               log_a_events = NA_real_, b_pointcloud = NA_real_,
                               b_events = NA_real_, concavity = NA_real_,
                               n_recession_events = NA_real_,
                               alpha_linear = NA_real_)
```

**Step 2: Add whole-record alpha collection vector**

After `all_recession_events <- list()` (line 143), add:

```r
  all_alpha_linear <- numeric(0)
```

**Step 3: Add per-year alpha list at top of year loop**

Inside the `for (yr in years)` loop (line 145), after `event_concavities <- numeric(0)` (line 158), add:

```r
    year_alpha_linear <- numeric(0)
```

**Step 4: Add alpha computation inside the event loop**

Inside the `for (event in recession_events)` loop, right after `params <- fit_recession_event(Q_event, remove_first_day = TRUE)` (line 165), BEFORE the `if (!is.na(params$log_a)` check (line 167), add:

```r
      # Compute discrete recession constant alpha = Q_{i+1}/Q_i (b=1 linear reservoir)
      # INDEPENDENT of power-law fit — depends only on raw Q pairs
      # Remove first day (storm peak) consistent with point-cloud fitting
      if (length(Q_event) > 2) {
        Q_alpha <- Q_event[-1]  # Remove first day
        for (j in seq_len(length(Q_alpha) - 1)) {
          if (Q_alpha[j] > 0 && !is.na(Q_alpha[j]) && !is.na(Q_alpha[j + 1])) {
            a_i <- Q_alpha[j + 1] / Q_alpha[j]
            if (a_i > 0 && a_i < 1) {
              year_alpha_linear <- c(year_alpha_linear, a_i)
              all_alpha_linear <- c(all_alpha_linear, a_i)
            }
          }
        }
      }
```

**Step 5: Store per-year alpha_linear median**

After the concavity block (after line 241 `annual_metrics$concavity[...] <- mean(event_concavities...)`), add:

```r
    # Per-year alpha_linear (linear reservoir assumption)
    if (length(year_alpha_linear) > 10) {
      annual_metrics$alpha_linear[annual_metrics$water_year == yr] <- median(year_alpha_linear, na.rm = TRUE)
    }
```

**Step 6: Add alpha_linear to signatures_with_stats**

Find line 118 where `signatures_with_stats` is defined. Add `"alpha_linear"`:

```r
  signatures_with_stats <- c("log_a_pointcloud", "log_a_events",
                             "b_pointcloud", "b_events", "concavity",
                             "alpha_linear")
```

**Step 7: Add whole-record scalar + update NA template**

Right before the min_events check (line ~256, `if (length(all_recession_events) < min_events)`), after `n_events_stats`, add:

```r
  # Whole-record recession alpha (scalar — independent of min_events gate)
  recession_alpha_scalar <- if (length(all_alpha_linear) > 10) {
    median(all_alpha_linear, na.rm = TRUE)
  } else {
    NA_real_
  }
```

Update `make_na_result()` (line 128) to include the scalar:

```r
  make_na_result <- function() {
    out <- unlist(lapply(signatures_with_stats, empty_stats), recursive = FALSE)
    out <- c(out, unlist(empty_stats("n_recession_events"), recursive = FALSE))
    for (s in seasonality_sigs) out[[s]] <- NA_real_
    out[["recession_alpha_point_cloud_linear_reservoir"]] <- NA_real_
    out
  }
```

In the `< min_events` NA-return path (around line 258), add before `return(na_result)`:
```r
    na_result[["recession_alpha_point_cloud_linear_reservoir"]] <- recession_alpha_scalar
```

At the very end of the function (line 301, before `result`), add:
```r
  result[["recession_alpha_point_cloud_linear_reservoir"]] <- recession_alpha_scalar
```

**Step 8: Commit**

```bash
git add rpkg/R/recession.R
git commit -m "feat(rpkg): add alpha_linear + recession_alpha scalar to recession analysis"
```

---

### Task 5: Add analyze_baseflow_indices_with_parameters() to rpkg baseflow.R

**Files:**
- Modify: `rpkg/R/baseflow.R` (append after line 132)

**Context:** Mirrors Task 2 (Python) but in R. Uses existing `eckhardt_filter()` and `lyne_hollick_filter()` from the same file.

**Step 1: Add the new function**

Append after `analyze_baseflow_indices()` (after line 132):

```r
#' Analyze baseflow indices with recession-derived parameters
#'
#' Same as \code{analyze_baseflow_indices()} but uses recession-derived alpha
#' instead of fixed defaults. BFImax remains fixed at 0.8.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @param alpha Recession-derived filter parameter (from
#'   recession_alpha_point_cloud_linear_reservoir). Must be in (0, 1).
#' @param BFImax Maximum baseflow index for Eckhardt filter (default from config).
#' @param passes Number of Lyne-Hollick filter passes (default from config).
#' @param trend_completeness Forwarded to generate_stats().
#' @param decade_completeness Forwarded to generate_stats().
#' @return Named list of 8 statistics per BFI metric (2 x 8 = 16 values).
#' @export
analyze_baseflow_indices_with_parameters <- function(
    streamflow_data,
    alpha,
    BFImax = pkg_env$eckhardt_bfimax,
    passes = pkg_env$lyne_hollick_passes,
    trend_completeness = NULL,
    decade_completeness = NULL
) {
  metrics <- c("BFI_Eckhardt_param", "BFI_LyneHollick_param")

  # Validate alpha
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
    result <- unlist(lapply(metrics, empty_stats), recursive = FALSE)
    return(result)
  }

  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  bfi_by_year <- data.frame(water_year = years, BFI_Eckhardt_param = NA_real_,
                            BFI_LyneHollick_param = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    year_data <- year_data[order(year_data$dowy), ]
    Q <- year_data$Q

    # Apply filters with recession-derived alpha
    bf_eck <- eckhardt_filter(Q, BFImax = BFImax, a = alpha)
    bf_lh <- lyne_hollick_filter(Q, alpha = alpha, passes = passes)

    # Eckhardt: paired masking
    eck_valid <- !is.na(Q) & !is.na(bf_eck)
    total_eck <- sum(Q[eck_valid])
    if (total_eck > 0) {
      bfi_eck <- sum(bf_eck[eck_valid]) / total_eck
      bfi_eck <- max(0, min(1, bfi_eck))
      bfi_by_year$BFI_Eckhardt_param[bfi_by_year$water_year == yr] <- bfi_eck
    }

    # Lyne-Hollick: paired masking
    lh_valid <- !is.na(Q) & !is.na(bf_lh)
    total_lh <- sum(Q[lh_valid])
    if (total_lh > 0) {
      bfi_lh <- sum(bf_lh[lh_valid]) / total_lh
      bfi_lh <- max(0, min(1, bfi_lh))
      bfi_by_year$BFI_LyneHollick_param[bfi_by_year$water_year == yr] <- bfi_lh
    }
  }

  generate_stats(bfi_by_year, value_cols = metrics, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness)
}
```

**Step 2: Commit**

```bash
git add rpkg/R/baseflow.R
git commit -m "feat(rpkg): add analyze_baseflow_indices_with_parameters()"
```

---

### Task 6: Wire rpkg orchestration + exports

**Files:**
- Modify: `rpkg/R/signatures.R:80-83`
- Modify: `rpkg/NAMESPACE`

**Step 1: Wire recession alpha extraction + parameterized BFI call in signatures.R**

Replace lines 80-83:
```r
  # Recession: event-based (inherently sparse), no trend completeness
  out <- safe_call(analyze_recession_parameters, "Recession parameters",
                   streamflow_data)
  if (!is.null(out)) results <- c(results, out)
```

With:
```r
  # Recession: event-based (inherently sparse), no trend completeness
  recession_alpha <- NaN
  out <- safe_call(analyze_recession_parameters, "Recession parameters",
                   streamflow_data)
  if (!is.null(out)) {
    # Extract scalar for parameterized BFI before merging
    if ("recession_alpha_point_cloud_linear_reservoir" %in% names(out)) {
      recession_alpha <- out[["recession_alpha_point_cloud_linear_reservoir"]]
    }
    results <- c(results, out)
  }

  # Parameterized BFI using recession-derived alpha (requires recession to have run first)
  # Uses trend completeness (same as fixed-parameter BFI)
  out <- safe_call(analyze_baseflow_indices_with_parameters, "Parameterized baseflow",
                   streamflow_data, alpha = recession_alpha,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness)
  if (!is.null(out)) results <- c(results, out)
```

**Step 2: Add export in NAMESPACE**

Add after `export(analyze_baseflow_indices)`:
```
export(analyze_baseflow_indices_with_parameters)
```

**Step 3: Commit**

```bash
git add rpkg/R/signatures.R rpkg/NAMESPACE
git commit -m "feat(rpkg): wire recession-parameterized BFI into orchestration + exports"
```

---

### Task 7: Run Python benchmark and compare vs Julia golden

**Files:**
- Run: `python docs/benchmarks/run_python_benchmark.py` (~70-130 min)
- Run: `python docs/benchmarks/compare_python_vs_golden_julia.py`
- Run: `python docs/benchmarks/build_julia_golden_dashboard.py python`

**Step 1: Run Python benchmark**

```bash
cd C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures
python docs/benchmarks/run_python_benchmark.py
```

Expected: 7,313 gages, 624 signature columns, ~70-130 min.

**Step 2: Run comparison vs Julia golden**

```bash
python docs/benchmarks/compare_python_vs_golden_julia.py
```

Expected: 624 common columns, ~619+ Perfect (R² >= 0.999), ~3-5 Good, <=3 Poor (same irreducible library differences as before).

**Step 3: Build dashboard**

```bash
python docs/benchmarks/build_julia_golden_dashboard.py python
```

**Step 4: Review results**

Check that the 25 new columns (alpha_linear×8, BFI_Eckhardt_param×8, BFI_LyneHollick_param×8, recession_alpha scalar) all show R² >= 0.999 vs Julia golden.

**Step 5: Commit benchmark results**

```bash
git add docs/benchmarks/python_vs_golden_julia_*.md docs/benchmarks/python_vs_golden_julia_*.csv
git commit -m "benchmark: Python vs Julia golden — 624 columns validated"
```

---

### Task 8: Run rpkg benchmark and compare vs Julia golden

**Files:**
- Run: `Rscript docs/benchmarks/run_rpkg_benchmark.R` (~2-3 hours)
- Run: `python docs/benchmarks/compare_rpkg_vs_golden_julia.py`
- Run: `python docs/benchmarks/build_julia_golden_dashboard.py rpkg`

**Step 1: Run rpkg benchmark**

```bash
cd C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures
Rscript docs/benchmarks/run_rpkg_benchmark.R
```

Expected: 7,313 gages, 624 signature columns, ~2-3 hours.

**Step 2: Run comparison vs Julia golden**

```bash
python docs/benchmarks/compare_rpkg_vs_golden_julia.py
```

Expected: 624 common columns, similar alignment profile to Python.

**Step 3: Build dashboard**

```bash
python docs/benchmarks/build_julia_golden_dashboard.py rpkg
```

**Step 4: Review results**

Check that the 25 new columns all show R² >= 0.999 vs Julia golden.

**Step 5: Commit**

```bash
git add docs/benchmarks/rpkg_vs_golden_julia_*.md docs/benchmarks/rpkg_vs_golden_julia_*.csv
git commit -m "benchmark: rpkg vs Julia golden — 624 columns validated"
```

---

### Task 9: Update documentation and CHANGELOG

**Files:**
- Modify: `CHANGELOG.md` — Update "Recession-Parameterized Baseflow Signatures" section to reflect Python/rpkg port completion
- Modify: `docs/CROSS_LANGUAGE_STATUS.md` — Update column counts and benchmark results
- Modify: `docs/DEVELOPMENT.md` — Update benchmark timing table with new results

**Step 1: Update CHANGELOG.md**

In the "Recession-Parameterized Baseflow Signatures" section, change:
```
Three new signatures using recession-derived filter parameters (Collischonn & Fan 2013, Eckhardt 2005). Julia only (canonical); Python and rpkg ports deferred.
```
To:
```
Three new signatures using recession-derived filter parameters (Collischonn & Fan 2013, Eckhardt 2005). Implemented in all three languages (Julia canonical, Python, rpkg).
```

Add Python/rpkg files to the "Files modified" list.

**Step 2: Update CROSS_LANGUAGE_STATUS.md**

Update the Python description to remove "port pending" note. Update column counts to 624 for all three implementations.

**Step 3: Update Planned section**

Remove "Port recession-parameterized BFI to Python and rpkg" from the Planned list.

**Step 4: Commit**

```bash
git add CHANGELOG.md docs/CROSS_LANGUAGE_STATUS.md docs/DEVELOPMENT.md
git commit -m "docs: update for Python/rpkg recession-parameterized BFI port completion"
```
