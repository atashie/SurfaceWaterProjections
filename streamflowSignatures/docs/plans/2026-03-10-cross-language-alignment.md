# Cross-Language Alignment: Fix R/Python/Julia Signature Divergences

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Bring Python and Julia streamflow signature outputs into >0.99 Spearman correlation with the canonical R implementation on all signature columns.

**Architecture:** Nine issues ranked by severity across three tiers. Each fix targets a specific root cause identified by the three-way comparison (`benchmarks/compare_three_way.py`). After each fix, re-run the comparison on the affected columns to verify improvement.

**Tech Stack:** Julia 1.9+ (julia/src/), Python 3.12+ (python/streamflow_signatures/), R (R/helperFunctions.R is canonical — do NOT modify unless noted). Shared config: `config/signatures_config.json`. Benchmark outputs in `benchmarks/`.

**Validation command** (run after each tier):
```bash
py benchmarks/compare_three_way.py
```

---

## Tier 1: Clear Bugs (Issues 1-3)

These are definitive bugs with known root causes and clear fixes.

---

### Task 1: Fix Julia D25_to_D75 — All NaN

**Root Cause:** `julia/src/timing.jl` lines 104-105 look up `D25_day` and `D75_day` from the row dictionary, but `CFG_D_PERCENTILES = [5,10,20,30,40,50,60,70,80,90,95]` does NOT include 25 or 75. The keys never exist, so `get()` always returns NaN.

**How R/Python solve it:** They compute D25/D75 directly from cumulative percentages (see `python/streamflow_signatures/timing.py` lines 102-109), independent of the D_PERCENTILES list.

**Files:**
- Modify: `julia/src/timing.jl:103-110`
- Test: `julia/test/runtests.jl` (Timing test set)

**Step 1: Write a failing test**

Add to `julia/test/runtests.jl` inside the `@testset "Timing"` block:

```julia
@testset "D25_to_D75 computed" begin
    result = analyze_flow_timing_trends(df)
    # D25_to_D75 must NOT be NaN — it should be a positive number
    @test haskey(result, "D25_to_D75_mean")
    @test !isnan(result["D25_to_D75_mean"])
    @test result["D25_to_D75_mean"] > 0
end
```

**Step 2: Run test to verify it fails**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `D25_to_D75_mean` is NaN.

**Step 3: Fix `timing.jl` — compute D25/D75 directly from cumulative percentages**

Replace lines 103-110 in `julia/src/timing.jl`:

```julia
        # D25 to D75 duration — compute directly from cumulative percentages
        # (NOT from D_PERCENTILES, which may not include 25 and 75)
        idx_25 = findfirst(cum_pct .>= 25)
        idx_75 = findfirst(cum_pct .>= 75)
        if idx_25 !== nothing && idx_75 !== nothing
            row["D25_to_D75"] = dowy_sorted[idx_75] - dowy_sorted[idx_25]
        else
            row["D25_to_D75"] = NaN
        end
```

**Step 4: Run test to verify it passes**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS — D25_to_D75_mean is a positive number.

**Step 5: Commit**

```bash
git add julia/src/timing.jl julia/test/runtests.jl
git commit -m "fix(julia): compute D25_to_D75 directly from cumulative percentages

D25 and D75 are not in CFG_D_PERCENTILES, so the dict lookup always
returned NaN. Now computed directly like R/Python implementations."
```

---

### Task 2: Fix Julia log_a_seasonality_minimum — 365-inverted

**Root Cause:** Julia's `fit_recession_seasonality()` (line 231) uses a different formula than R/Python to convert sinusoidal phase to minimum day.

- **R/Python:** Basis order `[sin, cos, 1]`, phase = `atan2(-cos_coeff, sin_coeff)`, minimum = `phase_days + 273.75`
- **Julia:** Basis order `[1, cos, sin]`, phase = `atan(sin_coeff, cos_coeff)`, minimum = `mod((-phase + pi)/omega, 365.25)`

Result: Julia computes the complement — `Julia ≈ 365 - R`. Confirmed by data: R_mean=157, Julia_mean=207, (R + Julia)/2 ≈ 182.

**Files:**
- Modify: `julia/src/recession.jl:211-231`
- Test: `julia/test/runtests.jl` (Recession test set)

**Step 1: Write a failing test**

Add to `julia/test/runtests.jl` inside the `@testset "Recession"` block:

```julia
@testset "Seasonality minimum matches R formula" begin
    # Known test case: sin curve with minimum at day 100
    # Model: log_a = 2.0 * cos(2pi/365 * (dowy - 100))
    # Minimum at dowy=100 + 182.5 = 282.5 (opposite phase)
    # Actually: minimum of cos(x - phi) is at x = phi + pi
    # So if we set up a signal with known minimum, we can verify the formula

    n = 1000
    dowy = rand(1:365, n)
    omega = 2 * pi / 365
    # Signal with minimum at day ~200
    log_a = -0.5 .* cos.(omega .* (dowy .- 200)) .+ randn(n) .* 0.01

    result = StreamflowSignatures.fit_recession_seasonality(
        Float64.(log_a), Int.(dowy)
    )

    # The minimum day should be near 200, NOT near 365-200=165
    @test abs(result["minimum_day"] - 200) < 30  # within 30 days
end
```

**Step 2: Run test to verify it fails**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — minimum_day is ~165 instead of ~200.

**Step 3: Fix `recession.jl` — match R/Python formula exactly**

Replace lines 211-231 in `julia/src/recession.jl`:

```julia
    # Convert dowy to radians (match R's 365, not 365.25)
    omega = 2 * pi / 365
    theta = omega .* dowy

    # Fit: log_a = B1*sin(theta) + B2*cos(theta) + C
    # Using least squares (same basis order as R/Python: sin, cos, intercept)
    n = length(log_a)
    X = hcat(sin.(theta), cos.(theta), ones(n))
    y = log_a

    # Solve normal equations
    coeffs = X \ y

    B1, B2 = coeffs[1], coeffs[2]  # sin coeff, cos coeff

    # Calculate amplitude
    amplitude = sqrt(B1^2 + B2^2)

    # Calculate phase (matching R's atan2(-B2, B1))
    phase_rad = atan(-B2, B1)
    phase_days = phase_rad * 365 / (2 * pi)

    # Ensure phase is between 0 and 365
    if phase_days < 0
        phase_days += 365
    end

    # Minimum occurs at phase + 273.75 days (3/4 of a cycle, matching R)
    minimum_day = phase_days + 273.75
    if minimum_day > 365
        minimum_day -= 365
    end

    return Dict("amplitude" => amplitude, "minimum_day" => minimum_day)
```

**Step 4: Run test to verify it passes**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS — minimum_day is near 200.

**Step 5: Commit**

```bash
git add julia/src/recession.jl julia/test/runtests.jl
git commit -m "fix(julia): match R formula for recession seasonality minimum day

Was using cos/sin basis with different phase formula, producing
365-complement of correct value. Now uses sin/cos basis and
atan2(-B2, B1) + 273.75 formula matching R/Python exactly."
```

---

### Task 3: Diagnose and Fix Python BFI_LyneHollick for Canadian Gages

**Root Cause:** Unknown — requires diagnosis. All top 20 outliers are Canadian gages (05CB007, 06GD001, etc.). Python BFI is 0.15-0.59 lower than R/Julia (which agree at rho=0.9983). The filter code appears identical across all three languages.

**Hypothesis:** The difference may stem from data preprocessing (water year assignment, dowy computation, or sort order) for Canadian gages in the Python benchmark, not the filter algorithm itself.

**Files:**
- Create: `benchmarks/diagnose_bfi_canadian.py`
- Possibly modify: `python/streamflow_signatures/baseflow.py` or `benchmarks/run_python_benchmark.py`
- Test: The diagnostic script IS the test

**Step 1: Write diagnostic script**

Create `benchmarks/diagnose_bfi_canadian.py`:

```python
"""
Diagnose BFI_LyneHollick discrepancy for Canadian gages.

Compares Python filter output against R/Julia for gage 05CB007.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
from streamflow_signatures.baseflow import lyne_hollick_filter, eckhardt_filter
from streamflow_signatures.config import (
    LYNE_HOLLICK_ALPHA, LYNE_HOLLICK_PASSES,
    BASEFLOW_MIN_DAYS, BASEFLOW_MAX_MISSING_FRAC,
    WATER_YEAR_START_MONTH,
)
from streamflow_signatures.io import add_water_year_columns

GAGE = "05CB007"  # Worst outlier: R=0.8796, Py=0.2925, Jl=0.8774

# Load raw data
parquet_path = Path(r"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet")
df = pd.read_parquet(parquet_path)
df.columns = [c.lower() for c in df.columns]
if "date" not in df.columns and "Date" in df.columns:
    df.rename(columns={"Date": "date"}, inplace=True)

# Filter to gage
gage_df = df[df["gage_id"].astype(str) == GAGE].copy()
print(f"Gage {GAGE}: {len(gage_df)} rows, date range {gage_df['date'].min()} to {gage_df['date'].max()}")

# Add water year columns
gage_df = add_water_year_columns(gage_df)
print(f"Water years: {sorted(gage_df['water_year'].unique())}")

# Check data quality per year
for yr in sorted(gage_df["water_year"].unique()):
    year_data = gage_df[gage_df["water_year"] == yr]
    Q = year_data["Q"].values.astype(float)
    n_days = len(Q)
    n_na = np.isnan(Q).sum()
    na_frac = n_na / n_days if n_days > 0 else 1.0
    print(f"  WY {yr}: {n_days} days, {n_na} NaN ({na_frac:.1%}), "
          f"Q range [{np.nanmin(Q):.4f}, {np.nanmax(Q):.4f}]")

# Apply filter year by year (matching analyze_baseflow_indices)
print(f"\nBFI per water year:")
for yr in sorted(gage_df["water_year"].unique()):
    year_data = gage_df[gage_df["water_year"] == yr].sort_values("dowy")
    Q = year_data["Q"].values.astype(float)

    if len(Q) < BASEFLOW_MIN_DAYS:
        print(f"  WY {yr}: SKIP (too few days: {len(Q)})")
        continue

    na_frac = np.isnan(Q).sum() / len(Q)
    if na_frac > BASEFLOW_MAX_MISSING_FRAC:
        print(f"  WY {yr}: SKIP (too many NaN: {na_frac:.1%})")
        continue

    bf = lyne_hollick_filter(Q, alpha=LYNE_HOLLICK_ALPHA, passes=LYNE_HOLLICK_PASSES)
    total_q = np.nansum(Q)
    total_bf = np.nansum(bf)
    bfi = total_bf / total_q if total_q > 0 else np.nan

    n_bf_nan = np.isnan(bf).sum()
    n_q_valid = (~np.isnan(Q)).sum()
    n_bf_valid = (~np.isnan(bf)).sum()

    print(f"  WY {yr}: BFI={bfi:.4f}, Q_valid={n_q_valid}, BF_valid={n_bf_valid}, "
          f"BF_NaN={n_bf_nan}, sum(Q)={total_q:.2f}, sum(BF)={total_bf:.2f}")

    # If there are filter-introduced NaNs, show where
    if n_bf_nan > np.isnan(Q).sum():
        q_valid_bf_nan = (~np.isnan(Q)) & np.isnan(bf)
        n_contaminated = q_valid_bf_nan.sum()
        print(f"    ** {n_contaminated} positions: Q valid but BF=NaN (filter contamination)")

# Show benchmark output for comparison
print(f"\nBenchmark CSV values:")
r_df = pd.read_csv(Path(__file__).parent / "r_signatures.csv", low_memory=False)
py_df = pd.read_csv(Path(__file__).parent / "python_signatures.csv", low_memory=False)
jl_df = pd.read_csv(Path(__file__).parent / "julia_signatures.csv", low_memory=False)

for label, df in [("R", r_df), ("Python", py_df), ("Julia", jl_df)]:
    row = df[df["gage_id"].astype(str) == GAGE]
    if len(row) > 0:
        bfi_lh = row["BFI_LyneHollick_mean"].values[0]
        bfi_eck = row["BFI_Eckhardt_mean"].values[0] if "BFI_Eckhardt_mean" in row.columns else "N/A"
        print(f"  {label}: BFI_LyneHollick_mean={bfi_lh}, BFI_Eckhardt_mean={bfi_eck}")
    else:
        print(f"  {label}: gage not found")
```

**Step 2: Run diagnostic**

```bash
py benchmarks/diagnose_bfi_canadian.py
```

Examine the output. Look for:
- Years where BF_NaN count is much higher than Q NaN count → filter contamination
- Differences in which years are included vs excluded
- Whether BFI values per year differ from expected

**Step 3: Root-cause and fix**

Based on diagnostic output, apply the fix. Most likely scenarios:

**(a) Filter NaN contamination is worse in Python than R/Julia:**
Fix the BFI calculation to use paired masking (like Julia does):
```python
# In analyze_baseflow_indices, replace np.nansum approach:
valid_mask = ~np.isnan(Q) & ~np.isnan(bf)
if valid_mask.sum() > 0:
    bfi = np.sum(bf[valid_mask]) / np.sum(Q[valid_mask])
```

**(b) Data preprocessing differs for Canadian gages:**
Fix the water year assignment or dowy computation in the benchmark script.

**(c) Sort order differs:**
Ensure `year_data.sort_values("dowy")` is applied consistently.

**Step 4: Verify fix**

```bash
py benchmarks/diagnose_bfi_canadian.py  # Should show BFI ~0.88 for 05CB007
```

Then spot-check 2-3 more Canadian gages from the top-20 list.

**Step 5: Commit**

```bash
git add python/streamflow_signatures/baseflow.py benchmarks/diagnose_bfi_canadian.py
git commit -m "fix(python): resolve BFI_LyneHollick discrepancy for Canadian gages

[Description of actual fix based on diagnosis]"
```

---

## Tier 2: Algorithmic Differences (Issues 4-6)

These require harmonizing naming conventions and statistical methods.

---

### Task 4: Fix FDC Column Naming Mismatch + Investigate Julia FDC Divergence

**Root Cause (naming):** R uses `FDCall/FDC90th/FDCmid`, Python/Julia use `FDC_all/FDC_90th/FDC_mid`. This means 24 columns are excluded from cross-language comparison.

**Root Cause (Julia divergence):** When manually compared (bypassing naming), R-Julia FDC slopes show rho=0.81-0.93, while R-Python is 0.97-0.99. The Julia FDC module uses a different exceedance range for mid-flows: Julia uses `20-70%` vs R's `20-80%`.

**Files:**
- Modify: `julia/src/fdc.jl:109-113` (fix mid-range)
- Modify: `python/streamflow_signatures/fdc.py` OR `R/helperFunctions.R` (harmonize naming — pick ONE to change)
- Test: `julia/test/runtests.jl`

**Step 1: Decide naming convention**

The naming should be harmonized. Since R is canonical, the cleanest fix is to rename Python/Julia columns to match R: `FDC_all` → `FDCall`, `FDC_90th` → `FDC90th`, `FDC_mid` → `FDCmid`. This keeps R unchanged.

**However**, if the R naming is considered inconsistent with other metrics (all others use underscores), it may be better to rename R's output. This is a **decision point for the user**. For this plan, we'll rename Python/Julia to match R.

**Step 2: Fix Julia FDC mid-range (20-70% → 20-80%)**

In `julia/src/fdc.jl`, find the exceedance range definitions (~lines 109-113). Change the FDC_mid range from `(0.2, 0.7)` to `(0.2, 0.8)` to match R's `fdc$exceedance >= 0.2 & fdc$exceedance <= 0.8`.

```julia
# Before:
ranges = Dict(
    "FDC_all" => (0.0, 1.0),
    "FDC_90th" => (0.9, 1.0),
    "FDC_mid" => (0.2, 0.7),   # <-- WRONG
)

# After:
ranges = Dict(
    "FDCall" => (0.0, 1.0),    # renamed to match R
    "FDC90th" => (0.9, 1.0),   # renamed to match R
    "FDCmid" => (0.2, 0.8),    # renamed + range fixed
)
```

Also update the metric names list and any references to the old names.

**Step 3: Rename Python FDC columns**

In `python/streamflow_signatures/fdc.py`, change the column renaming at the end of `analyze_fdc_trends()`:

```python
# Before:
rename_map = {"slp_all": "FDC_all", "slp_90th": "FDC_90th", "slp_mid": "FDC_mid"}
# After:
rename_map = {"slp_all": "FDCall", "slp_90th": "FDC90th", "slp_mid": "FDCmid"}
```

Also update `EXPECTED_SIGNATURE_BASES` in `python/streamflow_signatures/io.py`.

**Step 4: Update tests and config**

Update any references to `FDC_all` etc. in:
- `julia/test/runtests.jl`
- `python/tests/test_signatures.py`
- `config/signatures_config.json` (if FDC base names are listed)

**Step 5: Run tests and comparison**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
py -m pytest python/tests/ -v
py benchmarks/compare_three_way.py  # FDC columns should now appear in common set
```

**Step 6: Commit**

```bash
git add julia/src/fdc.jl python/streamflow_signatures/fdc.py python/streamflow_signatures/io.py julia/test/runtests.jl python/tests/test_signatures.py
git commit -m "fix: harmonize FDC column names to match R and fix Julia mid-range

Renamed FDC_all→FDCall, FDC_90th→FDC90th, FDC_mid→FDCmid in Python
and Julia to match canonical R naming. Fixed Julia FDC mid-range from
20-70% to 20-80% matching R's definition."
```

---

### Task 5: Fix Julia Spearman Correlation — ordinalrank vs tiedrank

**Root Cause:** Julia's `spearman_correlation()` (stats.jl:150-151) uses `ordinalrank()` which assigns sequential ranks to ties (1,2,3 for three tied values). R's `cor.test(method="spearman")` and Python's `scipy.stats.spearmanr()` use **average/fractional ranks** (2,2,2 for three ties). This matters enormously for integer-valued metrics like pulse counts, which have many ties.

**Impact:** `n_high_pulses_year_spearman_pval` R-Julia rho=0.50, despite `_mean` agreement of 0.9997. Also affects all other Spearman rho/pval columns (~100 columns with rho < 0.99).

**Files:**
- Modify: `julia/src/stats.jl:150-151`
- Test: `julia/test/runtests.jl`

**Step 1: Write a failing test exposing the tie-ranking issue**

Add to `julia/test/runtests.jl` inside the `@testset "Stats Module"` block:

```julia
@testset "Spearman handles ties with average rank" begin
    # Integer data with many ties — this is what pulse metrics look like
    x = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = Float64[3, 3, 3, 4, 4, 5, 5, 5, 6, 6]  # many ties

    rho, pval = StreamflowSignatures.spearman_correlation(x, y)

    # Reference: scipy.stats.spearmanr([1..10], [3,3,3,4,4,5,5,5,6,6])
    # Returns rho ≈ 0.9567, pval ≈ 0.0000726
    @test abs(rho - 0.9567) < 0.01
    @test pval < 0.001
end
```

**Step 2: Run test to verify it fails**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — rho differs due to ordinalrank producing different ranks for ties.

**Step 3: Fix stats.jl — change ordinalrank to tiedrank**

In `julia/src/stats.jl`, replace lines 150-151:

```julia
    # Calculate ranks (use tiedrank for average ranks, matching R/Python)
    rank_x = tiedrank(x_valid)
    rank_y = tiedrank(y_valid)
```

`tiedrank()` from StatsBase.jl assigns average ranks to tied values, matching R's `rank()` and Python's `scipy.stats.rankdata()`.

**Step 4: Run test to verify it passes**

```bash
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS — rho ≈ 0.9567.

**Step 5: Commit**

```bash
git add julia/src/stats.jl julia/test/runtests.jl
git commit -m "fix(julia): use tiedrank instead of ordinalrank for Spearman correlation

ordinalrank assigns sequential ranks to ties (1,2,3), while R and
Python use average ranks (2,2,2). This caused large divergences in
Spearman rho/pval for integer-valued metrics like pulse counts."
```

---

### Task 6: Investigate Julia Concavity — Extra Non-NA Gages

**Root Cause:** Julia produces concavity for 1,571 gages where R/Python return NA. Where all three have values, Python-Julia rho=1.0 (perfect). R-Python rho=0.94, suggesting R computes slightly different values.

**Hypothesis:** Julia's recession event identification has a less strict threshold, qualifying more gages. OR Julia counts recession events differently, reaching the `min_events=25` threshold for more gages.

**Files:**
- Investigate: `julia/src/recession.jl` — event identification and counting
- Compare: `R/helperFunctions.R` — same logic
- Possibly modify: `julia/src/recession.jl`

**Step 1: Diagnostic — compare event counts**

Add diagnostic output to understand event count differences. Create a small test script:

```python
# benchmarks/diagnose_concavity.py
"""Compare recession event counts across R/Python/Julia."""
import pandas as pd
import numpy as np

r = pd.read_csv("benchmarks/r_signatures.csv", low_memory=False)
py = pd.read_csv("benchmarks/python_signatures.csv", low_memory=False)
jl = pd.read_csv("benchmarks/julia_signatures.csv", low_memory=False)

for df in (r, py, jl):
    df["gage_id"] = df["gage_id"].astype(str).str.strip()

common = set(r.gage_id) & set(py.gage_id) & set(jl.gage_id)
r = r[r.gage_id.isin(common)].set_index("gage_id")
py = py[py.gage_id.isin(common)].set_index("gage_id")
jl = jl[jl.gage_id.isin(common)].set_index("gage_id")

# Gages where Julia has concavity but R doesn't
r_na = r["concavity_mean"].isna()
jl_na = jl["concavity_mean"].isna()
extra_jl = r_na & ~jl_na

print(f"Gages where Julia has concavity but R doesn't: {extra_jl.sum()}")
print(f"Gages where R has concavity but Julia doesn't: {(~r_na & jl_na).sum()}")

# For extra Julia gages, check if log_a_events and b_events also have values
extra_gages = extra_jl[extra_jl].index[:20]
print(f"\nSample extra Julia gages (first 20):")
for g in extra_gages:
    jl_conc = jl.loc[g, "concavity_mean"]
    jl_b = jl.loc[g, "b_events_mean"] if "b_events_mean" in jl.columns else "N/A"
    r_b = r.loc[g, "b_events_mean"] if "b_events_mean" in r.columns else "N/A"
    print(f"  {g}: Jl concavity={jl_conc:.4f}, Jl b_events={jl_b}, R b_events={r_b}")
```

**Step 2: Run diagnostic and analyze**

```bash
py benchmarks/diagnose_concavity.py
```

If Julia counts more recession events than R for the same data, check:
- Julia's `identify_recession_events()` look-ahead algorithm vs R's
- Whether Julia allows shorter recession events
- Whether Julia handles data gaps differently in event identification

**Step 3: Fix based on findings**

Most likely fix: tighten Julia's event identification to match R's exact behavior. The concavity itself is a derived metric (b_first_half - b_second_half), so if the event identification matches, concavity will match.

**Step 4: Verify**

```bash
py benchmarks/compare_three_way.py  # Check concavity NA counts align
```

**Step 5: Commit**

```bash
git add julia/src/recession.jl benchmarks/diagnose_concavity.py
git commit -m "fix(julia): align recession event identification with R

[Description based on actual findings]"
```

---

## Tier 3: Methodology and Filtering Differences (Issues 7-9)

These may require methodology discussions rather than pure bug fixes.

---

### Task 7: Investigate Runoff Ratio Linear Slope Divergence

**Root Cause:** `annual_runoff_ratio_linear_slp` R-Python rho=0.47. The `_mean` values agree well (0.9922), but trend slopes diverge wildly. This suggests the underlying annual time series has subtle differences that get amplified in slope computation.

**Impact:** Affects 8 columns (annual_runoff_ratio × 8 stats), plus seasonal runoff ratios (~40 columns total).

**Files:**
- Investigate: `python/streamflow_signatures/runoff_ratios.py`
- Compare: `R/helperFunctions.R` — `analyze_Q_PPT_relationships()`
- Investigate: `benchmarks/run_python_benchmark.py` — how climate data is merged

**Step 1: Diagnostic — compare annual runoff ratios per year for problem gages**

Create `benchmarks/diagnose_runoff_ratios.py`:

```python
"""Compare annual runoff ratio time series for divergent gages."""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
from streamflow_signatures.runoff_ratios import analyze_Q_PPT_relationships
from streamflow_signatures.io import add_water_year_columns

# Load data for a problem gage
GAGE = "11239300"  # Worst outlier for runoff ratio linear_slp

parquet_path = Path(r"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet")
climate_path = Path(r"D:\processedOuts_feb2026\daymet_1980_2023.parquet")

df = pd.read_parquet(parquet_path)
df.columns = [c.lower() for c in df.columns]
gage_df = df[df["gage_id"].astype(str) == GAGE].copy()

# Load and merge climate
climate = pd.read_parquet(climate_path)
climate.columns = [c.lower() for c in climate.columns]
climate = climate[climate["gage_id"].astype(str) == GAGE]

gage_df = gage_df.merge(climate[["gage_id", "date", "PPT"]], on=["gage_id", "date"], how="left")
gage_df = add_water_year_columns(gage_df)

print(f"Gage {GAGE}: {len(gage_df)} rows")
print(f"PPT available: {gage_df['PPT'].notna().sum()} / {len(gage_df)}")

# Show annual Q/P and ratio
for yr in sorted(gage_df["water_year"].unique()):
    yd = gage_df[gage_df["water_year"] == yr]
    q_sum = yd["Q"].sum()
    p_sum = yd["PPT"].sum() if "PPT" in yd.columns else np.nan
    ratio = q_sum / p_sum if p_sum > 10 else np.nan
    print(f"  WY {yr}: Q_sum={q_sum:.2f}, PPT_sum={p_sum:.2f}, ratio={ratio:.4f}")

# Run the signature function
result = analyze_Q_PPT_relationships(gage_df)
for k, v in sorted(result.items()):
    if "annual_runoff_ratio" in k:
        print(f"  {k}: {v}")
```

**Step 2: Run diagnostic and compare with R output**

```bash
py benchmarks/diagnose_runoff_ratios.py
```

Compare the annual time series to see where they diverge. Check:
- Is the climate data merge correct (matching dates)?
- Does Python use the same minimum PPT threshold (10mm annual)?
- Are the same years included/excluded?

**Step 3: Fix based on findings**

Likely scenarios:
- Climate data merge produces different PPT for some years → fix merge logic
- Different years qualify due to minimum PPT threshold → align thresholds
- Linear slope is sensitive to outlier years → consider whether this is acceptable

**Step 4: Commit**

```bash
git add [affected files]
git commit -m "fix: align runoff ratio computation across languages

[Description based on findings]"
```

---

### Task 8: Investigate avg_storage Trend Statistic Disagreement

**Root Cause:** All three implementations disagree on `avg_storage` trend statistics (rho 0.77-0.83 for p-values), though means agree well (0.99). This metric is inherently noisy (ignores ET, simplified water balance).

**Impact:** 8 columns (avg_storage × 8 stats), but primarily p-value and slope columns.

**Files:**
- Compare: `python/streamflow_signatures/storage.py`, `julia/src/storage.jl`, `R/helperFunctions.R`

**Step 1: Compare storage computation logic**

Check for differences in:
- How cumulative storage (cumsum(P - Q)) is computed
- How annual storage is interpolated at mean discharge
- Whether years are filtered differently (min_days, max_na_frac)

**Step 2: Diagnostic script**

Similar pattern to Tasks 3 and 7: pick a problem gage, compare year-by-year storage values.

**Step 3: Fix or document as acceptable**

If the `_mean` correlation is >0.99, the algorithm is fundamentally correct. The p-value/slope divergence may be acceptable noise from the simplified water balance model. In that case, document the expected divergence in SIGNATURES.md and the comparison report.

---

### Task 9: Investigate Extra Gages in Python/Julia vs R

**Root Cause:** R qualifies 5,707 gages vs Python/Julia's 7,636. Filtering parameters in `config/signatures_config.json` are shared, so the difference must be in how filters are applied.

**Impact:** Not a correctness issue per se, but understanding the difference helps validate that all implementations use the same data quality criteria.

**Files:**
- Investigate: `benchmarks/run_r_benchmark.R`, `benchmarks/run_python_benchmark.py`, `benchmarks/run_julia_benchmark.jl`

**Step 1: Compare qualifying-gage logic**

Check each benchmark script for:
- How many water years are counted per gage
- What constitutes a "good" year (min days, max NA fraction)
- Whether R applies additional filters (e.g., minimum flow threshold) that Python/Julia skip

**Step 2: Diagnostic**

```python
# Check which gages are in Python but not R
py_only = set(py_df.gage_id) - set(r_df.gage_id)
# For a sample of py_only gages, check why R rejected them
# Load raw parquet data and compute qualifying years using both criteria
```

**Step 3: Align or document**

If R's stricter criteria are intentional (e.g., it requires `min_Q_value_and_days` that Python/Julia skip), document the difference. If it's a bug, align the filtering.

---

## Verification

After all fixes, run:

```bash
# 1. Julia tests
cd julia && julia --project=. -e 'using Pkg; Pkg.test()'

# 2. Python tests
py -m pytest python/tests/ -v

# 3. Re-run benchmarks (if data is available on D: drive)
# This takes hours — only do if critical
# py benchmarks/run_python_benchmark.py
# julia --project=julia benchmarks/run_julia_benchmark.jl

# 4. Three-way comparison
py benchmarks/compare_three_way.py
```

**Success criteria:**
- All columns with rho < 0.99 in Tier 1 issues should improve to > 0.99
- D25_to_D75 should no longer be all-NaN in Julia
- log_a_seasonality_minimum R-Julia correlation should flip from -0.98 to +0.98
- BFI_LyneHollick R-Python correlation should improve from 0.987 to > 0.995
- Pulse metric Spearman p-values should improve from rho ~0.50 to > 0.95
- FDC columns should appear in common comparison set (no longer "only in R")
