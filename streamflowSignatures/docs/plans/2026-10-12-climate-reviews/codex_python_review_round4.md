**Round 4: all four round-3 findings PASS. No new findings.** Repository source files were not edited. Reproduction scripts and logs are under `/tmp/round4_*`; Julia configuration variants were evaluated in fresh modules.

1. **MINOR, resolved — preprocessing overrides.** `python/streamflow_signatures/io.py:598–610`.

   Expected: reject invalid limits before processing. On the original complete WY2001 fixture, each of all eight limits now rejects `-1`, `3.9`, `NaN`, and `True`: **32/32 Python ValueErrors**. Julia rejects negatives with `ArgumentError` and the other values with `TypeError`. The original 31-isolated-PPT-NaNs bypass now raises `ValueError`; default qualification remains `[]`. All five flags reject `0`, `1`, `"false"`, and `None`, and accept Boolean values. The proposed shared validation is implemented; no further fix proposed.

2. **MINOR, resolved — Arrow decimal observations.** `python/streamflow_signatures/climate.py:65–79`.

   Original WY2001 Decimal inputs: PPT=2, TMIN=1, TMAX=10. Expected and actual, for both object and Arrow decimal storage: **20 precipitation NaNs, 19 temperature NaNs, `valid_temp_years=[]`**. Annual precipitation, temperature maximum, and annual GDD are all NaN. The restricted fast path fixes the bypass. Valid NumPy, nullable, and Arrow real/Boolean unit-valued inputs retain annual precipitation/GDD **365/365** and temperature qualification `[2001]`. No further fix proposed.

3. **MINOR, resolved — forbidden canonical percentage endpoints.** `python/streamflow_signatures/config.py:382–385`; `julia/src/config.jl:223`; `rpkg/R/config.R:300`.

   Expected and actual: raw **1e-14** and **0.99999999999999** reject in all three languages before emitting annual metrics. Interior **0.500000000001 → 0.5** still passes; the original storms at dowy 1 and 200 give Python **`prcp_50_dowy=1`**. Valid endpoints **0.01/0.99** remain accepted in all three languages. The proposed canonical-range check is implemented; no further fix proposed.

4. **NOTE, resolved — Julia complex observations.** `julia/src/climate.jl:129–134`, used at `:296`, `:460–461`.

   Expected and actual: zero-imaginary and nonzero-imaginary complex temperature inputs produce **19 NaNs**, with **`valid_temp_years=[]`**, in both languages; Julia no longer throws `InexactError`. Complex precipitation produces **20 NaNs**. Julia Bool observations still give annual precipitation/GDD **365/365**; rational temperature inputs retain maximum **10.5**, GDD **2007.5**, and qualification `[2001]`. The Real-only conversion fixes the discrepancy; no further fix proposed.

**What I verified and found correct**

- **Numerical parity:** regenerated fixtures yield **15,737 annual entries**, with zero bit/NaN-pattern mismatches, and **30,576 summary comparisons**, all passing `atol=1e-12, rtol=4e-14`. Each of PPT/TMIN/TMAX has **19,363** exactly matching preprocessed observations. **93 tests pass.**
- **Independent formulas and accumulation:** the 100-mm single-storm fixture gives Gini **0.9972602739726028**, SI **1.8333333333333337**, SD **5.234239225902137**, and CV **19.1049731745428**, matching `364/365`, `11/6`, `100/√365`, and `√365` to rounding. Two 50-day 0.1-mm blocks give total/final cumulative **9.99999999999998**, crossing **50**. Exact logical counts satisfy the clarified exemption.
- **Calendars and missingness:** century/leap grid checks retain totals **730/732** at 2 mm/day; invalid grids fail. Daymet restoration gives dowy **92**, PPT **4**, TMIN **−1**. June 30 cutoffs remain **274/273** in WY2000/2001. Thirty isolated missing temperatures and three-day gaps qualify; 31 and four-day gaps fail. Original interpolation retains annual GDD **1.7250000000000003**.
- **Gate and integration:** seven/eight first-decade frost years retain slopes **NaN/1.0**. All **864** compared existing summary fields remain unchanged with temperature columns/calculations added. Temperature-only, unnormalized input gives maximum-temperature mean **10**, with no precipitation keys. The runner retains only **[2025]** from temperature-valid `[1992,2025,2026]` for WY1993–2025.
- **R supplementary checks:** original Q/PPT invalid overrides and malformed seasons reject; configured January-only winter completeness is now **1,1,1**, all true.

**Verdict: GO for using this Python port as a validated mirror of Julia under the reviewed contract.**
