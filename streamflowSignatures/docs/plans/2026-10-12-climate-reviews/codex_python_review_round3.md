**NO-GO for an unqualified validated mirror.** All five round-2 findings pass their original reproductions. Three additional boundary defects remain: incomplete validation of explicit preprocessing overrides, a dtype-dependent exception to the observation contract, and cumulative-fraction canonicalisation into forbidden endpoints. A separate canonical follow-up remains for direct complex input. No annual mismatch was found in the ordinary numeric fixtures under the shipped configuration.

Repository source files were not edited. Scripts, logs and this report are under `/tmp`. Julia configuration variants were loaded from source in fresh modules to avoid the package's documented precompile-cache behavior.

## Findings

1. **MINOR — Explicit Q/PPT/SWE preprocessing overrides bypass the new validation contract.**

   **Location:** `python/streamflow_signatures/io.py:574–595`, especially the temperature-only validation loop at `:592`; canonical `julia/src/io.jl:250–254`.

   **Julia/spec:** Explicit missing-count and interpolation-gap policies must be whole, finite, non-negative counts. Julia types all eight arguments as `Int` and checks all eight for non-negativity.

   **Port:** Only `max_raw_na_temp` and `max_gap_temp` are validated. The six corresponding Q/PPT/SWE settings accept negative values, fractions, NaN and booleans. An invalid negative override silently excludes valid years; NaN can disable a qualification limit.

   **Minimal reproduction:** Complete WY2001; finite valid Q, PPT and SWE = 2/day, TMIN = 1, TMAX = 10. Call `preprocess_daily_data(df, config=...)`.

   | Override | Expected / Julia | Actual Python |
   |---|---|---|
   | `{"max_gap_ppt": -1}` | Reject; Julia `ArgumentError` | Returns normally; `valid_years=[2001]`, `valid_climate_years=[]` despite no missing PPT |
   | `{"max_gap_days": -1}` | Reject; Julia `ArgumentError` | Returns normally; all four valid-year lists empty |
   | `{"max_gap_ppt": 3.9}` | Reject; Julia `TypeError` | Accepts; precipitation year `[2001]` |
   | `{"max_raw_na_ppt": NaN}` with 31 isolated PPT NaNs | Reject; Julia `TypeError` | Accepts, interpolates all 31, qualifies `[2001]`, `prcp_annual=730`; shipped policy correctly rejects the same input |

   The same invalid-value acceptance was reproduced for `max_raw_na_per_year`, `max_gap_days`, `max_raw_na_ppt`, `max_gap_ppt`, `max_raw_na_swe` and `max_gap_swe`, using `-1`, `3.9`, NaN and `True`. Both temperature settings reject all four values.

   **Proposed fix:** Apply the finite/whole/non-boolean/non-negative validator to all eight count/gap settings before processing any year. Normalize accepted whole values to integers. Add direct-override tests, including a NaN raw-count limit with 31 isolated missing observations.

2. **MINOR — Arrow decimal columns bypass the stipulated `numbers.Real` observation contract.**

   **Location:** `python/streamflow_signatures/climate.py:65–74`; temperature preprocessing calls this helper at `python/streamflow_signatures/io.py:798–799`.

   **Spec:** Python kernels and temperature preprocessing must accept `numbers.Real` observations, including ordinary booleans, and treat other types as missing. `decimal.Decimal` is not a `numbers.Real`.

   **Port:** The fast path treats every non-complex pandas numeric dtype as admissible. Pandas identifies Arrow decimal columns as numeric even though their scalar values are `Decimal`. The same values therefore change qualification and annual results solely when the column dtype changes.

   **Minimal reproduction:** Complete WY2001 and finite numeric Q. Construct columns with:

   ```python
   from decimal import Decimal
   import pandas as pd
   import pyarrow as pa

   dtype = pd.ArrowDtype(pa.decimal128(10, 2))
   df["PPT"] = pd.Series([Decimal("2.00")] * 365, dtype=dtype)
   df["TMIN"] = pd.Series([Decimal("1.00")] * 365, dtype=dtype)
   df["TMAX"] = pd.Series([Decimal("10.00")] * 365, dtype=dtype)
   ```

   **Expected:** Under the stated Python type rule, temperature-valid years `[]`, and both direct annual families entirely NaN.

   **Actual:** `valid_temp_years=[2001]`; direct `prcp_annual=730`, `temp_max_annual=10`, `gdd_cum_annual=2007.5`. Changing only `dtype` to `object` gives the expected missing treatment: `valid_temp_years=[]`, precipitation and temperature annual values NaN. This is a Python contract defect; no equivalent Julia Decimal representation is assumed.

   **Proposed fix:** Restrict the fast path to dtypes known to represent admissible real observations, or use the scalar contract for extension dtypes such as Arrow decimal. Test identical scalar values across object and extension-array storage.

3. **MINOR — Canonicalisation can create 0% or 100% after the open-interval validation passes.**

   **Location:** `python/streamflow_signatures/config.py:374–382`; shared defect at `julia/src/config.jl:210–222`.

   **Julia/spec:** The declared configuration domain is a whole percentage strictly inside `(0,1)`. The newly canonicalised computational fraction must also satisfy that domain.

   **Port and Julia:** Validate the raw value's range, accept a rounded whole percentage, then canonicalise without rechecking the range.

   **Minimal reproduction:** Set `climate.cumulative_fraction` as below. Use complete WY2001 with PPT = 1 at dowy 100 and 200 and zero elsewhere.

   | Raw setting | Expected | Actual in both languages |
   |---|---|---|
   | `1e-14` | Reject: canonical percentage is 0 | Accepts `CLIMATE_CUM_FRACTION=0.0`; emits `prcp_0_dowy=1`, before the first storm |
   | `0.99999999999999` | Reject: canonical percentage is 100 | Accepts fraction `1.0`; emits `prcp_100_dowy=200` |

   This does not invalidate the repaired `0.500000000001` case: that now correctly canonicalises to `0.5`.

   **Proposed fix:** After computing the canonical integer percentage, require `1 <= percentage <= 99`, then divide by 100. Retain the documented crossing tolerance.

4. **NOTE — Julia's direct kernels still disagree with Python's required complex-observation handling.**

   **Location:** `julia/src/climate.jl:279`, `:443–444`, calling `julia/src/stats.jl:614–624`; compare `python/streamflow_signatures/climate.py:65–75`.

   **Spec / Python:** Complex observations are missing under the stipulated Python observation contract. Python correctly returns an entirely NaN annual row. Both preprocessors reject an all-complex temperature year.

   **Canonical behavior:** The direct Julia climate kernels call `coalesce_q`, whose numeric branch accepts `Number`, including complex types, then attempts `Float64(val)`.

   **Minimal reproduction:** Complete WY2001; TMIN = `1+0im`, TMAX = `10+0im` on every day. Python's corresponding `1+0j` / `10+0j` direct temperature input gives **19 NaNs**. Julia gives **17 finite metrics**, including maximum **10** and annual GDD **2007.5**; only the two frost metrics are NaN. Both preprocessors return `valid_temp_years=[]`. Changing the imaginary parts to 2 causes Julia's direct call to raise **`InexactError`**, while Python treats those observations as missing.

   **Proposed fix:** Give Julia's climate kernels an explicit `Real` observation conversion matching the new preprocessing policy. Keep Python's prescribed behavior. Prefer a climate-local conversion unless changing the shared `coalesce_q` contract across existing families is intended.

## Round-2 findings: original reproductions rerun

All entries below are **PASS for the original reproduction**, not a claim that every related boundary is covered.

| Round-2 finding | Result | Recomputed expected and actual values |
|---|---|---|
| **1. Non-finite preprocessing** | **PASS** | Complete WY2001 with PPT = 2 except dowy 21: both `+Inf` and `-Inf` become missing, interpolate to **2**, qualify `[2001]`, and produce annual PPT **730** in Python and Julia. Four consecutive infinities and 31 isolated infinities give `valid_climate_years=[]`. SWE gives the corresponding valid/rejected lists and isolated fill **2**. Thirty-one infinite Q observations give `valid_years=[]` and `valid_temp_years=[]`. Internal Q infinities on the original ramp fill to **1.02**. Direct non-finite kernels still void their annual rows. |
| **2. Observation types** | **PASS** | All-string direct PPT/TMIN/TMAX input yields NaN annual families in both languages. All-string temperature preprocessing gives `valid_temp_years=[]`. All-complex TMIN/TMAX likewise gives `[]`. `Fraction(1,2)` / `Fraction(21,2)`, compared with Julia rationals, qualifies `[2001]`; annual GDD **2007.5**, maximum **10.5**. The newly found Arrow-decimal exception is finding 2 above. |
| **3. Fractional season month** | **PASS** | Winter `[12,1,2.9]` is rejected by Python `ValueError` and Julia configuration error; no **90** winter total is emitted. Python also rejects months **0**, **13** and `True`. |
| **4. Near-whole cumulative fraction** | **PASS** | Raw `0.500000000001` loads as exactly **0.5** in both languages. Original storms of 1 at dowy **1** and **200** produce **`prcp_50_dowy=1`**, not 200. |
| **5. Canonical scalar strictness** | **PASS** | `gdd_base_c=true`, `record_anchored_decade_gate=1`, and `max_raw_na_per_year_temp=true` are rejected in both languages. Threshold list `[true]` also rejects. An overflowing `1e309` GDD base rejects. |

## Earlier regressions rerun

| Round-1 finding | Result | Actual values / behavior |
|---|---|---|
| **1. Interpolation arithmetic** | **PASS** | Both languages fill `[-5,NaN,NaN,NaN,1.7]` to `[-5,-3.325,-1.65,0.025000000000000355,1.7]`; annual GDD **1.7250000000000003**. The original frost fixture fills dowy 13 to **+0.0** and reports last frost **12**. |
| **2. Signed zeros** | **PASS** | PPT/TMAX `[-0,+0,…]`, TMIN `[+0,-0,+0,…]`: all three extreme values are **+0.0**, bits **0x0000000000000000**, at dowy **1**, directly and after preprocessing in both languages. |
| **3. Zero temperature gap** | **PASS** | One internal missing TMIN, gap limit **0**: normal return, `valid_years=[2001]`, `valid_temp_years=[]`. Negative temperature limits reject. |
| **4. Original invalid settings** | **PASS** | Window **5.9**, frost month **6.9**, temperature raw limit **30.9**, gap **3.9**, flag `"false"`, literal JSON `NaN` dry threshold and `Infinity` GDD base all reject in both languages. |
| **5. Scientific-notation labels** | **PASS** | `1e-5` yields **`prcp_ndays_abv_1em05mm`** under the revised convention. Additional Python labels: `1e-6 → 1em06`, `0.0001 → 0p0001`, `0.5 → 0p5`. |
| **6. String temperatures** | **PASS** | TMIN `"1.0"`, TMAX `"10.0"`: `valid_temp_years=[]`; numeric Q retains `[2001]`. |
| **7. Infinite annual observations** | **PASS** | Original WY2000, otherwise-zero PPT with infinity at dowy 21: direct kernel yields all **20** NaNs; preprocessing repairs to zero and annual precipitation is **0**. Infinite direct temperature input voids all **19** values. |

## What I verified and found correct

- **Tests:** **93 passed**: 15 climate tests plus 78 tests covering phase-1 review fixes, area-normalization gating, schema validation, snow and annual collection. These tests do not catch the findings above.
- **Cross-language annual values:** Regenerated the shipped Julia fixture: **1,170** annual rows, zero bit/NaN-pattern mismatches, and **624** summary comparisons passing `atol=1e-12, rtol=4e-14`. Reran the prior adversarial fixtures: **10,433** annual rows and **19,968** summary fields, zero mismatches/tolerance failures. Reran the original independent fixture: **4,134** annual rows and **9,984** summary fields, also zero failures. These counts include intentionally invalid annual NaN rows. The 46-year random preprocessed subset contains **1,794** exactly matching annual values.
- **Interpolation parity:** All **19,363** compared preprocessed values in each of PPT, TMIN and TMAX matched Julia's bits or NaN pattern. The random-fixture maximum absolute summary difference was **4.547473508864641e-13**. The other adversarial fixture's maximum was **1.8189894035458565e-12** near an annual mean of **7305**, within the relative tolerance.
- **Crossing rule and sequential accumulation:** Two separated seven-day blocks of **0.1** give total/final cumulative **1.4000000000000001**, crossing **7**. Uniform **0.1** across WY2000 gives total **36.60000000000025**, crossing **183**. `[49,51,0,…]` crosses on **2**. Two fifty-day blocks give total/final cumulative **9.99999999999998**, crossing **50**; Python 3.12 builtin `sum` instead gives **10.0**.
- **Closed-form precipitation:** WY2001 at **2/day** yields annual **730**, five-day maximum **10**, D50/centre **183**, SD/Gini **0**, SI **0.022374429223744285**. One **100** storm gives Gini **0.9972602739726028**, SI **1.8333333333333337**, sample SD **5.234239225902137**, CV **19.1049731745428**, agreeing with `364/365`, `11/6`, `100/sqrt(365)` and `sqrt(365)` under sequential-double rounding.
- **Calendars and seasons:** WY1900/2001/2100 have **365** days; WY2000 has **366**. At 2/day totals are **730/732**, centres **183/183.5**. Shuffling preserves results; duplicated dowy and truncated grids void both families. A corrupted month column still gives WY2000 winter/spring/summer/fall precipitation **461/368/645/910**. Zero-PPT WY2000 gives dry length/start/end **366/1/366**, with NaN D50/centre/CV/Gini/SI.
- **Daymet and frost:** Removing December 31, 2000 fails the direct complete-grid guard. Preprocessing restores dowy **92**, PPT **4** between 2 and 6, TMIN **-1** between -4 and 2. June 30 is included at dowy **274** in WY2000 and **273** in WY2001; July 1 is excluded from last frost. Frost-free timing is NaN. Exact-threshold/next-double precipitation probes give strict-above counts **5/3/1** at 10/25/50. The final five-day `[1,2,3,4,5]` window gives **15**.
- **Temperature policy and gate:** Thirty isolated missing values qualify; 31 fail. Three-day internal gaps qualify; four-day and edge gaps fail. Both temperature series must pass; negative temperatures remain admissible; failed Q qualification suppresses temperature qualification. Seven/eight frost years in the first decade yield slope **NaN/1.0**, retaining means **11.764705882352942/11.61111111111111**. Helpers leave the central NaN in `[0,NaN,+Inf]` and `[0,NaN,-Inf]` unfilled in both languages.
- **Orchestrator and existing families:** Explicit climate/temperature frames remain required. Temperature runs with `has_climate=False`, `area_normalized=False`: maximum-temperature mean **10**, zero precipitation keys. All **864** pre-existing fields in the additivity fixture were unchanged by temperature-column presence or enabled temperature calculations; temperature adds **152** fields with changepoints disabled.
- **Runner:** Lowercase aliases normalize to TMIN/TMAX. Executing the runner's actual window-filter and temperature-frame statements for WY1993–2025 retains only **[2025]** from temperature-valid years `[1992,2025,2026]`.

**Verdict: NO-GO for using the port as an unqualified validated mirror of Julia.** The original reproductions and shipped numeric fixtures pass. The remaining failures concern explicitly required public-input/configuration boundaries; two are Python-only and one is shared with Julia. Direct complex inputs also expose a canonical follow-up. No failure of ordinary finite numeric annual parity under the shipped configuration was observed.
