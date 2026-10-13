**NO-GO.** I reproduced an annual GDD mismatch and a one-day frost-date mismatch under the shipped configuration. All 13 implementer tests and the shipped cross-language fixture nevertheless pass.

Repository files were not modified. Reproductions used complete, chronologically ordered water-year frames, valid Q, and the `AnnualCollector` to inspect values before statistical gating.

## Findings

1. **BLOCKER — Temperature interpolation changes annual values and frost dates.**

   **Location:** [python/streamflow_signatures/io.py:799](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/io.py:799); canonical `julia/src/io.jl:199`.

   **Julia/spec:** Annual values must be bit-identical. Julia interpolates each missing value using the scalar expression `left + (right-left)*offset/span`.

   **Port:** Pandas interpolation produces different rounding. These differences survive the otherwise matching climate calculations.

   **Minimal reproductions:** Use WY2001; positions below are **one-based dowy**.

   - Set both temperature series to −10°C except dowy 10–14:
     `[-5, NaN, NaN, NaN, 1.7]`.
     Julia fills dowy 13 with `0.025000000000000355`; Python produces `0.025000000000000133`.
     **Expected annual GDD:** `1.7250000000000003`.
     **Actual:** `1.725`.
   - Set TMIN to 1°C, TMAX to 10°C, then set TMIN dowy 10–14 to
     `[-0.30000000000000004, NaN, NaN, NaN, 0.1]`.
     Julia fills dowy 13 with `0.0`; Python produces `−2.7755575615628914e−17`, creating an additional frost day.
     **Expected `frost_last_dowy`: 12. Actual: 13.**

   **Fix:** Port Julia’s gap-filling loop and arithmetic evaluation order explicitly. Add both reproductions as exact cross-language regression tests. The existing interpolation test at `python/tests/test_climate_metrics.py:280` uses constant neighbours, which cannot expose this failure.

2. **MAJOR — Signed-zero extrema disagree with Julia; Julia itself violates the numerical first-tie rule.**

   **Location:** [python/streamflow_signatures/climate.py:246](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/climate.py:246), also lines 357–362; canonical `julia/src/climate.jl:279` and `434`.

   **Julia/spec:** Decisions require the first occurrence on ties. Julia’s extrema distinguish signed zeros, preferring `+0.0` for maxima and `−0.0` for minima.

   **Port:** NumPy treats these values as equal when finding their positions.

   **Minimal reproduction:** Complete WY2000:

   | Series | Input | Julia annual value / dowy | Python annual value / dowy |
   |---|---|---|---|
   | PPT | `[-0.0, +0.0, …]` | `+0.0 / 2` | `−0.0 / 1` |
   | TMAX | `[-0.0, +0.0, …]` | `+0.0 / 2` | `−0.0 / 1` |
   | TMIN | `[+0.0, −0.0, +0.0, …]` | `−0.0 / 2` | `+0.0 / 1` |

   These differences persist through preprocessing. Under numerical equality, the specification requires dowy 1; the canonical implementation therefore also needs correction.

   **Fix:** Normalize physical zero values consistently in both implementations before extrema calculations, then retain the first numerical tie. Test float bits as well as positions: ordinary `==` cannot distinguish signed zeros.

3. **MAJOR — A zero interpolation-gap setting crashes preprocessing when temperature data are present.**

   **Location:** [python/streamflow_signatures/config.py:228](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/config.py:228), [python/streamflow_signatures/io.py:799](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/io.py:799).

   **Julia/spec:** A maximum gap of zero permits no interpolation. Julia leaves missing temperatures unfilled, excludes the temperature year, and still returns streamflow results.

   **Port:** Config accepts zero, then passes `limit=0` to Pandas, which raises `ValueError: Limit must be greater than 0`.

   **Minimal reproduction:** Complete WY2001 with valid Q/PPT and one internal TMIN NaN; set `max_interpolation_gap_temp=0`.

   **Expected:** `valid_years=[2001]`, `valid_temp_years=[]`, normal return.
   **Actual:** Preprocessing aborts. Consequently, adding temperature columns can prevent existing families from running. The benchmark calls preprocessing without a surrounding recovery block.

   **Fix:** Handle zero by skipping interpolation. Reject negative limits and negative raw-NA thresholds during configuration validation.

4. **MAJOR — Config validation accepts settings Julia rejects, silently changing their meaning.**

   **Location:** [python/streamflow_signatures/config.py:227](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/config.py:227), lines 309–334, and `load_config()`.

   **Julia/spec:** Settings must retain their declared meaning and match canonical validation. Julia rejects nonintegral values passed to `Int`, string booleans, and the nonstandard JSON tokens tested below.

   **Port:** `int()` truncates; `bool()` accepts arbitrary truthy objects; Python’s JSON decoder accepts `NaN`/`Infinity`, and climate validation does not consistently require finite values.

   **Minimal reproductions:**

   | Setting | Expected Julia result | Actual Python result |
   |---|---|---|
   | `multi_day_max_window_days: 5.9` | Reject | Uses 5 |
   | `frost_last_end_month: 6.9` | Reject | Uses June |
   | `max_raw_na_per_year_temp: 30.9` | Reject | Uses 30 |
   | `max_interpolation_gap_temp: 3.9` | Reject | Uses 3 |
   | `record_anchored_decade_gate: "false"` | Reject | Enables gate |
   | `dry_day_threshold_mm: NaN` | Reject JSON | Accepts |
   | `gdd_base_c: Infinity` | Reject JSON | Accepts |

   With the accepted NaN dry threshold, an entirely dry 365-day year produces **`prcp_intermit_annual=0`**, with missing start/end, instead of being rejected at configuration loading.

   **Fix:** Validate types, integrality, finiteness and ranges before conversion. Require actual booleans and reject nonstandard JSON constants.

5. **MINOR — Scientific-notation thresholds generate different column names.**

   **Location:** [python/streamflow_signatures/climate.py:65](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/climate.py:65).

   **Julia/spec:** Accepted configuration must generate the same schema in both languages.

   **Port:** Native float string formatting differs from Julia’s.

   **Minimal reproduction:** Set `wet_day_thresholds_mm` to `[0.00001]`, accepted by both validators.

   **Expected Julia base:** `prcp_ndays_abv_1p0e-5mm`.
   **Actual Python base:** `prcp_ndays_abv_1e-05mm`.

   All associated summary columns and annual signature names differ.

   **Fix:** Specify and implement one language-independent threshold formatter, including exponent notation.

6. **MINOR — Numeric strings qualify temperature years that Julia rejects.**

   **Location:** [python/streamflow_signatures/io.py:794](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/io.py:794); canonical `julia/src/io.jl:580`.

   **Julia/spec:** Canonical preprocessing treats values that are not `Real` as missing, then applies the temperature NA policy.

   **Port:** `pd.to_numeric(..., errors="coerce")` converts numeric strings into observations.

   **Minimal reproduction:** Complete WY2001, valid numeric Q, TMIN containing the string `"1.0"` on every day, TMAX containing `"10.0"`.

   **Expected Julia:** `valid_temp_years=[]` because all temperatures are missing under its type rule.
   **Actual Python:** `valid_temp_years=[2001]`.

   **Fix:** Match Julia’s accepted input types before conversion, or explicitly change and validate the ingestion contract in both languages.

7. **NOTE — Both implementations allow infinite climate observations through their annual guards.**

   **Location:** [python/streamflow_signatures/climate.py:233](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/streamflow_signatures/climate.py:233), also line 352; PPT qualification at `python/streamflow_signatures/io.py:728`.

   **Julia/spec:** Temperature preprocessing correctly treats non-finite observations as missing. Annual guards, however, check only NaN. This is a shared robustness limitation, not a Python-only parity failure.

   **Minimal reproduction:** Complete WY2000, PPT zero except `+Inf` at dowy 21.

   **Actual in both:** The PPT year qualifies; annual precipitation and maxima are infinite, while `prcp_50_dowy=21` and wet-day counts are finite. Passing infinite temperatures directly also produces infinite temperature/GDD metrics instead of rejecting the annual row.

   **Expected robust behavior:** Treat non-finite observations as missing during qualification and reject any remaining non-finite annual input.

   **Fix:** Apply an explicit finite-value policy consistently to preprocessing and defensive annual guards. Document the direct API’s behavior.

## What I verified and found correct

- **Ordinary cross-language results:** The shipped fixture passed all 1,170 annual comparisons and 624 summary comparisons. My independent 46-year random fixture produced 1,794 bit-identical annual values after preprocessing. Its 624 summary comparisons passed `atol=1e-12`, `rtol=4e-14`; the largest absolute difference was `4.547473508864641e−13`.
- **Sequential accumulation:** For two separated blocks of fifty `0.1` values, the port returns total `9.99999999999998`, exactly equal to its final cumulative value, and D50=50. Python 3.12’s builtin sum instead gives `10.0`; the production climate summation helpers avoid it.
- **Uniform precipitation:** WY2001 at 2 mm/day gives annual precipitation 730, D50=183, centre of mass=183, five-day maximum=10, SD=0, Gini=0, and seasonality index `0.022374429223744285`. Uniform daily precipitation correctly does not imply equal monthly totals.
- **Single storm:** A 100 mm storm at dowy 100 gives Gini `364/365 = 0.9972602739726028`, SI approximately `11/6`, sample SD `100/√365`, and CV `√365`. The binding Gini formula correctly takes precedence over §A’s loose statement that a single wet day gives exactly 1.
- **Zero/all-dry year:** A zero-precipitation leap water year gives dry-spell length 366, start 1, end 366, and NaN cumulative timing, centre of mass, CV, Gini and SI.
- **Calendar handling:** Leap-year truncation is rejected. Seasons ignore a corrupted month column. For WY2000 with precipitation equal to calendar month, seasonal totals are winter 461, spring 368, summer 645 and fall 910.
- **Daymet omission:** Removing December 31, 2000 produces a rejected direct grid; preprocessing restores WY2001 dowy 92. With neighbouring PPT values 2 and 6, and TMIN values −4 and 2, the inserted values are respectively 4 and −1.
- **Frost/GDD boundaries:** June 30 in WY2000 is dowy 274 and is included; July 1 is excluded from last frost. Frost-free timing is NaN. Negative daily mean temperatures contribute zero GDD.
- **NA qualification:** Thirty isolated missing temperatures pass; 31 fail. Internal gaps of three pass; four fail. Boundary gaps fail. Both series must pass. Negative temperatures remain valid. A rejected Q year cannot become temperature-valid.
- **Record gate:** Eight qualifying frost years in the first ten-year window pass an 80% gate; seven fail. The latter suppresses trend statistics while retaining the mean.
- **Orchestrator/additivity:** Explicit input frames are required; precipitation runs without area normalization. Across the tested record, all 864 existing summary fields remained unchanged when temperature columns were added and when temperature calculations were enabled.
- **Runner:** Raw aliases normalize correctly. Executing the actual window-filter and temperature-frame statements for 1993–2025 retained only `[2025]` from temperature-valid years `[1992, 2025, 2026]`.

**Verdict: NO-GO for a validated Julia mirror.** The shipped-default interpolation failure alone breaks both exact annual parity and frost timing. The passing tests do not establish the required contract.