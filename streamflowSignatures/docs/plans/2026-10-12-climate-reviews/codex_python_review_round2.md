**NO-GO.** Six original findings pass their original reproductions. Finding 7 is only partially fixed: the annual finite guards work, but precipitation preprocessing still treats infinity as an observation. Additional input-type and configuration boundary failures remain. Repository source files were not modified; reproductions and this report are under `/tmp`.

## Findings

1. **MAJOR — Non-finite preprocessing remains inconsistent with the binding convention; the new annual guard discards otherwise recoverable precipitation years.**

   **Location:** `python/streamflow_signatures/io.py:635`, `:725`, `:754`; `python/streamflow_signatures/climate.py:239`. Canonical counterparts: `julia/src/io.jl:296`, `:448` and the SWE conversion block.

   **Julia/spec:** The round-2 convention makes non-finite observations missing before qualification and interpolation. An isolated internal missing value between finite neighbours should be filled; four consecutive missing values or 31 missing values should fail the shipped policy. Remaining non-finite annual input must void the annual row.

   **Actual:** Python and Julia still use NaN-only checks for Q/PPT/SWE preprocessing. Temperature alone converts infinity to missing. The new precipitation kernel rejects the surviving infinity, losing a recoverable year rather than returning its interpolated metrics.

   **Minimal reproduction:** Complete WY2001; finite valid Q; PPT = 2 mm/day except dowy 21.

   | Replacement | Expected preprocessing and annual result | Actual in Python and Julia |
   |---|---|---|
   | `+Inf` | Fill to 2; `valid_climate_years=[2001]`; `prcp_annual=730` | Leaves `Inf`; qualifies `[2001]`; all 20 precipitation metrics NaN |
   | `-Inf` | Same result as `+Inf` | Rejects precipitation year as negative; `valid_climate_years=[]` |
   | Four consecutive `+Inf` values | Reject precipitation year | Qualifies `[2001]`; kernel subsequently voids it |
   | 31 isolated `+Inf` values | Reject precipitation year | Qualifies `[2001]`; kernel subsequently voids it |

   The same qualification defect exists for SWE. With 31 isolated infinite Q observations, both implementations return `valid_years=[2001]` and `valid_temp_years=[2001]`; both should be empty under the stated convention. This matters to the climate family's streamflow-valid-year prerequisite and to existing families consuming the supposedly qualified frames.

   **Proposed fix:** Convert non-finite Q/PPT/SWE observations to missing before raw counts, run lengths, sign checks and interpolation, in both languages. Keep the annual finite guards. Test preprocessing plus family calculations, not only direct kernel rejection.

2. **MINOR — The observation-type contract still differs from Julia, including a regression in the numeric-string fix.**

   **Location:** `python/streamflow_signatures/climate.py:234`, `:353`; `python/streamflow_signatures/io.py:789–795`. Canonical: `julia/src/stats.jl:614` (`coalesce_q`) and `julia/src/io.jl:580`.

   **Julia/spec:** Julia's annual conversion makes strings missing. Temperature preprocessing accepts finite `Real` values and treats non-real values as missing.

   **Minimal reproductions:** Complete WY2001 with valid Q and otherwise valid inputs.

   | Input and entry point | Expected / Julia | Actual Python |
   |---|---|---|
   | Direct kernels: PPT all `"2.0"`, TMIN all `"1.0"`, TMAX all `"10.0"` | All 20 precipitation and 19 temperature annual values NaN | `prcp_annual=730`, `temp_max_annual=10`, `gdd_cum_annual=2007.5` |
   | Preprocessing: TMIN all `1+2j`, TMAX all `10+2j` | `valid_temp_years=[]` | Discards imaginary components, qualifies `[2001]`, GDD = `2007.5` |
   | Preprocessing: TMIN all `Fraction(1,2)`, TMAX all `Fraction(21,2)`; corresponding Julia `Rational` values | Qualifies `[2001]` | `valid_temp_years=[]` |

   The original numeric-string preprocessing reproduction now passes. Its replacement type check, however, enumerates a few numeric classes instead of implementing the stated real-number contract. Meanwhile the direct kernels retain permissive string coercion.

   **Proposed fix:** Use one explicit real-observation conversion policy in preprocessing and both kernels: accept supported real numeric values, reject complex/string values, and map missing/non-finite observations consistently. Add direct-API and object-dtype cases.

3. **MINOR — Fractional season months still truncate silently instead of failing config validation.**

   **Location:** `python/streamflow_signatures/config.py:243`; `python/streamflow_signatures/climate.py:62`. Canonical: `julia/src/config.jl:347`.

   **Julia/spec:** Months are integral calendar months; strict validation must preserve their meaning. Julia rejects a nonintegral month during `Int` conversion.

   **Minimal reproduction:** Set `na_handling.seasonal_completeness.season_definitions.winter` to `[12,1,2.9]`; use complete WY2001 with PPT = 1 mm/day.

   **Expected:** Reject the configuration. Julia raises `InexactError: Int64(2.9)`.

   **Actual:** Python imports successfully, converts the climate season to `{12,1,2}`, and returns `prcp_win=90`. Preprocessing still reads the original list, so the two paths do not even use the same month set.

   **Proposed fix:** Validate season definitions at config loading: numeric integral months in 1–12, with an explicit boolean policy. Supply the validated definitions to both preprocessing and climate kernels; do not truncate in `_season_months`.

4. **MINOR — Whole-percentage validation accepts a different fraction under the unchanged D50 label.**

   **Location:** `python/streamflow_signatures/config.py:366`; `python/streamflow_signatures/climate.py:77`, `:271`. Canonical: `julia/src/config.jl:190`.

   **Julia/spec:** The configuration requires a whole percentage, and the derived column name must describe the fraction used in the calculation.

   **Minimal reproduction:** Set `cumulative_fraction=0.500000000001`. Complete WY2001, PPT = 1 at dowy 1 and 200, zero elsewhere; total = 2.

   **Expected:** Reject the non-whole percentage, or explicitly normalize it to exactly `0.5` if tolerance is intentional. For a normalized `prcp_50_dowy`, the answer is 1.

   **Actual in both languages:** Config accepts the value, names the column `prcp_50_dowy`, but retains the original fraction for calculation. The threshold exceeds 1, so the reported day is **200**. Rounding only the validation/name leaves different semantics behind the same column.

   **Proposed fix:** Normalize an accepted percentage to the canonical integer percentage divided by 100, or reject values differing from that canonical representation. Test exact cumulative ties after config loading.

5. **NOTE — Python's strict scalar checks now exceed the canonical implementation; the stated shared validation convention is not yet shared.**

   **Location:** `julia/src/config.jl:160–169`, `:363`; compare `python/streamflow_signatures/config.py:57–79`.

   **Spec:** Boolean, numeric and integer config values must retain their types rather than be silently coerced.

   **Minimal reproductions and actual results:**

   | Config value | Expected / actual Python | Actual Julia |
   |---|---|---|
   | `gdd_base_c: true` | Reject / `ValueError` | Accepts as `1.0` |
   | `record_anchored_decade_gate: 1` | Reject / `ValueError` | Accepts as `true` |
   | `max_raw_na_per_year_temp: true` | Reject / `ValueError` | Accepts as `1` |

   **Proposed fix:** Add strict type checks before Julia's `Float64`, `Int` and `Bool` conversions. Retain Python's stricter checks. This is a canonical follow-up, not a reason to loosen the port.

## Round-1 findings: original reproductions rerun

| Original finding | Result | Recomputed expected and actual values |
|---|---|---|
| **1. Temperature interpolation** | **PASS** | Both languages fill `[-5,NaN,NaN,NaN,1.7]` as `[-5,-3.325,-1.65,0.025000000000000355,1.7]`; annual GDD is exactly `1.7250000000000003`. The frost reproduction fills dowy 13 with `+0.0`; `frost_last_dowy=12` in both. |
| **2. Signed-zero extrema** | **PASS** | Original WY2000 PPT/TMAX `[-0,+0,…]` and TMIN `[+0,-0,+0,…]`: both languages now return `+0.0` and dowy **1** for all three extrema, directly and after preprocessing. Values checked as float bits: `0x0000000000000000`. |
| **3. Zero interpolation-gap setting** | **PASS** | One internal TMIN NaN and gap limit 0: both return normally with `valid_years=[2001]`, `valid_temp_years=[]`. Negative gap/raw-NA limits are rejected at config loading. |
| **4. Config validation** | **PASS for every original input** | `5.9` window, `6.9` frost month, `30.9` raw-NA limit and `3.9` gap limit: Python `ValueError`, Julia rejection. String `"false"`: rejected by both. Literal JSON `NaN` dry threshold and `Infinity` GDD base: rejected by both. No invalid annual dry-spell result is produced. New validation gaps are findings 3–5 above. |
| **5. Scientific-notation labels** | **PASS under the revised convention** | Threshold `0.00001` produces **`prcp_ndays_abv_1em05mm`** in both languages. The former Julia spelling is intentionally superseded. Additional Python labels: `1e-6 → 1em06`, `0.0001 → 0p0001`, `0.5 → 0p5`. |
| **6. Numeric-string temperature preprocessing** | **PASS for the original reproduction** | TMIN all `"1.0"`, TMAX all `"10.0"`, numeric valid Q: both return `valid_temp_years=[]`; Python retains `valid_years=[2001]`. Direct kernels and other numeric types remain inconsistent: finding 2. |
| **7. Non-finite observations** | **FAIL — partially fixed** | Original WY2000, zero PPT except `+Inf` at dowy 21: direct precipitation kernel now returns NaN for all **20** annual metrics; direct infinite-temperature input voids all **19**. But preprocessing still qualifies PPT `[2000]`, leaves `Inf`, and consequently loses the annual row instead of filling that isolated missing observation to zero. |

## What I verified and found correct

- **Tests:** All **14** climate tests pass. Another **78** tests pass across phase-1 regression fixes, area-normalization gating, schema validation, snow metrics and annual collection. These passing tests do not cover the failures above.
- **Exact annual parity on finite numeric fixtures:** Regenerated the shipped Julia fixture: **1,170 annual values** pass a `uint64` bit comparison, and **624 summary fields** pass `atol=1e-12, rtol=4e-14`. Reran the original independent 46-year random fixture: all **1,794** preprocessed annual values match exactly; every interpolated PPT/TMIN/TMAX double also matches. Largest absolute summary difference: **4.547473508864641e-13**, for `prcp_annual_pettitt_post_mean` (`2529.81216967731` Julia, `2529.8121696773105` Python).
- **Additional cross-language probes:** **10,453** annual rows, including invalid-year NaNs, had no bit/NaN-pattern mismatches; **19,968** summary comparisons had no tolerance failures. Cases included repaired interpolation, infinities, signed zeros, three-year monotonic records, tied annual series, nearly constant series and cancellation-prone Pettitt inputs. Shared specification failures still exist despite this parity.
- **Sequential accumulation:** Two separated blocks of fifty `0.1` values give total and final cumulative value **9.99999999999998**, D50 **50**. Python 3.12 builtin `sum` gives **10.0** and is not used for these climate reductions.
- **Closed-form precipitation:** WY2001 at 2 mm/day gives annual **730**, five-day maximum **10**, D50/centre of mass **183**, SD/Gini **0**, SI **0.022374429223744285**. One 100 mm storm gives Gini **0.9972602739726028**, SI **1.8333333333333337**, sample SD **5.234239225902137**, CV **19.1049731745428**, consistent with `364/365`, `11/6`, `100/sqrt(365)` and `sqrt(365)` within the sequential-double contract.
- **Calendars and grid guards:** WY1900/2001/2100 have **365** days; WY2000 has **366**. At 2 mm/day, totals are **730/732**, with centre of mass **183/183.5**. Shuffled complete grids reproduce the values; duplicate interior dowy and truncated grids void both families. A deliberately corrupted month column still yields WY2000 seasonal PPT totals **461, 368, 645, 910**.
- **Daymet:** Omitting December 31, 2000 is rejected by the direct grid guard. Preprocessing restores dowy **92**, PPT **4** between 2 and 6, and TMIN **-1** between -4 and 2.
- **Dry/frost/threshold boundaries:** Zero-PPT WY2000 has dry spell **366**, start **1**, end **366**, with NaN cumulative timing/centre/CV/Gini/SI. June 30 is included at dowy **274** in WY2000 and **273** in WY2001; July 1 is excluded from last frost. Exact thresholds and their next representable doubles give strict-above counts **5,3,1** at 10/25/50 mm in the tested fixture. A final five-day `[1,2,3,4,5]` window gives **15**.
- **Temperature policy and gates:** Thirty isolated missing values qualify; 31 fail. Three-day gaps qualify; four-day and boundary gaps fail. Negative temperatures remain valid; both series must pass; rejected Q years cannot qualify. The original seven/eight frost-year gate gives respectively slope **NaN/1.0**, retaining means **11.764705882352942/11.61111111111111**. The record-anchor helper correctly uses **1/2 versus 2/2** when only two anchor years occur in the first decade.
- **Orchestrator and existing families:** Explicit precipitation/temperature frames remain required. Temperature runs with `has_climate=False` and `area_normalized=False`: `temp_max_annual_mean=10`, no precipitation keys. All **864** previously shared summary fields in the original additivity fixture remain unchanged with temperature columns present and with temperature calculations enabled; temperature adds **152** fields when changepoints are disabled.
- **Runner:** Lowercase temperature aliases normalize correctly. Executing the runner's actual window-filter and temperature-frame statements for WY1993–2025 retains only **[2025]** from temperature-valid years `[1992,2025,2026]`.

**Verdict: NO-GO for using this port as a validated mirror of Julia.** Default finite numeric fixtures now match, but non-finite preprocessing violates the binding convention, public API type handling still differs from Julia, and accepted configuration boundaries remain inconsistent or mislabeled.
