# Sonnet 5.5 blind review — Python port of the climate metrics (guidelines §3.9)

Reviewer: Sonnet 5.5, adversarial, blind (review brief + code only). Date: 2026-10-12.
Verdict (final snapshot below): **GO** — no BLOCKER, no MAJOR open; four MINOR (F1–F4) and five NOTE (F5–F9) items.

## 0. What was reviewed, on which bytes, and how

The tree moved three times while I was reviewing it (12:19–12:50, 13:08–13:12, 13:47). I
reproduced defects on the early snapshots (§2, all since fixed) and re-ran every check on the
final one. **Everything in §1 and §3 refers to this snapshot** (newest mtime 13:47:06, unchanged for
the last 90 minutes of testing; sha256, first 12 hex digits):

| file | sha256 |
|---|---|
| `python/streamflow_signatures/climate.py` | 5f9051c16a45 |
| `python/streamflow_signatures/io.py` | d2b04b145e99 |
| `python/streamflow_signatures/config.py` | 0e3b87dab8bb |
| `python/streamflow_signatures/signatures.py` | 64ea98d33cc5 |
| `python/tests/test_climate_metrics.py` | 82f4eb5dff6a |
| `docs/benchmarks/run_python_benchmark.py` | 47d5b2a2d37c |
| `config/signatures_config.json` (= bundled copies, `cmp` identical) | a9e8ed31fd37 |
| canonical: `julia/src/climate.jl` / `io.jl` / `config.jl` / `signatures.jl` | 9fe0b8197b68 / 5fa81db78599 / f4f938f32850 / e4d2bfe8029a |

Method (all scratch work under the session scratchpad; no repository file other than this report was touched):

1. **Bit-exact differential test against canonical Julia.** Inputs are hex-encoded IEEE doubles so no
   CSV parser can move a bit. A private snapshot copy of `julia/` (the repo's `Project.toml` lacked
   `Printf` between 12:47 and 13:08) was re-taken whenever the tree changed. Per case: preprocessor
   output (`valid_*_years`, every cleaned Q/PPT/TMIN/TMAX value), the annual values through the
   `AnnualCollector`, and all 624 summary values with the production gates (60/80 %, floor 20,
   Pettitt 1980–2024).
   Cases: 27 hand-built (leap/non-leap, Daymet 365-day calendar, ties, all-dry and constant years, subnormals,
   threshold boundaries, frost on Jun 30 / Jul 1 / Sep 30 / Oct 1, ±0.0, ±Inf, NaN runs of 1–5 days, exactly 30 and 31
   scattered NaN, decade-gate 7/8/9 of 10, record spans 9/10/12), 150 random "pre" fuzz cases, 60 random "direct-kernel" fuzz
   cases (shuffled, dropped and duplicated rows, Inf), and three config variants (window 7/3/366, fraction 25/80/99 %, wet thresholds
   incl. 0.5 / 12.5 / 100 mm, dry 0/1 mm, frost +2/−2.5 °C, cut-off month 12/3/9, GDD base 5/−1/10 °C, gate off) over 12 + 80 cases each,
   plus 80 cases with the stats floor removed.
2. **Independent naive reference** written from the specification text (pure Python, `math.fsum`, exact `Fraction`
   arithmetic for the 50 % crossing, double-sum definition of the Gini, date-based season and Jun 30 logic).
3. **Real Daymet data** (reprocessed 1980–2025 input, batch b1): 65 selected gages (frost-free, arid, wet, Canadian
   un-normalized, shortest records) with BOTH runners from the final tree, and the whole b1 batch (824 gages) with the Python
   runner against the Julia product the implementer produced (`additivity/with_climate`).
4. Config-validation matrix (43 variants, Python vs Julia), orchestrator-contract matrix (7 call patterns), additivity
   (TMIN/TMAX present vs absent, synthetic and real), a HEAD-vs-current comparison of the existing families (155 cases),
   calendar brute force (water years 1901–2200), and the unit suite.

## 1. Findings on the final snapshot

Ordered by severity. None changes a single annual value or summary value in any test I ran.

### F1 — MINOR — PPT and SWE blocks write the interpolated series back before the validity decision (frame differs from Julia in invalid years)
- **Where:** `python/streamflow_signatures/io.py:762-766` (PPT) and `:793-797` (SWE).
- **Julia / spec:** `io.jl` interpolates into a local vector and assigns `joined[!, :PPT] = P` / `:SWE = S` only when the year is valid. The new temperature block of the port (`io.py:810-826`) already does the same.
- **Port:** `yr_data["PPT"] = _interpolate_internal_gaps(...)` runs first and the residual-NaN rejection follows, so a year rejected for a leading/trailing NaN comes back with its interior gaps (and Daymet's absent Dec 31) filled.
- **Repro:** one water year, PPT = 2.0 mm every day except dowy 1 = NaN (year invalid) and dowy 49/50/51 = 1.0/NaN/3.0.
  Expected (Julia): `valid_climate_years = []`, PPT[dowy 50] = NaN. Actual (Python): `valid_climate_years = []`, PPT[dowy 50] = 2.0.
  Fuzz: 36 of 150 cases, 285 rows, none in a PPT-valid year; Q, TMIN, TMAX identical in all rows.
- **Impact:** no signature (every consumer filters to valid years; `swe_max_to_ppt` guards with `valid_climate_years`), but the returned `data` frames are not interchangeable between languages.
- **Fix:** compute `filled = _interpolate_internal_gaps(...)`, test `np.isnan(filled).any()`, assign to `yr_data` only when the year is valid (PPT and SWE).

### F2 — MINOR — the new public API is not exported from the package root
- **Where:** `python/streamflow_signatures/__init__.py:35-36, 77-83` (unchanged by this work).
- **Julia:** `StreamflowSignatures.jl` exports `calculate_precipitation_metrics`, `calculate_temperature_metrics`, `PRECIP_METRICS`, `TEMP_METRICS`; Python exports the snow and drought equivalents (`calculate_snow_metrics`, `SNOW_METRICS`, `calculate_drought_metrics`, …).
- **Repro:** `from streamflow_signatures import calculate_precipitation_metrics` → `ImportError: cannot import name ...` (works only as `from streamflow_signatures.climate import ...`).
- **Fix:** import the four names in `__init__.py` and add them to `__all__` under a "Climate (Daymet prcp/tmin/tmax)" group.

### F3 — MINOR — the round-3 strict validators reject NumPy integer scalars
- **Where:** `python/streamflow_signatures/io.py:593-610`.
- **Julia:** the keyword limits are typed `::Int`: an `Int64` is accepted, `Int32` and `3.0` raise a `TypeError`, a negative value an `ArgumentError` (checked).
- **Port:** `isinstance(_v, (int, float))` is false for `np.int64` / `np.int32` (they do not subclass `int`), so the exact analogue of Julia's accepted `Int64` is rejected, while `3.0` is accepted — a caller that takes the limits from a DataFrame or NumPy array fails. Flags given as `0/1` are rejected too (intended).
- **Repro:** `preprocess_daily_data(df, config={"max_gap_days": np.int64(3)})` → `ValueError: max_gap_days must be a whole number >= 0, got np.int64(3)`; `np.float64(3.0)`, `3`, `3.0` and `np.bool_(False)` pass. The shipped runner is unaffected (JSON ints).
- **Fix:** test with `numbers.Integral` / `numbers.Real` (NumPy registers its scalars with both) before the integrality check.

### F4 — MINOR (pre-existing, newly visible) — PPT/SWE residual-NaN rejection is conditional on `reject_residual_na`; Julia's is unconditional
- **Where:** `io.py:765` and `io.py:796` (`if reject_residual and ... > 0`); the port's temperature block (`io.py:820-822`) is unconditional, as in Julia.
- **Julia:** `io.jl` PPT/SWE/temperature blocks always reject a residual NaN; only the Q block honours `reject_residual_na`.
- **Repro:** WY 1990 with PPT NaN on dowy 1 and TMIN NaN on dowy 365, `config={"reject_residual_na": False}`: Python `valid_climate_years` contains 1990 and `valid_temp_years` does not (the two blocks disagree); Julia `preprocess_daily_data(...; reject_residual_na=false)` excludes 1990 from both.
- **Impact:** none with the shipped config (`reject_residual_na: true`); with the flag off the port would pass a year with a NaN into `prcp_*`, where the finite guard then makes the whole row NaN.
- **Fix:** drop the `reject_residual and` conjunct for PPT and SWE.

### F5 — NOTE — the 50 % crossing now carries a 1e-12 relative tolerance that is not in the specification or the decisions table
- **Where:** `climate.py:121-123` and `:288-291` (`thr = frac*total − 1e-12*total`); `julia/src/climate.jl` identical.
- **Spec / decisions:** "first dowy with cumulative P ≥ 0.5 ΣP" (brief §B). The tolerance is documented in `docs/SIGNATURES.md` (§ caveats) but not in the guidelines text or the decisions list.
- **Evidence it is harmless, and in one respect beneficial:** (a) P = 0.1 mm on every day of a 366-day year: exact arithmetic says day 183, the strict float test said 184, the tolerant test gives 183. On decimal-quantized series (0.1 mm / integer mm — my cases `c4_one_decimal`, `t2_boundaries`) the cumulative sum hits the half **exactly in decimal arithmetic** in 2 of the 60 water years of those two cases (WY 1999 → day 175, WY 1996 → day 181); the doubles are 1.6e-18 and 4.6e-19 (relative) below the half, so an exact-over-doubles reference says 176 / 182 and the tolerant implementation correctly says 175 / 181 (decimal-exact). Without the tolerance those years would be decided by representation error; (b) boundary identical in both languages: storms A (dowy 10) and B (dowy 200) with cum(A) = T/2 − δT cross at 10 for δ ≤ 9e-13 and at 200 for δ ≥ 1.1e-12; (c) on the real b1 batch (47,047 gage-years) the closest approach of a cumulative value to 50 % from below is 2.5e-7 of the annual total, so 0 years change.
- **Residual:** a constant subnormal series (5e-324 mm/day) still gives day 182 vs exact 183 because 1e-12·T underflows — irrelevant physically. The metric is intrinsically discontinuous at ties (δ from 9e-13 to 1.1e-12 moves the date by 190 days).
- **Fix / ask:** record it as a user decision (it silently refines a stated rule), and add a boundary test (the δ = 5e-13 / 2e-12 pair above) to the three unit suites; the current tests pin only the two-equal-blocks case.

### F6 — NOTE — runner memory: two more float64 columns, no batching in the Python runner
- **Where:** `docs/benchmarks/run_python_benchmark.py:81-83, 236-248`.
- **Measured:** the projected climate frame goes from 39.8 to 55.8 bytes/row (b1: 669 → 938 MB); for the 134.7 M-row reprocessed input that is 5.4 → 7.5 GB before `climate.groupby(...)` duplicates it into the per-gage dict (peak ≈ 2×). The docstring already records that the 8-column frame OOMs on 16 GB machines; a full-scale Python validation of the 1980–2025 input needs the batch inputs (`inputs/*_b1…b8`) or a bigger machine. Wall time was 1.0 gage/s on b1 (824 gages in ≈ 14 min), so batching is affordable.
- **Fix:** drive the Python runner per batch (the Julia batch harness' `batches.json`), or read the climate parquet by gage group instead of holding it whole; no code defect.

### F7 — NOTE — tests that guard the existing families around the interpolation refactor are skipped here, and the cross-language fixture cannot see this class of defect
- 9 tests in `tests/test_against_golden.py` skip (`test-data/sample_streamflow.parquet` missing), so the refactor that replaced `Series.interpolate` by `_interpolate_internal_gaps` for Q/PPT/SWE has no golden regression on this machine. Compensating evidence (measured on the 13:08 snapshot; the 13:47 edits do not touch interpolation): HEAD vs current Python on 155 cases (random Q gaps of 1–5 days): 11 cases change Q in the last bit (3-day gaps), 2 cases change 3 summary cells (≤ 1.0e-15 relative, `drought_deficit_fixed_p20_pettitt_*`) and 3 annual values; no valid-year list changes, no NaN flips.
- `docs/benchmarks/crosscheck_climate_fixture.*` builds 30 gap-free deterministic years (no NaN, no absent Dec 31, no ties), so it cannot detect interpolation, tie, signed-zero or gap-policy defects — which are exactly the classes that were found.
- **Fix:** give the fixture a few leap-year Dec 31 gaps (1-, 2- and 3-day), a tie year and a frost-free stretch; make the golden tests self-contained or fail loudly when their data are missing.

### F8 — NOTE — stale documentation and docstrings
`io.py:6-8, 48, 188, 338` still describe the 1,653-column product; `calculate_all_signatures`' `temp_data` docstring (`signatures.py:106-114`) says precipitation runs "whenever it carries PPT" but the code also requires `has_climate=True` (same as Julia); `python/README.md` has no `temp_data`/climate-metrics text and the root `README.md:72` ("the full 1,653-column summary additionally requires … snow_data … has_climate = true") omits `climate_data` and `temp_data`.
- **Fix:** update the four counts to 2,277 / 139, state in the `temp_data` docstring that precipitation needs `has_climate=True` and an explicit `climate_data`, and add `temp_data` / the climate section to `python/README.md` and the README recipe.

### F9 — NOTE — shared (Julia = Python) numerical / design properties worth a line in the co-author relay
- `prcp_gini` is bounded by (n−1)/n (0.99726 for one wet day of 365), not 1 as the guideline sentence says; implemented as documented.
- An all-dry year returns `prcp_max_dowy = 1` (first of 365 tied zeros) and `prcp_intermit_annual = 365`; `prcp_cm_dowy`, `prcp_50_dowy`, `prcp_gini`, `prcp_seasonality_index`, `prcp_annual_cv` are NaN. A max day for a zero maximum is not meaningful; decision #9 (first tie) makes it 1.
- Summary comparisons need an absolute floor: Julia's Spearman p-value is `2·(1 − cdf)`, which loses relative accuracy below ~1e-7 (observed 1.3e-9 relative on p ≈ 5e-8, absolute 6e-17), and a linear slope of a constant non-zero series is ±1e-18 (`prcp_ndays_abv_50mm_linear_slp`: 1.3e-18 vs 0.0 — "relative difference 1"). With a 1e-12 absolute floor the worst real-data difference is 5.2e-14.
- **Fix / ask:** none for the code; relay the first two bullets to the co-authors with the other §3.9 notes.

## 2. Defects I reproduced on earlier snapshots and verified fixed on the final one

These existed in the version first handed to me (12:19) and are the reason the early differential runs failed; they are listed so the fixes are credited to a second, independent reproduction.

| # | Sev. at the time | Reproduction (earlier snapshot) | State on final snapshot |
|---|---|---|---|
| R1 | BLOCKER for the bit-identity contract | TMIN/TMAX (and Q/PPT/SWE) used `Series.interpolate`; NumPy's `interp` kernel contracts `slope*(x−xj)+yj` to an FMA on arm64, so the **third day of a 3-day gap** differs from Julia's `left + (right−left)*m/span` in 27 % of draws (18,015 / 66,861; all at m = 3, 100 % explained by FMA). Case `i1_gap3`: 21 TMAX, 26 TMIN, 52 PPT cleaned values and 4 of 1,170 annual values differed. | fixed by `_interpolate_internal_gaps` (Julia's expression and order). `i1_gap3`, 3-day gaps in the 150 fuzz cases: 0 differences. |
| R2 | MAJOR | signed zeros: Julia `argmax([-0.0, 0.0]) = 2` / `argmin([0.0, -0.0]) = 2`, NumPy picks index 0; Julia `unique` separates ±0 so Mann-Kendall ties differ. Case `f1_frost_patterns` (WY 2005/2006 with TN = 0.0 / −0.0): `temp_min_annual_mk_rho` −0.33451 vs −0.33426 (7e-4 relative), `temp_min_sum_pettitt_post_mk_pval` 4.9e-3 relative. | fixed by `+ 0.0` before the kernels in both languages; `f2_frost_nosz` and the signed-zero fuzz: identical. |
| R3 | MAJOR | guard tested `isnan` only: +Inf in P made `prcp_*` Inf/NaN instead of a NaN row; engine-level Pettitt p-values NaN vs number (`x3_direct_nonfinite`). | fixed (`isfinite` guard; non-finite = missing in the preprocessor). |
| R4 | MAJOR | `prcp_50_dowy` decided by the last bit of the sums: P = 0.1 mm on a 366-day year → day 184 although the exact crossing is 183 (8 leap years in case `c3_tenth`; also `z1_zero_years` WY 1996). | tolerance added in both languages (F5). |
| R5 | MINOR (code reading, not executed) | wet-threshold label: Python `str(1e-5)` = `1e-05`, Julia `string(1.0e-5)` = `1.0e-5` → different column names for exponent-form thresholds. | shared `%g` formatter + uniqueness check; labels identical in both languages for 0.5, 12.5, 1e-5 (`1em05`), 1e15 (`1e+15`) and the colliding pair 10 / 10.0000001 is rejected in both. |
| R6 | MINOR | observed on the 12:50 snapshot: `record_anchored_decade_gate = 1` accepted by Julia, rejected by Python; by code reading of the 12:19 snapshot: `int(5.5)` = 5, string `"false"` → `True`, JSON `NaN` accepted. | strict readers in both; OK/ERR and derived names identical on all 43 variants (§3). |
| R7 | MINOR (code reading) | `pd.to_numeric(..., errors="coerce")` turned numeric strings into observations; Julia treats non-`Real` as missing. | `_real_values` (one observation contract; Arrow/decimal columns go through the scalar path). |

## 3. What I verified and found correct (final snapshot)

**Bit-identity with canonical Julia**
- 593 case-runs (27 base + 150 fuzz + 60 direct fuzz + 80 floor-free + 3 × 12 variant + 3 × 80 variant fuzz): **410,156 annual rows, 0 non-identical or missing** (380,589 finite); 368,560 summary values: 0 key-set mismatches, 0 NaN-pattern mismatches, 0 Pettitt `cp_year` flips in 6,277 finite pairs; worst relative difference with a 1e-12 absolute floor 2.0e-12 (`gdd_cum_annual_pettitt_delta_mean`, a difference of two means), everything else ≤ 4.2e-16. Cleaned Q, TMIN, TMAX columns identical in every row of every case (PPT differs only in invalid years, F1).
- Real data, both runners from this tree, 65 gages: 266,254 annual rows (139 signatures) merge 1:1; **75,254 finite climate annual pairs bit-identical**, 0 NaN-pattern mismatches; 2,276 CSV columns in identical order; 624 climate summary columns ≤ 1.0e-14 (gdd Pettitt delta); QA flags identical; no family failure logged.
- Real data, whole batch b1 (Python at the 12:50 snapshot vs the implementer's Julia product): 3,415,156 annual rows, 969,384 climate rows, **968,666 finite pairs bit-identical**, 0 NaN mismatches, 0 Pettitt flips in 31,673 pairs, 624 summary columns ≤ 5.2e-14 (abs floor 1e-12), QA flags identical. (The later edits are small; the 65-gage re-run above covers the final tree.)
- Real data, mid-record window (`STREAMFLOW_START/END_WATER_YEAR = 1999/2023`, 60 % rule; the preprocessor sees the full record, the result is windowed), both runners: 56 gages qualify in both, 186,747 annual rows merge 1:1, 52,816 finite climate pairs bit-identical, 0 NaN mismatches, 624 summary columns ≤ 8.2e-15, QA flags identical — so `valid_temp_years`, `temp_data` and the window filter are wired like Julia's.
- Additivity on real data (Python, final tree, same 65 gages with and without `tmin`/`tmax` in the climate input): 2,276 vs 1,972 columns; the 304 extra are exactly the temperature family; **all 1,972 shared columns are bitwise equal** (every Q-only family, snow, runoff ratios, elasticity, storage, drought, QA flags, metadata).

**Independent reference (naive, from the spec; run on the final tree)**: 15 cases, 20,054 annual values: dowy, counts and dry-spell start/end agree exactly except two `prcp_50_dowy` values that are explained in F5 (decimal ties: the port equals the decimal-exact day, the exact-over-doubles reference is one wet day later); continuous values ≤ 7e-15 relative (sequential vs `fsum`); the sorted-array Gini equals the double-sum definition Σ|yi−yj|/(2n²ȳ) to 1e-12; the 50 % day equals the exact-rational crossing in every other year except the subnormal series; the Jun 30 cut-off equals the date-based rule in leap and non-leap years (WY 2001 frost only on Jun 30 → 273; WY 2004 Jun 30 + Jul 1 → 274, Jul 1 excluded; WY 2000 Jul 1 only → first 275, last NaN; WY 2008 Sep 30 only → first 366, last NaN).

**Calendar and boundaries**: water years 1901–2200 (incl. 2000, 2100, 2200): day counts, month-of-day, complete-grid guard (full grid accepted, grid short by one rejected) and the cut-off dowy for all 12 possible `frost_last_end_month` values agree with an independent day-by-day count. Daymet's calendar confirmed on the real input: Dec 31 is absent in 2020 and 2024, Feb 29 present (12 of 45 climate years of the 1980–2025 product, 9 of 33 of 1993–2025, as docs state). Decade gate: 8/10 passes, 7/10 gates, 9/10 passes, record span 9 → no gate, span 10 → gate active — identical in both languages with and without the stats floor.

**Sequential-summation contract**: no builtin `sum`, `np.sum`, `np.cumsum`, `np.mean` on any float series in `climate.py` (only `_seqsum` / `_seqcumsum`; integer counts use `.sum()` on booleans). Python 3.12.13: builtin `sum([0.1]*10) == 1.0` (compensated) while `_seqsum` gives 0.9999999999999999 as in Julia. Operations that could differ by platform (`dowy*P`, `(2i−n−1)*y`, `(P−μ)²`, `np.abs`) are single element-wise IEEE operations; no FMA can form across them, so annual values do not depend on the BLAS/SIMD build.

**Config validation and column names**: 43 variants (negative/zero/string/duplicate/unsorted thresholds, window 1/2/366/367/5.5/"5", fractions 0/1/0.25/0.255/0.07/0.29/0.995/0.5000000004, end month 0/1/6.0/6.5/12/13, gate "false"/1, negative and fractional NA limits, absent/empty `climate`, absent `temperature_na_policy`): identical OK/ERR outcome and identical derived column names in Python and Julia for every variant.

**Orchestrator contract** (measured on the 12:50 snapshot; `signatures.py`/`signatures.jl` are unchanged since 12:19; 7 call patterns, counts 160 precipitation / 152 temperature keys with Pettitt off, identical in both languages): no `climate_data` ⇒ no `prcp_*` (no fallback to `gage_data`); `climate_data` with `has_climate=False` ⇒ none; `temp_data` alone ⇒ 152 temperature keys; `area_normalized=False` ⇒ both families present, Q-to-PPT absent; frame without PPT / without TMAX ⇒ family absent; 0-row `climate_data` ⇒ 160 NaN-stat keys; 0-row `temp_data` ⇒ 152 NaN keys.

**Existing families**: with vs without TMIN/TMAX the Python outputs are bitwise equal on all 1,972 shared summary columns of the 65 real gages (snow, QA flags and metadata included) and (12:50 snapshot) on all shared keys and annual rows of 5 synthetic cases; the +304 temperature columns are the only addition. `flagged_for_high_na` did not change in either test (the denominator grows by 624 columns by design).

**Runner**: window filter carries `valid_temp_years`; `temp_data` is built only when both columns exist; column projection accepts `tmin/tmax/TMIN/TMAX`; sorted column order identical to Julia's CSV; annual parquet has 139 signatures.

**Unit suite**: `pytest python/tests` 175 passed, 9 skipped (F7); `tests/test_climate_metrics.py` 15 passed also with `FutureWarning` and `DeprecationWarning` turned into errors. Claims I tried to falsify: no assertion is vacuous or circular; `test_record_anchored_gate` genuinely discriminates the anchored gate (the metric's own 10-year span is 100 % complete, so only the record-anchored rule can suppress the trends); `test_round2_parity_regressions` pins the 3-day-gap expression order. Not pinned by any test: the 1e-12 boundary (F5), NumPy-integer limits (F3), PPT/SWE write-back parity (F1), the `reject_residual_na=False` path (F4). Closed-form expectations in the tests re-derived by hand and correct: one 100 mm storm on dowy 100 (Jan 8) → Gini 364/365 = 0.997260, SI = 11/6 = 1.833333 (the spec's "1.83" maximum), longest dry run 265 days (dowy 101–365), `prcp_win` = 100, other seasons 0; 2 mm every day of a non-leap year → total 730, `prcp_50_dowy` 183 (cum 364 < 365 ≤ 366), centre of mass 183, Gini < 1e-12, `prcp_5d_max` 10; P = 0.1 / 0.1000001 / 0.0 / 0.1 on dowy 10–13 → dry run 12–13 (the 0.1000001 day is wet), counts above 10/25/50 mm of 4/2/1 for 10.0 / 10.0001 / 25 / 50 / 50.5 mm; TN 10 °C, TX 30 °C all year → GDD 20 × 365 = 7,300 and `gdd_cum_win` = 20 × 90; TN = month − 30, TX = month → `temp_max_annual` 12 on dowy 62 (Dec 1), `temp_min_annual` −29 on dowy 93 (Jan 1).

## 4. Verdict

**GO** for using the Python port, as of the snapshot hashed in §0, as a validated mirror of canonical Julia: the annual values are bit-identical on every input I could construct (593 synthetic case-runs, 65 real gages with both final runners, 824 real gages against the Julia product) and the 624 summary columns agree to the library tolerance with identical NaN patterns, key sets and column order. The specification and the 2026-10-12 decisions are met (definitions, dowy numbering, season membership, tie rules, Jun 30 cut-off, record-anchored gate, temperature NA policy, orchestrator contract).

Conditions for the freeze (none blocks the verdict): fix F1–F3 (a few lines each) or record them as accepted; record F5 as a user decision and pin its boundary in the three test suites; re-hash the §0 files before the production run — any later edit voids the numbers above.
