# Blind review (Sonnet 5.5) — rpkg port of the climate metrics (guidelines §3.9)

Reviewer inputs: `review_brief.md` (§A spec, §B decisions, §C global rules) and the code. I did not read the
implementer's narrative or any other reviewer's report. Julia (`julia/src/climate.jl`) is the reference; the
specification and decisions are the second reference. Review window 2026-10-12 12:20 – 2026-10-13 09:05 EDT.

## 0. Baseline, hashes, method

**Final code reviewed** (sha256, recorded 2026-10-13 09:01 EDT; unchanged at 09:05; the system-library build of
2026-10-13 08:58:25 EDT equals this source function for function, 33/33 functions compared with `deparse`):

| file | sha256 |
|---|---|
| `rpkg/R/climate.R` (mtime 10-12 15:22:28) | `e236f219fa41e6de534831f8c01277613fc2ad15f02504c7fff1868475c85d72` |
| `rpkg/R/io.R` (mtime 10-13 08:58:25, "round 5") | `9ca6882009370cf823993905a57027855936828ae12ec60f8315b02b1f57ef5d` |
| `rpkg/R/config.R` (10-12 15:22:28) | `2868a89ebb0217f62e88f64627accfeb15548006b24cc370e4950348c594c99a` |
| `rpkg/R/signatures.R` (10-12 12:23:34) | `32876bf809bb73c306ed48b88ffbe1b5ff8c21a3f27c577de1413890d329efb6` |
| `rpkg/R/stats.R` (10-12 12:51:25) | `b7e6f73a5b240dab623dee36c6df486889343feb4287b6e9dcd9b425dc33f4e7` |
| `rpkg/NAMESPACE` | `09b2a4b31f1e8771bcd403c0d9d49736852fe43d1307e3556f43814d38be95c9` |
| `docs/benchmarks/run_rpkg_benchmark.R` | `cd6be06fe6dcc1a55209a1bc51ea8f5a73a37698548b5f74fc64401075be40cb` |
| `rpkg/tests/testthat/test-climate_metrics.R` (10-13 08:58:25) | `7e61717d51619487b6974dbb39df0e30942dcc82b04c2dff2e211376cc20c5f8` |
| `rpkg/inst/config/signatures_config.json` = `config/signatures_config.json` (byte-identical, `cmp`) | `a9e8ed31fd3728a7ccaa6acfc3047dd357d479fb29243b8c685128f8bcca8dff` |
| Julia `climate.jl` (10-12 13:47:06) | `9fe0b8197b680f4ca252b21f9a60a6228d42853e1cbc6f5df391b7ab207fa9c1` |
| Julia `io.jl` / `config.jl` | `5fa81db78599c7d7f830656301aaad2da1317b5fc56575f83b7f758afa375aa7` / `f4f938f32850404d72fe452d464a9ae731cf734349b48a591269ba7a6f527caa` |

`climate.R`, `config.R`, `signatures.R`, `stats.R` are byte-identical to the build I tested from 15:22 on 10-12;
`io.R` changed once more (08:58, the whole-frame `.real_values` conversion at `io.R:238-245`) and I re-ran everything on it.

**The target moved.** R sources changed eight times (12:23, 12:27, 12:50, 12:51, 13:08, 13:13, 15:22, 08:58) and
Julia three times (12:47, 13:08, 13:47). Between ~12:47 and 13:08 the canonical Julia package did not precompile
(`using Printf` in `StreamflowSignatures.jl` without `Printf` in `julia/Project.toml`). To stay reproducible I
tested frozen private copies (Julia trees under `/tmp/rev/t2`, `/tmp/rev/t4`; R libraries `/tmp/rev/rlib*`), never the
system library, and re-verified at the end that the live files and the system build equal the last snapshot.
For all final-build numbers Julia was byte-identical to my reference tree `/tmp/rev/t4` (`cmp` on every `julia/src/*.jl` and `Project.toml`);
the pre-fix evidence in §2 was produced against the earlier Julia snapshots (12:47–13:08 sources).

**Method (all scratch, under `/tmp/rev`, nothing written to the repository except this report):**

1. Independent exact-rational reference (`ref_exact.py`, Python `fractions`, written from the spec/decisions, not from
   the code) for all 39 bases on the preprocessed frames: 13,595 annual cells.
2. Differential harness, installed R package vs canonical Julia, annual values compared **bit for bit**:
   `cases.csv` 16 gages (random, precipitation edge years, temperature edge years, NA policy, record-anchored gate
   variants, short records, full calendar, scattered NAs) = 13,595 cells; `cases2.csv` interpolation stress, 80 years
   with eight 2–3-day gaps per series per year = 3,120 cells; `cases3.csv` tie-rich quantised precipitation = 4,680
   cells; `cases4.csv` ±Inf/NaN in Q, PPT, TMIN, TMAX = 1,052 cells; direct-call fuzz, 400 water years with
   drop/duplicate/NaN/Inf/shuffle/truncate/signed-zero mutations = 15,600 cells.
3. Four config variants (dry 0/1 mm; wet [0.5, 12.5, 40] / [1e-5, 5]; window 3/7/366; fraction 0.25/0.9; frost −1.5/1;
   end month 12/2/10/9; GDD base 5/−5; gate on/off), Julia built per variant, R via `STREAMFLOW_SIGNATURES_CONFIG`.
4. The real runners end to end (`run_rpkg_benchmark.R` vs `run_julia_benchmark.jl`) on synthetic parquet and on 7 real
   gages cut from `daymet_1980_2023_rebuilt_10aug2026.parquet` (schema `site_id, Date, prcp, tmin, tmax, swe, vp, srad`)
   and `combined_streamflow_data_09feb2026.parquet`: 01011000 ME, 02296750 FL, 04085427 WI, 06614800 CO, 07263295 AR,
   09504000 AZ, 12010000 WA, WY 1981–2023, plus window/min-fraction gates and five degenerate-Daymet scenarios.
5. Orchestrator contract matrix (9 combinations), config-validation matrix (37 configs), additivity switch,
   TMIN/TMAX presence/absence, `R CMD check` on a copy, the full rpkg test suite, the author's cross-check script.

## 1. Findings still open on the final code

No BLOCKER and no MAJOR remains open. Severity words: MINOR = wrong or fragile but with no effect on the delivered
values today; NOTE = observation, decision or process item.

### 1. MINOR — the roxygen block of `preprocess_daily_data` is now attached to the new helper (latent NAMESPACE break)

- **Where:** `rpkg/R/io.R:96-126` (title, `@param`, `@return`, `@export` at :117, `@importFrom` at :118-119) is followed
  directly, with no blank line, by the helper's own roxygen text (:120-126) and `.fill_internal_gaps` (:127);
  `preprocess_daily_data` (:147) has no roxygen block any more.
- **Spec/Julia:** n/a (packaging). The exported API must keep `export(preprocess_daily_data)` (both runners call it).
- **Port does:** works today because `rpkg/NAMESPACE` is hand-synchronised. The next `devtools::document()` rewrites it.
- **Reproduction (copy of the final source, roxygen2 8.1.0):**
  `roxygen2::roxygenise("/tmp/rev/rox3/pkg", roclets = "namespace")` then diff against the repository NAMESPACE:
  `< export(preprocess_daily_data)` is lost, `> export(.fill_internal_gaps)` is added. Expected: unchanged. Actual: as stated.
- **Fix:** insert a blank line and move the helper (`io.R:120-145`) either above line 96 or below the closing brace of
  `preprocess_daily_data`, so each roxygen block precedes exactly one function. Re-run the NAMESPACE diff.

### 2. MINOR — an explicit JSON `null` in the climate / NA-policy sections is accepted by R and rejected by Julia

- **Where:** `rpkg/R/config.R:10-26` (`.cfg_int/.cfg_num/.cfg_bool` return the default for `NULL`), `:260-283` (climate
  section), and the raw-JSON screen `:287-310`, which only rejects arrays/booleans/strings (`is.list(rc[[key]])`).
- **Julia:** `julia/src/config.jl` strict readers: `ERROR: config frost_temp_c: expected a finite number, got nothing`
  (precompile fails). The 27 other invalid values I tried are rejected by R (§3); Julia's strict readers reject the same classes (by reading; `null` by running).
- **Reproduction:** config with `"climate": {..., "frost_temp_c": null, ...}`. Expected: load error in both.
  Actual: Julia fails; R loads with `frost_temp_c = 0` (also `gdd_base_c`, `frost_last_end_month`, ..., `temperature_na_policy`,
  `snow_na_policy` keys).
- **Fix:** in the raw screen add `for (key in names(rc)) if (is.null(rc[[key]])) stop(sprintf("climate.%s must not be null", key))`
  (`fromJSON(simplifyVector = FALSE)` keeps `null` entries in `names()`, checked), and the same for the two NA-policy lists.
  Invalid configs only; no effect on the shipped config.

### 3. MINOR — the shipped cross-check cannot see the defects found in this review, and its header and tolerance are stale

- **Where:** `docs/benchmarks/crosscheck_climate_fixture.jl:14-30` (periodic deterministic series: no gaps, no exact
  decimal tie, monotone/periodic annual series) and `crosscheck_climate_fixture.R:3-6,15,61` (`mk_rtol = 1e-6`, header
  says rpkg uses `Kendall::MannKendall`'s single-precision Fortran).
- **Evidence:** the fixture returned `CROSSCHECK PASS` against the 12:28 build while, on random data, that build had
  178 Spearman p-values (among 9,984 compared statistics) off by 1e-6 … 4.4 % (80.6 % at p ≈ 1e-7 on real Daymet), a Mann-Kendall tau off by
  1.07e-2 absolute (n = 20), `prcp_50_dowy` errors of 1 and 8 days on decimal ties (294 days in the minimal two-block case), and 212/163/75 last-bit
  differences in interpolated PPT/TMIN/TMAX. Since 12:51 `mann_kendall_test` is a double-precision port (final R vs Julia ≤ 1.1e-15), so
  `mk_rtol = 1e-6` would now hide a regression of 9 orders of magnitude.
- **Fix:** add to the fixture (i) a random-valued set with 2- and 3-day gaps in PPT/TMIN/TMAX, (ii) a 0.1 mm-quantised
  series with an exact 50 % tie (`P[1:7] <- 0.1; P[301:307] <- 0.1` is the minimal one), (iii) tie-free random annual
  series (exact vs approximate Spearman p-value matters only there); set `mk_rtol <- rtol`; rewrite the header. The generators
  in `/tmp/rev/gen_cases*.R` run in seconds and can be copied.

### 4. NOTE — definitional boundaries, identical in Julia and R (spec / decision level)

- `prcp_gini` (`climate.R:215`, `climate.jl:344`): for a constant series the sum cancels only up to rounding. 365 × 0.05 mm
  gives **−5.8667442997643451e-18** in both languages (1 of 577 Gini values in my suites), expected 0 ("0 is perfectly
  uniform"). The attainable maximum is (n−1)/n = 0.99726027397260275 for one wet day in 365 (the guideline text says
  "1 is all rain on a single day"). Both are consequences of the decided formula; clamp at 0 or document.
- `prcp_max_dowy` (`climate.R:195-197`) is 1 for an all-dry year (first of 365 equal maxima) while `prcp_50_dowy`, `prcp_cm_dowy`,
  `prcp_gini`, `prcp_seasonality_index`, `prcp_annual_cv` are NaN for ΣP ≤ 0 (edge_p WY2001). An arbitrary value then enters
  that metric's trend statistics. Negligible for basin means; say so in the docs or make it NaN for ΣP ≤ 0.
- Crossing tolerance edge: the rule admits a shortfall of 1e-12 × total below the target. For two storms of 1−d and 1+d mm
  the flip is at d = 2.0e-12 (cum/total = 0.5 − d/2), not 1e-12: d = 2.00e-12 → `prcp_50_dowy` 7, d = 2.01e-12 → 301, in both
  languages (the new boundary test uses 1e-12 and 4e-12, which sit safely on either side). Wording only.
- The Dec-31 filler (one interpolated day, dowy 92, in every water year after a Daymet leap calendar year, 11 of the 43 complete
  water years of 1980–2023) is a fabricated value counted against the 30-NA budget; it can turn a dry day wet and so split a dry spell.
  Identical in both languages (the 1-day fill is bit-exact under every formula). Already a known design choice; mention in the guideline relay.
- `fal` = Oct + Nov + Sep is not a contiguous season (decision #2/§C); `temp_*_fal`, `prcp_fal`, `gdd_cum_fal` inherit it.

### 5. NOTE — documentation and packaging state of the new exports

`R CMD check --no-tests` on a copy of the final source: 4 WARNING, 2 NOTE. Attributable to this change: no `man/*.Rd` for
`calculate_precipitation_metrics`, `calculate_temperature_metrics`, `precip_metrics`, `temp_metrics`; `temp_data` missing
from `man/calculate_all_signatures.Rd` (codoc mismatch); `Kendall` is an unused Import since `stats.R` stopped calling
`Kendall::MannKendall` (`DESCRIPTION:22`); the `@return` text at `signatures.R:30` still says ~1,653 elements; the
`preprocess_daily_data` docs omit `valid_temp_years`/`valid_swe_years`. Pre-existing: undocumented trend arguments in ~15 Rd
files, `high_na_denominator_columns`, `mann_kendall_test`, non-ASCII in `signatures.R`, `utils` unused. Regenerate the Rd
files after fixing finding 1.

### 6. NOTE — runner and process items

- Memory: `run_rpkg_benchmark.R:147-151` reads TMIN and TMAX on top of PPT/SWE: two more double columns on a ~98 M-row
  frame ≈ 1.6 GB resident, next to the 111 M-row streamflow table. Not measured at full scale (the Julia runner documents an
  OOM on the 16 GB machine at 8 columns); check before the full R run or run in gage batches.
- `.onLoad` reads `STREAMFLOW_SIGNATURES_CONFIG`, Julia reads `STREAMFLOW_CONFIG` (known deferred item). A climate-section
  variant experiment silently runs R with the shipped config unless both are set; I had to set them separately.
- Column order: the two runner CSVs differ in 1,941 of 1,978 column positions on identical data because R's `sort()` is
  locale-collated (`en_US.UTF-8`) and Julia's is byte order (pre-existing; includes the 624 new columns). Compare by name.
- Cost: precipitation family 72 ms + temperature family 53 ms per 45-year gage (≈ +18 % of the ~0.7 s per gage of the other families).

### 7. NOTE — what this review did not cover

All parity runs are on aarch64 macOS (R 4.6.1, Julia 1.12.6). There `capabilities("long.double")` is FALSE, so a stray `sum()`,
`mean()` or `cumsum()` on doubles (x87/long-double accumulation on x86-64) would be invisible empirically. I checked by reading:
`climate.R` contains only integer counts in `sum()` (:134, :136, :223) and every floating accumulation goes through
`seq_sum`/`seq_cumsum` (:36, :44). One x86-64 run of the battery is still owed. Not covered either: the ~6,000-gage scale,
real Daymet series with NaN gaps (the 7 real gages have none; gap handling was exercised synthetically), the Python port.

## 2. Findings I raised that the fixes resolved (dropped from the open list)

| # | Severity as found | Finding | Fix (where, when) | Verified on the final code |
|---|---|---|---|---|
| R1 | MAJOR | `prcp_50_dowy` used the strict `cum >= frac*total` while canonical Julia (13:08) uses `frac*total − 1e-12*total` (`climate.jl:82,336`). Minimal case: `P[1:7] <- 0.1; P[301:307] <- 0.1` gave **301** (strict; total 1.4000000000000001, 0.5·total 0.70000000000000007, cum[7] 0.69999999999999996) instead of 7. Decimal truth on a 0.01-quantised year (`gate_none` WY2015, 75,634 hundredths, exact half at dowy 199): R 200 vs Julia 199. Tie-rich set: 9 exact ties in 90 complete years, strict R wrong in 2 (WY2035 182 vs 181; WY2016 **145 vs 137**), tolerant rule right in 90/90 | `climate.R:88,210-211` (13:13) | 4,680/4,680 and 13,595/13,595 cells bit-identical; decimal truth 90/90; the boundary scan `d ∈ {0 … 4e-12, −1e-12, −4e-12}` identical in R and Julia |
| R2 | MINOR | Gap filling used `stats::approx`; on arm64 R's C code is FMA-contracted (`fma(r−l, k/span, l)`, matches 100 %), Julia is `l + (r−l)*k/span` evaluated left to right. 2-day gaps: only 69.2 % (k=1) / 56.5 % (k=2) of filled values bit-equal Julia's, 3-day gaps (k=3) 75.2 %; on x86-64 the unfused C formula still differs in 21–25 % of 2-day-gap values. Stress record: 212 PPT, 163 TMIN, 75 TMAX filled values differed → 5 of 3,120 annual cells (e.g. `gdd_cum_fal` 1022.1150000000004 vs 1022.1150000000001) | `.fill_internal_gaps` `io.R:127-145`, used for Q, PPT, SWE, TMIN/TMAX (13:08) | stress 3,120/3,120, all other suites 0 differences; filled frames bit-equal (cmp_frames) |
| R3 | MINOR as found, but repo-wide | Mann-Kendall through `Kendall::MannKendall` (single precision) and Spearman p-value through R's exact AS 89: 178 Spearman p-values (among 9,984 compared statistics) 1e-6 … 4.4 % off on tie-free series; real Daymet: Spearman p up to 80.6 % relative at p ≈ 1e-7, MK p up to 2.1 %; float32 near-ties (45.300000000000004 vs 45.299999999999997, the same decimal sum in another order) gave tau 0.0213 vs 0.0106 for n = 20, far above the documented ~1e-7 | double-precision port `stats.R:16-42`, `cor.test(..., exact = FALSE)` `stats.R:196-199` (12:51) | final R vs Julia on real data: `mk_rho` ≤ 5.0e-16, `mk_pval` ≤ 7.8e-16, `spearman_pval` ≤ 1.1e-15 absolute; synthetic ≤ 1e-12. **Consequence for the project:** every rpkg family's MK/Spearman values change relative to the 2026-08-27 rpkg reference CSV (re-baseline it; update `CROSS_LANGUAGE_STATUS.md` "irreducible library differences" and the brief's "~1e-7" statement) |
| R4 | NOTE (by code reading of the 12:28 build) | Season months hard-coded; `any(is.na(P))` instead of `all(isfinite)` (±Inf passed in R, NaN row in Julia); signed zero not normalised; wet-threshold label `gsub(".","p")` gave `prcp_ndays_abv_1e-05mm` vs Julia `1em05`; `as.integer(dowy)` truncation accepted fractional dowy; config readers truncated 5.9 → 5; no duplicate-column-name guard; temperature limits not overridable through `config` | `climate.R:28-32,54-58,96-100,181-183,297`, `config.R:10-26,40-52`, `io.R` (12:50-13:08, 15:22) | label/season/validation matrix §3; fuzz and non-finite suites identical |
| R5 | NOTE | Observation contract of Julia's `_real_vector` (strings, complex, factors → missing; Bool → value) and the SWE policy from `config` with typed readers; list-typed columns left a list in the cleaned frame | `config.R:31-38` `.real_values`, `io.R:238-245` (whole-frame conversion, round 5), SWE block `io.R` | character/factor PPT → NaN year, logical → 365, list → computed; `p$data$PPT` and `Q` are doubles; climate test blocks 5 and 6 pass |

## 3. What I verified and found correct

**Independent derivation (exact rationals, 13,595 cells, default config).** Every discrete cell (dowy, counts, run lengths,
NaN patterns) equals the reference; continuous cells agree to ≤ 1.5e-14 relative (`prcp_annual_sd/_cv` 1.46e-14, Gini 1.14e-14,
`prcp_cm_dowy` 9.3e-15, `prcp_annual` 9.3e-15, `gdd_cum_annual` 9.2e-15 — the price of the sequential-sum contract: 365 × 0.05 mm sums to
18.250000000000124). The one cell the strict-double reference rejects is the decimal-tie year of R1, where the code is the decimal-correct one. Under the four config
variants the reference also agrees (0 real mismatches).

**Hand-derivable values reproduced in both languages.** One 100 mm storm on dowy 100 of a 365-day year: Gini = 364/365 =
0.99726027397260275, seasonality index = 22/12 = 1.8333333333333337, `prcp_intermit_annual` 265 (dowy 101–365), 50 % day = centre
of mass = 100, `prcp_win` 100, `prcp_fal` 0. Uniform 2 mm in the 366-day WY2004: `prcp_annual` 732, `prcp_50_dowy` 183,
`prcp_cm_dowy` 183.5, Gini 0, seasonality index 7/366 = 0.0191256830601, `prcp_win` = `prcp_fal` = 182 (Dec+Jan+Feb = 91 d; Sep+Oct+Nov = 91 d).
Single 5 mm day on dowy 366 (leap): max/50 % day 366, dry spell 365 (1–365); on dowy 1: dry spell 364 (2–365). All-dry year:
`prcp_intermit_annual` 365 (1–365), five metrics NaN as specified. 0.05 mm every day: all days dry (≤ 0.1), spell 365.
Frost cut-off: Jun 30 = dowy 273 (non-leap) / 274 (leap), Jul 1 = 274 / 275; frost only on Jun 30 → first = last = 273/274; only on Jul 1
→ first = 274/275, last NaN; late-September frost only → last NaN; TMIN exactly 0.00 → no frost; Oct 1 frost → 1/1; Jun 29+Jun 30+Jul 1+Sep 30
→ first 272, last 273. Constant 10/20 °C: `gdd_cum_annual` = 15 × 365 = 5,475; TMIN −10, TMAX 10 with one 0 °C day: clamp gives 5. Ties: first
occurrence for `temp_min_annual_dowy` (40), `which.max/min`, longest dry run (first of two equal runs), 5-day windows inside the year.
End months 12 / 2 (leap and non-leap February) / 10 / 9 / 6 all match calendar arithmetic (the five I tested; other months use the same expression).

**R vs Julia, final build, annual values bit-identical:** 13,595 + 3,120 + 4,680 + 1,052 + 15,600 = **38,047 cells, 0 differences, 0
NaN-pattern differences**; valid-year lists (streamflow, PPT, temperature) identical for all 16 gages; summary statistics: 9,984 + 624 + 1,248 + 624 +
4,992 values, identical NaN patterns, worst relative difference 5.2e-14 (`pettitt_delta_mean`), 4.9e-14 (`linear_slp`), 3.9e-16 (`mean`); every
`senn_slp`, `spearman_*`, `mk_*`, `median`, `pettitt_*` within 1e-12 absolute. The physical gates agree: the record-anchored decade gate on
`frost_first_dowy/frost_last_dowy` (7/10 → skip, 8/10 → no skip, 9-year record → no gate, 10-year record → gate, 12 years, invalid middle years excluded
from the anchor), the 60 % / 80 % trend gates and the 20-value floor.

**Real Daymet, real runners.** 7 gages × WY 1981–2023 (Dec-31 filler present in 11 leap-year successors per gage): both runners qualify the same gages and emit the same
2,266 columns (624 climate); 11,232 climate annual rows bit-identical; 4,368 climate summary cells, NaN patterns identical, worst relative difference 6.4e-15
(`mean`); all 12 QA flags identical. Window `1993–2020` with `min_qualifying_fraction 0.6` selects the same 5 gages and gives 4,467/4,467 identical climate annual rows
(`valid_temp_years` window filter at `run_rpkg_benchmark.R:317-319` correct). Degenerate Daymet (5 gages): all-NaN Daymet rows, no Daymet rows, NaN precipitation
only, NaN temperature only, normal: non-NA cell counts and NaN patterns identical (0 vs 304/318 as expected; precipitation absent but temperature present and
vice versa follows the `gage_has_climate` / explicit `temp_data` contract). The lower-case Daymet names are renamed by the runner `nmap` (`:153-155`) and by the exported `read_parquet`.

**Orchestrator contract** (`signatures.R:196-217`; 9 combinations run in R, the 8 that exist in the Julia script match it except the one noted below): precipitation only on an explicit `climate_data` with
`has_climate`, never from `streamflow_data`; temperature only on an explicit `temp_data`; both independent of `area_normalized` (160 + 152 keys with 22 Q-PPT keys
when FALSE); a 0-row `temp_data` or `climate_data` yields the full NaN key set (304 / 320 keys with changepoint); a `temp_data` without TMAX yields nothing;
data.table and data.frame inputs identical. The one difference (0-row `climate_data`) is in the *existing* Q-PPT families only (R 49 vs Julia 90 keys); it is unreachable from the runner
(`climate_data` is built only when a PPT-valid year exists) and not part of this change.

**Existing families unchanged by this port's inputs:** with vs without TMIN/TMAX columns 21,184 non-climate values string-identical; with the `climate` section removed
1,332 keys, 0 climate keys, the other values identical (additivity switch works; `climate_enabled` false). The R-side behaviour changes for existing families come only
from R2/R3/R5 (interpolation formula, MK/Spearman port, ±Inf/strings missing), all toward Julia.

**Config.** 37 configs: 27 rejected at load (unsorted/duplicate/empty/zero/string/Boolean thresholds, window 1/367/5.5/"5"/[5], fraction 0/1/0.333, end month 0/13/6.5, negative dry
threshold, "yes"/0 for the gate, string GDD base, negative or fractional NA limits, season month 13/1.9/empty, near-duplicate labels), 10 loaded as intended (default; section absent
or `{}` → disabled; scalar threshold; `5.0`; `0.29`; absent gate/policy/season keys) plus the `null` divergence of finding 2. Config-derived names agree in both languages
(`prcp_3d_max`, `prcp_25_dowy`, `prcp_ndays_abv_0p5mm`, `…12p5mm`, `…1em05mm`, `prcp_366d_max`, `prcp_90_dowy`). Bundled rpkg config is byte-identical to the canonical one.

**Tests.** `test-climate_metrics.R`: 12 blocks, 0 failures; full rpkg suite on the final build: 98 blocks, 1,218 expectations, 0 failures, 0 errors; the author's
cross-check script passes. Each hand-derived claim in the tests holds (seasons, 273/274 cut-off, GDD clamp, the 7-day storm, the 1e-12 boundary). None of them contains a gap ≥ 2 days
or a decimal tie except the new ones, which is why R1/R2 survived until the differential harness.

## 4. Verdict

**GO** — use `rpkg` as a validated mirror of canonical Julia for the 39 climate metrics at the sources hashed in §0 (Julia 13:47:06, R as of 2026-10-13 08:58:25).
On the final build every annual value I could generate (38,047 synthetic cells and 11,232 real-Daymet cells, through the real runners) is bit-identical to Julia, summary
statistics agree to ≤ 5.2e-14, and the specification and decisions are met by an independent exact-arithmetic reference. The findings that mattered (R1–R3) were found
and fixed during the review; none is open. Conditions: (a) fix finding 1 before anyone runs `devtools::document()`; (b) any further edit to `climate.R`, `io.R`,
`config.R`, `stats.R` or `climate.jl` voids this verdict until the differential battery (cheap: ~10 minutes) is re-run — the R sources changed eight times in 20 hours;
(c) one x86-64 run of the battery (finding 7) and a memory check of the full-scale R runner (finding 6) before the production rerun; (d) re-baseline the rpkg reference CSV and
correct the "~1e-7 Mann-Kendall residual" and "irreducible Spearman" statements because R3 changed all rpkg families.

## Appendix — scratch harness (not in the repository)

`/tmp/rev/gen_cases.R, gen_cases2.R, gen_cases3.R, gen_cases4.R, gen_fuzz.R` (inputs, integer hundredths so both languages derive identical doubles),
`run_r.R`, `run_jl.jl` (preprocess → explicit `climate_data`/`temp_data` → orchestrator with collector, as the runners do), `compare.jl`, `cmp_frames.jl`,
`ref_exact.py` / `ref_exact2.py` (exact-rational reference, config-parametrised), `tie_truth.py`, `fuzz_r.R`, `fuzz_jl.jl`, `fuzz_cmp.jl`, `edge_r.R`, `orch_r.R`/`orch_jl.jl`,
`cfg_tests.py`, `boundary_r.R`/`boundary_jl.jl`, `runner/*.R`. They can be copied under `docs/benchmarks/` if wanted; I did not add them.
