# Plan: Streamflow Drought Signature Family (`drought_duration`, `drought_deficit`)

**Status**: IMPLEMENTED + UNIT-TESTED (Julia), 2026-07-27 — all design decisions
resolved; `julia/test/runtests.jl` green at **2,042 assertions, 0 failures** (drought
file: 670). Outstanding: `smoke_test.jl` (needs the `D:\` inputs), Codex review,
benchmark + additivity gate. See §14.

**Requested by**: user, 2026-07-27.
**Source**: Adelsperger et al. (in review), "A novel severity-based approach for
assessing streamflow drought characteristics and drivers". Verbatim method text in §2.0.
**Language**: Julia-first (canonical), per CLAUDE.md. Python/rpkg ports deferred to the
existing port queue (annual-values collector, b=1 recession alpha, snow metrics — this
becomes the fourth item).

---

## 0. Summary

Per gage, per water year, per threshold level:

| Metric | Definition | Units |
|---|---|---|
| `drought_duration_p{n}` | total number of days with 7-day-smoothed Q below the threshold | days |
| `drought_deficit_p{n}` | sum of flow departures below the threshold | mm |

**Scope as decided**: the **fixed** (whole-record) threshold method only — one threshold
per gage per level — at **all five levels (2, 5, 10, 20, 30 %)**. The paper's variable
(day-of-year) method is deliberately **out of scope** (§2.5, §13). That yields
**10 metric bases** (2 metrics × 5 levels) × (8 statistics + 8 Pettitt fields) =
**+165 columns** (160 stat/Pettitt fields + 5 threshold scalars; 1,488 → **1,653**) and
90 → **100** signature bases in the annual parquet (scalars are not annual series).

Both metrics are **dense annual series with meaningful zeros** (a wet year has 0 drought
days / 0 deficit), so they behave well under the existing trend gates and the 20-value
stats floor, and need no snow-style record-anchored gate (§3.5). Q-only — no climate
dependency, so no `area_normalized` Q-to-PPT gate (but see §3.4 for the deficit units
caveat).

### Resolved decisions (user, 2026-07-27)

| Decision | Resolution |
|---|---|
| Threshold method | **Fixed only** (whole-record). Variable/day-of-year dropped — it also removes the low-tail sample-size problem of §2.5 |
| Threshold levels | **All five**: 2, 5, 10, 20, 30 % (USDM D4→D0 ladder) |
| Annual grouping | **Water year (Oct 1 – Sep 30)** — framework-consistent; documented deviation from the paper's climate year (§3.2) |
| Below-plotting-range policy | **NaN** when `p < 1/(n+1)` — retained as a defensive guard; cannot fire for the fixed method (§2.5) |
| Smoothing alignment | **Centered** 7-day window |
| Column naming | **Explicit method infix**: `drought_duration_fixed_p10` (user, 2026-07-27) — a future variable method adds `_var_` non-breakingly |
| Threshold diagnostics | **Included**: five `drought_threshold_fixed_p{n}` per-gage scalars (§3.6) |
| Deficit reference series | **Smoothed** Q (the same series used for identification) — the internally consistent reading; confirmed with the user 2026-07-27 |

---

## 2. Method specification

### 2.0 Source text (verbatim)

> We converted 7-day smoothed daily observed streamflow to percentile values to enable
> drought identification via consistent percentile-based thresholds between streamgages.
> Streamflow percentiles were computed with the unbiased Weibull plotting position (e.g.,
> Laaha et al., 2017) for each smoothed daily flow using two approaches: (a) fixed—all
> flows in the period of record were used to calculate one fixed threshold (b) variable -
> unique thresholds were calculated for each day of the year using only the values for
> that day from all years of record (Figure 2). We set percentile thresholds of 2%, 5%,
> 10%, 20%, and 30% for drought identification (analogous to the percentiles currently
> used by the operational U.S. Drought Monitor), where the 10% flow equates to the flow
> value that is exceeded 90% of the time. Thresholds were calculated using all data at
> the site for each of the three periods of analysis. For each site, each climate year,
> and each threshold, we calculated:
> drought duration: the total number of days below the threshold (d)
> drought deficit: the sum of flow departures below the threshold (mm)

### 2.1 Pipeline (per gage)

1. **7-day centered smoothing** of daily Q, applied to the continuous date-indexed
   series — **not** per water year, so no artificial discontinuity at Oct 1.
   - The preprocessor hands us only qualifying years, so the frame can contain temporal
     gaps between non-adjacent valid years. The window must never blend across a calendar
     discontinuity: smooth within **contiguous runs of valid daily dates**, with a
     shrinking window at run edges requiring ≥ 4 of 7 days present (below that → NaN, and
     that day is excluded from both the threshold sample and the annual aggregates).
     **As implemented this loses no days in practice**: a centered 7-day window shrunk to
     a run edge still holds 4 values, so only runs SHORTER than 4 days (or sparse NaN
     regions inside a run) produce NaN — verified on a 25-year single-run gage where all
     9,131 days are smoothed.
   - Centered (±3 days) per decision; introduces no timing lag. Effect on annual
     aggregates is negligible either way (§12.2 notes the authors-confirmation item).
2. **Fixed thresholds via the unbiased Weibull plotting position** `p_i = i/(n+1)` on the
   ascending-sorted smoothed values pooled over the whole record — i.e. **Hyndman & Fan
   Definition 6**, which in Julia is `quantile(x, p; alpha=0.0, beta=0.0)`. This is NOT
   the default (`Type 7`, linear interpolation) used elsewhere in this repo, and the
   difference is largest exactly in the low tail we care about. The repo already uses
   Weibull plotting positions in `analyze_fdc_trends` (`p = i/(n+1)`), so the convention
   is not foreign.
   - One threshold vector of length 5 per gage (one per level). Sample size
     n ≈ 365 × years ≈ 7,300–17,000.
3. **Percentile levels**: 2, 5, 10, 20, 30 % — *magnitude* (non-exceedance) quantiles.
   The paper states the convention explicitly ("the 10% flow equates to the flow value
   that is exceeded 90% of the time"), which matches this repo's
   `Q10 = quantile(Q, 0.10)` convention (`julia/src/flow_volumes.jl:127`). USDM analogy
   for the docs: 30 % ≈ D0, 20 % ≈ D1, 10 % ≈ D2, 5 % ≈ D3, 2 % ≈ D4.
4. **Per water year, per level**:
   - `drought_duration_p{n} = count(Q_smooth < Q_thr)`
   - `drought_deficit_p{n}  = Σ max(0, Q_thr − Q_smooth)` (mm/day × 1 day = mm)
   - Departures are taken from the **smoothed** series, consistently with the
     identification step (§12.2 — a wording ambiguity in the source; smoothed is the
     internally consistent reading).
5. **No pooling, no minimum event duration** — the paper defines pure day counts and
   summed departures. The 7-day smoothing is itself the de-facto pooling mechanism (its
   standard role in threshold-level drought analysis).

### 2.2 Threshold reference window — resolved by the source text

"Thresholds were calculated using all data at the site for each of the three periods of
analysis" ⇒ thresholds are recomputed **per analysis period**, which maps exactly onto our
run-window convention: each run derives its own thresholds from its own window. So the
WY 1993–2025 and WY 1980–2025 standard products will legitimately carry **different**
thresholds and therefore different duration/deficit values.

⇒ Both metrics join the documented **record-dependent** signature list (`*_all` pulses,
elasticity, parameterized BFI): valid within a product, never mixed across products.

### 2.3 Threshold sample = the smoothed series over qualifying years
The threshold pool is the smoothed daily values of the gage's **qualifying** years only
(the preprocessor already excludes rejected years — automatic, no extra filtering), minus
the run-edge NaN days of §2.1. Require ≥ `min_years_for_threshold` (10) valid years, else
all drought metrics are NaN for that gage.

### 2.4 `<` vs `<=`
"Below the threshold" ⇒ strict `<`. Pinned by test. Matters for intermittent gages where
a low threshold can be exactly 0 (§3.3).

### 2.5 Why the variable method was dropped (rationale, retained for the record)

Weibull plotting positions span `1/(n+1)` … `n/(n+1)`, so a threshold at probability `p`
is only interpolable when `n ≥ 1/p − 1`:

| Level | Years needed (variable method) | WY 1993–2025 (n ≤ 33) | WY 1980–2025 (n ≤ 46) |
|---|---|---|---|
| 2 % | **49** | ✗ below `1/(n+1)` | ✗ below `1/(n+1)` |
| 5 % | 19 | marginal (20-yr floor) | ✓ |
| 10 % | 9 | ✓ | ✓ |
| 20 % | 4 | ✓ | ✓ |
| 30 % | 3 | ✓ | ✓ |

Under the paper's variable method the sample for each day-of-year is one value per year
(n = 20–46), so the 2 % threshold sits below the smallest plotting position at every gage
in both standard products, and 5 % is resolvable only just above our 20-year minimum.
`quantile` silently clamps in that regime, which would have shipped a "2 % threshold"
that really means "the lowest flow ever observed on that calendar day".

**The fixed method has no such problem** (n ≈ 7,300–17,000 ⇒ `1/(n+1) ≈ 0.0001`, three
orders of magnitude below the 2 % level). The `below_plotting_range_policy: "na"` guard
is implemented and tested anyway, as a defensive invariant that cannot silently degrade.

**Framing correction (Codex MAJOR-6, 2026-07-27)**: "unresolvable" overstated it. Standard
Hyndman-Fan type 6 *does* return a value below the first plotting position — it clamps to
the sample minimum. The real objections are (a) very large sampling uncertainty at the low
levels with n = 20–46, and (b) this project's deliberate refusal to report an endpoint
estimate (`below_plotting_range_policy: "na"`), which is a policy layered on top of type 6
rather than a property of it. Documentation should say "too uncertain under the selected
no-extrapolation policy", not "mathematically unresolvable".

---

## 3. Conventions & edge cases

### 3.1 Naming (frozen once shipped)
Ten bases:

```
drought_duration_p2   drought_duration_p5   drought_duration_p10   drought_duration_p20   drought_duration_p30
drought_deficit_p2    drought_deficit_p5    drought_deficit_p10    drought_deficit_p20    drought_deficit_p30
```

`p{n}` (not `Q{n}`) avoids collision with the existing magnitude-percentile flow columns.
No `_fixed` infix: the fixed method is the only one shipped, and a future variable method
would carry an explicit `_var_p{n}` marker — a non-breaking addition. (Override this if
you'd rather have the method in the name from day one; the CSV contract freezes on first
delivery.)

### 3.2 Water year vs the paper's climate year (documented deviation)
The paper aggregates over **climate years** (Apr 1 – Mar 31, the USGS low-flow
convention), specifically so the summer–autumn low-flow season is not split. Per decision
this implementation uses the **water year** (Oct 1 – Sep 30), consistent with CLAUDE.md
Critical Constraint #2, the guidelines doc ("All annual metrics should be calculated on
the water year"), the annual parquet's `water_year` key, the preprocessor's valid-year
gating, and manuscript §2.2.3.

**Consequence to document**: a drought peaking in September and persisting into October is
split across two annual values. The day-level indicator is identical under either
grouping — only the grouping differs, so record totals match; the annual series and hence
all trend statistics differ. This deviation belongs in `docs/SIGNATURES.md`, the
claude-skill, and the manuscript methods (§8).

### 3.3 Zero-flow / intermittent gages
On gages with many zero-flow days the low levels (2 %, 5 %) can have a threshold of
exactly **0**, so strict `Q < 0` never fires: duration ≡ 0, deficit ≡ 0 across all years,
and trend statistics on a constant series degenerate (Theil-Sen 0, Spearman NaN —
precedent: near-constant `D1_day`). Handling: leave as-is (honest — no sub-threshold flow
is possible) and expose the threshold value as a diagnostic (§3.6) so downstream users can
identify these gages rather than mistaking a structural zero for "no drought".

### 3.4 Un-normalized gages (`area_normalized = FALSE`, 28–37 Canadian gages)
No PPT involved ⇒ the Q-to-PPT gate does **not** apply; both metrics compute.
- `drought_duration_p{n}` — dimensionless (days), threshold is a percentile of the same
  series ⇒ **scale-invariant and valid**.
- `drought_deficit_p{n}` — **unit-carrying**: mm for normalized gages, m³/s·day for the
  un-normalized ones ⇒ add to the documented "filter on `area_normalized == TRUE` before
  cross-gage comparison of unit-carrying signatures" list.

### 3.5 Gates
- **Trend completeness** (60 % overall / 80 % first+last decade): applies; series is dense.
- **Stats floor** (20 non-NA annual values): applies — NOT exempt (dense, unlike
  recession/elasticity).
- **Record-anchored decade gate**: not applicable — the analogous "no drought this year"
  case emits a valid **0**, not NaN, which is exactly the signal a trend should see.
- **Changepoint (Pettitt)**: automatic via `generate_stats`.
- **Smoothing-induced NaN days** at valid-run edges (§2.1) slightly reduce the day
  denominator; quantified in the run report, not gated.

### 3.6 Proposed diagnostics (optional, 5 columns — §12.2)
`drought_threshold_p{2,5,10,20,30}` — the five per-gage threshold values (mm/day), as
non-8-stat scalars (same pattern as `runoff_ratio_high_count`). Makes the thresholds
auditable in the delivered CSV and exposes the §3.3 zero-threshold gages directly. Cheap
(the values are computed anyway) and recommended.

---

## 4. Relationship to existing columns

`calculate_pulse_metrics` already computes `n_low_pulses_all` / `dur_low_pulses_all` from
a fixed period-of-record 10th-percentile threshold on **raw** Q via `quantile` Type 7. The
drought family differs in three ways at once — 7-day smoothing, Weibull (Type 6) plotting
position, and five threshold levels — so `drought_duration_fixed_p10` measures a
**definitionally different** quantity, and `deficit` is magnitude-weighted, which nothing
in the current output is.

### MEASURED (2026-07-28, full WY 1993–2025 benchmark, 200,834 gage-years, 6,678 gages)

`docs/benchmarks/analyze_drought_redundancy.jl` compares the ANNUAL series — where
`n_low_pulses_all × dur_low_pulses_all` reconstructs that year's sub-threshold day count
(mathematically exact — `count_pulses` stores integer spell durations and returns their
arithmetic mean — subject only to floating representation),
since `dur` is the mean duration of that year's pulses. (The earlier summary-mean ratio
was invalid: a product of two means is not the mean of the product unless annual count
and duration are uncorrelated — Codex MAJOR-3.)

| stratum | n gage-years | Pearson r | Spearman ρ | med \|d\| | p90 \|d\| | max \|d\| | exact-equal |
|---|---|---|---|---|---|---|---|
| all | 200,834 | **0.979** | 0.971 | **2 d** | 9 d | 318 d | 32.5 % |
| intermittent (p10 threshold = 0) | 9,652 | 0.981 | **0.492** | 0 d | 0 d | 8 d | **99.87 %** (13 differ) |
| perennial | 191,182 | 0.979 | 0.970 | 2 d | 9 d | 318 d | 29.1 % |
| driest quartile (runoff ratio) | 43,093 | 0.974 | 0.954 | 1 d | 8 d | 318 d | 49.1 % |
| wettest quartile | 42,358 | 0.985 | 0.983 | 2 d | 8 d | 184 d | 25.1 % |

**Judged against the interannual signal, not the level** (Codex MAJOR-2 — a median
absolute difference is meaningless against a series *mean*; a trend analysis consumes
year-to-year variation). Within-gage over 6,304 gages with variable series:
**median r = 0.994**, 10th-percentile r = 0.971; RMSE of the difference divided by the
gage's own drought-duration SD: **median 0.117**, p90 0.253.

**Conclusion — `drought_duration_fixed_p10` IS largely redundant** with the existing pulse
pair, and the conclusion survives the correct framing: the two series track each other
year-to-year at r ≈ 0.99 within a gage, disagreeing by ~12 % of the interannual SD. The
7-day smoothing and Weibull-vs-type-7 choice shift the count slightly; they do not measure
a different phenomenon.

Two honest caveats:
- Redundancy is an **aggregate** statement. Only 32.5 % of gage-years agree exactly and the
  maximum disagreement is 318 days — there is a real tail where the smoothed Weibull
  threshold and the raw type-7 threshold diverge sharply (e.g. a mass of flow values sitting
  exactly at the raw threshold, which the strict `<` excludes but smoothing shifts below).
- On intermittent gages the earlier claim of "exactly zero disagreement" was **false**
  (Codex MAJOR-1): 99.87 % of those gage-years agree exactly, 13 of 9,652 differ, max 8 days.
  The Spearman ρ of 0.492 is precisely the signature of those 13 discordant pairs among
  otherwise-tied values — it cannot coexist with an identical-series claim, and reporting
  only median/p90 (both 0) hid that.

What the family **does** add, and the honest basis for keeping it:
1. **`drought_deficit_*` has no counterpart anywhere in the output** — it is the only
   magnitude-weighted low-flow measure, so it separates a long shallow drought from a
   short severe one, which no duration or count column can.
2. **Four of the five levels are MEASURABLY distinct** — not merely "different thresholds,
   therefore different" (Codex MINOR-3). Each level's annual duration vs the same p10
   pulse pair:

   | level | Pearson r | Spearman ρ | med \|d\| | mean(drought − pulse) |
   |---|---|---|---|---|
   | p2 | 0.712 | 0.680 | 18 d | −27.4 d |
   | p5 | 0.902 | 0.878 | 12 d | −17.0 d |
   | **p10** | **0.979** | **0.971** | **2 d** | +0.7 d |
   | p20 | 0.846 | 0.847 | 32 d | +36.7 d |
   | p30 | 0.731 | 0.740 | 70 d | +73.1 d |

   Only p10 collapses onto the existing signal; the other four carry systematically
   different severity information (the USDM D4/D3/D1/D0 rungs), which is the point of the
   source paper's severity ladder.
3. Alignment with a published, citable method (Adelsperger et al.) for the HISSS paper.

### DECISION (user, 2026-07-28): KEEP `drought_duration_fixed_p10`

The measured redundancy is **not** grounds for removal. Rationale: err toward abundance
rather than pruning a metric that is potentially redundant *across signature families* —
the same quantity computed by two documented methods (raw type-7 pulses vs 7-day-smoothed
Weibull drought) has independent value, and the drought family should be internally
complete at all five USDM rungs rather than carry a hole at D2 because a neighbouring
family happens to overlap there.

Downstream guidance (documentation, not exclusion): prefer `drought_deficit_*` and the
non-p10 levels for new analysis, and do not present `drought_duration_fixed_p10` and
`n_low_pulses_all × dur_low_pulses_all` as independent evidence for the same conclusion —
within a gage they track at median r = 0.994.

---

## 5. Code changes (Julia canonical)

| File | Change |
|---|---|
| `julia/src/drought.jl` | **NEW**. `smooth_flow_centered(dates, Q; window, min_valid)` (contiguous-run aware); `weibull_quantile(x, p)` (Definition 6 + the `p < 1/(n+1)` guard of §2.5); `drought_thresholds(Q_smooth; levels)`; `calculate_drought_metrics(df; trend_completeness, decade_completeness, min_values_for_stats, changepoint, collector)` → standard `Dict{String,Float64}`. Mirrors `snow.jl`: metric-name constants, `validate_columns`, `coalesce_q`, and a 0-row / missing-column path emitting the IDENTICAL key set via explicit `value_cols`. |
| `julia/src/StreamflowSignatures.jl` | `include("drought.jl")` (after `pulses.jl`) + `export calculate_drought_metrics` |
| `julia/src/signatures.jl` | Call in the non-climate block after `calculate_pulse_metrics`, inside the existing per-signature `try/catch`, passing `tc`, `dc`, `mv`, `cp`, `coll` (gated, not exempt) |
| `config/signatures_config.json` | New `drought` section (§6). rpkg's bundled copy updates only when the port lands |
| `julia/src/config.jl` | `CFG_DROUGHT_*` constants with absent-section defaults (same pattern as `_snow_config`) |
| `config.R` | `EXPECTED_SIGNATURE_BASES` += the 10 bases (+ a new exception constant for the §3.6 scalars if adopted) |
| `docs/benchmarks/validate_production_run.py` | Gate 4 hardcodes `ann.signature.nunique() == 90` (line 91) → **100** (or parameterize) |

`calculate_drought_metrics` needs a **`date` column** to detect contiguous runs for
smoothing — already present on every gage frame (`gage_id, date, Q, water_year, month,
dowy`), so no plumbing change. Untouched: `stats.jl`, `changepoint.jl`, the
`AnnualCollector`, the annual-parquet writer, `build_signature_explorer.py` (discovers
bases from the CSV), and the benchmark runner.

**Runtime**: one smoothing pass + 5 vectorized comparisons over ~12–17 k days per gage —
milliseconds. Benchmark wall-clock impact is negligible; the annual-parquet write grows
(§8).

---

## 6. Configuration (`config/signatures_config.json` → `drought`)

```jsonc
"drought": {
  "smoothing_window_days": 7,
  "smoothing_alignment": "center",
  "smoothing_min_valid_days": 4,
  "threshold_method": "fixed",
  "threshold_percentiles": [2, 5, 10, 20, 30],
  "plotting_position": "weibull",
  "below_plotting_range_policy": "na",
  "min_years_for_threshold": 10,
  "comment": "Streamflow drought duration/deficit after Adelsperger et al. (in review). threshold_percentiles are MAGNITUDE (non-exceedance) percentiles — 10 = the flow exceeded 90% of the time — matching the paper and this repo's Q{n} convention; USDM analogy 30/20/10/5/2 = D0/D1/D2/D3/D4. Thresholds are FIXED (whole-record, per run window => record-dependent, matching the paper's per-analysis-period thresholds) and computed with the unbiased Weibull plotting position i/(n+1) (Hyndman-Fan definition 6), NOT the Type-7 default used elsewhere. The paper's variable day-of-year method is intentionally not implemented (user decision 2026-07-27: unresolvable at the 2% level with 20-46 years of record). Aggregation is by WATER year, deviating from the paper's climate year (Apr-Mar). Absent section => family disabled."
}
```

Absent section ⇒ family disabled (no new columns), matching the conservative
`annual_values` / `stats_floor` fallback precedent. A `threshold_method` other than
`"fixed"` is an explicit error, not a silent fallback.

---

## 7. Statistics & export integration — all automatic

`generate_stats()` supplies the 8 statistics, 8 Pettitt fields, trend gates, and stats
floor; the opt-in `AnnualCollector` writes per-year rows into
`{prefix}_signatures_annual.parquet`. Both metrics are dense — every qualifying year gets
a row (0 is a value, not a placeholder).

---

## 8. Downstream impacts

- **Columns**: 1,488 → **1,653** (+165 = 10 bases × 16 fields + 5 threshold scalars; the
  scalars ship, they are not optional). Annual signature bases 90 → 100.
- **Annual parquet**: +10 dense series ≈ +11 % rows (WY 1993–2025: 16.9 M → ~18.8 M).
- **`flagged_for_high_na`**: its denominator counts all signature fields, so +165 was
  expected to flip the flag for some gages (April 2026 Pettitt / July 2026 snow
  precedent). **MEASURED 2026-07-28: it does NOT change at all** — the drought columns
  are 79.5–100 % populated, so no gage crossed the 30 % NA threshold. Zero pre-existing
  columns move (§14 additivity gate).
- **Additivity is a REQUIREMENT, not a verified result** (Codex MAJOR-4). Code inspection
  found no shared-state mutation, and the only intended existing-column change is
  `flagged_for_high_na`, but unit tests cannot establish production additivity — the
  benchmark diff (§11 step 4) is what does. Note also that the orchestrator's
  per-signature `try/catch` (`julia/src/signatures.jl:134-140`) would convert an
  unexpected drought failure on some production gage into *missing* new fields rather
  than a failed run, so the diff must check the new columns are fully populated too.
  Until the diff is run: every pre-existing column *should* be byte-identical, the family
  is additive by construction and
  touches no shared state ⇒ an unusually strong regression gate (§11).
- **Explorer**: 1,456 → 1,616 mapped signature-statistic combinations, HTML ~+11 %
  (the WY 1993–2025 explorer was already 64.6 MB — a standing watch item; the existing
  variable-list trim option may be needed).
- **Standard products**: both delivered products predate this family; carrying it means
  re-running both (~25–28 min each) plus explorers, comparison dashboards, and validation.
- **Manuscript**: §2.2.2's metric-family list gains drought metrics; the Adelsperger
  citation, the fixed-threshold-only scope, and the water-year deviation (§3.2) need
  stating → append to the queued manuscript edits in CHANGELOG → Manuscript
  Reconciliation Log.
- **Guidelines doc**: this family originates from a paper, not the governing Google Doc.
  Flag for the user to add the definitions there so the doc remains methodology ground
  truth (repo convention: doc leads, code implements).

---

## 9. Tests (`julia/test/test_drought_metrics.jl` + `runtests.jl` include)

Mirroring `test_snow_metrics.jl` rigor — exact hand-derived expectations:

1. **Weibull quantile**: hand-computed Definition 6 values on small samples; equality with
   `quantile(x, p; alpha=0, beta=0)`; **`p < 1/(n+1)` returns NaN** under the `na` policy
   (and the sample minimum under `clamp`, to pin the guard's contrast).
2. **Centered smoothing**: exact values on a ramp; shrinking window at run edges with the
   ≥ 4-of-7 floor; **no blending across a year gap** (two non-adjacent valid years ⇒ NaN
   at the seam, never a value averaged across the discontinuity).
3. **Exact arithmetic**: synthetic gage with a known threshold and a known set of
   sub-threshold days → duration and deficit match hand-computed values exactly.
4. **Wet year**: all smoothed Q above threshold → `0.0` / `0.0` (valid zeros, NOT NaN).
5. **Monotonicity in level**: duration and deficit are non-decreasing in the threshold
   percentile (p2 ≤ p5 ≤ p10 ≤ p20 ≤ p30) at every gage-year — a strong structural
   invariant that catches level/threshold mix-ups.
6. **Threshold definition**: on a long synthetic record, `drought_duration_p10` averages
   ≈ 10 % of the year by construction (the fixed-threshold self-check of §11.6).
7. **Strict comparison** (`<`, not `<=`) and the zero-threshold/intermittent case (§3.3).
8. **Threshold data floor**: fewer than `min_years_for_threshold` valid years → all-NaN
   key set.
9. **Contract invariants**: 0-row and missing-column inputs emit the IDENTICAL key set;
   with/without `collector` the statistics are identical; zero warnings on the happy path;
   no duplicate collector keys.
10. **Gate interaction**: stats floor (< 20 values → all NaN), trend completeness applied,
    changepoint fields present when enabled.
11. **`smoke_test.jl`**: presence + ranges (`0 ≤ duration ≤ 366`, `deficit ≥ 0`,
    `deficit == 0 ⟺ duration == 0`, monotone in level).

---

## 10. Docs

- `docs/SIGNATURES.md`: new category section (definitions, the fixed-threshold method,
  the five levels + USDM analogy, Weibull plotting position vs. the Type-7 default,
  7-day centered smoothing + the contiguous-run rule, units, record-dependence, the
  water-year deviation of §3.2, the dropped variable method and why, Adelsperger + Laaha
  references); summary table row; changepoint scope 90 → 100.
- `docs/DEVELOPMENT.md`: signature/column counts; annual-parquet row growth.
- `CHANGELOG.md`: dated entry + the record-dependent-signature note in the cross-product
  consistency guidance.
- `claude-skill/streamflow-signatures.md`: interpretation guidance — rising deficit with
  flat duration = deepening but not lengthening droughts; the five levels as a severity
  ladder (USDM analogy) and their strong mutual correlation; the deficit units caveat;
  never mix across the two standard products.
- `config.R` comment block alongside the new bases.

---

## 11. Validation & rollout

1. `julia --project=julia julia/test/runtests.jl` — full suite green (1,302 assertions + new file).
2. `smoke_test.jl` on the 10-gage subset — schema + ranges + level monotonicity.
3. Full Julia benchmark on ONE window (WY 1993–2025 @ 60 % first, ~28 min) into its own
   output folder per the experiment-folder convention.
4. **Additivity gate**: diff against the existing WY 1993–2025 product — every
   pre-existing column identical; the only delta is the 160 new columns plus the
   whitelisted `flagged_for_high_na` shift.
5. `validate_annual_values.py --floor 20` (drought NOT exempt) and
   `validate_production_run.py` with the 100-signature gate.
6. Distribution sanity: duration ∈ [0, 366]; deficit ≥ 0; monotone in level;
   **`drought_duration_p10` mean ≈ 36.5 days** (10 % of a year, by construction — a strong
   self-check that the threshold is right); correlation with
   `n_low_pulses_all × dur_low_pulses_all` high but < 1 (§4); arid/intermittent gages
   inspected against §3.3; smoothing-edge NaN-day count reported.
7. Codex adversarial review (plan → implementation → results), per project convention.
8. Then decide with the user whether to regenerate both standard products.

---

## 12. Decisions

### 12.1 Resolved (user, 2026-07-27)
See the table in §0: fixed method only, all five levels, water-year grouping, NaN
below-plotting-range policy, centered smoothing.

### 12.2 Resolved 2026-07-27 (second round)
- Diagnostics **included** (five `drought_threshold_fixed_p{n}` scalars).
- Column naming: **explicit `_fixed_` infix**.
- Deficits measured from the **smoothed** series.

### 12.3 Still open
- Confirm with the paper's authors: smoothing alignment (centered assumed; USGS
  WaterWatch convention would be trailing — the `smoothing_alignment` config knob
  implements both) and that deficits are indeed differenced against the smoothed series.
- Sequencing: land the family and batch the standard-product re-runs with the queued
  Python/rpkg ports, or re-run both products immediately after validation?

---

## 13. Out of scope

- **The paper's variable (day-of-year) threshold method** — dropped by decision (§2.5).
  The config and module are shaped to admit it later under `_var_p{n}` column names if a
  longer record or a ±window-pooled variant is ever wanted.
- **Climate-year (Apr–Mar) aggregation** — the paper's grouping; deviation documented
  (§3.2).
- Event-level drought metrics (`n_drought_events`, max duration, event
  severity/intensity, onset/termination timing, pooling by inter-event criteria) — the
  natural extension once the aggregate metrics are validated, and likely where the
  paper's "severity-based" novelty lives; not in the current request.
- Seasonal drought variants (winter/spring/summer/fall duration + deficit).
- Precipitation- or SWE-based drought (SPI/SSI-style) — this family is streamflow-only.
- Python/rpkg ports (deferred to the existing port queue).

---

## 14. Implementation record (2026-07-27)

**Delivered** (Julia):

| File | What |
|---|---|
| `julia/src/drought.jl` | NEW — `weibull_quantile`, `smooth_daily_flow`, `calculate_drought_metrics` |
| `julia/src/config.jl` | `CFG_DROUGHT_*` + fail-fast validation of unimplemented options (method / plotting position / alignment / below-range policy) |
| `julia/src/StreamflowSignatures.jl` | include + exports (`calculate_drought_metrics`, `weibull_quantile`, `smooth_daily_flow`, `CFG_DROUGHT_ENABLED`, `CFG_DROUGHT_PERCENTILES`) |
| `julia/src/signatures.jl` | orchestrator call after `calculate_pulse_metrics`, gated on `CFG_DROUGHT_ENABLED` |
| `config/signatures_config.json` | `drought` section |
| `config.R` | 10 bases + `EXPECTED_DROUGHT_THRESHOLDS`, wired into `validate_output_schema` |
| `julia/test/test_drought_metrics.jl` | NEW (+ `runtests.jl` include) |
| `julia/test/smoke_test.jl` | drought presence / range / level-monotonicity checks |
| `docs/benchmarks/validate_production_run.py` | signature-count gate 90 → 100 |
| Docs | `SIGNATURES.md` §13 + summary table + changepoint scope; `DEVELOPMENT.md`; `CHANGELOG.md`; claude-skill |

Also updated: `julia/test/test_annual_collector.jl` — `EXPECTED_DENSE_SIGNATURES` gained
the 10 drought bases. **That test asserts set EQUALITY of the collected signatures, so
ANY new dense family fails it until registered** — the one real failure the suite caught,
now added to the CLAUDE.md / DEVELOPMENT.md "Adding New Signatures" checklists.

**Verification status.** Julia 1.12.6 was installed locally (juliaup, `~/.juliaup/bin`)
and `julia/test/runtests.jl` runs **green: 2,042 assertions, 0 failures** — 670 of them
in `test_drought_metrics.jl`. `smoke_test.jl` is written but not exercised here (it
hardcodes the Windows `D:\processedOuts_feb2026\` inputs).

Expectations were first derived against a purpose-built pure-Python mirror of
`drought.jl`, then confirmed by the Julia run — both agree on —

- Weibull quantile hand values on `n=9` (p10→1.0, p25→2.5, p30→3.0, p05→NaN, clamp→1.0)
- centered smoothing of `Q = 1…10` → `[2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 8.5]`;
  trailing → first three NaN; three-run gap case with no cross-gap leakage; duplicate-date
  run break → `[NaN, 2.5, 3.0, 3.0, 3.0, 3.5]`; sparse-NaN windows
- 25-year distinct-valued gage: pool `n = 9,131`, all days smoothed, and for all five
  levels `Σ duration == count(pool < thr) == floor((n+1)p)` exactly (182 / 456 / 913 /
  1,826 / 2,739) with `Σ deficit` matching the pooled departure sum; monotone in level;
  `drought_duration_fixed_p10_mean = 36.52`
- wet/dry-contrast gage: wet years report valid `0.0` for duration AND deficit at all levels
- intermittent gage (300 zero-flow days/yr): every threshold `0.0`, every duration `0.0`
  (confirms the strict `<`)
- 9-year record: NaN thresholds and NaN metrics (threshold floor)

Mirror script (disposable): the session scratchpad `drought_ref.py`.

**Real-data validation (2026-07-27, 10 smoke gages × 45–46 WYs, feb2026 inputs from the
thumbdrive at `/Volumes/Untitled/`).** `smoke_test.jl` PASSED (all 10 drought bases
present, all ranges valid, level monotonicity holds). Mean duration tracks the
construction target `p × 365.25` to within 0.04 days — a **weak, near-circular
consistency check**, since the threshold is a percentile of the same series it is then
counted against (Codex MAJOR-2): it would also hold if both steps used the wrong series,
and plotting-position choice moves it by only a few observations. It rules out gross
mismatch, nothing finer:

| level | mean duration across gages | target | per-gage spread |
|---|---|---|---|
| p2 | 7.296 d/yr | 7.30 | 7.24 – 7.30 |
| p10 | 36.499 d/yr | 36.52 | 36.46 – 36.52 |
| p30 | 109.542 d/yr | 109.58 | 109.46 – 109.57 |

Fixed p10 thresholds ranged 0.128–0.610 mm/day, p10 deficits 1.88–4.39 mm/yr — plausible
for humid New England gages. A ratio of `duration_fixed_p10` to
`n_low_pulses_all × dur_low_pulses_all` came out at mean 0.992, range [0.811, 1.145],
but **that comparison is not valid evidence of non-redundancy** and the earlier claim
that it was has been withdrawn (Codex MAJOR-3): the denominator multiplies two summary
means, which equals mean annual sub-threshold days only when annual count and duration
are uncorrelated. See the revised §4.

**Adjacent pre-existing bug surfaced by the smoke run (NOT introduced here, not fixed):**
`calculate_negative_days` (`julia/src/pulses.jl:323`) applies `!isnan(x) && x < 0`
directly to a `Union{Missing,Float64}` Q column, so a `missing` value throws
`TypeError: non-boolean (Missing) used in boolean context`. The orchestrator's per-signature
`try/catch` swallows it and the gage silently loses all 8 `negative_ann` columns (gage
01073000 in the smoke set: 808 signature columns vs 816 for the other nine). Every other
signature routes Q through `coalesce_q`. Production is unaffected while
`use_legacy_filtering: false` (the preprocessor emits Float64 with NaN), so this only
bites the raw/legacy path. One-line fix available; logged in CHANGELOG → Known Issues.

**Benchmark run (2026-07-28, WY 1993–2025 @ 60 %, thumbdrive inputs, M1).** Wrapper
`docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl` → own folder
`/Volumes/Untitled/processedOuts_drought_28jul2026/`.

- **6.05 min**, 18.4 gages/s (the July x86 run of the same window took 27.6 min).
- **6,678 gages × 1,653 columns** — the gage set is IDENTICAL to standard product #1, so
  the family changes no gage qualification, and 1,653 confirms the corrected +165
  arithmetic (MAJOR-5).
- Annual parquet **18,898,406 rows / exactly 100 signatures**; drought contributes
  10 × 200,834 = 2,008,340 rows — perfectly dense (one row per gage-year per level),
  **zero NaN, zero duplicate keys**.
- **Two attempts emitted 1,642 columns**, silently missing the 11 GAGES-II/HYDAT
  human-interference columns: `metadata.gages_ii_dir` pointed at the Windows
  `D:/gagesMetadata` and the loader only warns. The first fix — an ENV-derived
  **constant** — did not work, because a const is baked at precompile time and the
  wrapper sets the ENV var at runtime (the §14 caching trap, third occurrence in one day).
  The working fix is a **runtime accessor** `gages_ii_dir()` used as
  `load_gages_ii_interference`'s default argument. A useful reminder that "column count
  differs" can be an input-path problem, not a signature problem — and that any ENV
  override consumed through a `const` is unreliable in this codebase.

- **Provenance re-run (2026-07-28)**: re-ran the same configuration after adding the
  provenance block. Output is **bitwise identical** to the audited run (both the
  signatures CSV and the annual parquet), which simultaneously (a) proves the provenance
  addition is inert with respect to results, (b) demonstrates same-machine reproducibility
  of the whole pipeline, and (c) means the Codex-audited values ARE the values in the
  final folder. The additivity gate was re-run against the drought-disabled baseline and
  passes again. Final artifacts + `RUN_NOTES.md`:
  `/Volumes/Untitled/processedOuts_drought_28jul2026/`.

### ✅ ADDITIVITY GATE: PASS (2026-07-28) — Codex MAJOR-4 closed

The rigorous test is a **same-machine, same-Julia baseline with the drought family
disabled** (a config copy with the `drought` section removed, via `STREAMFLOW_CONFIG`),
diffed against the drought-enabled run with EXACT equality:

```
new: 6678 rows x 1653 cols | reference: 6678 rows x 1488 cols
[PASS] 1 no reference column dropped;  added 165 (expected 165)
[PASS] 2 gage sets identical (6678 gages)
[PASS] 3 all 1487 shared columns unchanged (bitwise)
[PASS] 4 every added column populated (finite fraction min 0.795 / median 0.955 / max 1.000)
ADDITIVITY GATE: PASS
```

**Stronger than predicted: `flagged_for_high_na` did NOT shift.** Every earlier note
(here, in CHANGELOG, and in the snow-family precedent) expected it to move because its
denominator counts all signature fields. Direct per-gage recomputation found **zero
crossings** — 1,224 gages flagged before and after (independently confirmed by Codex).
The whitelist was not needed; zero pre-existing columns changed by any amount. **Do not
generalize**: the mechanism is per-gage (each gage's prior NA fraction vs its drought NA
pattern), not a general property of adding populated columns, and the closest gage sits
only ~6.1e-5 from the 30 % threshold — a different window could well flip some gages.

**Why ~20 % of some drought columns are NaN — two distinct mechanisms, not one**
(Codex MINOR-1; earlier text conflated them):
- **495 gages (7.4 %)** lose the six trend statistics to the standard trend-completeness
  gate (60 % overall / 80 % first-and-last decade). This is NOT drought-specific: the set
  is *identical* to the gages missing `Qann_senn_slp`, and their median completeness is
  0.839 vs 1.000 for those that pass.
- **A further ~300 gages** lose the rank statistics (Spearman/MK) because their annual
  series is constant — a rank correlation is undefined with zero variance. Theil-Sen
  still returns 0.
- The 0.795 column minimum is `drought_duration_fixed_p2_pettitt_post_mk_pval`
  specifically, and **93.6 % of its 1,368 NaNs are gages whose p2 series is constant
  all-zero** — i.e. never experienced an exceptional (D4) drought in the window. Expected
  and concentrated at the most extreme severity level.
- `_mean` / `_median` / the five threshold scalars are **100 %** populated.

**Cross-machine diff vs the delivered product #1** (kept as a secondary result): no column
dropped, 165 added, gage sets identical. 431 shared columns differ only as floating-point
noise (≤ 3.2e-06 relative). A further 66 differ *materially* on 1–38 rows out of 6,678 —
**all of them statistics that are discretely sensitive to a last-bit change**, and none
attributable to the drought family:
- `TQmean*` (a day-count threshold: one flipped day = 100/365 = 0.274 pp; the observed
  0.137 median step is exactly a half-day, i.e. a median between two years),
- `FDC90th*` (OLS on `log10(Q + 1e-10)` in the near-zero low-flow tail — already
  documented in this repo as FP-fragile: "28-gage NA mismatch from floating-point
  precision in near-zero regression"),
- 38 rank/p-value fields (`*_mk_pval`, `*_spearman_*`, `*_pettitt_*`) where a tie flip
  changes a rank statistic discretely.

The defensible statement (Codex MINOR-5 — the earlier wording overclaimed): these 66
differences are **consistent with platform/version floating-point effects and are
demonstrably unrelated to enabling drought** (the same-machine gate proves the latter).
Proving the former would need the July build reproduced with recorded input hashes and
software versions, which we do not have.

Lesson recorded: **a cross-machine diff cannot serve as an additivity gate** for this
codebase; only a same-machine, drought-disabled baseline can. `--tol` was added to the
tool so the distinction is visible rather than drowned in 431 noise rows.

**Independent verification by Codex (2026-07-28)**, beyond what was done here: it
confirmed `config_nodrought.json` is structurally identical to the committed config after
removing only `drought` (no key or numeric value differs), and compared the two runs'
annual parquets directly — **all 16,890,066 shared `(gage_id, water_year, signature)`
keys match with zero value mismatches**, so shared summary-column equality is not
concealing changed annual series.

**Outstanding step**
1. Standard-product decision: regenerate both delivered products to carry the family.

---

## 15. Review targets (claims to falsify)

Ranked by consequence. Each is a specific, checkable assertion — not "review the diff".

**Correctness of the metric definition**
1. `weibull_quantile` is genuinely Hyndman-Fan definition 6. Attack: compare against
   `quantile(x, p; alpha=0, beta=0)` on random samples INCLUDING ties and n where
   `(n+1)p` is exactly integral; check the `h >= n` upper clamp and the `h < 1` guard
   boundaries at `p == 1/(n+1)` exactly.
2. Thresholds are percentiles of the **smoothed** pool, not raw Q — and the deficit is
   differenced against the same smoothed series. A mismatch here would silently break
   the ≈`p`-fraction-of-days property while still producing plausible numbers.
3. The exact-count invariant used in the tests (`Σ duration == floor((n+1)p)`) holds only
   for distinct-valued samples. Attack: is the test gage actually distinct-valued, and
   does the invariant silently weaken (rather than fail) on a tied sample — i.e. is the
   test proving what it claims?

**Attribution and windowing**
4. Days are partitioned across water years with no double-count or loss. The tests assert
   this via the summed-duration identity; attack the case where the frame's first/last
   water year is partial, and where `water_year` is `Union{Missing,Int}`.
5. Smoothing never crosses a temporal discontinuity, and a valid year adjacent to a
   REJECTED year gets a shrunken-window value rather than a gap-contaminated one.
   Attack: is shrinking (vs NaN-ing) the right call at a rejected-year boundary — it
   mixes a 4-day mean with 7-day means inside the same threshold pool.
6. `min_years_for_threshold = 10` interacts with the 20-year gage floor and the 20-value
   stats floor. Is a 10–19-year gage able to emit thresholds and means but no trends, and
   is that coherent?

**Integration and contracts**
7. Purely additive: no pre-existing column can change except `flagged_for_high_na`.
   Attack via the benchmark diff, not by reading code.
8. The all-NaN early-return paths emit a key set byte-identical to the happy path
   (including Pettitt fields), so the CSV schema is stable across gages.
9. Config fail-fast: an unsupported `threshold_method` / `plotting_position` /
   `alignment` / below-range policy must error at module load, never silently fall back.
10. Registration completeness — `config.R` bases + scalar constant,
    `EXPECTED_DENSE_SIGNATURES`, the `validate_production_run.py` count gate. One of
    these was missed on the first pass and only the unit suite caught it; are there
    others (e.g. anything keyed off signature counts in the explorer or dashboards)?

**Interpretation claims made in the docs**
11. "`drought_duration_fixed_p10_mean` ≈ 36.5 days/yr" — verified on synthetic data and
    on the 10-gage smoke set; does it hold across the full 6,000+ gage benchmark, and is
    the claim safe for arid/intermittent gages where the threshold is 0?
12. The redundancy claim in §4 (drought vs `n_low_pulses_all × dur_low_pulses_all`):
    correlated but not equal. Quantify on real output rather than asserting it.

---

## 16. Codex adversarial review record (2026-07-27, codex-cli 0.145.0, read-only)

**Verdict: GO-WITH-FIXES.** No CRITICAL defect in the fixed-threshold arithmetic. Six
MAJOR and five MINOR findings; all resolved below. Codex independently re-ran the unit
suite (2,042 assertions at the time) and confirmed the exact-count invariant on the test
fixture. It could not re-run the mounted-data smoke test (its `Pkg.activate` needs to
write a manifest-usage lock, blocked by the read-only sandbox).

| # | Finding | Resolution |
|---|---|---|
| MAJOR-1 | The "record invariant" proves **conservation, not attribution** — shifting every `water_year` by one keeps the totals identical | **Fixed.** New testset "water-year attribution": (a) hand-derived boundary fixture — a 7-day low-flow block straddling Sep 30 → Oct 1 splits as exactly 6 days / 27.0 mm in WY 2000 and 7 days / 36.0 mm in WY 2001 at every level; (b) every per-year duration recomputed from the smoothed series against an INDEPENDENT `year + (month ≥ 10)` mapping. Docs reworded: the old invariant proves conservation |
| MAJOR-2 | The `p × 365.25` threshold check is **substantially circular** — a percentile of a series counted against that same series yields ≈p regardless of whether the series is the right one | **Fixed (docs).** Relabelled a "weak internal consistency check" in SIGNATURES.md, the skill, CHANGELOG and §14; "Threshold correctness confirmed" and "large departures indicate a threshold problem, not a signal" withdrawn; noted it can legitimately fail on tied/intermittent records |
| MAJOR-3 | The redundancy ratio is **not evidence of non-redundancy**: `mean(n) × mean(d) ≠ mean(n × d)` unless count and duration are uncorrelated, and 10 humid gages are unrepresentative | **MEASURED 2026-07-28, and the claim REVERSED.** Proper annual-series comparison over 200,834 gage-years (`analyze_drought_redundancy.jl`): r = 0.979, median disagreement 2 d of ~36, uniform across aridity quartiles, exactly 0 on intermittent gages ⇒ `drought_duration_fixed_p10` IS largely redundant with the pulse pair. §4 rewritten; the family's real contribution is `deficit` (no counterpart) + the four non-p10 levels. Docs and skill now steer users accordingly |
| MAJOR-4 | Additivity described as established while the gate is still outstanding; the orchestrator's `try/catch` could hide a drought failure as missing columns | **CLOSED 2026-07-28 by measurement.** Same-machine drought-disabled baseline vs drought-enabled run: 165 added, 1,487 shared columns **bitwise identical**, every added column populated. `flagged_for_high_na` did not even shift. New reusable tool `check_additivity.jl`; §14 records why a cross-machine diff cannot substitute |
| MAJOR-5 | Column arithmetic wrong — 160 + 5 shipped scalars = **+165 → 1,653**, not 1,648 | **Fixed** in the plan, CHANGELOG (annual base count of 100 was already right — scalars are not annual series) |
| MAJOR-6 | "Unresolvable" overstates it: type 6 clamps to the sample minimum below `1/(n+1)`; the NaN rule is a project policy | **Fixed (docs).** Reframed as large sampling uncertainty + a deliberate no-extrapolation policy, in SIGNATURES.md, the skill, CHANGELOG, config comment and §2.5 |
| MINOR-1 | Rows with missing/unparseable dates could still be smoothed with their neighbours | **Fixed (code).** `smooth_daily_flow` now ends the run *before* an invalid row; the row itself stays NaN. Regression test added |
| MINOR-2 | `min_years_for_threshold` counts distinct `water_year` labels, not complete years | **Documented** as an explicit caller contract in the `calculate_drought_metrics` docstring |
| MINOR-3 | The smoke test could PASS with every drought value NaN (WARN not FAIL; monotonicity vacuously true) | **Fixed.** All-NaN and partially-finite columns now FAIL on this 45+-year gage set, and an incomplete monotonicity sequence fails instead of being skipped |
| MINOR-4 | Config validation caught bad strings but not bad numbers (even "centered" window, `min_valid > window`, out-of-range/duplicate/unsorted percentiles, non-positive minimum years) | **Fixed (code).** All validated at module load in `config.jl` |
| MINOR-5 | The "Adding New Signatures" checklist still omitted total-column-count, additivity, and finite-value checks | **Fixed.** Added to both CLAUDE.md and DEVELOPMENT.md |

Post-fix suite: **2,700 assertions, 0 failures** (drought file 670 → 1,328).

Accepted as-is (Codex agreed these are defensible, disclosed decisions, not defects): the
water-year deviation from the paper's climate year, the fixed-thresholds-only scope, and
Weibull/type 6 for the fixed thresholds. Its caveat is adopted in the docs: these outputs
should not be described as direct reproductions of the paper's annual metrics.

---

## 17. Codex review #2 — RESULTS and ANALYSES (2026-07-28, read-only)

**Verdict: GO-WITH-FIXES, and the product is approved for promotion** to standard product
#1 once the documentation errors below are corrected. No CRITICAL. 4 MAJOR + 7 MINOR; all
resolved. Codex independently re-verified the config equivalence, the gage sets, and
**all 16,890,066 shared annual values** (zero mismatches), and confirmed the earlier
round's fixes are substantive rather than cosmetic (it specifically checked that the new
water-year attribution test would fail if the mapping were corrupted).

| # | Finding | Resolution |
|---|---|---|
| MAJOR-1 | **"Exactly zero disagreement on intermittent gage-years" is FALSE** — r = 0.981 with ρ = 0.492 cannot coexist with identical series; the script printed only median/p90, which rounding to 0 proves ≥90 % agree, not all | **Fixed.** Measured: **99.87 % exact-equal, 13 of 9,652 differ, max 8 days**. `report()` now prints exact-equal %, nonzero count and MAX \|diff\| alongside the quantiles; the false claim is corrected in §4, CHANGELOG, SIGNATURES.md and the skill |
| MAJOR-2 | The redundancy argument scaled a median absolute difference against the series **mean**; a trend consumes interannual variation | **Fixed.** Added within-gage analysis: median r = **0.994** (p10 0.971), RMSE/SD median **0.117** (p90 0.253) over 6,304 variable-series gages. The redundancy conclusion **survives** the correct framing; docs now cite these instead |
| MAJOR-3 | `check_additivity.jl`'s population gate passed a column with **one** finite value of 6,678 — it could not detect the very failure mode it existed for | **Fixed.** New `--min-finite-frac` (default 0.5) plus a by-suffix finite-fraction breakdown, so gating (trend/floor) is distinguishable from breakage. Re-run at 0.75: **PASS** |
| MAJOR-4 | Run construction strongly corroborated but **not provenance-complete** — timing JSON lacked input hashes, config hash, git revision, Julia version | **Fixed (tooling).** `run_julia_benchmark.jl` now records a `provenance` block: resolved paths + size/mtime for every input, SHA-256 for files < 50 MB (config, metadata; `STREAMFLOW_HASH_INPUTS=1` to hash the multi-GB parquets), git revision + dirty flag, Julia version, hostname, and the experiment ENV overrides. **The existing product predates this** — capturing it needs a 6-min re-run |
| MINOR-1 | The "~20 % lack trends" summary conflated two mechanisms | **Fixed.** Separated: 495 gages (7.4 %) lose 6 trend fields to the completeness gate — an **identical set** to those missing `Qann_senn_slp`, so not drought-specific; ~300 more lose rank stats to constant series; the 0.795 floor is `drought_duration_fixed_p2_pettitt_post_mk_pval`, 93.6 % of whose NaNs are constant all-zero p2 series |
| MINOR-2 | The `flagged_for_high_na` causal explanation was too loose | **Fixed.** Now reported as a direct per-gage recomputation (0 crossings, 1,224 before and after) with an explicit "do not generalize — closest gage is 6.1e-5 from the threshold" |
| MINOR-3 | The four non-p10 levels were argued distinct, not measured | **Fixed.** Measured table added to §4 (r = 0.712 / 0.902 / 0.846 / 0.731 for p2/p5/p20/p30) |
| MINOR-4 | "exactly" overstated the `n × dur` identity | **Fixed.** "Mathematically exact, subject to floating representation" |
| MINOR-5 | Calling the 66 cross-machine differences "platform FP noise" overclaims | **Fixed.** Now: "consistent with platform/version FP effects and demonstrably unrelated to enabling drought" — the same-machine gate proves only the latter |
| MINOR-6 | The config-cache workaround is correct but fragile; past experiments may be affected | **Fixed (documented).** DEVELOPMENT.md carries an audit: no recorded past result is invalidated (the July 60 % gate has empirical evidence it took effect; the snow gates involved source changes), but any future/historical run relying **only** on `STREAMFLOW_CONFIG` is untrustworthy. Durable fix (config as runtime data) logged as a follow-up |
| MINOR-7 | — | Confirmed the round-1 fixes are real, not cosmetic |

**Follow-ups deliberately NOT done here** (logged, not silently dropped): making config
values runtime data rather than precompiled constants; hashing the multi-GB inputs by
default; expected-eligibility masks per statistic in the additivity checker.

---

## 18. Promotion record (2026-08-10 / 08-11)

The drought family shipped to BOTH standard products; the July 2026 folders are superseded.

| | Standard product #1 | Standard product #2 |
|---|---|---|
| Window | WY 1993–2025 @ 60 % | WY 1980–2025 @ 60 % |
| Folder | `processedOuts_drought_28jul2026` | `processedOuts_1980_2025_11aug2026` |
| Promoted | 2026-08-10 | 2026-08-11 |
| Gages × cols | 6,678 × 1,653 | 6,250 × 1,653 |
| Annual parquet | 18,898,406 rows / 100 sigs | 24,366,487 rows / 100 sigs |
| Climate input | ORIGINAL daymet parquet | REBUILT from the annual CSVs |
| Additivity gate | PASS (28 Jul, same-machine, bitwise) | PASS (11 Aug, same-machine, bitwise) |
| Annual-values validation | PASS (612,046 pairs) | PASS (587,829 pairs) |
| Comparison vs superseded | 1,446/1,451 Perfect, min R² 0.9977 | 1,447/1,451 Perfect, min R² 0.9977 |

Both folders carry the full standard-product artifact set (explorer + `_annual/` sidecars,
comparison CSV/summary/dashboard, `validation_gates.txt`, additivity + control reports,
`RUN_NOTES.md`).

### Codex review of the promotion (2 rounds, read-only, codex-cli 0.145.0)

Round 1: **GO-WITH-FIXES on #1, NO-GO on #2** — 1 CRITICAL, 6 MAJOR, 3 MINOR. The CRITICAL
rejected "different machine" as an explanation for #2's 74 material shared-column
differences vs the July build, correctly noting several were **Q-only** metrics no climate
change could touch. Answered with three controls (input-only ≤ 3.4e-13; drought-only
bitwise; machine-only reproducing the pattern with the ORIGINAL climate input on both
sides) plus an annual-level comparison showing the underlying series agree to 5.68e-14.

Round 2 (delta): **GO-WITH-FIXES on both** — "the former NO-GO blocker is cleared".
It independently reproduced both annual-control numbers to full precision and confirmed the
drought-disabled twin proves itself by emitting 1,488 columns / 90 annual signatures. It
found one further overstatement, now corrected: the cross-machine divergence is **not**
"only rank-based statistics" — 14 of 67 material columns are `TQmean` / Pettitt / `FDC90th`
fields. The accurate framing is **discretely FP-sensitive** statistics (rank ties flip;
`TQmean` is a day count so one day = 0.274 pp; the Pettitt changepoint LOCATION jumps;
`FDC90th` is OLS on `log10(Q + 1e-10)` in the near-zero tail).

### Infrastructure discovered / built during promotion

- `docs/benchmarks/convert_daymet_csvs_to_parquet.py` — rebuilds the climate parquet from
  the 44 annual CSVs after the canonical one was found truncated. **Daymet uses a 365-day
  calendar** (leap years drop Dec 31, not Feb 29).
- `build_experiment_vs_julia_dashboard.py`'s `SIGNATURE_GROUPS` had never included the SNOW
  family, so every dashboard built before 2026-08-10 silently omitted all 14 snow bases.
  Fixed (selectable targets 623 → 820).
