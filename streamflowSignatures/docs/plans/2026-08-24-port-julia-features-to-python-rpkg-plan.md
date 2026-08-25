# Port the six Julia-only features to Python and rpkg — campaign plan

**Date**: 2026-08-24 · **Status**: PLANNED — revised after Codex adversarial review
(initial verdict NO-GO; all findings incorporated, see §8) · **Owner**: Arik + Claude

## 1. Objective & drivers

Bring `python/streamflow_signatures/` and `rpkg/` to full parity with the canonical
Julia implementation: the **1,653-column summary CSV** (100 time-series signature
bases + scalars, incl. Pettitt fields) **and the annual-values parquet**, validated
against a same-machine Julia reference.

Drivers:
- **CHANGELOG → Planned** six-feature port queue (standing since July 2026).
- **HISSS manuscript (synced 2026-08-24)** now claims: *"Because every signature is
  computed by an open, cross-language library (Julia, Python, and R), each metric is
  reproducible and extensible"* and *"the code for the full workflow is available in
  an open library in Julia, Python, and R."* Currently false — Python/rpkg reproduce
  only the April 2026 623-signature-column subset. Submission target Nov 9 2026.

## 2. Decisions (user, 2026-08-24)

1. **Scope = all six features** (incl. the annual-values collector + parquet export).
2. **Install R locally on this Mac** — both ports developed AND benchmarked here,
   same M1 machine as the Julia reference products (avoids the documented
   cross-machine FP-divergence problem in comparisons).
3. **Fix Julia's seasonal runoff-ratio flag-name bug FIRST** (Julia-first workflow),
   then produce a fresh Julia reference run for port comparison. Delivered standard
   products are untouched; the fix lands in future product regenerations.

## 3. Gap summary (code survey 2026-08-24; every claim below re-verified by Codex)

The six features have **no functional implementation or orchestration path** in
either port (Codex F20: "100% absent / zero grep hits" was rhetorically overstated —
config sections are visible to loaders and the words appear in docs — but the
functional gap is confirmed: neither orchestrator imports or calls any of it, and
neither `generate_stats` has the new kwargs).

| Feature | Julia source | Notes for the port |
|---|---|---|
| Pettitt changepoint | `changepoint.jl` (365 ln), `stats.jl:520-556`, `CP_SUFFIXES` `stats.jl:16-19` | +8 `_pettitt_*` suffixes per time-series signature; `changepoint` kwarg threaded to 14 orchestrator call sites. The BIC `detect_changepoint` path is **exported Julia-only API not used by the production pipeline — intentionally excluded from parity scope** (Codex F11). ⚠ `test_changepoint.jl` is NOT wired into `julia/test/runtests.jl` — fix in Phase 0 |
| 20-value stats floor | `stats.jl:334, 390-405`, `config.jl:311-315` | Orchestrator passes `min_values_for_stats` to 11 call sites, **deliberately NOT** recession/elasticity (`signatures.jl:66`). **This changes existing port columns** (below-floor fields → NaN) — its validation gate must be a floor-transition check, not additivity (Codex F1) |
| Annual collector | `stats.jl:235-283` (struct + `_collect_annual!`), collect **before** all gates (`stats.jl:352-353`) | `collector` kwarg through all 14 signature functions; runner drains to long-format parquet |
| b=1 recession alpha | `recession.jl:471, 490, 193-204` (`event_log_a_b1`) | Ports still compute `log_a_pointcloud`/`log_a_events` with **free-fit b** (semantic divergence behind existing columns); `alpha_linear` + `recession_alpha_point_cloud_linear_reservoir` already correct in all three. Seasonality sinusoid input must switch to per-event b=1 log_a values too |
| Snow (14 metrics) | `snow.jl` (324 ln) + preprocessor SWE plumbing (`io.jl:42-43, 192-194, 449-515, 554`) + record-anchored decade gate (`snow.jl:284-315` via `force_skip_trends`, `stats.jl:335, 417`) | Port preprocessors emit **no** `valid_swe_years` and have zero SWE awareness; orchestrator needs `snow_data`/`snow_climate_years` kwargs; both port runners merge only PPT — SWE never reaches the ports today |
| Drought (10 + 5 scalars) | `drought.jl` (315 ln): `weibull_quantile` (H-F type 6 + `below_plotting_range_policy: "na"`), `smooth_daily_flow` (7-day centered within contiguous date runs), `CFG_DROUGHT_ENABLED` gate | Do NOT reuse `fdc`'s Weibull exceedance code path — different convention |

Supporting gaps (survey + Codex-verified):
- **`generate_stats`**: Julia has 4 kwargs the ports lack (`changepoint`, `collector`,
  `min_values_for_stats`, `force_skip_trends`). `trend_completeness`/`decade_completeness`
  already ported.
- **Config**: Python silently ignores the `changepoint`/`stats_floor`/`snow`/`drought`/
  `annual_values`/`snow_na_policy`/`filtering` sections. **rpkg's bundled
  `inst/config/signatures_config.json` is stale**: 89 canonical leaf keys vs 60 bundled
  — 29 missing, 0 extra, 0 differing values (Codex-confirmed exactly). Latent fallback
  divergence: `config.R:94-102` defaults `na_reject_negative_flow` to `TRUE`
  (Python/JSON: false) — fix while there. **Python's config path is repo-relative**
  (`config.py:12-23`) — an installed wheel outside the source tree cannot find the
  JSON; package the config as package data with an explicit override/fallback.
- **Python has no `empty_stats`** (inlined NaN loops per call site); factor it out
  before snow/drought, whose all-NaN schema contracts depend on it.
- **Error visibility**: Python orchestrator swallows exceptions with bare
  `except Exception: pass` (`signatures.py:94-183`) and logs nothing; rpkg lacks the
  `gage_id` kwarg for warning context. Both must gain Julia-style per-family warnings
  BEFORE the feature ports. Note (Codex F13): pandas builds the output schema from the
  **union** of per-gage dicts, so a family failing on thousands of gages still creates
  its columns — per-gage key-set assertions are required, not column counts.
- **Benchmark runners**: neither port runner passes any of the six features; neither
  has the Julia runner's ENV overrides. **`run_python_benchmark.py` has NO bounded-
  memory climate strategy** — it loads the full ~98M-row climate parquet and merges it
  wholesale (`run_python_benchmark.py:88-118`), unlike Julia which projects to 4
  columns and releases raw frames (`run_julia_benchmark.jl:97-123, 268-280`). Adding
  SWE grows the payload; a concrete redesign is Phase-0 work (Codex F3, CRITICAL).
  `run_rpkg_benchmark.R` hardcodes all paths (`:49-58`), merges only PPT (`:156-169`),
  and silently drops whole-gage failures (`:143-148, 243-248`).
- **Comparison tooling blind spots** (Codex F5/F6): `compare_*_vs_golden_julia.py`
  compares only the column **intersection** and exits 0 with Julia-only columns
  missing; `build_julia_golden_dashboard.py`'s `SIGNATURE_GROUPS` contains neither
  snow nor drought; `validate_annual_values.py` is **self-referential** (port parquet
  vs the same port's summary). None of these is an acceptance gate as-is — see §4.
- **Tests**: 7 Julia feature test files (~2,025 lines) have no port counterparts;
  Python's `smoke_test.py` contains no `def test_` functions (pytest collects
  nothing); rpkg's golden test only checks function existence.
- **Stale docs**: `docs/CROSS_LANGUAGE_STATUS.md`, CLAUDE.md's "synced April 14-15",
  rpkg DESCRIPTION/README/`signatures.Rd` (~478/551 counts), `python/README.md`,
  comparison-tool hardcoded 623s.

**Corrected Known-Issue assumptions** (survey + Codex both confirm; CHANGELOG updated):
- The seasonal runoff-ratio completeness bug is **Julia-only**; ports mask correctly
  (Python builds `<configured-season>_complete` from the config's full names;
  rpkg likewise). Fixing Julia moves it *toward* the ports.
- `negative_ann` is NA-safe in both ports. rpkg nuance: `aggregate`'s formula
  interface may drop an all-NA water year entirely — pin with a parity test.

## 4. Reference, comparison & acceptance-gate strategy

- **Reference = a fresh Julia benchmark run on this M1**, produced in Phase 0 *after*
  the Julia RR fix, with standard-product #1 settings (WY 1993–2025 @ 60%, drought on,
  rebuilt climate parquet), in its own experiment folder per the one-folder convention.
- Port benchmarks run with the **same window/fraction/inputs**, always on the
  **non-legacy preprocessing path — runners must FAIL FAST if `use_legacy_filtering`
  is true** (Codex F14), and record provenance: package/runtime versions and input
  sizes/hashes (a multi-session campaign is not protected by "same machine" alone).

**Acceptance gates (used at every phase exit; diagnostic tools are NOT gates):**
1. **Exact schema equality** vs the expected column set for that phase (final target:
   exactly 1,653 columns, exactly 100 annual signature bases; FAIL on any
   reference-only or port-only signature). New small checker script — the existing
   comparison scripts' intersection-only behavior and exit-0-on-missing-columns
   disqualify them as gates (they remain diagnostics).
2. **Zero unexpected per-family warnings** + **per-gage expected-key-set assertions**
   (guards the union-schema silent-per-gage-loss mode) + exact expected gage count and
   zero dropped-gage errors in the runner log.
3. **Value gates split by field class** (Codex F12 + omitted-risk): *exact* equality
   for schemas, key sets, NA patterns, counts, discrete fields (`_pettitt_cp_year`),
   and within-language unchanged columns; *tolerance/identity-R² tiers* only for
   genuinely continuous cross-language outputs (slopes, p-values — different stat
   libraries legitimately differ in last bits). Pettitt cp-year parity is proven on
   **deterministic unit fixtures** (tied maxima, unsorted years, repeated values,
   boundary segment sizes, all-constant series), not benchmark tolerances.
4. **Additivity checks scoped correctly**: `check_additivity.jl` for genuinely
   additive phases; for phases that legitimately move existing columns (stats floor,
   `flagged_for_high_na` denominator shifts) the gate is a **transition check** — only
   the expected columns move, and their new values are independently recomputed and
   asserted (never a blanket `--allow-shift`). Per-family finite-coverage expectations
   replace the generic `--min-finite-frac 0.5` for gated/sparse families (snow gated
   metrics, Pettitt segment p-values, per-gage scalars).
- Per the **cross-language-alignment skill**: Julia is canonical — match exactly,
  never "improve"; one feature per phase; benchmark before moving on; revert anything
  that regresses other columns; all parameters from the shared config; NaN==NaN
  comparisons; three-way comparison to localize outliers; maintain the
  known-irreducible list (Spearman p small-n, FDC90th near-zero OLS, FP-tie-sensitive
  rank/Pettitt fields).

## 5. Phases

### Phase 0 — Environment, Julia fixes, reference, baseline
1. **Python env**: uv venv; `pip install -e python` + pytest; unit suite green.
   Package the config JSON as package data (or explicit fallback) so an installed
   wheel works outside the repo.
2. **R env**: install R (arm64 CRAN) + `data.table, arrow, zyp, Kendall, mblm,
   testthat, roxygen2, devtools`; **sync `rpkg/inst/config/signatures_config.json`**
   (29 missing keys) and align the `na_reject_negative_flow` fallback to `false`;
   install rpkg from source; testthat green **as installed package** (not only
   `load_all()`).
3. **Julia fixes** (Julia-first): (a) RR flag names `runoff_ratios.jl:114-117` →
   `win_/spr_/sum_/fal_complete` + a regression test asserting the masking fires;
   (b) wire `test_changepoint.jl` into `runtests.jl` (currently omitted — a green
   suite doesn't run it). Full Julia suite green.
4. **Fresh Julia reference run** (product #1 settings) → own folder;
   `check_additivity.jl` vs product #1 to attribute the delta. Expected changes:
   the four seasonal runoff-ratio stat/Pettitt families, possibly
   `flagged_for_high_na` — **and the seasonal runoff-ratio rows of the annual
   parquet** (masking precedes collection; compare and document those too, Codex F4).
   Annual runoff ratio + `runoff_ratio_high_count` must NOT move.
5. **Runner prep, both ports**: ENV overrides (window/fraction/paths/output dir);
   fail-fast on legacy filtering; provenance block (versions, input sizes/hashes);
   **Python bounded-memory climate redesign** — column-projected parquet read
   (site_id, Date, prcp, swe only), per-gage grouped storage, release raw frames
   (mirror `run_julia_benchmark.jl:97-123, 268-280`); R equivalents (select columns
   before merge; keep failure-drop behavior but make any drop a gate failure).
6. **Baseline port benchmarks** (unmodified port library code) → comparison report.
   Expected: ≈April-level agreement on the 623-column subset; whitelisted divergence
   = b≠1 log_a family (until Phase 2).

Exit gate: baseline comparison documented; both runners on non-legacy path with
provenance; environments reproducible.

### Phase 1 — Stats layer: Pettitt, stats floor, plumbing (split gates)
**1a — Pettitt (additive)**: port `pettitt_test` + `segment_differential_metrics` +
`_run_changepoint_block!` equivalent + `CP_SUFFIXES`; read the `changepoint` config
section; thread `changepoint` through the orchestrators. Python `empty_stats`
factored out; Python logged per-family warnings; rpkg `gage_id` kwarg. Port
`test_changepoint.jl` incl. **exact-equality fixtures** for cp-year determinism.
Gate: additivity vs Phase-0 baseline (only `_pettitt_*` columns added), Pettitt
fixture parity exact, benchmark comparison vs reference on the Pettitt fields.

**1b — Stats floor + `force_skip_trends` (transition, NOT additive)**: kwargs +
config + orchestrator threading with recession/elasticity exemption; direct unit
tests of the `force_skip_trends` mechanism now (its production consumer arrives in
Phase 4 — don't let it ride untested, Codex F10). Port `test_stats_floor.jl`.
Gate: **floor-transition check** — run `docs/benchmarks/apply_stats_floor_mask.py`
on the 1a output and require the 1b output to equal it exactly (independent
recomputation of exactly which fields must become NaN); above-floor fields
bit-identical to 1a; NaN pattern matches the Julia reference.

Exit: **Codex adversarial review of the stats layer** before building on it.

### Phase 2 — b=1 recession alpha
- `log_a_pointcloud` → `median(log(-dQ/dt) − log(Q))`; `log_a_events` via an
  `event_log_a_b1` equivalent; seasonality sinusoid fed per-event b=1 values.
  `b_*`/`concavity`/`alpha_linear` untouched.
- Port `test_recession_alpha_b1.jl`; update the ports' old free-fit tests.

Exit gate: the historically-whitelisted log_a divergence disappears (family
converges to reference); transition check — ONLY the log_a family moves.

### Phase 3 — Annual-values collector
- Python collector class; R environment-based accumulator (no O(n²) append).
  `collector` kwarg collects **before** all gates; thread through all signature
  functions + orchestrators; runners drain to parquet, gated on `annual_values.save`.
- Port `test_annual_collector.jl` (incl. with/without-collector stat identity).

Exit gate: (a) `validate_annual_values.py` per port (internal consistency —
necessary, NOT sufficient: it is self-referential, Codex F6); (b) **new
cross-language annual comparator** vs the Julia reference parquet: exact signature
set, exact `(gage_id, signature, water_year)` key set with documented exceptions
only, NA-pattern equality, row-level value agreement, zero duplicate keys. Base
count at this phase = the pre-snow/pre-drought set; the contract is re-asserted at
90 bases after Phase 4 and 100 after Phase 5 (Codex F10).

### Phase 4 — Snow family
- Preprocessor: `swe`→`SWE` normalization, per-year SWE policy, `valid_swe_years`
  return (additive key — verify consumers). `snow` module + record-anchored decade
  gate via `force_skip_trends`; orchestrator `snow_data`/`snow_climate_years`
  (explicit frame, no implicit fallback); runners select + merge SWE and build
  `snow_data` (neither does today).
- Port `test_snow_metrics.jl` + `test_snow_record_decade_gate.jl`; regenerate rpkg
  NAMESPACE/man + installed-package smoke (per-phase packaging rule, Codex F8).

Exit gate: snow columns match reference with **per-metric expected coverage ranges**
(gated metrics are legitimately sparse — no generic 50% finite rule); snow-bearing
gage count matches reference exactly; transition check — only `flagged_for_high_na`
may move, its new values independently recomputed (`refresh_qa_flags` approach);
annual comparator re-run at 90 bases.

### Phase 5 — Drought family
- `drought` module (`weibull_quantile` type 6 + NA-below-range,
  `smooth_daily_flow`, duration/deficit + threshold scalars); config fail-fast
  validation; `CFG_DROUGHT_ENABLED`-equivalent gate.
- Port `test_drought_metrics.jl` incl. the Sep 30 → Oct 1 boundary-attribution
  fixture and record invariants; rpkg NAMESPACE/man regeneration + installed smoke.

Exit gate: exactly +165 columns; transition check — only `flagged_for_high_na` may
move (denominator counts the new fields; Julia's zero-crossing outcome in July was
window-specific — do not assume it), independently recomputed; drought columns match
reference; annual comparator at 100 bases.

### Phase 6 — Full validation, adversarial review, docs
- Full-window benchmarks (both ports); **exact-schema gate: 1,653 columns, set
  equality with the reference, fail on any one-side-only signature**; three-way
  comparison; comparison dashboards as diagnostics (extend
  `build_julia_golden_dashboard.py`'s `SIGNATURE_GROUPS` with snow + drought first —
  same class of gap as the 2026-08-10 fix); final annual comparators.
- **Codex adversarial review** of the whole campaign.
- Docs sweep: rewrite `CROSS_LANGUAGE_STATUS.md`; CHANGELOG release entry; CLAUDE.md
  table + sync note; python/rpkg READMEs + `signatures.Rd` counts; claude-skill
  (drop "Julia only" markers); pyproject/DESCRIPTION version bumps; comparison-tool
  hardcoded counts.
- Follow-ups to file (not in-scope): Julia `negative_ann` coalesce one-liner + rpkg
  all-NA-year `aggregate` parity test; consider rpkg reading the canonical config
  path like Python.

## 6. Risks & watch items

- **Wall-clock**: Python full benchmark ~150 min, rpkg ~2–4 h; run unattended at
  phase exits, smoke-first. Smoke fixtures must use long-record gages (Pettitt needs
  ≥20 obs, ≥10 per segment — a short-record smoke "failure" is correct behavior).
- **exFAT thumbdrive**: verify input sizes + `PAR1` footers against provenance before
  every long run. Climate input = `daymet_1980_2023_rebuilt_10aug2026.parquet`.
- **16 GB RAM**: Python runner redesign is Phase-0 scope (see §5.0.5) — there is no
  existing chunked join; R runner needs column projection before merge.
- **Config-at-import**: Python freezes config into module constants at import — keep
  runner overrides runtime-read, as the Julia runner does.
- **R stats NA-swallowing** (`rpkg/R/stats.R:89-113` turns library errors into NA):
  port tests must distinguish legitimate degenerate-series NA from library failure.
- **Discretely FP-sensitive fields**: even same-machine, cross-language library
  differences flip ties exactly like cross-machine FP — hence the §4 split between
  exact-equality gates and tolerance tiers.
- **rpkg install**: arrow/data.table CRAN arm64 binaries exist; `nanoparquet` is the
  fallback if arrow compiles painfully (decide only if hit).

## 7. Effort estimate

~1,004 lines of genuinely new module code + ~250–300 lines of stats-layer additions +
~120 lines io + ~60 lines orchestrator, per language (rpkg ~40% of Julia's density,
Python ~85%); ~2,025 lines of Julia feature tests to port; plus the Phase-0 runner
redesigns and the new gate tooling (schema-equality checker, cross-language annual
comparator, floor-transition check). Realistic: 3–4 working sessions for Python,
2–3 for rpkg, plus benchmark time and the docs/review session.

## 8. Codex adversarial review record (2026-08-24, pre-implementation)

Read-only review of this plan against the code. **Initial verdict: NO-GO** — 3
CRITICAL, 11 MAJOR, 7 MINOR/verified + 5 omitted risks. All findings accepted and
incorporated above; dispositions:

| # | Sev | Finding | Resolution |
|---|---|---|---|
| 1 | CRIT | Phase-1 gate ("pre-existing columns unchanged") logically impossible once the stats floor lands — the floor NaNs existing below-floor fields | Phase 1 split into 1a (Pettitt, additive gate) and 1b (floor, transition gate verified against `apply_stats_floor_mask.py` output) |
| 2 | CRIT | Phase-5 "exactly +165, 0 changed" impossible — `flagged_for_high_na` denominator counts the new fields | Gate now allows only that flag to move AND independently recomputes it (blanket `--allow-shift` rejected) |
| 3 | CRIT | Plan claimed a chunked Python climate join exists — it does not (`run_python_benchmark.py:88-118` wholesale merge); SWE grows the payload | Concrete bounded-memory redesign added to Phase-0 scope |
| 4 | MAJ | RR-fix delta expectation incomplete — the annual parquet's seasonal RR rows change too (masking precedes collection) | Phase 0 §4 now compares CSV + annual rows; annual RR/`runoff_ratio_high_count` asserted unmoved |
| 5 | MAJ | Comparison scripts pass on column intersection, exit 0 with missing schemas; dashboard `SIGNATURE_GROUPS` lacks snow/drought | New exact-schema-equality gate (§4.1); dashboards demoted to diagnostics; groups extended in Phase 6 |
| 6 | MAJ | `validate_annual_values.py` is self-referential | New cross-language annual comparator required at Phases 3/4/5 (§5.3 gate) |
| 7 | MAJ | Smoke gates weak/uninterpretable (Pettitt min-obs; generic finite thresholds wrong for gated families) | Long-record smoke fixtures; per-family expected coverage ranges replace generic thresholds |
| 8 | MAJ | rpkg NAMESPACE/roxygen deferred to Phase 6; `load_all()` masks installed-package failures | Regeneration + installed-package smoke in every phase adding public R functions |
| 9 | MAJ | R runner: hardcoded paths, PPT-only merge, silent whole-gage failure drops | Phase-0 runner prep expanded; zero-error + exact gage count in every gate |
| 10 | MAJ | Collector contract must be re-asserted at 90/100 bases after snow/drought; `force_skip_trends` untested until Phase 4 | Annual-base recounts added to Phases 4/5; direct mechanism tests in 1b |
| 11 | MAJ | "BIC dead code" overstated (exported API); `test_changepoint.jl` not in `runtests.jl` | Rephrased ("exported Julia-only API, excluded from parity scope"); runtests wiring added to Phase 0 |
| 12 | MAJ | Pettitt cp-year cannot be tolerance-gated — ties flip the year and all segment fields | Exact-equality deterministic fixtures (ties, unsorted, constant, boundary sizes); tolerances only for continuous outputs (§4.3) |
| 13 | MAJ | Union-built pandas schema hides per-gage silent column loss | Per-gage expected-key-set assertions + zero-unexpected-warnings gate (§4.2) |
| 14 | MAJ | Legacy-filtering branches in port runners could silently bypass modern preprocessing | Runners fail fast unless non-legacy (§4, §5.0.5) |
| 15–19, 21 | — | Verified correct: `generate_stats` signatures, collector placement, recession characterization, `negative_ann` NA-safety (+ rpkg all-NA-year nuance), rpkg config staleness (89 vs 60 leaf keys exactly), 1,653/+165 arithmetic | No change |
| 20 | MIN | "Zero grep hits" rhetorically overstated | Rephrased in §3 |
| — | risk | Python config path repo-relative (installed wheel breaks); R stats NA-swallowing; exact-equality wrongly implied for cross-language p-values; no version/hash provenance; comparison exit codes not acceptance signals | All added: §5.0.1, §6, §4.3, §5.0.5, §4.1 |

Post-revision status: plan proceeds to implementation under the amended gates; a
follow-up Codex review is scheduled at the Phase-1 exit (stats layer) and at Phase 6.

## 9. Execution log

### Phase 0 (2026-08-24, first session)

Completed:
1. **Python env** — uv venv (py3.12) at repo `.venv` (gitignored),
   `pip install -e './python[dev]'`, 26 unit tests green.
2. **rpkg config** — bundled copy synced byte-identical to canonical (21
   sections); `config.R` `na_reject_negative_flow` fallback TRUE → FALSE.
3. **Julia fixes** — RR flag-name fix + `test_runoff_ratio_seasonal_mask.jl`
   (32 assertions, proves the mask fires); `test_changepoint.jl` wired into
   `runtests.jl`. Suite: **2,798 assertions, 0 failures**.
4. **Reference run** — new wrapper
   `run_julia_benchmark_portref_1993_2025_60pct.jl` →
   `/Volumes/Untitled/processedOuts_portref_24aug2026/` (6,678 × 1,653;
   annual 18,898,406 rows / 100 signatures; 10.9 min; inputs hashed).
   **Delta vs product #1 fully attributed** (RUN_NOTES.md + two validation
   files in the folder): exactly the documented 98-column climate-input
   residual (1e-18…6.3e-13; max 2.56e-13 over 729 of 18.9M annual rows) —
   product #1 predates the daymet rebuild. **The RR fix fires ZERO times in
   this window** (no qualifying year has an incomplete season; the >3-day-gap
   rejection makes that nearly impossible) — Codex F4's expectation of
   seasonal-RR annual shifts resolved as "none in this window", proven rather
   than assumed. Reference VALIDATED; ports share its climate input, so the
   residual never enters port validation.
5. **Runners rebuilt** (both): Julia-mirrored ENV names, runtime
   window/fraction gates (Julia semantics copied verbatim), legacy fail-fast
   (`STREAMFLOW_ALLOW_LEGACY=1` escape hatch), bounded-memory climate handling
   (Codex F3), provenance blocks with package versions, prefixed outputs;
   R runner exits nonzero on any errored gage (Codex F9). Python runner
   verified by synthetic end-to-end smoke (window, fraction, zero-padded
   metadata join, area gate). Two more schema gaps closed at runner level:
   rpkg gains the 11 GAGES-II/HYDAT interference columns +
   `human_interference_class` (previously absent); Python fills
   RHBN/REGULATED + hic from `metadata/canadian_hydat_interference.csv`
   (stale "not available from Python" removed). NOTE: the Python Phase-5
   metadata merge now joins on the leading-zero-normalized id (mirroring
   Julia) — the old raw-string join silently missed zero-padded gages'
   lat/lon/basin_area.
6. **Gate tooling** — `docs/benchmarks/check_schema_equality.py` (strict
   column-set + gage-set equality, explicit logged waivers; self-tested).

7. **Python baseline — COMPLETE**: 59.3 min (was ~150), 6,678 gages —
   gage set IDENTICAL to the reference (schema gate PASS, 0 either-only;
   the 997 missing columns = exactly 608 Pettitt + 224 snow + 165 drought);
   623 shared columns: 575 Perfect / 26 Good / 22 Extremely Low, the latter
   ENTIRELY the intentional b≠1 log_a family (Phase-2 target). Expected
   pre-port baseline exactly. Evidence in the reference folder's RUN_NOTES.

### Phase 1 — Python side (2026-08-24, same session)

**Implemented (1a Pettitt + 1b floor/force_skip_trends), 53/53 unit tests:**
- `python/streamflow_signatures/changepoint.py` — `pettitt_test` +
  `segment_differential_metrics` (BIC path excluded per plan §3). **Fixture
  cross-check vs Julia: 7 deterministic fixtures (clean step, tied-maximum
  "tent", constant≠0, constant=0, boundary n=20, n=19→NaN, NaN-interspersed)
  produce IDENTICAL cp_year everywhere — including the tie (first-max rule) —
  bit-identical p-values, means within 2 ulp.** Frozen into
  `python/tests/test_changepoint.py` (16 tests) with unsorted-input invariance,
  window filtering, independence-from-trend-gating, and below-min_rows CP-NaN
  contract.
- `stats.py`: `STAT_SUFFIXES`/`CP_SUFFIXES`, `empty_stats` (8 keys only —
  mirrors Julia; CP keys surface via union schema), `changepoint` +
  `min_values_for_stats` + `force_skip_trends` kwargs with Julia's exact gate
  ordering (floor shares the min_rows branch incl. CP NaN; forced skip goes
  through the trend-suppression path; CP independent of trend gates).
- Kwarg sweep: `changepoint` in all 13 module functions + 16 forwardings
  (matches Julia's per-module counts); `min_values_for_stats` in the 9
  non-exempt modules only — **recession/elasticity structurally exempt (their
  signatures don't accept the kwarg, pinned by TypeError tests)**, as in Julia.
- `config.py`: `changepoint` + `stats_floor` sections (Julia-matching
  fallbacks); orchestrator rewritten with per-family LOGGED warnings +
  `gage_id`/`changepoint`/`min_values_for_stats` kwargs; runner passes both.
- `python/tests/test_stats_floor.py` (11 tests) incl. the 07292500-shaped
  clustered-series regression and direct `force_skip_trends` mechanism tests
  (Codex F10).

**Gates prepared** (`scratchpad/phase1_gate/run_phase1_gates.py`): A schema
(+608 CP exactly, complete 8-suffix sets per base, gage-set equality); B floor
transition (phase-1 == `apply_stats_floor_mask.py`-masked baseline EXACTLY —
mask counts from the Julia reference annual parquet; 7,234 expected value→NaN
transitions; flags excluded, recomputed in-run); C cp_year exact agreement
≥99.9% + NA-pattern vs the reference. Note: `flagged_for_high_na` is NOT
reference-comparable until Phase 5 (denominator counts 1,653 vs 1,231
columns) — excluded from gates, reported informationally.

**History**: the first 1a-only benchmark was externally stopped mid-run; on
user go-ahead relaunched as ONE combined Phase-1 run (Pettitt + floor), with
both effects independently attributable (additivity handles the +608; the
mask script independently recomputes the floor's NaN set).

In progress / pending:
- Combined Phase-1 Python benchmark running → gates A/B/C + 6-tier comparison.
- R install requires one sudo command from the user (CRAN R 4.6.1 arm64
  downloaded); then `setup_r_env.sh` (deps, `R CMD INSTALL rpkg`,
  installed-package testthat), rpkg Phase-0 baseline, then rpkg Phases 1a/1b.
- Codex stats-layer review at the Phase-1 exit (after gates pass).
- rpkg gap carried to Phase 1: `na_cause_ice` preprocessor tracking behind
  `ice_affected_days_total`.
- R runner syntax/behavior untested until R exists (no R runtime on this
  machine yet).
