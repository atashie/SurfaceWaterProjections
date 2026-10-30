# Sonnet review — code and documentation currency after the default flip, area JULIA (2026-10-29, blind, read-only)

Brief: `review_brief_currency.md` (same folder). Verdict GO WITH CHANGES; every finding applied 2026-10-29 (dispositions in `../2026-10-21-qualifier-mask-fix-rerun-plan.md` §11, row "currency review"). Verbatim report:

# Julia review (currency after the default flip), HEAD 713fdbf, read-only

## Verdict: GO WITH CHANGES (no BLOCKER, no MAJOR; 4 MINOR doc/record items, 3 NOTES)
The Julia side does what the docs say. `julia/` is unchanged since 54f0483 (`git diff 54f0483..HEAD --stat -- julia/` is empty). The helper, the runner, all 12 wrappers, the guards and the input hashes check out. The unit suite and the self-test pass.

## Findings
| Rank | file:line | Evidence | Fix |
|---|---|---|---|
| MINOR | docs/benchmarks/run_batched_julia.py:465 | The stage record writes `os.environ.get("STREAMFLOW_CLIMATE_FRAME", "outer (default)")`. When `--climate-frame` is not passed and the wrapper is one of the ten `left`-pinned ones (e.g. the drought wrapper), the wrapper's `get!` makes Julia run `left`, but the stage record says "outer (default)". Line 273 `climate_frame_env` is None in that case. Claim 3, "records the resolved default", holds only for the `*_climate` wrappers. The per-batch timing JSON `provenance.climate_frame` is correct. | Say "unset (the wrapper pins or the helper default applies; see the batch timing JSON)", or read the frame from the batch timing JSON. |
| MINOR | .claude/skills/run-benchmark/SKILL.md:63-64 | The Launch block labels the `*_climate.jl` wrappers "product #1 VERSION 2 (climate family, reprocessed Daymet; local verified inputs)". They are now the version-3 wrappers (corrected input, prefix `..._climate_v3_26oct2026`, public name "version 2"), and #1/#2 in line 61-62 are the drought wrappers that pin `left`. | Relabel: "product #1 v3 (public 'version 2'; corrected input, frame outer; run through run_batched_julia.py)". |
| MINOR | docs/DEVELOPMENT.md (benchmark-files table, row for `prod_{1993,1980}_2025_60pct_climate.jl`) | The row says "the version-2 prefix and inputs until 2026-10-28". The wrappers now carry the v3 prefix and the corrected inputs; the sentence is ambiguous history in a current-state table. | Replace with "prefix `..._climate_v3_26oct2026`, inputs from `qualifier-rerun-26oct2026/input/`". |
| MINOR | julia/README.md:37 (also python/README.md:37, mirrored) | The usage example builds climate with `leftjoin(df, climate[...], on=[:gage_id,:date])`. That is the `left` frame, which the project calls a defect (a day absent from the streamflow record loses its Daymet value). There is no note beside the example. | Add one comment line: "this joins climate on the streamflow dates; the production runners use the outer frame (docs/benchmarks/climate_frame.jl)". |
| NOTE | docs/benchmarks/run_julia_benchmark.jl:342-343 | The legacy branch (`CFG_USE_LEGACY_FILTERING`) still uses a plain `leftjoin`, ignoring `STREAMFLOW_CLIMATE_FRAME`. The banner at line 74 and the timing JSON still report the resolved frame (outer by default). The rpkg runner handles this explicitly (`if use_legacy "left"`, run_rpkg_benchmark.R:289). Production uses the non-legacy branch (line 231-233, the helper). | Optional: record `climate_frame = "left"` in the provenance when `use_legacy`. |
| NOTE | wrappers `*_climate.jl` lines 69-70 / 70-71 | The guard fires only when `STREAMFLOW_OUTPUT_DIR` is unset. Explicitly setting it to the product folder is not refused, and an empty product folder passes. This is by design for the batched runner, which sets the directory per batch. | None. |
| NOTE | docs/benchmarks/audit_qualification.jl | It is Q-only (`preprocess_daily_data` on the streamflow parquet, no climate join). It has no dependence on the frame or the qualifier rule, and its frozen-cohort contract text is consistent with the docs. The default `SF` is still the Windows Feb-2026 path, but that is pre-existing and overridable by argument 4. | None. |

## Brief questions, Julia parts
- A (stale text): no current-state statement of the old behaviour in my scope. The sweeps for `left (unset)` and "default ... left" found only dated or historical hits (CHANGELOG.md:523 states the new default). `keep_codes`, `USGS_ACCEPTED_CODES`, `A e` and `P e` do not appear in `julia/src`, `julia/test` or `julia/README.md`. The library docstrings do not mention a climate join other than the preprocessor's own date-grid `leftjoin` (`io.jl:297`), which is independent of the frame. The stale items are the labels listed in the findings.
- B (code and docs agree): confirmed.
  - `climate_frame_mode()` (climate_frame.jl:37-43) defaults to "outer" and errors on any other value.
  - The runner calls it at run_julia_benchmark.jl:51 and uses the helper at 231-233.
  - The runner prints the frame at line 74 and records it at line 722 (`climate_frame`) and line 727 (`env_overrides`, only if the variable is set).
  - The self-test checks the default (selftest_climate_frame.jl:40-43).
- C (tests): see the verbatim outputs below.
- F: the verdict is above.

## Per-wrapper table
All 12 pin with `get!(ENV, "STREAMFLOW_CLIMATE_FRAME", ...)`. No wrapper hard-sets the frame. Every pin comes before the `include(.../run_julia_benchmark.jl)` line, and the include is the last statement.

| Wrapper | Frame pinned | get! | Pin before include | Docstring / comment matches env block |
|---|---|---|---|---|
| run_julia_benchmark_drought_1993_2025_60pct.jl (l.46) | left | yes | yes | yes (comment l.43-45) |
| run_julia_benchmark_portref_1993_2025_60pct.jl (l.46) | left | yes | yes | yes |
| run_julia_benchmark_prod_1980_2025_60pct.jl (l.24) | left | yes | yes | yes |
| run_julia_benchmark_prod_1980_2025_60pct_drought.jl (l.51) | left | yes | yes | yes |
| run_julia_benchmark_prod_1980_60pct.jl (l.21) | left | yes | yes | yes |
| run_julia_benchmark_prod_1993_2025_60pct.jl (l.27) | left | yes | yes | yes |
| run_julia_benchmark_prod_1993_60pct.jl (l.18) | left | yes | yes | yes |
| run_julia_benchmark_startIn1993.jl (l.15) | left | yes | yes | yes |
| run_julia_benchmark_startIn1993_60pct.jl (l.12) | left | yes | yes | yes |
| run_julia_benchmark_startIn1993_80pct.jl (l.12) | left | yes | yes | yes |
| run_julia_benchmark_prod_1993_2025_60pct_climate.jl (l.47) | outer | yes | yes | yes (docstring l.14, comment l.45-46) |
| run_julia_benchmark_prod_1980_2025_60pct_climate.jl (l.48) | outer | yes | yes | yes (docstring l.14) |

Ten pre-rerun wrappers are `left`, two `*_climate` wrappers are `outer`: this matches claim 3.

### The two `*_climate.jl` wrappers
- **Guard test.** I extracted the lines before `include(` into scratch files and ran them with `STREAMFLOW_OUTPUT_DIR` unset. Both stopped with: `ERROR: LoadError: refusing to run: /Users/arik/HISSS_data/qualifier-rerun-26oct2026/processedOuts_{1993,1980}_2025_climate_v3_26oct2026 holds the delivered product — set STREAMFLOW_OUTPUT_DIR to run elsewhere`. With `STREAMFLOW_OUTPUT_DIR` set to a scratch path, the guard does not fire and `STREAMFLOW_CLIMATE_FRAME` resolves to `outer`.
- **Prefix and docstrings.** The prefixes are `streamflow_{1993,1980}_2025_60pct_climate_v3_26oct2026`. The docstring gage counts are 6,678 and 6,250; the would-qualify counts are 31 and 29; the untouched-gage counts are 5,069 and 4,894.
- **Input paths.**
  - The corrected streamflow parquet exists. It is 628,450,987 bytes, and `shasum -a 256` gives `2d44c1afd9a88597b0f4037ba1bce17a98bc832764c396f934b3eeaded1ec7a0`. That equals `input/provenance.json` (`outputs.patched_parquet`; the docstring's "2d44c1af…" prefix is also right).
  - The patched metadata CSV exists. It is 2,130,953 bytes, and `shasum -a 256` gives `ca335b2d67fc77f98645ce9d3bf20792e628d41aaf09d985767ae4297b8ad36d`. That equals `outputs.metadata_patched`.
  - The Daymet parquet `daymet-processed-29sep2026/daymet_1980_2025_29sep2026.parquet` is 5,134,193,404 bytes. Its md5 is `c059076088cd8abd963c64094468e90f`, matching the wrapper comment.
  - The GAGES-II directory `daymet-replay-01oct2026/inputs/gagesMetadata` exists.
  - Both product folders exist.

## Verbatim test outputs
- **Julia unit suite.** `~/.juliaup/bin/julia --project=julia julia/test/runtests.jl` exited with no error and needed no cache purge. The 13 testsets, as printed, all pass with Pass = Total:
  - StreamflowSignatures.jl 78/78
  - Annual Collector 639/639
  - Recession alpha with fixed b = 1 60/60
  - Area-normalized gate on Q-to-PPT signatures 27/27
  - Snow metrics 377/377
  - Stats floor 98/98
  - Snow record-anchored decade gate 96/96
  - Drought signatures 1328/1328
  - Changepoint Detection 39/39
  - generate_stats() changepoint integration 27/27
  - Seasonal runoff-ratio completeness masking 32/32
  - flagged_for_high_na denominator 14/14
  - Climate metrics (guidelines 3.9) 369/369

  The sum is 3,184 passed, 0 failed, 0 errored. The suite prints no grand total, so this sum is mine.
- **Self-test.** `selftest_climate_frame.jl <scratch>` exited 0 and its last line is `41 passed, 0 failed`. The documented 41 is the real count. I could not check the Python and R counts (38 and 43) or 28 or 37; they are outside my scope.

## Verified versus not verified
**Verified**
- `julia/` is unchanged since 54f0483.
- The helper default is `outer`, and the runner resolves it through `climate_frame.jl`.
- The non-legacy path uses the helper.
- The runner provenance records `climate_frame` and `env_overrides`.
- The self-test passes 41/41.
- All 12 wrappers, with pins, `get!`, ordering and docstrings.
- The two product-folder guards fire.
- The input sizes, the shas and the md5 match.
- The whole-repo grep sweeps. The sweep for `"left"` was dominated by pandas and `merge` hits; the relevant ones are listed above. `docs/plans/` and `R/worktree/` were excluded from the sweep.
- audit_qualification.jl is frame-independent.

**Not verified**
- The Python and R code and test counts.
- The Python-side behaviour of `run_batched_julia.py` beyond the lines cited.
- Mirror contents.
- The staged deposit.
- The R ingestion chain.
- The row counts, which I did not recompute.
- The sweep of `docs/plans/` and CHANGELOG.md for stale statements (excluded; I only checked CHANGELOG for "left (unset)" and "default left").
