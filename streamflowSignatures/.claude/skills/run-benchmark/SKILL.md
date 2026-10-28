---
name: run-benchmark
description: Run and validate a full Julia signature-extraction benchmark, experiment, or standard-product run (WY 1993–2025 / WY 1980–2025 @ 60 %) — the precompile-cache and truncated-input gotchas, the runtime env overrides, the one-folder-per-run convention, provenance, and the post-run validation gates. Use when asked to run, rerun, replay, or validate a benchmark, an experiment, a config variant, or a production product in any of the three languages.
---

# Run and validate a benchmark

## Before launching

1. **Inputs live on an exFAT drive and can be silently truncated with an unchanged
   mtime.** Check byte sizes against the `provenance` block of the most recent timing
   JSON and check that each parquet still ends in `PAR1`:
   ```bash
   for f in combined_streamflow_data_09feb2026.parquet daymet_1980_2023_rebuilt_10aug2026.parquet; do
     p=/path/to/processedOuts_feb2026/$f
     echo "$f $(stat -c%s "$p" 2>/dev/null || stat -f%z "$p") footer=$(tail -c 4 "$p")"
   done
   ```
   The canonical-named `daymet_1980_2023.parquet` is TRUNCATED and left in place.
   **Since 2026-10-12 the climate input is the reprocessed
   `daymet_1980_2025_29sep2026.parquet`** (8,017 basins, calendar 1980–2025, six
   variables; 5,134,193,404 bytes, md5 `c059076088cd8abd963c64094468e90f`; verified local
   copy `~/HISSS_data/daymet-processed-29sep2026/`, md5-verified copy on the drive). The
   runner keeps SIX climate columns (`prcp, swe, tmin, tmax` + keys) for the climate
   metrics; `run_batched_julia.py split` writes the same six. The 2026 products used
   `daymet_1980_2023_rebuilt_10aug2026.parquet` (rebuild recipe:
   `docs/benchmarks/convert_daymet_csvs_to_parquet.py`; Daymet uses a 365-day calendar).
2. **A config variant does NOT take effect on a precompiled package.** Every `CFG_*`
   constant is evaluated at precompile time; an env-var change does not invalidate the
   cache and neither does `touch` (Julia checks file CONTENT). Purge, then probe:
   ```bash
   rm -rf ~/.julia/compiled/v1.12/StreamflowSignatures      # adjust the version dir
   STREAMFLOW_CONFIG=/path/variant.json julia --project=julia \
     -e 'using StreamflowSignatures; println(CFG_DROUGHT_ENABLED)'
   ```
   Purge again afterwards, or later "normal" runs silently keep the variant. A variant
   result without a purge, a probe, or an observed expected delta is untrustworthy.
3. **Standard products: clean working tree** (or retain the diff) and
   `STREAMFLOW_HASH_INPUTS=1`. Both delivered products recorded
   `git_working_tree_dirty = true` and cannot be tied to a source tree — do not repeat that.
4. Runtime env overrides, read inside `main()` and safe to set per run:
   `STREAMFLOW_DATA_PATH`, `STREAMFLOW_METADATA_PATH`, `STREAMFLOW_CLIMATE_PATH`,
   `STREAMFLOW_GAGES_II_DIR`, `STREAMFLOW_OUTPUT_DIR`, `STREAMFLOW_OUTPUT_PREFIX`,
   `STREAMFLOW_START_WATER_YEAR`, `STREAMFLOW_END_WATER_YEAR`,
   `STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION`, `STREAMFLOW_HASH_INPUTS`,
   `STREAMFLOW_CLIMATE_FRAME` (`left` | `outer`, 2026-10-26: `left` joins climate onto the
   streamflow dates only — every product before the qualifier rerun; `outer` also takes
   the gage's Daymet rows on absent streamflow dates inside the record span and the
   gage's water years — `docs/benchmarks/climate_frame.{jl,py,R}`; the resolved mode is
   in the timing JSON `provenance.climate_frame`; default `outer` since 2026-10-28 (user
   decision; `left` reproduces every pre-rerun product); every product wrapper pins its own
   frame with `get!`, so the batched runner's `--climate-frame` still wins).
   `STREAMFLOW_CONFIG` is precompile-time (item 2). rpkg reads
   `STREAMFLOW_SIGNATURES_CONFIG` instead (deferred fix — `docs/STATUS.md`).
5. rpkg: `R CMD INSTALL rpkg` first — the package bundles the config.

## Launch

```bash
julia --project=julia docs/benchmarks/run_julia_benchmark.jl                               # generic, ~27 min
julia --project=julia docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl       # standard product #1 wrapper
julia --project=julia docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct_drought.jl  # standard product #2 wrapper (rebuilt climate parquet by default)
julia --project=julia docs/benchmarks/run_julia_benchmark_prod_1993_2025_60pct_climate.jl  # product #1 VERSION 2 (climate family, reprocessed Daymet; local verified inputs)
julia --project=julia docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct_climate.jl  # product #2 VERSION 2
python docs/benchmarks/run_python_benchmark.py                                             # Python port
Rscript docs/benchmarks/run_rpkg_benchmark.R                                               # rpkg port
```

On the 16 GB laptop a full Julia run needs more memory than is usually free. Run it in gage
batches with `docs/benchmarks/run_batched_julia.py`, using a Python that has duckdb:
- `split` records the sources' sha256;
- `run` kills Julia's whole process group when memory runs low or a batch hangs;
- `merge` refuses an incomplete or inconsistent set of batches and orders the gages like
  `--reference`, or by gage id without one.
  Pass `--memory-limit 3GB` for the 46-year product (34.7 M annual rows; the 1GB default
  ran out of memory 2026-10-13) and read `merge.log` — a chain's echoed `$?` after `$(date)`
  is the date command's status.
Self-test: `docs/benchmarks/selftest_run_batched_julia.py`.

**Every artifact of a run goes in that run's OWN folder** (`processedOuts_<experiment>_<date>`):
signatures CSV, annual parquet, timing JSON (with its provenance block), run log, signature
explorer + `_annual/` sidecar, and every comparison dashboard/CSV/summary produced for it.
`docs/benchmarks/` keeps only the tools and the long-lived April reference CSVs.

## After the run

1. `python docs/benchmarks/validate_production_run.py` (see `--help`) — column count,
   signature count (139 since the climate family; `--expect-signatures 100` for an older
   run), gage set.
2. `python docs/benchmarks/validate_annual_values.py` — annual parquet vs summary CSV.
   Self-referential: necessary, not sufficient.
3. `python docs/benchmarks/check_signature_failures.py <run log>` — swallowed per-gage
   family exceptions resurface as ordinary NA, invisible to every column-based check.
4. When columns were added: `julia --project=julia docs/benchmarks/check_additivity.jl
   NEW.csv OLD.csv --expect-added N`. Only `flagged_for_high_na` may legitimately move.
   Cross-MACHINE diffs flip rank statistics on last-bit ties (`FDC90th_spearman_pval`
   moved 0.81 while the annual series agreed to 5.7e-14) — a bare FAIL there needs
   interpretation, not a fix.
5. Port runs: `check_schema_equality.py` (strict; waivers named on the command line),
   `check_annual_parquet_equality.py`, `run_rpkg_acceptance_gates.sh`;
   `compare_*_vs_golden_julia.py` and `build_julia_golden_dashboard.py` are diagnostics.
6. Explorer: `docs/benchmarks/build_signature_explorer.py`, output into the run folder.
6b. A rerun on a new climate input: `docs/benchmarks/compare_climate_rerun.py --new NEW.csv
   --old OLD.csv [--new-annual --old-annual]` classifies every column (Q-only must be
   identical; climate-dependent reported with R² and gained/lost cells; new columns
   populated + sniff ranges; annual values on shared gage-years).
7. **If ANY portion of the data was rerun**: regenerate `flagged_for_high_na` in both
   delivered products (`docs/benchmarks/recompute_high_na_flag.py --write`), update the
   HydroShare READMEs and dictionary row, and rerun the rpkg benchmark so its reference
   CSV picks up the constant-series Mann-Kendall fix (user decision 2026-09-04).
8. Record: CHANGELOG entry (severity, counts, folder, provenance) and the `docs/STATUS.md`
   product table if a product changed.
