---
name: add-signature
description: Checklist for adding or changing a streamflow signature or signature family — Julia canonical first, the registries that fail silently, the benchmark additivity proof, the Python and rpkg ports, and the docs to update. Use when asked to add, implement, rename, remove, or change the method of a signature or metric.
---

# Add or change a signature

Julia is canonical: every change starts in `julia/src/` and is then ported. The
guidelines doc (`docs/SIGNATURE_GUIDELINES.md`) is the methodology ground truth — if the
definition being implemented is not there yet, say so; the user relays wording to the
co-authors (`/sync-docs`).

## 1. Julia

1. Function in the appropriate `julia/src/*.jl` module, returning per-water-year annual
   values (DataFrame with `water_year` + metric columns). Use the preprocessed frame as
   given: no NA filling, no per-year min-days checks (rule `signatures-code`).
2. `generate_stats(annual_df, [:metric], :water_year; collector=collector)` — exactly
   the 8 statistics + 8 Pettitt fields. Thread the `AnnualCollector` through. A per-gage
   scalar that is not an 8-stat output is an exception and must be documented as one.
3. Register the call in `calculate_all_signatures()` (`julia/src/signatures.jl`) inside
   the per-family `try/catch` pattern; honour `area_normalized` if the metric uses PPT.
4. Config knobs go in `config/signatures_config.json` AND its two bundled copies
   (`python/streamflow_signatures/data/signatures_config.json`,
   `rpkg/inst/config/signatures_config.json`). `CFG_*` constants are baked at
   PRECOMPILE time — purge `~/.julia/compiled/<ver>/StreamflowSignatures` after a
   config edit (`/run-benchmark` explains the probe).

## 2. Registries — each one fails loudly only after the fact

- `EXPECTED_SIGNATURE_BASES` in `config.R` (8-stat bases); per-gage scalars need their
  own `EXPECTED_*` constant wired into `validate_output_schema()`.
- `EXPECTED_DENSE_SIGNATURES` in `julia/test/test_annual_collector.jl` — asserts SET
  EQUALITY of the collected annual series.
- The signature-count gate in `docs/benchmarks/validate_production_run.py`
  (`ann.signature.nunique() == N`).
- `qa_qc.high_na_denominator` in the config — register any new per-gage scalar (8-stat
  columns are matched by suffix).
- `docs/signature_categories.csv` — one row per output (category, function, module).
- Python `validate_schema()` and the rpkg schema tests, once ported.
- The expected TOTAL column count in the docs (1,653 today): 16 fields per annual
  signature PLUS every non-8-stat scalar (the drought family shipped documented as +160
  when it is +165).

## 3. Verify — a green unit suite proves nothing here

1. `julia --project=julia julia/test/runtests.jl`
2. Benchmark via `/run-benchmark` (~27 min) into a NEW experiment folder.
3. Prove additivity: `julia --project=julia docs/benchmarks/check_additivity.jl NEW.csv
   PREVIOUS.csv --expect-added N` — every pre-existing column unchanged (only
   `flagged_for_high_na` may legitimately shift) AND the new columns populated. The
   orchestrator's `try/catch` turns a failure on one production gage into silently
   missing columns; smoke checks must assert the new values are FINITE, not merely present.
4. `python docs/benchmarks/check_signature_failures.py <run log>` — no swallowed
   family exceptions.

## 4. Port

- Python: function in `python/streamflow_signatures/<module>.py`, called from
  `signatures.py`; fixture cross-check against Julia to ≤ 4e-14; tests in `python/tests/`.
- rpkg: function in `rpkg/R/<module>.R`, called from `signatures.R`, exported in
  `NAMESPACE`; tests in `rpkg/tests/testthat/`; `R CMD INSTALL rpkg` before any benchmark
  (the package bundles the config).
- Gates: `check_schema_equality.py`, `check_annual_parquet_equality.py`,
  `check_signature_failures.py` (`run_rpkg_acceptance_gates.sh` runs all four for rpkg).
  The `compare_*_vs_golden_julia.py` scripts are diagnostics only. Divergences:
  `/cross-language-alignment`.

## 5. Document

- `docs/SIGNATURES.md`: the family section (metrics table, method, units, caveats), the
  Summary Table row, and the Overview category counts.
- `README.md` Signature Categories table; `claude-skill/streamflow-signatures.md` (the
  user-facing interpretation skill — keep its category list and caveats current).
- `CHANGELOG.md` entry (rule `changelog`); `docs/STATUS.md` if a delivered product is
  affected.
- After the merge to main: `./publish_to_hisss.sh`.
