---
paths:
  - "docs/benchmarks/**"
  - "julia/src/config.jl"
  - "julia/src/io.jl"
---

# Benchmarks, experiments and standard products

- `docs/benchmarks/` holds TOOLS and the long-lived cross-language reference CSVs only.
  Run artifacts (signatures CSV, annual parquet, timing JSON, log, explorer, dashboards,
  comparison CSVs) go in the run's own folder, one folder per experiment (user
  convention, 2026-07-22).
- Two delivered standard products, both 1,653 columns at 60 % qualifying fraction:
  WY 1993–2025 (`processedOuts_drought_28jul2026`, 6,678 gages) and WY 1980–2025
  (`processedOuts_1980_2025_11aug2026`, 6,250). Neither is a subset of the other, and
  record-dependent signatures (drought thresholds, `*_all` pulses, elasticity,
  parameterized BFI) are never compared across them nor re-aggregated onto another
  window. Never rewrite a delivered product without an explicit user decision
  (`docs/STATUS.md` records the `flagged_for_high_na` decision).
- `STREAMFLOW_CONFIG` is read at PRECOMPILE time; the window, fraction, path and output
  overrides are read at runtime inside `main()`. A config-variant result obtained
  without a cache purge, a probe, or an observed expected delta is untrustworthy.
- Climate input is `daymet_1980_2023_rebuilt_10aug2026.parquet`; the canonical-named
  file is truncated. Verify parquet byte sizes and `PAR1` footers against the last timing
  JSON's `provenance` block before any long run.
- Gates vs diagnostics: `check_schema_equality.py`, `check_annual_parquet_equality.py`,
  `check_signature_failures.py` and `check_additivity.jl` are GATES (non-zero exit,
  waivers named on the command line). The `compare_*` scripts intersect columns and
  exit 0 on a missing family — diagnostics only. `validate_annual_values.py` is
  self-referential (necessary, not sufficient).
- Dashboards and explorers group by FUNCTION FAMILY (the 15 computing functions); the
  eight signature CATEGORIES come from `docs/signature_categories.csv`. Keep every
  `SIGNATURE_GROUPS` list complete when a family is added — the snow family was silently
  absent from every dashboard before 2026-08-10.
- Procedure, env vars and the post-run checklist: `/run-benchmark`.
