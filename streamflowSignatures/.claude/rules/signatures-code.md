---
paths:
  - "julia/src/**"
  - "julia/test/**"
  - "python/streamflow_signatures/**"
  - "python/tests/**"
  - "rpkg/R/**"
  - "rpkg/tests/**"
  - "config/signatures_config.json"
---

# Signature code — Julia canonical, Python and rpkg ports

- **Julia first, then port.** `julia/src/` is canonical. Python
  (`python/streamflow_signatures/`) and rpkg (`rpkg/R/`) mirror it module for module and
  are validated at full scale (1,653 columns, 6,678 gages — `docs/CROSS_LANGUAGE_STATUS.md`).
  A change is not done until all three agree; use `/cross-language-alignment` for
  divergences. `R/helperFunctions.R` is the DEPRECATED legacy shim — never extend it.
- **Method definitions come from `docs/SIGNATURE_GUIDELINES.md`** (synced Google Doc,
  ground truth); per-family reference in `docs/SIGNATURES.md`. If code and guidelines
  disagree, say so — the user decides which side moves.
- **Output contract.** Summary-CSV column names and order are frozen; every annual
  signature yields exactly 16 columns via `generate_stats()` (8 statistics + 8 Pettitt
  fields); per-gage scalars are the documented exceptions only. Never rename a shipped
  column (`recession_alpha_point_cloud_linear_reservoir`, `Q95_Q10`, ...).
- **NA handling is central.** `preprocess_daily_data()` runs once per gage before any
  signature: interpolates internal gaps ≤ 3 days, rejects years with > 30 raw NAs or a
  gap > 3 days, negative-Q rejection only if `reject_negative_flow`. Inside signature
  functions: NO `fillna(0)`, NO per-year `min_days` / `max_na_frac` checks, NO
  `min_Q_value_and_days` filter, NO implicit SWE use (snow runs only on an explicit
  `snow_data` frame). Constant-SD is a QA flag, never a rejection.
- **Gates live in the orchestrator, not in signature functions**: trend completeness
  (60 % overall, 80 % in the first and last decade; recession and elasticity exempt),
  the 20-value stats floor (same exemptions), the record-anchored snow decade gate, and
  the `area_normalized` gate (runoff ratios, elasticity, Q-P seasonality and storage are
  skipped when false).
- **Annual values.** Thread the `collector` kwarg into every `generate_stats()` call;
  with no collector the behaviour must stay byte-identical.
- **Config.** `config/signatures_config.json` is the source of truth; keep the bundled
  copies (`python/streamflow_signatures/data/`, `rpkg/inst/config/`) byte-identical.
  Julia bakes `CFG_*` at precompile time — purge the compiled cache after a config edit
  (`/run-benchmark`). rpkg needs `R CMD INSTALL rpkg` before tests or benchmarks see a
  change.
- **Numerics that bit before:** summation order (R `mean()` vs a sequential sum flipped
  whole drought plateaus through the strict `<`); Julia `unique` uses `isequal`, so
  −0.0 ≠ +0.0 (agreed fix: Julia moves to `==`); mid-event day-of-water-year off by one
  for even-length events; guards written `== 0` instead of `<= 0`; Mann-Kendall p-value
  is the continuity-corrected normal approximation with tie-corrected variance in ALL
  languages (scipy's default differs); constant series → NaN tau and p.
- **Tests.** Julia `julia --project=julia julia/test/runtests.jl`; Python `pytest
  python/tests`; rpkg `testthat` against the INSTALLED package. Cross-language fixtures
  agree to ≤ 4e-14. A green unit suite does not prove a production run — the
  orchestrator's per-family `try/catch` hides failures as missing columns
  (`/add-signature` §3).
