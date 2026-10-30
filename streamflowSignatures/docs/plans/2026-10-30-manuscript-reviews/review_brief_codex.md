# Codex brief — adversarial read of the manuscript claim ledger (2026-10-30)

Context: the HISSS data paper (Scientific Data) and the co-authors' signature guidelines doc are
being checked for methodological accuracy and for every "number of X" claim against the data
files as deposited on HydroShare. A census script measured the files; a ledger assigns every
claim a status. Your job is to break the ledger: wrong statuses, claims the ledger missed,
measurements that do not support the stated verdict, and suggested replacements that are
themselves wrong.

Read, in this order:
1. docs/plans/2026-10-30-manuscript-claim-ledger.md — the ledger (statuses ✔ ✘ ≈ ⚠ —).
2. docs/plans/2026-10-30-manuscript-claim-ledger.csv — the 201 measurements the ledger cites
   by id in brackets.
3. docs/benchmarks/manuscript_claim_census.py — the script that produced them (check that each
   measurement actually computes what its label says).
4. docs/MANUSCRIPT_DRAFT.md lines 31–350 — the manuscript text under review (snapshot of the
   live Google Doc, 2026-10-30).
5. docs/SIGNATURE_GUIDELINES.md — the guidelines doc (Parts 1–5).
6. For method claims: config/signatures_config.json, julia/src/io.jl, julia/src/signatures.jl,
   docs/benchmarks/run_julia_benchmark.jl (gage-inclusion rule, l. 250–275),
   docs/benchmarks/climate_frame.jl, R/usgs_qualifiers.R, EO_data_processing/README_DAYMET.md.

Questions to answer, with file:line evidence:
- Is any ✔ wrong? Is any ✘ or ≈ wrong, i.e. the manuscript was right and the ledger misread the
  files or the code?
- Row H4: the ledger says the 60 % qualification denominator is "window start to the gage's last
  water year in the window". Confirm or refute from run_julia_benchmark.jl and the products'
  num_water_years / start_water_year / end_water_year columns (the CSV ids wy1980_2025.min_fraction).
- Row H5: the ledger says the seasonal 80 % rule masked no valid water year in either product and
  that season_excluded_years_* counts rejected years too. Confirm or refute from io.jl.
- Row M2: is "storage" structurally NA for un-normalized gages or not?
- Which numeric or method claims in the manuscript lines 31–350 have NO ledger row?
- Are the suggested replacement wordings (the ✘ rows) correct and complete?
- Anything in the census script that would silently mis-measure (id-key handling, NaN vs null,
  column selection, date handling)?

Do not edit any file. Output: a ranked findings table (BLOCKER / MAJOR / MINOR), each with the
ledger row id, the evidence, and a one-line fix; then a one-paragraph verdict (GO / GO WITH
CHANGES / NO-GO) on using the ledger as the basis of the co-author relay note.
