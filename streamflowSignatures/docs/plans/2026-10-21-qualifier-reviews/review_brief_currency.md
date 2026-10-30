# Review brief — code and documentation currency after the default flip (2026-10-28, blind, read-only)

Repository: /Users/arik/github/SurfaceWaterProjections/streamflowSignatures (main, HEAD = `713fdbf`, clean tree).
Commit range under review: `54f0483` (tag `qualifier-rerun-baseline-2026-10-26`) → `161c56e` (Phases D–H of the
qualifier-mask + climate-frame rerun) → `713fdbf` (the default flip). `git diff 54f0483..HEAD` is the change set.
Run folder: /Users/arik/HISSS_data/qualifier-rerun-26oct2026 (**R**). Staged deposit (local only, nothing uploaded):
/Users/arik/Downloads/Signatures_v3_28oct2026 (**D**). Public mirror: https://github.com/CZ-Sync/HISSS (a snapshot
of the main tree minus `docs/plans/` and the other exclusions in `publish_to_hisss.sh`, published from `713fdbf`).

Hard rules: do NOT edit any file of the repository, the run folder, the staged deposit or the worktree
`R/worktree/`; no state-changing git (no checkout / stash / commit / reset); no production job; scratch writes only
under the folder named in your assignment; DuckDB `memory_limit` ≤ 1 GB (RAM is shared); never read or print
credentials. Reading, grep, small Python / R / Julia checks, the unit suites and the self-tests are welcome.
Python with duckdb + pandas + pyarrow: /Users/arik/HISSS_data/venvs/daymet/bin/python; the repo `.venv/bin/python`
has pandas / pyarrow. Julia: `~/.juliaup/bin/julia --project=julia`.

## The claims under review

1. **USGS qualifier rule (HIGH defect, fixed).** The February 2026 compile kept a USGS value only for the
   qualification strings `A`, `A e`, `P`, `P e`; approved days with a remark code (`A R`, `A <`, `A >`, `A [4]`)
   were NA. The rule now lives in `R/usgs_qualifiers.R` (status letter A or P + finite value, sentinel −999999
   rejected; remark codes never mask, carried verbatim in `flag`), is applied by `generate_streamflow_dt` in
   `R/helperFunctions.R`, sourced by `run_ingest_usgs_hydat.R`, tested by `R/tests/test_qualifier_codes.R`.
   `config.R` no longer carries the stale `USGS_ACCEPTED_CODES`. The splice tool
   `docs/benchmarks/splice_nwis_qualifiers.py` implements the same rule.
2. **Climate frame (MEDIUM defect, fixed).** The three runners (`docs/benchmarks/run_julia_benchmark.jl`,
   `run_python_benchmark.py`, `run_r_benchmark.R` / `run_rpkg_benchmark.R`) join Daymet through
   `docs/benchmarks/climate_frame.{jl,py,R}` under `STREAMFLOW_CLIMATE_FRAME = left | outer`.
   **DECISION (user, 2026-10-28): the default is `outer`**; `left` reproduces every product before the rerun and
   the S0 / S1 attribution baselines. The three signature libraries (`julia/src/`, `python/streamflow_signatures/`,
   `rpkg/`) are untouched by the rerun.
3. **Wrappers.** Every Julia wrapper in `docs/benchmarks/` pins its frame with `get!` (ten pre-rerun wrappers
   `left`; the two `*_climate.jl` product wrappers `outer`, rewritten for the version-3 inputs, prefix
   `streamflow_<w>_60pct_climate_v3_26oct2026`, with a guard that refuses to write into the delivered product
   folder when `STREAMFLOW_OUTPUT_DIR` is unset). The batch runner `run_batched_julia.py` records the resolved
   default. Self-tests `selftest_climate_frame.{jl,py,R}` check the default (41 / 38 / 43 checks),
   `selftest_run_batched_julia.py` 28 checks.
4. **Products.** Version-3 products (public name "version 2") in R: #1 WY 1993–2025, 6,678 gages, 2,277 columns,
   27,766,578 annual rows; #2 WY 1980–2025, 6,250 gages, 34,813,881 rows. Gage sets frozen at the version-2
   cohorts (no new watershed; 31 / 29 would-qualify gages listed, not added). Restored: 193,542 days compiled
   (987 gages); 1,177 / 950 water years valid again. Nothing deposited to HydroShare yet.
5. **Docs updated** for 1–4: `README.md`, `docs/STATUS.md`, `CHANGELOG.md` ([October 2026] + Known Issues),
   `docs/DEVELOPMENT.md`, `docs/DATA_SOURCES.md`, `docs/CROSS_LANGUAGE_STATUS.md`, `docs/SIGNATURES.md`,
   `claude-skill/streamflow-signatures.md`, `python/README.md`, `.claude/rules/benchmarks.md`,
   `.claude/skills/run-benchmark/SKILL.md`, `docs/QUALIFIER_RERUN_RESULTS.md`, plan §11,
   `docs/reconciliation/{guidelines_todos,manuscript_log}.md`, the staged READMEs / run notes / dictionary in D.

## Questions (answer each with evidence; name file:line; recompute rather than trust)

A. **Stale text.** Anywhere in the repository (code comments, docstrings, help strings, READMEs, docs, skills,
   rules, the staged deposit text) that still states the OLD behaviour as current: the four-string qualifier rule
   (`A`, `A e`, `P`, `P e` / `keep_codes` / `USGS_ACCEPTED_CODES`), "climate joins on the streamflow dates" as the
   current join, `left` as the default frame, "left (unset)", the version-2 row counts (27,602,161 / 34,681,331)
   presented as current, "version 3" leaking into public-facing text, or a claim that something is deposited.
   Distinguish historical statements (allowed, when dated) from current-state statements (must be right).
B. **Code ↔ docs agreement.** Does the code do what the docs say, and vice versa: the default resolution in each
   helper and each runner; which wrapper pins which frame; the product-folder guard; what `run_batched_julia.py`
   records; what the ingestion actually keeps (trace the R call chain from `run_ingest_usgs_hydat.R`); whether any
   file still references the removed `USGS_ACCEPTED_CODES`; whether every documented CLI flag exists.
C. **Tests.** Run the relevant suites and self-tests and report pass / fail counts verbatim: Julia
   `julia --project=julia julia/test/runtests.jl`, `pytest python/tests`, the rpkg `testthat` suite (against the
   installed package; skip if not installed and say so), `R/tests/test_qualifier_codes.R`, the three
   `selftest_climate_frame.*`, `selftest_run_batched_julia.py` (daymet venv). Are the documented check counts
   (41 / 38 / 43 / 28 / 37) the real ones?
D. **Numbers.** Every count, date, commit hash, folder name and file name in the docs listed under claim 5 is
   consistent across documents and with the artefacts in R and D (sizes, row counts via parquet metadata, the
   manifests, `input/provenance.json`, `hisss_staging_manifest_*.csv`).
E. **Public mirror.** Does the mirror match HEAD minus the exclusions (clone it read-only to your scratch folder
   and diff)? Do public-facing documents link to paths excluded from the mirror (`docs/plans/...`) or to the
   private repository, so a public reader hits a dead link?
F. Verdict GO / GO WITH CHANGES / NO-GO with findings ranked BLOCKER / MAJOR / MINOR / NOTE, each with the
   file:line, the evidence and the one-line fix. Say explicitly which claims you verified and which you could not.
