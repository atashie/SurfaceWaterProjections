# Streamflow Signatures — HISSS

@docs/STATUS.md

## What this is

Extraction of 100 annual hydrological signatures (+ 21 per-gage scalars → a 1,653-column
summary CSV and an annual-values parquet) from daily streamflow for ~8,000 US and Canadian
gages, for the HISSS data paper (Scientific Data, submission target 2026-11-09). Three
implementations: **Julia is canonical** (`julia/src/`); Python
(`python/streamflow_signatures/`) and R (`rpkg/`) are ports validated at full scale. The R
code at the repo root (`R/`, `run_*.R`, `config.R`) is the still-active raw-data INGESTION
path; `R/helperFunctions.R` is a deprecated shim. `EO_data_processing/` (Python) builds the
per-watershed MODIS and NLCD products and has its own CLAUDE.md.

## Ground truth and change flow

1. **Methodology** = the co-authors' guidelines Google Doc, snapshot
   `docs/SIGNATURE_GUIDELINES.md`. Code implements it; a disagreement is surfaced to the
   user, never resolved silently.
2. **Manuscript** (`docs/MANUSCRIPT_DRAFT.md`, read-only snapshot) must stay consistent
   with the code and the repo docs; corrections to it are relayed to the user, since the
   Google Doc cannot be edited from here.
3. **Julia first**, benchmark, then port to Python and rpkg. A change is done only when all
   three agree.
4. **At the start of every session run `/sync-docs`** — it fetches and diffs both Google
   Docs and is nearly free when nothing changed. Report "guidelines unchanged, manuscript
   unchanged" in one line when so.

## Critical constraints (always)

- **CSV output contract**: column names and order never change; every annual signature
  has exactly 16 columns (8 statistics + 8 Pettitt fields) from `generate_stats()`. The
  single-valued exceptions are enumerated in `docs/SIGNATURES.md` and `config.R` only.
- **Water year** Oct 1 – Sep 30. **Flow units** mm/day; 37 Canadian gages are raw m³/s
  (`area_normalized = FALSE`) and their Q-to-PPT signatures are NA by design.
- **Qualification**: 20+ water years per gage and ≥ 60 % of the window's years. Year
  rejection happens ONLY in `preprocess_daily_data()` (> 30 raw NAs or a gap > 3 days;
  negative Q only if configured). No NA filling and no per-year thresholds inside
  signature functions.
- **Config** `config/signatures_config.json` is the source of truth, with byte-identical
  bundled copies in `python/streamflow_signatures/data/` and `rpkg/inst/config/`. Julia
  bakes it at PRECOMPILE time — purge the compiled cache after a config edit.
- **Every artifact of a run lives in that run's own folder**; `docs/benchmarks/` holds
  tools and the long-lived reference CSVs only.
- **Delivered products are never rewritten** without an explicit user decision; STATUS.md
  lists the standing ones.
- **Inputs on the exFAT drive can be silently truncated** — verify sizes and `PAR1`
  footers before a long run, and use the rebuilt Daymet parquet.
- **Public mirror**: https://github.com/CZ-Sync/HISSS is a snapshot built by
  `./publish_to_hisss.sh` — run it after every merge to main; never commit there directly.
  Every external-facing repo pointer uses that URL.

## Signature statistics rule

Suffixes `_senn_slp, _linear_slp, _spearman_rho, _spearman_pval, _mk_rho, _mk_pval, _mean,
_median` plus `_pettitt_{cp_year, pval, pre_mean, post_mean, delta_mean, pct_change,
pre_mk_pval, post_mk_pval}`. Trend statistics require 60 % overall and 80 %
first-and-last-decade completeness and ≥ 20 annual values (recession and elasticity are
exempt). Eight signature CATEGORIES (`docs/signature_categories.csv`) group the 15
computing FUNCTION FAMILIES; "category" is the manuscript's and the guidelines' word.

## Procedures and scoped rules

- Skills (body loads when invoked): `/sync-docs`, `/add-signature`, `/run-benchmark`,
  `/cross-language-alignment`.
- Rules in `.claude/rules/` (load when matching files are touched): `signatures-code.md`
  for the Julia/Python/rpkg sources and config; `benchmarks.md` for `docs/benchmarks/`;
  `changelog.md` for CHANGELOG, STATUS and the reconciliation logs.
- `claude-skill/streamflow-signatures.md` is the USER-facing interpretation skill, not a
  Claude Code skill — keep it current when methodology, formats or validation change.
- Record user DECISIONS as such, with the date, in CHANGELOG and as a STATUS.md one-liner.

## Reference docs — read on demand, never `@`-import

| Read when you need | File |
|---|---|
| a signature's definition, method, units, caveats | `docs/SIGNATURES.md` (per family), then `docs/SIGNATURE_GUIDELINES.md` (ground truth) |
| architecture, data flow, NA pipeline, parquet inventory, HydroATLAS metadata, explorer builds, benchmark history | `docs/DEVELOPMENT.md` |
| what changed and why; open items in full | `CHANGELOG.md`; history in `changelog-old.md` and `docs/CHANGELOG_ARCHIVE.md` |
| cross-language residuals and gate results | `docs/CROSS_LANGUAGE_STATUS.md` |
| the 11 external data sources | `docs/DATA_SOURCES.md` |
| open guidelines / manuscript items | `docs/reconciliation/guidelines_todos.md`, `manuscript_log.md` |
| design records, run plans, manuscript edit lists | `docs/plans/` (excluded from the mirror) |
| EO products (MODIS LAI/LULC, Annual NLCD) | `EO_data_processing/README.md`, `README_NLCD.md` |

## Tests

`julia --project=julia julia/test/runtests.jl`; `pytest python/tests`; rpkg `testthat`
against the installed package (`R CMD INSTALL rpkg` first). A green unit suite does not
prove a production run — the orchestrator swallows per-family exceptions into missing
columns; see `/add-signature` §3.
