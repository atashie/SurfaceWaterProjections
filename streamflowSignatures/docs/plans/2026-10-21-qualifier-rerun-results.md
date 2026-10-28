# Qualifier-mask + climate-frame rerun — fixes, final artefacts and results (2026-10-28)

The consolidated record of the rerun that produced the **version-3 products** (public name "version 2").
Design and decisions: `2026-10-21-qualifier-mask-fix-rerun-plan.md` (v6; §10 the stage design, §10.5 the
census gates, §10.6 decisions 12–19); its §11 is the dated progress log with every run, defect and review
disposition. Reviews verbatim: `2026-10-21-qualifier-reviews/` (index there). Run folder:
`~/HISSS_data/qualifier-rerun-26oct2026/` (its `README.md` maps every subfolder, log and script).

## 1. What was wrong and what was fixed

| | Defect (every product before version 3) | Fix |
|---|---|---|
| HIGH | The February 2026 compile kept a USGS value only for the qualification strings `A`, `A e`, `P`, `P e`, so every approved day carrying a remark code (`A R` revised, `A <` / `A >` censored, `A [4]`) was NA: 193,542 days on 987 compiled gages; 1,177 / 950 water years of products #1 / #2 rejected by the per-year screening for that gap alone. | Rule: `R/usgs_qualifiers.R` (status letter A or P + finite value; remark codes never mask). Input: the 193,542 days re-fetched from NWIS (26 Oct 2026) and spliced into the compile — `input/combined_streamflow_data_09feb2026_qualifier_patched.parquet` + `correction_ledger.csv` + raw responses + audit B1–B8 PASS; metadata year counts recomputed in place. |
| MEDIUM | The three runners joined climate onto the streamflow dates, so a day absent from the streamflow record lost its Daymet value and PPT / TMIN / TMAX / SWE were interpolated across it. | `docs/benchmarks/climate_frame.{jl,py,R}` under `STREAMFLOW_CLIMATE_FRAME = left | outer`; the products run `outer` (Daymet's own calendar inside the record span and the gage's water years). The library is untouched. |

User decisions that shape the result: the gage sets are **frozen** at the version-2 cohorts (6,678 / 6,250;
no new watershed, 2026-10-23) — 31 / 29 gages that would qualify under the corrected record are listed, not
added; the two fixes go into one rerun (2026-10-26); internal "v3", public "version 2".

## 2. Final code

- Products built from the tagged baseline `qualifier-rerun-baseline-2026-10-26` (commit `54f0483`) in a clean
  worktree; the batch runner `docs/benchmarks/run_batched_julia.py` (`--gage-list`, `--climate-frame`,
  `--qa-schema`, `--stage`, `--require-clean-tree`) with the product wrappers
  `run_julia_benchmark_prod_{1993,1980}_2025_60pct_climate.jl`.
- Verification tools (main tree, post-tag patches recorded in plan §11): `splice_nwis_qualifiers.py` (fetch,
  splice, control, audit), `qualifier_year_oracle.py`, `climate_frame_inventory.py`,
  `preprocess_state_ledger.py` (daily delta ledgers + state tables), `assess_rerun_changes.py` (the change
  census, gates 5a–5g), `build_stage_cohorts.py`, `stage_attribution.py`, `check_annual_parquet_equality.py`
  and `check_summary_equality.py` (the E7 port gates), `validate_production_run.py`,
  `validate_annual_values.py`, `audit_qualification.jl`, `build_signature_explorer.py`,
  `check_html_basemap.py`, `stage_hydroshare_qualifier_v3.py` (+ `reconstruct_qualifier_patched_streamflow.py`).
- Driver scripts of the execution: run folder `scripts/` (`phase_d1*.sh`, `phase_d2*.sh`, `phase_e6.sh`,
  `phase_f.sh`, `merge_rpkg_batches.py`, the restage chains).

## 3. Final outputs

| Artefact | Where |
|---|---|
| Product #1 v3 — WY 1993–2025, 6,678 gages, 2,277 columns, 27,766,578 annual rows | `processedOuts_1993_2025_climate_v3_26oct2026/` (summary CSV, annual parquet, timing / merge JSON, explorer + `_annual/`, validation logs, v3-vs-v2 dashboard) |
| Product #2 v3 — WY 1980–2025, 6,250 gages, 2,277 columns, 34,813,881 annual rows | `processedOuts_1980_2025_climate_v3_26oct2026/` (same contents) |
| Corrected input + provenance | `input/` (patched parquet + metadata, ledger, `nwis_cache/`, control report, `oracle/`, audit), `manifests/` (selection, affected, expected-identical, stage cohorts) |
| rpkg reference on #1 v3 | `rpkg_reference_1993_2025_28oct2026/` (gates in `gate_*.txt`) |
| HydroShare staging (local, public name "version 2", nothing uploaded) | `~/Downloads/Signatures_v3_28oct2026/` R1 / R2 / R3 + `hisss_staging_manifest_*.csv` |
| Verification record per window | `census_<WU>/` (change census, gates JSON, classes), `e7_<WU>/` (port gates), `attribution_<WU>/`, `bprime_<WU>/` (delta ledgers, state tables), `S0_<WU>/`, `S1_<WU>/`, `S1prime_<WU>/` |
| Internal data review (user only, not for publication) | `review_<WU>/` (min / max / mean summaries) + the dashboards |

Superseded or failed artefacts kept for the record: `census_1993_2025_run1_bad_inventory/`,
`review_1993_2025_run1_bad_inventory/`, `bprime/inserted_rows_all_pandas_roundtrip_bad.csv`,
`bprime_<WU>/*_run1_bad_inventory.csv`, `logs_interrupted/`, `input/nwis_cache_to_2025-09-30/` (an earlier
request window), `preview_1993_2025/` (the pre-checkpoint preview).

## 4. Results (both windows, #1 / #2)

| Check | Result |
|---|---|
| S0 (February input, `left`) = version 2 | byte-identical summary CSVs, annual and summary gates bitwise — the tree reproduces its baseline (decision 18) |
| E1 production gates | 6 / 6 PASS |
| E2 annual values ↔ summary | PASS (`--floor 20`) |
| E3 swallowed failures | none in 16 batch logs |
| E4 qualification audit | emitted = selected ⊆ qualified; would-qualify kept out 31 / 29 |
| E6 explorers | built; basemap check PASS |
| B′ daily delta ledgers | COMPLETE: 2,702,639 / 2,217,736 changed daily values, every one attributed (restored day, inserted climate row, eligibility, secondary interpolation) |
| E8 change census vs version 2 | verdict **PASS**: 0 unexplained of 387,055 / 360,065 non-equal annual cells; 5,069 / 4,894 untouched gages identical at both layers; precipitation and GDD totals = daily sums; gained years = the oracle's, gage for gage; 0 years lost |
| Restored | 144,724 / 119,933 days inside the window on 899 / 757 gages; 1,177 / 950 water years valid again on 619 / 502 gages |
| E7 Python port on the 1,609 / 1,356 changed gages | annual key sets identical; 117 / 188 values beyond 1e-9, all threshold ties or float noise in 19 named signatures; summary 0 gated mismatches, 103 finite cells in the named waived columns and, separately, 8 / 7 gate-state cells accepted after recomputation reproduced both sides |
| Stage attribution | S1 = S0 (200 control gages), S1′ = S2 (197 / 199 join-only), S2 = S0 (class U) — all PASS at both layers; the interaction residual is reported, not gated (the elasticity family, non-additive by construction) |
| F rpkg on #1 v3 | 0 errors; schema and failure gates PASS; annual 1e-6 PASS with the named residuals (1 signed-zero key, 323 threshold ties); climate family bit-identical; R² 0.999922 |
| G reviews (2026-10-28) | Sonnet ×2 + Codex, blind: GO WITH CHANGES, no blocker; headline numbers reproduced independently; every finding applied |

PASS means agreement under the named allowances listed in the gate files, not unconditional equality.

## 5. Defects found during execution (none in the products; all in verification tooling, fixed post-tag)

1. The combined inserted-row inventory had been round-tripped through pandas' default CSV float parser
   (1–2 ulp off on 38 % of rows) → rebuilt by text concatenation; the reconstruction recipe reads the ledger
   with `float_precision="round_trip"` for the same reason.
2. The daily-diff tool expected an interpolate on a climate NA day of a year the variable is invalid in (the
   preprocessor leaves it NaN).
3. The census failed gate 5f and the `total_possible` freeze on the gages the replay does not cover, and left
   116 record-dependent keys emitted on old years unexplained → Qann fallback, scoped freeze, explanation
   branch (validated by the E7 key-set identity).
4. The E7 climate subset used the streamflow column names; the Python runner silently runs without climate
   when the file is absent → fixed with an assertion.
5. E2 launched without `--floor 20`; pandas 3 NaN ≠ NaN in the attribution's string comparison.
6. From the reviews: the attribution tool's class-U identity was documented but not run; interaction gages
   absent from S1 were silently intersected (now a coverage gate + explicit "S1 inferred = S0" labels); absent
   keys were conflated with present-but-NaN values; staged dictionaries carried stale rows; a stale manifest.

## 6. Open (user)

Done 2026-10-28: commit `161c56e` pushed, tag pushed, mirror published; DECISION (user): the
climate-frame default is `outer` (helpers, batch-runner record, wrappers, self-tests, docs updated
the same day — CHANGELOG [October 2026], Phase H). Open: the HydroShare deposit of the staged
folders (the user uploads, planned 2026-10-29); the co-author relays (`docs/reconciliation/`).
