# HydroShare deposit — Codex adversarial review, triage (2026-10-30)

**DECISION (user, 2026-10-30): none of these findings holds up the release of this version; all are noted for a
potential future release** (CHANGELOG [Unreleased] → Planned, with checkboxes). The deposit ships unchanged.

Five serial Codex runs (codex-cli 0.165.0, read-only, reasoning high), one per resource, 12:32–14:13 EDT, on the
byte-verified upload copies + the live page metadata of ~12:30 (abstracts since replaced by the user).
Reports: `codex_<resource>.md`; brief: `review_brief_shared.md`. Every item below marked VERIFIED was
re-checked by Claude (code read or recomputed); the rest are Codex's evidence, plausible, not re-run.

## A. Values in deposited files that a fix would change (user decision — delivered products)
| # | Item | Resources | Status |
|---|---|---|---|
| A1 | `Flow_Reversals_fall`: the water-year series is masked to Sep+Oct+Nov and passed to `count_flow_reversals` as one run, so Nov 30 → Sep 1 is treated as adjacent (`julia/src/pulses.jl` seasonal loop). Up to ±2 reversals per gage-year; other seasons are contiguous. | Signatures 1993-2025, 1980-2025 (+ explorers) | VERIFIED (code) |
| A2 | Mann–Kendall tie correction double-counts a tie group holding both −0.0 and +0.0 (`julia/src/stats.jl`: `unique` uses isequal, `count(==)` does not). 13 series (1993-2025) / 14 (1980-2025), almost all `elasticity_annual`; deposited `_mk_rho` reproduces the buggy ties exactly in every case; tau shifts ~0.0005–0.007; Codex: no p-value crosses 0.05. Same root cause as the STATUS storage −0.0 item. | Signatures both | VERIFIED (recomputed, `mk_signed_zero.py`) |
| A3 | HydroATLAS attributes of mis-snapped outlets: 36 gages ≤ 100 km² aggregated over > 100,000 km² (e.g. `07001910`, 3 km² → 1,828,066 km², Mississippi main stem); 166 rows > 100× and 2,020 rows `watershed_area_rel_diff` > 0.5 (median 61 km²: small basins below HydroBASINS L12 granularity). README documents the screening column but not the scale of the problem. | Watershed Geometry | VERIFIED (counts) |

## B. Deposited text that is wrong (fix = re-render + re-upload the file)
| # | Item | File | Status |
|---|---|---|---|
| B1 | Six attributes stored with HydroATLAS scale factors the dictionary omits: `lka_pc_use`, `dor_pc_pva` (% ×10), `pop_ct_usu` (thousands), `nli_ix_uav` (index ×100), `hft_ix_u93`, `hft_ix_u09` (index ×10). All other scaled fields document their factor. `gdp_ud_usu` = GDP total, `gdp_ud_sav` = GDP per capita (2011 international USD, PPP, 2015) — both described as "GDP (USD)". | HydroATLAS dictionary | VERIFIED (BasinATLAS catalog v10) |
| B2 | `outlet_ORDER` called "Strahler stream order"; values are HydroBASINS `ORDER_` (classical order, decreases downstream; 0 = coastal). | HydroATLAS dictionary | Codex evidence |
| B3 | Readme says all three explorers show the boundary layer's `low_confidence`; the NLCD explorer embeds the NLCD table's union flag (3 partial-coverage gages). | Vegetation readme | VERIFIED (builder) |
| B4 | `n_composites` "mean over the basin's tiles": the code takes the mean over ALL tiles of the month's run (`lai_pipeline.py:173`), one value per month for every basin (242/270 months single-valued). | Vegetation dictionary | VERIFIED (code) |
| B5 | "≥ 60 % of the window's years": the implemented denominator is window start → the gage's last in-window year (`run_julia_benchmark.jl` qualification block); 325 gages of 1980-2025 are below 60 % of the full 46-year window. Gage sets stay frozen — wording only. | Signatures readmes, input readme + dictionary | VERIFIED (code) |
| B6 | `_pettitt_pct_change` = 100 × delta / abs(pre_mean), NA if abs(pre_mean) ≤ 1e-10 — dictionary omits abs and the NA rule (sign differs for ~57k values). | Signatures dictionaries | VERIFIED (`changepoint.jl:348`) |
| B7 | Drought zero-threshold "reads 0 by construction" — false for 4 gages with negative flows. `avg_storage` = storage interpolated at mean Q, not the mean of cumsum(P−Q). `_pettitt_delta_mean` rows inherit bounded units (0–1, −1…+1, day of WY) that differences violate. Percentile-order flag written as strict `<`. `negative_ann_*` counts interpolated days. Recession seasonality minimum = continuous 365-day phase. | Signatures dictionaries | Codex evidence |
| B8 | `no_area_reference` = HydroBASINS fallback AND no published area (46 official outlines without area are not flagged); the 0.5 area and `coverage_sum < 4` thresholds are not stated. | Input dictionary | Codex evidence (matches the open STATUS follow-up) |
| B9 | Old expansion "Hydrologic Information Signatures and Summary Statistics" — one line in each of the five readmes (titles now "Hydroclimate Information, Signatures and Summary Statistics"). LAI explorer text "the older version of this idea". | 5 readmes, LAI explorer | VERIFIED |

## C. Completeness gaps (documentation to add)
20-value floor for ALL statistics (recession/elasticity exempt) undocumented (19,051 / 15,166 series have annual values but no summaries);
climate/snow signatures share the streamflow year mask and no Canadian gage has WY 2025 (HYDAT caveat absent from the signature readmes);
annual-Parquet columns and data types missing from the signature dictionaries; Daymet forcing-quality flags not referenced from the
signature readmes; agency flag-code legend (`Dis`, `Eqp`, `Fld`, `Ssn`, `Mnt`, `***`, `[4]`) absent; HydroATLAS class legends and the
boundary-layer columns lack a dictionary; 37 MODIS LULC class columns (LCCS1/2/3, prop) carry codes without names; LAI QA fractions'
zero-denominator NA rule; signature explorers pool the raw-m³/s gages into the colour scale and histogram.

## D. Repo documents out of step with the deposit (not deposited)
`EO_data_processing/geometry/README.md` (signature-bearing, zero-padding, 04HA001), `README_DAYMET.md` (53 large basins "omitted"),
`EO_data_processing/README.md` (~200 LULC class columns; deposit has 102), `R/aggregate_hydroatlas_metadata.R` unit placeholders.

## Verified clean (the reviews' positive results)
All populated means/medians reproduce from the annual values (895,157 / 842,497); all Pettitt splits, means and deltas reproduce; all
12 QA flags reproduce; selection lists, the 31/29 would-qualify rows and the qualifier ledger reproduce; Daymet complete (8,017 sites ×
16,790 dates, no nulls) and annual precipitation/GDD totals reproduce from daily Daymet to 1e-11; the outer climate frame confirmed;
every simplified polygon valid, areas reproduce in EPSG:6933; zero-strip joins leave no unexplained failure; explorer payloads match
the tables; LULC schemes and NLCD classes sum to 100; no "previous version" wording besides B9's phrase; the README's 17,166-id union
= 16,994 metadata ids + 172 HYDAT-interference-only stations (closes the "unverifiable" follow-up).
