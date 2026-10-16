# Pacific-coast precipitation review — findings and the streamflow-frame defect

**Date**: 2026-10-16 · **Status**: review complete; defect DOCUMENTED, fix ON HOLD
(user decision 2026-10-16: discuss with the research team first).
**Trigger**: the user found basin annual precipitation of ~4 m on the Pacific coast
implausible and asked for (1) a comparison with public data, (2) a review of the Daymet
aggregation, (3) a review of the annual totals computed from the basin series.
**Product reviewed**: WY 1993–2025 version-2 (`processedOuts_1993_2025_climate_12oct2026`,
staged copy `~/Downloads/Signatures_v2_12oct2026/resource1_signatures_wy1993_2025/`),
column `prcp_annual_mean` and the `prcp_annual` rows of the annual parquet; input
`daymet_1980_2025_29sep2026.parquet`.
**Evidence folder**: `~/HISSS_data/climate-rerun-12oct2026/precip_review_16oct2026/`
(scripts `pac_prcp.py`, `indep_compare.py`, `wc_vs_prism.py`, `daymet_check.py`,
`annual_vs_daymet.py`, `full_domain_check.py`, `gap_detail.py`; outputs
`prcp_subset_{1993,1980}.csv`, `compare_prism.csv`, `compare_worldclim.csv`,
`full_domain_prcp_check.parquet` — one row per product gage-year: product value, Daymet
day count, Daymet WY sum, expected value, difference). Python:
`~/HISSS_data/venvs/daymet/bin/python -I` (duckdb, pandas, pyarrow).

## 1. The product values are plausible

Only a handful of gages reach 4 m; the wettest are headwater basins behind famously wet
stations. WY 1993–2025 product, `prcp_annual_mean` (mm/yr), 6,678 gages: median 1,117,
p95 1,764, p99 2,948, max 4,462; 58 gages ≥ 3,000 mm, 3 gages ≥ 4,000 mm (WY 1980–2025
product: 52 and 1; max 4,555).

| Region | Gages | Median | p90 | Max |
|---|---|---|---|---|
| WA west of 121°W | 138 | 2,323 | 3,139 | 3,816 |
| OR west of 121°W | 153 | 1,643 | 2,466 | 3,350 |
| CA west of 120°W | 248 | 1,145 | 1,703 | 2,190 |
| BC and SE Alaska | 286 | 995 | 3,159 | 4,462 |

Wettest basins: 08MH141 Coquitlam R above Coquitlam Lake 4,462; 15088000 Sawmill Cr nr
Sitka 4,239; 15052000 Lemon Cr nr Juneau 4,122; 15101490 Greens Cr 3,994; 08HE006/08HE008
Zeballos R 3,874/3,860; 08HF005 Nimpkish 3,824; 12137290 SF Sultan 3,816; 12147600 SF Tolt
3,798; 08MH147 Stave 3,538; 08HB086 Tofino Cr 3,534; 08GA030 Seymour 3,450; 14138900 NF
Bull Run 3,350; 12036000/12035400 Wynoochee 3,243/3,219; 12040500 Queets 3,178; 12039500
Humptulips 3,138; 12041200 Hoh 2,947.

Seasonal totals sum to the annual total exactly (max |diff| 1.8e-12 mm). Runoff ratio
`Qann_mean / prcp_annual_mean`: median 0.34; 100 gages have Q > P (43 BC/AK, 16 WA, 4 OR,
2 CA, 35 elsewhere), e.g. Lemon Cr Q 5,531 vs P 4,122; Klaskish 08HE009 Q 6,288 vs P 3,278;
Tofino Cr Q 5,652 vs P 3,534 — Daymet UNDER-estimates windward coastal mountains (gauge-
based grid; plus glacier melt), it does not over-estimate.

**Why it looked like "often 4 m"**: the explorer's "Robust color scale" clips at the
domain-wide p5–p95 (`build_signature_explorer.py`), i.e. 1,764 mm, so every coastal basin
above that is painted the top colour; hover values are unclipped.

## 2. Independent comparisons

**GAGES-II PRISM 1971–2000 basin mean** (`PPTAVG_BASIN`, cm;
`~/HISSS_data/daymet-replay-01oct2026/inputs/gagesMetadata/conterm_climate.txt`;
n = 4,935 CONUS gages): median Daymet(WY 1993–2025) / PRISM ratio WA 1.01 (p10–p90
0.85–1.14), OR 1.03 (0.95–1.15), CA 0.99 (0.85–1.10), rest of CONUS 1.07 (0.94–1.14). The
wettest Olympic basins are well BELOW PRISM: Wynoochee nr Grisdale 0.59 (PRISM 5,480 mm),
Wynoochee above Save Cr 0.66 (4,917), Hoh 0.73, Humptulips 0.75, Queets 0.87 (3,635);
Sultan/Tolt 0.93–1.05 (PRISM SF Sultan 4,109; SF Tolt 3,819).

**HydroATLAS WorldClim 1970–2000** (`pre_mm_uyr`,
`resource4_geometry_attributes/hisss_watershed_hydroatlas_attributes.parquet`): Daymet /
WorldClim median 1.13–1.16 in CONUS, 1.53 in BC/AK (up to 3.1 for Coquitlam). But
WorldClim / PRISM on the same CONUS coastal gages is 0.84–0.93 — WorldClim is the low
product in the mountains, not evidence against Daymet.

**Station normals and basin estimates** (sources in §6):

| Reference (period) | Value (mm) | Nearest product basin (mm) |
|---|---|---|
| Forks WA (1991–2020) | 3,046 | Queets 3,178; Hoh 2,947 |
| Hoh Rain Forest station | 3,298 | Hoh 2,947 |
| Quinault ~3,485; USGS HCDN basin estimate, Quinault R at Quinault Lake | 4,724 | Humptulips 3,138; Wynoochee 3,219 / 3,243 |
| Bull Run watershed (Portland Water Bureau; 80–170 in range) | ≈ 3,300 | NF Bull Run 3,350 |
| Tofino airport (1981–2010) | 3,271 | Tofino Cr 3,534 |
| Coquitlam Lake station (1951–1981 normals) | 3,616 | Coquitlam above lake 4,462 |
| Metro Vancouver watersheds | ~3.5 m rain + 4.5 m snowpack | Seymour 3,450 |
| Prince Rupert | 2,619 | — |
| Ketchikan (1991–2020) | 3,798 | 15072000 3,358; 15085100 3,196 |
| Yakutat (1991–2020) | 3,564 | 15129600 3,191; 15129500 3,026 |
| Sitka airport / Hidden Falls (1991–2020) | 2,146 / 3,253 | Sawmill Cr 4,239 |
| Juneau airport / downtown (1991–2020) | 1,702 / 2,276 | Lemon Cr 4,122 (glacierized; Q 5,531) |
| Little Port Walter; Hucuktlis (Henderson) Lake | ~5,700; 6,903–7,296 | no gage |

Henn et al. (2018, J. Hydrology 556): gridded products differ by > 200 mm/yr at high
elevation (5–60 % relative), true biases likely larger, gauge-based grids underestimate
in complex terrain. **Conclusion**: 3–4.5 m basin totals on the windward Olympics,
Cascades, BC Coast Mountains and SE Alaska are consistent with or below the references.

## 3. Daymet aggregation (`EO_data_processing/daymet/`) — no defect

`daymet_weights.py`: exactextract coverage fraction × true ellipsoidal cell area on the
1 km LCC grid; checks ≥ 1 cell per polygon, coverage in (0, 1], representative point inside
own cells, `area_ratio` reported. `daymet_aggregate.py`: basin-day MEAN = Σ w·v / Σ w, fill
cells excluded and weights renormalised per day, fatal check when valid weight exceeds the
basin total (double counting), range check, deterministic; units stay mm/day; 365-day
calendar (Daymet drops Dec 31 in leap years). Validated (`README_DAYMET.md`): matches
exactextract's own weighted mean to ≤ 5e-11; co-authors' gdptools product on 5,965 basins
R² ≥ 0.9999999 daily, prcp annual totals within ±0.0011 % (p01–p99), largest daily
difference 0.015 mm; ORNL single-pixel API checks. That validation proves the aggregation,
not Daymet's accuracy.

## 4. Annual totals in Julia (`julia/src/climate.jl`) — correct as designed

`calculate_precipitation_metrics` requires a complete 365/366-day grid
(`_complete_year_grid`) and sums with `_seqsum`; seasons via `_grid_months`. Full-domain
check (`full_domain_check.py`): of 200,832 product gage-years, 199,913 equal the Daymet
water-year sum to |diff| ≤ 1e-6 mm, with the documented Dec 31 fill — Daymet drops Dec 31
of leap calendar years, so the following water year has 364 Daymet days and
`preprocess_daily_data` (`julia/src/io.jl`, full date grid → `leftjoin(grid, yr_df)` →
PPT NaN ≤ 30 → linear interpolation of internal gaps ≤ 3 days) fills it with
(Dec 30 + Jan 1)/2 (`docs/SIGNATURES.md` §14 caveats: mean 0.25 %, max 2.1 % of a year).

## 5. DEFECT — climate rides on the streamflow frame (MEDIUM; fix on hold)

**Mechanism.** The three benchmark runners join each gage's Daymet series onto that gage's
streamflow ROWS before the library builds the date grid:

| Runner | Line | Statement |
|---|---|---|
| `docs/benchmarks/run_julia_benchmark.jl` | 220, 332 | `leftjoin(gage_df, gage_climate, on=[:gage_id, :date])` |
| `docs/benchmarks/run_python_benchmark.py` | 319 | `gage_df.merge(cslice, on=["gage_id", "date"], how="left")` |
| `docs/benchmarks/run_rpkg_benchmark.R` | 276 | `merge(gage_data, climate_merge, by = "date", all.x = TRUE)` |

(`run_batched_julia.py` docstring line 7 describes the same.) A date with NO streamflow
row (as opposed to a row with NA Q) therefore reaches `preprocess_daily_data()` without
PPT/TMIN/TMAX/SWE, the grid left-join makes it NaN in every series, and step 4 linearly
interpolates it — climate that exists in the input is replaced by an interpolate. Runs
> 3 days reject the year for Q anyway, so only 1–3-day holes survive, and gages drop out
preferentially during storms: the interpolate is biased low. The library is correct for
the frame it is given; the defect is in the runners' frame construction.

**Evidence** (WY 1993–2025 version-2 product vs Daymet input; `full_domain_check.py`,
`gap_detail.py`):
- 919 of 200,832 `prcp_annual` gage-years (0.46 %, 679 gages) deviate from the true sum.
- Every deviating gage-year has ≥ 1 date absent from the streamflow record (0 counter-
  examples). 1,172 gage-years have absent dates; 919 of them deviate (the rest match by
  chance: dry or linear days).
- |diff| median 6 mm, p90 54, p99 327, max −606 mm.
- Worst: 02146285 (Charlotte, NC) WY 2024 — 22 absent dates in 1–2-day runs carrying
  713 mm (32.4 mm/day vs 2.4 mm/day on present days); product 942.9 mm, Daymet 1,548.9 mm.
  Its 22-year `prcp_annual_mean` is 892 vs 1,182 mm (−24.5 %).
- Per-gage bias of `prcp_annual_mean` > 1 % at 8 gages: 02146285 −24.5 %, 0214297160
  −6.1 %, 0214269560 −3.4 %, 02329342 −3.3 %, 02336030 −2.4 %, 11363930 −1.9 %, 11153650
  −1.8 %, 03250190 −1.1 %.
- Pacific coast (west of 120°W): 70 affected gages, max |bias| 1.9 % — NOT the cause of
  the coastal totals.

**Reach (by code reading, not measured).** The same frame (`cdata`) feeds the Q-to-PPT
family (`analyze_Q_PPT_relationships`: runoff ratios, elasticity, Q-P seasonality,
avg_storage), the snow family (`swe_max_to_ppt`; SWE itself), and the temperature metrics,
in every product: the 2026 deliveries (#1, #2), both version-2 reruns and the rpkg/Python
reference runs (parity holds because all three runners do the same). Long-standing; not
introduced by the climate family.

**Proposed fix** (when approved):
1. Runners, Julia first: outer-join the climate on `date` (Julia `outerjoin`, pandas
   `how="outer"`, R `all = TRUE`), or hand the climate to the preprocessor as its own
   frame. Added rows have NA Q — the grid normalisation already yields NaN for absent
   dates, so the Q diagnostics, year rejection and every Q-only signature must be
   bit-identical. Restrict added rows to the streamflow record span if memory matters
   (years with Q all-NaN are rejected anyway; decision 13 of the climate plan keeps
   climate years inside the streamflow-valid years).
2. Gate: `docs/benchmarks/check_additivity.jl NEW.csv V2.csv --expect-added 0` — zero
   change outside the PPT/TMIN/TMAX/SWE-derived families; changes confined to the 1,172
   gage-years with absent dates (list from `full_domain_prcp_check.parquet`, extended to
   the WY 1980–2025 window); regression test with a fixture that drops streamflow rows on
   wet days.
3. Ports (Python, rpkg) → `/cross-language-alignment`; then rerun both version-2 products
   (batched Julia) before the HydroShare deposit, regenerate explorers/dashboards, update
   the R1/R2 READMEs; CHANGELOG Fixed (MEDIUM) and move the Known Issue to resolved.

**Decision (user, 2026-10-16)**: document for future work; do not implement until the
research team has discussed it. Relay item: `docs/plans/2026-10-12-coauthor-relay.md` →
"For discussion with the research team".

## 6. Sources

Forks, WA — https://en.wikipedia.org/wiki/Forks,_Washington · Hoh Rainforest —
https://en.wikipedia.org/wiki/Hoh_Rainforest · Quinault, WA —
https://en.wikipedia.org/wiki/Quinault,_Washington · USGS HCDN 12039300 Quinault R at
Quinault Lake — https://pubs.usgs.gov/wri/wri934076/stations/12039300.html · Portland Water
Bureau, Bull Run — https://www.portland.gov/water/about-portlands-water-system/about-bull-run
· Tofino — https://en.wikipedia.org/wiki/Tofino · Coquitlam Lake —
https://en.wikipedia.org/wiki/Coquitlam_Lake · Metro Vancouver watersheds —
https://en.wikipedia.org/wiki/Metro_Vancouver_watersheds · Prince Rupert —
https://en.wikipedia.org/wiki/Prince_Rupert,_British_Columbia · Ketchikan —
https://en.wikipedia.org/wiki/Ketchikan,_Alaska · Yakutat —
https://en.wikipedia.org/wiki/Yakutat,_Alaska · Sitka — https://en.wikipedia.org/wiki/Sitka,_Alaska
· Juneau — https://en.wikipedia.org/wiki/Juneau,_Alaska · Hucuktlis Lake —
https://en.wikipedia.org/wiki/Hucuktlis_Lake · Henderson Lake rainfall —
https://newwestrecord.ca/local-news/this-unique-bc-place-receives-the-most-rainfall-in-all-of-north-america-3094242
· Port Walter — https://franco.wiki/fr/Port_Walter.html · Henn et al. 2018 —
https://impacts.ucar.edu/en/publications/an-assessment-of-differences-in-gridded-precipitation-datasets-in/
· Adam, Lettenmaier & Wood (AMS) — https://ams.confex.com/ams/84Annual/webprogram/Paper72487.html
· UCAR Climate Data Guide, Daymet — https://climatedataguide.ucar.edu/climate-data/daymet-daily-surface-weather-and-climatological-summaries
