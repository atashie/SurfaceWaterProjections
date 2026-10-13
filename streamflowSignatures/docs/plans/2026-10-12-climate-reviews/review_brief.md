# Blind review brief — climate metrics (guidelines §3.9)

## A. The specification as written by the co-authors (guidelines doc, synced 2026-10-12)

3.9 Climate Metrics

prcp_annual / _xSeason: total amount of precipitation in one year or season (depth, mm)

prcp_1d_max: maximum precipitation (depth, mm) in a water year

prcp_5d_max: maximum 5 consecutive day precipitation (depth, mm) total in a given water year

prcp_max_dowy: day of water year when maximum precipitation (i.e., prcp_1d_max) occurs

prcp_50_dowy: day of water year when precipitation hits 50% of total annual precipitation (i.e., prcp_annual)

prcp_cm_dowy: center of mass is calculated by summing the day of water year weighted rainfall (mm) and dividing by prcp_annual

prcp_ndays_abv_10mm: number of days in a water year with precipitation above 10mm

prcp_ndays_abv_25mm: number of days in a water year with precipitation above 25mm

prcp_ndays_abv_50mm: number of days in a water year with precipitation above 50mm

prcp_annual_sd: standard deviation for annual precipitation

prcp_annual_cv: coefficient of variability for annual precipitation

precip_intermit_annual: Annual precipitation intermittency; maximum number of consecutive days in water year without precipitation

precip_intermit_annual_dowy_start: day of water year corresponding to the start of the maximum number of consecutive days in water year without precipitation

precip_intermit_annual_dowy_end: day of water year corresponding to the end of the maximum number of consecutive days in water year without precipitation

prcp_gini : Gini index of precipitation, where 0 is perfectly uniform and 1 is all rain occurs on a single day (Sahbeni et al. 2023)

Where G is the Gini index, n is the total number of days, and yi is the precipitation for the ith day

prcp_seasonality_index: quantification of the variability of precipitation through the year. Hypothetically values can vary from 0-1.83, if every month of a year has the same magnitude of rainfall then the value is 0, if total annual rainfall occurs within one month then value is 1.83. Based on calculation from Sumner et al. 1988

where Ri is the total annual precipitation for the particular year under study and Xmm is the actual monthly precipitation for month n.

temp_max_annual: Maximum annual temperature (deg C)

temp_max_annual_dowy: day of water year on which temp_max_annual occurs

temp_min_annual: Minimum annual temp (deg C)

temp_min_annual_dowy: day of water year on which temp_min_annual occurs

temp_max_xSeason: Maximum temperature for each season

temp_min_xSeason: Minimum temperature for each season

frost_first_dowy: First day of water year below 0 degrees C

frost_last_dowy: Last day of water year below 0 degrees C between October and June

GDD_cum_annual: annual growing degree days are a metric of heat accumulation; calculated as the annual cumulative sum of daily growing degree days

GDD = ∑((Tmax+Tmin)/2)

GDD_cum_xSeason: growing degree days are a metric of heat accumulation; calculated as the cumulative sum of daily growing degree days for each season

GDD = ∑((Tmax+Tmin)/2)


## B. Decisions and definitions (user, 2026-10-12) — authoritative where the spec is silent

## 1. Decisions (user, 2026-10-12)

| # | Question | Decision |
|---|---|---|
| 1 | Rename shipped columns to the doc's `_annual` / `_xSeason` convention? | NO. Shipped names stay; the convention applies to new metrics only; relay a doc-side revert; a rename manifest only if the co-authors insist. |
| 2 | Season suffixes for new metrics | the Q seasonal totals' suffixes as in the code: `win`, `spr`, `sum`, `fal` (not the doc's `sumr`, `fall`). |
| 3 | Prefix normalisation | yes: `precip_intermit_*` → `prcp_intermit_*`; `GDD_cum_*` → `gdd_cum_*`; relay. |
| 4, 6, 11 | Snowfall metrics (§3.8 insert) | REMOVED from the doc by the user. Not implemented. |
| 5 | Growing degree days | base temperature 0 °C, clamped: Σ max(0, (Tmax+Tmin)/2 − 0); the user updated the doc. |
| 7 | Dry day | prcp ≤ 0.1 mm (not the ETCCDI 1 mm). |
| 8 | `prcp_annual_sd` / `_cv` | within-year spread of the daily values → annual signatures with the 16 statistics. |
| 9 | 5-day maximum, ties, dry spells | window inside the water year; ties → first occurrence; the first longest dry spell. |
| 10 | Frost | basin-mean tmin < 0 °C; first frost searched from Oct 1; last frost within Oct 1 – Jun 30; NA when no day qualifies. |
| 12 | Gini, seasonality index | Gini over all days of the water year, zeros included; Walsh–Lawler index on the 12 calendar months of the water year. |
| 13 | Year qualification | the snow-family precedent: climate years are evaluated only for streamflow-valid years (`preprocess_daily_data`), precipitation metrics on PPT-valid years, temperature metrics on temperature-valid years (new policy, same thresholds, no negative-value check). |
| 14 | Gating of NaN-prone timing metrics | the snow timing precedent: the record-anchored decade gate for `frost_first_dowy` and `frost_last_dowy`, anchored to the temperature-valid years. |
| 15 | Rerun | both windows (WY 1993–2025 @ 60 %, WY 1980–2025 @ 60 %), same config otherwise, as new product versions superseding #1 and #2; the delivered products are not rewritten. |
| 16 | "Slightly altered" acceptance | per-column R² ≥ 0.999 on shared gage-years, every maximum reported and explained. |
| 17 | Ports | Julia products first; Python and rpkg validated against them before publication. |
| 18 | Fold-ins | regenerate `flagged_for_high_na`; rerun the rpkg reference benchmark (constant-series Mann-Kendall fix). |
| 19 | Reviews | Codex + Sonnet 5.5 per language, blind (guidelines text + code, no implementation narrative), plus one numeric cross-language audit. |
| 20 | HydroShare | stage locally only; no upload before the user's review. |
| 21 | Relay list | yes, one line per item, incl. the manuscript season error. |

Open with the co-authors (relay): the published guidelines page still carries the snowfall
lines and the GDD formula without a base (the user's edits had not propagated at the
2026-10-12 re-sync); the Gini and seasonality equations are images.

## 2. The 39 metrics — definitions as implemented

Category **Climate** (ninth; the manuscript's §2 preamble already lists it). Two function
families in `julia/src/climate.jl`: `calculate_precipitation_metrics` (20 bases, PPT-valid
years) and `calculate_temperature_metrics` (19 bases, temperature-valid years). Every base
gets the 16 statistics; no per-gage scalar. Columns: 39 × 16 = **624** → the summary CSV
grows from 1,653 to **2,277**; annual signatures 100 → **139**.

Conventions shared by both families: a year is used only when its daily grid is complete
(365/366 rows, dowy 1..n, no NaN after preprocessing); seasons are the calendar-month sets
of `na_handling.seasonal_completeness.season_definitions` (DJF, MAM, JJA, SON) applied
WITHIN the water year, so `fal` = Oct + Nov of the water-year start plus Sep of its end,
exactly as `Qfal` and `fall_runoff_ratio`; dowy 1 = Oct 1; ties take the first occurrence;
no Q-derived seasonal mask is applied (these are pure climate series).

### Precipitation (daily P, mm; PPT-valid years)

| Base | Definition | Units | Year NaN when |
|---|---|---|---|
| `prcp_annual` | Σ P | mm | — |
| `prcp_win`, `prcp_spr`, `prcp_sum`, `prcp_fal` | Σ P over the season's months | mm | — |
| `prcp_1d_max` | max P | mm | — |
| `prcp_5d_max` | max over i of Σ P[i..i+4], window inside the water year | mm | — |
| `prcp_max_dowy` | dowy of the first max P | day | — |
| `prcp_50_dowy` | first dowy with cumulative P ≥ 0.5 Σ P | day | Σ P ≤ 0 |
| `prcp_cm_dowy` | Σ (dowy · P) / Σ P | day | Σ P ≤ 0 |
| `prcp_ndays_abv_10mm`, `_25mm`, `_50mm` | count(P > 10 / 25 / 50 mm) | days | — |
| `prcp_annual_sd` | sample sd (n − 1) of the daily P | mm | — |
| `prcp_annual_cv` | sd / mean of the daily P | – | mean ≤ 0 |
| `prcp_intermit_annual` | longest run of dry days (P ≤ 0.1 mm); 0 when every day is wet | days | — |
| `prcp_intermit_annual_dowy_start` / `_dowy_end` | first / last dowy of the first longest dry run | day | no dry day |
| `prcp_gini` | G = Σᵢ (2i − n − 1) y₍ᵢ₎ / (n Σ y), y ascending, all n days incl. zeros | – | Σ P ≤ 0 |
| `prcp_seasonality_index` | SI = (1/R) Σₘ₌₁¹² \|Xₘ − R/12\|, Xₘ monthly totals, R = Σ Xₘ (0 … 1.83) | – | R ≤ 0 |

### Temperature (daily tmin, tmax, °C; temperature-valid years)

| Base | Definition | Units | Year NaN when |
|---|---|---|---|
| `temp_max_annual`, `temp_max_annual_dowy` | max of daily tmax; dowy of its first occurrence | °C, day | — |
| `temp_min_annual`, `temp_min_annual_dowy` | min of daily tmin; dowy of its first occurrence | °C, day | — |
| `temp_max_win/spr/sum/fal` | max of daily tmax over the season's months | °C | — |
| `temp_min_win/spr/sum/fal` | min of daily tmin over the season's months | °C | — |
| `frost_first_dowy` | first dowy with tmin < 0 °C | day | no frost day (gated, §4) |
| `frost_last_dowy` | last dowy with tmin < 0 °C among days from Oct 1 to Jun 30 | day | no such day (gated, §4) |
| `gdd_cum_annual` | Σ max(0, (tmax + tmin)/2 − 0) | °C·day | — |
| `gdd_cum_win/spr/sum/fal` | the same over the season's months | °C·day | — |


## C. Part 1 global rules and the seasonal convention (guidelines)

### 1.1 Time conventions

All annual metrics are calculated on the water year, October 1 – September 30, naming convention `variableName`_annual

Seasonal periods are defined as: Winter: December–February; Spring: March–May; Summer: June–August; Fall: September–November. Naming conventions are `variableName`_xSeason, where xSeason is win, spr, sumr, fall

