# Relay list for the co-authors — 2026-10-12 sync (guidelines §3.9 + manuscript)

One line per item, "replace X with Y" where possible (user preference 2026-09-04). Sources:
`docs/reconciliation/guidelines_todos.md` and `manuscript_log.md` (2026-10-12 entries);
decisions: `docs/plans/2026-10-12-climate-metrics-plan.md` §1.

## Guidelines doc (Summary Documentation for Streamflow Signatures)

1. Part 1.1 naming convention: keep the shipped column names (Qann, Qwin/Qspr/Qsum/Qfal,
   Negative_ann, annual/winter/spring/summer/fall_runoff_ratio, Flow_Reversals_annual/
   _winter/…); apply `variableName_annual` / `variableName_<season>` to NEW metrics only,
   with the season suffixes the code already uses: `win`, `spr`, `sum`, `fal` (not `sumr`,
   `fall`). Revert the six renames in 3.1, 3.4 and 3.7 (Q_annual → Qann, Qsumr → Qsum,
   Qfall → Qfal, Negative_annual → Negative_ann, runoff_ratio_annual / _Xseason → the five
   `<season>_runoff_ratio` names, Flow_Reversals_xSeason → _winter/_spring/_summer/_fall).
2. §3.9 names as implemented: `precip_intermit_annual` (+ `_dowy_start`, `_dowy_end`) →
   `prcp_intermit_annual` (…); `GDD_cum_annual` / `GDD_cum_xSeason` → `gdd_cum_annual`,
   `gdd_cum_win/spr/sum/fal`; every `_xSeason` → `_win`, `_spr`, `_sum`, `_fal`.
3. §3.9 GDD: replace "GDD = ∑((Tmax+Tmin)/2)" with "GDD = Σ max(0, (Tmax + Tmin)/2 − 0 °C)"
   (base temperature 0 °C, negative days contribute 0) — confirm the published page shows it.
4. §3.9 dry day: add "a dry day has precipitation ≤ 0.1 mm" to the three
   `prcp_intermit_annual*` lines; add "strictly above" to the three `prcp_ndays_abv_*mm` lines.
5. §3.9 `prcp_5d_max`: add "the five-day window lies within the water year";
   `prcp_max_dowy` and the dry-spell dates: add "first occurrence on ties".
6. §3.9 `prcp_annual_sd` / `prcp_annual_cv`: replace "standard deviation for annual
   precipitation" with "sample standard deviation of the daily precipitation within the water
   year (mm)" and "coefficient of variation of the daily precipitation within the water year
   (sd / mean)".
7. §3.9 `prcp_gini`: type the formula — G = Σᵢ (2i − n − 1)·y₍ᵢ₎ / (n·Σy), y sorted ascending
   over all n days of the water year (zeros included); `prcp_seasonality_index`: SI = (1/R)·
   Σₘ₌₁¹² |Xₘ − R/12| over the water year's calendar months (Walsh & Lawler 1981; Sumner 1988);
   both are pasted images today.
8. §3.9 frost: `frost_first_dowy` → "first day of the water year with basin-mean tmin < 0 °C
   (NA when none)"; `frost_last_dowy` → "last day with tmin < 0 °C between Oct 1 and Jun 30
   (NA when none)". Temperature extremes: "maximum of daily tmax" / "minimum of daily tmin".
9. §3.9 header: add "Function: calculate_precipitation_metrics(), calculate_temperature_metrics().
   Requires Daymet prcp, tmin, tmax" and a Requirements block: precipitation metrics on
   PPT-qualified years, temperature metrics on years where both tmin and tmax pass the same
   screen (≤ 30 missing days, gaps ≤ 3 days interpolated; no sign check); frost dates carry
   the snow family's first/last-decade gate; 16 statistics per metric, no scalars.
10. §3.8: the five snowfall lines are deleted on the editing copy — confirm the published page
    follows; §3.9 is the ninth category "Climate" (precipitation + temperature families).
11. Part 1: the new "1.1 Streamflow data sources" duplicates the number of "1.1 Time
    conventions" → renumber 1.1–1.4; 3.1 "daicly" → "daily".
12. Part 6 "Resource Collections" still carries working notes ("Insert the table…", "Here is the
    original…") — replace with the resources table or remove; References is now Part 7.
13. Carried forward: 3.3 parameterized-BFI interval "[0, 1]" → "(0, 1)"; add Hecht et al. 2024
    to the references; "recession_alpha_point_cloud" → the shipped
    `recession_alpha_point_cloud_linear_reservoir`; Part 1.2 "(see 3.12)" → "(see 3.1)".
14. §3.9 `prcp_50_dowy` (from the blind reviews, 2026-10-12): add "the first day on which the
    cumulative total reaches 50 % of the water-year total, to within 1e-12 of the total (so an
    exact tie is not decided by floating-point rounding); NA when the total is 0" — an
    implementer decision confirmed by the user 2026-10-13 (plan §1, item 22).
15. §3.9 conventions worth one line each: `prcp_gini` ranges from 0 to (n − 1)/n, not to 1;
    `prcp_max_dowy` is day 1 in an all-zero year (none exists in the input); Daymet's 365-day
    calendar omits Dec 31 of leap years, which the preprocessor fills linearly, so the water
    year following a leap calendar year carries one interpolated day in every metric.

## Manuscript

1. §2.1.1: replace "Seasonal periods are defined as: Winter: January-March, Spring: April-June,
   Summer: July - September, Fall: October-December" with "Winter: December–February, Spring:
   March–May, Summer: June–August, Fall: September–November" (the code, the guidelines and the
   sentence deleted from §2.2.1 all use these).
2. §2 preamble "121 signatures across nine signature families (…, climate)" vs §2.2.1 "eight
   signature families … 121": once the climate metrics ship, both become "160 signatures
   (139 annually resolved + 21 diagnostics) across nine categories"; until then "eight".
3. §2.1.3: after the rerun on the reprocessed Daymet input, replace "covering 1980–2023" with
   "covering 1980–2025", "gdptools" with "exactextract-based area weighting over the
   full-resolution basin polygons" and add "for all 8,017 basins (every gage with a boundary
   plus the 53 basins > 100,000 km²)"; the four St. Johns River basins that were NaN now have
   series.
4. §2.2.1 / §2.2.2: restore "any WY-season with fewer than 80 % of raw observations passing QA
   has its seasonal metrics set to NA" and "we computed 12 automated data-quality flags", or
   state that they were dropped deliberately — both are still implemented.
5. §4: "median 698 km², 25 % > 2,700 km², 25 % < 180 km²" describes the WY 1980–2025 product
   (184 / 698 / 2,718 km²); the WY 1993–2025 product gives 166 / 645 / 2,500 km² — say which.
6. §5.1.2: "(see resource README for details))" has a double closing parenthesis.
7. Carried forward (2026-09-29): gap rule "those days were set to NA" → "the water year was
   excluded"; "at least 20 qualifying water years across its full period of record" → "within
   the analysis window"; 85,000 vs 100,000 km²; cross-references after the §2.1.1/§2.1.2 swap;
   title/abstract naming; citations (Albers 2026 vs 2017, Myneni DOI, Hatchett year, Condon
   2010 → 2020, "htpps://"); "Claude Code 0.145.0"; §5.1.2 "0103500".
