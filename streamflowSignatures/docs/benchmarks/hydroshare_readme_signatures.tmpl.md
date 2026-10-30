# HISSS Streamflow and Climate Signatures — Water Years ${W}–2025 (October 2026)

Part of the **Hydrologic Information Signatures and Summary Statistics (HISSS)** dataset:
hydrological signatures with trend and changepoint statistics for **${n_gages} gaged
watersheds** across the United States and Canada, computed over water years (WY)
**${W}–2025**. From quality-screened daily discharge (USGS National Water Information
System; Water Survey of Canada HYDAT) and basin-averaged Daymet climate/snow forcing,
${n_sigs} annually resolved signatures are computed in nine categories — flow volume, flow
duration, flow timing, flashiness, storage, drought, precipitation–streamflow, snow, and
climate (39 metrics of basin-mean precipitation and air temperature) — and each
signature's annual series is summarized with 16 statistics (eight trend/summary
statistics and eight Pettitt changepoint fields).

A companion resource provides the same product for WY ${companion}–2025 (DOI: to be assigned at
publication). Inputs, watershed geometry/attributes, and vegetation/land-cover
companions are published as separate resources in the HISSS collection (link below). The
computation code is at https://github.com/CZ-Sync/HISSS; the data paper will be linked at
publication.

## Files

| File | Description |
|---|---|
${files_table}

## Reading the data

```python
import pandas as pd
summary = pd.read_csv("hisss_signatures_${tag}.csv", dtype={"gage_id": str})
annual  = pd.read_parquet("hisss_signatures_annual_${tag}.parquet")
```

```r
summary <- readr::read_csv("hisss_signatures_${tag}.csv",
                           col_types = readr::cols(gage_id = readr::col_character()))
annual  <- arrow::read_parquet("hisss_signatures_annual_${tag}.parquet")
```

Always read `gage_id` as a **character string** — most US identifiers begin with a
leading zero that numeric parsing silently destroys. **Joining across HISSS resources:**
the tables of this resource and the daily streamflow and Daymet tables of the input-data
resource store US identifiers zero-padded, but the gage-metadata table of the input-data
resource stores every USGS identifier without its leading zeros (as retrieved, a number), and
the tables of the watershed-geometry and vegetation/land-cover resources store some long USGS
identifiers that way. Join on the zero-stripped
form of the all-digit (USGS) identifiers on both sides — R:
`ifelse(grepl("^[0-9]+$$", id), sub("^0+", "", id), id)`; Python:
`i.lstrip("0") if i.isdigit() else i` — and never re-pad, because US identifiers are not
all 8 digits. Canadian identifiers contain letters (e.g. `05BL015`) and are left as they are.

Missing values in the summary CSV appear as empty cells or as the literal `NaN`; both
mean "not computable" and both are read as missing by pandas and readr.

## Conventions

- **Water year**: 1 October – 30 September, labeled by the ending calendar year.
- **Column naming**: `{signature}_{statistic}`, e.g. `Qann_senn_slp` is the Theil–Sen
  slope of annual total streamflow. The eight summary statistics are `_mean`,
  `_median`, `_senn_slp`, `_linear_slp`, `_spearman_rho`, `_spearman_pval`, `_mk_rho`,
  `_mk_pval`; the eight Pettitt changepoint fields are `_pettitt_cp_year`,
  `_pettitt_pval`, `_pettitt_pre_mean`, `_pettitt_post_mean`, `_pettitt_delta_mean`,
  `_pettitt_pct_change`, `_pettitt_pre_mk_pval`, `_pettitt_post_mk_pval`. All
  definitions and units are in `hisss_data_dictionary.csv`.
- **Annual values parquet**: a missing row and a stored NaN are equivalent — both mean
  "not computable for that gage-year".
- **Units**: daily flow is mm/day and flow totals are mm **where
  `area_normalized = TRUE`**. For the small number of Canadian gages with no published
  drainage area (`area_normalized = FALSE`), daily-flow magnitudes and percentiles
  remain in native m³/s, and seasonal/annual totals are sums of daily m³/s values
  (m³/s·days; multiply by 86,400 for volume in m³). Filter to
  `area_normalized == TRUE` before any cross-gage comparison of unit-carrying columns.
- **Trend statistics may be NA where means are not.** Slopes, correlations, and their
  p-values are reported only when a signature's annual series has ≥ 20 valid values,
  ≥ 60% valid overall, and ≥ 80% valid in its first and last decade (recession and
  elasticity signatures are exempt). Mean, median, and Pettitt fields are still
  reported when only the trend gate fails.
- **Pettitt fields are computed over a fixed WY 1980–2024 changepoint window**${pettitt_note} A
  changepoint year is reported when at least 20 valid values fall in the window with
  at least 10 on each side of the candidate split; otherwise all Pettitt fields are
  NA. The reported year is not itself evidence of change — judge significance with
  `_pettitt_pval`, and with ~${n_sigs} signatures tested per gage, apply your own
  multiple-testing correction for per-gage analyses.
- **Streamflow input.** The daily discharge is the HISSS February 2026 compile after its
  qualifier correction (the input-data resource holds the corrected table and the
  correction ledger): the compile had set every USGS day whose qualification code
  carries an approved remark (`A R` revised, `A <` / `A >` censored, `A [4]`) to NA —
  ${n_led_all} approved days on ${n_led_gages_all} compiled gages; those days were re-retrieved from
  NWIS (26 October 2026) and restored before this product was computed. In this window
  ${n_restored_days} such days on ${n_restored_gages} gages are finite; ${n_gained_years} water years on ${n_gained_gages} gages are
  valid only because of them. Record ends 2013–2025 carry most of them.
- **Climate forcing.** Basin means of Daymet V4 R1 prcp, tmin, tmax and swe over the
  full-resolution basin polygons, calendar 1980–2025 (the input-data resource holds the
  table), so every gage has a climate series and climate-dependent signatures extend to
  WY 2025.${clim_note} **Precipitation-dependent signatures** (runoff ratios, elasticity,
  Q–P seasonality, storage) additionally require area-normalized discharge and are NA for
  the ${unnorm} gages without a published drainage area. **Snow signatures** and the
  **climate metrics** (`prcp_*`, `temp_*`, `frost_*`, `gdd_*`) need no discharge and exist
  for those gages too.
- **Climate frame.** Each gage's Daymet series enters the per-year screening on Daymet's
  own dates inside the gage's record span, so a day absent from the streamflow record
  keeps its climate values instead of being interpolated across (${n_J} gages have such
  days in valid years). Precipitation and degree-day totals equal the daily sums.
- **Climate metrics.** 20 precipitation and 19 temperature signatures of the forcing
  itself: totals (water year and DJF/MAM/JJA/SON seasons, where fall = Oct–Nov of the
  water-year start plus Sep of its end), 1- and 5-day maxima and the timing of the 1-day maximum, 50 % and
  centre-of-mass dates, days above 10/25/50 mm, within-year sd and cv, longest dry spell
  (≤ 0.1 mm days) with its start/end, Gini index, Walsh–Lawler seasonality index; annual
  (with timing) and seasonal TMAX maxima / TMIN minima, first and last frost (TMIN < 0 °C, last
  frost within Oct–Jun), growing degree days (base 0 °C). They describe the forcing, not
  the response;${wy1980_note} frost dates are NA in frost-free years and carry the snow
  family's first/last-decade trend gate.
- **Gage set.** A gage qualifies with ≥ 20 valid water years that are ≥ 60 % of the
  window's years. The ${n_gages} gages are those that qualified under the February 2026 compile;
  the qualifier correction restored days but did not change the gage set, so ${n_wq} gages that
  qualify only under the corrected table are not included. ${gage_set_files}
- **QA flags** (`flagged_for_*`) mark automated range and consistency checks; see the
  dictionary. `flagged_for_high_na` is TRUE for the ${n_high_na} gages where more than 30 % of
  the ${n_sig_cols} signature-derived columns are NA.
- **Known issues.** (i) `ice_affected_days_total` is 0 for every gage (the implementation
  reports 0 regardless of the record), so a zero does not establish their
  absence — use the daily `flag` column of the input-data resource instead. (ii)
  `flagged_for_qann_range` catches only ${n_qann_unnorm} of the ${unnorm} gages whose discharge is in
  native m³/s; filter on `area_normalized` rather than on that flag.

## Do not compare record-dependent signatures across windows

Drought metrics (fixed thresholds from this window's own record), period-of-record
pulse metrics (`*_all`), elasticity, and the recession-parameterized baseflow indices
(`BFI_*_param`) are computed from **this product's window**. Values are comparable
within this resource but **not** with the WY ${companion}–2025 companion resource, and the
annual values must not be re-aggregated onto a different analysis window. Note also
that `drought_duration_fixed_p10` is statistically largely redundant with
`n_low_pulses_all × dur_low_pulses_all` — do not present them as independent evidence.

## License and citation

License: CC-BY 4.0. Please cite the HISSS data paper and this resource (DOI: to be
assigned at publication; the citation will be added here). Source data credits: USGS
NWIS; Water Survey of Canada HYDAT; Daymet V4 R1 (ORNL DAAC).

This resource is part of the HISSS collection on HydroShare:
https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/
(collection DOI upon publication: https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768).
