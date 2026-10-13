#!/usr/bin/env python
"""Extend the HydroShare data dictionary and category table with the climate family (2026-10-12).

    python docs/benchmarks/extend_data_dictionary_climate.py \
        --dictionary ~/Downloads/Signatures/hisss_data_dictionary.csv \
        --categories ~/Downloads/Signatures/hisss_signature_categories.csv \
        --out-dir STAGING_DIR [--template-base D50_day]

Appends 39 bases x 16 statistics = 624 rows to the dictionary (category "Climate") and 39 rows
to the category table, reusing the 16 per-statistic templates (statistic label, unit pattern,
definition prefix) from an existing base's rows so the wording stays uniform. The new rows are
inserted after the last signature row and before the per-gage scalars / quality flags, in the
summary CSV's column order (alphabetical within the signature block). Writes
<out-dir>/hisss_data_dictionary.csv and <out-dir>/hisss_signature_categories.csv and prints
the row counts. Base definitions and units: docs/SIGNATURES.md §14 / guidelines §3.9.
"""
import argparse
import csv
import os
import sys

SUFFIX_ORDER = ["_linear_slp", "_mean", "_median", "_mk_pval", "_mk_rho", "_pettitt_cp_year",
                "_pettitt_delta_mean", "_pettitt_pct_change", "_pettitt_post_mean",
                "_pettitt_post_mk_pval", "_pettitt_pre_mean", "_pettitt_pre_mk_pval",
                "_pettitt_pval", "_senn_slp", "_spearman_pval", "_spearman_rho"]

SEASON_TXT = {"win": "winter (Dec-Feb)", "spr": "spring (Mar-May)", "sum": "summer (Jun-Aug)",
              "fal": "fall (Sep of the water year's end plus Oct-Nov of its start)"}

# base -> (units, definition)
BASES = {
    "prcp_annual": ("mm", "Water-year total precipitation (basin-mean Daymet prcp)."),
    "prcp_1d_max": ("mm", "Maximum daily precipitation in the water year."),
    "prcp_5d_max": ("mm", "Maximum total precipitation of any 5 consecutive days inside the water year."),
    "prcp_max_dowy": ("day of water year", "Day of the water year of the maximum daily precipitation (first day on ties; day 1 = Oct 1)."),
    "prcp_50_dowy": ("day of water year", "First day of the water year on which cumulative precipitation reaches 50% of the water-year total (to within 1e-12 of the total, so a mathematical tie is not decided by floating-point rounding); NA when the total is 0."),
    "prcp_cm_dowy": ("day of water year", "Precipitation centre of mass: sum(day-of-water-year x daily precipitation) / water-year total; NA when the total is 0."),
    "prcp_ndays_abv_10mm": ("days", "Number of days with daily precipitation strictly above 10 mm."),
    "prcp_ndays_abv_25mm": ("days", "Number of days with daily precipitation strictly above 25 mm."),
    "prcp_ndays_abv_50mm": ("days", "Number of days with daily precipitation strictly above 50 mm."),
    "prcp_annual_sd": ("mm", "Sample standard deviation (n-1) of the daily precipitation values within the water year."),
    "prcp_annual_cv": ("dimensionless", "Coefficient of variation of the daily precipitation values within the water year (sd / mean); NA when the mean is 0."),
    "prcp_intermit_annual": ("days", "Precipitation intermittency: longest run of consecutive dry days (daily precipitation <= 0.1 mm) in the water year; 0 when every day is wet."),
    "prcp_intermit_annual_dowy_start": ("day of water year", "First day of the longest dry run (first such run on ties); NA when the year has no dry day."),
    "prcp_intermit_annual_dowy_end": ("day of water year", "Last day of the longest dry run (first such run on ties); NA when the year has no dry day."),
    "prcp_gini": ("dimensionless", "Gini index of the daily precipitation values over every day of the water year (zeros included): 0 = uniform, (n-1)/n = all precipitation on one day (Sahbeni et al. 2023); NA when the total is 0."),
    "prcp_seasonality_index": ("dimensionless", "Walsh and Lawler (1981) seasonality index: (1/R) x sum over the water year's 12 calendar months of |monthly total - R/12|, R = water-year total; 0 = equal monthly totals, 1.83 = all precipitation in one month (Sumner 1988); NA when R is 0."),
    "temp_max_annual": ("deg C", "Maximum of the daily basin-mean maximum air temperature (Daymet tmax) in the water year."),
    "temp_max_annual_dowy": ("day of water year", "Day of the water year of the annual maximum temperature (first day on ties)."),
    "temp_min_annual": ("deg C", "Minimum of the daily basin-mean minimum air temperature (Daymet tmin) in the water year."),
    "temp_min_annual_dowy": ("day of water year", "Day of the water year of the annual minimum temperature (first day on ties)."),
    "frost_first_dowy": ("day of water year", "First day of the water year with basin-mean daily minimum temperature below 0 deg C; NA in frost-free years. Trend statistics are suppressed unless >= 80% of the first and last decade of the temperature record have a frost day."),
    "frost_last_dowy": ("day of water year", "Last day with basin-mean daily minimum temperature below 0 deg C between Oct 1 and Jun 30; NA when no such day. Trend statistics are suppressed unless >= 80% of the first and last decade of the temperature record have a frost day."),
    "gdd_cum_annual": ("deg C-days", "Growing degree days: water-year sum of max(0, (tmax + tmin)/2 - 0 deg C) over the daily basin-mean temperatures."),
}
for s, txt in SEASON_TXT.items():
    BASES[f"prcp_{s}"] = ("mm", f"Total precipitation of the {txt} season of the water year.")
    BASES[f"temp_max_{s}"] = ("deg C", f"Maximum of the daily maximum temperature over the {txt} season of the water year.")
    BASES[f"temp_min_{s}"] = ("deg C", f"Minimum of the daily minimum temperature over the {txt} season of the water year.")
    BASES[f"gdd_cum_{s}"] = ("deg C-days", f"Growing degree days (base 0 deg C) summed over the {txt} season of the water year.")
assert len(BASES) == 39, len(BASES)

COMMON = (" Computed on water years whose precipitation / temperature series passed the same screen as streamflow "
          "(<= 30 missing days, gaps <= 3 days interpolated; Daymet omits Dec 31 of leap calendar years, so the "
          "following water year has one linearly filled day); pure climate, so present for the un-normalized gages too.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dictionary", required=True); ap.add_argument("--categories", required=True)
    ap.add_argument("--out-dir", required=True); ap.add_argument("--template-base", default="D50_day")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with open(a.dictionary, newline="") as f:
        rows = list(csv.DictReader(f)); fields = list(rows[0].keys())
    tmpl = {}
    tb_unit = None
    for r in rows:
        if r["signature"] == a.template_base:
            suffix = r["column_name"][len(a.template_base):]
            tb_unit = tb_unit or ("day of water year" if a.template_base.endswith("_day") else None)
            prefix = r["definition"].split(" Signature: ")[0]
            tmpl[suffix] = dict(statistic=r["statistic"], units=r["units"], prefix=prefix)
    missing = [s for s in SUFFIX_ORDER if s not in tmpl]
    if missing or tb_unit is None:
        sys.exit(f"template base {a.template_base} lacks suffixes {missing} or has no known unit")
    existing = {r["column_name"] for r in rows}
    new_rows = []
    for base in sorted(BASES):
        unit, definition = BASES[base]
        for suffix in SUFFIX_ORDER:
            t = tmpl[suffix]
            units = t["units"].replace(tb_unit, unit)
            col = base + suffix
            if col in existing:
                sys.exit(f"{col} already in the dictionary")
            new_rows.append({"column_name": col, "category": "Climate", "signature": base,
                             "statistic": t["statistic"], "units": units,
                             "definition": f"{t['prefix']} Signature: {definition}{COMMON}"})
    # insert in the summary CSV's column order: the signature block is alphabetical, followed by
    # the per-gage scalars and the quality flags
    sig_rows = [r for r in rows if r["category"] not in ("Gage metadata", "Per-gage scalar", "Quality flag")]
    meta_rows = [r for r in rows if r["category"] == "Gage metadata"]
    tail_rows = [r for r in rows if r["category"] in ("Per-gage scalar", "Quality flag")]
    merged = sorted(sig_rows + new_rows, key=lambda r: r["column_name"])
    out_rows = meta_rows + merged + tail_rows
    assert len(out_rows) == len(rows) + 624, (len(out_rows), len(rows))
    with open(os.path.join(a.out_dir, "hisss_data_dictionary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(out_rows)

    with open(a.categories, newline="") as f:
        crows = list(csv.DictReader(f)); cfields = list(crows[0].keys())
    have = {r["signature"] for r in crows}
    for base in sorted(BASES):
        if base not in have:
            crows.append({"category": "Climate", "signature": base})
    with open(os.path.join(a.out_dir, "hisss_signature_categories.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cfields); w.writeheader(); w.writerows(crows)
    print(f"dictionary: {len(rows)} -> {len(out_rows)} rows; categories: {len(have)} -> {len(crows)} signatures; out: {a.out_dir}")


if __name__ == "__main__":
    main()
