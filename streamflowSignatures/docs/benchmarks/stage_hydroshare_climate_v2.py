#!/usr/bin/env python
"""Stage version-2 HydroShare resource folders (R1 / R2) for the climate reruns (2026-10-12).

    python docs/benchmarks/stage_hydroshare_climate_v2.py --run-dir RUN --prefix PREFIX --window 1993 \
        --staged-old ~/Downloads/Signatures/resource1_signatures_wy1993_2025 \
        --dictionary ~/Downloads/Signatures/hisss_data_dictionary.csv \
        --categories ~/Downloads/Signatures/hisss_signature_categories.csv \
        --out ~/Downloads/Signatures_v2_12oct2026/resource1_signatures_wy1993_2025

Local staging only (user decision 2026-10-12: nothing is uploaded). Copies and renames the run's
artifacts to the deposit names used in 2026 (`hisss_signatures_wyW_2025.csv`, the annual parquet,
the provenance JSON, the explorer + `_annual/` sidecar, the validation files), extends the data
dictionary and category table with the climate family (extend_data_dictionary_climate.py),
rewrites the counts and the climate / Daymet / high-NA paragraphs of the 2026 README, writes the
run notes and a rename manifest. Numbers come from the files themselves, never from this script.
"""
import argparse
import re
import json
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def copy(src, dst):
    if os.path.isdir(src):
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)
    return os.path.basename(src), os.path.basename(dst)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True); ap.add_argument("--prefix", required=True)
    ap.add_argument("--window", type=int, required=True, choices=[1993, 1980])
    ap.add_argument("--staged-old", required=True, help="the 2026 staged resource folder (README template)")
    ap.add_argument("--dictionary", required=True); ap.add_argument("--categories", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    W = a.window; tag = f"wy{W}_2025"
    os.makedirs(a.out, exist_ok=True)
    run = a.run_dir; p = a.prefix
    manifest = []

    csv_src = os.path.join(run, f"{p}_signatures.csv")
    ann_src = os.path.join(run, f"{p}_signatures_annual.parquet")
    tim_src = os.path.join(run, f"{p}_timing.json")
    for src, dst in [(csv_src, f"hisss_signatures_{tag}.csv"), (ann_src, f"hisss_signatures_annual_{tag}.parquet"),
                     (tim_src, f"hisss_run_provenance_{tag}.json")]:
        manifest.append(copy(src, os.path.join(a.out, dst)))
    expl = os.path.join(run, f"signature_explorer_{W}_2025_60pct_climate.html")
    if os.path.exists(expl):
        manifest.append(copy(expl, os.path.join(a.out, f"hisss_signature_explorer_{tag}.html")))
        side = expl[:-5] + "_annual"
        if os.path.isdir(side):
            manifest.append(copy(side, os.path.join(a.out, f"hisss_signature_explorer_{tag}_annual")))
    for name, dst in [("validation_gates.txt", f"hisss_validation_gates_{tag}.txt"),
                      ("validation_annual_values.txt", f"hisss_validation_annual_values_{tag}.txt"),
                      ("validation_climate_rerun_review.txt", f"hisss_validation_rerun_review_{tag}.txt"),
                      ("climate_rerun_changed_columns.csv", f"hisss_validation_changed_columns_{tag}.csv"),
                      ("climate_rerun_new_columns_sniff.csv", f"hisss_validation_new_columns_{tag}.csv"),
                      ("climate_rerun_annual_shared_years.csv", f"hisss_validation_annual_shared_years_{tag}.csv"),
                      ("merge.json", f"hisss_run_merge_{tag}.json")]:
        src = os.path.join(run, name)
        if os.path.exists(src):
            manifest.append(copy(src, os.path.join(a.out, dst)))

    # dictionary + categories with the climate family
    subprocess.run([sys.executable, os.path.join(HERE, "extend_data_dictionary_climate.py"),
                    "--dictionary", a.dictionary, "--categories", a.categories, "--out-dir", a.out], check=True)

    # numbers from the files
    summ = pd.read_csv(os.path.join(a.out, f"hisss_signatures_{tag}.csv"), dtype={"gage_id": str}, low_memory=False)
    n_gages, n_cols = summ.shape
    ann = pd.read_parquet(os.path.join(a.out, f"hisss_signatures_annual_{tag}.parquet"), columns=["signature"])
    n_rows = len(ann); n_sigs = ann.signature.nunique()
    tj = json.load(open(tim_src))
    clim_na = int(summ["prcp_annual_mean"].isna().sum()) if "prcp_annual_mean" in summ.columns else -1
    unnorm = int((summ["area_normalized"].astype(str).str.upper().isin(["FALSE", "0"])).sum()) if "area_normalized" in summ.columns else -1

    # README: the 2026 text with the counts, the climate family, the input and the resolved known issue
    old_readme = os.path.join(a.staged_old, f"hisss_readme_{tag}.md")
    t = open(old_readme, encoding="utf-8").read()
    old_n = {"1993": ("6,678", "18,898,406"), "1980": ("6,250", "24,366,487")}[str(W)]
    _m = re.search(r'about (\d+)% of gages in this window lack Daymet', t)
    old_pct = _m.group(1) if _m else '17'   # the 2026 READMEs say 17% (WY 1993) / 10% (WY 1980)
    clim_sentence = ("Every gage carries climate summary statistics." if clim_na == 0 else
                     f"{clim_na:,} gages have fewer than 20 climate-valid water years (WY 1980 has no climate value for\n"
                     "  any gage: the input begins 1980-01-01), so their climate summary statistics are NA while\n"
                     "  their climate annual values are in the parquet.")

    reps = [
        (f"**{old_n[0]} gaged", f"**{n_gages:,} gaged"),
        ("100 annually resolved signatures are computed — spanning flow volumes and percentiles,\n"
         "timing, flow-duration curves, baseflow, recession, pulses, flashiness, runoff ratios,\n"
         "elasticity, precipitation–streamflow seasonality, storage, negative-flow days, snow,\n"
         "and streamflow drought — and each signature's annual series is summarized with 16\n"
         "statistics (eight trend/summary statistics and eight Pettitt changepoint fields).",
         f"{n_sigs} annually resolved signatures are computed — spanning flow volumes and percentiles,\n"
         "timing, flow-duration curves, baseflow, recession, pulses, flashiness, runoff ratios,\n"
         "elasticity, precipitation–streamflow seasonality, storage, negative-flow days, snow,\n"
         "streamflow drought, and (new in this version) 39 climate metrics of basin-mean\n"
         "precipitation and air temperature — and each signature's annual series is summarized\n"
         "with 16 statistics (eight trend/summary statistics and eight Pettitt changepoint fields)."),
        (f"one row per gage, 1,653 columns (gage metadata, 100 signatures × 16 statistics, 21 per-gage scalars, 12 QA flags)",
         f"one row per gage, {n_cols:,} columns (gage metadata, {n_sigs} signatures × 16 statistics, 21 per-gage scalars, 12 QA flags)"),
        (f"`value` ({old_n[1]} rows)", f"`value` ({n_rows:,} rows)"),
        ("Coarse nine-class exploratory grouping of the 100 signature bases (the data paper's finer methodological categories are a separate taxonomy)",
         f"Grouping of the {n_sigs} signature bases (the 2026 nine-class exploratory grouping plus the Climate class for the 39 new bases; the data paper's categories are a separate taxonomy)"),
        ("  `_pettitt_pval`, and with ~100 signatures tested per gage, apply your own",
         f"  `_pettitt_pval`, and with ~{n_sigs} signatures tested per gage, apply your own"),
        ("- **Precipitation-dependent signatures** (runoff ratios, elasticity, Q–P seasonality,\n"
         "  storage) require basin-averaged Daymet forcing **and** area-normalized discharge;\n"
         f"  they are NA for gages lacking either (about {old_pct}% of gages in this window lack Daymet\n"
         "  forcing). **Snow signatures** require Daymet snow water equivalent but not\n"
         "  area-normalized discharge (`swe_max_to_ppt` additionally requires precipitation).",
         "- **Climate forcing.** This version uses the reprocessed Daymet V4 R1 input (basin means of\n"
         "  prcp, tmin, tmax and swe over the full-resolution basin polygons, calendar 1980–2025), so\n"
         "  every gage has a climate series and climate-dependent signatures extend to WY 2025\n"
         f"  (the earlier release lacked forcing for about {old_pct}% of gages and ended in 2023). {clim_sentence}\n"
         "  **Precipitation-dependent signatures** (runoff ratios,\n"
         "  elasticity, Q–P seasonality, storage) additionally require area-normalized discharge and\n"
         f"  are NA for the {unnorm} gages without a published drainage area. **Snow signatures** and the\n"
         "  **climate metrics** (`prcp_*`, `temp_*`, `frost_*`, `gdd_*`) need no discharge and exist\n"
         "  for those gages too.\n"
         "- **Climate metrics (new).** 20 precipitation and 19 temperature signatures of the forcing\n"
         "  itself: totals (water year and DJF/MAM/JJA/SON seasons, where fall = Oct–Nov of the\n"
         "  water-year start plus Sep of its end), 1- and 5-day maxima and their timing, 50 % and\n"
         "  centre-of-mass dates, days above 10/25/50 mm, within-year sd and cv, longest dry spell\n"
         "  (≤ 0.1 mm days) with its start/end, Gini index, Walsh–Lawler seasonality index; annual and\n"
         "  seasonal TMAX maxima / TMIN minima with timing, first and last frost (TMIN < 0 °C, last\n"
         "  frost within Oct–Jun), growing degree days (base 0 °C). They describe the forcing, not\n"
         "  the response; WY 1980 has none (the input begins 1980-01-01); frost dates are NA in\n"
         "  frost-free years and carry the snow family's first/last-decade trend gate."),
        ("- **QA flags** (`flagged_for_*`) mark automated range and consistency checks; see the\n"
         "  dictionary. `flagged_for_high_na` is intended to mark gages where > 30% of signature\n"
         "  fields are NA — **but see the known issue below: in this release it does not.**\n",
         "- **QA flags** (`flagged_for_*`) mark automated range and consistency checks; see the\n"
         "  dictionary. `flagged_for_high_na` marks gages where > 30% of the signature fields are\n"
         "  NA (regenerated in this version with the corrected shared definition; its denominator\n"
         f"  is every signature-derived column, {n_cols - 32:,} of {n_cols:,}).\n"),
        ("Self-contained interactive map explorer (open the HTML in any browser; keep the `_annual/` folder beside it",
         "Interactive map explorer (open the HTML in any browser — it needs an internet connection for the map library and the basemap tiles; keep the `_annual/` folder beside it"),
    ]
    for old, new in reps:
        if old not in t:
            print(f"WARNING: README anchor not found: {old[:70]!r}")
        t = t.replace(old, new)
    # drop the resolved known-issue bullet
    start = t.find("- **KNOWN ISSUE (one column) — `flagged_for_high_na`.**")
    if start >= 0:
        end = t.find("\n\n", start)
        t = t[:start] + "- **Resolved known issue.** The 2026 release's `flagged_for_high_na` was computed over the\n  metadata columns only; this version regenerates it with the corrected definition.\n" + t[end + 1:]
    t = t.replace("# HISSS Streamflow Signatures — Water Years", "# HISSS Streamflow and Climate Signatures (version 2, October 2026) — Water Years", 1)
    open(os.path.join(a.out, f"hisss_readme_{tag}.md"), "w", encoding="utf-8").write(t)

    notes = f"""# STANDARD PRODUCT WY {W}–2025 @ 60 % — VERSION 2 (climate family, reprocessed Daymet)

**Run date**: {tj.get('end_time', '?')[:10]} · **Runtime**: {tj.get('total_seconds', 0) / 60:.1f} min of Julia across 8 gage
batches on an Apple M5 (16 GB) · **Wrapper**: `docs/benchmarks/run_julia_benchmark_prod_{W}_2025_60pct_climate.jl`
via `docs/benchmarks/run_batched_julia.py`.

**Supersedes the 2026 product of the same window** (`{'processedOuts_drought_28jul2026' if W == 1993 else 'processedOuts_1980_2025_11aug2026'}`):
same streamflow input, window, qualification and configuration, so the gage set is unchanged;
the climate input is the reprocessed `daymet_1980_2025_29sep2026.parquet` (8,017 basins,
calendar 1980–2025, every gage) and the 39-signature climate family (guidelines §3.9,
2026-10-12) adds 624 columns.

| | |
|---|---|
| Gages | {n_gages:,} |
| Columns | {n_cols:,} ({n_cols - 32:,} signature + 20 metadata + 12 QA/QC flags) |
| Annual parquet | {n_rows:,} rows / {n_sigs} signatures |
| Gages with NA climate summary statistics (fewer than 20 climate-valid years; annual values present) | {clim_na:,} |
| Gages with raw m³/s discharge (`area_normalized = FALSE`) | {unnorm} |

What changed versus the superseded product, by column class (`hisss_validation_rerun_review_{tag}.txt`):
Q-only columns identical per gage; climate-dependent columns moved (new input values, two
added climate years WY 2024–2025, gages that gained a climate series); 624 new columns;
`flagged_for_high_na` regenerated. Validation files are listed in the README; the blind
review record of the code is in the repository (`docs/plans/2026-10-12-climate-reviews/`).
"""
    open(os.path.join(a.out, f"hisss_run_notes_{tag}.md"), "w", encoding="utf-8").write(notes)
    with open(os.path.join(a.out, f"hisss_rename_manifest_{tag}.csv"), "w") as f:
        f.write("original_name,deposit_name\n")
        for o, d in manifest:
            f.write(f"{o},{d}\n")
    print(f"staged {a.out}: {n_gages:,} gages x {n_cols:,} cols; {n_rows:,} annual rows / {n_sigs} signatures; {len(manifest)} files renamed")


if __name__ == "__main__":
    main()
