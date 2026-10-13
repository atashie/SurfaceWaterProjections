#!/usr/bin/env python
"""Output review of a climate rerun against the superseded standard product (2026-10-12).

    python docs/benchmarks/compare_climate_rerun.py --new NEW.csv --old OLD.csv
        [--new-annual NEW_annual.parquet --old-annual OLD_annual.parquet]
        [--categories docs/signature_categories.csv] [--out-dir DIR] [--r2-min 0.999]

Three column classes, three expectations (plan: docs/plans/2026-10-12-climate-metrics-plan.md §6):
  1. Q-ONLY columns (every family that never reads PPT/SWE/TMIN/TMAX) must be IDENTICAL
     per gage, NaN == NaN. Any difference is a gate failure.
  2. CLIMATE-DEPENDENT existing columns (runoff ratios, elasticity, Q-P seasonality,
     storage, snow + their scalars + flagged_for_high_na) are EXPECTED to move: new input
     values, two added climate years (WY 2024-2025), gages gaining a climate series. The
     summary statistics therefore change for every gage (information only, reported per
     column with R², max |diff|, gained / lost cells). The GATES (decision 16, 2026-10-12:
     R² >= --r2-min on shared GAGE-YEARS, maxima reported and explained) are:
       2a. annual values on shared gage-years, per signature: continuous signatures need
           identity R² >= --r2-min and a 99th-percentile relative difference <= --rel-p99;
           discrete signatures (days, counts, ssm) are reported as the fraction of shared
           gage-years that changed; elasticity is reported only (a ratio of year-to-year
           differences, unstable by construction).
       2b. the old product's _mean recomputed from the NEW annual values on the OLD year
           set, per signature: R² >= --r2-min — isolates the input effect from the added
           years at the summary level.
  3. NEW columns (the 39 climate bases x 16 fields) must be populated (finite fraction
     >= 0.5 per field, reported per base) and pass sniff ranges on their _mean.
With the annual parquets, the climate-dependent signatures are also compared on SHARED
gage-years (where the expectation is agreement at the Daymet-agreement level) and the
added gage-years are counted — the two effects the summary statistics confound.
Exit 1 when any gate fails.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
import pandas as pd

STAT_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
                 "_mean", "_median", "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
                 "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
                 "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
CLIMATE_FUNCTIONS = {"analyze_Q_PPT_relationships", "calculate_streamflow_elasticity", "calculate_qp_seasonality",
                     "calculate_average_storage", "calculate_snow_metrics",
                     "calculate_precipitation_metrics", "calculate_temperature_metrics"}
NEW_FUNCTIONS = {"calculate_precipitation_metrics", "calculate_temperature_metrics"}
METADATA = {"gage_id", "latitude", "longitude", "basin_area", "gage_type", "num_water_years",
            "start_water_year", "end_water_year", "area_normalized", "NDAMS_2009", "MAJ_DDENS_2009",
            "STOR_NID_2009", "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
            "CLASS", "RHBN", "REGULATED", "human_interference_class"}
# Sniff ranges on the _mean of the new bases (physical plausibility over North America)
SNIFF = {
    "prcp_annual": (50, 6000), "prcp_win": (0, 3000), "prcp_spr": (0, 3000), "prcp_sum": (0, 3000),
    "prcp_fal": (0, 3000), "prcp_1d_max": (5, 600), "prcp_5d_max": (10, 1500), "prcp_max_dowy": (1, 366),
    "prcp_50_dowy": (30, 330), "prcp_cm_dowy": (60, 300), "prcp_ndays_abv_10mm": (0, 200),
    "prcp_ndays_abv_25mm": (0, 120), "prcp_ndays_abv_50mm": (0, 60), "prcp_annual_sd": (0.5, 60),
    "prcp_annual_cv": (0.5, 8), "prcp_intermit_annual": (2, 250), "prcp_intermit_annual_dowy_start": (1, 366),
    "prcp_intermit_annual_dowy_end": (1, 366), "prcp_gini": (0.3, 1.0), "prcp_seasonality_index": (0.0, 1.5),
    "temp_max_annual": (0, 50), "temp_max_annual_dowy": (1, 366), "temp_min_annual": (-60, 15),
    "temp_min_annual_dowy": (1, 366), "temp_max_win": (-30, 40), "temp_max_spr": (-5, 50),
    "temp_max_sum": (5, 50), "temp_max_fal": (0, 50), "temp_min_win": (-60, 20), "temp_min_spr": (-50, 25),
    "temp_min_sum": (-25, 30), "temp_min_fal": (-50, 25), "frost_first_dowy": (1, 366),
    "frost_last_dowy": (1, 273), "gdd_cum_annual": (0, 12000), "gdd_cum_win": (0, 4000),
    "gdd_cum_spr": (0, 4000), "gdd_cum_sum": (0, 4000), "gdd_cum_fal": (0, 4000),
}


_SUFFIXES_LONGEST_FIRST = sorted(STAT_SUFFIXES, key=len, reverse=True)


def base_of(col):
    # longest suffix first: "_pettitt_delta_mean" must win over "_mean"
    for s in _SUFFIXES_LONGEST_FIRST:
        if col.endswith(s):
            return col[: -len(s)]
    return col


def load_categories(path):
    func_of = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            func_of[row["signature"]] = row["function"]
    return func_of


def identical(a, b):
    a = a.to_numpy(dtype=object); b = b.to_numpy(dtype=object)
    out = np.ones(len(a), dtype=bool)
    for i in range(len(a)):
        x, y = a[i], b[i]
        xn = x is None or (isinstance(x, float) and np.isnan(x)) or (isinstance(x, str) and x == "")
        yn = y is None or (isinstance(y, float) and np.isnan(y)) or (isinstance(y, str) and y == "")
        if xn or yn:
            out[i] = xn and yn
        else:
            try:
                out[i] = float(x) == float(y)
            except (TypeError, ValueError):
                out[i] = str(x) == str(y)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True); ap.add_argument("--old", required=True)
    ap.add_argument("--new-annual"); ap.add_argument("--old-annual")
    ap.add_argument("--categories", default=os.path.join(os.path.dirname(__file__), "..", "signature_categories.csv"))
    ap.add_argument("--out-dir", default=None); ap.add_argument("--r2-min", type=float, default=0.999)
    ap.add_argument("--allow-low", default="", help="comma-separated columns whose low R² is reported, not failed")
    ap.add_argument("--rel-p99", type=float, default=1e-3, help="gate 2a: 99th-percentile relative difference on shared gage-years")
    a = ap.parse_args()
    out_dir = a.out_dir or os.path.dirname(os.path.abspath(a.new))
    allow_low = set(x for x in a.allow_low.split(",") if x)
    func_of = load_categories(a.categories)

    new = pd.read_csv(a.new, dtype={"gage_id": str}, low_memory=False)
    old = pd.read_csv(a.old, dtype={"gage_id": str}, low_memory=False)
    results = []
    def gate(name, ok, detail):
        results.append((name, ok, detail)); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    # --- gage sets and column classes ---
    sn, so = set(new.gage_id), set(old.gage_id)
    gate("gage set identical", sn == so, f"new={len(sn)} old={len(so)} shared={len(sn & so)}")
    shared_gages = sorted(sn & so)
    new_i = new.set_index("gage_id").loc[shared_gages]; old_i = old.set_index("gage_id").loc[shared_gages]
    old_cols = [c for c in old.columns if c != "gage_id"]; new_cols = [c for c in new.columns if c != "gage_id"]
    added = [c for c in new_cols if c not in old.columns]; dropped = [c for c in old_cols if c not in new.columns]
    gate("no column dropped", not dropped, f"dropped={len(dropped)} {dropped[:5]}")
    gate("624 columns added (39 bases x 16)", len(added) == 624, f"added={len(added)}")

    def cls(col):
        if col in METADATA: return "metadata"
        if col.startswith("flagged_"): return "flag"
        f = func_of.get(base_of(col))
        if f is None: return "unregistered"
        if f in NEW_FUNCTIONS: return "new"
        if f in CLIMATE_FUNCTIONS: return "climate"
        return "qonly"
    classes = {c: cls(c) for c in new_cols}
    counts = pd.Series(classes).value_counts().to_dict()
    print("column classes:", counts)
    unregistered = [c for c, k in classes.items() if k == "unregistered"]
    gate("every column registered in signature_categories.csv", not unregistered, f"unregistered={unregistered[:8]}")

    # --- 1. Q-only (and metadata) columns identical ---
    qcols = [c for c in old_cols if classes.get(c) in ("qonly", "metadata")]
    bad_q = []
    for c in qcols:
        eq = identical(new_i[c], old_i[c])
        if not eq.all():
            bad_q.append((c, int((~eq).sum())))
    gate("Q-only + metadata columns identical per gage", not bad_q,
         f"{len(qcols)} columns; differing: {len(bad_q)} {bad_q[:8]}")

    # --- 2. climate-dependent existing columns: expected to move ---
    rows = []
    low = []
    ccols = [c for c in old_cols if classes.get(c) == "climate"]
    flag_cols = [c for c in old_cols if classes.get(c) == "flag"]
    for c in ccols + flag_cols:
        xo = pd.to_numeric(old_i[c], errors="coerce").astype(float); xn = pd.to_numeric(new_i[c], errors="coerce").astype(float)
        fo, fn = xo.notna(), xn.notna()
        both = fo & fn
        gained = int((fn & ~fo).sum()); lost = int((fo & ~fn).sum())
        n_both = int(both.sum())
        r2 = np.nan; maxabs = np.nan; n_changed = 0
        if n_both > 0:
            d = (xn[both] - xo[both]).abs()
            maxabs = float(d.max()); n_changed = int((d > 0).sum())
            if n_both >= 3 and xo[both].std() > 0:
                r2 = float(np.corrcoef(xo[both], xn[both])[0, 1] ** 2)
        rows.append(dict(column=c, base=base_of(c), n_shared_finite=n_both, n_changed=n_changed,
                         max_abs_diff=maxabs, r2=r2, gained=gained, lost=lost))
        if classes[c] == "climate" and n_both >= 30 and not np.isnan(r2) and r2 < a.r2_min and c not in allow_low:
            low.append((c, round(r2, 6), n_both))
    cd = pd.DataFrame(rows)
    cd.to_csv(os.path.join(out_dir, "climate_rerun_changed_columns.csv"), index=False)
    print(f"[INFO] climate-dependent summary columns (expected to move: 2 added climate years + gages gaining climate): "
          f"{len(ccols)} columns, {len(low)} with identity R² < {a.r2_min} on shared gages -> climate_rerun_changed_columns.csv")
    by_sfx = cd[cd.column.isin(ccols)].assign(suffix=[c[len(b):] for c, b in zip(cd[cd.column.isin(ccols)].column, cd[cd.column.isin(ccols)].base)]).groupby("suffix").r2.agg(["min", "median"])
    print(by_sfx.to_string())
    gained_tot = int(cd.loc[cd.column.isin(ccols), "gained"].sum()); lost_tot = int(cd.loc[cd.column.isin(ccols), "lost"].sum())
    print(f"       climate-dependent cells gained={gained_tot} lost={lost_tot}; per-column table -> climate_rerun_changed_columns.csv")
    for fc in flag_cols:
        r = cd[cd.column == fc].iloc[0]
        print(f"       {fc}: changed on {int(r.n_changed)} gages")

    # --- 3. new columns: populated + sniff ranges ---
    pop = []
    for c in added:
        xn = pd.to_numeric(new_i[c], errors="coerce")
        pop.append(dict(column=c, base=base_of(c), finite_frac=float(xn.notna().mean())))
    pop = pd.DataFrame(pop)
    by_suffix = pop.assign(suffix=[c[len(b):] for c, b in zip(pop.column, pop.base)]).groupby("suffix").finite_frac.agg(["min", "median"])
    print(by_suffix.to_string())
    gate("new columns populated (finite fraction >= 0.5 for every field)", bool((pop.finite_frac >= 0.5).all()),
         f"min finite fraction {pop.finite_frac.min():.3f} at {pop.loc[pop.finite_frac.idxmin(), 'column']}")
    sniff_bad = []
    sniff_rows = []
    for b, (lo, hi) in SNIFF.items():
        col = f"{b}_mean"
        if col not in new_i.columns:
            sniff_bad.append((col, "missing")); continue
        v = pd.to_numeric(new_i[col], errors="coerce").dropna()
        q = v.quantile([0.0, 0.01, 0.5, 0.99, 1.0]).to_list()
        out_of = int(((v < lo) | (v > hi)).sum())
        sniff_rows.append(dict(base=b, n=len(v), min=q[0], p01=q[1], median=q[2], p99=q[3], max=q[4], lo=lo, hi=hi, n_outside=out_of))
        if out_of > 0.01 * len(v):
            sniff_bad.append((col, out_of, q[0], q[4]))
    pd.DataFrame(sniff_rows).to_csv(os.path.join(out_dir, "climate_rerun_new_columns_sniff.csv"), index=False)
    gate("new columns sniff ranges (<= 1 % of gages outside the physical range per base)", not sniff_bad,
         f"{len(SNIFF)} bases; offending: {sniff_bad[:6]}")
    # internal identities on means: seasonal sums add up; frost_last <= 273; gini in [0,1]
    seas = new_i[[f"prcp_{s}_mean" for s in ["win", "spr", "sum", "fal"]]].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    ann = pd.to_numeric(new_i["prcp_annual_mean"], errors="coerce")
    ok_id = ((seas - ann).abs() <= 1e-6 * ann.abs() + 1e-9) | seas.isna() | ann.isna()
    gate("identity: seasonal prcp means add up to the annual mean", bool(ok_id.all()), f"violations={int((~ok_id).sum())}")

    # --- 4. annual-level comparison of climate-dependent signatures (optional) ---
    if a.new_annual and a.old_annual:
        na = pd.read_parquet(a.new_annual); oa = pd.read_parquet(a.old_annual)
        csigs = sorted({base_of(c) for c in ccols if base_of(c) in func_of})
        na = na[na.signature.isin(csigs)]; oa = oa[oa.signature.isin(csigs)]
        m = oa.merge(na, on=["gage_id", "signature", "water_year"], how="outer", suffixes=("_old", "_new"), indicator=True)
        arows = []
        for sig, g in m.groupby("signature"):
            both = g[(g._merge == "both") & g.value_old.notna() & g.value_new.notna()]
            d = (both.value_new - both.value_old).abs()
            scale = both.value_old.abs().replace(0, np.nan)
            rel = (d / scale).replace([np.inf], np.nan)
            arows.append(dict(signature=sig, shared_years=len(both), max_abs_diff=float(d.max()) if len(d) else np.nan,
                              p99_rel_diff=float(rel.quantile(0.99)) if rel.notna().any() else np.nan,
                              n_changed_gt_1e9=int((d > 1e-9).sum()),
                              added_gage_years=int(((g._merge == "right_only") | (g.value_old.isna() & g.value_new.notna())).sum()),
                              removed_gage_years=int(((g._merge == "left_only") | (g.value_new.isna() & g.value_old.notna())).sum())))
        ad = pd.DataFrame(arows).sort_values("signature")
        # identity R² on shared gage-years, per signature
        r2s = {}
        for sig, g in m.groupby("signature"):
            both = g[(g._merge == "both") & g.value_old.notna() & g.value_new.notna()]
            r2s[sig] = float(np.corrcoef(both.value_old, both.value_new)[0, 1] ** 2) if len(both) >= 3 and both.value_old.std() > 0 else np.nan
        ad["r2_shared_years"] = ad.signature.map(r2s)
        ad.to_csv(os.path.join(out_dir, "climate_rerun_annual_shared_years.csv"), index=False)
        print(ad.to_string(index=False))
        print("       (shared gage-years should agree at the Daymet-agreement level; added gage-years are WY 2024-2025 and the gages that gained climate)")
        DISCRETE = {"swe_max_dowy", "snow_on_dowy", "snow_off_dowy", "melt_season_days", "melt_com_dowy",
                    "snow_cover_days", "ssm"}
        UNSTABLE = {"elasticity_annual", "elasticity_rolling"}
        bad_cont = []
        for _, r in ad.iterrows():
            if r.signature in DISCRETE or r.signature in UNSTABLE or r.shared_years < 30:
                continue
            if (not np.isnan(r.r2_shared_years) and r.r2_shared_years < a.r2_min) or (not np.isnan(r.p99_rel_diff) and r.p99_rel_diff > a.rel_p99):
                bad_cont.append((r.signature, round(r.r2_shared_years, 6), r.p99_rel_diff))
        gate(f"2a annual values on shared gage-years: continuous signatures R² >= {a.r2_min} and p99 rel diff <= {a.rel_p99}",
             not bad_cont, f"{len(ad)} signatures; failing: {bad_cont}")
        disc = ad[ad.signature.isin(DISCRETE)]
        print("       discrete signatures, fraction of shared gage-years changed: " +
              ", ".join(f"{r.signature} {r.n_changed_gt_1e9 / max(r.shared_years, 1):.2%}" for _, r in disc.iterrows()))
        # 2b: the old _mean recomputed from the NEW annual values on the OLD year set
        old_years = oa[oa.value.notna()][["gage_id", "signature", "water_year"]]
        new_on_old = na.merge(old_years, on=["gage_id", "signature", "water_year"], how="inner")
        mean_new = new_on_old.groupby(["gage_id", "signature"]).value.mean().unstack()
        bad_mean = []; mrows = []
        for sig in mean_new.columns:
            col = f"{sig}_mean"
            if col not in old_i.columns:
                continue
            x_old = pd.to_numeric(old_i[col], errors="coerce").astype(float)
            x_new = mean_new[sig].reindex(old_i.index)
            both = x_old.notna() & x_new.notna()
            n_b = int(both.sum())
            r2 = float(np.corrcoef(x_old[both], x_new[both])[0, 1] ** 2) if n_b >= 3 and x_old[both].std() > 0 else np.nan
            mad = float((x_old[both] - x_new[both]).abs().max()) if n_b else np.nan
            mrows.append(dict(signature=sig, n_gages=n_b, r2_mean_on_old_years=r2, max_abs_diff=mad))
            if sig not in UNSTABLE and n_b >= 30 and not np.isnan(r2) and r2 < a.r2_min:
                bad_mean.append((sig, round(r2, 6), mad))
        md = pd.DataFrame(mrows).sort_values("signature")
        md.to_csv(os.path.join(out_dir, "climate_rerun_means_on_old_years.csv"), index=False)
        print(md.to_string(index=False))
        gate(f"2b old-product _mean reproduced from the NEW annual values on the OLD year set: R² >= {a.r2_min}",
             not bad_mean, f"{len(md)} signatures; failing: {bad_mean}")

    summary = {"gates": [{"name": n, "pass": bool(ok), "detail": d} for n, ok, d in results],
               "n_new_gages": len(sn), "n_old_gages": len(so), "added_columns": len(added), "dropped_columns": len(dropped),
               "column_classes": counts}
    with open(os.path.join(out_dir, "climate_rerun_review_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    ok = all(r[1] for r in results)
    print("\nCLIMATE RERUN REVIEW:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
