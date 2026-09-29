"""Compare newly aggregated Daymet basin series with the stale co-author product.

The stale input (`daymet_1980_2023_rebuilt_10aug2026.parquet`, = HydroShare Resource 3
`hisss_daymet_basin_daily.parquet`) was aggregated by the co-authors with gdptools on
their own polygon set. Agreement here is agreement with THAT run, not proof of
correctness (daymet_crosscheck.py checks the aggregation itself).

Usage:
    python daymet_validate.py --new-dir <dir with <var>_<year>.parquet> --old <stale parquet>
        --basins <daymet_basins.csv> --years 2023 --out-dir <dir> [--vars prcp,swe,...]

Joins on the canonical (zero-stripped) id and Date. Per basin and variable: n days,
old/new mean, bias (new - old), MAE, RMSE, max |diff|, Pearson r, identity R2
(1 - SSE/SST about the 1:1 line; NaN when the old series is constant), ratio of annual
sums, and the correlation at lags -1/0/+1 day (a date shift shows up as lag != 0).

Outputs in --out-dir:
    validation_<years>_per_basin.csv   one row per (basin, variable, year) with the QA columns
    validation_<years>_summary.json    distributions and counts per variable
"""
import argparse
import json
import os
import sys

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import DAYMET_VARS, canon, write_json  # noqa: E402

QS = [0.0, 0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1.0]


def _lagcorr(a, b, k):
    if k > 0:
        a, b = a[k:], b[:-k]
    elif k < 0:
        a, b = a[:k], b[-k:]
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def per_basin(m, var):
    rows = []
    for c, x in m.groupby("canon_id", sort=False):
        x = x.sort_values("Date")
        o, n = x[f"{var}_old"].to_numpy(), x[f"{var}_new"].to_numpy()
        ok = np.isfinite(o) & np.isfinite(n)
        rec = {"canon_id": c, "var": var, "n_days": len(x), "n_both_finite": int(ok.sum()),
               "n_old_nan": int(np.isnan(o).sum()), "n_new_nan": int(np.isnan(n).sum())}
        if ok.sum() >= 2:
            o, n = o[ok], n[ok]
            d = n - o
            sst = float(((o - o.mean()) ** 2).sum())
            rec.update({
                "old_mean": o.mean(), "new_mean": n.mean(), "bias": d.mean(), "mae": np.abs(d).mean(),
                "rmse": float(np.sqrt((d ** 2).mean())), "max_abs_diff": np.abs(d).max(),
                "r": _lagcorr(n, o, 0), "r_lag_m1": _lagcorr(n, o, -1), "r_lag_p1": _lagcorr(n, o, 1),
                "r2_identity": 1 - float((d ** 2).sum()) / sst if sst > 0 else np.nan,
                "sum_ratio": n.sum() / o.sum() if o.sum() != 0 else np.nan,
                "both_constant_equal": bool(sst == 0 and np.allclose(n, o, atol=1e-9)),
            })
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize(df, var):
    s = df[df["var"] == var]
    fin = s[s["n_both_finite"] > 0]
    r2 = fin["r2_identity"]
    out = {
        "basins_compared": int(len(s)), "basins_old_all_nan": int((s["n_old_nan"] == s["n_days"]).sum()),
        "basins_new_any_nan": int((s["n_new_nan"] > 0).sum()),
        "basins_constant_old_series": int(r2.isna().sum()),
        "of_which_new_equal": int(fin["both_constant_equal"].fillna(False).sum()),
        "r2_identity_quantiles": {str(q): round(float(v), 6) for q, v in r2.quantile(QS).items()},
        "share_r2_ge_0.999": round(float((r2 >= 0.999).sum() / max(1, r2.notna().sum())), 5),
        "share_r2_ge_0.99": round(float((r2 >= 0.99).sum() / max(1, r2.notna().sum())), 5),
        "mae_quantiles": {str(q): float(v) for q, v in fin["mae"].quantile(QS).items()},
        "max_abs_diff_quantiles": {str(q): float(v) for q, v in fin["max_abs_diff"].quantile(QS).items()},
        "bias_quantiles": {str(q): float(v) for q, v in fin["bias"].quantile(QS).items()},
        "sum_ratio_quantiles": {str(q): round(float(v), 6) for q, v in fin["sum_ratio"].quantile(QS).items()},
        "median_r_lag0": float(fin["r"].median()), "median_r_lag_m1": float(fin["r_lag_m1"].median()),
        "median_r_lag_p1": float(fin["r_lag_p1"].median()),
        "basins_lag_beats_lag0": int(((fin["r_lag_m1"] > fin["r"]) | (fin["r_lag_p1"] > fin["r"])).sum()),
    }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new-dir", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--basins", required=True)
    ap.add_argument("--years", required=True, help="e.g. 2023 or 1980-2023")
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    y = [int(v) for v in a.years.split("-")]
    years = list(range(y[0], y[-1] + 1))
    variables = [v for v in a.vars.split(",") if v]
    basins = pd.read_csv(a.basins, dtype={"site_id": str, "canon_id": str, "boundary_gage_id": str})

    all_rows, summary = [], {"years": years, "old": a.old, "new_dir": a.new_dir, "vars": {}}
    for year in years:
        cols = ", ".join(["site_id", "Date"] + variables)
        old = duckdb.sql(f"SELECT {cols} FROM read_parquet('{a.old}') WHERE Date >= DATE '{year}-01-01' "
                         f"AND Date <= DATE '{year}-12-31'").df()
        old["canon_id"] = old["site_id"].map(canon)
        old["Date"] = pd.to_datetime(old["Date"])
        for var in variables:
            p = os.path.join(a.new_dir, f"{var}_{year}.parquet")
            if not os.path.exists(p):
                print(f"missing {p}; skipped")
                continue
            new = pd.read_parquet(p)
            new["canon_id"] = new["site_id"].map(canon)
            new["Date"] = pd.to_datetime(new["Date"])
            m = new[["canon_id", "Date", var]].merge(old[["canon_id", "Date", var]], on=["canon_id", "Date"],
                                                     suffixes=("_new", "_old"), how="inner")
            join = {"old_sites": int(old["canon_id"].nunique()), "new_sites": int(new["canon_id"].nunique()),
                    "shared_sites": int(m["canon_id"].nunique()), "joined_rows": int(len(m)),
                    "old_rows": int(len(old)), "new_rows": int(len(new)),
                    "shared_sites_x_days": int(m["canon_id"].nunique()) * 365}
            pb = per_basin(m, var)
            pb["year"] = year
            all_rows.append(pb)
            summary["vars"].setdefault(var, {})[str(year)] = {"join": join, **summarize(pb, var)}
            s = summary["vars"][var][str(year)]
            print(f"{var} {year}: {join['shared_sites']} shared basins ({join['joined_rows']:,} rows); "
                  f"R2 min {s['r2_identity_quantiles']['0.0']}, p01 {s['r2_identity_quantiles']['0.01']}, "
                  f"median {s['r2_identity_quantiles']['0.5']}; share >= 0.999: {s['share_r2_ge_0.999']}; "
                  f"sum ratio p01-p99 {s['sum_ratio_quantiles']['0.01']}-{s['sum_ratio_quantiles']['0.99']}; "
                  f"old all-NaN {s['basins_old_all_nan']}, constant old {s['basins_constant_old_series']} "
                  f"({s['of_which_new_equal']} equal); lag beats lag0 in {s['basins_lag_beats_lag0']}", flush=True)
    df = pd.concat(all_rows, ignore_index=True)
    keep = ["canon_id", "site_id", "watershed_geom_source", "gage_type", "geom_area_km2", "geom_simplified",
            "low_confidence", "n_cells", "coverage_sum", "low_pixel_support"]
    df = df.merge(basins[keep], on="canon_id", how="left")
    tag = a.years.replace("-", "_")
    df.to_csv(os.path.join(a.out_dir, f"validation_{tag}_per_basin.csv"), index=False)
    write_json(os.path.join(a.out_dir, f"validation_{tag}_summary.json"), summary)
    print(f"wrote validation_{tag}_per_basin.csv ({len(df):,} rows) and validation_{tag}_summary.json")


if __name__ == "__main__":
    main()
