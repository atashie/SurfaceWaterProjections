"""Independent check of daymet_aggregate.py on real data.

For a few days of one variable-year file, reads each full-day slice, hands it to
exactextract as a raster with nodata = -9999, and lets exactextract compute the basin
mean itself for every polygon: "mean" (coverage-weighted) for --weight coverage, or
"weighted_mean" with a raster of true cell areas for --weight area. That result shares
no code with the chunk-aligned sparse accumulation in daymet_aggregate.py, so agreement
to rounding error shows the tiling, the time blocking, the fill handling, the weights
and the row/column bookkeeping are right, for every basin including coastal and
tile-straddling ones.

Usage:
    python daymet_crosscheck.py --file <nc> --polygons <boundaries> --weights-dir <dir>
        --agg-dir <dir with <var>_<year>.parquet> [--weight area|coverage]
        [--days 0,59,181,272,364] [--out <json>]
"""
import argparse
import os
import sys
import time

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import FILL, NX, NY, cell_true_area_km2, read_grid, write_json  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True)
    ap.add_argument("--polygons", required=True)
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--agg-dir", required=True)
    ap.add_argument("--weight", choices=["area", "coverage"], default="area")
    ap.add_argument("--days", default="0,59,181,272,364", help="0-based day indices")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    grid = read_grid(a.file)
    var, year = grid.var, int(str(grid.dates[0])[:4])
    basins = pd.read_csv(os.path.join(a.weights_dir, "daymet_basins.csv"),
                         dtype={"site_id": str, "boundary_gage_id": str})
    g = gpd.read_parquet(a.polygons) if a.polygons.endswith(".parquet") else gpd.read_file(a.polygons)
    g = g.reset_index(drop=True)
    if len(g) != len(basins) or not (g["gage_id"].to_numpy() == basins["boundary_gage_id"].to_numpy()).all():
        sys.exit("polygon layer order differs from the weights' basin table")
    gl = g[["geometry"]].to_crs(grid.crs)
    agg = pd.read_parquet(os.path.join(a.agg_dir, f"{var}_{year}.parquet"))
    agg["Date"] = pd.to_datetime(agg["Date"])
    xmin, ymin, xmax, ymax = grid.bounds
    wkt = grid.crs.to_wkt()
    wrast = None
    if a.weight == "area":
        t0 = time.time()
        area = cell_true_area_km2(np.arange(NX * NY, dtype=np.int64), grid.crs).reshape(NY, NX)
        wrast = NumPyRasterSource(area, xmin, ymin, xmax, ymax, srs_wkt=wkt, name="cell_area")
        print(f"true-area raster for {NX * NY:,} cells in {time.time() - t0:.0f} s")
    days = [int(d) for d in a.days.split(",")]
    res = {"file": os.path.basename(a.file), "var": var, "year": year, "weight": a.weight, "days": {}}
    with h5py.File(a.file, "r") as f:
        d = f[var]
        for k in days:
            t0 = time.time()
            rast = NumPyRasterSource(d[k, :, :], xmin, ymin, xmax, ymax, nodata=FILL, srs_wkt=wkt, name="v")
            if a.weight == "area":
                ee = exact_extract(rast, gl, ["weighted_mean"], weights=wrast, output="pandas")
            else:
                ee = exact_extract(rast, gl, ["mean"], output="pandas")
            ee = ee.iloc[:, -1].to_numpy(dtype=np.float64)
            date = pd.Timestamp(grid.dates[k])
            ours = (agg[agg["Date"] == date].set_index("site_id").loc[basins["site_id"], var]
                    .to_numpy(dtype=np.float64))
            ok = np.isfinite(ee) & np.isfinite(ours)
            diff = np.abs(ee[ok] - ours[ok])
            big = np.abs(ee[ok]) > 1e-6
            r = {"date": str(date.date()), "basins": len(ee), "compared": int(ok.sum()),
                 "both_nan": int((np.isnan(ee) & np.isnan(ours)).sum()),
                 "nan_mismatch": int((np.isnan(ee) ^ np.isnan(ours)).sum()),
                 "max_abs_diff": float(diff.max()) if diff.size else None,
                 "max_rel_diff": float((diff[big] / np.abs(ee[ok][big])).max()) if big.any() else 0.0,
                 "seconds": round(time.time() - t0, 1)}
            res["days"][str(k)] = r
            print(f"{var} {r['date']}: {r['compared']} basins, max |diff| {r['max_abs_diff']:.3e}, "
                  f"max rel {r['max_rel_diff']:.3e}, NaN mismatches {r['nan_mismatch']} ({r['seconds']} s)",
                  flush=True)
    worst = max(v["max_abs_diff"] or 0 for v in res["days"].values())
    res["verdict"] = "PASS" if worst < 1e-6 and all(v["nan_mismatch"] == 0 for v in res["days"].values()) else "FAIL"
    print(f"cross-check {var} {year} ({a.weight} weights): {res['verdict']} (worst max |diff| {worst:.3e})")
    if a.out:
        write_json(a.out, res)
    if res["verdict"] != "PASS":
        sys.exit(1)


if __name__ == "__main__":
    main()
