"""Georeferencing check against ORNL's Daymet Single Pixel Extraction API.

For a few lat/lon points, fetches the daily series ORNL itself serves for the pixel
containing the point (https://daymet.ornl.gov/single-pixel/api) and compares it with
the value this toolchain reads at the same point: lat/lon -> Daymet LCC (x, y) -> row,
col from the file's own grid. Agreement shows the CRS, the grid origin and the row/col
arithmetic are right -- something daymet_crosscheck.py cannot show, because it shares
that georeferencing with daymet_aggregate.py. Also useful for years with no stale
product to compare against (2024 onward).

Usage:
    python daymet_pixelcheck.py --file <daymet_v4_daily_na_<var>_<yyyy>.nc>
        [--points "35.96,-84.29;47.60,-120.66;..."] [--out <json>]
"""
import argparse
import io
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import requests
from pyproj import Transformer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import FILL, read_grid, write_json  # noqa: E402

API = "https://daymet.ornl.gov/single-pixel/api/data"
# Oak Ridge TN, Cascades WA, Colorado Rockies, Quebec, Florida, Yukon (default spread)
DEFAULT_POINTS = "35.96,-84.29;47.60,-120.66;39.60,-106.00;47.00,-72.50;29.00,-81.50;61.00,-135.00"
API_COL = {"prcp": "prcp (mm/day)", "swe": "swe (kg/m^2)", "tmin": "tmin (deg c)", "tmax": "tmax (deg c)",
           "vp": "vp (Pa)", "srad": "srad (W/m^2)", "dayl": "dayl (s)"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True)
    ap.add_argument("--points", default=DEFAULT_POINTS)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    grid = read_grid(a.file)
    var, year = grid.var, int(str(grid.dates[0])[:4])
    xmin, _, _, ymax = grid.bounds
    to_lcc = Transformer.from_crs("EPSG:4326", grid.crs, always_xy=True)
    res = {"file": os.path.basename(a.file), "var": var, "year": year, "points": []}
    with h5py.File(a.file, "r") as f:
        d = f[var]
        for p in a.points.split(";"):
            lat, lon = (float(v) for v in p.split(","))
            x, y = to_lcc.transform(lon, lat)
            col, row = int(np.floor((x - xmin) / 1000.0)), int(np.floor((ymax - y) / 1000.0))
            ours = d[:, row, col].astype(np.float64)
            txt = None
            for attempt in range(4):                 # the API times out now and then
                try:
                    r = requests.get(API, params={"lat": lat, "lon": lon, "vars": var, "start": f"{year}-01-01",
                                                  "end": f"{year}-12-31"}, timeout=120)
                    r.raise_for_status()
                    txt = r.text
                    break
                except requests.RequestException as e:
                    print(f"  API attempt {attempt + 1} failed for ({lat}, {lon}): {type(e).__name__}", flush=True)
                    time.sleep(15 * (attempt + 1))
            if txt is None:
                res["points"].append({"lat": lat, "lon": lon, "row": row, "col": col, "api_error": True})
                continue
            xy = txt.split("X & Y on Lambert Conformal Conic:")[1].split("\n")[0].split()
            body = txt[txt.index("year,yday"):]
            api = pd.read_csv(io.StringIO(body))
            theirs = api[API_COL[var]].to_numpy(dtype=np.float64)
            rec = {"lat": lat, "lon": lon, "row": row, "col": col, "n_api": len(theirs), "n_file": len(ours),
                   "file_fill_days": int((ours == FILL).sum()),
                   "api_x": float(xy[0]), "api_y": float(xy[1]), "dx_m": x - float(xy[0]), "dy_m": y - float(xy[1])}
            if len(theirs) == len(ours) and rec["file_fill_days"] == 0:
                diff = np.abs(ours - theirs)
                rec.update({"max_abs_diff": float(diff.max()), "mean_abs_diff": float(diff.mean()),
                            "api_mean": float(theirs.mean()),
                            "neighbour_max_abs_diff": float(np.abs(d[:, row, col + 1] - theirs).max())})
            res["points"].append(rec)
            print(f"{var} {year} ({lat}, {lon}) -> row {row} col {col}, LCC offset vs API "
                  f"({rec['dx_m']:.2f}, {rec['dy_m']:.2f}) m: " +
                  (f"max |file - API| {rec['max_abs_diff']:.3g} (API mean {rec['api_mean']:.3g}; "
                   f"east neighbour differs by up to {rec['neighbour_max_abs_diff']:.3g})"
                   if "max_abs_diff" in rec else f"not compared ({rec})"), flush=True)
            time.sleep(1)
    ok = [p for p in res["points"] if "max_abs_diff" in p]
    res["verdict"] = "PASS" if ok and all(p["max_abs_diff"] < 0.01 and abs(p["dx_m"]) < 1 and abs(p["dy_m"]) < 1
                                          for p in ok) else "CHECK"
    n_err = sum(1 for p in res["points"] if p.get("api_error"))
    print(f"pixel check {var} {year}: {res['verdict']} ({len(ok)} points compared"
          f"{f', {n_err} skipped after API errors' if n_err else ''})")
    if a.out:
        write_json(a.out, res)


if __name__ == "__main__":
    main()
