"""Check the assembled basin series against ORNL's Daymet Single Pixel API (no raw files).

For a few small basins, fetch the daily series of every weighted cell from ORNL's own
service (https://daymet.ornl.gov/single-pixel/api), form the run's weighted mean
sum(w * v) / sum(w) with the run's weights, and compare it with the basin's series in
the assembled parquet. The cell values do not pass through this toolchain, so the check
covers what the stale-product comparison cannot:
  * years after 2023;
  * basins the stale product lacks;
  * the assembled file itself.
It also catches a date shift or a misplaced cell: the query goes to each cell's centre,
and ORNL's projected coordinates for that point must equal the centre to < 1 m.
It shares the polygons and the weights with the run. The raw files are never read.

Usage:
    python daymet_outputcheck.py --parquet <daymet_....parquet> --weights-dir <dir>
        [--years 1980,2019,2020,2024,2025] [--sites 02BF013,12181200 | --auto 6]
        [--max-cells 6] [--vars prcp,tmin,tmax,swe,vp,srad] [--tol 1e-3] [--out <json>]
        [--cache-dir <dir>]
--auto picks basins of 2..--max-cells cells, half USGS and half Canadian, spread evenly
over the sorted ids. A cell the API serves no data for (fill) leaves its basin
"not compared".
"""
import argparse
import io
import json
import os
import sys
import time

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import requests
from pyproj import CRS, Transformer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import CELL, DAYMET_VARS, NX, X0_CENTRE, Y0_CENTRE, utc_now, write_json  # noqa: E402
from daymet_pixelcheck import API, API_COL  # noqa: E402


def api_series(lat, lon, years, variables, cache_dir):
    """ORNL's daily series at (lat, lon): (DataFrame year, yday, <vars>), (x, y) on LCC."""
    key = os.path.join(cache_dir, f"{lat:.7f}_{lon:.7f}_{'-'.join(map(str, years))}_{'-'.join(variables)}.csv")
    txt = open(key).read() if os.path.exists(key) else None
    for attempt in range(4):
        if txt is not None:
            break
        try:
            r = requests.get(API, params={"lat": f"{lat:.7f}", "lon": f"{lon:.7f}", "vars": ",".join(variables),
                                          "years": ",".join(map(str, years))}, timeout=90)
            r.raise_for_status()
            txt = r.text
            open(key, "w").write(txt)
        except requests.RequestException as e:
            print(f"  API attempt {attempt + 1} failed at ({lat:.4f}, {lon:.4f}): {type(e).__name__}", flush=True)
            time.sleep(15 * (attempt + 1))
        time.sleep(0.5)
    if txt is None or "year,yday" not in txt:
        return None, None
    xy = [float(v) for v in txt.split("X & Y on Lambert Conformal Conic:")[1].split("\n")[0].split()[:2]]
    df = pd.read_csv(io.StringIO(txt[txt.index("year,yday"):]))
    return df.rename(columns={API_COL[v]: v for v in variables}), xy


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--years", default="1980,2019,2020,2024,2025")
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    ap.add_argument("--sites", default=None)
    ap.add_argument("--auto", type=int, default=6)
    ap.add_argument("--max-cells", type=int, default=6)
    ap.add_argument("--weight", choices=["area", "coverage"], default="area")
    ap.add_argument("--tol", type=float, default=1e-3, help="max |file - API mean| per day (float32 rounding ~1e-4)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cache-dir", default=None)
    a = ap.parse_args()
    years = [int(y) for y in a.years.split(",")]
    variables = [v for v in a.vars.split(",") if v]
    cache = a.cache_dir or ((a.out or "outputcheck") + ".cache")
    os.makedirs(cache, exist_ok=True)

    basins = pd.read_csv(os.path.join(a.weights_dir, "daymet_basins.csv"), dtype={"site_id": str})
    if a.sites:
        chosen = basins[basins["site_id"].isin(a.sites.split(","))]
    else:
        cand = basins[basins["n_cells"].between(2, a.max_cells)].sort_values("site_id")
        us, ca = cand[cand["gage_type"] == "USGS"], cand[cand["gage_type"] != "USGS"]
        pick = []
        for sub, k in ((us, a.auto // 2), (ca, a.auto - a.auto // 2)):
            if k and len(sub):
                pick.append(sub.iloc[np.unique(np.linspace(0, len(sub) - 1, k).round().astype(int))])
        chosen = pd.concat(pick)
    if chosen.empty:
        sys.exit("no basins selected")
    wcol = {"area": "area_km2", "coverage": "coverage"}[a.weight]
    w = pq.read_table(os.path.join(a.weights_dir, "daymet_weights.parquet"),
                      filters=[("basin_idx", "in", chosen["basin_idx"].tolist())],
                      columns=["basin_idx", "cell_id", wcol]).to_pandas()
    meta = json.load(open(os.path.join(a.weights_dir, "daymet_weights_meta.json")))
    to_ll = Transformer.from_crs(CRS.from_wkt(meta["crs_wkt"]), "EPSG:4326", always_xy=True)

    con = duckdb.connect()
    con.execute("SET memory_limit='3GB'; SET threads=4")
    ids = ",".join(f"'{s}'" for s in chosen["site_id"])
    out = con.execute(f"SELECT site_id, Date, {', '.join(variables)} FROM read_parquet('{a.parquet}') "
                      f"WHERE site_id IN ({ids}) AND year(Date) IN ({', '.join(map(str, years))})").df()
    out["year"] = pd.to_datetime(out["Date"]).dt.year
    out["yday"] = (pd.to_datetime(out["Date"]) - pd.to_datetime(out["year"].astype(str) + "-01-01")).dt.days + 1

    res = {"parquet": os.path.abspath(a.parquet), "weights_md5": meta.get("weights_md5"), "years": years,
           "vars": variables, "weight": a.weight, "tol": a.tol, "checked_utc": utc_now(), "sites": []}
    for _, b in chosen.iterrows():
        cells = w[w["basin_idx"] == b["basin_idx"]]
        rec = {"site_id": b["site_id"], "gage_type": b["gage_type"], "n_cells": int(len(cells)),
               "geom_area_km2": float(b["geom_area_km2"])}
        num, den, xy_off, ok = None, 0.0, 0.0, True
        for cid, wt in zip(cells["cell_id"].to_numpy(), cells[wcol].to_numpy()):
            x, y = X0_CENTRE + (cid % NX) * CELL, Y0_CENTRE - (cid // NX) * CELL
            lon, lat = to_ll.transform(x, y)
            df, axy = api_series(lat, lon, years, variables, cache)
            if df is None or (df[variables] <= -9999).any().any():
                ok = False
                break
            xy_off = max(xy_off, abs(axy[0] - x), abs(axy[1] - y))
            df = df.set_index(["year", "yday"])[variables]
            num = df * wt if num is None else num + df * wt
            den += wt
        if not ok:
            rec["verdict"] = "not compared (no API data for a cell)"
            res["sites"].append(rec)
            continue
        api = (num / den).reset_index()
        m = out[out["site_id"] == b["site_id"]].merge(api, on=["year", "yday"], suffixes=("", "_api"))
        rec["days"] = int(len(m))
        rec["api_xy_offset_m"] = round(xy_off, 4)
        rec["max_abs_diff"], rec["max_abs_diff_shift1"] = {}, {}
        for v in variables:
            d = (m[v] - m[v + "_api"]).abs()
            rec["max_abs_diff"][v] = float(d.max())
            g = m.sort_values(["year", "yday"])
            rec["max_abs_diff_shift1"][v] = float((g[v].shift(1) - g[v + "_api"]).abs().max())   # what a 1-day shift would show
        rec["verdict"] = "PASS" if (rec["days"] == 365 * len(years) and xy_off < 1.0
                                    and all(x <= a.tol for x in rec["max_abs_diff"].values())) else "FAIL"
        res["sites"].append(rec)
        print(f"{b['site_id']} ({b['gage_type']}, {len(cells)} cells, {b['geom_area_km2']:.2f} km2): {rec['verdict']}; "
              f"days {rec['days']}, API xy offset {xy_off:.3f} m, max |file - API| "
              + ", ".join(f"{v} {rec['max_abs_diff'][v]:.2g}" for v in variables), flush=True)
    done = [s for s in res["sites"] if s["verdict"] in ("PASS", "FAIL")]
    res["verdict"] = "PASS" if done and all(s["verdict"] == "PASS" for s in done) else "FAIL"
    print(f"output check: {res['verdict']} ({len(done)} basins compared, "
          f"{len(res['sites']) - len(done)} not compared; years {years})")
    if a.out:
        write_json(a.out, res)
    sys.exit(0 if res["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
