"""Build a self-contained HTML dashboard: fresh Daymet basin series vs the original ones.

For one calendar year, compares the recomputed basin means (EO_data_processing/daymet) with
the stale co-author product: summary tiles, an original-vs-fresh scatter of every shared
basin's annual mean, a map of the difference, daily series (overlay + difference) for a
curated set of basins, and the basins that agree least. One or more "configurations"
(e.g. full-resolution vs published polygons) can be embedded and switched. Plotly loads
from cdnjs; all data is inline, so the page is a single file.

Usage:
    python build_daymet_comparison_dashboard.py --year 2023 --old <stale parquet>
        --polygons <polygon parquet, same basin order as the first config's weights>
        --config "fr=Full-resolution polygons,<out dir>,<weights dir>,<validation dir>"
        [--config "pub=Published Resource 4 polygons,<out dir>,<weights dir>,<validation dir>"]
        --out <page.html> [--fragment <body-only.html>] [--assets-dir <cache dir>]

<out dir> holds <var>_<year>.parquet; <weights dir> holds daymet_basins.csv; any
validation_*_per_basin.csv under <validation dir> (from daymet_validate.py) is read and
filtered to the year. Every configuration must list the same basins in the same order.
"""
import argparse
import base64
import glob
import json
import os
import sys
import time
import urllib.request

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import LineString, MultiLineString

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "daymet"))
from daymet_common import DAYMET_VARS, canon  # noqa: E402

LCC = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=m"
NAMES = {"prcp": "Precipitation", "swe": "Snow water equivalent", "tmin": "Minimum temperature",
         "tmax": "Maximum temperature", "vp": "Vapour pressure", "srad": "Shortwave radiation"}
SHORT = {"prcp": "Precipitation", "swe": "SWE", "tmin": "Tmin", "tmax": "Tmax", "vp": "Vapour pressure",
         "srad": "Shortwave"}
UNITS = {"prcp": "mm/day", "swe": "mm", "tmin": "°C", "tmax": "°C", "vp": "Pa", "srad": "W/m²"}
NE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
WORST_N = 10


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype="<f4").tobytes()).decode("ascii")


def fetch(url, cache):
    if not os.path.exists(cache):
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        urllib.request.urlretrieve(url, cache + ".part")
        os.replace(cache + ".part", cache)
    return cache


def outlines(cache_dir, to_lcc):
    """Country outlines (USA, CAN, MEX) + US state / Canadian province lines, in km."""
    c = gpd.read_file(fetch(NE + "ne_50m_admin_0_countries.geojson",
                            os.path.join(cache_dir, "ne_50m_admin_0_countries.geojson")))
    c = c[c["ADM0_A3"].isin(["USA", "CAN", "MEX"])]
    l1 = gpd.read_file(fetch(NE + "ne_50m_admin_1_states_provinces_lines.geojson",
                             os.path.join(cache_dir, "ne_50m_admin_1_states_provinces_lines.geojson")))
    col = "ADM0_A3" if "ADM0_A3" in l1.columns else "adm0_a3"
    l1 = l1[l1[col].isin(["USA", "CAN"])]
    xs, ys = [], []
    for geom in list(c.geometry.boundary) + list(l1.geometry):
        parts = geom.geoms if isinstance(geom, MultiLineString) else [geom]
        for ln in parts:
            if not isinstance(ln, LineString) or ln.length == 0:
                continue
            lon, lat = np.asarray(ln.coords).T[:2]
            if lat.max() < 14 or lon.min() > -50:            # far-off islands
                continue
            x, y = to_lcc.transform(lon, lat)
            s = LineString(np.c_[x, y] / 1000.0).simplify(4.0)   # km; 4 km tolerance
            cx, cy = np.asarray(s.coords).T
            xs += np.round(cx, 1).tolist() + [None]
            ys += np.round(cy, 1).tolist() + [None]
    return {"x": xs, "y": ys}


def parse_config(s):
    key, rest = s.split("=", 1)
    label, out_dir, weights_dir, val_dir = [p.strip() for p in rest.split(",")]
    return {"key": key.strip(), "label": label, "out": os.path.expanduser(out_dir),
            "weights": os.path.expanduser(weights_dir), "val": os.path.expanduser(val_dir)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--polygons", required=True)
    ap.add_argument("--config", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fragment", default=None)
    ap.add_argument("--assets-dir", default=os.path.expanduser("~/HISSS_data/dashboard_assets"))
    a = ap.parse_args()
    t0 = time.time()
    Y = a.year
    cfgs = [parse_config(s) for s in a.config]

    # ---- basins: first configuration defines the order --------------------------------
    tables = [pd.read_csv(os.path.join(c["weights"], "daymet_basins.csv"),
                          dtype={"site_id": str, "boundary_gage_id": str, "canon_id": str}) for c in cfgs]
    bf = tables[0]
    for c, t in zip(cfgs[1:], tables[1:]):
        if not t["canon_id"].equals(bf["canon_id"]):
            sys.exit(f"configuration {c['key']} lists basins in a different order")
    g = gpd.read_parquet(a.polygons)
    if not (g["gage_id"].map(canon).to_numpy() == bf["canon_id"].to_numpy()).all():
        sys.exit("polygon layer order differs from the basin table")
    lat = g["latitude"].to_numpy(dtype=float)
    lon = g["longitude"].to_numpy(dtype=float)
    miss = ~(np.isfinite(lat) & np.isfinite(lon))
    if miss.any():
        rp = g.geometry[miss].representative_point()
        lon[miss], lat[miss] = rp.x.to_numpy(), rp.y.to_numpy()
    to_lcc = Transformer.from_crs("EPSG:4326", LCC, always_xy=True)
    bx, by = to_lcc.transform(lon, lat)
    nb = len(bf)
    idx_of = {c: i for i, c in enumerate(bf["canon_id"])}
    simp = g["geom_simplified"].to_numpy(dtype=bool) if "geom_simplified" in g else np.zeros(nb, bool)

    # ---- per-basin validation ----------------------------------------------------------
    val = {}
    for c in cfgs:
        files = glob.glob(os.path.join(c["val"], "**", "validation_*_per_basin.csv"), recursive=True)
        d = pd.concat([pd.read_csv(f, dtype={"canon_id": str}) for f in files], ignore_index=True)
        d = d[d["year"] == Y].drop_duplicates(["canon_id", "var", "year"], keep="last")
        for v in DAYMET_VARS:
            dv = d[d["var"] == v].set_index("canon_id")
            if dv.empty:
                sys.exit(f"no {v} {Y} validation rows under {c['val']}")
            val[(c["key"], v)] = dv
    k0 = cfgs[0]["key"]
    shared = [c for c in val[(k0, "prcp")].index if c in idx_of]
    shared_idx = np.array([idx_of[c] for c in shared])
    old_nan = set(val[(k0, "prcp")].index[val[(k0, "prcp")]["n_both_finite"] == 0])
    status = np.full(nb, 2, dtype=int)                  # 2 = no original series (new coverage)
    status[shared_idx] = 0                               # 0 = compared
    for c in old_nan:
        status[idx_of[c]] = 1                            # 1 = original all NaN

    summary = {}
    for v in DAYMET_VARS:
        s = {"old_mean": val[(k0, v)].loc[shared, "old_mean"].to_numpy(dtype=float)}
        for c in cfgs:
            d = val[(c["key"], v)].reindex(shared)
            s[f"dmean_{c['key']}"] = d["bias"].to_numpy(dtype=float)          # mean(fresh - original)
            s[f"omr2_{c['key']}"] = 1.0 - d["r2_identity"].to_numpy(dtype=float)
            s[f"max_{c['key']}"] = d["max_abs_diff"].to_numpy(dtype=float)
        summary[v] = {name: b64(arr) for name, arr in s.items()}

    # ---- curated basins with daily series ----------------------------------------------
    reasons = {}

    def add(c, why):
        reasons.setdefault(c, [])
        if why not in reasons[c]:
            reasons[c].append(why)

    for v in DAYMET_VARS:
        for c in cfgs:
            d = val[(c["key"], v)]
            d = d[d["n_both_finite"] > 0].sort_values("r2_identity")
            for cid in d.index[:WORST_N]:
                add(cid, f"worst:{v}:{c['key']}")
    for cid in sorted(old_nan):
        add(cid, "orignan")
    base = val[(k0, "prcp")]
    base = base[base["n_both_finite"] > 0]
    rng = np.random.default_rng(Y)
    area = bf.set_index("canon_id")["geom_area_km2"]
    country = bf.set_index("canon_id")["gage_type"]
    for lo_, hi_ in [(0, 50), (50, 500), (500, 5000), (5000, 1e9)]:
        for ctry, n in [("USGS", 6), ("Canada", 3)]:
            pool = [cid for cid in base.index if lo_ <= area[cid] < hi_ and country[cid] == ctry]
            for cid in rng.choice(pool, size=min(n, len(pool)), replace=False):
                add(cid, "sample")
    for why, (v, asc) in {"wettest": ("prcp", False), "snowiest": ("swe", False),
                          "coldest": ("tmin", True), "hottest": ("tmax", False)}.items():
        d = val[(k0, v)]
        d = d[d["n_both_finite"] > 0].sort_values("old_mean", ascending=asc)
        add(d.index[0], why)
    cur = list(reasons)
    ncur = len(cur)
    cur_sites = set(bf.set_index("canon_id").loc[cur, "site_id"])

    stale_ids = duckdb.sql(f"SELECT DISTINCT site_id FROM read_parquet('{a.old}')").df()["site_id"]
    ids_sql = ",".join("'" + s.replace("'", "''") + "'" for s in stale_ids if canon(s) in reasons)
    old = duckdb.sql(f"SELECT site_id, Date, {', '.join(DAYMET_VARS)} FROM read_parquet('{a.old}') "
                     f"WHERE site_id IN ({ids_sql}) AND Date >= DATE '{Y}-01-01' AND Date <= DATE '{Y}-12-31'"
                     ).df() if ids_sql else pd.DataFrame(columns=["site_id", "Date"] + DAYMET_VARS)
    old["canon_id"] = old["site_id"].map(canon)
    old["Date"] = pd.to_datetime(old["Date"])
    fresh = {}
    for c in cfgs:
        for v in DAYMET_VARS:
            f = pd.read_parquet(os.path.join(c["out"], f"{v}_{Y}.parquet"))
            f = f[f["site_id"].isin(cur_sites)].copy()
            f["canon_id"] = f["site_id"].map(canon)
            f["Date"] = pd.to_datetime(f["Date"])
            fresh[(c["key"], v)] = f.pivot(index="Date", columns="canon_id", values=v)
    dates = fresh[(k0, "prcp")].index
    if len(dates) != 365 or dates[0] != pd.Timestamp(f"{Y}-01-01"):
        sys.exit("fresh series are not the 365-day Daymet calendar")
    old_missing = []
    for cid in cur:
        o = old[old["canon_id"] == cid]
        old_missing.append(bool(len(o) == 0 or o[DAYMET_VARS].isna().all().all()))
    series = {}
    for v in DAYMET_VARS:
        ov = old.pivot(index="Date", columns="canon_id", values=v) if len(old) else pd.DataFrame(index=dates)
        o_arr = np.full((ncur, 365), np.nan)
        d_arr = {c["key"]: np.full((ncur, 365), np.nan) for c in cfgs}
        for j, cid in enumerate(cur):
            if not old_missing[j]:
                o_arr[j] = ov[cid].reindex(dates).to_numpy(dtype=float)
            for c in cfgs:
                fr = fresh[(c["key"], v)][cid].reindex(dates).to_numpy(dtype=float)
                # with an original: the difference (float64, stored float32); without: the fresh values
                d_arr[c["key"]][j] = fr if old_missing[j] else fr - o_arr[j]
        series[v] = {"old": b64(o_arr.ravel()), **{f"d_{k}": b64(arr.ravel()) for k, arr in d_arr.items()}}

    data = {
        "meta": {"year": Y, "vars": DAYMET_VARS, "names": NAMES, "short": SHORT, "units": UNITS,
                 "cfgs": [{"key": c["key"], "label": c["label"]} for c in cfgs], "worstN": WORST_N,
                 "generated": time.strftime("%Y-%m-%d"), "n_polygons": int(nb), "n_shared": int(len(shared)),
                 "n_orignan": int(len(old_nan)), "n_new": int((status == 2).sum())},
        "basins": {"id": bf["site_id"].tolist(),
                   "country": ["US" if t == "USGS" else "CA" for t in bf["gage_type"].fillna("")],
                   "area": np.round(bf["geom_area_km2"].to_numpy(dtype=float), 1).tolist(),
                   "src": bf["watershed_geom_source"].fillna("").tolist(),
                   "simp": [bool(s) for s in simp],
                   "x": np.round(np.asarray(bx) / 1000.0, 1).tolist(),
                   "y": np.round(np.asarray(by) / 1000.0, 1).tolist(), "status": status.tolist()},
        "shared_idx": shared_idx.tolist(),
        "summary": summary,
        "curated": [{"i": idx_of[cid], "why": reasons[cid], "oldMissing": old_missing[j]} for j, cid in enumerate(cur)],
        "series": series,
        "outlines": outlines(a.assets_dir, to_lcc),
    }
    payload = json.dumps(data, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "daymet_dashboard_template.html")).read()
    page = tpl.replace("__DATA__", payload).replace("__YEAR__", str(Y))
    if a.fragment:
        with open(a.fragment, "w") as fh:
            fh.write(page)
    with open(a.out, "w") as fh:
        fh.write("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
                 "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
                 "</head>\n<body>\n" + page + "\n</body>\n</html>\n")
    print(f"{ncur} curated basins ({sum(old_missing)} without an original series); {len(shared)} shared "
          f"basins; {len(cfgs)} configuration(s); page {os.path.getsize(a.out) / 1e6:.1f} MB in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
