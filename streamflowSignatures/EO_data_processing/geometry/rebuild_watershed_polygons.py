#!/usr/bin/env python3
"""Rebuild the HISSS watershed polygon layer from the official sources.

This is the script that regenerated the delivered boundary layer (HydroShare Resource 4,
`hisss_watershed_boundaries.*`, 7,964 basins) after S3 access was lost (2026-08-25). It
consolidates the pipeline pieces in this folder (geometry_build_us / _ca / _hb_fallback,
build_v3) with an inlined residual-list derivation and the June post-step exclusion.
Committed 2026-09-29 with two opt-in variants for the Daymet reprocess; the defaults
reproduce Resource 4 exactly.

Stages: A US = GAGES-II boundaries; B Canada = ECCC MDA_ADP drainage basins; C+D residual
gages -> HydroBASINS lev-12 upstream union (outlet from Downstream_HB_ID or point-in-basin);
E merge + 4 inclusive-universe GAGES-II extras, zero-padded ids from the signature products,
exclusion, 200 m simplification (a basin stays full-resolution when simplifying moves its
area > 2 %), area QA flags.

Options:
  --no-simplify     keep every polygon at full resolution (geom_simplified = False).
  --include-large   keep the basins > 100,000 km² that the delivered layer drops; only the
                    known-bad polygons stay out (05KH009: a HydroBASINS fallback that drew a
                    328,000 km² river as 200 km²). Universe 8,018 -> 8,017 basins.

Usage:
  python rebuild_watershed_polygons.py [--source-dir ~/Downloads/geometry_rebuild]
      [--metadata <hisss_gage_metadata.csv>] [--golden <golden signatures csv>]
      [--signatures <R1 signatures csv>] --out-dir <dir> [--stamp rebuild_25aug2026]
      [--no-simplify] [--include-large] [--formats parquet,gpkg,csv]

Inputs (not in the repo): --source-dir holds basinAt_NorAm_polys.gpkg (HydroBASINS lev 12,
North America), official/gagesII_bnd/boundaries-shapefiles-by-aggeco/*.shp (USGS GAGES-II)
and official/ca/MDA_ADP_<nn>.gpkg (ECCC). The golden and R1 signature CSVs only supply the
zero-padded id spelling.

Validation targets for the defaults (June 27 delivery): 7,964 features; gagesii 6,164 /
wsc_eccc 1,771 / hydrobasins 29; 0 dups/invalid/empty; 141 basins kept full-res by the 2 %
guard; low_confidence = 29 HB + area outliers. (--no-simplify --include-large: 8,017; 30 HB;
57 low_confidence, i.e. 30 HB + 28 area outliers - 1 counted twice.)

A <stem>.provenance.json sidecar (added 2026-10-01) records the arguments, the md5 and
size of every input file, the outputs' md5, the git state of this folder and the library
versions.
"""
import argparse
import glob
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from collections import defaultdict

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import make_valid, union_all

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
LAYER = "DrainageBasin_BassinDeDrainage"
EQ = "EPSG:6933"
TOL = 0.002          # simplification tolerance, degrees (~200 m)
MAXREL = 0.02        # keep full resolution when simplifying changes the area by more than this
LARGE_KM2 = 100_000
KNOWN_BAD = {"05KH009"}


def canon(x):
    s = str(x).strip()
    return (str(int(s)) if s.isdigit() else s.upper())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", default=os.path.expanduser("~/Downloads/geometry_rebuild"))
    ap.add_argument("--metadata", default=os.path.expanduser(
        "~/Downloads/Signatures/resource3_input_data/hisss_gage_metadata.csv"))
    ap.add_argument("--golden", default=os.path.join(REPO, "golden-outputs", "streamflow_signatures_full_10feb2026.csv"))
    ap.add_argument("--signatures", default=os.path.expanduser(
        "~/Downloads/Signatures/resource1_signatures_wy1993_2025/hisss_signatures_wy1993_2025.csv"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stamp", default="rebuild_25aug2026")
    ap.add_argument("--no-simplify", action="store_true")
    ap.add_argument("--include-large", action="store_true")
    ap.add_argument("--formats", default="parquet,gpkg,csv")
    a = ap.parse_args()
    basin_gpkg = os.path.join(a.source_dir, "basinAt_NorAm_polys.gpkg")
    gages2 = os.path.join(a.source_dir, "official", "gagesII_bnd", "boundaries-shapefiles-by-aggeco")
    cadir = os.path.join(a.source_dir, "official", "ca")
    os.makedirs(a.out_dir, exist_ok=True)

    m = pd.read_csv(a.metadata, dtype=str, low_memory=False)
    m["basin_area"] = pd.to_numeric(m["basin_area"], errors="coerce")
    succ = m[m["processing_status"] == "success"].copy()
    ids = succ["gage_id"].astype(str)
    is_can = ids.str.match(r"^\d{2}[A-Za-z]")
    us_meta = succ[~is_can].copy()
    ca_meta = succ[is_can].copy()
    print(f"success gages: {len(succ)} (US-like {len(us_meta)}, CAN {len(ca_meta)})")

    # ---------- Stage A: US (GAGES-II) ----------
    us_canon = {canon(g): g for g in us_meta["gage_id"]}
    parts = []
    for f in sorted(glob.glob(os.path.join(gages2, "*.shp"))):
        g = gpd.read_file(f)
        parts.append(g[["GAGE_ID", "geometry"]])
    g2 = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=parts[0].crs)
    g2["geom_area_km2"] = g2.geometry.area / 1e6  # native Albers, equal-area
    g2["gid_canon"] = g2["GAGE_ID"].map(canon)
    g2u = g2[g2["gid_canon"].isin(us_canon)].copy()
    g2u = g2u.sort_values("geom_area_km2", ascending=False).drop_duplicates("gid_canon", keep="first")
    g2u["gage_id"] = g2u["gid_canon"].map(us_canon)
    us = g2u.to_crs(4326)[["gage_id", "geometry"]].copy()
    us["watershed_geom_source"] = "gagesii"
    us_matched = set(g2u["gage_id"])
    print(f"Stage A: US matched {len(us_matched)} / {len(us_meta)} (target 6,160)")

    # ---------- Stage B: Canada (ECCC MDA_ADP) ----------
    ca_set = set(ca_meta["gage_id"].astype(str))
    parts = []
    for nn in [f"{i:02d}" for i in range(1, 12)]:
        f = os.path.join(cadir, f"MDA_ADP_{nn}.gpkg")
        want = sorted(g for g in ca_set if g[:2] == nn)
        if not want:
            continue
        g = gpd.read_file(f, layer=LAYER)
        g = g[g["StationNum"].astype(str).isin(want)]
        parts.append(g)
        print(f"  MDA_{nn}: wanted {len(want)}, got {len(g)}")
    ca = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=parts[0].crs)
    ca["geom_area_km2"] = ca.geometry.area / 1e6  # Canada Albers, equal-area
    ca["gage_id"] = ca["StationNum"].astype(str).str.strip()
    ca = ca.sort_values("geom_area_km2", ascending=False).drop_duplicates("gage_id", keep="first")
    ca4 = ca.to_crs(4326)[["gage_id", "geometry"]].copy()
    ca4["watershed_geom_source"] = "wsc_eccc"
    ca_matched = set(ca4["gage_id"])
    print(f"Stage B: CAN matched {len(ca_matched)} / {len(ca_meta)} (June build: 1,823)")

    # ---------- Stage C+D: residual -> HydroBASINS fallback ----------
    resid = succ[~succ["gage_id"].isin(us_matched | ca_matched)].copy()
    resid["latitude"] = pd.to_numeric(resid["latitude"], errors="coerce")
    resid["longitude"] = pd.to_numeric(resid["longitude"], errors="coerce")
    print(f"Stage C: residual gages for HB fallback: {len(resid)} (June build: 31)")
    print("loading basinAt (167k lev12 basins)...")
    b = gpd.read_file(basin_gpkg)
    if b.crs is None or b.crs.to_epsg() != 4326:
        b = b.to_crs(4326)
    b["HYBAS_ID"] = b["HYBAS_ID"].astype("int64")
    b["NEXT_DOWN"] = b["NEXT_DOWN"].astype("int64")
    hyb = b["HYBAS_ID"].to_numpy()
    nxt = b["NEXT_DOWN"].to_numpy()
    row_of = {h: i for i, h in enumerate(hyb)}
    child_of = defaultdict(list)
    for h, n in zip(hyb, nxt):
        if n > 0:
            child_of[n].append(h)

    def bfs(outlet):
        if outlet not in row_of:
            return None
        seen = set()
        st = [outlet]
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x)
            st.extend(child_of.get(x, ()))
        return seen

    outlet_map = {}
    for _, r in resid.iterrows():
        hb_id = str(r["Downstream_HB_ID"])
        if hb_id.isdigit():
            outlet_map[r["gage_id"]] = int(hb_id)
    need_pt = resid[~resid["gage_id"].isin(outlet_map)].copy()
    if len(need_pt):
        pts = gpd.GeoDataFrame(need_pt, geometry=gpd.points_from_xy(need_pt.longitude, need_pt.latitude), crs=4326)
        j = gpd.sjoin(pts, b[["HYBAS_ID", "geometry"]], predicate="within", how="left")
        for _, r in j.iterrows():
            if pd.notna(r.get("HYBAS_ID")):
                outlet_map[r["gage_id"]] = int(r["HYBAS_ID"])
    print(f"  outlets resolved: {len(outlet_map)}/{len(resid)} (point-in-basin for {len(need_pt)})")
    rows, failed = [], []
    for _, r in resid.iterrows():
        g = r["gage_id"]
        o = outlet_map.get(g)
        if o is None:
            failed.append((g, "no outlet"))
            continue
        mem = bfs(o)
        if not mem:
            failed.append((g, "outlet not in basinAt"))
            continue
        idx = [row_of[h] for h in mem if h in row_of]
        rows.append({"gage_id": g, "watershed_geom_source": "hydrobasins",
                     "geometry": union_all(b.geometry.values[idx])})
    hb = gpd.GeoDataFrame(rows, crs=4326)
    print(f"Stage D: HB delineated {len(hb)} | failed {len(failed)} {failed}")

    # ---------- Stage E: merge, attrs, exclusion, simplify, flags ----------
    auth = set()
    for p in [a.golden, a.signatures]:
        auth |= set(pd.read_csv(p, usecols=[0], dtype=str).iloc[:, 0].astype(str).str.strip())
    amap = {}
    for x in auth:
        amap.setdefault(canon(x), set()).add(x)
    assert all(len(v) == 1 for v in amap.values()), "canon collision in signature ids!"
    amap = {k: next(iter(v)) for k, v in amap.items()}

    def padded(c):
        return amap.get(c, c.zfill(8) if c.isdigit() else c)

    new_ids = {"1591000", "1591400", "1591610", "1591700"}
    parts = []
    for f in glob.glob(os.path.join(gages2, "*.shp")):
        gg = gpd.read_file(f)
        gg = gg[gg["GAGE_ID"].map(canon).isin(new_ids)]
        if len(gg):
            parts.append(gg[["GAGE_ID", "geometry"]].rename(columns={"GAGE_ID": "gage_id"}))
    new = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=parts[0].crs).to_crs(4326)
    new["watershed_geom_source"] = "gagesii"
    print(f"Stage E: inclusive-universe extras from GAGES-II: {len(new)} (target 4)")

    final = gpd.GeoDataFrame(pd.concat([us, new, ca4, hb], ignore_index=True), crs="EPSG:4326")
    final["canon_id"] = final["gage_id"].map(canon)
    final["gage_id"] = final["canon_id"].map(padded)
    assert final["canon_id"].is_unique and final["gage_id"].is_unique
    m2 = m.copy()
    m2["canon"] = m2["gage_id"].map(canon)
    for c in ("latitude", "longitude"):
        m2[c] = pd.to_numeric(m2[c], errors="coerce")
    final = final.merge(m2[["canon", "basin_area", "latitude", "longitude", "gage_type"]],
                        left_on="canon_id", right_on="canon", how="left").drop(columns=["canon"])

    # full-res equal-area area, then the exclusion (the June delivery's post-step)
    af_all = final.geometry.to_crs(EQ).area / 1e6
    bad = final["gage_id"].isin(KNOWN_BAD)
    drop = bad if a.include_large else (af_all > LARGE_KM2) | bad
    print(f"exclusion: dropping {int(drop.sum())} basins ({'known-bad only' if a.include_large else 'target 54'}): "
          f"{final.loc[drop, 'watershed_geom_source'].value_counts().to_dict()}; "
          f"basins > {LARGE_KM2:,} km² kept: {int(((af_all > LARGE_KM2) & ~drop).sum())}")
    final = final[~drop].reset_index(drop=True)

    full = final.geometry
    if a.no_simplify:
        use_full = pd.Series(True, index=final.index)
    else:
        simp = full.simplify(TOL, preserve_topology=True)
        simp = gpd.GeoSeries([make_valid(s) if (s is not None and not s.is_empty and not s.is_valid) else s
                              for s in simp], crs=4326)
        af = full.to_crs(EQ).area / 1e6
        asi = simp.to_crs(EQ).area / 1e6
        rel = ((asi - af).abs() / af).replace([np.inf, -np.inf], np.nan)
        use_full = (rel > MAXREL) | asi.isna() | simp.is_empty
        final["geometry"] = gpd.GeoSeries([f if uf else s for uf, f, s in zip(use_full, full, simp)], crs=4326)
    final["geom_simplified"] = ~use_full.values
    final["geom_area_km2"] = final.geometry.to_crs(EQ).area / 1e6
    final["area_rel_diff"] = ((final["geom_area_km2"] - final["basin_area"]).abs()
                              / final["basin_area"]).replace([np.inf, -np.inf], np.nan)
    final["area_flag"] = final["area_rel_diff"] > 0.5
    final["low_confidence"] = (final["watershed_geom_source"] == "hydrobasins") | final["area_flag"].fillna(False)

    print("=== QA (rebuild) vs recorded June targets (defaults only) ===")
    print("rows:", len(final), "(target 7,964) | by source:", final["watershed_geom_source"].value_counts().to_dict(),
          "(target gagesii 6,164 / wsc_eccc 1,771 / hydrobasins 29)")
    print("dups:", int(final["gage_id"].duplicated().sum()), "| invalid:", int((~final.geometry.is_valid).sum()),
          "| empty:", int(final.geometry.is_empty.sum()), "(targets 0/0/0)")
    print("kept full-res:", int((~final["geom_simplified"]).sum()), "(June: 141) | simplified:",
          int(final["geom_simplified"].sum()))
    g = final["area_rel_diff"].dropna()
    print(f"area QA n={len(g)}: median {g.median():.4f} | <10% {(g < 0.1).mean() * 100:.1f}% | >50% {int((g > 0.5).sum())}")
    print("low_confidence:", int(final["low_confidence"].sum()))

    cols = ["gage_id", "canon_id", "watershed_geom_source", "geom_area_km2", "basin_area", "area_rel_diff",
            "area_flag", "geom_simplified", "low_confidence", "latitude", "longitude", "gage_type", "geometry"]
    stem = os.path.join(a.out_dir, f"watershed_polygons_{a.stamp}")
    formats = set(a.formats.split(","))
    if "gpkg" in formats:
        if os.path.exists(stem + ".gpkg"):
            os.remove(stem + ".gpkg")
        final[cols].to_file(stem + ".gpkg", driver="GPKG", layer="watersheds")
    if "parquet" in formats:
        final[cols].to_parquet(stem + ".parquet")
    if "csv" in formats:
        final[cols].drop(columns=["geometry"]).to_csv(stem + "_qa.csv", index=False)
    print(f"wrote {stem}.{{{','.join(sorted(formats))}}}")
    write_provenance(a, stem, final, [basin_gpkg, a.metadata, a.golden, a.signatures]
                     + sorted(glob.glob(os.path.join(gages2, "*")))
                     + sorted(glob.glob(os.path.join(cadir, "MDA_ADP_*.gpkg"))))


def write_provenance(a, stem, final, inputs):
    def md5(path):
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 24), b""):
                h.update(b)
        return h.hexdigest()

    def git(*args):
        r = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__))] + list(args),
                           capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    import pyogrio
    import pyproj
    import shapely
    outs = [p for p in (stem + ".parquet", stem + ".gpkg", stem + "_qa.csv") if os.path.exists(p)]
    prov = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "argv": sys.argv, "args": vars(a),
            "inputs": {p: {"bytes": os.path.getsize(p), "md5": md5(p)} for p in inputs if os.path.isfile(p)},
            "outputs": {os.path.basename(p): {"bytes": os.path.getsize(p), "md5": md5(p)} for p in outs},
            "rows": int(len(final)), "by_source": final["watershed_geom_source"].value_counts().to_dict(),
            "low_confidence": int(final["low_confidence"].sum()), "area_flag": int(final["area_flag"].fillna(False).sum()),
            "git": {"commit": git("rev-parse", "HEAD"),
                    "geometry_tools_dirty": bool(git("status", "--porcelain", "--", "."))},
            "software": {"python": sys.version.split()[0], "platform": platform.platform(),
                         "geopandas": gpd.__version__, "pandas": pd.__version__, "shapely": shapely.__version__,
                         "geos": shapely.geos_version_string, "pyogrio": pyogrio.__version__,
                         "gdal": getattr(pyogrio, "__gdal_version_string__", "unknown"),
                         "pyproj": pyproj.__version__, "proj": pyproj.proj_version_str}}
    with open(stem + ".provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2, default=str)
    print(f"wrote {stem}.provenance.json ({len(prov['inputs'])} inputs hashed)")


if __name__ == "__main__":
    main()
