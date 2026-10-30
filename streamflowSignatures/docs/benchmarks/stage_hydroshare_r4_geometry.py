#!/usr/bin/env python
"""Stage the HydroShare Resource 4 (watershed geometry and basin attributes) update of 2026-10-30.

DECISION (user, 2026-10-30): R4 keeps its seven deposited files unchanged except the README, and ADDS the
unsimplified boundary layer as a GeoPackage — the full-resolution outlines of the same 7,964 basins plus the
53 basins above 100,000 km² (8,017 polygons), i.e. the layer over which the Daymet forcing of Resource 3 was
aggregated (`EO_data_processing/README_DAYMET.md`, "Polygons"). The parquet of that layer was built on
2026-09-29 by `EO_data_processing/geometry/rebuild_watershed_polygons.py --no-simplify --include-large` and
rebuilt byte-identically on 2026-10-01 (`*.provenance_rebuild_2026-10-01.json` beside it).

What this tool does
  1. reads the full-resolution GeoParquet, writes it as a GeoPackage (layer `watersheds`, the same call the
     rebuild script uses for the deposited layer) and verifies the GeoPackage against its source: row count,
     ids, every attribute column, and every geometry byte-identical as WKB (polygons promoted to multipolygons,
     as the GeoPackage driver stores them);
  2. verifies the two layers against each other and against the attribute table, the frozen product gage sets
     and the Resource 3 Daymet flags — every number printed in the README is computed here;
  3. renders the README and the abstract, copies the six unchanged files from the August staging, and writes
     the staging and rename manifests; with --upload it also writes the upload folder (README + GeoPackage)
     and the `_METADATA_and_MD5.txt` note.

Usage
  python docs/benchmarks/stage_hydroshare_r4_geometry.py --out ~/Downloads/Signatures_v3_28oct2026/resource4_geometry_attributes \
      [--upload ~/Downloads/hisss_R4_upload_final_2026-10-30] [--skip-gpkg]   # --skip-gpkg re-renders the README only

Inputs (defaults below): the August R4 staging (byte-identical to the live resource, verified by MD5 on
2026-10-29), the full-resolution parquet, the patched gage metadata and the frozen selection lists of the
version-3 staging, and the Daymet flags table.
"""
import argparse
import csv
import datetime as dt
import hashlib
import os
import shutil
import statistics
import subprocess
import sys

HOME = os.path.expanduser("~")
DEF_AUG = f"{HOME}/Downloads/Signatures/resource4_geometry_attributes"
DEF_FULL = f"{HOME}/HISSS_data/daymet-processed-29sep2026/polygons/watershed_polygons_fullres_large_29sep2026.parquet"
DEF_R3 = f"{HOME}/Downloads/Signatures_v3_28oct2026/resource3_input_data"
FULL_NAME = "hisss_watershed_boundaries_unsimplified.gpkg"
UNCHANGED = ["hisss_hydroatlas_attributes_dictionary.csv", "hisss_watershed_boundaries_qa.csv", "hisss_watershed_boundaries.gpkg",
             "hisss_watershed_boundaries.parquet", "hisss_watershed_hydroatlas_attributes.csv", "hisss_watershed_hydroatlas_attributes.parquet"]
LARGE_KM2 = 100000
# From the Daymet build's 2023 gate (docs/plans/2026-09-29-daymet-reprocessing-action-plan.md §2, EO_data_processing/README_DAYMET.md):
# aggregating over the simplified Resource 4 polygons instead of the full-resolution ones moves prcp annual totals
# within ±0.1 % (1st–99th percentile across basins; up to 1.2 % in small high-relief basins) and swe annual sums
# within ±2.5 % (p01–p99; far more in trace-snow basins).
DAYMET_SHIFT_PRCP, DAYMET_SHIFT_PRCP_MAX, DAYMET_SHIFT_SWE = "0.1%", "1.2%", "2.5%"
FULL_BUILT, FULL_REBUILT = "29 September 2026", "1 October 2026"
COLLECTION = "f702201faa5d46069a5ee83ffa4c9768"


def sha(p, blk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(blk), b""):
            h.update(b)
    return h.hexdigest()


def md5(p, blk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(blk), b""):
            h.update(b)
    return h.hexdigest()


def copy(src, dst):
    if os.path.exists(dst):
        os.remove(dst)
    r = subprocess.run(["cp", "-p", "-c", src, dst], capture_output=True)
    if r.returncode != 0:
        shutil.copy2(src, dst)
    return os.path.basename(src), os.path.basename(dst)


def key(s):
    s = s.strip()
    return s.lstrip("0") if s.isdigit() else s.upper()


def readlist(p):
    with open(p) as f:
        return [l.strip() for l in f if l.strip()]


def to_multi(geoms):
    from shapely import MultiPolygon
    return [g if g.geom_type == "MultiPolygon" else MultiPolygon([g]) for g in geoms]


def build_gpkg(full_parquet, out_gpkg):
    """Write the GeoPackage from the GeoParquet and verify it byte-for-byte (attributes and WKB)."""
    import geopandas as gpd
    import numpy as np
    from shapely import to_wkb
    gdf = gpd.read_parquet(full_parquet)
    assert gdf.crs is not None and gdf.crs.to_epsg() == 4326, gdf.crs
    assert gdf.geometry.name == "geometry"
    cols = [c for c in gdf.columns if c != "geometry"]
    if os.path.exists(out_gpkg):
        os.remove(out_gpkg)
    gdf.to_file(out_gpkg, driver="GPKG", layer="watersheds")
    src_attr = gdf[cols].reset_index(drop=True)
    src_wkb = np.array([bytes(b) for b in to_wkb(np.array(to_multi(gdf.geometry.values), dtype=object))], dtype=object)
    n_src = len(gdf)
    del gdf
    back = gpd.read_file(out_gpkg, layer="watersheds")
    assert back.crs.to_epsg() == 4326, back.crs
    assert len(back) == n_src, (len(back), n_src)
    assert list(back.columns) == cols + ["geometry"], list(back.columns)
    for c in cols:
        a, b = src_attr[c].tolist(), back[c].tolist()
        if src_attr[c].dtype.kind == "f":
            ok = all((x == y) or (x != x and y != y) for x, y in zip(a, b))
        else:
            ok = a == b
        assert ok, f"attribute column {c} differs between the parquet and the GeoPackage"
    back_wkb = np.array([bytes(b) for b in to_wkb(np.array(to_multi(back.geometry.values), dtype=object))], dtype=object)
    n_diff = int(sum(1 for x, y in zip(src_wkb, back_wkb) if x != y))
    assert n_diff == 0, f"{n_diff} geometries differ as WKB between the parquet and the GeoPackage"
    n_poly_promoted = int(sum(1 for g in back.geometry.values if g.geom_type != "MultiPolygon"))
    valid = int(back.geometry.is_valid.sum())
    return {"n": n_src, "wkb_identical": n_src - n_diff, "n_not_multi": n_poly_promoted, "valid": valid}


def facts(a):
    """Every number the README states, computed from the files."""
    import geopandas as gpd
    import pyarrow.parquet as pq
    F = {}
    simp = gpd.read_parquet(os.path.join(a.aug_dir, "hisss_watershed_boundaries.parquet"))
    full = gpd.read_parquet(a.full_parquet)
    assert list(simp.columns) == list(full.columns), (list(simp.columns), list(full.columns))
    F["cols"] = [c for c in simp.columns if c != "geometry"]
    F["n_simp"], F["n_full"] = len(simp), len(full)
    assert simp["gage_id"].is_unique and full["gage_id"].is_unique
    s_ids, f_ids = set(simp["gage_id"]), set(full["gage_id"])
    assert s_ids <= f_ids, "the simplified layer holds ids absent from the full-resolution layer"
    extra = full[~full["gage_id"].isin(s_ids)]
    F["n_extra"] = len(extra)
    assert (extra["gage_type"] == "Canada").all() and (extra["geom_area_km2"] > LARGE_KM2).all(), "the extra basins are not all Canadian and > 100,000 km²"
    assert "05KH009" not in f_ids and "05KH009" not in s_ids
    assert (~full["geom_simplified"]).all(), "the full-resolution layer has simplified rows"
    F["extra_area_min"], F["extra_area_max"] = float(extra["geom_area_km2"].min()), float(extra["geom_area_km2"].max())
    for name, g in (("simp", simp), ("full", full)):
        vc = g["watershed_geom_source"].value_counts()
        F[f"src_{name}"] = {k: int(vc.get(k, 0)) for k in ("gagesii", "wsc_eccc", "hydrobasins")}
        F[f"lc_{name}"] = int(g["low_confidence"].sum())
        F[f"af_{name}"] = int(g["area_flag"].sum())
        tv = g["gage_type"].value_counts()
        F[f"us_{name}"], F[f"ca_{name}"] = int(tv.get("USGS", 0)), int(tv.get("Canada", 0))
    F["n_guard"] = int((~simp["geom_simplified"]).sum())
    # co-located gages can share one outline: ids are unique, geometries need not be (Codex review 2026-10-30)
    import collections
    import shapely
    for name, g in (("simp", simp), ("full", full)):
        cnt = collections.Counter(bytes(b) for b in shapely.to_wkb(g.geometry.values))
        groups = [k for k, v in cnt.items() if v > 1]
        F[f"dup_{name}"] = (len(groups), sum(cnt[k] - 1 for k in groups))
        assert g["gage_id"].is_unique and g.geometry.is_valid.all() and (~g.geometry.is_empty).all()
    F["pct_wsc_simp"] = 100 * F["src_simp"]["wsc_eccc"] / F["ca_simp"]
    F["pct_wsc_full"] = 100 * F["src_full"]["wsc_eccc"] / F["ca_full"]
    lc_only_full = set(full.loc[full["low_confidence"], "gage_id"]) - set(simp.loc[simp["low_confidence"], "gage_id"])
    assert lc_only_full <= set(extra["gage_id"]), lc_only_full
    F["lc_extra_ids"] = sorted(lc_only_full)
    # the 141 basins the guard kept at full resolution have the same vertices in both layers — the two builds ran in
    # different environments, so coordinates agree to floating-point precision (~1e-13°), not byte for byte; the other
    # shared basins differ by < 2 % in area (the guard)
    import shapely
    m = simp.merge(full[["gage_id", "geom_area_km2", "geometry"]], on="gage_id", suffixes=("", "_full"))
    guard = m[~m["geom_simplified"]]
    nv_equal = int((shapely.get_num_coordinates(guard["geometry"].values) == shapely.get_num_coordinates(guard["geometry_full"].values)).sum())
    hd = shapely.hausdorff_distance(guard["geometry"].values, guard["geometry_full"].values)
    F["guard"] = (nv_equal, len(guard), float(hd.max()))
    assert nv_equal == len(guard) and hd.max() < 1e-9, F["guard"]
    rel = ((m["geom_area_km2"] - m["geom_area_km2_full"]).abs() / m["geom_area_km2_full"])
    assert rel.max() < 0.02, float(rel.max())
    F["area_shift_max_pct"] = float(100 * rel.max())
    F["area_shift_median_pct"] = float(100 * rel.median())
    ard = simp["area_rel_diff"].dropna().abs()
    F["area_median"], F["area_within10"], F["n_area_ref"] = float(ard.median()), float(100 * (ard <= 0.10).mean()), int(len(ard))
    # un-padded ids: boundary ids whose agency form in the attribute table differs
    hp = pq.ParquetFile(os.path.join(a.aug_dir, "hisss_watershed_hydroatlas_attributes.parquet"))
    F["attr_cols"] = hp.schema_arrow.names
    attr_ids = hp.read(columns=["gage_id"]).column("gage_id").to_pylist()
    F["n_attr"] = len(attr_ids)
    lens = sorted({len(g) for g in attr_ids if g.isdigit()})
    F["usgs_id_lengths"] = ", ".join(str(x) for x in lens[:-1]) + f" or {lens[-1]}"
    F["ex_simp"] = sorted(simp.loc[simp["watershed_geom_source"] == "gagesii", "gage_id"])[0]
    F["ex_large"] = sorted(extra["gage_id"])[0]
    agency = {key(g): g for g in attr_ids}
    unp = sorted((g, agency[c]) for g, c in zip(simp["gage_id"], simp["canon_id"]) if c in agency and agency[c] != g)
    unp_full = sorted((g, agency[c]) for g, c in zip(full["gage_id"], full["canon_id"]) if c in agency and agency[c] != g)
    assert unp == unp_full, "the two boundary layers do not share the same un-padded ids"
    F["unpadded"] = unp
    assert all(key(g) == c for g, c in zip(full["gage_id"], full["canon_id"]))
    sb, sa = set(simp["canon_id"]), set(agency)
    F["overlap"], F["attr_only"], F["bound_only"] = len(sb & sa), len(sa - sb), sorted(sb - sa)
    F["full_overlap"] = len(set(full["canon_id"]) & sa)
    # the CSV is a decimal export of the parquet: same rows, columns and ids; floats agree to ~1e-9 (Codex review)
    import numpy as np
    import pandas as pd
    csvt = pd.read_csv(os.path.join(a.aug_dir, "hisss_watershed_hydroatlas_attributes.csv"), dtype={"gage_id": str}, float_precision="round_trip")
    pqt = hp.read().to_pandas()
    assert list(csvt.columns) == list(pqt.columns) and len(csvt) == len(pqt) == F["n_attr"]
    assert (csvt["gage_id"].tolist() == pqt["gage_id"].astype(str).tolist())
    mx = 0.0
    for c in pqt.columns:
        if pd.api.types.is_numeric_dtype(pqt[c]):
            x = pd.to_numeric(csvt[c], errors="coerce").to_numpy(dtype=float); y = pqt[c].to_numpy(dtype=float)
            assert (np.isnan(x) == np.isnan(y)).all(), c
            if (~np.isnan(x)).any():
                mx = max(mx, float(np.nanmax(np.abs(x - y))))
        else:
            assert (csvt[c].astype(str) == pqt[c].astype(str)).all(), c
    assert mx < 1e-8, mx
    F["csv_vs_parquet_max"] = mx
    with open(a.metadata, newline="") as f:
        meta = {key(r["gage_id"]): r for r in csv.DictReader(f)}
    F["bound_only_status"] = sorted({meta[g]["processing_status"] for g in F["bound_only"] if g in meta})
    F["bound_only_agency"] = [agency.get(g) or next((m for m in meta if key(m) == g), g) for g in F["bound_only"]]
    F["bound_only_agency"] = [meta[g]["gage_id"].zfill(8) if g in meta and meta[g]["gage_id"].isdigit() else g for g in F["bound_only"]]
    sel = {w: set(map(key, readlist(os.path.join(DEF_R3, f"hisss_gage_selection_wy{w}_2025.txt")))) for w in ("1993", "1980")}
    assert not any(g in s for g in F["bound_only"] for s in sel.values()), "a boundary-only gage is in a signature product"
    F["n_sel"] = {w: len(s) for w, s in sel.items()}
    F["sel_with_simp"] = {w: len(s & sb) for w, s in sel.items()}
    F["sel_with_full"] = {w: len(s & set(full["canon_id"])) for w, s in sel.items()}
    F["sel_exact_simp"] = {w: len(set(map(key, readlist(os.path.join(DEF_R3, f"hisss_gage_selection_wy{w}_2025.txt")))) & sb) - len(set(readlist(os.path.join(DEF_R3, f"hisss_gage_selection_wy{w}_2025.txt"))) & set(simp["gage_id"])) for w in sel}
    with open(a.daymet_flags, newline="") as f:
        fl = list(csv.DictReader(f))
    F["n_flags"] = len(fl)
    fl_ids = {key(r["canon_id"]) for r in fl}
    assert fl_ids <= set(full["canon_id"]), "Daymet flags name basins absent from the full-resolution layer"
    F["flags_lc"] = sum(1 for r in fl if r["low_confidence"].upper() == "TRUE")
    assert F["flags_lc"] == F["lc_full"], (F["flags_lc"], F["lc_full"])
    return F


def render_readme(F, gpkg_bytes):
    gb = gpkg_bytes / 1e9
    unp_stored = ", ".join(f"`{g}`" for g, _ in F["unpadded"])
    unp_agency = ", ".join(f"`{x}`" for _, x in F["unpadded"])
    bo = ", ".join(F["bound_only_agency"])
    ss, sf = F["src_simp"], F["src_full"]
    return f"""# HISSS Watershed Geometry and Basin Attributes

Part of the **Hydrologic Information Signatures and Summary Statistics (HISSS)**
dataset. This resource contains (1) the **watershed boundary polygons** for the gaged
watersheds underlying HISSS — a simplified layer and the same outlines at full
resolution, plus the {F['n_extra']} large basins the simplified layer omits — and (2) a
**watershed-scale physiographic attribute table** —
approximately 200 hydro-environmental descriptors per gage aggregated over each gage's
full upstream contributing area from HydroATLAS (BasinATLAS v10; Linke et al. 2019),
spanning hydrology, climate, physiography, land cover, soils and geology, and
anthropogenic influence. All layers are keyed by `gage_id` for joins to the other HISSS
resources (join on the zero-stripped form — see *Reading the data*); row coverage varies
by product (see below).

## Files

| File | Description |
|---|---|
| `hisss_watershed_boundaries.gpkg` | Watershed boundary polygons for {F['n_simp']:,} gaged watersheds, simplified at ~200 m (GeoPackage, WGS84 / EPSG:4326), one polygon per gage |
| `hisss_watershed_boundaries.parquet` | The same layer as GeoParquet |
| `hisss_watershed_boundaries_qa.csv` | Per-watershed QA companion of the simplified layer: geometry source, simplification and confidence flags, area diagnostics |
| `{FULL_NAME}` | The **unsimplified** layer: the same {F['n_simp']:,} outlines at full resolution plus the {F['n_extra']} basins larger than {LARGE_KM2:,} km² that the simplified layer omits — {F['n_full']:,} polygons (GeoPackage, WGS84, {gb:.1f} GB, same columns). The layer over which the Daymet forcing of the harmonized-input resource was aggregated |
| `hisss_watershed_hydroatlas_attributes.csv` / `.parquet` | HydroATLAS-derived attribute table: {F['n_attr']:,} gages × {len(F['attr_cols'])} columns (the same table in both formats; the CSV's decimal export agrees with the parquet's floating-point values to within {F['csv_vs_parquet_max']:.0e}) |
| `hisss_hydroatlas_attributes_dictionary.csv` | Abbreviated definition, HydroATLAS theme/unit codes, and aggregation method for every attribute column (full variable documentation and class-code legends: BasinATLAS v10 catalog at hydrosheds.org) |

## Watershed boundaries

Boundaries come from official agency products wherever they exist, with a documented
fallback: **GAGES-II** basin boundaries for all US gages (Falcone 2011), the Water
Survey of Canada **hydrometric basin polygons** for {F['pct_wsc_simp']:.0f}% of Canadian gages
(Environment and Climate Change Canada 2016), and **HydroBASINS level-12 upstream
aggregation** (Lehner and Grill 2013) for the Canadian gages with no official polygon
({ss['hydrobasins']} in the simplified layer, {sf['hydrobasins']} in the unsimplified one). Each feature records its origin in `watershed_geom_source`
(`gagesii`, `wsc_eccc`, or `hydrobasins`).

**Simplified layer** (`hisss_watershed_boundaries.gpkg` / `.parquet`, {F['n_simp']:,} polygons;
sources gagesii {ss['gagesii']:,} / wsc_eccc {ss['wsc_eccc']:,} / hydrobasins {ss['hydrobasins']}). Geometries are
topology-preserving simplified (~200 m tolerance; {F['n_guard']} small or narrow basins retained
at full resolution where simplification would shift area by more than 2%;
`geom_simplified` flags which), and `geom_area_km2` is recomputed from the delivered
geometry in an equal-area projection. `low_confidence` flags HydroBASINS-fallback basins
and area outliers; {F['lc_simp']} of the {F['n_simp']:,} boundaries carry the flag. The layer excludes
{F['n_extra']} basins larger than {LARGE_KM2:,} km² plus `05KH009` (a documented HydroBASINS fallback
mis-delineation) — {F['n_extra'] + 1} exclusions total — so it covers {F['n_simp']:,} of the
{F['n_attr'] + len(F['bound_only']):,} gages in the HISSS universe (the {F['n_attr']:,} compiled gages plus the four US gages
noted under *Basin attributes*).

**Unsimplified layer** (`{FULL_NAME}`, {F['n_full']:,} polygons; sources
gagesii {sf['gagesii']:,} / wsc_eccc {sf['wsc_eccc']:,} / hydrobasins {sf['hydrobasins']}). The same outlines before
simplification (`geom_simplified` is FALSE for every row; `geom_area_km2` is the
full-resolution area) plus the {F['n_extra']} Canadian basins above {LARGE_KM2:,} km²
({F['extra_area_min']:,.0f}–{F['extra_area_max']:,.0f} km²); `05KH009` is excluded from both layers.
{F['lc_full']} polygons carry `low_confidence` (the {F['lc_simp']} of the simplified layer plus one large
HydroBASINS-fallback basin, `{F['lc_extra_ids'][0]}`). Every gage of both signature products has a
polygon in this layer ({F['n_sel']['1993']:,} and {F['n_sel']['1980']:,} gages; the simplified layer lacks
{F['n_sel']['1993'] - F['sel_with_simp']['1993']} and {F['n_sel']['1980'] - F['sel_with_simp']['1980']} of them, all among the large basins). The {F['guard'][1]} basins that the
simplified layer kept at full resolution have the same vertices in both layers (coordinates
agree to within {F['guard'][2]:.1e} degrees — floating-point precision; the two layers were built in
different software environments); every other shared basin differs by less than 2% in area
(at most {F['area_shift_max_pct']:.3f}%, median {F['area_shift_median_pct']:.3f}%). **The Daymet forcing in the
harmonized-input resource** (`hisss_daymet_basin_daily.parquet`, {F['n_full']:,} basins) **was
aggregated over this layer** with fractional cell-coverage weights; aggregating over the
simplified layer instead moves annual precipitation totals by within about ±{DAYMET_SHIFT_PRCP}
(1st–99th percentile across basins; up to {DAYMET_SHIFT_PRCP_MAX} in small high-relief basins) and
annual snow water equivalent by within about ±{DAYMET_SHIFT_SWE} (1st–99th percentile; far more in
basins with only trace snow), as measured in the Daymet build's validation on calendar year 2023. The `hisss_daymet_basin_flags.csv` table of
that resource lists the {F['n_flags']} basins of this layer with low-confidence outlines or few Daymet
cells. The vegetation and land-cover products were computed on the simplified layer (its
original June 2026 build).

**Provenance — August 2026 rebuild.** The original June 2026 boundary layer was lost
in a storage-platform transition; the simplified layer deposited here was reconstructed
on 25 August 2026 from the project's committed workflow and the same source products — a
fresh USGS GAGES-II download; ECCC hydrometric basin polygons as served on
2024-09-05 (slightly newer than the June build's source); and HydroBASINS North
America level-12 — with the June 2026 build's previously uncommitted {F['n_extra'] + 1}-basin
exclusion step made explicit in the rebuild script. The rebuild reproduced every
recorded property of the June layer exactly:
{F['n_simp']:,} features; source counts gagesii {ss['gagesii']:,} / wsc_eccc {ss['wsc_eccc']:,} / hydrobasins {ss['hydrobasins']}; unique
ids and zero invalid or empty geometries ({F['dup_simp'][0]} groups of co-located gages share one
outline — {F['dup_simp'][1]} further gages; {F['dup_full'][0]} such groups in the unsimplified layer); the identical
{F['n_guard']} basins retained at full resolution by the 2% simplification guard; and identical
area-agreement statistics (median relative difference {F['area_median']:.3f}; {F['area_within10']:.1f}% of the
{F['n_area_ref']:,} polygons with a published drainage area within 10%). The vegetation and land-cover
companion products were computed on the original June layer; the rebuild is
equivalent by all recorded measures but not guaranteed byte-identical to it. The
unsimplified layer was built on {FULL_BUILT} by the same script with simplification
off and the large basins kept, and rebuilt byte-identically from the recorded inputs on
{FULL_REBUILT}; the GeoPackage deposited here was written from that build, and every
geometry is identical to it when compared as well-known binary (WKB), with single polygons
stored as multipolygons.

## Basin attributes

One row per compiled gage ({F['n_attr']:,}). The attribute and simplified-boundary layers overlap on
{F['overlap']:,} gages: the {F['attr_only']} basins excluded from the simplified layer retain attribute
rows (all {F['n_attr']:,} compiled gages have one), while four US gages whose streamflow compile did
not complete ({bo}; in neither signature product) have a boundary
but no attribute row. The attribute and unsimplified-boundary layers overlap on {F['full_overlap']:,} gages.
Attributes are drawn from BasinATLAS v10 over each gage's delineated upstream set of
HydroBASINS level-12 subbasins, using a hybrid aggregation:

- **Upstream-accumulated (`_u`) and pour-point (`_p`) attributes** (91) pass through
  from the gage's outlet subbasin — HydroATLAS has already accumulated these over the
  same upstream extent, so the outlet value is the rigorous watershed value.
- **Subbasin-scale (`_s`) attributes** are re-aggregated over the member subbasins:
  area-weighted means for continuous, percentage, and monthly-climate variables
  (HydroATLAS −9999 NoData masked), spatial minimum/maximum for elevation, and
  dominant classes for categorical variables — via argmax of the outlet's upstream
  (`_u`) class fractions for land cover (`glc`) and potential natural vegetation
  (`pnv`), and the SUB_AREA-weighted mode for the other nine categorical variables.

Diagnostics include `watershed_area_rel_diff` — the relative difference between the
summed member-subbasin area and the gage's published drainage area — for flagging
gages whose delineated watershed diverges from the reported drainage (for example,
mis-snapped outlets). The dictionary states the aggregation method used for every
column; its definitions are abbreviated, and categorical class codes (for example
land cover, potential natural vegetation, and wetland classes) are defined in the
BasinATLAS v10 catalog rather than restated here.

## Reading the data

```python
import geopandas as gpd, pandas as pd
basins = gpd.read_file("hisss_watershed_boundaries.gpkg")
attrs  = pd.read_csv("hisss_watershed_hydroatlas_attributes.csv", dtype={{"gage_id": str}})
# the unsimplified layer is {gb:.1f} GB — read a subset with an SQL filter rather than the whole file
# ({F['ex_simp']} is in both layers, {F['ex_large']} only in the unsimplified one):
full   = gpd.read_file("{FULL_NAME}", layer="watersheds",
                       where="gage_id IN ('{F['ex_simp']}', '{F['ex_large']}')")
```

```r
basins <- sf::st_read("hisss_watershed_boundaries.gpkg")
attrs  <- readr::read_csv("hisss_watershed_hydroatlas_attributes.csv",
                          col_types = readr::cols(gage_id = readr::col_character()))
full   <- sf::st_read("{FULL_NAME}",
                      query = "SELECT * FROM watersheds WHERE gage_id IN ('{F['ex_simp']}', '{F['ex_large']}')")
```

**Join rule (all HISSS tables).** Read `gage_id` as a character string, never
numeric. Join on the **zero-stripped** form — strip leading zeros from all-digit ids,
leave HYDAT ids unchanged — and do so on *both* sides of every join; HISSS station ids are
unique under this rule (17,166 distinct stations across all HISSS tables, zero collisions). **Do not re-pad**:
padding to 8 digits cannot recover 9- and 10-digit USGS site numbers that begin with a
zero (`021650905`, `0212433550`), and some HISSS tables store exactly those ids
un-padded. In Python: `key = lambda s: s.lstrip("0") if s.isdigit() else s.upper()`;
in R: `key <- function(x) ifelse(grepl("^[0-9]+$", x), sub("^0+", "", x), toupper(x))`.

In both boundary layers and the QA CSV, **{len(F['unpadded'])}** long USGS site numbers are stored in
`gage_id` without their leading zero ({unp_stored} for agency ids {unp_agency}); the `canon_id` column is the stripped key for every row and
joins correctly. The attribute table stores all {F['n_attr']:,} ids in agency form (USGS site numbers of
{F['usgs_id_lengths']} digits with their leading zeros — the dictionary's "8-digit" wording applies to
most but not all of them). The boundary layers are deposited with the ids exactly as built
rather than rewritten.

## Caveats

- **Subbasin granularity**: HydroBASINS level-12 subbasins average roughly 100 km²
  (ranging from a few to over 200 km²), so for gages whose watersheds are much
  smaller than their containing subbasin, attribute values describe the subbasin
  rather than the exact watershed. Use `watershed_area_rel_diff` to screen such
  cases.
- **Area diagnostics for US gages are a consistency check, not independent
  validation** — the published drainage areas are themselves GAGES-II-derived.
- The attribute table's source values are from HydroATLAS/BasinATLAS v10; consult the
  BasinATLAS v10 catalog (hydrosheds.org) for the underlying source datasets of each
  variable.

## License and citation

License: CC-BY 4.0. HydroATLAS/BasinATLAS is itself licensed CC-BY 4.0 — cite Linke
et al. (2019) when using the attribute table. Boundary sources: Falcone (2011,
GAGES-II); Environment and Climate Change Canada (2016); Lehner and Grill (2013,
HydroBASINS). Please also cite the HISSS data paper and this resource (DOI: to be
assigned at publication).

This resource is part of the HISSS collection on HydroShare:
https://www.hydroshare.org/resource/{COLLECTION}/
(collection DOI upon publication: https://doi.org/10.4211/hs.{COLLECTION}).
"""


def render_abstract(F):
    ss = F["src_simp"]
    return (
        "This resource contains the watershed geometry and basin-attribute layers of the Hydrologic Information Signatures and Summary "
        f"Statistics (HISSS) dataset. The boundary layer provides watershed polygons for {F['n_simp']:,} gaged watersheds across the United States "
        "and Canada (GeoPackage and GeoParquet, WGS84), assembled from official agency products: USGS GAGES-II basin boundaries for all US gages "
        f"and Water Survey of Canada hydrometric basin polygons for {F['pct_wsc_simp']:.0f}% of Canadian gages, with HydroBASINS level-12 upstream aggregation "
        f"for the {ss['hydrobasins']} Canadian gages lacking an official polygon. Each polygon records its source, simplification status, and area "
        "diagnostics, repeated in a companion QA table. A second GeoPackage holds the same outlines at full resolution, before the 200 m "
        f"simplification, together with the {F['n_extra']} basins larger than {LARGE_KM2:,} km2 that the simplified layer omits ({F['n_full']:,} polygons); "
        "it is the layer over which the Daymet forcing of the harmonized-input resource was aggregated. The attribute table covers all "
        f"{F['n_attr']:,} compiled HISSS gages with {len(F['attr_cols'])} columns ({F['attr_cols'].index('dis_m3_pyr')} identifiers and diagnostics plus "
        f"{len(F['attr_cols']) - F['attr_cols'].index('dis_m3_pyr')} HydroATLAS-derived attributes), aggregated over each gage's full upstream contributing area "
        "from BasinATLAS v10 and spanning hydrology, climate, physiography, land cover, soils and geology, and anthropogenic influence. "
        "Upstream-accumulated HydroATLAS attributes pass through from the outlet subbasin; subbasin-scale attributes are re-aggregated as "
        "area-weighted means, elevation extrema, or dominant classes. A dictionary lists every attribute's definition, HydroATLAS theme and unit "
        "codes, and aggregation method. The deposited simplified layer is an August 2026 reconstruction of the original June 2026 layer from the "
        "project's committed workflow and the same source products, validated to match every recorded property of the original (see README for "
        "provenance). All layers are keyed by gage_id for joins to the other HISSS resources (join on the zero-stripped identifier, as the README "
        f"describes); the attribute table and the simplified layer overlap on {F['overlap']:,} gages. Companion resources provide the "
        "streamflow-signature products, harmonized input data, and vegetation and land-cover records.\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aug-dir", default=DEF_AUG)
    ap.add_argument("--full-parquet", default=DEF_FULL)
    ap.add_argument("--metadata", default=os.path.join(DEF_R3, "hisss_gage_metadata.csv"))
    ap.add_argument("--daymet-flags", default=os.path.join(DEF_R3, "hisss_daymet_basin_flags.csv"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--upload", default=None)
    ap.add_argument("--skip-gpkg", action="store_true", help="reuse an existing GeoPackage in --out (README re-render)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    gpkg = os.path.join(a.out, FULL_NAME)
    if a.skip_gpkg:
        assert os.path.exists(gpkg), gpkg
        g = {"n": "(reused)", "wkb_identical": "(reused)", "n_not_multi": "?", "valid": "?"}
    else:
        g = build_gpkg(a.full_parquet, gpkg)
        print(f"GeoPackage written and verified: {g['n']} features, {g['wkb_identical']} geometries byte-identical to the parquet, "
              f"{g['valid']} valid, {g['n_not_multi']} non-multipolygon")
    F = facts(a)
    manifest = []
    for name in UNCHANGED:
        manifest.append(copy(os.path.join(a.aug_dir, name), os.path.join(a.out, name)))
    manifest.append((os.path.basename(a.full_parquet) + " (as GeoPackage)", FULL_NAME))
    readme = render_readme(F, os.path.getsize(gpkg))
    with open(os.path.join(a.out, "hisss_readme_geometry_attributes.md"), "w") as f:
        f.write(readme)
    manifest.append(("(rendered)", "hisss_readme_geometry_attributes.md"))
    with open(os.path.join(a.out, "hisss_abstract_resource4.txt"), "w") as f:
        f.write(render_abstract(F))
    bad = [w for w in ("version 2", "version-2", "previous build", "supersede", "regenerated", "earlier release", "this version")
           if w in readme.lower() or w in render_abstract(F).lower()]
    assert not bad, f"forbidden wording: {bad}"
    with open(os.path.join(a.out, "hisss_rename_manifest_geometry.csv"), "w") as f:
        f.write("original_name,deposit_name\n")
        for o, d in manifest:
            f.write(f"{o},{d}\n")
    rows = []
    for name in sorted(os.listdir(a.out)):
        p = os.path.join(a.out, name)
        if os.path.isfile(p) and not name.startswith("hisss_staging_manifest"):
            rows.append((name, os.path.getsize(p), sha(p)))
    with open(os.path.join(a.out, "hisss_staging_manifest_geometry.csv"), "w") as f:
        f.write("file,bytes,sha256\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]}\n")
    print(f"staged {len(rows)} files in {a.out}")
    print(f"facts: simplified {F['n_simp']:,} / unsimplified {F['n_full']:,} (+{F['n_extra']} large, {F['extra_area_min']:,.0f}–{F['extra_area_max']:,.0f} km²); "
          f"sources simp {F['src_simp']} full {F['src_full']}; low_confidence {F['lc_simp']} / {F['lc_full']} (+{F['lc_extra_ids']}); "
          f"guard {F['n_guard']} (same vertices, n, max Hausdorff deg) {F['guard']}; area shift max {F['area_shift_max_pct']:.2f}% median {F['area_shift_median_pct']:.3f}%; "
          f"attr {F['n_attr']:,}×{len(F['attr_cols'])}; overlap {F['overlap']:,} / attr-only {F['attr_only']} / bound-only {F['bound_only']} "
          f"status {F['bound_only_status']} agency {F['bound_only_agency']}; full overlap {F['full_overlap']:,}; unpadded {len(F['unpadded'])}; "
          f"products with boundary {F['sel_with_simp']} (full {F['sel_with_full']}); exact-join deficit {F['sel_exact_simp']}; Daymet flags {F['n_flags']} (lc {F['flags_lc']})")
    if a.upload:
        os.makedirs(a.upload, exist_ok=True)
        up = []
        for name in ("hisss_readme_geometry_attributes.md", FULL_NAME):
            copy(os.path.join(a.out, name), os.path.join(a.upload, name))
            p = os.path.join(a.upload, name)
            up.append((name, os.path.getsize(p), md5(p)))
        note = os.path.join(os.path.dirname(a.upload.rstrip("/")), os.path.basename(a.upload.rstrip("/")) + "_METADATA_and_MD5.txt")
        with open(note, "w") as f:
            f.write(f"HydroShare Resource 4 (watershed geometry and basin attributes) — update of {dt.date.today().isoformat()}\n")
            f.write("Resource: https://www.hydroshare.org/resource/1283981f7f6a451a91f745902256443e/\n\n")
            f.write("TITLE: keep as is.\n\nABSTRACT (replace):\n" + render_abstract(F) + "\n")
            f.write("KEYWORDS: keep. COVERAGE / LICENSE (CC-BY 4.0): keep.\n\n")
            f.write("FILES TO UPLOAD (from this folder):\n")
            for n, b, m in up:
                action = "REPLACE the existing file" if n.startswith("hisss_readme") else "ADD (new file; DECISION user 2026-10-30)"
                f.write(f"  {n}  {b:,} bytes  md5 {m}  — {action}\n")
            f.write("\nThe other six files of the resource stay as they are (byte-identical to the August staging, verified by MD5 on 2026-10-29).\n")
            f.write("After upload, compare the MD5 HydroShare shows for each file with the values above.\n")
        print(f"upload folder {a.upload} ({len(up)} files) + note {note}")


if __name__ == "__main__":
    main()
