#!/usr/bin/env python
"""Stage the documents (and the rebuilt explorers) of the HydroShare resource
"Hydrologic Information Signatures and Summary Statistics - Vegetation and Land Cover"
(https://www.hydroshare.org/resource/27da62083ba440ebbda8b4a805aa2b1d/).

The three data tables and the exclusions list are deposited byte-identical to the August 2026 staging
(`--aug-dir`); nothing in them depends on the streamflow record. This tool

  1. computes every number the README states from those files, the companion boundary layer and attribute
     table (`--geom-dir`, the unsimplified layer for the 53 large basins) and the two signature products
     (`--products-dir`), and asserts the invariants the README relies on;
  2. renders the README (first-public-version framing: no version labels) and patches the data dictionary
     (the join rule on the three `canon_id` rows and the two exclusion-list rows, which named an internal
     table, `combined_watershed_metadata`, instead of the deposited `hisss_gage_metadata.csv`);
  3. copies the unchanged files and, with `--explorers`, the rebuilt explorers (checked for the keyless
     basemap and the embedded gage count) into `--out`, with a staging manifest;
  4. with `--upload`, writes the upload folder (changed files only, physical copies) and the
     `_METADATA_and_MD5.txt` note beside it.

Run with the daymet venv (pandas, pyarrow):

    ~/HISSS_data/venvs/daymet/bin/python docs/benchmarks/stage_hydroshare_eo_docs.py \
        --out ~/Downloads/Signatures_v3_28oct2026/resource5_eo_products \
        --explorers ~/HISSS_data/eo-explorers-rebuild-30oct2026 \
        --upload ~/Downloads/hisss_vegetation_and_land_cover_upload_2026-10-30
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import sqlite3
import subprocess

import pandas as pd
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

HOME = os.path.expanduser("~")
DEF_AUG = f"{HOME}/Downloads/Signatures/resource5_eo_products"
DEF_GEOM = f"{HOME}/Downloads/Signatures/resource4_geometry_attributes"
DEF_UNSIMP = f"{HOME}/Downloads/Signatures_v3_28oct2026/resource4_geometry_attributes/hisss_watershed_boundaries_unsimplified.gpkg"
DEF_PRODUCTS = f"{HOME}/Downloads/Signatures_v3_28oct2026"
RESOURCE = "https://www.hydroshare.org/resource/27da62083ba440ebbda8b4a805aa2b1d/"
LAI, LULC, NLCD, EXCL = ("hisss_modis_lai_monthly.parquet", "hisss_modis_lulc_annual.parquet",
                         "hisss_nlcd_annual.parquet", "hisss_nlcd_excluded_gages.csv")
DICT, README, ABSTRACT = "hisss_eo_data_dictionary.csv", "hisss_readme_eo.md", "hisss_abstract_resource5.txt"
EXPLORERS = ("hisss_modis_lai_explorer.html", "hisss_modis_lulc_explorer.html", "hisss_nlcd_explorer.html")
BUILDERS = ("build_lai_explorer.py", "build_lulc_explorer.py", "build_nlcd_explorer.py")
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
UNCHANGED = (LAI, LULC, NLCD, EXCL)
LARGE_KM2 = 100000
# dictionary definitions to replace (old -> new); the new text is quoted when it needs to be
DICT_FIXES = {
    "zero-stripped canonical id; joins to combined_watershed_metadata":
        "zero-stripped canonical id (leading zeros removed from all-digit ids); the join key to every other HISSS table, "
        "e.g. hisss_gage_metadata.csv in the harmonized-input resource stores USGS ids in this form",
    "Zero-padded gage identifier of an excluded gage (joins to the signature products)":
        "Agency gage identifier of an excluded gage, as text (join to the other HISSS tables on canon_id)",
    "Zero-stripped canonical identifier (joins to the gage metadata table)":
        "Zero-stripped canonical identifier; the join key to every other HISSS table "
        "(e.g. hisss_gage_metadata.csv in the harmonized-input resource)",
}
DICT_SUBSTRINGS = {  # substring replacements inside a definition (reviewer B and Codex, 2026-10-30)
    "excluded to nlcd_out_of_footprint (see companion CSV)": "excluded (see hisss_nlcd_excluded_gages.csv)",
    # fill (250) is excluded from the NUMERATOR; the denominator is the whole polygon area (nlcd_pipeline.py, QA block)
    "valid pixels / total intersecting pixels (fill=250 excluded from the denominator); ~1 for CONUS basins":
        "valid (non-fill) classified area / total basin polygon area, both on the equal-area grid; ~1 for CONUS basins, "
        "< 1 where the basin extends past the raster footprint",
    "n_nlcd_pixels x 900 m2 (valid classified area; equal-area grid)":
        "n_nlcd_pixels x 900 m2 / 1e6 (valid classified area in km2; equal-area grid)",
    # lai_pipeline.py: the count of composites read for the month (mean over the basin's tiles, rounded); QA is per pixel
    "# QA-good 4-day MCD15A3H composites contributing to the month (month-level)":
        "# 4-day MCD15A3H composites read for the month (mean over the basin's tiles, rounded); per-pixel QA screening "
        "happens within each composite",
}
FORBIDDEN = ("version 2", "version-2", "previous build", "supersede", "regenerated", "earlier release", "this version")


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


def copy(src, dst, physical=False):
    if os.path.exists(dst):
        os.remove(dst)
    if physical:
        shutil.copy2(src, dst)
        return
    r = subprocess.run(["cp", "-p", "-c", src, dst], capture_output=True)
    if r.returncode != 0:
        shutil.copy2(src, dst)


def key(s):
    return s.lstrip("0") if s.isdigit() else s.upper()


def product_ids(path):
    t = pacsv.read_csv(path, convert_options=pacsv.ConvertOptions(include_columns=["gage_id"],
                                                                 column_types={"gage_id": "string"}))
    ids = t.column("gage_id").to_pylist()
    assert len(ids) == len(set(ids)), path
    return ids


def table_facts(path, period_cols):
    md = pq.read_metadata(path)
    t = pq.read_table(path, columns=["gage_id", "canon_id"] + period_cols)
    g = t.column("gage_id").to_pylist()
    c = t.column("canon_id").to_pylist()
    f = {"rows": md.num_rows, "cols": md.num_columns, "gages": len(set(g)), "canon": len(set(c)),
         "gage_ids": set(g), "canon_ids": set(c), "columns": pq.read_schema(path).names}
    assert f["gages"] == f["canon"], path
    keys = t.select(["gage_id"] + period_cols).to_pandas()
    assert not keys.duplicated().any(), f"{path}: duplicate (gage, period) keys"
    periods = keys[period_cols].drop_duplicates()
    assert len(keys) == f["gages"] * len(periods), f"{path}: the gage x period grid is not complete"
    f["periods"] = len(periods)
    for col in period_cols:
        mm = pc.min_max(t.column(col)).as_py()
        f[col] = (mm["min"], mm["max"])
    # the id forms are consistent within a table
    pairs = {(a, b) for a, b in zip(g, c)}
    assert len(pairs) == f["gages"], f"{path}: a gage_id maps to several canon_ids or vice versa"
    assert all(key(a) == b for a, b in pairs), f"{path}: canon_id is not the zero-stripped gage_id"
    return f


def facts(a):
    F = {}
    lai = table_facts(os.path.join(a.aug_dir, LAI), ["year", "month"])
    lulc = table_facts(os.path.join(a.aug_dir, LULC), ["year"])
    nlcd = table_facts(os.path.join(a.aug_dir, NLCD), ["year"])
    F["lai"], F["lulc"], F["nlcd"] = lai, lulc, nlcd
    # LAI: full gage x month grid, nine spatial statistics, all-NA gages
    t = pq.read_table(os.path.join(a.aug_dir, LAI), columns=["gage_id", "year", "month", "lai_mean"]).to_pandas()
    first = int(t.loc[t.year == lai["year"][0], "month"].min())
    last = int(t.loc[t.year == lai["year"][1], "month"].max())
    lai_months = (lai["year"][1] - lai["year"][0]) * 12 + (last - first) + 1
    assert lai["rows"] == lai["gages"] * lai_months == lai["gages"] * lai["periods"], (lai["rows"], lai["gages"], lai_months)
    F["lai_months"], F["lai_first"], F["lai_last"] = lai_months, (lai["year"][0], first), (lai["year"][1], last)
    F["lai_all_na"] = int(t.groupby("gage_id")["lai_mean"].apply(lambda s: s.isna().all()).sum())
    stats = [c for c in lai["columns"] if re.fullmatch(r"lai_(mean|std|min|max|q\d\d)", c)]
    assert len(stats) == 9, stats
    F["lai_stats"] = stats
    assert "__index_level_0__" in lai["columns"]
    del t
    # LULC: 8 schemes x class columns summing to 100
    pct = [c for c in lulc["columns"] if c.endswith("_pct")]
    schemes = sorted({c.split("_c")[0] for c in pct})
    assert lulc["rows"] == lulc["gages"] * (lulc["year"][1] - lulc["year"][0] + 1)
    assert lulc["gage_ids"] == lai["gage_ids"], "LAI and LULC gage sets differ"
    tl = pq.read_table(os.path.join(a.aug_dir, LULC), columns=pct).to_pandas()
    worst = 0.0
    for s in schemes:
        cols = [c for c in pct if c.startswith(s + "_c")]
        tot = tl[cols].sum(axis=1, min_count=1).dropna()
        worst = max(worst, float((tot - 100).abs().max()))
    assert worst < 1e-6, f"a LULC scheme does not sum to 100 (max deviation {worst})"
    F["lulc_pct_cols"], F["lulc_schemes"], F["lulc_sum_dev"] = len(pct), len(schemes), worst
    del tl
    # NLCD: 16 class columns + impervious, CONUS only
    ncls = [c for c in nlcd["columns"] if re.fullmatch(r"nlcd_c\d+_pct", c)]
    assert len(ncls) == 16 and "impervious_mean_pct" in nlcd["columns"], ncls
    assert nlcd["rows"] == nlcd["gages"] * (nlcd["year"][1] - nlcd["year"][0] + 1)
    assert nlcd["gage_ids"] <= lai["gage_ids"], "NLCD gages outside the MODIS universe"
    F["nlcd_classes"] = len(ncls)
    tn = pq.read_table(os.path.join(a.aug_dir, NLCD), columns=ncls + ["impervious_mean_pct"]).to_pandas()
    F["nlcd_sum_dev"] = float((tn[ncls].sum(axis=1) - 100).abs().max())
    assert F["nlcd_sum_dev"] < 1e-6 and tn[ncls].min().min() >= 0 and tn.impervious_mean_pct.between(0, 100).all(), "NLCD class sums or ranges"
    del tn
    F["nlcd_versions"] = sorted(set(pq.read_table(os.path.join(a.aug_dir, NLCD), columns=["nlcd_collection"])
                                    .column("nlcd_collection").to_pylist()))
    # exclusions: US gages of the EO universe absent from NLCD, all in Alaska
    ex = pd.read_csv(os.path.join(a.aug_dir, EXCL), dtype=str)
    assert list(ex.columns) == ["gage_id", "canon_id", "latitude", "longitude", "exclusion_reason"], list(ex.columns)
    assert ex.exclusion_reason.str.contains("Alaska").all() and (ex.latitude.astype(float) > 51).all() and (ex.longitude.astype(float) < -129).all()
    us_missing = {g for g in lai["gage_ids"] - nlcd["gage_ids"] if g.isdigit()}
    assert us_missing == set(ex.gage_id), "exclusion list != EO-universe US gages absent from NLCD"
    F["n_excl"] = len(ex)
    # un-padded ids: EO gage_id vs the agency form in the attribute table
    attr = pd.read_csv(os.path.join(a.geom_dir, "hisss_watershed_hydroatlas_attributes.csv"), usecols=["gage_id"], dtype=str)
    agency = {key(g): g for g in attr.gage_id}
    unp = sorted(g for g in lai["gage_ids"] if g.isdigit() and agency.get(key(g), g) != g)
    assert all(agency[key(g)].startswith("0") and len(agency[key(g)]) in (9, 10) for g in unp), "un-padded ids are not all 9-10 digit ids beginning with 0"
    F["n_unpadded"], F["unpadded"] = len(unp), unp
    # EO universe = the deposited simplified layer; the gages absent from the attribute table are the
    # boundary-only compiles (in the layer, in no product), stored 8-digit in the EO tables
    qa = pd.read_csv(os.path.join(a.geom_dir, "hisss_watershed_boundaries_qa.csv"), usecols=["gage_id", "canon_id", "low_confidence", "area_rel_diff"], dtype={"gage_id": str, "canon_id": str})
    con = sqlite3.connect(os.path.join(a.geom_dir, "hisss_watershed_boundaries.gpkg"))
    layer_ids = {r[0] for r in con.execute("SELECT canon_id FROM watersheds")}
    con.close()
    assert layer_ids == set(qa.canon_id) == lai["canon_ids"], "EO gage set != deposited simplified boundary layer (GeoPackage / QA table)"
    # low_confidence: the LAI table carries the June 2026 layer's geometry flag per row, the LULC table that flag OR
    # low_pixel_support (< 5 MODIS pixels); the explorers show the deposited layer's flag, which differs from the
    # June flag for one basin (04HA001)
    layer = qa.set_index("canon_id")["low_confidence"].astype(bool)
    tl = pq.read_table(os.path.join(a.aug_dir, LAI), columns=["canon_id", "low_confidence"]).to_pandas()
    lai_lc = tl.groupby("canon_id")["low_confidence"].agg(["any", "all"])
    assert (lai_lc["any"] == lai_lc["all"]).all(), "LAI low_confidence varies within a gage"
    tl = pq.read_table(os.path.join(a.aug_dir, LULC), columns=["canon_id", "low_confidence", "low_pixel_support"]).to_pandas()
    lulc_lc = tl.groupby("canon_id")[["low_confidence", "low_pixel_support"]].agg("max")
    lai_diff = sorted(k for k in lai_lc.index if bool(lai_lc.loc[k, "any"]) != bool(layer[k]))
    assert lai_diff == ["04HA001"] and lai_lc.loc["04HA001", "any"] and not layer["04HA001"], lai_diff
    expected = (lai_lc["any"] | lulc_lc["low_pixel_support"]).reindex(lulc_lc.index)
    assert (lulc_lc["low_confidence"] == expected).all(), "LULC low_confidence != LAI flag OR low_pixel_support"
    F["lc_lai"], F["lc_lulc"], F["lc_layer"] = int(lai_lc["any"].sum()), int(lulc_lc["low_confidence"].sum()), int(layer.sum())
    F["lc_tiny"] = int((lulc_lc["low_pixel_support"] & ~lai_lc["any"].reindex(lulc_lc.index)).sum())
    F["lc_diff_rel"] = float(qa.set_index("canon_id").loc["04HA001", "area_rel_diff"])
    assert F["lc_diff_rel"] < 0.5 and F["lc_lulc"] == F["lc_lai"] + F["lc_tiny"]
    not_in_attr = sorted(g for g in lai["gage_ids"] if key(g) not in agency)
    assert all(g.isdigit() and len(g) == 8 and g.startswith("0") for g in not_in_attr), not_in_attr
    F["n_layer"], F["not_in_attr"] = len(qa), not_in_attr
    # the 53 large basins from the unsimplified layer (attributes only)
    con = sqlite3.connect(a.unsimplified)
    large = {key(r[0]) for r in con.execute("SELECT gage_id FROM watersheds WHERE geom_area_km2 > ?", (LARGE_KM2,))}
    n_unsimp = con.execute("SELECT COUNT(*) FROM watersheds").fetchone()[0]
    con.close()
    F["n_large"], F["n_unsimp"] = len(large), n_unsimp
    assert n_unsimp == F["n_layer"] + len(large), (n_unsimp, F["n_layer"], len(large))
    con = sqlite3.connect(a.unsimplified)
    unsimp_ids = {key(r[0]) for r in con.execute("SELECT gage_id FROM watersheds")}
    con.close()
    assert "05KH009" in agency and "05KH009" not in unsimp_ids and "05KH009" not in lai["canon_ids"], "05KH009 placement"
    assert unsimp_ids == lai["canon_ids"] | large, "unsimplified layer != EO universe + the large basins"
    # product coverage
    excl_keys = set(ex.canon_id)
    for tag, rel in (("p1", "resource1_signatures_wy1993_2025/hisss_signatures_wy1993_2025.csv"),
                     ("p2", "resource2_signatures_wy1980_2025/hisss_signatures_wy1980_2025.csv")):
        ids = product_ids(os.path.join(a.products_dir, rel))
        ks = {key(i) for i in ids}
        modis = len(ks & lai["canon_ids"]); modis_exact = len(set(ids) & lai["gage_ids"])
        missing = ks - lai["canon_ids"]
        assert missing <= large, f"{tag}: product gages without MODIS not all among the large basins"
        nl = len(ks & nlcd["canon_ids"]); nl_exact = len(set(ids) & nlcd["gage_ids"])
        nmiss = ks - nlcd["canon_ids"]
        other = {k for k in nmiss if k.isdigit() and k not in excl_keys and k not in large}
        assert not other, f"{tag}: {len(other)} US gages without NLCD that are neither Alaskan nor large basins"
        F[tag] = {"n": len(ids), "modis": modis, "modis_deficit": modis - modis_exact, "missing": len(missing),
                  "nlcd": nl, "nlcd_deficit": nl - nl_exact}
    assert F["p1"]["modis_deficit"] == F["p1"]["nlcd_deficit"] and F["p2"]["modis_deficit"] == F["p2"]["nlcd_deficit"]
    # dictionary coverage
    d = pd.read_csv(os.path.join(a.aug_dir, DICT), dtype=str)
    assert list(d.columns) == ["table", "column", "units", "definition"]
    for name, f in ((LAI, lai), (LULC, lulc), (NLCD, nlcd)):
        rows = set(d.loc[d.table == name, "column"])
        caveats = {r for r in rows if r.startswith("CAVEAT_")}
        assert rows - caveats == set(f["columns"]), f"dictionary rows != columns for {name}"
    assert set(d.loc[d.table == EXCL, "column"]) == set(ex.columns)
    F["dict_rows"], F["dict_caveats"] = len(d), int(d.column.str.startswith("CAVEAT_").sum())
    # explorers
    if True:
        F["explorers"] = {}
        for name, builder, tf in zip(EXPLORERS, BUILDERS, (lai, lulc, nlcd)):
            p = os.path.join(a.explorers, name)
            txt = open(p, encoding="utf-8").read()
            data = embedded_json(txt, os.path.join(a.repo, "EO_data_processing", "viz", builder))
            ids = data["__POINTS__"]["gage_id"]
            assert "basemaps.cartocdn.com" not in txt and "arcgisonline" in txt and "openstreetmap.org" in txt, name
            assert len(ids) == tf["gages"] and set(ids) == tf["gage_ids"], (name, len(ids), tf["gages"])
            F["explorers"][name] = {"gages": len(ids), "vars": len(data["__VARS__"]), "bytes": os.path.getsize(p)}
        assert F["explorers"][EXPLORERS[1]]["vars"] == 27 and F["explorers"][EXPLORERS[2]]["vars"] == 32, F["explorers"]
    return F


def embedded_json(html, builder):
    """The JSON blobs a builder substituted into its TEMPLATE (`__POINTS__`, `__VARS__`, ...), located by the
    template text around each marker."""
    tpl = re.search(r'TEMPLATE\s*=\s*r"""(.*?)"""', open(builder, encoding="utf-8").read(), re.S).group(1)
    markers = [m.group(0) for m in re.finditer(r"__[A-Z_]+__", tpl)]
    segs = re.split("|".join(map(re.escape, dict.fromkeys(markers))), tpl)
    out, pos = {}, 0
    for i, name in enumerate(markers):
        a_end, b_start = segs[i][-60:], segs[i + 1][:60]
        s = html.index(a_end, pos) + len(a_end); e = html.index(b_start, s)
        out[name] = json.loads(html[s:e]); pos = e
    return out


def render_readme(F):
    lai, lulc, nlcd = F["lai"], F["lulc"], F["nlcd"]
    p1, p2 = F["p1"], F["p2"]
    mon = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    lai_span = f"{mon[F['lai_first'][1]-1]} {F['lai_first'][0]} – {mon[F['lai_last'][1]-1]} {F['lai_last'][0]}"
    nl_years = nlcd["year"][1] - nlcd["year"][0] + 1
    lu_years = lulc["year"][1] - lulc["year"][0] + 1
    n54 = F["n_large"] + 1
    return f"""# HISSS Vegetation and Land-Cover Records

Part of the **Hydrologic Information Signatures and Summary Statistics (HISSS)**
dataset. This resource contains three per-watershed Earth-observation products —
monthly MODIS leaf area index (LAI), annual MODIS land use / land cover (LULC), and
annual NLCD land cover + impervious surface — each aggregated over the full upstream
watershed of the HISSS gages and keyed to the same gage identifiers as the signature
products. All three tables carry both `gage_id` (the agency station identifier) and
`canon_id` (the same identifier with leading zeros stripped). **Join on the stripped
form** (see *Reading and joining*): {F['n_unpadded']} long USGS site numbers — 9–10 digits beginning
with a zero — are stored in `gage_id` WITHOUT their leading zero, so an exact-string
join on `gage_id` silently drops whichever of those {F['n_unpadded']} gages the other table stores padded
({p1['modis_deficit']} / {p2['modis_deficit']} gages against the two signature products); `canon_id`, or stripping
leading zeros on both sides, matches every gage's identifier.

## Files

| File | Description |
|---|---|
| `{LAI}` | Monthly basin LAI: {lai['rows']:,} rows = {lai['gages']:,} watersheds × {F['lai_months']} months ({lai_span}), from MODIS MCD15A3H v061 (500 m, 4-day composites). Nine spatial statistics of the basin's monthly-mean LAI field (`lai_mean`, `lai_std`, min/max, 5 quantiles) plus per-month QA fractions |
| `{LULC}` | Annual basin land-cover composition: {lulc['rows']:,} rows = {lulc['gages']:,} watersheds × {lu_years} years ({lulc['year'][0]}–{lulc['year'][1]}), from MODIS MCD12Q1 v061 (500 m). Per-class % coverage for all {F['lulc_schemes']} MODIS classification schemes ({F['lulc_pct_cols']} class columns; each scheme sums to 100 per row) |
| `{NLCD}` | Annual basin land cover + imperviousness, CONUS only: {nlcd['rows']:,} rows = {nlcd['gages']:,} watersheds × {nl_years} years ({nlcd['year'][0]}–{nlcd['year'][1]}), from USGS/MRLC Annual NLCD Collection 1 ({', '.join(F['nlcd_versions'])}, 30 m). {F['nlcd_classes']} land-cover class % columns (sum to 100) plus basin-mean fractional impervious surface, with per-row QA flags |
| `{EXCL}` | The {F['n_excl']} Alaska gages excluded from the NLCD table (outside the CONUS Annual NLCD footprint), with coordinates and reason |
| `{DICT}` | Definition and units for every column of the four tables above (including the NLCD interpretation caveats, recorded as `CAVEAT_*` rows) |
| `{EXPLORERS[0]}` | Interactive map explorer for the LAI product: {lai['gages']:,} gage points colored by QA/summary variables; click a gage for its full monthly LAI time series with quantile bands and a per-month coverage strip |
| `{EXPLORERS[1]}` | Interactive map explorer for the MODIS LULC product: {F['explorers'][EXPLORERS[1]]['vars']} map variables; click a gage for its stacked-area land-cover composition {lulc['year'][0]}–{lulc['year'][1]} with a switcher across all {F['lulc_schemes']} classification schemes |
| `{EXPLORERS[2]}` | Interactive map explorer for the NLCD product: {F['explorers'][EXPLORERS[2]]['vars']} map variables; click a gage for its stacked-area composition {nlcd['year'][0]}–{nlcd['year'][1]} with the impervious-surface trace |

The explorers are self-contained HTML files (all data embedded) — download and
double-click to open. They need an internet connection only to load the Leaflet map
library and the basemap tiles (an Esri light-grey canvas by default, OpenStreetMap as
a fallback; the layers control under the zoom buttons in the map's top-left corner
switches between them).

## Reading and joining

```python
import pandas as pd
lai  = pd.read_parquet("{LAI}")
lulc = pd.read_parquet("{LULC}")
nlcd = pd.read_parquet("{NLCD}")
# join to the HISSS signature products on the zero-stripped key (NOT on gage_id directly —
# {F['n_unpadded']} long USGS ids are stored un-padded in these tables):
key  = lambda s: s.lstrip("0") if s.isdigit() else s.upper()
sig  = pd.read_csv("hisss_signatures_wy1993_2025.csv", dtype={{"gage_id": str}})
sig["key"] = sig.gage_id.map(key)              # lulc["canon_id"] is already the stripped key
lulc2020 = lulc[lulc.year == 2020].drop(columns=["gage_id"])   # keep a single gage_id column
merged = sig.merge(lulc2020, left_on="key", right_on="canon_id", how="left")
```

```r
key  <- function(x) ifelse(grepl("^[0-9]+$", x), sub("^0+", "", x), toupper(x))
sig  <- readr::read_csv("hisss_signatures_wy1993_2025.csv",
                        col_types = readr::cols(gage_id = readr::col_character()))
lulc <- arrow::read_parquet("{LULC}")
sig$key <- key(sig$gage_id)
merged <- dplyr::left_join(sig, dplyr::filter(lulc, year == 2020), by = c("key" = "canon_id"))
```

**Join rule (all HISSS tables).** Read `gage_id` as a character string, never
numeric. Join on the **zero-stripped** form — strip leading zeros from all-digit ids,
leave HYDAT ids unchanged — and do so on *both* sides of every join; HISSS station ids are
unique under this rule (17,166 distinct stations, zero collisions). **Do not re-pad**:
padding to 8 digits cannot recover 9- and 10-digit USGS site numbers that begin with a
zero (`021650905`, `0212433550`), and some HISSS tables store exactly those ids
un-padded. In Python: `key = lambda s: s.lstrip("0") if s.isdigit() else s.upper()`;
in R: `key <- function(x) ifelse(grepl("^[0-9]+$", x), sub("^0+", "", x), toupper(x))`.
The un-padded ids in these tables are an artefact of an 8-digit padding rule in the
processing pipeline (found 2026-09-04); the files are deposited byte-identical to their
validated builds rather than rewritten. The gage-metadata table of the companion
harmonized-input resource (`hisss_gage_metadata.csv`) stores every USGS id in the
zero-stripped form, so it joins to `canon_id` directly.

## Product notes

**MODIS LAI (monthly).** Two-stage aggregation: per pixel, a day-weighted mean of the
month's quality-screened 4-day composites gives a monthly-mean LAI field; the table's
nine statistics then describe the spatial distribution of that field across the
watershed's pixels (coverage-weighted). So `lai_mean` is the basin mean monthly LAI and
`lai_std` / the quantiles measure spatial heterogeneity. `lai_min`/`lai_max` are
unweighted extrema over any intersecting pixel — a single boundary sliver can drive
them; prefer `lai_q05`/`lai_q95` for a robust range. `good_coverage_frac` is the
share of expected pixel-observations that were valid; `partial_month` flags months
with fewer than 6 composites. {F['lai_all_na']} heavily urban basins are permanently all-NA (MODIS
retrieves no LAI over built-up land); they can be listed as the gages whose `lai_mean`
is NA in every month. High-latitude winter NAs reflect snow/darkness and are expected.

**MODIS LULC (annual).** All {F['lulc_schemes']} MCD12Q1 classification schemes are retained (IGBP, UMD,
LAI/fPAR biome, BGC, PFT, and the 3 FAO-LCCS property layers); column names follow
`{{scheme}}_c{{code}}_pct`. Within each scheme the class percentages sum to 100 over valid
pixels. Interpretation caveats: MODIS labels much of the southeastern-US forest as
woody savanna; the data provider flags reduced land-cover confidence for 2021 onward.

**Annual NLCD (annual, CONUS only).** USGS/MRLC Annual NLCD Collection 1, version
{', '.join(F['nlcd_versions'])}, at 30 m — roughly {round((500 / 30) ** 2, -1):.0f}× the pixel density of the MODIS products and a
{lulc['year'][0] - nlcd['year'][0]}-year longer look-back than the MODIS land-cover record, but limited to the conterminous United States. Includes the {F['nlcd_classes']}-class
land-cover composition and `impervious_mean_pct`, the basin-mean fractional impervious
surface (the continuous urbanization signal; developed-class % and impervious % are
related but not additive). The {F['n_excl']} Alaska gages in the HISSS watershed universe fall
outside the product footprint and are excluded rather than reported as false zeros
(see `{EXCL}`). Interpretation caveats (also in the
dictionary): year-to-year class flicker — most visibly a shrub↔grassland swap
artifact — means single-year changes should not be over-read (compare multi-year
means); the most recent one to two years of any collection come from a lighter-weight
annual-update pathway (a methodological seam at the series tail); and re-releases of
the collection reprocess historical years, so values are tied to the {', '.join(F['nlcd_versions'])} version
recorded in the `nlcd_collection` column.

**Watershed universe and coverage.** The three products were computed on the {F['n_layer']:,}
HISSS watershed polygons of the June 2026 build, deposited as `hisss_watershed_boundaries.gpkg`
in the companion watershed-geometry resource (the deposited layer is the August 2026
reconstruction of that build; it reproduces the build's recorded properties, and the one
QA flag that differs is described under *Provenance* — see also its README). That
layer leaves out the {F['n_large']} basins larger than {LARGE_KM2:,} km² and one mis-delineated HydroBASINS
fallback basin (`05KH009`), so those {n54} gages have no record in any of the three tables.
The geometry resource also holds an unsimplified layer that includes the {F['n_large']} large basins;
no EO product was computed on it. Of the signature-product gages, {p1['modis']:,} of {p1['n']:,}
(WY 1993–2025 product) and {p2['modis']:,} of {p2['n']:,} (WY 1980–2025 product) have MODIS
LAI/LULC records — every missing gage is one of the {F['n_large']} large basins; {p1['nlcd']:,} and {p2['nlcd']:,}
respectively have NLCD records (the remainder are Canadian, Alaskan, or among the large
basins). These counts use the zero-stripped join; an exact-string join on `gage_id`
returns {p1['modis_deficit']} / {p2['modis_deficit']} fewer because of the {F['n_unpadded']} un-padded ids described above.

**A note on one column.** `{LAI}` carries a vestigial
`__index_level_0__` column (an internal pandas row index retained at write); it holds
no information and should be ignored. The file is deposited byte-identical to the
validated June 2026 build rather than rewritten.

## Provenance

The MODIS products were computed in June 2026 and the NLCD product in July 2026, each
passing automated validation gates (per-scheme sum-to-100, value ranges, key
uniqueness, row completeness) and an independent adversarial code/results review. The
deposited files were retrieved from the project archive backup and re-verified against
the documented row counts, watershed counts, and validation invariants before deposit.
`{EXCL}` is an exact reconstruction of the original exclusion
record (the excluded set is precisely the EO-universe US gages absent from the NLCD
table; all {F['n_excl']} verified to lie in Alaska). The per-granule MODIS input manifest from
the original runs was not retained in the backup; granule-level provenance is
documented in the HISSS code repository instead. The three explorers were built in
October 2026 from the deposited tables, with gage coordinates and polygon flags taken
from the deposited boundary layer. The explorers' `low_confidence` map variable is that
layer's flag ({F['lc_layer']} basins: HydroBASINS fallbacks and area outliers). The tables carry their
own `low_confidence` columns, built with the June 2026 layer: the LAI table's is the June
geometry flag ({F['lc_lai']} basins), the LULC table's that flag or `low_pixel_support` ({F['lc_lulc']} basins),
the NLCD table's the union the dictionary states. The June and deposited layers disagree
on one basin, `04HA001`: the June build flagged it (so the LAI and LULC tables still do);
in the deposited layer its published drainage area differs from the polygon area by {F['lc_diff_rel']*100:.0f} %,
within the 50 % rule of the layer's `area_flag`, so the explorers do not. The explorers'
gage search uses the tables' id form (the {F['n_unpadded']} long USGS ids without their leading zero).

## License and citation

License: CC-BY 4.0. Please cite the HISSS data paper and this resource (DOI: to be
assigned at publication). Source data credits: MODIS MCD15A3H v061 (LAI) and MCD12Q1
v061 (land cover), NASA LP DAAC; Annual NLCD Collection 1 ({', '.join(F['nlcd_versions'])}), USGS / Multi-
Resolution Land Characteristics Consortium; watershed boundaries from GAGES-II
(Falcone 2011), the Water Survey of Canada Hydrometric Network Basin Polygons, and
HydroBASINS (Lehner and Grill 2013).

This resource is part of the HISSS collection on HydroShare:
https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/
(collection DOI upon publication: https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768).
"""


def patch_dictionary(src, dst):
    import csv, io
    rows = list(csv.reader(io.StringIO(open(src, encoding="utf-8").read().replace("\r\n", "\n"))))
    assert rows[0] == ["table", "column", "units", "definition"] and all(len(r) == 4 for r in rows), "dictionary shape"
    hits = {old: 0 for old in DICT_FIXES}
    sub_hits = {old: 0 for old in DICT_SUBSTRINGS}
    for r in rows[1:]:
        if r[3] in DICT_FIXES:
            hits[r[3]] += 1; r[3] = DICT_FIXES[r[3]]
        for old, new in DICT_SUBSTRINGS.items():
            if old in r[3]:
                r[3] = r[3].replace(old, new); sub_hits[old] += 1
    assert hits[list(DICT_FIXES)[0]] == 3 and all(v == 1 for v in list(hits.values())[1:]), hits
    assert all(v == 1 for v in sub_hits.values()), sub_hits
    buf = io.StringIO(); w = csv.writer(buf, lineterminator="\n"); w.writerows(rows); text = buf.getvalue()
    assert "combined_watershed_metadata" not in text and "nlcd_out_of_footprint" not in text
    with open(dst, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    d = pd.read_csv(dst, dtype=str)
    assert list(d.columns) == ["table", "column", "units", "definition"] and d.definition.notna().all() and len(d) == len(rows) - 1
    return len(d)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aug-dir", default=DEF_AUG)
    ap.add_argument("--geom-dir", default=DEF_GEOM)
    ap.add_argument("--unsimplified", default=DEF_UNSIMP)
    ap.add_argument("--products-dir", default=DEF_PRODUCTS)
    ap.add_argument("--explorers", required=True, help="folder with the three rebuilt explorer HTML files (checked: keyless basemap, embedded gage set, map-variable counts)")
    ap.add_argument("--repo", default=REPO, help="repository root (the explorer builders supply the template anchors)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--upload", default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    F = facts(a)
    readme = render_readme(F)
    abstract = open(os.path.join(a.aug_dir, ABSTRACT), encoding="utf-8").read()
    with open(os.path.join(a.out, README), "w", encoding="utf-8", newline="\n") as f:
        f.write(readme)
    n_dict = patch_dictionary(os.path.join(a.aug_dir, DICT), os.path.join(a.out, DICT))
    with open(os.path.join(a.out, ABSTRACT), "w", encoding="utf-8", newline="\n") as f:
        f.write(abstract)
    dtext = open(os.path.join(a.out, DICT), encoding="utf-8").read().lower()
    bad = [w for w in FORBIDDEN if w in readme.lower() or w in abstract.lower() or w in dtext]
    assert not bad, f"forbidden wording: {bad}"
    status = {README: "replaced (rendered)", DICT: "replaced (patched)", ABSTRACT: "unchanged"}
    for name in UNCHANGED:
        copy(os.path.join(a.aug_dir, name), os.path.join(a.out, name))
        status[name] = "unchanged"
    for name in EXPLORERS:
        copy(os.path.join(a.explorers, name), os.path.join(a.out, name))
        status[name] = "rebuilt"
    rows = []
    for name in sorted(os.listdir(a.out)):
        p = os.path.join(a.out, name)
        if os.path.isfile(p) and not name.startswith("hisss_staging_manifest"):
            rows.append((name, os.path.getsize(p), md5(p), sha(p), status.get(name, "")))
    with open(os.path.join(a.out, "hisss_staging_manifest_eo.csv"), "w", newline="\n") as f:
        f.write("file,bytes,md5,sha256,status\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    lai, lulc, nlcd, p1, p2 = F["lai"], F["lulc"], F["nlcd"], F["p1"], F["p2"]
    print(f"staged {a.out}: {len(rows)} files; README {len(readme.encode()):,} B; dictionary {n_dict} rows")
    print(f"facts: LAI {lai['rows']:,} rows {lai['gages']:,} gages {F['lai_months']} months {F['lai_first']}..{F['lai_last']} all-NA {F['lai_all_na']}; "
          f"LULC {lulc['rows']:,} rows {lulc['year']} {F['lulc_pct_cols']} pct cols {F['lulc_schemes']} schemes (max sum dev {F['lulc_sum_dev']:.3g}); "
          f"NLCD {nlcd['rows']:,} rows {nlcd['gages']:,} gages {nlcd['year']} {F['nlcd_classes']} classes {F['nlcd_versions']}; excl {F['n_excl']}; "
          f"unpadded {F['n_unpadded']}; layer {F['n_layer']:,} (+{F['n_large']} large = {F['n_unsimp']:,}); "
          f"P1 modis {p1['modis']:,}/{p1['n']:,} (exact −{p1['modis_deficit']}, missing {p1['missing']}) nlcd {p1['nlcd']:,}; "
          f"P2 modis {p2['modis']:,}/{p2['n']:,} (exact −{p2['modis_deficit']}, missing {p2['missing']}) nlcd {p2['nlcd']:,}; "
          f"dictionary {F['dict_rows']} rows ({F['dict_caveats']} caveats); explorers {F['explorers']}")
    if a.upload:
        os.makedirs(a.upload, exist_ok=True)
        names = [README, DICT] + list(EXPLORERS)
        up = []
        for name in names:
            copy(os.path.join(a.out, name), os.path.join(a.upload, name), physical=True)
            p = os.path.join(a.upload, name)
            up.append((name, os.path.getsize(p), md5(p)))
        note = os.path.join(os.path.dirname(a.upload.rstrip("/")), os.path.basename(a.upload.rstrip("/")) + "_METADATA_and_MD5.txt")
        keep = [(n, os.path.getsize(os.path.join(a.out, n)), md5(os.path.join(a.out, n))) for n in UNCHANGED]
        with open(note, "w", encoding="utf-8") as f:
            f.write(f"HydroShare resource \"Hydrologic Information Signatures and Summary Statistics - Vegetation and Land Cover\" — update of {dt.date.today().isoformat()}\n")
            f.write(f"Resource: {RESOURCE}\n\n")
            f.write("TITLE: keep as is.\n\nABSTRACT: keep if the live text equals the one below; otherwise replace with it:\n" + abstract.strip() + "\n\n")
            f.write("KEYWORDS: keep. COVERAGE / LICENSE (CC-BY 4.0): keep.\n\n")
            f.write("FILES TO UPLOAD (from this folder):\n")
            for n, b, m in up:
                action = "REPLACE the existing file (ADD if it is no longer in the resource)" if n == README else "REPLACE the existing file"
                f.write(f"  {n}  {b:,} bytes  md5 {m}  — {action}\n")
            f.write("\nFILES THAT STAY AS THEY ARE (byte-identical to the August 2026 staging) — confirm each MD5 on the resource page:\n")
            for n, b, m in keep:
                f.write(f"  {n}  {b:,} bytes  md5 {m}\n")
            f.write("\nFILES THAT DO NOT BELONG TO THIS RESOURCE (uploaded by mistake 2026-10-30, removed the same day):\n"
                    "  hisss_readme_geometry_attributes.md, hisss_watershed_boundaries_unsimplified.gpkg — confirm both are gone.\n")
            f.write("\nAfter upload, compare the MD5 HydroShare shows for each file with the values above.\n")
        print(f"upload folder {a.upload} ({len(up)} files) + note {note}")


if __name__ == "__main__":
    main()
