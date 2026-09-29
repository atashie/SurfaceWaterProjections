"""Build the Daymet weight matrix for the watershed polygons (run once per polygon layer).

Each polygon is reprojected to the Daymet Lambert Conformal Conic grid (the raster is
never resampled) and exactextract returns, for every 1 km cell the polygon touches, the
cell's flat index (row * 7814 + col, row 0 = north) and the fraction of the cell the
polygon covers. Each weight is also stored as the TRUE area (km2) of the covered cell
piece: coverage x cell area on the ellipsoid (the LCC grid is conformal, so its cells
are 1 km2 only in projected units). daymet_aggregate.py then computes
sum(w * v) / sum(w) over the cells valid that day (fill -9999 excluded, weights
renormalised). Area weights are its default: they reproduce the co-authors' gdptools
product (2023 validation, 2026-09-29); coverage weights give exactextract's plain "mean".

Usage:
    python daymet_weights.py --polygons <boundaries.parquet|gpkg> --grid-file <Daymet NA .nc>
        --streamflow <streamflow parquet with gage_id> --out-dir <dir>

Outputs in --out-dir:
    daymet_weights.parquet    basin_idx int32, cell_id int64, coverage float64, area_km2 float64
    daymet_basins.csv         one row per polygon: basin_idx, site_id (the streamflow
                              parquet's spelling of the id -- the climate join key),
                              boundary_gage_id, canon_id, polygon attributes, n_cells,
                              n_cells_full, coverage_sum (projected cells), area_sum_km2
                              (true area of the weights), area_ratio (area_sum_km2 /
                              geom_area_km2), low_pixel_support (coverage_sum < 4)
    daymet_weights_meta.json  polygon file + md5, grid and CRS, software, timings, checks

Checks (fatal): every polygon gets >= 1 cell; coverage in (0, 1]; each polygon's
representative point falls in one of its own cells (catches row/column or order
errors); site ids map one-to-one onto the streamflow ids. area_ratio is reported: the
weights' true area against the polygon's own equal-area (EPSG:6933) area.
"""
import argparse
import os
import sys
import time

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import (NX, NY, canon, cell_true_area_km2, md5sum, read_grid,  # noqa: E402
                           software_versions, write_json)

LOW_PIXEL_SUPPORT = 4.0   # coverage_sum below this -> low_pixel_support (plan section 5.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons", required=True)
    ap.add_argument("--grid-file", required=True, help="any Daymet NA daily file (grid + CRS source)")
    ap.add_argument("--streamflow", required=True, help="streamflow parquet (gage_id column) for the id spelling")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    t_start = time.time()

    grid = read_grid(a.grid_file)
    g = gpd.read_parquet(a.polygons) if a.polygons.endswith(".parquet") else gpd.read_file(a.polygons)
    g = g.reset_index(drop=True)
    n = len(g)
    if g.crs is None:
        sys.exit("polygon layer has no CRS")
    if not g.geometry.is_valid.all() or g.geometry.is_empty.any():
        sys.exit("invalid or empty polygons in the layer")
    print(f"{n:,} polygons ({a.polygons}), CRS {g.crs.to_string()}; grid from {os.path.basename(a.grid_file)}")

    # --- site_id = the streamflow parquet's spelling (the Julia runner joins on the exact string)
    sf_ids = duckdb.sql(f"SELECT DISTINCT gage_id FROM read_parquet('{a.streamflow}')").df()["gage_id"]
    sf_map = {}
    for s in sf_ids:
        c = canon(s)
        if c in sf_map and sf_map[c] != s:
            sys.exit(f"streamflow ids collide on canon id {c}: {sf_map[c]} vs {s}")
        sf_map[c] = s
    g["canon_id"] = g["gage_id"].map(canon)
    if not g["canon_id"].is_unique:
        sys.exit("polygon canon ids are not unique")
    g["site_id"] = g["canon_id"].map(sf_map)
    unmatched = g["site_id"].isna()
    print(f"id mapping: {int((~unmatched).sum()):,} matched to streamflow ids, "
          f"{int(unmatched.sum())} unmatched (kept under their boundary id), "
          f"{int((g['site_id'].notna() & (g['site_id'] != g['gage_id'])).sum())} respelled")
    g.loc[unmatched, "site_id"] = g.loc[unmatched, "gage_id"]

    # --- reproject polygons to the Daymet grid CRS; build a dummy raster with the grid geometry
    t0 = time.time()
    gl = g[["geometry"]].to_crs(grid.crs)
    xmin, ymin, xmax, ymax = grid.bounds
    rast = NumPyRasterSource(np.ones((NY, NX), dtype=np.uint8), xmin, ymin, xmax, ymax,
                             srs_wkt=grid.crs.to_wkt())
    res = exact_extract(rast, gl, ["cell_id", "coverage"], output="pandas")
    t_ee = time.time() - t0
    if len(res) != n:
        sys.exit(f"exactextract returned {len(res)} rows for {n} polygons")
    counts = res["cell_id"].map(len).to_numpy()
    if (counts == 0).any():
        sys.exit(f"{int((counts == 0).sum())} polygons got no cells: "
                 f"{g.loc[counts == 0, 'gage_id'].tolist()[:20]}")
    basin_idx = np.repeat(np.arange(n, dtype=np.int32), counts)
    cell_id = np.concatenate(res["cell_id"].to_list()).astype(np.int64)
    cov = np.concatenate(res["coverage"].to_list()).astype(np.float64)
    if not ((cov > 0) & (cov <= 1 + 1e-9)).all():
        sys.exit("coverage outside (0, 1]")
    if cell_id.min() < 0 or cell_id.max() >= NX * NY:
        sys.exit("cell_id outside the grid")
    print(f"exactextract: {len(cov):,} (basin, cell) weights in {t_ee:.1f} s")

    # --- order/orientation check: each polygon's representative point lies in one of its cells
    rp = gl.geometry.representative_point()
    col = np.floor((rp.x.to_numpy() - xmin) / 1000.0).astype(np.int64)
    row = np.floor((ymax - rp.y.to_numpy()) / 1000.0).astype(np.int64)
    rp_cell = row * NX + col
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    ok = np.array([rp_cell[i] in set(cell_id[starts[i]:starts[i] + counts[i]]) for i in range(n)])
    if not ok.all():
        sys.exit(f"representative point outside its own cells for {int((~ok).sum())} polygons")
    print("orientation check: every representative point falls in one of its polygon's cells")

    # --- true area of each covered cell piece
    t0 = time.time()
    ucell, inv = np.unique(cell_id, return_inverse=True)
    area = cov * cell_true_area_km2(ucell, grid.crs)[inv]
    print(f"true cell areas for {len(ucell):,} distinct cells in {time.time() - t0:.1f} s")

    # --- per-basin QA
    cov_sum = np.bincount(basin_idx, weights=cov, minlength=n)
    area_sum = np.bincount(basin_idx, weights=area, minlength=n)
    qa = pd.DataFrame({
        "basin_idx": np.arange(n, dtype=np.int32), "site_id": g["site_id"],
        "boundary_gage_id": g["gage_id"], "canon_id": g["canon_id"],
        "in_streamflow": ~unmatched.to_numpy(),
        "watershed_geom_source": g.get("watershed_geom_source"), "gage_type": g.get("gage_type"),
        "geom_area_km2": g["geom_area_km2"], "geom_simplified": g.get("geom_simplified"),
        "low_confidence": g.get("low_confidence"),
        "n_cells": counts, "n_cells_full": np.bincount(basin_idx, weights=(cov >= 0.999), minlength=n).astype(int),
        "coverage_sum": cov_sum, "area_sum_km2": area_sum,
        "area_ratio": area_sum / g["geom_area_km2"].to_numpy(),
        "low_pixel_support": cov_sum < LOW_PIXEL_SUPPORT,
    })
    qa.to_csv(os.path.join(a.out_dir, "daymet_basins.csv"), index=False)
    r = qa["area_ratio"]
    print(f"weights' true area / polygon area: median {r.median():.5f}, "
          f"p01 {r.quantile(0.01):.4f}, p99 {r.quantile(0.99):.4f}, min {r.min():.4f}, max {r.max():.4f}")
    print(f"low_pixel_support (< {LOW_PIXEL_SUPPORT:g} cells): {int(qa['low_pixel_support'].sum())}")

    wpath = os.path.join(a.out_dir, "daymet_weights.parquet")
    pq.write_table(pa.table({"basin_idx": basin_idx, "cell_id": cell_id, "coverage": cov, "area_km2": area}),
                   wpath, compression="zstd")
    write_json(os.path.join(a.out_dir, "daymet_weights_meta.json"), {
        "polygons": os.path.abspath(a.polygons), "polygons_md5": md5sum(a.polygons),
        "n_polygons": n, "n_weights": int(len(cov)), "n_distinct_cells": int(len(ucell)),
        "grid_file": os.path.basename(a.grid_file), "grid_bounds": grid.bounds,
        "grid_shape": [NY, NX], "crs_wkt": grid.crs.to_wkt(),
        "low_pixel_support_threshold": LOW_PIXEL_SUPPORT,
        "area_ratio_quantiles": {str(q): float(v) for q, v in r.quantile([0, .01, .5, .99, 1]).items()},
        "id_mapping": {"matched": int((~unmatched).sum()), "unmatched": int(unmatched.sum()),
                       "unmatched_ids": g.loc[unmatched, "gage_id"].tolist()},
        "seconds_exactextract": round(t_ee, 1), "seconds_total": round(time.time() - t_start, 1),
        "weights_md5": md5sum(wpath), "software": software_versions(),
    })
    print(f"wrote {wpath} ({os.path.getsize(wpath) / 1e6:.0f} MB) in {time.time() - t_start:.0f} s total")


if __name__ == "__main__":
    main()
