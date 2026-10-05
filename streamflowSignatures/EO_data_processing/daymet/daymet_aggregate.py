"""Aggregate one Daymet NA variable-year file to per-basin daily area-weighted means.

Usage:
    python daymet_aggregate.py --file <daymet_v4_daily_na_<var>_<yyyy>.nc> --weights-dir <dir>
        --out-dir <dir> [--workers 8] [--weight area|coverage] [--source-sha256 <hex>]
        [--allow-out-of-range]

Method (chunk-aligned; each HDF5 chunk is decompressed exactly once):
    for every spatial tile of the file's own chunk grid that holds a weighted cell,
    read the tile for a chunk-aligned block of days, keep the weighted cells, and
    accumulate   num[basin, day] += sum(w * v)   and   den[basin, day] += sum(w)
    over the cells whose value is valid that day (not -9999, finite).
    mean = num / den; NaN where no weighted cell is valid (den == 0).
Weights (daymet_weights.py): --weight area (default) = true area of each covered cell piece
(coverage x cell area on the ellipsoid), the area-weighted mean the co-authors' gdptools
run produced; --weight coverage = coverage fraction alone, exactextract's plain "mean".
Either way fill is excluded and the weights are renormalised per day.

Outputs in --out-dir:
    <var>_<year>.parquet      site_id str, Date date32, <var> float64; 365 rows per basin
    <var>_<year>_qa.csv       per basin: valid_frac_min / _max (valid weight / total weight
                              over the year), n_nan_days, min, max, mean of the daily means
    <var>_<year>_timing.json  wall and CPU seconds, peak RSS, tiles read, bytes decompressed,
                              out-of-range counts, input file size and mtime, the input's
                              SHA-256 as verified by the caller (--source-sha256), the
                              weights md5 and the git state of these tools
    <var>_<year>.done         written last; carries the same SHA-256, weights md5 and commit,
                              so daymet_stream.py can refuse a .done built on other weights
The run is resumable at the file level: an existing <var>_<year>.done skips the file.
A cell value outside RANGES (daymet_common.py) is fatal unless --allow-out-of-range: the
outputs are written for inspection, but no .done.
Deterministic: partial sums are accumulated in a fixed task order, so for a FIXED weights
file, reruns and different --workers values give bit-identical output. A different weights
file (e.g. more basins) changes the tile order and so the last bits (<= 2e-13 relative).
"""
import argparse
import json
import os
import resource
import sys
import time
from multiprocessing import get_context

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import FILL, NX, NY, RANGES, git_state, read_grid, utc_now, write_json  # noqa: E402

READ_TARGET_BYTES = 160 * 2**20   # uncompressed bytes per read, rounded to whole time chunks
COMPUTE_DAYS = 64                 # days per sparse product (bounds float64 temporaries)

_F = None
_D = None


def _init(path, var):
    global _F, _D
    _F = h5py.File(path, "r", rdcc_nbytes=0)
    _D = _F[var]


def _work(task):
    """One spatial tile, all 365 days, read in chunk-aligned blocks of tb days."""
    (ty, tx, tb, y0, y1, x0, x1, rows_u, cols_u, wdata, windices, windptr, lo, hi) = task
    tic = time.process_time()
    W = csr_matrix((wdata, windices, windptr), shape=(len(rows_u), len(cols_u)))
    num = np.empty((len(rows_u), 365))
    den = np.empty((len(rows_u), 365))
    n_oor = n_valid = nbytes = 0
    vmin, vmax = np.inf, -np.inf
    for t0 in range(0, 365, tb):
        t1 = min(t0 + tb, 365)
        block = _D[t0:t1, y0:y1, x0:x1]                  # (days, h, w) float32, whole chunks
        nbytes += block.nbytes
        flat = block.reshape(t1 - t0, -1)
        for s0 in range(t0, t1, COMPUTE_DAYS):
            s1 = min(s0 + COMPUTE_DAYS, t1)
            V = np.ascontiguousarray(flat[s0 - t0:s1 - t0][:, cols_u].T)   # (cells, days) float32
            valid = (V != FILL) & np.isfinite(V)
            if valid.any():
                vv = V[valid]
                n_valid += vv.size
                n_oor += int(((vv < lo) | (vv > hi)).sum())
                vmin, vmax = min(vmin, float(vv.min())), max(vmax, float(vv.max()))
            num[:, s0:s1] = W @ np.where(valid, V, 0.0).astype(np.float64)
            v_all, v_any = valid.all(axis=1), valid.any(axis=1)
            if np.array_equal(v_all, v_any):                    # validity constant in time here
                den[:, s0:s1] = (W @ v_all.astype(np.float64))[:, None]
            else:
                den[:, s0:s1] = W @ valid.astype(np.float64)
        del block, flat
    return (ty, tx, 0, 365, rows_u, num, den, nbytes, n_valid, n_oor, vmin, vmax,
            time.process_time() - tic)


def plan_tasks(basin_idx, cell_id, cov, chunks, lo, hi):
    """Group the weights by the file's spatial chunk tiles (int32 keys, slices not copies)."""
    ct, cy, cx = chunks
    ntx = -(-NX // cx)
    cell_id = cell_id.astype(np.int32, copy=False)     # < 63.1 M cells, so int32 is exact
    rows = cell_id // NX
    cols = cell_id - rows * NX
    key = (rows // cy) * ntx + (cols // cx)
    order = np.argsort(key, kind="stable")
    key, rows, cols = key[order], rows[order], cols[order]
    b, w = basin_idx[order], cov[order]
    del order
    bounds = np.concatenate([[0], np.flatnonzero(np.diff(key)) + 1, [len(key)]])
    tiles = []
    for s0, s1 in zip(bounds[:-1], bounds[1:]):
        k = int(key[s0])
        ty, tx = divmod(k, ntx)
        y0, x0 = ty * cy, tx * cx
        y1, x1 = min(y0 + cy, NY), min(x0 + cx, NX)
        local = (rows[s0:s1] - y0) * (x1 - x0) + (cols[s0:s1] - x0)
        rows_u, inv_r = np.unique(b[s0:s1], return_inverse=True)
        cols_u, inv_c = np.unique(local, return_inverse=True)
        W = csr_matrix((w[s0:s1], (inv_r, inv_c)), shape=(len(rows_u), len(cols_u)))
        W.sum_duplicates()
        tiles.append((ty, tx, y0, y1, x0, x1, rows_u, cols_u, W))
    day_bytes = cy * cx * 4
    tb = max(ct, (READ_TARGET_BYTES // (day_bytes * ct)) * ct)
    # one task per tile (its weights are pickled once); the worker loops over the time blocks
    tasks = [(ty, tx, tb, y0, y1, x0, x1, rows_u, cols_u.astype(np.int64), W.data, W.indices, W.indptr, lo, hi)
             for (ty, tx, y0, y1, x0, x1, rows_u, cols_u, W) in tiles]
    tasks.sort(key=lambda t: -len(t[9]))        # biggest first, for load balance (stable)
    n_tiles_total = -(-NY // cy) * -(-NX // cx)
    return tasks, len(tiles), n_tiles_total, tb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True)
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--weight", choices=["area", "coverage"], default="area")
    ap.add_argument("--source-sha256", default=None,
                    help="the input file's SHA-256 as verified by the caller; recorded in the outputs")
    ap.add_argument("--allow-out-of-range", action="store_true",
                    help="write the .done even when a cell value falls outside RANGES")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    wall0, cpu0 = time.time(), time.process_time()
    grid = read_grid(a.file)
    var, year = grid.var, int(str(grid.dates[0])[:4])
    if str(grid.dates[-1])[:4] != str(year):
        sys.exit(f"{a.file}: time axis spans two years")
    stem = os.path.join(a.out_dir, f"{var}_{year}")
    if os.path.exists(stem + ".done"):
        print(f"skip {var} {year}: {stem}.done exists")
        return
    lo, hi = RANGES[var]

    wcol, scol = {"area": ("area_km2", "area_sum_km2"), "coverage": ("coverage", "coverage_sum")}[a.weight]
    t = pq.read_table(os.path.join(a.weights_dir, "daymet_weights.parquet"), columns=["basin_idx", "cell_id", wcol])
    w = {c: t.column(c).to_numpy() for c in ("basin_idx", "cell_id", wcol)}
    del t
    basins = pd.read_csv(os.path.join(a.weights_dir, "daymet_basins.csv"), dtype={"site_id": str})
    nb = len(basins)
    wmeta = os.path.join(a.weights_dir, "daymet_weights_meta.json")
    weights_md5 = json.load(open(wmeta)).get("weights_md5") if os.path.exists(wmeta) else None
    tasks, n_tiles, n_tiles_total, tb = plan_tasks(
        w["basin_idx"], w["cell_id"], w[wcol], grid.chunks, lo, hi)
    del w
    print(f"{var} {year}: chunks {grid.chunks}, {n_tiles}/{n_tiles_total} tiles hold weighted cells "
          f"({100 * n_tiles / n_tiles_total:.1f}%), {n_tiles * -(-365 // tb)} reads of {min(tb, 365)} days, "
          f"{a.workers} workers, "
          f"{a.weight} weights",
          flush=True)

    num = np.zeros((nb, 365))
    den = np.zeros((nb, 365))
    nbytes = n_valid = n_oor = 0
    vmin, vmax, worker_cpu = np.inf, -np.inf, 0.0
    t_pool = time.time()
    with get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.file, var)) as pool:
        # ordered imap: partial sums are added in the fixed task order, so the output is
        # bit-identical across reruns and across --workers values
        for i, r in enumerate(pool.imap(_work, tasks, chunksize=1), 1):
            ty, tx, t0, t1, rows_u, nm, dn, nby, nv, no, mn, mx, cpu = r
            num[rows_u, t0:t1] += nm
            den[rows_u, t0:t1] += dn
            nbytes += nby; n_valid += nv; n_oor += no; worker_cpu += cpu
            vmin, vmax = min(vmin, mn), max(vmax, mx)
            if i % max(1, len(tasks) // 10) == 0:
                print(f"  {i}/{len(tasks)} tiles, {time.time() - t_pool:.0f} s", flush=True)
    t_read = time.time() - t_pool

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(den > 0, num / den, np.nan)
    vf = den / basins[scol].to_numpy()[:, None]
    qa = pd.DataFrame({
        "site_id": basins["site_id"], "valid_frac_min": vf.min(axis=1), "valid_frac_max": vf.max(axis=1),
        "n_nan_days": np.isnan(mean).sum(axis=1), "min": np.nanmin(np.where(np.isnan(mean), np.inf, mean), axis=1),
        "max": np.nanmax(np.where(np.isnan(mean), -np.inf, mean), axis=1),
        "mean": np.nanmean(mean, axis=1) if np.isfinite(mean).any() else np.nan,
    })
    if (vf > 1 + 1e-9).any():
        sys.exit(f"{var} {year}: valid weight exceeds the basin's total weight (double counting)")

    dates = grid.dates.astype("datetime64[D]")
    table = pa.table({
        "site_id": pa.array(np.repeat(basins["site_id"].to_numpy(), 365), pa.string()),
        "Date": pa.array(np.tile(dates, nb), pa.date32()),
        var: pa.array(mean.reshape(-1), pa.float64()),
    })
    pq.write_table(table, stem + ".parquet.tmp", compression="zstd")
    os.replace(stem + ".parquet.tmp", stem + ".parquet")
    qa.to_csv(stem + "_qa.csv", index=False)

    ru_self = resource.getrusage(resource.RUSAGE_SELF)
    ru_kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    timing = {
        "file": os.path.basename(a.file), "file_bytes": os.path.getsize(a.file),
        "file_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(a.file))),
        "var": var, "year": year, "chunks": list(grid.chunks), "shuffle": grid.shuffle,
        "hdf5_fillvalue": grid.hdf5_fillvalue, "workers": a.workers, "weight": a.weight,
        "weights_dir": os.path.abspath(a.weights_dir),
        "tiles_read": n_tiles, "tiles_total": n_tiles_total, "reads": n_tiles * -(-365 // tb),
        "days_per_read": min(tb, 365),
        "bytes_decompressed": int(nbytes), "wall_seconds_total": round(time.time() - wall0, 2),
        "wall_seconds_read_aggregate": round(t_read, 2), "worker_cpu_seconds": round(worker_cpu, 1),
        "parent_cpu_seconds": round(time.process_time() - cpu0, 1),
        "peak_rss_parent_mb": round(ru_self.ru_maxrss / 2**20, 1),      # macOS: bytes
        "peak_rss_largest_worker_mb": round(ru_kids.ru_maxrss / 2**20, 1),
        "valid_cell_days": int(n_valid), "out_of_range_cell_days": int(n_oor),
        "valid_value_min": vmin, "valid_value_max": vmax, "range_checked": [lo, hi],
        "basins": nb, "basins_with_nan_days": int((qa["n_nan_days"] > 0).sum()),
        "output_bytes": os.path.getsize(stem + ".parquet"),
        "source_sha256": a.source_sha256, "weights_md5": weights_md5, "git": git_state(),
        "finished_utc": utc_now(),
    }
    write_json(stem + "_timing.json", timing)
    if n_oor and not a.allow_out_of_range:
        sys.exit(f"{var} {year}: {n_oor} cell-days outside {RANGES[var]} (values {vmin:.6g} .. {vmax:.6g}); "
                 f"outputs written for inspection, no .done (rerun with --allow-out-of-range to accept)")
    open(stem + ".done", "w").write(json.dumps({
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"), "finished_utc": timing["finished_utc"],
        "source_sha256": a.source_sha256, "weights_md5": weights_md5, "commit": timing["git"]["commit"]}))
    print(f"{var} {year}: {timing['wall_seconds_total']} s wall ({t_read:.0f} s read+aggregate, "
          f"{nbytes / 1e9:.1f} GB decompressed, worker CPU {worker_cpu:.0f} s), peak RSS parent "
          f"{timing['peak_rss_parent_mb']} MB / worker {timing['peak_rss_largest_worker_mb']} MB; "
          f"valid range [{vmin:.3g}, {vmax:.3g}], out-of-range cell-days {n_oor}; "
          f"basins with NaN days {timing['basins_with_nan_days']}")


if __name__ == "__main__":
    main()
