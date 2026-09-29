"""Probe one Daymet NA file: layout facts, and the unstored-chunk hazard (plan unknown U3).

Prints shape, chunking, compression, fill values, CRS, time axis and calendar, then
counts the HDF5 chunks actually stored against the chunks the shape implies. A chunk
that was never written reads back as the HDF5 dataset fill value -- 0.0 in the
1980-2019 files, NOT the -9999 `_FillValue` -- so a basin whose cells fall in such a
chunk would silently average valid-looking zeros. With --weights-dir, the probe lists
every unstored chunk that holds a weighted cell (must be none) and reports what
fraction of the grid's tiles the basins touch (plan unknown U4).

Usage:
    python daymet_probe.py --file <nc> [--weights-dir <dir>] [--out <json>]
"""
import argparse
import calendar
import os
import sys
import time

import h5py
import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import FILL, NX, NY, read_grid, write_json  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True)
    ap.add_argument("--weights-dir", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t0 = time.time()
    g = read_grid(a.file)
    year = int(str(g.dates[0])[:4])
    ct, cy, cx = g.chunks
    nct, ncy, ncx = -(-365 // ct), -(-NY // cy), -(-NX // cx)
    res = {"file": os.path.basename(a.file), "bytes": os.path.getsize(a.file), "var": g.var, "year": year,
           "shape": list(g.shape), "chunks": list(g.chunks), "compression": g.compression, "shuffle": g.shuffle,
           "hdf5_fillvalue": g.hdf5_fillvalue, "cf_fillvalue": FILL, "crs": g.crs.to_proj4(),
           "first_date": str(g.dates[0]), "last_date": str(g.dates[-1]), "leap_year": calendar.isleap(year),
           "feb29_present": bool((g.dates == np.datetime64(f"{year}-02-29")).any()) if calendar.isleap(year) else None,
           "dec31_present": bool((g.dates == np.datetime64(f"{year}-12-31")).any()),
           "dates_strictly_daily": bool((np.diff(g.dates.astype("int64")) == 1).all()),
           "chunks_possible": nct * ncy * ncx}
    with h5py.File(a.file, "r") as f:
        d = f[g.var]
        res["global_attrs"] = {k: (v.decode() if isinstance(v, bytes) else str(v)) for k, v in f.attrs.items()}
        n = d.id.get_num_chunks()
        res["chunks_stored"] = int(n)
        stored = set()
        if n != res["chunks_possible"]:
            for i in range(n):
                off = d.id.get_chunk_info(i).chunk_offset
                stored.add((off[0] // ct, off[1] // cy, off[2] // cx))
    missing = []
    if n != res["chunks_possible"]:
        for it in range(nct):
            for iy in range(ncy):
                for ix in range(ncx):
                    if (it, iy, ix) not in stored:
                        missing.append((it, iy, ix))
    res["chunks_missing"] = len(missing)
    res["missing_spatial_tiles"] = sorted({(iy, ix) for _, iy, ix in missing})
    if a.weights_dir:
        cid = pq.read_table(os.path.join(a.weights_dir, "daymet_weights.parquet"), columns=["cell_id"])["cell_id"].to_numpy()
        tiles = np.unique((cid // NX // cy) * ncx + (cid % NX) // cx)
        touched = {(int(k) // ncx, int(k) % ncx) for k in tiles}
        res["tiles_touched_by_basins"] = len(touched)
        res["tiles_total"] = ncy * ncx
        res["missing_chunks_in_touched_tiles"] = sorted(
            {(iy, ix) for (iy, ix) in res["missing_spatial_tiles"] if (iy, ix) in touched})
    res["seconds"] = round(time.time() - t0, 1)
    for k, v in res.items():
        if k != "global_attrs":
            print(f"{k}: {v}")
    if a.out:
        write_json(a.out, res)
    if res.get("missing_chunks_in_touched_tiles"):
        print("HAZARD: unstored chunks under basin weights -- they read as the HDF5 fill value, not -9999")
        sys.exit(2)


if __name__ == "__main__":
    main()
