#!/usr/bin/env python3
"""climate_frame_inventory.py — the inserted-row inventory of the climate-frame fix
(docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.5 item 1b, B′; 2026-10-26).

For every gage of --gages (one id per line; default: every gage of the streamflow input),
builds the production frame with docs/benchmarks/climate_frame.py in `left` and `outer`
mode from the streamflow and climate parquets (the batch inputs of run_batched_julia.py
split, or the full files) and writes every row the outer frame ADDS: gage_id, date,
water_year, PPT, SWE, TMIN, TMAX, and whether the date is a Daymet-missing Dec 31 of a leap
year (absent from the climate input itself — those are never added; they are listed from
the streamflow side as `daymet_missing_dec31` for the ledger's secondary-change bookkeeping).

    python docs/benchmarks/climate_frame_inventory.py --streamflow SF.parquet --climate C.parquet
        [--gages LIST] [--window 1993-2025] --out inserted_rows.csv

Output columns: gage_id, date, water_year, PPT, SWE, TMIN, TMAX, kind (inserted |
daymet_missing_dec31_neighbour), in_window. The streamflow parquet is streamed row group by
row group (rows contiguous per gage, any gage order); the climate parquet must be sorted by
site id (the batch files of run_batched_julia.py split are) and is read per gage through a
row-group index; memory stays at a few row groups.
"""
import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))
from climate_frame import frame_with_climate  # noqa: E402
from streamflow_signatures import add_water_year_columns  # noqa: E402

RENAME = {"Date": "date", "site_id": "gage_id", "prcp": "PPT", "swe": "SWE", "tmin": "TMIN", "tmax": "TMAX"}


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def iter_gages(parquet, id_col, columns):
    """Stream a parquet whose rows are contiguous per gage: yield (gage_id, frame)."""
    f = pq.ParquetFile(parquet)
    cur, parts = None, []
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=columns).to_pandas()
        t = t.rename(columns={k: v for k, v in RENAME.items() if k in t.columns})
        gid = t["gage_id"].astype(str).to_numpy(dtype=object)
        change = np.flatnonzero(gid[1:] != gid[:-1]) + 1
        for s, e in zip(np.concatenate(([0], change)), np.concatenate((change, [len(gid)]))):
            g = gid[s]
            if g != cur:
                if cur is not None:
                    yield cur, pd.concat(parts, ignore_index=True)
                cur, parts = g, []
            parts.append(t.iloc[s:e])
    if cur is not None:
        yield cur, pd.concat(parts, ignore_index=True)


class ClimateIndex:
    """Per-gage access to a climate parquet sorted by site id (the batch files of
    run_batched_julia.py split are): the row groups' id ranges are scanned once, then a gage's
    rows come from the one or two row groups covering it, with the last groups cached."""

    def __init__(self, parquet, id_col, columns):
        self.f = pq.ParquetFile(parquet); self.id_col = id_col; self.columns = columns
        self.ranges = []
        for i in range(self.f.num_row_groups):
            ids = self.f.read_row_group(i, columns=[id_col]).column(0)
            self.ranges.append((str(ids[0]), str(ids[-1])))
        self.cache = {}

    def get(self, gage):
        parts = []
        for i, (lo, hi) in enumerate(self.ranges):
            if lo <= gage <= hi:
                if i not in self.cache:
                    if len(self.cache) > 3:
                        self.cache.pop(next(iter(self.cache)))
                    t = self.f.read_row_group(i, columns=self.columns).to_pandas()
                    t = t.rename(columns={k: v for k, v in RENAME.items() if k in t.columns})
                    t["gage_id"] = t["gage_id"].astype(str)
                    self.cache[i] = t
                t = self.cache[i]
                parts.append(t[t["gage_id"] == gage])
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["gage_id", "date"])


def iter_pairs(streamflow, sf_id, climate, cl_id, wanted):
    """Stream the streamflow parquet gage by gage (any gage order) and look each gage up in the
    climate parquet; yield (gage, sf, cl) for the wanted gages that have climate rows."""
    cl_cols = [cl_id, "Date" if "Date" in pq.read_schema(climate).names else "date",
               *[c for c in ("prcp", "PPT", "swe", "SWE", "tmin", "TMIN", "tmax", "TMAX") if c in pq.read_schema(climate).names]]
    idx = ClimateIndex(climate, cl_id, cl_cols)
    for g, sf in iter_gages(streamflow, sf_id, [sf_id, "Date" if "Date" in pq.read_schema(streamflow).names else "date", "Q", "water_year"]):
        if wanted is not None and g not in wanted:
            continue
        cl = idx.get(g)
        if len(cl):
            yield g, sf, cl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--streamflow", required=True); ap.add_argument("--climate", required=True)
    ap.add_argument("--gages", default=None); ap.add_argument("--window", default=None, help="W0-W1: count only water years inside")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sf_schema = pq.read_schema(a.streamflow).names; cl_schema = pq.read_schema(a.climate).names
    sf_id = "gage_id" if "gage_id" in sf_schema else "site_id"
    cl_id = "site_id" if "site_id" in cl_schema else "gage_id"
    wanted = set(l.strip() for l in open(a.gages) if l.strip()) if a.gages else None
    w = tuple(int(x) for x in a.window.split("-")) if a.window else None
    rows = []
    t0 = time.time()
    n_gages = 0
    for i, (g, sf, cl) in enumerate(iter_pairs(a.streamflow, sf_id, a.climate, cl_id, wanted), 1):
        n_gages += 1
        sf = sf.copy(); cl = cl.copy()
        sf["gage_id"] = sf["gage_id"].astype(str); sf["date"] = pd.to_datetime(sf["date"])
        if "water_year" not in sf.columns or sf["water_year"].isna().any():
            sf = add_water_year_columns(sf)
        cl["gage_id"] = cl["gage_id"].astype(str); cl["date"] = pd.to_datetime(cl["date"])
        keep = [c for c in ("gage_id", "date", "PPT", "SWE", "TMIN", "TMAX") if c in cl.columns]
        cl = cl[keep]
        L = frame_with_climate(sf, cl, "left"); O = frame_with_climate(sf, cl, "outer")
        add = O.iloc[len(L):]
        for _, r in add.iterrows():
            rows.append({"gage_id": g, "date": r["date"].date().isoformat(), "water_year": int(r["water_year"]),
                         **{c: r.get(c, np.nan) for c in ("PPT", "SWE", "TMIN", "TMAX")}, "kind": "inserted"})
        # Daymet-missing Dec 31 (leap years) inside the record span whose Dec 30 or Jan 1 was an absent streamflow date
        added_dates = set(pd.to_datetime(add["date"]).dt.date) if len(add) else set()
        lo, hi = sf["date"].min(), sf["date"].max()
        for y in range(lo.year, hi.year + 1):
            if (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)):
                d31 = dt.date(y, 12, 31)
                if lo.date() <= d31 <= hi.date() and (dt.date(y, 12, 30) in added_dates or dt.date(y + 1, 1, 1) in added_dates):
                    rows.append({"gage_id": g, "date": d31.isoformat(), "water_year": y + 1, "PPT": np.nan, "SWE": np.nan,
                                 "TMIN": np.nan, "TMAX": np.nan, "kind": "daymet_missing_dec31_neighbour"})
        if i % 500 == 0:
            log(f"  {i} gages, {len(rows):,} rows ({time.time() - t0:.0f}s)")
    inv = pd.DataFrame(rows, columns=["gage_id", "date", "water_year", "PPT", "SWE", "TMIN", "TMAX", "kind"])
    if w:
        inv["in_window"] = (inv["water_year"] >= w[0]) & (inv["water_year"] <= w[1])
    inv.to_csv(a.out, index=False)
    ins = inv[inv["kind"] == "inserted"]
    log(f"{n_gages} gages: {len(ins):,} inserted rows on {ins['gage_id'].nunique()} gages"
        + (f" ({int(ins['in_window'].sum()):,} inside {a.window} on {ins.loc[ins['in_window'], 'gage_id'].nunique()} gages)" if w else "")
        + f"; Dec 31 neighbours {int((inv['kind'] != 'inserted').sum())} -> {a.out}")


if __name__ == "__main__":
    main()
