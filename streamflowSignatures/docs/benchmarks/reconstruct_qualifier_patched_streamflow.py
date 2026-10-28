#!/usr/bin/env python
"""Reconstruction recipe for the corrected HISSS daily streamflow table (deposited with resource 3).

    python hisss_streamflow_reconstruction_recipe.py --base hisss_streamflow_daily_as_compiled_feb2026.parquet \
        --corrections hisss_streamflow_qualifier_corrections.csv --expected hisss_streamflow_daily.parquet

Applies the correction ledger (one row per restored USGS day: `gage_id`, `date`, `Q_new`) to the February 2026
compile and compares the result with the deposited corrected table row group by row group (both files keep the
compile's 909 row groups, row order and schema; only the `Q` cells named in the ledger differ). Needs pyarrow and
pandas and about 2 GB of memory; prints the number of rows compared and exits 0 only when every row group is equal.
"""
import argparse, sys
import pandas as pd, pyarrow as pa, pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--base", required=True); ap.add_argument("--corrections", required=True); ap.add_argument("--expected", required=True)
a = ap.parse_args()
led = pd.read_csv(a.corrections, dtype={"gage_id": str, "date": str}, usecols=["gage_id", "date", "Q_new", "status"], float_precision="round_trip")  # the default parser is not round-trip
led = led[led["status"] == "restored"]
led["date"] = pd.to_datetime(led["date"]).dt.date
key = dict(zip(zip(led["gage_id"], led["date"]), led["Q_new"].astype(float)))
base, exp = pq.ParquetFile(a.base), pq.ParquetFile(a.expected)
if base.metadata.num_row_groups != exp.metadata.num_row_groups:
    sys.exit(f"row-group layouts differ ({base.metadata.num_row_groups} vs {exp.metadata.num_row_groups}); this recipe assumes the compile's layout")
n_rows = n_applied = 0; bad = []
corr = pd.Series(key); corr.index = pd.MultiIndex.from_tuples(corr.index)
for i in range(base.metadata.num_row_groups):
    b = base.read_row_group(i).to_pandas(); e = exp.read_row_group(i).to_pandas()
    if list(b.columns) != list(e.columns) or len(b) != len(e):
        bad.append((i, "shape")); continue
    idx = pd.MultiIndex.from_arrays([b["gage_id"].astype(str).to_numpy(), pd.to_datetime(b["Date"]).dt.date.to_numpy()])
    v = corr.reindex(idx).to_numpy(); hit = ~pd.isna(v)
    q = b["Q"].astype(float).to_numpy(copy=True)
    if hit.any() and not pd.isna(q[hit]).all():
        bad.append((i, "base Q is not NA on a ledger day")); continue
    q[hit] = v[hit]; n_applied += int(hit.sum()); b["Q"] = q
    other = [c for c in b.columns if c != "Q"]
    same = b[other].equals(e[other]) and bool(((b["Q"].isna() & e["Q"].isna()) | (b["Q"] == e["Q"])).all())
    if not same:
        bad.append((i, "values"))
    n_rows += len(b)
print(f"rows compared {n_rows:,}; ledger days applied {n_applied:,} of {len(key):,}; row groups differing {len(bad)} {bad[:5]}")
sys.exit(0 if not bad and n_applied == len(key) else 1)
