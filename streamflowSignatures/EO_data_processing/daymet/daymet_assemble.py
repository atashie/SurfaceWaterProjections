"""Assemble per-variable-year basin series into the pipeline's climate parquet.

Usage:
    python daymet_assemble.py --in-dir <dir with <var>_<year>.parquet> --years 2023 (or 1980-2025)
        --out <daymet_<years>_<date>.parquet> [--vars prcp,tmin,tmax,swe,vp,srad]

Output schema = the stale input's (convert_daymet_csvs_to_parquet.py):
    site_id string, Date date32, prcp tmin tmax swe vp srad float64
one row per site-day, 365 rows per site-year (Dec 31 absent in leap years, as Daymet
publishes it), sorted by site_id then Date. Written year by year, so memory stays at
one year of the six variables (~2.9 M rows).

Structural checks (fatal): every variable-year file present; identical (site_id, Date)
keys across variables; 365 rows per site-year; no duplicate keys; dates inside the year;
Dec 31 absent exactly in leap years; site set identical across years. NaN counts per
variable are reported (a NaN day = no valid Daymet cell under the basin that day).
A <out>.provenance.json sidecar records inputs, their md5s, row counts and NaN counts.
"""
import argparse
import calendar
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import DAYMET_VARS, md5sum, write_json  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--years", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    a = ap.parse_args()
    y = [int(v) for v in a.years.split("-")]
    years = list(range(y[0], y[-1] + 1))
    variables = [v for v in a.vars.split(",") if v]
    schema = pa.schema([("site_id", pa.string()), ("Date", pa.date32())] + [(v, pa.float64()) for v in variables])

    prov = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "years": years, "vars": variables,
            "inputs": {}, "rows": 0, "nan_days": {v: 0 for v in variables}}
    sites_ref = None
    tmp = a.out + ".tmp"
    writer = pq.ParquetWriter(tmp, schema, compression="zstd")
    try:
        for year in years:
            frame = None
            for v in variables:
                p = os.path.join(a.in_dir, f"{v}_{year}.parquet")
                if not os.path.exists(p):
                    sys.exit(f"missing {p}")
                prov["inputs"][os.path.basename(p)] = md5sum(p)
                t = pd.read_parquet(p)
                if list(t.columns) != ["site_id", "Date", v]:
                    sys.exit(f"{p}: unexpected columns {list(t.columns)}")
                if frame is None:
                    frame = t
                else:
                    if not (frame["site_id"].equals(t["site_id"]) and frame["Date"].equals(t["Date"])):
                        sys.exit(f"{p}: (site_id, Date) keys differ from {variables[0]}_{year}")
                    frame[v] = t[v].to_numpy()
            frame["Date"] = pd.to_datetime(frame["Date"])
            per_site = frame.groupby("site_id", sort=False).size()
            if (per_site != 365).any():
                sys.exit(f"{year}: {int((per_site != 365).sum())} sites without 365 rows")
            if frame.duplicated(["site_id", "Date"]).any():
                sys.exit(f"{year}: duplicate (site_id, Date)")
            if (frame["Date"].dt.year != year).any():
                sys.exit(f"{year}: dates outside the year")
            has_dec31 = ((frame["Date"].dt.month == 12) & (frame["Date"].dt.day == 31)).any()
            if has_dec31 == calendar.isleap(year):
                sys.exit(f"{year}: Dec 31 {'present' if has_dec31 else 'absent'} in a "
                         f"{'leap' if calendar.isleap(year) else 'common'} year")
            sites = set(per_site.index)
            if sites_ref is None:
                sites_ref = sites
            elif sites != sites_ref:
                sys.exit(f"{year}: site set differs from {years[0]}")
            frame = frame.sort_values(["site_id", "Date"], kind="mergesort")
            for v in variables:
                prov["nan_days"][v] += int(frame[v].isna().sum())
            writer.write_table(pa.Table.from_pandas(frame[["site_id", "Date"] + variables], schema=schema,
                                                    preserve_index=False))
            prov["rows"] += len(frame)
            print(f"  {year}: {len(frame):,} rows, {len(sites):,} sites", flush=True)
    finally:
        writer.close()
    os.replace(tmp, a.out)
    with open(a.out, "rb") as fh:
        fh.seek(-4, 2)
        if fh.read(4) != b"PAR1":
            sys.exit(f"{a.out}: no PAR1 footer")
    prov.update({"sites": len(sites_ref), "out": os.path.basename(a.out), "out_bytes": os.path.getsize(a.out),
                 "out_md5": md5sum(a.out)})
    write_json(a.out + ".provenance.json", prov)
    print(f"wrote {a.out}: {prov['rows']:,} rows, {prov['sites']:,} sites, {prov['out_bytes'] / 1e6:.1f} MB, "
          f"NaN days {prov['nan_days']}")


if __name__ == "__main__":
    main()
