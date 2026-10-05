"""Assemble per-variable-year basin series into the pipeline's climate parquet.

Usage:
    python daymet_assemble.py --in-dir <dir with <var>_<year>.parquet> --years 2023 (or 1980-2025)
        --out <daymet_<years>_<date>.parquet> [--vars prcp,tmin,tmax,swe,vp,srad]
        [--provenance-only] [--no-verify]

Output schema = the stale input's (convert_daymet_csvs_to_parquet.py):
    site_id string, Date date32, prcp tmin tmax swe vp srad float64
One row per site-day, 365 rows per site-year (Dec 31 absent in leap years, as Daymet
publishes it). Written year by year, so memory stays at one year of the six variables
(~2.9 M rows). The file is therefore ordered (year, site_id, Date): it is NOT contiguous
per site across years. Readers group by site_id; the Julia runner does.

Structural checks (fatal):
  * every variable-year file present;
  * identical (site_id, Date) keys across variables;
  * 365 rows per site-year, no duplicate keys, dates inside the year;
  * Dec 31 absent exactly in leap years;
  * site set identical across years.
NaN counts per variable are reported (a NaN day = no valid Daymet cell under the basin
that day).

Inputs: every variable-year needs its <var>_<year>.done (an unfinished or out-of-range
aggregation has none), and all must come from one weights file (the md5 in the .done or
the timing JSON, where recorded).

Verification (fatal; skipped with --no-verify): the output is written to <out>.tmp, re-read
and compared value by value with the per-variable-year inputs. That covers keys, row order
and every float bit, NaN and -0.0 included (NaN is stored as NaN, never as NULL). Only a
verified file replaces <out>; on failure the .tmp is deleted and an earlier <out> is left
as it was. daymet_validate.py compares the per-year files with the stale product, so this
check is what carries that validation over to the assembled file.

Provenance sidecar <out>.provenance.json:
  * the inputs and their md5s, row and NaN counts, and the output size and md5;
  * the verification result;
  * the git state and software of this step;
  * from the run folder: each run_meta_*.json (commit, dirty flag), the CMR record (size,
    SHA-256, concept id, revision) of every source granule, the stream_log.csv tallies,
    and the weights and polygons md5 from the weights folder named in the timing files.
--provenance-only rewrites that sidecar for an existing output without touching the
parquet. It requires the output's md5, and the input md5s, to equal an existing
sidecar's, and runs the verification. It is used to upgrade a run assembled before
2026-10-01.
"""
import argparse
import calendar
import csv
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import DAYMET_VARS, git_state, md5sum, software_versions, utc_now, write_json  # noqa: E402


def load_year(in_dir, year, variables, md5s=None):
    """The year's six variables joined on (site_id, Date), checked and sorted as written."""
    frame = None
    for v in variables:
        p = os.path.join(in_dir, f"{v}_{year}.parquet")
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        if md5s is not None:
            md5s[os.path.basename(p)] = md5sum(p)
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
    return frame.sort_values(["site_id", "Date"], kind="mergesort"), set(per_site.index)


def verify_output(out, in_dir, years, variables):
    """Re-read `out` and compare it, bit for bit and in row order, with the per-year inputs."""
    batches = pq.ParquetFile(out).iter_batches(batch_size=1 << 20)
    rest, checked = None, 0
    for year in years:
        exp, _ = load_year(in_dir, year, variables)
        need, parts, have = len(exp), ([rest] if rest is not None else []), (len(rest) if rest is not None else 0)
        while have < need:
            b = next(batches, None)
            if b is None:
                sys.exit(f"VERIFY: {out} ends inside {year} ({checked + have:,} rows)")
            parts.append(pa.Table.from_batches([b]))
            have += b.num_rows
        tab = pa.concat_tables(parts)
        got, rest = tab.slice(0, need), (tab.slice(need) if have > need else None)
        if not pc.all(pc.equal(got["site_id"], pa.array(exp["site_id"].to_numpy(), pa.string()))).as_py():
            sys.exit(f"VERIFY: {year}: site_id differs from the inputs")
        if not pc.all(pc.equal(got["Date"], pa.array(exp["Date"].to_numpy().astype("datetime64[D]"),
                                                     pa.date32()))).as_py():
            sys.exit(f"VERIFY: {year}: Date differs from the inputs")
        for v in variables:
            if got[v].null_count:
                sys.exit(f"VERIFY: {year} {v}: {got[v].null_count} NULLs; the inputs hold values or NaN there")
            a = got[v].to_numpy().view(np.int64)
            b = np.ascontiguousarray(exp[v].to_numpy(dtype=np.float64)).view(np.int64)
            if not np.array_equal(a, b):
                sys.exit(f"VERIFY: {year} {v}: {int((a != b).sum())} values differ from the inputs")
        checked += need
    if (rest is not None and len(rest)) or next(batches, None) is not None:
        sys.exit(f"VERIFY: {out} has rows beyond the last year")
    return {"rows_compared": checked, "values_compared": checked * len(variables), "mismatches": 0,
            "method": "re-read output vs per-variable-year inputs: keys, row order, float64 bits"}


def run_provenance(in_dir, years, variables):
    """What the run folder records about how the inputs were made (no recomputation)."""
    timing = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(in_dir, "*_timing.json")))]
    weights = {}
    for wd in sorted({t.get("weights_dir") for t in timing if t.get("weights_dir")}):
        mp = os.path.join(wd, "daymet_weights_meta.json")
        m = json.load(open(mp)) if os.path.exists(mp) else {}
        weights[wd] = {k: m.get(k) for k in ("polygons", "polygons_md5", "weights_md5", "n_polygons",
                                             "grid_file", "crs_wkt")}
    run_meta = []
    for p in sorted(glob.glob(os.path.join(in_dir, "run_meta_*.json"))):
        m = json.load(open(p))
        run_meta.append({"file": os.path.basename(p), "started": m.get("started"), "git": m.get("git"),
                         "argv": m.get("argv"), "planned": m.get("planned")})
    source, mp = {}, os.path.join(in_dir, "cmr_manifest.json")
    if os.path.exists(mp):
        man = json.load(open(mp))
        for v in variables:
            for y in years:
                n = f"daymet_v4_daily_na_{v}_{y}.nc"
                if n in man:
                    source[n] = {k: man[n].get(k) for k in ("bytes", "sha256", "algorithm", "concept_id",
                                                            "revision_id", "cmr_revision") if k in man[n]}
    for t in timing:
        if t.get("source_sha256") and t.get("file") in source:
            source[t["file"]]["sha256_recorded_by_aggregator"] = t["source_sha256"]
    log, lp = None, os.path.join(in_dir, "stream_log.csv")
    if os.path.exists(lp):
        rows = list(csv.DictReader(open(lp)))
        log = {"rows": len(rows), "sha256_ok": sum(r["sha256_ok"] == "True" for r in rows),
               "probe_ok": sum(r["probe_ok"] == "True" for r in rows),
               "sources": {s: sum(r["source"] == s for r in rows) for s in sorted({r["source"] for r in rows})}}
    return {"run_meta": run_meta, "weights": weights,
            "cmr_manifest_md5": md5sum(mp) if os.path.exists(mp) else None,
            "source_granules": source, "stream_log": log,
            "timing_files": len(timing),
            "aggregate_commits": sorted({(t.get("git") or {}).get("commit") for t in timing} - {None})}


def check_output(path, a, years, variables):
    """PAR1 at both ends, then (unless --no-verify) the bit-for-bit re-read against the inputs."""
    with open(path, "rb") as fh:
        head = fh.read(4)
        fh.seek(-4, 2)
        if head != b"PAR1" or fh.read(4) != b"PAR1":
            sys.exit(f"{path}: no PAR1 magic at both ends")
    if a.no_verify:
        return None
    t0 = time.time()
    ver = dict(verify_output(path, a.in_dir, years, variables), seconds=round(time.time() - t0, 1), at_utc=utc_now())
    print(f"verified: {ver['values_compared']:,} values equal the inputs ({ver['seconds']} s)", flush=True)
    return ver


def done_weights(in_dir, years, variables):
    """Every variable-year must have its .done (an out-of-range or unfinished aggregation has none),
    and all must come from one weights file: the md5 recorded in the .done or the timing JSON."""
    missing = [f"{v}_{y}" for y in years for v in variables if not os.path.exists(os.path.join(in_dir, f"{v}_{y}.done"))]
    if missing:
        sys.exit(f"{len(missing)} variable-years have no .done (unfinished, or out of range), e.g. {missing[:5]}")
    seen = set()
    for y in years:
        for v in variables:
            for f in (f"{v}_{y}.done", f"{v}_{y}_timing.json"):
                try:
                    w = json.load(open(os.path.join(in_dir, f))).get("weights_md5")
                except (OSError, ValueError):
                    w = None
                if w:
                    seen.add(w)
                    break
    if len(seen) > 1:
        sys.exit(f"the inputs were aggregated with {len(seen)} different weights files: {sorted(seen)}")
    return seen.pop() if seen else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--years", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    ap.add_argument("--provenance-only", action="store_true",
                    help="verify an existing output and rewrite its sidecar; the parquet is not touched")
    ap.add_argument("--no-verify", action="store_true")
    a = ap.parse_args()
    y = [int(v) for v in a.years.split("-")]
    years = list(range(y[0], y[-1] + 1))
    variables = [v for v in a.vars.split(",") if v]
    schema = pa.schema([("site_id", pa.string()), ("Date", pa.date32())] + [(v, pa.float64()) for v in variables])
    side = a.out + ".provenance.json"

    if a.provenance_only:
        if not (os.path.exists(a.out) and os.path.exists(side)):
            sys.exit(f"--provenance-only needs {a.out} and {side}")
        prov = json.load(open(side))
        md5 = md5sum(a.out)
        if md5 != prov.get("out_md5") or os.path.getsize(a.out) != prov.get("out_bytes"):
            sys.exit(f"{a.out}: md5/size differ from {side}; it is not the file the sidecar describes")
        now = {}
        for year in years:
            for v in variables:
                p = os.path.join(a.in_dir, f"{v}_{year}.parquet")
                now[os.path.basename(p)] = md5sum(p)
        if now != prov.get("inputs"):
            sys.exit(f"the per-variable-year inputs in {a.in_dir} differ from those {side} records")
        prov["sidecar_rewritten_utc"] = utc_now()
        prov["sidecar_rewritten_by"] = {"git": git_state(), "software": software_versions()}
    else:
        prov = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "created_utc": utc_now(), "years": years,
                "vars": variables, "inputs": {}, "rows": 0, "nan_days": {v: 0 for v in variables},
                "weights_md5": done_weights(a.in_dir, years, variables)}
        sites_ref = None
        tmp = a.out + ".tmp"
        writer, written = pq.ParquetWriter(tmp, schema, compression="zstd"), False
        try:
            for year in years:
                frame, sites = load_year(a.in_dir, year, variables, md5s=prov["inputs"])
                if sites_ref is None:
                    sites_ref = sites
                elif sites != sites_ref:
                    sys.exit(f"{year}: site set differs from {years[0]}")
                for v in variables:
                    prov["nan_days"][v] += int(frame[v].isna().sum())
                # NaN stays NaN (pandas conversion would store it as NULL)
                cols = [pa.array(frame["site_id"].to_numpy(), pa.string()),
                        pa.array(frame["Date"].to_numpy().astype("datetime64[D]"), pa.date32())]
                cols += [pa.array(frame[v].to_numpy(dtype=np.float64), pa.float64()) for v in variables]
                writer.write_table(pa.Table.from_arrays(cols, schema=schema))
                prov["rows"] += len(frame)
                print(f"  {year}: {len(frame):,} rows, {len(sites):,} sites", flush=True)
            written = True
        finally:
            writer.close()
            if not written and os.path.exists(tmp):
                os.remove(tmp)
        try:
            ver = check_output(tmp, a, years, variables)
        except BaseException:
            os.remove(tmp)      # an unverified file never takes the output's name
            raise
        os.replace(tmp, a.out)
        prov.update({"sites": len(sites_ref), "out": os.path.basename(a.out), "out_bytes": os.path.getsize(a.out),
                     "out_md5": md5sum(a.out), "git": git_state(), "software": software_versions()})
    if a.provenance_only:
        ver = check_output(a.out, a, years, variables)
    prov["row_order"] = "year, site_id, Date (written year by year)"
    if ver:
        prov["verification"] = ver
    prov["run"] = run_provenance(a.in_dir, years, variables)
    write_json(side, prov)
    print(f"{'sidecar rewritten for' if a.provenance_only else 'wrote'} {a.out}: {prov['rows']:,} rows, "
          f"{prov['sites']:,} sites, {prov['out_bytes'] / 1e6:.1f} MB, NaN days {prov['nan_days']}")


if __name__ == "__main__":
    main()
