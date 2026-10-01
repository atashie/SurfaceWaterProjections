#!/usr/bin/env python3
"""Run a Julia benchmark wrapper in gage batches to bound memory, then merge the outputs.

The runner (`run_julia_benchmark.jl`) loads the whole streamflow and climate parquets
before it does anything per gage. The full inputs need well over 8 GB of RAM, more than
the 16 GB laptop can spare while other work runs. Everything after loading is per gage:
each gage's climate is left-joined onto its own streamflow dates, then its years are
filtered, it is qualified, and its signatures are computed. So the run can be split by
gage into batches with their own small input files and merged afterwards.

The one cross-gage quantity is `flagged_for_high_na`: its denominator is the set of
signature columns PRESENT in the assembled table. `merge` therefore requires every batch
to write the identical header, so each batch saw the same column set and the flag
equals that of an unbatched run. If the headers differ, it stops and says so.

Steps:
  split  --streamflow SF.parquet --batches K --out-dir IN [--climate NAME=PATH[,max_date=YYYY-MM-DD]
         [,sites_with_data_in=PATH]] ...
      Gages (sorted ids) go to batch rank % K.
      Writes:
        * streamflow_b<k>.parquet: every column, rows in the original file order;
        * climate_<NAME>_b<k>.parquet: site_id, Date, prcp, swe, the four columns the
          runner keeps;
        * batches.json.
      Climate filters: max_date drops later rows; sites_with_data_in keeps only sites
      that have a non-NaN prcp value in another file.
  run    --wrapper W.jl --project JULIA_PROJECT --inputs IN --climate NAME --out-dir OUT
         --metadata META.csv [--gages-ii-dir DIR] [--min-avail-gb 1.5]
         [--heap-size-hint 2500M] [--julia ~/.juliaup/bin/julia]
      Runs the batches one at a time (STREAMFLOW_* env overrides, STREAMFLOW_HASH_INPUTS=1).
      Every 2 s it samples the process RSS and the system's available memory (free +
      inactive + speculative pages). If available memory falls below --min-avail-gb it
      terminates Julia and stops, protecting other work on the machine. Writes
      OUT/batch_runs.json (exit code, seconds, peak RSS, lowest available memory).
  merge  --out-dir OUT --prefix PREFIX [--reference REF.csv]
      Batch CSVs -> OUT/<PREFIX>_signatures.csv. Merging is text-level: values are not
      re-serialized and the header is checked identical. Rows follow the reference's gage
      order when given. The annual parquets are merged and sorted by (gage_id, signature,
      water_year) into OUT/<PREFIX>_signatures_annual.parquet. Writes OUT/merge.json and
      OUT/<PREFIX>_timing.json (summed counts, batch 1's provenance) for
      validate_production_run.py.
"""
import argparse
import csv
import glob
import json
import os
import re
import subprocess
import sys
import time

import duckdb
import pyarrow.parquet as pq

CLIMATE_COLS = ["site_id", "Date", "prcp", "swe"]


def con(mem="2GB", threads=4):
    c = duckdb.connect()
    c.execute(f"SET memory_limit='{mem}'; SET threads={threads}; SET preserve_insertion_order=true")
    return c


def write_query(c, query, path, params=None):
    """Stream a query result into a zstd parquet (pyarrow writer, bounded memory)."""
    reader = c.execute(query, params or []).to_arrow_reader(1_000_000)
    n = 0
    with pq.ParquetWriter(path + ".tmp", reader.schema, compression="zstd") as w:
        for b in reader:
            w.write_batch(b)
            n += b.num_rows
    os.replace(path + ".tmp", path)
    return n


def split(a):
    os.makedirs(a.out_dir, exist_ok=True)
    c = con()
    gages = [r[0] for r in c.execute(f"SELECT DISTINCT gage_id FROM read_parquet('{a.streamflow}') "
                                     f"ORDER BY gage_id").fetchall()]
    batches = [gages[k::a.batches] for k in range(a.batches)]
    meta = {"streamflow": os.path.abspath(a.streamflow), "batches": a.batches, "gages": len(gages),
            "assignment": "sorted gage_id, rank % K", "files": {}}
    climates = []
    for spec in a.climate or []:
        name, _, rest = spec.partition("=")
        parts = rest.split(",")
        opts = dict(p.split("=", 1) for p in parts[1:])
        climates.append((name, parts[0], opts))
    keep_sites = {}
    for name, path, opts in climates:
        if "sites_with_data_in" in opts:
            keep_sites[name] = [r[0] for r in c.execute(
                f"SELECT site_id FROM read_parquet('{opts['sites_with_data_in']}') GROUP BY site_id "
                f"HAVING bool_or(prcp IS NOT NULL AND NOT isnan(prcp))").fetchall()]
    for k, ids in enumerate(batches, 1):
        c.execute("CREATE OR REPLACE TEMP TABLE b AS SELECT unnest(?::VARCHAR[]) AS gage_id", [ids])
        p = os.path.join(a.out_dir, f"streamflow_b{k}.parquet")
        n = write_query(c, f"SELECT * EXCLUDE (file_row_number) FROM read_parquet('{a.streamflow}', "
                           f"file_row_number=true) WHERE gage_id IN (SELECT gage_id FROM b) "
                           f"ORDER BY file_row_number", p)
        meta["files"][os.path.basename(p)] = {"gages": len(ids), "rows": n}
        for name, path, opts in climates:
            where = ["site_id IN (SELECT gage_id FROM b)"]
            if "max_date" in opts:
                where.append(f"Date <= DATE '{opts['max_date']}'")
            if name in keep_sites:
                c.execute("CREATE OR REPLACE TEMP TABLE ks AS SELECT unnest(?::VARCHAR[]) AS site_id", [keep_sites[name]])
                where.append("site_id IN (SELECT site_id FROM ks)")
            q = (f"SELECT {', '.join(CLIMATE_COLS)} FROM read_parquet('{path}') WHERE {' AND '.join(where)} "
                 f"ORDER BY site_id, Date")
            cp = os.path.join(a.out_dir, f"climate_{name}_b{k}.parquet")
            n = write_query(c, q, cp)
            sites = c.execute(f"SELECT count(DISTINCT site_id) FROM read_parquet('{cp}')").fetchone()[0]
            meta["files"][os.path.basename(cp)] = {"source": os.path.abspath(path), "options": opts,
                                                   "rows": n, "sites": sites}
        print(f"batch {k}: {len(ids)} gages, {meta['files'][os.path.basename(p)]['rows']:,} streamflow rows"
              + "".join(f"; {nm} {meta['files'][f'climate_{nm}_b{k}.parquet']['sites']} sites"
                        for nm, _, _ in climates), flush=True)
    meta["batch_gages"] = batches
    with open(os.path.join(a.out_dir, "batches.json"), "w") as fh:
        json.dump(meta, fh, indent=1)


def available_gb():
    out = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = int(re.search(r"page size of (\d+) bytes", out).group(1))
    pages = {m.group(1): int(m.group(2)) for m in re.finditer(r"^Pages (\w[\w ]*?):\s+(\d+)\.", out, re.M)}
    return (pages.get("free", 0) + pages.get("inactive", 0) + pages.get("speculative", 0)) * page / 2**30


def rss_gb(pid):
    r = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    return int(r) / 2**20 if r else 0.0


def run(a):
    # Julia runs with its cwd in the checkout, so every path it receives must be absolute
    for k in ("wrapper", "project", "inputs", "out_dir", "metadata", "gages_ii_dir"):
        if getattr(a, k):
            setattr(a, k, os.path.abspath(os.path.expanduser(getattr(a, k))))
    os.makedirs(a.out_dir, exist_ok=True)
    k_all = sorted(int(m.group(1)) for f in glob.glob(os.path.join(a.inputs, "streamflow_b*.parquet"))
                   if (m := re.search(r"_b(\d+)\.parquet$", f)))
    log_path = os.path.join(a.out_dir, "batch_runs.json")
    runs = json.load(open(log_path)) if os.path.exists(log_path) else {}
    for k in k_all:
        if runs.get(str(k), {}).get("rc") == 0:
            print(f"batch {k}: done earlier, skipped", flush=True)
            continue
        out = os.path.join(a.out_dir, f"b{k}")
        os.makedirs(out, exist_ok=True)
        env = dict(os.environ, STREAMFLOW_DATA_PATH=os.path.join(a.inputs, f"streamflow_b{k}.parquet"),
                   STREAMFLOW_CLIMATE_PATH=os.path.join(a.inputs, f"climate_{a.climate}_b{k}.parquet"),
                   STREAMFLOW_METADATA_PATH=a.metadata, STREAMFLOW_OUTPUT_DIR=out, STREAMFLOW_HASH_INPUTS="1")
        if a.gages_ii_dir:
            env["STREAMFLOW_GAGES_II_DIR"] = a.gages_ii_dir
        cmd = [os.path.expanduser(a.julia), f"--project={a.project}", f"--heap-size-hint={a.heap_size_hint}", a.wrapper]
        t0, peak, low = time.time(), 0.0, available_gb()
        with open(os.path.join(out, "run.log"), "w") as log:
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=os.path.dirname(a.project))
            killed = False
            while p.poll() is None:
                peak, low = max(peak, rss_gb(p.pid)), min(low, available_gb())
                if low < a.min_avail_gb:
                    p.terminate()
                    try:
                        p.wait(30)
                    except subprocess.TimeoutExpired:
                        p.kill()
                    killed = True
                    break
                time.sleep(2)
        runs[str(k)] = {"rc": p.returncode if not killed else -999, "seconds": round(time.time() - t0, 1),
                        "peak_rss_gb": round(peak, 2), "min_available_gb": round(low, 2),
                        "killed_low_memory": killed, "cmd": cmd,
                        "env": {e: env[e] for e in env if e.startswith("STREAMFLOW_")}}
        with open(log_path, "w") as fh:
            json.dump(runs, fh, indent=1)
        r = runs[str(k)]
        print(f"batch {k}: rc {r['rc']}, {r['seconds']} s, peak RSS {r['peak_rss_gb']} GB, "
              f"lowest available {r['min_available_gb']} GB{' -- KILLED (low memory)' if killed else ''}", flush=True)
        if r["rc"] != 0:
            sys.exit(f"batch {k} failed; see {out}/run.log")


def merge(a):
    dirs = sorted((d for d in glob.glob(os.path.join(a.out_dir, "b*"))
                   if os.path.isdir(d) and re.fullmatch(r"b\d+", os.path.basename(d))),
                  key=lambda d: int(os.path.basename(d)[1:]))
    if not dirs:
        sys.exit(f"no batch folders b<k> in {a.out_dir}")
    header, rows, per_batch = None, {}, {}
    for d in dirs:
        f = glob.glob(os.path.join(d, "*_signatures.csv"))
        if len(f) != 1:
            sys.exit(f"{d}: expected one *_signatures.csv, found {f}")
        with open(f[0], newline="") as fh:
            h = fh.readline()
            if header is None:
                header = h
            elif h != header:
                hs, hb = next(csv.reader([header])), next(csv.reader([h]))
                sys.exit(f"{d}: header differs from batch 1 (only here: {sorted(set(hb) - set(hs))[:10]}; "
                         f"missing here: {sorted(set(hs) - set(hb))[:10]}) -- the flagged_for_high_na "
                         f"denominator would differ from an unbatched run")
            ncol = len(next(csv.reader([header])))
            n = 0
            for line in fh:
                fields = next(csv.reader([line]))
                if len(fields) != ncol:
                    sys.exit(f"{f[0]}: a line has {len(fields)} fields, the header {ncol}")
                if fields[0] in rows:
                    sys.exit(f"gage {fields[0]} appears in two batches")
                rows[fields[0]] = line if line.endswith("\n") else line + "\n"
                n += 1
            per_batch[os.path.basename(d)] = n
    order = list(rows)
    missing, extra = [], []
    if a.reference:
        with open(a.reference, newline="") as fh:
            ref_ids = [r[0] for r in csv.reader(fh)][1:]
        missing = [g for g in ref_ids if g not in rows]
        extra = sorted(set(rows) - set(ref_ids))
        order = [g for g in ref_ids if g in rows] + extra
    out_csv = os.path.join(a.out_dir, f"{a.prefix}_signatures.csv")
    with open(out_csv + ".tmp", "w", newline="") as fh:
        fh.write(header)
        for g in order:
            fh.write(rows[g])
    os.replace(out_csv + ".tmp", out_csv)
    ann = [f for d in dirs for f in sorted(glob.glob(os.path.join(d, "*_signatures_annual.parquet")))]
    out_ann, n_ann = None, 0
    if ann:
        out_ann = os.path.join(a.out_dir, f"{a.prefix}_signatures_annual.parquet")
        c = con()
        files = "[" + ", ".join(f"'{p}'" for p in ann) + "]"
        n_ann = write_query(c, f"SELECT * FROM read_parquet({files}) ORDER BY gage_id, signature, water_year", out_ann)
    timing = {}
    for d in dirs:
        t = glob.glob(os.path.join(d, "*_timing.json"))
        if t:
            j = json.load(open(t[0]))
            timing[os.path.basename(d)] = {"total_seconds": j.get("total_seconds"),
                                           "n_gages_processed": j.get("n_gages_processed"),
                                           "n_annual_rows": j.get("n_annual_rows"),
                                           "git_revision": j.get("provenance", {}).get("git_revision"),
                                           "git_working_tree_dirty": j.get("provenance", {}).get("git_working_tree_dirty")}
    first = glob.glob(os.path.join(dirs[0], "*_timing.json"))
    merged_timing = {"merged_from_batches": len(dirs), "n_annual_rows": n_ann,
                     "n_gages_processed": sum(v["n_gages_processed"] or 0 for v in timing.values()),
                     "total_seconds": round(sum(v["total_seconds"] or 0 for v in timing.values()), 1),
                     "provenance": json.load(open(first[0])).get("provenance") if first else None,
                     "batches": timing}
    with open(os.path.join(a.out_dir, f"{a.prefix}_timing.json"), "w") as fh:
        json.dump(merged_timing, fh, indent=1)
    res = {"batches": per_batch, "gages": len(rows), "columns": len(next(csv.reader([header]))),
           "reference": a.reference, "missing_vs_reference": missing, "extra_vs_reference": extra,
           "signatures_csv": out_csv, "annual_parquet": out_ann, "annual_rows": n_ann, "timing": timing}
    with open(os.path.join(a.out_dir, "merge.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"merged {len(dirs)} batches: {len(rows):,} gages x {res['columns']} columns -> {out_csv}; "
          f"annual rows {n_ann:,}; vs reference: {len(missing)} missing, {len(extra)} extra")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--streamflow", required=True)
    s.add_argument("--batches", type=int, required=True)
    s.add_argument("--out-dir", required=True)
    s.add_argument("--climate", action="append")
    r = sub.add_parser("run")
    r.add_argument("--wrapper", required=True)
    r.add_argument("--project", required=True)
    r.add_argument("--inputs", required=True)
    r.add_argument("--climate", required=True)
    r.add_argument("--out-dir", required=True)
    r.add_argument("--metadata", required=True)
    r.add_argument("--gages-ii-dir", default=None)
    r.add_argument("--min-avail-gb", type=float, default=1.5)
    r.add_argument("--heap-size-hint", default="2500M")
    r.add_argument("--julia", default="~/.juliaup/bin/julia")
    m = sub.add_parser("merge")
    m.add_argument("--out-dir", required=True)
    m.add_argument("--prefix", required=True)
    m.add_argument("--reference", default=None)
    a = ap.parse_args()
    {"split": split, "run": run, "merge": merge}[a.cmd](a)


if __name__ == "__main__":
    main()
