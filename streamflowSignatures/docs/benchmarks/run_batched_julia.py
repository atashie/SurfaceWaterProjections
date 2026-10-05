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

Use a wrapper that honours STREAMFLOW_OUTPUT_DIR, such as the standard-product wrappers
`run_julia_benchmark_drought_1993_2025_60pct.jl` and
`run_julia_benchmark_prod_1980_2025_60pct_drought.jl`. A wrapper that assigns it outright
writes elsewhere; `run` then finds no fresh outputs in the batch folder and stops.
Needs duckdb and pyarrow (e.g. ~/HISSS_data/venvs/daymet; the repo .venv has no duckdb).

Steps:
  split  --streamflow SF.parquet --batches K --out-dir IN [--climate NAME=PATH[,max_date=YYYY-MM-DD]
         [,sites_with_data_in=PATH]] ...
      Gages (sorted ids) go to batch rank % K.
      Writes:
        * streamflow_b<k>.parquet: every column, rows in the original file order;
        * climate_<NAME>_b<k>.parquet: site_id, Date, prcp, swe, the four columns the
          runner keeps;
        * batches.json: the assignment, plus size, mtime, PAR1 footer and sha256 of every
          source file. These become the merged run's input provenance.
      Climate filters: max_date drops later rows; sites_with_data_in keeps only sites
      that have a non-NaN prcp value in another file. Any other option is an error.
  run    --wrapper W.jl --project JULIA_PROJECT --inputs IN --climate NAME --out-dir OUT
         --metadata META.csv [--gages-ii-dir DIR] [--min-avail-gb 1.5] [--max-minutes 120]
         [--heap-size-hint 2500M] [--julia ~/.juliaup/bin/julia]
      Checks that every input exists (Julia would skip a missing climate or metadata file
      silently), then runs the K batches of IN/batches.json one at a time. It sets the
      STREAMFLOW_* env overrides and STREAMFLOW_HASH_INPUTS=1.
      Julia runs in its own process group. Every 2 s the guard samples its RSS and the
      system's available memory (free + inactive + speculative pages). It kills the whole
      group and stops when:
        * available memory falls below --min-avail-gb, or a sample fails;
        * a batch exceeds --max-minutes;
        * this script is interrupted or terminated.
      A batch is done only with exit code 0 and one fresh signatures CSV and timing JSON in
      its folder. A rerun skips the done batches whose command and environment are
      unchanged.
      Writes OUT/batch_runs.json (exit code, seconds, peak RSS, lowest available memory,
      earlier attempts) and appends every attempt to OUT/b<k>/run.log.
  merge  --out-dir OUT --prefix PREFIX [--reference REF.csv [--allow-missing]] [--memory-limit 1GB]
      Refuses an incomplete or inconsistent set of batches. It requires every batch of
      batches.json, found through OUT/batch_runs.json, to be done: folders b1..bK, exit
      code 0, one CSV and one timing JSON each, and either one annual parquet in every
      batch or none. Each batch's CSV gages must belong to that batch. Its rows must equal
      the timing's n_gages_processed, and its annual rows n_annual_rows. All batches must
      share git revision, working-tree state, Julia version, config, metadata and env
      overrides.
      Batch CSVs -> OUT/<PREFIX>_signatures.csv. Merging is text-level: values are not
      re-serialized and the header is checked identical. Rows follow the reference's gage
      order. A gage in only one of the reference and the batches is an error unless
      --allow-missing is given. Without a reference the rows are sorted by gage_id; an
      unbatched run's order is an artifact of its joins.
      The annual parquets are merged and sorted by (gage_id, signature, water_year) into
      OUT/<PREFIX>_signatures_annual.parquet. Also writes OUT/<PREFIX>_merge.json and
      OUT/<PREFIX>_timing.json for validate_production_run.py: summed counts, the shared
      provenance, and the source files recorded by split in place of the batch inputs.
"""
import argparse
import collections
import csv
import glob
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import time

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

CLIMATE_COLS = ["site_id", "Date", "prcp", "swe"]
CLIMATE_OPTIONS = {"max_date", "sites_with_data_in"}
# provenance every batch must share; the streamflow and climate entries are per batch
SHARED_PROVENANCE = ["git_revision", "git_working_tree_dirty", "julia_version", "hostname", "config", "metadata",
                     "env_overrides"]
SHARED_COUNTS = ["language", "n_signature_columns", "n_metadata_columns", "n_qaqc_flags"]


def con(mem="2GB", threads=4, temp_dir=None):
    c = duckdb.connect()
    c.execute(f"SET memory_limit='{mem}'; SET threads={threads}; SET preserve_insertion_order=true")
    if temp_dir:
        c.execute(f"SET temp_directory='{temp_dir}'")
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


def write_sorted_annual(files, path, mem, temp_dir, chunk_rows=2_000_000, group_rows=1_000_000):
    """Merge annual parquets sorted by (gage_id, signature, water_year) in bounded memory.

    Batches hold disjoint gages, so the gages are sorted in chunks of about `chunk_rows`
    rows, one chunk after another. The output is re-sliced into row groups of `group_rows`,
    the same file a single global ORDER BY would write."""
    c = con(mem, 2, temp_dir)
    flist = "[" + ", ".join(f"'{p}'" for p in files) + "]"
    counts = dict(c.execute(f"SELECT gage_id, count(*) FROM read_parquet({flist}) GROUP BY gage_id").fetchall())
    chunks, cur, n = [], [], 0
    for g in sorted(counts):
        cur.append(g)
        n += counts[g]
        if n >= chunk_rows:
            chunks.append(cur)
            cur, n = [], 0
    if cur:
        chunks.append(cur)
    writer, buf, nbuf, total = None, [], 0, 0
    try:
        for ch in chunks:
            c.execute("CREATE OR REPLACE TEMP TABLE g AS SELECT unnest(?::VARCHAR[]) AS gage_id", [ch])
            reader = c.execute(f"SELECT * FROM read_parquet({flist}) WHERE gage_id IN (SELECT gage_id FROM g) "
                               f"ORDER BY gage_id, signature, water_year").to_arrow_reader(group_rows)
            for b in reader:
                if writer is None:
                    writer = pq.ParquetWriter(path + ".tmp", b.schema, compression="zstd")
                buf.append(b)
                nbuf += b.num_rows
                while nbuf >= group_rows:
                    t = pa.Table.from_batches(buf)
                    writer.write_table(t.slice(0, group_rows))
                    rest = t.slice(group_rows)
                    buf, nbuf = rest.to_batches(), rest.num_rows
                    total += group_rows
        if nbuf:
            writer.write_table(pa.Table.from_batches(buf))
            total += nbuf
    finally:
        if writer is not None:
            writer.close()
        c.close()
    os.replace(path + ".tmp", path)
    return total


def file_provenance(path):
    """Size, mtime, sha256 and (for parquet) the PAR1 footer of a source file."""
    p = os.path.abspath(path)
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
        fh.seek(-4, os.SEEK_END)
        tail = fh.read(4)
    st = os.stat(p)
    return {"path": p, "bytes": st.st_size, "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": h.hexdigest(), "par1_footer": (tail == b"PAR1") if p.endswith(".parquet") else None}


def split(a):
    climates = []
    for spec in a.climate or []:
        name, _, rest = spec.partition("=")
        parts = rest.split(",")
        opts = {}
        for o in parts[1:]:
            key, eq, val = o.partition("=")
            if not eq or key not in CLIMATE_OPTIONS:
                sys.exit(f"--climate {spec}: unknown option '{o}' (allowed: {sorted(CLIMATE_OPTIONS)})")
            opts[key] = val
        if not name or not parts[0]:
            sys.exit(f"--climate {spec}: expected NAME=PATH[,option=value...]")
        climates.append((name, parts[0], opts))
    os.makedirs(a.out_dir, exist_ok=True)
    print("hashing the source files ...", flush=True)
    sources = {"streamflow": file_provenance(a.streamflow), "climate": {}}
    for name, path, opts in climates:
        sources["climate"][name] = dict(file_provenance(path), options=opts)
        if "sites_with_data_in" in opts:
            sources["climate"][name]["sites_with_data_in"] = file_provenance(opts["sites_with_data_in"])
    entries = [sources["streamflow"]] + [e for s in sources["climate"].values()
                                          for e in (s, s.get("sites_with_data_in")) if e]
    bad = [e["path"] for e in entries if e["par1_footer"] is False]
    if bad:
        sys.exit(f"no PAR1 footer, so the file is truncated or not parquet: {bad}")
    c = con()
    gages = [r[0] for r in c.execute(f"SELECT DISTINCT gage_id FROM read_parquet('{a.streamflow}') "
                                     f"ORDER BY gage_id").fetchall()]
    batches = [gages[k::a.batches] for k in range(a.batches)]
    meta = {"streamflow": os.path.abspath(a.streamflow), "batches": a.batches, "gages": len(gages),
            "assignment": "sorted gage_id, rank % K", "sources": sources, "files": {}}
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
    out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=10).stdout
    page = int(re.search(r"page size of (\d+) bytes", out).group(1))
    pages = {m.group(1): int(m.group(2)) for m in re.finditer(r"^Pages (\w[\w ]*?):\s+(\d+)\.", out, re.M)}
    return (pages.get("free", 0) + pages.get("inactive", 0) + pages.get("speculative", 0)) * page / 2**30


def rss_gb(pid):
    r = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    return int(r) / 2**20 if r else 0.0


def batch_outputs(d):
    return {kind: sorted(glob.glob(os.path.join(d, pat))) for kind, pat in
            (("csv", "*_signatures.csv"), ("annual", "*_signatures_annual.parquet"), ("timing", "*_timing.json"))}


def outputs_ok(d, since=None):
    """One signatures CSV and one timing JSON (written after `since`, if given), at most one annual parquet."""
    o = batch_outputs(d)
    if len(o["csv"]) != 1 or len(o["timing"]) != 1 or len(o["annual"]) > 1:
        return False
    return since is None or all(os.path.getmtime(f) >= since for f in o["csv"] + o["timing"])


def kill_group(p, grace=30):
    """SIGTERM the batch's whole process group, then SIGKILL whatever is left after `grace` seconds."""
    try:
        os.killpg(p.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        p.wait(grace)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(p.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    p.wait()


def on_signal(signum, frame):
    raise SystemExit(f"terminated by signal {signum}")


def run(a):
    # Julia runs with its cwd in the checkout, so every path it receives must be absolute
    for k in ("wrapper", "project", "inputs", "out_dir", "metadata", "gages_ii_dir", "julia"):
        if getattr(a, k):
            setattr(a, k, os.path.abspath(os.path.expanduser(getattr(a, k))))
    bj = os.path.join(a.inputs, "batches.json")
    if not os.path.isfile(bj):
        sys.exit(f"{bj} not found: run split first")
    n_batches = json.load(open(bj))["batches"]
    missing = [p for p in [a.wrapper, a.metadata, a.julia] + [
        os.path.join(a.inputs, f) for k in range(1, n_batches + 1)
        for f in (f"streamflow_b{k}.parquet", f"climate_{a.climate}_b{k}.parquet")] if not os.path.isfile(p)]
    missing += [p for p in [a.project] + ([a.gages_ii_dir] if a.gages_ii_dir else []) if not os.path.isdir(p)]
    if missing:
        sys.exit("missing inputs (Julia would skip a missing climate or metadata file silently):\n  "
                 + "\n  ".join(missing))
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGHUP, on_signal)
    os.makedirs(a.out_dir, exist_ok=True)
    log_path = os.path.join(a.out_dir, "batch_runs.json")
    runs = json.load(open(log_path)) if os.path.exists(log_path) else {}
    for k in range(1, n_batches + 1):
        out = os.path.join(a.out_dir, f"b{k}")
        env = dict(os.environ, STREAMFLOW_DATA_PATH=os.path.join(a.inputs, f"streamflow_b{k}.parquet"),
                   STREAMFLOW_CLIMATE_PATH=os.path.join(a.inputs, f"climate_{a.climate}_b{k}.parquet"),
                   STREAMFLOW_METADATA_PATH=a.metadata, STREAMFLOW_OUTPUT_DIR=out, STREAMFLOW_HASH_INPUTS="1")
        if a.gages_ii_dir:
            env["STREAMFLOW_GAGES_II_DIR"] = a.gages_ii_dir
        envrec = {e: env[e] for e in env if e.startswith("STREAMFLOW_")}
        cmd = [a.julia, f"--project={a.project}", f"--heap-size-hint={a.heap_size_hint}", a.wrapper]
        prev = runs.get(str(k))
        if prev and prev.get("rc") == 0 and prev.get("cmd") == cmd and prev.get("env") == envrec and outputs_ok(out):
            print(f"batch {k}: done earlier with the same command and environment, skipped", flush=True)
            continue
        os.makedirs(out, exist_ok=True)
        t0, peak, low, why, interrupted = time.time(), 0.0, float("inf"), None, None
        with open(os.path.join(out, "run.log"), "a") as log:
            log.write(f"\n===== attempt started {time.strftime('%Y-%m-%dT%H:%M:%S')}: {' '.join(cmd)}\n")
            log.flush()
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=os.path.dirname(a.project),
                                 start_new_session=True)
            try:
                while p.poll() is None:
                    peak = max(peak, rss_gb(p.pid))
                    try:
                        avail = available_gb()
                    except Exception as e:      # no reading, so assume the worst
                        avail, why = 0.0, f"memory sampling failed ({type(e).__name__}: {e})"
                    low = min(low, avail)
                    if why is None and avail < a.min_avail_gb:
                        why = f"available memory {avail:.2f} GB fell below --min-avail-gb {a.min_avail_gb}"
                    if why is None and a.max_minutes and time.time() - t0 > 60 * a.max_minutes:
                        why = f"exceeded --max-minutes {a.max_minutes}"
                    if why:
                        break
                    time.sleep(2)
            except BaseException as e:          # Ctrl-C, SIGTERM/SIGHUP or a bug: never leave Julia running
                interrupted, why = e, f"interrupted ({type(e).__name__}: {e})"
            finally:
                if p.poll() is None:
                    kill_group(p)
        rc = p.returncode if why is None else -999
        rec = {"rc": rc, "seconds": round(time.time() - t0, 1), "peak_rss_gb": round(peak, 2),
               "min_available_gb": round(low, 2) if low != float("inf") else None, "stopped": why,
               "cmd": cmd, "env": envrec}
        if rc == 0 and not outputs_ok(out, since=t0 - 1):
            rec["rc"] = -998
            rec["stopped"] = ("exit 0, but the batch folder holds no fresh signatures CSV + timing JSON "
                              "(does the wrapper honour STREAMFLOW_OUTPUT_DIR?)")
        if prev:
            rec["earlier_attempts"] = prev.pop("earlier_attempts", []) + [prev]
        runs[str(k)] = rec
        with open(log_path + ".tmp", "w") as fh:
            json.dump(runs, fh, indent=1)
        os.replace(log_path + ".tmp", log_path)
        print(f"batch {k}: rc {rec['rc']}, {rec['seconds']} s, peak RSS {rec['peak_rss_gb']} GB, lowest available "
              f"{rec['min_available_gb']} GB" + (f" -- STOPPED: {rec['stopped']}" if rec["stopped"] else ""), flush=True)
        if interrupted is not None:
            raise interrupted
        if rec["rc"] != 0:
            sys.exit(f"batch {k} failed; see {out}/run.log")


def merge(a):
    a.out_dir = os.path.abspath(a.out_dir)
    rp = os.path.join(a.out_dir, "batch_runs.json")
    if not os.path.isfile(rp):
        sys.exit(f"{rp} not found: merge reads the run log to find and check the batches")
    runs = json.load(open(rp))
    inputs = {os.path.dirname(r["env"]["STREAMFLOW_DATA_PATH"]) for r in runs.values()}
    climates = {os.path.basename(r["env"]["STREAMFLOW_CLIMATE_PATH"]).rsplit("_b", 1)[0] for r in runs.values()}
    if len(inputs) != 1 or len(climates) != 1:
        sys.exit(f"batch_runs.json mixes input folders {sorted(inputs)} or climates {sorted(climates)}")
    bj_path = os.path.join(inputs.pop(), "batches.json")
    bj = json.load(open(bj_path))
    n_batches, batch_gages = bj["batches"], bj["batch_gages"]
    climate_name = climates.pop()[len("climate_"):]
    dirs = [os.path.join(a.out_dir, f"b{k}") for k in range(1, n_batches + 1)]
    problems = []
    for k, d in enumerate(dirs, 1):
        r = runs.get(str(k))
        if not os.path.isdir(d):
            problems.append(f"b{k}: folder missing")
        elif not r or r.get("rc") != 0:
            problems.append(f"b{k}: not done (rc {r.get('rc') if r else 'none'})")
        elif not outputs_ok(d):
            problems.append(f"b{k}: expected one signatures CSV, one timing JSON and at most one annual parquet; "
                            f"found {batch_outputs(d)}")
    extra = sorted(d for d in glob.glob(os.path.join(a.out_dir, "b*"))
                   if os.path.isdir(d) and re.fullmatch(r"b\d+", os.path.basename(d)) and d not in dirs)
    if extra:
        problems.append(f"batch folders beyond the {n_batches} of batches.json: {extra}")
    if problems:
        sys.exit("incomplete or inconsistent batches:\n  " + "\n  ".join(problems))
    outs = [batch_outputs(d) for d in dirs]
    if len({len(o["annual"]) for o in outs}) > 1:
        sys.exit("some batches wrote an annual parquet and some did not")
    timings = [json.load(open(o["timing"][0])) for o in outs]
    p0 = timings[0].get("provenance") or {}
    for k, t in enumerate(timings[1:], 2):
        diff = [f for f in SHARED_PROVENANCE if (t.get("provenance") or {}).get(f) != p0.get(f)]
        if diff:
            problems.append(f"b{k}: provenance differs from b1 in {diff}")
    for f in SHARED_COUNTS:
        if len({json.dumps(t.get(f)) for t in timings}) > 1:
            problems.append(f"{f} differs between batches")
    if problems:
        sys.exit("batches from different code, config or inputs:\n  " + "\n  ".join(problems))

    header, rows, where, per_batch = None, {}, {}, {}
    for k, (o, t) in enumerate(zip(outs, timings), 1):
        allowed = set(batch_gages[k - 1])
        with open(o["csv"][0], newline="") as fh:
            h, body = fh.readline(), fh.readlines()
        if header is None:
            header = h
        elif h != header:
            hs, hb = next(csv.reader([header])), next(csv.reader([h]))
            empty = [f"b{j}" for j, n in ((1, per_batch.get("b1")), (k, len(body))) if n == 0]
            sys.exit(f"b{k}: header differs from b1 (only here: {sorted(set(hb) - set(hs))[:10]}; missing here: "
                     f"{sorted(set(hs) - set(hb))[:10]})"
                     + (f" -- {', '.join(empty)} wrote no gages, a header-only CSV" if empty else
                        " -- the flagged_for_high_na denominator would differ from an unbatched run"))
        ncol = len(next(csv.reader([header])))
        for line in body:
            fields = next(csv.reader([line]))
            g = fields[0]
            if len(fields) != ncol:
                sys.exit(f"{o['csv'][0]}: a line has {len(fields)} fields, the header {ncol}")
            if g in rows:
                sys.exit(f"gage {g} appears twice in b{k}" if where[g] == k else f"gage {g} appears in b{where[g]} and b{k}")
            if g not in allowed:
                sys.exit(f"b{k}: gage {g} is not one of its batch's gages")
            rows[g], where[g] = (line if line.endswith("\n") else line + "\n"), k
        per_batch[f"b{k}"] = len(body)
        if len(body) != t.get("n_gages_processed"):
            sys.exit(f"b{k}: {len(body)} CSV rows but n_gages_processed {t.get('n_gages_processed')} in its timing JSON")

    missing_ref, extra_ref = [], []
    if a.reference:
        with open(a.reference, newline="") as fh:
            ref_ids = [r[0] for r in csv.reader(fh)][1:]
        dup = sorted(g for g, n in collections.Counter(ref_ids).items() if n > 1)
        if dup:
            sys.exit(f"the reference lists {len(dup)} gages more than once, e.g. {dup[:5]}")
        missing_ref = [g for g in ref_ids if g not in rows]
        extra_ref = sorted(set(rows) - set(ref_ids))
        if (missing_ref or extra_ref) and not a.allow_missing:
            sys.exit(f"vs the reference: {len(missing_ref)} of its gages are not in the batches (e.g. "
                     f"{missing_ref[:5]}) and {len(extra_ref)} batch gages are not in it (e.g. {extra_ref[:5]}); "
                     f"pass --allow-missing to merge anyway")
        order = [g for g in ref_ids if g in rows] + extra_ref
    else:
        order = sorted(rows)
        print("no --reference: rows sorted by gage_id", flush=True)
    out_csv = os.path.join(a.out_dir, f"{a.prefix}_signatures.csv")
    with open(out_csv + ".tmp", "w", newline="") as fh:
        fh.write(header)
        for g in order:
            fh.write(rows[g])
    os.replace(out_csv + ".tmp", out_csv)

    out_ann, n_ann = None, 0
    if outs[0]["annual"]:
        ann = [o["annual"][0] for o in outs]
        for k, (f, t) in enumerate(zip(ann, timings), 1):
            nr = pq.ParquetFile(f).metadata.num_rows
            if nr != t.get("n_annual_rows"):
                sys.exit(f"b{k}: {nr} annual rows but n_annual_rows {t.get('n_annual_rows')} in its timing JSON")
        out_ann = os.path.join(a.out_dir, f"{a.prefix}_signatures_annual.parquet")
        tmp = os.path.join(a.out_dir, ".duckdb_tmp")
        n_ann = write_sorted_annual(ann, out_ann, a.memory_limit, tmp)
        try:
            os.rmdir(tmp)
        except OSError:
            pass
        if n_ann != sum(t["n_annual_rows"] for t in timings):
            sys.exit(f"merged {n_ann} annual rows, the batches report {sum(t['n_annual_rows'] for t in timings)}")

    summary = {f"b{k}": {"total_seconds": t.get("total_seconds"), "n_gages_processed": t.get("n_gages_processed"),
                         "n_annual_rows": t.get("n_annual_rows")} for k, t in enumerate(timings, 1)}
    sources = bj.get("sources") or {}
    note = "not recorded: batches.json was written before split recorded source provenance"
    prov = dict(p0)
    prov["streamflow"] = sources.get("streamflow") or {"path": bj.get("streamflow"), "note": note}
    prov["climate"] = (sources.get("climate") or {}).get(climate_name) or {
        "path": bj["files"].get(f"climate_{climate_name}_b1.parquet", {}).get("source"), "note": note}
    prov["batch_inputs"] = {f"b{k}": {f: (t.get("provenance") or {}).get(f) for f in ("streamflow", "climate")}
                            for k, t in enumerate(timings, 1)}
    phases = timings[0].get("phases") or {}
    merged_timing = {
        "language": timings[0].get("language"),
        "start_time": min((t.get("start_time") for t in timings if t.get("start_time")), default=None),
        "end_time": max((t.get("end_time") for t in timings if t.get("end_time")), default=None),
        "total_seconds": round(sum(t.get("total_seconds") or 0 for t in timings), 1),
        "phases": {ph: round(sum((t.get("phases") or {}).get(ph, 0) for t in timings), 3) for ph in phases},
        "n_gages_processed": len(rows), "n_annual_rows": n_ann,
        **{f: timings[0].get(f) for f in SHARED_COUNTS if f != "language"},
        "provenance": prov,
        "batched": {"batches": n_batches, "batches_json": bj_path, "climate": climate_name,
                    "note": "seconds and phases are summed over the batches, which ran one after another",
                    "per_batch": summary},
    }
    with open(os.path.join(a.out_dir, f"{a.prefix}_timing.json"), "w") as fh:
        json.dump(merged_timing, fh, indent=1)
    res = {"batches": per_batch, "gages": len(rows), "columns": len(next(csv.reader([header]))),
           "order": "reference" if a.reference else "sorted gage_id",
           "reference": os.path.abspath(a.reference) if a.reference else None,
           "missing_vs_reference": missing_ref, "extra_vs_reference": extra_ref,
           "signatures_csv": out_csv, "annual_parquet": out_ann, "annual_rows": n_ann, "timing": summary}
    with open(os.path.join(a.out_dir, f"{a.prefix}_merge.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"merged {n_batches} batches: {len(rows):,} gages x {res['columns']} columns -> {out_csv}; "
          f"annual rows {n_ann:,}; vs reference: {len(missing_ref)} missing, {len(extra_ref)} extra")


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
    r.add_argument("--max-minutes", type=float, default=120)
    r.add_argument("--heap-size-hint", default="2500M")
    r.add_argument("--julia", default="~/.juliaup/bin/julia")
    m = sub.add_parser("merge")
    m.add_argument("--out-dir", required=True)
    m.add_argument("--prefix", required=True)
    m.add_argument("--reference", default=None)
    m.add_argument("--allow-missing", action="store_true")
    m.add_argument("--memory-limit", default="1GB")
    a = ap.parse_args()
    {"split": split, "run": run, "merge": merge}[a.cmd](a)


if __name__ == "__main__":
    main()
