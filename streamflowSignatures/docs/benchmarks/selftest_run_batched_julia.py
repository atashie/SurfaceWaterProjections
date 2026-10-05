#!/usr/bin/env python3
"""Self-test of run_batched_julia.py on synthetic inputs with a fake julia: no Julia, no real data.

Covers:
  * split: its options and the source provenance it records;
  * run: resume, the input checks, a wrapper that ignores STREAMFLOW_OUTPUT_DIR, and the
    memory guard (low memory, a failed sample, --max-minutes, SIGTERM). Each guard must kill
    the whole process group, grandchild included;
  * merge: order, the reference, and refusal of incomplete or inconsistent batches.

Run with a Python that has duckdb and pyarrow:
    ~/HISSS_data/venvs/daymet/bin/python docs/benchmarks/selftest_run_batched_julia.py
Everything happens in a temporary folder; exit 1 if any check fails.
"""
import csv
import datetime
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time

import pyarrow as pa
import pyarrow.parquet as pq

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_batched_julia.py")
FAKE_JULIA = '''#!{python}
"""Fake julia: writes a run's outputs for the gages of STREAMFLOW_DATA_PATH."""
import json, os, subprocess, sys, time
import pyarrow as pa, pyarrow.parquet as pq
env = os.environ
if env.get("FAKE_GRANDCHILD"):
    child = subprocess.Popen(["sleep", "300"])
    with open(env["FAKE_GRANDCHILD"], "w") as fh:
        fh.write(str(child.pid))
time.sleep(float(env.get("FAKE_SLEEP", "0")))
out = env.get("FAKE_OUTDIR_OVERRIDE") or env["STREAMFLOW_OUTPUT_DIR"]
os.makedirs(out, exist_ok=True)
gages = sorted(set(pq.read_table(env["STREAMFLOW_DATA_PATH"], columns=["gage_id"]).column(0).to_pylist()))
with open(os.path.join(out, "fake_signatures.csv"), "w") as fh:
    fh.write("gage_id,a,b\\n")
    for g in gages:
        fh.write(f"{{g}},1.5,{{len(g)}}\\n")
rows = [(g, s, y) for g in gages for s in ("x", "y") for y in (2000, 2001)]
pq.write_table(pa.table({{"gage_id": [r[0] for r in rows], "signature": [r[1] for r in rows],
                         "water_year": [r[2] for r in rows], "value": [1.0] * len(rows)}}),
               os.path.join(out, "fake_signatures_annual.parquet"))
prov = {{"git_revision": "abc", "git_working_tree_dirty": False, "julia_version": "1.12.6", "hostname": "h",
        "config": {{"sha256": "c"}}, "metadata": {{"path": env["STREAMFLOW_METADATA_PATH"]}}, "env_overrides": {{}},
        "streamflow": {{"path": env["STREAMFLOW_DATA_PATH"]}}, "climate": {{"path": env["STREAMFLOW_CLIMATE_PATH"]}}}}
with open(os.path.join(out, "fake_timing.json"), "w") as fh:
    json.dump({{"language": "Julia", "start_time": "2026-10-01T00:00:00", "end_time": "2026-10-01T00:01:00",
               "total_seconds": 1.0, "phases": {{"p": 1.0}}, "n_gages_processed": len(gages),
               "n_annual_rows": len(rows), "n_signature_columns": 2, "n_metadata_columns": 1, "n_qaqc_flags": 0,
               "provenance": prov}}, fh)
'''

results = []


def expect(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def tool(args, env=None, timeout=120):
    return subprocess.run([sys.executable, TOOL] + args, capture_output=True, text=True, env=env, timeout=timeout)


def alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False


def main():
    T = tempfile.mkdtemp(prefix="selftest_batched_")
    try:
        # synthetic inputs: 7 gages in a shuffled row order; climate with an extra site and a 2001 row
        gages = ["20", "01", "12", "03", "11", "02", "10"]
        d0 = datetime.date(2000, 1, 1)
        sf = [(g, d0 + datetime.timedelta(days=i), float(i)) for i in range(3) for g in gages]
        pq.write_table(pa.table({"gage_id": [r[0] for r in sf], "date": [r[1] for r in sf], "Q": [r[2] for r in sf]}),
                       f"{T}/sf.parquet")
        cl = [(g, d, 1.0, 0.0, -1.0) for g in gages + ["99"] for d in (d0, datetime.date(2001, 1, 1))]
        pq.write_table(pa.table({"site_id": [r[0] for r in cl], "Date": [r[1] for r in cl], "prcp": [r[2] for r in cl],
                                 "swe": [r[3] for r in cl], "tmin": [r[4] for r in cl]}), f"{T}/cl.parquet")
        os.makedirs(f"{T}/repo/julia")
        for f in ("meta.csv", "wrapper.jl"):
            open(f"{T}/{f}", "w").write("x\n")
        with open(f"{T}/julia", "w") as fh:
            fh.write(FAKE_JULIA.format(python=sys.executable))
        os.chmod(f"{T}/julia", 0o755)
        base = ["--wrapper", f"{T}/wrapper.jl", "--project", f"{T}/repo/julia", "--inputs", f"{T}/in", "--climate", "c",
                "--metadata", f"{T}/meta.csv", "--julia", f"{T}/julia"]
        env0 = dict(os.environ)

        # ---- split
        r = tool(["split", "--streamflow", f"{T}/sf.parquet", "--batches", "3", "--out-dir", f"{T}/in",
                  "--climate", f"c={T}/cl.parquet,max_dat=2000-12-31"])
        expect("split rejects an unknown climate option", r.returncode != 0 and "unknown option" in r.stderr, r.stderr[-200:])
        r = tool(["split", "--streamflow", f"{T}/sf.parquet", "--batches", "3", "--out-dir", f"{T}/in",
                  "--climate", f"c={T}/cl.parquet,max_date=2000-12-31"])
        expect("split runs", r.returncode == 0, r.stderr[-300:])
        bj = json.load(open(f"{T}/in/batches.json"))
        sha = hashlib.sha256(open(f"{T}/sf.parquet", "rb").read()).hexdigest()
        expect("split records the source sha256 and PAR1 footer",
               bj["sources"]["streamflow"]["sha256"] == sha and bj["sources"]["streamflow"]["par1_footer"] is True
               and bj["sources"]["climate"]["c"]["options"] == {"max_date": "2000-12-31"})
        expect("split assigns sorted gages by rank % K", bj["batch_gages"] == [["01", "10", "20"], ["02", "11"], ["03", "12"]])
        c1 = pq.read_table(f"{T}/in/climate_c_b1.parquet").to_pydict()
        expect("split keeps four climate columns and applies max_date",
               list(c1) == ["site_id", "Date", "prcp", "swe"] and max(c1["Date"]) == d0 and set(c1["site_id"]) == {"01", "10", "20"})

        # ---- run, resume
        r = tool(["run"] + base + ["--out-dir", f"{T}/out"])
        runs = json.load(open(f"{T}/out/batch_runs.json"))
        expect("run completes three batches", r.returncode == 0 and [runs[k]["rc"] for k in "123"] == [0, 0, 0], r.stderr[-300:])
        r = tool(["run"] + base + ["--out-dir", f"{T}/out"])
        expect("a rerun skips the done batches", r.stdout.count("skipped") == 3, r.stdout[-300:])
        shutil.rmtree(f"{T}/out/b2")
        r = tool(["run"] + base + ["--out-dir", f"{T}/out"])
        runs = json.load(open(f"{T}/out/batch_runs.json"))
        expect("a deleted batch folder is run again", "batch 2: rc 0" in r.stdout and r.stdout.count("skipped") == 2
               and len(runs["2"].get("earlier_attempts", [])) == 1, r.stdout[-300:])
        r = tool(["run"] + [x if x != f"{T}/meta.csv" else f"{T}/nope.csv" for x in base] + ["--out-dir", f"{T}/out_x"])
        expect("run refuses a missing input before starting", r.returncode != 0 and "nope.csv" in r.stderr
               and not os.path.exists(f"{T}/out_x/b1"), r.stderr[-200:])
        r = tool(["run"] + base + ["--out-dir", f"{T}/out_wrongdir"], env=dict(env0, FAKE_OUTDIR_OVERRIDE=f"{T}/elsewhere"))
        rj = json.load(open(f"{T}/out_wrongdir/batch_runs.json"))
        expect("a wrapper that ignores STREAMFLOW_OUTPUT_DIR fails the batch", r.returncode != 0 and rj["1"]["rc"] == -998)

        # ---- merge
        r = tool(["merge", "--out-dir", f"{T}/out", "--prefix", "m"])
        ids = [row[0] for row in csv.reader(open(f"{T}/out/m_signatures.csv"))][1:]
        tj = json.load(open(f"{T}/out/m_timing.json"))
        expect("merge without a reference sorts by gage_id", r.returncode == 0 and ids == sorted(gages), r.stderr[-300:])
        expect("merged timing carries the source provenance and the totals",
               tj["provenance"]["streamflow"]["sha256"] == sha and tj["n_gages_processed"] == 7 and tj["n_annual_rows"] == 28
               and set(tj["provenance"]["batch_inputs"]) == {"b1", "b2", "b3"} and os.path.exists(f"{T}/out/m_merge.json"))
        ann = pq.read_table(f"{T}/out/m_signatures_annual.parquet").to_pydict()
        keys = list(zip(ann["gage_id"], ann["signature"], ann["water_year"]))
        expect("merged annual parquet is complete and sorted", len(keys) == 28 and keys == sorted(keys))
        ref = list(reversed(sorted(gages)))
        with open(f"{T}/ref.csv", "w") as fh:
            fh.write("gage_id,z\n" + "".join(f"{g},0\n" for g in ref))
        r = tool(["merge", "--out-dir", f"{T}/out", "--prefix", "m", "--reference", f"{T}/ref.csv"])
        ids = [row[0] for row in csv.reader(open(f"{T}/out/m_signatures.csv"))][1:]
        expect("merge follows the reference order", r.returncode == 0 and ids == ref, r.stderr[-300:])
        with open(f"{T}/ref_short.csv", "w") as fh:
            fh.write("gage_id,z\n" + "".join(f"{g},0\n" for g in ref[1:]))
        r = tool(["merge", "--out-dir", f"{T}/out", "--prefix", "m", "--reference", f"{T}/ref_short.csv"])
        expect("a gage absent from the reference is an error", r.returncode != 0 and "--allow-missing" in r.stderr)
        r = tool(["merge", "--out-dir", f"{T}/out", "--prefix", "m", "--reference", f"{T}/ref_short.csv", "--allow-missing"])
        ids = [row[0] for row in csv.reader(open(f"{T}/out/m_signatures.csv"))][1:]
        expect("--allow-missing appends it", r.returncode == 0 and ids == ref[1:] + [ref[0]], r.stderr[-300:])
        with open(f"{T}/ref_dup.csv", "w") as fh:
            fh.write("gage_id,z\n" + "".join(f"{g},0\n" for g in ref + [ref[0]]))
        r = tool(["merge", "--out-dir", f"{T}/out", "--prefix", "m", "--reference", f"{T}/ref_dup.csv"])
        expect("a duplicated reference gage is an error", r.returncode != 0 and "more than once" in r.stderr)

        def tampered(name, change, needle):
            shutil.copytree(f"{T}/out", f"{T}/t")
            change(f"{T}/t")
            r = tool(["merge", "--out-dir", f"{T}/t", "--prefix", "m"])
            expect(name, r.returncode != 0 and needle in r.stderr, r.stderr[-300:])
            shutil.rmtree(f"{T}/t")

        def edit_json(path, fn):
            j = json.load(open(path))
            fn(j)
            json.dump(j, open(path, "w"))

        tampered("merge refuses a missing batch folder", lambda d: shutil.rmtree(f"{d}/b3"), "b3: folder missing")
        tampered("merge refuses a failed batch", lambda d: edit_json(f"{d}/batch_runs.json", lambda j: j["2"].update(rc=1)),
                 "b2: not done")
        tampered("merge refuses a second annual parquet",
                 lambda d: shutil.copy(f"{d}/b1/fake_signatures_annual.parquet", f"{d}/b1/old_signatures_annual.parquet"),
                 "at most one annual parquet")
        tampered("merge refuses a batch folder beyond K", lambda d: shutil.copytree(f"{d}/b1", f"{d}/b9"), "beyond the 3")
        tampered("merge refuses mismatched annual rows",
                 lambda d: edit_json(f"{d}/b3/fake_timing.json", lambda j: j.update(n_annual_rows=5)), "annual rows")
        tampered("merge refuses mixed code revisions",
                 lambda d: edit_json(f"{d}/b2/fake_timing.json", lambda j: j["provenance"].update(git_revision="zzz")),
                 "git_revision")

        def move_gage(d):
            with open(f"{d}/b2/fake_signatures.csv", "a") as fh:
                fh.write("01,1.5,2\n")
            edit_json(f"{d}/b2/fake_timing.json", lambda j: j.update(n_gages_processed=3))
        tampered("merge refuses a gage in two batches", move_gage, "appears in b1 and b2")

        # ---- the guard: every stop must take the grandchild down too
        def guard(name, extra, env_extra, needle, signal_after=None, need_grandchild=False):
            # an instant stop may come before the fake has spawned its grandchild; the slower ones must not
            # leave it alive
            out, pidf = f"{T}/g_{name.split()[0]}", f"{T}/gc_{name.split()[0]}.pid"
            env = dict(env0, FAKE_SLEEP="60", FAKE_GRANDCHILD=pidf, **env_extra)
            t0 = time.time()
            p = subprocess.Popen([sys.executable, TOOL, "run"] + base + ["--out-dir", out] + extra, env=env,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if signal_after:
                time.sleep(signal_after)
                p.send_signal(signal.SIGTERM)
            _, err = p.communicate(timeout=90)
            time.sleep(0.5)
            gc = int(open(pidf).read()) if os.path.exists(pidf) else None
            rec = json.load(open(f"{out}/batch_runs.json"))["1"] if os.path.exists(f"{out}/batch_runs.json") else {}
            ok = (p.returncode != 0 and rec.get("rc") == -999 and needle in (rec.get("stopped") or "")
                  and (gc is not None or not need_grandchild) and (gc is None or not alive(gc))
                  and time.time() - t0 < 60)
            expect(name, ok, f"rc {p.returncode}, rec {rec.get('rc')} {rec.get('stopped')}, grandchild "
                             f"{gc} alive={gc is not None and alive(gc)}; {err[-200:]}")

        guard("low memory stops the batch", ["--min-avail-gb", "100000"], {}, "available memory")
        guard("max-minutes stops the batch", ["--max-minutes", "0.05"], {}, "max-minutes", need_grandchild=True)
        os.makedirs(f"{T}/fakebin")
        with open(f"{T}/fakebin/vm_stat", "w") as fh:
            fh.write("#!/bin/sh\nexit 1\n")
        os.chmod(f"{T}/fakebin/vm_stat", 0o755)
        guard("a failed memory sample stops the batch", [], {"PATH": f"{T}/fakebin:" + env0["PATH"]}, "sampling failed")
        guard("SIGTERM to the runner stops the batch", [], {}, "interrupted", signal_after=4, need_grandchild=True)
    finally:
        shutil.rmtree(T, ignore_errors=True)
    bad = [n for n, ok in results if not ok]
    print(f"{len(results) - len(bad)}/{len(results)} passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
