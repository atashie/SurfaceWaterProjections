"""Self-test of the Daymet tools on a synthetic file and a local HTTP server (no network).

Builds a Daymet-like 2024 file on the real NA grid. Only two chunk tiles are written; the
rest stays unallocated, so the file is small. It then checks:
  * read_grid: accepts the file; rejects a changed CRS, last cell centre or calendar;
  * daymet_probe.py and daymet_aggregate.py: basin means equal an independent numpy
    computation (fill days and a tile-crossing basin included); the timing file and the
    .done carry the source SHA-256 and the weights md5; an out-of-range value is fatal
    unless --allow-out-of-range;
  * daymet_assemble.py: the written file verifies; --provenance-only works; a tampered
    copy fails verification; a missing .done is refused and the earlier output kept;
  * copy_verify.py: changed content is refused without --replace (the drive copy and its
    manifest line stay as they were) and replaces the line with it; hidden folders are skipped;
  * daymet_stream.fetch against a local server:
      - the bearer token reaches the server but never curl's argv;
      - HTTP 401 and 404 end a source at once;
      - a partial file resumes;
      - a slow connection that still progresses is resumed without counting failures;
      - stop_downloads() terminates a running curl.
Usage: python selftest_daymet_tools.py   (prints PASS/FAIL per check; exit 1 on any FAIL)
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import daymet_common as dc  # noqa: E402

PY = sys.executable
RESULTS = []
CF = {"grid_mapping_name": "lambert_conformal_conic", "longitude_of_central_meridian": -100.0,
      "latitude_of_projection_origin": 42.5, "false_easting": 0.0, "false_northing": 0.0,
      "standard_parallel": np.array([25.0, 60.0]), "semi_major_axis": 6378137.0,
      "inverse_flattening": 298.257223563}
T0 = 27028.5            # 2024-01-01 at noon, days since 1950-01-01
ROWS, COLS = (3000, 3300), (3000, 3600)       # two 300 x 300 tiles of the (10, 300, 300) layout


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{('  -- ' + detail) if detail and not ok else ''}", flush=True)


def values(d0, d1):
    """Synthetic float32 values for days d0..d1-1 over the written block."""
    d = np.arange(d0, d1)[:, None, None]
    r = np.arange(*ROWS)[None, :, None] - ROWS[0]
    c = np.arange(*COLS)[None, None, :] - COLS[0]
    return np.round(1.0 + 0.5 * np.sin(d / 9.0) + 0.001 * r + 0.0003 * c, 2).astype(np.float32)


def make_file(path, cf=CF, t0=T0, x_step=None, data=True, fill_cells=(), poke=None):
    with h5py.File(path, "w") as f:
        x = dc.X0_CENTRE + np.arange(dc.NX) * (x_step or dc.CELL)
        f["x"] = x
        f["y"] = dc.Y0_CENTRE - np.arange(dc.NY) * dc.CELL
        t = f.create_dataset("time", data=t0 + np.arange(365.0))
        t.attrs["units"] = "days since 1950-01-01 00:00:00"
        g = f.create_dataset("lambert_conformal_conic", data=np.int32(0))
        for k, v in cf.items():
            g.attrs[k] = v
        d = f.create_dataset("prcp", shape=(365, dc.NY, dc.NX), dtype="f4", chunks=(10, 300, 300),
                             compression="gzip", shuffle=True, fillvalue=np.float32(-9999.0))
        d.attrs["_FillValue"] = np.float32(-9999.0)
        d.attrs["units"] = "mm/day"
        if data:
            for d0 in range(0, 365, 73):
                block = values(d0, d0 + 73)
                for (r, c, days) in fill_cells:
                    for dd in days:
                        if d0 <= dd < d0 + 73:
                            block[dd - d0, r - ROWS[0], c - COLS[0]] = -9999.0
                if poke and d0 <= poke[2] < d0 + 73:
                    block[poke[2] - d0, poke[0] - ROWS[0], poke[1] - COLS[0]] = poke[3]
                d[d0:d0 + 73, ROWS[0]:ROWS[1], COLS[0]:COLS[1]] = block


BASINS = {   # site_id -> [(row, col, coverage, area_km2)]
    "SYN00001": [(3001, 3001, 1.0, 0.90), (3001, 3002, 1.0, 0.95), (3002, 3001, 1.0, 1.00), (3002, 3002, 1.0, 1.05)],
    "SYN00002": [(3100, 3100, 0.5, 0.50), (3100, 3101, 0.25, 0.24), (3101, 3100, 1.0, 0.98)],
    "SYN00003": [(3005, 3005, 1.0, 1.00), (3005, 3006, 1.0, 1.00)],       # first cell is fill on days 0-9
    "SYN00004": [(3200, 3299, 1.0, 1.00), (3200, 3300, 1.0, 1.10)],       # crosses a tile boundary
}
FILL_CELLS = [(3005, 3005, range(10))]


def make_weights(wdir):
    os.makedirs(wdir, exist_ok=True)
    rows = []
    for i, (s, cells) in enumerate(BASINS.items()):
        for r, c, cov, area in cells:
            rows.append((i, r * dc.NX + c, cov, area))
    w = pd.DataFrame(rows, columns=["basin_idx", "cell_id", "coverage", "area_km2"])
    w["basin_idx"] = w["basin_idx"].astype("int32")
    p = os.path.join(wdir, "daymet_weights.parquet")
    w.to_parquet(p, index=False)
    b = w.groupby("basin_idx").agg(coverage_sum=("coverage", "sum"), area_sum_km2=("area_km2", "sum")).reset_index()
    b.insert(1, "site_id", list(BASINS))
    b.to_csv(os.path.join(wdir, "daymet_basins.csv"), index=False)
    dc.write_json(os.path.join(wdir, "daymet_weights_meta.json"), {"weights_md5": dc.md5sum(p)})
    return w


def expected_means(w):
    vals = values(0, 365)
    for r, c, days in FILL_CELLS:
        for dd in days:
            vals[dd, r - ROWS[0], c - COLS[0]] = -9999.0
    out = {}
    for i, s in enumerate(BASINS):
        sub = w[w["basin_idx"] == i]
        v = np.stack([vals[:, r // dc.NX - ROWS[0], r % dc.NX - COLS[0]] for r in sub["cell_id"]], axis=1)
        wt = sub["area_km2"].to_numpy()
        ok = v != -9999.0
        out[s] = (np.where(ok, v, 0).astype(np.float64) * wt).sum(1) / (ok * wt).sum(1)
    return out


def run(args, **kw):
    return subprocess.run([PY] + args, capture_output=True, text=True, **kw)


def test_files(tmp):
    good = os.path.join(tmp, "daymet_v4_daily_na_prcp_2024.nc")
    make_file(good, fill_cells=FILL_CELLS)
    g = dc.read_grid(good)
    check("read_grid accepts the synthetic 2024 file (Feb 29 present, Dec 31 absent)",
          str(g.dates[0]) == "2024-01-01" and str(g.dates[-1]) == "2024-12-30"
          and np.datetime64("2024-02-29") in g.dates)
    for label, kw, msg in [
        ("a changed standard parallel", {"cf": dict(CF, standard_parallel=np.array([30.0, 60.0]))}, "CRS"),
        ("a sphere instead of WGS84", {"cf": {k: v for k, v in CF.items() if k != "inverse_flattening"}}, "CRS"),
        ("a 0.5 mm/cell drift (passes the spacing check, ends 3.9 m off)", {"x_step": dc.CELL + 0.0005},
         "last cell"),
        ("a calendar starting Jan 2", {"t0": T0 + 1}, "time axis"),
    ]:
        p = os.path.join(tmp, "bad.nc")
        make_file(p, data=False, **kw)
        try:
            dc.read_grid(p)
            check(f"read_grid rejects {label}", False, "no error raised")
        except ValueError as e:
            check(f"read_grid rejects {label}", msg in str(e), str(e)[:200])
    return good


def test_aggregate(tmp, good):
    wdir, out = os.path.join(tmp, "weights"), os.path.join(tmp, "by_var_year")
    w = make_weights(wdir)
    r = run([os.path.join(HERE, "daymet_probe.py"), "--file", good, "--weights-dir", wdir,
             "--out", os.path.join(tmp, "probe.json")])
    check("probe passes (both touched tiles stored)", r.returncode == 0, r.stdout[-500:] + r.stderr[-500:])
    sha = dc.sha256sum(good)
    r = run([os.path.join(HERE, "daymet_aggregate.py"), "--file", good, "--weights-dir", wdir, "--out-dir", out,
             "--workers", "2", "--source-sha256", sha])
    check("aggregate runs", r.returncode == 0, r.stdout[-800:] + r.stderr[-800:])
    t = pd.read_parquet(os.path.join(out, "prcp_2024.parquet"))
    exp = expected_means(w)
    worst = max(np.abs(t[t["site_id"] == s]["prcp"].to_numpy() - exp[s]).max() for s in BASINS)
    check("basin means equal the independent computation (fill days, tile-crossing basin)", worst < 1e-12,
          f"max diff {worst}")
    check("365 dates per basin, Feb 29 present, Dec 31 absent",
          (t.groupby("site_id").size() == 365).all() and pd.Timestamp("2024-02-29") in set(pd.to_datetime(t["Date"]))
          and pd.Timestamp("2024-12-31") not in set(pd.to_datetime(t["Date"])))
    tim = json.load(open(os.path.join(out, "prcp_2024_timing.json")))
    done = json.load(open(os.path.join(out, "prcp_2024.done")))
    wmd5 = json.load(open(os.path.join(wdir, "daymet_weights_meta.json")))["weights_md5"]
    check("timing and .done record the source SHA-256, the weights md5 and the commit",
          tim["source_sha256"] == sha == done["source_sha256"] and tim["weights_md5"] == wmd5 == done["weights_md5"]
          and "commit" in tim["git"])
    # out-of-range: one weighted cell at 2000 mm/day on day 5
    bad = os.path.join(tmp, "oor", "daymet_v4_daily_na_prcp_2024.nc")
    os.makedirs(os.path.dirname(bad))
    make_file(bad, fill_cells=FILL_CELLS, poke=(3001, 3001, 5, 2000.0))
    o2 = os.path.join(tmp, "oor_out")
    r = run([os.path.join(HERE, "daymet_aggregate.py"), "--file", bad, "--weights-dir", wdir, "--out-dir", o2,
             "--workers", "2"])
    check("an out-of-range cell value is fatal and leaves no .done",
          r.returncode != 0 and not os.path.exists(os.path.join(o2, "prcp_2024.done")), r.stdout[-300:])
    r = run([os.path.join(HERE, "daymet_aggregate.py"), "--file", bad, "--weights-dir", wdir, "--out-dir", o2,
             "--workers", "2", "--allow-out-of-range"])
    check("--allow-out-of-range accepts it", r.returncode == 0 and os.path.exists(os.path.join(o2, "prcp_2024.done")))
    return out


def test_assemble(tmp, in_dir):
    sys.path.insert(0, HERE)
    import daymet_assemble as da
    out = os.path.join(tmp, "daymet_2024_test.parquet")
    r = run([os.path.join(HERE, "daymet_assemble.py"), "--in-dir", in_dir, "--years", "2024", "--vars", "prcp",
             "--out", out])
    prov = json.load(open(out + ".provenance.json")) if os.path.exists(out + ".provenance.json") else {}
    check("assemble writes and verifies the file", r.returncode == 0 and prov.get("verification", {}).get("mismatches") == 0,
          r.stdout[-500:] + r.stderr[-500:])
    before = hashlib.md5(open(out, "rb").read()).hexdigest()
    done = os.path.join(in_dir, "prcp_2024.done")
    os.rename(done, done + ".away")
    r = run([os.path.join(HERE, "daymet_assemble.py"), "--in-dir", in_dir, "--years", "2024", "--vars", "prcp",
             "--out", out])
    os.rename(done + ".away", done)
    check("assemble refuses a variable-year without .done and leaves the earlier output as it was",
          r.returncode != 0 and ".done" in r.stderr and hashlib.md5(open(out, "rb").read()).hexdigest() == before
          and not os.path.exists(out + ".tmp"), r.stderr[-300:])
    r = run([os.path.join(HERE, "daymet_assemble.py"), "--in-dir", in_dir, "--years", "2024", "--vars", "prcp",
             "--out", out, "--provenance-only"])
    prov2 = json.load(open(out + ".provenance.json"))
    check("--provenance-only re-verifies and rewrites the sidecar",
          r.returncode == 0 and "sidecar_rewritten_utc" in prov2 and prov2["out_md5"] == prov["out_md5"], r.stderr[-400:])
    t = pq.read_table(out).to_pandas()
    t.loc[17, "prcp"] = np.nextafter(t.loc[17, "prcp"], 10.0)          # one ulp
    bad = os.path.join(tmp, "tampered.parquet")
    pq.write_table(pa.Table.from_pandas(t, preserve_index=False), bad)
    try:
        da.verify_output(bad, in_dir, [2024], ["prcp"])
        check("verification catches a one-ulp change", False, "no error")
    except SystemExit as e:
        check("verification catches a one-ulp change", "VERIFY" in str(e), str(e))


def test_copy(tmp):
    src, dst = os.path.join(tmp, "src"), os.path.join(tmp, "dst")
    os.makedirs(src)
    for n in ("a.txt", "b.txt"):
        open(os.path.join(src, n), "w").write(n)
    cv = os.path.join(HERE, "copy_verify.py")
    os.makedirs(os.path.join(src, "sub", ".hidden"))
    open(os.path.join(src, "sub", "c.txt"), "w").write("c")
    open(os.path.join(src, "sub", ".hidden", "d.txt"), "w").write("d")
    r1 = run([cv, src, dst, "a.txt", "./b.txt", "sub"])
    man = os.path.join(dst, "MANIFEST_md5.txt")
    lines1 = open(man).read().splitlines()
    open(os.path.join(src, "a.txt"), "w").write("")          # a truncated source must not replace a good copy
    r2 = run([cv, src, dst, "a.txt"])
    kept = open(os.path.join(dst, "a.txt")).read() == "a.txt" and open(man).read().splitlines() == lines1
    check("copy_verify refuses changed content without --replace and keeps the drive copy",
          r1.returncode == 0 and r2.returncode != 0 and kept and len(lines1) == 3
          and not os.path.exists(os.path.join(dst, "sub", ".hidden")), "\n".join(lines1))
    open(os.path.join(src, "a.txt"), "w").write("changed")
    r3 = run([cv, "--replace", src, dst, "a.txt"])
    lines = open(man).read().splitlines()
    md5a = hashlib.md5(b"changed").hexdigest()
    check("copy_verify --replace replaces the file and its manifest line",
          r3.returncode == 0 and len(lines) == 3 and lines[0].startswith(md5a)
          and open(os.path.join(dst, "a.txt")).read() == "changed", "\n".join(lines))


class Handler(BaseHTTPRequestHandler):
    files, seen = {}, []

    def log_message(self, *a):
        pass

    def serve(self, head):
        Handler.seen.append((self.command, self.path, self.headers.get("Authorization")))
        kind, _, name = self.path.strip("/").partition("/")
        if kind == "auth" and self.headers.get("Authorization") != "Bearer TESTTOKEN":
            self.send_response(401)
            self.end_headers()
            return
        body = Handler.files.get(name)
        if body is None:
            self.send_response(404)
            self.end_headers()
            return
        rng = self.headers.get("Range")
        start = int(rng.split("=")[1].split("-")[0]) if rng else 0
        part = body[start:]
        self.send_response(206 if rng else 200)
        self.send_header("Content-Length", str(len(part)))
        self.send_header("Accept-Ranges", "bytes")
        if rng:
            self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
        self.end_headers()
        if head:
            return
        step = 4096 if kind == "slow" else len(part) or 1
        try:
            for i in range(0, len(part), step):
                self.wfile.write(part[i:i + step])
                self.wfile.flush()
                if kind == "slow":
                    time.sleep(0.08)               # ~50 KB/s
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_GET(self):
        self.serve(False)

    def do_HEAD(self):
        self.serve(True)


def test_fetch(tmp):
    import daymet_stream as ds
    rng = np.random.default_rng(1)
    Handler.files = {"ok.nc": rng.bytes(3_000_000), "slow.nc": rng.bytes(300_000), "big.nc": rng.bytes(5_000_000)}
    srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}/"
    raw = os.path.join(tmp, "raw")
    os.makedirs(raw)
    tok = os.path.join(tmp, "tok")
    open(tok, "w").write("TESTTOKEN\n")
    argv_seen, real_popen = [], subprocess.Popen

    def spy(args, *a, **kw):
        argv_seen.append(list(args))
        return real_popen(args, *a, **kw)
    ds.subprocess.Popen = spy
    auth = {"token_file": tok, "cookie": os.path.join(raw, ".c"), "ornl_ok": True}
    try:
        ds.ORNL, ds.MIRROR = base + "auth/", base + "nomirror/"
        src, sec, n, got = ds.fetch("ok.nc", 3_000_000, raw, auth, "ornl")
        same = open(os.path.join(raw, "ok.nc"), "rb").read() == Handler.files["ok.nc"]
        sent = [h for (m, p, h) in Handler.seen if p == "/auth/ok.nc" and m == "GET"]
        check("token reaches the server as a bearer header", same and sent and sent[-1] == "Bearer TESTTOKEN")
        check("token never appears in curl's argv", argv_seen and not any("TESTTOKEN" in " ".join(a) for a in argv_seen))
        t = time.time()
        try:
            ds.fetch("ok2.nc", 10, raw, auth, "ornl")
            check("HTTP 404 ends the source at once", False, "no error")
        except RuntimeError as e:
            check("HTTP 404 ends the source at once", "404" in str(e) and time.time() - t < 15, str(e))
        open(tok, "w").write("WRONG\n")
        t = time.time()
        try:
            ds.fetch("ok3.nc", 10, raw, auth, "ornl")
            check("HTTP 401 ends the source at once", False, "no error")
        except RuntimeError as e:
            check("HTTP 401 ends the source at once", "401" in str(e) and time.time() - t < 15, str(e))
        open(tok, "w").write("TESTTOKEN\n")
        ds.MIRROR = base + "data/"
        open(os.path.join(raw, "big.nc.part"), "wb").write(Handler.files["big.nc"][:2_000_000])
        src, sec, n, got = ds.fetch("big.nc", 5_000_000, raw, {"token_file": None, "cookie": "", "ornl_ok": False}, "mirror")
        check("a partial file resumes (only the rest is fetched)",
              open(os.path.join(raw, "big.nc"), "rb").read() == Handler.files["big.nc"] and got == 3_000_000,
              f"got {got}")
        ds.MIRROR = base + "slow/"
        src, sec, n, got = ds.fetch("slow.nc", 300_000, raw, {"token_file": None, "cookie": "", "ornl_ok": False},
                                    "mirror", min_bps=100_000, speed_time=1)
        check("a slow but progressing connection is resumed without counting failures",
              open(os.path.join(raw, "slow.nc"), "rb").read() == Handler.files["slow.nc"] and n > 1,
              f"attempts {n}, {sec:.1f} s")
        Handler.files["stop.nc"] = rng.bytes(2_000_000)
        stopped = []

        def slow_fetch():
            try:
                ds.fetch("stop.nc", 2_000_000, raw, {"token_file": None, "cookie": "", "ornl_ok": False},
                         "mirror", min_bps=1, speed_time=60)
            except RuntimeError as e:
                stopped.append(str(e))
        th = threading.Thread(target=slow_fetch)
        n0 = len(argv_seen)
        th.start()
        time.sleep(1.5)
        ds.stop_downloads()
        th.join(10)
        check("stop_downloads() terminates the running curl and ends the download",
              not th.is_alive() and len(argv_seen) > n0 and stopped and "stopped" in stopped[0], str(stopped))
    finally:
        ds.subprocess.Popen = real_popen
        srv.shutdown()


def main():
    tmp = tempfile.mkdtemp(prefix="daymet_selftest_")
    try:
        good = test_files(tmp)
        in_dir = test_aggregate(tmp, good)
        test_assemble(tmp, in_dir)
        test_copy(tmp)
        test_fetch(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
