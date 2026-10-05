"""Stream Daymet NA files through fetch -> verify -> probe -> aggregate -> delete.

Usage:
    python daymet_stream.py --years 1980-2025 --vars prcp,swe,tmin,tmax,vp,srad
        --weights-dir <dir> --out-dir <dir> --raw-dir <dir> [--workers 8] [--keep-raw]
        [--edl-token-file ~/.config/earthdata/edl_token] [--min-speed-mbps 20 --speed-time 30]
        [--dry-run [--check-mirror]]

Order: the ORNL-only files come first (calendar 2025; the mirror ends at 2024), while the
Earthdata token is still valid. After that the order is variable-major (every year of the
first variable, then the next), so whole variables complete early.

Before processing (every start):
  * NASA CMR is re-queried for the expected size and SHA-256 of every granule, cached in
    <out-dir>/cmr_manifest.json with the granule concept id and revision id. A record
    that changed since the cache was written stops the run, and so does a granule listed
    twice.
  * A <var>_<year>.done written with weights other than --weights-dir stops the run.

Per (variable, year):
  * Skip when <out-dir>/<var>_<year>.done exists (fully resumable; a .part is resumed).
  * Source: ORNL over HTTPS with an Earthdata Login when available (--prefer auto; ~3x
    faster than the mirror from the M5 laptop). Otherwise the anonymous NCAR GDEX mirror,
    when it holds the file at the CMR size. The other source is the fallback.
  * ORNL auth: a bearer token from --edl-token-file, handed to curl as a config line on
    stdin (never on a command line or on disk); or else ~/.netrc
    (machine urs.earthdata.nasa.gov) with cookies.
  * curl resumes byte ranges. A connection below --min-speed-mbps for --speed-time
    seconds is dropped and resumed at once, without counting as a failure while it still
    makes progress. HTTP 401/403/404 ends that source at once. Otherwise up to 30
    attempts per source, with back-off. A disk-space guard runs before each download.
  * Size == CMR, then SHA-256 == CMR (fatal; a mismatching file is renamed <name>.bad,
    so the next start downloads it again).
  * daymet_probe.py: an unstored chunk under the basin weights is fatal (plan U3), as is
    any grid, CRS or calendar change (daymet_common.read_grid).
  * daymet_aggregate.py (area weights by default) gets the verified SHA-256 and records
    it. The raw file is deleted afterwards unless --keep-raw. One line per file goes to
    <out-dir>/stream_log.csv.

A downloader thread runs one file ahead of the processing, so at most three raw files
(processing, waiting, downloading) are on disk at once (<= ~65 GB with vp). On any exit
(error, SIGTERM, Ctrl-C) the running curl is terminated and the downloader joined. Each
start writes <out-dir>/run_meta_<timestamp>.json (git commit, arguments, software, host).
"""
import argparse
import base64
import csv
import hashlib
import json
import os
import platform
import queue
import shutil
import signal
import subprocess
import sys
import threading
import time

import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import DAYMET_VARS, git_state, software_versions, utc_now, write_json  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
COLLECTION = "C2532426483-ORNL_CLOUD"          # Daymet_Daily_V4R1_2129
MIRROR = "https://tds.gdex.ucar.edu/thredds/fileServer/files/d682806/"
ORNL = "https://data.ornldaac.earthdata.nasa.gov/protected/daymet/Daymet_Daily_V4R1/data/"
MIRROR_LAST_YEAR = 2024
AUTH_OR_MISSING = {"401", "403", "404"}         # no retry: bad credentials or no such file
DISK_HEADROOM = 5e9                             # bytes kept free beside a download
STOP = threading.Event()                        # set on exit: no new curl starts
_CURL = {"proc": None}                          # the running curl, for stop_downloads()
_LOCK = threading.Lock()


def fname(var, year):
    return f"daymet_v4_daily_na_{var}_{year}.nc"


def cmr_items(v):
    """Every CMR granule of variable v in the Daymet V4 R1 daily collection (all pages)."""
    items, page = [], 1
    while True:
        r = requests.get(CMR, params={"collection_concept_id": COLLECTION, "page_size": 200, "page_num": page,
                                      "readable_granule_name": f"*_na_{v}_*",
                                      "options[readable_granule_name][pattern]": "true"}, timeout=120)
        r.raise_for_status()
        got = r.json()["items"]
        items += got
        hits = int(r.headers.get("CMR-Hits", len(items)))
        if len(items) >= hits or not got:
            return items, hits
        page += 1


def cmr_manifest(variables, path):
    """Expected size and checksum per granule from NASA CMR, re-queried on every start.

    A cached record whose size or checksum has changed stops the run: variable-years
    already processed used the old file. So does a granule listed twice, a granule
    without exactly one file record with a checksum, or fewer granules read than CMR
    reports.
    """
    old = json.load(open(path)) if os.path.exists(path) else {}
    man, drift = dict(old), []
    for v in variables:
        items, hits = cmr_items(v)
        if len(items) != hits:
            sys.exit(f"CMR returned {len(items)} of {hits} {v} granules")
        seen = set()
        for it in items:
            files = it["umm"]["DataGranule"]["ArchiveAndDistributionInformation"]
            if len(files) != 1 or "Checksum" not in files[0]:
                sys.exit(f"{it['umm']['GranuleUR']}: {len(files)} file records (expected one, with a checksum)")
            f = files[0]
            name = it["umm"]["GranuleUR"].split(".", 1)[1]
            if name in seen:
                sys.exit(f"CMR lists {name} twice")
            seen.add(name)
            rec = {"bytes": int(f["SizeInBytes"]), "sha256": f["Checksum"]["Value"],
                   "algorithm": f["Checksum"]["Algorithm"], "cmr_revision": it["meta"]["revision-date"],
                   "concept_id": it["meta"]["concept-id"], "revision_id": it["meta"]["revision-id"]}
            prev = old.get(name)
            if prev and (prev["bytes"], prev["sha256"], prev["algorithm"]) != (
                    rec["bytes"], rec["sha256"], rec["algorithm"]):
                drift.append(name)
            man[name] = rec
    if drift:
        sys.exit(f"CMR changed {len(drift)} granule record(s) since {path} was written ({drift[:5]}); "
                 f"variable-years already processed used the old files -- inspect before continuing")
    write_json(path, man)
    return man


def checksum(path, algorithm, block=1 << 24):
    h = {"SHA-256": hashlib.sha256, "MD5": hashlib.md5}[algorithm]()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def mirror_size(name, tries=3):
    """Size of `name` on the mirror; -1 when absent (HTTP 404) or unreachable after `tries` attempts."""
    for i in range(tries):
        try:
            r = requests.head(MIRROR + name, timeout=60, allow_redirects=True)
            if r.status_code == 200:
                return int(r.headers.get("Content-Length", -1))
            if r.status_code == 404:
                return -1
        except requests.RequestException:
            pass
        time.sleep(5 * (i + 1))
    return -1


def token_expiry(path):
    """Expiry (epoch seconds) from the token's JWT payload, or None; never prints the token."""
    try:
        payload = open(path).read().strip().split(".")[1]
        return int(json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))["exp"])
    except Exception:          # noqa: BLE001 -- not a JWT: no expiry to check
        return None


def run_curl(args, config=None):
    """Run curl; `config` (the auth header) goes on stdin, never into argv or a file.

    Returns (exit code, final HTTP status). The process is registered so that
    stop_downloads() can terminate it from the main thread.
    """
    with _LOCK:
        if STOP.is_set():
            return -1, ""
        p = subprocess.Popen(["curl"] + (["-K", "-"] if config else []) + args + ["-w", "%{http_code}"],
                             stdin=subprocess.PIPE if config else subprocess.DEVNULL,
                             stdout=subprocess.PIPE, text=True)
        _CURL["proc"] = p
    try:
        out, _ = p.communicate(input=config)
    finally:
        with _LOCK:
            _CURL["proc"] = None
    return p.returncode, (out or "").strip()[-3:]


def stop_downloads():
    """No new curl starts; the running one is terminated (its .part stays for a resume)."""
    STOP.set()
    with _LOCK:
        p = _CURL["proc"]
    if p is not None and p.poll() is None:
        p.terminate()
        try:
            p.wait(10)
        except subprocess.TimeoutExpired:
            p.kill()


def fetch(name, exp_bytes, raw_dir, auth, prefer, min_bps=20e6, speed_time=30, mirror_min_bps=10e6):
    """Download to raw_dir/name (resumable). Returns (source(s), seconds, attempts, bytes fetched now).

    Sources are tried in order of preference; a partial file is resumed across sources
    (the two serve identical bytes; SHA-256 is checked afterwards). Per source:
      * HTTP 401/403/404 ends the source at once;
      * a connection below `min_bps` (`mirror_min_bps` on the mirror, which delivers ~15-21 MB/s)
        for `speed_time` s (curl exit 28) that still made progress is resumed at once and does not
        count as a failure;
    The mirror is asked for the file's size only when it is about to be used.
      * any other failure counts, with back-off, up to 30 per source.
    """
    dest = os.path.join(raw_dir, name)
    if os.path.exists(dest) and os.path.getsize(dest) == exp_bytes:
        return "present", 0.0, 0, 0
    part = dest + ".part"
    size = lambda: os.path.getsize(part) if os.path.exists(part) else 0      # noqa: E731
    if size() > exp_bytes:
        os.remove(part)
    free, need = shutil.disk_usage(raw_dir).free, exp_bytes - size() + DISK_HEADROOM
    if free < need:
        raise RuntimeError(f"{name}: {free / 1e9:.1f} GB free in {raw_dir}, {need / 1e9:.1f} GB needed")
    order = ["ornl", "mirror"] if prefer == "ornl" else ["mirror", "ornl"]
    order = [s for s in order if s != "ornl" or auth["ornl_ok"]]
    base = ["-sS", "-L", "-f", "-C", "-", "--connect-timeout", "30", "--speed-time", str(int(speed_time))]
    t0, start, n, used, why = time.time(), size(), 0, [], []
    for src in order:
        if STOP.is_set():
            break
        if src == "mirror" and mirror_size(name) != exp_bytes:
            why.append("mirror lacks the file at the CMR size")
            continue
        url = (MIRROR if src == "mirror" else ORNL) + name
        floor = ["--speed-limit", str(int(mirror_min_bps if src == "mirror" else min_bps))]
        config, extra = None, []
        if src == "ornl":
            if auth["token_file"]:
                # curl drops the header on the redirect to another host (the signed S3 URL)
                config = 'header = "Authorization: Bearer ' + open(auth["token_file"]).read().strip() + '"\n'
            else:
                extra = ["-n", "-c", auth["cookie"], "-b", auth["cookie"]]
        used.append(src)
        fails = 0
        while size() < exp_bytes and fails < 30 and not STOP.is_set():
            n += 1
            before = size()
            rc, http = run_curl(base + floor + extra + ["-o", part, url], config)
            grew = size() > before
            if size() > exp_bytes:                     # never expected: start the file over
                os.remove(part)
                fails += 1
            elif rc == 0 and grew:
                continue
            elif http in AUTH_OR_MISSING:
                why.append(f"{src} HTTP {http}")
                break
            elif rc == 28 and grew:                    # below the speed floor, still progressing
                continue
            else:
                fails += 1
                STOP.wait(min(300, 10 * fails))
        if size() >= exp_bytes:
            break
    got = size()
    if got != exp_bytes:
        raise RuntimeError(f"{name}: {got} of {exp_bytes} bytes after {n} attempts from "
                           f"{'+'.join(used) or 'no usable source'}{'; ' + ', '.join(why) if why else ''}"
                           f"{'; stopped' if STOP.is_set() else ''}")
    os.replace(part, dest)
    return "+".join(used), time.time() - t0, n, got - start


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", required=True, help="e.g. 1980-2025")
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--workers", type=int, default=6,
                    help="aggregation workers; 6 fit the 16 GB laptop with the 8,017-basin weights "
                         "(5.0 GB parent + 1.46 GB per worker, 2026-09-30)")
    ap.add_argument("--weight", choices=["area", "coverage"], default="area")
    ap.add_argument("--edl-token-file", default=os.path.expanduser("~/.config/earthdata/edl_token"))
    ap.add_argument("--prefer", choices=["auto", "ornl", "mirror"], default="auto",
                    help="download source to try first; auto = ORNL when Earthdata auth is available "
                         "(measured 50-54 MB/s from the M5 laptop vs ~18 MB/s from the mirror, 2026-09-29)")
    ap.add_argument("--min-speed-mbps", type=float, default=20.0,
                    help="drop and resume a connection below this rate (MB/s) for --speed-time s; the "
                         "2026-09-29 run's 200 KB/s floor let 26 of 269 files crawl at 3-30 MB/s (6.8 h lost)")
    ap.add_argument("--mirror-min-speed-mbps", type=float, default=10.0,
                    help="the same floor on the NCAR mirror, which delivers ~15-21 MB/s")
    ap.add_argument("--speed-time", type=int, default=30)
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check-mirror", action="store_true",
                    help="with --dry-run: HEAD every planned file on the mirror and compare its size with CMR")
    a = ap.parse_args()
    y = [int(v) for v in a.years.split("-")]
    years = list(range(y[0], y[-1] + 1))
    variables = [v for v in a.vars.split(",") if v]
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.raw_dir, exist_ok=True)
    man = cmr_manifest(variables, os.path.join(a.out_dir, "cmr_manifest.json"))
    wmeta = os.path.join(a.weights_dir, "daymet_weights_meta.json")
    weights_md5 = json.load(open(wmeta)).get("weights_md5") if os.path.exists(wmeta) else None
    stale, legacy = [], 0
    for v in variables:
        for yr in years:
            dp = os.path.join(a.out_dir, f"{v}_{yr}.done")
            if os.path.exists(dp):
                try:
                    d = json.loads(open(dp).read() or "{}")
                except ValueError:
                    sys.exit(f"{dp} is not JSON; delete it (and its outputs) to redo {v} {yr}")
                if not d.get("weights_md5"):
                    legacy += 1
                elif weights_md5 and d["weights_md5"] != weights_md5:
                    stale.append(f"{v}_{yr}")
    if stale:
        sys.exit(f"{len(stale)} .done file(s) were built on other weights than {a.weights_dir} "
                 f"({stale[:5]}); delete them (and their outputs) to redo those variable-years")
    if legacy:
        print(f"note: {legacy} .done file(s) predate the weights record (2026-10-01) and are trusted as they are")
    plan = [(v, yr) for v in variables for yr in years
            if not os.path.exists(os.path.join(a.out_dir, f"{v}_{yr}.done"))]
    plan.sort(key=lambda t: 0 if t[1] > MIRROR_LAST_YEAR else 1)      # ORNL-only files first (stable)
    missing = [fname(v, yr) for v, yr in plan if fname(v, yr) not in man]
    if missing:
        sys.exit(f"not in CMR: {missing[:10]}")
    total = sum(man[fname(v, yr)]["bytes"] for v, yr in plan)
    need_ornl = [fname(v, yr) for v, yr in plan if yr > MIRROR_LAST_YEAR]
    tok = a.edl_token_file if (a.edl_token_file and os.path.exists(a.edl_token_file)
                               and os.path.getsize(a.edl_token_file) > 0) else None
    netrc = os.path.expanduser("~/.netrc")
    has_netrc = os.path.exists(netrc) and "urs.earthdata.nasa.gov" in open(netrc).read()
    auth_desc = "EDL token file" if tok else ("~/.netrc" if has_netrc else "NONE")
    print(f"{len(plan)} variable-years to process, {total / 1e12:.3f} TB to download "
          f"(~{total / 50e6 / 3600:.1f} h at ORNL's ~50 MB/s, ~{total / 19e6 / 3600:.1f} h at the mirror's ~19 MB/s); "
          f"{len(need_ornl)} exist only at ORNL (auth: {auth_desc})")
    exp = token_expiry(tok) if tok else None
    if exp:
        print(f"EDL token expires {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(exp))}"
              f" ({(exp - time.time()) / 86400:.1f} days from now)")
        if exp < time.time():
            sys.exit("the EDL token has expired; put a fresh one in the token file")
    if a.dry_run:
        if a.check_mirror:
            rows = []
            for v, yr in plan:
                name = fname(v, yr)
                ms = mirror_size(name)
                rows.append((name, ms, man[name]["bytes"]))
            same = [r for r in rows if r[1] == r[2]]
            absent = [r[0] for r in rows if r[1] < 0]
            differ = [r for r in rows if r[1] >= 0 and r[1] != r[2]]
            print(f"mirror vs CMR: {len(same)} same size, {len(absent)} absent from (or unreachable on) the mirror, "
                  f"{len(differ)} differ")
            for r in differ:
                print(f"  DIFFERS {r[0]}: mirror {r[1]} B, CMR {r[2]} B")
            if absent:
                print(f"  absent: {absent}")
            write_json(os.path.join(a.out_dir, "mirror_check.json"),
                       {"same": [r[0] for r in same], "absent": absent, "differ": differ})
        return
    if need_ornl and not (tok or has_netrc):
        sys.exit("ORNL files need Earthdata auth: a token in --edl-token-file, or ~/.netrc "
                 "(machine urs.earthdata.nasa.gov login <user> password <pw>, chmod 600)")
    stamp = time.strftime("%Y%m%dT%H%M%S")
    write_json(os.path.join(a.out_dir, f"run_meta_{stamp}.json"), {
        "started": stamp, "started_utc": utc_now(), "argv": sys.argv, "years": [years[0], years[-1]],
        "vars": variables, "planned": len(planned := [fname(v, yr) for v, yr in plan]), "first_files": planned[:8],
        "bytes_to_download": total, "auth": auth_desc, "token_expiry_utc": exp, "prefer": a.prefer,
        "min_speed_mbps": a.min_speed_mbps, "mirror_min_speed_mbps": a.mirror_min_speed_mbps,
        "speed_time": a.speed_time, "weights_md5": weights_md5,
        "git": git_state(), "host": platform.node(), "software": software_versions()})

    auth = {"token_file": tok, "cookie": os.path.join(a.raw_dir, ".edl_cookies"), "ornl_ok": bool(tok or has_netrc)}
    prefer = ("ornl" if auth["ornl_ok"] else "mirror") if a.prefer == "auto" else a.prefer
    q = queue.Queue(maxsize=1)

    def put(item):
        while not STOP.is_set():
            try:
                q.put(item, timeout=1)
                return True
            except queue.Full:
                pass
        return False

    def downloader():
        for v, yr in plan:
            if STOP.is_set():
                return
            name = fname(v, yr)
            try:
                res = fetch(name, man[name]["bytes"], a.raw_dir, auth, prefer,
                            a.min_speed_mbps * 1e6, a.speed_time, a.mirror_min_speed_mbps * 1e6)
            except Exception as e:            # noqa: BLE001 -- surfaced to the main thread
                put((v, yr, name, "ERROR", str(e), 0, 0))
                return
            if not put((v, yr, name) + res):
                return
        put(None)

    def next_item():
        while True:
            try:
                return q.get(timeout=5)
            except queue.Empty:
                if not th.is_alive() and q.empty():
                    sys.exit("the downloader stopped without reporting a result")

    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))     # run the cleanup below on kill
    th = threading.Thread(target=downloader, name="downloader")   # not a daemon: joined on every exit
    th.start()
    log = os.path.join(a.out_dir, "stream_log.csv")
    new = not os.path.exists(log)
    try:
        with open(log, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["file", "bytes", "source", "download_s", "MB_per_s", "attempts", "sha256_ok",
                            "probe_ok", "aggregate_s", "finished"])
            while (item := next_item()) is not None:
                process(item, a, man, w, fh)
    finally:
        stop_downloads()
        th.join(timeout=60)
    print("STREAM DONE")


def process(item, a, man, w, fh):
    """Verify, probe and aggregate one downloaded file (main thread); exits on any failure."""
    v, yr, name, src, dl_s, n, nbytes = item
    if src == "ERROR":
        sys.exit(f"download failed: {dl_s}")
    path = os.path.join(a.raw_dir, name)
    digest = checksum(path, man[name]["algorithm"])
    ok_sha = digest == man[name]["sha256"]
    if not ok_sha:
        os.replace(path, path + ".bad")
        sys.exit(f"{man[name]['algorithm']} mismatch for {name} (source {src}); kept as {name}.bad, "
                 f"the next start downloads it again")
    probe = subprocess.run([sys.executable, os.path.join(HERE, "daymet_probe.py"), "--file", path,
                            "--weights-dir", a.weights_dir,
                            "--out", os.path.join(a.out_dir, f"probe_{v}_{yr}.json")],
                           capture_output=True, text=True)
    if probe.returncode != 0:
        sys.exit(f"probe failed for {name}:\n{probe.stdout[-2000:]}{probe.stderr[-2000:]}")
    t0 = time.time()
    agg = subprocess.run([sys.executable, os.path.join(HERE, "daymet_aggregate.py"), "--file", path,
                          "--weights-dir", a.weights_dir, "--out-dir", a.out_dir,
                          "--workers", str(a.workers), "--weight", a.weight,
                          "--source-sha256", digest],
                         capture_output=True, text=True)
    if agg.returncode != 0 or not os.path.exists(os.path.join(a.out_dir, f"{v}_{yr}.done")):
        sys.exit(f"aggregate failed for {name}:\n{agg.stdout[-2000:]}{agg.stderr[-2000:]}")
    agg_s = time.time() - t0
    if not a.keep_raw:
        os.remove(path)
    mbps = nbytes / 1e6 / dl_s if dl_s else ""        # bytes fetched in this start only
    w.writerow([name, man[name]["bytes"], src, round(dl_s, 1), round(mbps, 2) if mbps else "", n,
                ok_sha, True, round(agg_s, 1), time.strftime("%Y-%m-%dT%H:%M:%S")])
    fh.flush()
    print(f"{name}: {src} {dl_s:.0f} s, checksum OK, probe OK, aggregated in {agg_s:.0f} s"
          f"{'' if a.keep_raw else ', raw deleted'}", flush=True)


if __name__ == "__main__":
    main()
