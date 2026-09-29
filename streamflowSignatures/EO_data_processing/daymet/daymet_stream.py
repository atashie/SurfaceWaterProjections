"""Stream Daymet NA files through fetch -> verify -> probe -> aggregate -> delete.

Usage:
    python daymet_stream.py --years 1980-2025 --vars prcp,swe,tmin,tmax,vp,srad
        --weights-dir <dir> --out-dir <dir> --raw-dir <dir> [--workers 8] [--keep-raw]
        [--edl-token-file ~/.config/earthdata/edl_token] [--dry-run [--check-mirror]]

Order: the ORNL-only files (calendar 2025; the mirror ends at 2024) come first, while the
Earthdata token is still valid; then variable-major (every year of the first variable,
then the next), so whole variables complete early. Per (variable, year):
  * skip when <out-dir>/<var>_<year>.done exists (fully resumable; a .part is resumed);
  * expected size and SHA-256 come from NASA CMR (cached in <out-dir>/cmr_manifest.json);
  * source: the NCAR GDEX mirror (anonymous) when it holds the file at the CMR size,
    otherwise ORNL over HTTPS with an Earthdata Login: a bearer token read from
    --edl-token-file (sent through a mode-600 header file, never on a command line), or
    else ~/.netrc (machine urs.earthdata.nasa.gov) with cookies;
  * curl range-resume with stall detection; size == CMR, then SHA-256 == CMR (fatal);
  * daymet_probe.py: an unstored chunk under the basin weights is fatal (plan U3);
  * daymet_aggregate.py (area weights by default); the raw file is deleted afterwards
    unless --keep-raw; one line per file goes to <out-dir>/stream_log.csv.
A downloader thread runs one file ahead of the processing, so at most three raw files
(processing, waiting, downloading) are on disk at once (<= ~65 GB with vp). Each start
writes <out-dir>/run_meta_<timestamp>.json (git commit, arguments, software, host).
"""
import argparse
import base64
import csv
import hashlib
import json
import os
import platform
import queue
import subprocess
import sys
import tempfile
import threading
import time

import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import DAYMET_VARS, software_versions, write_json  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
COLLECTION = "C2532426483-ORNL_CLOUD"          # Daymet_Daily_V4R1_2129
MIRROR = "https://tds.gdex.ucar.edu/thredds/fileServer/files/d682806/"
ORNL = "https://data.ornldaac.earthdata.nasa.gov/protected/daymet/Daymet_Daily_V4R1/data/"
MIRROR_LAST_YEAR = 2024


def fname(var, year):
    return f"daymet_v4_daily_na_{var}_{year}.nc"


def cmr_manifest(variables, path):
    man = json.load(open(path)) if os.path.exists(path) else {}
    for v in variables:
        if any(k.startswith(f"daymet_v4_daily_na_{v}_") for k in man):
            continue
        r = requests.get(CMR, params={"collection_concept_id": COLLECTION, "page_size": 200,
                                      "readable_granule_name": f"*_na_{v}_*",
                                      "options[readable_granule_name][pattern]": "true"}, timeout=120)
        r.raise_for_status()
        for it in r.json()["items"]:
            f = it["umm"]["DataGranule"]["ArchiveAndDistributionInformation"][0]
            name = it["umm"]["GranuleUR"].split(".", 1)[1]
            man[name] = {"bytes": int(f["SizeInBytes"]), "sha256": f["Checksum"]["Value"],
                         "algorithm": f["Checksum"]["Algorithm"], "cmr_revision": it["meta"]["revision-date"]}
    write_json(path, man)
    return man


def checksum(path, algorithm, block=1 << 24):
    h = {"SHA-256": hashlib.sha256, "MD5": hashlib.md5}[algorithm]()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def mirror_size(name):
    try:
        r = requests.head(MIRROR + name, timeout=60, allow_redirects=True)
        return int(r.headers.get("Content-Length", -1)) if r.status_code == 200 else -1
    except requests.RequestException:
        return -1


def token_expiry(path):
    """Expiry (epoch seconds) from the token's JWT payload, or None; never prints the token."""
    try:
        payload = open(path).read().strip().split(".")[1]
        return int(json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))["exp"])
    except Exception:          # noqa: BLE001 -- not a JWT: no expiry to check
        return None


def fetch(name, exp_bytes, raw_dir, auth):
    """Download to raw_dir/name (resumable). Returns (source, seconds, attempts)."""
    dest = os.path.join(raw_dir, name)
    if os.path.exists(dest) and os.path.getsize(dest) == exp_bytes:
        return "present", 0.0, 0
    src = "mirror" if mirror_size(name) == exp_bytes else "ornl"
    url = (MIRROR if src == "mirror" else ORNL) + name
    part = dest + ".part"
    hdr = None
    try:
        extra = []
        if src == "ornl":
            if auth["token_file"]:
                fd, hdr = tempfile.mkstemp(prefix=".edl_hdr_", dir=raw_dir)   # mode 600
                with os.fdopen(fd, "w") as fh:
                    fh.write("Authorization: Bearer " + open(auth["token_file"]).read().strip() + "\n")
                extra = ["-H", "@" + hdr]      # curl drops it on the redirect to another host
            else:
                extra = ["-n", "-c", auth["cookie"], "-b", auth["cookie"]]
        t0, n = time.time(), 0
        while (os.path.getsize(part) if os.path.exists(part) else 0) < exp_bytes and n < 60:
            n += 1
            cmd = ["curl", "-sS", "-L", "-f", "-C", "-", "--connect-timeout", "30", "--speed-limit", "200000",
                   "--speed-time", "60"] + extra + ["-o", part, url]
            if subprocess.run(cmd).returncode != 0:
                time.sleep(min(300, 10 * n))
    finally:
        if hdr and os.path.exists(hdr):
            os.remove(hdr)
    got = os.path.getsize(part) if os.path.exists(part) else 0
    if got != exp_bytes:
        raise RuntimeError(f"{name}: {got} of {exp_bytes} bytes after {n} attempts from {src}")
    os.replace(part, dest)
    return src, time.time() - t0, n


def git_state():
    def run(*args):
        r = subprocess.run(["git", "-C", HERE] + list(args), capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    return {"commit": run("rev-parse", "HEAD"),
            "daymet_tools_dirty": bool(run("status", "--porcelain", "--", HERE))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", required=True, help="e.g. 1980-2025")
    ap.add_argument("--vars", default=",".join(DAYMET_VARS))
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--weight", choices=["area", "coverage"], default="area")
    ap.add_argument("--edl-token-file", default=os.path.expanduser("~/.config/earthdata/edl_token"))
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
          f"(~{total / 19e6 / 3600:.1f} h at 19 MB/s); {len(need_ornl)} need ORNL (auth: {auth_desc})")
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
            print(f"mirror vs CMR: {len(same)} same size, {len(absent)} absent from the mirror, "
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
        "started": stamp, "argv": sys.argv, "years": [years[0], years[-1]], "vars": variables,
        "planned": len(planned := [fname(v, yr) for v, yr in plan]), "first_files": planned[:8],
        "bytes_to_download": total, "auth": auth_desc, "token_expiry_utc": exp,
        "git": git_state(), "host": platform.node(), "software": software_versions()})

    auth = {"token_file": tok, "cookie": os.path.join(a.raw_dir, ".edl_cookies")}
    q = queue.Queue(maxsize=1)

    def downloader():
        for v, yr in plan:
            name = fname(v, yr)
            try:
                q.put((v, yr, name) + fetch(name, man[name]["bytes"], a.raw_dir, auth))
            except Exception as e:            # noqa: BLE001 -- surfaced to the main thread
                q.put((v, yr, name, "ERROR", str(e), 0))
                return
        q.put(None)

    threading.Thread(target=downloader, daemon=True).start()
    log = os.path.join(a.out_dir, "stream_log.csv")
    new = not os.path.exists(log)
    with open(log, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["file", "bytes", "source", "download_s", "MB_per_s", "attempts", "sha256_ok",
                        "probe_ok", "aggregate_s", "finished"])
        while (item := q.get()) is not None:
            v, yr, name, src, dl_s, n = item
            if src == "ERROR":
                sys.exit(f"download failed: {dl_s}")
            path = os.path.join(a.raw_dir, name)
            ok_sha = checksum(path, man[name]["algorithm"]) == man[name]["sha256"]
            if not ok_sha:
                sys.exit(f"{man[name]['algorithm']} mismatch for {name} (source {src}); file kept for inspection")
            probe = subprocess.run([sys.executable, os.path.join(HERE, "daymet_probe.py"), "--file", path,
                                    "--weights-dir", a.weights_dir,
                                    "--out", os.path.join(a.out_dir, f"probe_{v}_{yr}.json")],
                                   capture_output=True, text=True)
            if probe.returncode != 0:
                sys.exit(f"probe failed for {name}:\n{probe.stdout[-2000:]}{probe.stderr[-2000:]}")
            t0 = time.time()
            agg = subprocess.run([sys.executable, os.path.join(HERE, "daymet_aggregate.py"), "--file", path,
                                  "--weights-dir", a.weights_dir, "--out-dir", a.out_dir,
                                  "--workers", str(a.workers), "--weight", a.weight],
                                 capture_output=True, text=True)
            if agg.returncode != 0 or not os.path.exists(os.path.join(a.out_dir, f"{v}_{yr}.done")):
                sys.exit(f"aggregate failed for {name}:\n{agg.stdout[-2000:]}{agg.stderr[-2000:]}")
            agg_s = time.time() - t0
            if not a.keep_raw:
                os.remove(path)
            mbps = man[name]["bytes"] / 1e6 / dl_s if dl_s else ""
            w.writerow([name, man[name]["bytes"], src, round(dl_s, 1), round(mbps, 2) if mbps else "", n,
                        ok_sha, True, round(agg_s, 1), time.strftime("%Y-%m-%dT%H:%M:%S")])
            fh.flush()
            print(f"{name}: {src} {dl_s:.0f} s, checksum OK, probe OK, aggregated in {agg_s:.0f} s"
                  f"{'' if a.keep_raw else ', raw deleted'}", flush=True)
    print("STREAM DONE")


if __name__ == "__main__":
    main()
