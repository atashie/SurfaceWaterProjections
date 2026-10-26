#!/usr/bin/env python3
"""splice_nwis_qualifiers.py — patch the February 2026 streamflow parquet for the USGS
qualifier-mask defect (docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md, Phase B;
decisions 1–2, 4 and 7 of 2026-10-26; census 2026-09-04 §10.1).

The February compile kept a USGS daily value only when its qualifier string was exactly
"A", "A e", "P" or "P e"; every other code set Q to NA while the code stayed in `flag`.
This tool re-fetches the affected gages from the NWIS daily-value service, restores the
APPROVED masked days (A R, A <, A >, A [4] — decision 1a–1c) from the October 2026
retrieval into a patched copy of the parquet, and writes the ledgers and audits the plan
asks for. Nothing else in the parquet changes. Every artefact goes under --out-dir
(convention: ~/HISSS_data/qualifier-rerun-<date>/input/).

Subcommands (run them in this order):
  list     --parquet P --metadata M --out-dir D
           Scan the parquet for US rows with Q NA and a qualifier code: the audit-fetch
           universe (every masked gage) and the splice universe (gages holding a splice
           code). Writes masked_rows.parquet, masked_days_by_code.csv,
           masked_days_by_gage_code.csv, audit_universe.csv, splice_universe.csv, list.json.
  fetch    --universe audit_universe.csv --cache-dir C [--start 1979-10-01 --end 2025-09-30]
           [--sleep 1.0] [--retries 4]
           Legacy NWIS DV JSON per gage (parameter 00060, statistic 00003, siteStatus=all),
           gzip-cached verbatim, resumable (a cached response is reused only for the identical
           request URL, window included), retried, >= --sleep s apart. Writes
           fetch_report.csv (one row per gage and series) next to the cache.
  splice   --parquet P --metadata M --masked masked_rows.parquet --cache-dir C --out-dir D
           Q_new = value_today * conversion on the masked splice-code days that carry a
           finite, non-sentinel value today; everything else untouched; identical schema and
           row order. Writes the patched parquet, correction_ledger.csv, provenance.json and
           the metadata CSV with the spliced gages' rows updated in place (decision 7).
  control  --parquet P --metadata M --cache-dir C --universe splice_universe.csv --out-dir D
           Old rule vs new rule on the October fetch, and the October values against the
           February non-NA values (B2 conversion identity / source-difference report).
  audit    --parquet P --patched P2 --metadata M --metadata-patched M2 --ledger L
           --masked masked_rows.parquet --cache-dir C --out-dir D
           Gates B1–B6 and B8 -> input_audit.md + CSVs (B7 is qualifier_year_oracle.py).

Needs pyarrow and pandas only (the repo .venv): `.venv/bin/python -I docs/benchmarks/splice_nwis_qualifiers.py ...`.
The acceptance rule it applies is the one of R/usgs_qualifiers.R (status letter A or P, finite
value, not the no-data sentinel); keep the two in step.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

KEPT_CODES_FEB = ("A", "A e", "P", "P e")          # the February rule
SPLICE_CODES = ("A R", "A <", "A >", "A [4]")        # decision 1a–1c
NO_DATA_VALUE = -999999.0                            # NWIS DV `noDataValue`
NWIS_DV_URL = "https://waterservices.usgs.gov/nwis/dv/"
USER_AGENT = "HISSS-qualifier-fix/1.0 (streamflow signatures; contact via github.com/CZ-Sync/HISSS)"
CFS_TO_MM_NUM = 60 * 60 * 24 * 1e6                  # conversion = 86400 / (km2 * 3280.84^3) * 1e6 (R ingestion)
FT_PER_KM_CUBED = 3280.84 ** 3


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def normalize_id(gid):
    s = str(gid).strip()
    if s.isdigit():
        return s.lstrip("0") or "0"
    return s


def status_letter(flag):
    """First token of a NWIS qualification string ('' when absent) — R/usgs_qualifiers.R."""
    if flag is None or (isinstance(flag, float) and math.isnan(flag)):
        return ""
    s = str(flag).strip()
    for sep in (" ", "\t", ":", ","):
        s = s.split(sep)[0]
    return s


def categorize(flag):
    """splice | approved_other | provisional_remark | kept | other."""
    if flag in SPLICE_CODES:
        return "splice"
    if flag in KEPT_CODES_FEB:
        return "kept"
    st = status_letter(flag)
    if st == "A":
        return "approved_other"
    if st == "P":
        return "provisional_remark"
    return "other"


def conversion_for(basin_area_km2):
    """cfs -> mm/day, exactly as R/helperFunctions.R::generate_streamflow_dt."""
    return 60 * 60 * 24 / (basin_area_km2 * FT_PER_KM_CUBED) * 1e6


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def par1_footer(path):
    with open(path, "rb") as f:
        f.seek(-4, os.SEEK_END)
        return f.read(4) == b"PAR1"


def file_info(path):
    return {"path": os.path.abspath(path), "bytes": os.path.getsize(path), "sha256": sha256_file(path),
            "par1_footer": par1_footer(path),
            "mtime": dt.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%dT%H:%M:%S")}


def read_metadata(path):
    m = pd.read_csv(path, dtype=str, keep_default_na=False)
    m["norm_id"] = m["gage_id"].map(normalize_id)
    if m["norm_id"].duplicated().any():
        dups = m.loc[m["norm_id"].duplicated(keep=False), "gage_id"].tolist()
        raise SystemExit(f"metadata: normalized ids are not unique: {dups[:10]}")
    return m


# --------------------------------------------------------------------------------------
# list
# --------------------------------------------------------------------------------------
def cmd_list(a):
    os.makedirs(a.out_dir, exist_ok=True)
    meta = read_metadata(a.metadata)
    us_ids = set(meta.loc[meta["gage_type"] == "USGS", "norm_id"])
    log(f"metadata: {len(meta)} rows, {len(us_ids)} USGS")
    f = pq.ParquetFile(a.parquet)
    parts = []
    n_rows = 0
    n_us_rows = 0
    n_nocode_novalue = 0
    gage_rows = {}
    t0 = time.time()
    for i in range(f.num_row_groups):
        rg = f.read_row_group(i, columns=["gage_id", "Date", "Q", "flag"])
        n_rows += rg.num_rows
        gid = rg.column("gage_id")
        is_us = pa.array([normalize_id(g) in us_ids for g in gid.to_pylist()])
        n_us_rows += pc.sum(is_us).as_py()
        qnull = pc.is_null(rg.column("Q"))
        fnull = pc.is_null(rg.column("flag"))
        n_nocode_novalue += pc.sum(pc.and_(pc.and_(qnull, fnull), is_us)).as_py()
        masked = rg.filter(pc.and_(pc.and_(qnull, pc.invert(fnull)), is_us))
        if masked.num_rows:
            parts.append(masked)
        for g, n in zip(*np.unique(np.array(gid.to_pylist()), return_counts=True)):
            gage_rows[g] = gage_rows.get(g, 0) + int(n)
        if (i + 1) % 100 == 0:
            log(f"  row group {i + 1}/{f.num_row_groups} ({n_rows:,} rows, {time.time() - t0:.0f}s)")
    masked = pa.concat_tables(parts).to_pandas()
    masked = masked.rename(columns={"Date": "date"})
    masked["date"] = pd.to_datetime(masked["date"]).dt.date
    masked["category"] = masked["flag"].map(categorize)
    masked = masked.sort_values(["gage_id", "date"]).reset_index(drop=True)
    kept_but_null = int((masked["category"] == "kept").sum())   # kept codes with a NULL value (sentinel / absent)
    masked_rows = masked[masked["category"] != "kept"].copy()
    log(f"scanned {n_rows:,} rows ({n_us_rows:,} US) in {time.time() - t0:.0f}s: {len(masked_rows):,} masked US rows "
        f"(Q NA with a non-kept code), {kept_but_null:,} kept-code rows with Q NA, {n_nocode_novalue:,} US rows with no code and no value")
    pq.write_table(pa.Table.from_pandas(masked_rows[["gage_id", "date", "flag", "category"]], preserve_index=False),
                   os.path.join(a.out_dir, "masked_rows.parquet"), compression="zstd")
    by_code = (masked_rows.groupby(["flag", "category"]).agg(n_days=("date", "size"), n_gages=("gage_id", "nunique"))
               .reset_index().sort_values("n_days", ascending=False))
    by_code.to_csv(os.path.join(a.out_dir, "masked_days_by_code.csv"), index=False)
    by_gage_code = (masked_rows.groupby(["gage_id", "flag", "category"]).agg(n_days=("date", "size"), first_date=("date", "min"),
                    last_date=("date", "max")).reset_index())
    by_gage_code.to_csv(os.path.join(a.out_dir, "masked_days_by_gage_code.csv"), index=False)
    piv = masked_rows.pivot_table(index="gage_id", columns="category", values="date", aggfunc="size", fill_value=0)
    for c in ("splice", "approved_other", "provisional_remark", "other"):
        if c not in piv.columns:
            piv[c] = 0
    uni = piv[["splice", "approved_other", "provisional_remark", "other"]].astype(int)
    uni["n_masked_days"] = uni.sum(axis=1)
    uni["first_masked"] = masked_rows.groupby("gage_id")["date"].min()
    uni["last_masked"] = masked_rows.groupby("gage_id")["date"].max()
    uni["n_rows_in_parquet"] = [gage_rows.get(g, 0) for g in uni.index]
    uni = uni.reset_index().rename(columns={"splice": "n_splice_days", "approved_other": "n_approved_other_days",
                                            "provisional_remark": "n_provisional_days", "other": "n_other_days"})
    uni.to_csv(os.path.join(a.out_dir, "audit_universe.csv"), index=False)
    spl = uni[uni["n_splice_days"] > 0]
    spl.to_csv(os.path.join(a.out_dir, "splice_universe.csv"), index=False)
    summary = {
        "parquet": file_info(a.parquet), "metadata": file_info(a.metadata),
        "rows": n_rows, "us_rows": n_us_rows, "us_rows_no_code_no_value": n_nocode_novalue,
        "masked_rows": int(len(masked_rows)), "kept_code_rows_with_null_q": kept_but_null,
        "audit_fetch_universe_gages": int(len(uni)), "splice_universe_gages": int(len(spl)),
        "days_by_category": masked_rows["category"].value_counts().to_dict(),
        "splice_days_by_code": masked_rows[masked_rows["category"] == "splice"]["flag"].value_counts().to_dict(),
        "approved_other_codes": masked_rows[masked_rows["category"] == "approved_other"]["flag"].value_counts().to_dict(),
        "other_codes": masked_rows[masked_rows["category"] == "other"]["flag"].value_counts().to_dict(),
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
    }
    json.dump(summary, open(os.path.join(a.out_dir, "list.json"), "w"), indent=1, default=str)
    log(f"audit-fetch universe {len(uni)} gages; splice universe {len(spl)} gages; "
        f"splice days {summary['days_by_category'].get('splice', 0):,}; provisional-remark days "
        f"{summary['days_by_category'].get('provisional_remark', 0):,}; approved-other codes {summary['approved_other_codes']}; "
        f"other codes {summary['other_codes']}")


# --------------------------------------------------------------------------------------
# fetch
# --------------------------------------------------------------------------------------
def nwis_url(site, start, end):
    q = {"format": "json", "sites": site, "parameterCd": "00060", "statCd": "00003",
         "startDT": start, "endDT": end, "siteStatus": "all"}
    return NWIS_DV_URL + "?" + urllib.parse.urlencode(q)


def http_get(url, timeout=180):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, r.read()


def parse_series(raw_json):
    """One record per time series of a cached DV response."""
    j = json.loads(raw_json)
    out = []
    for s in j.get("value", {}).get("timeSeries", []):
        var = s["variable"]
        stat = var.get("options", {}).get("option", [{}])[0].get("optionCode")
        nodata = float(var.get("noDataValue", NO_DATA_VALUE))
        for blk in s.get("values", []):
            vals = blk.get("value", [])
            method = (blk.get("method") or [{}])[0]
            codes = {}
            n_nodata = 0
            dates = []
            for v in vals:
                q = " ".join(v.get("qualifiers", []))
                codes[q] = codes.get(q, 0) + 1
                try:
                    if float(v["value"]) == nodata:
                        n_nodata += 1
                except (TypeError, ValueError):
                    n_nodata += 1
                dates.append(v["dateTime"][:10])
            out.append({
                "site": s["sourceInfo"]["siteCode"][0]["value"],
                "variable": var["variableCode"][0]["value"], "statistic": stat,
                "method_id": method.get("methodID"), "method": method.get("methodDescription"),
                "unit": var.get("unit", {}).get("unitCode"), "no_data_value": nodata,
                "n_values": len(vals), "n_nodata": n_nodata,
                "first_date": min(dates) if dates else None, "last_date": max(dates) if dates else None,
                "n_unique_dates": len(set(dates)),
                "codes_seen": ";".join(f"{k}={v}" for k, v in sorted(codes.items(), key=lambda kv: -kv[1])),
            })
    return out


def cache_path(cache_dir, site):
    return os.path.join(cache_dir, f"{site}.json.gz")


def cmd_fetch(a):
    os.makedirs(a.cache_dir, exist_ok=True)
    uni = pd.read_csv(a.universe, dtype=str)
    sites = uni["gage_id"].astype(str).tolist()
    log(f"fetch: {len(sites)} sites, window {a.start}..{a.end}, cache {a.cache_dir}")
    meta_rows = {}
    meta_path = os.path.join(a.cache_dir, "fetch_meta.json")
    if os.path.exists(meta_path):
        meta_rows = json.load(open(meta_path))
    last_req = 0.0
    n_new = n_cached = n_fail = 0
    for i, site in enumerate(sites):
        cp = cache_path(a.cache_dir, site)
        url = nwis_url(site, a.start, a.end)
        if os.path.exists(cp) and site in meta_rows and meta_rows[site].get("ok") and meta_rows[site].get("url") == url:
            n_cached += 1          # the same request (window included) is cached; a different window refetches
            continue
        ok = False
        err = None
        for attempt in range(a.retries):
            wait = max(0.0, a.sleep - (time.time() - last_req))
            if wait > 0:
                time.sleep(wait)
            t = time.time()
            try:
                last_req = time.time()
                status, raw = http_get(url)
                json.loads(raw)   # must parse
                with gzip.open(cp, "wb") as f:
                    f.write(raw)
                meta_rows[site] = {"ok": True, "http_status": status, "bytes": len(raw), "url": url,
                                   "retrieved": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                                   "seconds": round(time.time() - t, 2), "attempts": attempt + 1}
                ok = True
                break
            except urllib.error.HTTPError as e:
                err = f"HTTP {e.code}"
                if e.code == 404 or e.code == 400:
                    break
            except Exception as e:  # noqa: BLE001
                err = repr(e)[:200]
            time.sleep(2 ** attempt)
        if ok:
            n_new += 1
        else:
            n_fail += 1
            meta_rows[site] = {"ok": False, "error": err, "url": url,
                               "retrieved": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
        if (i + 1) % 25 == 0 or not ok:
            json.dump(meta_rows, open(meta_path, "w"), indent=0)
            log(f"  {i + 1}/{len(sites)}: new {n_new}, cached {n_cached}, failed {n_fail}" + (f" (last error {err})" if not ok else ""))
    json.dump(meta_rows, open(meta_path, "w"), indent=0)
    # the report is rebuilt from every cached response so a resumed run reports everything
    rows = []
    start_d, end_d = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    n_window_days = (end_d - start_d).days + 1
    for site in sites:
        m = meta_rows.get(site, {})
        cp = cache_path(a.cache_dir, site)
        if not (m.get("ok") and os.path.exists(cp)):
            rows.append({"gage_id": site, "ok": False, "error": m.get("error"), "n_series": 0})
            continue
        series = parse_series(gzip.open(cp, "rb").read())
        primary = [s for s in series if s["variable"] == "00060" and s["statistic"] == "00003"]
        for s in series:
            rows.append({"gage_id": site, "ok": True, "http_status": m.get("http_status"), "bytes": m.get("bytes"),
                         "retrieved": m.get("retrieved"), "n_series": len(series), "n_primary_series": len(primary),
                         **{k: s[k] for k in ("variable", "statistic", "method_id", "method", "unit", "no_data_value",
                                               "n_values", "n_nodata", "first_date", "last_date", "n_unique_dates", "codes_seen")},
                         "n_absent_in_window": n_window_days - s["n_unique_dates"] if s["n_values"] else n_window_days})
        if not series:
            rows.append({"gage_id": site, "ok": True, "http_status": m.get("http_status"), "bytes": m.get("bytes"),
                         "retrieved": m.get("retrieved"), "n_series": 0, "n_primary_series": 0})
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(a.cache_dir, "fetch_report.csv"), index=False)
    n_ok = int(rep.drop_duplicates("gage_id")["ok"].sum()) if len(rep) else 0
    log(f"fetch done: {n_ok}/{len(sites)} sites cached ({n_new} new this run, {n_fail} failed); report {os.path.join(a.cache_dir, 'fetch_report.csv')}")


# --------------------------------------------------------------------------------------
# the acceptance rule (mirrors R/usgs_qualifiers.R) and today's record per site
# --------------------------------------------------------------------------------------
def accept(value, qualifier, no_data=NO_DATA_VALUE):
    """(accepted, parsed value or None, reason): status letter A or P, finite, not the sentinel."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return False, None, "nonnumeric"
    if not math.isfinite(v):
        return False, None, "nonfinite"
    if v == no_data:
        return False, None, "sentinel"
    if qualifier is None or status_letter(qualifier) not in ("A", "P"):
        return False, v, "status"
    return True, v, "ok"


def load_today(cache_dir, site):
    """{date: (raw value string, qualifier string)} of the unqualified method block of the one
    primary 00060 / 00003 series of a cached response, plus a checks dict. Returns (None, info)
    when the record cannot be used: not cached, no or several primary series, several unqualified
    method blocks, a duplicate date, a site code other than the one requested, a unit other than
    ft3/s. Values are parsed by accept(), never here."""
    cp = cache_path(cache_dir, site)
    if not os.path.exists(cp):
        return None, {"error": "not_cached"}
    j = json.loads(gzip.open(cp, "rb").read())
    prim = []
    for s in j.get("value", {}).get("timeSeries", []):
        var = s["variable"]
        stat = var.get("options", {}).get("option", [{}])[0].get("optionCode")
        if var["variableCode"][0]["value"] == "00060" and stat == "00003":
            prim.append(s)
    if len(prim) != 1:
        return None, {"error": f"{len(prim)}_primary_series", "n_primary": len(prim)}
    s = prim[0]
    served_site = s["sourceInfo"]["siteCode"][0]["value"]
    if normalize_id(served_site) != normalize_id(site):
        return None, {"error": f"site_mismatch:{served_site}"}
    unit = s["variable"].get("unit", {}).get("unitCode")
    if unit != "ft3/s":
        return None, {"error": f"unit:{unit}"}
    nodata = float(s["variable"].get("noDataValue", NO_DATA_VALUE))
    blocks = s.get("values", [])
    if len(blocks) > 1:
        unq = [b for b in blocks if not (b.get("method") or [{}])[0].get("methodDescription")]
        if len(unq) != 1:
            return None, {"error": f"{len(blocks)}_method_blocks_{len(unq)}_unqualified", "n_blocks": len(blocks)}
        blocks = unq
    rec = {}
    for blk in blocks:
        for v in blk.get("value", []):
            d = v["dateTime"][:10]
            if d in rec:
                return None, {"error": f"duplicate_date:{d}"}
            rec[d] = (v.get("value"), " ".join(v.get("qualifiers", [])))
    return rec, {"unit": unit, "no_data_value": nodata, "n_dates": len(rec), "n_blocks": len(s.get("values", [])),
                 "method_id": (blocks[0].get("method") or [{}])[0].get("methodID") if blocks else None,
                 "served_site": served_site, "first_date": min(rec) if rec else None, "last_date": max(rec) if rec else None}


# --------------------------------------------------------------------------------------
# shared: per-gage access to the February parquet (streamed; gages contiguous)
# --------------------------------------------------------------------------------------
def iter_parquet_gages(parquet, columns=("gage_id", "Date", "Q", "flag"), only=None):
    """Yield (gage_id, pandas frame of its rows in file order); `only` restricts to a set."""
    f = pq.ParquetFile(parquet)
    cur, parts, seen = None, [], set()
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=list(columns)).to_pandas()
        gid = t["gage_id"].to_numpy(dtype=object)
        change = np.flatnonzero(gid[1:] != gid[:-1]) + 1
        for s, e in zip(np.concatenate(([0], change)), np.concatenate((change, [len(gid)]))):
            g = gid[s]
            if g != cur:
                if cur is not None and parts and (only is None or cur in only):
                    yield cur, pd.concat(parts, ignore_index=True)
                if g in seen:
                    raise SystemExit(f"gage {g} is not contiguous in the parquet")
                if cur is not None:
                    seen.add(cur)
                cur, parts = g, []
            if only is None or g in only:
                parts.append(t.iloc[s:e])
    if cur is not None and parts and (only is None or cur in only):
        yield cur, pd.concat(parts, ignore_index=True)


def compile_year_stats(dates, q, min_q=0.0001, min_days=30):
    """The February compile's metadata rule (R/helperFunctions.R process_gages_*): calendar
    years with >= min_days non-NA days > min_q; num_days = rows in those years."""
    years = pd.to_datetime(dates).year.to_numpy()
    qv = np.asarray(q, dtype=float)
    good = []
    for y in np.unique(years):
        sel = years == y
        if int(np.sum(qv[sel] > min_q)) >= min_days:      # NaN > x is False
            good.append(int(y))
    if not good:
        return {"num_years": 0, "start_year": None, "end_year": None, "num_days": 0}
    return {"num_years": len(good), "start_year": min(good), "end_year": max(good),
            "num_days": int(np.isin(years, good).sum())}


def read_ledger(path):
    """The correction ledger with its doubles round-tripped exactly."""
    led = pd.read_csv(path, dtype={"gage_id": str, "date": str}, float_precision="round_trip")
    led["date"] = pd.to_datetime(led["date"]).dt.date
    return led


def basin_areas(meta):
    out = {}
    for nid, ba in zip(meta["norm_id"], meta["basin_area"]):
        try:
            out[nid] = float(ba)
        except (TypeError, ValueError):
            out[nid] = float("nan")
    return out


# --------------------------------------------------------------------------------------
# splice
# --------------------------------------------------------------------------------------
def cmd_splice(a):
    os.makedirs(a.out_dir, exist_ok=True)
    meta = read_metadata(a.metadata)
    area = basin_areas(meta)
    masked = pq.read_table(a.masked).to_pandas()
    masked["date"] = pd.to_datetime(masked["date"]).dt.date
    masked["category"] = masked["flag"].map(categorize)           # never trust the stored category
    targets = masked[masked["category"] == "splice"].copy()
    if a.universe:
        uni = set(pd.read_csv(a.universe, dtype=str)["gage_id"].astype(str))
        targets = targets[targets["gage_id"].isin(uni)]
    gages = sorted(targets["gage_id"].unique())
    log(f"splice: {len(targets):,} target days on {len(gages)} gages; codes {targets['flag'].value_counts().to_dict()}")
    fetch_meta = json.load(open(os.path.join(a.cache_dir, "fetch_meta.json")))
    led_rows, series_info, unusable = [], {}, {}
    for g in gages:
        rec, info = load_today(a.cache_dir, g)
        series_info[g] = info
        if rec is None:
            unusable[g] = info.get("error")
        ba = area.get(normalize_id(g), float("nan"))
        if not (math.isfinite(ba) and ba > 0):
            raise SystemExit(f"gage {g}: no finite basin_area in the metadata — the February compile normalised every USGS gage")
        conv = conversion_for(ba)
        retrieved = fetch_meta.get(g, {}).get("retrieved")
        for _, r in targets[targets["gage_id"] == g].iterrows():
            d = r["date"].isoformat()
            if rec is None:
                acc, val, reason, qual, raw = False, None, "source_" + str(info.get("error")).replace(" ", "_")[:48], None, None
            elif d not in rec:
                acc, val, reason, qual, raw = False, None, "absent", None, None
            else:
                raw, qual = rec[d]
                acc, val, reason = accept(raw, qual, info["no_data_value"])
            led_rows.append({"gage_id": g, "date": d, "old_flag": r["flag"], "old_Q": "", "value_today_cfs": val if acc else (val if reason == "status" else None),
                             "value_today_raw": raw, "qualifiers_today": qual,
                             "status_today": status_letter(qual) if qual is not None else None,
                             "Q_new": (val * conv) if acc else None, "conversion": conv, "basin_area_km2": ba,
                             "unit_today": info.get("unit") if rec is not None else None,
                             "method_id": info.get("method_id") if rec is not None else None,
                             "retrieved": retrieved, "status": "restored" if acc else f"unavailable_{reason}"})
    led = pd.DataFrame(led_rows)
    led_path = os.path.join(a.out_dir, "correction_ledger.csv")
    led.to_csv(led_path, index=False, float_format="%.17g")
    n_rest = int((led["status"] == "restored").sum())
    log(f"ledger: {len(led):,} rows, {n_rest:,} restored, {len(led) - n_rest:,} unavailable "
        f"({led.loc[led['status'] != 'restored', 'status'].value_counts().to_dict()}); unusable source records {len(unusable)} {dict(list(unusable.items())[:5])}")

    # the patched parquet: identical schema, row groups and row order; Q replaced on restored rows
    rest = led[led["status"] == "restored"]
    new_q = {(g, d): v for g, d, v in zip(rest["gage_id"], pd.to_datetime(rest["date"]).dt.date, rest["Q_new"])}
    rest_gages = set(rest["gage_id"])
    f = pq.ParquetFile(a.parquet)
    schema = f.schema_arrow
    codec = f.metadata.row_group(0).column(0).compression.lower()
    out_pq = os.path.join(a.out_dir, a.parquet_out)
    n_replaced = 0
    spliced_series = {g: [] for g in rest_gages}     # (dates, Q) of the spliced gages for the metadata rule
    with pq.ParquetWriter(out_pq + ".tmp", schema, compression=codec) as w:
        for i in range(f.num_row_groups):
            t = f.read_row_group(i)
            gid = np.asarray(t.column("gage_id").to_pylist(), dtype=object)
            hit = np.isin(gid, list(rest_gages))
            if hit.any():
                q = t.column("Q").to_numpy(zero_copy_only=False).astype(float)
                qnull = np.asarray(t.column("Q").is_null().to_pylist())
                q = np.where(qnull, np.nan, q)
                dates = t.column("Date").to_pandas().to_numpy()
                flags = t.column("flag").to_pylist()
                for j in np.flatnonzero(hit):
                    v = new_q.get((gid[j], dates[j]))
                    if v is not None:
                        if not np.isnan(q[j]):
                            raise SystemExit(f"gage {gid[j]} {dates[j]}: the February row is not NA ({q[j]}) — the ledger targets masked days only")
                        if flags[j] not in SPLICE_CODES:
                            raise SystemExit(f"gage {gid[j]} {dates[j]}: the February flag is {flags[j]!r}, not a splice code")
                        q[j] = v; n_replaced += 1
                for g in set(gid[hit]):
                    sel = gid == g
                    spliced_series[g].append((dates[sel], q[sel]))
                qcol = pa.array(q, type=schema.field("Q").type, mask=np.isnan(q))
                t = t.set_column(t.schema.get_field_index("Q"), "Q", qcol)
            w.write_table(t, row_group_size=t.num_rows)
            if (i + 1) % 200 == 0:
                log(f"  row group {i + 1}/{f.num_row_groups}, {n_replaced:,} replaced so far")
    os.replace(out_pq + ".tmp", out_pq)
    if n_replaced != n_rest:
        raise SystemExit(f"replaced {n_replaced} cells but the ledger restores {n_rest}")
    log(f"patched parquet written: {out_pq} ({n_replaced:,} cells replaced)")

    # metadata rows of the spliced gages updated in place (decision 7): the compile's year rule
    # on the patched series; error_message cleared; every other row byte-identical
    meta_out = meta.drop(columns=["norm_id"]).copy()
    meta_idx = {nid: i for i, nid in enumerate(meta["norm_id"])}
    meta_changes = []
    for g, parts in spliced_series.items():
        dates = np.concatenate([p[0] for p in parts]); q = np.concatenate([p[1] for p in parts])
        st = compile_year_stats(dates, q)
        i = meta_idx.get(normalize_id(g))
        if i is None:
            raise SystemExit(f"gage {g}: no metadata row")
        old = {k: meta_out.at[i, k] for k in ("num_years", "start_year", "end_year", "num_days", "error_message")}
        for k in ("num_years", "start_year", "end_year", "num_days"):
            meta_out.at[i, k] = "" if st[k] is None else str(st[k])
        meta_out.at[i, "error_message"] = ""
        meta_changes.append({"gage_id": meta_out.at[i, "gage_id"], **{f"old_{k}": v for k, v in old.items()}, **{f"new_{k}": st[k] for k in st},
                             "values_changed": any(str(old[k]) != ("" if st[k] is None else str(st[k])) for k in ("num_years", "start_year", "end_year", "num_days"))})
    meta_out_path = os.path.join(a.out_dir, a.metadata_out)
    meta_out.to_csv(meta_out_path, index=False)
    mc = pd.DataFrame(meta_changes); mc.to_csv(os.path.join(a.out_dir, "metadata_changes.csv"), index=False)
    log(f"metadata: {len(mc)} rows updated in place ({int(mc['values_changed'].sum()) if len(mc) else 0} with changed year totals) -> {meta_out_path}")
    prov = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "policy": {"splice_codes": list(SPLICE_CODES), "rule": "status letter A or P, finite value, not the no-data sentinel (R/usgs_qualifiers.R); only the unqualified method block of the one 00060/00003 series, unit ft3/s",
                   "scope": "the masked splice-code days of the February parquet only; provisional-remark days out of scope (decision 1e)",
                   "conversion": "60*60*24/(basin_area_km2*3280.84^3)*1e6 (cfs -> mm/day, as the R ingestion)"},
        "inputs": {"base_parquet": file_info(a.parquet), "metadata": file_info(a.metadata), "masked_rows": file_info(a.masked),
                   "fetch_meta": file_info(os.path.join(a.cache_dir, "fetch_meta.json")), "cache_dir": os.path.abspath(a.cache_dir)},
        "outputs": {"patched_parquet": file_info(out_pq), "ledger": file_info(led_path), "metadata_patched": file_info(meta_out_path)},
        "counts": {"target_days": int(len(led)), "restored": n_rest, "by_status": led["status"].value_counts().to_dict(),
                   "restored_by_code": rest["old_flag"].value_counts().to_dict(), "gages_targeted": len(gages),
                   "gages_restored": len(rest_gages), "metadata_rows_updated": len(mc), "unusable_source_records": unusable},
        "series": {g: series_info[g] for g in gages},
        "tool_sha256": sha256_file(os.path.abspath(__file__)),
    }
    json.dump(prov, open(os.path.join(a.out_dir, "provenance.json"), "w"), indent=1, default=str)
    log("provenance.json written")


# --------------------------------------------------------------------------------------
# control: old rule vs new rule on the October fetch; October vs February (B2)
# --------------------------------------------------------------------------------------
def cmd_control(a):
    os.makedirs(a.out_dir, exist_ok=True)
    meta = read_metadata(a.metadata)
    area = basin_areas(meta)
    uni = sorted(set(pd.read_csv(a.universe, dtype=str)["gage_id"].astype(str)))
    rows, yrows = [], []
    t0 = time.time()
    seen = set()
    for g, df in iter_parquet_gages(a.parquet, only=set(uni)):
        seen.add(g)
        rec, info = load_today(a.cache_dir, g)
        if rec is None:
            rows.append({"gage_id": g, "error": info.get("error")}); continue
        conv = conversion_for(area[normalize_id(g)])
        feb_dates = pd.to_datetime(df["Date"]).dt.date
        feb_key = feb_dates.astype(str).to_numpy()
        feb_q = df["Q"].to_numpy(dtype=float)
        feb_flag = df["flag"].to_numpy(dtype=object)
        present = np.array([d in rec for d in feb_key])
        acc = np.zeros(len(df), bool); tv = np.full(len(df), np.nan); old_keep = np.zeros(len(df), bool)
        for i, d in enumerate(feb_key):
            if d in rec:
                raw, qual = rec[d]
                ok, v, _ = accept(raw, qual, info["no_data_value"])
                if v is not None:
                    tv[i] = v
                acc[i] = ok
                old_keep[i] = (qual in KEPT_CODES_FEB) and v is not None
        feb_has = np.isfinite(feb_q)
        tq = tv * conv
        both = feb_has & acc
        rel = np.abs(tq[both] - feb_q[both]) / np.maximum(np.abs(feb_q[both]), 1e-300)
        nz = both & (feb_q != 0)          # a zero flow equals its conversion under any area: no evidence
        rel_nz = np.abs(tq[nz] - feb_q[nz]) / np.abs(feb_q[nz])
        restored_by_rule = acc & ~old_keep & ~feb_has
        wy = np.where(pd.DatetimeIndex(pd.to_datetime(df["Date"])).month >= 10, pd.DatetimeIndex(pd.to_datetime(df["Date"])).year + 1, pd.DatetimeIndex(pd.to_datetime(df["Date"])).year)
        rows.append({"gage_id": g, "n_feb_rows": len(df), "n_feb_with_value": int(feb_has.sum()),
                     "n_feb_dates_in_today": int(present.sum()), "n_feb_dates_absent_today": int((~present).sum()),
                     "n_feb_value_today_not_accepted": int((feb_has & ~acc).sum()),
                     "n_compared": int(both.sum()), "n_equal_1e9": int((rel <= 1e-9).sum()), "n_differ": int((rel > 1e-9).sum()),
                     "n_compared_nonzero": int(nz.sum()), "n_equal_1e9_nonzero": int((rel_nz <= 1e-9).sum()),
                     "max_rel_diff": float(rel.max()) if len(rel) else 0.0,
                     "n_restored_by_rule": int(restored_by_rule.sum()),
                     "n_restored_by_rule_splice_codes": int((restored_by_rule & np.isin(feb_flag, list(SPLICE_CODES))).sum()),
                     "n_feb_na_old_rule_value_today": int((~feb_has & old_keep).sum()),
                     "n_today_dates_not_in_feb": int(len(set(rec) - set(feb_key))), "unit_today": info.get("unit")})
        for y in np.unique(wy):
            s = wy == y
            r_ = rel if False else None
            b_ = both & s
            rr = np.abs(tq[b_] - feb_q[b_]) / np.maximum(np.abs(feb_q[b_]), 1e-300)
            yrows.append({"gage_id": g, "water_year": int(y), "n_feb_with_value": int((feb_has & s).sum()), "n_compared": int(b_.sum()),
                          "n_equal_1e9": int((rr <= 1e-9).sum()), "n_differ": int((rr > 1e-9).sum()),
                          "n_restored_by_rule": int((restored_by_rule & s).sum())})
    missing = sorted(set(uni) - seen)
    ctl = pd.DataFrame(rows)
    ctl.to_csv(os.path.join(a.out_dir, "control_rule_vs_vintage.csv"), index=False)
    json.dump({"generated": dt.datetime.now().isoformat(timespec="seconds"), "parquet": file_info(a.parquet), "metadata": file_info(a.metadata),
               "fetch_meta": file_info(os.path.join(a.cache_dir, "fetch_meta.json")), "universe": file_info(a.universe),
               "control_csv": file_info(os.path.join(a.out_dir, "control_rule_vs_vintage.csv"))},
              open(os.path.join(a.out_dir, "control_provenance.json"), "w"), indent=1, default=str)
    pd.DataFrame(yrows).to_csv(os.path.join(a.out_dir, "control_rule_vs_vintage_years.csv"), index=False)
    ok = ctl[ctl["error"].isna()] if "error" in ctl.columns else ctl
    log(f"control: {len(ctl)} gages in {time.time() - t0:.0f}s (universe gages absent from the parquet: {len(missing)} {missing[:5]}); "
        f"compared {int(ok['n_compared'].sum()):,} February values with today: {int(ok['n_equal_1e9'].sum()):,} equal at 1e-9 "
        f"({100 * ok['n_equal_1e9'].sum() / max(1, ok['n_compared'].sum()):.3f} %), {int(ok['n_differ'].sum()):,} differ (agency revisions since "
        f"February; worst gage share {(ok['n_differ'] / ok['n_compared'].clip(lower=1)).max():.4f}); restored by the rule "
        f"{int(ok['n_restored_by_rule'].sum()):,} ({int(ok['n_restored_by_rule_splice_codes'].sum()):,} on splice codes); "
        f"February-NA days a kept code carries today {int(ok['n_feb_na_old_rule_value_today'].sum()):,}")
    if missing:
        sys.exit(1)


# --------------------------------------------------------------------------------------
# audit: gates B1–B6, B8 -> input_audit.md
# --------------------------------------------------------------------------------------
def cmd_audit(a):
    os.makedirs(a.out_dir, exist_ok=True)
    lines = ["# Input audit — qualifier-mask splice (plan §5, gates B1–B6 and B8)", "",
             f"Generated {dt.datetime.now().isoformat(timespec='seconds')} by splice_nwis_qualifiers.py audit.", ""]
    gates = {}; diags = []

    def gate(name, ok, detail):
        gates[name] = bool(ok); lines.append(f"- **{name}: {'PASS' if ok else 'FAIL'}** — {detail}"); log(f"{name}: {'PASS' if ok else 'FAIL'} — {detail}")

    meta_b = read_metadata(a.metadata); meta_p = read_metadata(a.metadata_patched)
    meta_idx_b = {nid: i for i, nid in enumerate(meta_b["norm_id"])}
    area = basin_areas(meta_b)
    led = read_ledger(a.ledger)
    rest = led[led["status"] == "restored"]
    rest_keys = {(g, d): (v, f) for g, d, v, f in zip(rest["gage_id"], rest["date"], rest["Q_new"], rest["old_flag"])}
    audit_uni = sorted(set(pd.read_csv(a.audit_universe, dtype=str)["gage_id"].astype(str))) if a.audit_universe else sorted(set(led["gage_id"]))
    fetch_meta = json.load(open(os.path.join(a.cache_dir, "fetch_meta.json")))
    # ---- B1: every cache reopened (never the report alone)
    records = {}; problems = []; multi = []; cover = {}
    for g in audit_uni:
        rec, info = load_today(a.cache_dir, g)
        if rec is None:
            problems.append((g, info.get("error")))
        else:
            records[g] = (rec, info)
            if info.get("n_blocks", 1) > 1:
                multi.append(g)
            cover[g] = (info["first_date"], info["last_date"])
    url_windows = {}
    for g in audit_uni:
        u = fetch_meta.get(g, {}).get("url", "")
        qs = urllib.parse.parse_qs(urllib.parse.urlparse(u).query)
        url_windows[(qs.get("startDT", [None])[0], qs.get("endDT", [None])[0])] = url_windows.get((qs.get("startDT", [None])[0], qs.get("endDT", [None])[0]), 0) + 1
    expect_w = tuple(a.expect_window.split(".."))
    gate("B1 fetch completeness", not problems and list(url_windows) == [expect_w],
         f"{len(audit_uni)} listed gages reopened from the cache; unusable records {len(problems)} {problems[:8]}; sites with several method "
         f"blocks {len(multi)} {multi} (resolved to the unqualified block); request windows {url_windows} (expected {expect_w}); served first dates "
         f"{pd.Series([c[0] for c in cover.values()]).value_counts().head(2).to_dict()}, last dates {pd.Series([c[1] for c in cover.values()]).value_counts().head(2).to_dict()}")
    # ---- B5 + the keyed identity: scan both parquets row group by row group
    fb, fp = pq.ParquetFile(a.parquet), pq.ParquetFile(a.patched)
    b5 = {"schema_equal": fb.schema_arrow.equals(fp.schema_arrow), "row_groups_equal": fb.num_row_groups == fp.num_row_groups,
          "rows_equal": fb.metadata.num_rows == fp.metadata.num_rows, "non_q_columns_identical": True, "row_order_equal": True,
          "q_changed_cells": 0, "q_changed_off_ledger": 0, "ledger_keys_seen": 0, "ledger_base_not_na": 0, "ledger_base_flag_not_splice": 0,
          "ledger_patched_mismatch": 0, "patched_par1": par1_footer(a.patched)}
    base_targets = {}          # (gage, date) -> flag for every masked splice-code day of the base parquet (recomputed, not trusted)
    parquet_gages = set(); gage_rows_b = {}
    patched_series = {g: [] for g in set(rest["gage_id"])}; base_series = {g: [] for g in set(rest["gage_id"])}
    if b5["row_groups_equal"]:
        for i in range(fb.num_row_groups):
            tb, tp = fb.read_row_group(i), fp.read_row_group(i)
            for c in tb.column_names:
                if c != "Q" and not tb.column(c).equals(tp.column(c)):
                    b5["non_q_columns_identical"] = False
                    if c == "gage_id":
                        b5["row_order_equal"] = False
            gid = tb.column("gage_id").to_pylist(); dts = tb.column("Date").to_pandas().to_numpy(); fl = tb.column("flag").to_pylist()
            qb = tb.column("Q").to_pandas().to_numpy(dtype=float); qp = tp.column("Q").to_pandas().to_numpy(dtype=float)
            parquet_gages.update(gid)
            for g, n in zip(*np.unique(np.array(gid, dtype=object), return_counts=True)):
                gage_rows_b[g] = gage_rows_b.get(g, 0) + int(n)
            for j in range(len(gid)):
                key = (gid[j], dts[j])
                if np.isnan(qb[j]) and fl[j] in SPLICE_CODES and normalize_id(gid[j]) in area:
                    base_targets[key] = fl[j]
                if key in rest_keys:
                    b5["ledger_keys_seen"] += 1
                    if not np.isnan(qb[j]):
                        b5["ledger_base_not_na"] += 1
                    if fl[j] not in SPLICE_CODES:
                        b5["ledger_base_flag_not_splice"] += 1
                    if not (qp[j] == rest_keys[key][0]):
                        b5["ledger_patched_mismatch"] += 1
                elif not ((np.isnan(qb[j]) and np.isnan(qp[j])) or qb[j] == qp[j]):
                    b5["q_changed_off_ledger"] += 1
            diff = ~((np.isnan(qb) & np.isnan(qp)) | (qb == qp)); b5["q_changed_cells"] += int(diff.sum())
            if patched_series:
                g_arr = np.array(gid, dtype=object)
                for g in set(g_arr) & set(patched_series):
                    sel = g_arr == g
                    patched_series[g].append((dts[sel], qp[sel])); base_series[g].append((dts[sel], qb[sel]))
    b5["patched_sha256"] = sha256_file(a.patched); b5["base_sha256"] = sha256_file(a.parquet)
    gate("B5 cell-level diff of the whole parquet",
         all(b5[k] for k in ("schema_equal", "row_groups_equal", "rows_equal", "non_q_columns_identical", "row_order_equal", "patched_par1"))
         and b5["q_changed_off_ledger"] == 0 and b5["ledger_keys_seen"] == len(rest_keys) and b5["ledger_base_not_na"] == 0
         and b5["ledger_base_flag_not_splice"] == 0 and b5["ledger_patched_mismatch"] == 0,
         f"schema {b5['schema_equal']}, row groups {b5['row_groups_equal']}, rows {b5['rows_equal']}, non-Q columns identical {b5['non_q_columns_identical']}, "
         f"row order {b5['row_order_equal']}, Q cells changed {b5['q_changed_cells']:,} (off-ledger {b5['q_changed_off_ledger']}); ledger keys found {b5['ledger_keys_seen']:,} "
         f"of {len(rest_keys):,} (base not NA {b5['ledger_base_not_na']}, base flag not a splice code {b5['ledger_base_flag_not_splice']}, patched Q ≠ ledger "
         f"{b5['ledger_patched_mismatch']}); PAR1 {b5['patched_par1']}; sha256 patched {b5['patched_sha256'][:16]}…, base {b5['base_sha256'][:16]}…")
    # ---- B2: conversion identity (control), coverage of every splice gage
    ctl = pd.read_csv(a.control, dtype={"gage_id": str})
    ctl_ok = ctl[ctl["error"].isna()] if "error" in ctl.columns else ctl
    # the control is bound to the metadata and the cache this audit sees
    cprov_path = os.path.join(os.path.dirname(os.path.abspath(a.control)), "control_provenance.json")
    cprov = json.load(open(cprov_path)) if os.path.isfile(cprov_path) else {}
    ctl_bound = (cprov.get("metadata", {}).get("sha256") == sha256_file(a.metadata)
                 and cprov.get("fetch_meta", {}).get("sha256") == sha256_file(os.path.join(a.cache_dir, "fetch_meta.json"))
                 and cprov.get("control_csv", {}).get("sha256") == sha256_file(a.control))
    splice_gages = sorted(set(led["gage_id"]))
    missing_ctl = sorted(set(splice_gages) - set(ctl_ok["gage_id"]))
    share = ctl_ok["n_equal_1e9"].sum() / max(1, ctl_ok["n_compared"].sum())
    worst = (ctl_ok["n_differ"] / ctl_ok["n_compared"].clip(lower=1))
    # per gage: the ±3-day neighbours of every restored day must equal value_today × conversion at 1e-9 —
    # the identity that proves the area and the conversion for THAT gage (a revised gage reports elsewhere)
    nb_bad = {}; nb_n = 0
    for g, df in rest.groupby("gage_id"):
        rec, info = records.get(g, (None, {}))
        if rec is None:
            continue
        conv = conversion_for(area[normalize_id(g)])
        bq = {}
        for p_ in base_series.get(g, []):
            for d_, q_ in zip(p_[0], p_[1]):
                bq[d_] = q_
        n_eq = n_cmp_ = 0
        for d in df["date"]:
            for k in (-3, -2, -1, 1, 2, 3):
                dd = d + dt.timedelta(days=k)
                qb_ = bq.get(dd)
                if qb_ is None or np.isnan(qb_):
                    continue
                raw, qual = rec.get(dd.isoformat(), (None, None)); ok, v, _ = accept(raw, qual, info["no_data_value"])
                if not ok:
                    continue
                n_cmp_ += 1
                if abs(v * conv - qb_) <= 1e-9 * abs(qb_):
                    n_eq += 1
        nb_n += n_cmp_
        if n_cmp_ and n_eq < n_cmp_:
            nb_bad[g] = (n_cmp_ - n_eq, n_cmp_)
    # the per-gage GATE proves the area and the conversion of every restored gage: a wrong basin area
    # makes every converted value differ, so the gage's share of compared non-zero February values
    # equal at 1e-9 must stay above --b2-gage-min-share (agency revisions since February lower it a
    # little and are reported, never applied); the neighbour counts above are the diagnostic
    ctl_g = ctl_ok.set_index("gage_id")
    gage_share = {g: (ctl_g.at[g, "n_equal_1e9_nonzero"] / max(1, ctl_g.at[g, "n_compared_nonzero"])) for g in splice_gages if g in ctl_g.index}
    no_nonzero = [g for g in splice_gages if g in ctl_g.index and ctl_g.at[g, "n_compared_nonzero"] == 0]
    low_share = {g: round(s, 4) for g, s in gage_share.items() if s < a.b2_gage_min_share}
    gate("B2 conversion identity on the February values", ctl_bound and not missing_ctl and share >= a.b2_min_share and not low_share and not no_nonzero,
         f"{int(ctl_ok['n_equal_1e9'].sum()):,} of {int(ctl_ok['n_compared'].sum()):,} stored values equal value_today × conversion at 1e-9 "
         f"({100 * share:.3f} %; the plan's 'overwhelming majority' operationalised as >= {100 * a.b2_min_share:.1f} %); splice gages without a "
         f"control row {len(missing_ctl)} {missing_ctl[:5]}; differing (agency revisions since February, reported, NOT applied): "
         f"{int(ctl_ok['n_differ'].sum()):,}, worst gage share {worst.max():.4f} ({ctl_ok.loc[worst.idxmax(), 'gage_id']}); per-gage: the ±3-day neighbours of every "
         f"restored day ({nb_n:,} compared) equal February at 1e-9 on every gage except {len(nb_bad)} {dict(list(nb_bad.items())[:5])} (diagnostic: revised neighbours); "
         f"per-gage GATE over NON-ZERO February values (a zero equals its conversion under any area): every restored gage's share equal at 1e-9 >= "
         f"{a.b2_gage_min_share:.2f} (min {min(gage_share.values()):.4f} at {min(gage_share, key=gage_share.get)}; below the floor: {low_share}; gages with no "
         f"non-zero evidence: {no_nonzero}); control bound to this metadata and cache: {ctl_bound}; per-year table control_rule_vs_vintage_years.csv")
    # ---- B3: restoration accounting per code over the audit universe (every masked day's state today)
    masked = pq.read_table(a.masked).to_pandas(); masked["date"] = pd.to_datetime(masked["date"]).dt.date
    masked["category"] = masked["flag"].map(categorize)
    b3 = []
    for g, df in masked.groupby("gage_id"):
        rec, info = records.get(g, (None, {}))
        for _, r in df.iterrows():
            d = r["date"].isoformat()
            if rec is None:
                state, st = "no_record", ""
            elif d not in rec:
                state, st = "absent", ""
            else:
                raw, qual = rec[d]; ok, v, reason = accept(raw, qual, info["no_data_value"])
                state = "finite_accepted" if ok else {"status": "finite_other_status", "sentinel": "sentinel", "nonnumeric": "nonnumeric", "nonfinite": "nonfinite"}[reason]
                st = status_letter(qual)
            b3.append((r["flag"], r["category"], state, st))
    b3 = pd.DataFrame(b3, columns=["flag", "category", "state_today", "status_today"])
    tab = b3.pivot_table(index=["category", "flag"], columns="state_today", values="status_today", aggfunc="size", fill_value=0)
    tab.to_csv(os.path.join(a.out_dir, "b3_restoration_accounting.csv"))
    b3.groupby(["category", "flag", "state_today", "status_today"]).size().rename("days").reset_index().to_csv(os.path.join(a.out_dir, "b3_restoration_accounting_by_status.csv"), index=False)
    n_splice_days = int((b3["category"] == "splice").sum()); n_prov = int((b3["category"] == "provisional_remark").sum())
    n_splice_acc = int(((b3["category"] == "splice") & (b3["state_today"] == "finite_accepted")).sum())
    n_prov_acc = int(((b3["category"] == "provisional_remark") & (b3["state_today"] == "finite_accepted")).sum())
    gate("B3 restoration accounting", n_splice_days == a.expect_splice_days and n_prov == a.expect_provisional_days,
         f"splice-code days {n_splice_days:,} (expected {a.expect_splice_days:,}; accepted today {n_splice_acc:,}; states "
         f"{b3[b3['category'] == 'splice']['state_today'].value_counts().to_dict()}), provisional-remark days {n_prov:,} (expected {a.expect_provisional_days:,}; "
         f"accepted today {n_prov_acc:,}, OUT OF SCOPE by decision 1e); tables b3_restoration_accounting*.csv")
    # ---- B4: every restored value re-derived from the cache and the metadata, exact; plausibility diagnostics
    mism = 0; implaus = []
    for g, df in rest.groupby("gage_id"):
        rec, info = records.get(g, (None, {}))
        conv = conversion_for(area[normalize_id(g)])
        # diagnostics against the FEBRUARY (base) series: its finite range widened by 10 % on the far side of
        # zero, and the ratio of each restored value to the median of its calendar-day ±3 neighbours against the
        # base series' own day-to-neighbour ratio distribution (p0.1–p99.9)
        bdates = np.concatenate([p[0] for p in base_series.get(g, [])]) if base_series.get(g) else np.array([])
        bq = np.concatenate([p[1] for p in base_series.get(g, [])]) if base_series.get(g) else np.array([])
        bmap = {d_: q_ for d_, q_ in zip(bdates, bq) if np.isfinite(q_)}
        fin = np.array([q_ for q_ in bmap.values()])
        if len(fin):
            lo = fin.min() * (1 / 1.1 if fin.min() > 0 else 1.1); hi = fin.max() * (1.1 if fin.max() > 0 else 1 / 1.1)
        else:
            lo, hi = np.nan, np.nan
        def nb_median(d_, src):
            nb = [src[d_ + dt.timedelta(days=k)] for k in (-3, -2, -1, 1, 2, 3) if (d_ + dt.timedelta(days=k)) in src and src[d_ + dt.timedelta(days=k)] > 0]
            return np.median(nb) if nb else np.nan
        base_ratios = np.array([q_ / nb_median(d_, bmap) for d_, q_ in list(bmap.items())[::7] if q_ > 0 and np.isfinite(nb_median(d_, bmap))])
        p01, p99 = (np.quantile(base_ratios, [0.001, 0.999]) if len(base_ratios) > 20 else (np.nan, np.nan))
        for d, qn in zip(df["date"], df["Q_new"]):
            if rec is None:
                mism += 1; continue
            raw, qual = rec.get(d.isoformat(), (None, None)); ok, v, _ = accept(raw, qual, info["no_data_value"])
            if not ok or v * conv != qn:
                mism += 1
            if len(fin) and not (lo <= qn <= hi):
                implaus.append((g, d, "outside the February range ×1.1"))
            m_ = nb_median(d, bmap)
            if qn > 0 and np.isfinite(m_) and np.isfinite(p01) and not (p01 <= qn / m_ <= p99):
                implaus.append((g, d, f"ratio to the February ±3-day neighbours {qn / m_:.3g} outside the base p0.1–p99.9 [{p01:.3g}, {p99:.3g}]"))
    gate("B4 restored values equal the cached record × the metadata conversion", mism == 0,
         f"{len(rest):,} restored cells re-derived from the cache and the metadata basin area: {mism} mismatches (the patched parquet's values are checked in B5)")
    pd.DataFrame(implaus, columns=["gage_id", "date", "diagnostic"]).to_csv(os.path.join(a.out_dir, "b4_plausibility.csv"), index=False)
    prov_today = rest[rest["status_today"] != "A"]
    diags.append(f"B4 plausibility diagnostics (not gates): {len(implaus)} of {len(rest):,} restored values flagged "
                 f"(range: {sum(1 for x in implaus if x[2].startswith('outside'))}, neighbour ratio: {sum(1 for x in implaus if x[2].startswith('ratio'))}); "
                 f"b4_plausibility.csv — disposition: the restored days are SELECTED by their February approved flags (A R, A <, A >, A [4]) and their values "
                 f"are the agency's record of the October 2026 retrieval under its status today: {len(rest) - len(prov_today):,} approved (A) and "
                 f"{len(prov_today):,} provisional today ({prov_today.groupby(['gage_id', 'status_today']).size().to_dict()}; the ledger's qualifiers_today "
                 f"carries each); extremes are listed for the reviewer, not applied.")
    # ---- B6: metadata
    same_shape = len(meta_b) == len(meta_p) and (meta_b["gage_id"].to_numpy() == meta_p["gage_id"].to_numpy()).all() and list(meta_b.columns) == list(meta_p.columns)
    rest_norm = {normalize_id(g) for g in rest["gage_id"]}
    b6_bad = []; n_unrelated_changed = 0
    if same_shape:
        cols = [c for c in meta_b.columns if c != "norm_id"]
        for i in range(len(meta_b)):
            nid = meta_b.at[i, "norm_id"]
            rb, rp = meta_b.loc[i, cols].astype(str), meta_p.loc[i, cols].astype(str)
            if nid in rest_norm:
                g = next(g for g in rest["gage_id"].unique() if normalize_id(g) == nid)
                dates = np.concatenate([p[0] for p in patched_series[g]]); q = np.concatenate([p[1] for p in patched_series[g]])
                st = compile_year_stats(dates, q)
                exp = {k: ("" if st[k] is None else str(st[k])) for k in ("num_years", "start_year", "end_year", "num_days")}
                got = {k: str(meta_p.at[i, k]) for k in exp}
                if got != exp or str(meta_p.at[i, "error_message"]) != "" or any(rb[c] != rp[c] for c in cols if c not in exp and c != "error_message"):
                    b6_bad.append((meta_b.at[i, "gage_id"], exp, got))
            elif (rb != rp).any():
                n_unrelated_changed += 1
    base_repro_bad = []
    for g in rest["gage_id"].unique():
        i = meta_idx_b.get(normalize_id(g))
        dates = np.concatenate([p[0] for p in base_series[g]]); q = np.concatenate([p[1] for p in base_series[g]])
        st = compile_year_stats(dates, q)
        exp = {k: ("" if st[k] is None else str(st[k])) for k in ("num_years", "start_year", "end_year", "num_days")}
        got = {k: str(meta_b.at[i, k]) for k in exp}
        if got != exp:
            base_repro_bad.append((g, exp, got))
    not_joining = sorted(g for g in parquet_gages if normalize_id(g) not in set(meta_b["norm_id"]))
    gate("B6 metadata updated in place", same_shape and not b6_bad and n_unrelated_changed == 0 and not not_joining and not base_repro_bad,
         f"same rows/order/columns {same_shape}; the compile rule applied to the BASE series reproduces the February metadata of every restored gage "
         f"(mismatches {len(base_repro_bad)} {base_repro_bad[:2]}); restored gages {len(rest_norm)} recomputed from the patched series — mismatches {len(b6_bad)} {b6_bad[:3]}; "
         f"unrelated rows changed {n_unrelated_changed}; parquet gages not joining the metadata {len(not_joining)} {not_joining[:5]}; normalized ids unique (asserted on read)")
    # ---- B8: source / ledger / output join over the recomputed base targets
    elig_not_restored = []; ledger_not_target = 0; ineligible_restored = 0; n_elig = 0
    led_keys = set(zip(led["gage_id"], led["date"])); led_gages = set(led["gage_id"])
    target_gages = {g for g, _ in base_targets} - (set() if not a.universe_scope else (set(g for g, _ in base_targets) - led_gages))
    no_source = sorted(g for g in target_gages if g not in records)
    for key, fl in base_targets.items():
        g, d = key
        if a.universe_scope and g not in led_gages:
            continue
        rec, info = records.get(g, (None, {}))
        if rec is None:
            continue            # counted in no_source: a coverage failure below, never "ineligible"
        ok = False
        if rec is not None and d.isoformat() in rec:
            raw, qual = rec[d.isoformat()]; ok, v, _ = accept(raw, qual, info["no_data_value"])
        if ok:
            n_elig += 1
            if key not in rest_keys:
                elig_not_restored.append(key)
        elif key in rest_keys:
            ineligible_restored += 1
    for key in led_keys:
        if key not in base_targets:
            ledger_not_target += 1
    gate("B8 no eligible finite source value discarded", not no_source and not elig_not_restored and ineligible_restored == 0 and ledger_not_target == 0 and (n_elig > 0 or not base_targets),
         f"base targets (masked splice-code US days, recomputed from the parquet) {len(base_targets):,} on {len(target_gages)} gages, source records opened for all but "
         f"{len(no_source)} {no_source[:5]} (a missing record is a coverage failure); eligible today {n_elig:,}; eligible but not restored "
         f"{len(elig_not_restored)} {elig_not_restored[:3]}; restored although ineligible {ineligible_restored}; ledger rows that are not base targets {ledger_not_target}")
    lines += [""] + [f"- {d}" for d in diags] + [""]
    verdict = all(gates.values())
    lines.append(f"**Verdict: {'ALL GATES PASS' if verdict else 'FAILED: ' + ', '.join(k for k, v in gates.items() if not v)}** "
                 "(B7, the year oracle vs the census on the patched parquet, is qualifier_year_oracle.py + qualification_census.jl).")
    open(os.path.join(a.out_dir, "input_audit.md"), "w").write("\n".join(lines) + "\n")
    json.dump({"gates": gates, "b5": b5, "diagnostics": diags}, open(os.path.join(a.out_dir, "input_audit.json"), "w"), indent=1, default=str)
    log("input_audit.md written: " + ("ALL GATES PASS" if verdict else "FAILURES"))
    sys.exit(0 if verdict else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list"); p.add_argument("--parquet", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--out-dir", required=True); p.set_defaults(func=cmd_list)
    p = sub.add_parser("fetch"); p.add_argument("--universe", required=True); p.add_argument("--cache-dir", required=True)
    p.add_argument("--start", default="1979-10-01"); p.add_argument("--end", default="2025-10-01")   # the February compile window (config.R DEFAULT_*_DATE, end inclusive)
    p.add_argument("--sleep", type=float, default=1.0); p.add_argument("--retries", type=int, default=4)
    p.set_defaults(func=cmd_fetch)
    p = sub.add_parser("splice"); p.add_argument("--parquet", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--masked", required=True); p.add_argument("--cache-dir", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--universe", default=None, help="splice_universe.csv (default: every gage with a splice-code day)")
    p.add_argument("--parquet-out", default="combined_streamflow_data_09feb2026_qualifier_patched.parquet")
    p.add_argument("--metadata-out", default="combined_watershed_metadata_09feb2026_qualifier_patched.csv")
    p.set_defaults(func=cmd_splice)
    p = sub.add_parser("control"); p.add_argument("--parquet", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--cache-dir", required=True); p.add_argument("--universe", required=True); p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_control)
    p = sub.add_parser("audit"); p.add_argument("--parquet", required=True); p.add_argument("--patched", required=True)
    p.add_argument("--metadata", required=True); p.add_argument("--metadata-patched", required=True); p.add_argument("--ledger", required=True)
    p.add_argument("--masked", required=True); p.add_argument("--cache-dir", required=True); p.add_argument("--control", required=True)
    p.add_argument("--audit-universe", default=None); p.add_argument("--out-dir", required=True)
    p.add_argument("--expect-splice-days", type=int, default=193542); p.add_argument("--expect-provisional-days", type=int, default=15257)
    p.add_argument("--b2-min-share", type=float, default=0.95)
    p.add_argument("--b2-gage-min-share", type=float, default=0.5, help="per gage: the minimum share of its compared February values equal to value_today × conversion at 1e-9 (a wrong basin area gives 0)")
    p.add_argument("--expect-window", default="1979-10-01..2025-10-01", help="the request window every cached response must carry (the compile window)")
    p.add_argument("--universe-scope", action="store_true", help="B8 over the ledger's gages only (a partial splice test)")
    p.set_defaults(func=cmd_audit)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
