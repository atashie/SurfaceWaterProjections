"""Re-fetch the daily records of the 80 `insufficient_years` candidate gages and replay the
ingestion gate (R/helperFunctions.R::process_gages_rawToRaw + generate_streamflow_dt) under
four variants of the per-year rule, to isolate the 0.0001 mm/day flow threshold from the
record-length pre-check and the USGS qualifier-code mask.

USGS: legacy NWIS daily-values service (what dataRetrieval::readNWISdv used in Feb 2026);
Canada: ECCC GeoMet OGC API (hydrometric-daily-mean, served from HYDAT). Raw responses are
cached under <out_dir>/fetch/raw/ so a rerun is offline.

Usage (from the repo root, ~10 min online):
  python docs/benchmarks/min_flow_rule_replay.py <out_dir>
Reads processed_streamflow_data_{AK,canada,conus}/watershed_metadata.csv (16,994 candidates),
writes <out_dir>/fetch/gage_summary.csv and year_counts.csv. Result + classification:
docs/filtering_min_flow_rule_2026-10-21.csv; census addendum
docs/plans/2026-09-04-filtering-stage-census.md §10.
"""
import csv, json, sys, time, urllib.request, urllib.parse, datetime as dt
from pathlib import Path
from collections import defaultdict

S = Path(sys.argv[1])            # output dir
ROOT = Path(__file__).resolve().parents[2]
META = S / "insufficient_years_80.csv"
if not META.exists():                # the 80 insufficient_years rows of the 16,994-candidate metadata
    import pandas as pd
    m = pd.concat([pd.read_csv(ROOT / f"processed_streamflow_data_{d}/watershed_metadata.csv", dtype={"gage_id": str})
                   for d in ("AK", "canada", "conus")], ignore_index=True)
    m = m[m.processing_status == "insufficient_years"].copy()
    m["n_years"] = m.error_message.str.extract(r"Only (\d+) years").astype(int)
    S.mkdir(parents=True, exist_ok=True); m.to_csv(META, index=False)
RAW = S / "fetch" / "raw"; RAW.mkdir(parents=True, exist_ok=True)
START, END = dt.date(1979, 10, 1), dt.date(2025, 10, 1)   # DEFAULT_START_DATE / DEFAULT_END_DATE
THR, MIN_DAYS, MIN_YEARS = 0.0001, 30, 20
KEEP_CODES = {"A", "A e", "P", "P e"}                        # USGS mask in generate_streamflow_dt

def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "hisss-census (urllib)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa
            print("  retry", k, e, flush=True); time.sleep(5 * (k + 1))
    raise RuntimeError(url)

def usgs_daily(site):
    f = RAW / f"usgs_{site}.json"
    if not f.exists():
        q = urllib.parse.urlencode({"format": "json", "sites": site, "parameterCd": "00060",
                                    "startDT": START.isoformat(), "endDT": END.isoformat(), "siteStatus": "all"})
        f.write_bytes(get("https://waterservices.usgs.gov/nwis/dv/?" + q))
    js = json.loads(f.read_bytes())
    series = [ts for ts in js["value"]["timeSeries"]
              if ts["variable"]["variableCode"][0]["value"] == "00060"
              and ts["variable"]["options"]["option"][0].get("optionCode") == "00003"]
    nodata = None
    rows = []
    if series:
        ts = series[0]                         # primary series (= `Flow` after renameNWISColumns)
        nodata = ts["variable"].get("noDataValue")
        for v in ts["values"][0]["value"]:
            d = dt.date.fromisoformat(v["dateTime"][:10])
            val = float(v["value"]) if v["value"] not in (None, "") else None
            if nodata is not None and val is not None and val == nodata:
                val = None
            code = " ".join(v.get("qualifiers") or [])
            rows.append((d, val, code))
    return rows, len(series)

def eccc_daily(stn):
    f = RAW / f"eccc_{stn}.json"
    if not f.exists():
        feats, offset = [], 0
        while True:
            q = urllib.parse.urlencode({"f": "json", "STATION_NUMBER": stn, "limit": 10000, "offset": offset,
                                        "datetime": f"{START.isoformat()}/{END.isoformat()}",
                                        "sortby": "DATE"})
            js = json.loads(get("https://api.weather.gc.ca/collections/hydrometric-daily-mean/items?" + q))
            feats += js["features"]
            if len(js["features"]) < 10000:
                break
            offset += 10000
        f.write_text(json.dumps({"features": feats}))
    js = json.loads(f.read_text())
    rows = []
    for ft in js["features"]:
        p = ft["properties"]
        if p.get("DISCHARGE") is None:          # level-only day: not a Flow row in HYDAT
            continue
        rows.append((dt.date.fromisoformat(p["DATE"][:10]), float(p["DISCHARGE"]), p.get("DISCHARGE_SYMBOL_EN") or ""))
    return rows, 1

def replay(rows, conv, usgs):
    """Return per-calendar-year counts and the gate decisions."""
    rows = [r for r in rows if START <= r[0] <= END]
    n_rows = len(rows)
    pre_gate = n_rows > 365 * MIN_YEARS and (max(r[0] for r in rows) > START if rows else False)
    per = defaultdict(lambda: [0, 0, 0, 0, 0, 0])   # rows, raw_nonNA, kept_nonNA, kept_pos, kept_thr, raw_thr
    for d, val, code in rows:
        y = per[d.year]; y[0] += 1
        if val is None:
            continue
        q = val * conv
        y[1] += 1
        if q > THR: y[5] += 1
        kept = (code in KEEP_CODES) if usgs else True
        if kept:
            y[2] += 1
            if q > 0: y[3] += 1
            if q > THR: y[4] += 1
    yrs = lambda i: sum(1 for y in per.values() if y[i] >= MIN_DAYS)
    return n_rows, pre_gate, per, {"actual": yrs(4), "nonNA_only": yrs(2), "positive_only": yrs(3), "no_flag_mask": yrs(5)}

meta = list(csv.DictReader(open(META)))
out_g = open(S / "fetch" / "gage_summary.csv", "w", newline=""); wg = csv.writer(out_g)
wg.writerow(["gage_id", "gage_type", "basin_area", "area_normalized", "recorded_n_years", "n_series", "n_rows",
             "pre_gate_pass", "yrs_actual", "yrs_nonNA_only", "yrs_positive_only", "yrs_no_flag_mask",
             "codes_seen", "max_q_mmday"])
out_y = open(S / "fetch" / "year_counts.csv", "w", newline=""); wy = csv.writer(out_y)
wy.writerow(["gage_id", "year", "rows", "raw_nonNA", "kept_nonNA", "kept_pos", "kept_thr", "raw_thr"])
for i, m in enumerate(meta, 1):
    gid, gtype = m["gage_id"], m["gage_type"]
    area = float(m["basin_area"]) if m["basin_area"] not in ("", "nan") else None
    print(f"[{i}/{len(meta)}] {gtype} {gid}", flush=True)
    try:
        if gtype == "USGS":
            rows, nser = usgs_daily(gid)
            conv = 1.0 if area is None else 86400 / (area * 3280.84 ** 3) * 1e6
        else:
            rows, nser = eccc_daily(gid)
            conv = 1.0 if area is None else 86400 * 1e9 / (area * 1e12)
        n_rows, pre, per, yrs = replay(rows, conv, gtype == "USGS")
        codes = sorted({r[2] for r in rows if r[1] is not None})
        maxq = max((r[1] * conv for r in rows if r[1] is not None), default=float("nan"))
        wg.writerow([gid, gtype, area, m["area_normalized"], m["n_years"], nser, n_rows, pre,
                     yrs["actual"], yrs["nonNA_only"], yrs["positive_only"], yrs["no_flag_mask"],
                     "|".join(codes), f"{maxq:.4g}"])
        for y in sorted(per):
            wy.writerow([gid, y, *per[y]])
    except Exception as e:  # noqa
        print("  FAILED", gid, e, flush=True)
        wg.writerow([gid, gtype, area, m["area_normalized"], m["n_years"], "", "", "", "", "", "", "", "FAILED", ""])
    out_g.flush(); out_y.flush()
print("done", flush=True)
