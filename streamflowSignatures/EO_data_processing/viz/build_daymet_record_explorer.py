"""Build a self-contained HTML explorer: the reprocessed 1980-2025 Daymet basin series against the
original co-author series, over the whole record.

The page makes no network request (plain JS, canvas and SVG; data inline). It shows
  * the agreement metrics of daymet_validate.py over 1980-2023 for every basin-year with data in
    both files: per-variable distributions of each metric, a quantile table, the one-day lag test
    and the headline numbers;
  * a map of every basin, coloured by coverage or by its record RMSE / largest daily |difference|;
  * daily series of all six variables, original and new side by side plus their difference, for
    3 x N embedded basins: N least matching, N random, N without an original series.

Selection (seeded, so the same inputs give the same page):
  least matching  the basins holding each variable's largest daily |new - original|, then, one
                  variable after another, the next basin by record RMSE (1980-2023)
  random          uniform among the compared basins not already chosen
  new coverage    one basin whose original is NaN on every day, one > 100,000 km2, one
                  low-confidence polygon, the rest uniform among the basins without an original
                  series (absent from the original file, or NaN on every day there)

Encoding. Per basin and variable, the new series is stored as int16 day-to-day steps of a fixed
quantum (0.01 mm, 0.01 degC, 0.1 Pa, 0.01 W/m2; coarser only when the series' range needs more
than 32,000 quanta). The difference new - original is stored as int16 scaled to its own largest
|value|. Both are byte-shuffled and gzipped, and the browser decodes them (DecompressionStream).
The original is drawn as new - difference, so the two panels differ by exactly the stored
difference. Before encoding, the builder checks every embedded basin-year against the validation
table: the per-year largest |difference| recomputed from the extracted series must equal it.

Usage:
    python build_daymet_record_explorer.py --new <daymet_1980_2025_29sep2026.parquet>
        --old <stale co-author parquet> --basins <weights/daymet_basins.csv>
        --polygons-qa <polygons/..._qa.csv> --validation <validation_1980_2023_per_basin.csv>
        --gages-ii-dir <gagesMetadata> [--flags <daymet_basin_flags.csv>]
        [--outputcheck <validation/outputcheck_*.json>] [--n 10] [--seed 20261001]
        --out <page.html> [--fragment <body-only.html>] [--assets-dir <Natural Earth cache>]
"""
import argparse
import base64
import gzip
import json
import os
import re
import sys
import time

import duckdb
import numpy as np
import pandas as pd
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "daymet"))
sys.path.insert(0, HERE)
from daymet_common import DAYMET_VARS, canon, git_state  # noqa: E402
from build_daymet_comparison_dashboard import LCC, NAMES, UNITS, outlines  # noqa: E402

VARS = ["prcp", "tmin", "tmax", "swe", "vp", "srad"]
SHORT = {"prcp": "Precipitation", "tmin": "Tmin", "tmax": "Tmax", "swe": "SWE", "vp": "Vapour pressure",
         "srad": "Shortwave"}
QUANTUM = {"prcp": 0.01, "tmin": 0.01, "tmax": 0.01, "swe": 0.01, "vp": 0.1, "srad": 0.01}
YEAR0, YEAR_OLD_END, YEAR_NEW_END = 1980, 2023, 2025
N_OLD = (YEAR_OLD_END - YEAR0 + 1) * 365
N_NEW = (YEAR_NEW_END - YEAR0 + 1) * 365
LARGE_KM2 = 100_000
QS = {"p01": 0.01, "p50": 0.5, "p99": 0.99, "p999": 0.999}
# metric key -> label, validation column expression, axis kind, which variables
METRICS = [
    ("max_abs", "Largest daily |Δ|", "max_abs_diff", "log", VARS,
     "The largest absolute daily difference |new − original| in each basin-year."),
    ("rmse", "RMSE", "rmse", "log", VARS, "Root-mean-square daily difference in each basin-year."),
    ("mae", "Mean |Δ|", "mae", "log", VARS, "Mean absolute daily difference in each basin-year."),
    ("bias", "Mean Δ", "bias", "lin", VARS,
     "Mean daily difference new − original in each basin-year (signed; a systematic offset would "
     "shift this away from 0)."),
    ("omr2", "1 − R²", "1 - r2_identity", "log", VARS,
     "1 − R² about the 1:1 line (identity R²), per basin-year. The seasonal cycle dominates R², so "
     "the absolute differences are the informative numbers. Basin-years whose original is constant "
     "(all-zero SWE) have no R² and are left out."),
    ("annual", "Annual total Δ", "100 * (sum_ratio - 1)", "lin", ["prcp", "swe", "vp", "srad"],
     "Relative difference of the annual sum, new vs original, in %. Not shown for temperatures: "
     "their annual sums cross zero, so a ratio means nothing."),
]


def shuffle16(a):
    b = np.ascontiguousarray(a, dtype="<i2").view(np.uint8).reshape(-1, 2)
    return np.concatenate([b[:, 0], b[:, 1]]).tobytes()


def unshuffle16(u8, off, n):
    return (u8[off:off + n].astype(np.uint16) | (u8[off + n:off + 2 * n].astype(np.uint16) << 8)).view(np.int16)


def check_block(b64, desc, new, old):
    """Decode a block the way the page does and compare it with the source series."""
    u8 = np.frombuffer(gzip.decompress(base64.b64decode(b64)), dtype=np.uint8)
    off = 0
    for var in VARS:
        d = desc[var]
        k = np.cumsum(unshuffle16(u8, off, N_NEW).astype(np.int64))
        off += 2 * N_NEW
        nv = d["lo"] + k * d["q"]
        if np.abs(nv - new[var].to_numpy(float)).max() > d["q"] / 2 + 1e-9:
            sys.exit(f"{var}: decoded new series off by more than half a quantum")
        if old is not None:
            e = unshuffle16(u8, off, N_OLD).astype(float) * d["ds"]
            off += 2 * N_OLD
            true = new[var].to_numpy(float)[:N_OLD] - old[var].to_numpy(float)
            if np.abs(e - true).max() > d["ds"] / 2 + 1e-15:
                sys.exit(f"{var}: decoded difference off by more than half its quantum")
    if off != len(u8):
        sys.exit("a block has unread bytes")


def daymet_dates(y0, y1):
    return pd.DatetimeIndex(np.concatenate([pd.date_range(f"{y}-01-01", periods=365).to_numpy()
                                            for y in range(y0, y1 + 1)]))


def sig(x, n=3):
    return float(f"{x:.{n}g}") if np.isfinite(x) else None


def histogram(vals, kind, nbins_lin=60, per_decade=8):
    """Bins for one metric and variable. Log: decades from the smallest positive value; exact
    zeros counted apart. Linear: symmetric about 0 out to the 0.1/99.9 percentiles; values beyond
    are counted apart, never clipped silently."""
    v = vals[np.isfinite(vals)]
    out = {"n": int(len(v)), "q": {k: sig(np.quantile(v, q), 4) for k, q in QS.items()},
           "min": sig(v.min(), 4), "max": sig(v.max(), 4)}
    if kind == "log":
        pos = v[v > 0]
        out["zeros"] = int((v == 0).sum())
        out["neg"] = int((v < 0).sum())
        lo, hi = np.floor(np.log10(pos.min())), np.ceil(np.log10(pos.max()))
        if hi == lo:
            hi += 1
        edges = np.linspace(lo, hi, int(hi - lo) * per_decade + 1)
        counts, _ = np.histogram(np.log10(pos), bins=edges)
        out.update({"lo": float(lo), "hi": float(hi), "counts": counts.tolist()})
    else:
        lim = float(max(abs(np.quantile(v, 0.001)), abs(np.quantile(v, 0.999))))
        lim = float(f"{lim * 1.05:.2g}") if lim > 0 else 1.0
        edges = np.linspace(-lim, lim, nbins_lin + 1)
        counts, _ = np.histogram(v[(v >= -lim) & (v <= lim)], bins=edges)
        out.update({"lo": -lim, "hi": lim, "counts": counts.tolist(),
                    "below": int((v < -lim).sum()), "above": int((v > lim).sum())})
    return out


def clean(x):
    return "" if pd.isna(x) else str(x).strip()


SMALL_WORDS = re.compile(r"\b(At|Near|Nr|Of|The|And|Above|Abv|Below|Blw|In|On|To|From|For|By)\b")


def tidy_name(raw):
    """Station names come upper-cased: title-case them, keep small words lower-case, ordinals
    lower-case (232nd) and a trailing state or province code upper-case."""
    t = clean(raw)
    if not t:
        return ""
    t = re.sub(r"(\d)([A-Za-z]+)", lambda m: m.group(1) + m.group(2).lower(), t.title())
    t = SMALL_WORDS.sub(lambda m: m.group(1).lower(), t)
    t = t[0].upper() + t[1:]
    return re.sub(r",\s*([A-Za-z]{2})\.?$", lambda m: ", " + m.group(1).upper(), t)


def station_names(gdir):
    names = {}
    for f in ["conterm_basinid.txt", "AKHIPR_basinid.txt"]:
        p = os.path.join(gdir, f)
        if os.path.exists(p):
            t = pd.read_csv(p, dtype=str, encoding="latin-1", usecols=["STAID", "STANAME", "STATE"])
            for s, n, st in zip(t["STAID"], t["STANAME"], t["STATE"]):
                names[canon(s)] = (tidy_name(n), clean(st))
    p = os.path.join(gdir, "Canadian_gages_goodones.csv")
    if os.path.exists(p):
        t = pd.read_csv(p, dtype=str, encoding="latin-1")
        for s, n, pr in zip(t["STATION_NUMBER"], t["STATION_NAME"], t["PROV_TERR_STATE_LOC"]):
            names.setdefault(canon(s), (tidy_name(n), clean(pr)))
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--basins", required=True)
    ap.add_argument("--polygons-qa", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--gages-ii-dir", required=True)
    ap.add_argument("--flags", default=None)
    ap.add_argument("--outputcheck", default=None)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--memory-limit", default="2GB")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fragment", default=None)
    ap.add_argument("--assets-dir", default=os.path.expanduser("~/HISSS_data/dashboard_assets"))
    a = ap.parse_args()
    t0 = time.time()
    if set(VARS) != set(DAYMET_VARS):
        sys.exit("variable list differs from daymet_common.DAYMET_VARS")
    con = duckdb.connect(config={"memory_limit": a.memory_limit, "threads": 4})

    # ---- basins ---------------------------------------------------------------------------
    b = pd.read_csv(a.basins, dtype={"site_id": str, "canon_id": str, "boundary_gage_id": str})
    q = pd.read_csv(a.polygons_qa, dtype={"gage_id": str, "canon_id": str})
    b = b.merge(q[["canon_id", "latitude", "longitude", "basin_area"]], on="canon_id", how="left",
                validate="one_to_one")
    if b[["latitude", "longitude"]].isna().any().any():
        sys.exit("a basin has no gage coordinates in the polygon QA table")
    names = station_names(a.gages_ii_dir)
    b["name"] = [names.get(c, ("", ""))[0] for c in b["canon_id"]]
    b["region"] = [names.get(c, ("", ""))[1] for c in b["canon_id"]]
    flags = {}
    if a.flags:
        f = pd.read_csv(a.flags, dtype={"site_id": str, "canon_id": str})
        flags = dict(zip(f["canon_id"], f["reason"]))
    old_ids = con.sql(f"""SELECT site_id, count(*) n,
        sum(CASE WHEN prcp IS NULL OR isnan(prcp) THEN 1 ELSE 0 END) n_nan
        FROM read_parquet('{a.old}') GROUP BY site_id""").df()
    old_ids["canon_id"] = old_ids["site_id"].map(canon)
    if old_ids["canon_id"].duplicated().any():
        sys.exit("two original site ids share a canonical id")
    old_spelling = dict(zip(old_ids["canon_id"], old_ids["site_id"]))
    old_allnan = set(old_ids.loc[old_ids["n_nan"] == old_ids["n"], "canon_id"])
    in_old = b["canon_id"].isin(set(old_ids["canon_id"]))
    status = np.where(~in_old, 2, np.where(b["canon_id"].isin(old_allnan), 1, 0))  # 0 both, 1 orig NaN, 2 new only
    b["status"] = status
    idx_of = {c: i for i, c in enumerate(b["canon_id"])}
    n_old_only = int((~old_ids["canon_id"].isin(set(b["canon_id"]))).sum())

    # ---- validation metrics ---------------------------------------------------------------
    src = (f"read_parquet('{a.validation}')" if a.validation.endswith(".parquet") else
           f"read_csv('{a.validation}', types={{'canon_id': 'VARCHAR', 'site_id': 'VARCHAR'}})")
    con.sql(f"""CREATE TEMP TABLE v AS SELECT canon_id, var, year, n_both_finite, bias, mae, rmse,
        max_abs_diff, r, r_lag_m1, r_lag_p1, r2_identity, sum_ratio FROM {src} WHERE n_both_finite > 0""")
    years = con.sql("SELECT min(year), max(year) FROM v").fetchone()
    if tuple(years) != (YEAR0, YEAR_OLD_END):
        sys.exit(f"validation covers {years}, expected {YEAR0}-{YEAR_OLD_END}")
    rec = con.sql("""SELECT canon_id, var, sum(n_both_finite) n, max(max_abs_diff) max_abs,
        arg_max(year, max_abs_diff) max_year, sqrt(sum(rmse * rmse * n_both_finite) / sum(n_both_finite)) rmse
        FROM v GROUP BY canon_id, var""").df()
    compared = sorted(set(rec["canon_id"]))
    if set(compared) != set(b.loc[b["status"] == 0, "canon_id"]):
        sys.exit("the validation table's basins differ from the basins with data in both files")

    metrics, headline = [], {"vars": {}}
    for key, label, expr, kind, vlist, desc in METRICS:
        m = {"key": key, "label": label, "kind": kind, "desc": desc, "vars": {}}
        for var in vlist:
            vals = con.sql(f"SELECT {expr} AS x FROM v WHERE var = '{var}'").df()["x"].to_numpy(dtype=float)
            m["vars"][var] = histogram(vals, kind)
        metrics.append(m)
    lag = {}
    for var in VARS:
        r = con.sql(f"""SELECT count(*) n, median(r_lag_m1) m1, median(r) r0, median(r_lag_p1) p1,
            sum(CASE WHEN r_lag_m1 > r OR r_lag_p1 > r THEN 1 ELSE 0 END) beats,
            sum(n_both_finite) n_day_total, count(DISTINCT canon_id) n_basin FROM v WHERE var = '{var}'""").fetchone()
        lag[var] = {"n": int(r[0]), "m1": sig(r[1], 6), "r0": sig(r[2], 12), "p1": sig(r[3], 6), "beats": int(r[4])}
        top = rec[rec["var"] == var].sort_values("max_abs", ascending=False).iloc[0]
        headline["vars"][var] = {"max_abs": sig(top["max_abs"], 3), "max_site": b.loc[idx_of[top["canon_id"]], "site_id"],
                                 "max_year": int(top["max_year"]), "basin_years": int(r[0]), "days": int(r[5]),
                                 "basins": int(r[6])}
    ann = next(m for m in metrics if m["key"] == "annual")["vars"]["prcp"]
    headline.update({"compared": len(compared), "shared": int(in_old.sum()), "orig_nan": int((status == 1).sum()),
                     "new_only": int((status == 2).sum()), "old_only": n_old_only, "n_new": int(len(b)),
                     "n_old": int(len(old_ids)), "prcp_annual_p01": ann["q"]["p01"], "prcp_annual_p99": ann["q"]["p99"],
                     "prcp_annual_min": ann["min"], "prcp_annual_max": ann["max"]})

    # ---- selection --------------------------------------------------------------------------
    rng = np.random.default_rng(a.seed)
    sel, why = {"worst": [], "random": [], "new": []}, {}

    def add(group, cid, reason):
        why.setdefault(cid, [])
        if reason not in why[cid]:
            why[cid].append(reason)
        if not any(cid in g for g in sel.values()):
            sel[group].append(cid)

    for var in VARS:
        top = rec[rec["var"] == var].sort_values("max_abs", ascending=False).iloc[0]
        add("worst", top["canon_id"], {"kind": "max", "var": var, "value": sig(top["max_abs"]), "year": int(top["max_year"])})
    ranked = {var: rec[rec["var"] == var].sort_values("rmse", ascending=False) for var in VARS}
    k = 0
    while len(sel["worst"]) < a.n:
        for var in VARS:
            if len(sel["worst"]) >= a.n:
                break
            r = ranked[var].iloc[k]
            add("worst", r["canon_id"], {"kind": "rmse", "var": var, "rank": k + 1, "value": sig(r["rmse"])})
        k += 1
    pool = [c for c in compared if c not in sel["worst"]]
    for cid in rng.choice(pool, size=a.n, replace=False):
        add("random", str(cid), {"kind": "random"})
    no_orig = sorted(b.loc[b["status"] > 0, "canon_id"])
    area = dict(zip(b["canon_id"], b["geom_area_km2"]))
    lowconf = dict(zip(b["canon_id"], b["low_confidence"].astype(bool)))
    picks = [("orig_nan", sorted(old_allnan & set(no_orig))),
             ("large", [c for c in no_orig if area[c] > LARGE_KM2]),
             ("low_confidence", [c for c in no_orig if lowconf[c]])]
    for kind, cands in picks:
        cands = [c for c in cands if c not in sel["new"]]
        if cands and len(sel["new"]) < a.n:
            add("new", str(rng.choice(cands)), {"kind": kind})
    rest = [c for c in no_orig if c not in sel["new"]]
    for cid in rng.choice(rest, size=a.n - len(sel["new"]), replace=False):
        add("new", str(cid), {"kind": "random"})

    # ---- daily series of the embedded basins ------------------------------------------------
    chosen = [c for g in ("worst", "random", "new") for c in sel[g]]
    new_sites = [b.loc[idx_of[c], "site_id"] for c in chosen]
    sql_in = ",".join("'" + s.replace("'", "''") + "'" for s in new_sites)
    cols = ", ".join(VARS)
    new = con.sql(f"SELECT site_id, Date, {cols} FROM read_parquet('{a.new}') WHERE site_id IN ({sql_in}) "
                  "ORDER BY site_id, Date").df()
    old_in = [old_spelling[c] for c in chosen if c in old_spelling]
    sql_old = ",".join("'" + s.replace("'", "''") + "'" for s in old_in) or "''"
    old = con.sql(f"SELECT site_id, Date, {cols} FROM read_parquet('{a.old}') WHERE site_id IN ({sql_old}) "
                  "ORDER BY site_id, Date").df()
    new["canon_id"], old["canon_id"] = new["site_id"].map(canon), old["site_id"].map(canon)
    new_dates, old_dates = daymet_dates(YEAR0, YEAR_NEW_END), daymet_dates(YEAR0, YEAR_OLD_END)
    vmax = con.sql("SELECT canon_id, var, year, max_abs_diff FROM v").df().set_index(["canon_id", "var", "year"])["max_abs_diff"]
    blocks, embedded, worst_check = {}, [], 0.0
    for cid in chosen:
        x = new[new["canon_id"] == cid]
        if not (pd.DatetimeIndex(x["Date"]) == new_dates).all() or len(x) != N_NEW:
            sys.exit(f"{cid}: the new series is not the 1980-2025 Daymet calendar")
        o = old[old["canon_id"] == cid]
        has_old = cid in old_spelling and cid not in old_allnan
        if has_old and (len(o) != N_OLD or not (pd.DatetimeIndex(o["Date"]) == old_dates).all()):
            sys.exit(f"{cid}: the original series is not the 1980-2023 Daymet calendar")
        raw, desc = b"", {}
        for var in VARS:
            xv = x[var].to_numpy(dtype=float)
            if not np.isfinite(xv).all():
                sys.exit(f"{cid} {var}: the new series has a NaN")
            lo = float(np.floor(xv.min() / QUANTUM[var]) * QUANTUM[var])
            qn = max(QUANTUM[var], (xv.max() - lo) / 32000)
            kq = np.rint((xv - lo) / qn).astype(np.int64)
            steps = np.diff(kq, prepend=0)
            if np.abs(steps).max() > 32767:
                sys.exit(f"{cid} {var}: a day-to-day step exceeds int16")
            raw += shuffle16(steps)
            d = {"lo": lo, "q": qn, "err": sig(np.abs(lo + kq * qn - xv).max(), 2)}
            if has_old:
                dv = xv[:N_OLD] - o[var].to_numpy(dtype=float)
                per_year = np.abs(dv).reshape(-1, 365).max(axis=1)
                ref = vmax.loc[cid, var].reindex(range(YEAR0, YEAR_OLD_END + 1)).to_numpy()
                worst_check = max(worst_check, float(np.nanmax(np.abs(per_year - ref))))
                m = float(np.abs(dv).max())
                s = m / 32767 if m > 0 else 1.0
                raw += shuffle16(np.rint(dv / s))
                d["ds"] = s
            desc[var] = d
        key = f"b{len(blocks)}"
        blocks[key] = base64.b64encode(gzip.compress(raw, 9, mtime=0)).decode("ascii")
        check_block(blocks[key], desc, x, o if has_old else None)
        group = next(g for g in sel if cid in sel[g])
        embedded.append({"i": idx_of[cid], "group": group, "why": why[cid], "block": key, "hasOld": bool(has_old),
                         "origNaN": cid in old_allnan, "vars": desc, "flag": flags.get(cid, "")})
    if worst_check > 1e-12:
        sys.exit(f"extracted series disagree with the validation table (per-year max |diff| off by {worst_check:.3g})")

    # ---- map geometry and per-basin record metrics ----------------------------------------------
    to_lcc = Transformer.from_crs("EPSG:4326", LCC, always_xy=True)
    bx, by = to_lcc.transform(b["longitude"].to_numpy(float), b["latitude"].to_numpy(float))
    cidx = np.array([idx_of[c] for c in compared])
    recmap = {}
    for metric in ("rmse", "max_abs"):
        piv = rec.pivot(index="canon_id", columns="var", values=metric).reindex(compared)
        recmap[metric] = {var: [sig(v, 3) for v in piv[var].to_numpy()] for var in VARS}

    check = None
    if a.outputcheck:
        oc = json.load(open(a.outputcheck))
        check = {"years": oc["years"], "n_sites": len(oc["sites"]), "verdict": oc["verdict"],
                 "checked": oc["checked_utc"][:10],
                 "max": {var: sig(max(s["max_abs_diff"][var] for s in oc["sites"]), 2) for var in VARS},
                 "sites": [s["site_id"] for s in oc["sites"]]}
    gs = git_state(HERE)
    prov = json.load(open(a.new + ".provenance.json")) if os.path.exists(a.new + ".provenance.json") else {}
    data = {
        "meta": {"generated": time.strftime("%Y-%m-%d"), "commit": (gs["commit"] or "")[:7], "n": a.n, "seed": a.seed,
                 "new_file": os.path.basename(a.new), "new_md5": prov.get("out_md5"),
                 "old_file": os.path.basename(a.old), "old_bytes": os.path.getsize(a.old),
                 "validation_file": os.path.basename(a.validation), "year0": YEAR0, "yearOldEnd": YEAR_OLD_END,
                 "yearNewEnd": YEAR_NEW_END, "nOld": N_OLD, "nNew": N_NEW, "largeKm2": LARGE_KM2},
        "vars": VARS, "names": NAMES, "short": SHORT, "units": UNITS,
        "headline": headline, "metrics": metrics, "lag": lag, "check": check,
        "basins": {"id": b["site_id"].tolist(), "name": b["name"].tolist(), "region": b["region"].tolist(),
                   "country": ["US" if t == "USGS" else "CA" for t in b["gage_type"].fillna("")],
                   "area": [sig(v, 3) for v in b["geom_area_km2"]],
                   "x": np.round(np.asarray(bx) / 1000.0, 1).tolist(), "y": np.round(np.asarray(by) / 1000.0, 1).tolist(),
                   "status": b["status"].astype(int).tolist(), "lowconf": b["low_confidence"].astype(bool).astype(int).tolist(),
                   "lowpix": b["low_pixel_support"].astype(bool).astype(int).tolist(),
                   "cells": b["n_cells"].astype(int).tolist()},
        "rec": {"idx": cidx.tolist(), **recmap},
        "embedded": embedded,
        "outlines": outlines(a.assets_dir, to_lcc),
    }
    payload = json.dumps(data, separators=(",", ":"), allow_nan=False, ensure_ascii=False).replace("</", "<\\/")
    bl = json.dumps(blocks, separators=(",", ":"))
    tpl = open(os.path.join(HERE, "daymet_record_explorer_template.html"), encoding="utf-8").read()
    page = tpl.replace("__DATA__", payload).replace("__BLOCKS__", bl)
    if a.fragment:
        with open(a.fragment, "w", encoding="utf-8") as fh:
            fh.write(page)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
                 "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
                 "</head>\n<body>\n" + page + "\n</body>\n</html>\n")
    print(f"embedded {len(embedded)} basins (worst {len(sel['worst'])}, random {len(sel['random'])}, "
          f"new {len(sel['new'])}); series check vs validation table: max |per-year max diff| gap {worst_check:.2g}; "
          f"metadata {len(payload) / 1e6:.2f} MB, series {len(bl) / 1e6:.2f} MB; page "
          f"{os.path.getsize(a.out) / 1e6:.1f} MB in {time.time() - t0:.0f} s")
    for g in ("worst", "random", "new"):
        print(f"  {g}: " + ", ".join(b.loc[idx_of[c], "site_id"] for c in sel[g]))


if __name__ == "__main__":
    main()
