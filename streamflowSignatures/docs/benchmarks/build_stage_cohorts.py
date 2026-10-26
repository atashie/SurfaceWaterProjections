#!/usr/bin/env python3
"""build_stage_cohorts.py — the executed-subset manifests of the attribution stages S1 and S1′
(docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.4, decision 14; 2026-10-26), composed
from the REALIZED ledgers, never from the forecast.

    python docs/benchmarks/build_stage_cohorts.py --window 1993-2025 --selection selection_1993-2025.txt
        --ledger correction_ledger.csv --inserted inserted_rows_all.csv --old-annual V2_annual.parquet
        --out-dir manifests [--control 200] [--join-only 200]

Per window, over the product's frozen cohort (the selection manifest):
  A-Q  = gages with >= 1 restored ledger day in a water year of the window;
  J    = gages with >= 1 inserted climate row in a water year that is Q-valid in version 2
         (holds annual rows) — the join-affected gages of §10.3;
  S1   = A-Q + a control sample of --control gages that are in neither set (evenly spaced
         over the sorted ids, deterministic);
  S1′  = the interaction cohort: A-Q ∩ J, plus every A-Q gage whose restored (newly valid or
         restored-day) years hold inserted rows (the realized overlap of §10.3), plus the
         eligibility-gain cases (a Q-valid v2 year with inserted rows and no climate annual
         row in v2 — a climate-invalid year that the frame can make valid), plus the Dec 31
         neighbour cases (inventory kind daymet_missing_dec31_neighbour in a Q-valid year),
         plus --join-only join-only gages evenly spaced.
Writes stage_S1_<window>.txt, stage_S1prime_<window>.txt and stage_cohorts_<window>.json
(membership by reason, counts, sha256).
"""
import argparse
import hashlib
import json
import os

import duckdb
import numpy as np
import pandas as pd


def wy_of(s):
    d = pd.to_datetime(s)
    return np.where(d.dt.month >= 10, d.dt.year + 1, d.dt.year)


def spaced(ids, n):
    ids = sorted(ids)
    if n >= len(ids):
        return ids
    pick = np.linspace(0, len(ids) - 1, n).round().astype(int)
    return [ids[i] for i in sorted(set(pick))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", required=True); ap.add_argument("--selection", required=True)
    ap.add_argument("--ledger", required=True); ap.add_argument("--inserted", required=True)
    ap.add_argument("--old-annual", required=True); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--control", type=int, default=200); ap.add_argument("--join-only", type=int, default=200)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    w0, w1 = (int(x) for x in a.window.split("-"))
    cohort = [l.strip() for l in open(a.selection) if l.strip()]; cset = set(cohort)
    led = pd.read_csv(a.ledger, dtype={"gage_id": str}); led = led[led["status"] == "restored"].copy()
    led["wy"] = wy_of(led["date"]); led = led[(led["wy"] >= w0) & (led["wy"] <= w1) & led["gage_id"].isin(cset)]
    restored_years = {g: set(d["wy"]) for g, d in led.groupby("gage_id")}
    inv = pd.read_csv(a.inserted, dtype={"gage_id": str}); inv = inv[(inv["water_year"] >= w0) & (inv["water_year"] <= w1) & inv["gage_id"].isin(cset)]
    con = duckdb.connect()
    vy = con.execute(f"SELECT CAST(gage_id AS VARCHAR) g, CAST(water_year AS INTEGER) wy, bool_or(signature = 'prcp_annual') has_clim "
                     f"FROM read_parquet('{a.old_annual}') WHERE water_year BETWEEN {w0} AND {w1} GROUP BY 1, 2").df()
    qvalid = {}; climvalid = {}
    for g, wy, hc in zip(vy["g"], vy["wy"], vy["has_clim"]):
        qvalid.setdefault(g, set()).add(int(wy))
        if hc:
            climvalid.setdefault(g, set()).add(int(wy))
    ins_years = {g: set(d["water_year"].astype(int)) for g, d in inv[inv["kind"] == "inserted"].groupby("gage_id")}
    nb_years = {g: set(d["water_year"].astype(int)) for g, d in inv[inv["kind"] != "inserted"].groupby("gage_id")}
    AQ = set(restored_years)
    J = {g for g, ys in ins_years.items() if ys & qvalid.get(g, set())}
    overlap = AQ & J
    restored_with_ins = {g for g in AQ if restored_years[g] & ins_years.get(g, set())}
    elig = {g for g, ys in ins_years.items() if (ys & qvalid.get(g, set())) - climvalid.get(g, set())}
    dec31 = {g for g, ys in nb_years.items() if ys & qvalid.get(g, set())}
    U = cset - AQ - J
    control = spaced(U, a.control)
    join_only = spaced(J - AQ, a.join_only)
    s1 = sorted(AQ) + control
    s1p = sorted(overlap | restored_with_ins | elig | dec31 | set(join_only))
    def write(name, ids):
        p = os.path.join(a.out_dir, f"{name}_{a.window}.txt"); open(p, "w").write("\n".join(ids) + "\n")
        return {"path": p, "n": len(ids), "sha256": hashlib.sha256(open(p, "rb").read()).hexdigest()}
    out = {"window": a.window, "cohort": len(cset), "A-Q": len(AQ), "J": len(J), "overlap_AQ_J": len(overlap),
           "AQ_restored_years_with_inserted_rows": len(restored_with_ins), "eligibility_gain_gages": sorted(elig),
           "dec31_neighbour_gages": len(dec31), "U": len(U), "control_sample": len(control), "join_only_sample": len(join_only),
           "S1": write("stage_S1", s1), "S1prime": write("stage_S1prime", s1p),
           "S1prime_members": {"overlap": sorted(overlap), "restored_with_inserted": sorted(restored_with_ins - overlap),
                               "eligibility": sorted(elig - overlap), "dec31": sorted(dec31 - overlap - elig)}}
    json.dump(out, open(os.path.join(a.out_dir, f"stage_cohorts_{a.window}.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "S1prime_members"}, indent=1))


if __name__ == "__main__":
    main()
