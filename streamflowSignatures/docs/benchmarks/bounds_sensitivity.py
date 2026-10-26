#!/usr/bin/env python3
"""bounds_sensitivity.py — decision 1b of the qualifier rerun plan: the censored `A <` / `A >` days are
restored with the reported BOUND as the value; this run quantifies what that choice does, on
exactly the gages that hold such days, by running the canonical Julia runner twice on the patched
input — bounds as values (the product policy) and bounds as NA (the alternative) — and comparing
qualification (which gages and water years qualify) and every signature (the change census).

    python docs/benchmarks/bounds_sensitivity.py --patched PATCHED.parquet --ledger correction_ledger.csv
        --climate-batches "~/HISSS_data/climate-rerun-12oct2026/inputs/climate_daymet2025_b*.parquet"
        --wrapper docs/benchmarks/run_julia_benchmark_prod_1993_2025_60pct_climate.jl --window 1993-2025
        --qa-schema V2.csv --metadata META_PATCHED.csv --gages-ii-dir DIR --out-dir DIR [--julia ~/.juliaup/bin/julia]

Writes DIR/bounds_value/ and DIR/bounds_na/ (the two runs), DIR/qualification_sensitivity.csv (per
gage: qualifies and num_water_years under each policy), DIR/sensitivity_census/ (assess_rerun_changes
on the shared gages with the bound days as the ledger) and DIR/bounds_sensitivity.md.
Needs duckdb + pyarrow + pandas (~/HISSS_data/venvs/daymet/bin/python).
"""
import argparse
import datetime as dt
import glob
import json
import os
import subprocess
import sys

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patched", required=True); ap.add_argument("--ledger", required=True); ap.add_argument("--climate-batches", required=True)
    ap.add_argument("--wrapper", required=True); ap.add_argument("--window", required=True); ap.add_argument("--qa-schema", required=True)
    ap.add_argument("--metadata", required=True); ap.add_argument("--gages-ii-dir", required=True); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--julia", default=os.path.expanduser("~/.juliaup/bin/julia")); ap.add_argument("--project", default="julia")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    led = pd.read_csv(a.ledger, dtype={"gage_id": str, "date": str}, float_precision="round_trip")
    bounds = led[(led["status"] == "restored") & led["old_flag"].isin(["A <", "A >"])].copy()
    bounds["date"] = pd.to_datetime(bounds["date"]).dt.date
    gages = sorted(bounds["gage_id"].unique())
    keys = set(zip(bounds["gage_id"], bounds["date"]))
    log(f"{len(bounds):,} censored-bound days on {len(gages)} gages ({bounds['old_flag'].value_counts().to_dict()})")
    bounds.to_csv(os.path.join(a.out_dir, "bound_days_ledger.csv"), index=False, float_format="%.17g")
    # the two streamflow inputs (the gages only; identical rows except Q on the bound days)
    f = pq.ParquetFile(a.patched); schema = f.schema_arrow
    wv = pq.ParquetWriter(os.path.join(a.out_dir, "bounds_value.parquet"), schema, compression="zstd")
    wn = pq.ParquetWriter(os.path.join(a.out_dir, "bounds_na.parquet"), schema, compression="zstd")
    n_rows = n_masked = 0
    gset = set(gages)
    for i in range(f.num_row_groups):
        t = f.read_row_group(i)
        gid = np.asarray(t.column("gage_id").to_pylist(), dtype=object)
        hit = np.isin(gid, list(gset))
        if not hit.any():
            continue
        t = t.filter(pa.array(hit)); n_rows += t.num_rows
        wv.write_table(t)
        q = t.column("Q").to_pandas().to_numpy(dtype=float).copy(); dts = t.column("Date").to_pandas().to_numpy(); g2 = t.column("gage_id").to_pylist()
        for j in range(len(g2)):
            if (g2[j], dts[j]) in keys:
                q[j] = np.nan; n_masked += 1
        t = t.set_column(t.schema.get_field_index("Q"), "Q", pa.array(q, type=schema.field("Q").type, mask=np.isnan(q)))
        wn.write_table(t)
    wv.close(); wn.close()
    if n_masked != len(keys):
        sys.exit(f"masked {n_masked} bound days but the ledger lists {len(keys)}")
    log(f"streamflow subsets written: {n_rows:,} rows, {n_masked:,} bound days set to NA in bounds_na.parquet")
    con = duckdb.connect(); con.execute("SET memory_limit='2GB'")
    files = sorted(glob.glob(os.path.expanduser(a.climate_batches)))
    con.execute("CREATE TEMP TABLE gs AS SELECT unnest(?::VARCHAR[]) AS site_id", [gages])
    n_cl = con.execute(f"COPY (SELECT * FROM read_parquet({files}) WHERE site_id IN (SELECT site_id FROM gs) ORDER BY site_id, Date) TO '{os.path.join(a.out_dir, 'climate_sub.parquet')}' (FORMAT PARQUET, COMPRESSION ZSTD)").fetchone()
    log(f"climate subset written from {len(files)} batch files")
    runs = {}
    for name in ("bounds_value", "bounds_na"):
        out = os.path.join(a.out_dir, name); os.makedirs(out, exist_ok=True)
        env = dict(os.environ, STREAMFLOW_DATA_PATH=os.path.join(a.out_dir, f"{name}.parquet"), STREAMFLOW_CLIMATE_PATH=os.path.join(a.out_dir, "climate_sub.parquet"),
                   STREAMFLOW_METADATA_PATH=a.metadata, STREAMFLOW_GAGES_II_DIR=a.gages_ii_dir, STREAMFLOW_OUTPUT_DIR=out, STREAMFLOW_OUTPUT_PREFIX=name,
                   STREAMFLOW_CLIMATE_FRAME="outer", STREAMFLOW_QA_SCHEMA=a.qa_schema, STREAMFLOW_HASH_INPUTS="1")
        log(f"running the Julia wrapper: {name} ...")
        with open(os.path.join(out, "run.log"), "w") as lg:
            rc = subprocess.run([a.julia, f"--project={a.project}", a.wrapper], env=env, stdout=lg, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            sys.exit(f"{name} run failed (rc {rc}); see {out}/run.log")
        runs[name] = {"csv": os.path.join(out, f"{name}_signatures.csv"), "annual": os.path.join(out, f"{name}_signatures_annual.parquet")}
    v = pd.read_csv(runs["bounds_value"]["csv"], dtype={"gage_id": str}, low_memory=False).set_index("gage_id")
    n_ = pd.read_csv(runs["bounds_na"]["csv"], dtype={"gage_id": str}, low_memory=False).set_index("gage_id")
    q = pd.DataFrame({"gage_id": gages})
    q["qualifies_bounds_value"] = q["gage_id"].isin(v.index); q["qualifies_bounds_na"] = q["gage_id"].isin(n_.index)
    q["num_water_years_value"] = q["gage_id"].map(v["num_water_years"]); q["num_water_years_na"] = q["gage_id"].map(n_["num_water_years"])
    q["n_bound_days"] = q["gage_id"].map(bounds.groupby("gage_id").size())
    q.to_csv(os.path.join(a.out_dir, "qualification_sensitivity.csv"), index=False)
    shared = sorted(set(v.index) & set(n_.index))
    open(os.path.join(a.out_dir, "shared_gages.txt"), "w").write("\n".join(shared) + "\n")
    cen = os.path.join(a.out_dir, "sensitivity_census")
    cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "assess_rerun_changes.py"),
           "--old-csv", runs["bounds_na"]["csv"], "--old-annual", runs["bounds_na"]["annual"], "--new-csv", runs["bounds_value"]["csv"],
           "--new-annual", runs["bounds_value"]["annual"], "--window", a.window, "--label", "bounds_value_vs_na", "--ledger",
           os.path.join(a.out_dir, "bound_days_ledger.csv"), "--gages", os.path.join(a.out_dir, "shared_gages.txt"), "--out-dir", cen, "--memory-limit", "2GB"]
    rc = subprocess.run(cmd).returncode
    gj = json.load(open(os.path.join(cen, f"gates_bounds_value_vs_na_{a.window}.json")))
    lines = [f"# Decision 1b sensitivity — censored bounds as values vs as NA ({a.window})", "",
             f"{len(bounds):,} `A <` / `A >` days on {len(gages)} gages (the patched input, Julia runner, outer frame, pinned schema).", "",
             f"Qualification: {int(q['qualifies_bounds_value'].sum())} gages qualify with the bounds as values, {int(q['qualifies_bounds_na'].sum())} with them as NA; "
             f"qualify ONLY with bounds as values: {sorted(q.loc[q.qualifies_bounds_value & ~q.qualifies_bounds_na, 'gage_id'])}; "
             f"only with NA: {sorted(q.loc[~q.qualifies_bounds_value & q.qualifies_bounds_na, 'gage_id'])}; "
             f"num_water_years differs on {int(((q.num_water_years_value != q.num_water_years_na) & q.qualifies_bounds_value & q.qualifies_bounds_na).sum())} of the "
             f"{int((q.qualifies_bounds_value & q.qualifies_bounds_na).sum())} gages that qualify under both policies (qualification_sensitivity.csv).", "",
             f"Signatures: change census over the {len(shared)} shared gages — annual {gj['annual']}, explanations {gj['explanations']}; "
             f"gates {gj['gates']} (verdict {gj['verdict']}: the census's oracle / E7 gates are NOT RUN here by design — this is a diagnostic of the policy, not a gate); "
             f"per-column magnitudes in sensitivity_census/summary_columns_bounds_value_vs_na_{a.window}.csv.", ""]
    open(os.path.join(a.out_dir, "bounds_sensitivity.md"), "w").write("\n".join(lines))
    log("bounds_sensitivity.md written")


if __name__ == "__main__":
    main()
