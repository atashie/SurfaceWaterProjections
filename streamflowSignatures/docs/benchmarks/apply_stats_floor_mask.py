#!/usr/bin/env python
"""Post-hoc application of the stats floor (min_values_for_stats) to a delivered
summary CSV.

For every (gage, signature) whose annual parquet series has fewer than --floor
non-NaN values, sets ALL 16 statistic columns (8 stats + 8 Pettitt fields) to NaN
in the summary CSV — mathematically identical to a rerun with the floor active,
since the floor only suppresses outputs and never changes computed values.
Recession and elasticity bases are exempt (same as the in-code orchestration).

The original CSV is preserved as <name>_prefloor.csv. Note the rewritten CSV's
float TEXT formatting is pandas' (round-trip lossless numerically, but not
byte-identical to the Julia writer).

Usage:
    python apply_stats_floor_mask.py --csv <signatures.csv> --annual <annual.parquet>
        [--floor 20]
"""
import argparse
import os
import shutil

import numpy as np
import pandas as pd

STAT_SUFFIXES = ["_mean", "_median", "_senn_slp", "_linear_slp",
                 "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval"]
CP_SUFFIXES = ["_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
               "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
               "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
FLOOR_EXEMPT = {"log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events",
                "concavity", "n_recession_events", "alpha_linear",
                "elasticity_rolling", "elasticity_annual"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--annual", required=True)
    ap.add_argument("--floor", type=int, default=20)
    args = ap.parse_args()

    ann = pd.read_parquet(args.annual, columns=["gage_id", "signature", "value"])
    ann["gage_id"] = ann["gage_id"].astype(str)
    counts = (ann.dropna(subset=["value"])
                 .groupby(["gage_id", "signature"]).size())
    wide = counts.unstack(fill_value=0)  # gage x signature non-NaN counts

    # round_trip: pandas' default fast float parser is ~1 ulp off on some values
    # (same defect class the Daymet CSV converter fixed, CHANGELOG Aug 2026) —
    # without it the rewritten CSV perturbs last bits of unmasked values, which
    # breaks exact-equality gates built on this tool's output (found 2026-08-24
    # by the port campaign's Phase-1 gate).
    df = pd.read_csv(args.csv, dtype={"gage_id": str}, low_memory=False,
                     float_precision="round_trip")
    backup = os.path.splitext(args.csv)[0] + "_prefloor.csv"
    if not os.path.exists(backup):
        shutil.copy2(args.csv, backup)
        print(f"Backup: {backup}")
    else:
        print(f"Backup already exists (not overwritten): {backup}")

    bases = sorted({c[:-5] for c in df.columns
                    if c.endswith("_mean") and "_pettitt_" not in c})
    masked_cells = 0
    per_base = {}
    for b in bases:
        if b in FLOOR_EXEMPT or b not in wide.columns:
            continue
        n = wide[b].reindex(df["gage_id"]).fillna(0).to_numpy()
        under = n < args.floor
        if not under.any():
            continue
        cols = [b + s for s in STAT_SUFFIXES + CP_SUFFIXES if (b + s) in df.columns]
        # only count cells that actually change (were non-NaN)
        changed = 0
        for c in cols:
            vals = pd.to_numeric(df[c], errors="coerce")
            was = vals.notna().to_numpy() & under
            changed += int(was.sum())
            vals[under] = np.nan
            df[c] = vals
        if changed:
            per_base[b] = (int(under.sum()), changed)
            masked_cells += changed

    df.to_csv(args.csv, index=False)
    print(f"Floor {args.floor}: masked {masked_cells} previously non-NaN cells "
          f"across {len(per_base)} signatures in {os.path.basename(args.csv)}")
    for b, (n_g, n_c) in sorted(per_base.items(), key=lambda kv: -kv[1][1])[:20]:
        print(f"  {b:32s} {n_g:5d} gages under floor, {n_c:6d} cells masked")
    if len(per_base) > 20:
        print(f"  ... and {len(per_base) - 20} more signatures")


if __name__ == "__main__":
    main()
