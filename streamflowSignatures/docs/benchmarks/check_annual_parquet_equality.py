#!/usr/bin/env python3
"""
Cross-language gate for the annual-values parquet: a port vs the Julia reference.

WHY A SEPARATE TOOL
-------------------
`validate_annual_values.py` checks a run against ITSELF — it recomputes mean/median
from the parquet and cross-checks that run's own summary CSV. That is necessary but
NOT sufficient: a port whose annual series are wrong in the same way in both artifacts
passes it cleanly (Codex finding F6). This tool is the sufficient check — it compares
the port's annual parquet against canonical Julia's, which is the only artifact that
can prove the per-year series themselves agree.

GATES (all must pass)
---------------------
  1. Signature-set equality        — exact, modulo explicit --waive-prefix
  2. Key-set equality              — (gage_id, signature, water_year), exact
  3. No duplicate keys             — in either side
  4. NA-pattern equality           — NaN-vs-finite must agree cell for cell
  5. Value agreement               — |diff| <= --tol on finite pairs (NaN==NaN ok)

Values are compared with a tolerance because different statistics libraries
legitimately differ in the last bits; the schema, key, duplicate and NA-pattern
gates are EXACT, because none of them has a legitimate reason to differ.

USAGE
-----
    python docs/benchmarks/check_annual_parquet_equality.py \
        --ref  REFERENCE_signatures_annual.parquet \
        --port PORT_signatures_annual.parquet \
        [--tol 1e-9] [--waive-prefix drought_ ...] [--max-show 20] \
        [--allow-key-diff N] [--allow-value-diff N]

QUOTE THE TOLERANCE WITH ANY COUNT THIS PRINTS. The same run reads differently at
different thresholds: Python vs Julia gives 267 over-tolerance rows at 1e-6 (all
discrete threshold-crossing metrics, where a last-bit tie moves a whole day) and
465 at 1e-9 (the extra 198 being qp_slope_sd / elasticity_annual at ~1e-8, ordinary
cross-library floating-point noise). Neither number means anything on its own.

--allow-key-diff / --allow-value-diff exist to record a residual that has already
been characterised in writing, NOT to quiet an unexplained one: the counts, the
families involved and the worst offending rows are printed either way.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ["gage_id", "signature", "water_year"]


def load(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_parquet(path)
    missing = [c for c in KEY + ["value"] if c not in df.columns]
    if missing:
        raise SystemExit(f"FAIL: {label} parquet missing columns: {missing}")
    df["gage_id"] = df["gage_id"].astype(str)
    df["signature"] = df["signature"].astype(str)
    df["water_year"] = df["water_year"].astype("int64")
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", type=Path, required=True, help="Julia reference parquet")
    ap.add_argument("--port", type=Path, required=True, help="Port's parquet")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--waive-prefix", action="append", default=[], metavar="PREFIX",
                    help="Waive reference signatures starting with PREFIX (repeatable). "
                         "Every waiver is printed, so a waived run still reports what "
                         "was excluded from the gate.")
    ap.add_argument("--max-show", type=int, default=20)
    ap.add_argument("--allow-key-diff", type=int, default=0, metavar="N",
                    help="Permit up to N rows present on only one side. Use ONLY for a "
                         "residual that has been characterised and written down; the "
                         "actual count and its signatures are always printed.")
    ap.add_argument("--allow-diff-signature", action="append", default=[],
                    metavar="SIG",
                    help="Restrict --allow-key-diff/--allow-value-diff to these "
                         "signatures (repeatable). A difference in ANY other signature "
                         "then fails regardless of count — so an allowance names the "
                         "residual it excuses instead of merely bounding it.")
    ap.add_argument("--allow-value-diff", type=int, default=0, metavar="N",
                    help="Permit up to N finite pairs over --tol. Same discipline: the "
                         "count, the families and the worst rows are always printed, so "
                         "an allowance documents a residual rather than hiding it.")
    args = ap.parse_args()

    # A bare count allowance is sized for one residual but silently excuses any
    # other (Codex delta review, 2026-08-26). Require the residual to be named.
    if (args.allow_key_diff or args.allow_value_diff) and not args.allow_diff_signature:
        print("FAIL: --allow-key-diff/--allow-value-diff require at least one "
              "--allow-diff-signature naming the residual being excused. A bare "
              "count would also excuse an unrelated difference of the same size.")
        return 2

    ref = load(args.ref, "reference")
    port = load(args.port, "port")
    failures: list[str] = []

    # Emptiness check FIRST. Without it two structurally-valid zero-row parquets
    # satisfy every downstream check vacuously and the gate prints PASS.
    if len(ref) == 0 or len(port) == 0:
        print(f"\nFAIL: empty artifact — reference has {len(ref):,} rows, "
              f"port has {len(port):,}. A gate cannot certify nothing.")
        return 1

    print(f"reference : {args.ref.name}  ({len(ref):,} rows)")
    print(f"port      : {args.port.name}  ({len(port):,} rows)")
    if args.waive_prefix:
        print(f"waived prefixes: {', '.join(args.waive_prefix)}")

    # ---- 1. signature-set equality -------------------------------------------
    r_sigs, p_sigs = set(ref.signature.unique()), set(port.signature.unique())
    waived = {s for s in r_sigs
              if any(s.startswith(w) for w in args.waive_prefix)}
    if waived:
        print(f"  waived signatures ({len(waived)}): "
              f"{', '.join(sorted(waived)[:8])}{' ...' if len(waived) > 8 else ''}")
        r_sigs -= waived
        ref = ref[~ref.signature.isin(waived)]
        port = port[~port.signature.isin(waived)]
        p_sigs -= waived

    only_ref, only_port = sorted(r_sigs - p_sigs), sorted(p_sigs - r_sigs)
    print(f"\n[1] signature sets: ref={len(r_sigs)} port={len(p_sigs)}")
    if only_ref or only_port:
        failures.append("signature-set mismatch")
        if only_ref:
            print(f"    REFERENCE-ONLY ({len(only_ref)}): {only_ref[:args.max_show]}")
        if only_port:
            print(f"    PORT-ONLY ({len(only_port)}): {only_port[:args.max_show]}")
    else:
        print("    OK — identical")

    # ---- 3. duplicate keys (checked before the merge, which they would corrupt)
    r_dup, p_dup = int(ref.duplicated(KEY).sum()), int(port.duplicated(KEY).sum())
    print(f"\n[3] duplicate keys: ref={r_dup} port={p_dup}")
    if r_dup or p_dup:
        failures.append("duplicate keys")
        print("    FAIL — duplicate keys make the row-level comparison ill-defined")
    else:
        print("    OK — none")

    # ---- 2. key-set equality --------------------------------------------------
    m = ref.merge(port, on=KEY, how="outer", suffixes=("_ref", "_port"), indicator=True)
    n_ref_only = int((m._merge == "left_only").sum())
    n_port_only = int((m._merge == "right_only").sum())
    n_both = int((m._merge == "both").sum())
    print(f"\n[2] key sets: shared={n_both:,} ref-only={n_ref_only:,} "
          f"port-only={n_port_only:,}")
    n_key_diff = n_ref_only + n_port_only
    unscoped_keys = set()
    if n_key_diff and args.allow_diff_signature:
        off = m[m._merge != "both"]
        unscoped_keys = set(off.signature.unique()) - set(args.allow_diff_signature)
    if n_key_diff > args.allow_key_diff or unscoped_keys:
        failures.append("key-set mismatch")
        if unscoped_keys:
            print(f"    NOT COVERED by --allow-diff-signature: "
                  f"{sorted(unscoped_keys)[:args.max_show or 10]}")
    elif n_key_diff:
        print(f"    ALLOWED: {n_key_diff} <= --allow-key-diff {args.allow_key_diff} "
              f"(documented residual)")
    if n_key_diff:
        for side, lbl in (("left_only", "REFERENCE-ONLY"), ("right_only", "PORT-ONLY")):
            sub = m[m._merge == side]
            if len(sub):
                by_sig = sub.signature.value_counts()
                print(f"    {lbl} rows by signature (top): "
                      f"{dict(by_sig.head(8))}")
    else:
        print("    OK — identical")

    both = m[m._merge == "both"]
    if len(both) == 0:
        print("\nFAIL: the two parquets share no keys at all — nothing to compare.")
        return 1
    vr = both.value_ref.to_numpy(dtype=float)
    vp = both.value_port.to_numpy(dtype=float)

    # ---- 4. NA-pattern equality ----------------------------------------------
    # NaN means "not computable that year" and is the only missing marker.
    # +/-Inf are real values: previously ~isfinite lumped them in with NaN, so a
    # reference +Inf against a port -Inf at the same key compared as "both NA".
    na_r, na_p = np.isnan(vr), np.isnan(vp)
    n_na_mismatch = int((na_r != na_p).sum())
    print(f"\n[4] NA patterns: mismatches={n_na_mismatch:,}")
    if n_na_mismatch:
        failures.append("NA-pattern mismatch")
        bad = both[na_r != na_p]
        print(f"    by signature (top): {dict(bad.signature.value_counts().head(8))}")
    else:
        print("    OK — identical")

    # ---- 5. value agreement ---------------------------------------------------
    present = ~na_r & ~na_p
    inf_pair = present & (np.isinf(vr) | np.isinf(vp))
    n_inf_bad = int((vr[inf_pair] != vp[inf_pair]).sum()) if inf_pair.any() else 0
    if inf_pair.any():
        print(f"\n[5a] infinities: {int(inf_pair.sum()):,} pairs where either side is "
              f"+/-Inf; sign/'value' mismatches={n_inf_bad:,}")
        if n_inf_bad:
            failures.append("infinity mismatch")
            print("    FAIL — an infinity must match exactly (same sign) to agree")
        else:
            print("    OK — every infinity matches exactly")

    finite = present & np.isfinite(vr) & np.isfinite(vp)
    if finite.any():
        diff = np.abs(vr[finite] - vp[finite])
        n_over = int((diff > args.tol).sum())
        print(f"\n[5] values (finite pairs={int(finite.sum()):,}, tol={args.tol:g}): "
              f"over-tolerance={n_over:,}  max|diff|={diff.max():.3e}")
        unscoped_vals = set()
        if n_over and args.allow_diff_signature:
            sigs_over = both[finite].signature.to_numpy()[diff > args.tol]
            unscoped_vals = set(np.unique(sigs_over)) - set(args.allow_diff_signature)
        if n_over > args.allow_value_diff or unscoped_vals:
            failures.append("value mismatch")
            if unscoped_vals:
                print(f"    NOT COVERED by --allow-diff-signature: "
                      f"{sorted(unscoped_vals)[:10]}")
        elif n_over:
            print(f"    ALLOWED: {n_over} <= --allow-value-diff "
                  f"{args.allow_value_diff} (documented residual)")
        if n_over:
            bad = both[finite].assign(d=diff)
            bad = bad[bad.d > args.tol].sort_values("d", ascending=False)
            print(f"    by signature (top): {dict(bad.signature.value_counts().head(8))}")
            for _, r in bad.head(args.max_show).iterrows():
                print(f"      {r.gage_id} {r.signature} {int(r.water_year)}: "
                      f"ref={r.value_ref!r} port={r.value_port!r} d={r.d:.3e}")
        else:
            print("    OK — all within tolerance")
    else:
        print("\n[5] values: no finite pairs to compare")

    print("\n" + "=" * 62)
    if failures:
        print(f"GATE FAILED: {', '.join(failures)}")
        return 1
    print("GATE PASSED: annual parquets agree "
          "(signatures, keys, duplicates, NA patterns, values).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
