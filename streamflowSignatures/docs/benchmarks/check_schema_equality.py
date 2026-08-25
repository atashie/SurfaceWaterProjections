#!/usr/bin/env python3
"""
Strict schema-equality gate for cross-language signature outputs.

The comparison scripts (compare_*_vs_golden_julia.py) intentionally compare only
the INTERSECTION of columns and exit 0 even when whole families are missing from
one side — they are diagnostics, not gates (Codex F5, port plan §4.1). This tool
is the acceptance gate: it FAILS (nonzero exit) unless

  * the candidate's column set equals the reference's, and
  * the candidate's gage_id set equals the reference's,

modulo waivers that must be spelled out on the command line so every acceptance
log shows exactly what was excused:

  --allow-candidate-missing  columns that may be absent from the candidate
                             (e.g. features not yet ported at this phase)
  --allow-candidate-extra    columns that may exist only in the candidate
  --allow-prefix-missing     waive whole families by prefix (comma list), e.g.
                             'drought_,swe_' before those phases land

Optional hard counts pin the target schema so a waiver can't hide drift:
  --expect-reference-columns N   fail unless the reference has exactly N columns
  --expect-candidate-columns N   fail unless the candidate has exactly N columns

Usage:
  python check_schema_equality.py REFERENCE.csv CANDIDATE.csv [options]

Only the header row and the gage_id column are read — cheap on 1,653-column CSVs.
"""

import argparse
import csv
import sys


def read_header(path):
    with open(path, newline="", encoding="utf-8") as f:
        return next(csv.reader(f))


def read_gage_ids(path):
    ids = set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "gage_id" not in (reader.fieldnames or []):
            return None
        for row in reader:
            ids.add(row["gage_id"])
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("--allow-candidate-missing", default="",
                    help="comma-separated columns allowed to be absent from the candidate")
    ap.add_argument("--allow-candidate-extra", default="",
                    help="comma-separated columns allowed to exist only in the candidate")
    ap.add_argument("--allow-prefix-missing", default="",
                    help="comma-separated column-name prefixes allowed to be absent from the candidate")
    ap.add_argument("--baseline-mode", action="store_true",
                    help="pre-port baseline: report (never fail on) reference columns missing "
                         "from the candidate; still fail on candidate-only columns and gage-set "
                         "mismatches. NEVER valid for a phase-exit or final acceptance gate.")
    ap.add_argument("--expect-reference-columns", type=int, default=None)
    ap.add_argument("--expect-candidate-columns", type=int, default=None)
    ap.add_argument("--skip-gage-set", action="store_true",
                    help="skip the gage_id set-equality check (columns only)")
    args = ap.parse_args()

    allow_missing = {c for c in args.allow_candidate_missing.split(",") if c}
    allow_extra = {c for c in args.allow_candidate_extra.split(",") if c}
    allow_prefixes = tuple(p for p in args.allow_prefix_missing.split(",") if p)

    ref_cols = read_header(args.reference)
    cand_cols = read_header(args.candidate)
    ref_set, cand_set = set(ref_cols), set(cand_cols)

    failures = []
    print(f"reference: {args.reference} ({len(ref_cols)} columns)")
    print(f"candidate: {args.candidate} ({len(cand_cols)} columns)")

    if len(ref_set) != len(ref_cols):
        failures.append(f"reference has {len(ref_cols) - len(ref_set)} duplicate column names")
    if len(cand_set) != len(cand_cols):
        failures.append(f"candidate has {len(cand_cols) - len(cand_set)} duplicate column names")

    if args.expect_reference_columns is not None and len(ref_cols) != args.expect_reference_columns:
        failures.append(f"reference column count {len(ref_cols)} != expected {args.expect_reference_columns}")
    if args.expect_candidate_columns is not None and len(cand_cols) != args.expect_candidate_columns:
        failures.append(f"candidate column count {len(cand_cols)} != expected {args.expect_candidate_columns}")

    missing = sorted(ref_set - cand_set)
    unwaived_missing = [c for c in missing
                        if c not in allow_missing and not c.startswith(allow_prefixes)]
    waived_missing = [c for c in missing if c not in unwaived_missing]
    if waived_missing:
        print(f"  waived missing ({len(waived_missing)}): "
              f"{', '.join(waived_missing[:8])}{' …' if len(waived_missing) > 8 else ''}")
    if unwaived_missing:
        if args.baseline_mode:
            print(f"  BASELINE MODE: {len(unwaived_missing)} reference columns missing from "
                  f"candidate (reported, not failed): "
                  f"{', '.join(unwaived_missing[:10])}{' …' if len(unwaived_missing) > 10 else ''}")
        else:
            failures.append(f"{len(unwaived_missing)} reference columns absent from candidate: "
                            f"{', '.join(unwaived_missing[:20])}"
                            f"{' …' if len(unwaived_missing) > 20 else ''}")

    extra = sorted(cand_set - ref_set)
    unwaived_extra = [c for c in extra if c not in allow_extra]
    if extra and not unwaived_extra:
        print(f"  waived extra ({len(extra)}): {', '.join(extra[:8])}")
    if unwaived_extra:
        failures.append(f"{len(unwaived_extra)} candidate-only columns: "
                        f"{', '.join(unwaived_extra[:20])}"
                        f"{' …' if len(unwaived_extra) > 20 else ''}")

    if not args.skip_gage_set:
        ref_ids = read_gage_ids(args.reference)
        cand_ids = read_gage_ids(args.candidate)
        if ref_ids is None or cand_ids is None:
            failures.append("gage_id column missing from one of the files")
        else:
            only_ref = ref_ids - cand_ids
            only_cand = cand_ids - ref_ids
            print(f"gage sets: reference={len(ref_ids)} candidate={len(cand_ids)} "
                  f"only-ref={len(only_ref)} only-cand={len(only_cand)}")
            if only_ref:
                failures.append(f"{len(only_ref)} gages only in reference "
                                f"(e.g. {sorted(only_ref)[:5]})")
            if only_cand:
                failures.append(f"{len(only_cand)} gages only in candidate "
                                f"(e.g. {sorted(only_cand)[:5]})")

    print()
    if failures:
        print("SCHEMA GATE: FAIL")
        for f_ in failures:
            print(f"  - {f_}")
        return 1
    print("SCHEMA GATE: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
