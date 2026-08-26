#!/usr/bin/env python3
"""
Gate a benchmark run log for per-gage signature-family failures.

WHY THIS EXISTS
---------------
`calculate_all_signatures()` wraps every signature family in a try/catch in all
three languages, so an exception on one family for one gage degrades to a
warning: that gage simply loses those columns. The final CSV is rectangular
(the writer unions keys across gages), so the loss reappears as ordinary NA and
is indistinguishable from a legitimately non-computable value. A green unit
suite proves nothing about it, and the column-set gate cannot see it either.

This is not hypothetical. The August 2026 port campaign hit the same failure
CLASS twice, both times as clean warning-only degradation:
  * rpkg's runner passed none of the new arguments, so a 19-hour run could not
    have produced the ported columns;
  * rpkg's drought family returned an all-NA family for every real gage because
    the preprocessor renames `date` -> `Date`.

So: treat ANY per-gage family failure as a run-level failure unless a human has
explicitly waived it.

WHAT IT READS
-------------
The runners emit a recognizable line per failure:
  Python : "Gage <id>: <family> signatures failed: <ExcType>: <msg>"
  rpkg   : "Gage <id>: <family> failed: <msg>"
  Julia  : "<family> failed for gage <id>"   (@warn, message on its own line)

NOTE ON R: R defers warnings and truncates at 50 ("There were 50 or more
warnings"), so a log from a runner that does NOT set `options(warn = 1)` can
hide failures. This tool FAILS LOUDLY when it sees that truncation banner
rather than reporting a clean pass it cannot substantiate.

USAGE
-----
    python docs/benchmarks/check_signature_failures.py RUN.log [--allow-family F]...
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

# Python/rpkg: "Gage 01013500: recession signatures failed: ValueError: ..."
RE_GAGE_FIRST = re.compile(r"Gage\s+(\S+?):\s+(.+?)\s+(?:signatures\s+)?failed(?::|\b)")
# Julia: "analyze_recession_parameters failed for gage 01013500"
RE_JULIA = re.compile(r"(\S+)\s+failed\s+for\s+gage\s+(\S+)")
# R's deferred-warning truncation banner — makes a clean result unprovable.
RE_R_TRUNCATED = re.compile(r"There were (\d+|50 or more) warnings")


def scan(path: Path):
    failures = []  # (gage, family)
    truncated = False
    with path.open(errors="replace") as fh:
        for line in fh:
            if RE_R_TRUNCATED.search(line):
                truncated = True
            m = RE_GAGE_FIRST.search(line)
            if m:
                failures.append((m.group(1), m.group(2).strip()))
                continue
            m = RE_JULIA.search(line)
            if m:
                failures.append((m.group(2), m.group(1).strip()))
    return failures, truncated


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path)
    ap.add_argument("--allow-family", action="append", default=[],
                    metavar="FAMILY",
                    help="Waive failures for this family (repeatable). Each waiver is "
                         "printed, so a waived run still reports what was ignored.")
    ap.add_argument("--max-show", type=int, default=15)
    args = ap.parse_args()

    if not args.log.exists():
        print(f"FAIL: log not found: {args.log}")
        return 2

    failures, truncated = scan(args.log)
    waived = {f.lower() for f in args.allow_family}

    hard = [(g, f) for g, f in failures if f.lower() not in waived]
    soft = [(g, f) for g, f in failures if f.lower() in waived]

    print(f"Signature-failure gate: {args.log}")
    print(f"  per-gage family failures found : {len(failures)}")
    if soft:
        by_fam = Counter(f for _, f in soft)
        print(f"  WAIVED (--allow-family)        : {len(soft)}")
        for fam, n in by_fam.most_common():
            print(f"      {fam}: {n}")

    if truncated:
        print("\nFAIL: the log contains R's deferred-warning truncation banner "
              "('There were N warnings'), so per-gage failures may be hidden.")
        print("      Re-run the R benchmark with options(warn = 1) so every warning "
              "prints immediately; a pass cannot be substantiated from this log.")
        return 1

    if not hard:
        print("\nPASS: no unwaived per-gage signature-family failures.")
        return 0

    by_fam = Counter(f for _, f in hard)
    by_gage = Counter(g for g, _ in hard)
    print(f"\nFAIL: {len(hard)} unwaived failures across "
          f"{len(by_gage)} gages and {len(by_fam)} families.")
    print("\n  by family:")
    for fam, n in by_fam.most_common():
        print(f"      {fam}: {n}")
    print("\n  sample (gage, family):")
    for g, f in hard[:args.max_show]:
        print(f"      {g}  {f}")
    if len(hard) > args.max_show:
        print(f"      ... and {len(hard) - args.max_show} more")
    print("\n  These gages silently lost those families' columns to NA in the "
          "output CSV. Fix the cause, or waive explicitly with --allow-family.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
