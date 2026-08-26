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

All three patterns were verified empirically (2026-08-26) against real emitted
output, not just read off the source: Python family failures were induced by
feeding the orchestrator a malformed frame, and all 8 resulting lines were
detected. That test also confirmed a real hazard is absent — the Python runner
never configures logging, but `logging.lastResort` still writes WARNING records
to stderr, so the messages do reach the log. (Had they not, this gate would
have reported a clean pass on a blind log.)

NOTE ON R: R defers warnings and truncates at 50 ("There were 50 or more
warnings"), so a log from a runner that does NOT set `options(warn = 1)` can
hide failures. This tool FAILS LOUDLY when it sees that truncation banner
rather than reporting a clean pass it cannot substantiate.

USAGE
-----
    python docs/benchmarks/check_signature_failures.py RUN.log \
        [--expect-prefix PREFIX] [--require-error-tally] [--allow-family F]...

--expect-prefix binds the log to the run whose artifacts are being gated: without
it an operator can gate a broken run's outputs while this gate certifies a
different run's clean log. --require-error-tally additionally demands the
"Gages errored: N" line the rpkg runner prints AFTER "BENCHMARK COMPLETE",
catching a footer truncated between the two. --expect-gages/--expect-columns
require the log's own footer counts to match the artifacts under test.

LIMIT OF THAT BINDING — state it rather than overclaim. Prefix and counts
identify a run's CONFIGURATION and its output SHAPE, not the run itself: two
executions of the same configuration produce the same prefix and the same
counts, so an earlier clean log can still be paired with a later run's artifacts
if someone supplies it. Together with the single-footer rule and the harness's
refusal to guess a log path, this makes ACCIDENTAL mispairing hard; it does not
make a deliberately or carelessly supplied stale log detectable.

The durable fix is a per-run identifier emitted by the runner into BOTH the log
header and the timing JSON, which the gate then requires to match — unforgeable
by copying or touching a file, and unique per execution. That is a runner change
and is queued with the other deferred rpkg fixes (CHANGELOG -> Planned).
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

# A log must PROVE it observed a completed benchmark. Without this, a zero-byte
# log, a stdout-only capture that dropped stderr, or a log cut off before
# processing started all contain no failure lines and would "pass".
#
# This is the ONLY true end-of-run marker. Earlier drafts accepted mid-run lines such as
# rpkg's "Signatures processed in" or "Writing annual values parquet", but those
# all precede output assembly, the metadata merge, and the CSV/timing writes — a
# log truncated just after one of them would have been certified (Codex delta
# review, 2026-08-26). All three runners print this AFTER the timing JSON.
RE_COMPLETE = re.compile(r"BENCHMARK COMPLETE")
# The rpkg footer reports the per-gage error tally — AFTER "BENCHMARK COMPLETE",
# so a log truncated between the two lines has a completion marker and no tally.
RE_ERRORED = re.compile(r"Gages errored:\s*(\d+)")
# Footer counts. These bind a log to an artifact SET by content rather than by
# name or timestamp: a stale log from a different run carries that run's numbers,
# and unlike an mtime they cannot be altered by copying or touching the file.
RE_N_GAGES = re.compile(r"Gages processed:\s*(\d+)")
RE_N_COLS = re.compile(r"Total columns:\s*(\d+)")
# Every runner prints this near the top; it binds a log to the run that produced it.
RE_PREFIX = re.compile(r"OUTPUT_PREFIX=\s*([^\s,]+)")


def scan(path: Path):
    failures = []  # (gage, family)
    truncated = False
    completed = False
    n_complete = 0
    errored_values: list[int] = []
    tally_after_complete = False
    prefixes: set[str] = set()
    footer_gages: list[int] = []
    footer_cols: list[int] = []
    n_lines = 0
    with path.open(errors="replace") as fh:
        for line in fh:
            n_lines += 1
            if RE_COMPLETE.search(line):
                # COUNT them: a set of prefixes collapses duplicates, so two runs
                # sharing a prefix would otherwise look like one (Codex delta 3).
                n_complete += 1
                completed = True
            me = RE_ERRORED.search(line)
            if me:
                # Collect ALL tallies rather than letting the last one win: a
                # concatenated log with "Gages errored: 7" followed by
                # "Gages errored: 0" must not pass on the second value.
                errored_values.append(int(me.group(1)))
                if completed:
                    tally_after_complete = True
            mp = RE_PREFIX.search(line)
            if mp:
                prefixes.add(mp.group(1))
            mg = RE_N_GAGES.search(line)
            if mg:
                footer_gages.append(int(mg.group(1)))
            mc = RE_N_COLS.search(line)
            if mc:
                footer_cols.append(int(mc.group(1)))
            if RE_R_TRUNCATED.search(line):
                truncated = True
            m = RE_GAGE_FIRST.search(line)
            if m:
                failures.append((m.group(1), m.group(2).strip()))
                continue
            m = RE_JULIA.search(line)
            if m:
                failures.append((m.group(2), m.group(1).strip()))
    return (failures, truncated, completed, n_lines,
            errored_values, tally_after_complete, prefixes, n_complete,
            footer_gages, footer_cols)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path)
    ap.add_argument("--allow-family", action="append", default=[],
                    metavar="FAMILY",
                    help="Waive failures for this family (repeatable). Each waiver is "
                         "printed, so a waived run still reports what was ignored.")
    ap.add_argument("--max-show", type=int, default=15)
    ap.add_argument("--expect-prefix", metavar="PREFIX",
                    help="Require the log to record OUTPUT_PREFIX=PREFIX (and only "
                         "that one). Binds the log to the run whose artifacts are "
                         "being gated.")
    ap.add_argument("--expect-gages", type=int, metavar="N",
                    help="Require the log's footer to report N processed gages, "
                         "matching the artifacts under test.")
    ap.add_argument("--expect-columns", type=int, metavar="N",
                    help="Require the log's footer to report N total columns.")
    ap.add_argument("--require-error-tally", action="store_true",
                    help="Require a 'Gages errored: N' line after the completion "
                         "marker (rpkg runner). Catches a footer truncated between "
                         "'BENCHMARK COMPLETE' and the tally.")
    ap.add_argument("--skip-completion-check", action="store_true",
                    help="Proceed even without an end-of-run marker. Only for a log "
                         "whose completeness you have established another way.")
    args = ap.parse_args()

    if not args.log.exists():
        print(f"FAIL: log not found: {args.log}")
        return 2

    (failures, truncated, completed, n_lines,
     errored_values, tally_after_complete, prefixes,
     n_complete, footer_gages, footer_cols) = scan(args.log)
    waived = {f.lower() for f in args.allow_family}

    hard = [(g, f) for g, f in failures if f.lower() not in waived]
    soft = [(g, f) for g, f in failures if f.lower() in waived]

    print(f"Signature-failure gate: {args.log}")
    print(f"  log lines read                 : {n_lines:,}")
    print(f"  BENCHMARK COMPLETE footer      : {'found' if completed else 'NOT FOUND'}")
    print(f"  log's footer counts            : "
          f"gages={footer_gages or 'n/a'} columns={footer_cols or 'n/a'}")
    print(f"  log's OUTPUT_PREFIX            : "
          f"{', '.join(sorted(prefixes)) if prefixes else 'not recorded'}")
    print(f"  runner-reported gage errors    : "
          f"{errored_values if errored_values else 'not reported'}")
    print(f"  per-gage family failures found : {len(failures)}")

    # Bind the log to the run under test. Without this an operator can gate one
    # run's artifacts while this gate certifies a different run's clean log.
    if args.expect_prefix:
        if args.expect_prefix not in prefixes:
            print(f"\nFAIL: this log does not belong to '{args.expect_prefix}'. "
                  f"It records OUTPUT_PREFIX="
                  f"{', '.join(sorted(prefixes)) if prefixes else '(none)'}.")
            return 1
        if n_complete > 1:
            print(f"\nFAIL: this log contains {n_complete} 'BENCHMARK COMPLETE' "
                  f"footers, so it covers more than one run. A clean stretch in it "
                  f"cannot vouch for the artifacts under test — even at the same "
                  f"OUTPUT_PREFIX, since a re-run overwrites them.")
            return 1
        if len(prefixes) > 1:
            print(f"\nFAIL: this log records more than one OUTPUT_PREFIX "
                  f"({', '.join(sorted(prefixes))}) — it covers multiple runs, so a "
                  f"clean stretch cannot vouch for the run under test.")
            return 1

    # Content binding. Name and mtime can both be faked by copying or touching a
    # file; the counts a run printed in its own footer cannot (Codex delta 4).
    for label, want, got in (("gages", args.expect_gages, footer_gages),
                             ("columns", args.expect_columns, footer_cols)):
        if want is None:
            continue
        if not got:
            print(f"\nFAIL: the log reports no {label} count, so it cannot be shown "
                  f"to describe the artifacts under test.")
            return 1
        if got[-1] != want:
            print(f"\nFAIL: the log reports {got[-1]} {label} but the artifacts under "
                  f"test have {want}. This log belongs to a different run.")
            return 1

    if any(v > 0 for v in errored_values):
        print(f"\nFAIL: the runner itself reported errored gages {errored_values}; "
              f"those gages are missing from the output entirely.")
        return 1

    if args.require_error_tally and not args.skip_completion_check:
        if not errored_values:
            print("\nFAIL: no 'Gages errored: N' tally in this log. The rpkg runner "
                  "prints it AFTER 'BENCHMARK COMPLETE', so its absence means the "
                  "footer is truncated and the run's error count is unknown.")
            return 1
        if completed and not tally_after_complete:
            print("\nFAIL: the error tally appears only BEFORE 'BENCHMARK COMPLETE', "
                  "so this log's footer is not the one belonging to the completed "
                  "run.")
            return 1

    if not completed and not args.skip_completion_check:
        print("\nFAIL: no end-of-run marker in this log, so it cannot be shown to "
              "cover a completed benchmark.")
        print("      An absent failure line proves nothing in a log that is empty, "
              "truncated, still running, or captured without stderr.")
        print("      Re-check the path, ensure the runner's stderr was redirected "
              "(2>&1), or pass --skip-completion-check if you have verified "
              "completeness another way.")
        return 1
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
