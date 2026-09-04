#!/usr/bin/env python3
"""
Recompute ONLY the `flagged_for_high_na` column of a delivered signatures CSV,
byte-preserving everything else.

Why a text-level tool: the delivered products were written by CSV.jl and every
other value must stay byte-identical (the flag is a boolean; no statistic changes).
Re-serializing through pandas / CSV.jl would reformat floats. This script parses
each line with the csv module for the NA count, but rewrites the line by replacing
one field in the trailing flag block, so nothing else is touched.

Definition (config qa_qc.high_na_denominator, shared by Julia/Python/rpkg since
2026-09-04): fraction of the SIGNATURE columns present in the table (names ending
in one of the 16 statistic suffixes, plus the per-gage scalar registry) whose value
is NA/NaN, compared with max_na_fraction using a strict ">".

Usage:
    recompute_high_na_flag.py SIGNATURES.csv            # dry run: report only
    recompute_high_na_flag.py SIGNATURES.csv --write    # rewrite in place (keeps a .bak)
Options: --config PATH (default: config/signatures_config.json next to this repo),
         --no-backup
"""
import argparse
import csv
import json
import os
import shutil
import sys
import tempfile
from collections import Counter

NA_TOKENS = {"", "NA", "NaN", "nan", "NAN", "missing", "None"}
BOOL_STYLES = {"true": ("true", "false"), "True": ("True", "False"), "TRUE": ("TRUE", "FALSE")}


def denominator_columns(header, suffixes, scalars):
    scal = set(scalars)
    return [c for c in header
            if not c.startswith("flagged_") and (c in scal or any(c.endswith(s) for s in suffixes))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_path")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "..", "..", "config", "signatures_config.json"))
    a = ap.parse_args()

    cfg = json.load(open(a.config))["qa_qc"]
    thr = float(cfg["max_na_fraction"])
    man = cfg["high_na_denominator"]
    suffixes, scalars = list(man["statistic_suffixes"]), list(man["per_gage_scalars"])

    with open(a.csv_path, newline="") as f:
        header_line = f.readline()
    newline = "\r\n" if header_line.endswith("\r\n") else "\n"
    header = next(csv.reader([header_line]))
    if "flagged_for_high_na" not in header:
        sys.exit("no flagged_for_high_na column")
    flag_idx = header.index("flagged_for_high_na")
    n_trailing = len(header) - flag_idx          # fields from the flag to the end of the line
    trailing = header[flag_idx:]
    if not all(c.startswith("flagged_") for c in trailing):
        sys.exit("flag block is not trailing; refusing text-level rewrite")
    den_idx = [header.index(c) for c in denominator_columns(header, suffixes, scalars)]
    gt_idx = header.index("gage_type") if "gage_type" in header else None
    print(f"{os.path.basename(a.csv_path)}: {len(header)} columns; denominator = {len(den_idx)} signature columns; threshold > {thr}")

    old_true = new_true = 0
    flips = Counter(); by_type = Counter(); style = None
    out_lines = []
    with open(a.csv_path, newline="") as f:
        lines = f.read().split(newline)
    if lines and lines[-1] == "":
        lines.pop(); trailing_newline = True
    else:
        trailing_newline = False
    out_lines.append(lines[0])
    for line in lines[1:]:
        fields = next(csv.reader([line]))
        if len(fields) != len(header):
            sys.exit(f"row with {len(fields)} fields vs {len(header)} in header — aborting")
        old = fields[flag_idx]
        if style is None:
            for k, v in BOOL_STYLES.items():
                if old in v:
                    style = v
            if style is None:
                sys.exit(f"unrecognised boolean token {old!r}")
        na = sum(1 for i in den_idx if fields[i] in NA_TOKENS)
        new_flag = (na / len(den_idx)) > thr
        new = style[0] if new_flag else style[1]
        old_true += old == style[0]; new_true += new_flag
        if old != new:
            flips[(old, new)] += 1
            if gt_idx is not None:
                by_type[(fields[gt_idx], old, new)] += 1
        # text-level replacement inside the trailing flag block only: the last
        # n_trailing fields are booleans (no commas/quotes), so splitting on the
        # last n_trailing commas isolates them without re-serializing anything else
        parts = line.rsplit(",", n_trailing)
        head, tail = parts[0], parts[1:]
        if len(tail) != n_trailing or tail[0] != old:
            sys.exit("trailing flag block did not split cleanly — aborting")
        tail[0] = new
        out_lines.append(head + "," + ",".join(tail))
    print(f"  flagged_for_high_na TRUE: {old_true} -> {new_true} ({sum(flips.values())} gages change)")
    for (o, n), k in sorted(flips.items()):
        print(f"    {o} -> {n}: {k}")
    if by_type:
        for (t, o, n), k in sorted(by_type.items()):
            print(f"    {t}: {o} -> {n}: {k}")
    if not a.write:
        print("  dry run — file not modified (pass --write to apply)")
        return
    if not a.no_backup:
        bak = a.csv_path + ".bak_highna_pre_recompute"
        if not os.path.exists(bak):
            shutil.copy2(a.csv_path, bak)
            print(f"  backup: {bak}")
    tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=os.path.dirname(a.csv_path) or ".", newline="")
    with tmp:
        tmp.write(newline.join(out_lines) + (newline if trailing_newline else ""))
    os.replace(tmp.name, a.csv_path)
    print(f"  rewrote {a.csv_path}")


if __name__ == "__main__":
    main()
