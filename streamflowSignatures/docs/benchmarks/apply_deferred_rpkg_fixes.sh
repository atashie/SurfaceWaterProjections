#!/usr/bin/env bash
# Apply the rpkg-side fixes still deferred after 2026-08-26.
#
# THREE ARE ALREADY APPLIED to run_rpkg_benchmark.R (options(warn=1),
# report_interval 500->50, and the streamflow column projection) — the
# benchmark they were queued behind was stopped externally, which made them
# safe to apply. Their steps below are idempotent no-ops now. What REMAINS is
# rpkg package source: the STREAMFLOW_CONFIG env-var name, the na_handling
# fallbacks, the SWE-merged-only-under-PPT branch, and the run manifest.
#
# All four change rpkg source or its runner, which must not happen while a
# benchmark is executing — an `R CMD INSTALL` mid-run can break lazy-load access
# for the running process. Run this ONLY after the benchmark has finished and its
# outputs have been gated.
#
# Each fix was verified INERT for the 2026-08-26 run before being deferred, so
# applying them does not invalidate that run's results:
#   1. options(warn = 1)          — logging only
#   2. STREAMFLOW_CONFIG          — that run set no config env var
#   3. na_handling fallbacks      — that run's bundled na_handling is byte-identical
#                                   to canonical, so no fallback was reached
#   4. SWE merged outside PPT     — that run's log records PPT and SWE both present
#
# Usage:  bash docs/benchmarks/apply_deferred_rpkg_fixes.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

if pgrep -f run_rpkg_benchmark >/dev/null 2>&1; then
  echo "REFUSING: an rpkg benchmark is still running (pgrep matched)."
  echo "These fixes reinstall the package and must not be applied mid-run."
  exit 1
fi

python3 - <<'PYEOF'
from pathlib import Path

# --- 1 + 4: the benchmark runner -----------------------------------------------
p = Path("docs/benchmarks/run_rpkg_benchmark.R"); s = p.read_text()

# (1) Print each warning as it happens. R otherwise defers warnings to the end of
# the script and truncates at 50, so per-gage family failures can be invisible and
# check_signature_failures.py cannot certify the log.
if "options(warn = 1)" not in s:
    marker = "#!/usr/bin/env Rscript"
    ins = (marker + "\n\n# Print warnings as they occur. R defaults to deferring them to the end of the\n"
           "# script and truncating at 50, which can hide per-gage signature-family\n"
           "# failures entirely (2026-08-26 Codex MAJOR).\noptions(warn = 1)")
    s = s.replace(marker, ins, 1) if s.startswith(marker) else "options(warn = 1)\n" + s
    print("  [1] options(warn = 1) added")

# (5) Report progress far more often. R is slow enough here that a 500-gage
# interval put the first marker over an hour away, leaving the run unobservable.
# NOT a buffering problem — cat() was verified to flush promptly on both the
# local SSD and the exFAT volume; the loop simply had not reached iteration 500.
if "report_interval <- 500" in s:
    s = s.replace("report_interval <- 500", "report_interval <- 50", 1)
    print("  [5] report_interval 500 -> 50")

# (6) Column-project the streamflow read. rpkg's read_parquet() takes no
# col_select, so call arrow directly and apply the same Date -> date rename.
# VERIFIED SAFE (2026-08-26): the parquet carries
# gage_id, Date, Q, year, month, doy, water_year, dowy, flag; the runner uses
# only the first three, because it calls add_water_year_columns() immediately
# after subsetting, which re-derives water_year/month/dowy from date. Nothing in
# the runner reads flag, doy or year. Dropping 6 of 9 columns removes roughly
# 4 GB of resident data on a 111.6M-row table.
old_read = "streamflow <- read_parquet(parquet_path)"
if old_read in s:
    s = s.replace(old_read, (
        'streamflow <- data.table::as.data.table(arrow::read_parquet(\n'
        '  parquet_path, col_select = c("gage_id", "Date", "Q")))\n'
        'data.table::setnames(streamflow, "Date", "date")'), 1)
    print("  [6] streamflow read column-projected to gage_id, Date, Q")

print("  [4] REVIEW MANUALLY: move the SWE merge out of the `if (\"PPT\" %in% clim_cols)`")
print("      branch near line 250 so a climate source with SWE but no PPT still")
print("      yields snow metrics (Python merges the slice independently).")
p.write_text(s)

# --- 2 + 3: rpkg config --------------------------------------------------------
p = Path("rpkg/R/config.R"); s = p.read_text()

# (2) Honor the canonical env var name, keeping the rpkg-specific one as a fallback.
if 'STREAMFLOW_CONFIG' not in s:
    s = s.replace('Sys.getenv("STREAMFLOW_SIGNATURES_CONFIG"',
                  'Sys.getenv("STREAMFLOW_CONFIG", Sys.getenv("STREAMFLOW_SIGNATURES_CONFIG"',
                  1)
    print("  [2] STREAMFLOW_CONFIG honored (verify the paren balance by hand)")

print("  [3] REVIEW MANUALLY: when the whole `na_handling` section is absent, rpkg")
print("      must default reject_negative_flow to FALSE (not TRUE) and must still")
print("      initialize na_max_raw_na_swe / na_max_gap_swe / na_reject_negative_swe,")
print("      matching julia/src/config.jl and python/streamflow_signatures/config.py.")
p.write_text(s)
PYEOF

echo
echo "Two fixes are marked REVIEW MANUALLY on purpose: both need judgement about"
echo "matching canonical semantics exactly, and a blind sed would risk diverging"
echo "from Julia in a way the unit suite would not catch."
echo
echo "After editing:  R CMD INSTALL rpkg && Rscript -e 'testthat::test_dir(\"rpkg/tests/testthat\")'"
