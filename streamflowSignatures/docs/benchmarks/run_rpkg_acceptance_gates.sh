#!/usr/bin/env bash
# Acceptance gates for the feature-complete rpkg benchmark, against the Julia
# reference produced in Phase 0 of the August 2026 port campaign.
#
# Steps 1-3 are GATES: each exits nonzero on failure and drives the consolidated
# verdict. Step 4 is a DIAGNOSTIC and is reported separately — the
# compare_*_vs_golden_julia.py scripts compare only the column INTERSECTION and
# exit 0 even when whole families are missing, so their exit status must never be
# read as acceptance. Gate 1 is what proves the column and gage sets are equal.
#
# Usage:  bash docs/benchmarks/run_rpkg_acceptance_gates.sh [PREFIX] [RUN_LOG]
#           PREFIX  defaults to rpkg_full_26aug2026
#           RUN_LOG defaults to <PREFIX>_run.log. If that does not exist the script
#                   FAILS rather than guessing — it never falls back to a
#                   conventionally-named log that may belong to a different run.
#
# The log must be tied to the prefix under test. Hard-coding one log path lets the
# schema and annual gates examine a retry while the failure gate scans a clean log
# from an earlier run, so the harness can report success on a broken retry
# (Codex 2026-08-26 MAJOR).
set -uo pipefail

DIR=/Volumes/Untitled/processedOuts_portref_24aug2026
REF=$DIR/streamflow_1993_2025_60pct_portref_24aug2026
PREFIX=${1:-rpkg_full_26aug2026}
PORT=$DIR/$PREFIX
PY=${PY:-.venv/bin/python}

# Resolve the run log: explicit arg > <PREFIX>_run.log > fail. Never guess.
if [ "${2:-}" != "" ]; then
  LOG=$2
elif [ -f "$DIR/${PREFIX}_run.log" ]; then
  LOG=$DIR/${PREFIX}_run.log
else
  echo "FAIL: no $DIR/${PREFIX}_run.log, and no log given as argument 2."
  echo "      Refusing to guess. Falling back to a conventionally-named log would"
  echo "      let this harness gate one run's artifacts against another run's log,"
  echo "      which is exactly the failure this argument exists to prevent."
  echo "      Pass the log explicitly:"
  echo "        bash $0 $PREFIX /path/to/that/run.log"
  exit 2
fi
echo "prefix : $PREFIX"
echo "log    : $LOG"
if [ ! -f "$LOG" ]; then echo "MISSING: $LOG"; exit 2; fi

# Prefix equality does NOT bind a log to an artifact set: a re-run under the same
# prefix overwrites the CSV/parquet while an earlier run's clean log still matches
# the prefix. Timestamps do not fix that either — a stale log can be copied or
# touched (Codex delta 4), and `stat` flags are not portable. Bind by CONTENT
# instead: the run's own timing JSON, written beside the CSV by the same run,
# carries the counts that run printed in its log footer.
EXPECT_ARGS=()
TIMING=${PORT}_timing.json
if [ ! -f "$TIMING" ]; then
  echo "FAIL: no $(basename "$TIMING"). Without the run's own timing JSON the"
  echo "      failure gate cannot be bound to these artifacts at all, and a log"
  echo "      from a different run could certify them. Refusing to proceed."
  exit 2
fi
n_g=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))['n_gages_processed'])" "$TIMING" 2>/dev/null || true)
n_c=$($PY -c "import json,sys;print(json.load(open(sys.argv[1])).get('n_columns',''))" "$TIMING" 2>/dev/null || true)
# Must be a positive integer. A JSON null prints as "None", which is non-empty
# and would sail past a bare -z test and then break the integer comparison with
# "integer expression expected" — and with no set -e the script would carry on
# unbound (Codex delta 6).
case "${n_g:-}" in
  "" | *[!0-9]*)
    echo "FAIL: $(basename "$TIMING") has no usable n_gages_processed"
    echo "      (got: '${n_g:-<empty>}'). Refusing to run an unbound failure gate."
    exit 2 ;;
esac
if [ "$n_g" -le 0 ]; then
  echo "FAIL: $(basename "$TIMING") reports n_gages_processed=$n_g."
  exit 2
fi
# The CSV and the timing JSON are NOT written atomically together: the runner
# writes the CSV, then collects provenance, then writes the timing JSON. A rerun
# that aborts in that interval leaves a NEW CSV beside the PREVIOUS run's timing
# JSON (Codex delta 5). Bind them by content: the CSV's own row count must equal
# the gage count the timing JSON claims.
# Count rows by PARSING the CSV, not with wc -l: a file without a trailing
# newline undercounts by one, and any quoted field containing a newline
# overcounts. Reading one column keeps this to a second or two.
csv_rows=$($PY -c "import pandas as pd,sys;print(len(pd.read_csv(sys.argv[1],usecols=['gage_id'],dtype=str)))" "$PORT"_signatures.csv 2>/dev/null || echo "")
if [ -z "$csv_rows" ]; then
  echo "FAIL: could not read ${PREFIX}_signatures.csv to count its rows."
  exit 2
fi
if [ "$csv_rows" -ne "$n_g" ]; then
  echo "FAIL: ${PREFIX}_signatures.csv has $csv_rows data rows but"
  echo "      $(basename "$TIMING") reports n_gages_processed=$n_g."
  echo "      The timing JSON does not describe this CSV — most likely a rerun"
  echo "      aborted between writing the CSV and writing the timing JSON, leaving"
  echo "      a stale JSON beside new artifacts. Re-run, or supply matching files."
  exit 2
fi
echo "binding: CSV rows ($csv_rows) == timing n_gages_processed ($n_g)"
EXPECT_ARGS+=(--expect-gages "$n_g")
case "${n_c:-}" in "" | *[!0-9]*) n_c="" ;; esac
[ -n "${n_c:-}" ] && EXPECT_ARGS+=(--expect-columns "$n_c")
echo "binding: timing.json reports n_gages_processed=$n_g n_columns=${n_c:-absent}"

fail=0
run_diag () {  # like run(), but its exit status does NOT affect the verdict
  local label=$1; shift
  echo; echo "═══════════════════════════════════════════════════════════════"
  echo "DIAGNOSTIC (not a gate): $label"
  echo "═══════════════════════════════════════════════════════════════"
  "$@" || echo ">>> $label reported a nonzero exit (informational only)"
}

run () {  # run <label> <cmd...>
  local label=$1; shift
  echo; echo "═══════════════════════════════════════════════════════════════"
  echo "GATE: $label"
  echo "═══════════════════════════════════════════════════════════════"
  "$@"
  local rc=$?
  if [ $rc -ne 0 ]; then echo ">>> $label FAILED (exit $rc)"; fail=1
  else echo ">>> $label passed"; fi
}

for f in "$REF"_signatures.csv "$PORT"_signatures.csv; do
  [ -f "$f" ] || { echo "MISSING: $f"; exit 2; }
done

# 1. Strict schema + gage-set equality (the acceptance gate proper).
run "schema equality (1,653 columns / 6,678 gages)" \
  $PY docs/benchmarks/check_schema_equality.py \
      "$REF"_signatures.csv "$PORT"_signatures.csv \
      --expect-reference-columns 1653 --expect-candidate-columns 1653

# 2. Swallowed per-gage family exceptions in the run log.
#    NOTE: this REFUSES to certify a log ending in R's deferred-warning
#    truncation banner. run_rpkg_benchmark.R sets options(warn = 1) at the top
#    so warnings print immediately; a log from an OLDER runner (or one whose
#    warn option was overridden) can still end in the banner and will be
#    refused — the correct outcome, not a green light.
# --require-error-tally is rpkg-specific: only the R runner prints
# "Gages errored: N". Dry-running this harness against a Python prefix therefore
# fails this gate by design, not because of a defect.
run "no swallowed signature-family failures" \
  $PY docs/benchmarks/check_signature_failures.py "$LOG" \
      --expect-prefix "$PREFIX" --require-error-tally "${EXPECT_ARGS[@]}"

# 3. Cross-language annual parquet (the sufficient annual check — unlike
#    validate_annual_values.py, which is self-referential).
if [ -f "$PORT"_signatures_annual.parquet ]; then
  run "annual parquet equality vs Julia (tol 1e-6)" \
    $PY docs/benchmarks/check_annual_parquet_equality.py \
        --ref "$REF"_signatures_annual.parquet \
        --port "$PORT"_signatures_annual.parquet --tol 1e-6 --max-show 10
else
  echo; echo ">>> annual parquet MISSING for $PREFIX — collector did not write"; fail=1
fi

# 4. DIAGNOSTIC, NOT A GATE. compare_rpkg_vs_golden_julia.py compares only the
#    column INTERSECTION and exits 0 even when whole families are missing, so its
#    exit status must never be read as acceptance — gate 1 is what proves the sets
#    are equal. It is run here for its R² tier breakdown, and its result is
#    reported separately from the verdict below.
run_diag "identity-R² tiers vs Julia golden (DIAGNOSTIC)" \
  env COMPARE_REFERENCE_CSV="$REF"_signatures.csv \
      COMPARE_CANDIDATE_CSV="$PORT"_signatures.csv \
      COMPARE_OUTPUT_DIR="$DIR" \
      $PY docs/benchmarks/compare_rpkg_vs_golden_julia.py

echo; echo "═══════════════════════════════════════════════════════════════"
if [ $fail -eq 0 ]; then echo "ALL rpkg ACCEPTANCE GATES PASSED"; else
  echo "ONE OR MORE rpkg ACCEPTANCE GATES FAILED — see above"; fi
echo "═══════════════════════════════════════════════════════════════"
exit $fail
