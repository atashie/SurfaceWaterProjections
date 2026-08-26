#!/usr/bin/env bash
# Acceptance gates for the feature-complete rpkg benchmark, against the Julia
# reference produced in Phase 0 of the August 2026 port campaign.
#
# These are GATES, not diagnostics: each exits nonzero on failure and the script
# reports a consolidated verdict. The compare_*_vs_golden_julia.py scripts are
# deliberately NOT gates — they compare only the column intersection and exit 0
# when whole families are missing — so one is run here for its identity-R² tiers
# only, after the strict schema gate has already proven the sets are equal.
#
# Usage:  bash docs/benchmarks/run_rpkg_acceptance_gates.sh [PREFIX]
#           PREFIX defaults to rpkg_full_26aug2026
set -uo pipefail

DIR=/Volumes/Untitled/processedOuts_portref_24aug2026
REF=$DIR/streamflow_1993_2025_60pct_portref_24aug2026
PREFIX=${1:-rpkg_full_26aug2026}
PORT=$DIR/$PREFIX
PY=${PY:-.venv/bin/python}
LOG=$DIR/rpkg_full_run.log

fail=0
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
#    truncation banner. Until run_rpkg_benchmark.R sets options(warn = 1),
#    expect this gate to report that it cannot substantiate a pass — which is
#    the correct outcome, not a green light.
run "no swallowed signature-family failures" \
  $PY docs/benchmarks/check_signature_failures.py "$LOG"

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

# 4. Identity-R² tiers (diagnostic detail, meaningful only after gate 1 passed).
#    This script is configured by env vars, not flags.
run "identity-R² comparison vs Julia golden" \
  env COMPARE_REFERENCE_CSV="$REF"_signatures.csv \
      COMPARE_CANDIDATE_CSV="$PORT"_signatures.csv \
      COMPARE_OUTPUT_DIR="$DIR" \
      $PY docs/benchmarks/compare_rpkg_vs_golden_julia.py

echo; echo "═══════════════════════════════════════════════════════════════"
if [ $fail -eq 0 ]; then echo "ALL rpkg ACCEPTANCE GATES PASSED"; else
  echo "ONE OR MORE rpkg ACCEPTANCE GATES FAILED — see above"; fi
echo "═══════════════════════════════════════════════════════════════"
exit $fail
