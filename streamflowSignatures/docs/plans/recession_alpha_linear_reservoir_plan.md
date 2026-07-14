# Design Note: Recession Alpha Under a Fixed Linear-Reservoir Assumption (b = 1)

**Status**: IMPLEMENTED (2026-07-14, Julia). Benchmark re-run pending — bundled with the
annual-values export feature (see `annual_values_export_plan.md`) into a single re-run.
**Requested by**: user (domain decision), 2026-07-14.
**Scope**: Julia canonical only; Python/rpkg ports deferred (consistent with the
annual-values feature). Column names unchanged (frozen CSV contract) — this is a
methodology change behind existing columns, like the April 2026 recession rewrite.

## 1. Problem

In the power-law recession model `-dQ/dt = a·Q^b`, the parameters are estimated from
`log(-dQ/dt) = log(a) + b·log(Q)` — so `log(a)` is the intercept of a regression whose
slope is `b`. Intercept and slope estimates are strongly negatively correlated
(convolved): any variation or uncertainty in `b` propagates directly into `log(a)`.
The previous implementation let `b` vary everywhere alpha was computed:

- `log_a_events`: per-event log_a recomputed against the **year's median b**
  (varies by year and gage) — `recession.jl` (old lines 441-457).
- `log_a_pointcloud`: `median(log(-dQ/dt) − b_pc·log(Q))` with `b_pc` = the year's
  free point-cloud OLS slope — old lines 419-435.
- `log_a_seasonality_*`: sinusoid fit to raw **per-event free-fit intercepts**
  (per-event b) — old line 354/376.

Result: alpha trends/seasonality partly reflect changes in b, not changes in
recession rate.

## 2. Decision (user spec)

- **Alpha is computed under a fixed linear-reservoir assumption, b ≡ 1, for all
  locations and all periods of time.** With b fixed at 1 there is no regression:
  `log_a = median(log(-dQ/dt) − log(Q))` over the relevant point set.
- **b and concavity analyses are unchanged** — `b_pointcloud`, `b_events`, and
  `concavity` keep their free fits exactly as before.

## 3. What changed (all in `julia/src/recession.jl`)

| Output | Before | After |
|---|---|---|
| `log_a_pointcloud` (per year) | `median(log_dQdt − b_pc·log_Q)`, gated on the b-fit succeeding (`var(log_Q) ≥ 1e-8`, `!isnan(b_pc)`) | `median(log_dQdt − log_Q)` whenever the year's point cloud has > 10 pairs — independent of the b fit (no regression → no singularity gate needed) |
| `log_a_events` (per year) | per-event `median(log_dQdt − median_b_year·log_Q)`, year median | per-event `median(log_dQdt − log_Q)` (new helper `event_log_a_b1`), year median; same event set (fit-success gating unchanged) |
| `log_a_seasonality_*` (6 scalars) | sinusoid over per-event **free-fit intercepts** | sinusoid over per-event **b=1 log_a** values (same events, same day-of-year keys) |
| `b_pointcloud`, `b_events`, `concavity` | free OLS fits | **unchanged** |
| `alpha_linear`, `recession_alpha_point_cloud_linear_reservoir`, param BFIs | already b=1 (discrete `Q_{i+1}/Q_i`) | **unchanged** |
| `n_recession_events`, event identification, min_events gate | — | **unchanged** (event acceptance still keyed to the free fit succeeding, which requires the same ≥ 2 valid (Q, −dQ) pairs the b=1 median needs, so event counts are identical) |

Point-set conventions preserved: first day of each event removed (storm peak), pairs
require `Q > 0` and `−dQ > 0`, natural log, medians.

## 4. Known mathematical relationships (intended, worth remembering)

- With the forward-difference discretization used here (`−dQ_j = Q_j − Q_{j+1}`,
  evaluated at `Q_j`), each b=1 point satisfies
  `log(−dQ) − log(Q) = log(1 − α_j)` where `α_j = Q_{j+1}/Q_j` is the discrete
  recession ratio. So the new `log_a_pointcloud` is the median of `log(1 − α_j)` —
  i.e., a monotone transform of the same pairs behind `alpha_linear`
  (`median(α_j)`). The two are now consistent views of the same linear-reservoir
  assumption rather than independent estimates; per-year values relate as
  `log_a_pointcloud ≈ log(1 − alpha_linear)` (exact for odd pair counts, since
  medians commute with monotone transforms). Event sets differ slightly
  (`alpha_linear` includes events where the free fit failed), so they are not
  strictly redundant.
- For a true linear reservoir (`Q_{t+1} = α·Q_t`): `log_a = log(1 − α)` exactly.
  The deterministic recession test gage (α = 0.95) yields `log_a = log(0.05) ≈
  −2.9957` under both the old free fit (because its b truly is 1) and the new
  b=1 computation — a built-in continuity check.

## 5. Behavioral edge cases (intentional)

- Years whose point cloud is nearly singular in `log(Q)` (constant-Q clouds) used
  to get NaN for BOTH `b_pointcloud` and `log_a_pointcloud`; now `log_a_pointcloud`
  is computed (it needs no regression) while `b_pointcloud` stays NaN. Rare, and
  arguably the correct reading of "alpha assumes b=1 everywhere".
- `log_a_events` remains gated on the same per-event fit success as before, so the
  event population behind alpha, b, seasonality, and the min-events gate is
  unchanged.

## 6. Validation

- `julia/test/test_recession_alpha_b1.jl`:
  - `event_log_a_b1` unit test — exact `log(0.1)` for a pure exponential event.
  - **Quadratic-reservoir gage** (`Q_{t+1} = Q_t − k·Q_t²`, so `−dQ = k·Q²` exactly):
    free fits must recover `b ≈ 2` (unchanged behavior) while `log_a_pointcloud`
    must equal an independent in-test b=1 recomputation over the same pooled pairs —
    pinning that alpha no longer tracks the fitted b.
  - **Linear-reservoir gage** (α = 0.95): `log_a_pointcloud ≈ log_a_events ≈
    log(0.05)`, `b ≈ 1`, seasonality-input continuity.
- Existing suite + annual-collector suite must stay green (the collector's
  recession-gage assertions double as continuity checks).
- Full-scale effect is quantified at the next benchmark run via the standard
  new-vs-golden Julia comparison (log_a columns and seasonality will diverge from
  golden intentionally; b/concavity/alpha_linear columns must NOT diverge).

## 7. Codex review record (2026-07-14, comprehensive review of both July features)

Initial verdict **NO-GO-as-is — "for documentation consistency rather than a code-path
defect."** All hard constraints verified by code trace: frozen output schema on every
early-return path; event acceptance/counts, min-events gating, `n_recession_events`,
`alpha_linear`, and the recession-alpha scalar behaviorally identical (fit success ⇒
the b=1 median is defined — same valid-pair mask); `all_log_a`/`all_b`/`all_dowy`/
`all_event_years` index-aligned (single guarded push block); no in-tree consumer
assumed `log_a_pointcloud` and `b_pointcloud` are jointly NaN. Feature-A re-scan clean.

| # | Severity | Finding | Fix applied |
|---|---|---|---|
| 1 | MAJOR | CHANGELOG / SIGNATURES.md / claude-skill overstated `log_a_pointcloud ≈ log(1 − alpha_linear)` as a same-sample identity ("by construction"); only this design note carried the event-pool caveat | All three docs reworded: relation is approximate — `alpha_linear` pools ALL identified events, the point cloud only fit-success events; medians commute with the monotone transform only at odd pair counts |
| 2 | MINOR | Tests didn't directly pin (a) seasonality-input alignment or (b) unchanged event counts / gate behavior | Added `b1_seasonal_gage` (alpha alternates 0.90/0.98 by season → strong seasonal log_a signal; pipeline seasonality outputs must exactly reproduce a sinusoid fit to independently recomputed (log_a, dowy) pairs — misalignment would change the fit) and `b1_sparse_gage` (3 events < 25 → all fitted metrics NaN while `n_recession_events` = 3/25 and the gate-independent alpha scalar = 0.9 exactly); plus exact `n_recession_events == 12.0` per-year assertions on the monthly-spike gage |

## 8. Out of scope

- Python/rpkg ports (deferred; sync after the Julia benchmark validates).
- Renaming columns or adding "convolved vs fixed-b" variants (user chose replace).
- The runoff-ratio seasonal-flag bug (tracked separately in CHANGELOG Known Issues).
