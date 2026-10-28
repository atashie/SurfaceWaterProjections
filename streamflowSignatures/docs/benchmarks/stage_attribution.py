#!/usr/bin/env python3
"""stage_attribution.py — the §10.4 attribution identities and the interaction residual of the qualifier-mask +
climate-frame rerun (docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.4, decision 14; 2026-10-26).

    python docs/benchmarks/stage_attribution.py --window 1993-2025 --cohorts stage_cohorts_1993-2025.json
        --s0-csv ... --s0-annual ... --s1-csv ... --s1-annual ... --s1p-csv ... --s1p-annual ... --s2-csv ... --s2-annual ...
        --ledger correction_ledger.csv --out-dir attribution_1993_2025 [--atol 1e-9] [--rtol 1e-6]

Stages (all on the tagged tree): S0 = February input, left; S1 = spliced input, left (executed on the splice-affected
gages + the control sample); S1′ = February input, outer (executed on the composed interaction cohort); S2 = spliced
input, outer (the product). The tool
  1. gates the identities on EXECUTED gages only: S1 = S0 on the control sample, S1′ = S2 on the join-only sample,
     S0 = S2 on class U (both layers, cell for cell, NaN-aware);
  2. defines the full-cohort S1 as the executed S1 rows ∪ the S0 rows of every other gage, and counts every cell so
     taken from S0 as INFERRED (per class);
  3. REPORTS the interaction residual (S2 − S1) − (S1′ − S0) per cell on the interaction cohort (overlap, restored-
     with-inserted, eligibility-gain and Dec 31-neighbour gages — every one executed in S1 and S1′) with the
     tolerance atol + rtol·max(|S0|, |S1|, |S1′|, |S2|); gained / lost cells (a key present in some stages only) are
     listed separately. The residual is not required to be zero (the two fixes interact through eligibility and
     interpolation neighbours); the report names every cell outside the tolerance.
Writes <out-dir>/attribution_<window>.md, identities_<window>.json, residual_cells_<window>.csv,
residual_keys_<window>.csv (gained / lost keys) and inferred_<window>.csv.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd, duckdb


def read_csv(p):
    return pd.read_csv(p, dtype={"gage_id": str}, low_memory=False).set_index("gage_id")


def cells_equal(a, b):
    """NaN-aware exact equality of two aligned numeric frames; non-numeric columns compared as strings."""
    out = pd.DataFrame(index=a.index, columns=a.columns, dtype=bool)
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            xf, yf = x.astype(float), y.astype(float)
            out[c] = ((xf == yf) | (xf.isna() & yf.isna())).to_numpy()
        else:
            out[c] = ((x.astype(str) == y.astype(str)) | (x.isna() & y.isna())).to_numpy()   # pandas 3 keeps missing strings missing
    return out


def summary_identity(csv_a, csv_b, gages, label):
    """Cell-for-cell identity of two summary CSVs on the listed gages (executed on both sides)."""
    a, b = read_csv(csv_a), read_csv(csv_b)
    ids = sorted(set(gages) & set(a.index) & set(b.index)); missing = sorted(set(gages) - set(ids))
    if list(a.columns) != list(b.columns):
        return {"label": label, "status": "FAIL", "reason": "header differs", "n_gages": len(ids), "missing": missing}
    eq = cells_equal(a.loc[ids], b.loc[ids]); bad = (~eq).sum()
    cols = bad[bad > 0].sort_values(ascending=False)
    return {"label": label, "status": "PASS" if not len(cols) and not missing and ids else "FAIL", "n_gages": len(ids),
            "n_cells": int(eq.size), "n_unequal": int((~eq).to_numpy().sum()), "columns_unequal": cols.head(25).to_dict(), "missing": missing}


def annual_identity(con, pq_a, pq_b, gages, label):
    """Keyed (gage_id, water_year, signature) identity of two annual parquets on the listed gages, NaN-aware."""
    con.execute("CREATE OR REPLACE TEMP TABLE gs AS SELECT unnest(?::VARCHAR[]) AS gage_id", [sorted(gages)])
    q = f"""
      WITH A AS (SELECT CAST(gage_id AS VARCHAR) AS gage_id, water_year, signature, "value" FROM read_parquet('{pq_a}') WHERE CAST(gage_id AS VARCHAR) IN (SELECT gage_id FROM gs)),
           B AS (SELECT CAST(gage_id AS VARCHAR) AS gage_id, water_year, signature, "value" FROM read_parquet('{pq_b}') WHERE CAST(gage_id AS VARCHAR) IN (SELECT gage_id FROM gs))
      SELECT count(*) FILTER (WHERE A.gage_id IS NOT NULL AND B.gage_id IS NOT NULL) AS shared,
             count(*) FILTER (WHERE B.gage_id IS NULL) AS only_a, count(*) FILTER (WHERE A.gage_id IS NULL) AS only_b,
             count(*) FILTER (WHERE A.gage_id IS NOT NULL AND B.gage_id IS NOT NULL AND NOT (A."value" = B."value" OR (isnan(A."value") AND isnan(B."value")) OR (A."value" IS NULL AND B."value" IS NULL))) AS unequal
      FROM A FULL OUTER JOIN B USING (gage_id, water_year, signature)"""
    shared, only_a, only_b, unequal = con.execute(q).fetchone()
    return {"label": label, "status": "PASS" if (shared and not only_a and not only_b and not unequal) else "FAIL",
            "shared": int(shared), "only_a": int(only_a), "only_b": int(only_b), "unequal": int(unequal)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", required=True); ap.add_argument("--cohorts", required=True); ap.add_argument("--ledger", required=True)
    for s in ("s0", "s1", "s1p", "s2"):
        ap.add_argument(f"--{s}-csv", required=True); ap.add_argument(f"--{s}-annual", required=True)
    ap.add_argument("--out-dir", required=True); ap.add_argument("--atol", type=float, default=1e-9); ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--memory-limit", default="3GB")
    ap.add_argument("--class-u", default=None, metavar="LIST", help="the realised class-U gages (untouched by both fixes; e.g. the selection manifest minus the affected manifest): the S0 = S2 identity runs on them at both layers; without it that identity is reported NOT RUN (2026-10-28)")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True); w0, w1 = map(int, a.window.split("-"))
    J = json.load(open(a.cohorts)); assert J["window"] == a.window, "cohort window mismatch"
    s1_gages = [l.strip() for l in open(J["S1"]["path"]) if l.strip()]; s1p_gages = [l.strip() for l in open(J["S1prime"]["path"]) if l.strip()]
    led = pd.read_csv(a.ledger, dtype={"gage_id": str, "date": str}, low_memory=False); led = led[led.status == "restored"]
    d = pd.to_datetime(led.date); led["wy"] = np.where(d.dt.month >= 10, d.dt.year + 1, d.dt.year)
    aq = set(led[(led.wy >= w0) & (led.wy <= w1)].gage_id)
    control = sorted(set(s1_gages) - aq)
    members = J["S1prime_members"]; interaction = sorted(set().union(*[set(v) for v in members.values()]))
    join_only = sorted(set(s1p_gages) - set(interaction))
    s0, s2 = read_csv(a.s0_csv), read_csv(a.s2_csv); s1, s1p = read_csv(a.s1_csv), read_csv(a.s1p_csv)
    cohort = list(s0.index)
    U = sorted(l.strip() for l in open(a.class_u) if l.strip()) if a.class_u else None
    if U is not None:
        bad_u = sorted(set(U) & aq); assert not bad_u, f"{len(bad_u)} class-U gages carry a restored day in the window: {bad_u[:5]}"
    con = duckdb.connect(); con.execute(f"SET memory_limit='{a.memory_limit}'"); con.execute("SET threads=4")
    # ---- 1. identities on executed gages
    ident = [summary_identity(a.s1_csv, a.s0_csv, control, "summary: S1 = S0 on the control sample"),
             annual_identity(con, a.s1_annual, a.s0_annual, control, "annual: S1 = S0 on the control sample"),
             summary_identity(a.s1p_csv, a.s2_csv, join_only, "summary: S1′ = S2 on the join-only sample"),
             annual_identity(con, a.s1p_annual, a.s2_annual, join_only, "annual: S1′ = S2 on the join-only sample")]
    if U is not None:
        ident += [summary_identity(a.s2_csv, a.s0_csv, U, "summary: S2 = S0 on class U"), annual_identity(con, a.s2_annual, a.s0_annual, U, "annual: S2 = S0 on class U")]
    else:
        ident += [{"label": "summary: S2 = S0 on class U", "status": "NOT RUN", "reason": "no --class-u list"}, {"label": "annual: S2 = S0 on class U", "status": "NOT RUN", "reason": "no --class-u list"}]
    exec_s1 = sorted(set(s1_gages) & set(s1.index)); exec_s1p = sorted(set(s1p_gages) & set(s1p.index))
    ident.append({"label": "execution: every S1 cohort gage present in the S1 output", "status": "PASS" if len(exec_s1) == len(s1_gages) else "FAIL", "n": len(exec_s1), "expected": len(s1_gages)})
    ident.append({"label": "execution: every S1′ cohort gage present in the S1′ output", "status": "PASS" if len(exec_s1p) == len(s1p_gages) else "FAIL", "n": len(exec_s1p), "expected": len(s1p_gages)})
    # a gage of the S1 cohort absent from the S1 output would be a qualification change by the splice alone (not expected: the gage sets are frozen)
    # ---- 2. the full-cohort S1: executed rows ∪ S0 rows, INFERRED cells counted per class
    inferred_gages = sorted(set(cohort) - set(exec_s1))
    inf_rows = {"gages_inferred": len(inferred_gages), "gages_executed": len(exec_s1), "summary_cells_inferred": int(len(inferred_gages) * s0.shape[1]),
                "annual_rows_inferred": int(con.execute(f"SELECT count(*) FROM read_parquet('{a.s0_annual}') WHERE CAST(gage_id AS VARCHAR) IN (SELECT unnest(?::VARCHAR[]))", [inferred_gages]).fetchone()[0]) if inferred_gages else 0,
                "note": "the full-cohort S1 = executed S1 rows ∪ S0 rows of every other gage; every cell so taken from S0 is INFERRED (equal to S0 on the control-sample evidence, not by execution)"}
    pd.DataFrame({"gage_id": cohort, "s1_source": ["executed" if g in set(exec_s1) else "inferred_from_S0" for g in cohort]}).to_csv(f"{a.out_dir}/inferred_{a.window}.csv", index=False)
    # ---- 3. the interaction residual on the interaction cohort
    # Coverage (Codex 2026-10-28): the S1 cohort is the splice-affected gages + the control sample, so an interaction gage of
    # the eligibility-gain or Dec 31 kind is NOT executed in S1. That is admissible only when the gage has no restored day
    # inside the window — then S1 = S0 by the rule the control-sample identity verifies — and its S1 cells are taken from
    # S0 and labelled INFERRED everywhere below; an interaction gage WITH a restored day absent from S1 fails the coverage gate.
    cat_of = {g: k for k, v in members.items() for g in v}
    s1_missing = sorted(set(interaction) - set(s1.index)); s1_missing_bad = sorted(set(s1_missing) & aq); s1p_missing = sorted(set(interaction) - set(s1p.index))
    ident.append({"label": "coverage: every interaction gage with a restored day inside the window executed in S1 (the others carry S1 inferred = S0)",
                  "status": "FAIL" if s1_missing_bad else "PASS", "interaction": len(interaction), "executed_in_S1": len(interaction) - len(s1_missing),
                  "s1_inferred_from_S0": len(s1_missing), "s1_inferred_gages": {g: cat_of.get(g) for g in s1_missing}, "with_restored_day_not_executed": s1_missing_bad})
    ident.append({"label": "coverage: every interaction gage executed in S1′", "status": "PASS" if not s1p_missing else "FAIL", "n": len(interaction) - len(s1p_missing), "expected": len(interaction), "missing": s1p_missing})
    con.execute("CREATE OR REPLACE TEMP TABLE ic AS SELECT unnest(?::VARCHAR[]) AS gage_id", [interaction])
    con.execute("CREATE OR REPLACE TEMP TABLE inf AS SELECT unnest(?::VARCHAR[]) AS gage_id", [s1_missing])
    def tbl(name, path):
        con.execute(f"""CREATE OR REPLACE TEMP TABLE {name} AS SELECT CAST(gage_id AS VARCHAR) AS gage_id, water_year, signature, "value" AS v, 1 AS p
                        FROM read_parquet('{path}') WHERE CAST(gage_id AS VARCHAR) IN (SELECT gage_id FROM ic)""")
    for n_, p_ in (("t0", a.s0_annual), ("t1", a.s1_annual), ("t1p", a.s1p_annual), ("t2", a.s2_annual)): tbl(n_, p_)
    con.execute("INSERT INTO t1 SELECT gage_id, water_year, signature, v, 2 AS p FROM t0 WHERE gage_id IN (SELECT gage_id FROM inf)")
    res = con.execute("""
      SELECT gage_id, water_year, signature, t0.v AS s0, t1.v AS s1, t1p.v AS s1p, t2.v AS s2, t0.p AS p0, t1.p AS p1, t1p.p AS p1p, t2.p AS p2
      FROM t0 FULL OUTER JOIN t1 USING (gage_id, water_year, signature) FULL OUTER JOIN t1p USING (gage_id, water_year, signature) FULL OUTER JOIN t2 USING (gage_id, water_year, signature)""").df()
    # key PRESENCE (the p markers) and value FINITENESS (a stored NaN / NULL) are separate masks (Codex 2026-10-28)
    pres = res[["p0", "p1", "p1p", "p2"]].notna(); pres.columns = ["s0", "s1", "s1p", "s2"]
    fin = res[["s0", "s1", "s1p", "s2"]].notna()
    res["s1_inferred"] = res["p1"] == 2
    all_present = pres.all(axis=1)
    keys = res[~all_present].copy()
    keys["present_in"] = [",".join(c for c in ("s0", "s1", "s1p", "s2") if r_[c]) for _, r_ in pres[~all_present].iterrows()]
    keys[["gage_id", "water_year", "signature", "present_in", "s1_inferred"]].to_csv(f"{a.out_dir}/residual_keys_{a.window}.csv", index=False)
    nan_all = all_present & ~fin.any(axis=1); nan_mixed = all_present & fin.any(axis=1) & ~fin.all(axis=1)
    n_nan_all = int(nan_all.sum()); n_nan_mixed = int(nan_mixed.sum())
    r = res[all_present & fin.all(axis=1)].copy()
    r["splice_effect"] = r.s1 - r.s0; r["join_effect"] = r.s1p - r.s0; r["residual"] = (r.s2 - r.s1) - (r.s1p - r.s0)
    r["tol"] = a.atol + a.rtol * r[["s0", "s1", "s1p", "s2"]].abs().max(axis=1)
    r["outside_tolerance"] = r.residual.abs() > r.tol
    r.sort_values("residual", key=lambda x: x.abs(), ascending=False).to_csv(f"{a.out_dir}/residual_cells_{a.window}.csv", index=False)
    out_tol = r[r.outside_tolerance]; n_r_exec = int((~r.s1_inferred).sum()); n_r_inf = int(r.s1_inferred.sum()); n_out_exec = int((~out_tol.s1_inferred).sum())
    per_sig = out_tol.groupby("signature").agg(cells=("residual", "size"), gages=("gage_id", "nunique"), max_abs_residual=("residual", lambda x: x.abs().max())).sort_values("cells", ascending=False)
    key_counts = keys.present_in.value_counts().to_dict()
    # summary layer on the interaction cohort: residual of the per-gage _mean columns and every other numeric column
    ids = sorted(set(interaction) & set(s0.index) & set(s1p.index) & set(s2.index)); ids_exec = [g for g in ids if g in s1.index]; ids_inf = [g for g in ids if g not in s1.index]
    num = [c for c in s0.columns if all(pd.api.types.is_numeric_dtype(x[c]) for x in (s0, s1, s1p, s2))]
    S0, S1P, S2 = (x.loc[ids, num].astype(float) for x in (s0, s1p, s2))
    S1 = pd.concat([s1.loc[ids_exec, num], s0.loc[ids_inf, num]]).loc[ids].astype(float)   # S1 of the inferred gages := S0
    RS = (S2 - S1) - (S1P - S0); TOL = a.atol + a.rtol * pd.concat([S0.abs(), S1.abs(), S1P.abs(), S2.abs()]).groupby(level=0).max()
    both = S0.notna() & S1.notna() & S1P.notna() & S2.notna()
    out_s = (RS.abs() > TOL) & both; col_counts = out_s.sum(); col_counts = col_counts[col_counts > 0].sort_values(ascending=False)
    state_changes = int((~both & (S0.notna() | S1.notna() | S1P.notna() | S2.notna())).to_numpy().sum())
    verdict = "FAIL" if any(i["status"] == "FAIL" for i in ident) else ("INCOMPLETE" if any(i["status"] == "NOT RUN" for i in ident) else "PASS")
    summary = {"window": a.window, "identities": ident, "verdict_identities": verdict, "inferred": inf_rows,
               "interaction_cohort": {"n_gages": len(interaction), "members": {k: len(v) for k, v in members.items()}, "executed_in_all_four_stages": len(interaction) - len(s1_missing),
                                      "s1_inferred_from_S0": {g: cat_of.get(g) for g in s1_missing}, "join_only_sample": len(join_only), "control_sample": len(control)},
               "residual_annual": {"cells_compared": int(len(r)), "cells_compared_executed": n_r_exec, "cells_compared_s1_inferred": n_r_inf, "outside_tolerance": int(len(out_tol)),
                                   "outside_tolerance_executed": n_out_exec, "max_abs_residual": float(r.residual.abs().max()) if len(r) else 0.0,
                                   "keys_not_in_every_stage": int(len(keys)), "keys_by_presence": key_counts, "present_in_all_four_NaN_in_all_four": n_nan_all,
                                   "present_in_all_four_NaN_state_differs": n_nan_mixed, "atol": a.atol, "rtol": a.rtol},
               "residual_summary": {"gages": len(ids), "gages_executed_in_all_four": len(ids_exec), "gages_s1_inferred": len(ids_inf), "cells_compared": int(both.to_numpy().sum()),
                                    "outside_tolerance": int(out_s.to_numpy().sum()), "outside_tolerance_executed": int(out_s.loc[ids_exec].to_numpy().sum()) if ids_exec else 0,
                                    "columns_outside": col_counts.head(40).to_dict(), "cells_with_a_NaN_state_difference": state_changes}}
    json.dump(summary, open(f"{a.out_dir}/identities_{a.window}.json", "w"), indent=1, default=str)
    L = [f"# Stage attribution — window {a.window} (plan §10.4)", "", f"Identities on executed gages: **{verdict}**", ""]
    for i in ident:
        L.append(f"- {i['status']}: {i['label']} — " + ", ".join(f"{k}={v}" for k, v in i.items() if k not in ("label", "status", "columns_unequal", "missing")) + (f"; unequal columns {i['columns_unequal']}" if i.get("columns_unequal") else "") + (f"; missing {i['missing'][:10]}" if i.get("missing") else ""))
    L += ["", f"Full-cohort S1 = executed S1 rows ∪ S0 rows: {inf_rows['gages_executed']} gages executed, {inf_rows['gages_inferred']} gages INFERRED from S0 "
          f"({inf_rows['summary_cells_inferred']:,} summary cells, {inf_rows['annual_rows_inferred']:,} annual rows) — `inferred_{a.window}.csv`.", "",
          f"## Interaction residual (S2 − S1) − (S1′ − S0) on the interaction cohort ({len(interaction)} gages: {', '.join(f'{k} {len(v)}' for k, v in members.items())})", "",
          f"Executed in all four stages: {len(interaction) - len(s1_missing)} gages; {len(s1_missing)} gages are not in the S1 cohort (no restored day inside the window: "
          f"{', '.join(f'{g} [{cat_of.get(g)}]' for g in s1_missing) or 'none'}) and carry S1 INFERRED = S0 in the residual below (`s1_inferred` column).", "",
          f"Annual layer: {len(r):,} cells present in all four stages with every value finite ({n_r_exec:,} on executed gages, {n_r_inf:,} with S1 inferred); "
          f"outside tolerance (atol {a.atol:g} + rtol {a.rtol:g}·max|value|): {len(out_tol):,} ({n_out_exec:,} on executed gages); max |residual| {summary['residual_annual']['max_abs_residual']:.4g}. "
          f"Keys absent from at least one stage: {len(keys):,} ({key_counts}); keys present in all four stages but NaN in all four: {n_nan_all:,}; present in all four with a NaN-state difference: {n_nan_mixed:,} — `residual_keys_{a.window}.csv`.", ""]
    if len(per_sig):
        L += ["| signature | cells outside tol | gages | max abs residual |", "|---|---|---|---|"] + [f"| {s} | {int(x.cells)} | {int(x.gages)} | {x.max_abs_residual:.4g} |" for s, x in per_sig.head(40).iterrows()] + [""]
    L += [f"Summary layer: {len(ids)} gages ({len(ids_exec)} executed in all four stages, {len(ids_inf)} with S1 inferred = S0), {summary['residual_summary']['cells_compared']:,} cells with every stage finite; outside tolerance {summary['residual_summary']['outside_tolerance']:,} ({summary['residual_summary']['outside_tolerance_executed']:,} on executed gages); "
          f"cells with a NaN-state difference between stages {state_changes:,}.", ""]
    if len(col_counts):
        L += ["| column | cells outside tol |", "|---|---|"] + [f"| {c} | {int(n)} |" for c, n in col_counts.head(40).items()] + [""]
    L += ["The residual is REPORTED, not gated (the two fixes interact through eligibility and interpolation neighbours); the identities above are the gates.", ""]
    open(f"{a.out_dir}/attribution_{a.window}.md", "w").write("\n".join(L)); print("\n".join(L[:12]))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
