#!/usr/bin/env python3
"""check_summary_equality.py — the numerical GATE on a summary CSV against a reference (plan E7 of
docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md, 2026-10-26): every summary statistic,
per-gage scalar and QA flag compared cell for cell.

    python docs/benchmarks/check_summary_equality.py --ref REF.csv --cand CAND.csv [--rtol 1e-6] [--atol 1e-12]
        [--ref-annual REF_annual.parquet --cand-annual CAND_annual.parquet] [--waive-column REGEX ...]
        [--gages LIST] [--max-show 20] [--out REPORT.json]

Gates (all must hold): the same gage set and header; every non-numeric cell identical (metadata,
flags, strings); every numeric cell within |d| <= atol + rtol*max(|a|,|b|) with GATE-STATE agreement
(a cell NaN in one file and finite in the other, or an infinity that is not the same infinity, is
a failure, never a residual). --waive-column REGEX excepts whole columns from the numeric gate
(reported; a gate-state difference is NEVER waived, not even there). With both annual parquets:
a (gage, base) whose reference annual series holds a tie or a near-tie (two values within 1e-9
relative) with annual values that differ within the annual tolerance, and a (gage, base) whose
annual series differs beyond --annual-rtol (an annual residual check_annual_parquet_equality.py
names on its own), have their summary cells gated against the statistics recomputed from each
side's OWN annual series with the port's generate_stats — the honest propagation of a named annual
residual passes and is reported, anything else is a gated mismatch. One library-level residual is
named explicitly: a series holding both −0.0 and +0.0, whose rank statistics Julia and the ports
compute with different tie counts (the standing Known Issue) — reported, not gated. Exit 1 on any
gated mismatch. The JSON written by --out carries the sha256 of both inputs and the compared-gage
list (consumed by assess_rerun_changes.py).
"""
import argparse
import json
import re
import sys

import numpy as np
import pandas as pd

SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval", "_mean", "_median",
            "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean", "_pettitt_delta_mean",
            "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
RANK = {s for s in SUFFIXES if s not in ("_senn_slp", "_linear_slp", "_mean", "_median")}
_SUF = sorted(SUFFIXES, key=len, reverse=True)


def split_col(c):
    for s in _SUF:
        if c.endswith(s):
            return c[:-len(s)], s
    return c, None


def annual_residual_pairs(ref_annual, cand_annual, gages, rtol=1e-9, atol=1e-12):
    """(gage, base) pairs with an annual value beyond the annual tolerance (or a NaN-state
    difference) between the two files: a summary difference of such a pair is the propagation
    of an annual residual that check_annual_parquet_equality.py already characterises."""
    import duckdb
    con = duckdb.connect()
    g_filter = ""
    if gages is not None:
        con.execute("CREATE TEMP TABLE gs AS SELECT unnest(?::VARCHAR[]) AS gage_id", [sorted(gages)]); g_filter = " WHERE gage_id IN (SELECT gage_id FROM gs)"
    q = f"""WITH r AS (SELECT CAST(gage_id AS VARCHAR) gage_id, CAST(signature AS VARCHAR) signature, CAST(water_year AS INTEGER) wy, CAST("value" AS DOUBLE) v FROM read_parquet('{ref_annual}'){g_filter}),
                 c AS (SELECT CAST(gage_id AS VARCHAR) gage_id, CAST(signature AS VARCHAR) signature, CAST(water_year AS INTEGER) wy, CAST("value" AS DOUBLE) v FROM read_parquet('{cand_annual}'){g_filter})
            SELECT DISTINCT coalesce(r.gage_id, c.gage_id), coalesce(r.signature, c.signature)
            FROM r FULL OUTER JOIN c USING (gage_id, signature, wy)
            WHERE r.gage_id IS NULL OR c.gage_id IS NULL
               OR ((r.v IS NULL OR isnan(r.v)) <> (c.v IS NULL OR isnan(c.v)))
               OR (NOT (r.v IS NULL OR isnan(r.v)) AND abs(r.v - c.v) > {atol} + {rtol} * greatest(abs(r.v), abs(c.v)))"""
    return set(con.execute(q).fetchall())


def recompute_pairs(ref_annual, cand_annual, pairs):
    """{(gage, base): (stats_from_ref_series, stats_from_cand_series)} — the 16 statistics recomputed with the
    Python port's generate_stats from each side's own annual series (the production gates; the family
    exemptions of recession / elasticity applied), so a summary cell can be checked against the honest
    propagation of its annual values instead of being exempted."""
    import duckdb
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).parent.parent.parent / "python"))
    from streamflow_signatures.stats import generate_stats
    from streamflow_signatures.config import (NA_TREND_MIN_FRACTION, NA_DECADE_MIN_FRACTION, MIN_VALUES_FOR_STATS,
                                              CP_START_WATER_YEAR, CP_END_WATER_YEAR, CP_MIN_TOTAL_OBS, CP_MIN_SEGMENT_OBS)
    if not pairs:
        return {}
    con = duckdb.connect()
    con.execute("CREATE TEMP TABLE pr AS SELECT * FROM (VALUES " + ", ".join(f"('{g}', '{b}')" for g, b in pairs) + ") t(gage_id, signature)")
    out = {}
    cp = {"start_year": CP_START_WATER_YEAR, "end_year": CP_END_WATER_YEAR, "min_total_obs": CP_MIN_TOTAL_OBS, "min_segment_obs": CP_MIN_SEGMENT_OBS}
    exempt = ("log_a_", "b_pointcloud", "b_events", "concavity", "n_recession_events", "alpha_linear", "elasticity_")
    series = {}
    for side, path in (("ref", ref_annual), ("cand", cand_annual)):
        q = con.execute(f"SELECT CAST(A.gage_id AS VARCHAR) g, CAST(A.signature AS VARCHAR) s, CAST(A.water_year AS INTEGER) wy, CAST(A.\"value\" AS DOUBLE) v "
                        f"FROM read_parquet('{path}') A JOIN pr ON CAST(A.gage_id AS VARCHAR) = pr.gage_id AND CAST(A.signature AS VARCHAR) = pr.signature ORDER BY 1, 2, 3").df()
        for (g, b), grp in q.groupby(["g", "s"]):
            series[(side, g, b)] = grp
    for g, b in pairs:
        res = []
        for side in ("ref", "cand"):
            grp = series.get((side, g, b))
            if grp is None or grp.empty:
                res.append({}); continue
            df = pd.DataFrame({"water_year": grp["wy"].astype(int).to_numpy(), b: grp["v"].to_numpy(dtype=float)})
            kw = dict(changepoint=cp)
            if not b.startswith(exempt):
                kw.update(trend_completeness=NA_TREND_MIN_FRACTION, decade_completeness=NA_DECADE_MIN_FRACTION, min_values_for_stats=MIN_VALUES_FOR_STATS)
            res.append(generate_stats(df, value_cols=[b], **kw))
        out[(g, b)] = tuple(res)
    return out


def signed_zero_pairs(ref_annual, gages):
    """(gage, base) pairs whose reference annual series holds a negative zero next to a zero: Julia's
    Mann-Kendall / Spearman tie correction counts −0.0 and +0.0 as distinct values (the standing
    Known Issue 'Julia unique treats −0.0 ≠ +0.0'), the ports as one tie, so the rank statistics
    of such a series legitimately differ between the languages — a named residual."""
    import duckdb
    con = duckdb.connect()
    g_filter = ""
    if gages is not None:
        con.execute("CREATE TEMP TABLE gs AS SELECT unnest(?::VARCHAR[]) AS gage_id", [sorted(gages)]); g_filter = " AND gage_id IN (SELECT gage_id FROM gs)"
    q = f"""SELECT CAST(gage_id AS VARCHAR), CAST(signature AS VARCHAR) FROM read_parquet('{ref_annual}')
            WHERE "value" = 0 {g_filter} GROUP BY 1, 2 HAVING bool_or(signbit(CAST("value" AS DOUBLE))) AND bool_or(NOT signbit(CAST("value" AS DOUBLE)))"""
    return set(con.execute(q).fetchall())


def tie_sensitive_pairs(ref_annual, cand_annual, gages):
    import duckdb
    con = duckdb.connect()
    g_filter = ""
    if gages is not None:
        con.execute("CREATE TEMP TABLE gs AS SELECT unnest(?::VARCHAR[]) AS gage_id", [sorted(gages)]); g_filter = " WHERE gage_id IN (SELECT gage_id FROM gs)"
    q = f"""WITH r AS (SELECT CAST(gage_id AS VARCHAR) gage_id, CAST(signature AS VARCHAR) signature, CAST(water_year AS INTEGER) wy, CAST("value" AS DOUBLE) v FROM read_parquet('{ref_annual}'){g_filter}),
                 c AS (SELECT CAST(gage_id AS VARCHAR) gage_id, CAST(signature AS VARCHAR) signature, CAST(water_year AS INTEGER) wy, CAST("value" AS DOUBLE) v FROM read_parquet('{cand_annual}'){g_filter}),
                 j AS (SELECT r.gage_id, r.signature, r.wy, r.v rv, c.v cv FROM r LEFT JOIN c USING (gage_id, signature, wy)),
                 ident AS (SELECT gage_id, signature, bool_and(coalesce(rv = cv, false) OR ((rv IS NULL OR isnan(rv)) AND (cv IS NULL OR isnan(cv)))) AS bit FROM j GROUP BY 1, 2),
                 within AS (SELECT gage_id, signature, bool_and(((rv IS NULL OR isnan(rv)) AND (cv IS NULL OR isnan(cv))) OR (rv IS NOT NULL AND cv IS NOT NULL AND abs(rv - cv) <= 1e-12 + 1e-9 * greatest(abs(rv), abs(cv)))) AS ok FROM j GROUP BY 1, 2),
                 s AS (SELECT gage_id, signature, v, lag(v) OVER (PARTITION BY gage_id, signature ORDER BY v) AS pv FROM r WHERE v IS NOT NULL AND NOT isnan(v)),
                 ties AS (SELECT gage_id, signature, bool_or(abs(v - pv) <= 1e-9 * greatest(abs(v), 1e-300)) AS tie FROM s WHERE pv IS NOT NULL GROUP BY 1, 2)
            SELECT ties.gage_id, ties.signature FROM ties JOIN ident USING (gage_id, signature) JOIN within USING (gage_id, signature)
            WHERE ties.tie AND NOT ident.bit AND within.ok"""
    return set(con.execute(q).fetchall())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True); ap.add_argument("--cand", required=True)
    ap.add_argument("--rtol", type=float, default=1e-6); ap.add_argument("--atol", type=float, default=1e-12)
    ap.add_argument("--ref-annual"); ap.add_argument("--cand-annual"); ap.add_argument("--waive-column", action="append", default=[])
    ap.add_argument("--gages"); ap.add_argument("--max-show", type=int, default=20); ap.add_argument("--out")
    ap.add_argument("--allow-gate-state-column", action="append", default=[], metavar="REGEX",
                    help="the ONE named exception to 'a gate-state difference is never waived': in a column matching REGEX, a NaN-vs-finite "
                         "cell of a tie-bearing (gage, base) pair is accepted only when each side's value is reproduced by the statistic "
                         "recomputed from that side's own annual series (the rank / Pettitt ambiguity of a zero-heavy series, where the "
                         "languages break a tie differently and one side's period is constant). Every such cell is printed and counted "
                         "(gate_state_tie_verified); the regex is recorded in the JSON and must be accepted by name downstream.")
    ap.add_argument("--annual-rtol", type=float, default=1e-9, help="annual tolerance above which a (gage, base) pair counts as an annual residual")
    a = ap.parse_args()
    gages = set(l.strip() for l in open(a.gages) if l.strip()) if a.gages else None
    ref = pd.read_csv(a.ref, dtype={"gage_id": str}, low_memory=False); cand = pd.read_csv(a.cand, dtype={"gage_id": str}, low_memory=False)
    if gages is not None:
        ref = ref[ref["gage_id"].isin(gages)]; cand = cand[cand["gage_id"].isin(gages)]
    ref = ref.set_index("gage_id").sort_index(); cand = cand.set_index("gage_id").sort_index()
    ok = True; rep = {}
    if list(ref.columns) != list(cand.columns):
        print(f"FAIL header: ref {ref.shape[1]} cols, cand {cand.shape[1]}; only-ref {sorted(set(ref.columns) - set(cand.columns))[:10]}, only-cand {sorted(set(cand.columns) - set(ref.columns))[:10]}"); ok = False
    if set(ref.index) != set(cand.index):
        print(f"FAIL gage set: ref {len(ref)}, cand {len(cand)}; only-ref {sorted(set(ref.index) - set(cand.index))[:10]}, only-cand {sorted(set(cand.index) - set(ref.index))[:10]}"); ok = False
    common = sorted(set(ref.index) & set(cand.index)); cols = [c for c in ref.columns if c in cand.columns]
    ref = ref.loc[common, cols]; cand = cand.loc[common, cols]
    waived = [c for c in cols if any(re.search(p, c) for p in a.waive_column)]
    ties, resid, recomputed, szero = set(), set(), {}, set()
    if a.ref_annual and a.cand_annual:
        ties = tie_sensitive_pairs(a.ref_annual, a.cand_annual, set(common)); print(f"tie-sensitive (gage, base) pairs: {len(ties)}")
        resid = annual_residual_pairs(a.ref_annual, a.cand_annual, set(common), a.annual_rtol, a.atol)
        print(f"(gage, base) pairs with an annual residual beyond {a.annual_rtol:g}: {len(resid)} {sorted(resid)[:8]} — their summary cells are gated against "
              "the statistics recomputed from each side's OWN annual series (the honest propagation of a named annual residual)")
        recomputed = recompute_pairs(a.ref_annual, a.cand_annual, resid | ties)
        szero = signed_zero_pairs(a.ref_annual, set(common))
        print(f"(gage, base) pairs with a signed-zero tie in the reference series (rank statistics reported as the named Known-Issue residual): {len(szero)} {sorted(szero)[:8]}")
    n_gated = n_bad = n_state = n_waived = n_tie = n_resid = n_szero = n_state_ok = 0; worst = (0.0, None, None); bad_cols = {}; shown = []
    for c in cols:
        r, k = ref[c], cand[c]
        base, suf = split_col(c)
        if pd.api.types.is_numeric_dtype(r) and pd.api.types.is_numeric_dtype(k):
            rn, kn = r.isna(), k.isna()
            state = rn != kn
            both = ~rn & ~kn
            rv, kv = r[both].astype(float), k[both].astype(float)
            # a non-finite value is a state of its own: equal only to the same infinity, never within a tolerance
            inf_r, inf_k = ~np.isfinite(rv), ~np.isfinite(kv)
            inf_mis = (inf_r | inf_k) & ~((inf_r & inf_k) & (np.sign(rv) == np.sign(kv)))
            state = state.copy(); state.loc[rv.index[inf_mis]] = True
            fin = ~inf_r & ~inf_k
            diff = (rv - kv).abs().where(fin, 0.0); scale = np.maximum(rv.abs(), kv.abs()).where(fin, 0.0)
            over = (diff > a.atol + a.rtol * scale) & fin
            rel = (diff / np.maximum(scale, 1e-300))
            if len(rel) and rel.max() > worst[0]:
                worst = (float(rel.max()), c, rel.idxmax())
            is_waived = c in waived
            tie_mask = pd.Series([(g, base) in ties for g in rv.index], index=rv.index) if (suf in RANK and ties) else pd.Series(False, index=rv.index)
            resid_mask = pd.Series([(g, base) in resid for g in rv.index], index=rv.index) if (suf is not None and resid) else pd.Series(False, index=rv.index)
            # a cell of a residual or tie pair is NOT exempt: it must equal the statistic recomputed from its own
            # side's annual series (within rtol); a cell that does not is a gated mismatch
            exempt_ok = pd.Series(False, index=rv.index)
            for g in rv.index[(tie_mask | resid_mask) & over]:
                rec_ = recomputed.get((g, base))
                if not rec_ or c not in rec_[0] or c not in rec_[1]:
                    continue
                er, ec = rec_[0][c], rec_[1][c]
                ok_r = (np.isnan(er) and np.isnan(rv[g])) or (np.isfinite(er) and abs(er - rv[g]) <= a.atol + a.rtol * max(abs(er), abs(rv[g])))
                ok_c = (np.isnan(ec) and np.isnan(kv[g])) or (np.isfinite(ec) and abs(ec - kv[g]) <= a.atol + a.rtol * max(abs(ec), abs(kv[g])))
                exempt_ok[g] = bool(ok_r and ok_c)
            sz_mask = pd.Series([(g, base) in szero for g in rv.index], index=rv.index) if (suf in RANK and szero) else pd.Series(False, index=rv.index)
            gated_over = over & ~((tie_mask | resid_mask) & exempt_ok) & ~sz_mask
            n_gated += int(both.sum()); n_state += int(state.sum()); n_tie += int((over & tie_mask & ~resid_mask & exempt_ok).sum()); n_resid += int((over & resid_mask & exempt_ok).sum())
            n_szero += int((over & sz_mask & ~((tie_mask | resid_mask) & exempt_ok)).sum())
            # a gate-state difference (NaN vs finite) is NEVER waived by column or residual class; the one named exception
            # (--allow-gate-state-column) is a cell of a tie-bearing pair whose two values are each reproduced by the
            # statistics recomputed from that side's own annual series — verified, printed and counted, never assumed
            state_ok = pd.Series(False, index=r.index)
            if state.any() and ties and any(re.search(p_, c) for p_ in a.allow_gate_state_column):
                for g in r.index[state]:
                    rec_ = recomputed.get((g, base)) if (g, base) in ties else None
                    if not rec_ or c not in rec_[0] or c not in rec_[1]:
                        continue
                    er, ec = rec_[0][c], rec_[1][c]
                    rg = np.nan if pd.isna(r.get(g)) else float(r.get(g)); kg = np.nan if pd.isna(k.get(g)) else float(k.get(g))
                    er = np.nan if er is None or pd.isna(er) else float(er); ec = np.nan if ec is None or pd.isna(ec) else float(ec)
                    ok_r = (np.isnan(er) and np.isnan(rg)) or (np.isfinite(er) and np.isfinite(rg) and abs(er - rg) <= a.atol + a.rtol * max(abs(er), abs(rg)))
                    ok_c = (np.isnan(ec) and np.isnan(kg)) or (np.isfinite(ec) and np.isfinite(kg) and abs(ec - kg) <= a.atol + a.rtol * max(abs(ec), abs(kg)))
                    state_ok[g] = bool(ok_r and ok_c)
            n_state_ok += int(state_ok.sum())
            for g in list(r.index[state_ok])[:3]:
                if len(shown) < a.max_show:
                    shown.append((c, g, r.get(g), k.get(g), "gate-state of a tie pair, both sides reproduced by recomputation (named allowance)"))
            state_bad = state & ~state_ok
            if state_bad.any():
                n_bad += int(state_bad.sum()); bad_cols[c] = bad_cols.get(c, 0) + int(state_bad.sum())
                for g in list(r.index[state_bad])[:2]:
                    if len(shown) < a.max_show:
                        shown.append((c, g, r.get(g), k.get(g), "gate-state"))
            if is_waived:
                n_waived += int(over.sum())
                continue
            nb = int(gated_over.sum())
            if nb:
                n_bad += nb; bad_cols[c] = nb
                for g in list(rv.index[gated_over])[:2] + list(r.index[state_bad])[:1]:
                    if len(shown) < a.max_show:
                        shown.append((c, g, r.get(g), k.get(g)))
        else:
            eq = (r.isna() & k.isna()) | (r.astype(str) == k.astype(str))
            n_gated += len(r)
            if (~eq).any():
                if c in waived:
                    n_waived += int((~eq).sum()); continue
                n_bad += int((~eq).sum()); bad_cols[c] = int((~eq).sum())
                for g in list(r.index[~eq])[:2]:
                    if len(shown) < a.max_show:
                        shown.append((c, g, r.get(g), k.get(g)))
    print(f"compared {len(common)} gages x {len(cols)} columns ({n_gated:,} numeric cells); worst relative difference {worst[0]:.3e} at {worst[1]} / {worst[2]}")
    print(f"gated mismatches {n_bad:,} in {len(bad_cols)} columns (gate-state differences {n_state:,}); waived by column {n_waived:,}; "
          f"tie-sensitive rank statistics verified by recomputation: {n_tie:,}; propagated annual residuals verified by recomputation: {n_resid:,}; "
          f"signed-zero tie rank statistics (Known Issue, reported): {n_szero:,}; gate-state cells of tie pairs verified by recomputation (named allowance): {n_state_ok:,}")
    for s_ in shown:
        print("  MISMATCH", s_)
    if bad_cols:
        print("  columns:", sorted(bad_cols.items(), key=lambda kv: -kv[1])[:30])
    ok = ok and n_bad == 0
    import hashlib
    def _sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 24), b""):
                h.update(blk)
        return h.hexdigest()
    gl = (a.out + ".gages.txt") if a.out else None
    if gl:
        open(gl, "w").write("\n".join(common) + "\n")
    rep = {"ref": a.ref, "cand": a.cand, "ref_sha256": _sha(a.ref), "cand_sha256": _sha(a.cand), "waive_column": a.waive_column,
           "gages_compared_file": gl, "n_gages_compared": len(common), "n_cells_compared": n_gated,
           "gages": len(common), "columns": len(cols), "rtol": a.rtol, "atol": a.atol, "gated_mismatches": n_bad,
           "gate_state_differences": n_state, "gate_state_tie_verified": n_state_ok, "allow_gate_state_column": a.allow_gate_state_column, "waived": n_waived, "tie_sensitive_reported": n_tie, "annual_residual_reported": n_resid, "signed_zero_tie_reported": n_szero, "worst_rel": worst[0], "worst_col": worst[1],
           "bad_columns": bad_cols, "verdict": "PASS" if ok else "FAIL"}
    if a.out:
        json.dump(rep, open(a.out, "w"), indent=1, default=str)
    print("SUMMARY GATE", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
