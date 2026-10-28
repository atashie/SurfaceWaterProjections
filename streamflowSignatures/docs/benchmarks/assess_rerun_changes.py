#!/usr/bin/env python3
"""assess_rerun_changes.py — the change census of the qualifier rerun: every difference between
two runs (a stage pair, or the product vs version 2) at the annual and the summary layer,
classified against the dependency registry and the input ledgers, with the gates of
docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §5 E5 and §10.5 (5a–5g) reported as
PASS / FAIL / NOT RUN. The verdict is PASS only when every required gate PASSED; a gate that
could not run (an input not given) makes the verdict INCOMPLETE (exit 2), never PASS.
2026-10-26, reworked after the Codex tooling review.

    python docs/benchmarks/assess_rerun_changes.py --old-csv OLD.csv --old-annual OLD_annual.parquet
        --new-csv NEW.csv --new-annual NEW_annual.parquet --window 1993-2025 --label S2_vs_v2 --out-dir DIR
        [--registry docs/qualifier_rerun_registry.csv] [--config config/signatures_config.json]
        [--ledger correction_ledger.csv] [--inserted inserted_rows.csv] [--oracle acceptance_oracle.csv]
        [--state-a state_S0.csv --state-b state_S2.csv --state-changes state_changes_S0_S2.csv]
        [--gage-a gage_S0.csv --gage-b gage_S2.csv] [--delta-ledger daily_delta_ledger_S0_S2.csv]
        [--daily-a daily_S0.parquet --daily-b daily_S2.parquet] [--e7-annual-json J] [--e7-summary-json J]
        [--expect-identical LIST] [--gages LIST] [--memory-limit 3GB]

Classes per gage inside the window (plan §5 E5): A-Q = a restored splice-ledger day in a
water year of the window, or a state-table change the splice caused; A-J = an inserted climate
row in a water year that is Q-valid in either run, or a state-table change the join caused;
A-QJ = both; U = neither. Q-valid years come from the state tables (q_valid) when given, else
from the `Qann` keys of the annual parquets (the always-emitted year-local base).

Evidence (an explanation without it is UNEXPLAINED): a changed Q-only value needs a restored day
or a splice-caused state change in that year (or a record-dependent base of a splice-affected
gage, validated by 5c); a changed climate-only or coupled value needs a changed daily value of a
variable the base depends on in that year (the delta ledger; PPT for prcp and the coupled bases,
TMIN/TMAX for temp / gdd / frost, SWE for the snow bases) or a record-dependent base of an
affected gage; a gained or withdrawn key on an old year needs the emission gate's counts to
cross the threshold in the right direction (recession: Σ n_recession_events vs min_events; Q–P
seasonality: usable years vs min_years) or, for a climate base, the variable's validity flip
in the state tables; a gained year needs the oracle (splice) or a validity flip (join).
Gates: freeze (header, gage set, immutable site attributes, total_possible), U identical at both
layers, no lost year / key, 5a (0 unexplained changes; Q-only summary columns identical on
join-only gages), 5b (prcp_annual and the four seasonal totals equal the scoped sums of the
daily PPT deltas on EVERY shared eligible year, unchanged cells included; gdd_cum_* equal the
thresholded daily deltas — needs the daily series), 5c (the E7 annual gate on the affected
gages: an absolute recomputation in a second language), 5d (gained years = the oracle's; join
eligibility gains = validity flips), 5e (expected-identical gages identical), 5f (valid-year
identity = the oracle's final sets, derived metadata = oracle), 5g (the E7 summary gate),
registry coverage. Needs duckdb + pyarrow + pandas (~/HISSS_data/venvs/daymet/bin/python).
"""
import argparse
import datetime as dt
import json
import os
import sys

import duckdb
import numpy as np
import pandas as pd

SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval", "_mean", "_median",
            "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean", "_pettitt_delta_mean",
            "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
_SUF = sorted(SUFFIXES, key=len, reverse=True)
METADATA = ["gage_id", "latitude", "longitude", "basin_area", "gage_type", "num_water_years", "start_water_year",
            "end_water_year", "area_normalized", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009", "IMPNLCD06",
            "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED", "human_interference_class"]
IMMUTABLE = [c for c in METADATA if c not in ("gage_id", "num_water_years", "start_water_year", "end_water_year")]
PVAL_SUF = ("_spearman_pval", "_mk_pval", "_pettitt_pval", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval")
PPT_TOTALS = {"prcp_annual": None, "prcp_win": "win", "prcp_spr": "spr", "prcp_sum": "sum", "prcp_fal": "fal"}
GDD_TOTALS = {"gdd_cum_annual": None, "gdd_cum_win": "win", "gdd_cum_spr": "spr", "gdd_cum_sum": "sum", "gdd_cum_fal": "fal"}
SEASON_KEY = {"win": "winter", "spr": "spring", "sum": "summer", "fal": "fall"}
ANNUAL_TOL = (1e-9, 1e-9)


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def split_col(c):
    for s in _SUF:
        if c.endswith(s):
            return c[:-len(s)], s
    return c, None


def wy_of_dates(s):
    d = pd.to_datetime(s)
    return np.where(d.dt.month >= 10, d.dt.year + 1, d.dt.year)


def load_csv(path, gages=None):
    df = pd.read_csv(path, dtype={"gage_id": str}, low_memory=False)
    if gages is not None:
        df = df[df["gage_id"].isin(gages)]
    return df.set_index("gage_id").sort_index()


def cells_equal(a, b):
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return (a.isna() & b.isna()) | (a == b)
    return (a.isna() & b.isna()) | (a.astype(str) == b.astype(str))


def dep_vars(base, cls):
    """The daily variables an annual base depends on (evidence for a changed value)."""
    if base.startswith(("prcp_",)):
        return {"PPT"}
    if base.startswith(("temp_", "gdd_", "frost_")):
        return {"TMIN", "TMAX"}
    if base.startswith(("swe_", "melt_", "snow_")) or base == "ssm":
        return {"SWE", "PPT"} if base == "swe_max_to_ppt" else {"SWE"}
    if cls == "coupled":
        return {"Q", "PPT"}
    return {"Q"}


def region_of(g):
    return ("HUC" + g[:2]) if g.isdigit() else ("prov_" + g[:2])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-csv", required=True); ap.add_argument("--old-annual", required=True)
    ap.add_argument("--new-csv", required=True); ap.add_argument("--new-annual", required=True)
    ap.add_argument("--window", required=True); ap.add_argument("--label", required=True)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--registry", default=os.path.join(here, "..", "qualifier_rerun_registry.csv"))
    ap.add_argument("--config", default=os.path.join(here, "..", "..", "config", "signatures_config.json"))
    ap.add_argument("--ledger"); ap.add_argument("--inserted"); ap.add_argument("--oracle")
    ap.add_argument("--state-a"); ap.add_argument("--state-b"); ap.add_argument("--state-changes")
    ap.add_argument("--gage-a"); ap.add_argument("--gage-b"); ap.add_argument("--delta-ledger")
    ap.add_argument("--daily-a"); ap.add_argument("--daily-b")
    ap.add_argument("--e7-annual-json"); ap.add_argument("--e7-summary-json"); ap.add_argument("--expect-identical")
    ap.add_argument("--gages", help="restrict both runs to these gages (executed subset)")
    ap.add_argument("--accept-waiver", action="append", help="a named E7 residual the verdict JSON may carry (signature, prefix, column regex, or allow_value_diff)")
    ap.add_argument("--out-dir", required=True); ap.add_argument("--memory-limit", default="3GB"); ap.add_argument("--max-list", type=int, default=25)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    w0, w1 = (int(x) for x in a.window.split("-"))
    tag = f"{a.label}_{a.window}"
    cfg = json.load(open(a.config))
    season_months = {k: set(v) for k, v in cfg["na_handling"]["seasonal_completeness"]["season_definitions"].items()}
    gdd_base = float(cfg["climate"]["gdd_base_c"]); rec_min = int(cfg["recession"]["min_events"]); qp_min = int(cfg["qp_seasonality"]["min_years"])
    gates = {}; notes = []
    lines = [f"# Change census {a.label}, window {a.window}", "",
             f"Generated {dt.datetime.now().isoformat(timespec='seconds')}; old `{a.old_csv}`, new `{a.new_csv}`.", ""]

    def gate(name, status, detail):
        gates[name] = status; lines.append(f"- **{name}: {status}** — {detail}"); log(f"{name}: {status} — {detail}")

    # ---------------- registry and summary CSVs ----------------
    reg = pd.read_csv(a.registry, comment="#", dtype=str).set_index("signature")
    cls = reg["climate_class"].to_dict(); dep = reg["annual_dependence"].to_dict(); egate = reg["emission_gate"].fillna("").to_dict()
    gages_sel = set(l.strip() for l in open(a.gages) if l.strip()) if a.gages else None
    old = load_csv(a.old_csv, gages_sel); new = load_csv(a.new_csv, gages_sel)
    gate("freeze: header identical", "PASS" if list(old.columns) == list(new.columns) else "FAIL",
         f"old {old.shape[1]} columns, new {new.shape[1]}; only-old {sorted(set(old.columns) - set(new.columns))[:5]}, only-new {sorted(set(new.columns) - set(old.columns))[:5]}")
    gate("freeze: gage set identical", "PASS" if set(old.index) == set(new.index) else "FAIL",
         f"old {len(old)} gages, new {len(new)}; only-old {sorted(set(old.index) - set(new.index))[:5]}, only-new {sorted(set(new.index) - set(old.index))[:5]}")
    common = sorted(set(old.index) & set(new.index)); cols = [c for c in old.columns if c in new.columns]
    old = old.loc[common, cols]; new = new.loc[common, cols]
    imm = [c for c in IMMUTABLE if c in cols]
    imm_bad = {c: int((~cells_equal(old[c], new[c])).sum()) for c in imm if (~cells_equal(old[c], new[c])).any()}
    gate("freeze: immutable site attributes equal", "PASS" if not imm_bad else "FAIL", f"{len(imm)} columns compared; differing {imm_bad}")
    sig_cols = [c for c in cols if c not in METADATA and not c.startswith("flagged_")]
    bases = {split_col(c)[0] for c in sig_cols}
    not_in_reg = sorted(b for b in bases if b not in reg.index)
    gate("registry coverage", "PASS" if not not_in_reg else "FAIL", f"{len(bases)} bases / scalars in the header; not in the registry {len(not_in_reg)} {not_in_reg[:10]}")

    # ---------------- exposure per gage ----------------
    restored_years, splice_state_years, join_state_years = {}, {}, {}
    if a.ledger:
        led = pd.read_csv(a.ledger, dtype={"gage_id": str, "date": str}, float_precision="round_trip"); led = led[led["status"] == "restored"]
        led["wy"] = wy_of_dates(led["date"]); led = led[(led["wy"] >= w0) & (led["wy"] <= w1)]
        for g, d in led.groupby("gage_id"):
            restored_years[g] = set(d["wy"].astype(int))
    inserted_years, neighbour_years = {}, {}
    if a.inserted:
        inv = pd.read_csv(a.inserted, dtype={"gage_id": str})
        inv = inv[(inv["water_year"] >= w0) & (inv["water_year"] <= w1)]
        for g, d in inv[inv["kind"] == "inserted"].groupby("gage_id"):
            inserted_years[g] = set(d["water_year"].astype(int))
        for g, d in inv[inv["kind"] != "inserted"].groupby("gage_id"):
            neighbour_years[g] = set(d["water_year"].astype(int))
    if a.state_changes:
        scg = pd.read_csv(a.state_changes, dtype=str, keep_default_na=False)
        scg["wy"] = scg["water_year"].astype(float).astype(int); scg = scg[(scg["wy"] >= w0) & (scg["wy"] <= w1)]
        for g, d in scg[scg["cause"].str.startswith("splice")].groupby("gage_id"):
            splice_state_years[g] = set(d["wy"])
        for g, d in scg[scg["cause"].str.startswith("inserted")].groupby("gage_id"):
            join_state_years[g] = set(d["wy"])
    state_a = state_b = None
    if a.state_a and a.state_b:
        state_a = pd.read_csv(a.state_a, dtype=str, keep_default_na=False); state_b = pd.read_csv(a.state_b, dtype=str, keep_default_na=False)
        for s in (state_a, state_b):
            s["wy"] = s["water_year"].astype(float).astype(int)
    oracle = None; oracle_gained = {}
    if a.oracle:
        orc = pd.read_csv(a.oracle, dtype={"gage_id": str}); orc = orc[orc["window"] == a.window].set_index("gage_id")
        oracle = orc
        for g in orc.index:
            gy = orc.at[g, "gained_years"]
            oracle_gained[g] = set(int(y) for y in str(gy).split()) if isinstance(gy, str) and gy.strip() else set()
    delta = None; delta_by_year = {}
    if a.delta_ledger:
        delta = pd.read_csv(a.delta_ledger, dtype={"gage_id": str, "date": str}, float_precision="round_trip")
        delta["wy"] = delta["water_year"].astype(int); delta = delta[(delta["wy"] >= w0) & (delta["wy"] <= w1)]
        for (g, y), d in delta.groupby(["gage_id", "wy"]):
            delta_by_year[(g, int(y))] = set(d["variable"])

    # ---------------- annual layer: keyed diff in duckdb ----------------
    con = duckdb.connect()
    con.execute(f"SET memory_limit='{a.memory_limit}'; SET threads=4; SET temp_directory='{os.path.join(a.out_dir, '.duckdb_tmp')}'")
    gfilter = ""
    if gages_sel is not None:
        con.execute("CREATE TEMP TABLE gsel AS SELECT unnest(?::VARCHAR[]) AS gage_id", [sorted(gages_sel)]); gfilter = " AND gage_id IN (SELECT gage_id FROM gsel)"
    V = '"value"'
    for name, path in (("o", a.old_annual), ("n", a.new_annual)):
        con.execute(f"CREATE TEMP VIEW {name} AS SELECT CAST(gage_id AS VARCHAR) gage_id, CAST(signature AS VARCHAR) signature, "
                    f"CAST(water_year AS INTEGER) water_year, CAST({V} AS DOUBLE) AS {V} "
                    f"FROM read_parquet('{path}') WHERE water_year BETWEEN {w0} AND {w1}{gfilter}")
    log("annual layer: full outer join of the two parquets ...")
    con.execute("""CREATE TEMP TABLE d AS
        SELECT coalesce(o.gage_id, n.gage_id) gage_id, coalesce(o.signature, n.signature) signature,
               coalesce(o.water_year, n.water_year) water_year, o."value" AS old_value, n."value" AS new_value,
               CASE WHEN o.gage_id IS NULL THEN 'gained' WHEN n.gage_id IS NULL THEN 'lost'
                    WHEN ((o."value" IS NULL OR isnan(o."value")) AND (n."value" IS NULL OR isnan(n."value"))) OR coalesce(o."value" = n."value", false) THEN 'equal'
                    WHEN (o."value" IS NULL OR isnan(o."value")) <> (n."value" IS NULL OR isnan(n."value")) THEN 'nan_state' ELSE 'changed' END AS status
        FROM o FULL OUTER JOIN n ON o.gage_id = n.gage_id AND o.signature = n.signature AND o.water_year = n.water_year""")
    tot = dict(con.execute("SELECT status, count(*) FROM d GROUP BY status").fetchall())
    n_old = con.execute("SELECT count(*) FROM o").fetchone()[0]; n_new = con.execute("SELECT count(*) FROM n").fetchone()[0]
    log(f"annual rows: old {n_old:,}, new {n_new:,}; {tot}")
    diff = con.execute("SELECT * FROM d WHERE status <> 'equal'").df()
    # Q-valid years per gage: the state tables, else the Qann keys
    valid_old = {g: set(int(y) for y in ys) for g, ys in con.execute("SELECT gage_id, list(DISTINCT water_year) FROM o WHERE signature = 'Qann' GROUP BY gage_id").fetchall()}
    valid_new = {g: set(int(y) for y in ys) for g, ys in con.execute("SELECT gage_id, list(DISTINCT water_year) FROM n WHERE signature = 'Qann' GROUP BY gage_id").fetchall()}
    valid_source = "Qann annual keys"; not_replayed = []
    if state_a is not None:
        # the state tables cover the replayed (affected) gages only (B′ runs on the affected manifest, 2026-10-27);
        # a gage absent from them keeps its Qann keys
        so = {g: set(int(y) for y in d.loc[d["q_valid"] == "True", "wy"]) for g, d in state_a[(state_a["wy"] >= w0) & (state_a["wy"] <= w1)].groupby("gage_id")}
        sn = {g: set(int(y) for y in d.loc[d["q_valid"] == "True", "wy"]) for g, d in state_b[(state_b["wy"] >= w0) & (state_b["wy"] <= w1)].groupby("gage_id")}
        valid_old.update(so); valid_new.update(sn)
        not_replayed = [g for g in common if g not in so or g not in sn]
        valid_source = f"state tables for {len(set(so) & set(sn))} replayed gages, Qann annual keys for the {len(not_replayed)} others"
    rec_old = {g: v for g, v in con.execute("SELECT gage_id, sum(\"value\") FROM o WHERE signature = 'n_recession_events' AND NOT isnan(\"value\") GROUP BY gage_id").fetchall()}
    rec_new = {g: v for g, v in con.execute("SELECT gage_id, sum(\"value\") FROM n WHERE signature = 'n_recession_events' AND NOT isnan(\"value\") GROUP BY gage_id").fetchall()}
    qp_old = {g: v for g, v in con.execute("SELECT gage_id, count(DISTINCT water_year) FROM o WHERE signature = 'qp_bimodality' AND NOT isnan(\"value\") GROUP BY gage_id").fetchall()}
    qp_new = {g: v for g, v in con.execute("SELECT gage_id, count(DISTINCT water_year) FROM n WHERE signature = 'qp_bimodality' AND NOT isnan(\"value\") GROUP BY gage_id").fetchall()}
    # variable validity per (gage, year) from the state tables
    def validity(state, g, y, var):
        if state is None:
            return None
        r = state[(state["gage_id"] == g) & (state["wy"] == y)]
        if r.empty:
            return None
        k = {"Q": "q_valid", "PPT": "ppt_valid", "SWE": "swe_valid", "TMIN": "temp_valid", "TMAX": "temp_valid"}[var]
        return r.iloc[0][k] == "True"
    va_cache = {}
    def valid_var(state, which, g, y, var):
        key = (which, g, y, var)
        if key not in va_cache:
            va_cache[key] = validity(state, g, y, var)
        return va_cache[key]

    # ---------------- classes ----------------
    rows = []
    for g in common:
        vo, vn = valid_old.get(g, set()), valid_new.get(g, set())
        aq = bool(restored_years.get(g)) or bool(splice_state_years.get(g))
        aj = bool(inserted_years.get(g, set()) & (vo | vn)) or bool(join_state_years.get(g, set()) & (vo | vn))
        c = "A-QJ" if (aq and aj) else "A-Q" if aq else "A-J" if aj else "U"
        rows.append({"gage_id": g, "class": c, "restored_years": " ".join(map(str, sorted(restored_years.get(g, set())))),
                     "inserted_valid_years": " ".join(map(str, sorted(inserted_years.get(g, set()) & (vo | vn)))),
                     "n_valid_old": len(vo), "n_valid_new": len(vn), "gained_years": " ".join(map(str, sorted(vn - vo))),
                     "lost_years": " ".join(map(str, sorted(vo - vn))), "oracle_gained": " ".join(map(str, sorted(oracle_gained.get(g, set())))) if oracle is not None else ""})
    classes = pd.DataFrame(rows).set_index("gage_id")
    ccount = classes["class"].value_counts().to_dict()
    lines += ["", f"Classes: {ccount} (Q-valid years from the {valid_source}).", ""]
    classes.to_csv(os.path.join(a.out_dir, f"classes_{tag}.csv"))

    # ---------------- annual classification with evidence ----------------
    def explain(r):
        g, sig, wy, st = r["gage_id"], r["signature"], int(r["water_year"]), r["status"]
        if sig not in reg.index:
            return "UNEXPLAINED: signature not in the registry"
        c = classes.at[g, "class"]; k = cls[sig]; d_ = dep[sig]; eg = egate[sig]
        vo = valid_old.get(g, set()); vn = valid_new.get(g, set())
        rest = wy in restored_years.get(g, set()) or wy in splice_state_years.get(g, set())
        dv = dep_vars(sig, k)
        changed_vars = delta_by_year.get((g, wy), set()) if delta is not None else None
        if c == "U":
            return "UNEXPLAINED: class U gage"
        if st == "lost":
            if wy not in vn and wy in vo:
                return "UNEXPLAINED: lost year"
            if eg.startswith("recession") and rec_old.get(g, 0) >= rec_min > rec_new.get(g, 0):
                return "gate-emitted old year withdrawn (recession min_events crossed downward)"
            if eg.startswith("qp") and qp_old.get(g, 0) >= qp_min > qp_new.get(g, 0):
                return "gate-emitted old year withdrawn (Q–P seasonality min_years crossed downward)"
            if k in ("Q_only", "coupled") and rest:
                return f"restored-day year ({'Q-only' if k == 'Q_only' else 'coupled'}; key withdrawn on an old year)"
            if d_ == "record_dependent" and (c in ("A-Q", "A-QJ") if k == "Q_only" else c != "U"):
                return "record-dependent base of an affected gage: key withdrawn on an old year (validated by the E7 key-set identity of 5c)"
            return "UNEXPLAINED: lost key"
        if st == "gained" and wy not in vo:
            if wy in oracle_gained.get(g, set()):
                return "newly eligible year (restored; oracle)"
            if k in ("climate_only", "coupled") and wy in vn and wy in vo:
                pass
            if k == "Q_only" and oracle is None and rest:
                return "newly eligible year (restored day; no oracle given)"
            return "UNEXPLAINED: gained year not in the oracle"
        if st == "gained":   # a new key on an old year
            if eg.startswith("recession") and rec_old.get(g, 0) < rec_min <= rec_new.get(g, 0):
                return "gate-emitted old year (recession min_events crossed upward)"
            if eg.startswith("qp") and qp_old.get(g, 0) < qp_min <= qp_new.get(g, 0):
                return "gate-emitted old year (Q–P seasonality min_years crossed upward)"
            if k in ("climate_only", "coupled"):
                for var in dv - {"Q"}:
                    if valid_var(state_a, "a", g, wy, var) is False and valid_var(state_b, "b", g, wy, var) is True:
                        return f"newly eligible climate year ({var} validity flipped; state tables)"
            if k in ("Q_only", "coupled") and rest:
                return f"restored-day year ({'Q-only' if k == 'Q_only' else 'coupled'}; key emitted on an old year)"
            if d_ == "record_dependent" and (c in ("A-Q", "A-QJ") if k == "Q_only" else c != "U"):
                return "record-dependent base of an affected gage: key emitted on an old year (validated by the E7 key-set identity and value recomputation of 5c)"
            if k in ("climate_only", "coupled") and state_a is None:
                return "UNEXPLAINED: gained key on an old year (no state tables to show a validity flip)"
            return "UNEXPLAINED: gained key on an old year"
        # changed / nan_state on a shared key
        if k == "Q_only":
            if rest:
                return "restored-day year (Q-only)"
            if d_ == "record_dependent" and c in ("A-Q", "A-QJ"):
                return "record-dependent base of a splice-affected gage (validated by 5c)"
            return "UNEXPLAINED: Q-only change" + (" on a join-only gage" if c == "A-J" else " on a year without a restored day")
        if changed_vars is None:
            return "EVIDENCE MISSING: climate-dependent change without the delta ledger"
        hit = sorted(changed_vars & dv)
        if hit:
            return f"changed daily {'/'.join(hit)} in the year (delta ledger)"
        if k == "coupled" and rest:
            return "restored-day year (coupled; Q delta)"
        if d_ == "record_dependent" and c != "U":
            return "record-dependent base of an affected gage (validated by 5c)"
        if k == "climate_only" and rest:
            return "UNEXPLAINED: climate-only value moved on a restored-day year without a climate delta (the splice moved the climate frame?)"
        return "UNEXPLAINED: no daily change of a dependent variable in the year"

    if len(diff):
        diff["explanation"] = diff.apply(explain, axis=1)
        diff["climate_class"] = diff["signature"].map(cls); diff["gage_class"] = diff["gage_id"].map(classes["class"])
        diff["abs_delta"] = (diff["new_value"] - diff["old_value"]).abs()
        # a changed cell within the annual tolerance (|Δ| <= 1e-12 + 1e-9·max) is a NUMERICAL difference (summation
        # order, last bits), never a hydrological effect; the report keeps the two apart (final review, 2026-10-26)
        scale = np.maximum(diff["new_value"].abs(), diff["old_value"].abs())
        diff["magnitude"] = np.where(diff["status"] != "changed", diff["status"], np.where(diff["abs_delta"] <= ANNUAL_TOL[0] + ANNUAL_TOL[1] * scale, "numerical", "substantive"))
    else:
        for c_ in ("explanation", "climate_class", "gage_class", "abs_delta", "magnitude"):
            diff[c_] = []
    diff.to_csv(os.path.join(a.out_dir, f"annual_changes_{tag}.csv"), index=False)
    unexpl = diff[diff["explanation"].str.startswith("UNEXPLAINED")]; missing_ev = diff[diff["explanation"].str.startswith("EVIDENCE MISSING")]
    unexpl.to_csv(os.path.join(a.out_dir, f"annual_unexplained_{tag}.csv"), index=False)
    expl_counts = diff["explanation"].value_counts().to_dict()
    mag = diff["magnitude"].value_counts().to_dict() if len(diff) else {}
    lines += ["## Annual layer", "", f"old rows {n_old:,}, new rows {n_new:,}; equal {tot.get('equal', 0):,}, changed {tot.get('changed', 0):,} "
              f"(substantive {mag.get('substantive', 0):,}, numerical-only within the annual tolerance {mag.get('numerical', 0):,}), "
              f"NaN-state {tot.get('nan_state', 0):,}, gained {tot.get('gained', 0):,}, lost {tot.get('lost', 0):,}.", "",
              "| explanation | cells |", "|---|---|"] + [f"| {k} | {v:,} |" for k, v in sorted(expl_counts.items(), key=lambda kv: -kv[1])] + [""]
    u_annual = diff[diff["gage_class"] == "U"]
    gate("class U identical (annual)", "PASS" if len(u_annual) == 0 else "FAIL", f"{len(u_annual)} non-equal annual cells on U gages ({u_annual['gage_id'].nunique()} gages)")
    gate("no lost year", "PASS" if not any(classes["lost_years"] != "") else "FAIL", f"gages losing a Q-valid year: {int((classes['lost_years'] != '').sum())}")
    if len(missing_ev):
        gate("5a: 0 unexplained annual changes", "NOT RUN", f"{len(missing_ev):,} climate-dependent changes need the delta ledger (--delta-ledger); unexplained so far {len(unexpl):,}")
    else:
        gate("5a: 0 unexplained annual changes", "PASS" if len(unexpl) == 0 else "FAIL",
             f"{len(unexpl):,} unexplained of {len(diff):,} non-equal cells; by reason {unexpl['explanation'].value_counts().to_dict() if len(unexpl) else {}}")
    # 5d / 5f: gained years and the valid-year identity against the oracle
    if oracle is not None:
        bad_d, bad_f = [], []
        for g in common:
            vo, vn = valid_old.get(g, set()), valid_new.get(g, set())
            og = oracle_gained.get(g, set())
            if (vn - vo) != og:
                bad_d.append((g, sorted(vn - vo), sorted(og)))
            if g in oracle.index:
                fin = set(int(y) for y in str(oracle.at[g, "final_valid_years"]).split()) if isinstance(oracle.at[g, "final_valid_years"], str) else set()
                base = set(int(y) for y in str(oracle.at[g, "baseline_valid_years"]).split()) if isinstance(oracle.at[g, "baseline_valid_years"], str) else set()
                nw = int(new.at[g, "num_water_years"]) if "num_water_years" in cols else None
                sw = int(new.at[g, "start_water_year"]) if "start_water_year" in cols else None; ew = int(new.at[g, "end_water_year"]) if "end_water_year" in cols else None
                if vn != fin or vo != base or (nw is not None and nw != len(fin)) or (sw is not None and fin and sw != min(fin)) or (ew is not None and fin and ew != max(fin)):
                    bad_f.append((g, len(vn), len(fin), nw, sw, ew))
            else:
                bad_f.append((g, "not in the oracle"))
        gate("5d: gained years = the oracle's gained years, gage for gage", "PASS" if not bad_d else "FAIL", f"{len(bad_d)} gages differ: {bad_d[:5]}")
        gate("5f: valid-year identity = the oracle (baseline and final sets exact; num/start/end_water_year = final)", "PASS" if not bad_f else "FAIL", f"{len(bad_f)} gages off: {bad_f[:5]}")
    else:
        gate("5d: gained years = the oracle's gained years, gage for gage", "NOT RUN", "no --oracle")
        gate("5f: valid-year identity = the oracle", "NOT RUN", "no --oracle")
    if a.gage_a and a.gage_b:
        ga = pd.read_csv(a.gage_a, dtype=str, keep_default_na=False).set_index("gage_id"); gb = pd.read_csv(a.gage_b, dtype=str, keep_default_na=False).set_index("gage_id")
        col = f"{a.window}_total_possible"
        missing_tp = [g for g in common if g not in ga.index or g not in gb.index]
        missing_not_u = [g for g in missing_tp if classes.at[g, "class"] != "U"]
        bad_tp = [(g, ga.at[g, col], gb.at[g, col]) for g in common if g in ga.index and g in gb.index and ga.at[g, col] != gb.at[g, col]]
        gate("freeze: total_possible identical (gage tables)", "PASS" if not bad_tp and not missing_not_u else "FAIL",
             f"{len(common) - len(missing_tp)} replayed gages compared; {len(bad_tp)} differ {bad_tp[:5]}; not replayed {len(missing_tp)}, "
             f"all class U: {not missing_not_u} (non-U absent {missing_not_u[:5]})")
    else:
        gate("freeze: total_possible identical (gage tables)", "NOT RUN", "no --gage-a/--gage-b")

    # ---------------- summary layer ----------------
    flag_cols = [c for c in cols if c.startswith("flagged_")]
    derived = [c for c in ("num_water_years", "start_water_year", "end_water_year") if c in cols]
    changed_bases = {(g, s) for g, s in zip(diff["gage_id"], diff["signature"])} if len(diff) else set()
    srows, gage_changes, sum_unexpl = [], pd.Series(0, index=common), []
    pop = []
    for c in sig_cols + flag_cols + derived:
        ao, an = old[c], new[c]
        numeric = pd.api.types.is_numeric_dtype(ao) and pd.api.types.is_numeric_dtype(an)
        ne = ~cells_equal(ao, an)
        gage_changes = gage_changes.add(ne.astype(int), fill_value=0)
        base, suf = split_col(c)
        kind = "flag" if c.startswith("flagged_") else "metadata" if c in derived else cls.get(base, "?")
        n_numerical = 0
        if numeric and ne.any():
            bo = ao.notna() & an.notna() & ne
            dd = (an[bo].astype(float) - ao[bo].astype(float)).abs(); sc = np.maximum(an[bo].astype(float).abs(), ao[bo].astype(float).abs())
            n_numerical = int((dd <= 1e-12 + 1e-9 * sc).sum())
        rec = {"column": c, "base": base, "suffix": suf, "climate_class": kind, "n_changed": int(ne.sum()), "n_changed_numerical_only": n_numerical,
               "n_gained": int((ao.isna() & an.notna()).sum()), "n_lost": int((ao.notna() & an.isna()).sum()),
               "n_changed_U": int(ne[classes["class"] == "U"].sum()), "n_changed_AJ": int(ne[classes["class"] == "A-J"].sum())}
        for g in ne.index[ne]:
            gc = classes.at[g, "class"]
            if gc == "U":
                sum_unexpl.append((c, g, "class U gage"))
            elif kind == "Q_only" and gc == "A-J":
                sum_unexpl.append((c, g, "Q-only column on a join-only gage"))
            elif kind in ("Q_only", "climate_only", "coupled") and suf is not None and (g, base) not in changed_bases:
                sum_unexpl.append((c, g, "summary changed with no annual change of the base for this gage"))
        if numeric and ne.any():
            both = ao.notna() & an.notna() & ne
            d_ = (an[both].astype(float) - ao[both].astype(float))
            if len(d_):
                rec["max_abs_delta"] = float(d_.abs().max()); rec["max_abs_delta_gage"] = d_.abs().idxmax()
                rec["p50_abs_delta"] = float(d_.abs().median()); rec["p99_abs_delta"] = float(d_.abs().quantile(0.99))
                nz = ao[both].astype(float) != 0
                rec["p99_rel_delta"] = float((d_[nz].abs() / ao[both][nz].astype(float).abs()).quantile(0.99)) if nz.any() else np.nan
                o_all = ao[ao.notna() & an.notna()].astype(float); n_all = an[ao.notna() & an.notna()].astype(float)
                ss = float(((o_all - o_all.mean()) ** 2).sum())
                rec["identity_r2"] = float(1 - ((n_all - o_all) ** 2).sum() / ss) if ss > 0 else np.nan
                sd = float(o_all.std()) if len(o_all) > 1 else np.nan
                rec["max_std_delta"] = float(d_.abs().max() / sd) if sd and sd > 0 else np.nan
                if suf in PVAL_SUF:
                    rec["cross_0.05"] = int(((ao[both] < 0.05) != (an[both] < 0.05)).sum()); rec["cross_0.01"] = int(((ao[both] < 0.01) != (an[both] < 0.01)).sum())
                if suf == "_pettitt_cp_year":
                    rec["cp_year_moved"] = int((d_ != 0).sum()); rec["cp_shift_p50"] = float(d_[d_ != 0].abs().median()) if (d_ != 0).any() else 0.0
                if suf in ("_senn_slp", "_linear_slp"):
                    rec["sign_changes"] = int((np.sign(ao[both].astype(float)) != np.sign(an[both].astype(float))).sum())
                if "dowy" in base and suf in ("_mean", "_median"):
                    rec["days_moved_p50"] = float(d_.abs().median())
        if c.startswith("flagged_") and ne.any():
            rec["flips_true"] = int(((ao.astype(str).str.lower() != "true") & (an.astype(str).str.lower() == "true")).sum())
            rec["flips_false"] = int(((ao.astype(str).str.lower() == "true") & (an.astype(str).str.lower() != "true")).sum())
        if numeric and suf == "_mean":
            pop.append({"column": c, "class": kind, "pop_mean_old": float(ao.astype(float).mean()), "pop_mean_new": float(an.astype(float).mean()),
                        "pop_median_old": float(ao.astype(float).median()), "pop_median_new": float(an.astype(float).median()), "n_changed": int(ne.sum())})
        srows.append(rec)
    summ = pd.DataFrame(srows); summ.to_csv(os.path.join(a.out_dir, f"summary_columns_{tag}.csv"), index=False)
    pd.DataFrame(pop).to_csv(os.path.join(a.out_dir, f"population_means_{tag}.csv"), index=False)
    gage_changes.rename("n_changed_cells").to_frame().join(classes[["class"]]).sort_values("n_changed_cells", ascending=False).to_csv(
        os.path.join(a.out_dir, f"summary_changes_by_gage_{tag}.csv"))
    su = pd.DataFrame(sum_unexpl, columns=["column", "gage_id", "reason"]); su.to_csv(os.path.join(a.out_dir, f"summary_unexplained_{tag}.csv"), index=False)
    n_u = int(summ["n_changed_U"].sum()); n_aj_q = int(summ.loc[summ["climate_class"] == "Q_only", "n_changed_AJ"].sum())
    gate("class U identical (summary)", "PASS" if n_u == 0 else "FAIL", f"{n_u} changed summary cells on U gages; columns {summ.loc[summ['n_changed_U'] > 0, 'column'].tolist()[:a.max_list]}")
    gate("5a: Q-only summary columns identical on join-only gages", "PASS" if n_aj_q == 0 else "FAIL",
         f"{n_aj_q} changed Q-only summary cells on A-J gages; columns {summ.loc[(summ['climate_class'] == 'Q_only') & (summ['n_changed_AJ'] > 0), 'column'].tolist()[:a.max_list]}")
    gate("5a: every changed summary cell traces to an annual change of its base (or a flag / scalar of an affected gage)", "PASS" if len(su) == 0 else "FAIL",
         f"{len(su)} untraced summary cells {su['reason'].value_counts().to_dict() if len(su) else {}}")
    # regional and standardized views
    try:
        reg_rows = []
        for c in [c for c in sig_cols if c.endswith("_mean") and cls.get(split_col(c)[0]) == "climate_only"]:
            ao, an = old[c].astype(float), new[c].astype(float); ne = ~cells_equal(old[c], new[c]) & ao.notna() & an.notna()
            if ne.any():
                dd = (an - ao)[ne]
                for reg_, grp in dd.groupby([region_of(g) for g in dd.index]):
                    reg_rows.append({"column": c, "region": reg_, "n_changed": len(grp), "median_delta": float(grp.median()), "max_abs_delta": float(grp.abs().max())})
        pd.DataFrame(reg_rows).to_csv(os.path.join(a.out_dir, f"regional_climate_deltas_{tag}.csv"), index=False)
    except Exception as e:  # noqa: BLE001
        notes.append(f"regional view not written: {e}")
    chg = summ[summ["n_changed"] > 0].sort_values("n_changed", ascending=False)
    lines += ["## Summary layer", "", f"{len(sig_cols)} signature columns, {len(flag_cols)} flags; columns with any change {len(chg)} of {len(summ)}; "
              f"cells changed {int(summ['n_changed'].sum()):,} (of which numerical-only within 1e-9: {int(summ['n_changed_numerical_only'].sum()):,}; "
              f"gained {int(summ['n_gained'].sum()):,}, lost {int(summ['n_lost'].sum()):,}).", "",
              "| column | class | changed | gained | lost | max abs Δ (gage) | max Δ / sd | p99 rel Δ | R² |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in chg.head(a.max_list).iterrows():
        lines.append(f"| {r['column']} | {r['climate_class']} | {r['n_changed']} | {r['n_gained']} | {r['n_lost']} | "
                     f"{r.get('max_abs_delta', float('nan')):.4g} ({r.get('max_abs_delta_gage', '')}) | {r.get('max_std_delta', float('nan')):.3g} | "
                     f"{r.get('p99_rel_delta', float('nan')):.3g} | {r.get('identity_r2', float('nan')):.6f} |")
    pv = summ[summ["suffix"].isin(PVAL_SUF)]
    lines += ["", f"p-value columns: {len(pv)}; crossings of 0.05: {int(pv['cross_0.05'].sum()) if 'cross_0.05' in pv else 0}, of 0.01: "
              f"{int(pv['cross_0.01'].sum()) if 'cross_0.01' in pv else 0}; Pettitt cp_year moved: {int(summ['cp_year_moved'].sum()) if 'cp_year_moved' in summ else 0}; "
              f"slope sign changes: {int(summ['sign_changes'].sum()) if 'sign_changes' in summ else 0}; flag flips to TRUE "
              f"{int(summ['flips_true'].sum()) if 'flips_true' in summ else 0} / to FALSE {int(summ['flips_false'].sum()) if 'flips_false' in summ else 0}; "
              f"population means / medians of every `_mean` before vs after in population_means_{tag}.csv; regional climate deltas in regional_climate_deltas_{tag}.csv.", ""]

    # ---------------- 5b: totals by summed daily deltas, over EVERY shared eligible year ----------------
    if delta is not None:
        ppt = delta[delta["variable"] == "PPT"].copy(); ppt["d"] = ppt["new_value"].astype(float).fillna(0) - ppt["old_value"].astype(float).fillna(0)
        ppt["month"] = pd.to_datetime(ppt["date"]).dt.month
        shared = con.execute("SELECT gage_id, signature, water_year, old_value, new_value FROM d WHERE status IN ('equal', 'changed') AND signature IN "
                             "('prcp_annual', 'prcp_win', 'prcp_spr', 'prcp_sum', 'prcp_fal')").df()
        shared = shared[shared["gage_id"].isin(set(classes.index[classes["class"] != "U"]) | set(ppt["gage_id"]))]
        sums = {}
        for sig, season in PPT_TOTALS.items():
            sub = ppt if season is None else ppt[ppt["month"].isin(season_months[SEASON_KEY[season]])]
            sums[sig] = sub.groupby(["gage_id", "wy"])["d"].sum().to_dict()
        bad_b = []; n_cmp = 0
        for r in shared.itertuples(index=False):
            exp = sums[r.signature].get((r.gage_id, int(r.water_year)), 0.0)
            act = (0.0 if np.isnan(r.new_value) or np.isnan(r.old_value) else r.new_value - r.old_value)
            n_cmp += 1
            if abs(act - exp) > ANNUAL_TOL[0] + ANNUAL_TOL[1] * max(abs(exp), abs(act)):
                bad_b.append((r.gage_id, r.signature, int(r.water_year), act, exp))
        pd.DataFrame(bad_b, columns=["gage_id", "signature", "water_year", "actual_delta", "ledger_sum"]).to_csv(os.path.join(a.out_dir, f"gate5b_ppt_mismatches_{tag}.csv"), index=False)
        gate("5b: precipitation totals (annual + 4 seasonal) = summed daily PPT deltas on every shared eligible year", "PASS" if (not bad_b and (n_cmp > 0 or len(ppt) == 0)) else "FAIL",
             f"{n_cmp:,} shared (gage, total, year) cells of affected gages compared (unchanged cells included); off the ledger sum {len(bad_b)} {bad_b[:3]}")
        if a.daily_a and a.daily_b:
            tdates = delta[delta["variable"].isin(["TMIN", "TMAX"])][["gage_id", "date", "wy"]].drop_duplicates()
            if len(tdates) and len(delta[delta["variable"].isin(["TMIN", "TMAX"])]):
                con.execute("CREATE OR REPLACE TEMP TABLE td AS SELECT * FROM (VALUES " + ", ".join(f"('{g}', DATE '{d}')" for g, d in zip(tdates["gage_id"], tdates["date"])) + ") t(gage_id, date)")
                tv = con.execute(f"""SELECT A.gage_id, A.date, A.water_year, A.TMIN_pre tn_a, A.TMAX_pre tx_a, B.TMIN_pre tn_b, B.TMAX_pre tx_b
                                     FROM read_parquet('{a.daily_a}') A JOIN read_parquet('{a.daily_b}') B ON A.gage_id = B.gage_id AND A.date = B.date
                                     JOIN td ON td.gage_id = A.gage_id AND td.date = CAST(A.date AS DATE)""").df()
                tv["month"] = pd.to_datetime(tv["date"]).dt.month
                def gdd(tx, tn):
                    return 0.0 if (np.isnan(tx) or np.isnan(tn)) else max(0.0, (tx + tn) / 2 - gdd_base)
                tv["dg"] = [gdd(x, n) - gdd(xa, na) for xa, na, x, n in zip(tv["tx_a"], tv["tn_a"], tv["tx_b"], tv["tn_b"])]
                gsums = {}
                for sig, season in GDD_TOTALS.items():
                    sub = tv if season is None else tv[tv["month"].isin(season_months[SEASON_KEY[season]])]
                    gsums[sig] = sub.groupby(["gage_id", "water_year"])["dg"].sum().to_dict()
                sharedg = con.execute("SELECT gage_id, signature, water_year, old_value, new_value FROM d WHERE status IN ('equal', 'changed') AND signature IN "
                                      "('gdd_cum_annual', 'gdd_cum_win', 'gdd_cum_spr', 'gdd_cum_sum', 'gdd_cum_fal')").df()
                sharedg = sharedg[sharedg["gage_id"].isin(set(classes.index[classes["class"] != "U"]))]
                bad_g = []
                for r in sharedg.itertuples(index=False):
                    exp = gsums[r.signature].get((r.gage_id, int(r.water_year)), 0.0)
                    act = (0.0 if np.isnan(r.new_value) or np.isnan(r.old_value) else r.new_value - r.old_value)
                    if abs(act - exp) > ANNUAL_TOL[0] + ANNUAL_TOL[1] * max(abs(exp), abs(act)):
                        bad_g.append((r.gage_id, r.signature, int(r.water_year), act, exp))
                n_expected = int(con.execute("SELECT count(*) FROM d WHERE status IN ('equal', 'changed') AND signature LIKE 'gdd_cum_%'").fetchone()[0])
                gate("5b: GDD totals = Σ thresholded daily degree-day deltas on every shared eligible year", "PASS" if (not bad_g and len(sharedg) > 0 and len(tv) == len(tdates)) else "FAIL",
                     f"{len(sharedg):,} cells compared ({n_expected:,} shared GDD cells in the window) on {len(tdates)} changed temperature dates ({len(tv)} found in the daily series); off {len(bad_g)} {bad_g[:3]}")
            else:
                n_ne = int(con.execute("SELECT count(*) FROM d WHERE status = 'changed' AND signature LIKE 'gdd_cum_%'").fetchone()[0])
                gate("5b: GDD totals = Σ thresholded daily degree-day deltas on every shared eligible year", "PASS" if n_ne == 0 else "FAIL",
                     f"no changed temperature date in the ledger, so every shared GDD cell must be unchanged: {n_ne} changed")
        else:
            gate("5b: GDD totals = Σ thresholded daily degree-day deltas on every shared eligible year", "NOT RUN", "no --daily-a/--daily-b")
    else:
        gate("5b: precipitation totals = summed daily PPT deltas", "NOT RUN", "no --delta-ledger")
        gate("5b: GDD totals = Σ thresholded daily degree-day deltas", "NOT RUN", "no --delta-ledger")
    # ---------------- 5c / 5g: the E7 gates in a second language ----------------
    # the E7 verdicts are bound to the files assessed here (sha256 of the reference side must be the
    # new run's file) and their waivers must be accepted by name on the command line
    import hashlib
    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 24), b""):
                h.update(blk)
        return h.hexdigest()
    accepted = set(a.accept_waiver or [])
    for nm, path, bound_to, waiver_keys in (
            ("5c: absolute recomputation of every changed or gained annual value (E7 annual gate, second language)", a.e7_annual_json, a.new_annual, ("allow_diff_signature", "waive_prefix")),
            ("5g: summary statistics, scalars and flags validated numerically (E7 summary gate)", a.e7_summary_json, a.new_csv, ("waive_column", "allow_gate_state_column"))):
        if not (path and os.path.isfile(path)):
            gate(nm, "NOT RUN", "no E7 verdict JSON given"); continue
        j = json.load(open(path))
        ref_sha = j.get("ref_sha256")
        bound = ref_sha == sha(bound_to) if ref_sha else False
        waivers = [w for k in waiver_keys for w in (j.get(k) or [])]
        unaccepted = [w for w in waivers if w not in accepted]
        affected = set(classes.index[classes["class"] != "U"])
        cov_file = j.get("gages_compared_file")
        covered = set(l.strip() for l in open(cov_file) if l.strip()) if cov_file and os.path.isfile(cov_file) else None
        coverage_ok = covered is not None and affected <= covered and j.get("n_gages_compared", 0) == len(covered)
        if not coverage_ok:
            unaccepted.append(f"coverage: {('no compared-gage list' if covered is None else f'{len(affected - covered)} affected gages not compared')}")
        if j.get("allow_value_diff", 0) and "allow_value_diff" not in accepted:
            unaccepted.append(f"allow_value_diff={j.get('allow_value_diff')}")
        status = "PASS" if (j.get("verdict") == "PASS" and bound and not unaccepted) else "FAIL"
        gate(nm, status, f"{path}: verdict {j.get('verdict')}; bound to the assessed file {bound}" + ("" if ref_sha else " (no ref_sha256 in the JSON)")
             + f"; compared gages {j.get('n_gages_compared')} cover the {len(affected)} affected gages: {coverage_ok}"
             + f"; waivers {waivers} (unaccepted {unaccepted}; accept with --accept-waiver); " + ", ".join(f"{k}={v}" for k, v in j.items() if k in ("failures", "gated_mismatches", "gate_state_differences", "waived", "tie_sensitive_reported", "annual_residual_reported")))
    # ---------------- 5e: expected-identical gages ----------------
    if a.expect_identical:
        listed = set(l.strip() for l in open(a.expect_identical) if l.strip())
        ids = listed & set(common); absent = sorted(listed - set(common))
        bad_e = sorted(set(diff["gage_id"]) & ids) if len(diff) else []
        bad_s = sorted(g for g in ids if gage_changes.get(g, 0) > 0)
        gate("5e: expected-identical gages identical at both layers", "PASS" if not bad_e and not bad_s and not absent and ids else "FAIL",
             f"{len(listed)} gages expected identical ({len(absent)} absent from both runs {absent[:5]}); annual differences on {len(bad_e)} {bad_e[:5]}, summary differences on {len(bad_s)} {bad_s[:5]}")
    else:
        gate("5e: expected-identical gages identical at both layers", "NOT RUN", "no --expect-identical")
    # ---------------- exposure: interpolated Q days in affected years ----------------
    if state_b is not None:
        ex = []
        for g in classes.index[classes["class"] != "U"]:
            yrs = restored_years.get(g, set()) | inserted_years.get(g, set())
            for y in sorted(yrs):
                r = state_b[(state_b["gage_id"] == g) & (state_b["wy"] == y)]
                if len(r):
                    ex.append({"gage_id": g, "water_year": y, "interpolated_q_days": len(r.iloc[0]["q_interp_dowys"].split()), "raw_na_q": r.iloc[0]["q_raw_na"], "q_valid": r.iloc[0]["q_valid"]})
        pd.DataFrame(ex).to_csv(os.path.join(a.out_dir, f"interpolated_q_exposure_{tag}.csv"), index=False)

    failed = [k for k, v in gates.items() if v == "FAIL"]; not_run = [k for k, v in gates.items() if v == "NOT RUN"]
    verdict = "FAIL" if failed else ("INCOMPLETE" if not_run else "PASS")
    lines += ["", "## Gates", "", "| gate | status |", "|---|---|"] + [f"| {k} | {v} |" for k, v in gates.items()] + [""]
    lines += [f"**Verdict: {verdict}**" + (f" — failed: {', '.join(failed)}" if failed else "") + (f" — not run: {', '.join(not_run)}" if not_run else ""), ""]
    if notes:
        lines += ["Notes: " + "; ".join(notes), ""]
    open(os.path.join(a.out_dir, f"change_census_{tag}.md"), "w").write("\n".join(lines) + "\n")
    json.dump({"gates": gates, "verdict": verdict, "classes": ccount, "annual": tot, "explanations": expl_counts, "label": a.label, "window": a.window},
              open(os.path.join(a.out_dir, f"gates_{tag}.json"), "w"), indent=1, default=str)
    log(f"change_census_{tag}.md written — verdict {verdict}")
    sys.exit(0 if verdict == "PASS" else 2 if verdict == "INCOMPLETE" else 1)


if __name__ == "__main__":
    main()
