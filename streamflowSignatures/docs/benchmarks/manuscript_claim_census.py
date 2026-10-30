#!/usr/bin/env python
"""Manuscript claim census (2026-10-30).

Recomputes, from the HydroShare files AS DEPOSITED (the staged upload folders that were
MD5-matched to the live resources on 2026-10-29/30), every count the HISSS manuscript and the
signature guidelines doc state — gages, rows, columns, years, variables, flags, overlaps.
Nothing is copied from a readme. Output: a CSV of (id, measure, value, source file) that the
claim ledger `docs/plans/2026-10-30-manuscript-accuracy-review-plan.md` §1 cites, plus a
printed summary. Plan: that file; results: `docs/plans/2026-10-30-manuscript-claim-ledger.md`.

Usage (repo .venv: pandas + pyarrow; sqlite3 for the GeoPackages):
  .venv/bin/python docs/benchmarks/manuscript_claim_census.py \
      --r1 ~/Downloads/hisss_R1_upload_2026-10-29 --r2 ~/Downloads/hisss_R2_upload_2026-10-29 \
      --r3 ~/Downloads/hisss_R3_upload_2026-10-29 \
      --r4-aug ~/Downloads/Signatures_v3_28oct2026/resource4_geometry_attributes \
      --r4-new ~/Downloads/hisss_R4_upload_final_2026-10-30 \
      --r5-eo ~/Downloads/Signatures_v3_28oct2026/resource5_eo_products \
      --lai ~/Downloads/hisss_modis_lai_monthly.parquet --excluded ~/Downloads/hisss_nlcd_excluded_gages.csv \
      --out docs/plans/2026-10-30-manuscript-claim-ledger.csv
Heavy scans (the 111.6 M-row streamflow table, the 134.6 M-row Daymet table) can be skipped
with --skip-heavy for a quick re-run of everything else.
"""
import argparse, csv, os, sqlite3, sys, time
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval", "_mean", "_median",
            "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean", "_pettitt_delta_mean",
            "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
WINDOW = {"wy1993_2025": (1993, 2025, 33), "wy1980_2025": (1980, 2025, 46)}
REC = []
T0 = time.time()


def key(s):
    s = str(s)
    return s.lstrip("0") if s.isdigit() else s.upper()


def rec(id_, measure, value, source):
    REC.append({"id": id_, "measure": measure, "value": value, "source": source})
    print(f"  [{id_}] {measure}: {value}", flush=True)


def section(name):
    print(f"\n=== {name}  (t={time.time()-T0:.0f}s)", flush=True)


# ---------------------------------------------------------------- signatures products (R1 / R2)
def census_product(folder, win):
    y0, y1, nyears = WINDOW[win]
    csvp = os.path.join(folder, f"hisss_signatures_{win}.csv")
    section(f"Product {win}: {csvp}")
    df = pd.read_csv(csvp, dtype={"gage_id": str}, low_memory=False)
    src = os.path.basename(csvp)
    cols = list(df.columns)
    rec(f"{win}.gages", "rows (gages)", len(df), src)
    rec(f"{win}.columns", "columns", len(cols), src)
    meta = cols[:20]
    flags = [c for c in cols if c.startswith("flagged_for_")]
    statcols = [c for c in cols if any(c.endswith(s) for s in SUFFIXES)]
    bases = sorted({c[: -len("_mean")] for c in cols if c.endswith("_mean") and not c.endswith("_pettitt_pre_mean") and not c.endswith("_pettitt_post_mean") and not c.endswith("_pettitt_delta_mean")})
    scalars = [c for c in cols if c not in meta and c not in flags and c not in statcols]
    rec(f"{win}.meta_cols", "metadata columns (first 20)", f"{len(meta)}: {meta}", src)
    rec(f"{win}.stat_cols", "signature × statistic columns", len(statcols), src)
    rec(f"{win}.bases", "annual signatures (bases with _mean)", len(bases), src)
    rec(f"{win}.bases_x16", "bases × 16 == stat columns", f"{len(bases)*16} == {len(statcols)}: {len(bases)*16 == len(statcols)}", src)
    rec(f"{win}.scalars", "per-gage scalar columns", f"{len(scalars)}: {scalars}", src)
    rec(f"{win}.flags", "QA flag columns", f"{len(flags)}: {flags}", src)
    rec(f"{win}.composition", "20 + bases×16 + scalars + flags", f"{len(meta)} + {len(statcols)} + {len(scalars)} + {len(flags)} = {len(meta)+len(statcols)+len(scalars)+len(flags)}", src)
    gt = df["gage_type"].value_counts(dropna=False).to_dict()
    rec(f"{win}.gage_type", "gage_type counts", gt, src)
    us = df["gage_id"].str.isdigit().sum()
    rec(f"{win}.us_ca_by_id", "US (all-digit id) / Canada", f"{us} / {len(df)-us}", src)
    an = df["area_normalized"].astype(str).str.upper()
    rec(f"{win}.area_normalized_false", "area_normalized FALSE", int((an == "FALSE").sum()), src)
    rec(f"{win}.area_normalized_false_by_type", "area_normalized FALSE by gage_type", df.loc[an == "FALSE", "gage_type"].value_counts().to_dict(), src)
    for f in flags:
        v = df[f].astype(str).str.upper()
        rec(f"{win}.{f}", f"{f} TRUE", int((v == "TRUE").sum()), src)
    qr = df.loc[an == "FALSE", "flagged_for_qann_range"].astype(str).str.upper()
    rec(f"{win}.qann_flag_on_unnormalized", "flagged_for_qann_range TRUE among un-normalized", int((qr == "TRUE").sum()), src)
    nwy = df["num_water_years"]
    rec(f"{win}.num_water_years", "num_water_years min / median / max", f"{nwy.min()} / {nwy.median()} / {nwy.max()}", src)
    rec(f"{win}.min_fraction", "min num_water_years / window years", f"{nwy.min()/nyears:.3f} (window {nyears})", src)
    rec(f"{win}.complete_records", f"gages with num_water_years == {nyears}", int((nwy == nyears).sum()), src)
    rec(f"{win}.ge30", "gages with num_water_years >= 30", int((nwy >= 30).sum()), src)
    rec(f"{win}.start_end", "start_water_year min / end_water_year max", f"{df['start_water_year'].min()} / {df['end_water_year'].max()}", src)
    ba = df["basin_area"].dropna()
    rec(f"{win}.basin_area_quartiles", "basin_area q25 / median / q75 (km2), n with area", f"{ba.quantile(.25):.0f} / {ba.median():.0f} / {ba.quantile(.75):.0f}, n={len(ba)}", src)
    # trend gate witness
    n_mean = {b: int(df[f"{b}_mean"].notna().sum()) for b in bases}
    n_senn = {b: int(df[f"{b}_senn_slp"].notna().sum()) for b in bases}
    n_cp = {b: int(df[f"{b}_pettitt_cp_year"].notna().sum()) for b in bases}
    gated = [b for b in bases if n_senn[b] < n_mean[b]]
    equal = [b for b in bases if n_senn[b] == n_mean[b]]
    rec(f"{win}.trend_gate_bases", "bases where trend count < mean count (gated) / equal (exempt or saturated)", f"{len(gated)} / {len(equal)}; equal: {equal}", src)
    rec(f"{win}.trend_gate_Qann", "Qann: mean / senn / pettitt_cp_year non-NA", f"{n_mean['Qann']} / {n_senn['Qann']} / {n_cp['Qann']}", src)
    rec(f"{win}.mean_count_max", "max non-NA _mean count over bases", max(n_mean.values()), src)
    # NA pattern on un-normalized gages by category
    cat = pd.read_csv(os.path.join(folder, "hisss_signature_categories.csv"))
    un = df[an == "FALSE"]
    for c in sorted(cat["category"].dropna().unique()):
        sigs = cat.loc[(cat["category"] == c) & (cat["kind"] == "annual"), "signature"].tolist()
        mcols = [f"{s}_mean" for s in sigs if f"{s}_mean" in df.columns]
        if not mcols:
            continue
        tot = un[mcols].notna().sum().sum(); allc = df[mcols].notna().sum().sum()
        rec(f"{win}.unnorm_na.{c}", f"non-NA _mean cells on un-normalized gages / all gages — category {c} ({len(mcols)} signatures)", f"{tot} of {len(un)*len(mcols)} / {allc} of {len(df)*len(mcols)}", src)
    sto = cat.loc[(cat["category"] == "Storage") & (cat["kind"] == "annual"), "signature"].tolist()
    rec(f"{win}.unnorm_na.storage_detail", "Storage signatures: non-NA _mean count on the un-normalized gages (of n un-normalized)", {s_: int(un[f"{s_}_mean"].notna().sum()) for s_ in sto if f"{s_}_mean" in df.columns} | {"n_unnormalized": len(un)}, src)
    fla = cat.loc[(cat["category"] == "Flashiness") & (cat["kind"] == "annual"), "signature"].tolist()
    rec(f"{win}.unnorm_na.flashiness_detail", "Flashiness signatures with any NA _mean on the un-normalized gages: non-NA count", {s_: int(un[f"{s_}_mean"].notna().sum()) for s_ in fla if f"{s_}_mean" in df.columns and un[f"{s_}_mean"].isna().any()}, src)
    sex = {c: int(df[c].sum()) for c in cols if c.startswith("season_excluded_years_")}
    rec(f"{win}.season_excluded_sum", "Σ season_excluded_years_* over gages (seasons failing the 80 % raw-observation rule → seasonal metrics NA)", sex | {"gages_with_any": int((df[[c for c in sex]].sum(axis=1) > 0).sum())}, src)
    negcols = [c for c in cols if "egative" in c or "constant" in c.lower()]
    rec(f"{win}.negative_constant_cols", "columns named Negative*/constant*", f"{len(negcols)}: {negcols[:4]}…", src)
    # dictionary + categories
    dic = pd.read_csv(os.path.join(folder, "hisss_data_dictionary.csv"))
    rec(f"{win}.dictionary_rows", "dictionary rows", len(dic), "hisss_data_dictionary.csv")
    rec(f"{win}.dictionary_matches_csv", "dictionary column_name == CSV columns (order)", list(dic["column_name"]) == cols, "hisss_data_dictionary.csv")
    rec(f"{win}.dictionary_categories", "dictionary category counts", dic["category"].value_counts().to_dict(), "hisss_data_dictionary.csv")
    rec(f"{win}.categories_rows", "categories file rows; kind counts", f"{len(cat)}; {cat['kind'].value_counts().to_dict()}", "hisss_signature_categories.csv")
    canon = pd.read_csv("docs/signature_categories.csv")
    key3 = lambda d: set(map(tuple, d[["signature", "kind", "category"]].astype(str).values))
    rec(f"{win}.categories_vs_canonical", "deposited categories (signature, kind, category) == docs/signature_categories.csv", f"{key3(cat) == key3(canon)}; canonical rows {len(canon)}", "hisss_signature_categories.csv")
    rec(f"{win}.categories_by_cat", "categories file: annual signatures per category", cat[cat.kind == "annual"]["category"].value_counts().to_dict(), "hisss_signature_categories.csv")
    rec(f"{win}.categories_scalars", "categories file: scalar per category", cat[cat.kind == "scalar"]["category"].value_counts(dropna=False).to_dict(), "hisss_signature_categories.csv")
    # annual parquet
    pqp = os.path.join(folder, f"hisss_signatures_annual_{win}.parquet")
    pf = pq.ParquetFile(pqp)
    rec(f"{win}.annual_rows", "annual parquet rows", pf.metadata.num_rows, os.path.basename(pqp))
    rec(f"{win}.annual_schema", "annual parquet columns", [f.name for f in pf.schema_arrow], os.path.basename(pqp))
    sigs, gages, wy_min, wy_max = set(), set(), 10**9, -1
    seas = Counter()
    for b in pf.iter_batches(columns=["gage_id", "signature", "water_year", "value"], batch_size=2_000_000):
        sigs.update(pc.unique(b.column("signature")).to_pylist())
        gages.update(pc.unique(b.column("gage_id")).to_pylist())
        wy_min = min(wy_min, pc.min(b.column("water_year")).as_py()); wy_max = max(wy_max, pc.max(b.column("water_year")).as_py())
        m = pc.is_in(b.column("signature"), value_set=pa.array(["Qann", "Qwin", "Qspr", "Qsum", "Qfal", "negative_ann"]))
        sub = b.filter(m)
        for s, n in zip(*np.unique(sub.column("signature").to_numpy(zero_copy_only=False), return_counts=True)):
            seas[s] += int(n)
        subv = sub.filter(pc.fill_null(pc.is_finite(sub.column("value")), False))
        for s, n in zip(*np.unique(subv.column("signature").to_numpy(zero_copy_only=False), return_counts=True)):
            seas[s + " (finite)"] += int(n)
    rec(f"{win}.annual_signatures", "annual parquet distinct signatures", len(sigs), os.path.basename(pqp))
    rec(f"{win}.annual_sigs_eq_bases", "annual parquet signatures == CSV bases", sorted(sigs) == bases, os.path.basename(pqp))
    rec(f"{win}.annual_gages", "annual parquet distinct gages == CSV gages", f"{len(gages)}; equal sets: {gages == set(df['gage_id'])}", os.path.basename(pqp))
    rec(f"{win}.annual_wy_range", "annual parquet water_year min / max", f"{wy_min} / {wy_max}", os.path.basename(pqp))
    keys = pq.read_table(pqp, columns=["gage_id", "signature", "water_year"], read_dictionary=["gage_id", "signature"]).to_pandas()
    ndup = int(keys.duplicated().sum())
    rec(f"{win}.annual_keys_unique", "annual parquet (gage_id, signature, water_year) keys unique (duplicate rows)", f"{ndup} duplicates of {len(keys)} rows: {ndup == 0}", os.path.basename(pqp))
    del keys
    rec(f"{win}.annual_seasonal_rows", "rows and finite rows per Qann / Qwin / Qspr / Qsum / Qfal / negative_ann (seasonal 80 % masking witness: seasonal finite < Qann finite)", dict(sorted(seas.items())), os.path.basename(pqp))
    return set(df["gage_id"]), df


# ---------------------------------------------------------------- inputs (R3)
def census_inputs(r3, p1, p2, skip_heavy):
    section(f"Inputs R3: {r3}")
    met = pd.read_csv(os.path.join(r3, "hisss_gage_metadata.csv"), dtype={"gage_id": str}, low_memory=False)
    src = "hisss_gage_metadata.csv"
    rec("r3.meta_rows", "metadata rows (candidate gages)", len(met), src)
    rec("r3.meta_cols", "metadata columns", list(met.columns), src)
    rec("r3.meta_by_type", "candidates by gage_type", met["gage_type"].value_counts(dropna=False).to_dict(), src)
    ps = met.groupby(["processing_status", "gage_type"]).size().unstack(fill_value=0)
    rec("r3.meta_status", "processing_status × gage_type", ps.to_dict(), src)
    succ = met[met["processing_status"] == "success"]
    rec("r3.compiled", "compiled (success) gages; by type", f"{len(succ)}; {succ['gage_type'].value_counts().to_dict()}", src)
    an = succ["area_normalized"].astype(str).str.upper()
    rec("r3.compiled_unnormalized", "compiled with area_normalized FALSE; by type", f"{int((an=='FALSE').sum())}; {succ.loc[an=='FALSE','gage_type'].value_counts().to_dict()}", src)
    rec("r3.all_unnormalized", "all candidates with area_normalized FALSE", int((met['area_normalized'].astype(str).str.upper() == 'FALSE').sum()), src)
    ba = succ["basin_area"]
    rec("r3.compiled_gt85k", "compiled with published basin_area > 85,000 km2 / > 100,000 km2 / no area", f"{int((ba > 85000).sum())} / {int((ba > 100000).sum())} / {int(ba.isna().sum())}", src)
    rec("r3.meta_years", "compiled: num_years min / max; start_year min; end_year max", f"{succ['num_years'].min()} / {succ['num_years'].max()}; {succ['start_year'].min()}; {succ['end_year'].max()}", src)
    rec("r3.meta_padding", "metadata USGS ids that start with 0 (should be 0 = zero-stripped)", int(met.loc[met.gage_id.str.isdigit(), "gage_id"].str.startswith("0").sum()), src)
    compiled_keys = set(succ["gage_id"].map(key))
    # ledger, flags, lists, would-qualify, interference, dictionary
    led = pd.read_csv(os.path.join(r3, "hisss_streamflow_qualifier_corrections.csv"), dtype={"gage_id": str}, low_memory=False)
    rec("r3.ledger_rows", "qualifier ledger rows / distinct gages", f"{len(led)} / {led['gage_id'].nunique()}", "hisss_streamflow_qualifier_corrections.csv")
    rec("r3.ledger_old_flag", "ledger old_flag counts", led["old_flag"].value_counts().to_dict(), "hisss_streamflow_qualifier_corrections.csv")
    rec("r3.ledger_status", "ledger status counts; status_today counts", f"{led['status'].value_counts().to_dict()}; {led['status_today'].value_counts(dropna=False).to_dict()}", "hisss_streamflow_qualifier_corrections.csv")
    rec("r3.ledger_dates", "ledger date min / max", f"{led['date'].min()} / {led['date'].max()}", "hisss_streamflow_qualifier_corrections.csv")
    fl = pd.read_csv(os.path.join(r3, "hisss_daymet_basin_flags.csv"), dtype={"site_id": str})
    rec("r3.flags_rows", "Daymet flag rows; low_confidence / low_pixel_support; in_p1 / in_p2", f"{len(fl)}; {int(fl['low_confidence'].sum())} / {int(fl['low_pixel_support'].sum())}; {int(fl['in_p1'].sum())} / {int(fl['in_p2'].sum())}", "hisss_daymet_basin_flags.csv")
    rec("r3.flags_reasons", "Daymet flag component counts", {c: int(fl[c].sum()) for c in ["hydrobasins_fallback", "area_mismatch", "no_area_reference", "low_pixel_support"]}, "hisss_daymet_basin_flags.csv")
    for win, pset in (("wy1993_2025", p1), ("wy1980_2025", p2)):
        lst = [l.strip() for l in open(os.path.join(r3, f"hisss_gage_selection_{win}.txt")) if l.strip()]
        rec(f"r3.list_{win}", f"gage list {win}: ids; == product gage set", f"{len(lst)}; {set(lst) == pset}", f"hisss_gage_selection_{win}.txt")
    wq = pd.read_csv(os.path.join(r3, "hisss_would_qualify_not_added.csv"), dtype={"gage_id": str})
    rec("r3.would_qualify", "would-qualify rows per window", wq["window"].value_counts().to_dict(), "hisss_would_qualify_not_added.csv")
    hi = pd.read_csv(os.path.join(r3, "hisss_hydat_interference.csv"), dtype=str)
    rec("r3.interference", "HYDAT interference rows; RHBN counts; REGULATED counts", f"{len(hi)}; {hi['RHBN'].value_counts(dropna=False).to_dict()}; {hi['REGULATED'].value_counts(dropna=False).to_dict()}", "hisss_hydat_interference.csv")
    idd = pd.read_csv(os.path.join(r3, "hisss_input_data_dictionary.csv"))
    rec("r3.input_dictionary", "input dictionary rows per table", idd["table"].value_counts().to_dict(), "hisss_input_data_dictionary.csv")
    # Daymet schema (cheap)
    dpf = pq.ParquetFile(os.path.join(r3, "hisss_daymet_basin_daily.parquet"))
    rec("r3.daymet_rows", "Daymet rows (file metadata)", dpf.metadata.num_rows, "hisss_daymet_basin_daily.parquet")
    rec("r3.daymet_schema", "Daymet columns", [f.name for f in dpf.schema_arrow], "hisss_daymet_basin_daily.parquet")
    spf = pq.ParquetFile(os.path.join(r3, "hisss_streamflow_daily.parquet"))
    rec("r3.q_rows", "streamflow rows (file metadata)", spf.metadata.num_rows, "hisss_streamflow_daily.parquet")
    rec("r3.q_schema", "streamflow columns", [f.name for f in spf.schema_arrow], "hisss_streamflow_daily.parquet")
    daymet_sites = None
    if skip_heavy:
        print("  (heavy scans skipped)")
        return met, compiled_keys, daymet_sites
    # Daymet scan
    section("Daymet scan (site_id, Date, prcp)")
    parts = []
    DV = ["prcp", "tmin", "tmax", "swe", "vp", "srad"]
    bad = {v: {"null": 0, "nan": 0, "inf": 0} for v in DV}
    for b in dpf.iter_batches(columns=["site_id", "Date"] + DV, batch_size=2_000_000):
        t = pa.Table.from_batches([b])
        for v in DV:
            c = t.column(v)
            bad[v]["null"] += pc.sum(pc.cast(pc.is_null(c), pa.int64())).as_py() or 0
            bad[v]["nan"] += pc.sum(pc.cast(pc.fill_null(pc.is_nan(c), False), pa.int64())).as_py() or 0
            bad[v]["inf"] += pc.sum(pc.cast(pc.fill_null(pc.is_inf(c), False), pa.int64())).as_py() or 0
        parts.append(t.select(["site_id", "Date"]).group_by("site_id").aggregate([("Date", "min"), ("Date", "max"), ("Date", "count"), ("Date", "count_distinct")]))
    g = pa.concat_tables(parts).group_by("site_id").aggregate([("Date_min", "min"), ("Date_max", "max"), ("Date_count", "sum"), ("Date_count_distinct", "sum")]).to_pandas()
    rec("r3.daymet_missing", "Daymet null / NaN / inf per variable (all six)", bad, "hisss_daymet_basin_daily.parquet")
    rec("r3.daymet_grid", "Daymet rows == distinct dates per site (summed over batches; equality ⇒ no duplicate dates within a batch) and rows per site == 46 × 365", f"{int((g['Date_count_sum'] == g['Date_count_distinct_sum']).sum())} of {len(g)} sites; all 16,790: {bool((g['Date_count_sum'] == 16790).all())}", "hisss_daymet_basin_daily.parquet")
    rec("r3.daymet_sites", "Daymet distinct site_id", len(g), "hisss_daymet_basin_daily.parquet")
    rec("r3.daymet_dates", "Daymet Date min / max; rows per site min / max", f"{g['Date_min_min'].min()} / {g['Date_max_max'].max()}; {g['Date_count_sum'].min()} / {g['Date_count_sum'].max()}", "hisss_daymet_basin_daily.parquet")
    dk = set(g["site_id"].map(key))
    rec("r3.daymet_vs_compiled", "compiled gages with Daymet / without; Daymet sites not compiled", f"{len(compiled_keys & dk)} / {sorted(compiled_keys - dk)}; {len(dk - compiled_keys)}", "hisss_daymet_basin_daily.parquet")
    daymet_sites = dk
    # streamflow scan
    section("Streamflow scan (gage_id, Date, Q)")
    parts = []
    for b in spf.iter_batches(columns=["gage_id", "Date", "Q"], batch_size=2_000_000):
        t = pa.Table.from_batches([b])
        q = t.column("Q")
        t = t.append_column("isnan", pc.cast(pc.is_nan(q), pa.int64())).append_column("isneg", pc.cast(pc.fill_null(pc.less(q, 0.0), False), pa.int64())).append_column("isnull", pc.cast(pc.is_null(q), pa.int64()))
        parts.append(t.group_by("gage_id").aggregate([("Date", "min"), ("Date", "max"), ("Date", "count"), ("Date", "count_distinct"), ("isnan", "sum"), ("isneg", "sum"), ("isnull", "sum")]))
    g = pa.concat_tables(parts).group_by("gage_id").aggregate([("Date_min", "min"), ("Date_max", "max"), ("Date_count", "sum"), ("Date_count_distinct", "sum"), ("isnan_sum", "sum"), ("isneg_sum", "sum"), ("isnull_sum", "sum")]).to_pandas()
    g.columns = ["gage_id", "dmin", "dmax", "n", "ndist", "nan", "neg", "null"]
    rec("r3.q_dates_unique", "streamflow gages whose row count == distinct Date count (batch-summed)", f"{int((g.n == g.ndist).sum())} of {len(g)}", "hisss_streamflow_daily.parquet")
    rec("r3.q_nonmissing", "streamflow rows with a discharge value (rows − null Q)", int(g.n.sum() - g['null'].sum()), "hisss_streamflow_daily.parquet")
    rec("r3.q_gages", "streamflow distinct gages; US (digit ids) / Canada", f"{len(g)}; {int(g.gage_id.str.isdigit().sum())} / {int((~g.gage_id.str.isdigit()).sum())}", "hisss_streamflow_daily.parquet")
    rec("r3.q_rows_sum", "streamflow rows (scan)", int(g.n.sum()), "hisss_streamflow_daily.parquet")
    rec("r3.q_dates", "streamflow Date min / max", f"{g.dmin.min()} / {g.dmax.max()}", "hisss_streamflow_daily.parquet")
    us = g[g.gage_id.str.isdigit()]; ca = g[~g.gage_id.str.isdigit()]
    rec("r3.q_dates_by_agency", "max Date US / Canada; gages ending 2025-10-01", f"{us.dmax.max()} / {ca.dmax.max()}; {int((g.dmax == pd.Timestamp('2025-10-01').date()).sum())}", "hisss_streamflow_daily.parquet")
    rec("r3.q_na_rows", "rows with NA/NaN Q (NaN / null); negative Q rows / gages", f"{int(g['nan'].sum())} / {int(g['null'].sum())}; {int(g.neg.sum())} / {int((g.neg > 0).sum())}", "hisss_streamflow_daily.parquet")
    rec("r3.q_padding", "streamflow USGS ids starting with 0", int(us.gage_id.str.startswith("0").sum()), "hisss_streamflow_daily.parquet")
    rec("r3.q_vs_meta", "streamflow gage set == compiled (success) metadata keys", set(g.gage_id.map(key)) == compiled_keys, "hisss_streamflow_daily.parquet")
    return met, compiled_keys, daymet_sites


# ---------------------------------------------------------------- geometry + attributes (R4)
def census_geometry(aug, new, p1, p2):
    section(f"Geometry R4: {aug} + {new}")
    out = {}
    for label, path in (("simplified", os.path.join(aug, "hisss_watershed_boundaries.gpkg")), ("unsimplified", os.path.join(new, "hisss_watershed_boundaries_unsimplified.gpkg"))):
        con = sqlite3.connect(path); src = os.path.basename(path)
        n = con.execute("select count(*) from watersheds").fetchone()[0]
        rec(f"r4.{label}.features", "polygons", n, src)
        rec(f"r4.{label}.sources", "watershed_geom_source counts", dict(con.execute("select watershed_geom_source, count(*) from watersheds group by 1").fetchall()), src)
        rec(f"r4.{label}.low_conf", "low_confidence TRUE; geom_simplified TRUE; area_flag TRUE", tuple(con.execute("select sum(low_confidence), sum(geom_simplified), sum(area_flag) from watersheds").fetchone()), src)
        rec(f"r4.{label}.areas", "geom_area_km2 > 85,000 / > 100,000; max; polygons 85,000–100,000", tuple(con.execute("select sum(geom_area_km2>85000), sum(geom_area_km2>100000), round(max(geom_area_km2)), sum(geom_area_km2>85000 and geom_area_km2<=100000) from watersheds").fetchone()), src)
        rec(f"r4.{label}.by_type", "gage_type counts", dict(con.execute("select gage_type, count(*) from watersheds group by 1").fetchall()), src)
        ids = [r[0] for r in con.execute("select gage_id from watersheds")]
        rec(f"r4.{label}.unique_ids", "distinct gage_id", len(set(ids)), src)
        ks = set(map(key, ids))
        rec(f"r4.{label}.product_overlap", "product gages with a polygon (zero-stripped join): #1 / #2", f"{len({key(x) for x in p1} & ks)} / {len({key(x) for x in p2} & ks)}", src)
        out[label] = ks
        con.close()
    rec("r4.layers_diff", "unsimplified minus simplified ids; simplified minus unsimplified", f"{len(out['unsimplified'] - out['simplified'])}; {sorted(out['simplified'] - out['unsimplified'])}", "both GeoPackages")
    gp = pq.ParquetFile(os.path.join(aug, "hisss_watershed_boundaries.parquet"))
    rec("r4.geoparquet", "GeoParquet rows / columns", f"{gp.metadata.num_rows} / {len(gp.schema_arrow.names)}", "hisss_watershed_boundaries.parquet")
    qa = pd.read_csv(os.path.join(aug, "hisss_watershed_boundaries_qa.csv"), dtype={"gage_id": str})
    rec("r4.qa_rows", "QA csv rows; low_confidence TRUE", f"{len(qa)}; {int(qa['low_confidence'].astype(str).str.upper().eq('TRUE').sum())}", "hisss_watershed_boundaries_qa.csv")
    at = pd.read_csv(os.path.join(aug, "hisss_watershed_hydroatlas_attributes.csv"), dtype={"gage_id": str}, low_memory=False)
    src = "hisss_watershed_hydroatlas_attributes.csv"
    rec("r4.attr_rows_cols", "attribute table rows × columns", f"{len(at)} × {len(at.columns)}", src)
    rec("r4.attr_has_downstream", "Downstream_HB_ID column present", "Downstream_HB_ID" in at.columns, src)
    rec("r4.attr_first_cols", "first 13 columns", list(at.columns[:13]), src)
    rec("r4.attr_product_overlap", "product gages with an attribute row: #1 / #2 (zero-stripped)", f"{len({key(x) for x in p1} & set(at.gage_id.map(key)))} / {len({key(x) for x in p2} & set(at.gage_id.map(key)))}", src)
    rec("r4.attr_overlap_layers", "attribute rows with a polygon: simplified / unsimplified", f"{len(set(at.gage_id.map(key)) & out['simplified'])} / {len(set(at.gage_id.map(key)) & out['unsimplified'])}", src)
    rec("r4.attr_padding", "attribute USGS ids starting with 0", int(at.loc[at.gage_id.str.isdigit(), 'gage_id'].str.startswith('0').sum()), src)
    if "n_member_basins" in at.columns and "watershed_area_sum_SUB_AREA" in at.columns:
        nm, sa = at["n_member_basins"].sum(), at["watershed_area_sum_SUB_AREA"].sum()
        rec("r4.level12_mean_area", "mean member level-12 subbasin area over the domain (Σ SUB_AREA / Σ n_member_basins), km2; median per-gage mean", f"{sa/nm:.1f}; {(at['watershed_area_sum_SUB_AREA']/at['n_member_basins']).median():.1f}", src)
    dd = pd.read_csv(os.path.join(aug, "hisss_hydroatlas_attributes_dictionary.csv"))
    rec("r4.attr_dictionary", "attribute dictionary rows; aggregation counts", f"{len(dd)}; {dd['aggregation'].value_counts().to_dict()}", "hisss_hydroatlas_attributes_dictionary.csv")
    return out


# ---------------------------------------------------------------- vegetation / land cover (R5)
def census_eo(eo, lai, excluded, p1, p2, daymet):
    section(f"EO R5: {eo}, {lai}")
    k1, k2 = {key(x) for x in p1}, {key(x) for x in p2}
    lp = pq.ParquetFile(lai)
    t = lp.read(columns=["gage_id", "canon_id", "year", "month"]).to_pandas()
    rec("r5.lai_rows", "LAI rows", len(t), os.path.basename(lai))
    rec("r5.lai_gages", "LAI distinct gage_id / canon_id", f"{t.gage_id.nunique()} / {t.canon_id.nunique()}", os.path.basename(lai))
    ym = t[["year", "month"]].drop_duplicates().sort_values(["year", "month"])
    rec("r5.lai_months", "LAI distinct (year, month); first / last", f"{len(ym)}; {ym.iloc[0].tolist()} / {ym.iloc[-1].tolist()}", os.path.basename(lai))
    rec("r5.lai_cols", "LAI columns", lp.schema_arrow.names, os.path.basename(lai))
    lk = set(t.canon_id.unique())
    rec("r5.lai_product_overlap", "product gages with LAI rows (canon_id join): #1 / #2; exact gage_id join: #1 / #2", f"{len(k1 & lk)} / {len(k2 & lk)}; {len(p1 & set(t.gage_id))} / {len(p2 & set(t.gage_id))}", os.path.basename(lai))
    lu = pq.read_table(os.path.join(eo, "hisss_modis_lulc_annual.parquet")).to_pandas()
    pct = [c for c in lu.columns if c.endswith("_pct")]
    rec("r5.lulc_rows", "LULC rows; gages; years min / max / n", f"{len(lu)}; {lu.gage_id.nunique()}; {lu.year.min()} / {lu.year.max()} / {lu.year.nunique()}", "hisss_modis_lulc_annual.parquet")
    rec("r5.lulc_schemes", "LULC class columns; by scheme", f"{len(pct)}; {Counter(c.split('_c')[0] for c in pct)}", "hisss_modis_lulc_annual.parquet")
    rec("r5.lulc_product_overlap", "product gages with LULC rows (canon_id): #1 / #2", f"{len(k1 & set(lu.canon_id))} / {len(k2 & set(lu.canon_id))}", "hisss_modis_lulc_annual.parquet")
    if daymet is not None:
        rec("r5.sig_daymet_modis", "product gages with Daymet AND MODIS LULC: #1 / #2; either", f"{len(k1 & set(lu.canon_id) & daymet)} / {len(k2 & set(lu.canon_id) & daymet)}; {len((k1 | k2) & set(lu.canon_id) & daymet)}", "R3 Daymet + R5 LULC")
    nl = pq.read_table(os.path.join(eo, "hisss_nlcd_annual.parquet")).to_pandas()
    ncls = [c for c in nl.columns if c.startswith("nlcd_c") and c.endswith("_pct")]
    rec("r5.nlcd_rows", "NLCD rows; gages; years min / max / n; class columns; collection", f"{len(nl)}; {nl.gage_id.nunique()}; {nl.year.min()} / {nl.year.max()} / {nl.year.nunique()}; {len(ncls)}; {nl.nlcd_collection.unique().tolist()}", "hisss_nlcd_annual.parquet")
    rec("r5.nlcd_product_overlap", "product gages with NLCD rows (canon_id): #1 / #2", f"{len(k1 & set(nl.canon_id))} / {len(k2 & set(nl.canon_id))}", "hisss_nlcd_annual.parquet")
    ex = pd.read_csv(excluded, dtype=str)
    rec("r5.excluded", "NLCD excluded gages; reasons", f"{len(ex)}; {ex['exclusion_reason'].value_counts().to_dict()}", os.path.basename(excluded))
    rec("r5.nlcd_plus_excluded", "NLCD gages + excluded == US gages in LULC universe", f"{nl.gage_id.nunique() + len(ex)} vs {int(lu.loc[lu.gage_id.str.isdigit(), 'gage_id'].nunique())}", "hisss_nlcd_annual.parquet")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("r1", "r2", "r3", "r4-aug", "r4-new", "r5-eo", "lai", "excluded", "out"):
        ap.add_argument(f"--{a}", required=True)
    ap.add_argument("--skip-heavy", action="store_true")
    a = ap.parse_args()
    ex = lambda p: os.path.expanduser(p)
    p1, _ = census_product(ex(a.r1), "wy1993_2025")
    p2, _ = census_product(ex(a.r2), "wy1980_2025")
    rec("products.overlap", "gages in both products / only #1 / only #2", f"{len(p1 & p2)} / {len(p1 - p2)} / {len(p2 - p1)}", "both CSVs")
    met, compiled, daymet = census_inputs(ex(a.r3), p1, p2, a.skip_heavy)
    rec("products.in_compiled", "product gages ⊂ compiled: #1 / #2", f"{p1 <= set(met.loc[met.processing_status=='success','gage_id'].map(key)) or len({key(x) for x in p1} - compiled)==0} / {len({key(x) for x in p2} - compiled)==0}", "metadata")
    census_geometry(ex(a.r4_aug), ex(a.r4_new), p1, p2)
    census_eo(ex(a.r5_eo), ex(a.lai), ex(a.excluded), p1, p2, daymet)
    with open(ex(a.out), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "measure", "value", "source"]); w.writeheader()
        for r in REC:
            w.writerow({k: str(v) for k, v in r.items()})
    print(f"\nwrote {len(REC)} measurements to {a.out}  (t={time.time()-T0:.0f}s)")


if __name__ == "__main__":
    main()
