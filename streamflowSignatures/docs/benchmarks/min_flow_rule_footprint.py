"""Footprint of the ingestion's per-year flow rule (>= 30 days with Q > 0.0001 mm/day per
calendar year; R/helperFunctions.R::process_gages_rawToRaw) among the 8,014 ADMITTED gages,
plus the USGS qualifier-code mask of the February 2026 compile (generate_streamflow_dt kept only
A, A e, P, P e until 2026-10-26; the rule is R/usgs_qualifiers.R since).

Per gage and calendar year: non-NA days, days > 0, days > 0.0001 (Q as stored: mm/day, raw
units for the 73 un-normalized gages). Tabulates the years with >= 30 non-NA days that the rule
would not count, the gages whose qualifying-year count depends on the threshold, and the days
whose Q was masked to NA by qualifier code. Admission is unaffected by construction (every
admitted gage has >= 20 counted years); the ingestion keeps every day of an admitted gage.

Usage: python docs/benchmarks/min_flow_rule_footprint.py <combined_streamflow_data.parquet> <out_dir>
Census addendum: docs/plans/2026-09-04-filtering-stage-census.md §10 (2026-10-21).
"""
import sys
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
with open(src, "rb") as f:
    f.seek(-4, 2); footer = f.read(4)
pf = pq.ParquetFile(src)
print("size", src.stat().st_size, "footer", footer, "rows", pf.metadata.num_rows, flush=True)
acc, flags = [], []
for i, b in enumerate(pf.iter_batches(batch_size=4_000_000, columns=["gage_id", "year", "Q", "flag"])):
    d = b.to_pandas()
    d["nonna"] = d.Q.notna(); d["pos"] = d.Q > 0; d["thr"] = d.Q > 0.0001
    g = d.groupby(["gage_id", "year"], sort=False)[["nonna", "pos", "thr"]].sum()
    g["rows"] = d.groupby(["gage_id", "year"], sort=False).size()
    acc.append(g.reset_index())
    fl = d.groupby([d.flag.fillna("<NA>"), d.Q.isna().rename("q_is_na")], sort=False).size().rename("n").reset_index()
    fl["n_gages"] = d.groupby([d.flag.fillna("<NA>"), d.Q.isna().rename("q_is_na")], sort=False).gage_id.nunique().values
    flags.append(fl)
    if i % 10 == 0: print("batch", i, flush=True)
a = pd.concat(acc).groupby(["gage_id", "year"], as_index=False).sum()
a.to_parquet(out / "admitted_year_counts.parquet", index=False)
a = a[(a.year >= 1979) & (a.year <= 2025)]
ok_nonna, ok_thr, ok_pos = a.nonna >= 30, a.thr >= 30, a.pos >= 30
print("gage-years (CY 1979-2025):", len(a))
print("  >=30 non-NA days:", ok_nonna.sum(), "| >=30 days > 0.0001:", ok_thr.sum())
print("  >=30 non-NA but <30 days > 0.0001 (rule-dropped years):", (ok_nonna & ~ok_thr).sum(),
      "in", a.loc[ok_nonna & ~ok_thr, "gage_id"].nunique(), "gages")
print("  >=30 non-NA but <30 days > 0 (zero-flow years):", (ok_nonna & ~ok_pos).sum(),
      "in", a.loc[ok_nonna & ~ok_pos, "gage_id"].nunique(), "gages")
print("  >=30 days > 0 but <30 days > 0.0001 (threshold-level only):", (ok_pos & ~ok_thr).sum(),
      "in", a.loc[ok_pos & ~ok_thr, "gage_id"].nunique(), "gages")
per_gage = a.groupby("gage_id").agg(yrs_nonna=("nonna", lambda s: (s >= 30).sum()),
                                    yrs_pos=("pos", lambda s: (s >= 30).sum()),
                                    yrs_thr=("thr", lambda s: (s >= 30).sum()))
per_gage.to_csv(out / "admitted_gage_years.csv")
print("gages:", len(per_gage), "| yrs_thr < 20:", (per_gage.yrs_thr < 20).sum(),
      "| yrs_thr < yrs_nonna:", (per_gage.yrs_thr < per_gage.yrs_nonna).sum(),
      "| yrs_pos < yrs_nonna:", (per_gage.yrs_pos < per_gage.yrs_nonna).sum())
print("qualifying-year gap (yrs_nonna - yrs_thr):")
print((per_gage.yrs_nonna - per_gage.yrs_thr).value_counts().sort_index().to_string())
f = pd.concat(flags).groupby(["flag", "q_is_na"], as_index=False).agg(n=("n", "sum"), n_gages=("n_gages", "max"))
f = f.sort_values("n", ascending=False)
f.to_csv(out / "flag_mask_tally.csv", index=False)
print("\nqualifier codes (n days, max gages per batch), Q masked to NA where q_is_na:")
print(f.head(40).to_string(index=False))

# ---- pass 2: the qualifier-code mask, per gage and water year, and the water years the
# preprocessor's two rules (> 30 NA days, NA run > 3 days) reject ONLY because of the `A R` mask
import pyarrow.parquet as pq2
parts = []
for b in pf.iter_batches(batch_size=4_000_000, columns=["gage_id", "water_year", "Q", "flag"]):
    d = b.to_pandas(); m = d[d.Q.isna() & d.flag.notna()]
    parts.append(m.groupby(["gage_id", "water_year", "flag"]).size().rename("n").reset_index())
mk = pd.concat(parts).groupby(["gage_id", "water_year", "flag"], as_index=False).n.sum()
mk.to_csv(out / "masked_days_by_gage_wy_flag.csv", index=False)
print("\nmasked days by qualifier code (distinct gages, gage-water-years):")
print(mk.groupby("flag").agg(days=("n", "sum"), gages=("gage_id", "nunique"), gage_wys=("n", "size"))
        .sort_values("days", ascending=False).to_string())
ar_ids = sorted(set(mk.loc[mk.flag == "A R", "gage_id"]))
t = pq2.read_table(src, columns=["gage_id", "Date", "water_year", "Q", "flag"],
                   filters=[("gage_id", "in", ar_ids)]).to_pandas().sort_values(["gage_id", "Date"])
t["na"] = t.Q.isna(); t["na_wo"] = t.na & ~(t.flag == "A R")
def maxrun(x):
    best = cur = 0
    for v in x:
        cur = cur + 1 if v else 0; best = max(best, cur)
    return best
rows = []
for (gid, wy), d in t.groupby(["gage_id", "water_year"]):
    if not (d.flag == "A R").any():
        continue
    rows.append((gid, wy, len(d), int(d.na.sum()), int(d.na_wo.sum()), maxrun(d.na.values), maxrun(d.na_wo.values)))
r = pd.DataFrame(rows, columns=["gage_id", "water_year", "n", "na_with", "na_without", "run_with", "run_without"])
rej_with = (r.na_with > 30) | (r.run_with > 3); rej_wo = (r.na_without > 30) | (r.run_without > 3)
r["lost_to_mask"] = rej_with & ~rej_wo
r.to_csv(out / "AR_gage_wy_budget.csv", index=False)
print(f"\n'A R' gages {len(ar_ids)} | gage-WYs with A R days {len(r)} | rejected WITH the mask {rej_with.sum()} | "
      f"WITHOUT {rej_wo.sum()} | water years lost to the mask {r.lost_to_mask.sum()} in {r[r.lost_to_mask].gage_id.nunique()} gages")
print(r[r.lost_to_mask].groupby("gage_id").water_year.agg(["count", "min", "max"]).sort_values("count", ascending=False).head(10).to_string())
