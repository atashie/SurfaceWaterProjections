"""Write the basin quality flags that travel with the Daymet climate file.

User decision 2026-10-01: the questionable polygons stay in the climate file and are
FLAGGED, not dropped. The products' signature CSV keeps its column contract, so the flags
are a companion table keyed by the same `site_id`. One row is written per basin with any
flag:
  hydrobasins_fallback   the polygon is a HydroBASINS lev-12 upstream union (no official
                         GAGES-II / ECCC outline exists for the gage)
  area_mismatch          polygon area > 50 % off the metadata drainage area
  no_area_reference      no metadata drainage area, so the area screen that caught 05KH009
                         cannot run (only flagged for HydroBASINS polygons; official outlines
                         without a metadata area are not flagged)
  low_pixel_support      the basin's area is under 4 Daymet cells (coverage sum < 4; it may
                         touch more cells partially, so `n_cells` can exceed 4)
`low_confidence` (from the polygon layer) = hydrobasins_fallback or area_mismatch. Which
of the two areas is wrong is NOT decided here: the table gives both.

Usage:
    python daymet_basin_flags.py --polygons-qa <..._qa.csv> --basins <daymet_basins.csv>
        --out <daymet_basin_flags.csv> [--product NAME=SIGNATURES.csv ...]
--product adds an in_<NAME> column and the product's `area_normalized` value (raw m³/s
gages: their Q-to-PPT signatures are NA by design).
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from daymet_common import canon  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons-qa", required=True)
    ap.add_argument("--basins", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--product", action="append", default=[])
    a = ap.parse_args()
    b = pd.read_csv(a.basins, dtype={"site_id": str, "boundary_gage_id": str, "canon_id": str})
    q = pd.read_csv(a.polygons_qa, dtype={"gage_id": str, "canon_id": str})
    q = q[["canon_id", "basin_area", "area_rel_diff", "area_flag"]]
    t = b.merge(q, on="canon_id", how="left", validate="one_to_one", indicator=True)
    miss = t.loc[t["_merge"] == "left_only", "site_id"].tolist()
    if miss:      # an unmatched basin would silently skip the area screen
        sys.exit(f"{len(miss)} basins have no polygon QA row (e.g. {miss[:5]}): check the canon_id "
                 f"spelling of the two files")
    t = t.drop(columns="_merge")
    t["hydrobasins_fallback"] = t["watershed_geom_source"] == "hydrobasins"
    t["area_mismatch"] = t["area_flag"].fillna(False).astype(bool)
    t["no_area_reference"] = t["hydrobasins_fallback"] & t["basin_area"].isna()
    t["low_pixel_support"] = t["low_pixel_support"].astype(bool)
    flags = ["hydrobasins_fallback", "area_mismatch", "no_area_reference", "low_pixel_support"]
    t = t[t[flags].any(axis=1)].copy()

    def km2(x):
        return f"{x:,.0f}" if x >= 100 else f"{x:.3g}"

    def reason(r):
        out = []
        if r["hydrobasins_fallback"]:
            out.append("HydroBASINS fallback polygon")
        if r["area_mismatch"]:
            out.append(f"polygon {km2(r['geom_area_km2'])} km2 vs metadata {km2(r['basin_area'])} km2 "
                       f"({100 * (r['geom_area_km2'] / r['basin_area'] - 1):+.0f} %)")
        if r["no_area_reference"]:
            out.append("no metadata drainage area to check against")
        if r["low_pixel_support"]:
            out.append(f"{r['coverage_sum']:.2f} Daymet cells")
        return "; ".join(out)
    t["reason"] = t.apply(reason, axis=1)
    cols = ["site_id", "canon_id", "gage_type", "watershed_geom_source", "geom_area_km2", "basin_area",
            "area_rel_diff", "n_cells", "coverage_sum", "low_confidence"] + flags
    t = t[cols + ["reason"]].rename(columns={"basin_area": "metadata_area_km2"})
    for spec in a.product:
        name, _, path = spec.partition("=")
        p = pd.read_csv(path, usecols=["gage_id", "area_normalized"], dtype={"gage_id": str})
        p["canon_id"] = p["gage_id"].map(canon)
        t[f"in_{name}"] = t["canon_id"].isin(set(p["canon_id"]))
        an = dict(zip(p["canon_id"], p["area_normalized"]))
        t[f"area_normalized_{name}"] = t["canon_id"].map(an)
    t = t.sort_values(["low_confidence", "site_id"], ascending=[False, True])
    t.to_csv(a.out, index=False)
    print(f"{len(t)} flagged basins -> {a.out}: " + ", ".join(f"{f} {int(t[f].sum())}" for f in flags)
          + f"; low_confidence {int(t['low_confidence'].sum())}"
          + "".join(f"; in {n.split('=')[0]}: {int(t['in_' + n.split('=')[0]].sum())} "
                    f"(low_confidence {int((t['in_' + n.split('=')[0]] & t['low_confidence']).sum())})"
                    for n in a.product))


if __name__ == "__main__":
    main()
