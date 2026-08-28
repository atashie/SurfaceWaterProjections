#!/usr/bin/env python
import os
import pandas as pd
import re

# Working dir for intermediate CSVs — ENV-overridable (was a session-specific scratchpad path).
SCRATCH = os.environ.get("EO_RECONCILE_SCRATCH", os.path.join(os.path.dirname(os.path.abspath(__file__)), "scratch_eo_task1"))
GOLDEN = os.environ.get(
    "EO_RECONCILE_GOLDEN",
    "/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/golden-outputs/streamflow_signatures_full_10feb2026.csv")

def load_ids(path, col="gage_id"):
    df = pd.read_csv(path, usecols=[col], dtype=str)
    s = df[col].dropna().astype(str).str.strip()
    return s

def is_numeric(x):
    return x.isdigit()

def is_canadian(x):
    # alphanumeric (contains a letter)
    return bool(re.search('[A-Za-z]', x))

def canon(x):
    """Canonical key: numeric ids -> strip leading zeros (int-form);
    alphanumeric (Canadian) -> uppercase, as-is."""
    x = str(x).strip()
    if x.isdigit():
        return str(int(x))   # strip leading zeros
    return x.upper()

def characterize(name, ids):
    ids = list(ids)
    n = len(ids)
    nu = len(set(ids))
    numeric = [x for x in ids if is_numeric(x)]
    canadian = [x for x in ids if is_canadian(x)]
    other = [x for x in ids if not is_numeric(x) and not is_canadian(x)]
    # zero-padded check among numeric
    zeropad = [x for x in numeric if x.startswith('0')]
    # length distribution of numeric
    from collections import Counter
    lens = Counter(len(x) for x in numeric)
    print(f"\n### {name}")
    print(f"  rows={n}  unique={nu}")
    print(f"  numeric ids={len(numeric)}  (zero-padded leading-0: {len(zeropad)})")
    print(f"  canadian/alpha ids={len(canadian)}  other={len(other)}")
    print(f"  numeric length dist: {dict(sorted(lens.items()))}")
    print(f"  numeric examples: {numeric[:4]}")
    print(f"  canadian examples: {canadian[:4]}")
    if other:
        print(f"  OTHER examples: {other[:6]}")
    return ids

# ---------- 1. INVENTORY ----------
print("="*70)
print("1. INVENTORY")
print("="*70)

# metadata
meta = pd.read_csv(f"{SCRATCH}/combined_watershed_metadata_09feb2026.csv", dtype=str)
print(f"\nmetadata total rows: {len(meta)}")
print(f"processing_status value counts:")
print(meta['processing_status'].value_counts(dropna=False).to_string())
meta_success = meta[meta['processing_status'] == 'success']['gage_id'].dropna().astype(str).str.strip()
characterize("metadata processing_status==success", meta_success)

# golden 10feb
gold = load_ids(GOLDEN)
characterize("golden streamflow_signatures_full_10feb2026", gold)

# other signature files
sig_files = {
    "JAN2026": f"{SCRATCH}/s3pull/streamflow_signatures_full_JAN2026.csv",
    "1993-min20yrs": f"{SCRATCH}/s3pull/streamflow_signatures_full_1993-to-2022-min-20yrs.csv",
    "1993-min20yrs_v2": f"{SCRATCH}/s3pull/streamflow_signatures_full_1993-to-2022-min-20yrs_v2.csv",
    "OCT2025_summary": f"{SCRATCH}/s3pull/streamflowSignature_summaryData_OCT2025.csv",
}
sig_ids = {}
for nm, p in sig_files.items():
    ids = load_ids(p)
    sig_ids[nm] = ids
    characterize(f"signatures {nm}", ids)

# QA geometry
geo = load_ids(f"{SCRATCH}/watershed_polygons_26jun2026_qa.csv")
characterize("watershed_polygons_26jun2026_qa (geometry)", geo)

# ---------- 2. LEADING-ZERO HYPOTHESIS ----------
print("\n" + "="*70)
print("2. LEADING-ZERO HYPOTHESIS: metadata-success vs golden-10feb")
print("="*70)

meta_raw = set(meta_success)
gold_raw = set(gold)
raw_int = meta_raw & gold_raw
print(f"\nmetadata-success: {len(meta_raw)} unique")
print(f"golden-10feb: {len(gold_raw)} unique")
print(f"RAW string intersection: {len(raw_int)}")

meta_canon = {canon(x): x for x in meta_success}
gold_canon = {canon(x): x for x in gold}
canon_int = set(meta_canon) & set(gold_canon)
print(f"CANONICAL intersection: {len(canon_int)}")
print(f"  gap closed by canon: {len(canon_int) - len(raw_int)} additional matches")

# example pairs that matched only under canon
new_matches = []
for k in canon_int:
    mv, gv = meta_canon[k], gold_canon[k]
    if mv != gv:
        new_matches.append((mv, gv, k))
print(f"\n  example canon-only match pairs (metadata_raw <-> golden_raw -> canon):")
for mv, gv, k in sorted(new_matches)[:8]:
    print(f"    {mv:>10}  <->  {gv:>10}   (canon={k})")

# What remains: in golden but NOT in metadata-success (canonical)
gold_only = set(gold_canon) - set(meta_canon)
meta_only = set(meta_canon) - set(gold_canon)
print(f"\n  golden-10feb NOT in metadata-success (canon): {len(gold_only)}")
print(f"  metadata-success NOT in golden-10feb (canon): {len(meta_only)}")
# characterize gold_only
go_can = [x for x in gold_only if is_canadian(x)]
go_num = [x for x in gold_only if not is_canadian(x)]
print(f"    golden_only breakdown: US/numeric={len(go_num)}  canadian={len(go_can)}")
print(f"    golden_only examples: {[gold_canon[k] for k in list(gold_only)[:6]]}")

# ---------- 3. MOST-INCLUSIVE UNION ----------
print("\n" + "="*70)
print("3. MOST-INCLUSIVE CANONICAL UNION (all passing/signature sources)")
print("="*70)

# Build canonical sets per source. Use canon key -> representative padded display
def to_canon_set(ids):
    return set(canon(x) for x in ids)

src_canon = {
    "metadata_success": to_canon_set(meta_success),
    "golden_10feb": to_canon_set(gold),
    "JAN2026": to_canon_set(sig_ids["JAN2026"]),
    "1993_min20": to_canon_set(sig_ids["1993-min20yrs"]),
    "1993_min20_v2": to_canon_set(sig_ids["1993-min20yrs_v2"]),
    "OCT2025": to_canon_set(sig_ids["OCT2025_summary"]),
}
for nm, s in src_canon.items():
    print(f"  {nm}: {len(s)} canonical ids")

# Union of EVERYTHING (most inclusive)
union_all = set()
for s in src_canon.values():
    union_all |= s
print(f"\n  UNION of ALL sources (canonical): {len(union_all)}")

# Union of just "passing" definitions the task cares about: metadata-success + signature files
# (metadata-success is the broadest single passing definition)
union_pass_sig = src_canon["metadata_success"] | src_canon["golden_10feb"] | src_canon["JAN2026"] | src_canon["1993_min20"] | src_canon["1993_min20_v2"] | src_canon["OCT2025"]
print(f"  UNION metadata-success + all signature files: {len(union_pass_sig)}")

# US vs Canadian split of union_all
u_can = [x for x in union_all if is_canadian(x)]
u_num = [x for x in union_all if not is_canadian(x)]
print(f"\n  Union split: US/numeric={len(u_num)}  Canadian/alpha={len(u_can)}")

# Contribution breakdown: each canon id - which sources have it
print("\n  Source contribution to union (canonical):")
ms = src_canon["metadata_success"]
# signature-union (any signature file)
sig_union = src_canon["golden_10feb"] | src_canon["JAN2026"] | src_canon["1993_min20"] | src_canon["1993_min20_v2"] | src_canon["OCT2025"]
in_both = ms & sig_union
in_meta_only = ms - sig_union
in_sig_only = sig_union - ms
print(f"    in metadata-success AND >=1 signature file: {len(in_both)}")
print(f"    in metadata-success ONLY (no signature file): {len(in_meta_only)}")
print(f"    in signature file(s) ONLY (not metadata-success): {len(in_sig_only)}")
print(f"    sig-only examples: {list(in_sig_only)[:8]}")

# does metadata-success cover everything? Is metadata the superset?
print(f"\n  Is metadata-success a superset of every signature file (canonical)?")
for nm in ["golden_10feb","JAN2026","1993_min20","1993_min20_v2","OCT2025"]:
    not_covered = src_canon[nm] - ms
    print(f"    {nm}: {len(not_covered)} ids NOT in metadata-success")
    if not_covered:
        print(f"       examples: {list(not_covered)[:6]}")

# ---------- 4. GEOMETRY COVERAGE ----------
print("\n" + "="*70)
print("4. GEOMETRY COVERAGE (watershed_polygons_26jun2026)")
print("="*70)
geo_canon = to_canon_set(geo)
print(f"  geometry ids (canonical): {len(geo_canon)}")
# Recommended inclusive universe = union_all (most inclusive defensible)
universe = union_all
covered = universe & geo_canon
missing = universe - geo_canon
extra_geo = geo_canon - universe
print(f"  recommended universe (union_all): {len(universe)}")
print(f"  covered by geometry: {len(covered)}")
print(f"  NEWLY needed (in universe, no geometry): {len(missing)}")
print(f"  geometry ids NOT in universe (extra polygons): {len(extra_geo)}")
m_can = [x for x in missing if is_canadian(x)]
m_num = [x for x in missing if not is_canadian(x)]
print(f"    missing breakdown: US/numeric={len(m_num)}  Canadian={len(m_can)}")
print(f"    missing examples: {list(missing)[:10]}")

# Also: coverage of metadata-success specifically (geometry was built from 8014 = metadata-success)
cov_ms = ms & geo_canon
print(f"\n  metadata-success covered by geometry: {len(cov_ms)} of {len(ms)}")
print(f"  metadata-success NOT covered: {len(ms - geo_canon)}")
print(f"  geometry NOT in metadata-success: {len(geo_canon - ms)}")

# Save universe + coverage flags
out = pd.DataFrame({"canon_id": sorted(universe)})
out["is_canadian"] = out["canon_id"].apply(is_canadian)
out["in_metadata_success"] = out["canon_id"].isin(ms)
out["in_signature_any"] = out["canon_id"].isin(sig_union)
out["in_golden_10feb"] = out["canon_id"].isin(src_canon["golden_10feb"])
out["has_geometry"] = out["canon_id"].isin(geo_canon)
out.to_csv(f"{SCRATCH}/canonical_universe.csv", index=False)
print(f"\n  wrote canonical_universe.csv ({len(out)} rows)")
print("\nDONE")
