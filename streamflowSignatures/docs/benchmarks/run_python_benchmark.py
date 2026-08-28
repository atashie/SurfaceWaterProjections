"""
Python Benchmark Runner for Streamflow Signatures

Runs full signature extraction on production data and captures timing.
Mirrors docs/benchmarks/run_julia_benchmark.jl (the canonical runner):
same ENV overrides, same window/qualifying-fraction semantics (read at runtime),
same bounded-memory strategy (column-projected climate read, per-gage join before
preprocessing, raw frames released after Phase 3, preprocess cache evicted per
gage), and the same provenance block in the timing JSON.

Rebuilt 2026-08-24 for the port campaign (Phase 0 —
docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md; Codex findings
F3/F9/F14): the previous version merged the full ~98M-row climate parquet
wholesale into the streamflow frame and had no window/fraction support.
"""

import sys
import os
import gc
import json
import time
import hashlib
import platform
import subprocess
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import pyarrow.parquet as pq_arrow

# Add project root for imports (docs/benchmarks/ -> project root)
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from streamflow_signatures import (
    read_parquet,
    add_water_year_columns,
    filter_qualifying_years,
    preprocess_daily_data,
    calculate_all_signatures,
)
from streamflow_signatures import AnnualCollector
from streamflow_signatures.qa_qc import compute_qa_flags
from streamflow_signatures.metadata import load_gages_ii_interference
from streamflow_signatures.config import (
    MIN_NUM_YEARS,
    MIN_FRAC_GOOD_DATA,
    MIN_Q_VALUE,
    MIN_DAYS_ABOVE_THRESHOLD,
    USE_LEGACY_FILTERING,
    NA_TREND_MIN_FRACTION,
    NA_DECADE_MIN_FRACTION,
    GAGES_II_DIR,
    CHANGEPOINT_ENABLED,
    CP_START_WATER_YEAR,
    CP_END_WATER_YEAR,
    CP_MIN_TOTAL_OBS,
    CP_MIN_SEGMENT_OBS,
    MIN_VALUES_FOR_STATS,
    SAVE_ANNUAL_VALUES,
)

# Configuration — ENV names mirror the Julia runner (STREAMFLOW_*); the pre-2026-08-24
# Python-only names are kept as fallbacks for compatibility.
_DATA_DIR = os.environ.get("STREAMFLOW_PARQUET_DIR", r"D:\processedOuts_feb2026")
STREAMFLOW_PATH = Path(os.environ.get("STREAMFLOW_DATA_PATH",
                       os.environ.get("STREAMFLOW_PATH",
                       os.path.join(_DATA_DIR, "combined_streamflow_data_09feb2026.parquet"))))
CLIMATE_PATH = Path(os.environ.get("STREAMFLOW_CLIMATE_PATH",
                    os.environ.get("CLIMATE_PATH",
                    os.path.join(_DATA_DIR, "daymet_1980_2023.parquet"))))
METADATA_PATH = Path(os.environ.get("STREAMFLOW_METADATA_PATH",
                     os.environ.get("METADATA_PATH",
                     os.path.join(_DATA_DIR, "combined_watershed_metadata_09feb2026.csv"))))
OUTPUT_DIR = Path(os.environ.get("STREAMFLOW_OUTPUT_DIR", str(Path(__file__).parent)))
OUTPUT_PREFIX = os.environ.get("STREAMFLOW_OUTPUT_PREFIX", "python")
GAGES_II_DIR_EFFECTIVE = os.environ.get("STREAMFLOW_GAGES_II_DIR", GAGES_II_DIR)

# Climate columns this runner consumes (post-normalization names): PPT for the
# climate signatures, SWE for snow. Extend here AND in the raw-name projection.
CLIMATE_KEEP = ["gage_id", "date", "PPT", "SWE"]
_CLIMATE_RAW_WANTED = ("gage_id", "site_id", "date", "Date", "PPT", "prcp", "SWE", "swe")


def _env_int(name):
    v = os.environ.get(name, "")
    return int(v) if v != "" else None


def _env_float(name):
    v = os.environ.get(name, "")
    return float(v) if v != "" else None


def _file_info(path, hash_big):
    """Provenance record for one input file (mirrors the Julia runner)."""
    path = Path(path)
    if not path.is_file():
        return {"path": str(path), "exists": False}
    info = {
        "path": str(path.resolve()),
        "exists": True,
        "bytes": path.stat().st_size,
        "mtime": datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%dT%H:%M:%S"),
    }
    if hash_big or path.stat().st_size < 50_000_000:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        info["sha256"] = h.hexdigest()
    return info


def _provenance():
    hash_big = os.environ.get("STREAMFLOW_HASH_INPUTS", "") == "1"
    here = Path(__file__).parent
    try:
        git_rev = subprocess.run(["git", "-C", str(here), "rev-parse", "HEAD"],
                                 capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        git_rev = "unavailable"
    try:
        dirty_out = subprocess.run(["git", "-C", str(here), "status", "--porcelain"],
                                   capture_output=True, text=True, check=True).stdout
        git_dirty = bool(dirty_out.strip())
    except Exception:
        git_dirty = None
    config_file = os.environ.get(
        "STREAMFLOW_CONFIG",
        str(here.parent.parent / "config" / "signatures_config.json"))
    from importlib.metadata import version as _pkg_version
    pkg_versions = {}
    for pkg in ("numpy", "pandas", "scipy", "pyarrow"):
        try:
            pkg_versions[pkg] = _pkg_version(pkg)
        except Exception:
            pkg_versions[pkg] = "unavailable"
    env_keys = ["STREAMFLOW_START_WATER_YEAR", "STREAMFLOW_END_WATER_YEAR",
                "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "STREAMFLOW_CONFIG",
                "STREAMFLOW_GAGES_II_DIR", "STREAMFLOW_OUTPUT_PREFIX",
                "STREAMFLOW_ALLOW_LEGACY"]
    return {
        "streamflow": _file_info(STREAMFLOW_PATH, hash_big),
        "climate": _file_info(CLIMATE_PATH, hash_big),
        "metadata": _file_info(METADATA_PATH, hash_big),
        "config": _file_info(config_file, hash_big),
        "git_revision": git_rev,
        "git_working_tree_dirty": git_dirty,
        "python_version": platform.python_version(),
        "package_versions": pkg_versions,
        "hostname": platform.node(),
        "env_overrides": {k: os.environ[k] for k in env_keys if k in os.environ},
    }


def main():
    # Experiment controls read at RUNTIME (not import) — same discipline as the
    # Julia runner (module constants are frozen at import there too).
    start_water_year = _env_int("STREAMFLOW_START_WATER_YEAR")
    end_water_year = _env_int("STREAMFLOW_END_WATER_YEAR")
    min_qualifying_frac = _env_float("STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION")

    print("=" * 70)
    print("PYTHON BENCHMARK - Streamflow Signatures")
    print("=" * 70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Config: MIN_NUM_YEARS={MIN_NUM_YEARS}, MIN_FRAC_GOOD_DATA={MIN_FRAC_GOOD_DATA}, "
          f"MIN_Q_VALUE={MIN_Q_VALUE}, MIN_DAYS_ABOVE={MIN_DAYS_ABOVE_THRESHOLD}")
    print(f"Changepoint: enabled={CHANGEPOINT_ENABLED}, WY={CP_START_WATER_YEAR}-{CP_END_WATER_YEAR}, "
          f"min_obs={CP_MIN_TOTAL_OBS}, min_seg={CP_MIN_SEGMENT_OBS}")
    print(f"Stats floor: min_values_for_stats={MIN_VALUES_FOR_STATS} (recession/elasticity exempt)")
    print(f"Annual values: save={SAVE_ANNUAL_VALUES}")
    print(f"Experiment: OUTPUT_PREFIX={OUTPUT_PREFIX}, START_WY={start_water_year}, "
          f"END_WY={end_water_year}, MIN_QUAL_FRAC={min_qualifying_frac}")
    print(f"Output dir: {OUTPUT_DIR}")
    print()

    # Changepoint analysis configuration (mirrors the Julia runner's cp_config)
    cp_config = None
    if CHANGEPOINT_ENABLED:
        cp_config = {
            "start_year": CP_START_WATER_YEAR,
            "end_year": CP_END_WATER_YEAR,
            "min_total_obs": CP_MIN_TOTAL_OBS,
            "min_segment_obs": CP_MIN_SEGMENT_OBS,
        }

    # Legacy filtering bypasses the modern preprocessor plumbing (no seasonal
    # flags, trend gates, or window semantics) — a reference/port-campaign run on
    # that path would compare apples to oranges. Fail fast (Codex F14).
    use_legacy = USE_LEGACY_FILTERING
    if use_legacy and os.environ.get("STREAMFLOW_ALLOW_LEGACY") != "1":
        print("ERROR: config has na_handling.use_legacy_filtering=true. This runner "
              "requires the modern preprocessing path; set STREAMFLOW_ALLOW_LEGACY=1 "
              "only for deliberate legacy-compat experiments.")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    timing = {
        "language": "Python",
        "start_time": datetime.now().isoformat(),
        "phases": {}
    }

    # Phase 1: Load streamflow data
    print("Phase 1: Loading streamflow data...")
    t0 = time.perf_counter()

    if not STREAMFLOW_PATH.exists():
        print(f"ERROR: Streamflow file not found: {STREAMFLOW_PATH}")
        return 1

    streamflow = read_parquet(str(STREAMFLOW_PATH))  # auto-normalizes columns

    if "water_year" not in streamflow.columns or "dowy" not in streamflow.columns:
        streamflow = add_water_year_columns(streamflow, date_col="date")

    t1 = time.perf_counter()
    timing["phases"]["load_streamflow"] = t1 - t0
    print(f"  Loaded {len(streamflow):,} rows in {t1-t0:.2f}s")
    print(f"  Unique gages: {streamflow['gage_id'].nunique()}")

    # Phase 2: Load climate data (column-projected; grouped per gage, NOT merged
    # wholesale — the full 8-column ~98M-row frame + a 111M-row merge OOMs on
    # 16 GB machines; mirrors the Julia runner's select!-then-group strategy)
    print("\nPhase 2: Loading climate data...")
    t0 = time.perf_counter()

    has_climate = False
    climate_by_gage = None
    if CLIMATE_PATH.exists():
        schema_names = set(pq_arrow.read_schema(str(CLIMATE_PATH)).names)
        raw_cols = [c for c in _CLIMATE_RAW_WANTED if c in schema_names]
        climate = read_parquet(str(CLIMATE_PATH), columns=raw_cols)  # normalizes names
        keep = [c for c in CLIMATE_KEEP if c in climate.columns]
        climate = climate[keep]

        t1 = time.perf_counter()
        print(f"  Loaded {len(climate):,} rows (columns: {list(climate.columns)}) in {t1-t0:.2f}s")
        print(f"  SWE column: {'present (snow metrics enabled)' if 'SWE' in climate.columns else 'absent (snow metrics skipped)'}")

        print("  Grouping climate data by gage...")
        climate_by_gage = {gid: frame for gid, frame in climate.groupby("gage_id", sort=False)}
        del climate
        gc.collect()

        has_climate = True
        t2 = time.perf_counter()
        timing["phases"]["load_climate"] = t2 - t0
        print(f"  Grouped {len(climate_by_gage)} climate gages in {t2-t1:.2f}s")
    else:
        print("  Climate file not found, skipping climate signatures")
        timing["phases"]["load_climate"] = 0

    # Phase 2b: Area-normalization lookup (Q-to-PPT gate).
    # Gages with no drainage area (HYDAT gap -- canals, dam outflows, channel
    # splits) carry Q in raw m3/s; Q-to-PPT signatures are skipped for them
    # because Q and PPT units don't match. Keys are leading-zero-stripped for
    # numeric ids (metadata stores USGS ids without leading zeros).
    def _normalize_join_id(gid):
        gid = str(gid)
        if gid.isdigit():
            return gid.lstrip("0") or "0"
        return gid

    area_norm_lookup = {}
    if METADATA_PATH.exists():
        # Callable usecols tolerates metadata files that predate the column
        _meta_an = pd.read_csv(
            METADATA_PATH, usecols=lambda c: c in ("gage_id", "area_normalized"))
        if {"gage_id", "area_normalized"} <= set(_meta_an.columns):
            for _gid, _an in zip(_meta_an["gage_id"], _meta_an["area_normalized"]):
                # Only an explicit false-like value gates (bool False, numeric 0,
                # "FALSE"/"0" strings — harmonized with the Julia/rpkg runners);
                # missing defaults to normalized
                if isinstance(_an, str):
                    explicit_false = _an.strip().upper() in ("FALSE", "0")
                else:
                    explicit_false = pd.notna(_an) and float(_an) == 0.0
                area_norm_lookup[_normalize_join_id(_gid)] = not explicit_false
            _n_unnorm = sum(1 for v in area_norm_lookup.values() if not v)
            print(f"  Area-normalization lookup: {len(area_norm_lookup)} gages, "
                  f"{_n_unnorm} not normalized (Q-to-PPT signatures will be skipped for these)")
        else:
            print("  Metadata has no area_normalized column -- assuming all gages area-normalized")

    # Phase 3: Per-year quality filtering
    print(f"\nPhase 3: Per-year quality filtering (legacy={use_legacy})...")
    t0 = time.perf_counter()

    grouped = streamflow.groupby("gage_id", sort=False)
    qualifying_entries = []  # list of (gage_id, qualifying water years)
    total_gages = 0

    # Cache preprocessing results for non-legacy path (evicted per gage in Phase 4)
    preprocess_cache = {}

    if use_legacy:
        for gage_id, gage_data in grouped:
            total_gages += 1
            qual_years, qualifies = filter_qualifying_years(gage_data)
            if qualifies:
                qualifying_entries.append((gage_id, qual_years))
    else:
        for gage_id, gage_data in grouped:
            total_gages += 1
            if total_gages % 500 == 0:
                print(f"  Preprocessing gage {total_gages}...")

            # Merge this gage's climate slice BEFORE preprocessing so PPT is
            # available for valid_climate_years (per-gage join: bounded memory)
            gage_df = gage_data
            if climate_by_gage is not None:
                cslice = climate_by_gage.get(gage_id)
                if cslice is not None:
                    gage_df = gage_df.merge(cslice, on=["gage_id", "date"], how="left")

            result = preprocess_daily_data(gage_df)

            # Apply water-year window filter if configured (start and/or end) —
            # semantics mirror the Julia runner: preprocess the FULL record, then
            # filter the result (never pre-filter raw rows; interpolation and
            # rejection must see the same series in every language)
            if start_water_year is not None or end_water_year is not None:
                wy_lo = start_water_year if start_water_year is not None else -(10**9)
                wy_hi = end_water_year if end_water_year is not None else 10**9
                rdata = result["data"]
                result = {
                    "data": rdata[(rdata["water_year"] >= wy_lo) & (rdata["water_year"] <= wy_hi)],
                    "valid_years": [yr for yr in result["valid_years"] if wy_lo <= yr <= wy_hi],
                    "valid_climate_years": [yr for yr in result["valid_climate_years"] if wy_lo <= yr <= wy_hi],
                    "rejected_years": result["rejected_years"],
                    "seasonal_flags": [d for d in result["seasonal_flags"] if wy_lo <= d["water_year"] <= wy_hi],
                    "diagnostics": [d for d in result["diagnostics"] if wy_lo <= d["water_year"] <= wy_hi],
                    # valid_swe_years is present whenever the preprocessor saw SWE
                    **({"valid_swe_years": [yr for yr in result["valid_swe_years"] if wy_lo <= yr <= wy_hi]}
                       if "valid_swe_years" in result else {}),
                }

            # Apply min qualifying data fraction filter if configured.
            # Denominator = per-gage possible years within the [start, end] window:
            # window start .. the gage's last in-window year (window-capped) —
            # identical to the Julia runner.
            if min_qualifying_frac is not None:
                all_wy = pd.unique(gage_df["water_year"].dropna())
                start_yr = start_water_year if start_water_year is not None else int(min(all_wy))
                end_yr = end_water_year if end_water_year is not None else 10**9
                wy_in_range = [int(yr) for yr in all_wy if start_yr <= yr <= end_yr]
                if wy_in_range:
                    total_possible = max(wy_in_range) - start_yr + 1
                    frac = len(result["valid_years"]) / total_possible
                    if frac < min_qualifying_frac:
                        continue  # Skip this gage

            if len(result["valid_years"]) >= MIN_NUM_YEARS:
                qualifying_entries.append((gage_id, result["valid_years"]))
                preprocess_cache[gage_id] = result

    print(f"  Total gages: {total_gages}")
    print(f"  Qualifying gages: {len(qualifying_entries)}")

    t1 = time.perf_counter()
    timing["phases"]["filter_gages"] = t1 - t0

    # Memory: the non-legacy Phase 4 consumes ONLY preprocess_cache — release the
    # raw frames + groupings before per-gage processing (mirrors the Julia runner)
    if not use_legacy:
        del grouped
        del streamflow
        climate_by_gage = None
        gc.collect()

    # Phase 4: Process signatures
    n_gages = len(qualifying_entries)
    print(f"\nPhase 4: Processing {n_gages} gages...")
    t0 = time.perf_counter()

    all_results = []

    # Annual-values accumulators (long format; one collector per gage, drained
    # here with the gage_id tag) — mirrors the Julia runner.
    save_annual = SAVE_ANNUAL_VALUES
    annual_gage_id, annual_signature, annual_water_year, annual_value = [], [], [], []

    if use_legacy:
        # Legacy path (deprecated compat mode; blocked above unless explicitly
        # allowed): raw data + per-year filter, no preprocessing semantics
        qual_map = dict(qualifying_entries)
        for i, (gage_id, gage_data) in enumerate(grouped):
            if gage_id not in qual_map:
                continue
            qual_years = qual_map[gage_id]
            gage_area_normalized = area_norm_lookup.get(_normalize_join_id(gage_id), True)
            if qual_years:
                gage_data = gage_data[gage_data["water_year"].isin(qual_years)]
            signatures = calculate_all_signatures(gage_data, has_climate,
                                                  area_normalized=gage_area_normalized)
            signatures["gage_id"] = gage_id
            signatures["start_water_year"] = min(qual_years) if qual_years else np.nan
            signatures["end_water_year"] = max(qual_years) if qual_years else np.nan
            signatures["num_water_years"] = len(qual_years)
            all_results.append(signatures)
    else:
        for i, (gage_id, qual_years) in enumerate(qualifying_entries):
            if (i + 1) % 100 == 0 or i == 0:
                elapsed = time.perf_counter() - t0
                rate = (i + 1) / elapsed if elapsed > 0 else 0
                eta = (n_gages - i - 1) / rate if rate > 0 else 0
                print(f"  [{i+1}/{n_gages}] Processing... ({rate:.1f} gages/s, ETA: {eta/60:.1f} min)")

            gage_area_normalized = area_norm_lookup.get(_normalize_join_id(gage_id), True)
            gage_collector = AnnualCollector() if save_annual else None

            # Evict as we go so the cache shrinks as Phase 4 progresses
            pp = preprocess_cache.pop(gage_id)
            pp_data = pp["data"]
            gage_has_climate = has_climate and len(pp.get("valid_climate_years", [])) > 0

            # Filter to valid_climate_years for climate signatures
            climate_data = None
            if gage_has_climate:
                climate_yr_set = set(pp["valid_climate_years"])
                climate_data = pp_data[pp_data["water_year"].isin(climate_yr_set)]

            # Snow data: explicit opt-in frame filtered to SWE-valid years.
            # Possibly 0 rows (gage has SWE but no valid SWE year -> full NaN
            # snow key set); None when the gage has no SWE column at all ->
            # snow keys absent. NO implicit fallback to pp_data in the
            # orchestrator (mirrors the Julia runner).
            snow_data = None
            if "SWE" in pp_data.columns:
                swe_yr_set = set(pp.get("valid_swe_years", []))
                snow_data = pp_data[pp_data["water_year"].isin(swe_yr_set)]
            signatures = calculate_all_signatures(
                pp_data, gage_has_climate,
                seasonal_flags=pp.get("seasonal_flags"),
                trend_completeness=NA_TREND_MIN_FRACTION,
                decade_completeness=NA_DECADE_MIN_FRACTION,
                climate_data=climate_data,
                area_normalized=gage_area_normalized,
                changepoint=cp_config,
                min_values_for_stats=MIN_VALUES_FOR_STATS,
                collector=gage_collector,
                snow_data=snow_data,
                snow_climate_years=pp.get("valid_climate_years"),
                gage_id=str(gage_id),
            )

            # Drain this gage's annual values into the global accumulators
            if gage_collector is not None and len(gage_collector) > 0:
                n_rows = len(gage_collector)
                annual_gage_id.extend([str(gage_id)] * n_rows)
                annual_signature.extend(gage_collector.signature)
                annual_water_year.extend(gage_collector.water_year)
                annual_value.extend(gage_collector.value)

            signatures["gage_id"] = gage_id
            signatures["start_water_year"] = min(qual_years) if qual_years else np.nan
            signatures["end_water_year"] = max(qual_years) if qual_years else np.nan
            signatures["num_water_years"] = len(qual_years)

            # Per-gage ice-affected day total from diagnostics
            signatures["ice_affected_days_total"] = sum(
                d.get("na_cause_ice", 0) for d in pp.get("diagnostics", []))

            all_results.append(signatures)

    t1 = time.perf_counter()
    timing["phases"]["process_signatures"] = t1 - t0
    print(f"  Processed {len(all_results)} gages in {t1-t0:.2f}s")

    # Zero-gage guard (mirrors the Julia runner)
    if not all_results:
        print("WARNING: Zero gages qualified. Writing empty output.")
        pd.DataFrame({"gage_id": []}).to_csv(OUTPUT_DIR / f"{OUTPUT_PREFIX}_signatures.csv", index=False)
        return 0

    # Phase 4b: Write annual values parquet (before the metadata merge, so a
    # metadata failure cannot lose it). Long format, deterministically sorted —
    # mirrors the Julia runner.
    if save_annual and annual_value:
        print("\nPhase 4b: Writing annual values parquet...")
        t0 = time.perf_counter()
        annual_df = pd.DataFrame({
            "gage_id": annual_gage_id,
            "signature": annual_signature,
            "water_year": np.asarray(annual_water_year, dtype="int32"),
            "value": annual_value,
        }).sort_values(["gage_id", "signature", "water_year"], kind="stable")
        annual_path = OUTPUT_DIR / f"{OUTPUT_PREFIX}_signatures_annual.parquet"
        annual_df.to_parquet(annual_path, index=False, compression="zstd")
        t1 = time.perf_counter()
        timing["phases"]["write_annual_values"] = t1 - t0
        timing["n_annual_rows"] = len(annual_df)
        print(f"  Saved {len(annual_df):,} rows ({annual_df['signature'].nunique()} signatures) "
              f"to {annual_path}")
        print(f"  Size: {annual_path.stat().st_size/1e6:.1f} MB, took {t1-t0:.2f}s")

    # Phase 5: Merge metadata and compute QA/QC flags
    print("\nPhase 5: Merging metadata and computing QA/QC flags...")
    t0 = time.perf_counter()

    results_df = pd.DataFrame(all_results)

    # Load and merge metadata
    if METADATA_PATH.exists():
        print("  Loading metadata...")
        metadata = pd.read_csv(METADATA_PATH)

        # Metadata stores USGS IDs without leading zeros; benchmark output keeps
        # them. Join on a normalized key (numeric ids leading-zero-stripped).
        results_df["gage_id"] = results_df["gage_id"].astype(str)
        metadata["gage_id"] = metadata["gage_id"].astype(str)
        results_df["_join_id"] = results_df["gage_id"].map(_normalize_join_id)
        metadata["_join_id"] = metadata["gage_id"].map(_normalize_join_id)

        # Rename basin_area_km2 to basin_area to match R output
        if "basin_area_km2" in metadata.columns:
            metadata = metadata.rename(columns={"basin_area_km2": "basin_area"})

        metadata_cols = [
            "_join_id",
            "latitude",
            "longitude",
            "basin_area",
            "gage_type",
            "area_normalized",
        ]
        available_meta_cols = [c for c in metadata_cols if c in metadata.columns]

        if len(available_meta_cols) > 1:  # More than just the join key
            print(f"  Merging {len(available_meta_cols)-1} metadata columns...")
            results_df = results_df.merge(
                metadata[available_meta_cols],
                on="_join_id",
                how="left"
            )
            n_matched = results_df["latitude"].notna().sum() if "latitude" in results_df.columns else 0
            print(f"  Matched {n_matched} / {len(results_df)} gages with metadata")
        else:
            print("  Warning: No metadata columns found to merge")
        results_df = results_df.drop(columns=["_join_id"])

        # Enrich with human interference metadata from GAGES-II
        print("  Loading GAGES-II interference metadata...")
        gages_ii = load_gages_ii_interference(GAGES_II_DIR_EFFECTIVE)
        if len(gages_ii) > 0 and "STAID" in gages_ii.columns:
            gages_ii["STAID"] = gages_ii["STAID"].astype(str)
            interference_cols = [
                "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                "HYDRO_DISTURB_INDX", "CLASS",
            ]
            avail_int_cols = [c for c in interference_cols if c in gages_ii.columns]
            gages_ii_subset = gages_ii[["STAID"] + avail_int_cols].rename(
                columns={"STAID": "gage_id"}
            )
            results_df = results_df.merge(gages_ii_subset, on="gage_id", how="left")
            n_matched = results_df[avail_int_cols[0]].notna().sum() if avail_int_cols else 0
            print(f"  Matched {n_matched} gages with GAGES-II data")

            # Compute human_interference_class from CLASS (USGS) and gage_type (Canadian)
            def classify_interference(row):
                if pd.notna(row.get("CLASS")):
                    cls = str(row["CLASS"]).strip()
                    if cls == "Ref":
                        return "reference"
                    elif cls == "Non-ref":
                        return "non-reference"
                return "unknown"

            results_df["human_interference_class"] = results_df.apply(
                classify_interference, axis=1
            )

            # Canadian HYDAT interference (RHBN/REGULATED) from the pre-exported
            # CSV — the same file the Julia runner reads via
            # load_canadian_interference() (the old "not available from Python"
            # behavior was stale). hic derivation mirrors julia/src/metadata.jl:
            # RHBN true -> reference, false -> non-reference, else unknown.
            results_df["RHBN"] = np.nan
            results_df["REGULATED"] = np.nan
            hydat_csv = Path(os.environ.get(
                "STREAMFLOW_HYDAT_PATH",
                str(Path(__file__).parent.parent.parent / "metadata" / "canadian_hydat_interference.csv")))
            if hydat_csv.is_file():
                can = pd.read_csv(hydat_csv)
                if "STATION_NUMBER" in can.columns:
                    can = can.rename(columns={"STATION_NUMBER": "gage_id"})
                    can["gage_id"] = can["gage_id"].astype(str)
                    can = can.drop_duplicates(subset="gage_id", keep="last").set_index("gage_id")
                    mask = results_df["gage_id"].isin(can.index)
                    idx = results_df.loc[mask, "gage_id"]
                    rhbn_vals = can.loc[idx, "RHBN"].values
                    results_df["RHBN"] = results_df["RHBN"].astype(object)
                    results_df["REGULATED"] = results_df["REGULATED"].astype(object)
                    results_df.loc[mask, "RHBN"] = rhbn_vals
                    results_df.loc[mask, "REGULATED"] = can.loc[idx, "REGULATED"].values
                    hic = np.where(rhbn_vals == True, "reference",  # noqa: E712
                          np.where(rhbn_vals == False, "non-reference", "unknown"))  # noqa: E712
                    results_df.loc[mask, "human_interference_class"] = hic
                    print(f"  Matched {int(mask.sum())} Canadian gages with HYDAT metadata")
            else:
                print(f"  Warning: Canadian HYDAT CSV not found: {hydat_csv}")
            print(f"  Interference columns added: {len(avail_int_cols) + 3}")
        else:
            print("  Warning: GAGES-II data not available")
    else:
        print(f"  Warning: Metadata file not found: {METADATA_PATH}")

    # Compute QA/QC flags
    print("  Computing QA/QC flags...")
    results_df = compute_qa_flags(results_df)

    # Organize columns: gage_id first, then metadata, then interference, then signatures, then flags
    metadata_order = [
        "gage_id", "latitude", "longitude", "basin_area",
        "gage_type", "num_water_years", "start_water_year", "end_water_year",
        "area_normalized",
        "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
        "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
        "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
        "human_interference_class",
    ]
    flag_cols = [c for c in results_df.columns if c.startswith("flagged_")]
    signature_cols = [c for c in results_df.columns
                      if c not in metadata_order and c not in flag_cols]

    final_cols = []
    for c in metadata_order:
        if c in results_df.columns:
            final_cols.append(c)
    final_cols.extend(sorted(signature_cols))
    final_cols.extend(sorted(flag_cols))

    results_df = results_df[final_cols]

    output_path = OUTPUT_DIR / f"{OUTPUT_PREFIX}_signatures.csv"
    results_df.to_csv(output_path, index=False)

    t1 = time.perf_counter()
    timing["phases"]["metadata_qaqc_save"] = t1 - t0
    print(f"  Saved to {output_path}")
    print(f"  Shape: {results_df.shape}")

    # Calculate totals
    timing["end_time"] = datetime.now().isoformat()
    timing["total_seconds"] = sum(timing["phases"].values())
    timing["n_gages_processed"] = len(all_results)
    n_meta = len([c for c in metadata_order if c in results_df.columns])
    n_flags = len([c for c in results_df.columns if c.startswith("flagged_")])
    timing["n_signature_columns"] = len(results_df.columns) - n_meta - n_flags
    timing["n_metadata_columns"] = n_meta
    timing["n_qaqc_flags"] = n_flags
    timing["provenance"] = _provenance()

    timing_path = OUTPUT_DIR / f"{OUTPUT_PREFIX}_timing.json"
    with open(timing_path, "w") as f:
        json.dump(timing, f, indent=2)
    print(f"  Timing saved to {timing_path}")

    # Summary
    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print(f"Total time: {timing['total_seconds']:.2f}s ({timing['total_seconds']/60:.2f} min)")
    print(f"Gages processed: {timing['n_gages_processed']}")
    print(f"Total columns: {len(results_df.columns)}")
    print(f"  Signature columns: {timing['n_signature_columns']}")
    print(f"  Metadata columns: {timing['n_metadata_columns']}")
    print(f"  QA/QC flag columns: {timing['n_qaqc_flags']}")
    print(f"Rate: {timing['n_gages_processed']/timing['total_seconds']:.2f} gages/s")

    return 0


if __name__ == "__main__":
    sys.exit(main())
