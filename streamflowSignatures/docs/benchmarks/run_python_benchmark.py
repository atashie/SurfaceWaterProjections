"""
Python Benchmark Runner for Streamflow Signatures

Runs full signature extraction on production data and captures timing.
Uses per-year quality filtering matching R's process_signatures_from_parquet().
"""

import sys
import json
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np

# Add project root for imports (docs/benchmarks/ -> project root)
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from streamflow_signatures import (
    read_parquet,
    add_water_year_columns,
    filter_qualifying_years,
    calculate_all_signatures,
)
from streamflow_signatures.qa_qc import compute_qa_flags
from streamflow_signatures.metadata import load_gages_ii_interference
from streamflow_signatures.config import (
    MIN_NUM_YEARS,
    MIN_FRAC_GOOD_DATA,
    MIN_Q_VALUE,
    MIN_DAYS_ABOVE_THRESHOLD,
    GAGES_II_DIR,
)

# Configuration — override with environment variables or edit defaults below
import os
_DATA_DIR = os.environ.get("STREAMFLOW_PARQUET_DIR", r"D:\processedOuts_feb2026")
STREAMFLOW_PATH = Path(os.environ.get("STREAMFLOW_PATH", os.path.join(_DATA_DIR, "combined_streamflow_data_09feb2026.parquet")))
CLIMATE_PATH = Path(os.environ.get("CLIMATE_PATH", os.path.join(_DATA_DIR, "daymet_1980_2023.parquet")))
METADATA_PATH = Path(os.environ.get("METADATA_PATH", os.path.join(_DATA_DIR, "combined_watershed_metadata_09feb2026.csv")))
OUTPUT_DIR = Path(__file__).parent



def main():
    print("=" * 70)
    print("PYTHON BENCHMARK - Streamflow Signatures")
    print("=" * 70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Config: MIN_NUM_YEARS={MIN_NUM_YEARS}, MIN_FRAC_GOOD_DATA={MIN_FRAC_GOOD_DATA}, "
          f"MIN_Q_VALUE={MIN_Q_VALUE}, MIN_DAYS_ABOVE={MIN_DAYS_ABOVE_THRESHOLD}")
    print()

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
        sys.exit(1)

    streamflow = read_parquet(str(STREAMFLOW_PATH))  # auto-normalizes columns

    # Check if water_year columns already exist
    if "water_year" not in streamflow.columns:
        streamflow = add_water_year_columns(streamflow, date_col="date")
    else:
        # Ensure dowy column exists (day of water year)
        if "dowy" not in streamflow.columns and "date" in streamflow.columns:
            streamflow = add_water_year_columns(streamflow, date_col="date")

    t1 = time.perf_counter()
    timing["phases"]["load_streamflow"] = t1 - t0
    print(f"  Loaded {len(streamflow):,} rows in {t1-t0:.2f}s")
    print(f"  Unique gages: {streamflow['gage_id'].nunique()}")

    # Phase 2: Load and merge climate data
    print("\nPhase 2: Loading climate data...")
    t0 = time.perf_counter()

    has_climate = False
    if CLIMATE_PATH.exists():
        climate = read_parquet(str(CLIMATE_PATH))  # auto-normalizes columns

        t1 = time.perf_counter()
        print(f"  Loaded {len(climate):,} rows in {t1-t0:.2f}s")

        # Pre-merge climate data with streamflow (much faster than per-gage merge)
        print("  Merging climate data with streamflow...")
        t_merge = time.perf_counter()

        # Only keep needed columns from climate
        climate_cols = climate[["gage_id", "date", "PPT"]].copy()

        # Merge
        streamflow = pd.merge(
            streamflow,
            climate_cols,
            on=["gage_id", "date"],
            how="left"
        )

        t2 = time.perf_counter()
        print(f"  Merged in {t2-t_merge:.2f}s")

        has_climate = True
        timing["phases"]["load_climate"] = t2 - t0
    else:
        print(f"  Climate file not found, skipping climate signatures")
        timing["phases"]["load_climate"] = 0

    # Phase 3: Per-year quality filtering (matching R's three-stage filter)
    print("\nPhase 3: Per-year quality filtering...")
    t0 = time.perf_counter()

    # Group by gage for per-gage, per-year filtering
    grouped = streamflow.groupby("gage_id", sort=False)
    qualifying_gages = []
    gage_qualifying_years = {}  # gage_id -> list of qualifying water years
    total_gages = 0

    for gage_id, gage_data in grouped:
        total_gages += 1
        qual_years, qualifies = filter_qualifying_years(gage_data)
        if qualifies:
            qualifying_gages.append(gage_id)
            gage_qualifying_years[gage_id] = qual_years

    print(f"  Total gages: {total_gages}")
    print(f"  Qualifying gages: {len(qualifying_gages)}")

    # Filter streamflow to only qualifying gages (reduces data size for groupby)
    streamflow = streamflow[streamflow["gage_id"].isin(qualifying_gages)]
    print(f"  Filtered streamflow rows: {len(streamflow):,}")

    t1 = time.perf_counter()
    timing["phases"]["filter_gages"] = t1 - t0

    # Phase 4: Process signatures using groupby (much faster than per-gage filtering)
    print(f"\nPhase 4: Processing {len(qualifying_gages)} gages...")
    t0 = time.perf_counter()

    all_results = []
    n_gages = len(qualifying_gages)

    # Use groupby for efficient iteration
    grouped = streamflow.groupby("gage_id", sort=False)

    for i, (gage_id, gage_data) in enumerate(grouped):
        if (i + 1) % 100 == 0 or i == 0:
            elapsed = time.perf_counter() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (n_gages - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{n_gages}] Processing... ({rate:.1f} gages/s, ETA: {eta/60:.1f} min)")

        # Filter to qualifying water years only (matching R's per-year filter)
        qual_years = gage_qualifying_years.get(gage_id, [])
        if qual_years:
            gage_data = gage_data[gage_data["water_year"].isin(qual_years)]

        # Calculate signatures
        signatures = calculate_all_signatures(gage_data, has_climate)
        signatures["gage_id"] = gage_id

        # Add computed metadata columns (matching R output)
        signatures["start_water_year"] = min(qual_years) if qual_years else np.nan
        signatures["end_water_year"] = max(qual_years) if qual_years else np.nan
        signatures["num_water_years"] = len(qual_years)

        all_results.append(signatures)

    t1 = time.perf_counter()
    timing["phases"]["process_signatures"] = t1 - t0
    print(f"  Processed {len(all_results)} gages in {t1-t0:.2f}s")

    # Phase 5: Merge metadata and compute QA/QC flags
    print("\nPhase 5: Merging metadata and computing QA/QC flags...")
    t0 = time.perf_counter()

    results_df = pd.DataFrame(all_results)

    # Load and merge metadata
    if METADATA_PATH.exists():
        print("  Loading metadata...")
        metadata = pd.read_csv(METADATA_PATH)

        # Ensure gage_id is string in both DataFrames for merging
        results_df["gage_id"] = results_df["gage_id"].astype(str)
        metadata["gage_id"] = metadata["gage_id"].astype(str)

        # Rename basin_area_km2 to basin_area to match R output
        if "basin_area_km2" in metadata.columns:
            metadata = metadata.rename(columns={"basin_area_km2": "basin_area"})

        # Select basic metadata columns to merge
        metadata_cols = [
            "gage_id",
            "latitude",
            "longitude",
            "basin_area",
            "gage_type",
            "area_normalized",
        ]

        # Only keep columns that exist in metadata
        available_meta_cols = [c for c in metadata_cols if c in metadata.columns]

        if len(available_meta_cols) > 1:  # More than just gage_id
            print(f"  Merging {len(available_meta_cols)-1} metadata columns...")
            results_df = results_df.merge(
                metadata[available_meta_cols],
                on="gage_id",
                how="left"
            )
        else:
            print("  Warning: No metadata columns found to merge")

        # Enrich with human interference metadata from GAGES-II
        print("  Loading GAGES-II interference metadata...")
        gages_ii = load_gages_ii_interference(GAGES_II_DIR)
        if len(gages_ii) > 0 and "STAID" in gages_ii.columns:
            gages_ii["STAID"] = gages_ii["STAID"].astype(str)
            # Match by gage_id (USGS gages have numeric IDs matching STAID)
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

            # Add empty RHBN/REGULATED columns (Canadian HYDAT not available from Python)
            results_df["RHBN"] = np.nan
            results_df["REGULATED"] = np.nan
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

    # Build final column order
    final_cols = []
    for c in metadata_order:
        if c in results_df.columns:
            final_cols.append(c)
    final_cols.extend(sorted(signature_cols))
    final_cols.extend(sorted(flag_cols))

    results_df = results_df[final_cols]

    output_path = OUTPUT_DIR / "python_signatures.csv"
    results_df.to_csv(output_path, index=False)

    t1 = time.perf_counter()
    timing["phases"]["metadata_qaqc_save"] = t1 - t0
    print(f"  Saved to {output_path}")
    print(f"  Shape: {results_df.shape}")

    # Calculate totals
    timing["end_time"] = datetime.now().isoformat()
    timing["total_seconds"] = sum(timing["phases"].values())
    timing["n_gages_processed"] = len(all_results)
    # Count columns: exclude gage_id, metadata cols, and flag cols
    n_meta = len([c for c in metadata_order if c in results_df.columns])
    n_flags = len([c for c in results_df.columns if c.startswith("flagged_")])
    timing["n_signature_columns"] = len(results_df.columns) - n_meta - n_flags
    timing["n_metadata_columns"] = n_meta
    timing["n_qaqc_flags"] = n_flags

    # Save timing
    timing_path = OUTPUT_DIR / "python_timing.json"
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
