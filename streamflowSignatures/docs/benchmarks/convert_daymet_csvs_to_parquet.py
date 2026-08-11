"""Rebuild daymet_1980_2023.parquet from the 44 annual Daymet CSVs.

Python replacement for `convert_daymet_zip_to_parquet()` (R,
`R/helperFunctions.R`), for machines without R. Written 2026-08-10 after the
climate parquet on the thumbdrive was found truncated (see CHANGELOG -> Known
Issues) and the user restored the source CSVs instead.

This is a STRICTER REIMPLEMENTATION of the same positional reconstruction -- not
an exact mirror (see "Differences from the R version" below). The CSVs carry
site_id, month, year but no day column, so:

    day  = row order within (site_id, year, month)   # data is site-major, day-sequential
    Date = year-month-day
    out  = site_id, Date, prcp, tmin, tmax, swe, vp, srad

`site_id` is read as a STRING and never coerced: the streamflow parquet keys on
zero-padded ids ("01011000"), and numeric coercion would silently break the
climate join for the 6,636 gages whose id has a leading zero. R's `fread()` may
infer these as numeric, which is one reason this is not a byte-for-byte mirror.

Differences from the R version:
  * reads a DIRECTORY of CSVs rather than a ZIP;
  * streams year by year instead of building a single globally sorted in-memory
    table -- ~98M rows x 8 columns does not fit in 16 GB. Each year's rows are
    sorted by (site_id, Date) and appended; pyarrow splits them into multiple row
    groups by size, so the artifact has ~132 row groups, NOT one per year. Row
    order is not semantically meaningful downstream: run_julia_benchmark.jl does
    `groupby(climate, :gage_id)` and joins per gage on (gage_id, date), and
    `preprocess_daily_data()` then re-sorts each gage onto a normalized daily
    grid. Verify end-to-end anyway by reproducing a known product (see below).
  * validation is FATAL here where R only warned (R logs a warning on unexpected
    per-site day counts and does not check per-month counts, input columns, year
    values, invalid dates, or duplicates at all).

DAYMET CALENDAR (verified against all 44 CSVs, 2026-08-10): Daymet publishes a
365-day year. Leap years keep Feb 29 and instead DROP Dec 31, so December has 30
rows and every year has exactly 365 rows per site. The positional reconstruction
therefore yields Dec 1-30 in leap years and no Dec 31 row at all -- a one-day
hole in the climate series every leap year, which the preprocessor absorbs as an
internal gap <= 3 days. This is a property of the source data reproduced
faithfully here, not an artifact of this script; the R version wrote the same
rows (its `expected_days` check assumed 366 and merely logged a warning).

Validation performed per year (all fatal unless noted):
  * required columns present (extra columns are ignored, not rejected)
  * per (site, month) row count == calendar days in that month, except December
    in leap years == 30
  * per-site row count == 365 (all years)
  * no unparseable dates, no duplicate (site_id, Date)
  * site set identical across years (warning only -- Daymet coverage may vary)

WHAT THAT DOES AND DOES NOT ESTABLISH (be precise about this -- a silent date
shift here would corrupt every climate and snow signature):
  * The per-(site, month) count check is a STRUCTURAL CARDINALITY check. It
    catches a month with the wrong number of rows -- e.g. a site short a day in
    February, which would otherwise shift every subsequent date in that month --
    but it CANNOT detect a permutation of the days WITHIN a month.
  * The chronology is instead corroborated END TO END: re-running a delivered
    product's config against the rebuilt input reproduced it with 0 columns
    added/dropped, an identical gage set, and <= 3.4e-13 on 98 of 1,653 columns.
    A within-month permutation could not survive that. Circumstantial support
    also comes from the data itself (high daily lag-1 autocorrelation in
    tmin/tmax/SWE; SWE increases never occurring on zero-precipitation days).

Usage:
    python convert_daymet_csvs_to_parquet.py --csv-dir <dir> --out <parquet>
        [--years 1980-2023] [--compression snappy]

After rebuilding, prove the input reproduces a delivered product before trusting
it for a new one, e.g. re-run the WY 1993-2025 standard config against it and
diff with docs/benchmarks/check_additivity.jl (expect 0 added / 0 dropped and
every column unchanged).
"""
import argparse
import calendar
import os
import re
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

VALUE_COLS = ["prcp", "tmin", "tmax", "swe", "vp", "srad"]
IN_COLS = ["site_id", "month", "year"] + VALUE_COLS
OUT_COLS = ["site_id", "Date"] + VALUE_COLS

SCHEMA = pa.schema(
    [("site_id", pa.string()), ("Date", pa.date32())]
    + [(c, pa.float64()) for c in VALUE_COLS]
)


def parse_years(spec):
    m = re.fullmatch(r"(\d{4})-(\d{4})", spec)
    if not m:
        raise argparse.ArgumentTypeError("--years must look like 1980-2023")
    lo, hi = int(m.group(1)), int(m.group(2))
    if hi < lo:
        raise argparse.ArgumentTypeError("--years end before start")
    return list(range(lo, hi + 1))


def convert_year(path, year):
    """Read one annual CSV and return the validated (site_id, Date, ...) table."""
    # float_precision="round_trip" is REQUIRED, not a nicety: pandas' default C
    # parser is ~1 ULP off on roughly 10% of these values. Measured 2026-08-10 --
    # with the default parser the rebuilt input reproduced the delivered
    # WY 1993-2025 product to 0 added / 0 dropped columns and an identical gage
    # set, but 137 climate-derived columns drifted (max |diff| 3.2e-07, on a
    # pct_change). round_trip cuts that to 98 columns at <= 3.4e-13 -- closer,
    # but NOT bitwise identical. This parse is verified correctly-rounded (0
    # mismatches vs strtod on 1.2M values); the remaining residual is consistent
    # with a last-bit ingestion/serialization difference in the original, and its
    # exact source cannot be isolated because that parquet is unreadable.
    df = pd.read_csv(path, dtype={"site_id": str}, float_precision="round_trip")

    missing = [c for c in IN_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{os.path.basename(path)}: missing columns {missing}")
    if not (df["year"] == year).all():
        bad = sorted(set(df["year"].unique()) - {year})
        raise ValueError(f"{os.path.basename(path)}: unexpected year values {bad}")

    # day = position within (site_id, year, month), preserving file row order.
    # sort=False keeps groups in order of appearance; cumcount is positional and
    # never reorders rows.
    df["day"] = df.groupby(["site_id", "year", "month"], sort=False).cumcount() + 1

    # Structural check: the reconstruction is only valid if each (site, month)
    # block holds exactly the days Daymet publishes for that month -- calendar
    # days, except December in a leap year (Dec 31 is dropped; see module docs).
    expect = {m: calendar.monthrange(year, m)[1] for m in range(1, 13)}
    if calendar.isleap(year):
        expect[12] = 30
    counts = df.groupby(["site_id", "month"], sort=False).size()
    wrong = counts[counts != counts.index.get_level_values("month").map(expect)]
    if len(wrong):
        raise ValueError(
            f"{os.path.basename(path)}: {len(wrong)} (site, month) blocks have a "
            f"row count != calendar days -- day reconstruction is unsafe. "
            f"First offenders:\n{wrong.head(10)}"
        )

    dates = pd.to_datetime(df[["year", "month", "day"]], errors="coerce")
    if dates.isna().any():
        raise ValueError(f"{os.path.basename(path)}: {int(dates.isna().sum())} unparseable dates")
    df["Date"] = dates

    expected_days = 365  # Daymet's 365-day calendar, leap years included
    per_site = df.groupby("site_id", sort=False).size()
    off = per_site[per_site != expected_days]
    if len(off):
        raise ValueError(
            f"{os.path.basename(path)}: {len(off)} sites have != {expected_days} rows. "
            f"First offenders:\n{off.head(10)}"
        )

    out = df[OUT_COLS].sort_values(["site_id", "Date"], kind="mergesort")
    if out.duplicated(["site_id", "Date"]).any():
        raise ValueError(f"{os.path.basename(path)}: duplicate (site_id, Date) rows")

    return out, set(per_site.index)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv-dir", required=True, help="directory of daymet_YYYY.csv files")
    ap.add_argument("--out", required=True, help="output parquet path")
    ap.add_argument("--years", type=parse_years, default=parse_years("1980-2023"))
    ap.add_argument("--compression", default="snappy",
                    help="parquet compression (default snappy, matching the R version)")
    args = ap.parse_args()

    paths = []
    for y in args.years:
        p = os.path.join(args.csv_dir, f"daymet_{y}.csv")
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing {p}")
        if os.path.getsize(p) == 0:
            sys.exit(f"ERROR: zero-byte file {p} (incomplete copy?)")
        paths.append((y, p))

    tmp = args.out + ".tmp"
    total_rows = 0
    site_sets = {}
    writer = pq.ParquetWriter(tmp, SCHEMA, compression=args.compression)
    try:
        for i, (year, path) in enumerate(paths, 1):
            out, sites = convert_year(path, year)
            writer.write_table(pa.Table.from_pandas(out, schema=SCHEMA,
                                                    preserve_index=False))
            total_rows += len(out)
            site_sets[year] = sites
            print(f"  [{i:2d}/{len(paths)}] {year}: {len(out):>9,} rows, "
                  f"{len(sites):,} sites", flush=True)
    finally:
        writer.close()

    base = site_sets[args.years[0]]
    for year, sites in site_sets.items():
        if sites != base:
            print(f"  WARNING: {year} site set differs from {args.years[0]} "
                  f"(+{len(sites - base)} / -{len(base - sites)})")

    os.replace(tmp, args.out)
    print(f"\nWrote {args.out}: {total_rows:,} rows, {len(base):,} sites, "
          f"{os.path.getsize(args.out) / 1024**3:.2f} GB")


if __name__ == "__main__":
    main()
