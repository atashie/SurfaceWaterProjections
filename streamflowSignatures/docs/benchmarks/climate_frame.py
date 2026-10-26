"""climate_frame.py — the production climate-frame helper of the Python runner
(Phase A′ of docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.2; decisions
12–13, user 2026-10-26). Mirrors climate_frame.jl (canonical Julia runner) and
climate_frame.R (rpkg runner); self-test: selftest_climate_frame.py.

The runner joins each gage's Daymet rows onto that gage's streamflow ROWS before the
library builds the water-year grid, so a date with no streamflow row reached
`preprocess_daily_data` without climate and was interpolated (CHANGELOG Known Issue
"climate rides on the streamflow dates", 2026-10-16). `frame_with_climate` builds the
frame under the runtime switch STREAMFLOW_CLIMATE_FRAME:
  left   today's left merge, bit for bit — the attribution baseline;
  outer  the left merge PLUS the gage's Daymet rows on dates with no streamflow row,
         restricted to the streamflow record span AND to water years that hold at least
         one streamflow row (whole interior water years without streamflow must not enter
         the preprocessor's year loop — Codex round 3).
Added rows carry gage_id, date, water_year / month / dowy (add_water_year_columns; year
and doy from the date when the frame has them), the climate columns, NaN for Q and
every other numeric column and None for the rest; they are appended after the merged
rows (the preprocessor re-grids every water year by date). Default mode `left` until
the rerun's products land (then `outer`, Phase H).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))
from streamflow_signatures import add_water_year_columns  # noqa: E402

CLIMATE_FRAME_MODES = ("left", "outer")


def climate_frame_mode():
    """STREAMFLOW_CLIMATE_FRAME at runtime; "left" when unset; any other value errors."""
    v = os.environ.get("STREAMFLOW_CLIMATE_FRAME", "left").strip().lower()
    if v not in CLIMATE_FRAME_MODES:
        raise ValueError(f"STREAMFLOW_CLIMATE_FRAME must be one of {CLIMATE_FRAME_MODES}; got {v!r}")
    return v


def _check_dates(df, name):
    """Both frames must carry tz-naive datetime64 dates at midnight: a time-of-day or a
    tz-aware column would match no climate row in the join and silently yield a
    climate-free frame (review 2026-10-26)."""
    if "date" not in df.columns:
        raise KeyError(f"frame_with_climate: the {name} frame has no date column")
    d = df["date"]
    if not pd.api.types.is_datetime64_any_dtype(d):
        raise TypeError(f"frame_with_climate: the {name} frame's date column is {d.dtype}, expected datetime64")
    if getattr(d.dt, "tz", None) is not None:
        raise TypeError(f"frame_with_climate: the {name} frame's date column is tz-aware")
    if len(d) and (d.isna().any() or not (d == d.dt.normalize()).all()):
        raise ValueError(f"frame_with_climate: the {name} frame's date column must be complete and at midnight")


def frame_with_climate(gage_df, cslice, mode):
    """gage_df: one gage's streamflow rows (gage_id, date, Q, water_year, ...); cslice:
    that gage's climate rows (gage_id, date, PPT[, SWE, TMIN, TMAX]); mode: left | outer."""
    _check_dates(gage_df, "streamflow")
    _check_dates(cslice, "climate")
    if mode == "left":
        return gage_df.merge(cslice, on=["gage_id", "date"], how="left")
    if mode != "outer":
        raise ValueError(f"frame_with_climate: unknown mode {mode!r}")
    joined = gage_df.merge(cslice, on=["gage_id", "date"], how="left")
    if len(gage_df) == 0 or len(cslice) == 0:
        return joined
    if "water_year" not in gage_df.columns:
        raise KeyError("frame_with_climate: the streamflow frame has no water_year column")

    sf_dates = gage_df["date"]
    lo, hi = sf_dates.min(), sf_dates.max()
    wy_set = set(gage_df["water_year"].dropna().astype(int).tolist())
    c_dates = cslice["date"]
    in_span = (c_dates >= lo) & (c_dates <= hi)
    if not in_span.any():
        return joined

    cand = cslice.loc[in_span].copy()
    # duplicate climate dates: identical duplicate rows collapse to one, conflicting ones
    # are an error — an appended duplicate would put two rows on one grid day
    if cand["date"].duplicated().any():
        cand = cand.drop_duplicates()
        if cand["date"].duplicated().any():
            bad = cand.loc[cand["date"].duplicated(), "date"].iloc[0]
            raise ValueError(f"frame_with_climate: conflicting duplicate climate rows for gage {gage_df['gage_id'].iloc[0]} on {bad.date()}")
    cand = add_water_year_columns(cand, date_col="date")  # water_year, month, dowy
    existing = set(sf_dates.tolist())
    keep = cand["water_year"].isin(wy_set).to_numpy() & ~cand["date"].isin(existing).to_numpy()
    if not keep.any():
        return joined
    add = cand.loc[keep].copy()

    add["gage_id"] = gage_df["gage_id"].iloc[0]
    if "year" in joined.columns:
        add["year"] = add["date"].dt.year
    if "doy" in joined.columns:
        add["doy"] = add["date"].dt.dayofyear
    for c in joined.columns:
        if c not in add.columns:
            if pd.api.types.is_numeric_dtype(joined[c]):
                add[c] = np.nan
            else:
                try:    # keep the column's dtype (e.g. pandas string) with its missing marker
                    add[c] = pd.array([pd.NA] * len(add), dtype=joined[c].dtype)
                except (TypeError, ValueError):
                    add[c] = None
    add = add[list(joined.columns)]
    for c in joined.columns:   # the added rows take the joined frame's dtypes where they can
        if add[c].dtype != joined[c].dtype:
            try:
                add[c] = add[c].astype(joined[c].dtype)
            except (TypeError, ValueError):
                pass
    return pd.concat([joined, add], ignore_index=True)
