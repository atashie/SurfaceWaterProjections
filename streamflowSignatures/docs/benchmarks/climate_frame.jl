# climate_frame.jl — the production climate-frame helper of the Julia runner
# (Phase A′ of docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.2; decisions
# 12–13, user 2026-10-26). Mirrors: climate_frame.py (Python runner), climate_frame.R
# (rpkg runner). Self-test: selftest_climate_frame.jl.
#
# The runner joins each gage's Daymet rows onto that gage's streamflow ROWS before the
# library builds the water-year grid. A date with no streamflow row therefore reached
# `preprocess_daily_data` without PPT / TMIN / TMAX / SWE, the grid made it NaN and the
# ≤ 3-day interpolation replaced climate that exists in the input (CHANGELOG Known Issue
# "climate rides on the streamflow dates", 2026-10-16). `frame_with_climate` builds the
# frame under a runtime switch, STREAMFLOW_CLIMATE_FRAME:
#   left   today's left join, bit for bit — the attribution baseline (stages S0, S1);
#   outer  the left join PLUS the gage's Daymet rows on dates that have no streamflow row,
#          restricted to dates inside the streamflow record span AND to water years that
#          hold at least one streamflow row (Codex round 3: the February records hold
#          whole interior water years with no streamflow at all; a span-only join would
#          put them into the preprocessor's year loop and move the Q-only
#          `season_excluded_years_*` scalars and `total_possible`).
# Both frames must carry a `Date` column `date` (asserted). Added rows carry gage_id, date, the temporal columns (water_year, month, dowy from
# add_water_year_columns; year and doy from the date when the frame has them), the
# climate columns, and `missing` for Q, flag and every other column; they are appended
# after the joined rows (the preprocessor re-grids every water year by date, so the row
# order is immaterial). The library is untouched (decision 12). The default mode is
# `left` until the qualifier rerun's products land (then `outer`, Phase H).

using DataFrames, Dates

const CLIMATE_FRAME_MODES = ("left", "outer")

"""
    climate_frame_mode() -> String

Read `STREAMFLOW_CLIMATE_FRAME` at runtime; "left" when unset; any other value errors.
"""
function climate_frame_mode()
    v = lowercase(strip(get(ENV, "STREAMFLOW_CLIMATE_FRAME", "left")))
    v in CLIMATE_FRAME_MODES ||
        error("STREAMFLOW_CLIMATE_FRAME must be one of $(join(CLIMATE_FRAME_MODES, ", ")); got \"$v\"")
    return v
end

"""
    frame_with_climate(gage_df, gage_climate, mode) -> DataFrame

`gage_df`: one gage's streamflow rows (gage_id, date, Q, water_year, ... as read from the
parquet); `gage_climate`: that gage's climate rows (gage_id, date, PPT[, SWE, TMIN, TMAX]);
`mode`: "left" or "outer" (see the file header).
"""
function frame_with_climate(gage_df::DataFrame, gage_climate::DataFrame, mode::AbstractString)
    # Both date columns must be Date: a DateTime streamflow column would match no climate
    # row in the join and silently yield a climate-free frame (review 2026-10-26)
    for (nm, df) in (("streamflow", gage_df), ("climate", gage_climate))
        "date" in names(df) || error("frame_with_climate: the $nm frame has no date column")
        nonmissingtype(eltype(df.date)) <: Date ||
            error("frame_with_climate: the $nm frame's date column is $(eltype(df.date)), expected Date")
    end
    if mode == "left"
        return leftjoin(gage_df, gage_climate, on=[:gage_id, :date])
    end
    mode == "outer" || error("frame_with_climate: unknown mode \"$mode\"")
    joined = leftjoin(gage_df, gage_climate, on=[:gage_id, :date])
    (nrow(gage_df) == 0 || nrow(gage_climate) == 0) && return joined
    "water_year" in names(gage_df) || error("frame_with_climate: the streamflow frame has no water_year column")

    sf_dates = Date.(gage_df.date)
    lo, hi = extrema(sf_dates)
    wy_set = Set(skipmissing(gage_df.water_year))
    c_dates = Date.(gage_climate.date)
    in_span = (c_dates .>= lo) .& (c_dates .<= hi)
    any(in_span) || return joined

    cand = gage_climate[in_span, :]
    # duplicate climate dates: identical duplicate rows collapse to one, conflicting ones are
    # an error — an appended duplicate would put two rows on one grid day (review 2026-10-26)
    if length(unique(cand.date)) != nrow(cand)
        cand = unique(cand)
        if length(unique(cand.date)) != nrow(cand)
            error("frame_with_climate: conflicting duplicate climate rows for gage $(gage_df.gage_id[1]) on $(first(cand.date[nonunique(cand, :date)]))")
        end
    end
    cand = add_water_year_columns(cand)          # water_year, month, dowy from the date
    existing = Set(sf_dates)
    keep = [(cand.water_year[i] in wy_set) && !(cand.date[i] in existing) for i in 1:nrow(cand)]
    any(keep) || return joined
    add = cand[keep, :]

    add.gage_id = fill(gage_df.gage_id[1], nrow(add))
    "year" in names(joined) && (add[!, :year] = year.(add.date))
    "doy" in names(joined) && (add[!, :doy] = dayofyear.(add.date))
    for c in names(joined)
        c in names(add) || (add[!, c] = fill(missing, nrow(add)))
    end
    return vcat(joined, add[:, names(joined)])
end
