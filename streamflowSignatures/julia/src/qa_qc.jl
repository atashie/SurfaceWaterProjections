"""
QA/QC flag computation for streamflow signatures.

Computes quality control flags based on expected ranges and consistency checks.
"""


"""
    compute_qa_flags(df::DataFrame) -> DataFrame

Compute QA/QC flags for a DataFrame of streamflow signatures.

Returns
-------
DataFrame
    Input DataFrame with 12 QA/QC flag columns added:
    - flagged_for_qann_range: Qann_mean outside 0-2000 mm
    - flagged_for_bfi_eckhardt_range: BFI_Eckhardt outside 0-1
    - flagged_for_bfi_lynehollick_range: BFI_LyneHollick outside 0-1
    - flagged_for_flashiness_range: flashinessRB outside 0-2
    - flagged_for_tqmean_range: TQmean outside 0-100%
    - flagged_for_d50_range: D50_day outside 1-366
    - flagged_for_elasticity_range: elasticity_static outside 0.1-5
    - flagged_for_runoff_ratio_range: annual_runoff_ratio outside 0.01-1.5
    - flagged_for_seasonal_sum: |seasonal_sum/annual - 1| > tolerance
    - flagged_for_percentile_order: Q5 < Q25 < Q50 < Q75 < Q95 violated
    - flagged_for_timing_order: D5 < D50 < D95 violated
    - flagged_for_high_na: NA/NaN fraction of the SIGNATURE columns > 30%
      (see high_na_denominator_columns)

A flag value of true indicates a potential quality issue.
"""
function compute_qa_flags(df::DataFrame)
    result = copy(df)
    n = nrow(result)

    # Initialize all flag columns to false
    flag_columns = [
        "flagged_for_qann_range",
        "flagged_for_bfi_eckhardt_range",
        "flagged_for_bfi_lynehollick_range",
        "flagged_for_flashiness_range",
        "flagged_for_tqmean_range",
        "flagged_for_d50_range",
        "flagged_for_elasticity_range",
        "flagged_for_runoff_ratio_range",
        "flagged_for_seasonal_sum",
        "flagged_for_percentile_order",
        "flagged_for_timing_order",
        "flagged_for_high_na",
    ]

    for col in flag_columns
        result[!, col] = falses(n)
    end

    # 1. Qann_mean outside 0-2000 mm
    if "Qann_mean" in names(result)
        qann = result.Qann_mean
        result.flagged_for_qann_range = coalesce.(
            (qann .< CFG_QAQC_QANN_RANGE[1]) .| (qann .> CFG_QAQC_QANN_RANGE[2]),
            false
        )
    end

    # 2. BFI_Eckhardt outside 0-1
    if "BFI_Eckhardt_mean" in names(result)
        bfi = result.BFI_Eckhardt_mean
        result.flagged_for_bfi_eckhardt_range = coalesce.(
            (bfi .< CFG_QAQC_BFI_RANGE[1]) .| (bfi .> CFG_QAQC_BFI_RANGE[2]),
            false
        )
    end

    # 3. BFI_LyneHollick outside 0-1
    if "BFI_LyneHollick_mean" in names(result)
        bfi = result.BFI_LyneHollick_mean
        result.flagged_for_bfi_lynehollick_range = coalesce.(
            (bfi .< CFG_QAQC_BFI_RANGE[1]) .| (bfi .> CFG_QAQC_BFI_RANGE[2]),
            false
        )
    end

    # 4. flashinessRB outside 0-2
    if "flashinessRB_mean" in names(result)
        flash = result.flashinessRB_mean
        result.flagged_for_flashiness_range = coalesce.(
            (flash .< CFG_QAQC_FLASHINESS_RANGE[1]) .| (flash .> CFG_QAQC_FLASHINESS_RANGE[2]),
            false
        )
    end

    # 5. TQmean outside 0-100%
    if "TQmean_mean" in names(result)
        tq = result.TQmean_mean
        result.flagged_for_tqmean_range = coalesce.(
            (tq .< CFG_QAQC_TQMEAN_RANGE[1]) .| (tq .> CFG_QAQC_TQMEAN_RANGE[2]),
            false
        )
    end

    # 6. D50_day outside 1-366
    if "D50_day_mean" in names(result)
        d50 = result.D50_day_mean
        result.flagged_for_d50_range = coalesce.(
            (d50 .< CFG_QAQC_D50_RANGE[1]) .| (d50 .> CFG_QAQC_D50_RANGE[2]),
            false
        )
    end

    # 7. elasticity_static outside 0.1-5
    if "elasticity_static" in names(result)
        elast = result.elasticity_static
        result.flagged_for_elasticity_range = coalesce.(
            (elast .< CFG_QAQC_ELASTICITY_RANGE[1]) .| (elast .> CFG_QAQC_ELASTICITY_RANGE[2]),
            false
        )
    end

    # 8. annual_runoff_ratio outside 0.01-1.5
    if "annual_runoff_ratio_mean" in names(result)
        rr = result.annual_runoff_ratio_mean
        result.flagged_for_runoff_ratio_range = coalesce.(
            (rr .< CFG_QAQC_RUNOFF_RATIO_RANGE[1]) .| (rr .> CFG_QAQC_RUNOFF_RATIO_RANGE[2]),
            false
        )
    end

    # 9. Seasonal sum deviation > tolerance
    seasonal_cols = ["Qwin_mean", "Qspr_mean", "Qsum_mean", "Qfal_mean"]
    if all(col -> col in names(result), seasonal_cols) && "Qann_mean" in names(result)
        seasonal_sum = result.Qwin_mean .+ result.Qspr_mean .+ result.Qsum_mean .+ result.Qfal_mean
        qann = result.Qann_mean

        ratio = map(eachindex(seasonal_sum)) do i
            if ismissing(seasonal_sum[i]) || ismissing(qann[i]) || qann[i] == 0
                return missing
            end
            return seasonal_sum[i] / qann[i]
        end

        deviation = abs.(coalesce.(ratio, 1.0) .- 1.0)
        result.flagged_for_seasonal_sum = coalesce.(deviation .> CFG_QAQC_SEASONAL_SUM_TOLERANCE, false)
    end

    # 10. Percentile order violation (Q5 < Q25 < Q50 < Q75 < Q95)
    pct_cols = ["Q5_mean", "Q25_mean", "Q50_mean", "Q75_mean", "Q95_mean"]
    if all(col -> col in names(result), pct_cols)
        q5 = result.Q5_mean
        q25 = result.Q25_mean
        q50 = result.Q50_mean
        q75 = result.Q75_mean
        q95 = result.Q95_mean

        result.flagged_for_percentile_order = coalesce.(
            (q5 .> q25) .| (q25 .> q50) .| (q50 .> q75) .| (q75 .> q95),
            false
        )
    end

    # 11. Timing order violation (D5 < D50 < D95)
    timing_cols = ["D5_day_mean", "D50_day_mean", "D95_day_mean"]
    if all(col -> col in names(result), timing_cols)
        d5 = result.D5_day_mean
        d50 = result.D50_day_mean
        d95 = result.D95_day_mean

        result.flagged_for_timing_order = coalesce.(
            (d5 .> d50) .| (d50 .> d95),
            false
        )
    end

    # 12. High NA fraction: > max_na_fraction of the gage's SIGNATURE columns are NA.
    # Denominator = the signature columns PRESENT in the table (16 statistic suffixes +
    # the per-gage scalar registry, config qa_qc.high_na_denominator — one manifest
    # shared by all three languages). Metadata, unregistered diagnostics and flag
    # columns never enter it; families NA-filled by the table union DO count (that
    # structural NA is what the flag is meant to surface). The NA test is value-level
    # so it works on untyped Vector{Any} columns: the benchmark runner builds them that
    # way, and the previous per-column eltype filter silently excluded every signature
    # column, leaving only the 16 numeric metadata columns (found 2026-09-04).
    signature_cols = high_na_denominator_columns(names(result))
    if !isempty(signature_cols)
        na_fraction = map(1:n) do i
            na_count = 0
            for col in signature_cols
                val = result[i, col]
                if ismissing(val) || (val isa Number && isnan(val))
                    na_count += 1
                end
            end
            return na_count / length(signature_cols)
        end
        result.flagged_for_high_na = na_fraction .> CFG_QAQC_MAX_NA_FRACTION
    end

    return result
end


"""
    high_na_denominator_columns(colnames) -> Vector{String}

The columns that enter the `flagged_for_high_na` denominator: every name ending in
one of the statistic suffixes (`CFG_QAQC_HIGH_NA_SUFFIXES`) plus the per-gage scalar
registry (`CFG_QAQC_HIGH_NA_SCALARS`), in input order; `flagged_*` columns excluded.
Shared definition with `python/streamflow_signatures/qa_qc.py` and `rpkg/R/qa_qc.R`.
"""
function high_na_denominator_columns(colnames)
    scalars = Set(CFG_QAQC_HIGH_NA_SCALARS)
    out = String[]
    for c in colnames
        name = String(c)
        startswith(name, "flagged_") && continue
        if name in scalars || any(sfx -> endswith(name, sfx), CFG_QAQC_HIGH_NA_SUFFIXES)
            push!(out, name)
        end
    end
    return out
end


"""
    get_flag_columns() -> Vector{String}

Return list of QA/QC flag column names.
"""
function get_flag_columns()
    return [
        "flagged_for_qann_range",
        "flagged_for_bfi_eckhardt_range",
        "flagged_for_bfi_lynehollick_range",
        "flagged_for_flashiness_range",
        "flagged_for_tqmean_range",
        "flagged_for_d50_range",
        "flagged_for_elasticity_range",
        "flagged_for_runoff_ratio_range",
        "flagged_for_seasonal_sum",
        "flagged_for_percentile_order",
        "flagged_for_timing_order",
        "flagged_for_high_na",
    ]
end
