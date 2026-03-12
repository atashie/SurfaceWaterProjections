"""
Average storage signature.

Catchment storage using simplified water balance (ignores ET).
"""

"""
    calculate_average_storage(df::DataFrame; min_years=10) -> Dict

Calculate average catchment storage trends.

Uses simplified water balance (dS = P - Q) without ET estimation.
Calculates storage amount corresponding to average discharge.

Parameters
----------
df : DataFrame
    Daily streamflow and climate data with columns: water_year, Q, PPT, dowy
min_years : Int
    Minimum years required

Returns
-------
Dict{String, Float64}
    Dictionary with avg_storage statistics (8 values)

Notes
-----
This calculation ignores evapotranspiration (ET), using only P - Q
for the water balance. This simplification may overestimate storage
in watersheds with significant ET losses.
"""
function calculate_average_storage(
    df::DataFrame;
    min_years::Int=10,
    min_days_per_year::Int=300,
    max_na_frac::Float64=0.1
)
    result = Dict{String, Float64}()

    # Check for required columns
    if !("PPT" in names(df)) || !("dowy" in names(df))
        return empty_stats("avg_storage")
    end

    if nrow(df) == 0
        return empty_stats("avg_storage")
    end

    years = unique(df.water_year)

    if length(years) < min_years
        return empty_stats("avg_storage")
    end

    annual_storage = Float64[]
    storage_years = Int[]

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_df = copy(df[year_mask, :])

        if nrow(year_df) < min_days_per_year
            continue
        end

        # Convert to clean vectors
        Q_clean = coalesce_q(year_df.Q)
        P_clean = coalesce_q(year_df.PPT)

        # Check NA fraction
        na_frac_Q = sum(isnan.(Q_clean)) / length(Q_clean)
        na_frac_P = sum(isnan.(P_clean)) / length(P_clean)
        if na_frac_Q > max_na_frac || na_frac_P > max_na_frac
            continue
        end

        # Replace NaN with 0 for cumsum
        Q_vals = copy(Q_clean)
        P_vals = copy(P_clean)
        Q_vals[isnan.(Q_vals)] .= 0.0
        P_vals[isnan.(P_vals)] .= 0.0

        # Sort by day of water year (handle Missing dowy)
        dowy_clean = [coalesce(d, 999) for d in year_df.dowy]
        perm = sortperm(dowy_clean)
        Q_vals = Q_vals[perm]
        P_vals = P_vals[perm]

        # Water balance: dS = P - Q
        dS = P_vals .- Q_vals

        # Cumulative storage (relative to start of year)
        S = cumsum(dS)

        # Mean Q for this year
        mean_Q = mean(Q_vals)

        # Sort by Q for interpolation
        sort_idx = sortperm(Q_vals)
        Q_sorted = Q_vals[sort_idx]
        S_sorted = S[sort_idx]

        # Remove NaN (shouldn't be any after replacement, but just in case)
        valid_mask = .!isnan.(Q_sorted) .& .!isnan.(S_sorted)
        Q_valid = Q_sorted[valid_mask]
        S_valid = S_sorted[valid_mask]

        if length(Q_valid) < 10
            continue
        end

        # Average S values at duplicate Q points to match R's approx(ties=mean).
        # Without this, picking a single S at tied Q points introduces noise
        # that corrupts trend statistics (slopes, p-values).
        Q_unique = unique(Q_valid)
        S_means = Float64[mean(S_valid[Q_valid .== q]) for q in Q_unique]

        if length(Q_unique) < 10
            continue
        end

        # Linear interpolation to find S at mean_Q (using de-duplicated arrays)
        idx_below = findlast(Q_unique .<= mean_Q)
        idx_above = findfirst(Q_unique .>= mean_Q)

        if idx_below === nothing || idx_above === nothing
            # mean_Q is outside the range
            if mean_Q < Q_unique[1]
                avg_storage = S_means[1]
            else
                avg_storage = S_means[end]
            end
        elseif idx_below == idx_above
            avg_storage = S_means[idx_below]
        else
            # Linear interpolation
            Q_lo, Q_hi = Q_unique[idx_below], Q_unique[idx_above]
            S_lo, S_hi = S_means[idx_below], S_means[idx_above]

            if abs(Q_hi - Q_lo) < 1e-10
                avg_storage = (S_lo + S_hi) / 2
            else
                t = (mean_Q - Q_lo) / (Q_hi - Q_lo)
                avg_storage = S_lo + t * (S_hi - S_lo)
            end
        end

        if !isnan(avg_storage)
            push!(annual_storage, avg_storage)
            push!(storage_years, yr)
        end
    end

    if length(annual_storage) < 5
        return empty_stats("avg_storage")
    end

    annual_df = DataFrame(
        water_year = storage_years,
        avg_storage = annual_storage
    )

    return generate_stats(annual_df; value_cols=["avg_storage"])
end
