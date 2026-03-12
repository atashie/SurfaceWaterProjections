"""
Q-P seasonality signatures.

Seasonality in cumulative streamflow vs precipitation relationship.
"""

# Configuration loaded from config.jl (CFG_QP_SLOPE_WINDOW_DAYS)


"""
    calculate_qp_seasonality(df::DataFrame; slope_window=30, min_years=10) -> Dict

Calculate Q-P seasonality metrics.

Quantifies how the relationship between cumulative streamflow and
cumulative precipitation varies seasonally.

Calculates:
- qp_slope_sd: Standard deviation of monthly Q-P slopes
- qp_bimodality: Bimodality coefficient (>0.555 suggests seasonal patterns)

Parameters
----------
df : DataFrame
    Daily streamflow and climate data with columns: water_year, Q, PPT, month, dowy
slope_window : Int
    Rolling window for slope calculation (days)
min_years : Int
    Minimum years required

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics
"""
function calculate_qp_seasonality(
    df::DataFrame;
    slope_window::Int=CFG_QP_SLOPE_WINDOW_DAYS,
    min_years::Int=10,
    min_days_per_year::Int=300,
    max_na_frac::Float64=0.1
)
    result = Dict{String, Float64}()
    metrics = ["qp_slope_sd", "qp_bimodality"]

    # Check for required columns
    if !("PPT" in names(df)) || !("month" in names(df))
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    if nrow(df) == 0
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    years = unique(df.water_year)

    if length(years) < min_years
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    annual_data = DataFrame()

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
        months = Int[coalesce(year_df.month[i], 0) for i in perm]

        n = length(Q_vals)
        if n < slope_window
            continue
        end

        # Cumulative sums
        cum_Q = cumsum(Q_vals)
        cum_P = cumsum(P_vals)

        # Calculate rolling slopes
        slopes = Float64[]
        slope_months = Int[]

        for end_idx in slope_window:n
            start_idx = end_idx - slope_window + 1

            delta_P = cum_P[end_idx] - cum_P[start_idx]
            delta_Q = cum_Q[end_idx] - cum_Q[start_idx]

            if abs(delta_P) >= 0.01
                push!(slopes, delta_Q / delta_P)
            else
                push!(slopes, NaN)
            end

            # Associate with middle of window (match R: end_idx - floor(window/2))
            mid_idx = end_idx - div(slope_window, 2)
            push!(slope_months, months[mid_idx])
        end

        # Calculate monthly mean slopes
        monthly_means = Dict{Int, Float64}()
        for m in 1:12
            month_slopes = [slopes[i] for i in 1:length(slopes)
                           if slope_months[i] == m && !isnan(slopes[i])]
            if !isempty(month_slopes)
                monthly_means[m] = mean(month_slopes)
            end
        end

        if length(monthly_means) < 6
            continue
        end

        # Metric 1: SD of monthly slopes
        qp_slope_sd = std(collect(values(monthly_means)))

        # Metric 2: Bimodality coefficient
        # B = (skewness^2 + 1) / kurtosis
        slopes_clean = filter(!isnan, slopes)
        if length(slopes_clean) >= 30
            m = mean(slopes_clean)
            s = std(slopes_clean)

            if s > 1e-10
                n_s = length(slopes_clean)
                skewness = sum((slopes_clean .- m).^3) / (n_s * s^3)
                kurtosis = sum((slopes_clean .- m).^4) / (n_s * s^4)

                if kurtosis > 1e-10
                    qp_bimodality = (skewness^2 + 1) / kurtosis
                else
                    qp_bimodality = NaN
                end
            else
                qp_bimodality = NaN
            end
        else
            qp_bimodality = NaN
        end

        push!(annual_data, (
            water_year = yr,
            qp_slope_sd = qp_slope_sd,
            qp_bimodality = qp_bimodality
        ))
    end

    if nrow(annual_data) < 5
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Generate statistics
    sd_stats = generate_stats(annual_data; value_cols=["qp_slope_sd"])
    merge!(result, sd_stats)

    # Bimodality - filter NaN
    bi_df = annual_data[.!isnan.(annual_data.qp_bimodality), :]
    if nrow(bi_df) >= 3
        bi_stats = generate_stats(bi_df; value_cols=["qp_bimodality"])
        merge!(result, bi_stats)
    else
        merge!(result, empty_stats("qp_bimodality"))
    end

    return result
end
