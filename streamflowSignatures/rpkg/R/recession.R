#' Identify recession events using position-level marking
#'
#' Marks each position where both recession criteria hold (Q decreasing AND
#' |dQ/dt| decreasing), then finds contiguous runs of valid positions >=
#' min_length. This captures the full extent of each recession event, unlike
#' the previous look-ahead algorithm which could truncate events.
#' @keywords internal
identify_recession_events <- function(Q_vector, min_length = NULL) {
  if (is.null(min_length)) min_length <- pkg_env$recession_min_length
  n <- length(Q_vector)
  if (n < min_length + 1) return(list())

  dQ_dt <- c(diff(Q_vector), NA)

  # Step 1: Mark each position where both recession criteria hold
  is_valid <- logical(n)
  for (i in seq_len(n - 1)) {
    if (!is.na(Q_vector[i]) && !is.na(Q_vector[i + 1]) &&
        !is.na(dQ_dt[i]) && !is.na(dQ_dt[i + 1])) {
      if (Q_vector[i + 1] < Q_vector[i] && abs(dQ_dt[i + 1]) < abs(dQ_dt[i])) {
        is_valid[i] <- TRUE
      }
    }
  }

  # Step 2: Find contiguous runs of valid positions >= min_length
  recession_events <- list()
  run_start <- 0L
  for (i in seq_len(n)) {
    if (is_valid[i]) {
      if (run_start == 0L) run_start <- i
    } else {
      if (run_start > 0L) {
        run_length <- i - run_start
        if (run_length >= min_length) {
          recession_events[[length(recession_events) + 1]] <- list(
            start = run_start, end = i, indices = run_start:i
          )
        }
        run_start <- 0L
      }
    }
  }
  # Handle run ending at data boundary
  if (run_start > 0L) {
    run_length <- n - run_start
    if (run_length >= min_length) {
      recession_events[[length(recession_events) + 1]] <- list(
        start = run_start, end = n, indices = run_start:n
      )
    }
  }

  recession_events
}

#' Fit recession power law for a single event
#' @keywords internal
fit_recession_event <- function(Q_values, remove_first_day = TRUE) {
  if (remove_first_day && length(Q_values) > 1) Q_values <- Q_values[-1]
  n <- length(Q_values)
  if (n < 3) return(list(log_a = NA_real_, b = NA_real_))

  dQ_dt <- -diff(Q_values)
  Q_sub <- Q_values[-length(Q_values)]
  valid <- which(Q_sub > 0 & dQ_dt > 0)
  if (length(valid) < 2) return(list(log_a = NA_real_, b = NA_real_))

  tryCatch({
    fit <- lm(log(dQ_dt[valid]) ~ log(Q_sub[valid]))
    list(log_a = unname(coef(fit)[1]), b = unname(coef(fit)[2]))
  }, error = function(e) list(log_a = NA_real_, b = NA_real_))
}

#' Fit sinusoidal model to recession log(a) vs day-of-year
#' @keywords internal
fit_sinusoidal_model <- function(doy_values, log_a_values) {
  valid <- which(!is.na(log_a_values) & !is.na(doy_values))
  if (length(valid) < 10) return(list(amplitude = NA_real_, minimum_doy = NA_real_))

  doy_clean   <- doy_values[valid]
  log_a_clean <- log_a_values[valid]

  tryCatch({
    X <- cbind(sin(2 * pi * doy_clean / 365),
               cos(2 * pi * doy_clean / 365), 1)
    fit <- lm(log_a_clean ~ X - 1)
    B1 <- coef(fit)[1]; B2 <- coef(fit)[2]

    amplitude <- sqrt(B1^2 + B2^2)
    phase_days <- atan2(-B2, B1) * 365 / (2 * pi)
    if (phase_days < 0) phase_days <- phase_days + 365
    minimum_doy <- phase_days + 273.75
    if (minimum_doy > 365) minimum_doy <- minimum_doy - 365

    list(amplitude = amplitude, minimum_doy = minimum_doy)
  }, error = function(e) list(amplitude = NA_real_, minimum_doy = NA_real_))
}

#' Per-event log(a) under a fixed linear reservoir (b = 1)
#'
#' Median of log(-dQ/dt) - log(Q) over the event's valid pairs — no regression.
#' Port of \code{event_log_a_b1} in julia/src/recession.jl (July 2026
#' convention): fixing b = 1 removes the slope-intercept convolution of the free
#' power-law fit. Conventions match the fitting path: first day removed (storm
#' peak), pairs require Q > 0 and -dQ > 0, natural log.
#'
#' @param Q_event Numeric vector of daily Q for one recession event.
#' @return Median log(a) under b = 1, or NA if there are no valid pairs.
#' @keywords internal
event_log_a_b1 <- function(Q_event) {
  if (length(Q_event) < 3) return(NA_real_)
  Q_work <- Q_event[-1]                 # Remove first day (storm peak)
  dQdt <- -diff(Q_work)
  Q_mid <- Q_work[-length(Q_work)]      # Q at start of each interval
  valid <- !is.na(Q_mid) & !is.na(dQdt) & Q_mid > 0 & dQdt > 0
  if (!any(valid)) return(NA_real_)
  median(log(dQdt[valid]) - log(Q_mid[valid]))
}


#' Analyze recession parameters
#'
#' Computes point cloud and event-based recession parameters (log_a, b),
#' concavity, and seasonality of recession rate. Requires a minimum of 25
#' recession events across the record.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list with 7 annual metrics x 8 stats (log_a_pointcloud,
#'   log_a_events, b_pointcloud, b_events, concavity, n_recession_events,
#'   alpha_linear) + 6 seasonality values + the per-gage
#'   recession_alpha_point_cloud_linear_reservoir scalar.
#' @export
analyze_recession_parameters <- function(streamflow_data,
                                        trend_completeness = NULL,
                                        decade_completeness = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  min_events <- pkg_env$recession_min_events

  signatures_with_stats <- c("log_a_pointcloud", "log_a_events",
                             "b_pointcloud", "b_events", "concavity",
                             "alpha_linear")
  seasonality_sigs <- c("log_a_seasonality_amplitude_all",
                        "log_a_seasonality_minimum_all",
                        "log_a_seasonality_amplitude_first_half",
                        "log_a_seasonality_minimum_first_half",
                        "log_a_seasonality_amplitude_last_half",
                        "log_a_seasonality_minimum_last_half")

  # Build NA result template
  make_na_result <- function() {
    out <- unlist(lapply(signatures_with_stats, empty_stats), recursive = FALSE)
    out <- c(out, unlist(empty_stats("n_recession_events"), recursive = FALSE))
    for (s in seasonality_sigs) out[[s]] <- NA_real_
    out[["recession_alpha_point_cloud_linear_reservoir"]] <- NA_real_
    out
  }

  years <- unique(streamflow_data$water_year)
  annual_metrics <- data.frame(water_year = years, log_a_pointcloud = NA_real_,
                               log_a_events = NA_real_, b_pointcloud = NA_real_,
                               b_events = NA_real_, concavity = NA_real_,
                               n_recession_events = NA_real_,
                               alpha_linear = NA_real_)

  streamflow_data <- streamflow_data[order(streamflow_data$water_year, streamflow_data$dowy), ]

  all_recession_events <- list()
  all_alpha_linear <- numeric(0)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    year_data <- year_data[order(year_data$dowy), ]

    recession_events <- identify_recession_events(year_data$Q)

    # Store event count for this year (INDEPENDENT of min_events gate)
    annual_metrics$n_recession_events[annual_metrics$water_year == yr] <- length(recession_events)

    all_Q <- numeric(0)
    all_dQ_dt <- numeric(0)
    event_log_a_values <- numeric(0)
    event_b_values <- numeric(0)
    event_concavities <- numeric(0)
    year_alpha_linear <- numeric(0)

    for (event in recession_events) {
      Q_event <- year_data$Q[event$indices]
      mid_idx <- event$indices[ceiling(length(event$indices) / 2)]
      event_dowy <- year_data$dowy[mid_idx]

      # Per-event FREE power-law fit — used for b only; fit success also
      # gates event acceptance (unchanged event counts)
      params <- fit_recession_event(Q_event, remove_first_day = TRUE)

      # Per-event alpha under fixed b = 1 (same valid-pair definition as the
      # free fit, so non-NA whenever the fit succeeds)
      log_a_b1 <- event_log_a_b1(Q_event)

      # Compute discrete recession constant alpha = Q_{i+1}/Q_i (b=1 linear reservoir)
      # INDEPENDENT of power-law fit — depends only on raw Q pairs
      # Remove first day (storm peak) consistent with point-cloud fitting
      if (length(Q_event) > 2) {
        Q_alpha <- Q_event[-1]  # Remove first day
        for (j in seq_len(length(Q_alpha) - 1)) {
          if (Q_alpha[j] > 0 && !is.na(Q_alpha[j]) && !is.na(Q_alpha[j + 1])) {
            a_i <- Q_alpha[j + 1] / Q_alpha[j]
            if (a_i > 0 && a_i < 1) {
              year_alpha_linear <- c(year_alpha_linear, a_i)
              all_alpha_linear <- c(all_alpha_linear, a_i)
            }
          }
        }
      }

      if (!is.na(params$log_a) && !is.na(params$b)) {
        event_log_a_values <- c(event_log_a_values, log_a_b1)
        event_b_values <- c(event_b_values, params$b)

        # log_a stored is the b=1 value — the seasonality sinusoid consumes
        # these (July 2026 convention)
        all_recession_events[[length(all_recession_events) + 1]] <- list(
          water_year = yr, dowy = event_dowy,
          log_a = log_a_b1, b = params$b
        )

        # Concavity
        if (length(Q_event) >= 6) {
          mid_point <- floor(length(Q_event) / 2)
          first_half  <- Q_event[1:mid_point]
          second_half <- Q_event[mid_point:length(Q_event)]

          p1 <- fit_recession_event(first_half,  remove_first_day = FALSE)
          p2 <- fit_recession_event(second_half, remove_first_day = FALSE)
          if (!is.na(p1$b) && !is.na(p2$b)) {
            event_concavities <- c(event_concavities, p2$b - p1$b)
          }
        }

        # Point cloud data
        if (length(Q_event) > 1) {
          Q_sub  <- Q_event[-1]
          dQ_sub <- -diff(Q_event[-1])
          valid  <- which(Q_sub[-length(Q_sub)] > 0 & dQ_sub > 0)
          if (length(valid) > 0) {
            all_Q     <- c(all_Q, Q_sub[-length(Q_sub)][valid])
            all_dQ_dt <- c(all_dQ_dt, dQ_sub[valid])
          }
        }
      }
    }

    median_b <- median(event_b_values, na.rm = TRUE)

    # Point cloud
    if (length(all_Q) > 10) {
      log_Q <- log(all_Q)
      log_dQdt <- log(all_dQ_dt)

      # log_a under fixed b = 1 (linear reservoir): log(a) = log(-dQ/dt) -
      # log(Q). No regression is involved, so this needs no singularity gate
      # and is independent of the free b fit below.
      annual_metrics$log_a_pointcloud[annual_metrics$water_year == yr] <-
        median(log_dQdt - log_Q, na.rm = TRUE)

      # b from the FREE point-cloud fit (unchanged)
      tryCatch({
        b_pc <- coef(lm(log_dQdt ~ log_Q))[2]
        if (!is.na(b_pc)) {
          annual_metrics$b_pointcloud[annual_metrics$water_year == yr] <- b_pc
        }
      }, error = function(e) NULL)
    }

    # Events. log_a_events uses the per-event fixed-b=1 values — the previous
    # median-b recalculation is gone: alpha is decoupled from b by assuming a
    # linear reservoir everywhere.
    if (length(event_b_values) > 0) {
      valid_b1 <- event_log_a_values[!is.na(event_log_a_values)]
      if (length(valid_b1) > 0) {
        annual_metrics$log_a_events[annual_metrics$water_year == yr] <- median(valid_b1)
      }
      annual_metrics$b_events[annual_metrics$water_year == yr] <- median_b
    }

    if (length(event_concavities) > 0) {
      annual_metrics$concavity[annual_metrics$water_year == yr] <- mean(event_concavities, na.rm = TRUE)
    }

    # Per-year alpha_linear (linear reservoir assumption)
    if (length(year_alpha_linear) > 10) {
      annual_metrics$alpha_linear[annual_metrics$water_year == yr] <- median(year_alpha_linear, na.rm = TRUE)
    }
  }

  # n_recession_events stats computed INDEPENDENTLY of min_events gate
  events_df <- annual_metrics[!is.na(annual_metrics$n_recession_events),
                              c("water_year", "n_recession_events")]
  if (nrow(events_df) >= 3) {
    n_events_stats <- generate_stats(events_df, value_cols = "n_recession_events",
                                     year_col = "water_year",
                                     trend_completeness = trend_completeness,
                                     decade_completeness = decade_completeness, changepoint = changepoint, collector = collector)
  } else {
    n_events_stats <- empty_stats("n_recession_events")
  }

  # Whole-record recession alpha (scalar — independent of min_events gate)
  recession_alpha_scalar <- if (length(all_alpha_linear) > 10) {
    median(all_alpha_linear, na.rm = TRUE)
  } else {
    NA_real_
  }

  # Check minimum events
  if (length(all_recession_events) < min_events) {
    na_result <- make_na_result()
    # Override n_recession_events with actual computed stats
    for (nm in names(n_events_stats)) na_result[[nm]] <- n_events_stats[[nm]]
    na_result[["recession_alpha_point_cloud_linear_reservoir"]] <- recession_alpha_scalar
    return(na_result)
  }

  result <- generate_stats(annual_metrics, value_cols = signatures_with_stats,
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, changepoint = changepoint, collector = collector)

  # Add n_recession_events stats (computed independently above)
  result <- c(result, n_events_stats)

  # Seasonality
  for (s in seasonality_sigs) result[[s]] <- NA_real_

  if (length(all_recession_events) >= 10) {
    event_dowys   <- vapply(all_recession_events, function(x) x$dowy, numeric(1))
    event_log_as  <- vapply(all_recession_events, function(x) x$log_a, numeric(1))
    event_wys     <- vapply(all_recession_events, function(x) x$water_year, numeric(1))

    seas_all <- fit_sinusoidal_model(event_dowys, event_log_as)
    result$log_a_seasonality_amplitude_all <- seas_all$amplitude
    result$log_a_seasonality_minimum_all   <- seas_all$minimum_doy

    med_wy <- median(unique(event_wys))
    first_idx <- which(event_wys <= med_wy)
    last_idx  <- which(event_wys > med_wy)

    if (length(first_idx) >= 10) {
      s1 <- fit_sinusoidal_model(event_dowys[first_idx], event_log_as[first_idx])
      result$log_a_seasonality_amplitude_first_half <- s1$amplitude
      result$log_a_seasonality_minimum_first_half   <- s1$minimum_doy
    }

    if (length(last_idx) >= 10) {
      s2 <- fit_sinusoidal_model(event_dowys[last_idx], event_log_as[last_idx])
      result$log_a_seasonality_amplitude_last_half <- s2$amplitude
      result$log_a_seasonality_minimum_last_half   <- s2$minimum_doy
    }
  }

  result[["recession_alpha_point_cloud_linear_reservoir"]] <- recession_alpha_scalar

  result
}
