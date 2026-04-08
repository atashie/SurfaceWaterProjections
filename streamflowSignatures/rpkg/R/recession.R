#' Identify recession events using look-ahead algorithm
#' @keywords internal
identify_recession_events <- function(Q_vector, min_length = NULL) {
  if (is.null(min_length)) min_length <- pkg_env$recession_min_length
  n <- length(Q_vector)
  if (n < min_length + 1) return(list())

  dQ_dt <- c(diff(Q_vector), NA)
  recession_events <- list()
  in_recession <- FALSE
  start_idx <- NA_integer_

  for (i in seq_len(n - min_length)) {
    if (is.na(Q_vector[i]) || is.na(dQ_dt[i])) {
      if (in_recession) {
        if (!is.na(start_idx) && (i - start_idx) >= min_length) {
          recession_events[[length(recession_events) + 1]] <- list(
            start = start_idx, end = i - 1L, indices = start_idx:(i - 1L)
          )
        }
        in_recession <- FALSE
        start_idx <- NA_integer_
      }
      next
    }

    if (i < n - 1) {
      is_rec <- TRUE
      for (j in 0:(min_length - 1)) {
        if (i + j + 1 > n || is.na(Q_vector[i + j]) || is.na(Q_vector[i + j + 1]) ||
            is.na(dQ_dt[i + j]) || is.na(dQ_dt[i + j + 1])) {
          is_rec <- FALSE; break
        }
        if (Q_vector[i + j + 1] >= Q_vector[i + j]) { is_rec <- FALSE; break }
        if (abs(dQ_dt[i + j + 1]) >= abs(dQ_dt[i + j])) { is_rec <- FALSE; break }
      }

      if (is_rec && !in_recession) {
        in_recession <- TRUE
        start_idx <- i
      } else if (!is_rec && in_recession) {
        if (!is.na(start_idx) && (i - start_idx) >= min_length) {
          recession_events[[length(recession_events) + 1]] <- list(
            start = start_idx, end = i - 1L, indices = start_idx:(i - 1L)
          )
        }
        in_recession <- FALSE
        start_idx <- NA_integer_
      }
    }
  }

  if (in_recession && !is.na(start_idx) && (n - start_idx) >= min_length) {
    recession_events[[length(recession_events) + 1]] <- list(
      start = start_idx, end = n, indices = start_idx:n
    )
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

#' Analyze recession parameters
#'
#' Computes point cloud and event-based recession parameters (log_a, b),
#' concavity, and seasonality of recession rate. Requires a minimum of 25
#' recession events across the record.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list with 5 metrics x 8 stats + 6 seasonality values.
#' @export
analyze_recession_parameters <- function(streamflow_data,
                                        trend_completeness = NULL,
                                        decade_completeness = NULL) {
  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  min_events <- pkg_env$recession_min_events

  signatures_with_stats <- c("log_a_pointcloud", "log_a_events",
                             "b_pointcloud", "b_events", "concavity")
  seasonality_sigs <- c("log_a_seasonality_amplitude_all",
                        "log_a_seasonality_minimum_all",
                        "log_a_seasonality_amplitude_first_half",
                        "log_a_seasonality_minimum_first_half",
                        "log_a_seasonality_amplitude_last_half",
                        "log_a_seasonality_minimum_last_half")

  # Build NA result template
  make_na_result <- function() {
    out <- unlist(lapply(signatures_with_stats, empty_stats), recursive = FALSE)
    for (s in seasonality_sigs) out[[s]] <- NA_real_
    out
  }

  years <- unique(streamflow_data$water_year)
  annual_metrics <- data.frame(water_year = years, log_a_pointcloud = NA_real_,
                               log_a_events = NA_real_, b_pointcloud = NA_real_,
                               b_events = NA_real_, concavity = NA_real_)

  streamflow_data <- streamflow_data[order(streamflow_data$water_year, streamflow_data$dowy), ]

  all_recession_events <- list()

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    year_data <- year_data[order(year_data$dowy), ]

    recession_events <- identify_recession_events(year_data$Q)

    all_Q <- numeric(0)
    all_dQ_dt <- numeric(0)
    event_log_a_values <- numeric(0)
    event_b_values <- numeric(0)
    event_concavities <- numeric(0)

    for (event in recession_events) {
      Q_event <- year_data$Q[event$indices]
      mid_idx <- event$indices[ceiling(length(event$indices) / 2)]
      event_dowy <- year_data$dowy[mid_idx]

      params <- fit_recession_event(Q_event, remove_first_day = TRUE)

      if (!is.na(params$log_a) && !is.na(params$b)) {
        event_log_a_values <- c(event_log_a_values, params$log_a)
        event_b_values <- c(event_b_values, params$b)

        all_recession_events[[length(all_recession_events) + 1]] <- list(
          water_year = yr, dowy = event_dowy,
          log_a = params$log_a, b = params$b
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
      tryCatch({
        pc_fit <- lm(log(all_dQ_dt) ~ log(all_Q))
        b_pc   <- coef(pc_fit)[2]
        log_a_pc <- median(log(all_dQ_dt) - b_pc * log(all_Q), na.rm = TRUE)
        annual_metrics$b_pointcloud[annual_metrics$water_year == yr]     <- b_pc
        annual_metrics$log_a_pointcloud[annual_metrics$water_year == yr] <- log_a_pc
      }, error = function(e) NULL)
    }

    # Events
    if (length(event_log_a_values) > 0) {
      log_a_recalc <- numeric(0)
      for (ei in seq_along(recession_events)) {
        if (ei <= length(event_b_values) && !is.na(event_b_values[ei])) {
          ev <- recession_events[[ei]]
          Qe <- year_data$Q[ev$indices]
          if (length(Qe) > 1) {
            Qs <- Qe[-1]
            dQs <- -diff(Qe[-1])
            vi <- which(Qs[-length(Qs)] > 0 & dQs > 0)
            if (length(vi) > 0) {
              la <- log(dQs[vi]) - median_b * log(Qs[-length(Qs)][vi])
              log_a_recalc <- c(log_a_recalc, median(la, na.rm = TRUE))
            }
          }
        }
      }
      if (length(log_a_recalc) > 0) {
        annual_metrics$log_a_events[annual_metrics$water_year == yr] <- median(log_a_recalc, na.rm = TRUE)
      }
      annual_metrics$b_events[annual_metrics$water_year == yr] <- median_b
    }

    if (length(event_concavities) > 0) {
      annual_metrics$concavity[annual_metrics$water_year == yr] <- mean(event_concavities, na.rm = TRUE)
    }
  }

  # Check minimum events
  if (length(all_recession_events) < min_events) return(make_na_result())

  result <- generate_stats(annual_metrics, value_cols = signatures_with_stats,
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness)

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

  result
}
