#' Eckhardt recursive digital filter for baseflow separation
#'
#' @param Q Numeric vector of daily streamflow.
#' @param BFImax Maximum baseflow index (default from config).
#' @param a Recession constant (default from config).
#' @return Numeric vector of baseflow values.
#' @keywords internal
eckhardt_filter <- function(Q, BFImax = pkg_env$eckhardt_bfimax,
                            a = pkg_env$eckhardt_alpha) {
  n <- length(Q)
  baseflow <- numeric(n)

  if (!is.na(Q[1]) && Q[1] > 0) {
    baseflow[1] <- min(Q[1] * BFImax, Q[1])
  } else {
    baseflow[1] <- 0
  }

  for (i in 2:n) {
    if (is.na(Q[i])) {
      baseflow[i] <- baseflow[i - 1]
      next
    }
    numerator   <- (1 - BFImax) * a * baseflow[i - 1] + (1 - a) * BFImax * Q[i]
    denominator <- 1 - a * BFImax
    baseflow[i] <- numerator / denominator
    baseflow[i] <- min(baseflow[i], Q[i])
    baseflow[i] <- max(baseflow[i], 0)
  }

  baseflow
}

#' Lyne-Hollick recursive digital filter for baseflow separation
#'
#' @param Q Numeric vector of daily streamflow.
#' @param alpha Filter parameter (default from config).
#' @param passes Number of forward-backward passes (default from config).
#' @return Numeric vector of baseflow values.
#' @keywords internal
lyne_hollick_filter <- function(Q, alpha = pkg_env$lyne_hollick_alpha,
                                passes = pkg_env$lyne_hollick_passes) {
  n <- length(Q)
  input_signal <- Q

  for (pass in seq_len(passes)) {
    qf_forward <- numeric(n)
    qf_forward[1] <- 0

    for (i in 2:n) {
      if (!is.na(input_signal[i]) && !is.na(input_signal[i - 1]) && !is.na(qf_forward[i - 1])) {
        qf_forward[i] <- alpha * qf_forward[i - 1] +
          ((1 + alpha) / 2) * (input_signal[i] - input_signal[i - 1])
        qf_forward[i] <- max(0, min(qf_forward[i], input_signal[i]))
      } else {
        qf_forward[i] <- NA
      }
    }

    qf_backward <- numeric(n)
    qf_backward[n] <- qf_forward[n]

    for (i in (n - 1):1) {
      if (!is.na(qf_forward[i]) && !is.na(qf_forward[i + 1]) && !is.na(qf_backward[i + 1])) {
        qf_backward[i] <- alpha * qf_backward[i + 1] +
          ((1 + alpha) / 2) * (qf_forward[i] - qf_forward[i + 1])
        qf_backward[i] <- max(0, min(qf_backward[i], input_signal[i]))
      } else {
        qf_backward[i] <- NA
      }
    }

    input_signal <- input_signal - qf_backward
    input_signal[input_signal < 0] <- 0
    input_signal[is.na(Q)] <- NA
  }

  input_signal
}

#' Analyze baseflow indices (Eckhardt and Lyne-Hollick)
#'
#' Calculates BFI_Eckhardt and BFI_LyneHollick for each water year, then
#' applies generate_stats() to produce trend statistics.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics per BFI metric (2 x 8 = 16 values).
#' @export
analyze_baseflow_indices <- function(streamflow_data,
                                    trend_completeness = NULL,
                                    decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  # Zero-row input is legitimate and canonical Julia guards it. Without this,
  # data.frame(water_year = <empty>, X = NA_real_) raises "arguments imply
  # differing number of rows", which safe_call() swallows — silently dropping
  # the whole family for that gage (the class of the snow defect, 2026-08-26).
  # Route through generate_stats so the key set is identical to the normal path.
  if (nrow(streamflow_data) == 0L || length(years) == 0L) {
    .empty <- data.frame(water_year = integer(0))
    for (.m in c("BFI_Eckhardt", "BFI_LyneHollick")) .empty[[.m]] <- numeric(0)
    return(generate_stats(.empty, value_cols = c("BFI_Eckhardt", "BFI_LyneHollick"),
                          year_col = "water_year",
                          trend_completeness = trend_completeness,
                          decade_completeness = decade_completeness,
                          min_values_for_stats = min_values_for_stats,
                          changepoint = changepoint, collector = collector))
  }

  bfi_by_year <- data.frame(water_year = years, BFI_Eckhardt = NA_real_,
                            BFI_LyneHollick = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    year_data <- year_data[order(year_data$dowy), ]
    Q <- year_data$Q

    bf_eck  <- eckhardt_filter(Q)
    bf_lyne <- lyne_hollick_filter(Q)

    # Eckhardt: paired masking — only sum where both Q and baseflow are valid
    eck_valid <- !is.na(Q) & !is.na(bf_eck)
    total_eck <- sum(Q[eck_valid])
    if (total_eck > 0) {
      bfi_eck <- sum(bf_eck[eck_valid]) / total_eck
      bfi_eck <- max(0, min(1, bfi_eck))
      bfi_by_year$BFI_Eckhardt[bfi_by_year$water_year == yr] <- bfi_eck
    }

    # Lyne-Hollick: paired masking
    lh_valid <- !is.na(Q) & !is.na(bf_lyne)
    total_lh <- sum(Q[lh_valid])
    if (total_lh > 0) {
      bfi_lh <- sum(bf_lyne[lh_valid]) / total_lh
      bfi_lh <- max(0, min(1, bfi_lh))
      bfi_by_year$BFI_LyneHollick[bfi_by_year$water_year == yr] <- bfi_lh
    }
  }

  # Julia APPENDS only computable years, so a non-computable year is ABSENT from
  # its frame. rpkg pre-allocates a row per year and assigns into it, so a skipped
  # year left an NA row behind — which the annual collector then exported as a row
  # Julia never emits (1,039 such rows across the sparse families, caught
  # 2026-08-26 by the cross-language annual gate). Drop them to match.
  bfi_by_year <- bfi_by_year[!is.na(bfi_by_year$BFI_Eckhardt) | !is.na(bfi_by_year$BFI_LyneHollick), , drop = FALSE]

  generate_stats(bfi_by_year, value_cols = c("BFI_Eckhardt", "BFI_LyneHollick"),
                 year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

#' Analyze baseflow indices with recession-derived parameters
#'
#' Same as \code{analyze_baseflow_indices()} but uses recession-derived alpha
#' instead of fixed defaults. BFImax remains fixed at 0.8.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @param alpha Recession-derived filter parameter (from
#'   recession_alpha_point_cloud_linear_reservoir). Must be in (0, 1).
#' @param BFImax Maximum baseflow index for Eckhardt filter (default from config).
#' @param passes Number of Lyne-Hollick filter passes (default from config).
#' @param trend_completeness Forwarded to generate_stats().
#' @param decade_completeness Forwarded to generate_stats().
#' @return Named list of 8 statistics per BFI metric (2 x 8 = 16 values).
#' @export
analyze_baseflow_indices_with_parameters <- function(
    streamflow_data,
    alpha,
    BFImax = pkg_env$eckhardt_bfimax,
    passes = pkg_env$lyne_hollick_passes,
    trend_completeness = NULL,
    decade_completeness = NULL,
    min_values_for_stats = NULL, changepoint = NULL, collector = NULL
) {
  metrics <- c("BFI_Eckhardt_param", "BFI_LyneHollick_param")

  # Validate alpha
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
    result <- unlist(lapply(metrics, empty_stats), recursive = FALSE)
    return(result)
  }

  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  # Zero-row input is legitimate and canonical Julia guards it. Without this,
  # data.frame(water_year = <empty>, X = NA_real_) raises "arguments imply
  # differing number of rows", which safe_call() swallows — silently dropping
  # the whole family for that gage (the class of the snow defect, 2026-08-26).
  # Route through generate_stats so the key set is identical to the normal path.
  if (nrow(streamflow_data) == 0L || length(years) == 0L) {
    .empty <- data.frame(water_year = integer(0))
    for (.m in c("BFI_Eckhardt_param", "BFI_LyneHollick_param")) .empty[[.m]] <- numeric(0)
    return(generate_stats(.empty, value_cols = c("BFI_Eckhardt_param", "BFI_LyneHollick_param"),
                          year_col = "water_year",
                          trend_completeness = trend_completeness,
                          decade_completeness = decade_completeness,
                          min_values_for_stats = min_values_for_stats,
                          changepoint = changepoint, collector = collector))
  }

  bfi_by_year <- data.frame(water_year = years, BFI_Eckhardt_param = NA_real_,
                            BFI_LyneHollick_param = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    year_data <- year_data[order(year_data$dowy), ]
    Q <- year_data$Q

    # Apply filters with recession-derived alpha
    bf_eck <- eckhardt_filter(Q, BFImax = BFImax, a = alpha)
    bf_lh <- lyne_hollick_filter(Q, alpha = alpha, passes = passes)

    # Eckhardt: paired masking
    eck_valid <- !is.na(Q) & !is.na(bf_eck)
    total_eck <- sum(Q[eck_valid])
    if (total_eck > 0) {
      bfi_eck <- sum(bf_eck[eck_valid]) / total_eck
      bfi_eck <- max(0, min(1, bfi_eck))
      bfi_by_year$BFI_Eckhardt_param[bfi_by_year$water_year == yr] <- bfi_eck
    }

    # Lyne-Hollick: paired masking
    lh_valid <- !is.na(Q) & !is.na(bf_lh)
    total_lh <- sum(Q[lh_valid])
    if (total_lh > 0) {
      bfi_lh <- sum(bf_lh[lh_valid]) / total_lh
      bfi_lh <- max(0, min(1, bfi_lh))
      bfi_by_year$BFI_LyneHollick_param[bfi_by_year$water_year == yr] <- bfi_lh
    }
  }

  # Julia APPENDS only computable years, so a non-computable year is ABSENT from
  # its frame. rpkg pre-allocates a row per year and assigns into it, so a skipped
  # year left an NA row behind — which the annual collector then exported as a row
  # Julia never emits (1,039 such rows across the sparse families, caught
  # 2026-08-26 by the cross-language annual gate). Drop them to match.
  bfi_by_year <- bfi_by_year[rowSums(!is.na(bfi_by_year[, metrics, drop = FALSE])) > 0, , drop = FALSE]

  generate_stats(bfi_by_year, value_cols = metrics, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}
