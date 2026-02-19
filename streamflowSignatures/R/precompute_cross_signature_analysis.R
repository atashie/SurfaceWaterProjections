# precompute_cross_signature_analysis.R
#
# Offline script to compute pairwise Theil-Sen slopes between z-scored
# signature means and z-scored signature trends across all gages.
#
# Output: cross_signature_analysis.csv (uploaded to S3)
# Each row = one ordered pair (A -> B) with:
#   - mean_slope:  Theil-Sen slope of z(A_mean) vs z(B_mean)
#   - trend_slope: Theil-Sen slope of z(A_senn_slp) vs z(B_senn_slp)
#
# Usage:
#   Rscript R/precompute_cross_signature_analysis.R
#
# Prerequisites:
#   - AWS credentials configured in .Renviron
#   - Packages: data.table, aws.s3, zyp

library(data.table)
library(aws.s3)

# Source app helper functions (for read_csv_from_s3_direct, fast_theil_sen_slope)
source("streamflowAndClimateVisualizationApp/helperFunctions.R")

# === CONFIGURATION ===
S3_BUCKET <- "climate-ai-data-science-shiny-app-data"
SIGNATURE_CSV_KEY <- "streamflow/streamflow_signatures_full_10feb2026.csv"
OUTPUT_FILE <- "cross_signature_analysis.csv"
OUTPUT_S3_KEY <- "streamflow/cross_signature_analysis.csv"

# === SIGNATURE CATEGORIES ===
# 13 categories mapping base metric names to their grouping
SIGNATURE_CATEGORIES <- list(
  "Flow Volume" = c("Qann", "Qwin", "Qspr", "Qsum", "Qfal"),
  "Flow Percentiles" = c("Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40",
                          "Q50", "Q60", "Q70", "Q75", "Q80", "Q90", "Q95",
                          "Q99", "Q95.Q10"),
  "FDC" = c("FDC90th", "FDCall", "FDCmid"),
  "Baseflow" = c("BFI_Eckhardt", "BFI_LyneHollick"),
  "Recession" = c("log_a_pointcloud", "log_a_events", "b_pointcloud",
                   "b_events", "concavity"),
  "Pulse Metrics" = c("n_high_pulses_year", "n_low_pulses_year",
                       "n_high_pulses_all", "n_low_pulses_all",
                       "dur_high_pulses_year", "dur_low_pulses_year",
                       "dur_high_pulses_all", "dur_low_pulses_all",
                       "TQmean"),
  "Flow Reversals" = c("Flow_Reversals_annual", "Flow_Reversals_winter",
                        "Flow_Reversals_spring", "Flow_Reversals_summer",
                        "Flow_Reversals_fall"),
  "Flashiness" = c("flashinessRB"),
  "Flow Timing" = c("D5_day", "D10_day", "D20_day", "D30_day", "D40_day",
                     "D50_day", "D60_day", "D70_day", "D80_day", "D90_day",
                     "D95_day", "D25_to_D75", "Dmax"),
  "Runoff Ratios" = c("annual_runoff_ratio", "winter_runoff_ratio",
                       "spring_runoff_ratio", "summer_runoff_ratio",
                       "fall_runoff_ratio"),
  "Elasticity" = c("elasticity"),
  "Q-P Seasonality" = c("qp_slope_sd", "qp_bimodality"),
  "Average Storage" = c("avg_storage")
)

# Build reverse lookup: metric -> category
category_lookup <- character(0)
for (cat_name in names(SIGNATURE_CATEGORIES)) {
  for (metric in SIGNATURE_CATEGORIES[[cat_name]]) {
    category_lookup[metric] <- cat_name
  }
}

# === LOAD DATA ===
message("Loading signature data from S3...")
readRenviron("streamflowAndClimateVisualizationApp/.Renviron")
sig_data <- read_csv_from_s3_direct(bucket = S3_BUCKET, object_key = SIGNATURE_CSV_KEY)
message(paste("Loaded", nrow(sig_data), "gages"))

# Filter to area-normalized gages with valid coordinates
if ("area_normalized" %in% names(sig_data)) {
  sig_data <- sig_data[is.na(area_normalized) | area_normalized == TRUE]
}
sig_data <- sig_data[!is.na(latitude) & !is.na(longitude)]
message(paste("After filtering:", nrow(sig_data), "area-normalized gages with coordinates"))

# === IDENTIFY VALID METRICS ===
# A metric is valid if it has both _mean and _senn_slp columns
all_cols <- names(sig_data)
mean_cols <- grep("_mean$", all_cols, value = TRUE)
metric_bases <- gsub("_mean$", "", mean_cols)

# Keep only metrics that have both _mean and _senn_slp
valid_metrics <- character(0)
for (m in metric_bases) {
  mean_col <- paste0(m, "_mean")
  slope_col <- paste0(m, "_senn_slp")
  if (mean_col %in% all_cols && slope_col %in% all_cols) {
    # Also require at least 10 non-NA values in each
    n_mean <- sum(!is.na(sig_data[[mean_col]]))
    n_slope <- sum(!is.na(sig_data[[slope_col]]))
    if (n_mean >= 10 && n_slope >= 10) {
      valid_metrics <- c(valid_metrics, m)
    }
  }
}
message(paste("Found", length(valid_metrics), "metrics with both _mean and _senn_slp columns"))

# === Z-SCORE EACH METRIC COLUMN ===
message("Z-scoring metric columns...")
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# Create z-scored columns
for (m in valid_metrics) {
  mean_col <- paste0(m, "_mean")
  slope_col <- paste0(m, "_senn_slp")
  set(sig_data, j = paste0("z_mean_", m), value = zscore(sig_data[[mean_col]]))
  set(sig_data, j = paste0("z_slope_", m), value = zscore(sig_data[[slope_col]]))
}

# === COMPUTE PAIRWISE SLOPES ===
n_metrics <- length(valid_metrics)
n_pairs <- n_metrics * (n_metrics - 1)
message(paste("Computing", n_pairs, "pairwise Theil-Sen slopes..."))

results <- vector("list", n_pairs)
idx <- 0
t_start <- Sys.time()

for (i in seq_along(valid_metrics)) {
  from_metric <- valid_metrics[i]
  z_mean_from <- sig_data[[paste0("z_mean_", from_metric)]]
  z_slope_from <- sig_data[[paste0("z_slope_", from_metric)]]

  for (j in seq_along(valid_metrics)) {
    if (i == j) next
    to_metric <- valid_metrics[j]
    z_mean_to <- sig_data[[paste0("z_mean_", to_metric)]]
    z_slope_to <- sig_data[[paste0("z_slope_", to_metric)]]

    # Theil-Sen slope: z(from) vs z(to) for means
    mean_slope <- fast_theil_sen_slope(z_mean_from, z_mean_to)
    n_gages_mean <- sum(complete.cases(z_mean_from, z_mean_to))

    # Theil-Sen slope: z(from) vs z(to) for trends
    trend_slope <- fast_theil_sen_slope(z_slope_from, z_slope_to)
    n_gages_trend <- sum(complete.cases(z_slope_from, z_slope_to))

    idx <- idx + 1
    results[[idx]] <- data.table(
      from_metric = from_metric,
      to_metric = to_metric,
      mean_slope = mean_slope,
      trend_slope = trend_slope,
      from_category = ifelse(from_metric %in% names(category_lookup),
                             category_lookup[from_metric], "Other"),
      to_category = ifelse(to_metric %in% names(category_lookup),
                           category_lookup[to_metric], "Other"),
      n_gages_mean = n_gages_mean,
      n_gages_trend = n_gages_trend
    )
  }

  # Progress
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  pct <- round(100 * i / n_metrics, 1)
  rate <- if (elapsed > 0) round(idx / elapsed, 1) else 0
  message(sprintf("  [%s%%] %d/%d metrics done (%d pairs, %.1f pairs/sec)",
                  pct, i, n_metrics, idx, rate))
}

result_dt <- rbindlist(results)
message(paste("Computed", nrow(result_dt), "pairwise slopes"))

# === WRITE AND UPLOAD ===
output_path <- file.path("streamflowAndClimateVisualizationApp", OUTPUT_FILE)
fwrite(result_dt, output_path)
message(paste("Wrote", output_path))

# Upload to S3
message(paste("Uploading to S3:", OUTPUT_S3_KEY))
put_object(
  file = output_path,
  object = OUTPUT_S3_KEY,
  bucket = S3_BUCKET,
  region = Sys.getenv("AWS_DEFAULT_REGION")
)
message("Upload complete!")

# === SUMMARY ===
message("\n=== SUMMARY ===")
message(paste("Gages used:", nrow(sig_data)))
message(paste("Metrics:", length(valid_metrics)))
message(paste("Pairs:", nrow(result_dt)))
message(paste("Mean slope range:", round(min(result_dt$mean_slope, na.rm = TRUE), 4),
              "to", round(max(result_dt$mean_slope, na.rm = TRUE), 4)))
message(paste("Trend slope range:", round(min(result_dt$trend_slope, na.rm = TRUE), 4),
              "to", round(max(result_dt$trend_slope, na.rm = TRUE), 4)))
message(paste("NA mean slopes:", sum(is.na(result_dt$mean_slope)),
              "| NA trend slopes:", sum(is.na(result_dt$trend_slope))))

cat_counts <- result_dt[, .N, by = from_category][order(-N)]
message("\nPairs by from_category:")
for (r in seq_len(nrow(cat_counts))) {
  message(paste(" ", cat_counts$from_category[r], ":", cat_counts$N[r]))
}
