################################################################################
# Comparison Script: Restricted (1993-2022) vs Baseline Analysis
#
# Compares signature outputs between:
#   - Baseline: D:/processedOuts_feb2026/streamflow_signatures_full_10feb2026.csv
#   - Restricted: D:/processedOuts_feb2026/streamflow_signatures_full_1993-to-2022-min-20yrs.csv
#
# Outputs:
#   - Text summary: D:/processedOuts_feb2026/comparison_summary_1993-to-2022.txt
#   - Visualization: D:/processedOuts_feb2026/comparison_plots_1993-to-2022.pdf
################################################################################

# Clear environment
rm(list = ls())

# Load required libraries
library(data.table)
library(ggplot2)

cat("========== COMPARISON: RESTRICTED vs BASELINE ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Define file paths
baseline_path <- "D:/processedOuts_feb2026/streamflow_signatures_full_10feb2026.csv"
restricted_path <- "D:/processedOuts_feb2026/streamflow_signatures_full_1993-to-2022-min-20yrs_v2.csv"
output_summary <- "D:/processedOuts_feb2026/comparison_summary_1993-to-2022_v2.txt"
output_plots <- "D:/processedOuts_feb2026/comparison_plots_1993-to-2022_v2.pdf"

# Verify files exist
if (!file.exists(baseline_path)) {
  stop("Baseline file not found: ", baseline_path)
}
if (!file.exists(restricted_path)) {
  stop("Restricted file not found: ", restricted_path)
}

# Load data
cat("Loading baseline data...\n")
baseline <- fread(baseline_path, colClasses = list(character = "gage_id"))
baseline[, gage_id := as.character(gage_id)]

cat("Loading restricted data...\n")
restricted <- fread(restricted_path, colClasses = list(character = "gage_id"))
restricted[, gage_id := as.character(gage_id)]

cat("\n")

# ============== A. GAGE COUNT SUMMARY ==============

cat("========== A. GAGE COUNT SUMMARY ==========\n")

n_baseline <- nrow(baseline)
n_restricted <- nrow(restricted)

gages_baseline <- baseline$gage_id
gages_restricted <- restricted$gage_id

gages_both <- intersect(gages_baseline, gages_restricted)
gages_baseline_only <- setdiff(gages_baseline, gages_restricted)
gages_restricted_only <- setdiff(gages_restricted, gages_baseline)

cat("Total gages in baseline:", n_baseline, "\n")
cat("Total gages in restricted (1993-2022):", n_restricted, "\n")
cat("Gages in both:", length(gages_both), "\n")
cat("Gages in baseline only (excluded by restriction):", length(gages_baseline_only), "\n")
cat("Gages in restricted only (gained):", length(gages_restricted_only), "\n\n")

# Breakdown by gage type
if ("gage_type" %in% names(baseline)) {
  cat("By gage type:\n")
  cat("\nBaseline:\n")
  print(table(baseline$gage_type))
  cat("\nRestricted:\n")
  print(table(restricted$gage_type))
}

cat("\n")

# ============== B. SUMMARY STATISTICS ==============

cat("========== B. SUMMARY STATISTICS ==========\n")

# Key signatures to compare
key_signatures <- c("Qann", "BFI_Eckhardt", "flashinessRB", "D50_day")
key_stats <- c("_mean", "_senn_slp")

stats_comparison <- data.table()

for (sig in key_signatures) {
  for (stat in key_stats) {
    col_name <- paste0(sig, stat)
    if (col_name %in% names(baseline) && col_name %in% names(restricted)) {
      baseline_vals <- baseline[[col_name]]
      restricted_vals <- restricted[[col_name]]

      comparison_row <- data.table(
        signature = sig,
        statistic = gsub("^_", "", stat),
        baseline_mean = mean(baseline_vals, na.rm = TRUE),
        baseline_median = median(baseline_vals, na.rm = TRUE),
        baseline_sd = sd(baseline_vals, na.rm = TRUE),
        restricted_mean = mean(restricted_vals, na.rm = TRUE),
        restricted_median = median(restricted_vals, na.rm = TRUE),
        restricted_sd = sd(restricted_vals, na.rm = TRUE),
        diff_mean = mean(restricted_vals, na.rm = TRUE) - mean(baseline_vals, na.rm = TRUE),
        diff_median = median(restricted_vals, na.rm = TRUE) - median(baseline_vals, na.rm = TRUE)
      )
      stats_comparison <- rbind(stats_comparison, comparison_row)
    }
  }
}

# Climate signature if available
if ("elasticity_static" %in% names(baseline) && "elasticity_static" %in% names(restricted)) {
  baseline_vals <- baseline$elasticity_static
  restricted_vals <- restricted$elasticity_static

  comparison_row <- data.table(
    signature = "elasticity_static",
    statistic = "value",
    baseline_mean = mean(baseline_vals, na.rm = TRUE),
    baseline_median = median(baseline_vals, na.rm = TRUE),
    baseline_sd = sd(baseline_vals, na.rm = TRUE),
    restricted_mean = mean(restricted_vals, na.rm = TRUE),
    restricted_median = median(restricted_vals, na.rm = TRUE),
    restricted_sd = sd(restricted_vals, na.rm = TRUE),
    diff_mean = mean(restricted_vals, na.rm = TRUE) - mean(baseline_vals, na.rm = TRUE),
    diff_median = median(restricted_vals, na.rm = TRUE) - median(baseline_vals, na.rm = TRUE)
  )
  stats_comparison <- rbind(stats_comparison, comparison_row)
}

print(stats_comparison, digits = 3)
cat("\n")

# ============== C. TREND COMPARISON (Common Gages Only) ==============

cat("========== C. TREND COMPARISON ==========\n")

# Merge common gages
baseline_common <- baseline[gage_id %in% gages_both]
restricted_common <- restricted[gage_id %in% gages_both]

# Ensure same order
setkey(baseline_common, gage_id)
setkey(restricted_common, gage_id)

# Compare Theil-Sen slopes for key signatures
trend_cols <- paste0(c("Qann", "BFI_Eckhardt", "flashinessRB", "D50_day"), "_senn_slp")

trend_comparison <- data.table()

for (col in trend_cols) {
  if (col %in% names(baseline_common) && col %in% names(restricted_common)) {
    baseline_trend <- baseline_common[[col]]
    restricted_trend <- restricted_common[[col]]

    # Correlation between trends
    valid_idx <- !is.na(baseline_trend) & !is.na(restricted_trend)
    if (sum(valid_idx) > 10) {
      cor_val <- cor(baseline_trend[valid_idx], restricted_trend[valid_idx], method = "spearman")

      trend_row <- data.table(
        metric = col,
        n_valid = sum(valid_idx),
        spearman_cor = cor_val,
        baseline_mean_trend = mean(baseline_trend, na.rm = TRUE),
        restricted_mean_trend = mean(restricted_trend, na.rm = TRUE),
        trend_diff_mean = mean(restricted_trend[valid_idx] - baseline_trend[valid_idx])
      )
      trend_comparison <- rbind(trend_comparison, trend_row)
    }
  }
}

cat("Trend correlation between baseline and restricted (common gages):\n")
print(trend_comparison, digits = 3)
cat("\n")

# Identify gages with largest trend differences
cat("Gages with largest Qann trend differences:\n")
if ("Qann_senn_slp" %in% names(baseline_common) && "Qann_senn_slp" %in% names(restricted_common)) {
  trend_diff <- data.table(
    gage_id = baseline_common$gage_id,
    baseline_trend = baseline_common$Qann_senn_slp,
    restricted_trend = restricted_common$Qann_senn_slp,
    diff = restricted_common$Qann_senn_slp - baseline_common$Qann_senn_slp
  )
  trend_diff <- trend_diff[!is.na(diff)]
  trend_diff <- trend_diff[order(-abs(diff))]
  print(head(trend_diff, 10))
}
cat("\n")

# ============== D. SPATIAL COVERAGE MAP ==============

cat("========== D. GENERATING PLOTS ==========\n")

# Create PDF with plots
pdf(output_plots, width = 11, height = 8.5)

# Plot 1: Spatial coverage map
if ("latitude" %in% names(baseline) && "longitude" %in% names(baseline)) {

  # Create combined dataset with membership info
  spatial_data <- data.table(
    gage_id = c(gages_baseline, gages_restricted_only),
    status = c(
      ifelse(gages_baseline %in% gages_both, "Both", "Baseline only"),
      rep("Restricted only", length(gages_restricted_only))
    )
  )

  # Merge with coordinates from baseline (or restricted for new gages)
  coords_baseline <- baseline[, .(gage_id, latitude, longitude)]
  coords_restricted <- restricted[, .(gage_id, latitude, longitude)]

  spatial_data <- merge(spatial_data,
                        rbind(coords_baseline,
                              coords_restricted[!gage_id %in% coords_baseline$gage_id]),
                        by = "gage_id", all.x = TRUE)

  # Remove rows without coordinates
  spatial_data <- spatial_data[!is.na(latitude) & !is.na(longitude)]

  p1 <- ggplot(spatial_data, aes(x = longitude, y = latitude, color = status)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Both" = "blue", "Baseline only" = "red", "Restricted only" = "green")) +
    labs(
      title = "Spatial Coverage: Baseline vs Restricted (1993-2022)",
      subtitle = paste0("Blue = Both (", length(gages_both), "), ",
                        "Red = Baseline only (", length(gages_baseline_only), "), ",
                        "Green = Restricted only (", length(gages_restricted_only), ")"),
      x = "Longitude",
      y = "Latitude",
      color = "Coverage"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p1)
}

# Plot 2: Scatter plot of Qann trends
if ("Qann_senn_slp" %in% names(baseline_common) && "Qann_senn_slp" %in% names(restricted_common)) {

  scatter_data <- data.table(
    baseline_trend = baseline_common$Qann_senn_slp,
    restricted_trend = restricted_common$Qann_senn_slp
  )
  scatter_data <- scatter_data[!is.na(baseline_trend) & !is.na(restricted_trend)]

  p2 <- ggplot(scatter_data, aes(x = baseline_trend, y = restricted_trend)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = "Qann Theil-Sen Slope: Baseline vs Restricted",
      subtitle = paste0("N = ", nrow(scatter_data), " common gages, ",
                        "r = ", round(cor(scatter_data$baseline_trend, scatter_data$restricted_trend, method = "spearman"), 3)),
      x = "Baseline Qann Trend (mm/year/year)",
      y = "Restricted (1993-2022) Qann Trend (mm/year/year)"
    ) +
    theme_minimal()

  print(p2)
}

# Plot 3: Distribution of num_water_years
if ("num_water_years" %in% names(baseline) && "num_water_years" %in% names(restricted)) {

  dist_data <- rbind(
    data.table(source = "Baseline", num_water_years = baseline$num_water_years),
    data.table(source = "Restricted (1993-2022)", num_water_years = restricted$num_water_years)
  )

  p3 <- ggplot(dist_data, aes(x = num_water_years, fill = source)) +
    geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
    labs(
      title = "Distribution of Qualifying Water Years",
      x = "Number of Water Years",
      y = "Count",
      fill = "Dataset"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p3)
}

# Plot 4: BFI comparison
if ("BFI_Eckhardt_mean" %in% names(baseline_common) && "BFI_Eckhardt_mean" %in% names(restricted_common)) {

  bfi_data <- data.table(
    baseline_bfi = baseline_common$BFI_Eckhardt_mean,
    restricted_bfi = restricted_common$BFI_Eckhardt_mean
  )
  bfi_data <- bfi_data[!is.na(baseline_bfi) & !is.na(restricted_bfi)]

  p4 <- ggplot(bfi_data, aes(x = baseline_bfi, y = restricted_bfi)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = "BFI Eckhardt Mean: Baseline vs Restricted",
      subtitle = paste0("N = ", nrow(bfi_data), " common gages"),
      x = "Baseline BFI Eckhardt (mean)",
      y = "Restricted (1993-2022) BFI Eckhardt (mean)"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_minimal()

  print(p4)
}

dev.off()
cat("Plots saved to:", output_plots, "\n")

# ============== WRITE SUMMARY FILE ==============

cat("\n========== WRITING SUMMARY FILE ==========\n")

sink(output_summary)
cat("COMPARISON SUMMARY: Restricted (1993-2022) vs Baseline\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================\n\n")

cat("A. GAGE COUNTS\n")
cat("--------------\n")
cat("Baseline gages:", n_baseline, "\n")
cat("Restricted gages:", n_restricted, "\n")
cat("Common gages:", length(gages_both), "\n")
cat("Baseline only (excluded):", length(gages_baseline_only), "\n")
cat("Restricted only (gained):", length(gages_restricted_only), "\n\n")

cat("B. SUMMARY STATISTICS\n")
cat("---------------------\n")
print(stats_comparison, digits = 3)
cat("\n")

cat("C. TREND COMPARISON (Common Gages)\n")
cat("----------------------------------\n")
print(trend_comparison, digits = 3)
cat("\n")

cat("D. OUTPUTS\n")
cat("----------\n")
cat("Summary file:", output_summary, "\n")
cat("Plots file:", output_plots, "\n")

sink()

cat("Summary saved to:", output_summary, "\n")
cat("\n========== COMPARISON COMPLETE ==========\n")
