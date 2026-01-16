################################################################################
# Visual Assessment for Streamflow Signatures QA/QC
#
# Generates diagnostic plots for hydrological validation
# Output: data_out/qa_plots/ folder with PNG files
################################################################################

# Clear environment
rm(list = ls())

# Set working directory
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  main_dir <- dirname(main_dir)
  setwd(main_dir)
}

# Load libraries
library(data.table)
library(ggplot2)
library(gridExtra)
library(viridis)

cat("========== STREAMFLOW SIGNATURES VISUAL ASSESSMENT ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Define paths
input_path <- "data_out/streamflow_signatures_full_JAN2026.csv"
plot_dir <- "data_out/qa_plots"

# Check if input exists
if (!file.exists(input_path)) {
  stop("Signature file not found: ", input_path, "\n  Run run_full_processing.R first")
}

# Create output directory
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading signature data...\n")
signatures <- fread(input_path)
cat("  Loaded", nrow(signatures), "gages\n\n")

# Set theme
theme_set(theme_minimal(base_size = 12))

################################################################################
# 1. DISTRIBUTION PLOTS
################################################################################
cat("Generating distribution plots...\n")

# 1.1 Annual mean flow distribution
if ("Qann_mean" %in% names(signatures)) {
  p1 <- ggplot(signatures[!is.na(Qann_mean)], aes(x = Qann_mean)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = median(signatures$Qann_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", size = 1) +
    scale_x_log10() +
    labs(title = "Annual Mean Streamflow (Qann)",
         subtitle = paste("n =", sum(!is.na(signatures$Qann_mean)),
                         "| Median =", round(median(signatures$Qann_mean, na.rm = TRUE), 2), "mm/day"),
         x = "Qann (mm/day, log scale)", y = "Count") +
    annotation_logticks(sides = "b")

  ggsave(file.path(plot_dir, "01_qann_distribution.png"), p1, width = 10, height = 6, dpi = 150)
}

# 1.2 BFI distribution (both methods)
if (all(c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean") %in% names(signatures))) {
  bfi_long <- melt(signatures[, .(BFI_Eckhardt_mean, BFI_LyneHollick_mean)],
                   measure.vars = c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean"),
                   variable.name = "Method", value.name = "BFI")
  bfi_long <- bfi_long[!is.na(BFI)]
  bfi_long[, Method := gsub("_mean", "", Method)]

  p2 <- ggplot(bfi_long, aes(x = BFI, fill = Method)) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    scale_fill_viridis_d(option = "D") +
    labs(title = "Baseflow Index Distribution",
         subtitle = "Comparison of Eckhardt and Lyne-Hollick methods",
         x = "Baseflow Index (0-1)", y = "Count") +
    theme(legend.position = "top")

  ggsave(file.path(plot_dir, "02_bfi_distribution.png"), p2, width = 10, height = 6, dpi = 150)
}

# 1.3 Flashiness distribution
if ("flashinessRB_mean" %in% names(signatures)) {
  p3 <- ggplot(signatures[!is.na(flashinessRB_mean)], aes(x = flashinessRB_mean)) +
    geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
    geom_vline(xintercept = median(signatures$flashinessRB_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Richards-Baker Flashiness Index",
         subtitle = paste("n =", sum(!is.na(signatures$flashinessRB_mean)),
                         "| Median =", round(median(signatures$flashinessRB_mean, na.rm = TRUE), 3)),
         x = "R-B Flashiness Index", y = "Count")

  ggsave(file.path(plot_dir, "03_flashiness_distribution.png"), p3, width = 10, height = 6, dpi = 150)
}

# 1.4 Elasticity distribution
if ("elasticity_static" %in% names(signatures)) {
  p4 <- ggplot(signatures[!is.na(elasticity_static) & elasticity_static > 0 & elasticity_static < 5],
               aes(x = elasticity_static)) +
    geom_histogram(bins = 40, fill = "darkgreen", alpha = 0.7) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median(signatures$elasticity_static, na.rm = TRUE),
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Streamflow Elasticity Distribution",
         subtitle = "Red line = proportional response (E=1), Blue = median",
         x = "Elasticity (dQ/dP / Q/P)", y = "Count")

  ggsave(file.path(plot_dir, "04_elasticity_distribution.png"), p4, width = 10, height = 6, dpi = 150)
}

################################################################################
# 2. CORRELATION PLOTS
################################################################################
cat("Generating correlation plots...\n")

# 2.1 BFI method comparison
if (all(c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean") %in% names(signatures))) {
  valid <- signatures[!is.na(BFI_Eckhardt_mean) & !is.na(BFI_LyneHollick_mean)]
  r_val <- cor(valid$BFI_Eckhardt_mean, valid$BFI_LyneHollick_mean)

  p5 <- ggplot(valid, aes(x = BFI_Eckhardt_mean, y = BFI_LyneHollick_mean)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "BFI Method Comparison",
         subtitle = paste("r =", round(r_val, 3), "| n =", nrow(valid)),
         x = "BFI Eckhardt", y = "BFI Lyne-Hollick") +
    coord_fixed()

  ggsave(file.path(plot_dir, "05_bfi_comparison.png"), p5, width = 8, height = 8, dpi = 150)
}

# 2.2 Theil-Sen vs Linear slope comparison
senn_col <- "Qann_senn_slp"
linear_col <- "Qann_linear_slp"
if (all(c(senn_col, linear_col) %in% names(signatures))) {
  valid <- signatures[!is.na(get(senn_col)) & !is.na(get(linear_col))]
  setnames(valid, c(senn_col, linear_col), c("senn", "linear"), skip_absent = TRUE)
  r_val <- cor(valid$senn, valid$linear)

  p6 <- ggplot(valid, aes(x = senn, y = linear)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "Trend Method Comparison: Qann",
         subtitle = paste("r =", round(r_val, 3), "| Theil-Sen vs Linear regression"),
         x = "Theil-Sen Slope", y = "Linear Regression Slope")

  ggsave(file.path(plot_dir, "06_slope_comparison.png"), p6, width = 8, height = 8, dpi = 150)
}

# 2.3 BFI vs Flashiness (expect negative correlation)
if (all(c("BFI_Eckhardt_mean", "flashinessRB_mean") %in% names(signatures))) {
  valid <- signatures[!is.na(BFI_Eckhardt_mean) & !is.na(flashinessRB_mean)]
  r_val <- cor(valid$BFI_Eckhardt_mean, valid$flashinessRB_mean)

  p7 <- ggplot(valid, aes(x = BFI_Eckhardt_mean, y = flashinessRB_mean)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = "BFI vs Flashiness",
         subtitle = paste("r =", round(r_val, 3), "(expect negative - baseflow dampens flashiness)"),
         x = "Baseflow Index (Eckhardt)", y = "R-B Flashiness Index")

  ggsave(file.path(plot_dir, "07_bfi_vs_flashiness.png"), p7, width = 10, height = 6, dpi = 150)
}

################################################################################
# 3. SPATIAL PLOTS (if coordinates available)
################################################################################
cat("Generating spatial plots...\n")

if (all(c("longitude", "latitude") %in% names(signatures))) {
  valid_coords <- signatures[!is.na(longitude) & !is.na(latitude) &
                              longitude > -180 & longitude < -50 &
                              latitude > 20 & latitude < 80]

  # 3.1 Map of Qann
  if ("Qann_mean" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    p8 <- ggplot(valid_coords[!is.na(Qann_mean)],
                 aes(x = longitude, y = latitude, color = log10(Qann_mean))) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "C", name = "log10(Qann)") +
      labs(title = "Spatial Distribution of Annual Mean Flow",
           subtitle = paste("n =", sum(!is.na(valid_coords$Qann_mean))),
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "08_map_qann.png"), p8, width = 12, height = 8, dpi = 150)
  }

  # 3.2 Map of BFI
  if ("BFI_Eckhardt_mean" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    p9 <- ggplot(valid_coords[!is.na(BFI_Eckhardt_mean)],
                 aes(x = longitude, y = latitude, color = BFI_Eckhardt_mean)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "D", name = "BFI") +
      labs(title = "Spatial Distribution of Baseflow Index",
           subtitle = "Higher values = more groundwater-dominated",
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "09_map_bfi.png"), p9, width = 12, height = 8, dpi = 150)
  }

  # 3.3 Map of Elasticity
  if ("elasticity_static" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    p10 <- ggplot(valid_coords[!is.na(elasticity_static) &
                               elasticity_static > 0 & elasticity_static < 5],
                  aes(x = longitude, y = latitude, color = elasticity_static)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "B", name = "Elasticity") +
      labs(title = "Spatial Distribution of Streamflow Elasticity",
           subtitle = "Higher values = more sensitive to precipitation",
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "10_map_elasticity.png"), p10, width = 12, height = 8, dpi = 150)
  }
}

################################################################################
# 4. REGIONAL COMPARISONS
################################################################################
cat("Generating regional comparison plots...\n")

# 4.1 Boxplots by gage type
if ("gage_type" %in% names(signatures)) {

  # Qann by gage type
  if ("Qann_mean" %in% names(signatures)) {
    p11 <- ggplot(signatures[!is.na(Qann_mean) & !is.na(gage_type)],
                  aes(x = gage_type, y = Qann_mean, fill = gage_type)) +
      geom_boxplot(outlier.alpha = 0.3) +
      scale_y_log10() +
      scale_fill_viridis_d() +
      labs(title = "Annual Mean Flow by Gage Type",
           x = "Gage Type", y = "Qann (mm/day, log scale)") +
      theme(legend.position = "none")

    ggsave(file.path(plot_dir, "11_qann_by_type.png"), p11, width = 8, height = 6, dpi = 150)
  }

  # BFI by gage type
  if ("BFI_Eckhardt_mean" %in% names(signatures)) {
    p12 <- ggplot(signatures[!is.na(BFI_Eckhardt_mean) & !is.na(gage_type)],
                  aes(x = gage_type, y = BFI_Eckhardt_mean, fill = gage_type)) +
      geom_boxplot(outlier.alpha = 0.3) +
      scale_fill_viridis_d() +
      labs(title = "Baseflow Index by Gage Type",
           x = "Gage Type", y = "BFI (Eckhardt)") +
      theme(legend.position = "none")

    ggsave(file.path(plot_dir, "12_bfi_by_type.png"), p12, width = 8, height = 6, dpi = 150)
  }
}

# 4.2 Latitude bands (if coordinates available)
if (all(c("latitude", "Qann_mean") %in% names(signatures))) {
  lat_data <- signatures[!is.na(latitude) & !is.na(Qann_mean)]
  lat_data[, lat_band := cut(latitude, breaks = seq(25, 70, by = 5),
                             labels = paste0(seq(25, 65, by = 5), "-",
                                           seq(30, 70, by = 5), "N"))]

  if (nrow(lat_data[!is.na(lat_band)]) > 0) {
    p13 <- ggplot(lat_data[!is.na(lat_band)], aes(x = lat_band, y = Qann_mean, fill = lat_band)) +
      geom_boxplot(outlier.alpha = 0.3) +
      scale_y_log10() +
      scale_fill_viridis_d() +
      labs(title = "Annual Mean Flow by Latitude Band",
           x = "Latitude Band", y = "Qann (mm/day, log scale)") +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(plot_dir, "13_qann_by_latitude.png"), p13, width = 10, height = 6, dpi = 150)
  }
}

################################################################################
# 5. SUMMARY DASHBOARD
################################################################################
cat("Generating summary dashboard...\n")

# Create a multi-panel summary plot
if (all(c("Qann_mean", "BFI_Eckhardt_mean", "flashinessRB_mean") %in% names(signatures))) {

  # Panel 1: Qann histogram
  p_sum1 <- ggplot(signatures[!is.na(Qann_mean)], aes(x = Qann_mean)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    scale_x_log10() +
    labs(title = "Qann Distribution", x = "mm/day (log)", y = "Count") +
    theme(plot.title = element_text(size = 10))

  # Panel 2: BFI histogram
  p_sum2 <- ggplot(signatures[!is.na(BFI_Eckhardt_mean)], aes(x = BFI_Eckhardt_mean)) +
    geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
    labs(title = "BFI Distribution", x = "BFI", y = "Count") +
    theme(plot.title = element_text(size = 10))

  # Panel 3: Flashiness histogram
  p_sum3 <- ggplot(signatures[!is.na(flashinessRB_mean)], aes(x = flashinessRB_mean)) +
    geom_histogram(bins = 30, fill = "darkorange", alpha = 0.7) +
    labs(title = "Flashiness Distribution", x = "R-B Index", y = "Count") +
    theme(plot.title = element_text(size = 10))

  # Panel 4: BFI vs Flashiness scatter
  valid <- signatures[!is.na(BFI_Eckhardt_mean) & !is.na(flashinessRB_mean)]
  p_sum4 <- ggplot(valid, aes(x = BFI_Eckhardt_mean, y = flashinessRB_mean)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = "BFI vs Flashiness", x = "BFI", y = "Flashiness") +
    theme(plot.title = element_text(size = 10))

  # Combine panels
  dashboard <- grid.arrange(p_sum1, p_sum2, p_sum3, p_sum4, ncol = 2,
                            top = "Streamflow Signatures Summary Dashboard")

  ggsave(file.path(plot_dir, "00_summary_dashboard.png"), dashboard,
         width = 12, height = 10, dpi = 150)
}

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n========== VISUALIZATION COMPLETE ==========\n")
plots_generated <- length(list.files(plot_dir, pattern = "\\.png$"))
cat("Generated", plots_generated, "plots in:", plot_dir, "\n")

# List all plots
cat("\nPlots generated:\n")
for (f in sort(list.files(plot_dir, pattern = "\\.png$"))) {
  cat("  ", f, "\n")
}

cat("\n========== DONE ==========\n")
