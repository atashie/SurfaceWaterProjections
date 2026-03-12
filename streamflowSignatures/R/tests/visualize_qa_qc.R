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
  if (!file.exists(file.path(main_dir, "config.R"))) {
    main_dir <- dirname(main_dir)
  }
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

# Set theme with white background
theme_set(theme_minimal(base_size = 12) +
            theme(plot.background = element_rect(fill = "white", color = NA),
                  panel.background = element_rect(fill = "white", color = NA)))

################################################################################
# HELPER FUNCTION: Filter flagged data
################################################################################

# Filter out flagged gages and return data with counts
filter_flagged <- function(data, value_col, flag_col) {
  # If flag column doesn't exist, return all non-NA data
  if (!flag_col %in% names(data)) {
    filtered <- data[!is.na(get(value_col))]
    return(list(data = filtered, n_flagged = 0, n_total = nrow(filtered)))
  }

  # Count flagged observations
  n_flagged <- sum(data[[flag_col]] == TRUE, na.rm = TRUE)

  # Filter to keep only non-flagged, non-NA observations
  filtered <- data[get(flag_col) == FALSE & !is.na(get(value_col))]

  list(data = filtered, n_flagged = n_flagged, n_total = nrow(filtered))
}

################################################################################
# 1. DISTRIBUTION PLOTS
################################################################################
cat("Generating distribution plots...\n")

# 1.1 Annual mean runoff distribution
if ("Qann_mean" %in% names(signatures)) {
  filt <- filter_flagged(signatures, "Qann_mean", "flagged_for_qann_range")

  p1 <- ggplot(filt$data, aes(x = Qann_mean)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = median(filt$data$Qann_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", size = 1) +
    scale_x_log10() +
    labs(title = "Annual Mean Runoff (Qann)",
         subtitle = paste("n =", filt$n_total,
                         "| Excluded:", filt$n_flagged, "flagged",
                         "| Median =", round(median(filt$data$Qann_mean, na.rm = TRUE), 1), "mm/yr"),
         x = "Qann (mm/yr, log scale)", y = "Count") +
    annotation_logticks(sides = "b")

  ggsave(file.path(plot_dir, "01a_qann_distribution.png"), p1, width = 10, height = 6, dpi = 150, bg = "white")
}

# 1.2 BFI distribution (both methods)
if (all(c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean") %in% names(signatures))) {
  # Filter for Eckhardt method flagging (primary BFI method)
  filt <- filter_flagged(signatures, "BFI_Eckhardt_mean", "flagged_for_bfi_eckhardt_range")

  # Keep only rows that also have Lyne-Hollick values and aren't flagged for it
  bfi_data <- filt$data
  if ("flagged_for_bfi_lynehollick_range" %in% names(bfi_data)) {
    bfi_data <- bfi_data[flagged_for_bfi_lynehollick_range == FALSE]
  }

  bfi_long <- melt(bfi_data[, .(BFI_Eckhardt_mean, BFI_LyneHollick_mean)],
                   measure.vars = c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean"),
                   variable.name = "Method", value.name = "BFI")
  bfi_long <- bfi_long[!is.na(BFI)]
  bfi_long[, Method := gsub("_mean", "", Method)]

  p2 <- ggplot(bfi_long, aes(x = BFI, fill = Method)) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    scale_fill_viridis_d(option = "D") +
    labs(title = "Baseflow Index Distribution",
         subtitle = paste("Comparison of Eckhardt and Lyne-Hollick methods | n =", nrow(bfi_data),
                         "| Excluded:", filt$n_flagged, "flagged"),
         x = "Baseflow Index (0-1)", y = "Count") +
    theme(legend.position = "top")

  ggsave(file.path(plot_dir, "02a_bfi_distribution.png"), p2, width = 10, height = 6, dpi = 150, bg = "white")
}

# 1.3 Flashiness distribution
if ("flashinessRB_mean" %in% names(signatures)) {
  filt <- filter_flagged(signatures, "flashinessRB_mean", "flagged_for_flashiness_range")

  p3 <- ggplot(filt$data, aes(x = flashinessRB_mean)) +
    geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
    geom_vline(xintercept = median(filt$data$flashinessRB_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Richards-Baker Flashiness Index",
         subtitle = paste("n =", filt$n_total,
                         "| Excluded:", filt$n_flagged, "flagged",
                         "| Median =", round(median(filt$data$flashinessRB_mean, na.rm = TRUE), 3)),
         x = "R-B Flashiness Index", y = "Count")

  ggsave(file.path(plot_dir, "03a_flashiness_distribution.png"), p3, width = 10, height = 6, dpi = 150, bg = "white")
}

# 1.4 Elasticity distribution
if ("elasticity_static" %in% names(signatures)) {
  filt <- filter_flagged(signatures, "elasticity_static", "flagged_for_elasticity_range")

  p4 <- ggplot(filt$data, aes(x = elasticity_static)) +
    geom_histogram(bins = 40, fill = "darkgreen", alpha = 0.7) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median(filt$data$elasticity_static, na.rm = TRUE),
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Streamflow Elasticity Distribution",
         subtitle = paste("n =", filt$n_total,
                         "| Excluded:", filt$n_flagged, "flagged",
                         "| Red = E=1 (proportional), Blue = median"),
         x = "Elasticity (dQ/dP / Q/P)", y = "Count")

  ggsave(file.path(plot_dir, "04a_elasticity_distribution.png"), p4, width = 10, height = 6, dpi = 150, bg = "white")
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

  ggsave(file.path(plot_dir, "05_bfi_comparison.png"), p5, width = 8, height = 8, dpi = 150, bg = "white")
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

  ggsave(file.path(plot_dir, "06_slope_comparison.png"), p6, width = 8, height = 8, dpi = 150, bg = "white")
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

  ggsave(file.path(plot_dir, "07_bfi_vs_flashiness.png"), p7, width = 10, height = 6, dpi = 150, bg = "white")
}

################################################################################
# 3. SPATIAL PLOTS (paired with distributions)
################################################################################
cat("Generating spatial plots...\n")

if (all(c("longitude", "latitude") %in% names(signatures))) {
  valid_coords <- signatures[!is.na(longitude) & !is.na(latitude) &
                              longitude > -180 & longitude < -50 &
                              latitude > 20 & latitude < 80]

  # 3.1 Map of Qann (paired with 01a)
  if ("Qann_mean" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    filt <- filter_flagged(valid_coords, "Qann_mean", "flagged_for_qann_range")

    p_qann_map <- ggplot(filt$data,
                 aes(x = longitude, y = latitude, color = log10(Qann_mean))) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "C", name = "log10(Qann)") +
      labs(title = "Spatial Distribution of Annual Mean Runoff",
           subtitle = paste("n =", filt$n_total, "| Excluded:", filt$n_flagged, "flagged"),
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "01b_qann_map.png"), p_qann_map, width = 12, height = 8, dpi = 150, bg = "white")
  }

  # 3.2 Map of BFI (paired with 02a)
  if ("BFI_Eckhardt_mean" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    filt <- filter_flagged(valid_coords, "BFI_Eckhardt_mean", "flagged_for_bfi_eckhardt_range")

    p_bfi_map <- ggplot(filt$data,
                 aes(x = longitude, y = latitude, color = BFI_Eckhardt_mean)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "D", name = "BFI") +
      labs(title = "Spatial Distribution of Baseflow Index",
           subtitle = paste("n =", filt$n_total, "| Excluded:", filt$n_flagged, "flagged",
                          "| Higher = groundwater-dominated"),
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "02b_bfi_map.png"), p_bfi_map, width = 12, height = 8, dpi = 150, bg = "white")
  }

  # 3.3 Map of Flashiness (paired with 03a) - NEW
  if ("flashinessRB_mean" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    filt <- filter_flagged(valid_coords, "flashinessRB_mean", "flagged_for_flashiness_range")

    p_flash_map <- ggplot(filt$data,
                 aes(x = longitude, y = latitude, color = flashinessRB_mean)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "A", name = "Flashiness") +
      labs(title = "Spatial Distribution of Flashiness",
           subtitle = paste("n =", filt$n_total, "| Excluded:", filt$n_flagged, "flagged",
                          "| Higher = more variable flow"),
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "03b_flashiness_map.png"), p_flash_map, width = 12, height = 8, dpi = 150, bg = "white")
  }

  # 3.4 Map of Elasticity (paired with 04a)
  if ("elasticity_static" %in% names(valid_coords) && nrow(valid_coords) > 0) {
    filt <- filter_flagged(valid_coords, "elasticity_static", "flagged_for_elasticity_range")

    p_elast_map <- ggplot(filt$data,
                  aes(x = longitude, y = latitude, color = elasticity_static)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis(option = "B", name = "Elasticity") +
      labs(title = "Spatial Distribution of Streamflow Elasticity",
           subtitle = paste("n =", filt$n_total, "| Excluded:", filt$n_flagged, "flagged",
                          "| Higher = more sensitive to precipitation"),
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "04b_elasticity_map.png"), p_elast_map, width = 12, height = 8, dpi = 150, bg = "white")
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

    ggsave(file.path(plot_dir, "08_qann_by_type.png"), p11, width = 8, height = 6, dpi = 150, bg = "white")
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

    ggsave(file.path(plot_dir, "09_bfi_by_type.png"), p12, width = 8, height = 6, dpi = 150, bg = "white")
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

    ggsave(file.path(plot_dir, "10_qann_by_latitude.png"), p13, width = 10, height = 6, dpi = 150, bg = "white")
  }
}

################################################################################
# 5. ADDITIONAL METRICS (Seasonal, Bimodality, Runoff Ratios, Flow Reversals)
################################################################################
cat("Generating additional metric plots...\n")

# 5.1 Seasonal Flow Proportions (calculated from Qseason / Qann)
seasonal_cols <- c("Qwin_mean", "Qspr_mean", "Qsum_mean", "Qfal_mean", "Qann_mean")
if (all(seasonal_cols %in% names(signatures))) {
  # Calculate seasonal proportions
  prop_data <- signatures[!is.na(Qann_mean) & Qann_mean > 0,
                          .(gage_id, Qann_mean, Qwin_mean, Qspr_mean, Qsum_mean, Qfal_mean,
                            longitude, latitude)]
  prop_data[, `:=`(
    Winter = Qwin_mean / Qann_mean,
    Spring = Qspr_mean / Qann_mean,
    Summer = Qsum_mean / Qann_mean,
    Fall = Qfal_mean / Qann_mean
  )]

  # Melt for plotting
  prop_long <- melt(prop_data[, .(gage_id, Winter, Spring, Summer, Fall, longitude, latitude)],
                    id.vars = c("gage_id", "longitude", "latitude"),
                    variable.name = "Season", value.name = "Proportion")
  prop_long <- prop_long[!is.na(Proportion) & Proportion >= 0 & Proportion <= 1]

  # 11a: Seasonal proportions distribution (faceted density)
  p_seasonal <- ggplot(prop_long, aes(x = Proportion, fill = Season)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("Winter" = "#4575b4", "Spring" = "#91cf60",
                                 "Summer" = "#fc8d59", "Fall" = "#d73027")) +
    labs(title = "Seasonal Flow Proportions (Qseason / Qann)",
         subtitle = paste("n =", nrow(prop_data), "gages"),
         x = "Proportion of Annual Flow", y = "Density") +
    theme(legend.position = "top") +
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "gray50", alpha = 0.5)

  ggsave(file.path(plot_dir, "11a_seasonal_proportions_distribution.png"), p_seasonal,
         width = 10, height = 6, dpi = 150, bg = "white")

  # 11b: Map of dominant season
  prop_data[, dominant_season := {
    seasons <- c("Winter", "Spring", "Summer", "Fall")
    props <- c(Winter, Spring, Summer, Fall)
    seasons[which.max(props)]
  }, by = gage_id]

  if (all(c("longitude", "latitude") %in% names(prop_data))) {
    map_data <- prop_data[!is.na(longitude) & !is.na(latitude) &
                          longitude > -180 & longitude < -50 &
                          latitude > 20 & latitude < 80]
    if (nrow(map_data) > 0) {
      p_dom_map <- ggplot(map_data, aes(x = longitude, y = latitude, color = dominant_season)) +
        geom_point(alpha = 0.6, size = 1) +
        scale_color_manual(values = c("Winter" = "#4575b4", "Spring" = "#91cf60",
                                      "Summer" = "#fc8d59", "Fall" = "#d73027"),
                          name = "Dominant\nSeason") +
        labs(title = "Spatial Distribution of Dominant Flow Season",
             subtitle = paste("n =", nrow(map_data), "gages | Season with highest proportion of annual flow"),
             x = "Longitude", y = "Latitude") +
        coord_quickmap() +
        theme(legend.position = "right")

      ggsave(file.path(plot_dir, "11b_seasonal_dominant_map.png"), p_dom_map,
             width = 12, height = 8, dpi = 150, bg = "white")
    }
  }
}

# 5.2 Q-P Bimodality
if ("qp_bimodality_mean" %in% names(signatures)) {
  filt_bimod <- signatures[!is.na(qp_bimodality_mean)]
  n_total <- nrow(filt_bimod)

  # 12a: Bimodality distribution
  p_bimod <- ggplot(filt_bimod, aes(x = qp_bimodality_mean)) +
    geom_histogram(bins = 40, fill = "purple4", alpha = 0.7) +
    geom_vline(xintercept = median(filt_bimod$qp_bimodality_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Q-P Bimodality Distribution",
         subtitle = paste("n =", n_total,
                         "| Median =", round(median(filt_bimod$qp_bimodality_mean, na.rm = TRUE), 3)),
         x = "Bimodality Index", y = "Count")

  ggsave(file.path(plot_dir, "12a_bimodality_distribution.png"), p_bimod,
         width = 10, height = 6, dpi = 150, bg = "white")

  # 12b: Bimodality map
  if (all(c("longitude", "latitude") %in% names(filt_bimod))) {
    map_bimod <- filt_bimod[!is.na(longitude) & !is.na(latitude) &
                            longitude > -180 & longitude < -50 &
                            latitude > 20 & latitude < 80]
    if (nrow(map_bimod) > 0) {
      p_bimod_map <- ggplot(map_bimod, aes(x = longitude, y = latitude, color = qp_bimodality_mean)) +
        geom_point(alpha = 0.6, size = 1) +
        scale_color_viridis(option = "E", name = "Bimodality") +
        labs(title = "Spatial Distribution of Q-P Bimodality",
             subtitle = paste("n =", nrow(map_bimod), "| Higher = more bimodal flow-precip relationship"),
             x = "Longitude", y = "Latitude") +
        coord_quickmap() +
        theme(legend.position = "right")

      ggsave(file.path(plot_dir, "12b_bimodality_map.png"), p_bimod_map,
             width = 12, height = 8, dpi = 150, bg = "white")
    }
  }
}

# 5.3 Runoff Ratios (Annual and Seasonal)
if ("annual_runoff_ratio_mean" %in% names(signatures)) {
  filt_rr <- filter_flagged(signatures, "annual_runoff_ratio_mean", "flagged_for_runoff_ratio_range")

  # 13a: Annual runoff ratio distribution
  p_rr_ann <- ggplot(filt_rr$data, aes(x = annual_runoff_ratio_mean)) +
    geom_histogram(bins = 50, fill = "darkblue", alpha = 0.7) +
    geom_vline(xintercept = median(filt_rr$data$annual_runoff_ratio_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = 1, color = "orange", linetype = "dotted", linewidth = 1) +
    labs(title = "Annual Runoff Ratio (Q/P)",
         subtitle = paste("n =", filt_rr$n_total,
                         "| Excluded:", filt_rr$n_flagged, "flagged",
                         "| Median =", round(median(filt_rr$data$annual_runoff_ratio_mean, na.rm = TRUE), 3),
                         "| Orange = 1.0 (all precip becomes runoff)"),
         x = "Runoff Ratio (Q/P)", y = "Count")

  ggsave(file.path(plot_dir, "13a_runoff_ratio_annual_distribution.png"), p_rr_ann,
         width = 10, height = 6, dpi = 150, bg = "white")

  # 13b: Annual runoff ratio map
  if (all(c("longitude", "latitude") %in% names(filt_rr$data))) {
    map_rr <- filt_rr$data[!is.na(longitude) & !is.na(latitude) &
                           longitude > -180 & longitude < -50 &
                           latitude > 20 & latitude < 80]
    if (nrow(map_rr) > 0) {
      p_rr_map <- ggplot(map_rr, aes(x = longitude, y = latitude, color = annual_runoff_ratio_mean)) +
        geom_point(alpha = 0.6, size = 1) +
        scale_color_viridis(option = "C", name = "Q/P", limits = c(0, 1.5), oob = scales::squish) +
        labs(title = "Spatial Distribution of Annual Runoff Ratio",
             subtitle = paste("n =", nrow(map_rr), "| Excluded:", filt_rr$n_flagged, "flagged",
                            "| Higher = more efficient runoff generation"),
             x = "Longitude", y = "Latitude") +
        coord_quickmap() +
        theme(legend.position = "right")

      ggsave(file.path(plot_dir, "13b_runoff_ratio_annual_map.png"), p_rr_map,
             width = 12, height = 8, dpi = 150, bg = "white")
    }
  }
}

# 14: Seasonal runoff ratios comparison
seasonal_rr_cols <- c("winter_runoff_ratio_mean", "spring_runoff_ratio_mean",
                      "summer_runoff_ratio_mean", "fall_runoff_ratio_mean")
if (all(seasonal_rr_cols %in% names(signatures))) {
  rr_long <- melt(signatures[, c("gage_id", seasonal_rr_cols), with = FALSE],
                  id.vars = "gage_id",
                  variable.name = "Season", value.name = "RunoffRatio")
  rr_long <- rr_long[!is.na(RunoffRatio) & RunoffRatio >= 0 & RunoffRatio <= 2]
  rr_long[, Season := gsub("_runoff_ratio_mean", "", Season)]
  rr_long[, Season := tools::toTitleCase(Season)]
  rr_long[, Season := factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))]

  p_rr_seas <- ggplot(rr_long, aes(x = Season, y = RunoffRatio, fill = Season)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Winter" = "#4575b4", "Spring" = "#91cf60",
                                 "Summer" = "#fc8d59", "Fall" = "#d73027")) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    labs(title = "Seasonal Runoff Ratios (Q/P)",
         subtitle = paste("n =", length(unique(rr_long$gage_id)), "gages | Dashed line = Q/P = 1"),
         x = "Season", y = "Runoff Ratio (Q/P)") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 2))

  ggsave(file.path(plot_dir, "14_runoff_ratio_seasonal_boxplot.png"), p_rr_seas,
         width = 8, height = 6, dpi = 150, bg = "white")
}

# 5.4 Flow Reversals (Annual and Seasonal)
if ("Flow_Reversals_annual_mean" %in% names(signatures)) {
  filt_fr <- signatures[!is.na(Flow_Reversals_annual_mean)]

  # 15a: Annual flow reversals distribution
  p_fr_ann <- ggplot(filt_fr, aes(x = Flow_Reversals_annual_mean)) +
    geom_histogram(bins = 50, fill = "darkcyan", alpha = 0.7) +
    geom_vline(xintercept = median(filt_fr$Flow_Reversals_annual_mean, na.rm = TRUE),
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Annual Flow Reversals",
         subtitle = paste("n =", nrow(filt_fr),
                         "| Median =", round(median(filt_fr$Flow_Reversals_annual_mean, na.rm = TRUE), 1),
                         "reversals/year"),
         x = "Flow Reversals per Year", y = "Count")

  ggsave(file.path(plot_dir, "15a_flow_reversals_annual_distribution.png"), p_fr_ann,
         width = 10, height = 6, dpi = 150, bg = "white")

  # 15b: Annual flow reversals map
  if (all(c("longitude", "latitude") %in% names(filt_fr))) {
    map_fr <- filt_fr[!is.na(longitude) & !is.na(latitude) &
                      longitude > -180 & longitude < -50 &
                      latitude > 20 & latitude < 80]
    if (nrow(map_fr) > 0) {
      p_fr_map <- ggplot(map_fr, aes(x = longitude, y = latitude, color = Flow_Reversals_annual_mean)) +
        geom_point(alpha = 0.6, size = 1) +
        scale_color_viridis(option = "D", name = "Reversals\n/year") +
        labs(title = "Spatial Distribution of Annual Flow Reversals",
             subtitle = paste("n =", nrow(map_fr), "| Higher = more variable day-to-day flow direction"),
             x = "Longitude", y = "Latitude") +
        coord_quickmap() +
        theme(legend.position = "right")

      ggsave(file.path(plot_dir, "15b_flow_reversals_annual_map.png"), p_fr_map,
             width = 12, height = 8, dpi = 150, bg = "white")
    }
  }
}

# 16: Seasonal flow reversals comparison
seasonal_fr_cols <- c("Flow_Reversals_winter_mean", "Flow_Reversals_spring_mean",
                      "Flow_Reversals_summer_mean", "Flow_Reversals_fall_mean")
if (all(seasonal_fr_cols %in% names(signatures))) {
  fr_long <- melt(signatures[, c("gage_id", seasonal_fr_cols), with = FALSE],
                  id.vars = "gage_id",
                  variable.name = "Season", value.name = "Reversals")
  fr_long <- fr_long[!is.na(Reversals)]
  fr_long[, Season := gsub("Flow_Reversals_|_mean", "", Season)]
  fr_long[, Season := tools::toTitleCase(Season)]
  fr_long[, Season := factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))]

  p_fr_seas <- ggplot(fr_long, aes(x = Season, y = Reversals, fill = Season)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Winter" = "#4575b4", "Spring" = "#91cf60",
                                 "Summer" = "#fc8d59", "Fall" = "#d73027")) +
    labs(title = "Seasonal Flow Reversals",
         subtitle = paste("n =", length(unique(fr_long$gage_id)), "gages"),
         x = "Season", y = "Flow Reversals per Season") +
    theme(legend.position = "none")

  ggsave(file.path(plot_dir, "16_flow_reversals_seasonal_boxplot.png"), p_fr_seas,
         width = 8, height = 6, dpi = 150, bg = "white")
}

################################################################################
# 6. SUMMARY DASHBOARD
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
         width = 12, height = 10, dpi = 150, bg = "white")
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
