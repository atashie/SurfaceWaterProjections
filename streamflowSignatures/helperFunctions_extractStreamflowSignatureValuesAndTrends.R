require(data.table)     # for fread
require(lubridate)		# for dates
require(dataRetrieval)	# for USGS hydro data
require(tidyhydat)		# for canadian hydro data
require(lfstat)			# for the baseflow() function for calculating DFIs
require(hddtools) 		# for grdc catalogue and some summary data; actual grdc daily data must be downloaded separately and stored locally; need to check again later
require(segmented)		# for breakpoint analysis
require(mblm) 			# for theil sen regression
#require(hydrostats)		# for LH13 baseflow separation
#require(FlowScreen)		# for bf_oneparam, bf_eckhardt, and bf_boughton
require(magrittr)
require(zyp)
require(dplyr)
require(tibble)
require(sf)
sf::sf_use_s2(FALSE)
require(ecmwfr)  # For accessing ERA5 data
require(terra)   # For spatial operations
require(ncdf4)




#########################################################################################
#########################################################################################
# SECTION 0: FINAL USER FACING WRAPPER FUNCTION
#########################################################################################
#########################################################################################

# Function to process streamflow signatures from parquet files (long format)
process_signatures_from_parquet <- function(
    parquet_file_path = "combined_streamflow_output/combined_daily_data.parquet",
    metadata_file_path = "combined_streamflow_output/combined_watershed_metadata.csv",
    output_file = "processed_signatures_from_parquet.csv",
    min_Q_value_and_days = c(0.0001, 30),
    min_num_years = 20,
    min_frac_good_data = 0.9  # NEW PARAMETER for filter 2b
) {
  
  # Load required libraries
  require(arrow)
  require(data.table)
  require(lubridate)
  
  cat("Reading parquet file and metadata...\n")
  
  # Read the combined parquet file
  streamflow_data <- arrow::read_parquet(parquet_file_path)
  streamflow_data <- as.data.table(streamflow_data)
  
  # Check the structure of the data
  cat("Data structure:\n")
  cat("Columns:", paste(names(streamflow_data), collapse = ", "), "\n")
  cat("Number of rows:", nrow(streamflow_data), "\n")
  
  # Read metadata - ensure gage_id is character
  metadata <- fread(metadata_file_path, colClasses = list(character = "gage_id"))
  metadata[, gage_id := as.character(gage_id)]
  cat("Number of gages in metadata:", nrow(metadata), "\n")
  
  # Initialize output data table
  summary_output <- data.table()
  
  # Get list of unique gage IDs from the parquet data
  unique_gages <- unique(as.character(streamflow_data$gage_id))
  cat("Found", length(unique_gages), "unique gages in parquet file\n")
  
  # Create a lookup table for metadata with multiple possible formats
  # important for USGS data where IDs often have leading 0s that may be dropped
  metadata_lookup <- data.table()
  for (i in 1:nrow(metadata)) {
    meta_row <- metadata[i]
    base_id <- as.character(meta_row$gage_id)
    
    # Store the original metadata row
    metadata_lookup <- rbind(metadata_lookup, meta_row, fill = TRUE)
    
    # Also store versions with leading zeros removed (if they exist)
    stripped_id <- gsub("^0+", "", base_id)
    if (stripped_id != base_id && stripped_id != "") {
      alt_row <- copy(meta_row)
      alt_row$gage_id <- stripped_id
      metadata_lookup <- rbind(metadata_lookup, alt_row, fill = TRUE)
    }
  }
  
  # Function to find metadata for a gage ID, trying different formats
  find_metadata <- function(gage_id, metadata_lookup) {
    gage_id_to_find <- as.character(gage_id)  # Rename to avoid collision
    
    # First try exact match
    meta <- metadata_lookup[gage_id == gage_id_to_find]  # Now this works correctly!
    if (nrow(meta) > 0) return(meta[1])
    
    # Try adding leading zeros (up to 4)
    for (num_zeros in 1:4) {
      padded_id <- sprintf(paste0("%0", nchar(gage_id_to_find) + num_zeros, "d"), 
                           as.numeric(gage_id_to_find))
      meta <- metadata_lookup[gage_id == padded_id]
      if (nrow(meta) > 0) {
        cat("  Found match for", gage_id_to_find, "as", padded_id, "\n")
        return(meta[1])
      }
    }
    
    # Try removing leading zeros
    stripped_id <- gsub("^0+", "", gage_id_to_find)
    if (stripped_id != gage_id_to_find && stripped_id != "") {
      meta <- metadata_lookup[gage_id == stripped_id]
      if (nrow(meta) > 0) {
        cat("  Found match for", gage_id_to_find, "as", stripped_id, "\n")
        return(meta[1])
      }
    }
    
    return(NULL)
  }
  
  # Track matching statistics
  matched_gages <- 0
  unmatched_gages <- character()
  
  # Process each unique gage
  for (i in seq_along(unique_gages)) {
    gage_id <- unique_gages[i]
    
    if (i %% 100 == 0) {
      cat("Processing gage", i, "of", length(unique_gages), "\n")
    }
    
    tryCatch({
      # Find metadata for this gage
      gage_meta <- find_metadata(gage_id, metadata_lookup)
      
      if (is.null(gage_meta)) {
        unmatched_gages <- c(unmatched_gages, gage_id)
        next
      }
      
      matched_gages <- matched_gages + 1
      
      # Extract streamflow data for this gage
      gage_flow <- streamflow_data[gage_id == unique_gages[i], ]
      
      if (nrow(gage_flow) == 0) {
        cat("  No data for gage", gage_id, "\n")
        next
      }
      
      # Ensure required columns exist
      if (!"Date" %in% names(gage_flow)) {
        cat("  Missing Date column for gage", gage_id, "\n")
        next
      }
      
      if (!"Q" %in% names(gage_flow)) {
        cat("  Missing Q column for gage", gage_id, "\n")
        next
      }
      
      # Remove NA values in Q
      gage_flow <- gage_flow[!is.na(Q)]
      
      if (nrow(gage_flow) == 0) {
        cat("  No valid Q data for gage", gage_id, "\n")
        next
      }
      
      # Check if temporal columns exist, if not create them
      if (!"year" %in% names(gage_flow)) {
        gage_flow[, year := year(Date)]
      }
      if (!"month" %in% names(gage_flow)) {
        gage_flow[, month := month(Date)]
      }
      if (!"doy" %in% names(gage_flow)) {
        gage_flow[, doy := yday(Date)]
      }
      
      # Add water year information if missing
      if (!all(c("water_year", "dowy") %in% names(gage_flow))) {
        wy_info <- calculate_water_year_info(gage_flow$Date)
        gage_flow[, `:=`(
          water_year = wy_info$water_year,
          dowy = wy_info$dowy
        )]
      }
      
      # ============= APPLY WATER YEAR BASED FILTERS =============
      # Filter 2: Check each water year for minimum days above threshold
      # Filter 2b: Check each water year for minimum fraction of non-NA data
      
      water_years_to_use <- NULL
      
      for (this_wy in unique(gage_flow$water_year)) {
        test_wy <- gage_flow[water_year == this_wy]
        
        # Filter 2: Check days with Q > threshold
        nonzero_rows <- sum(test_wy$Q > min_Q_value_and_days[1], na.rm = TRUE)
        if (nonzero_rows < min_Q_value_and_days[2]) {
          cat("  WY", this_wy, "rejected: only", nonzero_rows, 
              "days > threshold (need", min_Q_value_and_days[2], ")\n")
          next
        }
        
        # Filter 2b: Check fraction of good (non-NA) data
        # Determine expected days in water year (account for leap years)
        # Water year runs Oct 1 to Sep 30
        # If this_wy is 2020, it runs from Oct 1, 2019 to Sep 30, 2020
        # Leap year affects the calendar year that ends the water year
        is_leap <- ((this_wy %% 4 == 0) & (this_wy %% 100 != 0)) | (this_wy %% 400 == 0)
        expected_days <- ifelse(is_leap, 366, 365)
        min_required_days <- floor(expected_days * min_frac_good_data)
        
        n_good_days <- nrow(test_wy)  # Already filtered for non-NA Q above
        
        if (n_good_days < min_required_days) {
          cat("  WY", this_wy, "rejected: only", n_good_days, 
              "valid days (need", min_required_days, "for", 
              sprintf("%.1f%%", min_frac_good_data * 100), "coverage)\n")
          next
        }
        
        # This water year passes both filters
        water_years_to_use <- c(water_years_to_use, this_wy)
      }
      
      # Filter 3: Check minimum number of qualifying water years
      if (length(water_years_to_use) < min_num_years) {
        cat("  Insufficient qualifying water years for gage", gage_id, 
            "(", length(water_years_to_use), "water years)\n")
        next
      }
      
      # Filter to qualifying water years
      streamflow_data_filtered <- gage_flow[water_year %in% water_years_to_use]
      
      # Initialize metric results
      metrics_list <- list()
      
      # Calculate each metric group with error handling
      tryCatch({
        metrics_list$flow_vols <- calculate_flow_vols_by_year(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$fdc_trends <- analyze_fdc_trends_from_streamflow(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$flashiness <- analyze_flashiness_trends(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$flow_timing <- analyze_flow_timing_trends(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$pulses <- calculate_pulse_metrics(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$baseflow <- analyze_baseflow_indices(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      tryCatch({
        metrics_list$recession <- analyze_recession_parameters(streamflow_data_filtered)
      }, error = function(e) NULL)
      
      # Create output row
      gage_row <- data.table(
        gage_id = gage_id,
        gage_id_metadata = gage_meta$gage_id,  # Store the metadata version too
        latitude = gage_meta$latitude,
        longitude = gage_meta$longitude, 
        basin_area = gage_meta$basin_area,
        gage_type = gage_meta$gage_type,
        num_water_years = length(water_years_to_use),  # Changed from num_years
        start_water_year = min(water_years_to_use),     # Changed from start_year
        end_water_year = max(water_years_to_use),       # Changed from end_year
        num_upstream_basins = gage_meta$num_upstream_basins
      )
      
      # Add calculated metrics
      for (metric_name in names(metrics_list)) {
        if (!is.null(metrics_list[[metric_name]])) {
          gage_row <- cbind(gage_row, as.data.table(metrics_list[[metric_name]]))
        }
      }
      
      # Add to summary output
      summary_output <- rbind(summary_output, gage_row, fill = TRUE)
      
    }, error = function(e) {
      cat("  Error processing gage", gage_id, ":", e$message, "\n")
    })
    
    # Periodic save
    if (i %% 500 == 0 && nrow(summary_output) > 0) {
      fwrite(summary_output, output_file)
      cat("Saved intermediate results (", nrow(summary_output), "gages processed)\n")
    }
  }
  
  # Final save
  if (nrow(summary_output) > 0) {
    fwrite(summary_output, output_file)
    cat("\n========== PROCESSING COMPLETE ==========\n")
    cat("Total unique gages in parquet:", length(unique_gages), "\n")
    cat("Successfully matched to metadata:", matched_gages, "\n")
    cat("Could not match to metadata:", length(unmatched_gages), "\n")
    cat("Successfully processed:", nrow(summary_output), "\n")
    cat("Results saved to:", output_file, "\n")
    
    # Save list of unmatched gages for debugging
    if (length(unmatched_gages) > 0) {
      unmatched_file <- gsub("\\.csv$", "_unmatched_gages.txt", output_file)
      writeLines(unmatched_gages, unmatched_file)
      cat("List of unmatched gages saved to:", unmatched_file, "\n")
      cat("Sample of unmatched gages:", head(unmatched_gages, 10), "\n")
    }
  } else {
    cat("\n========== NO DATA PROCESSED ==========\n")
    cat("No gages were successfully processed.\n")
  }
  
  return(summary_output)
}




#########################################################################################
#########################################################################################
# SECTION 1: STATISTICAL PROCESSING FUNCTION
             # THE STATS / SIGNALS / TRENDS DEFINED IN THIS FUNCITON WILL BE APPLIED TO ALL SIGNATURES
#########################################################################################
#########################################################################################

  # helper function that receives a time series and outputs summary metrics (trends, averages, etc.)
generate_stats <- function(data, value_cols = NULL, year_col = "year", min_rows = 3) {
  # Check if zyp package is available
  if (!requireNamespace("zyp", quietly = TRUE)) {
    stop("Package 'zyp' is needed for this function. Please install it with install.packages('zyp')")
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame or data.table")
  }
  
  if (!year_col %in% colnames(data)) {
    stop(paste("Year column '", year_col, "' not found in data"))
  }
  
  # If value_cols not specified, use all numeric columns except year
  if (is.null(value_cols)) {
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    value_cols <- setdiff(numeric_cols, year_col)
  }
  
  # Initialize results list
  results <- list()
  
  # Process each value column
  for (col in value_cols) {
    if (!col %in% colnames(data)) {
      warning(paste("Column", col, "not found in data. Skipping."))
      next
    }
    
    # Create working data with only non-NA values
    working_data <- data.frame(
      year = data[[year_col]],
      value = data[[col]]
    )
    working_data <- working_data[!is.na(working_data$value), ]
    
    # Check if we have enough data
    if (nrow(working_data) < min_rows) {
      # Return NAs for all stats
      results[[paste0(col, "_slp")]] <- NA
      results[[paste0(col, "_rho")]] <- NA
      results[[paste0(col, "_pval")]] <- NA
      results[[paste0(col, "_mean")]] <- NA
      results[[paste0(col, "_median")]] <- NA
      next
    }
    
    # Calculate Theil-Sen slope
    sen_result <- try(zyp::zyp.sen(value ~ year, data = working_data), silent = TRUE)
    if (inherits(sen_result, "try-error")) {
      sen_slope <- NA
    } else {
      sen_slope <- sen_result$coefficients[2]
    }
    
    # Calculate Spearman correlation
    spearman_result <- try(cor.test(working_data$year, working_data$value, 
                                    method = "spearman"), silent = TRUE)
    if (inherits(spearman_result, "try-error")) {
      rho <- NA
      pval <- NA
    } else {
      rho <- spearman_result$estimate
      pval <- spearman_result$p.value
    }
    
    # Calculate mean and median
    mean_val <- mean(working_data$value, na.rm = TRUE)
    median_val <- median(working_data$value, na.rm = TRUE)
    
    # Store results
    results[[paste0(col, "_slp")]] <- sen_slope
    results[[paste0(col, "_rho")]] <- rho
    results[[paste0(col, "_pval")]] <- pval
    results[[paste0(col, "_mean")]] <- mean_val
    results[[paste0(col, "_median")]] <- median_val
  }
  
  # Convert to data frame
  return(as.data.frame(results))
}




#########################################################################################
#########################################################################################
# SECTION 2: STREAMFLOW SIGNATURES
             # DEFINE EACH SET OF SIGNATURES THAT WE WISH TO EXTRACT
#########################################################################################
#########################################################################################


#################################################################
# start streamflow signature analysis


calculate_flow_vols_by_year = function(streamflow_data){
  # Ensure required columns exist
  required_cols <- c("water_year", "Q", "month", "dowy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Calculate annual means
  annual_means <- aggregate(Q ~ water_year, data=streamflow_data, FUN=sum, na.rm=TRUE)
  
  # Check if we have any valid annual data
  if (nrow(annual_means) == 0) {
    # Return empty result with correct structure
    return(generate_stats(data.frame(water_year = numeric()), 
                          value_cols = character(), 
                          year_col = "water_year"))
  }
  
  names(annual_means)[2] <- "Qann"
  
  # Calculate seasonal means
  winter <- streamflow_data[streamflow_data$month %in% c(12, 1, 2), ]
  spring <- streamflow_data[streamflow_data$month %in% c(3, 4, 5), ]
  summer <- streamflow_data[streamflow_data$month %in% c(6, 7, 8), ]
  fall <- streamflow_data[streamflow_data$month %in% c(9, 10, 11), ]
  
  winter_means <- aggregate(Q ~ water_year, data=winter, FUN=sum, na.rm=TRUE)
  spring_means <- aggregate(Q ~ water_year, data=spring, FUN=sum, na.rm=TRUE)
  summer_means <- aggregate(Q ~ water_year, data=summer, FUN=sum, na.rm=TRUE)
  fall_means <- aggregate(Q ~ water_year, data=fall, FUN=sum, na.rm=TRUE)
  
  # Rename columns for seasonal means
  if (nrow(winter_means) > 0) names(winter_means)[2] <- "Qwin"
  if (nrow(spring_means) > 0) names(spring_means)[2] <- "Qspr"
  if (nrow(summer_means) > 0) names(summer_means)[2] <- "Qsum"
  if (nrow(fall_means) > 0) names(fall_means)[2] <- "Qfal"
  
  # After seasonal aggregations, check if any are empty
  if (nrow(winter_means) == 0) winter_means <- data.frame(water_year=numeric(), Qwin=numeric())
  if (nrow(spring_means) == 0) spring_means <- data.frame(water_year=numeric(), Qspr=numeric())
  if (nrow(summer_means) == 0) summer_means <- data.frame(water_year=numeric(), Qsum=numeric())
  if (nrow(fall_means) == 0) fall_means <- data.frame(water_year=numeric(), Qfal=numeric())
  
  # Calculate flow percentiles by year with error handling
  calculate_percentile <- function(data, percentile) {
    if (nrow(data) == 0 || all(is.na(data$Q))) {
      return(data.frame(water_year = numeric(), Q = numeric()))
    }
    
    agg <- aggregate(Q ~ water_year, data=data, 
                     FUN=function(x) {
                       if (all(is.na(x))) return(NA)
                       quantile(x, probs=percentile/100, na.rm=TRUE)
                     })
    
    if (nrow(agg) > 0) {
      names(agg)[2] <- paste0("Q", percentile)
    }
    return(agg)
  }
  
  # Calculate all percentiles with error handling
  q1 <- calculate_percentile(streamflow_data, 1)
  q5 <- calculate_percentile(streamflow_data, 5)
  q10 <- calculate_percentile(streamflow_data, 10)
  q20 <- calculate_percentile(streamflow_data, 20)
  q25 <- calculate_percentile(streamflow_data, 25)
  q30 <- calculate_percentile(streamflow_data, 30)
  q40 <- calculate_percentile(streamflow_data, 40)
  q50 <- calculate_percentile(streamflow_data, 50)
  q60 <- calculate_percentile(streamflow_data, 60)
  q70 <- calculate_percentile(streamflow_data, 70)
  q75 <- calculate_percentile(streamflow_data, 75)
  q80 <- calculate_percentile(streamflow_data, 80)
  q90 <- calculate_percentile(streamflow_data, 90)
  q95 <- calculate_percentile(streamflow_data, 95)
  q99 <- calculate_percentile(streamflow_data, 99)
  
  # Calculate Q95-Q10 difference by year with error handling
  if (nrow(q95) > 0 && nrow(q10) > 0) {
    q95_q10 <- merge(q95, q10, by="water_year", all=TRUE)
    q95_q10$`Q95-Q10` <- q95_q10$Q95 - q95_q10$Q10
    q95_q10 <- q95_q10[, c("water_year", "Q95-Q10")]
  } else {
    q95_q10 <- data.frame(water_year = annual_means$water_year, `Q95-Q10` = NA)
  }
  
  # Merge all metrics into a single data frame
  # Start with annual means
  all_metrics <- annual_means
  
  # Helper function for safe merging
  safe_merge <- function(df1, df2, by_col) {
    if (nrow(df2) == 0) {
      # If df2 is empty, just return df1
      return(df1)
    }
    return(merge(df1, df2, by=by_col, all.x=TRUE))
  }
  
  # Merge seasonal means
  all_metrics <- safe_merge(all_metrics, winter_means, "water_year")
  all_metrics <- safe_merge(all_metrics, spring_means, "water_year")
  all_metrics <- safe_merge(all_metrics, summer_means, "water_year")
  all_metrics <- safe_merge(all_metrics, fall_means, "water_year")
  
  # Merge percentiles
  all_metrics <- safe_merge(all_metrics, q1, "water_year")
  all_metrics <- safe_merge(all_metrics, q5, "water_year")
  if (nrow(q10) > 0) {
    all_metrics <- safe_merge(all_metrics, q10[, c("water_year", "Q10")], "water_year")
  } else {
    all_metrics$Q10 <- NA
  }
  all_metrics <- safe_merge(all_metrics, q20, "water_year")
  all_metrics <- safe_merge(all_metrics, q25, "water_year")
  all_metrics <- safe_merge(all_metrics, q30, "water_year")
  all_metrics <- safe_merge(all_metrics, q40, "water_year")
  all_metrics <- safe_merge(all_metrics, q50, "water_year")
  all_metrics <- safe_merge(all_metrics, q60, "water_year")
  all_metrics <- safe_merge(all_metrics, q70, "water_year")
  all_metrics <- safe_merge(all_metrics, q75, "water_year")
  all_metrics <- safe_merge(all_metrics, q80, "water_year")
  all_metrics <- safe_merge(all_metrics, q90, "water_year")
  if (nrow(q95) > 0) {
    all_metrics <- safe_merge(all_metrics, q95[, c("water_year", "Q95")], "water_year")
  } else {
    all_metrics$Q95 <- NA
  }
  all_metrics <- safe_merge(all_metrics, q99, "water_year")
  all_metrics <- safe_merge(all_metrics, q95_q10, "water_year")
  
  # Get list of metric columns (all except water_year)
  metric_columns <- setdiff(names(all_metrics), "water_year")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(all_metrics, value_cols = metric_columns, year_col = "water_year")
  
  return(result)
}



analyze_fdc_trends_from_streamflow <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Calculate FDC characteristics by year
  years <- unique(streamflow_data$water_year)
  
  # Check if we have any valid years
  if (length(years) == 0 || all(is.na(years))) {
    # Return empty result with correct structure
    return(data.frame(
      FDCall_slp = NA, FDCall_rho = NA, FDCall_pval = NA, FDCall_mean = NA, FDCall_median = NA,
      FDC90_slp = NA, FDC90th_rho = NA, FDC90th_pval = NA, FDC90th_mean = NA, FDC90th_median = NA,
      FDCmid_slp = NA, FDCmid_rho = NA, FDCmid_pval = NA, FDCmid_mean = NA, FDCmid_median = NA
    ))
  }
  
  # Initialize FDC_byYear data frame
  FDC_byYear <- data.frame(
    water_year = years,
    slp_all = NA,
    slp_90th = NA,
    slp_mid = NA
  )
  
  # For each year, calculate FDC slopes
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    # Need sufficient data points for the year
    if (nrow(year_data) < 30) {
      next
    }
    
    # Remove NA values from Q
    Q_values <- year_data$Q[!is.na(year_data$Q)]
    
    if (length(Q_values) < 30) {
      next
    }
    
    # Sort flows in descending order
    sorted_flows <- sort(Q_values, decreasing = TRUE)
    n <- length(sorted_flows)
    
    # Calculate exceedance probabilities
    exceedance <- (1:n) / (n + 1)
    
    # Create FDC data frame
    fdc <- data.frame(
      exceedance = exceedance,
      flow = sorted_flows
    )
    
    # Calculate slopes for different segments of the FDC
    if (n >= 10) {
      # Use log-transformed flow for better fit
      fdc$log_flow <- log10(fdc$flow + 1e-10)  # Add small constant to handle zeros
      
      # Overall slope
      all_model <- try(lm(log_flow ~ exceedance, data=fdc), silent=TRUE)
      if (!inherits(all_model, "try-error") && !is.na(coef(all_model)[2])) {
        FDC_byYear$slp_all[FDC_byYear$water_year == yr] <- coef(all_model)[2]
      }
      
      # Slope for 90th percentile and above (low flows)
      low_flow_data <- fdc[fdc$exceedance >= 0.9, ]
      if (nrow(low_flow_data) >= 3) {
        low_flow_model <- try(lm(log_flow ~ exceedance, data=low_flow_data), silent=TRUE)
        if (!inherits(low_flow_model, "try-error") && !is.na(coef(low_flow_model)[2])) {
          FDC_byYear$slp_90th[FDC_byYear$water_year == yr] <- coef(low_flow_model)[2]
        }
      }
      
      # Slope for mid-range flows (20th to 80th percentile)
      mid_flow_data <- fdc[fdc$exceedance >= 0.2 & fdc$exceedance <= 0.8, ]
      if (nrow(mid_flow_data) >= 3) {
        mid_flow_model <- try(lm(log_flow ~ exceedance, data=mid_flow_data), silent=TRUE)
        if (!inherits(mid_flow_model, "try-error") && !is.na(coef(mid_flow_model)[2])) {
          FDC_byYear$slp_mid[FDC_byYear$water_year == yr] <- coef(mid_flow_model)[2]
        }
      }
    }
  }
  
  # Use generate_stats for all three FDC metrics
  stats_result <- generate_stats(FDC_byYear, value_cols = c("slp_all", "slp_90th", "slp_mid"), year_col = "water_year")
  
  # Rename columns to match expected output
  names(stats_result) <- gsub("slp_all", "FDCall", names(stats_result))
  names(stats_result) <- gsub("slp_90th", "FDC90th", names(stats_result))
  names(stats_result) <- gsub("slp_mid", "FDCmid", names(stats_result))
  
  # Fix the one inconsistent column name
  names(stats_result) <- gsub("FDC90th_slp", "FDC90_slp", names(stats_result))
  
  # Convert to data frame
  result <- as.data.frame(stats_result)
  
  # Add FDC_byYear as an attribute to the result
  attr(result, "FDC_byYear") <- FDC_byYear
  
  return(result)
}




analyze_flashiness_trends <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Calculate Richards-Baker flashiness index by year
  years <- unique(streamflow_data$water_year)
  
  # Initialize flashiness_byYear data frame
  flashiness_byYear <- data.frame(
    water_year = years,
    RB_index = NA
  )
  
  
  # For each year, calculate the RB flashiness index
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    
    # Need sufficient data points for the year
    if (nrow(year_data) < 30) {
      next
    }
    
    # Sort by day to ensure chronological order
    if ("dowy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$dowy), ]
    }
    
    
    # Calculate RB index: sum of absolute day-to-day changes divided by total flow
    q_values <- year_data$Q
    
    # Check for missing values
    if (sum(is.na(q_values)) > 0) {
      # Skip if too many missing values (more than 20%)
      if (sum(is.na(q_values)) / length(q_values) > 0.2) {
        next
      }
      # Otherwise, interpolate missing values
      q_values <- approx(1:length(q_values), q_values, 1:length(q_values), rule=2)$y
    }
    
    # Calculate absolute day-to-day changes
    q_diff <- abs(diff(q_values))
    
    # Calculate RB index
    rb_index <- sum(q_diff, na.rm=TRUE) / sum(q_values, na.rm=TRUE)
    
    # Store in flashiness_byYear
    flashiness_byYear$RB_index[flashiness_byYear$water_year == yr] <- rb_index
  }
  
  # Use generate_stats to calculate all statistics
  result <- generate_stats(flashiness_byYear, value_cols = "RB_index", year_col = "water_year")
  
  # Rename columns to match expected output
  names(result) <- gsub("RB_index", "flashinessRB", names(result))
  
  # Add flashiness_byYear as an attribute to the result
  attr(result, "flashiness_byYear") <- flashiness_byYear
  
  return(result)
}




analyze_flow_timing_trends <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q", "dowy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Create a data frame to store Julian days when cumulative flow reaches each percentile
  years <- unique(streamflow_data$water_year)
  julday_max <- data.frame(water_year = years)
  
  # Define percentiles
  percentiles <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  
  # Initialize columns for all metrics
  for (p in percentiles) {
    julday_max[[paste0("D", p, "_day")]] <- NA
  }
  julday_max$D25_to_D75 <- NA
  julday_max$Dmax <- NA
  
  # For each year, find the day when cumulative flow reaches each percentile threshold
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    
    # Skip years with insufficient data
    if (nrow(year_data) < 300) {
      next
    }
    
    # Sort by day of year to ensure chronological order
    year_data <- year_data[order(year_data$dowy), ]
    
    # Calculate total annual flow
    total_flow <- sum(year_data$Q, na.rm = TRUE)
    
    # Skip years with zero or NA total flow
    if (total_flow <= 0 || is.na(total_flow)) {
      next
    }
    
    # Calculate cumulative flow for each day
    year_data$cum_flow <- cumsum(year_data$Q)
    
    # Calculate cumulative flow as percentage of total
    year_data$cum_pct <- (year_data$cum_flow / total_flow) * 100
    
    # For each percentile, find the first day when cumulative flow exceeds the threshold
    for (p in percentiles) {
      # Find days where cumulative percentage exceeds the threshold
      above_threshold <- which(year_data$cum_pct >= p)
      
      # If there are days above threshold, take the first one
      if (length(above_threshold) > 0) {
        julday_max[julday_max$water_year == yr, paste0("D", p, "_day")] <- 
          year_data$dowy[above_threshold[1]]
      }
    }
    
    # Calculate D25_to_D75 (days between 25% and 75% cumulative flow)
    # Find days for 25% and 75%
    above_25 <- which(year_data$cum_pct >= 25)
    above_75 <- which(year_data$cum_pct >= 75)
    
    if (length(above_25) > 0 && length(above_75) > 0) {
      day_25 <- year_data$dowy[above_25[1]]
      day_75 <- year_data$dowy[above_75[1]]
      julday_max[julday_max$water_year == yr, "D25_to_D75"] <- day_75 - day_25
    }
    
    # Calculate Dmax (day of maximum discharge)
    max_Q_idx <- which.max(year_data$Q)
    if (length(max_Q_idx) > 0) {
      julday_max[julday_max$water_year == yr, "Dmax"] <- year_data$dowy[max_Q_idx]
    }
  }
  
  # Define which columns to calculate statistics for
  metric_columns <- c(paste0("D", percentiles, "_day"), "D25_to_D75", "Dmax")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(julday_max, value_cols = metric_columns, year_col = "water_year")
  
  # Add julday_max as an attribute to the result
  attr(result, "julday_max") <- julday_max
  
  return(result)
}




calculate_pulse_metrics <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q", "dowy", "month")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Calculate overall 90th and 10th percentiles for entire period
  q90_all <- quantile(streamflow_data$Q, probs = 0.90, na.rm = TRUE)
  q10_all <- quantile(streamflow_data$Q, probs = 0.10, na.rm = TRUE)
  
  # Initialize data frame to store annual pulse metrics
  years <- unique(streamflow_data$water_year)
  pulse_metrics <- data.frame(
    water_year = years,
    n_high_pulses_year = NA,
    n_low_pulses_year = NA,
    n_high_pulses_all = NA,
    n_low_pulses_all = NA,
    dur_high_pulses_year = NA,
    dur_low_pulses_year = NA,
    dur_high_pulses_all = NA,
    dur_low_pulses_all = NA,
    TQmean = NA,
    Flow_Reversals_annual = NA,
    Flow_Reversals_winter = NA,
    Flow_Reversals_spring = NA,
    Flow_Reversals_summer = NA,
    Flow_Reversals_fall = NA
  )
  
  # Function to identify pulses (consecutive days above/below threshold)
  identify_pulses <- function(flow_vector, threshold, above = TRUE) {
    if (above) {
      exceeds <- flow_vector > threshold
    } else {
      exceeds <- flow_vector < threshold
    }
    
    # Handle NAs by treating them as FALSE
    exceeds[is.na(exceeds)] <- FALSE
    
    # Find runs of consecutive TRUE values
    runs <- rle(exceeds)
    
    # Extract pulses (runs where value is TRUE and length >= 1)
    pulse_lengths <- runs$lengths[runs$values == TRUE]
    
    # Return number of pulses and their durations
    if (length(pulse_lengths) > 0) {
      return(list(
        n_pulses = length(pulse_lengths),
        durations = pulse_lengths,
        mean_duration = mean(pulse_lengths)
      ))
    } else {
      return(list(
        n_pulses = 0,
        durations = numeric(0),
        mean_duration = NA
      ))
    }
  }
  
  # Function to count flow reversals
  count_flow_reversals <- function(flow_vector, threshold_pct = 0.02) {
    # Remove NAs
    flow_clean <- flow_vector[!is.na(flow_vector)]
    n <- length(flow_clean)
    
    if (n < 3) return(0)
    
    reversal_count <- 0
    
    for (i in 2:(n-1)) {
      # Calculate changes
      prev_change <- flow_clean[i] - flow_clean[i-1]
      next_change <- flow_clean[i+1] - flow_clean[i]
      
      # Check if change exceeds threshold (2% of current flow)
      threshold <- abs(threshold_pct * flow_clean[i])
      
      # Check for reversal: increasing to decreasing or decreasing to increasing
      if (abs(next_change) > threshold) {
        # Was increasing (or flat), now decreasing
        if (prev_change >= 0 && next_change < 0) {
          reversal_count <- reversal_count + 1
        }
        # Was decreasing (or flat), now increasing
        else if (prev_change <= 0 && next_change > 0) {
          reversal_count <- reversal_count + 1
        }
      }
    }
    
    return(reversal_count)
  }
  
  # Process each year
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    # Sort by day of year to ensure chronological order
    year_data <- year_data[order(year_data$dowy), ]
    
    # Calculate year-specific thresholds
    q90_year <- quantile(year_data$Q, probs = 0.90, na.rm = TRUE)
    q10_year <- quantile(year_data$Q, probs = 0.10, na.rm = TRUE)
    
    # Skip if thresholds can't be calculated
    if (is.na(q90_year) || is.na(q10_year)) {
      next
    }
    
    # Analyze pulses for year-specific thresholds
    high_pulses_year <- identify_pulses(year_data$Q, q90_year, above = TRUE)
    low_pulses_year <- identify_pulses(year_data$Q, q10_year, above = FALSE)
    
    # Analyze pulses for period-of-record thresholds
    high_pulses_all <- identify_pulses(year_data$Q, q90_all, above = TRUE)
    low_pulses_all <- identify_pulses(year_data$Q, q10_all, above = FALSE)
    
    # Calculate TQmean
    annual_mean_flow <- mean(year_data$Q, na.rm = TRUE)
    days_above_mean <- sum(year_data$Q > annual_mean_flow, na.rm = TRUE)
    total_days <- sum(!is.na(year_data$Q))
    tqmean_pct <- (days_above_mean / total_days) * 100
    
    # Calculate annual flow reversals
    annual_reversals <- count_flow_reversals(year_data$Q)
    
    # Calculate seasonal flow reversals
    winter_data <- year_data[year_data$month %in% c(12, 1, 2), ]
    spring_data <- year_data[year_data$month %in% c(3, 4, 5), ]
    summer_data <- year_data[year_data$month %in% c(6, 7, 8), ]
    fall_data <- year_data[year_data$month %in% c(9, 10, 11), ]
    
    winter_reversals <- if(nrow(winter_data) >= 30) count_flow_reversals(winter_data$Q) else NA
    spring_reversals <- if(nrow(spring_data) >= 30) count_flow_reversals(spring_data$Q) else NA
    summer_reversals <- if(nrow(summer_data) >= 30) count_flow_reversals(summer_data$Q) else NA
    fall_reversals <- if(nrow(fall_data) >= 30) count_flow_reversals(fall_data$Q) else NA
    
    # Store results
    idx <- which(pulse_metrics$water_year == yr)
    pulse_metrics$n_high_pulses_year[idx] <- high_pulses_year$n_pulses
    pulse_metrics$n_low_pulses_year[idx] <- low_pulses_year$n_pulses
    pulse_metrics$n_high_pulses_all[idx] <- high_pulses_all$n_pulses
    pulse_metrics$n_low_pulses_all[idx] <- low_pulses_all$n_pulses
    pulse_metrics$dur_high_pulses_year[idx] <- high_pulses_year$mean_duration
    pulse_metrics$dur_low_pulses_year[idx] <- low_pulses_year$mean_duration
    pulse_metrics$dur_high_pulses_all[idx] <- high_pulses_all$mean_duration
    pulse_metrics$dur_low_pulses_all[idx] <- low_pulses_all$mean_duration
    pulse_metrics$TQmean[idx] <- tqmean_pct
    pulse_metrics$Flow_Reversals_annual[idx] <- annual_reversals
    pulse_metrics$Flow_Reversals_winter[idx] <- winter_reversals
    pulse_metrics$Flow_Reversals_spring[idx] <- spring_reversals
    pulse_metrics$Flow_Reversals_summer[idx] <- summer_reversals
    pulse_metrics$Flow_Reversals_fall[idx] <- fall_reversals
  }
  
  # Define which columns to calculate statistics for (all except year)
  metric_columns <- setdiff(names(pulse_metrics), "water_year")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(pulse_metrics, value_cols = metric_columns, year_col = "water_year")
  
  
  # Add pulse_metrics as an attribute to the result
  attr(result, "pulse_metrics") <- pulse_metrics
  
  return(result)
}





analyze_Q_PPT_relationships <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q", "PPT", "month")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Calculate annual totals
  annual_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=streamflow_data, FUN=sum, na.rm=TRUE)
  # Calculate annual runoff ratio, handling cases where PPT is zero or very small
  annual_totals$annual_runoff_ratio <- ifelse(annual_totals$PPT > 0.001, 
                                              annual_totals$Q / annual_totals$PPT, 
                                              NA)
  
  # Calculate seasonal totals and ratios
  # Winter (December, January, February)
  winter <- streamflow_data[streamflow_data$month %in% c(12, 1, 2), ]
  winter_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=winter, FUN=sum, na.rm=TRUE)
  winter_totals$winter_runoff_ratio <- ifelse(winter_totals$PPT > 0.001, 
                                              winter_totals$Q / winter_totals$PPT, 
                                              NA)
  
  # Spring (March, April, May)
  spring <- streamflow_data[streamflow_data$month %in% c(3, 4, 5), ]
  spring_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=spring, FUN=sum, na.rm=TRUE)
  spring_totals$spring_runoff_ratio <- ifelse(spring_totals$PPT > 0.001, 
                                              spring_totals$Q / spring_totals$PPT, 
                                              NA)
  
  # Summer (June, July, August)
  summer <- streamflow_data[streamflow_data$month %in% c(6, 7, 8), ]
  summer_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=summer, FUN=sum, na.rm=TRUE)
  summer_totals$summer_runoff_ratio <- ifelse(summer_totals$PPT > 0.001, 
                                              summer_totals$Q / summer_totals$PPT, 
                                              NA)
  
  # Fall (September, October, November)
  fall <- streamflow_data[streamflow_data$month %in% c(9, 10, 11), ]
  fall_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=fall, FUN=sum, na.rm=TRUE)
  fall_totals$fall_runoff_ratio <- ifelse(fall_totals$PPT > 0.001, 
                                          fall_totals$Q / fall_totals$PPT, 
                                          NA)
  
  # Combine all runoff ratios into a single dataframe by year
  # Before merging, ensure all seasonal totals have matching structure:
  all_years <- unique(c(annual_totals$water_year, winter_totals$water_year, 
                        spring_totals$water_year, summer_totals$water_year, 
                        fall_totals$water_year))
  
  # Create base data frame with all years
  all_ratios <- data.frame(water_year = all_years)
  
  # Merge annual ratios first
  all_ratios <- merge(all_ratios, annual_totals[, c("water_year", "annual_runoff_ratio")], 
                      by = "water_year", all.x = TRUE)
  
  # Then merge seasonal ratios
  all_ratios <- merge(all_ratios, winter_totals[, c("water_year", "winter_runoff_ratio")], 
                      by = "water_year", all.x = TRUE)
  all_ratios <- merge(all_ratios, spring_totals[, c("water_year", "spring_runoff_ratio")], 
                      by = "water_year", all.x = TRUE)
  all_ratios <- merge(all_ratios, summer_totals[, c("water_year", "summer_runoff_ratio")], 
                      by = "water_year", all.x = TRUE)
  all_ratios <- merge(all_ratios, fall_totals[, c("water_year", "fall_runoff_ratio")], 
                      by = "water_year", all.x = TRUE)
  
  
  # Define which columns to calculate statistics for (all except year)
  metric_columns <- setdiff(names(all_ratios), "water_year")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(all_ratios, value_cols = metric_columns, year_col = "water_year")
  
  # Add the annual data as an attribute for reference
  attr(result, "runoff_ratios_by_year") <- list(
    annual = annual_totals,
    winter = winter_totals,
    spring = spring_totals,
    summer = summer_totals,
    fall = fall_totals
  )
  
  return(result)
}





analyze_baseflow_indices <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q", "doy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Function to apply Eckhardt filter
  eckhardt_filter <- function(Q, BFImax = 0.8, a = 0.98) {
    n <- length(Q)
    baseflow <- numeric(n)
    baseflow[1] <- Q[1] * BFImax  # Initialize with fraction of first flow
    
    for (i in 2:n) {
      if (!is.na(Q[i]) && !is.na(Q[i-1]) && !is.na(baseflow[i-1])) {
        # Eckhardt filter equation
        numerator <- (1 - BFImax) * a * baseflow[i-1] + (1 - a) * BFImax * Q[i]
        denominator <- 1 - a * BFImax
        baseflow[i] <- numerator / denominator
        
        # Baseflow cannot exceed total flow
        baseflow[i] <- min(baseflow[i], Q[i])
        # Baseflow cannot be negative
        baseflow[i] <- max(baseflow[i], 0)
      } else {
        baseflow[i] <- NA
      }
    }
    
    return(baseflow)
  }
  
  # Function to apply Lyne-Hollick filter
  lyne_hollick_filter <- function(Q, alpha = 0.925, passes = 2) {
    n <- length(Q)
    
    # Apply filter in forward and backward directions for specified passes
    quickflow <- Q
    
    for (pass in 1:passes) {
      # Forward pass
      qf_forward <- numeric(n)
      qf_forward[1] <- 0
      
      for (i in 2:n) {
        if (!is.na(Q[i]) && !is.na(Q[i-1]) && !is.na(qf_forward[i-1])) {
          qf_forward[i] <- alpha * qf_forward[i-1] + ((1 + alpha) / 2) * (Q[i] - Q[i-1])
          # Constrain quickflow
          qf_forward[i] <- max(0, min(qf_forward[i], Q[i]))
        } else {
          qf_forward[i] <- NA
        }
      }
      
      # Backward pass
      qf_backward <- numeric(n)
      qf_backward[n] <- qf_forward[n]
      
      for (i in (n-1):1) {
        if (!is.na(qf_forward[i]) && !is.na(qf_forward[i+1]) && !is.na(qf_backward[i+1])) {
          qf_backward[i] <- alpha * qf_backward[i+1] + ((1 + alpha) / 2) * (qf_forward[i] - qf_forward[i+1])
          # Constrain quickflow
          qf_backward[i] <- max(0, min(qf_backward[i], Q[i]))
        } else {
          qf_backward[i] <- NA
        }
      }
      
      quickflow <- qf_backward
    }
    
    # Calculate baseflow as total flow minus quickflow
    baseflow <- Q - quickflow
    baseflow[baseflow < 0] <- 0
    baseflow[is.na(Q)] <- NA
    
    return(baseflow)
  }
  
  # Initialize data frame to store annual BFI values
  years <- unique(streamflow_data$water_year)
  bfi_by_year <- data.frame(
    water_year = years,
    BFI_Eckhardt = NA,
    BFI_LyneHollick = NA
  )
  
  # Process each year
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    # Skip years with insufficient data
    if (nrow(year_data) < 250) {
      next
    }
    
    # Sort by day of year to ensure chronological order
    year_data <- year_data[order(year_data$dowy), ]
    
    # Get streamflow values
    Q <- year_data$Q
    
    # Skip if too many missing values
    if (sum(is.na(Q)) / length(Q) > 0.2) {
      next
    }
    
    # Apply Eckhardt filter
    baseflow_eckhardt <- eckhardt_filter(Q)
    
    # Apply Lyne-Hollick filter
    baseflow_lyne <- lyne_hollick_filter(Q)
    
    # Calculate annual totals and BFI
    total_flow <- sum(Q, na.rm = TRUE)
    
    if (total_flow > 0) {
      # Eckhardt BFI
      total_baseflow_eckhardt <- sum(baseflow_eckhardt, na.rm = TRUE)
      bfi_eckhardt <- total_baseflow_eckhardt / total_flow
      bfi_by_year$BFI_Eckhardt[bfi_by_year$water_year == yr] <- bfi_eckhardt
      
      # Lyne-Hollick BFI
      total_baseflow_lyne <- sum(baseflow_lyne, na.rm = TRUE)
      bfi_lyne <- total_baseflow_lyne / total_flow
      bfi_by_year$BFI_LyneHollick[bfi_by_year$water_year == yr] <- bfi_lyne
    }
  }
  
  # Define which columns to calculate statistics for (all except year)
  metric_columns <- setdiff(names(bfi_by_year), "water_year")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(bfi_by_year, value_cols = metric_columns, year_col = "water_year")
  
  
  # Add bfi_by_year as an attribute to the result
  attr(result, "bfi_by_year") <- bfi_by_year
  
  return(result)
}



analyze_recession_parameters <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("water_year", "Q", "dowy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Define signatures that will have standard stats
  signatures_with_stats <- c("log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity")
  
  # Define seasonality signatures (these will only have values, no stats)
  seasonality_signatures <- c("log_a_seasonality_amplitude_all", "log_a_seasonality_minimum_all",
                              "log_a_seasonality_amplitude_first_half", "log_a_seasonality_minimum_first_half",
                              "log_a_seasonality_amplitude_last_half", "log_a_seasonality_minimum_last_half")
  
  # Function to identify recession events (unchanged)
  identify_recession_events <- function(Q_vector, min_length = 5) {
    n <- length(Q_vector)
    if (n < min_length + 1) return(list())
    
    # Calculate dQ/dt (using forward difference)
    dQ_dt <- c(diff(Q_vector), NA)
    
    # Initialize list to store recession events
    recession_events <- list()
    
    # Track current recession
    in_recession <- FALSE
    start_idx <- NA
    
    for (i in 1:(n-min_length)) {
      # Check if we have valid data
      if (is.na(Q_vector[i]) || is.na(dQ_dt[i])) {
        if (in_recession) {
          # End current recession if we hit NA
          if (!is.na(start_idx) && (i - start_idx) >= min_length) {
            recession_events[[length(recession_events) + 1]] <- list(
              start = start_idx,
              end = i - 1,
              indices = start_idx:(i-1)
            )
          }
          in_recession <- FALSE
          start_idx <- NA
        }
        next
      }
      
      # Check for monotonic decrease in both Q and |dQ/dt|
      if (i < n - 1) {
        # Need at least min_length consecutive days
        is_recession <- TRUE
        
        # Check next min_length days
        for (j in 0:(min_length-2)) {
          if (i+j+1 > n || is.na(Q_vector[i+j]) || is.na(Q_vector[i+j+1]) || 
              is.na(dQ_dt[i+j]) || is.na(dQ_dt[i+j+1])) {
            is_recession <- FALSE
            break
          }
          
          # Check if Q is decreasing
          if (Q_vector[i+j+1] >= Q_vector[i+j]) {
            is_recession <- FALSE
            break
          }
          
          # Check if |dQ/dt| is decreasing (becoming less negative)
          if (abs(dQ_dt[i+j+1]) >= abs(dQ_dt[i+j])) {
            is_recession <- FALSE
            break
          }
        }
        
        if (is_recession && !in_recession) {
          # Start new recession
          in_recession <- TRUE
          start_idx <- i
        } else if (!is_recession && in_recession) {
          # End current recession
          if (!is.na(start_idx) && (i - start_idx) >= min_length) {
            recession_events[[length(recession_events) + 1]] <- list(
              start = start_idx,
              end = i - 1,
              indices = start_idx:(i-1)
            )
          }
          in_recession <- FALSE
          start_idx <- NA
        }
      }
    }
    
    # Check if we ended in a recession
    if (in_recession && !is.na(start_idx) && (n - start_idx) >= min_length) {
      recession_events[[length(recession_events) + 1]] <- list(
        start = start_idx,
        end = n,
        indices = start_idx:n
      )
    }
    
    return(recession_events)
  }
  
  # Function to fit recession parameters for a single event (now returning log_a)
  fit_recession_event <- function(Q_values, remove_first_day = TRUE) {
    if (remove_first_day && length(Q_values) > 1) {
      Q_values <- Q_values[-1]
    }
    
    n <- length(Q_values)
    if (n < 3) return(list(log_a = NA, b = NA))
    
    # Calculate -dQ/dt
    dQ_dt <- -diff(Q_values)
    Q_subset <- Q_values[-length(Q_values)]
    
    # Remove any non-positive values for log transformation
    valid_idx <- which(Q_subset > 0 & dQ_dt > 0)
    if (length(valid_idx) < 2) return(list(log_a = NA, b = NA))
    
    Q_valid <- Q_subset[valid_idx]
    dQ_dt_valid <- dQ_dt[valid_idx]
    
    # Fit in log-log space: log(-dQ/dt) = log(a) + b*log(Q)
    tryCatch({
      fit <- lm(log(dQ_dt_valid) ~ log(Q_valid))
      b_est <- coef(fit)[2]
      log_a_est <- coef(fit)[1]  # This is already log(a)
      
      return(list(log_a = log_a_est, b = b_est))
    }, error = function(e) {
      return(list(log_a = NA, b = NA))
    })
  }
  
  # Function to fit sinusoidal model to log(a) values
  fit_sinusoidal_model <- function(doy_values, log_a_values) {
    # Remove NA values
    valid_idx <- which(!is.na(log_a_values) & !is.na(doy_values))
    if (length(valid_idx) < 10) {
      return(list(amplitude = NA, minimum_doy = NA))
    }
    
    doy_clean <- doy_values[valid_idx]
    log_a_clean <- log_a_values[valid_idx]
    
    # Fit sinusoidal model: log(a) = A * sin(2π/365 * (doy - φ)) + C
    tryCatch({
      # Create design matrix
      X <- cbind(
        sin(2 * pi * doy_clean / 365),
        cos(2 * pi * doy_clean / 365),
        1  # intercept
      )
      
      # Fit linear model
      fit <- lm(log_a_clean ~ X - 1)  # -1 because X already includes intercept
      
      B1 <- coef(fit)[1]
      B2 <- coef(fit)[2]
      C <- coef(fit)[3]
      
      # Calculate amplitude and phase
      amplitude <- sqrt(B1^2 + B2^2)
      
      # Calculate phase (in days)
      phase_rad <- atan2(-B2, B1)
      phase_days <- phase_rad * 365 / (2 * pi)
      
      # Ensure phase is between 0 and 365
      if (phase_days < 0) phase_days <- phase_days + 365
      
      # Minimum occurs at phase + 273.75 days
      minimum_doy <- phase_days + 273.75
      if (minimum_doy > 365) minimum_doy <- minimum_doy - 365
      
      return(list(amplitude = amplitude, minimum_doy = minimum_doy))
      
    }, error = function(e) {
      return(list(amplitude = NA, minimum_doy = NA))
    })
  }
  
  # Initialize storage for annual metrics
  years <- unique(streamflow_data$water_year)
  annual_metrics <- data.frame(
    water_year = years,
    log_a_pointcloud = NA,
    log_a_events = NA,
    b_pointcloud = NA,
    b_events = NA,
    concavity = NA
  )
  
  # Store all recession events with their timing
  all_recession_events <- list()
  
  # Process entire dataset to identify all recession events
  # Sort by year and doy
  streamflow_data <- streamflow_data[order(streamflow_data$water_year, streamflow_data$dowy), ]
  
  # Process each year
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    
    # Sort by day of year
    year_data <- year_data[order(year_data$dowy), ]
    
    # Identify recession events
    recession_events <- identify_recession_events(year_data$Q)
    
    # Collect all recession data for point cloud analysis
    all_Q <- numeric()
    all_dQ_dt <- numeric()
    
    # Store individual event parameters
    event_log_a_values <- numeric()
    event_b_values <- numeric()
    event_concavities <- numeric()
    
    # Process each recession event
    for (event in recession_events) {
      Q_event <- year_data$Q[event$indices]
      
      # Get the middle day of the recession event for timing
      mid_idx <- event$indices[ceiling(length(event$indices)/2)]
      event_dowy <- year_data$dowy[mid_idx]
      
      # Fit parameters for this event
      event_params <- fit_recession_event(Q_event, remove_first_day = TRUE)
      
      if (!is.na(event_params$log_a) && !is.na(event_params$b)) {
        event_log_a_values <- c(event_log_a_values, event_params$log_a)
        event_b_values <- c(event_b_values, event_params$b)
        
        # Store event with its timing
        all_recession_events[[length(all_recession_events) + 1]] <- list(
          water_year = yr,
          dowy = event_dowy,
          log_a = event_params$log_a,
          b = event_params$b
        )
        
        # Calculate concavity (difference in b between first and second half)
        if (length(Q_event) >= 6) {  # Need at least 6 points
          mid_point <- floor(length(Q_event) / 2)
          first_half <- Q_event[1:mid_point]
          second_half <- Q_event[mid_point:length(Q_event)]
          
          params_first <- fit_recession_event(first_half, remove_first_day = FALSE)
          params_second <- fit_recession_event(second_half, remove_first_day = FALSE)
          
          if (!is.na(params_first$b) && !is.na(params_second$b)) {
            concavity <- params_second$b - params_first$b
            event_concavities <- c(event_concavities, concavity)
          }
        }
        
        # Add to point cloud data (removing first day)
        if (length(Q_event) > 1) {
          Q_subset <- Q_event[-1]
          dQ_subset <- -diff(Q_event[-1])
          
          valid_idx <- which(Q_subset[-length(Q_subset)] > 0 & dQ_subset > 0)
          if (length(valid_idx) > 0) {
            all_Q <- c(all_Q, Q_subset[-length(Q_subset)][valid_idx])
            all_dQ_dt <- c(all_dQ_dt, dQ_subset[valid_idx])
          }
        }
      }
    }
    
    # Calculate median b for use in calculating log(a)
    median_b <- median(event_b_values, na.rm = TRUE)
    
    # Point cloud analysis
    if (length(all_Q) > 10) {
      tryCatch({
        # Fit b using point cloud
        pc_fit <- lm(log(all_dQ_dt) ~ log(all_Q))
        b_pointcloud <- coef(pc_fit)[2]
        
        # Calculate log(a) using median b
        # For each point: log(a) = log(-dQ/dt) - b*log(Q)
        log_a_values_pc <- log(all_dQ_dt) - median_b * log(all_Q)
        log_a_pointcloud <- median(log_a_values_pc, na.rm = TRUE)
        
        annual_metrics$b_pointcloud[annual_metrics$water_year == yr] <- b_pointcloud
        annual_metrics$log_a_pointcloud[annual_metrics$water_year == yr] <- log_a_pointcloud
      }, error = function(e) {
        # Leave as NA
      })
    }
    
    # Individual events analysis
    if (length(event_log_a_values) > 0) {
      # Calculate log(a) for events using median b
      log_a_events_recalc <- numeric()
      for (i in 1:length(recession_events)) {
        if (i <= length(event_b_values) && !is.na(event_b_values[i])) {
          event <- recession_events[[i]]
          Q_event <- year_data$Q[event$indices]
          if (length(Q_event) > 1) {
            Q_subset <- Q_event[-1]  # Remove first day
            dQ_subset <- -diff(Q_event[-1])
            valid_idx <- which(Q_subset[-length(Q_subset)] > 0 & dQ_subset > 0)
            if (length(valid_idx) > 0) {
              log_a_vals <- log(dQ_subset[valid_idx]) - median_b * log(Q_subset[-length(Q_subset)][valid_idx])
              log_a_events_recalc <- c(log_a_events_recalc, median(log_a_vals, na.rm = TRUE))
            }
          }
        }
      }
      
      if (length(log_a_events_recalc) > 0) {
        annual_metrics$log_a_events[annual_metrics$water_year == yr] <- median(log_a_events_recalc, na.rm = TRUE)
      }
      
      annual_metrics$b_events[annual_metrics$water_year == yr] <- median(event_b_values, na.rm = TRUE)
    }
    
    # Concavity
    if (length(event_concavities) > 0) {
      annual_metrics$concavity[annual_metrics$water_year == yr] <- mean(event_concavities, na.rm = TRUE)
    }
  }
  
  # Check if we have enough data overall
  total_valid_years <- sum(!is.na(annual_metrics$b_events))
  if (total_valid_years < 3) {
    # Not enough data, return all NAs
    # Initialize results data frame with correct structure
    n_cols <- length(signatures_with_stats) * 5 + length(seasonality_signatures)
    result <- data.frame(matrix(NA, nrow=1, ncol=n_cols))
    
    # Set column names
    col_names <- c()
    for (sig in signatures_with_stats) {
      for (stat in c("slp", "rho", "pval", "mean", "median")) {
        col_names <- c(col_names, paste0(sig, "_", stat))
      }
    }
    col_names <- c(col_names, seasonality_signatures)
    colnames(result) <- col_names
    
    return(result)
  }
  
  # Use generate_stats to calculate statistics for signatures with trends
  stats_result <- generate_stats(annual_metrics, 
                                 value_cols = signatures_with_stats, 
                                 year_col = "water_year")
  
  
  # Initialize the final result with stats_result
  result <- stats_result
  
  # Add seasonality signatures (which don't have stats)
  for (sig in seasonality_signatures) {
    result[[sig]] <- NA
  }
  
  # Calculate seasonality of recession parameter log(a)
  if (length(all_recession_events) >= 10) {
    # Extract DOY and log(a) values from all events
    event_dowys <- sapply(all_recession_events, function(x) x$dowy)
    event_log_a_values <- sapply(all_recession_events, function(x) x$log_a)
    event_water_years <- sapply(all_recession_events, function(x) x$water_year)
    
    # Fit sinusoidal model to all data
    seasonality_all <- fit_sinusoidal_model(event_dowys, event_log_a_values)
    result$log_a_seasonality_amplitude_all <- seasonality_all$amplitude
    result$log_a_seasonality_minimum_all <- seasonality_all$minimum_doy
    
    # Split into first and last half of years
    median_water_year <- median(unique(event_water_years))
    first_half_idx <- which(event_water_years <= median_water_year)
    last_half_idx <- which(event_water_years > median_water_year)
    
    # Fit sinusoidal model to first half
    if (length(first_half_idx) >= 10) {
      seasonality_first <- fit_sinusoidal_model(event_dowys[first_half_idx], 
                                                event_log_a_values[first_half_idx])
      result$log_a_seasonality_amplitude_first_half <- seasonality_first$amplitude
      result$log_a_seasonality_minimum_first_half <- seasonality_first$minimum_doy
    }
    
    # Fit sinusoidal model to last half
    if (length(last_half_idx) >= 10) {
      seasonality_last <- fit_sinusoidal_model(event_dowys[last_half_idx], 
                                               event_log_a_values[last_half_idx])
      result$log_a_seasonality_amplitude_last_half <- seasonality_last$amplitude
      result$log_a_seasonality_minimum_last_half <- seasonality_last$minimum_doy
    }
  }
  
  # Add annual_metrics as an attribute
  attr(result, "recession_metrics_by_water_year") <- annual_metrics
  attr(result, "recession_events") <- all_recession_events
  
  return(result)
}






#########################################################################################
#########################################################################################
# SECTION 3: HELPER FUNCITONS FOR READING AND PLOTTING DATA

#########################################################################################
#########################################################################################

library(arrow)
library(ggplot2)
library(data.table)
library(viridis)
library(lubridate)

# Function to read parquet data for specific watersheds
read_streamflow_data <- function(output_dir, gage_ids = NULL, max_gages_to_plot = 10) {
  # Get list of parquet files
  parquet_files <- list.files(output_dir, 
                              pattern = "daily_data_chunk_.*\\.parquet", 
                              full.names = TRUE)
  
  if (length(parquet_files) == 0) {
    stop("No parquet files found in ", output_dir)
  }
  
  # Read and combine data
  all_data <- data.table()
  
  for (file in parquet_files) {
    cat("Reading", basename(file), "...\n")
    chunk <- read_parquet(file)
    
    # Filter for specific gage_ids if provided
    if (!is.null(gage_ids)) {
      chunk <- chunk[chunk$gage_id %in% gage_ids, ]
    }
    
    if (nrow(chunk) > 0) {
      all_data <- rbind(all_data, chunk, fill = TRUE)
    }
  }
  
  # If no specific gages requested, sample some for plotting
  if (is.null(gage_ids) && length(unique(all_data$gage_id)) > max_gages_to_plot) {
    cat("Sampling", max_gages_to_plot, "gages for plotting...\n")
    gage_ids <- sample(unique(all_data$gage_id), max_gages_to_plot)
    all_data <- all_data[gage_id %in% gage_ids]
  }
  
  return(all_data)
}

# Function to plot time series
plot_streamflow_timeseries <- function(streamflow_data, 
                                       title = "Streamflow Time Series by Gage",
                                       date_range = NULL,
                                       log_scale = FALSE) {
  
  # Ensure Date is in Date format
  streamflow_data$Date <- as.Date(streamflow_data$Date)
  
  # Filter by date range if specified
  if (!is.null(date_range)) {
    streamflow_data <- streamflow_data[Date >= date_range[1] & Date <= date_range[2]]
  }
  
  # Create the plot
  p <- ggplot(streamflow_data, aes(x = Date, y = Q, color = gage_id)) +
    geom_line(alpha = 0.7, size = 0.8) +
    scale_color_viridis_d(name = "Gage ID") +
    labs(title = title,
         x = "Date",
         y = "Streamflow (mm/day)") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10(labels = scales::comma) +
      annotation_logticks(sides = "l")
  }
  
  return(p)
}

# Function to create faceted plot for many gages
plot_streamflow_faceted <- function(streamflow_data, 
                                    ncol = 2,
                                    date_range = NULL,
                                    scales = "free_y") {
  
  # Ensure Date is in Date format
  streamflow_data$Date <- as.Date(streamflow_data$Date)
  
  # Filter by date range if specified
  if (!is.null(date_range)) {
    streamflow_data <- streamflow_data[Date >= date_range[1] & Date <= date_range[2]]
  }
  
  # Create faceted plot
  p <- ggplot(streamflow_data, aes(x = Date, y = Q)) +
    geom_line(color = "steelblue", alpha = 0.8) +
    facet_wrap(~ gage_id, ncol = ncol, scales = scales) +
    labs(title = "Streamflow Time Series by Gage (Faceted)",
         x = "Date",
         y = "Streamflow (mm/day)") +
    theme_minimal() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  return(p)
}

# Function to plot annual patterns
plot_annual_patterns <- function(streamflow_data, 
                                 aggregate_fun = "mean") {
  
  # Calculate day of year average across all years
  if (aggregate_fun == "mean") {
    daily_avg <- streamflow_data[, .(Q_avg = mean(Q, na.rm = TRUE)), 
                                 by = .(gage_id, dowy)]
  } else if (aggregate_fun == "median") {
    daily_avg <- streamflow_data[, .(Q_avg = median(Q, na.rm = TRUE)), 
                                 by = .(gage_id, dowy)]
  }
  
  # Create plot
  p <- ggplot(daily_avg, aes(x = dowy, y = Q_avg, color = gage_id)) +
    geom_line(size = 1.2, alpha = 0.8) +
    scale_color_viridis_d(name = "Gage ID") +
    labs(title = paste("Annual Streamflow Pattern (", aggregate_fun, ")", sep = ""),
         x = "Day of Water Year",
         y = paste(aggregate_fun, "Streamflow (mm/day)")) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"))
  
  return(p)
}





calculate_water_year_info <- function(dates) {
  # For Northern Hemisphere: water year starts October 1
  # TODO: calculate water year for alternative regions
  
  years <- year(dates)
  months <- month(dates)
  
  # Water year is the calendar year + 1 if month >= 10
  water_years <- ifelse(months >= 10, years + 1, years)
  
  # Calculate day of water year (dowy)
  # October 1 is day 1 of the water year
  wy_start <- as.Date(paste0(ifelse(months >= 10, years, years - 1), "-10-01"))
  dowy <- as.numeric(dates - wy_start + 1)
  
  return(list(water_year = water_years, dowy = dowy))
}






