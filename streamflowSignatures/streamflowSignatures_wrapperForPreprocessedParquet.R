


####################################################################################################################
# Final user-facing function

# Function to process streamflow signatures from parquet files (long format)
process_signatures_from_parquet <- function(
    parquet_file_path = "combined_streamflow_output/combined_daily_data.parquet",
    metadata_file_path = "combined_streamflow_output/combined_watershed_metadata.csv",
    output_file = "processed_signatures_from_parquet.csv",
    min_Q_value_and_days = c(0.0001, 30),
    min_num_years = 20,
    min_frac_good_data = 0.9  
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
    gage_id <- as.character(gage_id)
    
    # First try exact match
    meta <- metadata_lookup[gage_id == gage_id]
    if (nrow(meta) > 0) return(meta[1])
    
    # Try adding leading zeros (up to 4)
    for (num_zeros in 1:4) {
      padded_id <- sprintf(paste0("%0", nchar(gage_id) + num_zeros, "d"), 
                           as.numeric(gage_id))
      meta <- metadata_lookup[gage_id == padded_id]
      if (nrow(meta) > 0) {
        cat("  Found match for", gage_id, "as", padded_id, "\n")
        return(meta[1])
      }
    }
    
    # Try removing leading zeros
    stripped_id <- gsub("^0+", "", gage_id)
    if (stripped_id != gage_id && stripped_id != "") {
      meta <- metadata_lookup[gage_id == stripped_id]
      if (nrow(meta) > 0) {
        cat("  Found match for", gage_id, "as", stripped_id, "\n")
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
