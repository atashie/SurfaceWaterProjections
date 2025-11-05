#########################################################################################
#########################################################################################
# SECTION 0: PROCESS RAW DATA VIA API AND SAVE AS PARQUET FILE
#            THESE FUNCTIONS ARE ONLY HERE FOR REFERENCE
#########################################################################################
#########################################################################################

#########################################################################################

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


# ORIGINAL WRAPPER FOR INGESTING RAW DATA VIA API AND OUTPUTING SIGNATURE TRENDS AND METRICS DIRECTLY
process_gages_rawData <- function(gages_df, gage_type, min_num_years, start_date, end_date, 
                                  min_Q_value_and_days, basinAt_NorAm_polys, HB_dt, 
                                  upstream_hydrobasins, output_file) {
  
  # Check if output file exists; if so, read it with gage_id forced to character.
  if (file.exists(output_file)) {
    summary_output <- fread(output_file, colClasses = list(character = "gage_id"))
    # Also force conversion in case some rows were stored as integer64.
    summary_output[, gage_id := as.character(gage_id)]
    cat("Loaded existing summary data with", nrow(summary_output), "rows\n")
  } else {
    summary_output <- data.table(
      gage_id = character(),
      latitude = numeric(),
      longitude = numeric(),
      basin_area = numeric(),
      gage_type = character()
    )
    cat("Created new summary data table\n")
    fwrite(summary_output, output_file)
  }
  
  # Process each gage
  for (i in 1:nrow(gages_df)) {
    current_gage <- gages_df[i, ]
    
    # Extract gage ID and coordinate info based on gage type
    if (gage_type == "USGS") {
      gage_id <- as.character(current_gage$STAID)
      latitude <- as.numeric(current_gage$LAT_GAGE)
      longitude <- as.numeric(current_gage$LNG_GAGE)
      basin_area <- as.numeric(current_gage$DRAIN_SQKM)
    } else if (gage_type %in% c("Canada", "CANADIAN")) {
      # Rename Canadian coordinate columns to match USGS convention
      names(current_gage)[names(current_gage) == 'LATITUDE'] <- 'LAT_GAGE'
      names(current_gage)[names(current_gage) == 'LONGITUDE'] <- 'LNG_GAGE'
      gage_id <- as.character(current_gage$STATION_NUMBER)
      latitude <- as.numeric(current_gage$LAT_GAGE)
      longitude <- as.numeric(current_gage$LNG_GAGE)
      basin_area <- NA
    } else {
      stop("Unsupported gage type")
    }
    
    if (gage_id %in% summary_output$gage_id) {
      cat("Skipping gage", gage_id, "as it's already in the output\n")
      next
    }
    
    cat("Processing gage", gage_id, "(", i, "of", nrow(gages_df), ")\n")
    
    tryCatch({
      streamflow_data <- generate_streamflow_dt(current_gage, gage_type, 
                                                min_num_years, start_date, end_date)
      if (is.null(streamflow_data) || identical(streamflow_data, NA) ||
          (is.data.frame(streamflow_data) && nrow(streamflow_data) == 0)) {
        cat("No valid streamflow data for gage", gage_id, "\n")
        next
      }
      if (!is.data.frame(streamflow_data)) {
        cat("Invalid streamflow data format for gage", gage_id, "\n")
        next
      }
      
      # Updated required columns to include water_year and dowy
      required_cols <- c("year", "Q", "doy", "month", "water_year", "dowy")
      if (!all(required_cols %in% colnames(streamflow_data))) {
        if (!"year" %in% colnames(streamflow_data) && "Date" %in% colnames(streamflow_data)) {
          streamflow_data$year <- year(streamflow_data$Date)
        }
        if (!"doy" %in% colnames(streamflow_data) && "Date" %in% colnames(streamflow_data)) {
          streamflow_data$doy <- yday(streamflow_data$Date)
        }
        if (!"month" %in% colnames(streamflow_data) && "Date" %in% colnames(streamflow_data)) {
          streamflow_data$month <- month(streamflow_data$Date)
        }
        # Add water year information if missing
        if (!all(c("water_year", "dowy") %in% colnames(streamflow_data)) && "Date" %in% colnames(streamflow_data)) {
          wy_info <- calculate_water_year_info(streamflow_data$Date)
          streamflow_data$water_year <- wy_info$water_year
          streamflow_data$dowy <- wy_info$dowy
        }
        if (!all(required_cols %in% colnames(streamflow_data))) {
          cat("Missing required columns in streamflow data for gage", gage_id, "\n")
          next
        }
      }
      
      # Apply min_Q_value_and_days filter
      years_to_use <- NULL
      for (this_year in unique(streamflow_data$year)) {
        test_year <- subset(streamflow_data, year == this_year)
        # Ensure Q is numeric for comparison
        if(!is.numeric(test_year$Q)) test_year$Q <- as.numeric(test_year$Q)
        
        # Handle potential NAs in Q before comparison
        valid_q_values <- test_year$Q[!is.na(test_year$Q)]
        nonzero_rows <- which(valid_q_values > min_Q_value_and_days[1])
        
        if (length(nonzero_rows) >= min_Q_value_and_days[2]) {
          years_to_use <- c(years_to_use, this_year)
        }
      }
      
      if (length(years_to_use) < min_num_years) {
        cat("Insufficient years with valid data for gage", gage_id, "\n")
        next
      }
      
      streamflow_data_filtered <- streamflow_data[streamflow_data$year %in% years_to_use, ]
      
      if (nrow(streamflow_data_filtered) == 0) {
        cat("No data remaining after year filtering for gage", gage_id, "\n")
        next
      }
      
      # [Rest of the function remains the same...]
      # Find upstream basins
      upstream_basins <- NULL
      num_upstream_basins <- NA
      tryCatch({
        upstream_basins <- find_upstream_hydrobasins(
          current_gage = current_gage,
          basinAt_NorAm_polys = basinAt_NorAm_polys,
          HB_dt = HB_dt,
          upstream_hydrobasins = upstream_hydrobasins,
          save_path = file.path(dirname(output_file), "upstream_hydrobasins.RData")
        )
        if (!is.null(upstream_basins)) {
          num_upstream_basins <- length(upstream_basins)
        }
      }, error = function(e) {
        cat("Error finding upstream basins for gage", gage_id, ":", e$message, "\n")
      })
      
      # Calculate all streamflow metrics
      cat("Calculating metrics for gage", gage_id, "...\n")
      
      # Initialize empty metric results in case of errors
      metrics_flow_vols <- NULL
      metrics_fdc_trends <- NULL
      metrics_flashiness <- NULL
      metrics_flow_timing <- NULL
      metrics_pulses <- NULL
      metrics_baseflow <- NULL
      metrics_recession <- NULL
      
      # Calculate each metric with error handling
      tryCatch({
        metrics_flow_vols <- calculate_flow_vols_by_year(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating flow volumes for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_fdc_trends <- analyze_fdc_trends_from_streamflow(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating FDC trends for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_flashiness <- analyze_flashiness_trends(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating flashiness for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_flow_timing <- analyze_flow_timing_trends(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating flow timing for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_pulses <- calculate_pulse_metrics(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating pulse metrics for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_baseflow <- analyze_baseflow_indices(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating baseflow indices for gage", gage_id, ":", e$message, "\n")
      })
      
      tryCatch({
        metrics_recession <- analyze_recession_parameters(streamflow_data_filtered)
      }, error = function(e) {
        cat("Error calculating recession indices for gage", gage_id, ":", e$message, "\n")
      })
      
      # Create a row for this gage with base information
      gage_row <- data.table(
        gage_id = gage_id,
        latitude = latitude,
        longitude = longitude,
        basin_area = basin_area,
        gage_type = gage_type,
        num_years = length(years_to_use),
        start_year = min(years_to_use),
        end_year = max(years_to_use),
        num_upstream_basins = num_upstream_basins
      )
      
      # Combine base gage info with calculated metrics
      # Only add metrics that were successfully calculated
      if (!is.null(metrics_flow_vols)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_flow_vols))
      }
      if (!is.null(metrics_fdc_trends)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_fdc_trends))
      }
      if (!is.null(metrics_flashiness)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_flashiness))
      }
      if (!is.null(metrics_flow_timing)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_flow_timing))
      }
      if (!is.null(metrics_pulses)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_pulses))
      }
      if (!is.null(metrics_baseflow)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_baseflow))
      }
      if (!is.null(metrics_recession)) {
        gage_row <- cbind(gage_row, as.data.table(metrics_recession))
      }
      
      # Add NA columns for Q-PPT metrics that we can't calculate
      q_ppt_cols <- c("annual_runoff_ratio_slp", "annual_runoff_ratio_rho", 
                      "annual_runoff_ratio_pval", "annual_runoff_ratio_mean", 
                      "annual_runoff_ratio_median",
                      "winter_runoff_ratio_slp", "winter_runoff_ratio_rho", 
                      "winter_runoff_ratio_pval", "winter_runoff_ratio_mean", 
                      "winter_runoff_ratio_median",
                      "spring_runoff_ratio_slp", "spring_runoff_ratio_rho", 
                      "spring_runoff_ratio_pval", "spring_runoff_ratio_mean", 
                      "spring_runoff_ratio_median",
                      "summer_runoff_ratio_slp", "summer_runoff_ratio_rho", 
                      "summer_runoff_ratio_pval", "summer_runoff_ratio_mean", 
                      "summer_runoff_ratio_median",
                      "fall_runoff_ratio_slp", "fall_runoff_ratio_rho", 
                      "fall_runoff_ratio_pval", "fall_runoff_ratio_mean", 
                      "fall_runoff_ratio_median")
      
      for (col in q_ppt_cols) {
        gage_row[[col]] <- NA_real_
      }
      
      # Append to summary output
      summary_output <- rbind(summary_output, gage_row, fill = TRUE)
      
      # Save after each successful processing
      fwrite(summary_output, output_file)
      cat("Successfully processed gage", gage_id, "\n")
      
    }, error = function(e) {
      cat("Error processing gage", gage_id, ":", e$message, "\n")
    })
  }
  
  cat("Finished processing all gages. Final summary has", nrow(summary_output), "rows\n")
  return(summary_output)
}


# MAIN WRAPPER FUNCTION FOR PROCESSING RAW STREAMFLOW DATA VIA API AND SAVING AS BUNDLED DATA (E.G., PARQUET)
process_gages_rawToRaw <- function(gages_df, gage_type, min_num_years, start_date, end_date, 
                                   min_Q_value_and_days, output_dir, 
                                   storage_format = c("parquet", "csv", "rds", "feather"),
                                   chunk_size = 1000) {
  
  storage_format <- match.arg(storage_format)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize counters
  processed_gages <- 0
  failed_gages <- 0
  
  # Create metadata file
  metadata_file <- file.path(output_dir, "watershed_metadata.csv")
  metadata <- data.table(
    gage_id = character(),
    latitude = numeric(),
    longitude = numeric(),
    basin_area = numeric(),
    gage_type = character(),
    num_years = integer(),
    start_year = integer(),
    end_year = integer(),
    num_days = integer(),
    processing_status = character(),
    error_message = character()
  )
  
  # Process in chunks for memory efficiency
  chunk_num <- 1
  daily_data_chunk <- data.table()
  
  for (i in 1:nrow(gages_df)) {
    current_gage <- gages_df[i, ]
    
    # Extract gage ID and coordinate info based on gage type
    if (gage_type == "USGS") {
      gage_id <- as.character(current_gage$STAID)
      latitude <- as.numeric(current_gage$LAT_GAGE)
      longitude <- as.numeric(current_gage$LNG_GAGE)
      basin_area <- as.numeric(current_gage$DRAIN_SQKM)
    } else if (gage_type %in% c("Canada", "CANADIAN")) {
      names(current_gage)[names(current_gage) == 'LATITUDE'] <- 'LAT_GAGE'
      names(current_gage)[names(current_gage) == 'LONGITUDE'] <- 'LNG_GAGE'
      gage_id <- as.character(current_gage$STATION_NUMBER)
      latitude <- as.numeric(current_gage$LAT_GAGE)
      longitude <- as.numeric(current_gage$LNG_GAGE)
      basin_area <- NA
    } else {
      stop("Unsupported gage type")
    }
    
    cat("Processing gage", gage_id, "(", i, "of", nrow(gages_df), ")\n")
    
    # Initialize metadata row
    meta_row <- data.table(
      gage_id = gage_id,
      latitude = latitude,
      longitude = longitude,
      basin_area = basin_area,
      gage_type = gage_type,
      num_years = NA_integer_,
      start_year = NA_integer_,
      end_year = NA_integer_,
      num_days = NA_integer_,
      processing_status = "processing",
      error_message = NA_character_
    )
    
    tryCatch({
      # Get streamflow data
      streamflow_data <- generate_streamflow_dt(current_gage, gage_type, 
                                                min_num_years, start_date, end_date)
      
      if (is.null(streamflow_data) || identical(streamflow_data, NA) ||
          (is.data.frame(streamflow_data) && nrow(streamflow_data) == 0)) {
        meta_row$processing_status <- "no_data"
        meta_row$error_message <- "No valid streamflow data retrieved"
        metadata <- rbind(metadata, meta_row)
        failed_gages <- failed_gages + 1
        next
      }
      
      # Gate based on qualifying years, but do not drop non-qualifying years once accepted
      years_to_use <- NULL
      for (this_year in unique(streamflow_data$year)) {
        test_year <- subset(streamflow_data, year == this_year)
        if (!is.numeric(test_year$Q)) test_year$Q <- as.numeric(test_year$Q)
        valid_q_values <- test_year$Q[!is.na(test_year$Q)]
        nonzero_rows <- which(valid_q_values > min_Q_value_and_days[1])
        if (length(nonzero_rows) >= min_Q_value_and_days[2]) {
          years_to_use <- c(years_to_use, this_year)
        }
      }
      
      if (length(years_to_use) < min_num_years) {
        meta_row$processing_status <- "insufficient_years"
        meta_row$error_message <- paste("Only", length(years_to_use), "years with valid data")
        metadata <- rbind(metadata, meta_row)
        failed_gages <- failed_gages + 1
        next
      }
      
      # IMPORTANT: Keep ALL data within the requested window once the site qualifies
      streamflow_data_kept <- streamflow_data[
        streamflow_data$Date >= start_date & streamflow_data$Date <= end_date, ]
      
      # Metadata semantics unchanged:
      # - num_years/start_year/end_year = qualifying years
      # - num_days = number of days in those qualifying years (as before)
      meta_row$num_years  <- length(years_to_use)
      meta_row$start_year <- min(years_to_use)
      meta_row$end_year   <- max(years_to_use)
      meta_row$num_days   <- nrow(streamflow_data[streamflow_data$year %in% years_to_use, ])
      meta_row$processing_status <- "success"
      
      # Add gage_id and keep the flag in the saved output
      streamflow_data_kept$gage_id <- gage_id
      
      cols_to_keep <- c("gage_id", "Date", "Q", "year", "month", "doy", "water_year", "dowy", "flag")
      if ("PPT" %in% names(streamflow_data_kept)) cols_to_keep <- c(cols_to_keep, "PPT")
      if ("SWE" %in% names(streamflow_data_kept)) cols_to_keep <- c(cols_to_keep, "SWE")
      
      streamflow_data_kept <- as.data.table(streamflow_data_kept)[, ..cols_to_keep]
      
      # Add to chunk
      daily_data_chunk <- rbind(daily_data_chunk, streamflow_data_kept, fill = TRUE)
      processed_gages <- processed_gages + 1
      
      # Save chunk when it reaches chunk_size
      if (processed_gages %% chunk_size == 0) {
        save_chunk(daily_data_chunk, output_dir, chunk_num, storage_format)
        daily_data_chunk <- data.table()
        chunk_num <- chunk_num + 1
        
        fwrite(metadata, metadata_file)
        cat("\nSaved chunk", chunk_num - 1, "with", chunk_size, "watersheds\n")
      }
      
    }, error = function(e) {
      meta_row$processing_status <- "error"
      meta_row$error_message <- as.character(e$message)
      failed_gages <- failed_gages + 1
      cat("Error processing gage", gage_id, ":", e$message, "\n")
    })
    
    # Add metadata row
    metadata <- rbind(metadata, meta_row)
  }
  
  # Save final chunk if it exists
  if (nrow(daily_data_chunk) > 0) {
    save_chunk(daily_data_chunk, output_dir, chunk_num, storage_format)
  }
  
  # Save final metadata
  fwrite(metadata, metadata_file)
  
  # Create summary report
  summary_stats <- list(
    total_gages = nrow(gages_df),
    processed_successfully = processed_gages,
    failed_gages = failed_gages,
    output_directory = output_dir,
    storage_format = storage_format,
    chunk_size = chunk_size,
    num_chunks = chunk_num
  )
  
  saveRDS(summary_stats, file.path(output_dir, "processing_summary.rds"))
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Processing Complete!\n")
  cat("Total gages:", summary_stats$total_gages, "\n")
  cat("Successfully processed:", summary_stats$processed_successfully, "\n")
  cat("Failed:", summary_stats$failed_gages, "\n")
  cat("Output saved to:", output_dir, "\n")
  cat("Storage format:", storage_format, "\n")
  cat("Number of chunks:", summary_stats$num_chunks, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  return(metadata)
}



# Helper function to save data chunks in different formats
save_chunk <- function(data_chunk, output_dir, chunk_num, storage_format) {
  chunk_file <- file.path(output_dir, paste0("daily_data_chunk_", 
                                             sprintf("%04d", chunk_num)))
  
  if (storage_format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required for parquet format. Install with: install.packages('arrow')")
    }
    arrow::write_parquet(data_chunk, paste0(chunk_file, ".parquet"))
    
  } else if (storage_format == "feather") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required for feather format. Install with: install.packages('arrow')")
    }
    arrow::write_feather(data_chunk, paste0(chunk_file, ".feather"))
    
  } else if (storage_format == "csv") {
    fwrite(data_chunk, paste0(chunk_file, ".csv"))
    
  } else if (storage_format == "rds") {
    saveRDS(data_chunk, paste0(chunk_file, ".rds"))
  }
}

# Helper function to read all chunks back into memory (use with caution for large datasets)
read_all_chunks <- function(output_dir, storage_format = "parquet") {
  # Get list of chunk files
  pattern <- paste0("daily_data_chunk_.*\\.", storage_format)
  chunk_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  
  if (length(chunk_files) == 0) {
    stop("No chunk files found in ", output_dir)
  }
  
  # Read and combine all chunks
  all_data <- data.table()
  
  for (file in chunk_files) {
    cat("Reading", basename(file), "...\n")
    
    if (storage_format == "parquet") {
      chunk <- arrow::read_parquet(file)
    } else if (storage_format == "feather") {
      chunk <- arrow::read_feather(file)
    } else if (storage_format == "csv") {
      chunk <- fread(file)
    } else if (storage_format == "rds") {
      chunk <- readRDS(file)
    }
    
    all_data <- rbind(all_data, chunk, fill = TRUE)
  }
  
  return(all_data)
}

# Helper function to query specific watersheds efficiently
query_watersheds <- function(output_dir, gage_ids, storage_format = "parquet") {
  pattern <- paste0("daily_data_chunk_.*\\.", storage_format)
  chunk_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  
  result <- data.table()
  
  for (file in chunk_files) {
    if (storage_format == "parquet") {
      # Efficient filtering for parquet
      chunk <- arrow::read_parquet(file, 
                                   col_select = everything(),
                                   as_data_frame = TRUE)
    } else if (storage_format == "csv") {
      chunk <- fread(file)
    } else if (storage_format == "rds") {
      chunk <- readRDS(file)
    }
    
    # Filter for requested gage_ids
    chunk_filtered <- chunk[chunk$gage_id %in% gage_ids, ]
    
    if (nrow(chunk_filtered) > 0) {
      result <- rbind(result, chunk_filtered)
    }
  }
  
  return(result)
}



# CONVERT RAW DATA TO STANDARDIZED DATA TABLE
generate_streamflow_dt <- function(dt, data_origin, 
                                   min_num_years = 20, 
                                   start_date = as.Date("1900-01-01"), 
                                   end_date = as.Date("2024-12-31")) {
  if (!data_origin %in% c("USGS", "Canada")) {
    warning("Invalid data_origin provided. It must be either 'USGS' or 'Canada'. Returning NA.")
    return(NULL)
  }
  if (!inherits(dt, "data.table")) {
    dt <- as.data.table(dt)
  }
  output <- NA
  
  if (data_origin == "USGS") {
    gage_data <- readNWISdv(siteNumber = dt$STAID,
                            parameterCd = "00060",
                            startDate = "1900-01-01",
                            endDate = as.character(end_date))
    # Enforce the requested window
    gage_data <- subset(gage_data, Date >= start_date & Date <= end_date)
    if (nrow(gage_data) == 0) return(NA)
    
    # Try to standardize column names and find value/code
    gd2 <- tryCatch(dataRetrieval::renameNWISColumns(gage_data), error = function(e) gage_data)
    if ("Flow" %in% names(gd2) && "Flow_cd" %in% names(gd2)) {
      val_col <- "Flow"
      code_col <- "Flow_cd"
      gage_data <- gd2
    } else {
      # Fallback to raw names
      candidates_val <- grep("^X_00060_00003$", names(gage_data), value = TRUE)
      candidates_cd  <- grep("^X_00060_00003_cd$", names(gage_data), value = TRUE)
      if (length(candidates_val) == 1 && length(candidates_cd) == 1) {
        val_col <- candidates_val
        code_col <- candidates_cd
      } else {
        # Last resort: original assumption (4th is value, 5th is code)
        val_col <- names(gage_data)[4]
        code_col <- names(gage_data)[5]
      }
    }
    
    if (nrow(gage_data) > 365 * min_num_years && max(gage_data$Date) > start_date) {
      streamy <- gage_data[, c("Date", val_col, code_col)]
      names(streamy) <- c("Date", "Q_rawUnits", "flag")
      
      # Mask Q as before by acceptable codes, but keep the flag column
      keep_codes <- c("A", "A e", "P", "P e")
      streamy$Q_rawUnits[!(streamy$flag %in% keep_codes)] <- NA_real_
      
      # Convert to mm/day using drainage area (sqkm) from metadata row
      gage_id <- gage_data$site_no[1]
      sqkm <- dt$DRAIN_SQKM[dt$STAID == gage_id]
      conversion <- 60 * 60 * 24 / (sqkm * 3280.84^3) * 1e6
      streamy$Q <- as.numeric(streamy$Q_rawUnits) * conversion
      
      streamy$year  <- lubridate::year(streamy$Date)
      streamy$month <- lubridate::month(streamy$Date)
      streamy$doy   <- lubridate::yday(streamy$Date)
      
      wy_info <- calculate_water_year_info(streamy$Date)
      streamy$water_year <- wy_info$water_year
      streamy$dowy       <- wy_info$dowy
      
      output <- streamy
    } else {
      message("Insufficient Data to Process")
      output <- NA
    }
  }
  
  if (data_origin == "Canada") {
    can_stream <- hy_daily(station_number = paste(dt$STATION_NUMBER))
    can_stream_only <- subset(can_stream, Parameter == "Flow")
    # Enforce the requested window
    can_stream_only <- can_stream_only[
      can_stream_only$Date >= start_date & can_stream_only$Date <= end_date, ]
    
    if ("Flow" %in% can_stream$Parameter &&
        nrow(can_stream_only) > 365 * min_num_years) {
      streamy <- data.frame(
        Date = as.Date(can_stream_only$Date),
        Q_rawUnits = can_stream_only$Value,
        flag = can_stream_only$Symbol,   # HYDAT qualifier
        stringsAsFactors = FALSE
      )
      
      # Convert from m^3/s to mm/day
      sqkm <- hy_stations(paste(dt$STATION_NUMBER))$DRAINAGE_AREA_GROSS
      conversion <- ifelse(is.na(sqkm), 99999, 60 * 60 * 24 * 1e9 / (sqkm * 1e12))
      streamy$Q <- as.numeric(streamy$Q_rawUnits) * conversion
      
      streamy$year  <- lubridate::year(streamy$Date)
      streamy$month <- lubridate::month(streamy$Date)
      streamy$doy   <- lubridate::yday(streamy$Date)
      
      wy_info <- calculate_water_year_info(streamy$Date)
      streamy$water_year <- wy_info$water_year
      streamy$dowy       <- wy_info$dowy
      
      output <- streamy
    } else {
      message("Insufficient Data to Process")
      output <- NA
    }
  }
  
  return(output)
}




# DELINEATE BASINS VIA HYDROBASINS AND SAVE METADATA
# recursive search for all hydrobasins upstream of a a selected basin; only needed for delineating basins as defined by hydrobasins
find_upstream_hydrobasins <- function(current_gage, basinAt_NorAm_polys, HB_dt, upstream_hydrobasins = list(), save_path = NULL) {
  if (missing(current_gage) || missing(basinAt_NorAm_polys) || missing(HB_dt)) {
    stop("Required inputs missing: current_gage, basinAt_NorAm_polys, or HB_dt")
  }
  
  required_cols <- c("LNG_GAGE", "LAT_GAGE")
  if (!all(required_cols %in% colnames(current_gage))) {
    stop("current_gage must contain LNG_GAGE and LAT_GAGE columns")
  }
  
  current_gage$LNG_GAGE <- as.numeric(as.character(current_gage$LNG_GAGE))
  current_gage$LAT_GAGE <- as.numeric(as.character(current_gage$LAT_GAGE))
  
  # Ensure the basin polygon IDs are characters.
  basinAt_NorAm_polys$HYBAS_ID <- as.character(basinAt_NorAm_polys$HYBAS_ID)
  
  tryCatch({
    gage_point <- sf::st_point(c(current_gage$LNG_GAGE, current_gage$LAT_GAGE))
    gage_point_sfc <- sf::st_sfc(gage_point, crs = 4326)
    
    intersect_result <- sf::st_intersects(gage_point_sfc, basinAt_NorAm_polys, prepared = FALSE)
    if (length(intersect_result[[1]]) == 0) {
      warning("No hydrobasin found containing the gage location at: ", current_gage$LNG_GAGE, ", ", current_gage$LAT_GAGE)
      return(NULL)
    }
    
    gage_basin <- as.character(basinAt_NorAm_polys$HYBAS_ID[intersect_result[[1]][1]])
    if (gage_basin %in% names(upstream_hydrobasins)) {
      return(upstream_hydrobasins[[gage_basin]])
    }
    
    these_hydro_basins <- gage_basin
    HB_copy <- copy(HB_dt)
    HB_copy$HYBAS_ID <- as.character(HB_copy$HYBAS_ID)
    HB_copy$NEXT_DOWN <- as.character(HB_copy$NEXT_DOWN)
    
    HB_remaining <- HB_copy[HB_copy$HYBAS_ID != gage_basin]
    data.table::setindex(HB_remaining, NEXT_DOWN)
    basins_to_check <- these_hydro_basins
    
    while (length(basins_to_check) > 0) {
      new_basins <- character(0)
      for (basin in basins_to_check) {
        upstream_rows <- which(HB_remaining$NEXT_DOWN == basin)
        if (length(upstream_rows) > 0) {
          upstream_ids <- HB_remaining$HYBAS_ID[upstream_rows]
          new_basins <- c(new_basins, upstream_ids)
          HB_remaining <- HB_remaining[-upstream_rows]
        }
      }
      these_hydro_basins <- c(these_hydro_basins, new_basins)
      basins_to_check <- new_basins
      if (length(new_basins) == 0) break
    }
    
    if (!is.null(save_path)) {
      tryCatch({
        if (file.exists(save_path)) {
          existing_basins <- readRDS(save_path)
          existing_basins[[gage_basin]] <- these_hydro_basins
          saveRDS(existing_basins, save_path)
        } else {
          new_list <- list()
          new_list[[gage_basin]] <- these_hydro_basins
          saveRDS(new_list, save_path)
        }
      }, error = function(e) {
        warning("Failed to save upstream_hydrobasins to ", save_path, ": ", e$message)
      })
    }
    
    return(these_hydro_basins)
  }, error = function(e) {
    warning("Error in find_upstream_hydrobasins: ", e$message)
    return(NULL)
  })
}







################################################################################
# read and merge all parquet files

# Load required libraries
library(duckdb)
library(arrow)
library(data.table)

# Function to concatenate multiple parquet directories into one
concatenate_parquet_directories <- function(input_dirs, 
                                            output_file = "combined_streamflow_data.parquet",
                                            method = c("duckdb", "arrow")) {
  
  method <- match.arg(method)
  
  # Get all parquet files from all directories
  all_parquet_files <- character()
  
  for (dir in input_dirs) {
    files <- list.files(dir, 
                        pattern = "daily_data_chunk_.*\\.parquet$", 
                        full.names = TRUE)
    all_parquet_files <- c(all_parquet_files, files)
    cat("Found", length(files), "parquet files in", dir, "\n")
  }
  
  cat("Total files to concatenate:", length(all_parquet_files), "\n\n")
  
  if (method == "duckdb") {
    concatenate_with_duckdb(all_parquet_files, output_file)
  } else {
    concatenate_with_arrow(all_parquet_files, output_file)
  }
}

# Method 1: Using DuckDB (recommended for large files)
concatenate_with_duckdb <- function(parquet_files, output_file) {
  cat("Using DuckDB method...\n")
  
  # Create a DuckDB connection
  con <- dbConnect(duckdb::duckdb())
  
  tryCatch({
    # Create the UNION ALL query
    union_parts <- character()
    
    for (i in seq_along(parquet_files)) {
      file_path <- parquet_files[i]
      # Escape single quotes in file path
      file_path <- gsub("'", "''", file_path)
      
      if (i == 1) {
        union_parts <- paste0("SELECT * FROM read_parquet('", file_path, "')")
      } else {
        union_parts <- paste0(union_parts, " UNION ALL SELECT * FROM read_parquet('", file_path, "')")
      }
      
      # Process in batches to avoid overly long SQL queries
      if (i %% 50 == 0 || i == length(parquet_files)) {
        cat("Processing files", max(1, i-49), "to", i, "of", length(parquet_files), "...\n")
      }
    }
    
    # Execute the query and write to parquet
    cat("\nCombining all files and writing to", output_file, "...\n")
    query <- paste0("COPY (", union_parts, ") TO '", output_file, "' (FORMAT PARQUET, COMPRESSION 'SNAPPY')")
    
    dbExecute(con, query)
    
    cat("Successfully created:", output_file, "\n")
    
    # Get some statistics about the combined file
    stats_query <- paste0("SELECT COUNT(*) as total_rows, 
                          COUNT(DISTINCT gage_id) as n_gages,
                          MIN(Date) as start_date,
                          MAX(Date) as end_date
                          FROM read_parquet('", output_file, "')")
    
    stats <- dbGetQuery(con, stats_query)
    cat("\nCombined file statistics:\n")
    cat("Total rows:", format(stats$total_rows, big.mark = ","), "\n")
    cat("Number of gages:", stats$n_gages, "\n")
    cat("Date range:", as.character(stats$start_date), "to", as.character(stats$end_date), "\n")
    
  }, finally = {
    dbDisconnect(con, shutdown = TRUE)
  })
}

# Method 2: Using Arrow (alternative method)
concatenate_with_arrow <- function(parquet_files, output_file) {
  cat("Using Arrow method...\n")
  
  # Open dataset from multiple files
  dataset <- open_dataset(parquet_files, format = "parquet")
  
  # Write combined dataset to single file
  cat("Writing combined file to", output_file, "...\n")
  write_parquet(dataset, output_file)
  
  cat("Successfully created:", output_file, "\n")
}

# Function to concatenate and create metadata summary
concatenate_with_metadata <- function(input_dirs, output_dir = "combined_streamflow_data") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Combine metadata files first
  cat("Combining metadata files...\n")
  all_metadata <- data.table()
  
  for (dir in input_dirs) {
    metadata_file <- file.path(dir, "watershed_metadata.csv")
    if (file.exists(metadata_file)) {
      metadata <- fread(metadata_file)
      metadata$source_dir <- basename(dir)
      all_metadata <- rbind(all_metadata, metadata, fill = TRUE)
    }
  }
  
  # Save combined metadata
  fwrite(all_metadata, file.path(output_dir, "combined_watershed_metadata.csv"))
  cat("Combined metadata saved with", nrow(all_metadata), "watersheds\n")
  
  # Print summary statistics
  cat("\nMetadata summary:\n")
  cat("Total watersheds:", nrow(all_metadata), "\n")
  cat("Successful:", sum(all_metadata$processing_status == "success"), "\n")
  cat("Failed:", sum(all_metadata$processing_status != "success"), "\n")
  cat("By source:\n")
  print(table(all_metadata$source_dir, all_metadata$processing_status))
  
  # Concatenate parquet files
  output_parquet <- file.path(output_dir, "combined_streamflow_data.parquet")
  concatenate_parquet_directories(input_dirs, output_parquet, method = "duckdb")
  
  return(list(
    metadata = all_metadata,
    output_file = output_parquet
  ))
}




add_downstream_basin_ids <- function(metadata_file_path, 
                                     basinAt_NorAm_polys, 
                                     HB_dt, 
                                     upstream_hydrobasins_path = "upstream_hydrobasins.rds") {
  
  # Load the combined metadata
  cat("Loading metadata from:", metadata_file_path, "\n")
  metadata <- fread(metadata_file_path)
  
  # Add Downstream_HB_ID column if it doesn't exist
  if (!"Downstream_HB_ID" %in% names(metadata)) {
    metadata[, Downstream_HB_ID := NA_character_]
  }
  
  # Initialize upstream_hydrobasins
  upstream_hydrobasins <- list()
  
  # Handle both .rds and .RData extensions
  if (grepl("\\.RData$", upstream_hydrobasins_path)) {
    # User provided .RData path, convert to .rds
    rds_path <- gsub("\\.RData$", ".rds", upstream_hydrobasins_path)
    cat("Note: Converting to RDS format. Will use:", rds_path, "\n")
  } else {
    rds_path <- upstream_hydrobasins_path
  }
  
  # Try to load existing data
  if (file.exists(rds_path)) {
    tryCatch({
      upstream_hydrobasins <- readRDS(rds_path)
      cat("Loaded existing upstream_hydrobasins with", length(upstream_hydrobasins), "basins\n")
    }, error = function(e) {
      cat("Could not load existing RDS file. Starting fresh.\n")
      upstream_hydrobasins <- list()
    })
  } else if (file.exists(upstream_hydrobasins_path) && upstream_hydrobasins_path != rds_path) {
    # Try to load the incorrectly saved .RData file (which is actually RDS format)
    tryCatch({
      # Try readRDS first since find_upstream_hydrobasins uses saveRDS
      upstream_hydrobasins <- readRDS(upstream_hydrobasins_path)
      # Save with correct extension
      saveRDS(upstream_hydrobasins, rds_path)
      cat("Loaded existing data from", upstream_hydrobasins_path, "and converted to", rds_path, "\n")
      cat("Found", length(upstream_hydrobasins), "basins\n")
    }, error = function(e) {
      cat("Could not load existing file. Starting fresh.\n")
      upstream_hydrobasins <- list()
    })
  }
  
  # Filter for successful watersheds
  successful_watersheds <- metadata[processing_status == "success"]
  cat("Processing", nrow(successful_watersheds), "successful watersheds\n\n")
  
  # Process each successful watershed
  for (i in 1:nrow(successful_watersheds)) {
    current_watershed <- successful_watersheds[i,]
    
    # Create a gage object
    current_gage <- data.frame(
      LAT_GAGE = current_watershed$latitude,
      LNG_GAGE = current_watershed$longitude
    )
    
    cat("Processing", i, "of", nrow(successful_watersheds), 
        "- Gage:", current_watershed$gage_id, 
        "(", round(current_watershed$longitude, 3), ",", round(current_watershed$latitude, 3), ") ...")
    
    tryCatch({
      # Find the basin containing this gage
      gage_sf <- st_as_sf(current_gage, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326)
      gage_sf_transformed <- st_transform(gage_sf, st_crs(basinAt_NorAm_polys))
      
      intersection <- st_intersects(gage_sf_transformed, basinAt_NorAm_polys)
      
      if (length(intersection[[1]]) > 0) {
        basin_index <- intersection[[1]][1]
        hydro_id <- basinAt_NorAm_polys$HYBAS_ID[basin_index]
        
        if (!is.null(hydro_id) && length(hydro_id) > 0 && !is.na(hydro_id)) {
          downstream_basin_id <- as.character(hydro_id)
          
          if (length(downstream_basin_id) == 1 && !is.na(downstream_basin_id)) {
            metadata[gage_id == current_watershed$gage_id, Downstream_HB_ID := downstream_basin_id]
            
            cat(" Basin ID:", downstream_basin_id)
            
            # Check if upstream basins need to be calculated
            if (!downstream_basin_id %in% names(upstream_hydrobasins)) {
              cat(" - Finding upstream basins...")
              
              # Call find_upstream_hydrobasins with the RDS path
              upstream_basins <- tryCatch({
                find_upstream_hydrobasins(
                  current_gage = current_gage,
                  basinAt_NorAm_polys = basinAt_NorAm_polys,
                  HB_dt = HB_dt,
                  upstream_hydrobasins = upstream_hydrobasins,
                  save_path = rds_path  # Use RDS path since function uses saveRDS
                )
              }, error = function(e) {
                cat(" Error in find_upstream_hydrobasins:", e$message)
                NULL
              })
              
              if (!is.null(upstream_basins)) {
                # Update our local copy
                upstream_hydrobasins[[downstream_basin_id]] <- upstream_basins
                cat(" - Found", length(upstream_basins), "upstream basins")
                
                # Reload to get the updated version
                if (file.exists(rds_path)) {
                  tryCatch({
                    upstream_hydrobasins <- readRDS(rds_path)
                  }, error = function(e) {
                    # Continue with local copy
                  })
                }
              }
            } else {
              cat(" - Upstream basins already calculated")
            }
          } else {
            metadata[gage_id == current_watershed$gage_id, Downstream_HB_ID := "INVALID_HYDROID"]
          }
        } else {
          metadata[gage_id == current_watershed$gage_id, Downstream_HB_ID := "NO_HYDROID"]
        }
        
        cat("\n")
        
      } else {
        cat(" - No basin intersection found\n")
        metadata[gage_id == current_watershed$gage_id, Downstream_HB_ID := "NO_BASIN_FOUND"]
      }
      
    }, error = function(e) {
      error_msg <- as.character(e$message)
      cat(" - Error:", substr(error_msg, 1, 100), "\n")
      metadata[gage_id == current_watershed$gage_id, 
               Downstream_HB_ID := paste0("ERROR: ", substr(error_msg, 1, 50))]
    })
    
    # Save progress periodically
    if (i %% 100 == 0) {
      fwrite(metadata, metadata_file_path)
      cat("Saved progress at", i, "watersheds\n")
    }
  }
  
  # Save final metadata
  fwrite(metadata, metadata_file_path)
  
  cat("\nCompleted! Updated metadata saved to:", metadata_file_path, "\n")
  cat("Upstream basins data saved to:", rds_path, "\n")
  
  # Print summary
  cat("\nSummary:\n")
  cat("Total watersheds:", nrow(metadata), "\n")
  cat("Successful watersheds processed:", nrow(successful_watersheds), "\n")
  
  basin_count <- sum(!is.na(metadata$Downstream_HB_ID) & 
                       !metadata$Downstream_HB_ID %in% c("NO_BASIN_FOUND", "NO_HYDROID", "INVALID_HYDROID") & 
                       !grepl("ERROR:", metadata$Downstream_HB_ID), na.rm = TRUE)
  
  cat("Watersheds with basin ID:", basin_count, "\n")
  cat("Watersheds without basin:", sum(metadata$Downstream_HB_ID == "NO_BASIN_FOUND", na.rm = TRUE), "\n")
  cat("Unique basins with upstream data:", length(upstream_hydrobasins), "\n")
  
  return(metadata)
}


