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

# Helper to wrap metric calculations with error handling
# Returns the function result on success, NULL on error (with logging)
safe_calculate <- function(fn, data, gage_id, label, ctx, log_level = "warn") {
  tryCatch({
    fn(data)
  }, error = function(e) {
    if (log_level == "warn") {
      log_warn(label, "failed for gage", gage_id, ":", e$message, context = ctx)
    } else {
      log_debug(label, "failed for gage", gage_id, ":", e$message, context = ctx)
    }
    NULL
  })
}

##### Streamflow signatures
# tools:
#	BFI:
#		DFI				Stoelzle  et al. 2020, partially relying on baseflow() in lfstat
#		WMO				by World Meteorological Organisation WMO using baseflow() in lfstat
#		NM90			Nathan and McMahon 1990 using BaseflowSeparation() in EcoHydRology
#		LH13			Ladson et al 2013 using baseflows() in hydrostats; Ladson, A. R., R. Brown, B. Neal and R. Nathan (2013) A standard approach to baseflow separation using the Lyne and Hollick filter. Australian Journal of Water Resources 17(1): 173-18 and Lynne, V., Hollick, M. (1979) Stochastic time-variable rainfall-runoff modelling. In: pp. 89-93 Institute of Engineers Australia National Conference. Perth.
#		E05			Eckhardt 2005 1 parameter recursive digital filter using bf_oneparam() in FlowScreen
#		E12			Eckhardt 2012 2 paramater recursive digital filter using bf_eckhardt() in FlowScreen
#		B93			Boughton 1993 2 parameter recursive digital filter using bf_boughton() in FlowScreen
#	monthly
#	seasonsally
#	Qxx



#########################################################################################
# process gage data

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


# original function - processes api data directly, outputs metrics and trends
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
      # Fetch drainage area from tidyhydat database
      basin_area <- tryCatch({
        stn_info <- hy_stations(gage_id)
        as.numeric(stn_info$DRAINAGE_AREA_GROSS)
      }, error = function(e) NA_real_)
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

      ctx_raw <- "process_gages_rawData"
      metric_specs_raw <- list(
        list(name = "flow_vols",  fn = calculate_flow_vols_by_year,        label = "Flow volumes"),
        list(name = "fdc_trends", fn = analyze_fdc_trends_from_streamflow, label = "FDC trends"),
        list(name = "flashiness", fn = analyze_flashiness_trends,          label = "Flashiness"),
        list(name = "flow_timing",fn = analyze_flow_timing_trends,         label = "Flow timing"),
        list(name = "pulses",     fn = calculate_pulse_metrics,            label = "Pulse metrics"),
        list(name = "baseflow",   fn = analyze_baseflow_indices,           label = "Baseflow indices"),
        list(name = "recession",  fn = analyze_recession_parameters,       label = "Recession parameters")
      )

      metrics_list_raw <- list()
      for (spec in metric_specs_raw) {
        metrics_list_raw[[spec$name]] <- safe_calculate(spec$fn, streamflow_data_filtered,
                                                         gage_id, spec$label, ctx_raw)
      }

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

      # Combine base gage info with successfully calculated metrics
      for (metric_name in names(metrics_list_raw)) {
        if (!is.null(metrics_list_raw[[metric_name]])) {
          gage_row <- cbind(gage_row, as.data.table(metrics_list_raw[[metric_name]]))
        }
      }
      
      # Add NA columns for Q-PPT metrics that we can't calculate
      q_ppt_bases <- c("annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
                       "summer_runoff_ratio", "fall_runoff_ratio")
      q_ppt_suffixes <- c("_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                          "_mk_rho", "_mk_pval", "_mean", "_median")
      q_ppt_cols <- as.vector(outer(q_ppt_bases, q_ppt_suffixes, paste0))
      
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


  # processes api data and outputs timeseries of streamflow and an associated metadata file
process_gages_rawToRaw <- function(gages_df, gage_type, min_num_years, start_date, end_date, 
                                   min_Q_value_and_days, output_dir, 
                                   storage_format = c("parquet", "csv", "rds", "feather"),
                                   chunk_size = 1000) {
  
  storage_format <- match.arg(storage_format)
  force(start_date)
  force(end_date)

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
    error_message = character(),
    area_normalized = logical()
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
      # Fetch drainage area from tidyhydat database
      basin_area <- tryCatch({
        stn_info <- hy_stations(gage_id)
        as.numeric(stn_info$DRAINAGE_AREA_GROSS)
      }, error = function(e) NA_real_)
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
      error_message = NA_character_,
      area_normalized = !is.na(basin_area)
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








  # helper funciton that generates a data table for holding streamflow data; bespoke per data source
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
      if (length(sqkm) == 0 || is.na(sqkm)) {
        log_warn("Missing basin area for USGS gage", gage_id,
                 "- Q kept in raw units (cfs), not area-normalized", context = "process_gage")
        conversion <- 1
      } else {
        conversion <- 60 * 60 * 24 / (sqkm * 3280.84^3) * 1e6
      }
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
      if (length(sqkm) == 0 || is.na(sqkm)) {
        log_warn("Missing basin area for Canadian gage", dt$STATION_NUMBER,
                 "- Q kept in raw units (m3/s), not area-normalized", context = "process_gage")
        conversion <- 1
      } else {
        conversion <- 60 * 60 * 24 * 1e9 / (sqkm * 1e12)
      }
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


################################################################
# statistical processing functions

  # helper function that receives a time series and outputs summary metrics (trends, averages, etc.)
  # Returns 8 statistics per metric:
  #   _senn_slp     = Theil-Sen slope (robust non-parametric trend)
  #   _linear_slp   = Linear regression slope (parametric trend)
  #   _spearman_rho = Spearman's rank correlation coefficient
  #   _spearman_pval= Spearman's p-value for trend significance
  #   _mk_rho       = Mann-Kendall tau (non-parametric trend correlation)
  #   _mk_pval      = Mann-Kendall p-value for trend significance
  #   _mean         = Arithmetic mean across water years
  #   _median       = Median across water years
generate_stats <- function(data, value_cols = NULL, year_col = "year", min_rows = 3) {
  # Check if required packages are available
  if (!requireNamespace("zyp", quietly = TRUE)) {
    stop("Package 'zyp' is needed for this function. Please install it with install.packages('zyp')")
  }
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    stop("Package 'Kendall' is needed for this function. Please install it with install.packages('Kendall')")
  }

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame or data.table")
  }

  if (!year_col %in% colnames(data)) {
    stop(paste("Year column '", year_col, "' not found in data"))
  }

  # Ensure data is sorted by year for trend calculations
  data <- data[order(data[[year_col]]), ]

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
      # Return NAs for all 8 stats
      results[[paste0(col, "_senn_slp")]] <- NA
      results[[paste0(col, "_linear_slp")]] <- NA
      results[[paste0(col, "_spearman_rho")]] <- NA
      results[[paste0(col, "_spearman_pval")]] <- NA
      results[[paste0(col, "_mk_rho")]] <- NA
      results[[paste0(col, "_mk_pval")]] <- NA
      results[[paste0(col, "_mean")]] <- NA
      results[[paste0(col, "_median")]] <- NA
      next
    }

    # Calculate Theil-Sen slope (robust non-parametric)
    sen_result <- try(zyp::zyp.sen(value ~ year, data = working_data), silent = TRUE)
    if (inherits(sen_result, "try-error")) {
      sen_slope <- NA
    } else {
      sen_slope <- sen_result$coefficients[2]
    }

    # Calculate Linear regression slope (parametric)
    lm_result <- try(lm(value ~ year, data = working_data), silent = TRUE)
    if (inherits(lm_result, "try-error")) {
      linear_slope <- NA
    } else {
      linear_slope <- coef(lm_result)[2]
    }

    # Calculate Spearman correlation
    spearman_result <- try(cor.test(working_data$year, working_data$value,
                                    method = "spearman"), silent = TRUE)
    if (inherits(spearman_result, "try-error")) {
      spearman_rho <- NA
      spearman_pval <- NA
    } else {
      spearman_rho <- spearman_result$estimate
      spearman_pval <- spearman_result$p.value
    }

    # Calculate Mann-Kendall trend test
    mk_result <- try(Kendall::MannKendall(working_data$value), silent = TRUE)
    if (inherits(mk_result, "try-error")) {
      mk_tau <- NA
      mk_pval <- NA
    } else {
      mk_tau <- as.numeric(mk_result$tau)
      mk_pval <- as.numeric(mk_result$sl)  # sl = significance level (p-value)
    }

    # Calculate mean and median
    mean_val <- mean(working_data$value, na.rm = TRUE)
    median_val <- median(working_data$value, na.rm = TRUE)

    # Store results with new naming convention
    results[[paste0(col, "_senn_slp")]] <- sen_slope
    results[[paste0(col, "_linear_slp")]] <- linear_slope
    results[[paste0(col, "_spearman_rho")]] <- spearman_rho
    results[[paste0(col, "_spearman_pval")]] <- spearman_pval
    results[[paste0(col, "_mk_rho")]] <- mk_tau
    results[[paste0(col, "_mk_pval")]] <- mk_pval
    results[[paste0(col, "_mean")]] <- mean_val
    results[[paste0(col, "_median")]] <- median_val
  }

  # Convert to data frame
  return(as.data.frame(results))
}






#################################################################
# start streamflow signature analysis


calculate_flow_vols_by_year = function(streamflow_data){
  # Ensure required columns exist
  required_cols <- c("water_year", "Q", "month", "dowy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }

  # Filter out water years with insufficient data (< MIN_NONA_DAYS_ANNUAL non-NA days)
  years <- unique(streamflow_data$water_year)
  valid_years <- c()
  for (yr in years) {
    yr_data <- streamflow_data[streamflow_data$water_year == yr, ]
    if (sum(!is.na(yr_data$Q)) >= MIN_NONA_DAYS_ANNUAL) {
      valid_years <- c(valid_years, yr)
    }
  }
  streamflow_data <- streamflow_data[streamflow_data$water_year %in% valid_years, ]

  # Calculate annual total flow (mm) by summing daily Q (mm/day).
  # na.rm=TRUE sums available days rather than discarding years with minor gaps.
  # Safe here because the 250-day filter above guarantees substantial valid data per year.
  # Trade-off: if a year has a few missing days, the total is slightly underestimated
  # rather than lost entirely — standard practice for hydrological time series.
  annual_totals <- aggregate(Q ~ water_year, data=streamflow_data, FUN=sum, na.rm=TRUE)

  # Check if we have any valid annual data
  if (nrow(annual_totals) == 0) {
    # Return empty result with correct structure
    return(generate_stats(data.frame(water_year = numeric()),
                          value_cols = character(),
                          year_col = "water_year"))
  }

  names(annual_totals)[2] <- "Qann"

  # Calculate seasonal total flow (mm)
  winter <- streamflow_data[streamflow_data$month %in% c(12, 1, 2), ]
  spring <- streamflow_data[streamflow_data$month %in% c(3, 4, 5), ]
  summer <- streamflow_data[streamflow_data$month %in% c(6, 7, 8), ]
  fall <- streamflow_data[streamflow_data$month %in% c(9, 10, 11), ]

  # Seasonal sums also use na.rm=TRUE. Edge case: a season could have all NAs
  # while the year passes the 250-day filter (e.g., all winter days missing).
  # In that rare case sum returns 0 rather than NA — acceptable given the annual
  # completeness guarantee and standard hydrological practice.
  winter_totals <- aggregate(Q ~ water_year, data=winter, FUN=sum, na.rm=TRUE)
  spring_totals <- aggregate(Q ~ water_year, data=spring, FUN=sum, na.rm=TRUE)
  summer_totals <- aggregate(Q ~ water_year, data=summer, FUN=sum, na.rm=TRUE)
  fall_totals <- aggregate(Q ~ water_year, data=fall, FUN=sum, na.rm=TRUE)

  # Rename columns for seasonal totals
  if (nrow(winter_totals) > 0) names(winter_totals)[2] <- "Qwin"
  if (nrow(spring_totals) > 0) names(spring_totals)[2] <- "Qspr"
  if (nrow(summer_totals) > 0) names(summer_totals)[2] <- "Qsum"
  if (nrow(fall_totals) > 0) names(fall_totals)[2] <- "Qfal"

  # After seasonal aggregations, check if any are empty
  if (nrow(winter_totals) == 0) winter_totals <- data.frame(water_year=numeric(), Qwin=numeric())
  if (nrow(spring_totals) == 0) spring_totals <- data.frame(water_year=numeric(), Qspr=numeric())
  if (nrow(summer_totals) == 0) summer_totals <- data.frame(water_year=numeric(), Qsum=numeric())
  if (nrow(fall_totals) == 0) fall_totals <- data.frame(water_year=numeric(), Qfal=numeric())
  
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
    q95_q10 <- data.frame(water_year = annual_totals$water_year, `Q95-Q10` = NA)
  }

  # Merge all metrics into a single data frame
  # Start with annual totals
  all_metrics <- annual_totals

  # Helper function for safe merging
  safe_merge <- function(df1, df2, by_col) {
    if (nrow(df2) == 0) {
      # If df2 is empty, just return df1
      return(df1)
    }
    return(merge(df1, df2, by=by_col, all.x=TRUE))
  }
  
  # Merge seasonal totals
  all_metrics <- safe_merge(all_metrics, winter_totals, "water_year")
  all_metrics <- safe_merge(all_metrics, spring_totals, "water_year")
  all_metrics <- safe_merge(all_metrics, summer_totals, "water_year")
  all_metrics <- safe_merge(all_metrics, fall_totals, "water_year")
  
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
    # Return empty result with correct structure (8 stats per metric)
    return(data.frame(
      FDCall_senn_slp = NA, FDCall_linear_slp = NA, FDCall_spearman_rho = NA, FDCall_spearman_pval = NA,
      FDCall_mk_rho = NA, FDCall_mk_pval = NA, FDCall_mean = NA, FDCall_median = NA,
      FDC90th_senn_slp = NA, FDC90th_linear_slp = NA, FDC90th_spearman_rho = NA, FDC90th_spearman_pval = NA,
      FDC90th_mk_rho = NA, FDC90th_mk_pval = NA, FDC90th_mean = NA, FDC90th_median = NA,
      FDCmid_senn_slp = NA, FDCmid_linear_slp = NA, FDCmid_spearman_rho = NA, FDCmid_spearman_pval = NA,
      FDCmid_mk_rho = NA, FDCmid_mk_pval = NA, FDCmid_mean = NA, FDCmid_median = NA
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
      fdc$log_flow <- log10(fdc$flow + FDC_FLOW_FLOOR)  # Small constant to handle zeros
      
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
    
    
    # Need sufficient non-NA data points for the year
    if (sum(!is.na(year_data$Q)) < 30) {
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

    # Guard against division by zero (all flows are zero)
    total_q <- sum(q_values, na.rm=TRUE)
    if (total_q == 0) next

    # Calculate RB index
    rb_index <- sum(q_diff, na.rm=TRUE) / total_q
    
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
  
  # Create a data frame to store day of water year when cumulative flow reaches each percentile
  years <- unique(streamflow_data$water_year)
  timing_by_year <- data.frame(water_year = years)
  
  # Define percentiles
  percentiles <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  
  # Initialize columns for all metrics
  for (p in percentiles) {
    timing_by_year[[paste0("D", p, "_day")]] <- NA
  }
  timing_by_year$D25_to_D75 <- NA
  timing_by_year$Dmax <- NA
  
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
    
    # Replace NAs with 0 before cumsum (conservative: treats missing as zero flow)
    # R's cumsum() propagates NAs, which would corrupt all subsequent values
    q_for_cumsum <- year_data$Q
    q_for_cumsum[is.na(q_for_cumsum)] <- 0
    year_data$cum_flow <- cumsum(q_for_cumsum)

    # Calculate cumulative flow as percentage of total
    year_data$cum_pct <- (year_data$cum_flow / total_flow) * 100
    
    # For each percentile, find the first day when cumulative flow exceeds the threshold
    for (p in percentiles) {
      # Find days where cumulative percentage exceeds the threshold
      above_threshold <- which(year_data$cum_pct >= p)
      
      # If there are days above threshold, take the first one
      if (length(above_threshold) > 0) {
        timing_by_year[timing_by_year$water_year == yr, paste0("D", p, "_day")] <- 
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
      timing_by_year[timing_by_year$water_year == yr, "D25_to_D75"] <- day_75 - day_25
    }
    
    # Calculate Dmax (day of maximum discharge)
    max_Q_idx <- which.max(year_data$Q)
    if (length(max_Q_idx) > 0) {
      timing_by_year[timing_by_year$water_year == yr, "Dmax"] <- year_data$dowy[max_Q_idx]
    }
  }
  
  # Define which columns to calculate statistics for
  metric_columns <- c(paste0("D", percentiles, "_day"), "D25_to_D75", "Dmax")
  
  # Use generate_stats to calculate all statistics at once
  result <- generate_stats(timing_by_year, value_cols = metric_columns, year_col = "water_year")
  
  # Add timing_by_year as an attribute to the result
  attr(result, "timing_by_year") <- timing_by_year
  
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
    # Interpolate small NA gaps instead of removing them,
    # to avoid creating false adjacencies between non-adjacent days
    if (any(is.na(flow_vector))) {
      non_na_idx <- which(!is.na(flow_vector))
      if (length(non_na_idx) < 3) return(0)
      flow_clean <- approx(non_na_idx, flow_vector[non_na_idx],
                            xout = seq_len(length(flow_vector)), rule = 2)$y
    } else {
      flow_clean <- flow_vector
    }
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

    # Skip years with insufficient valid data (minimum days from config.R)
    if (sum(!is.na(year_data$Q)) < MIN_NONA_DAYS_ANNUAL) {
      next
    }

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
  # Calculate annual runoff ratio; require minimum 10mm annual PPT (consistent with elasticity)
  annual_totals$annual_runoff_ratio <- ifelse(annual_totals$PPT > ELASTICITY_MIN_ANNUAL_PPT,
                                              annual_totals$Q / annual_totals$PPT,
                                              NA)

  # Calculate seasonal totals and ratios
  # Seasonal minimum PPT threshold: 1mm (avoids extreme ratios from near-zero PPT)
  seasonal_min_ppt <- 1

  # Winter (December, January, February)
  winter <- streamflow_data[streamflow_data$month %in% c(12, 1, 2), ]
  winter_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=winter, FUN=sum, na.rm=TRUE)
  winter_totals$winter_runoff_ratio <- ifelse(winter_totals$PPT > seasonal_min_ppt,
                                              winter_totals$Q / winter_totals$PPT,
                                              NA)

  # Spring (March, April, May)
  spring <- streamflow_data[streamflow_data$month %in% c(3, 4, 5), ]
  spring_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=spring, FUN=sum, na.rm=TRUE)
  spring_totals$spring_runoff_ratio <- ifelse(spring_totals$PPT > seasonal_min_ppt,
                                              spring_totals$Q / spring_totals$PPT,
                                              NA)

  # Summer (June, July, August)
  summer <- streamflow_data[streamflow_data$month %in% c(6, 7, 8), ]
  summer_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=summer, FUN=sum, na.rm=TRUE)
  summer_totals$summer_runoff_ratio <- ifelse(summer_totals$PPT > seasonal_min_ppt,
                                              summer_totals$Q / summer_totals$PPT,
                                              NA)

  # Fall (September, October, November)
  fall <- streamflow_data[streamflow_data$month %in% c(9, 10, 11), ]
  fall_totals <- aggregate(cbind(Q, PPT) ~ water_year, data=fall, FUN=sum, na.rm=TRUE)
  fall_totals$fall_runoff_ratio <- ifelse(fall_totals$PPT > seasonal_min_ppt,
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
  # Check if required columns exist (dowy is used for sorting, not doy)
  required_cols <- c("water_year", "Q", "dowy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Function to apply Eckhardt filter
  eckhardt_filter <- function(Q, BFImax = 0.8, a = 0.98) {
    n <- length(Q)
    baseflow <- numeric(n)

    # Initialize matching Python/Julia: min(BFImax * Q[1], Q[1]) if valid and > 0
    if (!is.na(Q[1]) && Q[1] > 0) {
      baseflow[1] <- min(Q[1] * BFImax, Q[1])
    } else {
      baseflow[1] <- 0
    }

    for (i in 2:n) {
      if (is.na(Q[i])) {
        # Forward-fill: carry previous baseflow through gap.
        # Baseflow (groundwater component) changes slowly, so this is
        # hydrologically defensible for short gaps. Prevents NaN cascade
        # that would lose all post-gap baseflow, creating a numerator/
        # denominator mismatch in BFI calculation.
        # Matches Python and Julia implementations (Round 5 fix).
        baseflow[i] <- baseflow[i - 1]
        next
      }
      # Normal Eckhardt filter equation — baseflow[i-1] is always valid
      # because forward-fill guarantees no NA propagation
      numerator <- (1 - BFImax) * a * baseflow[i-1] + (1 - a) * BFImax * Q[i]
      denominator <- 1 - a * BFImax
      baseflow[i] <- numerator / denominator
      baseflow[i] <- min(baseflow[i], Q[i])
      baseflow[i] <- max(baseflow[i], 0)
    }

    return(baseflow)
  }
  
  # Function to apply Lyne-Hollick filter
  lyne_hollick_filter <- function(Q, alpha = 0.925, passes = 2) {
    n <- length(Q)

    # The input signal for each pass; first pass uses original Q,
    # subsequent passes use baseflow from the previous pass
    input_signal <- Q

    for (pass in 1:passes) {
      # Forward pass — filter the current input signal
      qf_forward <- numeric(n)
      qf_forward[1] <- 0

      for (i in 2:n) {
        if (!is.na(input_signal[i]) && !is.na(input_signal[i-1]) && !is.na(qf_forward[i-1])) {
          qf_forward[i] <- alpha * qf_forward[i-1] + ((1 + alpha) / 2) * (input_signal[i] - input_signal[i-1])
          # Constrain quickflow to [0, input_signal]
          qf_forward[i] <- max(0, min(qf_forward[i], input_signal[i]))
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
          # Constrain quickflow to [0, input_signal]
          qf_backward[i] <- max(0, min(qf_backward[i], input_signal[i]))
        } else {
          qf_backward[i] <- NA
        }
      }

      # Baseflow from this pass becomes input for the next pass
      input_signal <- input_signal - qf_backward
      input_signal[input_signal < 0] <- 0
      input_signal[is.na(Q)] <- NA
    }

    # After all passes, input_signal holds the final baseflow
    return(input_signal)
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
        
        # Check next min_length consecutive decreases (requires min_length+1 points)
        for (j in 0:(min_length-1)) {
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
  
  # Fit sinusoidal model to recession log(a) values as a function of day-of-year.
  #
  # Model: log(a) = B1*sin(2*pi*doy/365) + B2*cos(2*pi*doy/365) + C
  #   Amplitude = sqrt(B1^2 + B2^2)
  #   Phase (days) = atan2(-B2, B1) * 365 / (2*pi)
  #   Minimum doy = phase + 273.75 days (3/4 cycle offset from maximum)
  #
  # @param doy_values Numeric vector of day-of-year values

  # @param log_a_values Numeric vector of log(a) recession rate values
  # @return List with amplitude (dimensionless) and minimum_doy (1-365)
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

        # Calculate log(a) using b_pointcloud (not median_b from per-event fits)
        # For each point: log(a) = log(-dQ/dt) - b*log(Q)
        log_a_values_pc <- log(all_dQ_dt) - b_pointcloud * log(all_Q)
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
  
  # Check if we have enough recession events overall (minimum 25 per RECESSION_MIN_EVENTS in config.R)
  total_events <- length(all_recession_events)
  if (total_events < RECESSION_MIN_EVENTS) {
    # Not enough data, return all NAs with correct 8-statistic schema
    col_names <- c()
    for (sig in signatures_with_stats) {
      for (stat in STAT_SUFFIXES) {
        col_names <- c(col_names, paste0(sig, stat))
      }
    }
    col_names <- c(col_names, seasonality_signatures)
    result <- data.frame(matrix(NA, nrow=1, ncol=length(col_names)))
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






# -----------------------------------------------------------------------------
# Helper function to generate streamflow data.table from Caravan NetCDF files
# -----------------------------------------------------------------------------
generate_streamflow_dt_caravan <- function(nc_file_path, 
                                           min_num_years_data = 20,
                                           start_date_filter = as.Date("1900-01-01"), 
                                           end_date_filter = as.Date("2024-12-31")) {
  
  if (!file.exists(nc_file_path)) {
    warning("NetCDF file not found: ", nc_file_path)
    return(NULL)
  }
  
  tryCatch({
    nc_data <- nc_open(nc_file_path)
    
    # Extract streamflow data
    streamflow_raw <- ncvar_get(nc_data, "streamflow")
    ppt_raw <- ncvar_get(nc_data, "total_precipitation_sum")
    swe_raw <- ncvar_get(nc_data, "snow_depth_water_equivalent_mean")
    
    # Extract date dimension
    time_raw <- ncvar_get(nc_data, "date")
    time_units <- ncatt_get(nc_data, "date", "units")$value
    
    nc_close(nc_data)
    
    if (is.null(streamflow_raw) || is.null(time_raw)) {
      warning("Streamflow or date variable not found in: ", nc_file_path)
      return(NULL)
    }
    
    # Convert time to Date objects
    origin_date_str <- sub("days since ", "", time_units)
    origin_date <- as.Date(origin_date_str, format="%Y-%m-%d %H:%M:%S")
    if (is.na(origin_date)) {
      origin_date <- as.Date(sub("days since ", "", time_units))
    }
    if (is.na(origin_date)) {
      warning("Could not parse origin date from units: ", time_units, " in file: ", nc_file_path)
      return(NULL)
    }
    
    dates <- as.Date(time_raw, origin = origin_date)
    
    stream_dt <- data.table(Date = dates,
                            Q = as.numeric(streamflow_raw),
                            PPT = as.numeric(ppt_raw),
                            SWE = as.numeric(swe_raw))
    
    # Filter by overall start and end date
    stream_dt <- stream_dt[Date >= start_date_filter & Date <= end_date_filter]
    
    if (nrow(stream_dt) == 0) {
      message("No data within the specified date range for: ", basename(nc_file_path))
      return(NULL)
    }
    
    # Add year, month, doy
    stream_dt[, year := year(Date)]
    stream_dt[, month := month(Date)]
    stream_dt[, doy := yday(Date)]
    
    # Add water year information
    wy_info <- calculate_water_year_info(stream_dt$Date)
    stream_dt[, water_year := wy_info$water_year]
    stream_dt[, dowy := wy_info$dowy]
    
    # Check for minimum number of years of data
    if (length(unique(stream_dt$year)) < min_num_years_data) {
      message("Insufficient years of data (", length(unique(stream_dt$year)), 
              " years) after date filtering for: ", basename(nc_file_path), 
              " (min_num_years_data: ", min_num_years_data, ")")
      return(NULL)
    }
    
    stream_dt[, Q := as.numeric(Q)]
    
    return(stream_dt)
    
  }, error = function(e) {
    warning("Error processing NetCDF file ", nc_file_path, ": ", e$message)
    return(NULL)
  })
}


# -----------------------------------------------------------------------------
# Main wrapper function to process Caravan data for a given data_project
# -----------------------------------------------------------------------------
process_caravan_gages <- function(data_project_arg, caravan_base_dir, 
                                  min_num_years, start_date, end_date, 
                                  min_Q_value_and_days, output_file) {
  
  # Construct path to the data_project directory
  project_data_path <- file.path(caravan_base_dir, "timeseries", "netcdf", data_project_arg)
  
  if (!dir.exists(project_data_path)) {
    stop("Caravan project data directory not found: ", project_data_path)
  }
  
  # List all NetCDF files for the data_project
  nc_files <- list.files(project_data_path, pattern = paste0("^", data_project_arg, "_.*\\.nc$"), full.names = TRUE)
  
  if (length(nc_files) == 0) {
    cat("No NetCDF files found for data_project '", data_project_arg, "' in ", project_data_path, "\n")
    return(data.table()) # Return empty data.table
  }
  
  # Check if output file exists; if so, read it.
  if (file.exists(output_file)) {
    summary_output <- fread(output_file, colClasses = list(character = c("watershed_id", "data_project")))
    # Ensure types after reading, especially if integer64 might occur
    if ("watershed_id" %in% names(summary_output)) summary_output[, watershed_id := as.character(watershed_id)]
    if ("data_project" %in% names(summary_output)) summary_output[, data_project := as.character(data_project)]
    cat("Loaded existing summary data with", nrow(summary_output), "rows from '", output_file, "'\n")
  } else {
    summary_output <- data.table(
      watershed_id = character(),
      data_project = character(),
      latitude = numeric(),      # Will be NA for Caravan unless metadata is sourced elsewhere
      longitude = numeric(),     # Will be NA for Caravan
      basin_area = numeric()     # Will be NA for Caravan
      # Metric columns will be added by rbind with fill=TRUE
    )
    cat("Created new summary data table for '", output_file, "'\n")
    # fwrite(summary_output, output_file) # Write header only if file is new
  }
  
  # Create a unique identifier for already processed watersheds in the current output file
  if (nrow(summary_output) > 0 && all(c("watershed_id", "data_project") %in% names(summary_output))) {
    processed_identifiers <- paste(summary_output$watershed_id, summary_output$data_project, sep = "_")
  } else {
    processed_identifiers <- character(0)
  }
  
  # Process each NetCDF file (watershed)
  for (i in 1:length(nc_files)) {
    nc_file_path <- nc_files[i]
    filename_base <- tools::file_path_sans_ext(basename(nc_file_path))
    
    # Extract watershed_id from filename (e.g., "camels_01013500" -> "01013500")
    watershed_id <- sub(paste0("^", data_project_arg, "_"), "", filename_base)
    
    current_identifier <- paste(watershed_id, data_project_arg, sep = "_")
    if (current_identifier %in% processed_identifiers) {
      cat("Skipping watershed", watershed_id, "(data_project:", data_project_arg, ") as it's already in the output\n")
      # note: there are redundancies between the camels and hysets datasets, so redundancy is expected
      next
    }
    
    cat("Processing watershed", watershed_id, "(data_project:", data_project_arg, ") (", i, "of", length(nc_files), ")\n")
    
    tryCatch({
      streamflow_data <- generate_streamflow_dt_caravan(
        nc_file_path = nc_file_path,
        min_num_years_data = min_num_years, # Pass the overall min_num_years for the internal check
        start_date_filter = start_date,
        end_date_filter = end_date
      )
      
      if (is.null(streamflow_data) || !is.data.frame(streamflow_data) || nrow(streamflow_data) == 0) {
        cat("No valid streamflow data for watershed", watershed_id, "(data_project:", data_project_arg, ")\n")
        next
      }
      
      # Ensure required columns from generate_streamflow_dt_caravan are present
      required_cols_from_generation <- c("Date", "Q", "year", "month", "doy", "water_year", "dowy")
      if (!all(required_cols_from_generation %in% colnames(streamflow_data))) {
        # Add water year information if missing
        if ("Date" %in% colnames(streamflow_data) && !all(c("water_year", "dowy") %in% colnames(streamflow_data))) {
          wy_info <- calculate_water_year_info(streamflow_data$Date)
          streamflow_data$water_year <- wy_info$water_year
          streamflow_data$dowy <- wy_info$dowy
        }
        
        if (!all(required_cols_from_generation %in% colnames(streamflow_data))) {
          cat("Missing required columns from generate_streamflow_dt_caravan for watershed", watershed_id, "\n")
          next
        }
      }
      
      # Apply min_Q_value_and_days filter (similar to process_gages)
      years_to_use <- NULL
      for (this_year in unique(streamflow_data$year)) {
        test_year <- subset(streamflow_data, year == this_year)
        # Ensure Q is numeric for comparison
        if(!is.numeric(test_year$Q)) test_year$Q <- as.numeric(test_year$Q)
        
        # Handle potential NAs in Q before comparison
        valid_q_values <- test_year$Q[!is.na(test_year$Q)]
        nonzero_rows <- which(valid_q_values > min_Q_value_and_days[1])
        
        if (length(nonzero_rows) >= min_Q_value_and_days[2]) { # Use >= to match intent
          years_to_use <- c(years_to_use, this_year)
        }
      }
      
      # This check for min_num_years is applied *after* the min_Q_value_and_days filter per year
      if (length(years_to_use) < min_num_years) { # Use <, if 20 years are required, 19 is not enough
        cat("Insufficient years (", length(years_to_use), ") with valid data after Q filter for watershed", watershed_id, 
            "(data_project:", data_project_arg, "). Required:", min_num_years, "\n")
        next
      }
      
      streamflow_data_filtered <- streamflow_data[streamflow_data$year %in% years_to_use, ]
      
      if (nrow(streamflow_data_filtered) == 0) {
        cat("No data remaining after year filtering for watershed", watershed_id, "\n")
        next
      }
      
      # Calculate metrics using existing functions
      # Ensure streamflow_data_filtered has the columns expected by these functions
      # (Date, Q, year, month, doy)
      
      cat("Calculating metrics for watershed", watershed_id, "...\n")
      metrics_flow_vols <- calculate_flow_vols_by_year(streamflow_data_filtered)
      metrics_fdc_trends <- analyze_fdc_trends_from_streamflow(streamflow_data_filtered)
      metrics_flashiness <- analyze_flashiness_trends(streamflow_data_filtered)
      metrics_flow_timing <- analyze_flow_timing_trends(streamflow_data_filtered)
      metrics_pulses <- calculate_pulse_metrics(streamflow_data_filtered)
      metrics_QtoPPT <- analyze_Q_PPT_relationships(streamflow_data_filtered)
      metrics_baseflow <- analyze_baseflow_indices(streamflow_data_filtered)
      metrics_recession <- analyze_recession_parameters(streamflow_data_filtered)
      # Add other metric calculations here if needed
      
      # Create a row for this watershed
      watershed_row <- data.table(
        watershed_id = watershed_id,
        data_project = data_project_arg,
        latitude = NA_real_,  # Caravan .nc files for timeseries don't typically store this
        longitude = NA_real_, # Caravan .nc files for timeseries don't typically store this
        basin_area = NA_real_,# Caravan .nc files for timeseries don't typically store this
        num_years = length(years_to_use),
        start_year = min(years_to_use),
        end_year = max(years_to_use)
      )
      
      # Combine base watershed info with calculated metrics
      # Ensure metric results are data.tables or can be coerced
      watershed_row <- cbind(watershed_row, 
                             as.data.table(metrics_flow_vols),
                             as.data.table(metrics_fdc_trends),
                             as.data.table(metrics_flashiness),
                             as.data.table(metrics_flow_timing),
                             as.data.table(metrics_pulses),
                             as.data.table(metrics_QtoPPT),
                             as.data.table(metrics_baseflow),
                             as.data.table(metrics_recession))
      
      summary_output <- rbind(summary_output, watershed_row, fill = TRUE)
      fwrite(summary_output, output_file) # Write after each successful processing
      cat("Successfully processed and saved watershed", watershed_id, "(data_project:", data_project_arg, ")\n")
      
    }, error = function(e) {
      cat("Error processing watershed", watershed_id, "(data_project:", data_project_arg, "):", e$message, "\n")
    })
  }
  
  cat("Finished processing all files for data_project '", data_project_arg, "'.\n")
  return(summary_output)
}




####################################################################
## not yet implemented

get_era5land_for_basins <- function(basin_ids, basin_polygons, 
                                    start_year = 1980, end_year = 2021,
                                    variables = c("2m_temperature", "total_precipitation"),
                                    user_email = Sys.getenv("ECMWF_USER"),
                                    api_key = Sys.getenv("ECMWF_KEY")) {
  
  # Check credentials
  if (user_email == "" || api_key == "") {
    stop("ECMWF credentials not found. Please set them using wf_set_key() or environment variables.")
  }
  
  # 1. Create a spatial union of all basins
  basin_subset <- basin_polygons[basin_polygons$HYBAS_ID %in% basin_ids,]
  if (nrow(basin_subset) == 0) {
    stop("No matching basins found")
  }
  
  # Union all basins into a single polygon
  basin_union <- st_union(basin_subset)
  
  # 2. Get the bounding box with a small buffer
  bbox <- st_bbox(basin_union)
  bbox <- bbox + c(-0.2, -0.2, 0.2, 0.2)  # Add buffer
  
  # 3. Set up ERA5-Land request
  request <- list(
    product_type = "reanalysis",
    format = "netcdf",
    variable = variables,
    year = as.character(start_year:end_year),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    time = "00:00",
    area = c(bbox["ymax"], bbox["xmin"], bbox["ymin"], bbox["xmax"]),
    dataset_short_name = "reanalysis-era5-land",
    target = "era5land_data.nc"
  )
  
  # 4. Submit request to ECMWF
  file <- wf_request(
    user = user_email,
    request = request,
    transfer = TRUE,
    path = "."
  )
  
  # 5. Process the downloaded NetCDF file
  era5_raster <- rast(file)
  
  # 6. Extract data for each basin and calculate basin-average values
  basin_climate <- data.table()
  
  for (id in basin_ids) {
    # Get the basin polygon
    basin_poly <- basin_subset[basin_subset$HYBAS_ID == id,]
    
    # Extract values for this basin
    basin_values <- terra::extract(era5_raster, vect(basin_poly), fun = mean, na.rm = TRUE)
    
    # Convert to data.table with basin ID
    basin_dt <- as.data.table(basin_values)
    basin_dt[, HYBAS_ID := id]
    
    # Append to main results
    basin_climate <- rbind(basin_climate, basin_dt)
  }
  
  # 7. Post-process the data (convert units, calculate derived variables)
  # This will depend on the specific variables you requested
  
  return(basin_climate)
}


#################################################################################
# read in and plot full time series data

# Load required libraries
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


################################################################################
# DAYMET CLIMATE DATA INTEGRATION
################################################################################

#' Convert Daymet ZIP of CSVs to a single Parquet file with proper Date column
#'
#' One-time conversion that:
#' - Reads all yearly CSVs from the ZIP
#' - Reconstructs the Date column (CSVs only have month/year, not day)
#' - Combines into a single parquet with snappy compression
#'
#' @param daymet_zip_path Path to daymet_1980_2023.zip
#' @param output_parquet_path Path for output parquet file
#' @param years Vector of years to process (default: all 1980-2023)
#' @param max_sites Maximum number of sites to process (NULL for all, useful for testing)
#' @return Invisible path to output parquet
convert_daymet_zip_to_parquet <- function(daymet_zip_path,
                                           output_parquet_path,
                                           years = 1980:2023,
                                           max_sites = NULL) {

  # Validate inputs
  validate_file_exists(daymet_zip_path, "daymet_zip_path", required_ext = "zip",
                       context = "convert_daymet_zip_to_parquet")

  log_info("Starting Daymet ZIP to parquet conversion", context = "convert_daymet")
  log_info("Input:", daymet_zip_path, context = "convert_daymet")
  log_info("Output:", output_parquet_path, context = "convert_daymet")
  log_info("Years:", min(years), "to", max(years), context = "convert_daymet")

  # List files in the ZIP
  zip_contents <- unzip(daymet_zip_path, list = TRUE)
  csv_files <- zip_contents$Name[grepl("\\.csv$", zip_contents$Name)]

  # Filter to requested years
  year_pattern <- paste0("daymet_(", paste(years, collapse = "|"), ")\\.csv$")
  csv_files <- csv_files[grepl(year_pattern, csv_files)]

  if (length(csv_files) == 0) {
    stop("No matching CSV files found in ZIP for requested years")
  }

  log_info("Found", length(csv_files), "CSV files to process", context = "convert_daymet")

  # Process each year and combine
  all_data <- list()

  for (i in seq_along(csv_files)) {
    csv_name <- csv_files[i]
    year <- as.integer(gsub(".*daymet_(\\d{4})\\.csv$", "\\1", csv_name))

    log_info("Processing", csv_name, "(", i, "/", length(csv_files), ")",
             context = "convert_daymet")

    # Extract CSV from ZIP to temp directory, then read with fread
    # (fread doesn't work directly with unz() connections)
    temp_dir <- tempdir()
    dt <- tryCatch({
      # Extract single file from ZIP
      unzip(daymet_zip_path, files = csv_name, exdir = temp_dir, junkpaths = TRUE)
      extracted_path <- file.path(temp_dir, basename(csv_name))

      # Read with fread (fast)
      result <- fread(extracted_path, showProgress = FALSE)

      # Clean up temp file
      unlink(extracted_path)

      result
    }, error = function(e) {
      stop(paste0("Error reading ", csv_name, ": ", e$message))
    })

    # Limit sites if requested (for testing)
    if (!is.null(max_sites)) {
      unique_sites <- unique(dt$site_id)
      if (length(unique_sites) > max_sites) {
        sites_to_keep <- unique_sites[1:max_sites]
        dt <- dt[site_id %in% sites_to_keep]
      }
    }

    # Reconstruct the day column
    # Data is organized as: all days for site 1, then all days for site 2, etc.
    # Within each site/year/month, rows are sequential days
    dt[, day := seq_len(.N), by = .(site_id, year, month)]

    # Create proper Date column
    dt[, Date := as.Date(paste(year, month, day, sep = "-"))]

    # Validate: check we have expected number of days per site
    days_per_site <- dt[, .N, by = site_id]
    expected_days <- ifelse(year %% 4 == 0 & (year %% 100 != 0 | year %% 400 == 0), 366, 365)

    bad_sites <- days_per_site[N != expected_days]
    if (nrow(bad_sites) > 0) {
      log_warn("Year", year, ":", nrow(bad_sites), "sites have unexpected day count",
               context = "convert_daymet")
    }

    # Select and order columns for output
    dt <- dt[, .(site_id, Date, prcp, tmin, tmax, swe, vp, srad)]

    all_data[[i]] <- dt

    # Clear memory periodically
    if (i %% 10 == 0) {
      gc(verbose = FALSE)
    }
  }

  # Combine all years
  log_info("Combining all years...", context = "convert_daymet")
  combined <- rbindlist(all_data, use.names = TRUE)

  # Sort by site_id and Date
  setorder(combined, site_id, Date)

  # Check for duplicates
  n_rows <- nrow(combined)
  n_unique <- nrow(unique(combined[, .(site_id, Date)]))
  if (n_rows != n_unique) {
    log_warn("Found", n_rows - n_unique, "duplicate site_id + Date combinations",
             context = "convert_daymet")
  }

  # Summary statistics
  n_sites <- length(unique(combined$site_id))
  date_range <- range(combined$Date)
  log_info("Combined data:", n_rows, "rows,", n_sites, "sites",
           context = "convert_daymet")
  log_info("Date range:", as.character(date_range[1]), "to", as.character(date_range[2]),
           context = "convert_daymet")

  # Ensure output directory exists
  output_dir <- dirname(output_parquet_path)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write to parquet with compression
  log_info("Writing parquet file...", context = "convert_daymet")
  arrow::write_parquet(combined, output_parquet_path, compression = "snappy")

  # Report file size
  file_size_mb <- file.info(output_parquet_path)$size / 1024 / 1024
  log_info("Parquet file created:", round(file_size_mb, 1), "MB",
           context = "convert_daymet")

  invisible(output_parquet_path)
}


#' Load Daymet climate data for a specific gage
#'
#' Reads from the pre-processed parquet file and filters to a single gage.
#'
#' @param gage_id USGS gage ID (e.g., "01011000")
#' @param daymet_parquet_path Path to daymet_1980_2023.parquet
#' @param start_year First year to include (default: 1980)
#' @param end_year Last year to include (default: 2023)
#' @return data.table with Date, prcp, tmin, tmax, swe, vp, srad
load_daymet_for_gage <- function(gage_id,
                                  daymet_parquet_path,
                                  start_year = 1980,
                                  end_year = 2023) {

  validate_file_exists(daymet_parquet_path, "daymet_parquet_path",
                       required_ext = "parquet", context = "load_daymet_for_gage")

  log_debug("Loading Daymet data for gage:", gage_id, context = "load_daymet_for_gage")

  # Read parquet with filtering
  # Use arrow's filter pushdown for efficiency
  df <- arrow::read_parquet(
    daymet_parquet_path,
    as_data_frame = FALSE
  )

  # Filter by site_id and year range
  start_date <- as.Date(paste0(start_year, "-01-01"))
  end_date <- as.Date(paste0(end_year, "-12-31"))

  result <- df |>
    dplyr::filter(site_id == gage_id) |>
    dplyr::filter(Date >= start_date & Date <= end_date) |>
    dplyr::select(-site_id) |>
    dplyr::collect()

  result <- as.data.table(result)

  if (nrow(result) == 0) {
    log_warn("No Daymet data found for gage:", gage_id, context = "load_daymet_for_gage")
  } else {
    log_debug("Loaded", nrow(result), "days of climate data for gage:", gage_id,
              context = "load_daymet_for_gage")
  }

  return(result)
}


#' Integrate Daymet climate data with streamflow data
#'
#' Joins Daymet climate variables to a streamflow data.table based on Date.
#'
#' @param streamflow_data data.table with at least a Date column
#' @param gage_id USGS gage ID for fetching Daymet data
#' @param daymet_parquet_path Path to daymet parquet file
#' @return streamflow_data with added prcp, tmin, tmax columns
integrate_daymet_with_streamflow <- function(streamflow_data,
                                              gage_id,
                                              daymet_parquet_path) {

  ctx <- "integrate_daymet_with_streamflow"

  # Validate streamflow data has Date column
  if (!"Date" %in% names(streamflow_data)) {
    stop("streamflow_data must have a 'Date' column")
  }

  # Get year range from streamflow data
  start_year <- year(min(streamflow_data$Date))
  end_year <- year(max(streamflow_data$Date))

  # Load Daymet data
  daymet_data <- load_daymet_for_gage(
    gage_id = gage_id,
    daymet_parquet_path = daymet_parquet_path,
    start_year = start_year,
    end_year = end_year
  )

  if (nrow(daymet_data) == 0) {
    log_warn("No Daymet data available for gage", gage_id, "- returning original data",
             context = ctx)
    return(streamflow_data)
  }

  # Join on Date
  result <- merge(streamflow_data, daymet_data, by = "Date", all.x = TRUE)

  # Report coverage
  n_matched <- sum(!is.na(result$prcp))
  n_total <- nrow(result)
  coverage <- n_matched / n_total * 100

  if (coverage < 95) {
    log_warn("Daymet coverage for gage", gage_id, "is only", round(coverage, 1), "%",
             context = ctx)
  } else {
    log_debug("Daymet coverage:", round(coverage, 1), "% for gage", gage_id,
              context = ctx)
  }

  # Rename prcp to PPT for compatibility with analyze_Q_PPT_relationships()
  if ("prcp" %in% names(result)) {
    setnames(result, "prcp", "PPT")
  }

  return(result)
}


################################################################################
# CLIMATE-DEPENDENT STREAMFLOW SIGNATURES
################################################################################

#' Calculate streamflow elasticity (Sawicz et al. 2011)
#'
#' Elasticity measures how sensitive streamflow is to precipitation changes.
#' E = median((dQ/dP) / (Q_mean/P_mean))
#'
#' @param streamflow_data data.table with water_year, Q, prcp columns
#' @param rolling_window Number of years for rolling window (NULL for static only)
#' @return List with elasticity statistics
calculate_streamflow_elasticity <- function(streamflow_data,
                                             rolling_window = ELASTICITY_WINDOW_YEARS) {

  ctx <- "calculate_streamflow_elasticity"

  # Validate required columns
  required_cols <- c("water_year", "Q", "PPT")
  validate_columns(streamflow_data, required_cols, "streamflow_data", context = ctx)

  # Aggregate to annual totals with data quality check
  annual <- streamflow_data[, .(
    Q_annual = sum(Q, na.rm = TRUE),
    P_annual = sum(PPT, na.rm = TRUE),
    n_valid_days = sum(!is.na(Q) & !is.na(PPT))
  ), by = water_year]

  # Calculate expected days per water year (accounting for leap years)
  annual[, expected_days := ifelse(
    ((water_year %% 4 == 0) & (water_year %% 100 != 0)) | (water_year %% 400 == 0),
    366, 365
  )]

  # Filter: remove years with >10% missing data (threshold from config.R)
  annual <- annual[n_valid_days >= ELASTICITY_MIN_DATA_COMPLETENESS * expected_days]

  # Remove years with zero or very low precipitation (threshold from config.R)
  annual <- annual[P_annual > ELASTICITY_MIN_ANNUAL_PPT]

  if (nrow(annual) < MIN_YEARS_ELASTICITY) {
    log_warn("Insufficient years (", nrow(annual), ") for elasticity calculation",
             context = ctx)
    return(list(
      elasticity_static = NA_real_,
      elasticity_senn_slp = NA_real_,
      elasticity_linear_slp = NA_real_,
      elasticity_spearman_rho = NA_real_,
      elasticity_spearman_pval = NA_real_,
      elasticity_mk_rho = NA_real_,
      elasticity_mk_pval = NA_real_,
      elasticity_mean = NA_real_,
      elasticity_median = NA_real_
    ))
  }

  # Calculate long-term means
  Q_mean <- mean(annual$Q_annual, na.rm = TRUE)
  P_mean <- mean(annual$P_annual, na.rm = TRUE)

  # Calculate annual elasticity values
  # E_i = (dQ_i/dP_i) / (Q_mean/P_mean)
  # where dQ_i = Q_i - Q_mean, dP_i = P_i - P_mean
  annual[, dQ := Q_annual - Q_mean]
  annual[, dP := P_annual - P_mean]

  # Avoid division by zero
  annual[, elasticity := ifelse(abs(dP) > 0.1,
                                 (dQ / dP) / (Q_mean / P_mean),
                                 NA_real_)]

  # Static elasticity is the median
  elasticity_static <- median(annual$elasticity, na.rm = TRUE)

  # Calculate rolling window elasticity if requested
  if (!is.null(rolling_window) && nrow(annual) >= rolling_window) {
    # Calculate elasticity for each rolling window
    setorder(annual, water_year)
    n <- nrow(annual)

    rolling_elasticity <- data.table(
      water_year = annual$water_year[(rolling_window):n],
      elasticity_rolling = sapply((rolling_window):n, function(end_idx) {
        start_idx <- end_idx - rolling_window + 1
        window_data <- annual[start_idx:end_idx]

        Q_mean_w <- mean(window_data$Q_annual, na.rm = TRUE)
        P_mean_w <- mean(window_data$P_annual, na.rm = TRUE)

        window_data[, dQ_w := Q_annual - Q_mean_w]
        window_data[, dP_w := P_annual - P_mean_w]
        window_data[, e_w := ifelse(abs(dP_w) > 0.1,
                                     (dQ_w / dP_w) / (Q_mean_w / P_mean_w),
                                     NA_real_)]
        median(window_data$e_w, na.rm = TRUE)
      })
    )

    # Calculate trend statistics on rolling elasticity
    trend_stats <- generate_stats(rolling_elasticity,
                                   value_cols = "elasticity_rolling",
                                   year_col = "water_year")

    result <- list(
      elasticity_static = elasticity_static,
      elasticity_senn_slp = trend_stats$elasticity_rolling_senn_slp,
      elasticity_linear_slp = trend_stats$elasticity_rolling_linear_slp,
      elasticity_spearman_rho = trend_stats$elasticity_rolling_spearman_rho,
      elasticity_spearman_pval = trend_stats$elasticity_rolling_spearman_pval,
      elasticity_mk_rho = trend_stats$elasticity_rolling_mk_rho,
      elasticity_mk_pval = trend_stats$elasticity_rolling_mk_pval,
      elasticity_mean = trend_stats$elasticity_rolling_mean,
      elasticity_median = trend_stats$elasticity_rolling_median
    )
  } else {
    # No rolling window - use annual values for stats
    trend_stats <- generate_stats(annual[!is.na(elasticity)],
                                   value_cols = "elasticity",
                                   year_col = "water_year")

    result <- list(
      elasticity_static = elasticity_static,
      elasticity_senn_slp = trend_stats$elasticity_senn_slp,
      elasticity_linear_slp = trend_stats$elasticity_linear_slp,
      elasticity_spearman_rho = trend_stats$elasticity_spearman_rho,
      elasticity_spearman_pval = trend_stats$elasticity_spearman_pval,
      elasticity_mk_rho = trend_stats$elasticity_mk_rho,
      elasticity_mk_pval = trend_stats$elasticity_mk_pval,
      elasticity_mean = trend_stats$elasticity_mean,
      elasticity_median = trend_stats$elasticity_median
    )
  }

  return(result)
}


#' Calculate Q-P seasonality metrics (Wrede et al. 2015)
#'
#' Quantifies seasonality in the cumulative Q vs cumulative P relationship.
#' Uses two metrics:
#' - qp_slope_sd: Standard deviation of monthly Q-P slopes
#' - qp_bimodality: Bimodality coefficient of slope distribution
#'
#' @param streamflow_data data.table with water_year, Q, prcp, month, dowy columns
#' @param slope_window_days Window size for rolling slope calculation
#' @return List with seasonality statistics
calculate_qp_seasonality <- function(streamflow_data,
                                      slope_window_days = QP_SLOPE_WINDOW_DAYS) {


  ctx <- "calculate_qp_seasonality"

  # Validate required columns
  required_cols <- c("water_year", "Q", "PPT", "month", "dowy")
  validate_columns(streamflow_data, required_cols, "streamflow_data", context = ctx)

  years <- unique(streamflow_data$water_year)

  if (length(years) < 10) {
    log_warn("Insufficient years for Q-P seasonality", context = ctx)
    return(list(
      qp_slope_sd_senn_slp = NA_real_, qp_slope_sd_linear_slp = NA_real_,
      qp_slope_sd_spearman_rho = NA_real_, qp_slope_sd_spearman_pval = NA_real_,
      qp_slope_sd_mk_rho = NA_real_, qp_slope_sd_mk_pval = NA_real_,
      qp_slope_sd_mean = NA_real_, qp_slope_sd_median = NA_real_,
      qp_bimodality_senn_slp = NA_real_, qp_bimodality_linear_slp = NA_real_,
      qp_bimodality_spearman_rho = NA_real_, qp_bimodality_spearman_pval = NA_real_,
      qp_bimodality_mk_rho = NA_real_, qp_bimodality_mk_pval = NA_real_,
      qp_bimodality_mean = NA_real_, qp_bimodality_median = NA_real_
    ))
  }

  # Calculate annual metrics
  annual_metrics <- rbindlist(lapply(years, function(yr) {
    year_data <- copy(streamflow_data[water_year == yr])

    if (nrow(year_data) < 300) return(NULL)

    # Check for too many NA values in Q or PPT
    na_frac_Q <- sum(is.na(year_data$Q)) / nrow(year_data)
    na_frac_PPT <- sum(is.na(year_data$PPT)) / nrow(year_data)
    if (na_frac_Q > 0.1 || na_frac_PPT > 0.1) return(NULL)

    # Replace remaining NA with 0 for cumsum (small gaps only)
    year_data[is.na(Q), Q := 0]
    year_data[is.na(PPT), PPT := 0]

    # Sort by day of water year
    setorder(year_data, dowy)

    # Calculate cumulative Q and P
    year_data[, cum_Q := cumsum(Q)]
    year_data[, cum_P := cumsum(PPT)]

    # Calculate rolling slope of cum_Q vs cum_P
    n <- nrow(year_data)
    if (n < slope_window_days) return(NULL)

    # Calculate slope for each window
    slopes <- sapply((slope_window_days):n, function(end_idx) {
      start_idx <- end_idx - slope_window_days + 1
      window <- year_data[start_idx:end_idx]

      # Simple linear regression slope
      delta_P <- window$cum_P[nrow(window)] - window$cum_P[1]
      delta_Q <- window$cum_Q[nrow(window)] - window$cum_Q[1]

      if (abs(delta_P) < 0.01) return(NA_real_)
      delta_Q / delta_P
    })

    # Associate slopes with months (use middle of window, not end)
    mid_offsets <- (slope_window_days):n - floor(slope_window_days / 2)
    slope_months <- year_data$month[mid_offsets]

    # Calculate monthly mean slopes
    slope_dt <- data.table(month = slope_months, slope = slopes)
    monthly_slopes <- slope_dt[!is.na(slope), .(mean_slope = mean(slope)), by = month]

    if (nrow(monthly_slopes) < 6) return(NULL)

    # Metric 1: SD of monthly slopes
    qp_slope_sd <- sd(monthly_slopes$mean_slope, na.rm = TRUE)

    # Metric 2: Bimodality coefficient
    # B = (skewness^2 + 1) / kurtosis
    # Values > 5/9 (0.555) suggest bimodality
    slopes_clean <- slopes[!is.na(slopes)]
    if (length(slopes_clean) < 30) {
      qp_bimodality <- NA_real_
    } else {
      n_s <- length(slopes_clean)
      m <- mean(slopes_clean)
      s <- sd(slopes_clean)

      if (s < 1e-10) {
        qp_bimodality <- NA_real_
      } else {
        skewness <- sum((slopes_clean - m)^3) / (n_s * s^3)
        kurtosis <- sum((slopes_clean - m)^4) / (n_s * s^4)

        if (kurtosis < 1e-10) {
          qp_bimodality <- NA_real_
        } else {
          qp_bimodality <- (skewness^2 + 1) / kurtosis
        }
      }
    }

    data.table(
      water_year = yr,
      qp_slope_sd = qp_slope_sd,
      qp_bimodality = qp_bimodality
    )
  }))

  if (is.null(annual_metrics) || nrow(annual_metrics) < 5) {
    log_warn("Could not calculate Q-P seasonality for enough years", context = ctx)
    return(list(
      qp_slope_sd_senn_slp = NA_real_, qp_slope_sd_linear_slp = NA_real_,
      qp_slope_sd_spearman_rho = NA_real_, qp_slope_sd_spearman_pval = NA_real_,
      qp_slope_sd_mk_rho = NA_real_, qp_slope_sd_mk_pval = NA_real_,
      qp_slope_sd_mean = NA_real_, qp_slope_sd_median = NA_real_,
      qp_bimodality_senn_slp = NA_real_, qp_bimodality_linear_slp = NA_real_,
      qp_bimodality_spearman_rho = NA_real_, qp_bimodality_spearman_pval = NA_real_,
      qp_bimodality_mk_rho = NA_real_, qp_bimodality_mk_pval = NA_real_,
      qp_bimodality_mean = NA_real_, qp_bimodality_median = NA_real_
    ))
  }

  # Calculate trend statistics
  sd_stats <- generate_stats(annual_metrics, value_cols = "qp_slope_sd",
                              year_col = "water_year")
  bi_stats <- generate_stats(annual_metrics[!is.na(qp_bimodality)],
                              value_cols = "qp_bimodality", year_col = "water_year")

  result <- list(
    qp_slope_sd_senn_slp = sd_stats$qp_slope_sd_senn_slp,
    qp_slope_sd_linear_slp = sd_stats$qp_slope_sd_linear_slp,
    qp_slope_sd_spearman_rho = sd_stats$qp_slope_sd_spearman_rho,
    qp_slope_sd_spearman_pval = sd_stats$qp_slope_sd_spearman_pval,
    qp_slope_sd_mk_rho = sd_stats$qp_slope_sd_mk_rho,
    qp_slope_sd_mk_pval = sd_stats$qp_slope_sd_mk_pval,
    qp_slope_sd_mean = sd_stats$qp_slope_sd_mean,
    qp_slope_sd_median = sd_stats$qp_slope_sd_median,
    qp_bimodality_senn_slp = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_senn_slp else NA_real_,
    qp_bimodality_linear_slp = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_linear_slp else NA_real_,
    qp_bimodality_spearman_rho = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_spearman_rho else NA_real_,
    qp_bimodality_spearman_pval = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_spearman_pval else NA_real_,
    qp_bimodality_mk_rho = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_mk_rho else NA_real_,
    qp_bimodality_mk_pval = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_mk_pval else NA_real_,
    qp_bimodality_mean = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_mean else NA_real_,
    qp_bimodality_median = if(nrow(bi_stats) > 0) bi_stats$qp_bimodality_median else NA_real_
  )

  return(result)
}


#' Calculate average storage metric (Peters & Aulenbach 2011)
#'
#' Uses simplified water balance: dS = P - Q (no ET estimation)
#' Calculates storage amount corresponding to average discharge.
#'
#' @param streamflow_data data.table with water_year, Q, prcp, dowy columns
#' @return List with storage statistics
calculate_average_storage <- function(streamflow_data) {

  ctx <- "calculate_average_storage"

  # Validate required columns
  required_cols <- c("water_year", "Q", "PPT", "dowy")
  validate_columns(streamflow_data, required_cols, "streamflow_data", context = ctx)

  years <- unique(streamflow_data$water_year)

  if (length(years) < 10) {
    log_warn("Insufficient years for storage calculation", context = ctx)
    return(list(
      avg_storage_senn_slp = NA_real_,
      avg_storage_linear_slp = NA_real_,
      avg_storage_spearman_rho = NA_real_,
      avg_storage_spearman_pval = NA_real_,
      avg_storage_mk_rho = NA_real_,
      avg_storage_mk_pval = NA_real_,
      avg_storage_mean = NA_real_,
      avg_storage_median = NA_real_
    ))
  }

  # Calculate storage for each year
  annual_storage <- rbindlist(lapply(years, function(yr) {
    year_data <- copy(streamflow_data[water_year == yr])

    if (nrow(year_data) < 300) return(NULL)

    # Check for too many NA values in Q or PPT
    na_frac_Q <- sum(is.na(year_data$Q)) / nrow(year_data)
    na_frac_PPT <- sum(is.na(year_data$PPT)) / nrow(year_data)
    if (na_frac_Q > 0.1 || na_frac_PPT > 0.1) return(NULL)

    # Replace remaining NA with 0 for cumsum (small gaps only)
    year_data[is.na(Q), Q := 0]
    year_data[is.na(PPT), PPT := 0]

    # Sort by day of water year
    setorder(year_data, dowy)

    # Water balance: dS = P - Q
    year_data[, dS := PPT - Q]

    # Cumulative storage (relative to start of year)
    year_data[, S := cumsum(dS)]

    # Calculate mean Q for this year
    mean_Q <- mean(year_data$Q, na.rm = TRUE)

    # Find storage-discharge relationship
    # Interpolate to find S at mean Q
    # Sort by Q for interpolation
    sq_data <- year_data[!is.na(Q) & !is.na(S)]
    if (nrow(sq_data) < 10) return(NULL)

    setorder(sq_data, Q)

    # Linear interpolation to find S at mean_Q
    avg_storage <- tryCatch({
      approx(sq_data$Q, sq_data$S, xout = mean_Q, rule = 2)$y
    }, error = function(e) NA_real_)

    data.table(
      water_year = yr,
      avg_storage = avg_storage
    )
  }))

  if (is.null(annual_storage) || nrow(annual_storage) < 5) {
    log_warn("Could not calculate storage for enough years", context = ctx)
    return(list(
      avg_storage_senn_slp = NA_real_,
      avg_storage_linear_slp = NA_real_,
      avg_storage_spearman_rho = NA_real_,
      avg_storage_spearman_pval = NA_real_,
      avg_storage_mk_rho = NA_real_,
      avg_storage_mk_pval = NA_real_,
      avg_storage_mean = NA_real_,
      avg_storage_median = NA_real_
    ))
  }

  # Calculate trend statistics
  storage_stats <- generate_stats(annual_storage[!is.na(avg_storage)],
                                   value_cols = "avg_storage",
                                   year_col = "water_year")

  result <- list(
    avg_storage_senn_slp = storage_stats$avg_storage_senn_slp,
    avg_storage_linear_slp = storage_stats$avg_storage_linear_slp,
    avg_storage_spearman_rho = storage_stats$avg_storage_spearman_rho,
    avg_storage_spearman_pval = storage_stats$avg_storage_spearman_pval,
    avg_storage_mk_rho = storage_stats$avg_storage_mk_rho,
    avg_storage_mk_pval = storage_stats$avg_storage_mk_pval,
    avg_storage_mean = storage_stats$avg_storage_mean,
    avg_storage_median = storage_stats$avg_storage_median
  )

  return(result)
}


################################################################################
# temp: convert caravan to unique csvs
process_caravan_to_annual <- function(caravan_directory, 
                                      data_project = "camels",
                                      min_num_years_data = 20,
                                      start_date_filter = as.Date("1900-01-01"), 
                                      end_date_filter = as.Date("2024-12-31"),
                                      output_dir = "annualized_caravan_data",
                                      min_num_days = 328
) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get list of NetCDF files
  pattern <- paste0(data_project, "_.*\\.nc$")
  nc_files <- list.files(caravan_directory, pattern = pattern, full.names = TRUE)
  
  cat("Found", length(nc_files), "NetCDF files for", data_project, "\n")
  
  # Process each watershed
  for (nc_file in nc_files) {
    # Extract watershed ID from filename
    watershed_id <- gsub(paste0(data_project, "_"), "", basename(nc_file))
    watershed_id <- gsub("\\.nc$", "", watershed_id)
    
    cat("Processing watershed", watershed_id, "...")
    
    tryCatch({
      # Read data from NetCDF
      nc_data <- nc_open(nc_file)
      
      # Extract all variables including Temperature
      streamflow_raw <- ncvar_get(nc_data, "streamflow")
      ppt_raw <- ncvar_get(nc_data, "total_precipitation_sum")
      swe_raw <- ncvar_get(nc_data, "snow_depth_water_equivalent_mean")
      temp_raw <- ncvar_get(nc_data, "temperature_2m_mean")
      
      # Extract date dimension
      time_raw <- ncvar_get(nc_data, "date")
      time_units <- ncatt_get(nc_data, "date", "units")$value
      nc_close(nc_data)
      
      # Convert time to Date objects
      origin_date <- as.Date(sub("days since ", "", time_units))
      dates <- as.Date(time_raw, origin = origin_date)
      
      # Create data.table
      stream_dt <- data.table(
        Date = dates,
        Q = as.numeric(streamflow_raw),
        PPT = as.numeric(ppt_raw),
        SWE = as.numeric(swe_raw),
        Temperature = as.numeric(temp_raw)
      )
      
      # Filter by date range
      stream_dt <- stream_dt[Date >= start_date_filter & Date <= end_date_filter]
      
      # Add water year information
      wy_info <- calculate_water_year_info(stream_dt$Date)
      stream_dt[, water_year := wy_info$water_year]
      
      # Calculate annual aggregations by water year with valid day counts
      annual_data <- stream_dt[, .(
        Q_annual_sum = sum(Q, na.rm = TRUE),
        Q_valid_days = sum(!is.na(Q)),
        PPT_annual_sum = sum(PPT, na.rm = TRUE),
        Temperature_annual_mean = mean(Temperature, na.rm = TRUE),
        SWE_annual_mean = mean(SWE, na.rm = TRUE),
        SWE_annual_max = max(SWE, na.rm = TRUE)
      ), by = water_year]
      
      # Order by water year
      setorder(annual_data, water_year)
      
      # Apply filtering rules for streamflow
      # Identify water years with at least 328 valid Q days
      annual_data[, sufficient_Q_data := Q_valid_days >= min_num_days]
      
      # Find first and last water years with sufficient data
      sufficient_years <- annual_data[sufficient_Q_data == TRUE, water_year]
      
      if (length(sufficient_years) > 0) {
        first_good_year <- min(sufficient_years)
        last_good_year <- max(sufficient_years)
        
        # Filter to keep only years from first to last good year
        annual_data <- annual_data[water_year >= first_good_year & water_year <= last_good_year]
        
        # For years within the range that have insufficient data, set Q to 0
        annual_data[sufficient_Q_data == FALSE, Q_annual_sum := 0]
      } else {
        # No years with sufficient data
        cat(" No years with sufficient data\n")
        next
      }
      
      # Remove helper columns
      annual_data[, c("Q_valid_days", "sufficient_Q_data") := NULL]
      
      # Save to CSV
      output_file <- file.path(output_dir, paste0(watershed_id, "_annualized.csv"))
      fwrite(annual_data, output_file)
      
      cat(" Done (", nrow(annual_data), " years)\n", sep = "")
      
    }, error = function(e) {
      cat(" Error:", e$message, "\n")
    })
  }
  
  cat("Processing complete. Output saved to:", output_dir, "\n")
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
  combined_metadata_path <- file.path(output_dir, "combined_watershed_metadata.csv")
  fwrite(all_metadata, combined_metadata_path)
  cat("Combined metadata saved with", nrow(all_metadata), "watersheds\n")

  # Enrich with human interference metadata (GAGES-II for USGS, HYDAT for Canadian)
  cat("\nEnriching metadata with human interference data...\n")
  tryCatch({
    all_metadata <- enrich_metadata_with_interference(
      metadata_file_path = combined_metadata_path,
      gages_ii_dir = GAGES_II_DIR
    )
    cat("Human interference metadata enrichment complete\n")
  }, error = function(e) {
    log_warn("Human interference enrichment failed:", e$message,
             context = "concatenate_with_metadata")
    cat("Warning: Human interference enrichment failed:", e$message, "\n")
    cat("Continuing without enrichment...\n")
  })

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




















########################################################################################################################
## older versions of functions that are no longer used



# old; no longer used
process_gages <- function(gages_df, gage_type, min_num_years, start_date, end_date, 
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
      # Fetch drainage area from tidyhydat database
      basin_area <- tryCatch({
        stn_info <- hy_stations(gage_id)
        as.numeric(stn_info$DRAINAGE_AREA_GROSS)
      }, error = function(e) NA_real_)
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
      
      # Try to add required columns if missing
      required_cols <- c("year", "Q", "doy")
      if (!all(required_cols %in% colnames(streamflow_data))) {
        if (!"year" %in% colnames(streamflow_data) && "Date" %in% colnames(streamflow_data)) {
          streamflow_data$year <- year(streamflow_data$Date)
        }
        if (!"doy" %in% colnames(streamflow_data) && "Date" %in% colnames(streamflow_data)) {
          streamflow_data$doy <- yday(streamflow_data$Date)
        }
        if (!all(required_cols %in% colnames(streamflow_data))) {
          cat("Missing required columns in streamflow data for gage", gage_id, "\n")
          next
        }
      }
      
      years_to_use <- NULL
      for (this_year in unique(streamflow_data$year)) {
        test_year <- subset(streamflow_data, year == this_year)
        nonzero_rows <- which(test_year$Q > min_Q_value_and_days[1])
        if (length(nonzero_rows) > min_Q_value_and_days[2]) {
          years_to_use <- c(years_to_use, this_year)
        }
      }
      
      if (length(years_to_use) <= min_num_years) {
        cat("Insufficient years with valid data for gage", gage_id, "\n")
        next
      }
      
      streamflow_data <- streamflow_data[streamflow_data$year %in% years_to_use, ]
      
      upstream_basins <- NULL
      tryCatch({
        upstream_basins <- find_upstream_hydrobasins(
          current_gage = current_gage,
          basinAt_NorAm_polys = basinAt_NorAm_polys,
          HB_dt = HB_dt,
          upstream_hydrobasins = upstream_hydrobasins,
          save_path = file.path(dirname(output_file), "upstream_hydrobasins.RData")
        )
      }, error = function(e) {
        cat("Error finding upstream basins for gage", gage_id, ":", e$message, "\n")
      })
      
      # [Calculate metrics...]
      # (Your existing calls to calculate_flow_vols_by_year, etc., remain unchanged.)
      
      # Create a row for this gage
      gage_row <- data.table(
        gage_id = gage_id,
        latitude = latitude,
        longitude = longitude,
        basin_area = basin_area,
        gage_type = gage_type,
        num_years = length(years_to_use),
        start_year = min(years_to_use),
        end_year = max(years_to_use)
      )
      
      # Append metrics (if available) to gage_row...
      if (!is.null(upstream_basins)) {
        gage_row$num_upstream_basins <- length(upstream_basins)
      } else {
        gage_row$num_upstream_basins <- NA
      }
      
      summary_output <- rbind(summary_output, gage_row, fill = TRUE)
      fwrite(summary_output, output_file)
      cat("Successfully processed gage", gage_id, "\n")
      
    }, error = function(e) {
      cat("Error processing gage", gage_id, ":", e$message, "\n")
    })
  }
  
  return(summary_output)
}

 # generate_streamflow_dt_og() — moved to archive/deprecated_generate_streamflow_dt_og.R
 # Superseded by the parquet-based pipeline in process_signatures_from_parquet()




####################################################################################################################
# Final user-facing function

# Function to process streamflow signatures from parquet files (long format)
process_signatures_from_parquet <- function(
    parquet_file_path = "combined_streamflow_output/combined_daily_data.parquet",
    metadata_file_path = "combined_streamflow_output/combined_watershed_metadata.csv",
    output_file = "processed_signatures_from_parquet.csv",
    min_Q_value_and_days = c(0.0001, 30),
    min_num_years = 20,
    min_frac_good_data = 0.95,
    daymet_parquet_path = NULL,  # Path to Daymet parquet for climate signatures
    start_water_year = NULL,     # Earliest water year to include (e.g., 1993)
    end_water_year = NULL        # Latest water year to include (e.g., 2022)
) {

  ctx <- "process_signatures_from_parquet"

  # Load required libraries
  require(arrow)
  require(data.table)
  require(lubridate)

  log_info("Starting signature extraction from parquet", context = ctx)

  # ==== INPUT VALIDATION ====
  log_info("Validating inputs...", context = ctx)
  validate_file_exists(parquet_file_path, "parquet_file_path",
                       required_ext = "parquet", context = ctx)
  validate_file_exists(metadata_file_path, "metadata_file_path",
                       required_ext = "csv", context = ctx)
  validate_numeric(min_Q_value_and_days, "min_Q_value_and_days",
                   min_val = 0, context = ctx)
  validate_numeric(min_num_years, "min_num_years",
                   min_val = 1, context = ctx)
  validate_numeric(min_frac_good_data, "min_frac_good_data",
                   min_val = 0, max_val = 1, context = ctx)

  # Validate water year range parameters
  if (!is.null(start_water_year)) {
    validate_numeric(start_water_year, "start_water_year", min_val = 1900, max_val = 2100, context = ctx)
    log_info("Water year range filter: start_water_year =", start_water_year, context = ctx)
  }
  if (!is.null(end_water_year)) {
    validate_numeric(end_water_year, "end_water_year", min_val = 1900, max_val = 2100, context = ctx)
    log_info("Water year range filter: end_water_year =", end_water_year, context = ctx)
  }
  if (!is.null(start_water_year) && !is.null(end_water_year)) {
    if (start_water_year > end_water_year) {
      log_error("start_water_year (", start_water_year, ") > end_water_year (", end_water_year, ")", context = ctx)
      stop("start_water_year must be <= end_water_year")
    }
    potential_years <- end_water_year - start_water_year + 1
    if (potential_years < min_num_years) {
      log_warn("Water year range (", potential_years, " years) is less than min_num_years (", min_num_years, ")", context = ctx)
    }
  }

  # Validate Daymet parquet if provided
  has_daymet <- FALSE
  daymet_dataset <- NULL  # lazy Arrow dataset for batch loading
  if (!is.null(daymet_parquet_path)) {
    if (file.exists(daymet_parquet_path)) {
      validate_file_exists(daymet_parquet_path, "daymet_parquet_path",
                           required_ext = "parquet", context = ctx)
      has_daymet <- TRUE
      daymet_dataset <- arrow::open_dataset(daymet_parquet_path, format = "parquet")
      log_info("Daymet data available - climate signatures will be calculated",
               context = ctx)
    } else {
      log_warn("Daymet parquet not found:", daymet_parquet_path,
               "- climate signatures will be skipped", context = ctx)
    }
  } else {
    log_info("No Daymet path provided - climate signatures will be skipped",
             context = ctx)
  }

  # Validate output directory exists (create if needed)
  output_dir <- dirname(output_file)
  if (output_dir != "." && output_dir != "") {
    validate_directory(output_dir, "output_directory", create = TRUE, context = ctx)
  }

  log_info("Reading parquet file and metadata...", context = ctx)

  # Open parquet lazily via Arrow to avoid loading entire dataset into memory
  # (111M+ rows can exceed available RAM on 16 GB systems)
  streamflow_dataset <- tryCatch({
    arrow::open_dataset(parquet_file_path, format = "parquet")
  }, error = function(e) {
    log_error("Failed to open parquet file:", e$message, context = ctx)
    stop(paste0("Failed to open parquet file: ", e$message))
  })

  # Check the structure of the data
  parquet_cols <- names(streamflow_dataset)
  log_info("Data columns:", paste(parquet_cols, collapse = ", "), context = ctx)

  # Validate streamflow data structure
  required_cols <- c("gage_id", "Date", "Q")
  missing_cols <- setdiff(required_cols, parquet_cols)
  if (length(missing_cols) > 0) {
    stop(paste0("streamflow_data missing required columns: ",
                paste(missing_cols, collapse = ", ")))
  }

  # Read metadata - ensure gage_id is character
  metadata <- tryCatch({
    fread(metadata_file_path, colClasses = list(character = "gage_id"))
  }, error = function(e) {
    log_error("Failed to read metadata file:", e$message, context = ctx)
    stop(paste0("Failed to read metadata file: ", e$message))
  })
  metadata[, gage_id := as.character(gage_id)]
  log_info("Number of gages in metadata:", nrow(metadata), context = ctx)

  # Initialize output data table
  summary_output <- data.table()

  # Get list of unique gage IDs from the parquet data (lightweight Arrow query)
  all_streamflow_gages <- streamflow_dataset |>
    dplyr::distinct(gage_id) |>
    dplyr::collect()
  all_streamflow_gages <- as.character(all_streamflow_gages$gage_id)
  log_info("Found", length(all_streamflow_gages), "unique gages in streamflow parquet", context = ctx)

  # If Daymet is available, iterate over Daymet site_ids (guarantees climate data)
  # Build a mapping: streamflow_gage_id -> daymet_site_id
  daymet_id_for_gage <- character(0)  # named vector: names=streamflow_id, values=daymet_id

  if (has_daymet && !is.null(daymet_dataset)) {
    daymet_site_ids <- daymet_dataset |>
      dplyr::distinct(site_id) |>
      dplyr::collect()
    daymet_site_ids <- as.character(daymet_site_ids$site_id)
    log_info("Found", length(daymet_site_ids), "unique sites in Daymet parquet", context = ctx)

    # For each Daymet site_id, find the matching streamflow gage_id
    matched_streamflow_ids <- character(0)
    unmatched_daymet_ids <- character(0)

    for (d_id in daymet_site_ids) {
      # Try exact match first
      if (d_id %in% all_streamflow_gages) {
        matched_streamflow_ids <- c(matched_streamflow_ids, d_id)
        daymet_id_for_gage[d_id] <- d_id
        next
      }
      # Try adding leading zeros (1-4)
      found <- FALSE
      for (num_zeros in 1:4) {
        padded_id <- paste0(strrep("0", num_zeros), d_id)
        if (padded_id %in% all_streamflow_gages) {
          matched_streamflow_ids <- c(matched_streamflow_ids, padded_id)
          daymet_id_for_gage[padded_id] <- d_id
          found <- TRUE
          break
        }
      }
      if (found) next
      # Try stripping leading zeros
      stripped_id <- gsub("^0+", "", d_id)
      if (stripped_id != d_id && stripped_id != "" && stripped_id %in% all_streamflow_gages) {
        matched_streamflow_ids <- c(matched_streamflow_ids, stripped_id)
        daymet_id_for_gage[stripped_id] <- d_id
        next
      }
      unmatched_daymet_ids <- c(unmatched_daymet_ids, d_id)
    }

    unique_gages <- matched_streamflow_ids
    log_info("Matched", length(unique_gages), "of", length(daymet_site_ids),
             "Daymet sites to streamflow gages", context = ctx)
    if (length(unmatched_daymet_ids) > 0) {
      log_warn(length(unmatched_daymet_ids), "Daymet sites had no matching streamflow gage",
               context = ctx)
      if (length(unmatched_daymet_ids) <= 20) {
        log_debug("Unmatched Daymet IDs:", paste(unmatched_daymet_ids, collapse = ", "),
                  context = ctx)
      }
    }
  } else {
    # No Daymet: fall back to streamflow-only iteration
    unique_gages <- all_streamflow_gages
    log_info("No Daymet data - iterating over", length(unique_gages), "streamflow gages",
             context = ctx)
  }
  
  # Create a keyed lookup table for metadata with multiple possible ID formats
  # Important for USGS data where IDs often have leading 0s that may be dropped.
  # Vectorized construction (replaces O(n²) row-by-row rbind loop).
  metadata_lookup <- copy(metadata)
  metadata_lookup[, gage_id := as.character(gage_id)]

  # Build alternate-ID rows for entries with leading zeros
  stripped_ids <- gsub("^0+", "", metadata_lookup$gage_id)
  has_leading_zeros <- stripped_ids != metadata_lookup$gage_id & stripped_ids != ""
  if (any(has_leading_zeros)) {
    alt_rows <- copy(metadata_lookup[has_leading_zeros])
    alt_rows[, gage_id := stripped_ids[has_leading_zeros]]
    metadata_lookup <- rbind(metadata_lookup, alt_rows, fill = TRUE)
  }

  # Key the table for O(1) lookup
  setkey(metadata_lookup, gage_id)

  # Function to find metadata for a gage ID, trying different formats
  # NOTE: Parameter named target_gage_id to avoid data.table scoping collision
  # with the 'gage_id' column in metadata_lookup
  find_metadata <- function(target_gage_id, metadata_lookup) {
    target_gage_id <- as.character(target_gage_id)

    # O(1) keyed lookup — exact match
    meta <- metadata_lookup[.(target_gage_id), nomatch = NULL]
    if (nrow(meta) > 0) return(meta[1])

    # Try adding leading zeros (up to 4)
    for (num_zeros in 1:4) {
      padded_id <- paste0(strrep("0", num_zeros), target_gage_id)
      meta <- metadata_lookup[.(padded_id), nomatch = NULL]
      if (nrow(meta) > 0) {
        log_debug("Found match for", target_gage_id, "as", padded_id, context = ctx)
        return(meta[1])
      }
    }

    # Try removing leading zeros
    stripped_id <- gsub("^0+", "", target_gage_id)
    if (stripped_id != target_gage_id && stripped_id != "") {
      meta <- metadata_lookup[.(stripped_id), nomatch = NULL]
      if (nrow(meta) > 0) {
        log_debug("Found match for", target_gage_id, "as", stripped_id, context = ctx)
        return(meta[1])
      }
    }

    return(NULL)
  }
  
  # Track matching statistics
  matched_gages <- 0
  unmatched_gages <- character()

  # Per-metric success counters (for R5 coverage tracking)
  all_metric_names <- c("flow_vols", "fdc_trends", "flashiness", "flow_timing",
                        "pulses", "baseflow", "recession",
                        "qtoppt", "elasticity", "qp_seasonality", "avg_storage")
  metric_success <- setNames(integer(length(all_metric_names)), all_metric_names)

  # Batch loading: load groups of gages from parquet at a time to balance
  # memory usage (~1-2 GB per batch) vs speed (fast data.table lookups within batch)
  BATCH_SIZE <- PROCESSING_BATCH_SIZE  # From config.R
  n_batches <- ceiling(length(unique_gages) / BATCH_SIZE)
  log_info("Processing in", n_batches, "batches of up to", BATCH_SIZE, "gages", context = ctx)
  streamflow_batch <- NULL  # current batch data
  daymet_batch <- NULL       # current batch of Daymet climate data

  # Process each unique gage
  for (i in seq_along(unique_gages)) {
    # Use current_gage_id to avoid shadowing the data.table column name 'gage_id'
    current_gage_id <- unique_gages[i]

    if (i %% 100 == 0) {
      log_info("Processing gage", i, "of", length(unique_gages), context = ctx)
    }

    # Load a new batch of streamflow + Daymet data when needed
    if ((i - 1) %% BATCH_SIZE == 0) {
      batch_start <- i
      batch_end <- min(i + BATCH_SIZE - 1, length(unique_gages))
      batch_ids <- unique_gages[batch_start:batch_end]
      batch_num <- ceiling(i / BATCH_SIZE)
      log_info("Loading streamflow batch", batch_num, "of", n_batches,
               "(gages", batch_start, "-", batch_end, ")", context = ctx)

      # Free previous batches before loading new ones
      if (!is.null(streamflow_batch)) {
        rm(streamflow_batch)
        if (!is.null(daymet_batch)) rm(daymet_batch)
      }

      streamflow_batch <- streamflow_dataset |>
        dplyr::filter(gage_id %in% batch_ids) |>
        dplyr::collect() |>
        as.data.table()
      setkey(streamflow_batch, gage_id)

      # Also batch-load Daymet data for these gages (translate to Daymet IDs)
      if (has_daymet && !is.null(daymet_dataset)) {
        batch_daymet_ids <- unname(daymet_id_for_gage[batch_ids])
        batch_daymet_ids <- batch_daymet_ids[!is.na(batch_daymet_ids)]
        daymet_batch <- tryCatch({
          daymet_dataset |>
            dplyr::filter(site_id %in% batch_daymet_ids) |>
            dplyr::collect() |>
            as.data.table()
        }, error = function(e) {
          log_warn("Failed to batch-load Daymet data:", e$message, context = ctx)
          NULL
        })
        if (!is.null(daymet_batch) && nrow(daymet_batch) > 0) {
          setkey(daymet_batch, site_id)
          log_info("Loaded Daymet data for batch:", nrow(daymet_batch), "rows", context = ctx)
        } else {
          daymet_batch <- NULL
        }
      }
    }

    tryCatch({
      # Find metadata for this gage
      gage_meta <- find_metadata(current_gage_id, metadata_lookup)

      if (is.null(gage_meta)) {
        unmatched_gages <- c(unmatched_gages, current_gage_id)
        next
      }

      # If basin_area is NA and this looks like a Canadian gage (format: ##XX###),
      # try to fetch from tidyhydat
      if (is.na(gage_meta$basin_area) && grepl("^[0-9]{2}[A-Z]{2}", current_gage_id)) {
        gage_meta$basin_area <- tryCatch({
          stn_info <- hy_stations(current_gage_id)
          as.numeric(stn_info$DRAINAGE_AREA_GROSS)
        }, error = function(e) NA_real_)
      }

      matched_gages <- matched_gages + 1

      # Extract streamflow data for this gage from the current batch
      gage_flow <- streamflow_batch[gage_id == current_gage_id]

      if (nrow(gage_flow) == 0) {
        log_debug("No data for gage", current_gage_id, context = ctx)
        next
      }

      # Ensure required columns exist
      if (!"Date" %in% names(gage_flow)) {
        log_warn("Missing Date column for gage", current_gage_id, context = ctx)
        next
      }

      if (!"Q" %in% names(gage_flow)) {
        log_warn("Missing Q column for gage", current_gage_id, context = ctx)
        next
      }

      # Remove NA values in Q
      gage_flow <- gage_flow[!is.na(Q)]

      if (nrow(gage_flow) == 0) {
        log_debug("No valid Q data for gage", current_gage_id, context = ctx)
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

      # ============= INTEGRATE CLIMATE DATA =============
      gage_has_climate <- FALSE
      if (has_daymet && !is.null(daymet_batch)) {
        tryCatch({
          # Use pre-loaded batch data instead of per-gage parquet scan
          # Look up the Daymet site_id that corresponds to this streamflow gage_id
          daymet_id <- daymet_id_for_gage[current_gage_id]
          gage_climate <- if (!is.na(daymet_id)) daymet_batch[site_id == daymet_id] else data.table()
          if (nrow(gage_climate) > 0) {
            gage_climate <- copy(gage_climate)  # avoid modifying batch by reference
            gage_climate[, site_id := NULL]  # drop site_id before merge
            gage_flow <- merge(gage_flow, gage_climate, by = "Date", all.x = TRUE)
            # Rename prcp to PPT for compatibility
            if ("prcp" %in% names(gage_flow)) {
              setnames(gage_flow, "prcp", "PPT")
            }
            gage_has_climate <- "PPT" %in% names(gage_flow) &&
                                sum(!is.na(gage_flow$PPT)) > 0
            if (gage_has_climate) {
              log_debug("Climate data integrated for gage", current_gage_id, context = ctx)
            }
          } else {
            log_warn("No Daymet data found for gage:", current_gage_id,
                     context = "integrate_daymet_batch")
          }
        }, error = function(e) {
          log_warn("Could not integrate Daymet for gage", current_gage_id, ":",
                   e$message, context = ctx)
        })
      }

      # ============= APPLY WATER YEAR BASED FILTERS =============

      # Filter 1: Restrict to specified water year range (if provided)
      if (!is.null(start_water_year)) {
        n_before <- nrow(gage_flow)
        gage_flow <- gage_flow[water_year >= start_water_year]
        log_debug("Filtered to water_year >=", start_water_year, ":",
                  n_before, "->", nrow(gage_flow), "rows",
                  context = paste0(ctx, ":", current_gage_id))
      }
      if (!is.null(end_water_year)) {
        n_before <- nrow(gage_flow)
        gage_flow <- gage_flow[water_year <= end_water_year]
        log_debug("Filtered to water_year <=", end_water_year, ":",
                  n_before, "->", nrow(gage_flow), "rows",
                  context = paste0(ctx, ":", current_gage_id))
      }

      # Skip gage if no data remains after water year filtering
      if (nrow(gage_flow) == 0) {
        log_debug("Gage", current_gage_id, "has no data in specified water year range", context = ctx)
        next
      }

      # Filter 2: Check each water year for minimum days above threshold
      # Filter 2b: Check each water year for minimum fraction of non-NA data

      water_years_to_use <- NULL
      
      for (this_wy in unique(gage_flow$water_year)) {
        test_wy <- gage_flow[water_year == this_wy]
        
        # Filter 2: Check days with Q > threshold
        nonzero_rows <- sum(test_wy$Q > min_Q_value_and_days[1], na.rm = TRUE)
        if (nonzero_rows < min_Q_value_and_days[2]) {
          log_debug("WY", this_wy, "rejected: only", nonzero_rows,
                    "days > threshold (need", min_Q_value_and_days[2], ")",
                    context = paste0(ctx, ":", current_gage_id))
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
          log_debug("WY", this_wy, "rejected: only", n_good_days,
                    "valid days (need", min_required_days, "for",
                    sprintf("%.1f%%", min_frac_good_data * 100), "coverage)",
                    context = paste0(ctx, ":", current_gage_id))
          next
        }

        # This water year passes both filters
        water_years_to_use <- c(water_years_to_use, this_wy)
      }

      # Filter 3: Check minimum number of qualifying water years
      if (length(water_years_to_use) < min_num_years) {
        log_debug("Insufficient qualifying water years for gage", current_gage_id,
                  "(", length(water_years_to_use), "water years)", context = ctx)
        next
      }
      
      # Filter to qualifying water years
      streamflow_data_filtered <- gage_flow[water_year %in% water_years_to_use]
      
      # Initialize metric results
      metrics_list <- list()

      # Non-climate metric specifications: name, function, label
      non_climate_specs <- list(
        list(name = "flow_vols",  fn = calculate_flow_vols_by_year,            label = "Flow volumes"),
        list(name = "fdc_trends", fn = analyze_fdc_trends_from_streamflow,     label = "FDC trends"),
        list(name = "flashiness", fn = analyze_flashiness_trends,              label = "Flashiness"),
        list(name = "flow_timing",fn = analyze_flow_timing_trends,             label = "Flow timing"),
        list(name = "pulses",     fn = calculate_pulse_metrics,                label = "Pulse metrics"),
        list(name = "baseflow",   fn = analyze_baseflow_indices,               label = "Baseflow indices"),
        list(name = "recession",  fn = analyze_recession_parameters,           label = "Recession parameters")
      )

      for (spec in non_climate_specs) {
        result <- safe_calculate(spec$fn, streamflow_data_filtered, current_gage_id, spec$label, ctx, "warn")
        if (!is.null(result)) {
          metrics_list[[spec$name]] <- result
          metric_success[[spec$name]] <- metric_success[[spec$name]] + 1L
        }
      }

      # ============= CLIMATE-DEPENDENT SIGNATURES =============
      if (gage_has_climate) {
        climate_specs <- list(
          list(name = "qtoppt",        fn = analyze_Q_PPT_relationships,   label = "Q-PPT analysis"),
          list(name = "elasticity",    fn = calculate_streamflow_elasticity,label = "Elasticity calc"),
          list(name = "qp_seasonality",fn = calculate_qp_seasonality,      label = "Q-P seasonality calc"),
          list(name = "avg_storage",   fn = calculate_average_storage,     label = "Storage calc")
        )

        for (spec in climate_specs) {
          result <- safe_calculate(spec$fn, streamflow_data_filtered, current_gage_id, spec$label, ctx, "debug")
          if (!is.null(result)) {
            metrics_list[[spec$name]] <- result
            metric_success[[spec$name]] <- metric_success[[spec$name]] + 1L
          }
        }
      }

      # Create output row
      gage_row <- data.table(
        gage_id = current_gage_id,
        gage_id_metadata = gage_meta$gage_id,  # Store the metadata version too
        latitude = gage_meta$latitude,
        longitude = gage_meta$longitude,
        basin_area = gage_meta$basin_area,
        gage_type = gage_meta$gage_type,
        num_water_years = length(water_years_to_use),  # Changed from num_years
        start_water_year = min(water_years_to_use),     # Changed from start_year
        end_water_year = max(water_years_to_use),       # Changed from end_year
        analysis_start_water_year = start_water_year,    # Requested analysis range start (or NA if unrestricted)
        analysis_end_water_year = end_water_year,        # Requested analysis range end (or NA if unrestricted)
        num_upstream_basins = gage_meta$num_upstream_basins,
        area_normalized = gage_meta$area_normalized,
        # Human interference metadata (from GAGES-II / HYDAT)
        NDAMS_2009 = if ("NDAMS_2009" %in% names(gage_meta)) gage_meta$NDAMS_2009 else NA_real_,
        MAJ_DDENS_2009 = if ("MAJ_DDENS_2009" %in% names(gage_meta)) gage_meta$MAJ_DDENS_2009 else NA_real_,
        STOR_NID_2009 = if ("STOR_NID_2009" %in% names(gage_meta)) gage_meta$STOR_NID_2009 else NA_real_,
        IMPNLCD06 = if ("IMPNLCD06" %in% names(gage_meta)) gage_meta$IMPNLCD06 else NA_real_,
        DEVNLCD06 = if ("DEVNLCD06" %in% names(gage_meta)) gage_meta$DEVNLCD06 else NA_real_,
        FRESHW_WITHDRAWAL = if ("FRESHW_WITHDRAWAL" %in% names(gage_meta)) gage_meta$FRESHW_WITHDRAWAL else NA_real_,
        HYDRO_DISTURB_INDX = if ("HYDRO_DISTURB_INDX" %in% names(gage_meta)) gage_meta$HYDRO_DISTURB_INDX else NA_real_,
        CLASS = if ("CLASS" %in% names(gage_meta)) gage_meta$CLASS else NA_character_,
        RHBN = if ("RHBN" %in% names(gage_meta)) gage_meta$RHBN else NA,
        REGULATED = if ("REGULATED" %in% names(gage_meta)) gage_meta$REGULATED else NA,
        human_interference_class = if ("human_interference_class" %in% names(gage_meta)) gage_meta$human_interference_class else NA_character_
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
      log_error("Error processing gage", current_gage_id, ":", e$message, context = ctx)
    })

    # Periodic save
    if (i %% 500 == 0 && nrow(summary_output) > 0) {
      fwrite(summary_output, output_file)
      log_info("Saved intermediate results (", nrow(summary_output), "gages processed)",
               context = ctx)
    }
  }

  # Final save and validation
  if (nrow(summary_output) > 0) {
    fwrite(summary_output, output_file)

    log_info("========== PROCESSING COMPLETE ==========", context = ctx)
    log_info("Total unique gages in parquet:", length(unique_gages), context = ctx)
    log_info("Successfully matched to metadata:", matched_gages, context = ctx)
    log_info("Could not match to metadata:", length(unmatched_gages), context = ctx)
    log_info("Successfully processed:", nrow(summary_output), context = ctx)
    log_info("Results saved to:", output_file, context = ctx)

    # Per-metric success summary
    log_info("--- Per-metric success counts ---", context = ctx)
    for (mn in names(metric_success)) {
      log_info(sprintf("  %-20s %d / %d gages", mn, metric_success[[mn]], nrow(summary_output)),
               context = ctx)
    }

    # Validate output schema
    validation_result <- validate_output_schema(summary_output, strict = FALSE,
                                                 context = ctx)

    # Save list of unmatched gages for debugging
    if (length(unmatched_gages) > 0) {
      unmatched_file <- gsub("\\.csv$", "_unmatched_gages.txt", output_file)
      writeLines(unmatched_gages, unmatched_file)
      log_info("List of unmatched gages saved to:", unmatched_file, context = ctx)
      log_debug("Sample of unmatched gages:",
                paste(head(unmatched_gages, 10), collapse = ", "), context = ctx)
    }
  } else {
    log_warn("========== NO DATA PROCESSED ==========", context = ctx)
    log_warn("No gages were successfully processed.", context = ctx)
  }

  log_info("Signature extraction complete", context = ctx)
  return(summary_output)
}


# ==============================================================================
# HUMAN INTERFERENCE METADATA FUNCTIONS
# ==============================================================================

#' Load GAGES-II files for a specific region (CONUS or AKHIPR)
#'
#' Internal helper function to read and merge GAGES-II text files.
#'
#' @param gages_dir Path to directory containing GAGES-II files
#' @param file_list Named list of file names to load
#' @param region_name Region identifier for logging
#' @return data.table with STAID and extracted columns
load_gages_ii_region <- function(gages_dir, file_list, region_name) {
  ctx <- paste0("load_gages_ii_region:", region_name)

  result <- NULL

  for (file_type in names(file_list)) {
    file_path <- file.path(gages_dir, file_list[[file_type]])

    if (!file.exists(file_path)) {
      log_debug("GAGES-II file not found:", file_path, context = ctx)
      next
    }

    cols_to_keep <- GAGES_II_COLUMNS[[file_type]]
    if (is.null(cols_to_keep)) {
      log_warn("No column specification for file type:", file_type, context = ctx)
      next
    }

    tryCatch({
      dt <- fread(file_path, colClasses = c("STAID" = "character"))

      # Keep only specified columns that exist
      available_cols <- intersect(cols_to_keep, names(dt))
      if (!"STAID" %in% available_cols) {
        log_warn("STAID column missing in:", file_path, context = ctx)
        next
      }

      dt <- dt[, ..available_cols]

      if (is.null(result)) {
        result <- dt
      } else {
        result <- merge(result, dt, by = "STAID", all = TRUE)
      }

      log_debug("Loaded", nrow(dt), "rows from", basename(file_path), context = ctx)

    }, error = function(e) {
      log_warn("Failed to load", file_path, ":", e$message, context = ctx)
    })
  }

  if (is.null(result)) {
    result <- data.table(STAID = character())
  }

  return(result)
}


#' Load GAGES-II Human Interference Metadata
#'
#' Loads and combines human interference columns from GAGES-II text files
#' for CONUS and Alaska/Hawaii/PR gages.
#'
#' @param gages_dir Path to directory containing GAGES-II text files.
#'   Defaults to GAGES_II_DIR from config.R
#' @return data.table with STAID and interference columns:
#'   NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009, IMPNLCD06, DEVNLCD06,
#'   FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX, CLASS
load_gages_ii_interference <- function(gages_dir = GAGES_II_DIR) {
  ctx <- "load_gages_ii_interference"
  log_info("Loading GAGES-II human interference metadata from:", gages_dir, context = ctx)

  validate_directory(gages_dir, "gages_dir", context = ctx)

  # Load CONUS files
  conus_data <- load_gages_ii_region(gages_dir, GAGES_II_FILES_CONUS, "CONUS")
  log_info("Loaded CONUS GAGES-II data:", nrow(conus_data), "gages", context = ctx)

  # Load AKHIPR files
  akhipr_data <- load_gages_ii_region(gages_dir, GAGES_II_FILES_AKHIPR, "AKHIPR")
  log_info("Loaded AKHIPR GAGES-II data:", nrow(akhipr_data), "gages", context = ctx)

  # Combine (AKHIPR may not have all columns)
  combined <- rbind(conus_data, akhipr_data, fill = TRUE)

  # Replace sentinel missing value (-999) with NA
  numeric_cols <- names(combined)[sapply(combined, is.numeric)]
  for (col in numeric_cols) {
    combined[get(col) == GAGES_II_MISSING_VALUE, (col) := NA]
  }

  # Standardize CLASS values
  if ("CLASS" %in% names(combined)) {
    combined[, CLASS := trimws(CLASS)]
    combined[CLASS == "", CLASS := NA_character_]
  }

  log_info("Loaded GAGES-II interference data for", nrow(combined), "total gages", context = ctx)

  return(combined)
}


#' Load Canadian HYDAT Human Interference Metadata
#'
#' Extracts RHBN designation and regulation status from tidyhydat database.
#'
#' @param station_numbers Character vector of Canadian station numbers (optional).
#'   If NULL, loads data for all stations.
#' @return data.table with gage_id and interference columns:
#'   RHBN (Reference Hydrometric Basin Network flag),
#'   REGULATED (station regulation status),
#'   human_interference_class (unified classification)
load_canadian_interference <- function(station_numbers = NULL) {
  ctx <- "load_canadian_interference"

  if (!requireNamespace("tidyhydat", quietly = TRUE)) {
    log_warn("tidyhydat package not available - skipping Canadian interference metadata",
             context = ctx)
    return(data.table(
      gage_id = character(),
      RHBN = logical(),
      REGULATED = logical(),
      human_interference_class = character()
    ))
  }

  log_info("Loading Canadian HYDAT interference metadata", context = ctx)

  # Get station information (includes RHBN flag)
  stations <- tryCatch({
    if (is.null(station_numbers)) {
      as.data.table(tidyhydat::hy_stations())
    } else {
      as.data.table(tidyhydat::hy_stations(station_numbers))
    }
  }, error = function(e) {
    log_error("Failed to load hy_stations:", e$message, context = ctx)
    return(NULL)
  })

  if (is.null(stations) || nrow(stations) == 0) {
    return(data.table(
      gage_id = character(),
      RHBN = logical(),
      REGULATED = logical(),
      human_interference_class = character()
    ))
  }

  # Get regulation information
  regulation <- tryCatch({
    if (is.null(station_numbers)) {
      as.data.table(tidyhydat::hy_stn_regulation())
    } else {
      as.data.table(tidyhydat::hy_stn_regulation(station_numbers))
    }
  }, error = function(e) {
    log_warn("Failed to load hy_stn_regulation:", e$message, context = ctx)
    data.table(STATION_NUMBER = character(), REGULATED = logical())
  })

  # Select relevant columns from stations
  station_cols <- c("STATION_NUMBER", "RHBN")
  available_cols <- intersect(station_cols, names(stations))
  stations <- stations[, ..available_cols]

  # Merge with regulation status
  if ("STATION_NUMBER" %in% names(regulation) && "REGULATED" %in% names(regulation)) {
    # Keep most recent regulation status per station (any TRUE = regulated)
    regulation <- regulation[, .(REGULATED = any(REGULATED, na.rm = TRUE)),
                             by = STATION_NUMBER]
    stations <- merge(stations, regulation, by = "STATION_NUMBER", all.x = TRUE)
  } else {
    stations[, REGULATED := NA]
  }

  # Ensure RHBN column exists
  if (!"RHBN" %in% names(stations)) {
    stations[, RHBN := NA]
  }

  # Create human interference classification based on RHBN only
  # RHBN = TRUE means reference (minimal human impact)
  # RHBN = FALSE means non-reference
  # RHBN = NA means unknown
  stations[, human_interference_class := fcase(
    RHBN == TRUE, "reference",
    RHBN == FALSE, "non-reference",
    default = "unknown"
  )]

  # Rename for consistency with pipeline
  setnames(stations, "STATION_NUMBER", "gage_id")

  log_info("Loaded Canadian interference data for", nrow(stations), "stations", context = ctx)

  return(stations)
}


#' Enrich Combined Watershed Metadata with Human Interference Data
#'
#' Adds human interference columns from GAGES-II (USGS) and Canadian HYDAT
#' to existing combined_watershed_metadata.csv file.
#'
#' @param metadata_file_path Path to combined_watershed_metadata.csv
#' @param gages_ii_dir Path to GAGES-II data directory.
#'   Defaults to GAGES_II_DIR from config.R
#' @param output_file_path Path for enriched output. If NULL, overwrites input file.
#' @return data.table with enriched metadata including interference columns
enrich_metadata_with_interference <- function(
    metadata_file_path,
    gages_ii_dir = GAGES_II_DIR,
    output_file_path = NULL
) {
  ctx <- "enrich_metadata_with_interference"

  log_info("========== ENRICHING METADATA WITH HUMAN INTERFERENCE DATA ==========",
           context = ctx)

  # Load existing metadata
  validate_file_exists(metadata_file_path, "metadata_file_path",
                       required_ext = "csv", context = ctx)
  metadata <- fread(metadata_file_path, colClasses = c("gage_id" = "character"))
  log_info("Loaded metadata with", nrow(metadata), "gages", context = ctx)

  # Initialize interference columns with NA
  interference_cols <- c(
    "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
    "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
    "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
    "human_interference_class"
  )
  for (col in interference_cols) {
    if (!col %in% names(metadata)) {
      if (col %in% c("RHBN", "REGULATED")) {
        metadata[, (col) := NA]
      } else if (col %in% c("CLASS", "human_interference_class")) {
        metadata[, (col) := NA_character_]
      } else {
        metadata[, (col) := NA_real_]
      }
    }
  }

  # Load GAGES-II data
  gages_ii_data <- load_gages_ii_interference(gages_ii_dir)

  # Create lookup with both original and stripped gage IDs
  if (nrow(gages_ii_data) > 0) {
    gages_ii_data[, gage_id_stripped := gsub("^0+", "", STAID)]
  }

  # Load Canadian data
  canadian_gages <- metadata[gage_type %in% c("Canada", "CANADIAN"), gage_id]
  canadian_data <- load_canadian_interference(canadian_gages)

  # Join GAGES-II data to USGS gages
  usgs_idx <- which(metadata$gage_type == "USGS")
  if (length(usgs_idx) > 0 && nrow(gages_ii_data) > 0) {
    log_info("Joining GAGES-II data to", length(usgs_idx), "USGS gages", context = ctx)

    usgs_matched <- 0
    for (i in usgs_idx) {
      gage_id_orig <- metadata$gage_id[i]
      gage_id_stripped <- gsub("^0+", "", gage_id_orig)

      # Try exact match first, then stripped
      match_row <- gages_ii_data[STAID == gage_id_orig]
      if (nrow(match_row) == 0) {
        match_row <- gages_ii_data[gage_id_stripped == gage_id_stripped]
      }

      if (nrow(match_row) > 0) {
        match_row <- match_row[1]  # Take first if multiple
        usgs_matched <- usgs_matched + 1

        # Copy GAGES-II columns
        for (col in c("NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                      "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                      "HYDRO_DISTURB_INDX", "CLASS")) {
          if (col %in% names(match_row)) {
            set(metadata, i, col, match_row[[col]])
          }
        }

        # Compute human_interference_class for USGS based on CLASS column
        # CLASS = "Ref" -> reference, CLASS = "Non-ref" -> non-reference, else -> unknown
        class_val <- match_row$CLASS
        if (!is.na(class_val) && nchar(class_val) > 0) {
          interference_class <- ifelse(class_val == "Ref", "reference", "non-reference")
        } else {
          interference_class <- "unknown"
        }
        set(metadata, i, "human_interference_class", interference_class)
      }
    }

    log_info("GAGES-II match rate:", usgs_matched, "/", length(usgs_idx),
             sprintf("(%.1f%%)", 100 * usgs_matched / length(usgs_idx)), context = ctx)
  }

  # Join Canadian data
  canadian_idx <- which(metadata$gage_type %in% c("Canada", "CANADIAN"))
  if (length(canadian_idx) > 0 && nrow(canadian_data) > 0) {
    log_info("Joining HYDAT data to", length(canadian_idx), "Canadian gages", context = ctx)

    canadian_matched <- 0
    for (i in canadian_idx) {
      gage_id_orig <- metadata$gage_id[i]
      match_row <- canadian_data[gage_id == gage_id_orig]

      if (nrow(match_row) > 0) {
        match_row <- match_row[1]
        canadian_matched <- canadian_matched + 1

        for (col in c("RHBN", "REGULATED", "human_interference_class")) {
          if (col %in% names(match_row)) {
            set(metadata, i, col, match_row[[col]])
          }
        }
      }
    }

    log_info("Canadian HYDAT match rate:", canadian_matched, "/", length(canadian_idx),
             sprintf("(%.1f%%)", 100 * canadian_matched / length(canadian_idx)), context = ctx)
  }

  # Save enriched metadata
  if (is.null(output_file_path)) {
    output_file_path <- metadata_file_path
  }
  fwrite(metadata, output_file_path)
  log_info("Saved enriched metadata to:", output_file_path, context = ctx)

  # Summary statistics
  log_info("========== ENRICHMENT SUMMARY ==========", context = ctx)
  log_info("Total gages:", nrow(metadata), context = ctx)
  if ("human_interference_class" %in% names(metadata)) {
    class_counts <- table(metadata$human_interference_class, useNA = "ifany")
    for (class_name in names(class_counts)) {
      log_info("  ", ifelse(is.na(class_name), "NA", class_name), ":",
               class_counts[class_name], context = ctx)
    }
  }

  return(metadata)
}
