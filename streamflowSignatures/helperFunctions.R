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



generate_streamflow_dt <- function(dt, data_origin, 
                                   min_num_years = 20, 
                                   start_date = as.Date("1900-01-01"), 
                                   end_date = as.Date("2024-12-31")) {
  # Check that data_origin is valid; if not, warn and return NA.
  if (!data_origin %in% c("USGS", "Canada")) {
    warning("Invalid data_origin provided. It must be either 'USGS' or 'Canada'. Returning NA.")
    return(NULL)
  }
  
  # Ensure dt is a data.table for consistency.
  if (!inherits(dt, "data.table")) {
    dt <- as.data.table(dt)
  }
  
  output <- NA  # Default output
  
  if (data_origin == "USGS") {
    gage_data <- readNWISdv(siteNumber = dt$STAID,
                            parameterCd = "00060",
                            startDate = "1900-01-01", endDate = as.character(end_date))
    gage_data <- subset(gage_data, Date > start_date)
    gage_id <- gage_data$site_no[1]
    
    if (nrow(gage_data) > 365 * min_num_years & last(gage_data$Date) > start_date) {
      names(gage_data)[4] <- "Q_rawUnits"
      # Remove flagged data (assumes column 5 holds flags)
      gage_data$Q_rawUnits[-which(gage_data[, 5] %in% c("A", "A e", "P", "P e"))] <- NA
      
      streamy <- gage_data[, c("Date", "Q_rawUnits")]
      
      # Convert to mm/day using drainage area from dt
      sqkm <- dt$DRAIN_SQKM[dt$STAID == gage_id]
      conversion <- 60 * 60 * 24 / (sqkm * 3280.84^3) * 1e6
      streamy$Q <- as.numeric(streamy$Q_rawUnits) * conversion
      streamy$year = year(streamy$Date)
      streamy$month = month(streamy$Date)
      streamy$doy = yday(streamy$Date)
      
      output <- streamy
    } else {
      message("Insufficient Data to Process")
      output <- NA
    }
  }
  
  if (data_origin == "Canada") {
    can_stream <- hy_daily(station_number = paste(dt$STATION_NUMBER))
    can_stream_only <- subset(can_stream, Parameter == "Flow")
    
    if ("Flow" %in% can_stream$Parameter & last(can_stream_only$Date) > start_date & 
        nrow(can_stream_only) > 365 * min_num_years) {
      stream_all <- cbind.data.frame(as.Date(can_stream_only$Date), can_stream_only$Value)
      colnames(stream_all) <- c("Date", "Q_rawUnits")
      streamy <- subset(stream_all, Date > start_date)
      
      # Converting to mm/day for m^3/s
      sqkm <- hy_stations(paste(dt$STATION_NUMBER))$DRAINAGE_AREA_GROSS
      conversion <- ifelse(is.na(sqkm), 99999, 60 * 60 * 24 * 1e9 / (sqkm * 1e12))
      streamy$Q <- as.numeric(streamy$Q_rawUnits) * conversion
      
      output <- streamy
    } else {
      message("Insufficient Data to Process")
      output <- NA
    }
  }
  
  return(output)
}


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





calculate_flow_vols_by_year = function(streamflow_data){
    # Ensure required columns exist
    required_cols <- c("year", "Q", "month", "doy")
    if (!all(required_cols %in% colnames(streamflow_data))) {
      missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
      stop(paste("Missing required columns:", paste(missing, collapse=", ")))
    }
    
    # Initialize results data frame
    result <- data.frame(matrix(NA, nrow=1, ncol=36))
    
    # Set column names
    periods <- c("Qann", "Qwin", "Qspr", "Qsum", "Qfal", 
                 "Q5", "Q10", "Q20", "Q50", "Q70", "Q90", "Q95")
    stats <- c("slp", "rho", "pval")
    colnames(result) <- paste0(rep(periods, each=3), "_", rep(stats, times=length(periods)))
    
    # Calculate annual means
    annual_means <- aggregate(Q ~ year, data=streamflow_data, FUN=mean, na.rm=TRUE)
    
    # Calculate seasonal means
    winter <- streamflow_data[streamflow_data$month %in% c(12, 1, 2), ]
    spring <- streamflow_data[streamflow_data$month %in% c(3, 4, 5), ]
    summer <- streamflow_data[streamflow_data$month %in% c(6, 7, 8), ]
    fall <- streamflow_data[streamflow_data$month %in% c(9, 10, 11), ]
    
    winter_means <- aggregate(Q ~ year, data=winter, FUN=mean, na.rm=TRUE)
    spring_means <- aggregate(Q ~ year, data=spring, FUN=mean, na.rm=TRUE)
    summer_means <- aggregate(Q ~ year, data=summer, FUN=mean, na.rm=TRUE)
    fall_means <- aggregate(Q ~ year, data=fall, FUN=mean, na.rm=TRUE)
    
    # Calculate flow percentiles by year
    calculate_percentile <- function(data, percentile) {
      agg <- aggregate(Q ~ year, data=data, 
                       FUN=function(x) quantile(x, probs=percentile/100, na.rm=TRUE))
      return(agg)
    }
    
    q5 <- calculate_percentile(streamflow_data, 5)
    q10 <- calculate_percentile(streamflow_data, 10)
    q20 <- calculate_percentile(streamflow_data, 20)
    q50 <- calculate_percentile(streamflow_data, 50)
    q70 <- calculate_percentile(streamflow_data, 70)
    q90 <- calculate_percentile(streamflow_data, 90)
    q95 <- calculate_percentile(streamflow_data, 95)
    
    # Function to calculate trend statistics
    calculate_trend_stats <- function(data) {
      if (nrow(data) < 3) {
        return(list(coef = c(NA, NA), estimate = NA, p.value = NA))
      }
      
      # Linear regression
      mod <- try(lm(Q ~ year, data=data), silent=TRUE)
      if (inherits(mod, "try-error")) {
        mod_coef <- c(NA, NA)
      } else {
        mod_coef <- coef(mod)
        if (length(mod_coef) < 2) mod_coef <- c(mod_coef, NA)
      }
      
      # Spearman correlation
      spearmans <- try(cor.test(data$year, data$Q, method="spearman"), silent=TRUE)
      if (inherits(spearmans, "try-error")) {
        spearmans <- list(estimate = NA, p.value = NA)
      }
      
      return(list(
        coef = mod_coef,
        estimate = spearmans$estimate,
        p.value = spearmans$p.value
      ))
    }
    
    # Calculate all trend statistics
    Qann_stats <- calculate_trend_stats(annual_means)
    Qwin_stats <- calculate_trend_stats(winter_means)
    Qspr_stats <- calculate_trend_stats(spring_means)
    Qsum_stats <- calculate_trend_stats(summer_means)
    Qfal_stats <- calculate_trend_stats(fall_means)
    Q5_stats <- calculate_trend_stats(q5)
    Q10_stats <- calculate_trend_stats(q10)
    Q20_stats <- calculate_trend_stats(q20)
    Q50_stats <- calculate_trend_stats(q50)
    Q70_stats <- calculate_trend_stats(q70)
    Q90_stats <- calculate_trend_stats(q90)
    Q95_stats <- calculate_trend_stats(q95)
    
    # Populate results
    result$Qann_slp <- Qann_stats$coef[2]
    result$Qann_rho <- Qann_stats$estimate
    result$Qann_pval <- Qann_stats$p.value
    
    result$Qwin_slp <- Qwin_stats$coef[2]
    result$Qwin_rho <- Qwin_stats$estimate
    result$Qwin_pval <- Qwin_stats$p.value
    
    result$Qspr_slp <- Qspr_stats$coef[2]
    result$Qspr_rho <- Qspr_stats$estimate
    result$Qspr_pval <- Qspr_stats$p.value
    
    result$Qsum_slp <- Qsum_stats$coef[2]
    result$Qsum_rho <- Qsum_stats$estimate
    result$Qsum_pval <- Qsum_stats$p.value
    
    result$Qfal_slp <- Qfal_stats$coef[2]
    result$Qfal_rho <- Qfal_stats$estimate
    result$Qfal_pval <- Qfal_stats$p.value
    
    result$Q5_slp <- Q5_stats$coef[2]
    result$Q5_rho <- Q5_stats$estimate
    result$Q5_pval <- Q5_stats$p.value
    
    result$Q10_slp <- Q10_stats$coef[2]
    result$Q10_rho <- Q10_stats$estimate
    result$Q10_pval <- Q10_stats$p.value
    
    result$Q20_slp <- Q20_stats$coef[2]
    result$Q20_rho <- Q20_stats$estimate
    result$Q20_pval <- Q20_stats$p.value
    
    result$Q50_slp <- Q50_stats$coef[2]
    result$Q50_rho <- Q50_stats$estimate
    result$Q50_pval <- Q50_stats$p.value
    
    result$Q70_slp <- Q70_stats$coef[2]
    result$Q70_rho <- Q70_stats$estimate
    result$Q70_pval <- Q70_stats$p.value
    
    result$Q90_slp <- Q90_stats$coef[2]
    result$Q90_rho <- Q90_stats$estimate
    result$Q90_pval <- Q90_stats$p.value
    
    result$Q95_slp <- Q95_stats$coef[2]
    result$Q95_rho <- Q95_stats$estimate
    result$Q95_pval <- Q95_stats$p.value
    
    return(result)
}


analyze_fdc_trends_from_streamflow <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("year", "Q")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Check if zyp package is available
  if (!requireNamespace("zyp", quietly = TRUE)) {
    stop("Package 'zyp' is needed for this function. Please install it with install.packages('zyp')")
  }
  
  # Initialize results data frame
  result <- data.frame(matrix(NA, nrow=1, ncol=9))
  
  # Set column names
  colnames(result) <- c(
    "FDCall_slp", "FDC90_slp", "FDCmid_slp",
    "FDCall_rho", "FDCall_pval",
    "FDC90th_rho", "FDC90th_pval",
    "FDCmid_rho", "FDCmid_pval"
  )
  
  # Calculate FDC characteristics by year
  years <- unique(streamflow_data$year)
  
  # Initialize FDC_byYear data frame
  FDC_byYear <- data.frame(
    year = years,
    slp_all = NA,
    slp_90th = NA,
    slp_mid = NA
  )
  
  # For each year, calculate FDC slopes
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$year == yr, ]
    
    # Need sufficient data points for the year
    if (nrow(year_data) < 30) {
      next
    }
    
    # Sort flows in descending order
    sorted_flows <- sort(year_data$Q, decreasing = TRUE)
    n <- length(sorted_flows)
    
    # Calculate exceedance probabilities
    exceedance <- (1:n) / (n + 1)
    
    # Create FDC data frame
    fdc <- data.frame(
      exceedance = exceedance,
      flow = sorted_flows
    )
    
    # Calculate slopes for different segments of the FDC
    # 1. Overall slope (all data)
    if (n >= 10) {
      # Use log-transformed flow for better fit
      log_flow <- log10(fdc$flow + 1e-10)  # Add small constant to handle zeros
      
      # Overall slope
      all_model <- try(lm(log_flow ~ exceedance, data=fdc), silent=TRUE)
      if (!inherits(all_model, "try-error")) {
        FDC_byYear$slp_all[FDC_byYear$year == yr] <- coef(all_model)[2]
      }
      
      # 2. Slope for 90th percentile and above (low flows)
      low_flow_data <- fdc[fdc$exceedance >= 0.9, ]
      if (nrow(low_flow_data) >= 3) {
        low_flow_model <- try(lm(log10(low_flow_data$flow + 1e-10) ~ low_flow_data$exceedance), silent=TRUE)
        if (!inherits(low_flow_model, "try-error")) {
          FDC_byYear$slp_90th[FDC_byYear$year == yr] <- coef(low_flow_model)[2]
        }
      }
      
      # 3. Slope for mid-range flows (20th to 80th percentile)
      mid_flow_data <- fdc[fdc$exceedance >= 0.2 & fdc$exceedance <= 0.8, ]
      if (nrow(mid_flow_data) >= 3) {
        mid_flow_model <- try(lm(log10(mid_flow_data$flow + 1e-10) ~ mid_flow_data$exceedance), silent=TRUE)
        if (!inherits(mid_flow_model, "try-error")) {
          FDC_byYear$slp_mid[FDC_byYear$year == yr] <- coef(mid_flow_model)[2]
        }
      }
    }
  }
  
  # Calculate Sen's slopes
  # For all FDC
  all_data <- subset(FDC_byYear, !is.na(FDC_byYear$slp_all))
  if (nrow(all_data) > 2) {
    sen_all <- zyp::zyp.sen(slp_all ~ year, all_data)
    result$FDCall_slp <- sen_all$coef[2]
  }
  
  # For 90th percentile FDC
  p90_data <- subset(FDC_byYear, !is.na(FDC_byYear$slp_90th))
  if (nrow(p90_data) > 2) {
    sen_90 <- zyp::zyp.sen(slp_90th ~ year, p90_data)
    result$FDC90_slp <- sen_90$coef[2]
  }
  
  # For mid-range FDC
  mid_data <- subset(FDC_byYear, !is.na(FDC_byYear$slp_mid))
  if (nrow(mid_data) > 2) {
    sen_mid <- zyp::zyp.sen(slp_mid ~ year, mid_data)
    result$FDCmid_slp <- sen_mid$coef[2]
  }
  
  # Calculate Spearman correlations
  # For all FDC
  if (nrow(all_data) > 2) {
    FDC_all_cor <- cor.test(all_data$year, all_data$slp_all, method='spearman')
    result$FDCall_rho <- FDC_all_cor$estimate
    result$FDCall_pval <- FDC_all_cor$p.value
  }
  
  # For 90th percentile FDC
  if (nrow(p90_data) > 2) {
    FDC_90th_cor <- cor.test(p90_data$year, p90_data$slp_90th, method='spearman')
    result$FDC90th_rho <- FDC_90th_cor$estimate
    result$FDC90th_pval <- FDC_90th_cor$p.value
  }
  
  # For mid-range FDC
  if (nrow(mid_data) > 2) {
    FDC_mid_cor <- cor.test(mid_data$year, mid_data$slp_mid, method='spearman')
    result$FDCmid_rho <- FDC_mid_cor$estimate
    result$FDCmid_pval <- FDC_mid_cor$p.value
  }
  
  # Add FDC_byYear as an attribute to the result
  attr(result, "FDC_byYear") <- FDC_byYear
  
  return(result)
}


analyze_flashiness_trends <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("year", "Q")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Check if zyp package is available
  if (!requireNamespace("zyp", quietly = TRUE)) {
    stop("Package 'zyp' is needed for this function. Please install it with install.packages('zyp')")
  }
  
  # Initialize results data frame
  result <- data.frame(matrix(NA, nrow=1, ncol=3))
  
  # Set column names
  colnames(result) <- c(
    "flashinessRB_slp", 
    "flashinessRB_rho", 
    "flashinessRB_pval"
  )
  
  # Calculate Richards-Baker flashiness index by year
  years <- unique(streamflow_data$year)
  
  # Initialize flashiness_byYear data frame
  flashiness_byYear <- data.frame(
    year = years,
    RB_index = NA
  )
  
  # For each year, calculate the RB flashiness index
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$year == yr, ]
    
    # Need sufficient data points for the year
    if (nrow(year_data) < 30) {
      next
    }
    
    # Sort by day to ensure chronological order
    if ("doy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$doy), ]
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
    flashiness_byYear$RB_index[flashiness_byYear$year == yr] <- rb_index
  }
  
  # Calculate Sen's slope for RB index trend
  rb_data <- subset(flashiness_byYear, !is.na(flashiness_byYear$RB_index))
  
  if (nrow(rb_data) > 2) {
    # Calculate Sen's slope
    sen_rb <- zyp::zyp.sen(RB_index ~ year, rb_data)
    result$flashinessRB_slp <- sen_rb$coef[2]
    
    # Calculate Spearman correlation
    flashiness_RB_cor <- cor.test(rb_data$year, rb_data$RB_index, method='spearman')
    result$flashinessRB_rho <- flashiness_RB_cor$estimate
    result$flashinessRB_pval <- flashiness_RB_cor$p.value
  }
  
  # Add flashiness_byYear as an attribute to the result
  attr(result, "flashiness_byYear") <- flashiness_byYear
  
  return(result)
}



analyze_flow_timing_trends <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("year", "Q", "doy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Check if zyp package is available
  if (!requireNamespace("zyp", quietly = TRUE)) {
    stop("Package 'zyp' is needed for this function. Please install it with install.packages('zyp')")
  }
  
  # Initialize results data frame
  percentiles <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  
  # Add historical_avg to col_types
  col_types <- c("slp", "rho", "pval", "historical_avg")
  result <- data.frame(matrix(NA, nrow=1, ncol=length(percentiles)*length(col_types)))
  
  # Set column names
  col_names <- c()
  for (p in percentiles) {
    for (type in col_types) {
      col_names <- c(col_names, paste0("D", p, "_day_", type))
    }
  }
  colnames(result) <- col_names
  
  # Create a data frame to store Julian days when cumulative flow reaches each percentile
  years <- unique(streamflow_data$year)
  julday_max <- data.frame(year=years)
  
  # For each year, find the day when cumulative flow reaches each percentile threshold
  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$year == yr, ]
    
    # Skip years with insufficient data
    if (nrow(year_data) < 300) {
      next
    }
    
    # Sort by day of year to ensure chronological order
    year_data <- year_data[order(year_data$doy), ]
    
    # Calculate total annual flow
    total_flow <- sum(year_data$Q, na.rm=TRUE)
    
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
        julday_max[julday_max$year == yr, paste0("D", p, "_day")] <- 
          year_data$doy[above_threshold[1]]
      } else {
        julday_max[julday_max$year == yr, paste0("D", p, "_day")] <- NA
      }
    }
  }
  
  # Calculate trends and historical averages for each percentile
  for (p in percentiles) {
    day_col <- paste0("D", p, "_day")
    
    # Subset to non-NA values
    valid_data <- subset(julday_max, !is.na(julday_max[[day_col]]))
    
    if (nrow(valid_data) > 2) {
      # Calculate Sen's slope
      formula <- as.formula(paste(day_col, "~ year"))
      day_mod <- zyp::zyp.sen(formula, valid_data)
      
      # Calculate Spearman correlation
      day_spearmans <- cor.test(valid_data$year, valid_data[[day_col]], method='spearman')
      
      # Calculate historical average (mean rounded to nearest integer)
      historical_avg <- round(mean(valid_data[[day_col]], na.rm=TRUE))
      
      # Store results
      result[[paste0("D", p, "_day_slp")]] <- day_mod$coef[2]
      result[[paste0("D", p, "_day_rho")]] <- day_spearmans$estimate
      result[[paste0("D", p, "_day_pval")]] <- day_spearmans$p.value
      result[[paste0("D", p, "_day_historical_avg")]] <- historical_avg
    }
  }
  
  # Add julday_max as an attribute to the result
  attr(result, "julday_max") <- julday_max
  
  return(result)
}




#############################################################
### as yet unimplemented code

calculate_median_flow_dates <- function(streamflow_data) {
  # Check if required columns exist
  required_cols <- c("year", "Q", "doy")
  if (!all(required_cols %in% colnames(streamflow_data))) {
    missing <- required_cols[!required_cols %in% colnames(streamflow_data)]
    stop(paste("Missing required columns:", paste(missing, collapse=", ")))
  }
  
  # Initialize results data frame
  result <- data.frame(matrix(NA, nrow=1, ncol=11))
  
  # Set column names for percentiles 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95
  percentiles <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  colnames(result) <- paste0("D", percentiles, "_med_date")
  
  # Calculate flow percentiles for the entire dataset
  flow_percentiles <- quantile(streamflow_data$Q, probs=percentiles/100, na.rm=TRUE)
  
  # Create a data frame to store Julian days when flow drops below each percentile
  julday_max <- data.frame(year=unique(streamflow_data$year))
  
  # For each year, find the day when flow drops below each percentile threshold
  for (yr in unique(streamflow_data$year)) {
    year_data <- streamflow_data[streamflow_data$year == yr, ]
    
    # Sort by day of year to ensure chronological order
    year_data <- year_data[order(year_data$doy), ]
    
    # For each percentile, find the first day when flow drops below the threshold
    for (p in percentiles) {
      threshold <- flow_percentiles[as.character(p/100)]
      
      # Find days where flow is below the threshold
      below_threshold <- which(year_data$Q <= threshold)
      
      # If there are days below threshold, take the first one
      if (length(below_threshold) > 0) {
        julday_max[julday_max$year == yr, paste0("D", p, "_day")] <- 
          year_data$doy[below_threshold[1]]
      } else {
        julday_max[julday_max$year == yr, paste0("D", p, "_day")] <- NA
      }
    }
  }
  
  # Calculate median dates for each percentile
  for (p in percentiles) {
    col_name <- paste0("D", p, "_day")
    result_col <- paste0("D", p, "_med_date")
    
    result[1, result_col] <- median(julday_max[[col_name]], na.rm=TRUE)
  }
  
  return(result)
}





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

