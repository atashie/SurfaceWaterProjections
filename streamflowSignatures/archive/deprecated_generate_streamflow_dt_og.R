# DEPRECATED: Original version of generate_streamflow_dt()
# Moved from helperFunctions.R during code review (March 2026)
# Superseded by the parquet-based pipeline in process_signatures_from_parquet()
# Contains the 99999 bug for Canadian gages (see CHANGELOG.md H5)

generate_streamflow_dt_og <- function(dt, data_origin,
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

      # Add water year information
      wy_info <- calculate_water_year_info(streamy$Date)
      streamy$water_year <- wy_info$water_year
      streamy$dowy <- wy_info$dowy

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

      # Add temporal information
      streamy$year = year(streamy$Date)
      streamy$month = month(streamy$Date)
      streamy$doy = yday(streamy$Date)

      # Add water year information
      wy_info <- calculate_water_year_info(streamy$Date)
      streamy$water_year <- wy_info$water_year
      streamy$dowy <- wy_info$dowy

      output <- streamy
    } else {
      message("Insufficient Data to Process")
      output <- NA
    }
  }

  return(output)
}
