# Set the main directory where data and results will be stored
# NOTE: This script expects to be run from the streamflowSignatures directory
# Paths are configured in config.R
source("config.R")
source("helperFunctions.R")
metadata_dir <- "D://gagesMetadata"#METADATA_DATA_DIR
output_dir <- "D://processedOuts_feb2026"

# Load configuration and helper functions




# Output file for results
output_file = file.path(output_dir, "summary_data_07feb2026.csv")

################################################################################################
# LOAD REQUIRED DATASETS
################################################################################################

# Check if required directories exist
if (!dir.exists(file.path(metadata_dir, ""))) {
  stop("Metadata directory not found. Please ensure the following path exists: ", 
       file.path(metadata_dir, "metadata"))
}

# Load USGS CONUS reference gages
tryCatch({
  conus_gages_raw = fread(file.path(metadata_dir, "", "conterm_bas_classif.txt"),
                          colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )#[CLASS=="Ref"]
  conus_basinid = fread(file.path(metadata_dir, "", "conterm_basinid.txt"), 
                        colClasses = c("STAID" = "character"))
  conus_gages_raw = fread(file.path(metadata_dir, "", "conterm_bas_classif.txt"),
                          colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )#[-c(1:8860),]#[CLASS=="Ref"]
  conus_basinid = fread(file.path(metadata_dir, "", "conterm_basinid.txt"), 
                        colClasses = c("STAID" = "character"))#[-c(1:8860),]
  
  conus_gages = merge(conus_gages_raw, conus_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(conus_gages), "CONUS reference gages\n")
}, error = function(e) {
  stop("Error loading CONUS gage data: ", e$message)
})

# Load USGS Alaska reference gages
tryCatch({
  AK_gages_all = fread(file.path(metadata_dir, "", "AKHIPR_bas_classif.txt"),
                       colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )[AGGECOREGION == 'Alaska']# & CLASS == 'Ref']
  AK_basinid = fread(file.path(metadata_dir, "", "AKHIPR_basinid.txt"), 
                     colClasses = c("STAID" = "character"))
  AK_gages = merge(AK_gages_all, AK_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(AK_gages), "Alaska reference gages\n")
}, error = function(e) {
  stop("Error loading Alaska gage data: ", e$message)
})

# Load Canadian reference gages
tryCatch({
  canadian_gages_goodData = fread(file.path(metadata_dir, "", "Canadian_gages_goodones.csv"))
  regulation_info = as.data.table(hy_stn_regulation(canadian_gages_goodData$STATION_NUMBER))
  canadian_gages = merge(canadian_gages_goodData, regulation_info, 
                         by = "STATION_NUMBER", all.x = TRUE)#[REGULATED != TRUE]
  cat("Loaded", nrow(canadian_gages), "Canadian reference gages\n")
}, error = function(e) {
  stop("Error loading Canadian gage data: ", e$message, 
       "\nNote: You need to install and configure the tidyhydat package for Canadian data.")
})




# Process gages
metadata <- process_gages_rawToRaw(
  gages_df = canadian_gages,#conus_gages,#AK_gages,
  gage_type = "Canada",
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  output_dir = file.path(getwd(), "processed_streamflow_data_canada"),
  storage_format = "parquet",
  chunk_size = 1000#,
#  resume=TRUE
)

# Process gages
metadata <- process_gages_rawToRaw(
  gages_df = AK_gages,
  gage_type = "USGS",
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  output_dir = file.path(getwd(), "processed_streamflow_data_AK"),
  chunk_size = 1000#,
  #  resume=TRUE
)

metadata <- process_gages_rawToRaw(
  gages_df = conus_gages,#AK_gages,
  gage_type = "USGS",
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  output_dir = file.path(getwd(), "processed_streamflow_data_conus"),
  storage_format = "parquet",
  chunk_size = 1000#,
  #  resume=TRUE
)




# Query specific watersheds
selected_data <- query_watersheds(
  output_dir = "processed_streamflow_data",
  gage_ids = c("01010070", "01010500"),
  storage_format = "parquet"
)




# Freading in downloaded data
# Set the output directory where parquet files are stored
output_dir <- "processed_streamflow_data_conus"
output_dir <- "processed_streamflow_data_canada"
output_dir <- "processed_streamflow_data_AK"

# Option 1: Read specific gages
#specific_gages <- AK_gages$STAID[1:10]# c("01010070", "01010500", "01011000")
#data_specific <- read_streamflow_data(output_dir, gage_ids = specific_gages)
data_specific <- read_streamflow_data(output_dir)

# Create time series plot
p1 <- plot_streamflow_timeseries(data_specific, 
                                 title = "Streamflow Time Series for Selected Gages",
                                 log_scale = TRUE)
print(p1)

# Option 2: Read random sample of gages
data_sample <- read_streamflow_data(output_dir, max_gages_to_plot = 6)

# Create faceted plot for better visibility with many gages
p2 <- plot_streamflow_faceted(data_sample, 
                              ncol = 3,
                              date_range = c(as.Date("2010-01-01"), 
                                             as.Date("2020-12-31")))
print(p2)

# Option 3: Plot annual patterns
p3 <- plot_annual_patterns(data_sample, aggregate_fun = "median")
print(p3)

# Save plots
#ggsave("streamflow_timeseries.png", p1, width = 12, height = 6, dpi = 300)
#ggsave("streamflow_faceted.png", p2, width = 14, height = 10, dpi = 300)
#ggsave("streamflow_annual_pattern.png", p3, width = 10, height = 6, dpi = 300)

# Additional analysis: Summary statistics by gage
summary_stats <- data_specific[, .(
  n_days = .N,
  mean_Q = mean(Q, na.rm = TRUE),
  median_Q = median(Q, na.rm = TRUE),
  min_Q = min(Q, na.rm = TRUE),
  max_Q = max(Q, na.rm = TRUE),
  sd_Q = sd(Q, na.rm = TRUE),
  start_date = min(Date),
  end_date = max(Date)
), by = gage_id]

print(summary_stats)















#######################
# concatenate files



# Usage example:
input_directories <- c(
  "processed_streamflow_data_AK",
  "processed_streamflow_data_canada", 
  "processed_streamflow_data_conus"
)

# Method 1: Simple concatenation
#concatenate_parquet_directories(
#  input_dirs = input_directories,
#  output_file = "all_streamflow_data_combined.parquet",
#  method = "duckdb"  # or "arrow"
#)

# Method 2: Concatenation with metadata handling
result <- concatenate_with_metadata(
  input_dirs = input_directories,
  output_dir = "combined_streamflow_output"
)


basinAt_NorAm_polys = st_read(file.path(metadata_dir, "geospatial_derivedData/basinAt_NorAm_polys.gpkg"))
basinAt_NorAm_strip = basinAt_NorAm_polys
st_geometry(basinAt_NorAm_strip) = NULL
HB_dt = data.table(basinAt_NorAm_strip)

# Identifying associated hydroBasins ids


# Usage - note the .rds extension:
updated_metadata <- add_downstream_basin_ids(
  metadata_file_path = "combined_streamflow_output/combined_watershed_metadata.csv",
  basinAt_NorAm_polys = basinAt_NorAm_polys,
  HB_dt = HB_dt,
  upstream_hydrobasins_path = "upstream_hydrobasins.rds"  # Use .rds extension
)

# To load the saved data later:
upstream_hydrobasins <- readRDS("upstream_hydrobasins.rds")























######################################################################################
######################################################################################
######################################################################################
## Aggregating subwatersheds to basins
# Load required libraries
library(terra)
library(sf)
library(data.table)
library(arrow)
library(dplyr)

# Function to process Daymet NetCDF files for watersheds
process_daymet_from_files <- function(
    watersheds_file,
    downloads_dir = DOWNLOADS_DIR,
    output_dir = "daymet_processed_data",
    variable = "prcp",  # Can be: prcp, tmin, tmax, srad, vp, swe, dayl
    max_area_sqkm = 10000,
    start_year = 1980,
    end_year = 2023,
    delete_after_processing = TRUE,
    verbose = TRUE) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load watersheds
  if (verbose) cat("Loading watersheds...\n")
  watersheds <- st_read(watersheds_file, quiet = !verbose)
  
  # Filter by size
  watersheds_filtered <- watersheds[watersheds$area_sqkm < max_area_sqkm, ]
  watersheds_filtered <- st_make_valid(watersheds_filtered)
  
  if (verbose) {
    cat(sprintf("Processing %d watersheds < %d sq km\n", 
                nrow(watersheds_filtered), max_area_sqkm))
  }
  
  # Initialize list to store all watershed data
  all_watershed_data <- list()
  
  # Process each watershed
  for (w_idx in seq_len(nrow(watersheds_filtered))) {
    watershed <- watersheds_filtered[w_idx, ]
    gage_id <- watershed$gage_id
    
    if (verbose && w_idx %% 10 == 0) {
      cat(sprintf("Processing watershed %d/%d: %s\n", 
                  w_idx, nrow(watersheds_filtered), gage_id))
    }
    
    # Transform watershed to Daymet projection
    daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
    watershed_proj <- st_transform(watershed, daymet_crs)
    
    # Initialize list for this watershed's data
    watershed_yearly_data <- list()
    
    # Process each year
    for (year in start_year:end_year) {
      
      # Construct filename
      nc_file <- file.path(downloads_dir, 
                           sprintf("daymet_v4_daily_na_%s_%d.nc", variable, year))
      
      if (!file.exists(nc_file)) {
        if (verbose) cat(sprintf("  Warning: File not found for year %d\n", year))
        next
      }
      
      tryCatch({
        # Read the NetCDF file as raster
        r <- terra::rast(nc_file)
        
        # Extract data for watershed
        extracted <- terra::extract(r, vect(watershed_proj), 
                                    fun = "mean", 
                                    weights = TRUE,
                                    touches = TRUE)
        
        # Get dates for this year
        if (year %% 4 == 0) {  # Leap year
          dates <- seq(as.Date(sprintf("%d-01-01", year)), 
                       as.Date(sprintf("%d-12-31", year)), 
                       by = "day")
        } else {
          dates <- seq(as.Date(sprintf("%d-01-01", year)), 
                       as.Date(sprintf("%d-12-31", year)), 
                       by = "day")
        }
        
        # Create data table for this year
        if (ncol(extracted) > 1) {
          values <- as.numeric(extracted[1, -1])  # Remove ID column
          
          year_data <- data.table(
            gage_id = gage_id,
            Date = dates[1:length(values)],
            value = values
          )
          
          watershed_yearly_data[[as.character(year)]] <- year_data
        }
        
      }, error = function(e) {
        if (verbose) {
          cat(sprintf("  Error processing year %d: %s\n", year, e$message))
        }
      })
    }
    
    # Combine all years for this watershed
    if (length(watershed_yearly_data) > 0) {
      watershed_data <- rbindlist(watershed_yearly_data)
      all_watershed_data[[gage_id]] <- watershed_data
    }
    
    # Periodic garbage collection
    if (w_idx %% 50 == 0) gc()
  }
  
  # Combine all watershed data
  if (verbose) cat("\nCombining all watershed data...\n")
  
  if (length(all_watershed_data) > 0) {
    combined_data <- rbindlist(all_watershed_data)
    
    # Pivot to wide format to match streamflow data structure
    wide_data <- dcast(combined_data, Date ~ gage_id, value.var = "value")
    
    # Sort by date
    setorder(wide_data, Date)
    
    # Save as parquet
    output_file <- file.path(output_dir, sprintf("combined_%s_data.parquet", variable))
    arrow::write_parquet(wide_data, output_file)
    
    if (verbose) {
      cat(sprintf("\nSaved combined data to: %s\n", output_file))
      cat(sprintf("Data dimensions: %d dates x %d watersheds\n", 
                  nrow(wide_data), ncol(wide_data) - 1))
    }
    
    # Delete NetCDF files if requested
    if (delete_after_processing) {
      if (verbose) cat("\nDeleting processed NetCDF files...\n")
      
      for (year in start_year:end_year) {
        nc_file <- file.path(downloads_dir, 
                             sprintf("daymet_v4_daily_na_%s_%d.nc", variable, year))
        if (file.exists(nc_file)) {
          unlink(nc_file)
        }
      }
    }
    
    return(output_file)
    
  } else {
    warning("No data was extracted")
    return(NULL)
  }
}

# Function to read and inspect the processed data
inspect_daymet_output <- function(output_dir, variable = "prcp") {
  
  file_path <- file.path(output_dir, sprintf("combined_%s_data.parquet", variable))
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the parquet file
  data <- arrow::read_parquet(file_path)
  
  cat(sprintf("\nData for %s:\n", variable))
  cat(sprintf("Date range: %s to %s\n", min(data$Date), max(data$Date)))
  cat(sprintf("Number of watersheds: %d\n", ncol(data) - 1))
  cat(sprintf("Number of days: %d\n", nrow(data)))
  
  # Check for missing values
  n_missing <- sum(is.na(data[, -1]))
  n_total <- (ncol(data) - 1) * nrow(data)
  cat(sprintf("Missing values: %d (%.2f%%)\n", n_missing, 100 * n_missing / n_total))
  
  # Sample data
  cat("\nFirst 5 dates and 5 watersheds:\n")
  print(data[1:5, 1:min(6, ncol(data))])
  
  return(data)
}

# Execute processing for precipitation
result_prcp <- process_daymet_from_files(
  watersheds_file = file.path(metadata_dir, "geospatial_derivedData/unified_watersheds_simplified.gpkg"),
  downloads_dir = DOWNLOADS_DIR,
  output_dir = file.path(metadata_dir, "daymet_processed_data"),
  variable = "prcp",
  max_area_sqkm = 10000,
  start_year = 1980,
  end_year = 2023,
  delete_after_processing = FALSE,  # Set to TRUE when ready
  verbose = TRUE
)


# Inspect the results
prcp_data <- inspect_daymet_output(
  output_dir = file.path(metadata_dir, "daymet_processed_data"),
  variable = "prcp"
)

