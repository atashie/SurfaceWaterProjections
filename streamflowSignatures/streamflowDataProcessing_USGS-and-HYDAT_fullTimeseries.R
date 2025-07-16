# Set the main directory where data and results will be stored
main_dir = "C:\\Users\\18033\\Documents\\GitHub\\SurfaceWaterProjections\\streamflowSignatures"#"C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"
metadata_dir = "D:/"
# Load helper functions
source(file.path(main_dir, "helperFunctions.R"), local=FALSE)

################################################################################################
# USER CONFIGURATION - MODIFY THESE SETTINGS AS NEEDED
################################################################################################


# Analysis period
start_date = as.Date("1979-10-01")  # Modern data period
end_date = as.Date("2025-06-01")    # End of analysis period

# Data quality thresholds
min_num_years = 30    # Minimum number of years required for analysis
min_nona_days = floor(365*.9)   # Minimum number of non-NA days per year
min_Q_value_and_days = c(0.0001, 30)  # Min flow value (mm) and days above this value

# Output file for results
output_file = file.path(metadata_dir, "processedOuts/summary_data.csv")

################################################################################################
# LOAD REQUIRED DATASETS
################################################################################################

# Check if required directories exist
if (!dir.exists(file.path(main_dir, "metadata"))) {
  stop("Metadata directory not found. Please ensure the following path exists: ", 
       file.path(main_dir, "metadata"))
}

# Load USGS CONUS reference gages
tryCatch({
  conus_gages_raw = fread(file.path(metadata_dir, "gagesMetadata", "conterm_bas_classif.txt"),
                          colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )#[CLASS=="Ref"]
  conus_basinid = fread(file.path(metadata_dir, "gagesMetadata", "conterm_basinid.txt"), 
                        colClasses = c("STAID" = "character"))
  conus_gages_raw = fread(file.path(metadata_dir, "gagesMetadata", "conterm_bas_classif.txt"),
                          colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )#[-c(1:8860),]#[CLASS=="Ref"]
  conus_basinid = fread(file.path(metadata_dir, "gagesMetadata", "conterm_basinid.txt"), 
                        colClasses = c("STAID" = "character"))#[-c(1:8860),]
  
  conus_gages = merge(conus_gages_raw, conus_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(conus_gages), "CONUS reference gages\n")
}, error = function(e) {
  stop("Error loading CONUS gage data: ", e$message)
})

# Load USGS Alaska reference gages
tryCatch({
  AK_gages_all = fread(file.path(metadata_dir, "gagesMetadata", "AKHIPR_bas_classif.txt"),
                       colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )[AGGECOREGION == 'Alaska']# & CLASS == 'Ref']
  AK_basinid = fread(file.path(metadata_dir, "gagesMetadata", "AKHIPR_basinid.txt"), 
                     colClasses = c("STAID" = "character"))
  AK_gages = merge(AK_gages_all, AK_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(AK_gages), "Alaska reference gages\n")
}, error = function(e) {
  stop("Error loading Alaska gage data: ", e$message)
})

# Load Canadian reference gages
tryCatch({
  canadian_gages_goodData = fread(file.path(metadata_dir, "gagesMetadata", "Canadian_gages_goodones.csv"))
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
  min_num_years = 30,
  start_date = as.Date("1979-10-01"),
  end_date = as.Date("2025-06-01"),
  min_Q_value_and_days = c(0.001, 328),
  output_dir = "processed_streamflow_data_canada",
  storage_format = "parquet",
  chunk_size = 1000,
  resume=TRUE
)

metadata <- process_gages_rawToRaw(
  gages_df = conus_gages,#AK_gages,
  gage_type = "USGS",
  min_num_years = 30,
  start_date = as.Date("1979-10-01"),
  end_date = as.Date("2025-06-01"),
  min_Q_value_and_days = c(0.001, 328),
  output_dir = "processed_streamflow_data_conus",
  storage_format = "parquet",
  chunk_size = 1000,
  resume=TRUE
)

# Process gages
metadata <- process_gages_rawToRaw(
  gages_df = AK_gages,
  gage_type = "USGS",
  min_num_years = 30,
  start_date = as.Date("1979-10-01"),
  end_date = as.Date("2025-06-01"),
  min_Q_value_and_days = c(0.001, 328),
  output_dir = "processed_streamflow_data_AK",
  storage_format = "parquet",
  chunk_size = 1000
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
library(sf)
library(data.table)
library(rmapshaper)
library(dplyr)

# Function to create unified watersheds with adaptive simplification
create_unified_watersheds <- function(upstream_hydrobasins_path,
                                      basinAt_NorAm_polys,
                                      output_file = "unified_watersheds.gpkg",
                                      verbose = TRUE) {
  
  # Load the upstream hydrobasins mapping
  if (verbose) cat("Loading upstream hydrobasins data...\n")
  upstream_hydrobasins <- readRDS(upstream_hydrobasins_path)
  
  # Initialize list to store processed watersheds
  all_watersheds <- list()
  
  # Get unique gage IDs
  gage_ids <- names(upstream_hydrobasins)
  
  if (verbose) cat(sprintf("Processing %d unique gage IDs...\n", length(gage_ids)))
  
  # Process each gage ID
  for (i in seq_along(gage_ids)) {
    gage_id <- gage_ids[i]
    
    if (verbose && i %% 100 == 0) {
      cat(sprintf("Processing gage %d of %d (%.1f%%)...\n", 
                  i, length(gage_ids), 100 * i / length(gage_ids)))
    }
    
    # Get basin IDs for this gage
    basin_ids <- upstream_hydrobasins[[gage_id]]
    
    # Skip if no basins associated
    if (length(basin_ids) == 0) next
    
    # Extract polygons for these basin IDs
    gage_polygons <- basinAt_NorAm_polys[basinAt_NorAm_polys$HYBAS_ID %in% basin_ids, ]
    
    # Skip if no polygons found
    if (nrow(gage_polygons) == 0) next
    
    # Merge polygons into single multipolygon
    merged_watershed <- st_union(gage_polygons)
    
    # Apply simplification if more than 2 basins
    n_basins <- length(basin_ids)
    if (n_basins > 2) {
      # Calculate simplification factor
      simplify_factor <- n_basins / 500
      
      # Use rmapshaper for topology-preserving simplification
      # The 'keep' parameter represents the proportion of vertices to keep
      # So we need to convert our simplification factor to a keep proportion
      keep_proportion <- min(1, 500 / n_basins)  # Inverse of simplify_factor
      
      tryCatch({
        merged_watershed <- ms_simplify(merged_watershed, 
                                        keep = keep_proportion,
                                        keep_shapes = TRUE,
                                        method = "dp")  # Douglas-Peucker algorithm
      }, error = function(e) {
        if (verbose) cat(sprintf("Warning: Simplification failed for gage %s: %s\n", 
                                 gage_id, e$message))
      })
    }
    
    # Create sf object with attributes
    watershed_sf <- st_sf(
      gage_id = gage_id,
      n_basins = n_basins,
      basin_ids = paste(basin_ids, collapse = ";"),
      area_sqkm = as.numeric(st_area(merged_watershed)) / 1e6,
      simplified = n_basins > 2,
      simplify_factor = ifelse(n_basins > 2, n_basins / 500, NA),
      geometry = merged_watershed
    )
    
    # Add to list
    all_watersheds[[i]] <- watershed_sf
  }
  
  # Combine all watersheds
  if (verbose) cat("Combining all watersheds...\n")
  combined_watersheds <- do.call(rbind, all_watersheds[!sapply(all_watersheds, is.null)])
  
  # Sort by polygon size (area)
  if (verbose) cat("Sorting by polygon size...\n")
  combined_watersheds <- combined_watersheds[order(combined_watersheds$area_sqkm, decreasing = TRUE), ]
  
  # Add rank column
  combined_watersheds$size_rank <- seq_len(nrow(combined_watersheds))
  
  # Calculate some summary statistics
  if (verbose) {
    cat("\nSummary statistics:\n")
    cat(sprintf("Total watersheds processed: %d\n", nrow(combined_watersheds)))
    cat(sprintf("Watersheds simplified: %d\n", sum(combined_watersheds$simplified, na.rm = TRUE)))
    cat(sprintf("Area range: %.2f - %.2f sq km\n", 
                min(combined_watersheds$area_sqkm), 
                max(combined_watersheds$area_sqkm)))
    cat(sprintf("Basins per watershed range: %d - %d\n",
                min(combined_watersheds$n_basins),
                max(combined_watersheds$n_basins)))
  }
  
  # Save to file
  if (verbose) cat(sprintf("\nSaving to %s...\n", output_file))
  st_write(combined_watersheds, output_file, delete_dsn = TRUE, quiet = !verbose)
  
  return(combined_watersheds)
}

# Execute the function
unified_watersheds <- create_unified_watersheds(
  upstream_hydrobasins_path = file.path(main_dir, "upstream_hydrobasins.rds"),
  basinAt_NorAm_polys = basinAt_NorAm_polys,
  output_file = file.path(metadata_dir, "geospatial_derivedData/unified_watersheds_simplified.gpkg"),
  verbose = TRUE
)

# Optional: Create a quick visualization of the largest watersheds
library(ggplot2)
library(viridis)

# Plot top 20 largest watersheds
top_watersheds <- unified_watersheds[1:min(20, nrow(unified_watersheds)), ]

p <- ggplot() +
  geom_sf(data = top_watersheds, 
          aes(fill = log10(area_sqkm)), 
          color = "black", 
          size = 0.2) +
  scale_fill_viridis_c(name = "Log10 Area\n(sq km)") +
  theme_minimal() +
  labs(title = "Top 20 Largest Unified Watersheds",
       subtitle = sprintf("Colored by area; Total watersheds: %d", 
                          nrow(unified_watersheds))) +
  theme(legend.position = "right")
p


# Create summary table
summary_table <- unified_watersheds %>%
  st_drop_geometry() %>%
  group_by(simplified) %>%
  summarise(
    count = n(),
    mean_area = mean(area_sqkm),
    median_area = median(area_sqkm),
    mean_n_basins = mean(n_basins),
    max_n_basins = max(n_basins)
  )

print(summary_table)



## END: Aggregating subwatersheds to basins
######################################################################################
######################################################################################
######################################################################################








######################################################################################
######################################################################################
######################################################################################
## Aggregating daymet climate data
# Test function to check Daymet data availability and debug issues
test_daymet_availability <- function(location, start_year = 2020, end_year = 2024, var = "tmin") {
  cat("\nTesting Daymet data availability...\n")
  cat(sprintf("Location: %.4f, %.4f, %.4f, %.4f\n", location[1], location[2], location[3], location[4]))
  
  # Test each year individually to find which ones work
  for (year in start_year:end_year) {
    cat(sprintf("\nTesting year %d...\n", year))
    
    temp_dir <- tempdir()
    test_file <- file.path(temp_dir, sprintf("test_%d.nc", year))
    
    tryCatch({
      result <- download_daymet_ncss(
        location = location,
        start = year,
        end = year,
        param = var,
        path = temp_dir,
        silent = FALSE
      )
      
      # Check if file exists and size
      if (file.exists(test_file)) {
        size_mb <- file.size(test_file) / 1024^2
        cat(sprintf("  Success! File size: %.2f MB\n", size_mb))
        unlink(test_file)
      }
      
    }, error = function(e) {
      cat(sprintf("  Failed: %s\n", e$message))
    })
  }
}

# Corrected function with proper year handling
aggregate_daymet_for_watersheds_corrected <- function(
    watersheds_file,
    output_dir = "daymet_watershed_timeseries",
    max_area_sqkm = 10000,
    daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
    start_year = 1980,
    end_year = 2023,  # Changed from 2024 to 2023 - Daymet 2024 may not be available
    batch_degrees = 0.25,
    verbose = TRUE) {
  
  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "by_variable"), showWarnings = FALSE)
  
  # Create temporary directory
  temp_nc_dir <- file.path(tempdir(), paste0("daymet_nc_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_nc_dir, showWarnings = FALSE)
  
  # Load watersheds
  if (verbose) cat("Loading watersheds...\n")
  watersheds <- st_read(watersheds_file, quiet = !verbose)
  
  # Filter watersheds by size
  watersheds_filtered <- watersheds[watersheds$area_sqkm < max_area_sqkm, ]
  
  if (verbose) {
    cat(sprintf("Processing %d watersheds (out of %d total) smaller than %d sq km\n",
                nrow(watersheds_filtered), nrow(watersheds), max_area_sqkm))
  }
  
  # Transform to lat/lon
  watersheds_ll <- st_transform(watersheds_filtered, 4326)
  
  # Initialize storage
  all_results <- list()
  for (var in daymet_variables) {
    all_results[[var]] <- list()
  }
  
  # First, let's test with one watershed to check data availability
  if (verbose) {
    cat("\nTesting data availability with first watershed...\n")
    test_bbox <- st_bbox(watersheds_ll[1, ])
    test_location <- c(test_bbox["ymax"], test_bbox["xmin"], test_bbox["ymin"], test_bbox["xmax"])
    
    # Test what years are actually available
    available_years <- c()
    for (test_year in c(2020:2024)) {
      tryCatch({
        test_dir <- file.path(temp_nc_dir, "availability_test")
        dir.create(test_dir, showWarnings = FALSE)
        
        download_daymet_ncss(
          location = test_location,
          start = test_year,
          end = test_year,
          param = "tmin",
          path = test_dir,
          silent = TRUE
        )
        
        available_years <- c(available_years, test_year)
        cat(sprintf("  Year %d: Available\n", test_year))
        
      }, error = function(e) {
        cat(sprintf("  Year %d: Not available\n", test_year))
      }, finally = {
        unlink(test_dir, recursive = TRUE)
      })
    }
    
    # Adjust end_year based on availability
    if (length(available_years) > 0) {
      actual_end_year <- min(end_year, max(available_years))
      if (actual_end_year < end_year) {
        cat(sprintf("\nAdjusting end year from %d to %d based on data availability\n", 
                    end_year, actual_end_year))
        end_year <- actual_end_year
      }
    }
  }
  
  # Process each watershed
  for (i in seq_len(nrow(watersheds_ll))) {
    gage_id <- watersheds_ll$gage_id[i]
    area_sqkm <- watersheds_ll$area_sqkm[i]
    
    if (verbose) {
      cat(sprintf("\nProcessing watershed %d/%d: %s (%.1f sq km)...\n",
                  i, nrow(watersheds_ll), gage_id, area_sqkm))
    }
    
    # Get bounding box
    bbox <- st_bbox(watersheds_ll[i, ])
    lat_range <- bbox["ymax"] - bbox["ymin"]
    lon_range <- bbox["xmax"] - bbox["xmin"]
    
    # Calculate spatial batches
    lat_batches <- ceiling(lat_range / batch_degrees)
    lon_batches <- ceiling(lon_range / batch_degrees)
    total_batches <- lat_batches * lon_batches
    
    if (verbose && total_batches > 1) {
      cat(sprintf("  Watershed extent: %.3f x %.3f degrees\n", lat_range, lon_range))
      cat(sprintf("  Will process in %d spatial batches (%d x %d)\n",
                  total_batches, lat_batches, lon_batches))
    }
    
    # Process each variable
    for (var_idx in seq_along(daymet_variables)) {
      var <- daymet_variables[var_idx]
      
      if (verbose) cat(sprintf("  Variable %d/%d: %s\n", var_idx, length(daymet_variables), var))
      
      # Collect data for all spatial batches
      var_batch_data <- list()
      
      # Process spatial batches
      batch_num <- 0
      for (lat_batch in seq_len(lat_batches)) {
        for (lon_batch in seq_len(lon_batches)) {
          batch_num <- batch_num + 1
          
          # Calculate batch bounds
          lat_min <- bbox["ymin"] + (lat_batch - 1) * batch_degrees
          lat_max <- min(bbox["ymax"], bbox["ymin"] + lat_batch * batch_degrees)
          lon_min <- bbox["xmin"] + (lon_batch - 1) * batch_degrees
          lon_max <- min(bbox["xmax"], bbox["xmin"] + lon_batch * batch_degrees)
          
          if (verbose && total_batches > 1) {
            cat(sprintf("    Batch %d/%d: [%.4f,%.4f] x [%.4f,%.4f]\n", 
                        batch_num, total_batches, lat_min, lat_max, lon_min, lon_max))
          }
          
          # Create unique directory
          batch_dir <- file.path(temp_nc_dir, sprintf("%s_%s_b%d", gage_id, var, batch_num))
          dir.create(batch_dir, showWarnings = FALSE)
          
          tryCatch({
            # Small delay between requests
            if (batch_num > 1) Sys.sleep(1)
            
            # Download full period
            download_daymet_ncss(
              location = c(lat_max, lon_min, lat_min, lon_max),
              start = start_year,
              end = end_year,
              param = var,
              frequency = "daily",
              path = batch_dir,
              silent = TRUE
            )
            
            # Find downloaded file - it will be named with the variable and years
            nc_pattern <- sprintf("%s_.*\\.nc$", var)
            nc_files <- list.files(batch_dir, pattern = nc_pattern, full.names = TRUE)
            
            if (length(nc_files) > 0) {
              # Use the largest file (in case multiple were created)
              file_sizes <- file.info(nc_files)$size
              nc_file <- nc_files[which.max(file_sizes)]
              
              if (verbose) {
                cat(sprintf("      Downloaded: %.2f MB\n", file.size(nc_file) / 1024^2))
              }
              
              # Read and process
              r <- terra::rast(nc_file)
              
              # Transform watershed to Daymet projection
              daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
              watershed_proj <- st_transform(watersheds_ll[i, ], daymet_crs)
              
              # Extract mean values
              extracted <- terra::extract(r, vect(watershed_proj), 
                                          fun = "mean", 
                                          weights = TRUE,
                                          touches = TRUE)
              
              # Get dates
              dates <- time(r)
              if (is.null(dates) || length(dates) == 0) {
                # Generate full date sequence
                n_days <- nlyr(r)
                dates <- seq(as.Date(paste0(start_year, "-01-01")), 
                             by = "day", 
                             length.out = n_days)
              }
              
              # Create data table
              if (ncol(extracted) > 1) {
                values <- as.numeric(extracted[1, -1])
                batch_data <- data.table(
                  date = dates[1:length(values)],
                  value = values
                )
                var_batch_data[[length(var_batch_data) + 1]] <- batch_data
              }
            }
            
          }, error = function(e) {
            if (verbose) {
              cat(sprintf("      Error: %s\n", e$message))
            }
          }, finally = {
            # Clean up
            unlink(batch_dir, recursive = TRUE)
            if (exists("r")) rm(r)
            gc(verbose = FALSE)
          })
        }
      }
      
      # Combine spatial batches
      if (length(var_batch_data) > 0) {
        var_data <- rbindlist(var_batch_data)
        # Average overlapping dates
        var_data <- var_data[, .(value = mean(value, na.rm = TRUE)), by = date]
        var_data[, gage_id := gage_id]
        all_results[[var]][[gage_id]] <- var_data
      }
    }
    
    # Periodic cleanup
    if (i %% 5 == 0) gc(verbose = FALSE)
  }
  
  # Save results
  if (verbose) cat("\nSaving results to parquet files...\n")
  
  for (var in names(all_results)) {
    if (length(all_results[[var]]) > 0) {
      # Combine all watersheds
      combined_data <- rbindlist(all_results[[var]])
      
      # Pivot to wide format
      wide_data <- dcast(combined_data, date ~ gage_id, value.var = "value")
      
      # Sort by date
      setorder(wide_data, date)
      
      # Save
      output_file <- file.path(output_dir, "by_variable", sprintf("daymet_%s_daily.parquet", var))
      arrow::write_parquet(wide_data, output_file)
      
      if (verbose) {
        cat(sprintf("  Saved %s: %d dates x %d watersheds (%.1f MB)\n", 
                    var, nrow(wide_data), ncol(wide_data) - 1,
                    file.size(output_file) / 1024^2))
      }
    }
  }
  
  # Clean up
  unlink(temp_nc_dir, recursive = TRUE)
  
  if (verbose) cat("\nProcessing complete!\n")
  
  return(output_dir)
}

# Execute with corrected parameters
result <- aggregate_daymet_for_watersheds_corrected(
  watersheds_file = file.path(metadata_dir, "geospatial_derivedData/unified_watersheds_simplified.gpkg"),
  output_dir = file.path(metadata_dir, "daymet_watershed_timeseries"),
  max_area_sqkm = 10000,
  daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
  start_year = 1980,
  end_year = 2023,  # Avoiding 2024 which may not be available
  batch_degrees = 0.25,
  verbose = TRUE
)










































































# Load required libraries
library(daymetr)
library(sf)
library(terra)
library(ncdf4)
library(arrow)
library(data.table)
library(dplyr)

# Function to aggregate Daymet data for watersheds with improved file handling
# Full implementation for processing Daymet data
aggregate_daymet_for_watersheds_full <- function(
    watersheds_file,
    output_dir = "daymet_aggregated_data",
    max_area_sqkm = 10000,
    daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
    start_year = 1980,
    end_year = 2024,
    batch_degrees = 0.25,
    verbose = TRUE) {
  
  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "by_variable"), showWarnings = FALSE)
  
  # Create temporary directory for NetCDF files
  temp_nc_dir <- file.path(tempdir(), paste0("daymet_nc_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_nc_dir, showWarnings = FALSE)
  
  # Load watersheds
  if (verbose) cat("Loading watersheds...\n")
  watersheds <- st_read(watersheds_file, quiet = !verbose)
  
  # Filter watersheds by size
  watersheds_filtered <- watersheds[watersheds$area_sqkm < max_area_sqkm, ]
  
  if (verbose) {
    cat(sprintf("Processing %d watersheds (out of %d total) smaller than %d sq km\n",
                nrow(watersheds_filtered), nrow(watersheds), max_area_sqkm))
  }
  
  # Transform to lat/lon for Daymet queries
  watersheds_ll <- st_transform(watersheds_filtered, 4326)
  
  # Initialize storage for all variables
  all_results <- list()
  for (var in daymet_variables) {
    all_results[[var]] <- list()
  }
  
  # Process each watershed
  for (i in seq_len(nrow(watersheds_ll))) {
    gage_id <- watersheds_ll$gage_id[i]
    area_sqkm <- watersheds_ll$area_sqkm[i]
    
    if (verbose) {
      cat(sprintf("\nProcessing watershed %d/%d: %s (%.1f sq km)...\n",
                  i, nrow(watersheds_ll), gage_id, area_sqkm))
    }
    
    # Get bounding box
    bbox <- st_bbox(watersheds_ll[i, ])
    lat_range <- bbox["ymax"] - bbox["ymin"]
    lon_range <- bbox["xmax"] - bbox["xmin"]
    
    # Calculate number of batches needed
    lat_batches <- ceiling(lat_range / batch_degrees)
    lon_batches <- ceiling(lon_range / batch_degrees)
    total_batches <- lat_batches * lon_batches
    
    if (verbose && total_batches > 1) {
      cat(sprintf("  Large watershed - will download in %d batches (%d x %d)\n",
                  total_batches, lat_batches, lon_batches))
    }
    
    # Process each variable
    for (var_idx in seq_along(daymet_variables)) {
      var <- daymet_variables[var_idx]
      
      if (verbose) cat(sprintf("  Variable %d/%d: %s\n", var_idx, length(daymet_variables), var))
      
      # Collect data for all years
      var_data_all_years <- list()
      
      # Process each year
      for (year in start_year:end_year) {
        
        # Show progress every 5 years
        if (verbose && (year - start_year) %% 5 == 0) {
          cat(sprintf("    Years %d-%d...\n", year, min(year + 4, end_year)))
        }
        
        # Collect data for all batches in this year
        year_batch_data <- list()
        
        # Process spatial batches
        batch_num <- 0
        for (lat_batch in seq_len(lat_batches)) {
          for (lon_batch in seq_len(lon_batches)) {
            batch_num <- batch_num + 1
            
            # Calculate batch bounds
            lat_min <- bbox["ymin"] + (lat_batch - 1) * batch_degrees
            lat_max <- min(bbox["ymax"], bbox["ymin"] + lat_batch * batch_degrees)
            lon_min <- bbox["xmin"] + (lon_batch - 1) * batch_degrees
            lon_max <- min(bbox["xmax"], bbox["xmin"] + lon_batch * batch_degrees)
            
            # Create unique directory for this batch
            batch_dir <- file.path(temp_nc_dir, sprintf("%s_%s_%d_b%d", gage_id, var, year, batch_num))
            dir.create(batch_dir, showWarnings = FALSE)
            
            tryCatch({
              # Small delay to avoid overwhelming server
              if (batch_num > 1) Sys.sleep(0.5)
              
              # Download data
              download_daymet_ncss(
                location = c(lat_max, lon_min, lat_min, lon_max),
                start = year,
                end = year,
                param = var,
                frequency = "daily",
                path = batch_dir,
                silent = TRUE
              )
              
              # Find downloaded file
              nc_files <- list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE)
              
              if (length(nc_files) > 0 && file.exists(nc_files[1])) {
                # Read and process
                r <- terra::rast(nc_files[1])
                
                # Transform watershed to Daymet projection
                daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
                watershed_proj <- st_transform(watersheds_ll[i, ], daymet_crs)
                
                # Extract mean values for watershed
                extracted <- terra::extract(r, vect(watershed_proj), 
                                            fun = "mean", 
                                            weights = TRUE,
                                            touches = TRUE)
                
                # Get dates
                dates <- time(r)
                if (is.null(dates) || length(dates) == 0) {
                  # Generate dates for the year
                  n_days <- nlyr(r)
                  dates <- seq(as.Date(paste0(year, "-01-01")), 
                               by = "day", 
                               length.out = n_days)
                }
                
                # Create data table
                if (ncol(extracted) > 1) {
                  values <- as.numeric(extracted[1, -1])
                  batch_data <- data.table(
                    date = dates[1:length(values)],
                    value = values
                  )
                  year_batch_data[[length(year_batch_data) + 1]] <- batch_data
                }
              }
              
            }, error = function(e) {
              # Silent fail for individual batches
            }, finally = {
              # Always clean up
              unlink(batch_dir, recursive = TRUE)
              if (exists("r")) rm(r)
              gc(verbose = FALSE)
            })
          }
        }
        
        # Combine batches for this year
        if (length(year_batch_data) > 0) {
          year_data <- rbindlist(year_batch_data)
          # Average overlapping dates
          year_data <- year_data[, .(value = mean(value, na.rm = TRUE)), by = date]
          var_data_all_years[[length(var_data_all_years) + 1]] <- year_data
        }
      }
      
      # Combine all years for this variable
      if (length(var_data_all_years) > 0) {
        var_data_complete <- rbindlist(var_data_all_years)
        var_data_complete[, gage_id := gage_id]
        all_results[[var]][[gage_id]] <- var_data_complete
      }
    }
    
    # Periodic cleanup
    if (i %% 5 == 0) gc(verbose = FALSE)
  }
  
  # Save results by variable
  if (verbose) cat("\nSaving results to parquet files...\n")
  
  for (var in names(all_results)) {
    if (length(all_results[[var]]) > 0) {
      # Combine all watersheds for this variable
      combined_data <- rbindlist(all_results[[var]])
      
      # Pivot to wide format
      wide_data <- dcast(combined_data, date ~ gage_id, value.var = "value")
      
      # Sort by date
      setorder(wide_data, date)
      
      # Save
      output_file <- file.path(output_dir, "by_variable", sprintf("daymet_%s_daily.parquet", var))
      arrow::write_parquet(wide_data, output_file)
      
      if (verbose) {
        cat(sprintf("  Saved %s: %d dates x %d watersheds (%.1f MB)\n", 
                    var, nrow(wide_data), ncol(wide_data) - 1,
                    file.size(output_file) / 1024^2))
      }
    }
  }
  
  # Clean up
  unlink(temp_nc_dir, recursive = TRUE)
  
  # Create summary
  create_processing_summary_full(output_dir, watersheds_filtered, daymet_variables)
  
  if (verbose) cat("\nProcessing complete!\n")
  
  return(output_dir)
}

# Enhanced summary function
create_processing_summary_full <- function(output_dir, watersheds, variables) {
  
  # Get saved files
  var_files <- list.files(file.path(output_dir, "by_variable"), 
                          pattern = "daymet_.*_daily\\.parquet$", 
                          full.names = TRUE)
  
  # Create summary of each variable file
  file_summary <- data.frame()
  
  for (file in var_files) {
    var_name <- gsub(".*daymet_(.*)_daily\\.parquet", "\\1", basename(file))
    data <- arrow::read_parquet(file)
    
    file_info <- data.frame(
      variable = var_name,
      n_watersheds = ncol(data) - 1,
      n_dates = nrow(data),
      date_range = paste(min(data$date), "to", max(data$date)),
      file_size_mb = round(file.size(file) / 1024^2, 2),
      stringsAsFactors = FALSE
    )
    
    file_summary <- rbind(file_summary, file_info)
  }
  
  # Save summary
  write.csv(file_summary, file.path(output_dir, "processing_summary.csv"), row.names = FALSE)
  
  # Create detailed metadata
  metadata <- list(
    processing_date = Sys.time(),
    total_watersheds_input = nrow(watersheds),
    watersheds_processed = length(unique(unlist(lapply(var_files, function(f) {
      names(arrow::read_parquet(f))[-1]
    })))),
    variables = variables,
    file_summary = file_summary
  )
  
  saveRDS(metadata, file.path(output_dir, "processing_metadata.rds"))
  
  return(metadata)
}

# Execute the full processing
result <- aggregate_daymet_for_watersheds_full(
  watersheds_file = file.path(metadata_dir, "geospatial_derivedData/unified_watersheds_simplified.gpkg"),
  output_dir = file.path(metadata_dir, "daymet_watershed_timeseries"),
  max_area_sqkm = 10000,
  daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
  start_year = 1980,
  end_year = 2024,
  batch_degrees = 0.25,
  verbose = TRUE
)

# Function to read and inspect results
inspect_daymet_results <- function(output_dir, variable = "prcp") {
  
  file_path <- file.path(output_dir, "by_variable", sprintf("daymet_%s_daily.parquet", variable))
  
  if (!file.exists(file_path)) {
    cat("Available variables:\n")
    files <- list.files(file.path(output_dir, "by_variable"), pattern = "daymet_.*_daily\\.parquet$")
    vars <- gsub("daymet_(.*)_daily\\.parquet", "\\1", files)
    cat(paste(" -", vars, collapse = "\n"))
    return(NULL)
  }
  
  data <- arrow::read_parquet(file_path)
  
  cat(sprintf("\nData for %s:\n", variable))
  cat(sprintf("Date range: %s to %s\n", min(data$date), max(data$date)))
  cat(sprintf("Number of watersheds: %d\n", ncol(data) - 1))
  cat(sprintf("Number of days: %d\n", nrow(data)))
  cat(sprintf("File size: %.1f MB\n", file.size(file_path) / 1024^2))
  
  # Sample data
  cat("\nFirst few dates and watersheds:\n")
  print(data[1:5, 1:min(6, ncol(data))])
  
  return(data)
}

# Inspect results
prcp_data <- inspect_daymet_results(result, "prcp")















































# Updated helper function to handle year-specific processing
# Load required libraries
library(daymetr)
library(sf)
library(terra)
library(ncdf4)
library(arrow)
library(data.table)
library(dplyr)

# Modified function that processes in smaller year chunks
aggregate_daymet_for_watersheds_chunked <- function(
    watersheds_file,
    output_dir = "daymet_watershed_timeseries",
    max_area_sqkm = 10000,
    daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
    start_year = 1980,
    end_year = 2024,
    years_per_chunk = 10,  # Process 10 years at a time
    batch_degrees = 0.25,
    verbose = TRUE) {
  
  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "by_variable"), showWarnings = FALSE)
  
  # Create temporary directory for NetCDF files
  temp_nc_dir <- file.path(tempdir(), paste0("daymet_nc_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_nc_dir, showWarnings = FALSE)
  
  # Load watersheds
  if (verbose) cat("Loading watersheds...\n")
  watersheds <- st_read(watersheds_file, quiet = !verbose)
  
  # Filter watersheds by size
  watersheds_filtered <- watersheds[watersheds$area_sqkm < max_area_sqkm, ]
  
  if (verbose) {
    cat(sprintf("Processing %d watersheds (out of %d total) smaller than %d sq km\n",
                nrow(watersheds_filtered), nrow(watersheds), max_area_sqkm))
  }
  
  # Transform to lat/lon for Daymet queries
  watersheds_ll <- st_transform(watersheds_filtered, 4326)
  
  # Initialize storage for all variables
  all_results <- list()
  for (var in daymet_variables) {
    all_results[[var]] <- list()
  }
  
  # Create year chunks
  year_chunks <- split(start_year:end_year, 
                       ceiling(seq_along(start_year:end_year) / years_per_chunk))
  
  # Process each watershed
  for (i in seq_len(nrow(watersheds_ll))) {
    gage_id <- watersheds_ll$gage_id[i]
    area_sqkm <- watersheds_ll$area_sqkm[i]
    
    if (verbose) {
      cat(sprintf("\nProcessing watershed %d/%d: %s (%.1f sq km)...\n",
                  i, nrow(watersheds_ll), gage_id, area_sqkm))
    }
    
    # Get bounding box
    bbox <- st_bbox(watersheds_ll[i, ])
    lat_range <- bbox["ymax"] - bbox["ymin"]
    lon_range <- bbox["xmax"] - bbox["xmin"]
    
    # Calculate number of spatial batches needed
    lat_batches <- ceiling(lat_range / batch_degrees)
    lon_batches <- ceiling(lon_range / batch_degrees)
    total_batches <- lat_batches * lon_batches
    
    if (verbose && total_batches > 1) {
      cat(sprintf("  Large watershed - will download in %d spatial batches (%d x %d)\n",
                  total_batches, lat_batches, lon_batches))
    }
    
    # Process each variable
    for (var_idx in seq_along(daymet_variables)) {
      var <- daymet_variables[var_idx]
      
      if (verbose) cat(sprintf("  Variable %d/%d: %s\n", var_idx, length(daymet_variables), var))
      
      # Collect data for all year chunks
      var_all_data <- list()
      
      # Process each year chunk
      for (chunk_idx in seq_along(year_chunks)) {
        chunk_years <- year_chunks[[chunk_idx]]
        chunk_start <- min(chunk_years)
        chunk_end <- max(chunk_years)
        
        if (verbose) {
          cat(sprintf("    Years %d-%d...\n", chunk_start, chunk_end))
        }
        
        # Collect data for all spatial batches in this chunk
        chunk_batch_data <- list()
        
        # Process spatial batches
        batch_num <- 0
        for (lat_batch in seq_len(lat_batches)) {
          for (lon_batch in seq_len(lon_batches)) {
            batch_num <- batch_num + 1
            
            # Calculate batch bounds
            lat_min <- bbox["ymin"] + (lat_batch - 1) * batch_degrees
            lat_max <- min(bbox["ymax"], bbox["ymin"] + lat_batch * batch_degrees)
            lon_min <- bbox["xmin"] + (lon_batch - 1) * batch_degrees
            lon_max <- min(bbox["xmax"], bbox["xmin"] + lon_batch * batch_degrees)
            
            # Create unique directory for this batch
            batch_dir <- file.path(temp_nc_dir, 
                                   sprintf("%s_%s_y%d-%d_b%d", 
                                           gage_id, var, chunk_start, chunk_end, batch_num))
            dir.create(batch_dir, showWarnings = FALSE)
            
            tryCatch({
              # Small delay to avoid overwhelming server
              if (batch_num > 1) Sys.sleep(1)
              
              # Download data for this year chunk
              download_daymet_ncss(
                location = c(lat_max, lon_min, lat_min, lon_max),
                start = chunk_start,
                end = chunk_end,
                param = var,
                frequency = "daily",
                path = batch_dir,
                silent = TRUE
              )
              
              # Find downloaded file
              nc_files <- list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE)
              
              if (length(nc_files) > 0 && file.exists(nc_files[1])) {
                # Read and process
                r <- terra::rast(nc_files[1])
                
                # Transform watershed to Daymet projection
                daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
                watershed_proj <- st_transform(watersheds_ll[i, ], daymet_crs)
                
                # Extract mean values for watershed
                extracted <- terra::extract(r, vect(watershed_proj), 
                                            fun = "mean", 
                                            weights = TRUE,
                                            touches = TRUE)
                
                # Get dates
                dates <- time(r)
                if (is.null(dates) || length(dates) == 0) {
                  # Generate dates for this chunk
                  n_days <- nlyr(r)
                  dates <- seq(as.Date(paste0(chunk_start, "-01-01")), 
                               by = "day", 
                               length.out = n_days)
                }
                
                # Create data table
                if (ncol(extracted) > 1) {
                  values <- as.numeric(extracted[1, -1])
                  batch_data <- data.table(
                    date = dates[1:length(values)],
                    value = values
                  )
                  chunk_batch_data[[length(chunk_batch_data) + 1]] <- batch_data
                }
              }
              
            }, error = function(e) {
              if (verbose) {
                cat(sprintf("      Error in spatial batch %d: %s\n", batch_num, e$message))
              }
            }, finally = {
              # Always clean up
              unlink(batch_dir, recursive = TRUE)
              if (exists("r")) rm(r)
              gc(verbose = FALSE)
            })
          }
        }
        
        # Combine spatial batches for this year chunk
        if (length(chunk_batch_data) > 0) {
          chunk_data <- rbindlist(chunk_batch_data)
          # Average overlapping dates from different spatial batches
          chunk_data <- chunk_data[, .(value = mean(value, na.rm = TRUE)), by = date]
          var_all_data[[length(var_all_data) + 1]] <- chunk_data
        }
      }
      
      # Combine all year chunks for this variable
      if (length(var_all_data) > 0) {
        var_data <- rbindlist(var_all_data)
        var_data[, gage_id := gage_id]
        all_results[[var]][[gage_id]] <- var_data
      }
    }
    
    # Periodic cleanup
    if (i %% 5 == 0) gc(verbose = FALSE)
  }
  
  # Save results by variable
  if (verbose) cat("\nSaving results to parquet files...\n")
  
  for (var in names(all_results)) {
    if (length(all_results[[var]]) > 0) {
      # Combine all watersheds for this variable
      combined_data <- rbindlist(all_results[[var]])
      
      # Pivot to wide format (dates as rows, gage_ids as columns)
      wide_data <- dcast(combined_data, date ~ gage_id, value.var = "value")
      
      # Sort by date
      setorder(wide_data, date)
      
      # Save
      output_file <- file.path(output_dir, "by_variable", sprintf("daymet_%s_daily.parquet", var))
      arrow::write_parquet(wide_data, output_file)
      
      if (verbose) {
        cat(sprintf("  Saved %s: %d dates x %d watersheds (%.1f MB)\n", 
                    var, nrow(wide_data), ncol(wide_data) - 1,
                    file.size(output_file) / 1024^2))
      }
    }
  }
  
  # Clean up
  unlink(temp_nc_dir, recursive = TRUE)
  
  if (verbose) cat("\nProcessing complete!\n")
  
  return(output_dir)
}

# Execute with year chunking
result <- aggregate_daymet_for_watersheds_chunked(
  watersheds_file = file.path(metadata_dir, "geospatial_derivedData/unified_watersheds_simplified.gpkg"),
  output_dir = file.path(metadata_dir, "daymet_watershed_timeseries"),
  max_area_sqkm = 10000,
  daymet_variables = c("tmin", "tmax", "prcp", "srad", "vp", "swe", "dayl"),
  start_year = 1980,
  end_year = 2024,
  years_per_chunk = 10,  # Process 10 years at a time to stay under 6GB
  batch_degrees = 0.25,
  verbose = TRUE
)




# Summary function
create_processing_summary <- function(output_dir, watersheds, variables, start_year, end_year) {
  
  # Get saved files
  var_files <- list.files(file.path(output_dir, "by_variable"), 
                          pattern = "daymet_.*_daily\\.parquet$", 
                          full.names = TRUE)
  
  # Create summary of each variable file
  file_summary <- data.frame()
  
  for (file in var_files) {
    var_name <- gsub(".*daymet_(.*)_daily\\.parquet", "\\1", basename(file))
    data <- arrow::read_parquet(file)
    
    file_info <- data.frame(
      variable = var_name,
      n_watersheds = ncol(data) - 1,
      n_dates = nrow(data),
      date_range = paste(min(data$date), "to", max(data$date)),
      file_size_mb = round(file.size(file) / 1024^2, 2),
      stringsAsFactors = FALSE
    )
    
    file_summary <- rbind(file_summary, file_info)
  }
  
  # Save summary
  write.csv(file_summary, file.path(output_dir, "processing_summary.csv"), row.names = FALSE)
  
  # Create detailed metadata
  metadata <- list(
    processing_date = Sys.time(),
    total_watersheds_input = nrow(watersheds),
    watersheds_processed = length(unique(unlist(lapply(var_files, function(f) {
      names(arrow::read_parquet(f))[-1]
    })))),
    variables = variables,
    period = paste(start_year, "to", end_year),
    file_summary = file_summary
  )
  
  saveRDS(metadata, file.path(output_dir, "processing_metadata.rds"))
  
  return(metadata)
}


# Function to inspect results
inspect_daymet_results <- function(output_dir, variable = "prcp") {
  
  file_path <- file.path(output_dir, "by_variable", sprintf("daymet_%s_daily.parquet", variable))
  
  if (!file.exists(file_path)) {
    cat("Available variables:\n")
    files <- list.files(file.path(output_dir, "by_variable"), pattern = "daymet_.*_daily\\.parquet$")
    vars <- gsub("daymet_(.*)_daily\\.parquet", "\\1", files)
    cat(paste(" -", vars, collapse = "\n"))
    return(NULL)
  }
  
  data <- arrow::read_parquet(file_path)
  
  cat(sprintf("\nData for %s:\n", variable))
  cat(sprintf("Date range: %s to %s (%d years)\n", 
              min(data$date), max(data$date),
              as.numeric(difftime(max(data$date), min(data$date), units = "days")) / 365.25))
  cat(sprintf("Number of watersheds: %d\n", ncol(data) - 1))
  cat(sprintf("Number of days: %d\n", nrow(data)))
  cat(sprintf("File size: %.1f MB\n", file.size(file_path) / 1024^2))
  
  # Check for missing data
  n_missing <- sum(is.na(data[, -1]))
  n_total <- (ncol(data) - 1) * nrow(data)
  cat(sprintf("Missing values: %d (%.2f%%)\n", n_missing, 100 * n_missing / n_total))
  
  return(data)
}

# Inspect results
prcp_data <- inspect_daymet_results(result, "prcp")

