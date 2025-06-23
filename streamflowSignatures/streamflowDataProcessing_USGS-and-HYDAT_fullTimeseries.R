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
