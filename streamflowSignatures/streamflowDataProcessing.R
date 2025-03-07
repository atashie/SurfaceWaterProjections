################################################################################################
# Streamflow Signatures Analysis
# 
# This script processes streamflow data from USGS and Canadian gages to calculate various
# hydrological signatures and trends. It analyzes flow volumes, flow duration curves,
# flashiness, and flow timing.
#
# Last updated: 4-MAR-2025
################################################################################################

# Load helper functions
source(file.path(main_dir, "helperFunctions.R"), local=FALSE)

################################################################################################
# USER CONFIGURATION - MODIFY THESE SETTINGS AS NEEDED
################################################################################################

# Set the main directory where data and results will be stored
main_dir = "C:/Users/arik/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"

# Analysis period
start_date = as.Date("1973-01-01")  # Modern data period
end_date = as.Date("2024-12-31")    # End of analysis period

# Data quality thresholds
min_num_years = 20    # Minimum number of years required for analysis
min_nona_days = 250   # Minimum number of non-NA days per year
min_Q_value_and_days = c(0.0001, 30)  # Min flow value (mm) and days above this value

# Output file for results
output_file = file.path(main_dir, "summary_data.csv")

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
  conus_gages_raw = fread(file.path(main_dir, "metadata", "conterm_bas_classif.txt"),
                           colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )[CLASS=="Ref"]
  conus_basinid = fread(file.path(main_dir, "metadata", "conterm_basinid.txt"), 
                         colClasses = c("STAID" = "character"))
  conus_gages = merge(conus_gages_raw, conus_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(conus_gages), "CONUS reference gages\n")
}, error = function(e) {
  stop("Error loading CONUS gage data: ", e$message)
})

# Load USGS Alaska reference gages
tryCatch({
  AK_gages_all = fread(file.path(main_dir, "metadata", "AKHIPR_bas_classif.txt"),
                        colClasses = c("STAID" = "character", "AGGECOREGION" = "character")
  )[AGGECOREGION == 'Alaska' & CLASS == 'Ref']
  AK_basinid = fread(file.path(main_dir, "metadata", "AKHIPR_basinid.txt"), 
                      colClasses = c("STAID" = "character"))
  AK_gages = merge(AK_gages_all, AK_basinid, by="STAID", all.x=TRUE)
  cat("Loaded", nrow(AK_gages), "Alaska reference gages\n")
}, error = function(e) {
  stop("Error loading Alaska gage data: ", e$message)
})

# Load Canadian reference gages
tryCatch({
  canadian_gages_goodData = fread(file.path(main_dir, "metadata", "Canadian_gages_goodones.csv"))
  regulation_info = as.data.table(hy_stn_regulation(canadian_gages_goodData$STATION_NUMBER))
  canadian_gages = merge(canadian_gages_goodData, regulation_info, 
                          by = "STATION_NUMBER", all.x = TRUE)[REGULATED != TRUE]
  cat("Loaded", nrow(canadian_gages), "Canadian reference gages\n")
}, error = function(e) {
  stop("Error loading Canadian gage data: ", e$message, 
       "\nNote: You need to install and configure the tidyhydat package for Canadian data.")
})

# Load watershed boundaries (HydroBASINS)
tryCatch({
  basinAt_NorAm_polys = st_read(file.path(main_dir, "basinAt_NorAm_polys.gpkg"))
  basinAt_NorAm_strip = basinAt_NorAm_polys
  st_geometry(basinAt_NorAm_strip) = NULL
  HB_dt = data.table(basinAt_NorAm_strip)
  cat("Loaded HydroBASINS watershed boundaries\n")
}, error = function(e) {
  stop("Error loading HydroBASINS data: ", e$message, 
       "\nPlease ensure basinAt_NorAm_polys.gpkg is in the main directory.")
})

# Load or initialize upstream hydrobasins cache
upstream_hydrobasins_file = file.path(main_dir, "upstream_hydrobasins.RData")
if (file.exists(upstream_hydrobasins_file)) {
  upstream_hydrobasins = readRDS(upstream_hydrobasins_file)
  cat("Loaded cached upstream basin relationships\n")
} else {
  upstream_hydrobasins = list()
  cat("Initialized new upstream basin cache\n")
}

################################################################################################
# PROCESS STREAMFLOW DATA
################################################################################################

# Create output directory if it doesn't exist
output_dir = dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Process USGS Alaska gages
cat("\n========== PROCESSING ALASKA GAGES ==========\n")
summary_output = process_gages(
  gages_df = AK_gages,
  gage_type = "USGS",
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  basinAt_NorAm_polys = basinAt_NorAm_polys,
  HB_dt = HB_dt,
  upstream_hydrobasins = upstream_hydrobasins,
  output_file = output_file
)

# Process USGS CONUS gages
cat("\n========== PROCESSING CONUS GAGES ==========\n")
summary_output = process_gages(
  gages_df = conus_gages,
  gage_type = "USGS",
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  basinAt_NorAm_polys = basinAt_NorAm_polys,
  HB_dt = HB_dt,
  upstream_hydrobasins = upstream_hydrobasins,
  output_file = output_file
)

# Process Canadian gages
cat("\n========== PROCESSING CANADIAN GAGES ==========\n")
summary_output = process_gages(
  gages_df = canadian_gages,
  gage_type = "Canada",  
  min_num_years = min_num_years,
  start_date = start_date,
  end_date = end_date,
  min_Q_value_and_days = min_Q_value_and_days,
  basinAt_NorAm_polys = basinAt_NorAm_polys,
  HB_dt = HB_dt,
  upstream_hydrobasins = upstream_hydrobasins,
  output_file = output_file
)

# Final summary
cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Completed processing all gages. Final summary data has", nrow(summary_output), "rows\n")
cat("Results saved to:", output_file, "\n")
