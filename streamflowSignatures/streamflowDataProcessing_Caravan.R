# Set the main directory
main_dir = "C:\\Users\\18033\\Documents\\GitHub\\SurfaceWaterProjections\\streamflowSignatures"#"C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"

source(file.path(main_dir, "helperFunctions.R"), local=FALSE)


# Define Caravan base directory
caravan_main_dir = "J:/Downloads/Caravan-nc/usr/local/google/home/kratzert/Data/Caravan-Jan25-nc" 
caravan_main_dir = "D:/Caravan-nc/usr/local/google/home/kratzert/Data/Caravan-Jan25-nc" 
# (Adjust this path to your actual Caravan base directory)

# Define output file for a specific Caravan project
caravan_output_file_camels = file.path(main_dir, "summary_data_caravan_hysets.csv") 

# User configuration 
start_date_analysis = as.Date("1970-01-01") 
end_date_analysis = as.Date("2024-12-31")   
min_num_years_analysis = 20   # minimum period of record (years) required for inclusion in analysis
min_Q_val_days = c(0.0001, 30) # minimum flow (mm) required before setting to NA; minimum number of non-NA days required to count a year 

# Example: Process the "camels" data_project from Caravan
cat("\n========== PROCESSING CARAVAN CAMELS DATA ==========\n")
summary_camels_caravan <- process_caravan_gages(
  data_project_arg = "hysets", #e.g., "camels", "hysets", etc.
  caravan_base_dir = caravan_main_dir,
  min_num_years = min_num_years_analysis,
  start_date = start_date_analysis,
  end_date = end_date_analysis,
  min_Q_value_and_days = min_Q_val_days,
  output_file = caravan_output_file_camels
)


cat("\n========== CARAVAN ANALYSIS COMPLETE ==========\n")
if (exists("summary_camels_caravan") && nrow(summary_camels_caravan) > 0) {
  cat("Completed processing Caravan 'camels' project. Final summary data has", nrow(summary_camels_caravan), "rows\n")
  cat("Results saved to:", caravan_output_file_camels, "\n")
}


