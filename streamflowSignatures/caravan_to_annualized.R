# Set the main directory
main_dir = "C:/Users/18033/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"

# Load configuration and helper functions
source(file.path(main_dir, "config.R"))           # Centralized configuration
source(file.path(main_dir, "helperFunctions.R"), local=FALSE)


# Define Caravan base directory
caravan_main_dir = "J:/Downloads/Caravan-nc/usr/local/google/home/kratzert/Data/Caravan-Jan25-nc" 
caravan_main_dir = "D:/Caravan-nc/usr/local/google/home/kratzert/Data/Caravan-Jan25-nc" 
# (Adjust this path to your actual Caravan base directory)

# Define output file for a specific Caravan project
caravan_output_file_camels = file.path(main_dir, "summary_data_caravan_camels.csv") 
caravan_input_files_camels = file.path(caravan_main_dir, "timeseries/netcdf/camels") 

# User configuration 
start_date_analysis = as.Date("1979-09-01") 
end_date_analysis = as.Date("2025-06-01")   
min_num_years_analysis = 30   # minimum period of record (years) required for inclusion in analysis
min_Q_val_days = c(0.00001, 30) # minimum flow (mm) required before setting to NA; minimum number of non-NA days required to count a year 



process_caravan_to_annual(caravan_directory = caravan_input_files_camels, 
                            data_project = "camels",
                            min_num_years_data = 30,
                            start_date_filter = as.Date("1979-09-01"), 
                            end_date_filter = as.Date("2025-06-01"),
                            output_dir = "annualized_caravan_data"
                          )

