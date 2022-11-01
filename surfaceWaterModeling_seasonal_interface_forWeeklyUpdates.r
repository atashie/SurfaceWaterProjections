	# reading in the ncdfs to ensure timing / structure 
library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
library(data.table)		# for data.table
	
	# example input data
theseBasins = data.frame(basinSymbol = 	c('CLE',		'ISB',		'EXC',		'ORO',		'SHA',		'FOL',		'PNF',		'MIL',		'TRM'),
								Lon = 	c(-122.765, 	-118.479,	-120.264,	-121.480,	-122.417,	-121.155,	-119.318,	-119.6545,	-118.998),
								Lat = 	c(40.8225,		35.650,		37.591,		39.540,		40.720,		38.710,		36.845,		37.0425,	36.416),
						incldStorage = 	c(TRUE, 		TRUE,		TRUE,		TRUE,		TRUE,		TRUE,		TRUE,		FALSE,		FALSE),
						maxStorage = 	c(2447650,		568000,		1024600,	3537577,	4552000,	977000,		1000000,	520500,		185600))

yesterdaysDate = '2022-10-26'	# historic data is released every day for the day prior
forecastDate = '27OCT022'		# 
waterYearStart = as.Date('2022-10-01')

	# these file locations remain the same for all watersheds of a region
mainFolder = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
seas5DataNCDF = paste0(mainFolder, 'NuveenNorCal-testing-seas5.nc') 				# SoCal
#cfsDataNCDF =   'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\NuveenNorCal-testing-cfs.nc'			 	# SoCal
era5DataNCDF =   paste0(mainFolder, 'NuveenNorCal-testing-recent-era.nc')			# SoCal
basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'

	# correct dates must be manually selected for now
#cfsStartDate = as.Date('2022-02-28') #  + ncvar_get(ncin_cfs, 'time')/24 for calculating actual dates
era5StartDate =  as.Date('2020-08-01') # + ncvar_get(ncin_era5, 'time') for calculating actual dates 
#era5RecentStartDate =  as.Date('2001-07-01') # + ncvar_get(ncin_recentEra5, 'time') for calculating actual dates 
seas5StartDate = as.Date('2022-10-01') # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates

		### this section is purely for inspecting data before running models
#########################################################################
#ncin_cfs = nc_open(cfsDataNCDF)		;	ncin_cfs
ncin_era5 = nc_open(era5DataNCDF)	;	ncin_era5
#ncin_recentEra5 = nc_open(era5DataRecentNCDF)	;	ncin_recentEra5
ncin_seas5 = nc_open(seas5DataNCDF)	;	ncin_seas5

#ncatt_get(ncin_cfs, 'time','units')$value	# says hours but seems to be in days
ncatt_get(ncin_era5, 'time','units')$value
#ncatt_get(ncin_recentEra5, 'time','units')$value
ncatt_get(ncin_seas5, 'valid_time','units')$value
	#########################################################################

	# initializing dataframe for storing outputs
allForecasts = data.frame(Date=NA, Reservoir = NA, Units = NA, Lat = NA, Lon = NA,
		Pred_Q05 = NA, Pred_Q25 = NA, Pred_Q50 = NA, Pred_Q75 = NA, Pred_Q95 = NA,
		Clim_Q05 = NA, Clim_Q25 = NA, Clim_Q50 = NA, Clim_Q75 = NA, Clim_Q95 = NA,
		WaterYear_Curr = NA, WaterYear_1YrAgo = NA, WaterYear_2YrAgo = NA, WaterYear_3YrAgo = NA,
		MaxStorage = NA, StorageOrStreamflow = NA)


for(thisBasin in 1:nrow(theseBasins))	{

	basinSymbol = theseBasins$basinSymbol[thisBasin]
	basinName = basinSymbol # paste0(basinSymbol, '_atOutlet')

	gageLonLat = c(theseBasins$Lon[thisBasin], theseBasins$Lat[thisBasin])
	infOrFnf = 76 #8 for fnf, 76 for inflow

	historicStreamflowFileLoc =   paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=", infOrFnf, "&dur_code=D&Start=1900-01-01&End=", yesterdaysDate)
	historicReservoirFileLoc = paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=15&dur_code=D&Start=1900-01-01&End=", yesterdaysDate)

		# defining pathways to basin-specific files
	dataOut_location = paste0('J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\', basinName, '\\')
		
	
	#########################################################################################################
	####	This section is for running and making projections with a previously calibrated model
	####
	#########################################################################################################
						
	#########################################################################################################
		# step 8
		# import and convert projection climate data


	climateInputConversion_f(
		basinName = basinName,
		climateDataNCDF = seas5DataNCDF,	####!!!!! change btw cfs and seas5
		tempConversionFactor = NA,
		pptConversionFactor = NA,
		avgTempGiven = FALSE, 
		startDate = seas5StartDate, 	# when does the clock of the netcdf start?
		timeToDaysConversion = 1,	# convert time increments to days if necessary
		dataOut_location = dataOut_location,
		optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
		variableOrderOption = 'seas5', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
		precipName = 'tp')	# other options include: tp, tp_sum	



	climateInputConversion_f(
		basinName = basinName,
		climateDataNCDF = era5DataNCDF,	####!!!!! change btw cfs and seas5
		tempConversionFactor = NA,
		pptConversionFactor = NA,
		avgTempGiven = FALSE, 
		startDate = era5StartDate, 	# when does the clock of the netcdf start?
		timeToDaysConversion = 1,	# convert time increments to days if necessary
		dataOut_location = dataOut_location,
		optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
		variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
		precipName = 'tp_sum')	# other options include: tp, tp_sum	
		


		# step 9a 
		# run the model with forecasting data
		## Running the Model for Seasonal Forecasts 
	seasonalStreamflowForecast_f(
		basinName = basinName,
		historicStreamflowFileLoc = historicStreamflowFileLoc,
		dataOut_location = dataOut_location,
		dataSource = 1,							# 1 for cal.gov,
		waterYearStart = waterYearStart,
		forecastDate = forecastDate,
		gageLonLat = gageLonLat,
		biasCorrection = TRUE,
		uploadToGCS = TRUE)

	dataOut_fileLoc = paste0(dataOut_location, 'strmflwForecastsFor_', forecastDate)
	thisForecast = fread(paste0(dataOut_fileLoc, "forecastStrmCumSum_", basinName, '_', forecastDate, ".csv"))
	thisForecast$MaxStorage = NA
	thisForecast$StorageOrStreamflow = 'Streamflow'
	allForecasts = rbind(allForecasts, thisForecast)

		# step 9b 
		# run the model with forecasting data
		## Running the Model for Seasonal Forecasts 
	seasonalStorageForecast_f(
		basinName = basinName,
		historicStreamflowFileLoc = historicStreamflowFileLoc,
		historicReservoirFileLoc = historicReservoirFileLoc,
		dataOut_location = dataOut_location,
		dataSource = 1,							# 1 for cal.gov,
		waterYearStart = waterYearStart,
		forecastDate = forecastDate,
		gageLonLat = gageLonLat,
		biasCorrection = TRUE,
		uploadToGCS = TRUE,
		incldStorage = TRUE)
		
	dataOut_fileLoc = paste0(dataOut_location, 'storageForecastsFor_', forecastDate)
	thisForecast = fread(paste0(dataOut_fileLoc, "forecastStorage_", basinName, '_', forecastDate, ".csv"))
	thisForecast$MaxStorage = theseBasins$maxStorage[thisBasin]
	thisForecast$StorageOrStreamflow = 'Storage'
	allForecasts = rbind(allForecasts, thisForecast)
}
#allForecasts$Date = seq(as.Date('2022-10-01'), as.Date('2023-09-30'), by = 1)
fwrite(allForecasts, paste0(mainFolder, "allForecasts_", forecastDate, ".csv"))
	



