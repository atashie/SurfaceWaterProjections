	# example input data
basinSymbol = 'EXC'
basinName = basinSymbol # paste0(basinSymbol, '_atOutlet')
gageLonLat =       c(-120.264, 37.591)
infOrFnf = 76 #8 for fnf, 76 for inflow
	# list of basins by symbol and lon / lat
	# webpage to search for Cali reservoirs: https://cdec.water.ca.gov/dynamicapp/wsSensorData
				# for ENG: c(-121.270278, 39.240278)
				# for ISB: c(-118.479, 35.650) 
				# for PNF: c(-119.318, 36.845)
				# for TRM: c(-118.998, 36.416) 
				# for EXC: c(-120.264, 37.591) 
				# for ORO: c(-121.480, 39.540) 
				# for CLE: c(-122.765, 40.8225)
				# for MIL: c(-119.6545, 37.0425)
				# for SHA: c(-122.417, 40.720)		
				# for FOL: c(-121.155, 38.710)
					# these reservoirs are problematic, e.g. are immediately downstream of another reservoir
				# for SJF: c(-119.697,37.004) 
				# for YRS: c(-121.292500, 39.223611) 
				# for LEW: c(-122.800, 40.732)
					# this loc only has fnf, no inflow 
				# for BND: c(-122.185556, 40.288611) 

yesterdaysDate = '2022-09-13'	# historic data is released every day for the day prior
historicStreamflowFileLoc =   paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=", infOrFnf, "&dur_code=D&Start=1900-01-01&End=", yesterdaysDate)
historicReservoirFileLoc = paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=15&dur_code=D&Start=1900-01-01&End=", yesterdaysDate)

	# defining pathways to basin-specific files
dataOut_location = paste0('J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\', basinName, '\\')
forecastDate = '7SEP2022'
waterYearStart = as.Date('2022-10-01')

	# these file locations remain the same for all watersheds
seas5DataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-recent-seas5.nc'
cfsDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-cfs.nc'
	# no longer separating historical from recent era5... may revisit
#era5DataHistoricalNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-hist-era5.nc'
era5DataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-era.nc'
#seas5MultiDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\testing-multiple-forecasts-seas5-wy2019.nc'
basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'

	
	
	# reading in the ncdfs to ensure timing / structure 
library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
library(data.table)		# for data.table

	### this section is purely for inspecting data before running models
#########################################################################
ncin_cfs = nc_open(cfsDataNCDF)		;	ncin_cfs
ncin_era5 = nc_open(era5DataNCDF)	;	ncin_era5
#ncin_recentEra5 = nc_open(era5DataRecentNCDF)	;	ncin_recentEra5
ncin_seas5 = nc_open(seas5DataNCDF)	;	ncin_seas5

ncatt_get(ncin_cfs, 'time','units')$value	# says hours but seems to be in days
ncatt_get(ncin_era5, 'time','units')$value
#ncatt_get(ncin_recentEra5, 'time','units')$value
ncatt_get(ncin_seas5, 'valid_time','units')$value
#########################################################################

	# correct dates must be manually selected for now
cfsStartDate = as.Date('2022-02-28') #  + ncvar_get(ncin_cfs, 'time')/24 for calculating actual dates
era5StartDate =  as.Date('2001-07-01') # + ncvar_get(ncin_era5, 'time') for calculating actual dates 
#era5RecentStartDate =  as.Date('2001-07-01') # + ncvar_get(ncin_recentEra5, 'time') for calculating actual dates 
	# seas5 is incorrectly showing the second of the month, but should be the first
seas5StartDate = as.Date('2022-08-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates
seas5MultiStartDate = as.Date('1993-01-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates


#########################################################################################################
	# step 1
	# delineate a new basin
basinDelineation_f(
	gageLonLat = gageLonLat ,
	basinATLAS_locationAndFile = basinATLAS_locationAndFile,
	dataOut_location = dataOut_location,
	basinName = basinName)


#########################################################################################################
	# step 2
	# import and convert historical climate data
climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = era5DataNCDF,
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	startDate = era5StartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	
	
#########################################################################################################
	# step 3 
	# calibrate model
modelCalibration_f(
	historicStreamflowFileLoc = historicStreamflowFileLoc,	#'https://someplace.gov',
	dataOut_location = dataOut_location,					#'save_file_location',
	dataSource = 1,											# 1 for FNF from cal.gov, 
	numberOfRuns = 200000,
	targetMetric = 1, 										# 1 = KGE; other options not included yey
	targetMetricValue = 0.81,								# threshold for considering a value good
	minGoodRuns = 200,										# number of 'good' calibrations before the routine stops
	sfcf = c(runif(50000, .2, 1), runif(50000, 1, 3)),			#snowfall correction factor [-]
	tr   = runif(100000, -6, 5),								#solid and liquid precipitation threshold temperature [C]
	tt   = runif(10000, -5, 6),								#melt temperature [C]
	fm   = c(runif(50000, .2, 1.5), (runif(50000, 1.5, 8))),	#snowmelt factor [mm/C]
	fi   = c(runif(50000, .2, 1.5), (runif(50000, 1.5, 10))),	#icemelt factor [mm/C]
	fic  = runif(100000, 2, 10),								#debris-covered icemelt factor [mm/C]
	fc   = c(runif(50000, 25, 150), (runif(50000, 150, 1200))),	#field capacity
	lp   = runif(100000, .2, 0.999),								#parameter to actual ET
	beta_soils = runif(100000, 1, 3),							#beta - exponential value for nonlinear relations between soil storage and runoff
	k0   = c(runif(50000, .05, .5), (runif(50000, .5, 0.999))),		#top bucket drainage
	k1   = c(runif(50000, .005, .09), (runif(50000, .09, .5))),	#middle bucket drainage
	k2   = c(runif(50000, .0001, .01), (runif(50000, .01, .1))),#bottom bucket drainage	
	uz1  = c(runif(50000, .22, 10), (runif(50000, 10, 40))),	#max flux rate from STZ to SUZ in mm/d
	perc = c(runif(50000, .1, .5), (runif(50000, .5, 20)))	# max flux rate from SUZ to SLZ in mm/d
)	



#########################################################################################################
####	This section is for validating model performance 
####
#########################################################################################################

	# step 4
	# validate historical data and generate plots
validationHistoricalOutput = validationAndPlotGeneration_f(
	basinName = basinName,
#	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,


	# step 5
	# import and convert multi projection climate data for validationAndPlotGeneration_f
allWYs = 2002:2021
for(thisWY in allWYs)	{
	seas5MultiDataNCDF = paste0('J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\testing-multiple-forecasts-seas5_wy', thisWY, '.nc')
	climateInputConversion_f(
		basinName = basinName,
		climateDataNCDF = seas5MultiDataNCDF,
		tempConversionFactor = NA,
		pptConversionFactor = NA,
		avgTempGiven = FALSE, 
		startDate = seas5MultiStartDate, 	# when does the clock of the netcdf start?
		timeToDaysConversion = 1,	# convert time increments to days if necessary
		dataOut_location = dataOut_location,
		optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
		variableOrderOption = 'seas5Multi', 
		precipName = 'tp_sum',
		limitedModels = 25)	# other options include: tp, tp_sum	
}


	# step 6
	# validate streamflow projection data and generate plots
projectionValidationAndPlotGeneration_f(
	basinName = basinName,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1,							# 1 for cal.gov,
	biasCorrection = TRUE)							


#########################################################################################################
####	This section is for linking streamflow models to reservoirs
####
#########################################################################################################
	# step 7
	# ML for predicting reservoir outflows
#reservoirOutflowCalVal_f(
#	dataOut_location = dataOut_location,
#	historicReservoirFileLoc = historicReservoirFileLoc,
#	historicStreamflowFileLoc = historicStreamflowFileLoc,
#	basinName = basinName,
#	dataSource = 1, 			# 1 for cal.gov,
#	nTreeModel=500,
##	modelMetric = 'Rsquared',
#	modelMethod = 'rf',
#	metric="Rsquared",
#	numFolds = 5,
#	numRepeats = 30)

	# step 8
	# validation of combined model projections of total storage
projectedStorageValidationAndPlotGeneration_f(
	basinName = basinName,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	historicReservoirFileLoc = historicReservoirFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1,
	biasCorrection = TRUE)							# 1 for FNF from cal.gov,
	
	
#########################################################################################################
####	This section is for running and making projections with a previously calibrated model
####
#########################################################################################################
					
#########################################################################################################
	# step 9
	# import and convert projection climate data
climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
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
	

#cfsClimateDataframe = climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
#	basinName = basinName,
#	climateDataNCDF = cfsDataNCDF,	####!!!!! change btw cfs and seas5
#	tempConversionFactor = NA,
#	pptConversionFactor = NA,
#	avgTempGiven = FALSE, 
#	multipleModels = FALSE,	# are there multiple models that need to be stored?
#	startDate = cfsStartDate, 	# when does the clock of the netcdf start?
#	timeToDaysConversion = 24,	# convert time increments to days if necessary
#	dataOut_location = dataOut_location,
#	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
#	variableOrderOption = 'cfs', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
#	precipName = 'tp')	# other options include: tp, tp_sum	

#saveRDS(cfsClimateDataframe, paste0(climateInputsFileLoc, 'CFS.RData'))

climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
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
	precipName = 'tp_sum')	# other options include: tp, tp_sum	


	# step 10 
	# run the model with forecasting data
	## Running the Model for Seasonal Forecasts 
seasonalStreamflowForecast_f(
	basinName = basinName,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1,							# 1 for FNF from cal.gov,
	waterYearStart = waterYearStart,
	forecastDate = forecastDate,
	gageLonLat = gageLonLat)



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# works for EXC, SHA, CLE, TRM, PNF, ISB, ORO

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'EXC_atOutlet',
	basinSymbol = 'EXC',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'TRM_atOutlet',
	basinSymbol = 'TRM',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'SHA_atOutlet',
	basinSymbol = 'SHA',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'CLE_atOutlet',
	basinSymbol = 'CLE',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'PNF_atOutlet',
	basinSymbol = 'PNF',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ISB_atOutlet',
	basinSymbol = 'ISB',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ORO_atOutlet',
	basinSymbol = 'ORO',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


# works for FOL, MIL, LEW, ENG
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'FOL_atOutlet',
	basinSymbol = 'FOL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'MIL_atOutlet',
	basinSymbol = 'MIL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'LEW_atOutlet',
	basinSymbol = 'LEW',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ENG_atOutlet',
	basinSymbol = 'ENG',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


	
	
# works for FOL, MIL, LEW, ENG
#FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor =  as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

#FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
#FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = stor$VALUE
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = stor[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

plot(allDat$stor, allDat$outflow)
plot(allDat$stor, allDat$inflow)
#plot(allDat$inflow, allDat$fnf)
plot(allDat$inflow, allDat$outflow)
plot(allDat$Date, allDat$outflow, type='l', col='red', ylim=c(0,10000))
lines(allDat$Date, allDat$inflow, type='l', col='blue')


# BND has no dam, YRS has no dam, SJF has no dam




############ rf modeling of storage
library(data.table)
library(lubridate)
library(CAST)
library(caret)
library(sf)
library(zoo)

basinSymbol = 'TRM'
basinName = paste0(basinSymbol, '_atOutlet')
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
#pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
#pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
infOrFNF = 'fnf' #'inflow'
if(infOrFNF == 'fnf')	{
	inflow = as.data.table(read.csv(paste0(
		"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}	else	{
	if(infOrFNF == 'inflow')	{
	inflow = as.data.table(read.csv(paste0(
				"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}}
stor = read.csv(paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = inflow[stor[c('Date','stor')], on = 'Date']
allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$resLoss = c(diff(allDat$stor) - allDat$AcFtPDay_inflow[-1], NA) * -1


#######################################################################################################################################
### for monthly estimations


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, resLoss = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$resLoss, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$resLoss = as.numeric(monthDat$resLoss)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(6, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_inflow, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,3,4)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,3,4)], max(lagLength))

resSubset = monthDat[, c('ymonth','resLoss','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_inflow, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$resLoss,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=2000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$resLoss))
plot(resValid$resLoss, predValid)



# old versions


basinName = 'EXC_atOutlet'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")

	
#testing works for EXC, SHA, CLE, TRM, PNF, ISB, EXC
FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05")
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9

plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1])
lines(diff(allDat$stor))
plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1], diff(allDat$stor))










### for daily estimations
	# removing nas; need a better method moving forward
while(any(is.na(allDat$AcFtPDay_fnf)))	{
	allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))] = allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))+1]
}
while(any(is.na(allDat$AcFtPDay_inflow)))	{
	allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))] = allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))+1]
}
while(any(is.na(allDat$AcFtPDay_outflow)))	{
	allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))] = allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))+1]
}
while(any(is.na(allDat$stor)))	{
	allDat$stor[which(is.na(allDat$stor))] = allDat$stor[which(is.na(allDat$stor))+1]
}

recLength = nrow(allDat)
lagLength = seq(1, 365*2, 7)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(allDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(allDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,3,6,10)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,3,6,10)], max(lagLength))

resSubset = allDat[, c('Date','AcFtPDay_outflow')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(allDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(allDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resSubset$doy = yday(resSubset$Date)
resSubset$month = month(resSubset$Date)

resComplete = resSubset[complete.cases(resSubset)]
resTest = subset(resComplete, Date < as.Date("2016-12-31"))

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars

resValid = subset(resComplete, Date >= as.Date("2016-12-31"))
predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)


#######################################################################################################################################
### for monthly estimations

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, AcFtPDay_outflow = NA, AcFtPDay_fnf = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_outflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_fnf, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$AcFtPDay_outflow = as.numeric(monthDat$AcFtPDay_outflow)
monthDat$AcFtPDay_fnf = as.numeric(monthDat$AcFtPDay_fnf)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(1, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,4,6)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,4,6)], max(lagLength))

resSubset = monthDat[, c('ymonth','AcFtPDay_outflow','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)

=======

	# example input data
basinSymbol = 'PNF'
basinName = paste0(basinSymbol, '_atOutlet')
infOrFNF = 'fnf' # fnf or inflow	
gageLonLat = 	 c(-119.318, 36.845)
				# for ENG: c(-121.270278, 39.240278)
				# for ISB: c(-118.479, 35.650) 
				# for PNF: c(-119.318, 36.845)
				# for TRM: c(-118.998, 36.416) 
				# for EXC: c(-120.264, 37.591) 
				# for ORO: c(-121.480, 39.540) 
				# for CLE: c(-122.765, 40.8225)
				# for MIL: c(-119.6545, 37.0425)
				# for BND: c(-122.185556, 40.288611) 
				# for SHA: c(-122.417, 40.720)
				# for FOL: c(-121.155, 38.710)

				# for SJF: c(-119.697,37.004) 
				# for YRS: c(-121.292500, 39.223611) 
				# for LEW: c(-122.800, 40.732)
historicStreamflowFileLoc = "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=PNF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ENG: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ORO: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for EXC: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for TRM: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=TRM&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for PNF: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=PNF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ISB: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ISB&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for BND: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=BND&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for SHA: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SHA&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for CLE: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=CLE&SensorNums=8&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for MIL: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=MIL&SensorNums=76&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for FOL: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=FOL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
		# baddies
	# for YRS: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=YRS&SensorNums=8&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for SJF: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SJF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for LEW: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=LEW&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
		# to search for Cali reservoirs: https://cdec.water.ca.gov/dynamicapp/wsSensorData

basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
seas5DataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-seas5.nc'
cfsDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-cfs.nc'
era5DataHistoricalNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-hist-era5.nc'
era5DataRecentNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-recent-era5.nc'
seas5MultiDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\testing-multiple-forecasts-seas5-wy2019.nc'
climateInputsFileLoc = paste0(dataOut_location, basinName, '_processedClimateData')
forecastDate = '28JUL2022'
	
	
	# reading in the ncdfs to ensure timing / structure 
library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
library(data.table)		# for data.table

	### this section is purely for inspecting data before running models
#########################################################################
ncin_cfs = nc_open(cfsDataNCDF)		;	ncin_cfs
ncin_era5 = nc_open(era5DataHistoricalNCDF)	;	ncin_era5
ncin_recentEra5 = nc_open(era5DataRecentNCDF)	;	ncin_recentEra5
ncin_seas5 = nc_open(seas5DataNCDF)	;	ncin_seas5

ncatt_get(ncin_cfs, 'time','units')$value	# says hours but seems to be in days
ncatt_get(ncin_era5, 'time','units')$value
ncatt_get(ncin_recentEra5, 'time','units')$value
ncatt_get(ncin_seas5, 'valid_time','units')$value
#########################################################################

	# correct dates must be manually selected for now
cfsStartDate = as.Date('2022-02-28') #  + ncvar_get(ncin_cfs, 'time')/24 for calculating actual dates
era5StartDate =  as.Date('2001-06-01') # + ncvar_get(ncin_era5, 'time') for calculating actual dates 
era5RecentStartDate =  as.Date('2020-06-01') # + ncvar_get(ncin_recentEra5, 'time') for calculating actual dates 
	# seas5 is incorrectly showing the second of the month, but should be the first
seas5StartDate = as.Date('2022-06-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates
seas5MultiStartDate = as.Date('1993-01-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates


#########################################################################################################
	# step 1
	# delineate a new basin
basinDelineation_f(
	gageLonLat = gageLonLat ,
	basinATLAS_locationAndFile = basinATLAS_locationAndFile,
	dataOut_location = dataOut_location,
	basinName = basinName)


#########################################################################################################
	# step 2
	# import and convert historical climate data
climateDataframe =	climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = era5DataHistoricalNCDF,
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = era5StartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	
	
saveRDS(climateDataframe, paste0(climateInputsFileLoc, 'Historical.RData'))

#########################################################################################################
	# step 3 
	# calibrate model
modelCalibration_f(
	climateInputsFileLoc = climateInputsFileLoc,			#'file_location_and_name.RData',
	historicStreamflowFileLoc = historicStreamflowFileLoc,	#'https://someplace.gov',
	pathToWatershedsGPKG = pathToWatershedsGPKG,			#'file_location_and_name.gpkg',
	dataOut_location = dataOut_location,					#'save_file_location',
	dataSource = 1,											# 1 for FNF from cal.gov, 
	numberOfRuns = 200000,
	targetMetric = 1, 										# 1 = KGE, 2 = NSE, 3 = MAE, 4 = RMSE, 5 = bias
	targetMetricValue = 0.81,								# threshold for considering a value good
	minGoodRuns = 100,										# number of 'good' calibrations before the routine stops
	sfcf = c(runif(5000, .2, 1), runif(5000, 1, 3)),			#snowfall correction factor [-]
	tr   = runif(10000, -6, 5),								#solid and liquid precipitation threshold temperature [C]
	tt   = runif(10000, -5, 6),								#melt temperature [C]
	fm   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 8))),	#snowmelt factor [mm/C]
	fi   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 10))),	#icemelt factor [mm/C]
	fic  = runif(10000, 2, 10),								#debris-covered icemelt factor [mm/C]
	fc   = c(runif(5000, 25, 150), (runif(5000, 150, 1200))),	#field capacity
	lp   = runif(10000, .2, 0.999),								#parameter to actual ET
	beta_soils = runif(10000, 1, 3),							#beta - exponential value for nonlinear relations between soil storage and runoff
	k0   = c(runif(5000, .05, .5), (runif(5000, .5, 0.999))),		#top bucket drainage
	k1   = c(runif(5000, .005, .09), (runif(5000, .09, .5))),	#middle bucket drainage
	k2   = c(runif(5000, .0001, .01), (runif(5000, .01, .1))),#bottom bucket drainage	
	uz1  = c(runif(5000, .22, 10), (runif(5000, 10, 40))),	#max flux rate from STZ to SUZ in mm/d
	perc = c(runif(5000, .1, .5), (runif(5000, .5, 20)))
)	




#########################################################################################################
####	This section is for running and making projections with a previously calibrated model
####
#########################################################################################################
					
#########################################################################################################
	# step 4
	# import and convert projection climate data
recentEra5ClimateDataframe = climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = era5DataRecentNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = era5RecentStartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	
	
saveRDS(recentEra5ClimateDataframe, paste0(climateInputsFileLoc, 'Recent.RData'))

#cfsClimateDataframe = climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
#	basinName = basinName,
#	climateDataNCDF = cfsDataNCDF,	####!!!!! change btw cfs and seas5
#	tempConversionFactor = NA,
#	pptConversionFactor = NA,
#	avgTempGiven = FALSE, 
#	multipleModels = FALSE,	# are there multiple models that need to be stored?
#	startDate = cfsStartDate, 	# when does the clock of the netcdf start?
#	timeToDaysConversion = 24,	# convert time increments to days if necessary
#	dataOut_location = dataOut_location,
#	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
#	variableOrderOption = 'cfs', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
#	precipName = 'tp')	# other options include: tp, tp_sum	

#saveRDS(cfsClimateDataframe, paste0(climateInputsFileLoc, 'CFS.RData'))

seas5ClimateDataframe = climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = seas5DataNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = seas5StartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'seas5', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	

saveRDS(seas5ClimateDataframe, paste0(climateInputsFileLoc, 'SEAS5.RData'))


	# step 5 
	# run the model with forecasting data
	## Running the Model for Seasonal Forecasts 
allForecastsOutput = seasonalForecast_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,

saveRDS(allForecastsOutput, paste0(climateInputsFileLoc, 'ForecastsFor_', forecastDate, '.RData'))





#########################################################################################################
####	This section is for validating model performance 
####
#########################################################################################################

	# step 6
	# validate historical data and generate plots
validationHistoricalOutput = validationAndPlotGeneration_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,


	# step 7
	# import and convert multi projection climate data for validationAndPlotGeneration_f
climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = seas5MultiDataNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = seas5MultiStartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'seas5Multi', 
	precipName = 'tp_sum')	# other options include: tp, tp_sum	


	# step 8
	# validate streamflow projection data and generate plots
validationProjectionOutput = projectionValidationAndPlotGeneration_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,



#########################################################################################################
####	This section is for linking streamflow models to reservoirs
####
#########################################################################################################
	# step 9
	# ML for predicting reservoir outflows
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = basinName,
	basinSymbol = basinSymbol,
	infOrFNF = infOrFNF,
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=500,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 5,
	numRepeats = 30)

	# step 10
	# validation of combined model projections of total storage
validationProjectedStorageOutput = storageProjectionValidationAndPlotGeneration_f(
	basinName = basinName,
	infOrFNF = infOrFNF,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# works for EXC, SHA, CLE, TRM, PNF, ISB, ORO

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'EXC_atOutlet',
	basinSymbol = 'EXC',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'TRM_atOutlet',
	basinSymbol = 'TRM',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'SHA_atOutlet',
	basinSymbol = 'SHA',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'CLE_atOutlet',
	basinSymbol = 'CLE',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'PNF_atOutlet',
	basinSymbol = 'PNF',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ISB_atOutlet',
	basinSymbol = 'ISB',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ORO_atOutlet',
	basinSymbol = 'ORO',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


# works for FOL, MIL, LEW, ENG
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'FOL_atOutlet',
	basinSymbol = 'FOL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'MIL_atOutlet',
	basinSymbol = 'MIL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'LEW_atOutlet',
	basinSymbol = 'LEW',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ENG_atOutlet',
	basinSymbol = 'ENG',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


	
	
# works for FOL, MIL, LEW, ENG
#FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor =  as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

#FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
#FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = stor$VALUE
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = stor[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

plot(allDat$stor, allDat$outflow)
plot(allDat$stor, allDat$inflow)
#plot(allDat$inflow, allDat$fnf)
plot(allDat$inflow, allDat$outflow)
plot(allDat$Date, allDat$outflow, type='l', col='red', ylim=c(0,10000))
lines(allDat$Date, allDat$inflow, type='l', col='blue')


# BND has no dam, YRS has no dam, SJF has no dam




############ rf modeling of storage
library(data.table)
library(lubridate)
library(CAST)
library(caret)
library(sf)
library(zoo)

basinSymbol = 'TRM'
basinName = paste0(basinSymbol, '_atOutlet')
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
#pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
#pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
infOrFNF = 'fnf' #'inflow'
if(infOrFNF == 'fnf')	{
	inflow = as.data.table(read.csv(paste0(
		"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}	else	{
	if(infOrFNF == 'inflow')	{
	inflow = as.data.table(read.csv(paste0(
				"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}}
stor = read.csv(paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = inflow[stor[c('Date','stor')], on = 'Date']
allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$resLoss = c(diff(allDat$stor) - allDat$AcFtPDay_inflow[-1], NA) * -1


#######################################################################################################################################
### for monthly estimations


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, resLoss = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$resLoss, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$resLoss = as.numeric(monthDat$resLoss)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(6, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_inflow, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,3,4)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,3,4)], max(lagLength))

resSubset = monthDat[, c('ymonth','resLoss','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_inflow, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$resLoss,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=2000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$resLoss))
plot(resValid$resLoss, predValid)



# old versions


basinName = 'EXC_atOutlet'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")

	
#testing works for EXC, SHA, CLE, TRM, PNF, ISB, EXC
FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05")
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9

plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1])
lines(diff(allDat$stor))
plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1], diff(allDat$stor))










### for daily estimations
	# removing nas; need a better method moving forward
while(any(is.na(allDat$AcFtPDay_fnf)))	{
	allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))] = allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))+1]
}
while(any(is.na(allDat$AcFtPDay_inflow)))	{
	allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))] = allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))+1]
}
while(any(is.na(allDat$AcFtPDay_outflow)))	{
	allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))] = allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))+1]
}
while(any(is.na(allDat$stor)))	{
	allDat$stor[which(is.na(allDat$stor))] = allDat$stor[which(is.na(allDat$stor))+1]
}

recLength = nrow(allDat)
lagLength = seq(1, 365*2, 7)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(allDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(allDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,3,6,10)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,3,6,10)], max(lagLength))

resSubset = allDat[, c('Date','AcFtPDay_outflow')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(allDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(allDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resSubset$doy = yday(resSubset$Date)
resSubset$month = month(resSubset$Date)

resComplete = resSubset[complete.cases(resSubset)]
resTest = subset(resComplete, Date < as.Date("2016-12-31"))

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars

resValid = subset(resComplete, Date >= as.Date("2016-12-31"))
predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)


#######################################################################################################################################
### for monthly estimations

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, AcFtPDay_outflow = NA, AcFtPDay_fnf = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_outflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_fnf, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$AcFtPDay_outflow = as.numeric(monthDat$AcFtPDay_outflow)
monthDat$AcFtPDay_fnf = as.numeric(monthDat$AcFtPDay_fnf)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(1, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,4,6)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,4,6)], max(lagLength))

resSubset = monthDat[, c('ymonth','AcFtPDay_outflow','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)

>>>>>>> 1f901d54a7a7a166b2fc66def008c84eb38db122
=======

	# example input data
basinSymbol = 'PNF'
basinName = paste0(basinSymbol, '_atOutlet')
infOrFNF = 'fnf' # fnf or inflow	
gageLonLat = 	 c(-119.318, 36.845)
				# for ENG: c(-121.270278, 39.240278)
				# for ISB: c(-118.479, 35.650) 
				# for PNF: c(-119.318, 36.845)
				# for TRM: c(-118.998, 36.416) 
				# for EXC: c(-120.264, 37.591) 
				# for ORO: c(-121.480, 39.540) 
				# for CLE: c(-122.765, 40.8225)
				# for MIL: c(-119.6545, 37.0425)
				# for BND: c(-122.185556, 40.288611) 
				# for SHA: c(-122.417, 40.720)
				# for FOL: c(-121.155, 38.710)

				# for SJF: c(-119.697,37.004) 
				# for YRS: c(-121.292500, 39.223611) 
				# for LEW: c(-122.800, 40.732)
historicStreamflowFileLoc = "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=PNF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ENG: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ORO: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for EXC: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for TRM: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=TRM&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for PNF: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=PNF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ISB: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ISB&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for BND: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=BND&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for SHA: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SHA&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for CLE: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=CLE&SensorNums=8&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for MIL: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=MIL&SensorNums=76&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for FOL: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=FOL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
		# baddies
	# for YRS: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=YRS&SensorNums=8&dur_code=D&Start=1906-01-01&End=2022-08-01"
	# for SJF: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SJF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for LEW: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=LEW&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
		# to search for Cali reservoirs: https://cdec.water.ca.gov/dynamicapp/wsSensorData

basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
seas5DataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-seas5.nc'
cfsDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-cfs.nc'
era5DataHistoricalNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-hist-era5.nc'
era5DataRecentNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-recent-era5.nc'
seas5MultiDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\testing-multiple-forecasts-seas5-wy2019.nc'
climateInputsFileLoc = paste0(dataOut_location, basinName, '_processedClimateData')
forecastDate = '28JUL2022'
	
	
	# reading in the ncdfs to ensure timing / structure 
library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
library(data.table)		# for data.table

	### this section is purely for inspecting data before running models
#########################################################################
ncin_cfs = nc_open(cfsDataNCDF)		;	ncin_cfs
ncin_era5 = nc_open(era5DataHistoricalNCDF)	;	ncin_era5
ncin_recentEra5 = nc_open(era5DataRecentNCDF)	;	ncin_recentEra5
ncin_seas5 = nc_open(seas5DataNCDF)	;	ncin_seas5

ncatt_get(ncin_cfs, 'time','units')$value	# says hours but seems to be in days
ncatt_get(ncin_era5, 'time','units')$value
ncatt_get(ncin_recentEra5, 'time','units')$value
ncatt_get(ncin_seas5, 'valid_time','units')$value
#########################################################################

	# correct dates must be manually selected for now
cfsStartDate = as.Date('2022-02-28') #  + ncvar_get(ncin_cfs, 'time')/24 for calculating actual dates
era5StartDate =  as.Date('2001-06-01') # + ncvar_get(ncin_era5, 'time') for calculating actual dates 
era5RecentStartDate =  as.Date('2020-06-01') # + ncvar_get(ncin_recentEra5, 'time') for calculating actual dates 
	# seas5 is incorrectly showing the second of the month, but should be the first
seas5StartDate = as.Date('2022-06-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates
seas5MultiStartDate = as.Date('1993-01-02') - 2 # + ncvar_get(ncin_seas5, 'lead_time') for calculating actual dates


#########################################################################################################
	# step 1
	# delineate a new basin
basinDelineation_f(
	gageLonLat = gageLonLat ,
	basinATLAS_locationAndFile = basinATLAS_locationAndFile,
	dataOut_location = dataOut_location,
	basinName = basinName)


#########################################################################################################
	# step 2
	# import and convert historical climate data
climateDataframe =	climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = era5DataHistoricalNCDF,
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = era5StartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	
	
saveRDS(climateDataframe, paste0(climateInputsFileLoc, 'Historical.RData'))

#########################################################################################################
	# step 3 
	# calibrate model
modelCalibration_f(
	climateInputsFileLoc = climateInputsFileLoc,			#'file_location_and_name.RData',
	historicStreamflowFileLoc = historicStreamflowFileLoc,	#'https://someplace.gov',
	pathToWatershedsGPKG = pathToWatershedsGPKG,			#'file_location_and_name.gpkg',
	dataOut_location = dataOut_location,					#'save_file_location',
	dataSource = 1,											# 1 for FNF from cal.gov, 
	numberOfRuns = 200000,
	targetMetric = 1, 										# 1 = KGE, 2 = NSE, 3 = MAE, 4 = RMSE, 5 = bias
	targetMetricValue = 0.81,								# threshold for considering a value good
	minGoodRuns = 100,										# number of 'good' calibrations before the routine stops
	sfcf = c(runif(5000, .2, 1), runif(5000, 1, 3)),			#snowfall correction factor [-]
	tr   = runif(10000, -6, 5),								#solid and liquid precipitation threshold temperature [C]
	tt   = runif(10000, -5, 6),								#melt temperature [C]
	fm   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 8))),	#snowmelt factor [mm/C]
	fi   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 10))),	#icemelt factor [mm/C]
	fic  = runif(10000, 2, 10),								#debris-covered icemelt factor [mm/C]
	fc   = c(runif(5000, 25, 150), (runif(5000, 150, 1200))),	#field capacity
	lp   = runif(10000, .2, 0.999),								#parameter to actual ET
	beta_soils = runif(10000, 1, 3),							#beta - exponential value for nonlinear relations between soil storage and runoff
	k0   = c(runif(5000, .05, .5), (runif(5000, .5, 0.999))),		#top bucket drainage
	k1   = c(runif(5000, .005, .09), (runif(5000, .09, .5))),	#middle bucket drainage
	k2   = c(runif(5000, .0001, .01), (runif(5000, .01, .1))),#bottom bucket drainage	
	uz1  = c(runif(5000, .22, 10), (runif(5000, 10, 40))),	#max flux rate from STZ to SUZ in mm/d
	perc = c(runif(5000, .1, .5), (runif(5000, .5, 20)))
)	




#########################################################################################################
####	This section is for running and making projections with a previously calibrated model
####
#########################################################################################################
					
#########################################################################################################
	# step 4
	# import and convert projection climate data
recentEra5ClimateDataframe = climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = era5DataRecentNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = era5RecentStartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'era5', # # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	
	
saveRDS(recentEra5ClimateDataframe, paste0(climateInputsFileLoc, 'Recent.RData'))

#cfsClimateDataframe = climateInputConversion_f(
#	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
#	pathToWatershedsGPKG = pathToWatershedsGPKG,
#	basinName = basinName,
#	climateDataNCDF = cfsDataNCDF,	####!!!!! change btw cfs and seas5
#	tempConversionFactor = NA,
#	pptConversionFactor = NA,
#	avgTempGiven = FALSE, 
#	multipleModels = FALSE,	# are there multiple models that need to be stored?
#	startDate = cfsStartDate, 	# when does the clock of the netcdf start?
#	timeToDaysConversion = 24,	# convert time increments to days if necessary
#	dataOut_location = dataOut_location,
#	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
#	variableOrderOption = 'cfs', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
#	precipName = 'tp')	# other options include: tp, tp_sum	

#saveRDS(cfsClimateDataframe, paste0(climateInputsFileLoc, 'CFS.RData'))

seas5ClimateDataframe = climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = seas5DataNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = seas5StartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'seas5', # 'era5' = [longitude,latitude,time]; 'cfs' = [longitude, latitude, member, step]; 'seas5' = [longitude, latitude, member, lead_time] for tmax and tmin but [lead_time, longitude, latitude, member] for tp_sum]
	precipName = 'tp_sum')	# other options include: tp, tp_sum	

saveRDS(seas5ClimateDataframe, paste0(climateInputsFileLoc, 'SEAS5.RData'))


	# step 5 
	# run the model with forecasting data
	## Running the Model for Seasonal Forecasts 
allForecastsOutput = seasonalForecast_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,

saveRDS(allForecastsOutput, paste0(climateInputsFileLoc, 'ForecastsFor_', forecastDate, '.RData'))





#########################################################################################################
####	This section is for validating model performance 
####
#########################################################################################################

	# step 6
	# validate historical data and generate plots
validationHistoricalOutput = validationAndPlotGeneration_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,


	# step 7
	# import and convert multi projection climate data for validationAndPlotGeneration_f
climateInputConversion_f(
	pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	basinName = basinName,
	climateDataNCDF = seas5MultiDataNCDF,	####!!!!! change btw cfs and seas5
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = seas5MultiStartDate, 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = dataOut_location,
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = 'seas5Multi', 
	precipName = 'tp_sum')	# other options include: tp, tp_sum	


	# step 8
	# validate streamflow projection data and generate plots
validationProjectionOutput = projectionValidationAndPlotGeneration_f(
	basinName = basinName,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,



#########################################################################################################
####	This section is for linking streamflow models to reservoirs
####
#########################################################################################################
	# step 9
	# ML for predicting reservoir outflows
	
	# note: plans to remove infOrFNF and only use 'inflow' going forward
	infOrFNF = 'inflow' # fnf or inflow	

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = basinName,
	basinSymbol = basinSymbol,
	infOrFNF = infOrFNF,
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=500,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 5,
	numRepeats = 30)

	# step 10
	# validation of combined model projections of total storage
validationProjectedStorageOutput = storageProjectionValidationAndPlotGeneration_f(
	basinName = basinName,
	infOrFNF = infOrFNF,
	climateInputsFileLoc = climateInputsFileLoc,	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = pathToWatershedsGPKG,
	historicStreamflowFileLoc = historicStreamflowFileLoc,
	dataOut_location = dataOut_location,
	dataSource = 1)							# 1 for FNF from cal.gov,

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# works for EXC, SHA, CLE, TRM, PNF, ISB, ORO

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'EXC_atOutlet',
	basinSymbol = 'EXC',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'TRM_atOutlet',
	basinSymbol = 'TRM',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)

reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'SHA_atOutlet',
	basinSymbol = 'SHA',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'CLE_atOutlet',
	basinSymbol = 'CLE',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'PNF_atOutlet',
	basinSymbol = 'PNF',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ISB_atOutlet',
	basinSymbol = 'ISB',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ORO_atOutlet',
	basinSymbol = 'ORO',
	infOrFNF = 'fnf',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


# works for FOL, MIL, LEW, ENG
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'FOL_atOutlet',
	basinSymbol = 'FOL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'MIL_atOutlet',
	basinSymbol = 'MIL',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'LEW_atOutlet',
	basinSymbol = 'LEW',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)
reservoirOutflowCalVal = reservoirOutflowCalVal_f(
	dataOut_location = dataOut_location,
	basinName = 'ENG_atOutlet',
	basinSymbol = 'ENG',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=2000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 4,
	numRepeats = 20)


	
	
# works for FOL, MIL, LEW, ENG
#FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor =  as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ENG&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

#FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
#FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = stor$VALUE
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = stor[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

plot(allDat$stor, allDat$outflow)
plot(allDat$stor, allDat$inflow)
#plot(allDat$inflow, allDat$fnf)
plot(allDat$inflow, allDat$outflow)
plot(allDat$Date, allDat$outflow, type='l', col='red', ylim=c(0,10000))
lines(allDat$Date, allDat$inflow, type='l', col='blue')


# BND has no dam, YRS has no dam, SJF has no dam




############ rf modeling of storage
library(data.table)
library(lubridate)
library(CAST)
library(caret)
library(sf)
library(zoo)

basinSymbol = 'TRM'
basinName = paste0(basinSymbol, '_atOutlet')
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
#pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
#pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
infOrFNF = 'fnf' #'inflow'
if(infOrFNF == 'fnf')	{
	inflow = as.data.table(read.csv(paste0(
		"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}	else	{
	if(infOrFNF == 'inflow')	{
	inflow = as.data.table(read.csv(paste0(
				"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")))
	}}
stor = read.csv(paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", basinSymbol, "&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05"))
#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = inflow[stor[c('Date','stor')], on = 'Date']
allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$resLoss = c(diff(allDat$stor) - allDat$AcFtPDay_inflow[-1], NA) * -1


#######################################################################################################################################
### for monthly estimations


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, resLoss = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$resLoss, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$resLoss = as.numeric(monthDat$resLoss)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(6, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_inflow, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$resLoss,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,3,4)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,3,4)], max(lagLength))

resSubset = monthDat[, c('ymonth','resLoss','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_inflow, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$resLoss,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=2000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$resLoss))
plot(resValid$resLoss, predValid)



# old versions


basinName = 'EXC_atOutlet'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")

	
#testing works for EXC, SHA, CLE, TRM, PNF, ISB, EXC
FNF = as.data.table(read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-05"))
stor = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=15&dur_code=D&Start=1900-01-01&End=2022-08-05")
outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")
inflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-05")

FNF$Date = ymd(unlist(strsplit(FNF$DATE.TIME, " "))[seq(1,nrow(FNF)*2,2)])
FNF$fnf = FNF$VALUE
stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
stor$stor = as.numeric(stor$VALUE)
outflow$Date = ymd(unlist(strsplit(outflow$DATE.TIME, " "))[seq(1,nrow(outflow)*2,2)])
outflow$outflow = outflow$VALUE
inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
inflow$inflow = inflow$VALUE

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9

plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1])
lines(diff(allDat$stor))
plot(allDat$AcFtPDay_inflow[-1] - allDat$AcFtPDay_outflow[-1], diff(allDat$stor))










### for daily estimations
	# removing nas; need a better method moving forward
while(any(is.na(allDat$AcFtPDay_fnf)))	{
	allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))] = allDat$AcFtPDay_fnf[which(is.na(allDat$AcFtPDay_fnf))+1]
}
while(any(is.na(allDat$AcFtPDay_inflow)))	{
	allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))] = allDat$AcFtPDay_inflow[which(is.na(allDat$AcFtPDay_inflow))+1]
}
while(any(is.na(allDat$AcFtPDay_outflow)))	{
	allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))] = allDat$AcFtPDay_outflow[which(is.na(allDat$AcFtPDay_outflow))+1]
}
while(any(is.na(allDat$stor)))	{
	allDat$stor[which(is.na(allDat$stor))] = allDat$stor[which(is.na(allDat$stor))+1]
}

recLength = nrow(allDat)
lagLength = seq(1, 365*2, 7)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(allDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(allDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		allDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,3,6,10)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,3,6,10)], max(lagLength))

resSubset = allDat[, c('Date','AcFtPDay_outflow')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(allDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(allDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resSubset$doy = yday(resSubset$Date)
resSubset$month = month(resSubset$Date)

resComplete = resSubset[complete.cases(resSubset)]
resTest = subset(resComplete, Date < as.Date("2016-12-31"))

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars

resValid = subset(resComplete, Date >= as.Date("2016-12-31"))
predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)


#######################################################################################################################################
### for monthly estimations

allDat = FNF[stor[c('Date','stor')], on = 'Date']
allDat = allDat[outflow[c('Date','outflow')], on='Date']
allDat = allDat[inflow[c('Date','inflow')], on='Date']

allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
allDat$AcFtPDay_outflow = as.numeric(allDat$outflow) * (60*60*24) / 43559.9
allDat$AcFtPDay_fnf = 	  as.numeric(allDat$fnf)     * (60*60*24) / 43559.9


allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, AcFtPDay_outflow = NA, AcFtPDay_fnf = NA)
k = 0
for(i in unique(allDat$ymonth))	{
	k = k + 1
	subsetDat = subset(allDat, ymonth == i)
 	monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
		mean(subsetDat$stor, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_inflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_outflow, na.rm=TRUE),
		mean(subsetDat$AcFtPDay_fnf, na.rm=TRUE)))
}
monthDat$stor = as.numeric(monthDat$stor)
monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
monthDat$AcFtPDay_outflow = as.numeric(monthDat$AcFtPDay_outflow)
monthDat$AcFtPDay_fnf = as.numeric(monthDat$AcFtPDay_fnf)

plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1])
lines(diff(monthDat$stor))
plot(monthDat$AcFtPDay_inflow[-1] - monthDat$AcFtPDay_outflow[-1], diff(monthDat$stor))



recLength = nrow(monthDat)
lagLength = seq(1, 12*2, 1)

corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
k=0
for(i in lagLength)	{
	k=k+1
	corResults[k, 1] = cor.test(
		c(rep(NA, i), rollmean(monthDat$stor, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
	corResults[k, 2] = cor.test(
		c(rep(NA, i), rollmean(monthDat$AcFtPDay_fnf, i)[-c((recLength - (i-1)): recLength)]),
		monthDat$AcFtPDay_outflow,
		method='kendall')$p.value
}

bestInfLags = c(order(corResults[,1])[c(1,2,4,6)], max(lagLength))
bestStorLags = c(order(corResults[,2])[c(1,2,4,6)], max(lagLength))

resSubset = monthDat[, c('ymonth','AcFtPDay_outflow','month')]

for(j in 1:length(bestInfLags))	{
	infLag = bestInfLags[j]
	storLag = bestStorLags[j]

	resSubset[ , paste0('stor_mn_', storLag)] = c(rep(NA, storLag), rollmean(monthDat$stor, storLag)[-c((recLength - (storLag-1)): recLength)])
	resSubset[ , paste0('inf_mn_', infLag)] = c(rep(NA, infLag), rollmean(monthDat$AcFtPDay_fnf, infLag)[-c((recLength - (infLag-1)): recLength)])
}

resComplete = resSubset[complete.cases(resSubset), ]
resTest = resComplete[1:(floor(nrow(resComplete) / (4/3))), ]
resValid = resComplete[(floor(nrow(resComplete) / (4/3)) + 1):nrow(resComplete), ]

model_LLO = train(resTest[,-c(1,2)],
					resTest$AcFtPDay_outflow,
					metric="Rsquared",
					method='rf',
					#tuneLength=?,
					imporance=TRUE,
					ntree=5000,
					trControl = trainControl(
						method='repeatedcv',
						number=5,
						repeats=8))
						#method='cv'))
	#					index = indices$index))
	#					p=.75))

model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
varImp(model_LLO)$importance	
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


predValid = predict(model_LLO, newdata = resValid)

summary(lm(predValid ~ resValid$AcFtPDay_outflow))
plot(resValid$AcFtPDay_outflow, predValid)

>>>>>>> 90a879e5efdf100422fb9d469711758a65331831
