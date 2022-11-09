
##########################################################################################################
	## basin delineation
	
basinDelineation_f = function(
	gageLonLat = c(-79,36),
	basinATLAS_locationAndFile = 'basinATLAS_location.gdb',
	dataOut_location = 'save_file_location',
	basinName = 'inster_basin_or_outlet_name')
	{

		# create the folder for storing outputs
	if(!file.exists(paste0(dataOut_location)))	{
			dir.create(file.path(paste0(dataOut_location)))
	}
	# check to make sure the basin has not already been delineated
	if(file.exists(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")))	{
		print("This file already exists you big ole dummy!")
	}	else	{

		# load relevant libraries
	library(sf)				# for geospatial data
	sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat; this does not matter for the scale we are interested in
	library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
	
	
	##################################################
	#### reading basinATLAS data and intersecting with gage location
		#BasinATLAS available via hydrosheds.org
	basinAt12 = st_read(dsn=basinATLAS_locationAndFile, layer="BasinATLAS_v10_lev12")
	gageLocation_asSF = st_sfc(st_point(gageLonLat))
	st_crs(gageLocation_asSF) = 4326
	gageIntersection = as.numeric(st_intersects(gageLocation_asSF, st_buffer(basinAt12,0))) # keeps points
	
		# only the adjacent upstream watersheds is identified in the database, so first we identify the most downstream basin according to the intersection above
	basinIDs = basinAt12$HYBAS_ID[gageIntersection]
	basinAt12Remaining = basinAt12[-gageIntersection,]
	
		# then we search the database iteratively to identify each watershed upstream of the newly identified watershed, until there are not more upstream watersheds
	while(any(basinIDs %in% basinAt12Remaining$NEXT_DOWN))	{
		for(trueys in which(basinIDs %in% basinAt12Remaining$NEXT_DOWN))	{
			newIDRows = which(basinAt12Remaining$NEXT_DOWN == basinIDs[trueys])
			newIDs = basinAt12Remaining$HYBAS_ID[newIDRows]
			basinIDs = c(basinIDs, newIDs)
			print(basinIDs)	# this takes a while esp for large basins, so printing to ensure progress is being made
			basinAt12Remaining = basinAt12Remaining[-newIDRows,]
		}
	}

	rm(basinAt12Remaining)
	basinWatersheds = subset(basinAt12, HYBAS_ID %in% basinIDs)
	rm(basinAt12)

	# creating new polygon based on merging all the upstream hydrobasins
	delineatedBasin = st_sf(st_union(basinWatersheds))
	
		# save the watersheds both indivudually and as a single basin (st_union above) 
	st_write(delineatedBasin, paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg"), append=FALSE)
	st_write(basinWatersheds, paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"), append=FALSE)
	}
}	






#################################################################################################################
	# function for selecting climate inputs and converting to appropriate units

climateInputConversion_f = function(
	basinName = 'inster_basin_or_outlet_name',
	climateDataNCDF = 'file_location_and_name.nc',
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	startDate = as.Date("1990-01-01"), 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = 'save_file_location',
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = NA, # 'seas5' or 'era5' or 'seas5Multi'
	precipName = 'tp',	# other options include: tp, tp_sum	
	limitedModels = 25)	# the full seas5 multi only goes back to ~ 2017, so for hindcasting we may want to limit analysis to the first 25 models (which seem to be present for the whole record)
	{
	# load relevant libraries
	library(sf)				# for geospatial data
	sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat
	library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
	library(data.table)		# for data.frame and fread
	if(optionForPET == 1) {library(EcoHydRology)}	# for PET_fromTemp

		# creating the folder if it does not exist
	if(!file.exists(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))) {
		print('Gotta delineate the basin first!!')
	} else {
		
			# read in the previously delineated basins
		basinWatersheds = st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))

			# identifying area of each subbasin for rescaling the 
		basinArea = sum(basinWatersheds$SUB_AREA)

		######################################################################################################################
			# routine for seas5 monthly forecast data
		if(variableOrderOption == 'seas5')	{
				# reading in the ncdfs
			ncin = nc_open(climateDataNCDF)

				# current post processing of ncdfs results in some nas, so must remove these before proceeding
					# na removal by mdeian value needs to be revisited... a linear interp would probably be better
			ncTmin = ncvar_get(ncin, 't2m_min'); if(any(is.na(ncTmin)))	{ncTmin[is.na(ncTmin)] = median(ncTmin, na.rm=TRUE)} # [longitude,latitude,lead_time,member] 
			ncTmax = ncvar_get(ncin, 't2m_max'); if(any(is.na(ncTmax)))	{ncTmax[is.na(ncTmax)] = median(ncTmax, na.rm=TRUE)} # [longitude,latitude,lead_time,member] 
			ncPPT = ncvar_get(ncin, precipName); if(any(is.na(ncPPT)))	{ncPPT[is.na(ncPPT)] = median(ncPPT, na.rm=TRUE)}	 # [longitude,latitude,lead_time,member] 
				# current CAi climate data pipeline does not include tavg, so must estimate using tmin and tmax
			if(avgTempGiven)	{ncTavg = ncvar_get(ncin, 't2m')
			}	else {ncTavg = (ncTmin + ncTmax) / 2}	

				# assuming all lat / lon structures are the same
			nc_lat = ncvar_get(ncin, 'latitude')
			nc_lon = ncvar_get(ncin, 'longitude')
			nc_date = startDate + timeToDaysConversion * ncvar_get(ncin, 'lead_time') # time is days after jan 1 1990

				# list for storing time series of climate data
			allClimateData = list()
			for(numModels in 1:length(ncvar_get(ncin, 'member')))	{
				allPPT = 0	; allTmin = 0	; allTmax = 0	; allPET = 0
				for(numberOfWatersheds in 1:nrow(basinWatersheds))	{
						#identify climate data closest to the centroid of the subbasin of interest
					thisLonLat = st_coordinates(st_centroid(basinWatersheds[numberOfWatersheds, ]))
					nearestLon = which.min(abs(thisLonLat[1] - nc_lon))
					nearestLat = which.min(abs(thisLonLat[2] - nc_lat))
				
						# Tmax is sometimes < Tmin, so need to adjust before running for PET
					theseTmin = ncTmin[nearestLon, nearestLat, , numModels]	# [longitude,latitude,lead_time,member]
					theseTmax = ncTmax[nearestLon, nearestLat, , numModels]	# [longitude,latitude,lead_time,member]
					if(any(theseTmin >= theseTmax))	{
						theseTmin[theseTmin >= theseTmax] = theseTmax[theseTmin >= theseTmax] - 0.1
					}
					
						# summing all climate variables; they normalized by basin area below
							# ultimately, subbasins should be disaggregated and the model run on subcomponents, but this will take time to implement
					allTmin = allTmin +  theseTmin *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					allTmax = allTmax + theseTmax *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					allPPT = allPPT + ncPPT[nearestLon, nearestLat, , numModels] *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# [longitude,latitude,lead_time,member]
						# since we don't currently ingest PET data, we arecalculating from a penman monteith eq
					if(optionForPET == 1)	{
						allPET = allPET + PET_fromTemp(yday(nc_date), theseTmax , theseTmin,
							lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
					}	else	{
						print('need to figure this one out later')
					}
				}

				# normalizing data by basin area
				avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
				avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
				avgTavg = (avgTmin + avgTmax) / 2
				avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
				avgPET = allPET / basinArea

					# all climate data to be saved
				allClimateData[[numModels]] = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
				saveRDS(allClimateData, paste0(dataOut_location, 'SEAS5_', basinName, '.RData'))
			}
		}

		######################################################################################################################
			# routine for era5 data
		if(variableOrderOption == 'era5')	{
				# reading in the ncdfs
			ncin = nc_open(climateDataNCDF)

				# current post processing of ncdfs results in some nas, so must remove these before proceeding
			ncTmin = ncvar_get(ncin, 't2m_min'); if(any(is.na(ncTmin)))	{ncTmin[is.na(ncTmin)] = median(ncTmin, na.rm=TRUE)} 	#[longitude,latitude,time] 
			ncTmax = ncvar_get(ncin, 't2m_max'); if(any(is.na(ncTmax)))	{ncTmax[is.na(ncTmax)] = median(ncTmax, na.rm=TRUE)}	#[longitude,latitude,time] 
			ncPPT = ncvar_get(ncin, precipName); if(any(is.na(ncPPT)))	{ncPPT[is.na(ncPPT)] = median(ncPPT, na.rm=TRUE)} 		#[longitude,latitude,time] 
					# current CAi climate data pipeline does not include tavg, so must estimate using tmin and tmax
			if(avgTempGiven)	{ncTavg = ncvar_get(ncin, 't2m_avg')
			}	else {ncTavg = (ncTmin + ncTmax) / 2}	

				# assuming all lat / lon structures are the same
			nc_lat = ncvar_get(ncin, 'latitude')
			nc_lon = ncvar_get(ncin, 'longitude')
			nc_date = startDate + timeToDaysConversion * ncvar_get(ncin, 'time') # time is days after jan 1 1990

			allPPT = 0	; allTmin = 0	; allTmax = 0	; allPET = 0
			for(numberOfWatersheds in 1:nrow(basinWatersheds))	{
				thisLonLat = st_coordinates(st_centroid(basinWatersheds[numberOfWatersheds, ]))
				nearestLon = which.min(abs(thisLonLat[1] - nc_lon))
				nearestLat = which.min(abs(thisLonLat[2] - nc_lat))
						
					# summing all climate variables; they normalized by basin area below
						# ultimately, subbasins should be disaggregated and the model run on subcomponents, but this will take time to implement
				allTmin = allTmin + ncTmin[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
				allTmax = allTmax + ncTmax[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
				allPPT = allPPT + ncPPT[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					# since we don't currently ingest PET data, we arecalculating from a penman monteith eq
				if(optionForPET == 1)	{
					allPET = allPET + PET_fromTemp(yday(nc_date), ncTmax[nearestLon, nearestLat, ], ncTmin[nearestLon, nearestLat, ],
						lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
				}	else	{
					print('need to figure this one out later')
				}
			}
			
				# normalizing data by basin area
			avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
			avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
			avgTavg = (avgTmin + avgTmax) / 2
			avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
			avgPET = allPET / basinArea
			
				# save data output
			allClimateData = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
			saveRDS(allClimateData, paste0(dataOut_location, 'ERA5_', basinName, '.RData'))
		}
		
		######################################################################################################################
			# for long record of seas5 data for hindcasting
		if(variableOrderOption == 'seas5Multi')	{
			# creating the folder to hold the output
				if(!file.exists(paste0(dataOut_location, 'seas5MultiOutput')))	{
				dir.create(file.path(paste0(dataOut_location, 'seas5MultiOutput')))
			}

			
			# reading in the ncdfs
			ncin = nc_open(climateDataNCDF) # for tmax and tmin [longitude,latitude,member,lead_time,init_time] 
											# for tp_sum [lead_time,longitude,latitude,member,init_time]

			ncTmin = ncvar_get(ncin, 't2m_min'); if(any(is.na(ncTmin)))	{ncTmin[is.na(ncTmin)] = median(ncTmin, na.rm=TRUE)} 
			ncTmax = ncvar_get(ncin, 't2m_max'); if(any(is.na(ncTmax)))	{ncTmax[is.na(ncTmax)] = median(ncTmax, na.rm=TRUE)} 
			ncPPT = ncvar_get(ncin, precipName); if(any(is.na(ncPPT)))	{ncPPT[is.na(ncPPT)] = median(ncPPT, na.rm=TRUE)} 
			if(avgTempGiven)	{ncTavg = ncvar_get(ncin, 't2m_avg')
			}	else {ncTavg = (ncTmin + ncTmax) / 2}	

				# assuming all lat / lon structures are the same
			nc_lat = ncvar_get(ncin, 'latitude')
			nc_lon = ncvar_get(ncin, 'longitude')

			initTimes = ncvar_get(ncin, 'init_time')
			leadTimes = ncvar_get(ncin, 'lead_time')
			for(thisInitTime in 1:length(initTimes))	{
				
				nc_date = startDate + timeToDaysConversion * leadTimes + timeToDaysConversion * initTimes[thisInitTime] # time is days after jan 1 

				allClimateData = list()
				
					# the full seas5 multi only goes back to ~ 2017, so for hindcasting we may want to limit analysis to the first 25 models (which seem to be present for the whole record) 
				if(limitedModels) 	{
					totModels = 1:limitedModels
				} else	{
					totModels = 1:length(ncvar_get(ncin, 'member'))
				}
				
				for(numModels in totModels)	{
					allPPT = 0	; allTmin = 0	; allTmax = 0	; allPET = 0
					for(numberOfWatersheds in 1:nrow(basinWatersheds))	{
							#identify climate data closest to the centroid of the subbasin of interest
						thisLonLat = st_coordinates(st_centroid(basinWatersheds[numberOfWatersheds, ]))
						nearestLon = which.min(abs(thisLonLat[1] - nc_lon))
						nearestLat = which.min(abs(thisLonLat[2] - nc_lat))
					
							# Tmax is sometimes < Tmin, so need to adjust before running for PET
							# this step can be removed when the bias correction / post processing is standardized
						theseTmin = ncTmin[nearestLon, nearestLat, numModels, , thisInitTime]
						theseTmax = ncTmax[nearestLon, nearestLat, numModels, , thisInitTime]
						if(any(theseTmin >= theseTmax))	{
							theseTmin[theseTmin >= theseTmax] = theseTmax[theseTmin >= theseTmax] - 0.1
						}
						
							# summing all climate variables; they normalized by basin area below
							# ultimately, subbasins should be disaggregated and the model run on subcomponents, but this will take time to implement
						allTmin = allTmin +  theseTmin *  basinWatersheds$SUB_AREA[numberOfWatersheds]
						allTmax = allTmax + theseTmax *  basinWatersheds$SUB_AREA[numberOfWatersheds]
						allPPT = allPPT + ncPPT[ , nearestLon, nearestLat, numModels, thisInitTime] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
							# since we don't currently ingest PET data, we arecalculating from a penman monteith eq
						if(optionForPET == 1)	{
							allPET = allPET + PET_fromTemp(yday(nc_date), theseTmax , theseTmin,
								lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
						}	else	{
							print('need to figure this one out later')
						}
					}
					
					# normalizing data by basin area
					avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
					avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
					avgTavg = (avgTmin + avgTmax) / 2
					avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
					avgPET = allPET / basinArea

						# all climate data to be saved
					allClimateData[[numModels]] = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
				}
				saveRDS(allClimateData, paste0(dataOut_location, 'seas5MultiOutput\\SEAS5_startDate_', nc_date[1], '.RData'))
			}
		}
	}
}	
	
	
	
##########################################################################################################
	## Function to run HBV
runHBV_f = function(
	climateInput = climateInput,
	sfcf = 1,	#snowfall correction factor [-]
	tr   = 0,	#solid and liquid precipitation threshold temperature [C]
	tt   = 0,	#melt temperature [C]
	fm   = 2,	#snowmelt factor [mm/C]
	fi   = 2,	#icemelt factor [mm/C]
	fic  = 3,	#debris-covered icemelt factor [mm/C]
	fc   = 200, # field capacity 
	lp   = 0.75, # parameter for PET --> AET
	beta_soils = 1.5, #soil moisture drainage exponent
	k0   = 0.5,	# storage constant of top bucket
	k1   = 0.1,	# storage constant of middle bucket
	k2   = 0.01,	# storage constant of bottom bucket	
	uz1  = 5, #max flux rate from STZ to SUZ in mm/d
	perc = 2)  # max flux rate from SUZ to SLZ in mm/d
{
	library(HBV.IANIGLA) # for running HBV model
        
    pptModuleOutput =	cbind(climateInput,	
                            SnowGlacier_HBV(model = 1,
								inputData = cbind(climateInput$Tavg, climateInput$PPT),
                                initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                param = c(sfcf, tr, tt, fm, fi, fic)
							)	)        

    rechModuleOutput =	cbind(pptModuleOutput,
                            Soil_HBV(
                                model = 1,
                                inputData = cbind(pptModuleOutput$Total, pptModuleOutput$PET),
                                initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                                param = c(fc, lp, beta_soils)
							)	)
        
        
    allHBVoutput =		cbind(rechModuleOutput,
                            Routing_HBV(
								model = 1,	# model=1 gives three stores with routing for each
                                lake = FALSE,
                                inputData = cbind(rechModuleOutput$Rech),	# recharge time series
                                initCond = c(10,10,10),	# initial storage in each reservoir in mm
                                param = c(k0, k1, k2, uz1, perc)
							)	)
    return(allHBVoutput)
}    
    
	

##########################################################################################################
	## Calibration Function
modelCalibration_f = function(
	historicStreamflowFileLoc = 'https://someplace.gov',
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	dataOut_location = 'save_file_location',
	dataSource = 1,											# 1 for FNF from cal.gov, 
	numberOfRuns = 100000,									# maximum number of runs
	targetMetric = 1, 										# 1 = KGE, 2 = NSE, 3 = MAE, 4 = RMSE, 5 = bias
	targetMetricValue = 0.81,								# threshold for considering a value good
	minGoodRuns = 200,										# number of 'good' calibrations before the routine stops
	sfcf = c(runif(5000, .2, 1), runif(5000, 1, 3)),			#snowfall correction factor [-]
	tr   = runif(10000, -6, 5),								#solid and liquid precipitation threshold temperature [C]
	tt   = runif(10000, -5, 6),								#melt temperature [C]
	fm   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 8))),	#snowmelt factor [mm/C]
	fi   = c(runif(5000, .2, 1.5), (runif(5000, 1.5, 10))),	#icemelt factor [mm/C]
	fic  = runif(10000, 2, 10),								#debris-covered icemelt factor [mm/C]
	fc   = c(runif(5000, 25, 150), (runif(5000, 150, 1200))),	#field capacity
	lp   = runif(10000, .2, 1),								#parameter to actual ET
	beta_soils = runif(10000, 1, 3),							#beta - exponential value for nonlinear relations between soil storage and runoff
	k0   = c(runif(5000, .05, .5), (runif(5000, .5, 0.999))),		#top bucket drainage
	k1   = c(runif(5000, .005, .09), (runif(5000, .09, .5))),	#middle bucket drainage
	k2   = c(runif(5000, .0001, .01), (runif(5000, .01, .1))),#bottom bucket drainage	
	uz1  = c(runif(5000, .22, 10), (runif(5000, 10, 40))),	#max flux rate from STZ to SUZ in mm/d
	perc = c(runif(5000, .1, .5), (runif(5000, .5, 20))))		#max flux rate from SUZ to SLZ in mm/d
	{

	if(file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		print("This file already exists you big ole dummy!")
	}	else	{


		library(sf)
		library(data.table)
		library(lubridate)
		library(hydroGOF)		# for nse / kge calculations
		
		
			# reading in previously reanalyzed climate data
		climateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

			# reading in and reformatting streamflow input 
		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
				# removing negative streamflow
			if(any(historicStreamflow$historicQinOriginalUnits < 0))	{historicStreamflow$historicQinOriginalUnits[historicStreamflow$historicQinOriginalUnits < 0] = NA}
			basinArea = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
				
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}
		
		
			# merging historic streamflow record onto climate inputs data
		climateAndStreamflowInput = historicStreamflow[climateInput, on='Date']
			# apportioning data for calibration and validation
		lengthOfCalibrationInput = nrow(climateAndStreamflowInput)
		calRows = c(1:(floor(lengthOfCalibrationInput / 3)), ceiling(lengthOfCalibrationInput/(3/2)):lengthOfCalibrationInput)
		valRows = c(ceiling(lengthOfCalibrationInput / 3):floor(lengthOfCalibrationInput/(3/2)))
		
			# dataframe for capturing calibration metrics
		cal_out = data.frame(
			sfcf =	rep(NA,minGoodRuns),
			tr = 	rep(NA,minGoodRuns),
			tt = 	rep(NA,minGoodRuns),
			fm = 	rep(NA,minGoodRuns),
			fi = 	rep(NA,minGoodRuns),
			fic =	rep(NA,minGoodRuns),
			fc =	rep(NA,minGoodRuns),
			lp =	rep(NA,minGoodRuns),
			beta_soils = rep(NA,minGoodRuns),
			k0 = 	rep(NA,minGoodRuns),
			k1 = 	rep(NA,minGoodRuns),
			k2 = 	rep(NA,minGoodRuns),
			uz1 = 	rep(NA,minGoodRuns),
			perc = 	rep(NA,minGoodRuns),
			kgeCalibration = rep(NA,minGoodRuns),
			nseCalibration = rep(NA,minGoodRuns),
			maeCalibration = rep(NA,minGoodRuns),
			rmseCalibration = rep(NA,minGoodRuns),
			biasCalibration = rep(NA,minGoodRuns),
			kgeValidation = rep(NA,minGoodRuns),
			nseValidation = rep(NA,minGoodRuns),
			maeValidation = rep(NA,minGoodRuns),
			rmseValidation = rep(NA,minGoodRuns),
			biasValidation = rep(NA,minGoodRuns),
			kgeAll = rep(NA,minGoodRuns),
			nseAll = rep(NA,minGoodRuns),
			maeAll = rep(NA,minGoodRuns),
			rmseAll = rep(NA,minGoodRuns),
			biasAll = rep(NA,minGoodRuns),
			mnthSumAbsBias = rep(NA,minGoodRuns),
			mnthBias_1 = rep(NA,minGoodRuns),
			mnthBias_2 = rep(NA,minGoodRuns),
			mnthBias_3 = rep(NA,minGoodRuns),
			mnthBias_4 = rep(NA,minGoodRuns),
			mnthBias_5 = rep(NA,minGoodRuns),
			mnthBias_6 = rep(NA,minGoodRuns),
			mnthBias_7 = rep(NA,minGoodRuns),
			mnthBias_8 = rep(NA,minGoodRuns),
			mnthBias_9 = rep(NA,minGoodRuns),
			mnthBias_10 = rep(NA,minGoodRuns),
			mnthBias_11 = rep(NA,minGoodRuns),
			mnthBias_12 = rep(NA,minGoodRuns)
		)

		iter = 0
		while(iter < minGoodRuns)	{
			jj = 0

				# sampling parameter values for calibration
			sfcf =	sample(sfcf, numberOfRuns, replace=TRUE)
			tr = 	sample(tr, numberOfRuns, replace=TRUE)
			tt = 	sample(tt, numberOfRuns, replace=TRUE)
			fm = 	sample(fm, numberOfRuns, replace=TRUE)
			fi = 	sample(fi, numberOfRuns, replace=TRUE)
			fic =	sample(fic, numberOfRuns, replace=TRUE)
			fc =	sample(fc, numberOfRuns, replace=TRUE)
			lp =	sample(lp, numberOfRuns, replace=TRUE)
			beta_soils = sample(beta_soils, numberOfRuns, replace=TRUE)
			k0 = 	sample(k0, numberOfRuns, replace=TRUE)
			k1 = 	sample(k1, numberOfRuns, replace=TRUE)
			k2 = 	sample(k2, numberOfRuns, replace=TRUE)
			uz1 = 	sample(uz1, numberOfRuns, replace=TRUE)
			perc = 	sample(perc, numberOfRuns, replace=TRUE)



			  # since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing 'random' sampling
			if(any(k1 > k0))	{
				k1[which(k1 > k0)] = k0[which(k1 > k0)] * .99
			}
			if(any(k2 > k1))	{
				k2[which(k2 > k1)] = k1[which(k2 > k1)] * .99
			}
			if(any(uz1 < perc))	{
				uz1[which(uz1 < perc)] = perc[which(uz1 < perc)] * 1.01
			}


				# incrementally decreasing the target metric value every n runs
			targetMetricValue = targetMetricValue - 0.01
			print(targetMetricValue)
			print(iter)
			while(jj < numberOfRuns & iter < minGoodRuns) {
				jj = jj+1
					# running HBV
				HBVoutput = runHBV_f(
					climateInput = climateAndStreamflowInput,
					sfcf[jj],	#snowfall correction factor [-]
					tr[jj],	#solid and liquid precipitation threshold temperature [C]
					tt[jj],	#melt temperature [C]
					fm[jj],	#snowmelt factor [mm/C]
					fi[jj],	#icemelt factor [mm/C]
					fic[jj],	#debris-covered icemelt factor [mm/C]
					fc[jj], # field capacity 
					lp[jj], # parameter for PET --> AET
					beta_soils[jj], #soil moisture drainage exponent
					k0[jj],	# storage constant of top bucket
					k1[jj],	# storage constant of middle bucket
					k2[jj],	# storage constant of bottom bucket	
					uz1[jj], #max flux rate from STZ to SUZ in mm/d
					perc[jj])  # max flux rate from SUZ to SLZ in mm/d
			
			
					# identifying if a parameterization meets the criteria to be saved
				if(KGE(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows]) > targetMetricValue &
					KGE(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows]) > targetMetricValue)	{
					
					iter = iter + 1
					plot(HBVoutput$Date[-c(1:1000)], HBVoutput$historicQinmm[-c(1:1000)], main=paste0(jj, ' runs'))
					lines(HBVoutput$Date[-c(1:1000)], HBVoutput$Qg[-c(1:1000)], col='red')
					cal_out$sfcf[iter] = sfcf[jj]
					cal_out$tr[iter] = tr[jj]
					cal_out$tt[iter] = tt[jj]
					cal_out$fm[iter] = fm[jj]
					cal_out$fi[iter] = fi[jj]
					cal_out$fic[iter] = fic[jj]
					cal_out$fc[iter] = fc[jj]
					cal_out$lp[iter] = lp[jj]
					cal_out$beta_soils[iter] = beta_soils[jj]
					cal_out$k0[iter] = k0[jj]
					cal_out$k1[iter] = k1[jj]
					cal_out$k2[iter] = k2[jj]
					cal_out$uz1[iter] = uz1[jj]
					cal_out$perc[iter] = perc[jj]
					cal_out$kgeCalibration[iter] = KGE(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows])
					cal_out$nseCalibration[iter] = NSE(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows])
					cal_out$maeCalibration[iter] = mae(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows])
					cal_out$rmseCalibration[iter] = rmse(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows])
					cal_out$biasCalibration[iter] = pbias(HBVoutput$Qg[calRows], HBVoutput$historicQinmm[calRows])
					cal_out$kgeValidation[iter] = KGE(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows])
					cal_out$nseValidation[iter] = NSE(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows])
					cal_out$maeValidation[iter] = mae(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows])
					cal_out$rmseValidation[iter] = rmse(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows])
					cal_out$biasValidation[iter] = pbias(HBVoutput$Qg[valRows], HBVoutput$historicQinmm[valRows])
					cal_out$kgeAll[iter] = KGE(HBVoutput$Qg, HBVoutput$historicQinmm)
					cal_out$nseAll[iter] = NSE(HBVoutput$Qg, HBVoutput$historicQinmm)
					cal_out$maeAll[iter] = mae(HBVoutput$Qg, HBVoutput$historicQinmm)
					cal_out$rmseAll[iter] = rmse(HBVoutput$Qg, HBVoutput$historicQinmm)
					cal_out$biasAll[iter] = pbias(HBVoutput$Qg, HBVoutput$historicQinmm)
						# calculating monthly biases
					HBVoutput$month = month(HBVoutput$Date)
					HBVoutputMonth = subset(HBVoutput, month == 1)
					cal_out$mnthBias_1[iter] = pbias(HBVoutputMonth$Qg, HBVoutputMonth$historicQinmm)
#					cal_out$mnthSumAbsBias[iter] = abs(cal_out$mnthBias_1[iter])
					for(thisMonth in 2:12)	{
						HBVoutputMonth = subset(HBVoutput, month == thisMonth)
						cal_out[iter, paste0('mnthBias_', thisMonth)] = pbias(HBVoutputMonth$Qg, HBVoutputMonth$historicQinmm)
					}
					mnthBiasCols = which(names(cal_out) == 'mnthBias_1'):which(names(cal_out) == 'mnthBias_12')
					cal_out$mnthSumAbsBias[iter] = apply(abs(cal_out[iter, mnthBiasCols]), 1, mean)
					
					fwrite(cal_out, paste0(dataOut_location, "calibration_", basinName, ".csv"), append=FALSE)
				}
			}
		}
	}
}





##########################################################################################################
	## Running the Model for Seasonal Forecasts 
seasonalStreamflowForecast_f = function(
	basinName = 'inster_basin_or_outlet_name',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1,							# 1 for FNF from cal.gov,
	waterYearStart = as.Date('yyyy-mm-dd'),
	forecastDate = as.Date('yyyy-mm-dd'),
	gageLonLat = c(1,1),
	biasCorrection = TRUE,
	uploadToGCS = TRUE)
	{

	if(file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		# read in the necessary libraries
		library(data.table)
		library(lubridate)
		library(zoo)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in basin specific data that has been previously loaded
		seas5ClimateInput = readRDS(paste0(dataOut_location, 'SEAS5_', basinName, '.RData')) # reads in a list with each [[i]] being output from a model
		era5ClimateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

	
			# reading in calibrated parameterizations and selecting a subsample based on minimum monthly bias
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = head(calibratedVarsAll[order(abs(calibratedVarsAll$mnthSumAbsBias)), ], 40)
		rm(calibratedVarsAll)

			# reading in and reformatting streamflow input 
		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
				# removing negative streamflow and zero streamflow for later log transformation
			if(any(historicStreamflow$historicQinOriginalUnits < 0))	{historicStreamflow$historicQinOriginalUnits[historicStreamflow$historicQinOriginalUnits < 0] = NA}
			basinArea = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
				# converting cfs to af / day
			historicStreamflow$AcFtPDay_inflow =  as.numeric(historicStreamflow$historicQinOriginalUnits)  * (60*60*24) / 43559.9
				# interpolating nas if necessary)
			if(any(is.na(historicStreamflow$AcFtPDay_inflow)))	{
				historicStreamflow$AcFtPDay_inflow[1] = historicStreamflow$AcFtPDay_inflow[which(!is.na(historicStreamflow$AcFtPDay_inflow))][1]
				historicStreamflow$AcFtPDay_inflow[nrow(historicStreamflow)] = last(historicStreamflow$AcFtPDay_inflow[which(!is.na(historicStreamflow$AcFtPDay_inflow))])
				historicStreamflow$AcFtPDay_inflow = na.approx(historicStreamflow$AcFtPDay_inflow)
			}
		
			
				# routine for calculating climatolology of streamflow
			historicStreamflow$year = year(historicStreamflow$Date)
			historicStreamflow$doy = yday(historicStreamflow$Date)
			startWYs = which(historicStreamflow$doy == yday(waterYearStart))
			
			climatologyDF = data.frame(dayOfWY = 1:365)
			for(thisWY in startWYs[1:(length(startWYs) - 1)])	{
				histSubset = historicStreamflow[thisWY:(thisWY+364)]
				if(histSubset$Date[365] - histSubset$Date[1] <= 366)	{	# ensuring we aren't using years with large data gaps
					# backfilling for missing data
					if(any(is.na(histSubset$AcFtPDay_inflow)))	{
						histSubset$AcFtPDay_inflow[is.na(histSubset$AcFtPDay_inflow)] = mean(histSubset$AcFtPDay_inflow, na.rm=TRUE)
					}
					if(histSubset$AcFtPDay_inflow[1] < 1)	{histSubset$AcFtPDay_inflow[1] = 1}
					climatologyDF = cbind(climatologyDF, cumsum(histSubset$AcFtPDay_inflow))
				}
			}
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}
		
			# building a record of historical streamflow for comparison

	
			# identifying which time periods will use historic data before replacing with forecast data
		lastHistData = last(era5ClimateInput$Date)
		seas5Rows = (1 + which(as.character(seas5ClimateInput[[1]]$Date) == as.character(lastHistData))):length(seas5ClimateInput[[1]]$Date)

		allForecastsOutput = data.frame(Date = c(era5ClimateInput$Date, seas5ClimateInput[[1]]$Date[seas5Rows]))

		iter = 0
		for(numModels in 1:length(seas5ClimateInput))	{
			for(numCalibs in 1:nrow(calibratedVars))	{	
				iter = iter + 1
				climateInput = rbind(era5ClimateInput,seas5ClimateInput[[numModels]][seas5Rows,])
				

				allHBVoutput = runHBV_f(
					climateInput = climateInput,
					sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
					tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
					tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
					fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
					fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
					fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
					fc   = calibratedVars$fc[numCalibs], # field capacity 
					lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
					beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
					k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
					k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
					k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
					uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
					perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d
				if(dataSource == 1)	{
					allHBVoutput$projectedQinOriginalUnits = allHBVoutput$Qg * flowUnitConversion * areaUnitConversion * basinArea
					allHBVoutput$AcFtPDay_proj = allHBVoutput$projectedQinOriginalUnits  * (60*60*24) / 43559.9
				}	else	{ 
					return("we need to figure out how to read in and normalize this streamflow data")
				}
				
				#optional bias correction by month
				if(biasCorrection)	{
					allHBVoutput$month =  month(allHBVoutput$Date)
					for(ii in 1:12)	{
						biasCol = which(names(calibratedVars) == paste0('mnthBias_', ii))
						debiasVal = calibratedVars[numCalibs, ..biasCol]
						allHBVoutput$AcFtPDay_proj[allHBVoutput$month == ii] = allHBVoutput$AcFtPDay_proj[allHBVoutput$month == ii] *  as.numeric(100 - debiasVal) / 100
					}
				}
	
	
					# compiling all forecasts into a list
				allForecastsOutput = cbind(allForecastsOutput, allHBVoutput$AcFtPDay_proj)
	


			}
		}
	}	else	{print("Ya gotta calibrate the model first you big ole dummy!")}
		####			prevRecord = which(climateAndStreamflowOutput$Date < waterYearStart)
	forecastThisWY = subset(allForecastsOutput, Date >= waterYearStart)
	if(last(historicStreamflow$Date) >= waterYearStart)	{
		whichRecentStreamflow = which(historicStreamflow$Date >= waterYearStart)
		forecastThisWY[1:length(whichRecentStreamflow), -1] = historicStreamflow$AcFtPDay_inflow[whichRecentStreamflow]
	}
	forecastThisWY[ , -1] = apply(forecastThisWY[ , -1], 2, cumsum)
		
	
		# cleaning hist data by interpolation
	if(is.na(historicStreamflow$AcFtPDay_inflow[1]) | is.na(last(historicStreamflow$AcFtPDay_inflow)))	{
		historicStreamflow$AcFtPDay_inflow[1] = historicStreamflow$AcFtPDay_inflow[!is.na(historicStreamflow$AcFtPDay_inflow)][1]
		historicStreamflow$AcFtPDay_inflow[nrow(historicStreamflow)] = last(historicStreamflow$AcFtPDay_inflow[!is.na(historicStreamflow$AcFtPDay_inflow)])
	}
	historicStreamflow$strmIntrp = na.approx(historicStreamflow$AcFtPDay_inflow)
	currWY = subset(historicStreamflow, Date >= waterYearStart)
	prevWY = subset(historicStreamflow, Date < waterYearStart & Date >= waterYearStart %m-% years(1))
	prev2WY = subset(historicStreamflow, Date < waterYearStart %m-% years(1) & Date >= waterYearStart %m-% years(2))
	prev3WY = subset(historicStreamflow, Date < waterYearStart %m-% years(2) & Date >= waterYearStart %m-% years(3))

	lnForecast = nrow(forecastThisWY)
	lnNAs = 365 - lnForecast
	currWYNAs = max(365 - nrow(currWY), 0)
	prevWYNAs = max(365 - nrow(prevWY), 0)
	forecastOutput = data.frame(Date = seq(as.Date(waterYearStart), by='days', length.out=365),
		Reservoir = basinName,
		Units = NA,
		Lat = gageLonLat[2], Lon = gageLonLat[1],
		Pred_Q05 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.05), rep(NA, lnNAs)),
		Pred_Q25 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.25), rep(NA, lnNAs)),
		Pred_Q50 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.50), rep(NA, lnNAs)),
		Pred_Q75 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.75), rep(NA, lnNAs)),
		Pred_Q95 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.95), rep(NA, lnNAs)),
		Clim_Q05 = apply(climatologyDF[,-1], 1, quantile, probs=0.05),
		Clim_Q25 = apply(climatologyDF[,-1], 1, quantile, probs=0.25),
		Clim_Q50 = apply(climatologyDF[,-1], 1, quantile, probs=0.50),
		Clim_Q75 = apply(climatologyDF[,-1], 1, quantile, probs=0.75),
		Clim_Q95 = apply(climatologyDF[,-1], 1, quantile, probs=0.95),
		WaterYear_Curr = c(cumsum(currWY$strmIntrp), rep(NA, currWYNAs)),
		WaterYear_1YrAgo = c(cumsum(prevWY$strmIntrp), rep(NA, prevWYNAs)),
		WaterYear_2YrAgo = cumsum(prev2WY$strmIntrp[-366]),
		WaterYear_3YrAgo = cumsum(prev3WY$strmIntrp[-366]))


	if(dataSource == 1)	{
		forecastOutput$Units = 'Acre Feet'
	} else {
		return("we need to figure out how to read in and normalize this streamflow data")
	}


	# create the folder for storing outputs
	dataOut_fileLoc = paste0(dataOut_location, 'strmflwForecastsFor_', forecastDate)
	if(!file.exists(dataOut_fileLoc))	{
			dir.create(file.path(paste0(dataOut_fileLoc, '\\')))
	}



		#saving output
	fwrite(forecastOutput, paste0(dataOut_fileLoc, "forecastStrmCumSum_", basinName, '_', forecastDate, ".csv"))
	
		#cumulative streamflow forecast figure
	fstOfMnths = forecastOutput$Date[which(mday(forecastOutput$Date) == 1)]
	mnthNms = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
	png(paste0(dataOut_fileLoc, '\\', 'streamflowForecast.png'), width=1920, height=960)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
	plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(0,max(forecastOutput$Clim_Q95)*1.05),
		type='l', lwd=1, col='white', xaxt = 'n', #log='y',
		main='', ylab='Cumulative Streamflow (Acre Feet)', xlab='',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(1, at = fstOfMnths,col.lab='#1A232F', col.axis='#666D74', 
		labels = mnthNms)
	abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
	lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
		col='#FDB600', lwd=4, lty = 2)
	lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
		col='#FDB600', lwd=4, lty = 3)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
		col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
		col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
		col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
	forecastCompletes = forecastOutput[!is.na(forecastOutput$Pred_Q50),]
	polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q05, rev(forecastCompletes$Pred_Q95)),
		col=adjustcolor('#0098B2', alpha.f=0.1), border=NA)
	polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q25, rev(forecastCompletes$Pred_Q75)),
		col=adjustcolor('#0098B2', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Pred_Q50, 
		col='#0098B2', lwd=5)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.95,
		paste0('Forecast for ', basinName),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.87,
		paste0('Issued on ', forecastDate),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	dev.off()
	
			# log of cumulative streamflow forecast figure 
	png(paste0(dataOut_fileLoc, '\\', 'logStreamflowForecast.png'), width=1920, height=960)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
	plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(min(forecastOutput$Clim_Q05),max(forecastOutput$Clim_Q95)*1.05),
		type='l', lwd=1, col='white', log='y', yaxt='n', xaxt='n',
		main='', ylab='Cumulative Streamflow (Acre Feet)', xlab='',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
	axis(2, at = c(1, 10, 100, 1000, 10000,100000, 1000000, 10000000), col.lab='#1A232F', col.axis='#666D74', 
		labels = c('1', '10', '100', '1,000', '10,000', '100,000', '1,000,000', '10,000,000'))
	axis(1, at = fstOfMnths, col.lab='#1A232F', col.axis='#666D74', 
		labels = mnthNms)
	lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
		col='#FDB600', lwd=4, lty = 2)
	lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
		col='#FDB600', lwd=4, lty = 3)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
		col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
		col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
		col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
	polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q05, rev(forecastCompletes$Pred_Q95)),
		col=adjustcolor('#0098B2', alpha.f=0.1), border=NA)
	polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q25, rev(forecastCompletes$Pred_Q75)),
		col=adjustcolor('#0098B2', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Pred_Q50, 
		col='#0098B2', lwd=5)
	text(x=forecastOutput$Date[20], y=max(forecastOutput$Clim_Q95)*0.95,
		paste0('Forecast for ', basinName, ' Issued on ', forecastDate),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	dev.off()

	######################################################################################
	## prevWY figures
	if(Sys.Date() > waterYearStart)	{
		forecastOutput$WaterYear_3YrAgo = forecastOutput$WaterYear_2YrAgo
		forecastOutput$WaterYear_2YrAgo = forecastOutput$WaterYear_1YrAgo
		forecastOutput$WaterYear_1YrAgo = forecastOutput$WaterYear_Curr
	}
		#cumulative streamflow WY record figure
	fstOfMnths = forecastOutput$Date[which(mday(forecastOutput$Date) == 1)]
	mnthNms = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
	png(paste0(dataOut_fileLoc, '\\','streamflowWYcurrent.png'), width=1920, height=960)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
	plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(0,max(forecastOutput$Clim_Q95)*1.05),
		type='l', lwd=1, col='white', xaxt = 'n', #log='y',
		main='', ylab='Cumulative Streamflow (Acre Feet)', xlab='',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(1, at = fstOfMnths,col.lab='#1A232F', col.axis='#666D74', 
		labels = mnthNms)
	abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
	abline(v=forecastOutput$Date[length(which(!is.na(forecastOutput$WaterYear_1YrAgo)))], col='black')
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
		col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
		col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
		col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
	lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
		col='#0098B2', lwd=5, lty = 1)
	lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
		col='#FDB600', lwd=4, lty = 2)
	lines(forecastOutput$Date, forecastOutput$WaterYear_3YrAgo, 
		col='#FDB600', lwd=4, lty = 3)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.95,
		paste0('Water Year Record for ', basinName),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.87,
		paste0('Issued on ', forecastDate),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	lastRecData = last(which(!is.na(forecastOutput$WaterYear_1YrAgo)))
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.79,
		paste0('Currently ',
			signif(forecastOutput$WaterYear_1YrAgo[lastRecData] / forecastOutput$Clim_Q50[lastRecData], 2) * 100,
			'% of Average'),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	dev.off()
	
			# log of cumulative streamflow WY record figure 
	png(paste0(dataOut_fileLoc, '\\', 'logStreamflowWYcurrent.png'), width=1920, height=960)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
	plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(min(forecastOutput$Clim_Q05),max(forecastOutput$Clim_Q95)*1.05),
		type='l', lwd=1, col='white', log='y', yaxt='n', xaxt='n',
		main='', ylab='Cumulative Streamflow (Acre Feet)', xlab='',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at = c(1, 10, 100, 1000, 10000,100000, 1000000, 10000000), col.lab='#1A232F', col.axis='#666D74', 
		labels = c('1', '10', '100', '1,000', '10,000', '100,000', '1,000,000', '10,000,000'))
	axis(1, at = fstOfMnths, col.lab='#1A232F', col.axis='#666D74', 
		labels = mnthNms)
	abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
	abline(v=forecastOutput$Date[length(which(!is.na(forecastOutput$WaterYear_1YrAgo)))], col='black')
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
		col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
		col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
		col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
	lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
		col='#0098B2', lwd=5, lty = 1)
	lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
		col='#FDB600', lwd=4, lty = 2)
	lines(forecastOutput$Date, forecastOutput$WaterYear_3YrAgo, 
		col='#FDB600', lwd=4, lty = 3)
	text(x=forecastOutput$Date[20], y=min(forecastOutput$Clim_Q05)*1.05,
		paste0('Water Year Record for ', basinName, ';   Issued on ', forecastDate,
			';   Currently ',
				signif(forecastOutput$WaterYear_1YrAgo[lastRecData] / forecastOutput$Clim_Q50[lastRecData], 2) * 100,
				'% of Average'),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	dev.off()
		
	# uploading figs to GCStorage
	if(uploadToGCS == TRUE)	{
		setwd(paste0(dataOut_fileLoc))
		command = paste0('gsutil -m cp ', dataOut_fileLoc, '/*.png gs://seasonal_surface_water-looker_img_hosting/', basinName)
		system(command)
	}
}






##########################################################################################################
	## Running the Model for Seasonal Forecasts for storage
seasonalStorageForecast_f = function(
	basinName = 'inster_basin_or_outlet_name',
	historicStreamflowFileLoc = 'https://someplace.gov',
	historicReservoirFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1,							# 1 for FNF from cal.gov,
	waterYearStart = as.Date('yyyy-mm-dd'),
	forecastDate = as.Date('yyyy-mm-dd'),
	gageLonLat = c(1,1),
	biasCorrection = TRUE,
	uploadToGCS = TRUE,
	incldStorage = TRUE)
	{

	if(file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		# read in the necessary libraries
		library(data.table)
		library(lubridate)
		library(zoo)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in basin specific data that has been previously loaded
		seas5ClimateInput = readRDS(paste0(dataOut_location, 'SEAS5_', basinName, '.RData')) # reads in a list with each [[i]] being output from a model
		era5ClimateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

	
			# reading in calibrated parameterizations and selecting a subsample based on minimum monthly bias
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = head(calibratedVarsAll[order(abs(calibratedVarsAll$mnthSumAbsBias)), ], 40)
		rm(calibratedVarsAll)

			# reading in and reformatting streamflow and storage input s
		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		if(dataSource == 1)	{
			inflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			stor = read.csv(paste0(historicReservoirFileLoc))
			#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

			stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
			stor$stor = as.numeric(stor$VALUE)
			inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
			inflow$inflow = inflow$VALUE
			allDat = inflow[stor[c('Date','stor')], on = 'Date']
			whichFirstNoNA = which(!is.na(allDat$stor) & !is.na(allDat$inflow))[1]
			allDat = allDat[-c(1:(whichFirstNoNA-1)), ]
			allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
				# resLoss = Outflow + AET + Withdrawals = dS - Inflow
			allDat$storIntrp = na.approx(allDat$stor)
			allDat$infIntrp = allDat$AcFtPDay_inflow
			allDat$infIntrp[is.na(allDat$infIntrp)] = mean(allDat$infIntrp, na.rm=TRUE)
			
			allDat$resLoss = c(diff(allDat$storIntrp) - allDat$infIntrp[-nrow(allDat)], NA) * -1
			if(any(allDat$storIntrp <= 1)) {allDat$storIntrp[allDat$storIntrp <= 1] = 1}	# using a power law model so eliminating 0s
			if(any(allDat$resLoss <= 1)) {allDat$resLoss[allDat$resLoss <= 1] = 1}	# using a power law model so eliminating 0s
			if(any(allDat$infIntrp <= 1)) {allDat$infIntrp[allDat$infIntrp <= 1] = 1}	# using a power law model so eliminating 0s
				# debiasing resLoss
			allDat$resLoss = allDat$resLoss * (1 - (mean(allDat$resLoss, na.rm=TRUE) - mean(allDat$infIntrp)) / mean(allDat$resLoss, na.rm=TRUE))

				# basic models for relating storage and inflows to outflows
#			plot(allDat$infIntrp, allDat$resLoss)
			infModelLg = lm(allDat$resLoss ~ log(allDat$infIntrp))
			infModelLn = lm(allDat$resLoss ~ allDat$infIntrp)
#			infModelEx = lm(log(allDat$resLoss) ~ allDat$infIntrp)
#			plot(allDat$storIntrp, allDat$resLoss)
			storModelLg = lm(allDat$resLoss ~ log(allDat$storIntrp))
			storModelLn = lm(allDat$resLoss ~ allDat$storIntrp)
#			storModelEx = lm(log(allDat$resLoss) ~ allDat$storIntrp)
#			plot(yday(allDat$Date), allDat$resLoss)
			allDat$doy = yday(allDat$Date)
			allDat$doyEstResLoss = mean(allDat$resLoss, na.rm = TRUE)
			for(ydays in 1:366)	{
				allDat$doyEstResLoss[allDat$doy == ydays] = mean(allDat$resLoss[which(allDat$doy == ydays)], na.rm=TRUE)
			}
			ydayModelLg = lm(allDat$resLoss ~ log(allDat$doyEstResLoss))
			ydayModelLn = lm(allDat$resLoss ~ allDat$doyEstResLoss)
#			ydayModelEx = lm(log(allDat$resLoss) ~ allDat$doyEstResLoss)
				# normalizing for reweighting models for use in projections
			totR2 = summary(infModelLn)$adj + summary(infModelLg)$adj + summary(storModelLg)$adj + summary(storModelLn)$adj +summary(ydayModelLg)$adj + summary(ydayModelLn)$adj
			infModelLg_wt = summary(infModelLg)$adj / totR2
			infModelLn_wt = summary(infModelLn)$adj / totR2
#			infModelLn_wt = summary(infModelEx)$adj / totR2
			storModelLg_wt = summary(storModelLg)$adj / totR2
			storModelLn_wt = summary(storModelLn)$adj / totR2
#			storModelEx_wt = summary(storModelEx)$adj / totR2
			ydayModelLg_wt = summary(ydayModelLg)$adj / totR2
			ydayModelLn_wt = summary(ydayModelLn)$adj / totR2
#			ydayModelEx_wt = summary(ydayModelEx)$adj / totR2
		
				# calculating factor to convert mm to acre feet for storage calculations
			basinAreaInKm = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			acrePerSqKm = 247.105
			mmPerFoot = 304.8
			mmToAcrFt = basinAreaInKm * acrePerSqKm / mmPerFoot
			
				# building a record of historical storage for comparison
				
			allDat$yday = yday(allDat$Date)
			climatologyDF = data.frame(dayOfWY = 1:365, Clim_Q05=NA, Clim_Q25=NA, Clim_Q50=NA, Clim_Q75=NA, Clim_Q95=NA)
			ydayWY = yday(seq(waterYearStart, waterYearStart + 364, by = 1))
			wyYday = 0
			for(thisDay in ydayWY)	{
				wyYday = wyYday + 1
				histSubset = subset(allDat, yday == thisDay)
				climatologyDF[wyYday, -1] = quantile(histSubset$storIntrp, c(0.05, 0.25, 0.50, 0.75, 0.95))
			}
	
				
		} else	{return("need to figure out how to handle these data")}

	
			# identifying which time periods will use historic data before replacing with forecast data
		lastHistData = last(era5ClimateInput$Date)
		seas5Rows = (1 + which(as.character(seas5ClimateInput[[1]]$Date) == as.character(lastHistData))):length(seas5ClimateInput[[1]]$Date)

		allForecastsOutput = data.frame(Date = c(era5ClimateInput$Date, seas5ClimateInput[[1]]$Date[seas5Rows]))

		iter = 0
		for(numModels in 1:length(seas5ClimateInput))	{
			for(numCalibs in 1:nrow(calibratedVars))	{	
				iter = iter + 1
				climateInput = rbind(era5ClimateInput,seas5ClimateInput[[numModels]][seas5Rows,])
				

				allHBVoutput = runHBV_f(
					climateInput = climateInput,
					sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
					tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
					tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
					fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
					fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
					fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
					fc   = calibratedVars$fc[numCalibs], # field capacity 
					lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
					beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
					k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
					k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
					k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
					uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
					perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d
	



					# merging historic streamflow record onto climate inputs data
				climateAndStreamflowOutput = allDat[allHBVoutput, on='Date']#[projDates, ]  
				for(ydays in 1:366)	{
						climateAndStreamflowOutput$doyEstResLoss[climateAndStreamflowOutput$doy == ydays] = mean(climateAndStreamflowOutput$resLoss[which(climateAndStreamflowOutput$doy == ydays)], na.rm=TRUE)
					}
				climateAndStreamflowOutput$inflowProj = climateAndStreamflowOutput$Qg * mmToAcrFt
				climateAndStreamflowOutput$month = month(climateAndStreamflowOutput$Date)
							#identifying the last non-na storage value
				lastHistStor = last(which(!is.na(climateAndStreamflowOutput$stor)))
				modelOut = climateAndStreamflowOutput[-c(1:(lastHistStor - 1)), ]
	
			
				#optional bias correction by month
				if(biasCorrection)	{
					for(ii in 1:12)	{
						biasCol = which(names(calibratedVars) == paste0('mnthBias_', ii))
						debiasVal = calibratedVars[numCalibs, ..biasCol]
						modelOut$inflowProj[modelOut$month == ii] = modelOut$inflowProj[modelOut$month == ii] *  as.numeric(100 - debiasVal) / 100
					}
				}


					# estimating reservoir losses as a function of inflows and storage
				if(any(modelOut$inflowProj <=1))	{modelOut$inflowProj[modelOut$inflowProj <= 1] = 1}	# climate debiasing generates some negative precip, which generates negative streamflow, which must be eliminated for log transformation
				resLossEstInfYday = 
					infModelLg_wt * (infModelLg$coef[1] + log(modelOut$inflowProj) * infModelLg$coef[2]) +
					infModelLn_wt * (infModelLn$coef[1] + modelOut$inflowProj * infModelLn$coef[2]) +
					ydayModelLg_wt * (ydayModelLg$coef[1] + log(modelOut$doyEstResLoss) * ydayModelLg$coef[2]) +
					ydayModelLn_wt * (ydayModelLn$coef[1] + modelOut$doyEstResLoss * ydayModelLn$coef[2])
				if(any(resLossEstInfYday < 0))	{resLossEstInfYday[resLossEstInfYday < 0] = 0}
					
				storEstInf = as.numeric(climateAndStreamflowOutput[(projDates[1] - 1), 'storIntrp']) + cumsum(modelOut$inflowProj) - cumsum(resLossEstInfYday / (1 - (storModelLg_wt + storModelLn_wt) / 1))
				if(any(storEstInf <= 1))	{storEstInf[storEstInf <= 1] = 1}

					# smoothing initial estimates of storage before running final function of reservoir loss ~ storage
				storSmoother = 0
				while(storSmoother < 10)	{
					resLossEstStor = 
						storModelLg_wt * (storModelLg$coef[1] + log(storEstInf) * storModelLg$coef[2]) +
						storModelLn_wt * (storModelLn$coef[1] + storEstInf * storModelLn$coef[2])
					resLossEstTot = resLossEstInfYday + resLossEstStor
					if(any(resLossEstTot < 0))	{resLossEstTot[resLossEstTot < 0] = 0}
					storEstInfStor = as.numeric(climateAndStreamflowOutput[(projDates[1] - 1), 'storIntrp']) + cumsum(modelOut$inflowProj) - 
						cumsum(resLossEstTot)
					if(any(storEstInfStor < 1))	{storEstInfStor[storEstInfStor < 1] = 1}
					storSmoother = storSmoother + 1
				}
					
						
				modelOut$storEst = storEstInfStor
	
					# compiling all forecasts into a list
				allForecastsOutput = cbind(allForecastsOutput, c(climateAndStreamflowOutput$stor[1:(lastHistStor - 1)] , modelOut$storEst))

			}
		}
	}	else	{print("Ya gotta calibrate the model first you big ole dummy!")}
	




	
		# combining and printing storage forecasts
	forecastThisWY = subset(allForecastsOutput, Date >= waterYearStart)
	currWY = subset(allDat, Date >= waterYearStart)
	prevWY = subset(allDat, Date < waterYearStart & Date >= waterYearStart %m-% years(1))
	prev2WY = subset(allDat, Date < waterYearStart %m-% years(1) & Date >= waterYearStart %m-% years(2))
	prev3WY = subset(allDat, Date < waterYearStart %m-% years(2) & Date >= waterYearStart %m-% years(3))

	lnForecast = nrow(forecastThisWY)
	lnNAs = max(365 - lnForecast, 0)
	currWYNAs = max(365 - nrow(currWY), 0)
	prevWYNAs = max(365 - nrow(prevWY), 0)

	forecastOutput = data.frame(Date = seq(as.Date(waterYearStart), by='days', length.out=365),
		Reservoir = basinName,
		Units = NA,
		Lat = gageLonLat[2], Lon = gageLonLat[1],
		Pred_Q05 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.05), rep(NA, lnNAs)),
		Pred_Q25 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.25), rep(NA, lnNAs)),
		Pred_Q50 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.50), rep(NA, lnNAs)),
		Pred_Q75 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.75), rep(NA, lnNAs)),
		Pred_Q95 = c(apply(forecastThisWY[,-1], 1, quantile, probs=0.95), rep(NA, lnNAs)),
		Clim_Q05 = climatologyDF$Clim_Q05,
		Clim_Q25 = climatologyDF$Clim_Q25,
		Clim_Q50 = climatologyDF$Clim_Q50,
		Clim_Q75 = climatologyDF$Clim_Q75,
		Clim_Q95 = climatologyDF$Clim_Q95,
		WaterYear_Curr = c(currWY$storIntrp, rep(NA, currWYNAs)),
		WaterYear_1YrAgo = c(prevWY$storIntrp, rep(NA, prevWYNAs)),
		WaterYear_2YrAgo = prev2WY$storIntrp[-366],
		WaterYear_3YrAgo = prev3WY$storIntrp[-366])
		

	if(dataSource == 1)	{
		forecastOutput$Units = 'Acre Feet'
	} else {
		return("we need to figure out how to read in and normalize this streamflow data")
	}


	# create the folder for storing outputs
	dataOut_fileLoc = paste0(dataOut_location, 'storageForecastsFor_', forecastDate)
	if(!file.exists(dataOut_fileLoc))	{
			dir.create(file.path(paste0(dataOut_fileLoc, '\\')))
	}



		#saving output
	fwrite(forecastOutput, paste0(dataOut_fileLoc, "forecastStorage_", basinName, '_', forecastDate, ".csv"))
	
	if(incldStorage)	{
			#storage forecast figure
		fstOfMnths = forecastOutput$Date[which(mday(forecastOutput$Date) == 1)]
		mnthNms = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
		png(paste0(dataOut_fileLoc, '\\', 'storageForecast.png'), width=1920, height=960)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
		plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(0,max(forecastOutput$Clim_Q95)*1.05),
			type='l', lwd=1, col='white', xaxt = 'n', #log='y',
			main='', ylab='Storage (Acre Feet)', xlab='',
			col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
			family='A')
		axis(1, at = fstOfMnths,col.lab='#1A232F', col.axis='#666D74', 
			labels = mnthNms)
		abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
		lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
			col='#FDB600', lwd=4, lty = 2)
		lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
			col='#FDB600', lwd=4, lty = 3)
		polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
			col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
		polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
			col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
		lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
			col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
		forecastCompletes = forecastOutput[!is.na(forecastOutput$Pred_Q50),]
		polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q05, rev(forecastCompletes$Pred_Q95)),
			col=adjustcolor('#0098B2', alpha.f=0.1), border=NA)
		polygon(x=c(forecastCompletes$Date, rev(forecastCompletes$Date)), y=c(forecastCompletes$Pred_Q25, rev(forecastCompletes$Pred_Q75)),
			col=adjustcolor('#0098B2', alpha.f=0.2), border=NA)
		lines(forecastOutput$Date, forecastOutput$Pred_Q50, 
			col='#0098B2', lwd=5)
		text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.95,
			paste0('Forecast for ', basinName),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
		text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.87,
			paste0('Issued on ', forecastDate),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
		dev.off()
	}	else	{
		png(paste0(dataOut_fileLoc, '\\', 'storageForecast.png'), width=1920, height=960)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
		plot(0:1, 0:1, xaxt = 'n', col='white',
			main='', ylab='', xlab='',axes=FALSE)
		text(x=0.05, y=0.55,
			paste0('Storage Forecast Currently Unavailable'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=5*1.3)
		text(x=0.16, y=0.45,
			paste0('modeled reservoir releases too uncertain for reliable forecasting'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)	
			
		dev.off()
	}
	
	######################################################################################
	## prevWY figures
	if(Sys.Date() > waterYearStart)	{
		forecastOutput$WaterYear_3YrAgo = forecastOutput$WaterYear_2YrAgo
		forecastOutput$WaterYear_2YrAgo = forecastOutput$WaterYear_1YrAgo
		forecastOutput$WaterYear_1YrAgo = forecastOutput$WaterYear_Curr
	}
		#storage WY record figure
	fstOfMnths = forecastOutput$Date[which(mday(forecastOutput$Date) == 1)]
	mnthNms = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
	png(paste0(dataOut_fileLoc, '\\', 'storageWYcurrent.png'), width=1920, height=960)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
	plot(forecastOutput$Date, forecastOutput$Clim_Q95, ylim = c(0,max(forecastOutput$Clim_Q95)*1.05),
		type='l', lwd=1, col='white', xaxt = 'n', #log='y',
		main='', ylab='Storage(Acre Feet)', xlab='',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(1, at = fstOfMnths,col.lab='#1A232F', col.axis='#666D74', 
		labels = mnthNms)
	abline(v=fstOfMnths, lwd=1, col=adjustcolor('#666D74', alpha.f=0.1))
	abline(v=forecastOutput$Date[length(which(!is.na(forecastOutput$WaterYear_1YrAgo)))], col='black')
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q95, rev(forecastOutput$Clim_Q05)),
		col=adjustcolor('#666D74', alpha.f=0.1), border=NA)
	polygon(x=c(forecastOutput$Date, rev(forecastOutput$Date)), y=c(forecastOutput$Clim_Q25, rev(forecastOutput$Clim_Q75)),
		col=adjustcolor('#666D74', alpha.f=0.2), border=NA)
	lines(forecastOutput$Date, forecastOutput$Clim_Q50, 
		col=adjustcolor('#666D74', alpha.f=0.2), lwd=2.5)
	lines(forecastOutput$Date, forecastOutput$WaterYear_1YrAgo, 
		col='#0098B2', lwd=5, lty = 1)
	lines(forecastOutput$Date, forecastOutput$WaterYear_2YrAgo, 
		col='#FDB600', lwd=4, lty = 2)
	lines(forecastOutput$Date, forecastOutput$WaterYear_3YrAgo, 
		col='#FDB600', lwd=4, lty = 3)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.95,
		paste0('Water Year Record for ', basinName),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.87,
		paste0('Issued on ', forecastDate),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	lastRecData = last(which(!is.na(forecastOutput$WaterYear_1YrAgo)))
	text(x=forecastOutput$Date[1], y=max(forecastOutput$Clim_Q95)*0.79,
		paste0('Currently ',
			signif(forecastOutput$WaterYear_1YrAgo[lastRecData] / forecastOutput$Clim_Q50[lastRecData], 2) * 100,
			'% of Average'),
		adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.3)
	dev.off()
	
	
	# uploading figs to GCStorage
	if(uploadToGCS == TRUE)	{
		setwd(paste0(dataOut_fileLoc))
		command = paste0('gsutil cp ', dataOut_fileLoc, '/*.png gs://seasonal_surface_water-looker_img_hosting/', basinName)
		system(command)
	}
}































###################################################################################################################
	## Validation and plot generation
validationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1)							# 1 for FNF from cal.gov,
	{

	if(file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		# read in the necessary libraries
		library(data.table)
		library(scales)
		library(lubridate)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in basin specific data that has been previously loaded
		climateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

			# reading in calibrated parameterizations and selecting a subsample based on minimum monthly bias
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = head(calibratedVarsAll[order(abs(calibratedVarsAll$mnthSumAbsBias)), ], 40)
		#calibratedVars = tail(calibratedVarsAll[order(abs(calibratedVarsAll$kgeAll)), ], 40)
		rm(calibratedVarsAll)

			# reading in and reformatting streamflow input 
		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
				# removing negative streamflow
			if(any(historicStreamflow$historicQinOriginalUnits < 0))	{historicStreamflow$historicQinOriginalUnits[historicStreamflow$historicQinOriginalUnits < 0] = NA}
			basinArea = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}

			# merging climate and streamflow data
		climateAndStreamflowInput = historicStreamflow[climateInput, on='Date']
		allStreamflows = data.frame(Date = climateAndStreamflowInput$Date, GageQ = climateAndStreamflowInput$historicQinmm)
		iter = 0
		for(numCalibs in 1:nrow(calibratedVars))	{
			iter = iter + 1

			allHBVoutput = runHBV_f(
				climateInput = climateInput,
				sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
				tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
				tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
				fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
				fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
				fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
				fc   = calibratedVars$fc[numCalibs], # field capacity 
				lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
				beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
				k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
				k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
				k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
				uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
				perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d

			# merging historic streamflow record onto climate inputs data
			allStreamflows[ , paste0('model_', numCalibs)] = allHBVoutput$Qg
		}

		allStreamflows[, 'medianModelQ'] = apply(allStreamflows[ , -c(1:2)], 1, median)
		


			# period of record overall model performance
		png(paste0(dataOut_location, basinName, '_hydrograph_historicalValidation.png'), width=1920, height=960)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=2*c(5,5,2,2), mgp=2*c(3,1.3,0), font.lab=2, bty='l', cex.lab=2*1.8, cex.axis=2*1.4, cex.main=2*1.8, col='#1A232F')
		nonaAllQs = subset(allStreamflows[-c(1:365), ], !is.na(GageQ))
		plot(nonaAllQs$Date[-c(1:365)], nonaAllQs$GageQ[-c(1:365)], type='l', lwd=2, col='white',
			main='', ylab='Streamflow (mm)', xlab='Date', ylim = c(0, max(nonaAllQs$GageQ[-c(1:365)])*1.2),
			col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
			family='A')
		lines(smooth.spline(nonaAllQs$Date, nonaAllQs$GageQ, spar=0.000001), lwd=2*2, col='grey10', lty=1, type='l')
		for(numRows in 1:nrow(calibratedVars))	{
			lines(nonaAllQs$Date, nonaAllQs[, paste0('model_', numRows)], col=alpha('#0098B2', 0.05), lwd=2*0.5)
		}
		lines(nonaAllQs$Date, nonaAllQs$medianModelQ, col='#196CE1', lwd=2*1.2)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (9/10)),
			paste0('KGE = ', round(median(calibratedVars$kgeAll), 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (8.3/10)),
			paste0('NSE = ', round(median(calibratedVars$nseAll), 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (7.6/10)),
			paste0('r  = ', round(summary(lm(nonaAllQs$GageQ ~ nonaAllQs$medianModelQ))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (7.8/10)),
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (6.9/10)),
			paste0('MAE = ', round(median(calibratedVars$maeAll), 2), ' mm'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=nonaAllQs$Date[370], y=max(nonaAllQs$GageQ[-c(1:365)] * (6.2/10)),
			paste0('% Bias = ', round(median(calibratedVars$biasAll), 2), '%'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=max(nonaAllQs$Date), y=max(nonaAllQs$GageQ[-c(1:365)] * (8.3/10)),
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=2*1.9)
		dev.off()




			# annual average overall performance
		png(paste0(dataOut_location, basinName, '_annualDepth_historicalValidation.png'), width=960, height=960)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
		nonaAllQs$year = year(nonaAllQs$Date + 92)	# converting to water year by subtracting 92 days 
		yearQs = nonaAllQs[1,-1]
		thisRow = 0
		for(numYears in unique(nonaAllQs$year))	{
			thisRow = thisRow + 1
			thisYear = which(nonaAllQs$year == numYears)
			yearQs[thisRow, -ncol(yearQs)] = apply(nonaAllQs[thisYear, -c(1,ncol(nonaAllQs))], 2, sum)
		}
		yearQs$year = unique(nonaAllQs$year)
		plot(yearQs$GageQ, yearQs$medianModelQ, pch=1, lwd=2, col='white',
			main='', ylab='Annual Modeled Streamflow (mm)', xlab='Annual Historic Streamflow (mm)',
			col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
			family='A',
			ylim=c(0,max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2)), xlim=c(0,max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2)), cex=2)
		abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
		abline(lm(yearQs$medianModelQ ~ yearQs$GageQ), lwd=2, col=alpha('#F06000', 1), cex=2)
		for(numRows in 1:nrow(calibratedVars))	{
			points(yearQs$GageQ, yearQs[, paste0('model_', numRows)], pch=15, lwd=1.5,  col=alpha('#0098B2', 0.05), cex=2.5)
		}
		points(yearQs$GageQ, yearQs$medianModelQ, pch=22, lwd=2.5, bg=alpha('#0098B2', 0.05), col='#196CE1', cex=3)
		text(x=0, y=max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2) * (9.2/10),
			paste0('r  = ', round(summary(lm(yearQs$medianModelQ ~ yearQs$GageQ))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=0, y=max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2) * (9.4/10),
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
		text(x=0, y=max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2) * (8.4/10),
			paste0('slope  = ', round(summary(lm(yearQs$medianModelQ ~ yearQs$GageQ))$coef[2], 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=max(c(yearQs$GageQ,yearQs$medianModelQ)*1.2) , y=0,
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
		dev.off()


			# monthly average overall performance
		nonaAllQs$month = month(nonaAllQs$Date)
		monthsQs = nonaAllQs[1:12,-1]
		monthsQs = monthsQs[, -which(names(monthsQs) == 'year')]
		for(numMonths in 1:12)	{
			theseMonths = which(nonaAllQs$month == numMonths)
			monthsQs[numMonths, 1:(ncol(monthsQs)-1)] = apply(nonaAllQs[theseMonths, 2:(ncol(nonaAllQs) - 2)], 2, mean)*30
			monthsQs$month[numMonths] = numMonths
		}
		png(paste0(dataOut_location, basinName, '_monthlyDepth_historicalValidation.png'), width=960, height=960)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
		plot(monthsQs$GageQ, monthsQs$medianModelQ, pch=1, lwd=2, col='white',
			main='', ylab='Monthly Modeled Streamflow (mm)', xlab='Monthly Historic Streamflow (mm)',
			col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
			family='A',
			ylim=c(0,max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2)), xlim=c(0,max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2)), cex=2)
		abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
		abline(lm(monthsQs$medianModelQ ~ monthsQs$GageQ), lwd=2, col=alpha('#F06000', 1), cex=2)
		for(numRows in 1:nrow(calibratedVars))	{
			points(monthsQs$GageQ, monthsQs[, paste0('model_', numRows)], pch=15, lwd=1.5,  col=alpha('#0098B2', 0.05), cex=2.5)
		}
		points(monthsQs$GageQ, monthsQs$medianModelQ, pch=22, lwd=2.5, bg=alpha('#0098B2', 0.05), col='#196CE1', cex=3)
		text(x=0, y=max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2) * (9.2/10),
			paste0('r  = ', round(summary(lm(monthsQs$medianModelQ ~ monthsQs$GageQ))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=0, y=max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2) * (9.4/10),
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
		text(x=0, y=max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2) * (8.4/10),
			paste0('slope  = ', round(summary(lm(monthsQs$medianModelQ ~ monthsQs$GageQ))$coef[2], 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=max(c(monthsQs$GageQ,monthsQs$medianModelQ)*1.2) , y=0,
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.6*1.9)
		dev.off()



			# monthly average overall performance by %
		nonaAllQs$month = month(nonaAllQs$Date)
		monthsQsBias = nonaAllQs[1:12,-1]
		monthsQsBias = monthsQsBias[, -which(names(monthsQsBias) == 'year')]
		for(numMonths in 1:12)	{
			theseMonths = which(nonaAllQs$month == numMonths)
			monthsQsBias[numMonths, 2:(ncol(monthsQsBias) - 1)] = 100 * 
				(((apply(nonaAllQs[theseMonths, 3:(ncol(nonaAllQs) - 2)], 2, sum)) - sum(nonaAllQs[theseMonths, 2])) / sum(nonaAllQs[theseMonths, 2]))
			monthsQsBias$month[numMonths] = numMonths
		}
		monthlyBP = data.frame(modelError = NA, month=NA)
		for(numRows in 1:nrow(calibratedVars))	{
			newBP = data.frame(modelError = monthsQsBias[, paste0('model_', numRows)], month = 1:12)
			monthlyBP = rbind(monthlyBP, newBP)
		}
		png(paste0(dataOut_location, basinName, '_monthlyPercentage_historicalValidation.png'), width = 1200, height = 720)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
		boxplot(modelError ~ month, data=monthlyBP, pch=1, lwd=1, col='#0098B2', border='#666D74', cex=1.5, 
			main='', ylab='% Error', xlab='Month', xaxt='n',
			col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
			family='A')
		axis(1, at=1:12, labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
			col.lab='#1A232F', col.axis='#666D74')
		abline(h=0, lwd=2, col='#1A232F', cex=2)
		text(x=1, y=180,
			paste0(basinName),
			adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.9)
		dev.off()

	historicalValidation = list()
	historicalValidation[[1]] = nonaAllQs
	historicalValidation[[2]] = monthsQsBias
	saveRDS(historicalValidation, paste0(dataOut_location, 'historicalValidationData_', basinName, '.RData'))
	}	else	{print("Ya gotta calibrate the model first you big ole dummy!")}
}





projectionValidationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1,
	biasCorrection = TRUE)							# 1 for FNF from cal.gov,
	{
	
	if(!file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		print('gotta calibrate the model first you big ole dummy!!')
	}	else	{

		if(!file.exists(paste0(dataOut_location, basinName, '_projectedOutputFigures')))	{
			dir.create(file.path(paste0(dataOut_location, basinName, '_projectedOutputFigures')))
		}

			# read in the necessary libraries
		library(data.table)
		library(lubridate)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# read in climate data
		seas5ClimateInputFileNames = list.files(paste0(dataOut_location, 'seas5MultiOutput'))
		era5ClimateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

			# reading in calibrated parameterizations and selecting a subsample based on minimum monthly bias
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = head(calibratedVarsAll[order(abs(calibratedVarsAll$mnthSumAbsBias)), ], 40)
		#calibratedVars = tail(calibratedVarsAll[order(abs(calibratedVarsAll$kgeAll)), ], 40)
		rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
				# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
				# removing negative streamflow
			if(any(historicStreamflow$historicQinOriginalUnits < 0))	{historicStreamflow$historicQinOriginalUnits[historicStreamflow$historicQinOriginalUnits < 0] = NA}
			basinArea = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}
		
			# building a record of historical streamflow for comparison
		historicStreamflow$year = year(historicStreamflow$Date)
		historicStreamflow$month = month(historicStreamflow$Date)
		climatologyDF = data.frame(month = NA, year = NA, totQ = NA)
		climIter = 0
		for(thisYear in unique(historicStreamflow$year))	{
			for(thisMonth in 1:12)	{
				climIter = climIter + 1
				histSubset = subset(historicStreamflow, year == thisYear & month == thisMonth)
				climatologyDF[climIter, ] = c(thisMonth, thisYear, mean(histSubset$historicQinmm, na.rm=TRUE) * nrow(histSubset))
			}
		}
		

		save_iter = 0
		validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predQ = NA, actQ = NA)

		for(thisSeas5 in seas5ClimateInputFileNames)	{
			seas5ClimateInput = readRDS(paste0(dataOut_location, 'seas5MultiOutput\\', thisSeas5))

			# removing monthly forecasts that are not a full month
			partialMonths = FALSE
			if(last(mday(seas5ClimateInput[[1]]$Date + 1)) != 1) 	{
				partialMonths = TRUE
				rmSEAS5DateRm = length(last(which(mday(seas5ClimateInput[[1]]$Date) == 1)):nrow(seas5ClimateInput[[1]]$Date)) - 1
			}
			
			lastHistData = which(era5ClimateInput$Date == seas5ClimateInput[[1]]$Date[1])
			allForecastsOutput = list()
			iter = 0
			for(numModels in 1:length(seas5ClimateInput))	{
				for(numCalibs in 1:nrow(calibratedVars))	{	
					iter = iter + 1
					
					climateInput = rbind(era5ClimateInput[1:lastHistData,], seas5ClimateInput[[numModels]])
					if(partialMonths) 	{climateInput = climateInput[-c((nrow(climateInput)-rmSEAS5DateRm):nrow(climateInput)), ]}
					
					projDates = (1 + lastHistData):nrow(climateInput)



					allHBVoutput = runHBV_f(
						climateInput = climateInput,
						sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
						tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
						tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
						fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
						fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
						fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
						fc   = calibratedVars$fc[numCalibs], # field capacity 
						lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
						beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
						k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
						k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
						k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
						uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
						perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d

					# merging historic streamflow record onto climate inputs data
					climateAndStreamflowOutput = historicStreamflow[allHBVoutput, on='Date'][projDates, ]  
					climateAndStreamflowOutput$month = month(climateAndStreamflowOutput$Date)

					#optional bias correction by month
					if(biasCorrection)	{
						for(ii in 1:12)	{
							biasCol = which(names(calibratedVars) == paste0('mnthBias_', ii))
							debiasVal = calibratedVars[numCalibs, ..biasCol]
							climateAndStreamflowOutput$Qg[climateAndStreamflowOutput$month == ii] = climateAndStreamflowOutput$Qg[climateAndStreamflowOutput$month == ii] *  as.numeric(100 - debiasVal) / 100
						}
					}
					

					climateAndStreamflowOutput$year = year(climateAndStreamflowOutput$Date)
					monthsOut = 0
					predCumsum = 0
					actCumsum = 0 
					
					for(thisMonth in unique(climateAndStreamflowOutput$month))	{
						monthsOut = monthsOut + 1
						outputSubset = subset(climateAndStreamflowOutput, month == thisMonth)
						numDays = nrow(outputSubset)
						predCumsum = predCumsum + mean(outputSubset$Qg) * numDays
						actCumsum = actCumsum + mean(outputSubset$historicQinmm, na.rm=TRUE) * numDays
						validDF = rbind(validDF,
							c(numModels, numCalibs, outputSubset$Date[1] + 14, monthsOut,
								predCumsum, actCumsum))
					}
					
					
						
					
					
					if(nrow(validDF) > 20000)	{
						save_iter = save_iter + 1
						fwrite(validDF, paste0(dataOut_location, "temp_out_", save_iter, ".csv"))
						validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predQ = NA, actQ = NA)
					}
				}
			}
		}

			# combining saved iterations and returning the dataframe with formatted climatology
		save_iter = save_iter + 1
		fwrite(validDF, paste0(dataOut_location, "temp_out_", save_iter, ".csv"))
			
		validDF = fread(paste0(dataOut_location, "temp_out_", 1, ".csv"))
		for(files in 2:save_iter)	{
			validDF = rbind(validDF, fread(paste0(dataOut_location, "temp_out_", files, ".csv")))
		}
	
		validDF = validDF[complete.cases(validDF), ]
		validDF$Date = validDF$days + as.Date("1970-01-01")
		validDF$startDate = validDF$Date %m-% (months(validDF$monthsOut - 1)) 
		
		startDates = unique(validDF$startDate) 
			
			# datatable for assessing error / performance over entire period for which we have data
		sumModPerf = data.frame(month = NA, leadTime = NA, predTerc = NA, actTerc = NA, 
			predMAE25 = NA, climMAE25 = NA, predMAE50 = NA, climMAE50 = NA, predMAE75 = NA, climMAE75 = NA)
		
		for(thisStartDate in startDates)	{
			validSubset = subset(validDF, startDate == thisStartDate)
			
				# parsing the climatology df for plotting
			theseMonths = unique(month(validSubset$Date))
			numMonths = length(theseMonths)
			climDF = data.frame(month = NA, cumsumQ = NA, plotMonth = NA)
			for(startYear in climatologyDF$year[1]:(max(climatologyDF$year)-1))	{
				climSubset = subset(climatologyDF, year >= startYear)
				climStart = which(climSubset$month == month(validSubset$startDate)[1])[1]
				
				climDF = rbind(climDF,
					data.table(month = theseMonths, cumsumQ = cumsum(climSubset$totQ[climStart:(climStart+numMonths-1)]), plotMonth = c(1:7)))
			}	
		
			
			
			png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\streamflow_', validSubset$Date[1] - 14, '.png'), width = 1200, height = 720)
			windowsFonts(A = windowsFont("Roboto"))
			par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
			#validSubset$plotDate = rep(1:7, nrow(validSubset)/7) - .2
			boxplot(validSubset$predQ ~ validSubset$Date, pch='.', lwd=1, 
				#col=colorRampPalette(c('#F06000', '#B91863'))(length(unique(validDF$monthsOut))),
				ylim = c(0, max(validDF$predQ, na.rm=TRUE)),
				at = c(1:7) - 0.2,
				boxwex = 0.4,
				col = '#A9BF2C',
				border='#666D74', cex=1.5, 
				main='', ylab='Cumulative Streamflow (mm)', xlab='Month', xaxt='n',
				col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
				family='A')
			#climDF$plotDate =climDF$plotMonth + 0.2
			boxplot(cumsumQ ~ plotMonth, data = climDF, pch='.',
				col = '#FDB600', 
				at = c(1:7) + 0.2,
				boxwex = 0.4, 
				border='#666D74', cex=1.5, 
				main='', ylab='Cumulative Streamflow (mm)', xlab='Month', xaxt='n',
				add=TRUE)
			lines(1:length(unique(validSubset$Date)), unique(validSubset$actQ), type='l', col='#196CE1', lwd=1)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actQ), pch = '---', col='#F2F3F3', lwd=2, cex=9)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actQ), pch = '---', col='#196CE1', lwd=2, cex=7)
			axis(1, at=1:length(unique(validSubset$Date)),
				labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')[month(unique(validSubset$Date))],
				col.lab='#1A232F', col.axis='#666D74')
			text(x=0.5, y=max(validDF$predQ, na.rm=TRUE) * 0.95,
				paste0(basinName),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predQ, na.rm=TRUE)* 0.87,
				paste0('Forecast Beginning on ', validSubset$Date[1] - 14),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predQ, na.rm=TRUE) * 0.75,
				paste0('CAi Forecast'),
				adj = c(0,0), font=2, col='#A9BF2C', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predQ, na.rm=TRUE) * 0.7,
				paste0('Climatology'),
				adj = c(0,0), font=2, col='#FDB600', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predQ, na.rm=TRUE) * 0.65,
				paste0('Actual'),
				adj = c(0,0), font=2, col='#196CE1', family='A', cex=1.6*1.3)
			dev.off()

			monthsLead = 0
			for(thisLeadTime in unique(validSubset$Date)) {
				monthsLead = monthsLead + 1
				validSubsetLd = subset(validSubset, Date == thisLeadTime)
				climDFLd = subset(climDF, plotMonth == monthsLead)
				climQuant = quantile(climDFLd$cumsumQ, c(0.333, 0.5, 0.667), na.rm=TRUE) 
				predQuant = quantile(validSubsetLd$predQ, c(0.333, 0.5, 0.667)) 
				
				predTerc = 1
				if(predQuant[2] > climQuant[1]) {predTerc = 2}
				if(predQuant[2] > climQuant[3])	{predTerc = 3}
				
				actVal = validSubsetLd$actQ[1]
				actTerc = 1
				if(actVal > climQuant[1]) 		{actTerc = 2}
				if(actVal > climQuant[3])		{actTerc = 3}

				sumModPerf = rbind(sumModPerf, c(
					as.numeric(climDFLd$month[1]), monthsLead, predTerc, actTerc, 
					as.numeric(predQuant[1] - actVal), as.numeric(climQuant[1] - actVal),
					as.numeric(predQuant[2] - actVal), as.numeric(climQuant[2] - actVal),
					as.numeric(predQuant[3] - actVal), as.numeric(climQuant[3] - actVal)))
			}
		}
	}
	fwrite(sumModPerf, paste0(dataOut_location, basinName, '_summaryOfModelPerformance.csv'))

#	sumModPerf = data.frame(month = NA, leadTime = NA, predTerc = NA, actTerc = NA, 
#		predMAE25 = NA, climMAE25 = NA, predMAE50 = NA, climMAE50 = NA, predMAE75 = NA, climMAE75 = NA)

	highTercHeatMap = matrix(nrow=7, ncol=12)
	lowTercHeatMap  = matrix(nrow=7, ncol=12)
	medTercHeatMap  = matrix(nrow=7, ncol=12)
	tercPredHighBarplot = matrix(nrow=3, ncol=7)
	tercPredMedBarplot  = matrix(nrow=3, ncol=7)
	tercPredLowBarplot  = matrix(nrow=3, ncol=7)
	tercActHighBarplot = matrix(nrow=3, ncol=7)
	tercActMedBarplot  = matrix(nrow=3, ncol=7)
	tercActLowBarplot  = matrix(nrow=3, ncol=7)
	absMae25HeatMap = matrix(nrow=7, ncol=12)
	absMae50HeatMap = matrix(nrow=7, ncol=12)
	absMae75HeatMap = matrix(nrow=7, ncol=12)
	relMae25HeatMap = matrix(nrow=7, ncol=12)
	relMae50HeatMap = matrix(nrow=7, ncol=12)
	relMae75HeatMap = matrix(nrow=7, ncol=12)
	
	
	for(jj in 1:7)	{
		highPredTerc  = subset(sumModPerf, predTerc == 3 & leadTime == jj)
		medPredTerc   = subset(sumModPerf, predTerc == 2 & leadTime == jj)
		lowPredTerc   = subset(sumModPerf, predTerc == 1 & leadTime == jj)
		
		highActTerc = subset(sumModPerf, actTerc == 3 & leadTime == jj)
		medActTerc  = subset(sumModPerf, actTerc == 2 & leadTime == jj)
		lowActTerc  = subset(sumModPerf, actTerc == 1 & leadTime == jj)
		
		tercPredHighBarplot[1:3,jj] = c(
			length(which(highPredTerc$actTerc == 3)) / nrow(highPredTerc), 
			length(which(highPredTerc$actTerc == 2)) / nrow(highPredTerc),
			length(which(highPredTerc$actTerc == 1)) / nrow(highPredTerc))
		
		tercPredMedBarplot[1:3,jj] = c(
			length(which(medPredTerc$actTerc == 3)) / nrow(medPredTerc), 
			length(which(medPredTerc$actTerc == 2)) / nrow(medPredTerc),
			length(which(medPredTerc$actTerc == 1)) / nrow(medPredTerc))

		tercPredLowBarplot[1:3,jj] = c(
			length(which(lowPredTerc$actTerc == 3)) / nrow(lowPredTerc), 
			length(which(lowPredTerc$actTerc == 2)) / nrow(lowPredTerc),
			length(which(lowPredTerc$actTerc == 1)) / nrow(lowPredTerc))

		tercActHighBarplot[1:3,jj] = c(
			length(which(highActTerc$predTerc == 3)) / nrow(highActTerc), 
			length(which(highActTerc$predTerc == 2)) / nrow(highActTerc),
			length(which(highActTerc$predTerc == 1)) / nrow(highActTerc))
		
		tercActMedBarplot[1:3,jj] = c(
			length(which(medActTerc$predTerc == 3)) / nrow(medActTerc), 
			length(which(medActTerc$predTerc == 2)) / nrow(medActTerc),
			length(which(medActTerc$predTerc == 1)) / nrow(medActTerc))

		tercActLowBarplot[1:3,jj] = c(
			length(which(lowActTerc$predTerc == 3)) / nrow(lowActTerc), 
			length(which(lowActTerc$predTerc == 2)) / nrow(lowActTerc),
			length(which(lowActTerc$predTerc == 1)) / nrow(lowActTerc))


		for(ii in 1:12)	{
			subModPerf = subset(sumModPerf, month == ii & leadTime == jj)
			
			absMae25HeatMap[jj,ii] = median(subModPerf$predMAE25)
			absMae50HeatMap[jj,ii] = median(subModPerf$predMAE50)
			absMae75HeatMap[jj,ii] = median(subModPerf$predMAE75)

			relMae25HeatMap[jj,ii] = median((abs(subModPerf$predMAE25) - abs(subModPerf$climMAE25)) / abs(subModPerf$climMAE25)) * 100
			relMae50HeatMap[jj,ii] = median((abs(subModPerf$predMAE50) - abs(subModPerf$climMAE50)) / abs(subModPerf$climMAE50)) * 100 
			relMae75HeatMap[jj,ii] = median((abs(subModPerf$predMAE75) - abs(subModPerf$climMAE75)) / abs(subModPerf$climMAE75)) * 100

			subHighModPerf = subset(subModPerf, actTerc == 3)
			subMedModPerf = subset(subModPerf, actTerc == 2)
			subLowModPerf = subset(subModPerf, actTerc == 1)
			highTercHeatMap[jj,ii] = length(which(subHighModPerf$predTerc == 3)) / nrow(subHighModPerf)
			medTercHeatMap[jj,ii] = length(which(subMedModPerf$predTerc == 2)) / nrow(subMedModPerf)
			lowTercHeatMap[jj,ii] = length(which(subLowModPerf$predTerc == 1)) / nrow(subLowModPerf)
		}
	}
	row.names(tercActHighBarplot) = c('high', 'med', 'low')
	row.names(tercActMedBarplot) = c('high', 'med', 'low')
	row.names(tercActLowBarplot) = c('high', 'med', 'low')
	row.names(tercPredHighBarplot) = c('high', 'med', 'low')
	row.names(tercPredMedBarplot) = c('high', 'med', 'low')
	row.names(tercPredLowBarplot) = c('high', 'med', 'low')
	
	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11) # colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = seq(-100, 100, length.out = 12)
	
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\SummaryPlot_relMAE50.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	if(any(relMae50HeatMap > 100))	{relMae50HeatMap[which(relMae50HeatMap > 100)] = 99.9}
	image(t(relMae50HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae50HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()
			
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\SummaryPlot_relMAE25.png'), width = 720, height = 720)
	if(any(relMae25HeatMap > 100))	{relMae25HeatMap[which(relMae25HeatMap > 100)] = 99.9}
	image(t(relMae25HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae25HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\SummaryPlot_relMAE75.png'), width = 720, height = 720)
	if(any(relMae75HeatMap > 100))	{relMae75HeatMap[which(relMae75HeatMap > 100)] = 99.9}
	image(t(relMae75HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae75HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileHeatMap_high.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(highTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% High Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(highTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileHeatMap_med.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(medTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% Med Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(medTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileHeatMap_low.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(lowTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% Low Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(lowTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()
	
	
		# terciles based on prediction values
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercilePredBarPlot_high.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredHighBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='High Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
	
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercilePredBarPlot_med.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredMedBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Med Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()

	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercilePredBarPlot_low.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredLowBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Low Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()


		# terciles based on actual values
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileActBarPlot_high.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActHighBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='High Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
	
	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileActBarPlot_med.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActMedBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Med Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()

	png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\tercileActBarPlot_low.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActLowBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Low Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
}




















projectedStorageValidationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	historicStreamflowFileLoc = 'https://someplace.gov',
	historicReservoirFileLoc = 'insert_file_loc_here.com',
	dataOut_location = 'save_file_location',
	dataSource = 1,
	biasCorrection = TRUE)							# 1 for FNF from cal.gov,
	{
	
		# Confirming the hydro model and storage model have both been calibrated
	if(!file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		return('gotta calibrate the model first you big ole dummy!!')
	}	else	{

		if(!file.exists(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures')))	{
			dir.create(file.path(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures')))
		}

		############ loading libraries
		library(data.table)
		library(lubridate)
		library(CAST)
		library(caret)
		library(sf)
		library(zoo)

			# reading in data on reservoir inflows and storage
		if(dataSource == 1)	{
			inflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			stor = read.csv(paste0(historicReservoirFileLoc))
			#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

			stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
			stor$stor = as.numeric(stor$VALUE)
			inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
			inflow$inflow = inflow$VALUE
			allDat = inflow[stor[c('Date','stor')], on = 'Date']
			whichFirstNoNA = which(!is.na(allDat$stor) & !is.na(allDat$inflow))[1]
			allDat = allDat[-c(1:(whichFirstNoNA-1)), ]
			allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
				# resLoss = Outflow + AET + Withdrawals = dS - Inflow
			allDat$storIntrp = na.approx(allDat$stor)
			allDat$infIntrp = allDat$AcFtPDay_inflow
			allDat$infIntrp = na.approx(allDat$AcFtPDay_inflow)
			
			allDat$resLoss = c(diff(allDat$storIntrp) - allDat$infIntrp[-nrow(allDat)], NA) * -1
			if(any(allDat$storIntrp <= 1)) {allDat$storIntrp[allDat$storIntrp <= 1] = 1}	# using a power law model so eliminating 0s
			if(any(allDat$resLoss <= 1)) {allDat$resLoss[allDat$resLoss <= 1] = 1}	# using a power law model so eliminating 0s
			if(any(allDat$infIntrp <= 1)) {allDat$infIntrp[allDat$infIntrp <= 1] = 1}	# using a power law model so eliminating 0s
				# debiasing resLoss
			allDat$resLoss = allDat$resLoss * (1 - (mean(allDat$resLoss, na.rm=TRUE) - mean(allDat$infIntrp)) / mean(allDat$resLoss, na.rm=TRUE))

				# basic models for relating storage and inflows to outflows
#			plot(allDat$infIntrp, allDat$resLoss)
			infModelLg = lm(allDat$resLoss ~ log(allDat$infIntrp))
			infModelLn = lm(allDat$resLoss ~ allDat$infIntrp)
#			infModelEx = lm(log(allDat$resLoss) ~ allDat$infIntrp)
#			plot(allDat$storIntrp, allDat$resLoss)
			storModelLg = lm(allDat$resLoss ~ log(allDat$storIntrp))
			storModelLn = lm(allDat$resLoss ~ allDat$storIntrp)
#			storModelEx = lm(log(allDat$resLoss) ~ allDat$storIntrp)
#			plot(yday(allDat$Date), allDat$resLoss)
			allDat$doy = yday(allDat$Date)
			allDat$doyEstResLoss = mean(allDat$resLoss, na.rm = TRUE)
			for(ydays in 1:366)	{
				allDat$doyEstResLoss[allDat$doy == ydays] = mean(allDat$resLoss[which(allDat$doy == ydays)], na.rm=TRUE)
			}
			ydayModelLg = lm(allDat$resLoss ~ log(allDat$doyEstResLoss))
			ydayModelLn = lm(allDat$resLoss ~ allDat$doyEstResLoss)
#			ydayModelEx = lm(log(allDat$resLoss) ~ allDat$doyEstResLoss)
				# normalizing for reweighting models for use in projections
			totR2 = summary(infModelLn)$adj + summary(infModelLg)$adj + summary(storModelLg)$adj + summary(storModelLn)$adj +summary(ydayModelLg)$adj + summary(ydayModelLn)$adj
			infModelLg_wt = summary(infModelLg)$adj / totR2
			infModelLn_wt = summary(infModelLn)$adj / totR2
#			infModelLn_wt = summary(infModelEx)$adj / totR2
			storModelLg_wt = summary(storModelLg)$adj / totR2
			storModelLn_wt = summary(storModelLn)$adj / totR2
#			storModelEx_wt = summary(storModelEx)$adj / totR2
			ydayModelLg_wt = summary(ydayModelLg)$adj / totR2
			ydayModelLn_wt = summary(ydayModelLn)$adj / totR2
#			ydayModelEx_wt = summary(ydayModelEx)$adj / totR2
		
		
				# calculating factor to convert mm to acre feet for storage calculations
			basinAreaInKm = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			acrePerSqKm = 247.105
			mmPerFoot = 304.8
			mmToAcrFt = basinAreaInKm * acrePerSqKm / mmPerFoot

		} else	{return("need to figure out how to handle these data")}

			# read in climate data
			# read in climate data
		seas5ClimateInputFileNames = list.files(paste0(dataOut_location, 'seas5MultiOutput'))
		era5ClimateInput = as.data.table(readRDS(paste0(dataOut_location, 'ERA5_', basinName, '.RData')))

			# reading in calibrated parameterizations and selecting a subsample based on minimum monthly bias
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = head(calibratedVarsAll[order(abs(calibratedVarsAll$mnthSumAbsBias)), ], 40)
		#calibratedVars = tail(calibratedVarsAll[order(abs(calibratedVarsAll$kgeAll)), ], 40)
		rm(calibratedVarsAll)

			
			# building a record of historical streamflow for comparison
		allDat$year = year(allDat$Date)
		allDat$month = month(allDat$Date)
		climatologyDF = data.frame(month = NA, year = NA, Stor = NA)
		climIter = 0
		for(thisYear in unique(allDat$year)[-length(unique(allDat$year))])	{
			for(thisMonth in unique(month(subset(allDat, year == thisYear)$Date)))	{
				climIter = climIter + 1
				histSubset = subset(allDat, year == thisYear & month == thisMonth)
				climatologyDF[climIter, ] = c(thisMonth, thisYear, last(histSubset$storIntrp))
			}
		}
	

			
			

		save_iter = 0
		validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predStor = NA, actStor = NA)

		for(thisSeas5 in seas5ClimateInputFileNames)	{
			seas5ClimateInput = readRDS(paste0(dataOut_location, 'seas5MultiOutput\\', thisSeas5))

			# removing monthly forecasts that are not a full month
			partialMonths = FALSE
			if(last(mday(seas5ClimateInput[[1]]$Date + 1)) != 1) 	{
				partialMonths = TRUE
				rmSEAS5DateRm = length(last(which(mday(seas5ClimateInput[[1]]$Date) == 1)):nrow(seas5ClimateInput[[1]]$Date)) - 1
			}
				
			lastHistData = which(era5ClimateInput$Date == seas5ClimateInput[[1]]$Date[1])
			allForecastsOutput = list()
			iter = 0
			for(numModels in 1:length(seas5ClimateInput))	{
				for(numCalibs in 1:nrow(calibratedVars))	{	
					iter = iter + 1
						
					climateInput = rbind(era5ClimateInput[1:lastHistData,], seas5ClimateInput[[numModels]])
					if(partialMonths) 	{climateInput = climateInput[-c((nrow(climateInput)-rmSEAS5DateRm):nrow(climateInput)), ]}
						
					projDates = (1 + lastHistData):nrow(climateInput)



					allHBVoutput = runHBV_f(
						climateInput = climateInput,
						sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
						tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
						tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
						fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
						fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
						fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
						fc   = calibratedVars$fc[numCalibs], # field capacity 
						lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
						beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
						k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
						k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
						k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
						uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
						perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d

					# merging historic streamflow record onto climate inputs data
					climateAndStreamflowOutput = allDat[allHBVoutput, on='Date']#[projDates, ]  

					for(ydays in 1:366)	{
						climateAndStreamflowOutput$doyEstResLoss[climateAndStreamflowOutput$doy == ydays] = mean(climateAndStreamflowOutput$resLoss[which(climateAndStreamflowOutput$doy == ydays)], na.rm=TRUE)
					}
					modelOut = climateAndStreamflowOutput[projDates, ]
					modelOut$inflowProj = modelOut$Qg * mmToAcrFt
	
					#optional bias correction by month
					if(biasCorrection)	{
						for(ii in 1:12)	{
							biasCol = which(names(calibratedVars) == paste0('mnthBias_', ii))
							debiasVal = calibratedVars[numCalibs, ..biasCol]
							modelOut$inflowProj[modelOut$month == ii] = modelOut$inflowProj[modelOut$month == ii] *  as.numeric(100 - debiasVal) / 100
						}
					}



						# estimating reservoir losses as a function of inflows and storage
					if(any(modelOut$inflowProj <=1))	{modelOut$inflowProj[modelOut$inflowProj <= 1] = 1}	# climate debiasing generates some negative precip, which generates negative streamflow, which must be eliminated for log transformation
					resLossEstInfYday = 
						infModelLg_wt * (infModelLg$coef[1] + log(modelOut$inflowProj) * infModelLg$coef[2]) +
						infModelLn_wt * (infModelLn$coef[1] + modelOut$inflowProj * infModelLn$coef[2]) +
						ydayModelLg_wt * (ydayModelLg$coef[1] + log(modelOut$doyEstResLoss) * ydayModelLg$coef[2]) +
						ydayModelLn_wt * (ydayModelLn$coef[1] + modelOut$doyEstResLoss * ydayModelLn$coef[2])
					if(any(resLossEstInfYday < 0))	{resLossEstInfYday[resLossEstInfYday < 0] = 0}
					
					storEstInf = as.numeric(climateAndStreamflowOutput[(projDates[1] - 1), 'storIntrp']) + cumsum(modelOut$inflowProj) - cumsum(resLossEstInfYday / (1 - (storModelLg_wt + storModelLn_wt) / 1))
					if(any(storEstInf <= 1))	{storEstInf[storEstInf <= 1] = 1}

						# smoothing initial estimates of storage before running final function of reservoir loss ~ storage
					storSmoother = 0
					while(storSmoother < 10)	{
						resLossEstStor = 
							storModelLg_wt * (storModelLg$coef[1] + log(storEstInf) * storModelLg$coef[2]) +
							storModelLn_wt * (storModelLn$coef[1] + storEstInf * storModelLn$coef[2])
						resLossEstTot = resLossEstInfYday + resLossEstStor
						if(any(resLossEstTot < 0))	{resLossEstTot[resLossEstTot < 0] = 0}
						storEstInfStor = as.numeric(climateAndStreamflowOutput[(projDates[1] - 1), 'storIntrp']) + cumsum(modelOut$inflowProj) - 
							cumsum(resLossEstTot)
						if(any(storEstInfStor < 1))	{storEstInfStor[storEstInfStor < 1] = 1}
						storSmoother = storSmoother + 1
					}
					
						
					modelOut$storEst = storEstInfStor
					
					

					monthsOut = 0
					for(thisMonth in unique(modelOut$month))	{
						monthsOut = monthsOut + 1
						outputSubset = subset(modelOut, month == thisMonth)
						numDays = nrow(outputSubset)
						validDF = rbind(validDF,
							c(numModels, numCalibs, outputSubset$Date[1] + 14, monthsOut,
								mean(outputSubset$storEst), mean(outputSubset$storIntrp)))
					}
					
						
					if(nrow(validDF) > 20000)	{
						save_iter = save_iter + 1
						fwrite(validDF, paste0(dataOut_location, "temp_out_stor", save_iter, ".csv"))
						validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predStor = NA, actStor = NA)
					}
				}
			}
		}	

			# combining saved iterations and returning the dataframe with formatted climatology
		save_iter = save_iter + 1
		fwrite(validDF, paste0(dataOut_location, "temp_out_stor", save_iter, ".csv"))
			
		validDF = fread(paste0(dataOut_location, "temp_out_stor", 1, ".csv"))
		for(files in 2:save_iter)	{
			validDF = rbind(validDF, fread(paste0(dataOut_location, "temp_out_stor", files, ".csv")))
		}
	
		validDF = validDF[complete.cases(validDF), ]
		validDF$Date = validDF$days + as.Date("1970-01-01")
		validDF$startDate = validDF$Date %m-% (months(validDF$monthsOut - 1)) 
		
		startDates = unique(validDF$startDate)
			
			# datatable for assessing error / performance over entire period for which we have data
		sumModPerf = data.frame(month = NA, leadTime = NA, predTerc = NA, actTerc = NA, 
			predMAE25 = NA, climMAE25 = NA, predMAE50 = NA, climMAE50 = NA, predMAE75 = NA, climMAE75 = NA)
		
		for(thisStartDate in startDates)	{
			validSubset = subset(validDF, startDate == thisStartDate)
			
				# parsing the climatology df for plotting
			theseMonths = unique(month(validSubset$Date))
			numMonths = length(theseMonths)
			
			climDF = subset(climatologyDF, month %in% theseMonths)
			
			
			png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\streamflow_', validSubset$Date[1] - 14, '.png'), width = 1200, height = 720)
			windowsFonts(A = windowsFont("Roboto"))
			par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
			#validSubset$plotDate = rep(1:7, nrow(validSubset)/7) - .2
			boxplot(validSubset$predStor ~ validSubset$Date, pch='.', lwd=1, 
				#col=colorRampPalette(c('#F06000', '#B91863'))(length(unique(validDF$monthsOut))),
				ylim = c(0, max(validDF$predStor, na.rm=TRUE)),
				at = c(1:7) - 0.2,
				boxwex = 0.4,
				col = '#A9BF2C',
				border='#666D74', cex=1.5, 
				main='', ylab='Storage (Acre Feet)', xlab='Month', xaxt='n',
				col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
				family='A')
			#climDF$plotDate =climDF$plotMonth + 0.2
			boxplot(Stor ~ month, data = climDF, pch='.',
				col = '#FDB600', 
				at = c(1:7)[order(theseMonths)] + 0.2,
				boxwex = 0.4, 
				border='#666D74', cex=1.5, 
				main='', ylab='Storage (Acre Feet)', xlab='Month', xaxt='n',
				add=TRUE)
			lines(1:length(unique(validSubset$Date)), unique(validSubset$actStor), type='l', col='#196CE1', lwd=1)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actStor), pch = '---', col='#F2F3F3', lwd=2, cex=9)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actStor), pch = '---', col='#196CE1', lwd=2, cex=7)
			axis(1, at=1:length(unique(validSubset$Date)),
				labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')[month(unique(validSubset$Date))],
				col.lab='#1A232F', col.axis='#666D74')
			text(x=0.5, y=max(validDF$predStor, na.rm=TRUE) * 0.95,
				paste0(basinName),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predStor, na.rm=TRUE)* 0.87,
				paste0('Forecast Beginning on ', validSubset$Date[1] - 14),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predStor, na.rm=TRUE) * 0.75,
				paste0('CAi Forecast'),
				adj = c(0,0), font=2, col='#A9BF2C', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predStor, na.rm=TRUE) * 0.7,
				paste0('Climatology'),
				adj = c(0,0), font=2, col='#FDB600', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predStor, na.rm=TRUE) * 0.65,
				paste0('Actual'),
				adj = c(0,0), font=2, col='#196CE1', family='A', cex=1.6*1.3)
			dev.off()

			monthsLead = 0
			for(thisLeadTime in unique(validSubset$Date)) {
				monthsLead = monthsLead + 1
				validSubsetLd = subset(validSubset, Date == thisLeadTime)
				climDFLd = subset(climatologyDF, month == month(validSubsetLd$Date[1]))
				climQuant = quantile(climDFLd$Stor, c(0.333, 0.5, 0.667), na.rm=TRUE) 
				predQuant = quantile(validSubsetLd$predStor, c(0.333, 0.5, 0.667)) 
				
				predTerc = 1
				if(predQuant[2] > climQuant[1]) {predTerc = 2}
				if(predQuant[2] > climQuant[3])	{predTerc = 3}
				
				actVal = validSubsetLd$actStor[1]
				actTerc = 1
				if(actVal > climQuant[1]) 		{actTerc = 2}
				if(actVal > climQuant[3])		{actTerc = 3}

				sumModPerf = rbind(sumModPerf, c(
					as.numeric(month(validSubsetLd$Date[1])), monthsLead, predTerc, actTerc, 
					as.numeric(predQuant[1] - actVal), as.numeric(climQuant[1] - actVal),
					as.numeric(predQuant[2] - actVal), as.numeric(climQuant[2] - actVal),
					as.numeric(predQuant[3] - actVal), as.numeric(climQuant[3] - actVal)))
			}
		}
	}
	fwrite(sumModPerf, paste0(dataOut_location, basinName, '_summaryOfStorageModelPerformance.csv'))

#	sumModPerf = data.frame(month = NA, leadTime = NA, predTerc = NA, actTerc = NA, 
#		predMAE25 = NA, climMAE25 = NA, predMAE50 = NA, climMAE50 = NA, predMAE75 = NA, climMAE75 = NA)

	highTercHeatMap = matrix(nrow=7, ncol=12)
	lowTercHeatMap  = matrix(nrow=7, ncol=12)
	medTercHeatMap  = matrix(nrow=7, ncol=12)
	tercPredHighBarplot = matrix(nrow=3, ncol=7)
	tercPredMedBarplot  = matrix(nrow=3, ncol=7)
	tercPredLowBarplot  = matrix(nrow=3, ncol=7)
	tercActHighBarplot = matrix(nrow=3, ncol=7)
	tercActMedBarplot  = matrix(nrow=3, ncol=7)
	tercActLowBarplot  = matrix(nrow=3, ncol=7)
	absMae25HeatMap = matrix(nrow=7, ncol=12)
	absMae50HeatMap = matrix(nrow=7, ncol=12)
	absMae75HeatMap = matrix(nrow=7, ncol=12)
	relMae25HeatMap = matrix(nrow=7, ncol=12)
	relMae50HeatMap = matrix(nrow=7, ncol=12)
	relMae75HeatMap = matrix(nrow=7, ncol=12)
	
	
	for(jj in 1:7)	{
		highPredTerc  = subset(sumModPerf, predTerc == 3 & leadTime == jj)
		medPredTerc   = subset(sumModPerf, predTerc == 2 & leadTime == jj)
		lowPredTerc   = subset(sumModPerf, predTerc == 1 & leadTime == jj)
		
		highActTerc = subset(sumModPerf, actTerc == 3 & leadTime == jj)
		medActTerc  = subset(sumModPerf, actTerc == 2 & leadTime == jj)
		lowActTerc  = subset(sumModPerf, actTerc == 1 & leadTime == jj)
		
		tercPredHighBarplot[1:3,jj] = c(
			length(which(highPredTerc$actTerc == 3)) / nrow(highPredTerc), 
			length(which(highPredTerc$actTerc == 2)) / nrow(highPredTerc),
			length(which(highPredTerc$actTerc == 1)) / nrow(highPredTerc))
		
		tercPredMedBarplot[1:3,jj] = c(
			length(which(medPredTerc$actTerc == 3)) / nrow(medPredTerc), 
			length(which(medPredTerc$actTerc == 2)) / nrow(medPredTerc),
			length(which(medPredTerc$actTerc == 1)) / nrow(medPredTerc))

		tercPredLowBarplot[1:3,jj] = c(
			length(which(lowPredTerc$actTerc == 3)) / nrow(lowPredTerc), 
			length(which(lowPredTerc$actTerc == 2)) / nrow(lowPredTerc),
			length(which(lowPredTerc$actTerc == 1)) / nrow(lowPredTerc))

		tercActHighBarplot[1:3,jj] = c(
			length(which(highActTerc$predTerc == 3)) / nrow(highActTerc), 
			length(which(highActTerc$predTerc == 2)) / nrow(highActTerc),
			length(which(highActTerc$predTerc == 1)) / nrow(highActTerc))
		
		tercActMedBarplot[1:3,jj] = c(
			length(which(medActTerc$predTerc == 3)) / nrow(medActTerc), 
			length(which(medActTerc$predTerc == 2)) / nrow(medActTerc),
			length(which(medActTerc$predTerc == 1)) / nrow(medActTerc))

		tercActLowBarplot[1:3,jj] = c(
			length(which(lowActTerc$predTerc == 3)) / nrow(lowActTerc), 
			length(which(lowActTerc$predTerc == 2)) / nrow(lowActTerc),
			length(which(lowActTerc$predTerc == 1)) / nrow(lowActTerc))


		for(ii in 1:12)	{
			subModPerf = subset(sumModPerf, month == ii & leadTime == jj)
			
			absMae25HeatMap[jj,ii] = median(subModPerf$predMAE25)
			absMae50HeatMap[jj,ii] = median(subModPerf$predMAE50)
			absMae75HeatMap[jj,ii] = median(subModPerf$predMAE75)

			relMae25HeatMap[jj,ii] = median((abs(subModPerf$predMAE25) - abs(subModPerf$climMAE25)) / abs(subModPerf$climMAE25)) * 100
			relMae50HeatMap[jj,ii] = median((abs(subModPerf$predMAE50) - abs(subModPerf$climMAE50)) / abs(subModPerf$climMAE50)) * 100 
			relMae75HeatMap[jj,ii] = median((abs(subModPerf$predMAE75) - abs(subModPerf$climMAE75)) / abs(subModPerf$climMAE75)) * 100

			subHighModPerf = subset(subModPerf, actTerc == 3)
			subMedModPerf = subset(subModPerf, actTerc == 2)
			subLowModPerf = subset(subModPerf, actTerc == 1)
			highTercHeatMap[jj,ii] = length(which(subHighModPerf$predTerc == 3)) / nrow(subHighModPerf)
			medTercHeatMap[jj,ii] = length(which(subMedModPerf$predTerc == 2)) / nrow(subMedModPerf)
			lowTercHeatMap[jj,ii] = length(which(subLowModPerf$predTerc == 1)) / nrow(subLowModPerf)
		}
	}
	row.names(tercActHighBarplot) = c('high', 'med', 'low')
	row.names(tercActMedBarplot) = c('high', 'med', 'low')
	row.names(tercActLowBarplot) = c('high', 'med', 'low')
	row.names(tercPredHighBarplot) = c('high', 'med', 'low')
	row.names(tercPredMedBarplot) = c('high', 'med', 'low')
	row.names(tercPredLowBarplot) = c('high', 'med', 'low')
	
	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11) # colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = seq(-100, 100, length.out = 12)
	
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\SummaryPlot_relMAE50.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	if(any(relMae50HeatMap > 100))	{relMae50HeatMap[which(relMae50HeatMap > 100)] = 99.9}
	image(t(relMae50HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae50HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()
			
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\SummaryPlot_relMAE25.png'), width = 720, height = 720)
	if(any(relMae25HeatMap > 100))	{relMae25HeatMap[which(relMae25HeatMap > 100)] = 99.9}
	image(t(relMae25HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae25HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\SummaryPlot_relMAE75.png'), width = 720, height = 720)
	if(any(relMae75HeatMap > 100))	{relMae75HeatMap[which(relMae75HeatMap > 100)] = 99.9}
	image(t(relMae75HeatMap),# Rowv = NA, Colv = NA, scale='none',
		col = rev(rampcols), breaks = rampbreaks,
		main='Relative MAE (% change)', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(relMae75HeatMap, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileHeatMap_high.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(highTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% High Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(highTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileHeatMap_med.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(medTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% Med Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(medTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()


	rampcols = colorRampPalette(c('#B91863','white', '#196CE1'))(11)# colors are 'blue' to 'red' but in CAi scheme
	rampbreaks = c(seq(0,0.30, length.out = 6), seq(0.34,1,length.out=6))
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileHeatMap_low.png'), width = 720, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	image(t(lowTercHeatMap),
		col = rampcols, breaks = rampbreaks,
		main='% Low Tercile Correct', ylab='Months Out', xlab='Month', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	text(x=rep(1:12,each=7)/11 - 0.09, y=rep(1:7,12)/6 - 0.18, labels = paste0(round(lowTercHeatMap * 100, 0), '%'))
	axis(1, at=1:12 / 11 - 0.09,
		labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
		col.lab='#1A232F', col.axis='#666D74')
	axis(2, at = 1:7 / 6 - 0.18,
		labels =1:7,
		col.lab='#1A232F', col.axis='#666D74')
	dev.off()
	
	
		# terciles based on prediction values
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercilePredBarPlot_high.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredHighBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='High Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
	
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercilePredBarPlot_med.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredMedBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Med Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()

	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercilePredBarPlot_low.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercPredLowBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Low Tercile', ylab='% Tercile Correct by Forecast', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()


		# terciles based on actual values
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileActBarPlot_high.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActHighBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='High Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
	
	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileActBarPlot_med.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActMedBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Med Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()

	png(paste0(dataOut_location, basinName, '_projectedStorageOutputFigures\\tercileActBarPlot_low.png'), width = 900, height = 720)
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	barplot(tercActLowBarplot * 100, beside=TRUE, ylim=c(0,100),
		col = c('#B91863', '#A9BF2C', '#039CE2'), 
		main='Low Tercile', ylab='% Tercile Correct by Actual', xlab='Lead Time', xaxt='n', yaxt='n',
		col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
		family='A')
	axis(2, at=seq(0,100,20) ,
		#labels = 3*(1:7),
		col.lab='#1A232F', col.axis='#666D74')
	axis(1, at=seq(3,(7*4), 4) -1.5 ,
		labels = 1:7,
		col.lab='#1A232F', col.axis='#666D74')
	abline(h=100/3, col = '#666D74', lwd = 3, cex=1.5)
	dev.off()
}









#############################################################################################################################
## everything under this line is still in progress














reservoirOutflowCalVal_f = function(
	dataOut_location = 'save_file_location',
	historicReservoirFileLoc = 'insert_file_loc_here.com',
	historicStreamflowFileLoc = 'insert_file_loc_here.com',
	basinName = 'inster_basin_or_outlet_name',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=5000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 5,
	numRepeats = 12)
	{
	

		# create the folder for storing outputs
	if(!file.exists(paste0(dataOut_location, basinName, '_reservoirModel')))	{
			dir.create(file.path(paste0(dataOut_location, basinName, '_reservoirModel')))
	}

	############ loading libraries
	library(data.table)
	library(lubridate)
	library(CAST)
	library(caret)
	library(sf)
	library(zoo)

		# reading in data on reservoir inflows and storage
	if(dataSource == 1)	{
		inflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		stor = read.csv(paste0(historicReservoirFileLoc))
		#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

		stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
		stor$stor = as.numeric(stor$VALUE)
		inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
		inflow$inflow = inflow$VALUE
		allDat = inflow[stor[c('Date','stor')], on = 'Date']
		whichFirstNoNA = which(!is.na(allDat$stor) & !is.na(allDat$inflow))[1]
		allDat = allDat[-c(1:(whichFirstNoNA-1)), ]
		allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
			# resLoss = Outflow + AET + Withdrawals = dS - Inflow
		allDat$storIntrp = na.approx(allDat$stor)
		allDat$infIntrp = na.approx(allDat$AcFtPDay_inflow)
		allDat$resLoss = c(diff(allDat$storIntrp) - allDat$infIntrp[-nrow(allDat)], NA) * -1
		allDat$resLoss[allDat$resLoss < 0] = 0
	} else	{print("need to figure out how to handle these data")}
	

	resComplete = allDat[-nrow(allDat), c('Date')]
		# rounding to help reduce overfitting
	resComplete = cbind(resComplete, signif(allDat[-nrow(allDat), c('resLoss', 'storIntrp', 'infIntrp')], 4))
	resComplete$month = as.factor(month(resComplete$Date))
	
	#resComplete$yday = as.factor(yday(resComplete$Date))
	
	resTest = resComplete[(floor(nrow(resComplete) * .33) + 1):nrow(resComplete), ]
	resValid = resComplete[1:(floor(nrow(resComplete) * .33)), ]
	
	
	resStorModel = train(resTest[,-c(1,2)],
						resTest$resLoss,
						metric=modelMetric,
						method=modelMethod,
						#tuneLength=?,
						imporance=TRUE,
						ntree=nTreeModel,
						trControl = trainControl(
							method = 'repeatedcv',
							number = numFolds,
							repeats = numRepeats))
							#method='cv'))
		#					index = indices$index))
		#					p=.75))

			#saving the model
	saveRDS(resStorModel, paste0(dataOut_location, basinName, '_reservoirModel\\dailyProj.RData'))

	predValid = predict(resStorModel, newdata = resValid)
	predCalib = predict(resStorModel, newdata = resTest)
			
		# save ML model performance
	png(paste0(dataOut_location, basinName, '_reservoirModel\\', basinName, '_MLcalib_dailyProj.png'), width = 720, height = 720)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	plot(resTest$resLoss, predCalib, 
		main='', xlab = 'Act. Reservoir Loss (AcFt / Month)', ylab = 'Est. Reservoir Loss (AcFt / Month)',
		pch=2, lwd=1.15,  col=alpha('#0098B2', 0.5), cex=1.75)
		
	abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
	abline(lm(predCalib ~ resTest$resLoss),
		lwd=2, col=alpha('#F06000', 1), cex=2)
	text(x=min(resTest$resLoss, na.rm=TRUE), 
			y=max(predCalib, na.rm=TRUE) * 0.9,
			paste0('r  = ', round(summary(lm(predCalib ~ resTest$resLoss))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
	text(x=min(resTest$resLoss, na.rm=TRUE), 
			y=max(predCalib, na.rm=TRUE) * 0.92,
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
	text(x=max(resTest$resLoss, na.rm=TRUE) ,
			y=min(predCalib, na.rm=TRUE) * .9,
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
	text(x=min(resTest$resLoss, na.rm=TRUE) ,
			y=max(predCalib, na.rm=TRUE) * .82,
			paste0('Calib. Pred. Daily'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.5*1.9)
	dev.off()


		# save validation of ML model
	png(paste0(dataOut_location, basinName, '_reservoirModel\\', basinName, '_MLvalid_dailyProj.png'), width = 720, height = 720)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	plot(resValid$resLoss, predValid, 
		main='', xlab = 'Act. Reservoir Loss (AcFt / Month)', ylab = 'Est. Reservoir Loss (AcFt / Month)',
		pch=2, lwd=1.15,  col=alpha('#0098B2', 0.5), cex=1.75)
	
	abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
	abline(lm(predValid ~ resValid$resLoss),
		lwd=2, col=alpha('#F06000', 1), cex=2)
	text(x=min(resValid$resLoss, na.rm=TRUE), 
			y=max(predValid, na.rm=TRUE) * 0.9,
			paste0('r  = ', round(summary(lm(predValid ~ resValid$resLoss))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
	text(x=min(resValid$resLoss, na.rm=TRUE), 
			y=max(predValid, na.rm=TRUE) * 0.92,
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
	text(x=max(resValid$resLoss, na.rm=TRUE) ,
			y=min(predValid, na.rm=TRUE) * .9,
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
	text(x=min(resValid$resLoss, na.rm=TRUE) ,
			y=max(predValid, na.rm=TRUE) * .82,
			paste0('Valid. Pred. Daily'),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.5*1.9)
	dev.off()
}
























######################
# under construction

































###################
## !!! todo
storageProjectionValidationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	climateInputsFileLoc = 'file_location_and_name',	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	infOrFNF = 'inflow',
	dataSource = 1)							# 1 for FNF from cal.gov,
	{
	
	if(!file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		print('gotta calibrate the model first you big ole dummy!!')
	}	else	{
			# creating directory for holding figures
		if(!file.exists(paste0(dataOut_location, basinName, '_projectedOutputFigures')))	{
			dir.create(file.path(paste0(dataOut_location, basinName, '_projectedOutputFigures')))
		}


		# read in the necessary libraries
		library(data.table)
		library(lubridate)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in climate data
		seas5ClimateInputFileNames = list.files(paste0(dataOut_location, basinName, '_seas5MultiOutput'))
		era5ClimateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Historical.RData')))

			# reading in hydro model calibration values
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)


		###############
		## reading in data and models for reservoirs
		# reading in historic reservoir data
				# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			basinAreaInKm = sum(st_read(paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"))$SUB_AREA)
			acrePerSqKm = 247.105
			mmPerFoot = 304.8
			mmToAcrFt = basinAreaInKm * acrePerSqKm / mmPerFoot

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
			
	
				# building a record of historical sorage for comparison
			allDat$year = year(allDat$Date)
			allDat$month = month(allDat$Date)
			climatologyDF = data.frame(month = NA, year = NA, stor = NA)
			climIter = 0
			for(thisYear in unique(allDat$year))	{
				for(thisMonth in 1:12)	{
					climIter = climIter + 1
					histSubset = subset(allDat, year == thisYear & month == thisMonth)
					climatologyDF[climIter, ] = c(thisMonth, thisYear, mean(histSubset$stor, na.rm=TRUE))
				}
			}
			
			allDat$ymonth = paste0(year(allDat$Date), month(allDat$Date))
			monthDat = data.frame(ymonth = NA, month = NA, year = NA, stor = NA, AcFtPDay_inflow = NA, resLoss = NA)
			k = 0
		for(i in unique(allDat$ymonth))	{
				k = k + 1
				subsetDat = subset(allDat, ymonth == i)
				monthDat = rbind(monthDat, c(i, month(subsetDat$Date[1]), year(subsetDat$Date[1]), 
					last(subsetDat$stor[!is.na(subsetDat$stor)]),
					mean(subsetDat$AcFtPDay_inflow * nrow(subsetDat), na.rm=TRUE),	# getting a sum even for months with missing data
					mean(subsetDat$resLoss * nrow(subsetDat), na.rm=TRUE)))  		# getting a sum even for months with missing data
			}
			monthDat$stor = as.numeric(monthDat$stor)
			monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
			monthDat$resLoss = as.numeric(monthDat$resLoss)
			
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}
			#backfilling missing reservoir data, as ML can't handle NAs
		monthDat$stor[is.na(monthDat$stor)] = mean(monthDat$stor, na.rm=TRUE)
		monthDat$AcFtPDay_inflow[is.na(monthDat$AcFtPDay_inflow)] = mean(monthDat$AcFtPDay_inflow, na.rm=TRUE)
		monthDat$resLoss[is.na(monthDat$resLoss)] = mean(monthDat$resLoss, na.rm=TRUE)
		monthAvgResLoss = NULL
		for(kk in 1:12)	{
			monthAvgResLoss = c(monthAvgResLoss, mean(subset(monthDat, month==kk)$resLoss))
		}
		monthAvgResLoss = c(monthAvgResLoss, monthAvgResLoss)
		
###############!!!!!!!!!!


		reservoirVals = list()
		reservoirR2s = NULL
			# reading in and rebuilding the ML model for estimating reservoir outflow
		for(thisMLmod in 1:6)	{
			model_LLO = readRDS(paste0(dataOut_location, basinName, '_reservoirModel\\', thisMLmod, '_monthsOut.RData'))
			bestLags = model_LLO[[length(model_LLO)]]
			
			reservoirVals[[thisMLmod]] = monthDat[, c('year', 'month', 'ymonth','stor','resLoss')]
	
			for(j in 1:nrow(bestLags))	{
				storLag = bestLags$bestStorLags[j]
				storRolMn = rollmean(monthDat$stor, storLag)
				storRolMnLn = length(storRolMn)
				
				infLag = bestLags$bestInfLags[j]
				infRolMn = rollmean(monthDat$AcFtPDay_inflow, infLag)
				infRolMnLn = length(infRolMn)
			
				reservoirVals[[thisMLmod]][ , paste0('stor_mn_', storLag)] = 
					c(rep(NA, thisMLmod + (storLag - 1)), 
						storRolMn[-c((storRolMnLn - (thisMLmod-1)):storRolMnLn)])
				reservoirVals[[thisMLmod]][ , paste0('inf_mn_', infLag)] = 
					c(rep(NA, thisMLmod + (infLag - 1)), 
						infRolMn[-c((infRolMnLn - (thisMLmod-1)):infRolMnLn)])
			}


			# predicting reservoir outflows
			reservoirVals[[thisMLmod]] = reservoirVals[[thisMLmod]][complete.cases(reservoirVals[[thisMLmod]]), ]	# temporary method for backfilling NAs when data are missing
			reservoirVals[[thisMLmod]]$Predictions = predict(model_LLO, newdata = reservoirVals[[thisMLmod]])
			reservoirR2s = c(reservoirR2s, max(model_LLO$results$Rsquared))
		}



		save_iter = 0
		validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predStor = NA, actStor = NA)

		for(thisSeas5 in seas5ClimateInputFileNames)	{
			seas5ClimateInput = readRDS(paste0(dataOut_location, 'seas5MultiOutput\\', thisSeas5))

			# removing monthly forecasts that are not a full month
			partialMonths = FALSE
			if(last(mday(seas5ClimateInput[[1]]$Date + 1)) != 1) 	{
				partialMonths = TRUE
				rmSEAS5DateRm = length(last(which(mday(seas5ClimateInput[[1]]$Date) == 1)):nrow(seas5ClimateInput[[1]]$Date)) - 1
			}
			
			lastHistData = which(era5ClimateInput$Date == seas5ClimateInput[[1]]$Date[1])
			allForecastsOutput = list()
			iter = 0
			for(numModels in 1:length(seas5ClimateInput))	{
				for(numCalibs in 1:nrow(calibratedVars))	{	
					iter = iter + 1
					
					climateInput = rbind(era5ClimateInput[1:lastHistData,], seas5ClimateInput[[numModels]])
					if(partialMonths) 	{climateInput = climateInput[-c((nrow(climateInput)-rmSEAS5DateRm):nrow(climateInput)), ]}
					
					projDates = (1 + lastHistData):nrow(climateInput)



					allHBVoutput = runHBV_f(
						climateInput = climateInput,
						sfcf = calibratedVars$sfcf[numCalibs],	#snowfall correction factor [-]
						tr   = calibratedVars$tr[numCalibs],	#solid and liquid precipitation threshold temperature [C]
						tt   = calibratedVars$tt[numCalibs],	#melt temperature [C]
						fm   = calibratedVars$fm[numCalibs],	#snowmelt factor [mm/C]
						fi   = calibratedVars$fi[numCalibs],	#icemelt factor [mm/C]
						fic  = calibratedVars$fic[numCalibs],	#debris-covered icemelt factor [mm/C]
						fc   = calibratedVars$fc[numCalibs], # field capacity 
						lp   = calibratedVars$lp[numCalibs], # parameter for PET --> AET
						beta_soils = calibratedVars$beta_soils[numCalibs], #soil moisture drainage exponent
						k0   = calibratedVars$k0[numCalibs],	# storage constant of top bucket
						k1   = calibratedVars$k1[numCalibs],	# storage constant of middle bucket
						k2   = calibratedVars$k2[numCalibs],	# storage constant of bottom bucket	
						uz1  = calibratedVars$uz1[numCalibs], #max flux rate from STZ to SUZ in mm/d
						perc = calibratedVars$perc[numCalibs])  # max flux rate from SUZ to SLZ in mm/d

					# merging historic streamflow record onto climate inputs data
#					climateAndStreamflowOutput = historicStreamflow[allHBVoutput, on='Date'][projDates, ]  
#					plotRows = (lastHistData-365):nrow(climateAndStreamflowOutput)
#					plot(climateAndStreamflowOutput$Date[plotRows], climateAndStreamflowOutput$historicQinmm[plotRows], type='l')
#					lines(climateAndStreamflowOutput$Date[plotRows], climateAndStreamflowOutput$Qg[plotRows], col='red2')
#					abline(v=climateAndStreamflowOutput$Date[lastHistData])
					projectedOutput = allHBVoutput[projDates,]
					projectedOutput$QinAcFt = projectedOutput$Qg * mmToAcrFt

					projectedOutput$month = month(projectedOutput$Date)
					projectedOutput$year = year(projectedOutput$Date)
					monthsOut = 0
					predCumsum = 0
					outflowCumsum = 0
#					medResLoss = median(reservoirVals[[1]]$resLoss)
					for(thisMonth in unique(projectedOutput$month))	{
						monthsOut = monthsOut + 1
						if(monthsOut <=6)	{
							outputSubset = subset(projectedOutput, month == thisMonth)
							initialStorRow = which(reservoirVals[[monthsOut]]$year == outputSubset$year[1] & 
								reservoirVals[[monthsOut]]$month == outputSubset$month[1]) - 1
							initialStor = reservoirVals[[monthsOut]]$stor[initialStorRow]
							futureStor = reservoirVals[[monthsOut]]$stor[initialStorRow + monthsOut]
#							totalOutflow = sum(reservoirVals[[monthsOut]]$Predictions[initialStorRow:(initialStorRow + monthsOut - 1)])
#							predRatio = 1 * (1 / monthsOut)
							predRatio = reservoirR2s[monthsOut]
							numDays = nrow(outputSubset)
							predCumsum = predCumsum + mean(outputSubset$QinAcFt) * numDays

							outflowCumsum = min(
								outflowCumsum + monthAvgResLoss[thisMonth + monthsOut - 1] * (1-predRatio) + 
									predRatio * reservoirVals[[monthsOut]]$Predictions[initialStorRow + monthsOut - 1],
								initialStor + predCumsum)

							predStor = initialStor + predCumsum - outflowCumsum
							actStor = reservoirVals[[monthsOut]]$stor[initialStorRow + 1]

							validDF = rbind(validDF,
								c(numModels, numCalibs, outputSubset$Date[1] + 14, monthsOut,
									predStor, actStor))
						}
					}
					if(nrow(validDF) > 20000)	{
						save_iter = save_iter + 1
						fwrite(validDF, paste0(dataOut_location, "temp_out_", save_iter, ".csv"))
						validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predStor = NA, actStor = NA)
					}
				}
			}
		}

			# combining saved iterations and returning the dataframe with formatted climatology
		save_iter = save_iter + 1
		fwrite(validDF, paste0(dataOut_location, "temp_out_", save_iter, ".csv"))
			
		validDF = fread(paste0(dataOut_location, "temp_out_", 1, ".csv"))
		for(files in 2:save_iter)	{
			validDF = rbind(validDF, fread(paste0(dataOut_location, "temp_out_", files, ".csv")))
		}
	
		validDF = validDF[complete.cases(validDF), ]
		validDF$Date = validDF$days + as.Date("1970-01-01")
		validDF$startDate = validDF$Date %m-% (months(validDF$monthsOut - 1)) 
		
		startDates = unique(validDF$startDate) 
			
		
		for(thisStartDate in startDates)	{
			validSubset = subset(validDF, startDate == thisStartDate)
			
				# parsing the climatology df for plotting
			theseMonths = unique(month(validSubset$Date))
			numMonths = length(theseMonths)
			climDF = data.frame(month = NA, stor = NA, plotMonth = NA)
			for(startYear in climatologyDF$year[1]:(max(climatologyDF$year)-1))	{
				climSubset = subset(climatologyDF, year >= startYear)
				climStart = which(climSubset$month == month(validSubset$startDate)[1])[1]
				
				climDF = rbind(climDF,
					data.frame(month = theseMonths, stor = climSubset$stor[climStart:(climStart+numMonths-1)], plotMonth = c(1:6)))
			}	
		
			
			
			png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\storage_', validSubset$Date[1] - 14, '.png'), width = 1200, height = 720)
				windowsFonts(A = windowsFont("Roboto"))
			par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
			#validSubset$plotDate = rep(1:7, nrow(validSubset)/7) - .2
			boxplot(validSubset$predStor ~ validSubset$Date, pch='.', lwd=1, 
				#col=colorRampPalette(c('#F06000', '#B91863'))(length(unique(validDF$monthsOut))),
				ylim = c(0, max(validDF$predStor)),
				at = c(1:6) - 0.2,
				boxwex = 0.4,
				col = '#A9BF2C',
				border='#666D74', cex=1.5, 
				main='', ylab='Total Storage (AcFt)', xlab='Month', xaxt='n',
				col.lab='#1A232F', col.axis='#666D74', col.main='#1A232F',
				family='A')
			#climDF$plotDate =climDF$plotMonth + 0.2
			boxplot(stor ~ plotMonth, data = climDF, pch='.',
				col = '#FDB600', 
				at = c(1:6) + 0.2,
				boxwex = 0.4, 
				border='#666D74', cex=1.5, 
				main='', ylab='Total Storage (AcFt)', xlab='Month', xaxt='n',
				add=TRUE)
			lines(1:length(unique(validSubset$Date)), unique(validSubset$actStor), type='l', col='#196CE1', lwd=1)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actStor), pch = '---', col='#F2F3F3', lwd=2, cex=9)
			points(1:length(unique(validSubset$Date)), unique(validSubset$actStor), pch = '---', col='#196CE1', lwd=2, cex=7)
			axis(1, at=1:length(unique(validSubset$Date)),
				labels =c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')[month(unique(validSubset$Date))],
				col.lab='#1A232F', col.axis='#666D74')
			text(x=0.5, y=max(validDF$predStor) * 0.95,
				paste0(basinName),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predStor) * 0.87,
				paste0('Forecast Beginning on ', validSubset$Date[1] - 14),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predStor) * 0.75,
				paste0('CAi Forecast'),
				adj = c(0,0), font=2, col='#A9BF2C', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predStor) * 0.7,
				paste0('Climatology'),
				adj = c(0,0), font=2, col='#FDB600', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predStor) * 0.65,
				paste0('Actual'),
				adj = c(0,0), font=2, col='#196CE1', family='A', cex=1.6*1.3)
			dev.off()
		}

	}
}

































































####################################
## old functions no longer in use


#original reservoir model outflow model
OLD_reservoirOutflowCalVal_f = function(
	dataOut_location = 'save_file_location',
	historicReservoirFileLoc = 'insert_file_loc_here.com',
	historicStreamflowFileLoc = 'insert_file_loc_here.com',
	basinName = 'inster_basin_or_outlet_name',
	basinSymbol = 'EXC',
	infOrFNF = 'inflow',
	dataSource = 1, 			# 1 for FNF or inflow from cal.gov,
	nTreeModel=5000,
	modelMetric = 'Rsquared',
	modelMethod = 'rf',
	metric="Rsquared",
	numFolds = 5,
	numRepeats = 12)
	{
	
	############ loading libraries
	library(data.table)
	library(lubridate)
	library(CAST)
	library(caret)
	library(sf)
	library(zoo)

	if(dataSource == 1)	{
		inflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
		stor = read.csv(paste0(historicStreamflowFileLoc))
		#outflow = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=23&dur_code=D&Start=1900-01-01&End=2022-08-05")

		stor$Date = ymd(unlist(strsplit(stor$DATE.TIME, " "))[seq(1,nrow(stor)*2,2)])
		stor$stor = as.numeric(stor$VALUE)
		inflow$Date = ymd(unlist(strsplit(inflow$DATE.TIME, " "))[seq(1,nrow(inflow)*2,2)])
		inflow$inflow = inflow$VALUE
		allDat = inflow[stor[c('Date','stor')], on = 'Date']
		allDat$AcFtPDay_inflow =  as.numeric(allDat$inflow)  * (60*60*24) / 43559.9
			# resLoss = Outflow + AET + Withdrawals = dS - Inflow
		allDat$resLoss = c(diff(allDat$stor) - allDat$AcFtPDay_inflow[-1], NA) * -1
	} else	{print("need to figure out how to handle these data")}
	
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
			mean(subsetDat$AcFtPDay_inflow * 30.4375, na.rm=TRUE),
			mean(subsetDat$resLoss * 30.4375, na.rm=TRUE)))
	}
	monthDat$stor = as.numeric(monthDat$stor)
	monthDat$AcFtPDay_inflow = as.numeric(monthDat$AcFtPDay_inflow)
	monthDat$resLoss = as.numeric(monthDat$resLoss)

	if(!file.exists(paste0(dataOut_location, basinName, '_reservoirModel')))	{
		dir.create(file.path(paste0(dataOut_location, basinName, '_reservoirModel')))
	}

		# saving uncertainty in I + O = dS
	png(paste0(dataOut_location, basinName, '_reservoirModel\\', basinName, '_monthlyReservoirDstorVsLoss_historical.png'), width = 720, height = 720)
	windowsFonts(A = windowsFont("Roboto"))
	par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
	plot(diff(monthDat$stor), (monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1]), 
		main='', xlab = 'Delta Reservoir Storage (AcFt / Month)', ylab = 'Est. Inflow - Outflow (AcFt / Month)',
		pch=2, lwd=1.15,  col=alpha('#0098B2', 0.5), cex=1.75)
	
	abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
	abline(lm((monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1]) ~ diff(monthDat$stor)),
		lwd=2, col=alpha('#F06000', 1), cex=2)
	text(x=min(diff(monthDat$stor), na.rm=TRUE), 
			y=max(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], na.rm=TRUE) * 0.9,
			paste0('r  = ', round(summary(lm((30 * (monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1])) ~ diff(monthDat$stor)))$r.squared, 2)),
			adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
	text(x=min(diff(monthDat$stor), na.rm=TRUE), 
			y=max(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], na.rm=TRUE) * 0.96,
			'  2',
			adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
	text(x=max(diff(monthDat$stor), na.rm=TRUE) ,
			y=min(monthDat$AcFtPDay_inflow[-1] - monthDat$resLoss[-1], na.rm=TRUE) * .9,
			paste0(basinName),
			adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
	dev.off()

	
		# initializing model predictions
	for(h in 1:6)	{
		recLength = nrow(monthDat)
		lagLength = seq(1, (12*2 - h), 1)

		corResults = matrix(data = NA, nrow = length(lagLength), ncol = 2)
		k=0
		for(i in lagLength)	{
			k=k+1
			storRolMn = rollmean(monthDat$stor, i)
			storRolMnLn = length(storRolMn)
			corResults[k, 1] = cor.test(
				c(rep(NA, h + (i - 1)), storRolMn[-c((storRolMnLn - (h-1)):storRolMnLn)]),
				monthDat$resLoss,
				method='kendall')$p.value
			infRolMn = rollmean(monthDat$AcFtPDay_inflow, i)
			infRolMnLn = length(infRolMn)
			corResults[k, 2] = cor.test(
				c(rep(NA, h + (i - 1)), infRolMn[-c((infRolMnLn - (h-1)):infRolMnLn)]),
				monthDat$resLoss,
				method='kendall')$p.value
		}
	

		bestStorLags = c(order(corResults[,1])[c(1:5)], max(lagLength))
		bestInfLags = c(order(corResults[,2])[c(1:5)], max(lagLength))

		resSubset = monthDat[, c('ymonth','resLoss','month')]

		for(j in 1:length(bestInfLags))	{
			storLag = bestStorLags[j]
			storRolMn = rollmean(monthDat$stor, storLag)
			storRolMnLn = length(storRolMn)
			
			infLag = bestInfLags[j]
			infRolMn = rollmean(monthDat$AcFtPDay_inflow, infLag)
			infRolMnLn = length(infRolMn)
		
			resSubset[ , paste0('stor_mn_', storLag)] = 
				c(rep(NA, h + (storLag - 1)), 
					storRolMn[-c((storRolMnLn - (h-1)):storRolMnLn)])
			resSubset[ , paste0('inf_mn_', infLag)] = 
				c(rep(NA, h + (infLag - 1)), 
					infRolMn[-c((infRolMnLn - (h-1)):infRolMnLn)])
		}

		resComplete = resSubset[complete.cases(resSubset), ]
		resTest = resComplete[(floor(nrow(resComplete) * .33) + 1):nrow(resComplete), ]
		resValid = resComplete[1:(floor(nrow(resComplete) * .33)), ]

		model_LLO = train(resTest[,-c(1,2)],
							resTest$resLoss,
							metric=modelMetric,
							method=modelMethod,
							#tuneLength=?,
							imporance=TRUE,
							ntree=nTreeModel,
							trControl = trainControl(
								method='repeatedcv',
								number=numFolds,
								repeats=numRepeats))
								#method='cv'))
			#					index = indices$index))
			#					p=.75))

		#model_LLO
		#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
		#varImp(model_LLO)$importance	
		#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
		#	sortedVars

			# appending the best time lags to the model for rerunning
		model_LLO[[length(model_LLO) + 1]] = data.frame(bestStorLags = bestStorLags, bestInfLags = bestInfLags)


			# save ML model
		model_LLO[[length(model_LLO)+1]] = data.frame(bestStorLags = bestStorLags, bestInfLags = bestInfLags)
		saveRDS(model_LLO, paste0(dataOut_location, basinName, '_reservoirModel\\', h, '_monthsOut.RData'))

		predValid = predict(model_LLO, newdata = resValid)
		predCalib = predict(model_LLO, newdata = resTest)
			
			# save ML model performance
		png(paste0(dataOut_location, basinName, '_reservoirModel\\', basinName, '_MLcalib_', h, '_monthsOut.png'), width = 720, height = 720)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
		plot(resTest$resLoss, predCalib, 
			main='', xlab = 'Act. Reservoir Loss (AcFt / Month)', ylab = 'Est. Reservoir Loss (AcFt / Month)',
			pch=2, lwd=1.15,  col=alpha('#0098B2', 0.5), cex=1.75)
		
		abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
		abline(lm(predCalib ~ resTest$resLoss),
			lwd=2, col=alpha('#F06000', 1), cex=2)
		text(x=min(resTest$resLoss, na.rm=TRUE), 
				y=max(predCalib, na.rm=TRUE) * 0.9,
				paste0('r  = ', round(summary(lm(predCalib ~ resTest$resLoss))$r.squared, 2)),
				adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=min(resTest$resLoss, na.rm=TRUE), 
				y=max(predCalib, na.rm=TRUE) * 0.92,
				'  2',
				adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
		text(x=max(resTest$resLoss, na.rm=TRUE) ,
				y=min(predCalib, na.rm=TRUE) * .9,
				paste0(basinName),
				adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
		text(x=min(resTest$resLoss, na.rm=TRUE) ,
				y=max(predCalib, na.rm=TRUE) * .82,
				paste0('Calib. Pred. ', h, ' Months Out'),
				adj = c(0,0), font=2, col='#F06000', family='A', cex=1.5*1.9)
		dev.off()


			# save validation of ML model
		png(paste0(dataOut_location, basinName, '_reservoirModel\\', basinName, '_MLpreds_', h, '_monthsOut.png'), width = 720, height = 720)
		windowsFonts(A = windowsFont("Roboto"))
		par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
		plot(resValid$resLoss, predValid, 
			main='', xlab = 'Act. Reservoir Loss (AcFt / Month)', ylab = 'Est. Reservoir Loss (AcFt / Month)',
			pch=2, lwd=1.15,  col=alpha('#0098B2', 0.5), cex=1.75)
		
		abline(a=0, b=1, lwd=2, col='#1A232F', cex=2)
		abline(lm(predValid ~ resValid$resLoss),
			lwd=2, col=alpha('#F06000', 1), cex=2)
		text(x=min(resValid$resLoss, na.rm=TRUE), 
				y=max(predValid, na.rm=TRUE) * 0.9,
				paste0('r  = ', round(summary(lm(predValid ~ resValid$resLoss))$r.squared, 2)),
				adj = c(0,0), font=2, col='#F06000', family='A', cex=2*1.6)
		text(x=min(resValid$resLoss, na.rm=TRUE), 
				y=max(predValid, na.rm=TRUE) * 0.92,
				'  2',
				adj = c(0,0), font=2, col='#F06000', family='A', cex=1.2*1.6)
		text(x=max(resValid$resLoss, na.rm=TRUE) ,
				y=min(predValid, na.rm=TRUE) * .9,
				paste0(basinName),
				adj = c(1,0), font=2, col='#1A232F', family='A', cex=1.5*1.9)
		text(x=min(resValid$resLoss, na.rm=TRUE) ,
				y=max(predValid, na.rm=TRUE) * .82,
				paste0('Valid. Pred. ', h, ' Months Out'),
				adj = c(0,0), font=2, col='#F06000', family='A', cex=1.5*1.9)
		dev.off()


	}
}

