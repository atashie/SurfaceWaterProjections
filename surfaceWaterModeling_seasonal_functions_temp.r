
##########################################################################################################
	## basin delineation
	
basinDelineation_f = function(
	gageLonLat = c(-79,36),
	basinATLAS_locationAndFile = 'basinATLAS_location.gdb',
	dataOut_location = 'save_file_location',
	basinName = 'inster_basin_or_outlet_name')
	{
	# check to make sure the basin has not already been delineated
	if(file.exists(paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")))	{
		print("This file already exists you big ole dummy!")
	}	else	{

		# load relevant libraries
	library(sf)				# for geospatial data
	sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat
	library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
	
	
	##################################################
	#### reading basinATLAS data and intersecting with gage location
	basinAt12 = st_read(dsn=basinATLAS_locationAndFile, layer="BasinATLAS_v10_lev12")
	gageLocation_asSF = st_sfc(st_point(gageLonLat))
	st_crs(gageLocation_asSF) = 4326
	gageIntersection = as.numeric(st_intersects(gageLocation_asSF, st_buffer(basinAt12,0))) # keeps points
	
	basinIDs = basinAt12$HYBAS_ID[gageIntersection]
	basinAt12Remaining = basinAt12[-gageIntersection,]
	
		# identify all upstream watersheds
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
	
	st_write(delineatedBasin, paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg"), append=FALSE)
	st_write(basinWatersheds, paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg"), append=FALSE)
	}
}	






#################################################################################################################
	# function for selecting climate inputs and converting to appropriate units

climateInputConversion_f = function(
	pathToBasinBoundaryGPKG = 'file_location_and_name.gpkg',
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	basinName = 'inster_basin_or_outlet_name',
	climateDataNCDF = 'file_location_and_name.nc',
	tempConversionFactor = NA,
	pptConversionFactor = NA,
	avgTempGiven = FALSE, 
	multipleModels = FALSE,	# are there multiple models that need to be stored?
	startDate = as.Date("1990-01-01"), 	# when does the clock of the netcdf start?
	timeToDaysConversion = 1,	# convert time increments to days if necessary
	dataOut_location = 'save_file_location',
	optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
	variableOrderOption = NA, # 'seas5' or 'era5' or 'seas5Multi'
	precipName = 'tp')	# other options include: tp, tp_sum	
	{
	# load relevant libraries
	library(sf)				# for geospatial data
	sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat
	library(ncdf4)			# for loading netcdf that contains lat-lons of climate data
	library(data.table)	# for data.table
	if(optionForPET == 1) {library(EcoHydRology)}	# for PET_fromTemp

	
	delineatedBasin = st_read(pathToBasinBoundaryGPKG)
	basinWatersheds = st_read(pathToWatershedsGPKG)

	basinArea = sum(basinWatersheds$SUB_AREA)

	if(variableOrderOption == 'seas5')	{
		# reading in the ncdfs
		ncin = nc_open(climateDataNCDF)

		ncTmin = ncvar_get(ncin, 't2m_min'); if(any(is.na(ncTmin)))	{ncTmin[is.na(ncTmin)] = median(ncTmin, na.rm=TRUE)} 
		ncTmax = ncvar_get(ncin, 't2m_max'); if(any(is.na(ncTmax)))	{ncTmax[is.na(ncTmax)] = median(ncTmax, na.rm=TRUE)} 
		ncPPT = ncvar_get(ncin, precipName); if(any(is.na(ncPPT)))	{ncPPT[is.na(ncPPT)] = median(ncPPT, na.rm=TRUE)} 
		if(avgTempGiven)	{ncTavg = ncvar_get(ncin, 't2m_avg')
		}	else {ncTavg = (ncTmin + ncTmax) / 2}	

			# assuming all lat / lon structures are the same
		nc_lat = ncvar_get(ncin, 'latitude')
		nc_lon = ncvar_get(ncin, 'longitude')
		nc_date = startDate + timeToDaysConversion * ncvar_get(ncin, 'lead_time') # time is days after jan 1 1990

		allClimateData = list()
		for(numModels in 1:length(ncvar_get(ncin, 'member')))	{
			allPPT = 0	; allTmin = 0	; allTmax = 0	; allPET = 0
			for(numberOfWatersheds in 1:nrow(basinWatersheds))	{
				thisLonLat = st_coordinates(st_centroid(basinWatersheds[numberOfWatersheds, ]))
				nearestLon = which.min(abs(thisLonLat[1] - nc_lon))
				nearestLat = which.min(abs(thisLonLat[2] - nc_lat))
			
					# Tmax is sometimes < Tmin, so need to adjust before running for PET
				theseTmin = ncTmin[nearestLon, nearestLat, numModels, ]
				theseTmax = ncTmax[nearestLon, nearestLat, numModels, ]
				if(any(theseTmin >= theseTmax))	{
					theseTmin[theseTmin >= theseTmax] = theseTmax[theseTmin >= theseTmax] - 0.1
				}
				
				allTmin = allTmin +  theseTmin *  basinWatersheds$SUB_AREA[numberOfWatersheds]
				allTmax = allTmax + theseTmax *  basinWatersheds$SUB_AREA[numberOfWatersheds]
				allPPT = allPPT + ncPPT[ , nearestLon, nearestLat, numModels] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
				if(optionForPET == 1)	{
					allPET = allPET + PET_fromTemp(yday(nc_date), theseTmax , theseTmin,
						lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
				}	else	{
					print('need to figure this one out later')
				}
			}
			avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
			avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
			avgTavg = (avgTmin + avgTmax) / 2
			avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
			avgPET = allPET / basinArea

			allClimateData[[numModels]] = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
		}
	}

	if(variableOrderOption == 'era5')	{
			# reading in the ncdfs
		ncin = nc_open(climateDataNCDF)

		ncTmin = ncvar_get(ncin, 't2m_min'); if(any(is.na(ncTmin)))	{ncTmin[is.na(ncTmin)] = median(ncTmin, na.rm=TRUE)} 
		ncTmax = ncvar_get(ncin, 't2m_max'); if(any(is.na(ncTmax)))	{ncTmax[is.na(ncTmax)] = median(ncTmax, na.rm=TRUE)} 
		ncPPT = ncvar_get(ncin, precipName); if(any(is.na(ncPPT)))	{ncPPT[is.na(ncPPT)] = median(ncPPT, na.rm=TRUE)} 
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
					
			allTmin = allTmin + ncTmin[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
			allTmax = allTmax + ncTmax[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
			allPPT = allPPT + ncPPT[nearestLon, nearestLat, ] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
			if(optionForPET == 1)	{
				allPET = allPET + PET_fromTemp(yday(nc_date), ncTmax[nearestLon, nearestLat, ], ncTmin[nearestLon, nearestLat, ],
					lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
			}	else	{
				print('need to figure this one out later')
			}
		}
		
		avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
		avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
		avgTavg = (avgTmin + avgTmax) / 2
		avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
		avgPET = allPET / basinArea
		
		# save data output
		allClimateData = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
	}
	
	if(variableOrderOption == 'seas5Multi')	{
		# creating the folder to hold the output
			if(!file.exists(paste0(dataOut_location, basinName, '_seas5MultiOutput')))	{
			dir.create(file.path(paste0(dataOut_location, basinName, '_seas5MultiOutput')))
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
			for(numModels in 1:length(ncvar_get(ncin, 'member')))	{
				allPPT = 0	; allTmin = 0	; allTmax = 0	; allPET = 0
				for(numberOfWatersheds in 1:nrow(basinWatersheds))	{
					thisLonLat = st_coordinates(st_centroid(basinWatersheds[numberOfWatersheds, ]))
					nearestLon = which.min(abs(thisLonLat[1] - nc_lon))
					nearestLat = which.min(abs(thisLonLat[2] - nc_lat))
				
						# Tmax is sometimes < Tmin, so need to adjust before running for PET
					theseTmin = ncTmin[nearestLon, nearestLat, numModels, , thisInitTime]
					theseTmax = ncTmax[nearestLon, nearestLat, numModels, , thisInitTime]
					if(any(theseTmin >= theseTmax))	{
						theseTmin[theseTmin >= theseTmax] = theseTmax[theseTmin >= theseTmax] - 0.1
					}
					
					allTmin = allTmin +  theseTmin *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					allTmax = allTmax + theseTmax *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					allPPT = allPPT + ncPPT[ , nearestLon, nearestLat, numModels, thisInitTime] *  basinWatersheds$SUB_AREA[numberOfWatersheds]
					if(optionForPET == 1)	{
						allPET = allPET + PET_fromTemp(yday(nc_date), theseTmax , theseTmin,
							lat_radians =  min((thisLonLat[1,2]*pi/180), 1.1)) * 1000  *  basinWatersheds$SUB_AREA[numberOfWatersheds]	# output in m, convert to mm
					}	else	{
						print('need to figure this one out later')
					}
				}
				avgTmin	= if(is.na(tempConversionFactor))	{allTmin / basinArea}	else	{tempConversionFactor + allTmin / basinArea}
				avgTmax = if(is.na(tempConversionFactor)) 	{allTmax / basinArea}	else	{tempConversionFactor + allTmax / basinArea}
				avgTavg = (avgTmin + avgTmax) / 2
				avgPPT = if(is.na(pptConversionFactor)) 	{allPPT / basinArea}	else	{pptConversionFactor * allPPT / basinArea}
				avgPET = allPET / basinArea

				allClimateData[[numModels]] = data.frame(Date = nc_date, Tmin = avgTmin, Tmax = avgTmax, Tavg = avgTavg, PPT = avgPPT, PET = avgPET)
			}
			saveRDS(allClimateData, paste0(dataOut_location, basinName, '_seas5MultiOutput\\SEAS5_startDate_', nc_date[1], '.RData'))
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
	climateInputsFileLoc = 'file_location_and_name.RData',
	historicStreamflowFileLoc = 'https://someplace.gov',
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	dataOut_location = 'save_file_location',
	dataSource = 1,											# 1 for FNF from cal.gov, 
	numberOfRuns = 100000,									# maximum number of runs
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
		library(hydroGOF)		# for nse calculations
		
		
		climateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Historical.RData')))
		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (as.numeric(historicStreamflow$VALUE) / flowUnitConversion) / (areaUnitConversion * basinArea)
			
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}
		
		
			# merging historic streamflow record onto climate inputs data
		climateAndStreamflowInput = historicStreamflow[climateInput, on='Date']
		lengthOfCalibrationInput = nrow(climateAndStreamflowInput)
		calRows = c(1:(floor(lengthOfCalibrationInput / 3)), ceiling(lengthOfCalibrationInput/(3/2)):lengthOfCalibrationInput)
		valRows = c(ceiling(lengthOfCalibrationInput / 3):floor(lengthOfCalibrationInput/(3/2)))
		
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
			biasAll = rep(NA,minGoodRuns)
		)


		  # since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
		if(any(k1 > k0))	{
			k1[which(k1 > k0)] = k0[which(k1 > k0)] * .99
		}
		if(any(k2 > k1))	{
			k2[which(k2 > k1)] = k1[which(k2 > k1)] * .99
		}
		if(any(uz1 < perc))	{
			uz1[which(uz1 < perc)] = perc[which(uz1 < perc)] * 1.01
		}


		iter = 0
		while(iter < minGoodRuns)	{
			jj = 0
			targetMetricValue = targetMetricValue - 0.01
			while(jj < numberOfRuns & iter < minGoodRuns) {
				jj = jj+1	; print(jj)
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
					
					fwrite(cal_out, paste0(dataOut_location, "calibration_", basinName, ".csv"), append=FALSE)
				}
			}
		}
	}
}





##########################################################################################################
	## Running the Model for Seasonal Forecasts 
seasonalForecast_f = function(
	basinName = 'inster_basin_or_outlet_name',
	climateInputsFileLoc = 'file_location_and_name',	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1)							# 1 for FNF from cal.gov,
	{

	if(file.exists(paste0(dataOut_location, "calibration_", basinName, ".csv")))	{
		# read in the necessary libraries
		library(data.table)
		library(lubridate)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in basin specific data that has been previously loaded
		seas5ClimateInput = readRDS(paste0(climateInputsFileLoc, 'SEAS5.RData')) # reads in a list with each [[i]] being output from a model
		era5ClimateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Recent.RData')))
		
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}

	

		lastHistData = last(era5ClimateInput$Date)
		seas5Rows = (1 + which(as.character(seas5ClimateInput[[1]]$Date) == as.character(lastHistData))):length(seas5ClimateInput[[1]]$Date)
		allForecastsOutput = list()
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
				climateAndStreamflowOutput = historicStreamflow[allHBVoutput, on='Date']
				plot(climateAndStreamflowOutput$Date, climateAndStreamflowOutput$historicQinmm, type='l')
				lines(climateAndStreamflowOutput$Date, climateAndStreamflowOutput$Qg, col='red2')

				if(dataSource == 1)	{
					climateAndStreamflowOutput$historicQinOriginalUnits
					climateAndStreamflowOutput$projectedQinOriginalUnits = climateAndStreamflowOutput$Qg * flowUnitConversion * areaUnitConversion * basinArea
				}	else	{ 
					return("we need to figure out how to read in and normalize this streamflow data")
				}
				
				allForecastsOutput[[iter]] = climateAndStreamflowOutput
			}
		}
	}	else	{print("Ya gotta calibrate the model first you big ole dummy!")}
	return(allForecastsOutput)
}








###################################################################################################################
	## Validation and plot generation
validationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	climateInputsFileLoc = 'file_location_and_name',	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
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
		climateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Historical.RData')))

		
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}

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
	}	else	{print("Ya gotta calibrate the model first you big ole dummy!")}
	return(historicalValidation)
}





reservoirOutflowCalVal_f = function(
	dataOut_location = 'save_file_location',
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









projectionValidationAndPlotGeneration_f = function(
	basinName = 'inster_basin_or_outlet_name',
	climateInputsFileLoc = 'file_location_and_name',	# seas5 / cfs / era5 / Recent .RData is appended in the function
	pathToWatershedsGPKG = 'file_location_and_name.gpkg',
	historicStreamflowFileLoc = 'https://someplace.gov',
	dataOut_location = 'save_file_location',
	dataSource = 1)							# 1 for FNF from cal.gov,
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

		seas5ClimateInputFileNames = list.files(paste0(dataOut_location, basinName, '_seas5MultiOutput'))

		# reading in basin specific data that has been previously loaded
		era5ClimateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Historical.RData')))
			
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
				# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
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
			seas5ClimateInput = readRDS(paste0(dataOut_location, basinName, '_seas5MultiOutput\\', thisSeas5))

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
					plotRows = (lastHistData-365):nrow(climateAndStreamflowOutput)
#					plot(climateAndStreamflowOutput$Date[plotRows], climateAndStreamflowOutput$historicQinmm[plotRows], type='l')
#					lines(climateAndStreamflowOutput$Date[plotRows], climateAndStreamflowOutput$Qg[plotRows], col='red2')
#					abline(v=climateAndStreamflowOutput$Date[lastHistData])

					climateAndStreamflowOutput$month = month(climateAndStreamflowOutput$Date)
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
					data.frame(month = theseMonths, cumsumQ = cumsum(climSubset$totQ[climStart:(climStart+numMonths-1)]), plotMonth = c(1:7)))
			}	
		
			
			
			png(paste0(dataOut_location, basinName, '_projectedOutputFigures\\streamflow_', validSubset$Date[1] - 14, '.png'), width = 1200, height = 720)
			windowsFonts(A = windowsFont("Roboto"))
			par(mar=1.6*c(5,5,2,2), mgp=1.5*c(3,1.3,0), font.lab=2, bty='l', cex.lab=1.5*1.8, cex.axis=1.5*1.4, cex.main=1.5*1.8, col='#1A232F')
			#validSubset$plotDate = rep(1:7, nrow(validSubset)/7) - .2
			boxplot(validSubset$predQ ~ validSubset$Date, pch='.', lwd=1, 
				#col=colorRampPalette(c('#F06000', '#B91863'))(length(unique(validDF$monthsOut))),
				ylim = c(0, max(validDF$predQ)),
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
			text(x=0.5, y=max(validDF$predQ) * 0.95,
				paste0(basinName),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predQ) * 0.87,
				paste0('Forecast Beginning on ', validSubset$Date[1] - 14),
				adj = c(0,0), font=2, col='#1A232F', family='A', cex=1.6*1.6)
			text(x=0.5, y=max(validDF$predQ) * 0.75,
				paste0('CAi Forecast'),
				adj = c(0,0), font=2, col='#A9BF2C', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predQ) * 0.7,
				paste0('Climatology'),
				adj = c(0,0), font=2, col='#FDB600', family='A', cex=1.6*1.3)
			text(x=0.5, y=max(validDF$predQ) * 0.65,
				paste0('Actual'),
				adj = c(0,0), font=2, col='#196CE1', family='A', cex=1.6*1.3)
			dev.off()
		}

	}
}



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
			basinAreaInKm = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
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
					mean(subset(subsetDat, !is.na(stor))$stor),
					mean(subsetDat$AcFtPDay_inflow * 30.4375, na.rm=TRUE),
					mean(subsetDat$resLoss * 30.4375, na.rm=TRUE)))
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
		}



		save_iter = 0
		validDF = data.frame(modelClm = NA, modelHyd = NA, days = NA, monthsOut = NA, predStor = NA, actStor = NA)

		for(thisSeas5 in seas5ClimateInputFileNames)	{
			seas5ClimateInput = readRDS(paste0(dataOut_location, basinName, '_seas5MultiOutput\\', thisSeas5))

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
							predRatio = 1 * (1 / monthsOut)
							outflowCumsum = outflowCumsum + monthAvgResLoss[thisMonth + monthsOut - 1] * (1-predRatio) + 
								predRatio * reservoirVals[[monthsOut]]$Predictions[initialStorRow + monthsOut - 1]

							numDays = nrow(outputSubset)
							predCumsum = predCumsum + mean(outputSubset$QinAcFt) * numDays
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





































		# read in the necessary libraries
		library(data.table)
		library(scales)
		library(lubridate)
		library(sf)				# for geospatial data
		sf::sf_use_s2(FALSE)	# for suppressing issue with intersecting spherical w flat

			# reading in basin specific data that has been previously loaded
		climateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Historical.RData')))

		
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}

		climateAndStreamflowInput = historicStreamflow[climateInput, on='Date']



era5ClimateInput = as.data.table(readRDS(paste0(climateInputsFileLoc, 'Recent.RData')))
		
		calibratedVarsAll = subset(fread(paste0(dataOut_location, "calibration_", basinName, ".csv")), !is.na(kgeAll))
		calibratedVars = tail(calibratedVarsAll[order(calibratedVarsAll$kgeAll), ], 20)
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(abs(calibratedVarsAll$biasAll)), ], 10))
		calibratedVars = rbind(calibratedVars, tail(calibratedVarsAll[order(calibratedVarsAll$nseAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$maeAll), ], 4))
		calibratedVars = rbind(calibratedVars, head(calibratedVarsAll[order(calibratedVarsAll$rmseAll), ], 2))	; rm(calibratedVarsAll)

		historicStreamflow = as.data.table(read.csv(paste0(historicStreamflowFileLoc)))
			# reading in and reformatting streamflow input 
		if(dataSource == 1)	{
			historicStreamflow$Date = ymd(unlist(strsplit(historicStreamflow$DATE.TIME, " "))[seq(1,nrow(historicStreamflow)*2,2)])
			historicStreamflow$historicQinOriginalUnits = as.numeric(historicStreamflow$VALUE)
			basinArea = sum(st_read(paste0(pathToWatershedsGPKG))$SUB_AREA)
			flowUnitConversion = 4.08735e-13 # cubic mm / day in cfs
			areaUnitConversion = (1000000)^2     # sq mm per sq km
			historicStreamflow$historicQinmm = (historicStreamflow$historicQinOriginalUnits / flowUnitConversion) / (areaUnitConversion * basinArea)
		}	else	{ 
			return("we need to figure out how to read in and normalize this streamflow data")
		}

	

		lastHistData = last(era5ClimateInput$Date)
		seas5Rows = (1 + which(as.character(seas5ClimateInput[[1]]$Date) == as.character(lastHistData))):length(seas5ClimateInput[[1]]$Date)
		allForecastsOutput = list()
		iter = 0
		for(numModels in 1:length(seas5ClimateInput))	{
			for(numCalibs in 1:nrow(calibratedVars))	{	
				iter = iter + 1
				climateInput = rbind(era5ClimateInput,seas5ClimateInput[[numModels]][seas5Rows,])


	











	# example input data
basinName = 'FOL_atOutlet'
gageLonLat = 	# for SJF: c(-119.697,37.004) 
				c(-121.155, 38.710)
				# for YRS: c(-121.292500, 39.223611) 
				# for BND: c(-122.185556, 40.288611) 
				# for ISB: c(-118.479, 35.650) 
				# for PNF: c(-119.318, 36.845)
				# for TRM: c(-118.998, 36.416) 
				# for EXC: c(-120.264, 37.591) 
				# for ORO: c(-121.480, 39.540) 
				# for LEW: c(-122.800, 40.732)
				# for SHA: c(-122.417, 40.720)
basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
seas5DataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-seas5.nc'
cfsDataNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-cfs.nc'
era5DataHistoricalNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-hist-era5.nc'
era5DataRecentNCDF = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\cali-recent-era5.nc'
	# to search for Cali reservoirs: https://cdec.water.ca.gov/dynamicapp/wsSensorData
historicStreamflowFileLoc = 
	# for SJF: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SJF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=FOL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for SHA: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=SHA&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for LEW: 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=LEW&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-08-01'
	# for ORO: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for EXC: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for TRM: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=TRM&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for PNF: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=PNF&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for ISB: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ISB&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for BND: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=BND&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-08-01"
	# for YRS: "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=YRS&SensorNums=8&dur_code=D&Start=1906-01-01&End=2022-08-01"
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






