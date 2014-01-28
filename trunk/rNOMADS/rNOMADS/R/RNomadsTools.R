#Functions for doing specific useful tasks with rNOMADS
AtmosphericProfile <- function(variables, lon, lat, forecast.date, levels = NULL, model.date = "latest", spatial.average = TRUE, temporal.average = TRUE, verbose = TRUE) {
   #Get an atmospheric profile using all available pressure levels, utilizing nearest neighbor interpolation and time weighted averaging if requested
   #INPUTS
   #    VARIABLES - Model variables to get for profile
   #    LON - Longitude
   #    LAT - Latitude
   #    FORECAST.DATE - What date you want the profile for, as a date/time object, in GMT
   #    LEVELS - If not NULL, try and get data for the requested levels only.  If NULL, get all pressure levels
   #    MODEL.DATE - Which model run to use, in YYYYMMDDHH, where HH is 00, 06, 12, 18.  Defaults to "latest", which means get the latest model available.
   #        If MODEL.DATE is not "latest", rNOMADS will scan the entire directory, and this will take time. 
   #    SPATIAL.AVERAGE - Perform nearest neighbor interpolation for 4 grid nodes to get average profile at a specific point.  Default TRUE.  If FALSE, get data from nearest grid node.
   #    TEMPORAL.AVERAGE - Do a weighted average between forecasts to approximate the profile at a specific time. If FALSE, use closest forecast run.   
   #OUTPUTS
   #    PROFILE - Table of pressures and requested variables
   #    FORECAST.USED - Which forecast(s) were used to calculate the profile
   #    NEAREST.MODEL.GRID - The location of the grid node nearest to the requested point

   grid.resolution <- 0.75 #Set this in code since we are using only GFS, may make this an option in the future

   #Get atmosphere
   if(is.null(levels)) {
       pressure <- c(1, 2, 3, 5, 7, 10, 20, 30, 50, 70,
           seq(100, 1000, by = 25))
       levels <- paste(pressure, " mb", sep = "")
   }

   model.to.get <- GetClosestGFSForecasts(forecast.date, model.date, verbose = verbose)

   
   model.domain <- c(lon, lon, lat, lat) + c(-1, 1, 1, -1) * grid.resolution + c(-1, 1, 1, -1) * 0.1 * grid.resolution

   if(!temporal.average) {
      if(abs(model.to.get$back.hr) < abs(model.to.get$fore.hr)) {
          pred <- model.to.get$back.forecast
      } else {
          pred <- model.to.get$fore.forecast
      }
      grib.info <- GribGrab(model.to.get$model.url, pred, levels, variables, model.domain = model.domain)
      grib.data <- ReadGrib(file.path(grib.info$local.dir, grib.info$file.name), levels, variables)
      gridded.data <- ModelGrid(grib.data)
   }

## UNFINISHED

GetClosestGFSForecasts <- function(forecast.date, model.date = "latest", verbose = TRUE) {
   #Figure out the closest GFS forecasts to a given date, returns both the closest forecast behind and the closest forecast ahead, as well as how far beind and how far ahead
   #INPUTS
   #    FORECAST.DATE - What date you want a forecast for, as a date/time object, in GMT
   #    MODEL.DATE - Which model run to use, in YYYYMMDDHH, where HH is 00, 06, 12, 18.  Defaults to "latest", which means get the latest model available.
   #        If MODEL.DATE is not "latest", rNOMADS will scan the entire directory, and this will take time. 
   #    VERBOSE - Give a blow-by-blow account of progress.  Defaults to TRUE.
   #OUTPUTS
   #    MODEL.URL - Which model URL to use
   #    MODEL.RUN.DATE - When the model was run
   #    BACK.FORECAST - Nearest forecast string behind forecast date
   #    FORE.FORECAST - Nearest forecast past forecast date
   #    BACK.HR - How many hours before the requested time the back forecast is
   #    FORE.HR - How many hours after the requested time the forward forecast is

   if(verbose){
       print("Finding model run dates...")
   }
   if(model.date == "latest")  {
       urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 2, verbose = verbose)
       model.parameters <- ParseModelPage(urls.out[1])
       if(length(model.parameters$pred) == 0) { #No data in latest model, try next latest
           model.parameters <- ParseModelPage(urls.out[2]) 
           if(length(model.parameters$pred) == 0) { #Nothing here either, maybe the server is down?
               stop("No data found for most recent GFS model run.  Perhaps the NOMADS server is down?")
           } else {
               url.to.use <- urls.out[2]
           }
       } else {
           url.to.use <- urls.out[1]
       }
    } else {
         urls.out <- CrawlModels(abbrev = "gfs0.5", verbose = verbose)
         model.run.dates <- unlist(stringr::str_extract_all(urls.out, "\\d{10}")) 
         if(model.date %in% model.run.dates) {
              url.to.use <- urls.out[which(model.date == model.run.dates)]
         } else {
             stop(paste0("The model run date ", model.date, " does not appear in the list of model runs on the NOMADS server."))
         }
    }

    if(verbose) {
        print("Determining closest forecasts...")
    }
   
    #Get model run date, convert to POSIX date 
    run.date <- str_match_all(url.to.use, "\\d{10}")[[1]][1,1]
    d.vec <- strsplit(run.date, split = "")[[1]]
    nice.run.date <- strftime(paste0(paste(d.vec[1:4], collapse = ""),
        "-", paste(d.vec[5:6], collapse = ""), "-", paste(d.vec[7:8], collapse = ""),
        " ", paste(d.vec[9:10], collapse = ""), ":00:00", sep = "")) 

   #Figure out time difference between forecast date and model date 
   hr.shift <- as.numeric(difftime(forecast.date, as.POSIXlt(nice.run.date, tz = "GMT")), units = "hours")

   #Figure out the forward (in the future) forecast and back (in the past) forecast using hr.shift
   pred.hrs <- as.numeric(unlist(str_match_all(model.parameters$pred, "\\d{2,3}$")))
   hr.diff <- pred.hrs - hr.shift

   #Get forecasts
   back.forecast <- model.parameters$pred[which(max(hr.diff[which(hr.diff <=0)]) == hr.diff)]
   fore.forecast <- model.parameters$pred[which(min(hr.diff[which(hr.diff > 0)]) == hr.diff)]
   back.hr <- hr.diff[which(max(hr.diff[which(hr.diff <=0)]) == hr.diff)]
   fore.hr <- hr.diff[which(min(hr.diff[which(hr.diff > 0)]) == hr.diff)] 

   if(verbose) {
       print("Finished determining forecast times.")
   }
   return(list(model.url = url.to.use, model.run.date = nice.run.date, back.forecast = back.forecast, fore.forecast = fore.forecast, back.hr = back.hr, fore.hr = fore.hr))
}

ModelGrid <- function(model.data, levels = NULL, variables = NULL, model.domain = NULL) {
    #Transform model data array into a grid with dimensions levels x variables x lon range x lat range
    #This should reduce the size of the returned data by removing redundant information
    #INPUTS
    #    MODEL.DATA - Data returned by ReadGrib
    #    VARIABLES - variables to include in grid, if NULL, include all of them
    #    LEVELS - levels to include in grid, if NULL, include all of them
    #    MODEL.DOMAIN - vector c(LEFT LON, RIGHT LON, TOP LAT, BOTTOM LAT) of region to include in output. If NULL, include everything.
    #
    #OUTPUTS
    #   FCST.GRID - A list with elements:
    #       $Z An array of dimensions levels x variables x lon x lat; each level x variable contains the model grid of data from that variable and level
    #       $X Vector of longitudes
    #       $Y Vector of latitudes
    #       $VARIABLES - the variables contained in the grid
    #       $LEVELS - the levels in the grid
    #       $MODEL.RUN.DATE - when the forecast model was run
    #       $FCST.DATE - what date the forecast is for
  
    model.run.date <- unique(model.data[,1])

    lat.grid <- unique(round(diff(as.numeric(sort(unique(model.data[,6]))))))
    lon.grid <- unique(round(diff(as.numeric(sort(unique(model.data[,5])))))) 

    if(length(model.run.date) > 1) {
        warning("There appears to be more than one model run date in your model grid!")
    }

    fcst.date <- unique(model.data[,2])

    if(length(fcst.date) > 1) {
        warning("There appears to be more than one model run date in your model grid!")
    }

    data.grid <- matrix(as.numeric(model.data[,5:7]), nrow = nrow(model.data))

    if(is.null(variables)) {
        variables <- unique(model.data[,3]) 
    }
  
    nomatch.ind <- is.na(match(variables, unique(model.data[,3])))
    if(sum(nomatch.ind) > 0) {
        warning(paste("The following variables are NOT present in the model data:", paste(variables[nomatch.ind], collapse = " ")))
        variables <- variables[!nomatch.ind]
    }

 
    if(is.null(levels)) {
        levels <- unique(model.data[,4])
    }

    nomatch.ind <- is.na(match(levels, unique(model.data[,4])))
    if(sum(nomatch.ind) > 0) {
        warning(paste("The following levels are NOT present in the model data:", paste(levels[nomatch.ind], collapse = " ")))
        levels <- levels[!nomatch.ind]
    }


    if(is.null(model.domain)) {
        model.domain <- c(min(data.grid[,1]), max(data.grid[,1]), max(data.grid[,2]), min(data.grid[,2]))
    }

    #Build grid

    lons <- as.numeric(sort(unique(model.data[,5])))
    lats <- as.numeric(sort(unique(model.data[,6])))
    grid <- list(x = lons, y = lats)
    
    fcst.grid <- list(z = array(rep(NA, length(lons) * length(lats) * length(variables) * length(levels)),
        dim = c(length(levels), length(variables), length(lons), length(lats))), 
        x = sort(lons), y = sort(lats), variables = variables, levels = levels, 
        model.run.date = model.run.date, fcst.date = fcst.date)

    #Put variables and levels into a series of layered images
    for(lvl in levels) {
        for(var in variables) {
             mi <- which(var == model.data[,3] & lvl == model.data[,4] &
                   data.grid[,1] >= model.domain[1] & data.grid[,1] <= model.domain[2] &
                  data.grid[,2] <= model.domain[3] & data.grid[,2] >= model.domain[4])
             if(length(mi) > 0) {
                 fcst.grid$z[which(lvl == fcst.grid$levels), which(var == fcst.grid$variables),,] <- fields::as.image(
                     array(
                     data.grid[mi,3],
                     dim = c(length(lons), length(lats))),
                     grid = grid,
                     x = cbind(data.grid[,1], data.grid[,2]))$z
              }
        }
    }

    return(fcst.grid)
}
