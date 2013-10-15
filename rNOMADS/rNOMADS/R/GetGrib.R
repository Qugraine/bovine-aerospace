CrawlModels <- function(abbrev = NULL, url = NULL, depth = 1000) {
   #A simple web crawler that looks at the specified model directory online and gets information on all runs of the specified model.
   #See the NOMADSList function for available models.
   #Alternatively, pass CrawlModels a URL to get a model that I have not included yet.
   #INPUTS
   #    ABBREV - Model abbreviation as defined in NOMADSList().  
   #        If NULL, use the url you provided, if you did not provide one, throw error.
   #    URL - Use your own URL and attempt to get model data from it.  
   #        This is in case NOMADS updates its system before I have a chance to update rNOMADS
   #    DEPTH - How many links to return; this prevents infinite loops if something goes wrong
   #OUTPUTS
   #    URLS.OUT is a list of available models from the given ABBREV or URL
    
}
GetModelRunHour <- function(model.date = Sys.time(), fcst.date = Sys.time(),
    url.to.check = c("http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl?dir=%2Fgfs.", "%2Fmaster"), attempts = 10)
{
    #Checks for and returns the date for the latest model run time for the requested date.
    #By default, check for the system time, and get the closest forecast.
    #If the input is not the current time, get the model forecast closest behind the requested date.
    #
    #INPUTS
    #    MODEL.DATE - a date in POSIXct format saying which GFS model date to use (i.e. when the model was run).  Defaults to the current system time
    #    FCST.DATE = a date in POSIXct format saying what date the forecast should be for.  Defaults to the current system time.
    #    URL.TO.CHECK - what URL to append the GFS formatted model run time to
    #    We use this URL to check if the model has run yet.
    #    Perhaps, if the default fails, the user can modify this argument to make things work
    #    ATTEMPTS - number of model runs to check for before giving up
    #
    #OUTPUTS as list
    #    MODEL.HOUR - The model run hour to download
    #    MODEL.DATE - the full date of the model run
    #    URL.TESTED - The url that was tested to determine the model hour
    #    FCST.TDIFF - Time difference between model date and forecast date (i.e. how far in the future the forecast is from the model run that's available) in hours 
    #    FCST.BACK - The model forecast run immediately before the requested forecast date, in hours, in case that grib file is desired
    #    FCST.FORE - The model forecast run immediately after the requested forecast date, in hours, in case that grib file is desired
    

    model.hour <- seq(0, 18, by = 6)
    fcst.hour <- seq(0, 192, by = 3)
    #Convert to GMT
    model.date <- as.POSIXlt(model.date, tz = "GMT")
    fcst.date <- as.POSIXlt(fcst.date, tz = "GMT")
    
    c = 1
     
    while (1)
    {
       yr <- model.date$year + 1900
       mo <- model.date$mo + 1
       mday <- model.date$mday
       hr <- model.date$hour
       
       hr.diff <- model.hour - hr
       latest.model.run <- model.hour[hr.diff == max(hr.diff[hr.diff <= 0])]
       model.date <- as.POSIXlt(ISOdatetime(yr, mo, mday, latest.model.run, 0, 0, tz = "GMT"))   
       fcst.url <- paste(url.to.check[1], yr, sprintf("%02d", mo), sprintf("%02d", mday), sprintf("%02d", latest.model.run), url.to.check[2], sep = "")
       test <- suppressWarnings(tryCatch(url(fcst.url, open = "rb"), error = NoModelRun))
       
       if(test == "Failure") {
           model.date = as.POSIXlt(model.date - 3600 * 6) #Subtract 6 hours and try again
       }
       else {
           close(test)
           fcst.tdiff <- as.numeric(difftime(fcst.date, model.date, units = "hours"))
           fcst.hour.diff <- fcst.hour - fcst.tdiff
           fcst.back <- fcst.hour[fcst.hour.diff == max(fcst.hour.diff[fcst.hour.diff <=0])]
           fcst.fore <- fcst.hour[fcst.hour.diff == min(fcst.hour.diff[fcst.hour.diff >=0])]
           break
       }

      if (c > attempts) {
          model.hour <- NA
          break
      } 
      c <- c + 1
   } 
   return (list(model.hour = latest.model.run, model.date = as.POSIXlt(ISOdatetime(yr, mo, mday, latest.model.run, 0, 0, tz = "GMT")), 
        url.tested = fcst.url, fcst.tdiff = fcst.tdiff, fcst.back = fcst.back, fcst.fore = fcst.fore))
}

GribGrab <- function(levels, variables, which.fcst = "back", local.dir = ".", file.name = "fcst.grb", model.date = Sys.time(), fcst.date = Sys.time(), 
    model.domain = NULL, tidy = FALSE, verbose = TRUE)
{
    #Get grib file from the GFS forecast repository
    #INPUTS
    #    LEVELS is the vertical region to return data for,  as vector
    #    VARIABLES is the data to return, as vector
    #    WHICH.FCST determines if you want the forecast BEFORE the requested time ("back") or AFTER the requested time ("forward"), defaults to "back"
    #    LOCAL.DIR is the directory to save the files in
    #    FILE.NAME is the directory path and file name to save the grib file on disk, defaults to "fcst.grb" in current directory
    #    MODEL.DATE is the date and time of the requested model run, in GMT and POSIXlt format, defaults to current system time
    #    FCST.DATE is the requested forecast date, defaults to current system time
    #    MODEL.DOMAIN is a vector of latitudes and longitudes that specify the area to return a forecast for
    #    This is a rectangle with elements: west longitude, east longitude, north latitude, south latitude
    #    Defaults to entire planet
    #    LEVELS is the vertical region to return data for
    #    VARIABLES is the data to return
    #    TIDY asks whether to delete all grib files in the directory specified in FILE.NAME, default FALSE.
    #    This is useful to clear out previous model runs.
    #    It looks for all files named '.grb' and removes them.
    #    VERBOSE gives a blow by blow account of the download. Default TRUE.
    #OUTPUTS
    #    FILE.NAME is the name and location of the grib file that has been downloaded 

   if(tidy) {
        unlink(list.files(local.dir, pattern = "*\\.grb$"))
   }
   levels.str <- paste(gsub(" ", "_", levels), collapse = "=on&lev_")
   variables.str <- paste(variables, collapse = "=on&var_")

   #Convert dates to GMT
   
    model.date <- as.POSIXlt(model.date, tz = "GMT")
    fcst.date <- as.POSIXlt(fcst.date, tz = "GMT")

   #Check for latest model run date
   model.params <- GetModelRunHour(model.date = model.date, fcst.date = fcst.date) 
   if(is.na(model.params$model.hour)) {
       stop("Could not find the latest model run date.  Make sure you have a working Internet connection.  If you do and this code is still not working, it may be that the NOMADS website is down.
           Give it a an hour or so, and try again.")
   }
   if(model.params$model.date > fcst.date) {
      stop("The reqested model date is after the requested forecast date! This means you are trying to access a forecast from a model that has not been run yet.")
   }
   if(which.fcst == "back") {
       fcst.date <- as.POSIXlt(model.params$model.date + 3600 * model.params$fcst.back)
       grb.name <- paste("gfs.t", sprintf("%02d", model.params$model.hour), "z.mastergrb2f",
       sprintf("%02d", model.params$fcst.back), "&", sep = "")
   } else if (which.fcst == "forward") { 
       fcst.date <- as.POSIXlt(model.params$model.date + 3600 * model.params$fcst.fore)
       grb.name <- paste("gfs.t", sprintf("%02d", model.params$model.hour), "z.mastergrb2f",
       sprintf("%02d", model.params$fcst.fore), "&", sep = "")
   } else {
       stop(paste("Did not recognize the forecast designation ", 
       which.fcst, 
       ".  Please use either \"back\" for the nearest forecast time BEFORE the requested time, 
       or \"forward\" for the nearest forecast time AFTER the requested time.", sep = ""))
   }
   grb.dir <- paste("dir=", strsplit(model.params$url.tested, split = "dir=")[[1]][2], sep = "")

   if(!is.null(model.domain)) {
       subregion.str <- paste( "=on&subregion=",
       "&leftlon=", model.domain[1],
       "&rightlon=", model.domain[2],
       "&toplat=", model.domain[3],
       "&bottomlat=", model.domain[4],
       "&", sep = "")
    } else {
       subregion.str <- "=on&" 
    }

   grb.url <- paste("http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl?file=",
       grb.name,
       "lev_",
       levels.str,
       "=on&var_",
       variables.str,
       subregion.str,
       grb.dir,
       sep = "") 

      #now write download logic

   download.file(grb.url, paste(local.dir,file.name, sep = "/"), mode = "wb", quiet = !verbose)
   
   return(list(file.name = file.name, model.date = model.params$model.date, fcst.date = fcst.date))
}

NOMADSList <- function(abbrev = NULL) {
    #Returns a list of model abbreviations, a short description, and URL for each model offered by the NOMADS server
    #If a specific model abbreviation is requested, the abbreviation is checked against the model list.
    #If a match is found, information is returned about that model; otherwise an error occurs
    #INPUTS
    #    ABBREV is the model abbreviation that rNOMADS uses to figure out which model you want.
    #    if NULL, returns information on all models
    #OUTPUTS
    #    MODEL.LIST - a list of model metadata with elements
    #        $ABBREV - the abbrevation used to call the model in rNOMADS
    #        $NAME - the name of the model
    #        $URL - the location of the model on the NOMADS website

    abbrevs <- c(
    "fnl",
    "gfs1.0",
    "gfs0.5",
    "gfs2.5",
    "gfse_highres",
    "gfse_precip_biasc",
    "gfse_highres_biasc", 
    "gfse_ndgdres_biasc",
    "naefs_hires_biasc",
    "naefs_ndgdres_biasc",
    "ngac2d",
    "ngac3d",
    "ngac_aod",
    "aqm_dm",
    "aqm_hso",
    "hires_ak",
    "hires_econus",
    "hires_guam",
    "hires_hi",
    "hires_pr",
    "hires_wconus",
    "nam12_ak",
    "nam12_conus",
    "nam12_na",
    "nam12_carib_ca",
    "nam12_pa",
    "nam_nest_ak",
    "nam_nest_conus",
    "nam_nest_hi",
    "nam_nest_pr",
    "rtma_ak",
    "rtma_conus",
    "rtma_conus2.5",
    "rtma_guam",
    "rtma_hi",
    "rtma_pr",
    "rap",
    "rap_na32",
    "narre",
    "sref_conus",
    "sref_conus_bc",
    "sref_na32",
    "sref_na16",
    "rtofs_at",
    "rtofs_at_hires",
    "sea_ice",
    "wave",
    "gl_wave",
    "wave_mgrd",
    "wave_hur",
    "wave_nfc",
    "estofs",
    "cmc_en",
    "fnmoc_en")

    names <- c(
    "Final Operational Global Forecast System ",
    "Global Forecast System 1x1 Degree ",
    "Global Forecast System 0.5x0.5 Degree ",
    "Global Forecast System 2.5x2.5 Degree ",
    "Global Forecast System Ensemble",
    "Global Forecast System Ensemble Precipitation Bias Corrected ",
    "Global Forecast System Ensemble Bias Corrected ",
    "Global Forecast System Ensemble National Digital Guidance Database Bias Corrected ",
    "North American Ensemble Forecast System Bias Corrected ",
    "North American Ensemble Forecast System National Digital Guidance Database Bias Corrected ",
    "NOAA Environmental ing System Global Forecast System Aerosol Component 2D",
    "NOAA Environmental ing System Global Forecast System Aerosol Component 3D",
    "NOAA Environmental ing System Global Forecast System Aerosol Optical Depth",
    "Air Quality Daily Maximum",
    "Air Quality Hourly Surface Ozone",
    "High Res Window Alaska",
    "High Res Window - East Continental United States",
    "High Res Window - Guam",
    "High Res Window - Hawaii",
    "High Res Window - Puerto Rico",
    "High Res Window - West Continental United States", 
    "North American Mesoscale 12 km - Alaska",
    "North American Mesoscale 12 km - Continental United States",
    "North American Mesoscale 12 km - North America",
    "North American Mesoscale 12 km - Caribbean and Central America", 
    "North American Mesoscale 12 km - Pacific",
    "North American Mesoscale Nest - Alaska",
    "North American Mesoscale Nest - Continental United States",
    "North American Mesoscale Nest - Hawaii",
    "North American Mesoscale Nest - Puerto Rico",
    "Real-Time Mesoscale Analysis - Alaska",
    "Real Time Mesoscale Analysis - Continental United States",
    "Real Time Mesoscale Analysis - Continental United States 2.5 km Resolution",
    "Real Time Mesoscale Analysis - Guam",
    "Real Time Mesoscale Analysis - Hawaii",
    "Real Time Mesoscale Analysis - Puerto Rico",
    "Rapid Refresh Weather Prediction System",
    "Rapid Refresh Weather Prediction System - 32 km Resolution",
    "North American Rapid Refresh Ensemble",
    "Short Range Ensemble Forecast - Continental United States 40 km",
    "Short Range Ensemble Forecast - Continental United States 40 km Bias Corrected",
    "Short Range Ensemble Forecast - Continental United States 32 km",
    "Short Range Ensemble Forecast - Continental United States 16 km",
    "Real Time Ocean Forecast System - Atlantic",
    "Real Time Ocean Forecast System - Atlantic High Resolution",
    "Sea Ice",
    "Operational Ocean Wave Predictions",
    "Opreational Ocean Wave Predictions - Great Lakes",
    "Multi-grid Wave",
    "Hurricane Wave",
    "Combined Wave Ensemble",
    "Extratropical Surge and Tide Operational Forecast System",
    "Canadian Meterological Center Global Ensemble",
    "Fleet Numerical Meteorology and Oceanography Ensemble Forecast System") 

    urls <- c(
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_fnl.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_2p5.pl", 
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gens.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gensbc_precip.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gensbc.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gensbc_ndgd.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_naefsbc.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_naefsbc_ndgd.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_ngac_a2d.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_ngac_a3d.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_ngac_aod.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_aqm_daily.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_aqm_ozone_1hr.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hiresak.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hireseast.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hiresguam.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hireshi.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hirespr.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hireswest.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_ak.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_na.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_crb.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_pac.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_alaskanest.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl", 
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_hawaiinest.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_priconest.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_akrtma.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_rtma.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_rtma2p5.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_gurtma.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hirtma.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_prrtma.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_rap.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_rap32.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_narre.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_sref.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_srefbc.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_sref_na.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_sref_132.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_ofs.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_rtofs_hires.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_seaice.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_wave.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_glw.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_wave_multi.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_hurricane_wave.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_nfcens.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_estofs.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_cmcens.pl",
    "http://nomads.ncep.noaa.gov/cgi-bin/filter_fens.pl")

    if(!is.null(abbrev)) {
        i <- which(abbrevs == abbrev)
        if(length(i) == 0) {
            stop(paste("The model you searched for:\"", abbrev, "\"is not included in rNOMADS.  Sorry!"))
        } else {
             return(list(abbrev = abbrev, name = names[i], url = urls[i]))
        }
    }
    
    if(display) {
        cat(paste(abbrevs
    return(list(abbrevs = abbrevs, names = names, urls = urls))
}
NoModelRun <- function(e)
{
    #Called when code in GetModelRunDate tries to ping a GFS model that has not been run yet
    return ("Failure")
}

RecursiveWebCrawler <- function(url, url.list = c(), start.depth = 0, max.depth = 1000) {
#    This function recursively searches for links in the given url and follows every single link, to a maximum of DEPTH.
#    It returns a list of the final (dead end) URLs.
#    Many thanks to users David F and Adam Smith on stackoverflow for the link parser:
#    http://stackoverflow.com/questions/3746256/extract-links-from-webpage-using-r/3746290#3746290
#    INPUTS
#        URL is the url to start looking in
#        START.DEPTH is the current number of links that have been searched (useful for recursion)
#        MAX.DEPTH is the maximum number of links to follow; stops infinite loops
#    OUTPUTS
#        URLS.OUT are the URLs at the end of the road
#        DEPTH is the number of URLS that have been found so far

    doc <- htmlParse(url)
    links <- xpathSApply(doc, "//a/@href")
    free(doc)
    if(is.null(links)) {
        print(url)
        url.list <- append(url.list, url)
    } else {
        for(link in links) {
            url <- RecursiveWebCrawler(link, url.list, start.depth = start.depth, max.depth = max.depth)
        }
    }
}

