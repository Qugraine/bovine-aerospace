#This code is a prototype for the atmospheric profile generator that
#furiousblue will use to track balloon motion
#It will work like this
#1.  Determine forecast time.
#2.  Download the model forecast immediately before and immediately after the forecast time.
#3.  Project from lat/lon into cartesian with origin at the center of the grid.
#4.  Regrid both these forecasts to the resolution required by the balloon model using spatial interpolation
#5.  Determine which grid cell the balloon is located in.
#6.  Create atmospheric profile for that location linearly interpolated in time between the behind and foreward forecasts
#7.  Get data for vertical balloon location on that profile.
#8.  Iterate.

library(rNOMADS) #Interface with GFS forecast
library(MBA) #Spatial interpolation routines
library(GEOmap) #Map projection

levels <- c("1000 mb")
variables <- c("TMP")
hrs.ahead <- 0
fcst.date <- Sys.time() + 3600 * hrs.ahead
model.domain <- c(-80, -79, 36.5, 35.5) #Research triangle region

#Get back forecast
backinfo <- GribGrab(levels, variables, which.fcst = "back", 
   fcst.date = fcst.date, file.name = "fcst_back.grb",
   model.domain = model.domain)

#Get forward forecast
foreinfo <- GribGrab(levels, variables, which.fcst = "forward",
   fcst.date = fcst.date, file.name = "fcst_fore.grb",
   model.domain = model.domain)
