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

variables <- c("TMP", "HGT")
hrs.ahead <- 0
model.domain <- c(-80, -79, 36.5, 35.5) #Research triangle region
lat.res <- 0.5 #Resolution in degrees latitude
lon.res <- 0.5 #Resolution in degrees longitude
cart.res <- 2000 #Cartesian resolution of interpolated data, in meters
center.point <- c(mean(model.domain[1:2]), mean(model.domain[3:4])) #Projection point

point.coords <- c(-79.39, 35.51) #Initial coordinates of point of interest

#Get GFS 0.5 forecast data
urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 1)
model.parameters <- ParseModelPage(urls.out[1])

#Get all pressure levels
levels <- gsub("_", " ", model.parameters$levels[grep("\\d+_mb$", model.parameters$levels)])

#Get model run date, convert to POSIX date 
run.date <- str_match_all(urls.out[1], "\\d{10}")[[1]][1,1]
d.vec <- strsplit(run.date, split = "")[[1]]
nice.run.date <- strftime(paste0(paste(d.vec[1:4], collapse = ""), "-", paste(d.vec[5:6], collapse = ""),
    "-", paste(d.vec[7:8], collapse = ""), " ", paste(d.vec[9:10], collapse = ""), ":00:00", sep = ""))

#Figure out time difference between now and model run date
hr.shift <- as.numeric(difftime(as.POSIXlt(Sys.time(), tz = "GMT"), as.POSIXlt(nice.run.date, tz = "GMT")))

#Figure out the forward (in the future) forecast and back (in the past) forecast using hr.shift
pred.hrs <- as.numeric(unlist(str_match_all(model.parameters$pred, "\\d{2,3}$")))
hr.diff <- pred.hrs - hr.shift

#Get back forecast
back.pred <- model.parameters$pred[which(max(hr.diff[which(hr.diff <=0)]) == hr.diff)]
fore.pred <- model.parameters$pred[which(min(hr.diff[which(hr.diff > 0)]) == hr.diff)]

#Get back forecast
back.file <- GribGrab(urls.out[1], back.pred, levels, variables, 
   file.name = "fcst_back.grb", model.domain = model.domain)

#Get forward forecast
fore.file <- GribGrab(urls.out[1], fore.pred, levels, variables,
   file.name = "fcst_fore.grb", model.domain = model.domain)

## Make into model grid
back.data <- ReadGrib(back.file, levels, variables)
fore.data <- ReadGrib(fore.file, levels, variables)

back.grd <- ModelGrid(back.data)
fore.grd <- ModelGrid(fore.data)

## REGRID FOR HIGHER RESOLUTION
source("HighResLayers.R")
back.grd.int <- MakeLayerGrids(back.grd, cart.res, center.point)
fore.grd.int <- MakeLayerGrids(fore.grd, cart.res, center.point)

#PROJECT LOCATION OF POINT OF INTEREST
#Note - need to provide cartesian grid here to figure out where object is
point.xy <- GLOB.XY(point.coords[2], point.coords[1], back.grd.int$projection)

##GET PROFILE FOR POINT OF INTEREST

#Find which node our point is located in by minimizing absolute distance

abs.x.dist <- abs(point.xy$x - back.grd.int$x)
abs.y.dist <- abs(point.xy$y - back.grd.int$y)

#Find index of that node
x.point.ind <- which(abs.x.dist == min(abs.x.dist))
y.point.ind <- which(abs.y.dist == min(abs.y.dist))

#Get back and fore data for that point
back.profile <- back.grd.int$z[, , x.point.ind, y.point.ind]
fore.profile <- fore.grd.int$z[, , x.point.ind, y.point.ind]
