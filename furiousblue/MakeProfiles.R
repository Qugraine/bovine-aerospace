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

levels <- c("1000 mb", "100 mb")
variables <- c("TMP", "UGRD")
hrs.ahead <- 0
fcst.date <- Sys.time() + 3600 * hrs.ahead
model.domain <- c(-80, -79, 36.5, 35.5) #Research triangle region
lat.res <- 0.5 #Resolution in degrees latitude
lon.res <- 0.5 #Resolution in degrees longitude
cart.res <- 2000 #Cartesian resolution of interpolated data, in meters
center.point <- c(mean(model.domain[1:2]), mean(model.domain[3:4])) #Projection point

point.coords <- c(-79.39, 35.51) #Initial coordinates of point of interest

#Get back forecast
backinfo <- GribGrab(levels, variables, which.fcst = "back", 
   fcst.date = fcst.date, file.name = "fcst_back.grb",
   model.domain = model.domain)

#Get forward forecast
foreinfo <- GribGrab(levels, variables, which.fcst = "forward",
   fcst.date = fcst.date, file.name = "fcst_fore.grb",
   model.domain = model.domain)

## Make into model grid
back.data <- ReadGrib(backinfo$file.name, variables, levels)
fore.data <- ReadGrib(foreinfo$file.name, variables, levels)

back.grd <- ModelGrid(back.data, lat.res, lon.res)
fore.grd <- ModelGrid(fore.data, lat.res, lon.res)

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
