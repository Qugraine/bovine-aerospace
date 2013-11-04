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

variables <- c("TMP", "HGT", "UGRD", "VGRD")
hrs.ahead <- 0
model.domain <- c(-80, -79, 36.5, 35.5) #Research triangle region
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

#Get forecasts
back.pred <- model.parameters$pred[which(max(hr.diff[which(hr.diff <=0)]) == hr.diff)]
fore.pred <- model.parameters$pred[which(min(hr.diff[which(hr.diff > 0)]) == hr.diff)]

#Get time weighted average for profiles
weight.avg <- abs(hr.shift - c(as.numeric(str_match_all(back.pred, "\\d+")[[1]][3]),
    as.numeric(str_match_all(fore.pred, "\\d+")[[1]][3])))

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

#Make weighted profile by date
combined.profile <- CombinedProfile(point.coords, back.grd.int, fore.grd.int, weight.avg) 

#Interpolate
i.p   <- splinefun(combined.profile[,1], rev(as.numeric(unlist(str_match_all(levels, "\\d+")))), method = "natural")
i.tmp <- splinefun(combined.profile[,1], combined.profile[,2], method = "natural")
i.wu  <- splinefun(combined.profile[,1], combined.profile[,3], method = "natural")
i.wv  <- splinefun(combined.profile[,1], combined.profile[,4], method = "natural")

#BEGIN THE TRACKING LOOP
s <- 0 #Elapsed seconds
e <- 1000 #Elevation (m)
deltat <- 1 #Time step

point.xy <- GLOB.XY(point.coords[2], point.coords[1], back.grd.int$projection)

p.pos <- c(point.xy$x, point.xy$y) #Position of particle
for(k in seq_len(600)) {
    s <- s + 1
    p.pos[1] <- p.pos[1] + i.wu(e) * deltat
    p.pos[2] <- p.pos[2] + i.wv(e) * deltat
}

point.latlon <- XY.GLOB(p.pos[1]/1000, p.pos[2]/1000, back.grd.int$projection)
