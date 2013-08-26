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

pressure <- c(1, 2, 3, 5, 7,
   10, 20, 30, 50, 70,
   seq(100, 1000, by = 25))

levels <- paste(pressure, " mb", sep = "")
variables <- c("VGRD", "UGRD", "HGT")
hrs.ahead <- 0
start.date <- Sys.time() + 10
model.step <- 600 #Seconds to step
flight.time <- 12 * 3600 #How long to fly
model.domain <- c(-80, -79, 36.5, 35.5) #Research triangle region

#Start off
#Get back forecast
backinfo <- GribGrab(levels, variables, which.fcst = "back", 
   fcst.date = start.date, file.name = "fcst_back.grb",
   model.domain = model.domain)
back.data <- ReadGrib(backinfo$file.name, variables, levels)
back.array <- ModelGrid(back.data, 0.5, 0.5)

#Get forward forecast
foreinfo <- GribGrab(levels, variables, which.fcst = "forward",
   fcst.date = start.date, file.name = "fcst_fore.grb",
   model.domain = model.domain)
fore.data <-ReadGrib(foreinfo$file.name, variables, levels)
fore.array <- ModelGrid(fore.data, 0.5, 0.5)

plot(c(0, 150), c(0, 50000), type = "n", xlab = "Wind Speed (km/hr)", ylab = "Elevation (m)")
tt <- seq(0, 50000, by = 100)

model.date <- start.date
c <- 1
while (TRUE) {
    
    model.date <- model.date + model.step    

    if (model.date > foreinfo$fcst.date) {
        back.array <- fore.array
        foreinfo <- GribGrab(levels, variables, which.fcst = "forward", 
        fcst.date = model.date, file.name = "fcst_fore.grb",
        model.domain = model.domain)
        fore.data <- ReadGrib(foreinfo$file.name, variables, levels)
        fore.array <- ModelGrid(fore.data, 0.5, 0.5)
    }

    tmp.wind <- splinefun(fore.array$z[,1,1,1], sqrt(fore.array$z[,2,1,1]^2 + fore.array$z[,3,1,1]^2), method = "natural")
    synth.wind <- tmp.wind(tt)
    lines(synth.wind * 3.6, tt, col = c)
    c <- c + 1
    if(model.date > start.date + flight.time) {
        print(start.date + flight.time)
        break
    }
}
