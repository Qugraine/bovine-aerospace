## MAKE TEMP AND WIND PROFILES FOR A SPECIFIC POINT AT A SPECIFIC TIME

library(rNOMADS)

save.fig <- TRUE #If TRUE, save as postscript, if FALSE, display

#Get the date and time for right now
#Need to convert it to UTM/GMT

forecast.date <- as.POSIXlt(Sys.time(), tz = "GMT")

#Now figure out what forecasts are closest to this time
#This only works for the Global Forecast System
#but others can be added if users request them

forecasts <- GetClosestGFSForecasts(forecast.date, model.date = "latest")

#We get data from before and after, then do a weighted average
#to determine temperature and wind speed

#Get levels
pressure <- c(1, 2, 3, 5, 7,
    10, 20, 30, 50, 70,
    seq(100, 1000, by = 25))
levels <- paste(pressure, " mb", sep = "")

#Variables - temperature and height only
variables <- c("TMP", "HGT", "UGRD", "VGRD")

#Location - Sakura-jima volcano, Japan
#If it erupts right now, where will the ash go?
#Might be nice to know!

lon <- c(130.655115)
lat <- c(31.590929)
grid.type <- "latlon"
resolution <- c(0.5, 0.5)

#Get profiles
back.profile <- RTModelProfile(forecasts$model.url, forecasts$back.forecast,
    levels, variables, lon, lat, resolution = resolution, grid.type = grid.type)
fore.profile <- RTModelProfile(forecasts$model.url, forecasts$fore.forecast,
    levels, variables, lon, lat, resolution = resolution, grid.type = grid.type)

#Get weighted averages
total.time <- abs(forecasts$fore.hr) + abs(forecasts$back.hr)
fore.wt <- abs(forecasts$fore.hr)/total.time
back.wt <- abs(forecasts$back.hr)/total.time

#Final profile
profile.array <- (back.profile$profile.data[[1]] * back.wt + fore.profile$profile.data[[1]] * fore.wt)

#Make nice spline curves
tmp.spline <- splinefun(profile.array[,1], profile.array[,2] - 273.15, method = "natural")
zone.spline <- splinefun(profile.array[,1], profile.array[,3] * 3.6, method = "natural")
merid.spline <- splinefun(profile.array[,1], profile.array[,4] * 3.6, method = "natural")
synth.hgt <- seq(min(profile.array[,1]), max(profile.array[,1]), length.out = 1000)
synth.tmp <- tmp.spline(synth.hgt)
synth.zone <- zone.spline(synth.hgt)
synth.merid <- merid.spline(synth.hgt)

#Make plots
#Temperature

if(save.fig) {
   postscript("tmp_profile.eps")
}

plot(synth.tmp, synth.hgt, type = "l",
   xlab = "Temperature (C)", ylab = "Geopotential Height (m)", 
   main = paste("Temperature Profile at Sakura-Jima volcano,", forecast.date, "GMT"))
points(profile.array[,2] - 273.15, profile.array[,1], col = "red")

if(save.fig) {
   dev.off()
}

#Zonal/Meridional wind speed
wind.range <- c(min(cbind(synth.zone, synth.merid)),
   max(cbind(synth.zone, synth.merid)))

if(save.fig) {
   postscript("wind_profile.eps")
}

plot(wind.range, c(min(synth.hgt), max(synth.hgt)), type = "n",
   xlab = "Wind Speed (km/hr)", ylab = "Geopotential Height (m)",
   main = paste("Zonal and Meridional Wind Speeds at Sakura-Jima volcano,", forecast.date, "GMT"))

points(profile.array[,3] * 3.6, profile.array[,1], col = "red")
lines(synth.zone, synth.hgt, col = "red")
points(profile.array[,4] * 3.6, profile.array[,1], pch = 2, col = "blue")
lines(synth.merid, synth.hgt, col = "blue")

legend("topright", pch = c(1, 2), col = c("red", "blue"), lty = c(1, 1),
   legend = c("Zonal Wind", "Meridional Wind"))

if(save.fig) {
   dev.off()
}

