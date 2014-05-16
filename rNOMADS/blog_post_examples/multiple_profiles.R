library(rNOMADS)
save.fig <- TRUE #If TRUE, save as postscript, if FALSE, display

#Download list of USArray station locations
download.file("http://www.unc.edu/~haksaeng/rNOMADS/myTA.RDATA", destfile = "myTA.RDATA")
load("myTA.RDATA")

#Get variables and levels
variables <- c("TMP", "UGRD", "VGRD", "HGT")
pressure <- c(1, 2, 3, 5, 7,
    10, 20, 30, 50, 70,
    seq(100, 1000, by = 25))
levels <- paste(pressure, " mb", sep = "")

#Get analysis forecast for latest model run
urls.out <- CrawlModels(abbrev = "gfs_hd", depth = 1)
model.parameters <- ParseModelPage(urls.out[1])
pred <- model.parameters$pred[1]

#Get profiles for every station in the network
grid.type <- "latlon"
resolution <- c(0.5, 0.5)

profiles <- RTModelProfile(urls.out[1], pred,
   levels, variables, myTA$lon, myTA$lat, resolution = resolution,
   grid.type = grid.type)

if(save.fig) {
   postscript("usarray_tmp_profiles.eps")
}

##Plot all temperatures
plot(c(-100, 50), c(0, 50000), type = "n", 
    xlab = "Temperature (C)", ylab = "Geopotential Height (m)",
    main = paste("Temperature Profiles above USArray", profiles$model.date))
for(k in seq_len(length(myTA$lon))) {
   synth.hgt <- seq(min(profiles$profile.data[[k]][,1]), 
       max(profiles$profile.data[[k]][,1]), length.out = 1000)
   tmp.spline <- splinefun(profiles$profile.data[[k]][,1], 
       profiles$profile.data[[k]][,2] - 273.15, method = "natural")
   synth.tmp <- tmp.spline(synth.hgt)
   lines(synth.tmp, synth.hgt)
}

if(save.fig) {
   dev.off()
}

if(save.fig) {
   postscript("usarray_wind_profiles.eps")
}

## Plot all winds
plot(c(-250, 250), c(0, 50000), type = "n",
    xlab = "Wind Speed (km/hr)", ylab = "Geopotential Height (m)",
    main = paste("Wind Profiles above USArray", profiles$model.date))
for(k in seq_len(length(myTA$lon))) {
   synth.hgt <- seq(min(profiles$profile.data[[k]][,1]), 
       max(profiles$profile.data[[k]][,1]), length.out = 1000)
   zone.spline <- splinefun(profiles$profile.data[[k]][,1], 
       profiles$profile.data[[k]][,3] * 3.6, method = "natural")
   synth.zone <- zone.spline(synth.hgt)
   lines(synth.zone, synth.hgt, col = "red")
   merid.spline <- splinefun(profiles$profile.data[[k]][,1],
       profiles$profile.data[[k]][,4] * 3.6, method = "natural")
   synth.merid <- merid.spline(synth.hgt)
   lines(synth.merid, synth.hgt, col = "blue")
}

legend("topleft", col = c("red", "blue"),
   pch = c(NA, NA), lty = c(1, 1), 
   legend = c("Zonal Wind", "Meridional Wind"))
if(save.fig) {
   dev.off()
}

