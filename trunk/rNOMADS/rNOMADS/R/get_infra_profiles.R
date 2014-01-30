source("ReadGrib.R")
source("GetGrib.R")
source("RNomadsTools.R")

library(MBA)
library(stringr)
library(GEOmap)

load("myTA.RDATA")

variables <- c("TMP", "UGRD", "VGRD", "HGT")
lon <- -79.05228
lat <- 35.9079
forecast.date <- as.POSIXlt(Sys.time(), tz = "GMT")

for(k in seq_len(length(myTA$lat))) {
    prof <- AtmosphericProfile(variables, myTA$lon[k], myTA$lat[k], forecast.date, 
        spatial.average = FALSE, temporal.average = FALSE)
    gfs.data <- cbind(as.numeric(unlist(str_extract_all(prof$levels, "\\d+"))) * 100, prof$profile.data)
    atmos.prof <- list(data = gfs.data, 
        variables = c("Temp (K)", "Zonal Wind (m/s)", "Meridional Wind (m/s)", "Geopotential Height (m)"), forecast.date = forecast.date, lat = lat, lon = lon)
    save(file = paste0(myTA$name[k], ".RData"), atmos.prof)
}
