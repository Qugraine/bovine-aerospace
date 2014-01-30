source("ReadGrib.R")
source("GetGrib.R")
source("RNomadsTools.R")

library(MBA)
library(stringr)
library(GEOmap)

variables <- NULL 
lon <- -79.052029
lat <- 35.907475
forecast.date <- as.POSIXlt(Sys.time(), tz = "GMT")

variables <- unique(append(variables, "HGT", "UGRD", "VGRD"))
AtmosphericProfile(variables, lon, lat, forecast.date, levels = NULL, model.date = "latest", spatial.average = TRUE, temporal.average = TRUE, verbose = TRUE)
