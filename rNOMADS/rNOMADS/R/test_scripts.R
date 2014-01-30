source("ReadGrib.R")
source("GetGrib.R")
source("RNomadsTools.R")

library(MBA)
library(stringr)
library(GEOmap)

variables <- c("TMP", "HGT")
lon <- -79.05228
lat <- 35.9079
forecast.date <- as.POSIXlt(Sys.time(), tz = "GMT")

plot(c(200, 300), c(0, 40000), type = "n", xlab = "Temperature (K)", ylab = "Geopotential Height (m)")

prof.nos.not <- AtmosphericProfile(variables, lon, lat, forecast.date, spatial.average = FALSE, temporal.average = FALSE)

points(prof.nos.not$profile.data[,2:1])

prof.s.not <- AtmosphericProfile(variables, lon, lat, forecast.date, spatial.average = TRUE, temporal.average = FALSE)

points(prof.s.not$profile.data[,2:1], pch = 2, col = "red")

prof.nos.t <- AtmosphericProfile(variables, lon, lat, forecast.date, spatial.average = FALSE, temporal.average = TRUE)

points(prof.nos.t$profile.data[,2:1], pch = 3, col = "blue")

prof.s.t <- AtmosphericProfile(variables, lon, lat, forecast.date, spatial.average = TRUE, temporal.average = TRUE)

points(prof.s.t$profile.data[,2:1], pch = 4, col = "green")
