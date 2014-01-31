source("ReadGrib.R")
source("GetGrib.R")
source("RNomadsTools.R")

library(MBA)
library(stringr)
library(GEOmap)

if(FALSE) {
variables <- c()
lon <- -79.052029
lat <- 35.907475
forecast.date <- as.POSIXlt(Sys.time(), tz = "GMT")

variables <- unique(append(variables, c("HGT", "UGRD", "VGRD")))
profile <- AtmosphericProfile(variables, lon, lat, forecast.date, levels = NULL, model.date = "latest", spatial.average = TRUE, temporal.average = TRUE, verbose = TRUE)
}

#Test imager
urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 2, verbose = TRUE)
levels <- c("2 m above ground")
variables <- c("TMP")
model.parameters <- ParseModelPage(urls.out[1])
pred <- model.parameters$pred[1]
#model.domain <- c(-80, -79, 36.5, 35.5)
model.domain <- NULL
grb.info <- GribGrab(urls.out, pred, levels, variables,model.domain = model.domain)
model.data <- ReadGrib(grb.info$file.name, levels, variables)
gridded.data <- ModelGrid(model.data, c(0.5, 0.5), grid.type = "latlon")
image(gridded.data$z[1,1,,])

