#Get the GFS forecasts for launch day when we lost Jake 8
library(rNOMADS)

#GET GFS

#Launch date and time, in GMT
forecast.date <- "2014-05-21 10:00:00 UTC"
model.date <- 2014052106
depth <- NULL
forecasts <- GetClosestGFSForecasts(forecast.date, model.date, depth)

model.domain <- c(-81.685181, -74.807739, 36.615528, 33.155948)
pressure <- c(1, 2, 3, 5, 7,
10, 20, 30, 50, 70,
seq(100, 1000, by = 25), "10 m above ground", "80 m above ground", "100 m above ground")
levels <- paste(pressure, " mb", sep = "")
variables <- c("HGT", "TMP", "UGRD", "VGRD", "RH")

#Download back forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$back.forecast,
   levels, variables, file.name = "back_forecast.grb")

#Download fore forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$fore.forecast,
   levels, variables, file.name = "fore_forecast.grb")

