#Get the GFS forecasts for launch day when we lost Jake 8
library(rNOMADS)
library(XML)
library(RCurl)

#GET GFS AT LAUNCH

#Launch date and time, in GMT
forecast.date <- "2014-05-21 10:00:00 UTC"
model.date <- 2014052106
depth <- NULL
forecasts <- GetClosestGFSForecasts(forecast.date, model.date, depth)

model.domain <- c(-81.685181, -74.807739, 36.615528, 33.155948)
pressure <- c(1, 2, 3, 5, 7,
10, 20, 30, 50, 70,
seq(100, 1000, by = 25))
levels <- c(paste(pressure, " mb", sep = ""), "10 m above ground", "80 m above ground", "100 m above ground")
variables <- c("HGT", "TMP", "UGRD", "VGRD", "RH")

#Download back forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$back.forecast,
   levels, variables, model.domain = model.domain, file.name = "launch_back_forecast.grb")

#Download fore forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$fore.forecast,
   levels, variables, model.domain = model.domain, file.name = "launch_fore_forecast.grb")

#GET GFS AT POSSIBLE IMPACT SITE WHEN BALLOON WAS SEEN

forecast.date <- "2014-05-21 23:00:00 UTC" #Assume around 7 PM EST 
model.date <- 2014052118
depth <- NULL
forecasts <- GetClosestGFSForecasts(forecast.date, model.date, depth)

model.domain <- c(-79.6, -76.9, 36.6, 34.4)
pressure <- c(1, 2, 3, 5, 7,
10, 20, 30, 50, 70,
seq(100, 1000, by = 25))
levels <- c(paste(pressure, " mb", sep = ""), "10 m above ground", "80 m above ground", "100 m above ground")
variables <- c("HGT", "TMP", "UGRD", "VGRD", "RH")

#Download back forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$back.forecast,
   levels, variables, model.domain = model.domain, file.name = "gfs_witness_back_forecast.grb")

#Download fore forecast
grib.info <- GribGrab(forecasts$model.url, forecasts$fore.forecast,
   levels, variables, model.domain = model.domain, file.name = "gfs_witness_fore_forecast.grb")


#Get FNL
urls.out <- CrawlModels(abbrev = "fnl", depth = NULL)
url.to.get <- urls.out[which(grepl("20140521", urls.out))]
levels <- c("10 m above ground")
#Download back forecast (at 5 PM EST)
grib.info <- GribGrab(url.to.get, "gdas1.t18z.pgrbf03.grib2",
   levels, variables, model.domain = model.domain, file.name = "fnl_witness_fore_forecast.grb")
#Download fore forecast (at 8 PM EST)
grib.info <- GribGrab(url.to.get, "gdas1.t18z.pgrbf09.grib2",
   levels, variables, model.domain = model.domain, file.name = "fnl_witness_back_forecast.grb")

#Get NAM
urls.out <- CrawlModels(abbrev = "nam", depth = NULL)
url.to.get <- urls.out[which(grepl("20140521", urls.out))]
levels <- c("1000 m above ground", "10 m above ground", "80 m above ground")
#Download back forecast (at 6 PM EST)
grib.info <- GribGrab(url.to.get, "nam.t18z.awphys04.grb2.tm00",
   levels, variables, model.domain = model.domain, file.name = "nam_witness_fore_forecast.grb")
#Download fore forecast (at 8 PM EST)
grib.info <- GribGrab(url.to.get, "nam.t18z.awphys06.grb2.tm00",
   levels, variables, model.domain = model.domain, file.name = "nam_witness_back_forecast.grb")

