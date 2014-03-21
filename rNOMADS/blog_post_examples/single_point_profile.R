## MAKE TEMP AND WIND PROFILES FOR A SPECIFIC POINT AT A SPECIFIC TIME

library(rNOMADS)

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
