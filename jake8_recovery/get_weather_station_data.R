library(weatherData)
data <- getWeatherForDate("KHRJ", "2014-05-21", opt_detailed=TRUE, 
    opt_all_columns=TRUE)
ind <- 49:61
write("WEATHER DATA FOR COATS, NC, MAY 21 2014\n", file = "coats_weather.txt")
write("TIME (EDT)\tWIND DIR (deg)\tWIND SPEED (mph)\n", file = "coats_weather.txt", append = TRUE)
for(k in ind) {
   write(paste0(data$TimeEDT[k], "\t", 360 - data$WindDirDegrees[k], "\t", data$Wind_SpeedMPH[k], "\n"), file = "coats_weather.txt", append = TRUE)
}
