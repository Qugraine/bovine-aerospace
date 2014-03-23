library(rNOMADS)
library(GEOmap)
library(aqfig)

save.fig <- TRUE #If TRUE, save as postscript, if FALSE, display

#Get latest wave model run
urls.out <- CrawlModels(abbrev = "wave", depth = 1)

#Get a prediction
model.parameters <- ParseModelPage(urls.out[1])
pred <- model.parameters$pred[grepl("wna", model.parameters$pred)][1]

#Get height at the surface of the sea
levels <- c("surface")
variables <- c("WVHGT")

#Get data
grib.info <- GribGrab(urls.out[1], pred, levels, variables)

#Read data
grib.data <- ReadGrib(grib.info$file.name, levels, variables)

#Use most recent forecast date
wi <- which(grib.data$forecast.date == min(grib.data$forecast.date))

#Subset (maybe I should make a function for this!)

grib.data.sub <- NULL
for(n in names(grib.data)) {
    grib.data.sub[[n]] <- grib.data[[n]][wi]
}

#Format data into array
resolution <- c(0.25, 0.25)
ocean <- ModelGrid(grib.data.sub, resolution)

#Set up display
li <- which(ocean$levels == "surface")
vi <- which(ocean$variables == "WVHGT")

colormap <- rev(rainbow(500, start = 0 , end = 5/6))

if(save.fig) {
   postscript("atlantic_wave_height.eps")
}

image(ocean$x + 180, sort(ocean$y), ocean$z[li,vi,,], col = colormap,
    xlab = "Longitude", ylab = "Latitude",
    main = paste("Atlantic Wave Height (m):", ocean$fcst.date, "GMT"))

#Plot coastlines
plotGEOmap(coastmap, border = "black", add = TRUE,
    MAPcol = "black", shiftlon = 180)

#Plot legend, convert to Celsius
vertical.image.legend(col=colormap, zlim = range(ocean$z[li, vi, ,][!is.na(ocean$z[li,vi,,])]))

if(save.fig) {
   dev.off()
}
