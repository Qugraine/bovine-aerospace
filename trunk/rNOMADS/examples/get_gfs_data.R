library(rNOMADS)
urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 2)
model.parameters <- ParseModelPage(urls.out[1])

levels <- c("10 m above ground")
variables <- c("VGRD", "UGRD")
model.domain <- c(-180,180,85,-85)

for(i in 0:24) { #Get the first 24 predictions
    pred <- model.parameters$pred[grep(paste0(sprintf("%02d", i * 3), "$"), model.parameters$pred)] #Get prediction that matches requested forecast
    file.name <- paste(i,'.test',sep="")
    file.name <- GribGrab(urls.out[1], pred, levels, variables,file.name = file.name,model.domain = model.domain)
    model.data <- ReadGrib(file.name, levels, variables) 
    #Order data by variable, then by lon, then by lat
    model.data <- model.data[order(model.data[,3], model.data[,6], model.data[,5]), ]
    var.ind <- dim(model.data)[1] / 2 #Since there are 2 variables, the model.data array is 2x the size of the grid (since each grid point has 2 values)
    lats <- as.numeric(model.data[1:var.ind, 6])
    lons <- as.numeric(model.data[1:var.ind, 5])
    winddir <- 90 - (180 / pi) * atan2(as.numeric(model.data[(var.ind + 1):(dim(model.data)[1]), 7]) , as.numeric(model.data[1:var.ind, 7])) #Wind angle
    winddir[winddir < 0] <-  winddir[winddir < 0] + 360 #Convert to AZIMUTH (degrees from North)
    windspeed <- sqrt(as.numeric(model.data[(var.ind + 1):(dim(model.data)[1]), 7])^2 + as.numeric(model.data[1:var.ind, 7])^2) #Wind speed (m/s)
    timemark <- model.data[1:var.ind,2]
    write.table(cbind(seq_len(var.ind), lats, lons, round(winddir, 2), round(windspeed, 2), timemark), file = paste0(pred, ".txt"), sep = "\t",  #Write the table
        row.names = FALSE, col.names = c("meteoID", "latitude", "longitude", "winddir", "windspeed", "timemark"))
} 
