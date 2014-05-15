#Test script for DODS
source("RNomadsTools.R")
source("GetDODS.R")
source("Models.R")

variable <- "tcdcclm"
models.out <- GetDODSDates("gfs_hd")
model.url <- tail(models.out$url, 1)
urls.out <- GetDODSModelRuns(model.url)
model.info <- GetDODSModelRunInfo(model.url, tail(urls.out$model.run, 1))
model.data <- DODSGrab(model.url, tail(urls.out$model.run, 1), variable, 
    c(0, 0), c(0, 719), c(0, 360))
atmos <- ModelGrid(model.data, c(0.5, 0.5), "latlon")
image(atmos$z[1,1,,])
