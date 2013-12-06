library(rNOMADS) #Interface with GFS forecast
library(MBA) #Spatial interpolation routines
library(GEOmap) #Map projection
library(R2G2)

BuildAtmosphere <- function(model.date, model.center.point, model.span, model.res, variables, levels) {
    #This function prepares the model grids for interpolation.
    #It will always download a new model when called - so call it sparingly!
    #INPUTS
    #    MODEL.DATE - Date around which to build model
    #    MODEL.CENTER.POINT - Center of model
    #    MODEL.SPAN - Degrees latitude and longitude away from center point
    #    MODEL.RES - Resolution of interpolated grid, in meters
    #    VARIABLES - Variables to get from NOMADS
    #    LEVELS - Levels to get from NOMADS

    #Set up model domain

    model.domain <- c(model.center.point[1] - model.span/2, 
        model.center.point[1] + model.span/2,
        model.center.point[2] + model.span/2,
        model.center.point[2] - model.span/2)
    #Read NOMADS server 
    urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 1)
    model.parameters <- ParseModelPage(urls.out[1])
  
    #Get model run date, convert to POSIX date 
    run.date <- str_match_all(urls.out[1], "\\d{10}")[[1]][1,1]
    d.vec <- strsplit(run.date, split = "")[[1]]
    nice.run.date <- strftime(paste0(paste(d.vec[1:4], collapse = ""), 
        "-", paste(d.vec[5:6], collapse = ""), "-", paste(d.vec[7:8], collapse = ""), 
        " ", paste(d.vec[9:10], collapse = ""), ":00:00", sep = ""))

   
   #Figure out time difference between forecast date and balloon model date 
   hr.shift <- as.numeric(difftime(model.date, as.POSIXlt(nice.run.date, tz = "GMT")))

   #Figure out the forward (in the future) forecast and back (in the past) forecast using hr.shift
   pred.hrs <- as.numeric(unlist(str_match_all(model.parameters$pred, "\\d{2,3}$")))
   hr.diff <- pred.hrs - hr.shift

   #Get forecasts
   back.pred <- model.parameters$pred[which(max(hr.diff[which(hr.diff <=0)]) == hr.diff)]
   fore.pred <- model.parameters$pred[which(min(hr.diff[which(hr.diff > 0)]) == hr.diff)]

   #Get back forecast
   back.file <- GribGrab(urls.out[1], back.pred, levels, variables,
       file.name = "fcst_back.grb", model.domain = model.domain)

   #Get forward forecast
   fore.file <- GribGrab(urls.out[1], fore.pred, levels, variables,
      file.name = "fcst_fore.grb", model.domain = model.domain)

   ## Make into model grid
   back.data <- ReadGrib(back.file, levels, variables)
   fore.data <- ReadGrib(fore.file, levels, variables)

   back.grd <- ModelGrid(back.data)
   fore.grd <- ModelGrid(fore.data)
   back.grd$lon <- back.grd$x
   back.grd$lat <- back.grd$y
   fore.grd$lon <- fore.grd$x
   fore.grd$lat <- fore.grd$y
 
   back.grd.int <- MakeLayerGrids(back.grd, model.res, model.center.point)
   fore.grd.int <- MakeLayerGrids(fore.grd, model.res, model.center.point)
   invisible(list(back.grd.int = back.grd.int, fore.grd.int = fore.grd.int)) 
}

BuildProfile <- function(object.coords, grd.int, weight.avg) {
#This function takes the high resolution back and foreward interpolated grids
   #and produces a single profile (interpolated in space and time) for a given
   #particle position
   #INPUTS
   #    POINT.COORDS - lat/lon coords of point of interest
   #    BACK.GRD.INT - the model data before the particle time
   #    FORE.GRD.INT - the model data after the particle time
   #    WEIGHT.AVG - The weights (time gap) between models
   #
   #OUTPUTS
   #    ATMOS.PROFILE - An atmospheric profile with rows as layers and columns as variables

   #PROJECT LOCATION OF POINT OF INTEREST
   #Note - need to provide cartesian grid here to figure out where object is
   point.xy <- GLOB.XY(object.coords[2], object.coords[1], grd.int[[1]]$projection)

    ##GET PROFILE FOR POINT OF INTEREST

    #Find which node our point is located in by minimizing absolute distance

    abs.x.dist <- abs(point.xy$x - grd.int[[1]]$x)
    abs.y.dist <- abs(point.xy$y - grd.int[[1]]$y)

    #Find index of that node
    x.point.ind <- which(abs.x.dist == min(abs.x.dist))
    y.point.ind <- which(abs.y.dist == min(abs.y.dist))

    #Get back and fore data for that point
    back.profile <- grd.int[[1]]$z[, , x.point.ind, y.point.ind]
    fore.profile <- grd.int[[2]]$z[, , x.point.ind, y.point.ind]

    #Make weighted profile by date
    atmos.profile <- (back.profile * weight.avg[1] + fore.profile * weight.avg[2])/sum(weight.avg)

    atmos.profile <- atmos.profile[ order(atmos.profile[,1]), ]
    invisible(atmos.profile)
}


MakeLayerGrids <- function(fcst.grid, resolution, center.point) {
    #This function takes the ModelGrid output and regrids it at a given resolution projected from a given center.point
    #INPUTS
    #    FCST.GRID - data structure from ModelGrid
    #    RESOLUTION - x/y resolution of new grid, in meters
    #    CENTER.POINT - center point of projection used to convert lat/lon to meters, as c(LON, LAT)
    #OUTPUTS
    #    INTERP.GRID is a similar structure as FCST.GRID, except it is in cartesian coordinates, contains projection information, and has regridded data
    #Set the projection
    proj <- setPROJ(type = 2, LAT0 = center.point[2], LON0 = center.point[1])

    #Convert lat/lon grid to cartesian

    xygrd <- GLOB.XY(fcst.grid$y, fcst.grid$x, proj)

    x.cells <- ceiling(((max(xygrd$y) - min(xygrd$y)) * 1000) / resolution)
    y.cells <- ceiling(((max(xygrd$x) - min(xygrd$x)) * 1000) / resolution)

    x.vec <- vector()
    y.vec <- vector()

    for(k in seq_len(length(xygrd$x))) {
        x.vec <- append(x.vec, rep(xygrd$x[k], 3))
    }

    y.vec <- append(y.vec, rep(xygrd$y, 3))
   #Define new spatially interpolated model structure
    
    interp.grid <- fcst.grid
    interp.grid$center.point <- center.point
    interp.grid$projection <- proj
    interp.grid$z <- array(rep(0, x.cells * y.cells * length(fcst.grid$variables) * length(fcst.grid$levels)),
        dim = c(length(fcst.grid$levels), length(fcst.grid$variables), x.cells, y.cells))

    for(k in seq_len(length(interp.grid$levels))) {
        for(j in seq_len(length(interp.grid$variables))) {
            grid.est <- mba.surf(cbind(x.vec, y.vec, as.vector(fcst.grid$z[k, j, , ])), y.cells, x.cells)$xyz.est
            interp.grid$z[k, j, , ] <- t(array(grid.est$z, dim = c(y.cells, x.cells)))
        }
     }

     interp.grid$x <- grid.est$x
     interp.grid$y <- grid.est$y
     invisible(interp.grid)
}

#Let's fly something! 

#Define model parameters 
model.span <- 1 #degrees to span
model.res <- 2000 #Resolution of interpolated grid
grid.tol <- 10000 #Distance from center point that triggers rebuilding of model (m)
profile.tol <- c(600, 2000) #Rebuild atmospheric profile every X seconds or Y meters of drift from  point
variables <- c("TMP", "HGT", "UGRD", "VGRD")
levels <- paste(c(1, 2, 3, 5, 7, 10, 20, 30, 50, 70, seq(100, 1000, by = 25)), "mb")

#Define launch date and location
model.date <- as.POSIXlt(Sys.time() + 5, tz = "GMT") #Get data for this date
object.coords <- c(-79.39, 35.51, 1000) #Initial coordinates of point of interest
time.limit <- 3600 * 3 #How many seconds to fly

#Get this party started

#BALLOON MODEL PARAMETERS
t <- 0
deltat <- 60 #Time step (seconds)
grid.tol.tmp <- Inf 
tdiff.back <- Inf
tdiff.fore <- Inf
balloon <- list(lat = c(), lon = c(), elev = c(), time = c())
while(t < time.limit) { #Time limit

    #If the particle is at the edge of the model, or we are exiting the time domain
    #Project to max/min of model to test for tolerance, not of "center point"
    if(grid.tol.tmp > grid.tol | (tdiff.back + tdiff.fore) > 3) { 
        print("Downloading model...")
        grd.int <- BuildAtmosphere(model.date, object.coords, model.span, model.res, variables, levels)
        profile.tol.tmp <- c(Inf, Inf)  #Force profile rebuild
        cart.pos <- c(0, 0, object.coords[3])
    }
    
    #Rebuild the atmospheric profiles if necessary
    if(profile.tol[1] < profile.tol.tmp[1] | profile.tol[2] < profile.tol.tmp[2]) {
        print("Reprofiling...")
        profile.tol.tmp <- c(0, 0) #Reset profile
        tdiff.back <-  abs(as.numeric(difftime(as.POSIXlt(grd.int[[1]]$fcst.date, tz = "GMT"), model.date, units = "hours")))
        tdiff.fore <- abs(as.numeric(difftime(as.POSIXlt(grd.int[[2]]$fcst.date, tz = "GMT"), model.date, units = "hours")))
        atmos.profile <- BuildProfile(object.coords, grd.int, c(tdiff.back, tdiff.fore))
        #Interpolate
        i.p   <- splinefun(atmos.profile[,1], rev(as.numeric(unlist(str_match_all(levels, "\\d+")))), method = "natural")
        i.tmp <- splinefun(atmos.profile[,1], atmos.profile[,2], method = "natural")
        i.wu  <- splinefun(atmos.profile[,1], atmos.profile[,3], method = "natural")
        i.wv  <- splinefun(atmos.profile[,1], atmos.profile[,4], method = "natural")
        profile.cart.pos <- cart.pos #So we can determine how far the balloon has gone to trigger reprofiling if necessary
    }
   
   #Distance from center of profile
   cart.pos[1] <- cart.pos[1] + i.wu(cart.pos[3]) * deltat #East - west movement
   cart.pos[2] <- cart.pos[2] + i.wv(cart.pos[3]) * deltat #North - south movement 
   cart.pos[3] <- cart.pos[3] + 10 * t #Elevation gain or loss

   #Latitude and longitude of object
   obj.tmp <- XY.GLOB((cart.pos[1] + i.wu(cart.pos[3]) * deltat)/1000, (cart.pos[2] + i.wv(cart.pos[3]) * deltat) / 1000, grd.int[[1]]$projection)
   object.coords[1] <- obj.tmp$lon
   object.coords[2] <- obj.tmp$lat
   object.coords[3] <- cart.pos[3]

   t <- t + deltat
   model.date <- model.date + deltat
   profile.tol.tmp <- c(profile.tol.tmp[1] + deltat, sqrt((profile.cart.pos[1] - cart.pos[1])^2 + (profile.cart.pos[2] - cart.pos[2])^2))

   #Make sure we are still in the model domain
   model.xy <- GLOB.XY(object.coords[2], object.coords[1], grd.int[[1]]$projection)
   grid.tol.tmp <- sqrt(model.xy$x^2 + model.xy$y^2)

   #Append data to balloon list
    balloon$lat <- append(balloon$lat, object.coords[2])
    balloon$lon <- append(balloon$lon, object.coords[1])
    balloon$elev <- append(balloon$elev, object.coords[3])
    balloon$time <- append(balloon$time, t)
    print(object.coords)
    coords <- cbind(balloon$lon - 360, balloon$lat)
    Dots2GE(coords, balloon$time, goo = "test_trajectory.kml")
}

coords <- cbind(balloon$lon - 360, balloon$lat)
Dots2GE(coords, balloon$time, goo = "test_trajectory.kml")
