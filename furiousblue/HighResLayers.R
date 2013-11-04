library(MBA)
library(GEOmap)

CombinedProfile <- function(point.coords, back.grd.int, fore.grd.int, weight.avg) {
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
   #    COMBINED.PROFILE - An atmospheric profile with rows as layers and columns as variables

   #PROJECT LOCATION OF POINT OF INTEREST
   #Note - need to provide cartesian grid here to figure out where object is
   point.xy <- GLOB.XY(point.coords[2], point.coords[1], back.grd.int$projection)

    ##GET PROFILE FOR POINT OF INTEREST

    #Find which node our point is located in by minimizing absolute distance

    abs.x.dist <- abs(point.xy$x - back.grd.int$x)
    abs.y.dist <- abs(point.xy$y - back.grd.int$y)

    #Find index of that node
    x.point.ind <- which(abs.x.dist == min(abs.x.dist))
    y.point.ind <- which(abs.y.dist == min(abs.y.dist))

    #Get back and fore data for that point
    back.profile <- back.grd.int$z[, , x.point.ind, y.point.ind]
    fore.profile <- fore.grd.int$z[, , x.point.ind, y.point.ind]

    #Make weighted profile by date
    combined.profile <- (back.profile * weight.avg[1] + fore.profile * weight.avg[2])/sum(weight.avg)
   
    combined.profile <- combined.profile[ order(combined.profile[,1]), ] 
    invisible(combined.profile)
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
