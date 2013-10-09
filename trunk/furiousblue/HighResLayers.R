library(MBA)
library(GEOmap)
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
