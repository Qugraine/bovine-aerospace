library(rNOMADS) #Interface with GFS forecast
library(MBA) #Spatial interpolation routines
library(GEOmap) #Map projection
library(R2G2) #Write KML

AscentCd <- function(Re = NULL, formula = "gallice") {
    #Determine the coefficient of drag of an ascending balloon.
    #This can be done in many ways, some of which take Reynold's number into account.
    #INPUTS
    #    RE is Reynold's number
    #    formula is the formula to use to calculate the drag coefficient (see notes below)
    #OUTPUTS
    #    CD is the coefficient of drag    

    #smoothsphere taken from http://www.engineeringtoolbox.com/drag-coefficient-d_627.html

    #gallice taken from Gallice et al (2011) "Modeling the ascent of sounding balloons: derivation of the vertical air motion" Atmospheric Measuring Techniques (4) 2235-2253
  #They note it is only strictly valid for a certain type of weather balloon at night

    #mikhailov taken from Mikhailov and Freire (2013) "The drag coefficient of a sphere: An approximation using Shanks transform", Powder Technology 237 432-435

    Cd <- NULL
    if(formula == "smoothsphere") { #From Engineering Toolbox website
        Cd <- 0.5
    }

    if(formula == "gallice")  {
        Cd <- 0.04808 * (log(Re))^(2) - 1.406 * log(Re) + 10.49
    }

    if(formula == "mikhailov") {
        Cd <- (777/(646 * Re)) * ((669806/875) + (114976/1155) * Re + (707/1380)*Re^2)/((32869/952) + (924/643) * Re + (1/385718) * Re^2)
    }

    if(formula == "smoothsphere" | is.null(Cd)) { #From Engineering Toolbox website
        Cd <- 0.5
    }
    
    invisible(Cd)
}

AscentRe <- function(ma, pa, Vb, T, va, r) {
    #Calculate Reynold's number for ascending balloon
    #INPUTS
    #   MA is molar mass of air (kg) 
    #   PA is air pressure (pa)
    #   T is temperature (K)
    #   VB is the volume of the balloon (m3)
    #   VA is ascent velocity (m/s)
    #   R is balloon radius (m)
    #OUTPUTS
    #   RE is the Reynold's number of the balloon

    #Calculate dynamic viscosity of air using Sutherland's formula

    mu0 <- 0.00001827 #Reference viscosity (pa s)
    To <- 291.15 #Reference temperature (K)
    C <- 120 #Sutherland's constant for air
    R <- 8.3145 #Universal gas constant
 
    dm <- (pa/(R*T))*ma #Air density

    mu <- mu0 * ((To + C)/(T + C)) * (T/To)^(1.5)

    invisible(dm*va*r*2.0/mu) #Reynold's number
}

AscentVelocity <- function(mi, Mp, mp, ma, pa, po, acd, max.radius, T.ambient, T.diff = 0, T.interior = NULL, sealed = TRUE) {
    #Calculate the ascent velocity of the balloon.
    #INPUTS
    #    MI is mass of instrument package, parachute, etc (kg)
    #    Mp is mols of lift gas
    #    MP is molar mass of lift gas (kg)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    PA is ambient air pressure (pascals)
    #    PO is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
    #    ACD is coefficient of drag (unitless)
    #    MAX.RADIUS is the radius of a non elastic envelope or the burst radius of an elastic envelope
    #    T.AMBIENT is temperature outside of balloon (kelvin)
    #    T.DIFF is temperature difference between balloon and outside (kelvin)
    #    T.INTERIOR is absolute interior temperature (overrides T.DIFF) (kelvin)
    #    SEALED denotes if the envelope is open to the atmosphere (and thus can lose gas) or is sealed shut (gas is retained)
    #OUTPUTS
    #    VA is ascent velocity (m/s)

    if(T.diff != 0 & !is.null(T.interior)) {
         warning("A temperature difference from ambient and a fixed internal temperature are both defined for the balloon.  The fixed interior temperature shall override the temperature difference.")
    }
    
    R <- 8.3145 #Universal gas constant
    g <- 9.81 #Gravitational acceleration

    if(is.null(T.interior)) {
        rv <- BalloonSphere(Mp, pa, po, T.ambient + T.diff)
    } else {
        rv <- BalloonSphere(Mp, pa, po, T.interior)
    }

    if(max.radius < rv$r) {
        rv$r <- max.radius
        rv$V <- (4/3)* pi * max.radius ^ 3
    }

    disp.mass=((pa*rv$V)/(R*T.ambient))*ma #Mass of displaced fluid

    if(sealed) {
        prop.mass=Mp*mp #Mass of lift gas
    } else {
        if(is.null(T.interior)) {
            prop.mass <- ((pa*rv$V)/(R*(T.ambient + T.diff)))*mp
        } else {
            prop.mass <- ((pa*rv$V)/(R*(T.interior)))*mp
        }
    }

    L <- (disp.mass-prop.mass-mi)*g #Net lift force

    A <- pi*rv$r^2 #Cross sectional area of balloon

    dm <- disp.mass/rv$V #Air density
    va <- sqrt((2*L)/(acd*dm*A))
    print(c(Mp * mp, prop.mass, disp.mass, L, acd, dm, A, va)) 
    invisible(va) #Ascent velocity
}

BalloonSphere <- function(Mp, pa, po, T) {#Get balloon dimensions
   #Calculate radius and volume of balloon if it is elastic
   #INPUTS
   #   MP is mols of lift gas
   #   PA is ambient air pressure (pascals)
   #   po is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
   #   T is internal temperature (kelvin)
   #OUTPUTS
   #   A list with elementds
   #   r - Balloon radius (m)
   #   V - Balloon volume (m3)

   R <- 8.3145 #Universal gas constant

   V <- (Mp*R*T)/(pa*po) #Balloon volume (m3)

    numerator <- 3.0*Mp*R*T

    denominator <- 4*pi*(pa*po)



    r <- ((3 * V)/(4 * pi)) ^(1/3)
    print(r)
    invisible(list(r = r, V = V))
}

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
    model.ind <- 1
    model.domain <- c(model.center.point[1] - model.span/2, 
        model.center.point[1] + model.span/2,
        model.center.point[2] + model.span/2,
        model.center.point[2] - model.span/2)
    #Read NOMADS server 
    urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 2)
    model.parameters <- ParseModelPage(urls.out[1])
    if(length(model.parameters$pred) == 0) { #If the GFS model has no data yet, go to the previous one
        model.parameters <- ParseModelPage(urls.out[2])
        model.ind <- 2
    }
  
    #Get model run date, convert to POSIX date 
    run.date <- str_match_all(urls.out[model.ind], "\\d{10}")[[1]][1,1]
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
   back.file <- GribGrab(urls.out[model.ind], back.pred, levels, variables,
       file.name = "fcst_back.grb", model.domain = model.domain)

   #Get forward forecast
   fore.file <- GribGrab(urls.out[model.ind], fore.pred, levels, variables,
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
    back.profile <- grd.int[[1]]$z[, , y.point.ind, x.point.ind]
    fore.profile <- grd.int[[2]]$z[, , y.point.ind, x.point.ind]

    #Make weighted profile by date
    atmos.profile <- (back.profile * weight.avg[1] + fore.profile * weight.avg[2])/sum(weight.avg)

    atmos.profile <- atmos.profile[ order(atmos.profile[,1]), ]
    invisible(atmos.profile)
}

DescentVelocity <- function(mi, ma, pa, rp, T, dcd) {
    #Calculate velocity of a (non spherical) descending object
    #INPUTS
    #    MI is mass of instrument package, parachute, etc (kg)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    PA is ambient air pressure (pascals)
    #    RP is parachute radius (assume circular) (m)
    #    T is temperature (Kelvin)
    #    DCD is coefficient of drag (unitless)

    #CODE VERIFIED BY NASA ROCKET PARACHUTE DESCENT MODEL

    R <- 8.3145 #Universal gas constant
    g <- 9.81 #Gravitational acceleration

    dm <- pa*ma/(R*T) #Air density
    A <- pi*rp*rp #Cross sectional area of parachute

    vd <- sqrt(2*mi*g/(dcd*dm*A)) #Descent velocity
    invisible(vd)
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

    xygrd <- GLOB.XY(fcst.grid$lat, fcst.grid$lon, proj)

    x.cells <- ceiling(((max(xygrd$y) - min(xygrd$y)) * 1000) / resolution)
    y.cells <- ceiling(((max(xygrd$x) - min(xygrd$x)) * 1000) / resolution)

    x.vec <- vector()
    y.vec <- vector()

    for(k in seq_len(length(xygrd$x))) {
        x.vec <- append(x.vec, rep(xygrd$x[k], length(fcst.grid$variables)))
    }

    y.vec <- append(y.vec, rep(xygrd$y, length(fcst.grid$variables)))
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

MolsFromLift <- function(mb, L, po, ma, mp) {
    #Calculate how many mols of lift gas you added to your balloon to get the measured lift at the ground
    #This assumes ISOTHERMAL CONDITIONS...so the balloon is the same temperature as the ambient air
    #INPUTS
    #    MB is the mass of the balloon envelope (kg)
    #    L is the amount of lift measured in the field (kg)
    #    PO is the overpressure (how much higher the pressure is inside the balloon, as a dimensionless ratio balloon pressure/ambient air pressure)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    MP is the molar mass of the propellant
    #OUTPUTS
    #MP  how many mols of lift gas you have added to get measured lift

    total.mass <- L + mb

    Mp <- total.mass/((ma/po)-mp)

    invisible(Mp) #How many mols of lift gas you added to your balloon to get the measured lift
}

MolsFromVolume <- function(V, pa, T, T.diff = 0, T.interior = NULL) {
   #How many mols of gas are inside the balloon, assuming the balloon is open to the atmosphere
   #For example, how many mols of air are inside my fully inflated solar balloon

   if(T.diff != 0 & !is.null(T.interior)) {
         warning("A temperature difference from ambient and a fixed internal temperature are both defined for the balloon.  The fixed interior temperature shall override the temperature difference.")
    }

   R <- 8.3145 #Universal gas constant
 
   if(is.null(T.interior)) {
       Mp <- pa * V / (R * (T + T.diff))
   } else {
       Mp <- pa * V / (R * T.interior)
   }

   invisible(Mp)
}

#Let's fly something! 

#Define model parameters 
model.span <- 2 #degrees to span
model.res <- 2000 #Resolution of interpolated grid
model.tol <- 0.49 #Distance from center point that triggers rebuilding of model (degrees) 
profile.tol <- c(600, 2000) #Rebuild atmospheric profile every X seconds or Y meters of drift from  point
variables <- c("TMP", "HGT", "UGRD", "VGRD")
levels <- paste(c(1, 2, 3, 5, 7, 10, 20, 30, 50, 70, seq(100, 1000, by = 25)), "mb")

#Define launch date and location
model.date <- as.POSIXlt(Sys.time() + 5, tz = "GMT") #Get data for this date
object.coords <- c(-106.912295, 34.064383, 1702) #Initial coordinates of point of interest
time.limit <- 3600 * 24 #How many seconds to fly

#Get this party started

#MODEL PARAMETERS
t <- 0
deltat <- 60 #Time step (seconds)
reload.model <- TRUE
tdiff.back <- Inf
tdiff.fore <- Inf

#BALLOON PARAMETERS
mb<-0.8 #Mass of envelope (kg)
L<- 2.2 #Net lift before instrument package is added (kg)
po<-1.0 #Overpressure (unitless)
ma<-0.02897 #Molar mass of dry air
mp<-0.004002602 #Molar mass of propellant gas
mi<- 2 #Mass of payload
dcd<-1.5 #Parachute coefficient of drag
rp<-0.762 #Parachute radius
max.radius<-3.5 #Balloon burst radius
T.diff <- 0 #Difference in temperature between balloon and surroundings (K)
sealed <- TRUE #Balloon does not vent to atmosphere

##### FIRST LET'S TEST TO MAKE SURE WE GET REASONABLE BALLOON AND PARACHUTE VALUES #####

Mp <- MolsFromLift(mb, L, po, ma, mp)
pa <- 100000
T <- 300
acd <- 0.5
va <- 5
for (k in seq_len(25)) {
        rv=BalloonSphere(Mp, pa, 1.1, T)
        va=AscentVelocity(mi, Mp, mp, ma, pa, po, acd, max.radius, T)
        Re=AscentRe(ma, pa, rv$V, T, va, rv$r)
        acd=AscentCd(Re)
}

balloon <- list(lat = c(), lon = c(), elev = c(), time = c())
if(1 == 1) {
while(t < time.limit) { #Time limit

    #If the particle is at the edge of the model, or we are exiting the time domain
    if(reload.model | (tdiff.back + tdiff.fore) > 3) { 
        print("Downloading model...")
        grd.int <- BuildAtmosphere(model.date, object.coords, model.span, model.res, variables, levels)
        profile.tol.tmp <- c(Inf, Inf)  #Force profile rebuild
        cart.pos <- c(0, 0, object.coords[3])
        reload.model <- FALSE
    }
    
    #Rebuild the atmospheric profiles if necessary
    if(profile.tol[1] < profile.tol.tmp[1] | profile.tol[2] < profile.tol.tmp[2]) {
        print("Reprofiling...")
        profile.tol.tmp <- c(0, 0) #Reset profile
        tdiff.back <-  abs(as.numeric(difftime(as.POSIXlt(grd.int[[1]]$fcst.date, tz = "GMT"), model.date, units = "hours")))
        tdiff.fore <- abs(as.numeric(difftime(as.POSIXlt(grd.int[[2]]$fcst.date, tz = "GMT"), model.date, units = "hours")))
        atmos.profile <- BuildProfile(object.coords, grd.int, c(tdiff.back, tdiff.fore))
        #Interpolate
        i.p   <- splinefun(atmos.profile[,1], 100 * rev(as.numeric(unlist(str_match_all(levels, "\\d+")))), method = "natural")
        i.tmp <- splinefun(atmos.profile[,1], atmos.profile[,2], method = "natural")
        i.wu  <- splinefun(atmos.profile[,1], atmos.profile[,3], method = "natural")
        i.wv  <- splinefun(atmos.profile[,1], atmos.profile[,4], method = "natural")
        profile.cart.pos <- cart.pos #So we can determine how far the balloon has gone to trigger reprofiling if necessary
    }
   
   #Distance from center of profile
   cart.pos[1] <- cart.pos[1] + i.wu(cart.pos[3]) * deltat #East - west movement
   cart.pos[2] <- cart.pos[2] + i.wv(cart.pos[3]) * deltat #North - south movement 

   #Ascent rate calcs
   T.ambient  <- i.tmp(cart.pos[3])
   pa <- i.p(cart.pos[3])
   rv <- BalloonSphere(Mp, pa, po, T.ambient) #Calculate balloon size
   if (rv$r > max.radius) { #Balloon pops if it goes too high
       break
   }
   Re <- AscentRe(ma, pa, rv$V, T.ambient, va, rv$r) #Calculate Reynold's number
   acd <- AscentCd(formula = "smoothsphere") #Calculate ascent coefficient of drag
   va <- AscentVelocity(mi, Mp, mp, ma, pa, po, acd, max.radius, T.ambient, T.diff = T.diff) #Calculate ascent rate
   cart.pos[3] <- cart.pos[3] + va * deltat #Elevation gain or loss

   #Latitude and longitude of object
   obj.tmp <- XY.GLOB((cart.pos[1] + i.wu(cart.pos[3]) * deltat)/1000, (cart.pos[2] + i.wv(cart.pos[3]) * deltat) / 1000, grd.int[[1]]$projection)
   object.coords[1] <- obj.tmp$lon
   object.coords[2] <- obj.tmp$lat
   object.coords[3] <- cart.pos[3]

   t <- t + deltat
   model.date <- model.date + deltat
   profile.tol.tmp <- c(profile.tol.tmp[1] + deltat, sqrt((profile.cart.pos[1] - cart.pos[1])^2 + (profile.cart.pos[2] - cart.pos[2])^2))
   #Correct for west longitude
   if(object.coords[1] > 180) {
       west.lon <- object.coords[1] - 360
   } else {
       west.lon <- object.coords[1] 
   }

   #Make sure we are still in the model domain
   if(abs(west.lon - mean(grd.int[[2]]$lon)) > model.tol | abs(object.coords[2] - mean(grd.int[[2]]$lat)) > model.tol) {
       reload.model <- TRUE
   }

   #Append data to balloon list
    balloon$lat <- append(balloon$lat, object.coords[2])
    balloon$lon <- append(balloon$lon, object.coords[1])
    balloon$elev <- append(balloon$elev, object.coords[3])
    balloon$time <- append(balloon$time, t)
    coords <- cbind(balloon$lon - 360, balloon$lat, balloon$elev)
    PolyLines2GE(coords, goo = "test_trajectory.kml")
}
}
