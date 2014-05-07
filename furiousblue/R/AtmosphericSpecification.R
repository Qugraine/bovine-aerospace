#Code to get atmospheric data 

RNomads2FBProfile <- function(rnomads.profile, profile.index = 1) {
    #This function takes output from the BuildProfile function in rNOMADS
    #and converts it into a profile structure for use in furiousblue.
    #INPUTS
    #    RNOMADS.PROFILE - Output of the BuildProfile function in rNOMADS
    #    PROFILE.INDEX - Which profile to use in RNOMADS.PROFILE (as BuildProfile can return a list of profiles).
    #OUTPUTS
    #    PROFILE.ARRAY - A 5 x n array of atmospheric layers

    values <- rnomads.profile$profile.data[[1]]
    variables <- rnomads.profile$variables
    pressures <- as.numeric(stringr::str_replace(levels, " mb", "")) * 100

    profile.array <- array(rep(0, 5 * length(pressures)), dim = c(length(pressures), 5))
    colnames(profile.array) <- c("Elevation", "Pressure", "Temperature", "Zonal Wind", "Meridional Wind")
    profile.array[,1] <- values[, variables == "HGT"]
    profile.array[,2] <- rev(pressures)
    profile.array[,3] <- values[, variables == "TMP"]
    profile.array[,4] <- values[, variables == "UGRD"]
    profile.array[,5] <- values[, variables == "VGRD"]

    return(profile.array)
}

AtmosphericProfile2Spline  <- function(profile.array, method = "natural") {
   #Convert an atmospheric profile to spline functions with respect to elevation.
   #INPUTS
   #    PROFILE.ARRAY - Atmospheric profile in form specified by furiousblue (see RNomads2FBProfile)
   #    METHOD - Method for spline function, see R help for splinefun
   #OUTPUTS 
   #   PROFILE.SF - List of spline functions for pressure, temperature, zonal wind, and meridional wind, all with respect to elevation

   profile.sf <- list(Prs = splinefun(profile.array[, 1], profile.array[, 2], method = method),
       Tmp = splinefun(profile.array[, 1], profile.array[, 3], method = method),
       W.z = splinefun(profile.array[, 1], profile.array[, 4], method = method),
       W.m = splinefun(profile.array[, 1], profile.array[, 5], method = method))

   return(profile.sf)
}
