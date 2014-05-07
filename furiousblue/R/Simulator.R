# Flight simulation code

StaticProfileFlightSimulation <- function(balloon.config, atmospheric.profile) {
    #Simulates the balloon specified in BALLOON.CONFIG with a user-supplied atmospheric profile. 
    #This allows test atmospheres to be inputted.
    #INPUTS
    #    BALLOON.CONFIG - Configuration data structure for balloon flight
    #    ATMOSPHERIC.PROFILE - A n x 5 array, where n is the number of atmospheric layers.
    #        Column 1 is geopotential height (m), column 2 is pressure (Pa), column 3 is temperature (K),
    #        column 4 is zonal wind (m/s) and column 5 is meridional wind (m/s).
    #OUTPUTS
    #   TRAJECTORY - Data on balloon flight, as array. 
}
