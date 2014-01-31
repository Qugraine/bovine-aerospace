ReadGrib <- function(file.name, levels, variables) {
    #This is a function to read forecast data from a Grib file
    #INPUTS
    #    FILE.NAME - Grib file name
    #    VARIABLES - data to extract
    #    LEVELS - which levels to extract data from
    #OUTPUTS
    #    MODEL.DATA - the grib model as an array, with columns for the model run date (when the model was run)
    #       the forecast (when the model was for), the variable (what kind of data), the level (where in the atmosphere or the Earth, vertically)
    #       the longitude, the latitude, and the value of the variable.

    #Get specified data from grib file

    match.str <- ' -match "('
    for(var in variables) {
        match.str <- paste(match.str, var, "|", sep = "")
    }

    match.str.lst <- strsplit(match.str, split = "")[[1]]
    match.str <- paste(match.str.lst[1:(length(match.str.lst) - 1)], collapse = "")

    if(length(levels) > 0 & !is.null(levels)) {
        match.str <- paste(match.str, "):(", sep = "")
        for(lvl in levels) {
            match.str <- paste(match.str, lvl, "|", sep = "")
        }
    } else {
        match.str <- paste0(match.str, ")")
   }

    match.str.lst <- strsplit(match.str, split = "")[[1]]
    match.str <- paste(match.str, '"', sep = "")
    match.str <- paste(match.str.lst[1:(length(match.str.lst) - 1)], collapse = "")
    match.str <- paste(match.str, ")\"", sep = "")

    wg2.str <- paste('wgrib2 ',     
        file.name, ' -inv my.inv -csv - -no_header', 
        match.str, sep = "")
    
    #Get the data from the grib file in CSV format
    csv.str <- system(wg2.str, intern = TRUE)

    #HERE IS THE EXTRACTION
    model.data.vector <- strsplit(paste(gsub("\"", "", csv.str), collapse = ","), split = ",")[[1]]
    chunk.inds <- seq(1, length(model.data.vector) - 6, by = 7)
    model.data <- list(model.run.date = model.data.vector[chunk.inds],
        forecast.date = model.data.vector[chunk.inds + 1],
        variables = model.data.vector[chunk.inds + 2],
        levels = model.data.vector[chunk.inds + 3],
        lon = as.numeric(model.data.vector[chunk.inds + 4]),
        lat = as.numeric(model.data.vector[chunk.inds + 5]),
        value = model.data.vector[chunk.inds + 6]
        )
    return(model.data)
}
