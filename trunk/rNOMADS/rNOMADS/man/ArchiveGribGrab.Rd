\name{ArchiveGribGrab}
\alias{ArchiveGribGrab}
\title{
Download archived model data from the NOMADS server.
}
\description{
This function interfaces with the programming API at \url{http://nomads.ncdc.noaa.gov/} to download archived NOMADS model data.
The available models can be viewed by calling \code{\link{NOMADSArchiveList}} without arguments.
The data arrives in grib (gridded binary) format that can be read with \code{\link{ReadGrib}}.
Some of these files are in GRIB format, others are in GRIB2 format; select the appropriate file type when calling \code{\link{ReadGrib}}.
}
\usage{
GribGrab(abbrev, model.date, model.run, pred,
    local.dir = ".", file.name = "fcst.grb",
    model.domain = NULL, tidy = FALSE, verbose = TRUE,
    check.url = TRUE, download.method = NULL)
}
\arguments{
  \item{abbrev}{
      Model abbreviation per \code{\link{NOMADSArchiveList}}.
}
  \item{model.date}{
      The year, month, and day of the model run, in YYYYMMDD format
}
   \item{model.run}{
       Which hour the model was run (i.e. 00, 06, 12, 18 for GFS)
}
  \item{pred}{
      Which prediction to get (analysis is 00)
  }
  \item{local.dir}{
  Where to save the grib file, defaults to the current directory.
}
  \item{file.name}{
  What to name the grib file, defaults to \code{"fcst.grb"}.
}
  \item{tidy}{
  If \code{TRUE}, remove all files with the suffix ".grb" from \code{local.dir} prior to downloading a new grib file.
}
  \item{verbose}{
  If \code{TRUE}, give information on connection status.  Default \code{TRUE}
}
 \item{download.method}{
  Allows the user to set the download method used by \code{download.file}: \code{"internal"}, \code{"wget"} \code{"curl"}, \code{"lynx"}.
  If \code{NULL} (the default), let R decide.
}
}
\value{
\item{grib.info$local.dir}{The absolute path to the grib file that was downloaded.}
\item{grib.info$file.name}{The name of the grib file that was downloaded.}
\item{grib.info$url}{The URL that the grib file was downloaded from}
}
\references{
\url{http://nomads.ncdc.noaa.gov/}
}
\author{
Daniel Bowman \email{daniel.bowman@unc.edu}
}

\seealso{
\code{\link{NOMADSArchiveList}}, \code{\link{ReadGrib}}
}
\examples{
Need to edit these still.
#An example for the Global Forecast System 0.5 degree model

#Get the latest model url
\dontrun{urls.out <- CrawlModels(abbrev = "gfs0.5", depth = 1)}

#Get a list of forecasts, variables and levels
\dontrun{model.parameters <- ParseModelPage(urls.out[1])}

#Figure out which one is the 6 hour forecast
#provided by the latest model run
#(will be the forecast from 6-12 hours from the current date) 

\dontrun{my.pred <- model.parameters$pred[grep("06$", model.parameters$pred)]}

#What region of the atmosphere to get data for
levels <- c("2 m above ground", "800 mb") 

#What data to return
variables <- c("TMP", "RH") #Temperature and relative humidity

#Get the data
\dontrun{grib.info <- GribGrab(urls.out[1], my.pred, levels, variables)}

#Extract the data
\dontrun{model.data <- ReadGrib(grib.info$file.name, levels, variables)}

#Reformat it
\dontrun{model.grid <- ModelGrid(model.data)}

#Show an image of world temperature at ground level
\dontrun{image(model.grid$z[2, 1,,])}

}
\keyword{ connection }
