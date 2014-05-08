#Use the GrADS-DODS capability of NOMADS to get ascii data

GetDODSDates <- function(abbrev) {
    #Checks the real time GrADS data server to see what dates and model subsets are available for model specified by ABBREV.
    #INPUTS
    #    ABBREV - Model abbreviation for real time NOMADS data
    #OUTPUTS
    #    AVAILABLE.DATES - A list of model URLS and dates
    #        $ABBREV - Model abbreviation
    #        $DATE - Model run date, YYYYMMDD
    #        $URL - Model url

    date.pattern <- "[1-2]\\d{3}[0-1]\\d{1}[0-3]\\d{1}$"

    top.url <- NOMADSRealTimeList("dods", abbrev)$url

   if(!RCurl::url.exists(top.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", top.url, " but did not succeed..."))
   }

    html.tmp <- XML::htmlParse(top.url)
    top.links <- XML::xpathSApply(html.tmp, "//a/@href")
    XML::free(html.tmp)
 
    #List entries in html.tmp, see if they are dates or not
    date.links.lind <- grepl(date.pattern, top.links)
    if(sum(date.links.lind) == 0) { #If there do not appear to be dates here, go down one more directory level
        urls.tmp <- c()
        for(k in seq(4, length(top.links) - 3)) {
               if(!RCurl::url.exists(top.links[k])) {
                   stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
                       "Details:  Attempted to access ", top.links[k], " but did not succeed..."))
                }
            html.tmp.low <- XML::htmlParse(top.links[k])
            links.low <- XML::xpathSApply(html.tmp.low, "//a/@href")
            XML::free(html.tmp.low)
            urls.tmp <- append(urls.tmp, links.low[grepl(date.pattern, links.low)])
        }
    } else {
        urls.tmp <- top.links[date.links.lind]
    }

    urls.tmp <- as.character(urls.tmp)

    return(list(model = abbrev, date = stringr::str_extract(urls.tmp, date.pattern), url = urls.tmp)) 
}

GetDODSModelRuns <- function(model.url) {
   #Given a URL of a certain model date, determine which model runs are available for that date.
   #This is useful to check to see if a certain model run is on the website before attempting to get information or data from it.
   #The model url probably comes from GetDODSDates.
   #INPUTS 
   #    MODEL.URL - A URL pointing to the DODS model page for a certain date; get this URL from GetDODSDates.
   #OUTPUTS
   #    AVAILABLE.MODEL.RUNS - Which model runs are available for that date
   #        $MODEL.RUN - The model run
   #        $MODEL.RUN.INFO - Info about the model run, hence the name

      if(!RCurl::url.exists(model.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", model.url, " but did not succeed..."))
   }

   html.tmp <- XML::htmlParse(model.url)
   model.runs <- XML::xpathSApply(html.tmp, '//b', xmlValue) 
   XML::free(html.tmp)
   html.txt <- readLines(model.url)

   model.run.str <- NULL
   model.info.str <- NULL
   for(k in seq_len(length(model.runs))) {
       model.run.tmp <-  paste0(stringr::str_extract_all(stringr::str_replace(model.runs[k], "^\\d+:", ""), "[a-zA-Z0-9_.]")[[1]], collapse = "")
       model.info.tmp <- stringr::str_replace(html.txt[grepl(paste0(model.run.tmp, ":"), html.txt)], "</b>&nbsp;", " ")
       model.run.str <- append(model.run.str, model.run.tmp)
       model.info.str <- append(model.info.str, model.info.tmp)
   } 

   return(list(model.run = model.run.str, model.run.info = model.info.str))
}

GetDODSModelRunInfo <- function(model.url, model.run) {
   #Get description of the model run,  documentation (if present), longitude and latitude covered, time span, levels (if present), and variable list
   #INPUTS
   #    MODEL.URL is a URL pointing to a certain model date, probably from GetDODSDates.
   #    MODEL.RUN is a specified model run, probably from GetDODSModelRuns.
   #OUTPUTS
   #    MODEL.PARAMETERS - List of variables and model coverage information

   info.url <- paste0(model.url, "/", model.run, ".info")

   if(!RCurl::url.exists(info.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", info.url, " but did not succeed..."))
   }
            
   info.table <- readHTMLTable(info.url)[[2]]
   info.arr <- cbind(
       as.vector(info.table[,1]),
       as.vector(info.table[,2]),
       as.vector(info.table[,3]))

   info.arr[which(is.na(info.arr), arr.ind=TRUE)] <- ""
   model.parameters <- apply(stringr::str_replace_all(info.arr, "Â", ""), 1, paste, collapse = " ")

   return(model.parameters)
}

DODSGrab <- function(model.url, model.run, variable, time, lon, lat, levels = NULL) {
   #Get data from DODS.  Note that this is slower than GribGrab but will work on all operating systems.
   #Also, we can only do one variable at a time.
   #The output of this function will be the same as the output of ReadGrib in order to maintain consistency across rNOMADS.
   #ALL INDICES START FROM ZERO 
   #INPUTS
   #    MODEL.URL is a URL pointing to a certain model date, probably from GetDODSDates. 
   #    MODEL.RUN is a specified model run, probably from GetDODSModelRuns.
   #    VARIABLE is a variable from the list returned by GetDODSModelRunInfo
   #    TIME is an **index list** of times per info from GetDODSModelRunInfo, "c(x,y)"
   #    LON is an **index list** of longitudes per info from GetDODSModelRunInfo "c(x,y)"
   #    LAT is an **index list** of latitudes per info from GetDODSModelRunInfo "c(x,y)"
   #    LEVELS is an **index list** of levels per info from GetDODSModelRunInfo "c(x,y)"
   #         if not NULL, try to request the variable at a certain level.  Will fail if the variable does not have associated levels.
   #OUTPUTS
   #    MODEL.DATA - the model as an array, with columns for the model run date (when the model was run)
   #       the forecast (when the model was for), the variable (what kind of data), the level (where in the atmosphere or the Earth, vertically)
   #       the longitude, the latitude, and the value of the variable.

   preamble <- paste0(model.url, "/", model.run, ".ascii?", variable)
   time.str <- paste0("[", paste0(time, collapse = ":"), "]")

   l.ind <- !is.null(levels)

   if(l.ind) {
       level.str <- paste0("[", paste0(levels, collapse = ":"), "]")
   } else {
       level.str <- ""
   }
   lat.str <- paste0("[", paste0(lat, collapse = ":"), "]")
   lon.str <- paste0("[", paste0(lon, collapse = ":"), "]")
  
   data.url <- paste0(preamble, time.str, level.str, lat.str, lon.str)  

   #RCurl needs to be loaded for this to work I think
   data.txt <- readLines(data.url)
   
   lons <- as.numeric(unlist(strsplit(data.txt[grep("^lon,", data.txt) + 1], split = ",")))
   lats <- as.numeric(unlist(strsplit(data.txt[grep("^lat,", data.txt) + 1], split = ",")))
   if(l.ind) {
       levels <- as.numeric(unlist(strsplit(data.txt[grep("^lev,", data.txt) + 1], split = ",")))
   }
   t.ind <- grep("^time,", data.txt)
   times <- as.numeric(unlist(strsplit(data.txt[t.ind + 1], split = ",")))

   #Extract data values 

   val.txt <- data.txt[2:(t.ind - 1)]    
   val.txt <- val.txt[val.txt !=""] 
   val.txt <- stringr::str_replace_all(val.txt, "\\]\\[", ",")
   val.txt <- stringr::str_replace_all(val.txt, c("\\]|\\["), "")

   model.run.date <- paste0(stringr::str_extract(model.url, "[1-2]\\d{3}[0-1]\\d{1}[0-3]\\d{1}$"), model.run)
   
   row.num <- (stringr::str_count(val.txt[1], ",") - 3 + l.ind) * length(val.txt)
   model.data.tmp <- array(rep("", row.num), dim = c(row.num, 7))

   c <- 1 #Counter
   for(k in seq_len(length(val.txt))) {
       val.tmp <- sapply(strsplit(val.txt[k], split = ","), as.numeric)
       for(j in seq_len(length(val.tmp) - 2 - l.ind)) {
           
           model.data.tmp[c, 1] <- model.run.date
           model.data.tmp[c, 2] <- times[val.tmp[1] + 1]
           model.data.tmp[c, 3] <- variable
           if(l.ind) {
               model.data.tmp[c, 4] <- levels[val.tmp[2] + 1]
           }
           model.data.tmp[c, 5] <- lats[val.tmp[2 + l.ind] + 1]
           model.data.tmp[c, 6] <- lons[j]
           model.data.tmp[c, 7] <- val.tmp[2 + j + l.ind]
           c <- c + 1
       }     
   }

   model.data <- list(
       model.run.date = model.data.tmp[,1],
       forecast.date  = model.data.tmp[,2],
       variables      = model.data.tmp[,3],
       levels         = model.data.tmp[,4],
       lon            = as.numeric(model.data.tmp[,5]),
       lat            = as.numeric(model.data.tmp[,6]),
       value          = model.data.tmp[,7])
 
   return(model.data)
}

