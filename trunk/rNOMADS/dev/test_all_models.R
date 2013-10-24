source("../rNOMADS/R/ReadGrib.R")
source("../rNOMADS/R/GetGrib.R")
library(stringr)
library(scrapeR)
library(fields)

model.list <- NOMADSList()
cat("bad models\n", file = "badmodels.txt")
for(m in model.list$abbrevs) {
   urls.out <- CrawlModels(abbrev = m, depth = 2, verbose = FALSE)
   print(paste("Model", m))
   model.parameters <- ParseModelPage(urls.out[length(urls.out)])
   if(length(model.parameters$pred) < 1  | length(model.parameters$variables) <1 | length(model.parameters$levels) <1) {
      cat(paste(m, "\t", urls.out[length(urls.out)], "\n", sep = ""), file = "badmodels.txt", append=TRUE)
   } else {
       file.name <- GribGrab(urls.out[length(urls.out)], model.parameters$pred[1], model.parameters$levels[1], model.parameters$variables[1], file.name = paste0(m,".grb"))
       ind <-0 
       while (file.info(file.name)$size == 0) {
           if(ind > length(model.parameters$variables)) {
               stop(paste("No data returned for model", m))
           }
           ind <- ind + 1
           for(k in seq_len(length(model.parameters$levels))) {
                file.name <- GribGrab(urls.out[length(urls.out)], model.parameters$pred[1], model.parameters$levels[k], model.parameters$variables[ind], file.name = paste0(m,".grb"))
               }
            if(file.info(file.name)$size > 0) {
               break
           }
       }

   }
}



