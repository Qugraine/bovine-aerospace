library(XML)
url <- "http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl"

WebCrawler <- function(url) {
    doc <- htmlParse(url)
    links <- xpathSApply(doc, "//a/@href")
    free(doc)
    if(is.null(links)) {
        return(url)
    } else {
        ret <- vector("list", length = length(links))
        for(link in links) {
           print(link)
           ret[[link]] <- WebCrawler(link)
        }
        return(ret)
    }
}

foo <- WebCrawler(url)
model.urls <- unlist(foo, recursive = TRUE, use.names = FALSE)
