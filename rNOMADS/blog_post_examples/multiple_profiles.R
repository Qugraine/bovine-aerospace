library(rNOMADS)

#Download list of USArray station locations
download.file("http://www.unc.edu/~haksaeng/rNOMADS/myTA.RDATA", destfile = "myTA.RDATA")
load("myTA.RDATA")


variables <- c("TMP", "UGRD", "VGRD", "HGT")
