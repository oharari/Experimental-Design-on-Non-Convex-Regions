wd <- getwd()
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
setwd(wd)



#** reading the physical locations of the locations on the glacier
cands <- read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary.txt", sep=''))

low <- min(cands)
high <- max(cands)


#** rescaling the locations so that they are contained in the unit cube
cands.scaled <- (cands-low)/(high-low)
write.table(cands.scaled, paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt"))

