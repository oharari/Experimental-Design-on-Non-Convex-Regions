#***********************************************************************************************
#***** Before proceeding, make sure the following packages are installed on your machine: ******
#***********************************************************************************************
library(fields) #includes the function transformx()
library(rgl)
#***********************************************************************************************
#***********************************************************************************************



wd <- getwd()
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
setwd(wd)
source("Pre_Optimization_Functions.R")


N.Integration.pts <- 1000 # Number of integration points

#** reading the embedded locations of the previous step
data <- read.table(paste(wd.aux, "/Embedded_Integration_Pts_Glacier_Scaled_all.txt", sep=''))

#** calculating the (approximate) geodesic distance matrix
GeoD <- as.matrix(dist(data))

#** choosing a subset of the desired size according to the maximin geodesic distance criterion
#** to allow for integration with respect to a uniform measure over the embedded region
D <- MaximinDesign(data, N.Integration.pts, GeoD, nruns = 5, max.loop = 12000)$design

#writing the resultant sample into file
write.table(D, paste(wd.aux, "/Glacier_Uniform_Integration_Pts.txt", sep=''))

