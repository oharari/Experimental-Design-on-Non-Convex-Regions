library(proxy)

wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/Glacier"
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
wd.data <- paste(wd, "/Glacier_2009_Measurements", sep='')
setwd(wd)

source("Pre_Optimization_Functions.R")


delta <- 2e-5 #size of the neighborhood about each location
k <- 4  #embedding dimension
year <- 2009 #year of the data


#** reading the glacier data from the specified year
data <- read.table(paste(wd.data, "/Glacier_Data_", year, ".txt", sep=''))

#** reading the (original and scaled) locations on the glacier
origpts <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt", sep='')))
cands <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary.txt", sep='')))
low <- min(cands)
high <- max(cands)

#** reading the geodesic distance matrix
geoD <- read.table(paste(wd.aux, "/Glacier_Geodesic_Distances.txt", sep=''))

#** performing MDS
mdsG <- mds(geoD)

#** rescaling the glacier observations so that they're contained in the unit cube
newpts <- data[,1:3]
newpts.scaled <- (newpts-low)/(high-low)

#** retaining the desired number of eigenvalues and eigenvectors
evals <- mdsG$evals[1:k]
evecs <- mdsG$evecs[,1:k]

#** embedding the glacier observations in the 4D space
Embed_Data_Pts <- mdsemapping(geoD,origpts,newpts.scaled,delta,evals,evecs)
Embed_Data_Pts <- cbind(Embed_Data_Pts, data[,4])

#** writing the observations in the new coordinates (in 4D)
write.table(Embed_Data_Pts, paste(wd.data, "/Embedded_", year, "_Data_Glacier_Scaled.txt", sep=''))

library(rgl)
plot3d((evecs%*%diag(evals^.5))[,1:3])
points3d(Embed_Data_Pts[,1:3], col = 2, size=9)