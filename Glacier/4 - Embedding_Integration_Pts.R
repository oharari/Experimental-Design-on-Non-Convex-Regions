library(proxy)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
setwd(wd)


source("Pre_Optimization_Functions.R")

delta <- 2e-5 #size of the neighborhood about each location
k <- 4 #embedding dimension


#** reading the (original and scaled) locations on the glacier
origpts <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt", sep='')))
cands <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary.txt", sep='')))
low <- min(cands)
high <- max(cands)


#** reading the geodesic distance matrix
geoD <- read.table(paste(wd.aux, "/Glacier_Geodesic_Distances.txt", sep=''))

#** performing MDS 
mdsG <- mds(geoD)

#** reading another 5922 locations on the glacier 
newpts <- as.matrix(read.table(paste(wd.aux, "/Glacier_Integration_Points.txt", sep=''), header = T))

#** rescaling the locations previously read so that they're contained in the unit cube 
newpts.scaled <- (newpts-low)/(high-low)


#** retaining the desired number of eigenvectors and eigenvalues
evals <- mdsG$evals[1:k]
evecs <- mdsG$evecs[,1:k]

#** embedding the new locations in 4D using Nystrom's method
Embed_Integ_Pts <- mdsemapping(geoD,origpts,newpts.scaled,delta,evals,evecs)

#** Writing the embedded locations into file
write.table(Embed_Integ_Pts, paste(wd.aux, "/Embedded_Integration_Pts_Glacier_Scaled_all.txt", sep=''))



