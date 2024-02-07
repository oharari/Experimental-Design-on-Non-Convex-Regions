#***************************************************************************************************
#***************************************************************************************************
#********** This code creates some useful tables for later use in the Glacier Example     **********
#***************************************************************************************************
#***************************************************************************************************



#***************************************************************************************************
#***************************************************************************************************
#***** Before proceeding, make sure you have the following packages installed on your machine: *****
#***************************************************************************************************
#***************************************************************************************************
library(fOptions) #includes the runif.sobol function
library(proxy) #includes the cross-distance dist() function
#***************************************************************************************************
#***************************************************************************************************




wd <- getwd()
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
setwd(wd)



source("Pre_Optimization_Functions.R")




#** reading the list of scaled locations on the glacier
candidates <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt", sep=''), header = T))

#** calculating the Euclidean distance matrix between locations
distance <- as.matrix(dist(candidates))

#** size of the neighborhood about each candidate (to determine which pairs are connected by edges)
epsilon <- 2e-5

#** Dimension of the embedding space (number of eigenvalues/eigenvectors to be retained)
numapproxevals <- 4

#** performing Isomap on the glacier (represented by the set of scaled locations) 
Iso <- Isomap(candidates, epsilon, k=numapproxevals)

#** writing the embedded, 4-dimensional locations in their new coordinate system into file
ma <- Iso$EmbedPts
write.table(ma, paste(wd.aux, "/Glacier_Embeddings_", as.character(numapproxevals), "D.txt", sep = ''))

#** writing the geodesic distance matrix into file
GeoDist <- Iso$GeoDist
write.table(GeoDist, paste(wd.aux, "/Glacier_Geodesic_Distances.txt", sep = ''))








