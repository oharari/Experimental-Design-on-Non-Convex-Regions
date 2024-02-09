#***************************************************************************************************
#***************************************************************************************************
#********** This code creates some useful tables for later use in the Annulr Ring Example **********
#***************************************************************************************************
#***************************************************************************************************



#***************************************************************************************************
#***************************************************************************************************
#***** Before proceeding, make sure you have the following packages installed on your machine: *****
#***************************************************************************************************
#***************************************************************************************************
library(fOptions) #contains the runif.sobol function
library(fields) #includes the function transformx()
library(rgl)
library(proxy)
#***************************************************************************************************
#***************************************************************************************************





wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
wd.write <- paste(wd, "/Simulation/Auxiliary_Files", sep='')


source("Pre_Optimization_Functions.R")




epsilon <- 0.05 # epsilon distance between grid points
numapproxevals <- 4 # Embedding dimension

#** writing the list of candidates within the horseshoe region into file
dilution(epsilon)


#** reading the list of candidates previously written into file
candidates <- read.table(paste("Horseshoe_Candidates_" , as.character(epsilon), ".txt", sep = ''))
#plot(candidates)

#** performing Isomap
Iso <- Isomap(candidates, epsilon)

#** extracting the (n-dimensional) embeddings and the geodesic distances and writing them to files
ma <- Iso$EmbedPts
GeoDist <- Iso$GeoDist
write.table(ma, paste(wd.write, "/Horseshoe_Embedded_Candidates_Full_", as.character(epsilon), ".txt", sep = ''))
write.table(GeoDist, paste(wd.write, "/Horseshoe_Geodesic_Distances_", as.character(epsilon), ".txt", sep = ''))


#Focusing on the fisrt 4 dimensions of the embedding (and writing to file)
ma <- Iso$EmbedPts[,1:numapproxevals]
write.table(ma, paste(wd.write, "/Horseshoe_Embedded_Candidates_", as.character(numapproxevals), "Dim_", as.character(epsilon), ".txt", sep = ''))



#** plotting the neighborhood graph made of the horseshoe candidates (unnecessary)
G <- makeGdelta(GeoDist,sqrt(2)*epsilon)
plot3d(ma, size = 8, col = 4, xlab = '', ylab = '', zlab='')
for(i in 1:nrow(G))
{
	for(j in 1:nrow(G))
	{
		if(G[i,j]) lines3d(c(ma[i,1],ma[j,1]),c(ma[i,2],ma[j,2]),c(ma[i,3],ma[j,3]))
	}
}
rgl.snapshot("Horseshoe_Embedding_3D.png")




#** choosing locations for numerical integration (to be further diluted later)
pts=as.matrix(runif.sobol(4096,2))  #gives us points to use in computing the approx integral
pts=horseshoe(pts, epsilon) 
pts=pts$newcands

#** embedding the integration locations in the 4D region
temp <- Embed_Integ_Pts(pts, candidates, Iso, epsilon) 
temp <- temp$Embed_Integ_Pts

#** choosing a subset of size 1000 according to the maximin distance criterion
#** to make sure that the numerical integration is done with respect to a uniform
#** measure over the embedded region
D_Integ <- as.matrix(dist(temp))
Emb_Integ <- MaximinDesign(temp, 1000, D_Integ, nruns = 5, max.loop = 200)

#** writing the resultant set of location for numerical integration into file
write.table(Emb_Integ$design, paste(wd.write, "/Horseshoe_Embedded_Integration_", as.character(numapproxevals), "Dim_", as.character(epsilon), ".txt", sep = ''))






