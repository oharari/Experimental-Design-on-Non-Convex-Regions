#***********************************************************************************************
#***** Before proceeding, make sure the following packages are installed on your machine: ******
#***********************************************************************************************
library(fields) #includes the function transformx()
library(rgl)
#***********************************************************************************************
#***********************************************************************************************


wd <- "C:/Users/Ofir/Dropbox/Research/SFU/Non_Convex_Regions/Glacier"
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
wd.maximin <- paste(wd, "/Glacier_Maximin_Designs", sep='')
setwd(wd)
source("Pre_Optimization_Functions.R")


#** The desired number of stakes 
design.size <- 9

#** reading the least of (scaled) locations on the glacier and the geodesic distance matrix
cands <- read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt", sep=''))
GeoDist <- read.table(paste(wd.aux, "/Glacier_Geodesic_Distances.txt", sep=''))


#** choosing a subset of the desired size according to the maximin geodesic distance criterion
D <- MaximinDesign(cands, design.size, GeoDist, nruns = 1, max.loop = 200)$design

#** reading the least of embedded candidates in their 4D representation
cands.embedded <- read.table(paste(wd.aux, "/Glacier_Embeddings_4D.txt", sep=''))


#** plotting the design in 2D 
tit <- substitute(paste('Size ', nn, ' Maximin Design', sep = ''), list(nn=design.size))

boundary2D <- rbind(cands[1:247,1:2], cands[1,1:2])
pdf(paste(wd.maximin, "/Glacier_Size_", design.size, "_Maximin_Design.pdf", sep=''))
par(mar = c(5,5,2,2))
plot(cands[,1:2], col = 'gray', xlab = "Latitude", ylab = 'Longitude', cex.lab = 1.5)
title(main = paste("Size ", design.size, " Maximin Design", sep=''), cex.main = 1.5)
for(i in 1:247)
{
	lines(c(boundary2D[i,1],boundary2D[i+1,1]),c(boundary2D[i,2],boundary2D[i+1,2]), lwd = 3)
}
points(D[,1:2], cex = 1.5, pch = 19, col=4)
dev.off()


#** Writing the resultant design into file
write.table(D, paste(wd.maximin, "/Glacier_Size_", design.size, "_Maximin_Design.txt", sep = ''))

#** Writing the embedded (4D) resultant design into file
write.table(cands.embedded[rownames(D),], paste(wd.maximin, "/Glacier_Size_", design.size, "_Embedded_Maximin_Design.txt", sep = ''))



