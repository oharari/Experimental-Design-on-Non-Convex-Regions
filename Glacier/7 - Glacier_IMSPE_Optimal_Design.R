library(proxy)
library(DiceKriging)
library(digest)


year <- 2009 #year of observations


wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
wd.aux <- paste(wd, "/Auxiliary_Files", sep='')
wd.maximin <- paste(wd, "/Glacier_Maximin_Designs", sep='')
wd.imspe <- paste(wd, "/Glacier_IMSPE_Optimal_Designs/", year, "_Designs", sep='')
wd.data <- paste(wd, "/Glacier_2009_Measurements", sep='')

setwd(wd)
source("IMSPE_Optimal_Design_Functions.R")

design.size <- 14 #your desired number of runs
N.Iter <- 1 #(unused parameter)
maxit <- 3 #Number of external loops for the optimization algorithm


#** reading the embedded glacier observations from the year of interest
data <- as.matrix(read.table(paste(wd.data, "/Embedded_", year, "_Data_Glacier_Scaled.txt", sep=''), header = T))
p <- dim(data)[2]-1
X <- data[, 1:p]
y <- data[,p+1]




#***************************************************************************************************
#**		                          importing pre-prepared tables            			 			  **
#***************************************************************************************************

#**importing a maximin design of identical size, to be used as an initial guess
Init_Design <- as.matrix(read.table(paste(wd.maximin, "/Glacier_Size_", design.size, "_Embedded_Maximin_Design.txt", sep = '')))


#** importing the scaled locations on the glacier
Orig_Cands <- as.matrix(read.table(paste(wd.aux, "/Glacier_Interior_and_Boundary_Scaled.txt", sep='')))

#importing the embedded set of (4D) candidates
Embed_Cands <- as.matrix(read.table(paste(wd.aux, "/Glacier_Embeddings_4D.txt", sep='')))

#** importing the embedded set of (4D) points for numerical integration                          
Embed_Integ_Pts <- as.matrix(read.table(paste(wd.aux, "/Glacier_Uniform_Integration_Pts.txt", sep='')))  
#***************************************************************************************************


#** generating an IMSPE-optimal design on the glacier, based on the MLEs of the GP parameters
#** as estimated from the data
Embed_Opt_Design <- Non_Convex_IMSPE_Optim_Estimated(X,y,Init_Design,Embed_Cands,Embed_Integ_Pts,
													 N.Iter,maxit=maxit)


#** writing the embedded (4D) design to file													 
write.table(Embed_Opt_Design$design, paste(wd.imspe, "/Glacier_Size_", design.size, "_Embedded_IMSPE_Design_", year, ".txt", sep = ''))


#** finding the corresponding original (3D) coordinates of the optimal design						 
D <- matrix(0, nrow = design.size, ncol = 3)
for(i in 1:design.size)
{
	for(j in 1:nrow(Orig_Cands))
	if(sum(Embed_Opt_Design$design[i,]-Embed_Cands[j,])==0) D[i,] <- Orig_Cands[j,]
}

IMSPE <- Embed_Opt_Design$IMSPE


#** plotting the resultant design in 2D
tit <- substitute(paste('IMSPE=', I, sep = ''), list(I=round(IMSPE,4)))

boundary2D <- rbind(Orig_Cands[1:247,1:2], Orig_Cands[1,1:2])
pdf(paste(wd.imspe, "/Glacier_Size_", design.size, "_IMSPE_Design_", year, ".pdf", sep=''))
par(mar = c(6,5,2,2))
plot(Orig_Cands[,1:2], col = 'gray', xlab = "Latitude", ylab = 'Longitude', cex.lab = 1.5)
title(main = paste("Size ", design.size, " IMSPE-Optimal Design", sep=''), cex.main = 1.5)
title(sub=tit, cex.sub=1.2, line=5)
for(i in 1:247)
{
	lines(c(boundary2D[i,1],boundary2D[i+1,1]),c(boundary2D[i,2],boundary2D[i+1,2]), lwd = 3)
}
points(D[,1:2], cex = 1.5, pch = 19, col=4)
dev.off()


#** writing the resultant design into file
write.table(D, paste(wd.imspe, "/Glacier_Size_", design.size, "_IMSPE_Design_", year, ".txt", sep = ''))
