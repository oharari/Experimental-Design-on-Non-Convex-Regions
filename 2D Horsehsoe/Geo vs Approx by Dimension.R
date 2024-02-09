wd <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(wd)
source("pre_optimization_functions.r")

GeoD <- as.matrix(read.table("Horseshoe_Geodesic_Distances_0.05.txt"))
mdsG <- mds(GeoD)
V <- mdsG$embedding

K <- dim(V)[2]

error.k <- function(k)
{
	tempD <- V[,1:k]
	d <- as.matrix(dist(tempD))
	sqrt(sum((GeoD-d)^2))
}



e <- sapply(1:15,error.k)
m <- which.min(e)

pdf("Horseshoe_Geo_Dist_vs_Dimension.pdf", width = 9, height = 6)
par(mar = c(5,5,2,2))
plot(e, xlab = expression(k), ylab = expression(RMSE[k]), cex.lab = 1.6, type = 'l', lwd = 3)
lines(c(m,m), c(e[m], e[1]), lty = 3, lwd = 2)
dev.off()



