#***********************************************************************************************
#***** Before proceeding, make sure the following packages are installed on your machine: ******
#***********************************************************************************************
library(fields) #includes the function transformx()
#***********************************************************************************************
#***********************************************************************************************




wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/2D Horsehsoe"
setwd(wd)
wd.read <- paste(wd, "/Simulation/Auxiliary_Files", sep='')
wd.write <- paste(wd, "/Simulation/Euclidean_Maximin_Designs", sep='')

source("Pre_Optimization_Functions.R")


epsilon <- .05 # coarseness parameter of the grid
sizes <- c(8,12,16) # desired size of your designs

#** reading the list of candidates within the horseshoe region
cands <-  read.table(paste("Horseshoe_Candidates_", epsilon,".txt", sep = ''))

#** calculating the Euclidean distance matrix
Dist <- as.matrix(dist(cands))

#** reading the embedded (4D) candidates
cands.embedded <- read.table(paste(wd.read, "/Horseshoe_Embedded_Candidates_4Dim_",epsilon,".txt",sep =''))



#** generating optimal designs according to the maximin Euclidean distance criterion
#** plotting them and writing them (and their 4D embeddings) into files
for(i in 1:length(sizes))
{
	design.size <- sizes[i]

	D <- MaximinDesign(cands, design.size, Dist, nruns = 10, max.loop = 1000)$design


	tit <- substitute(paste('Size ', nn, ' Euclidean Maximin Design', sep = ''), list(nn=design.size))
	pdf(paste(wd.write, "/Horseshoe_Size_", design.size, "_Euclidean_Maximin_Design.pdf", sep = ''), 	width = 6.2, height = 6)
	par(mar=c(5,5,3,3))
	plot.cands=as.matrix(expand.grid(seq(0,1,.05),seq(0,1,.05)))
	hs=horseshoe(plot.cands)$newcands
	plot(plot.cands,xlab=expression(X[1]), ylab=expression(X[2]), col="white", xlim = c(0,1),
  	   ylim = c(0,1), cex.lab = 1.6)
	title(main =  tit, cex.main = 1.5)

	xt=seq(0.4,0,-.001)
	yt=sqrt(.4^2-(0.4-xt)^2)
	xt=c(xt,rev(xt))
	yt=c(yt+.5,rev(0.5-yt))
	xt=c(xt,1,1,.4)
	yt=c(yt,0.1,.45,.45)
	xt2=seq(0.4,0.35,-.001)
	yt2=sqrt(.05^2-(0.4-xt2)^2)
	xt2=c(xt2,rev(xt2))
	yt2=c(0.5-yt2,rev(0.5+yt2))
	xt=c(xt,xt2,.4,1,1,.4)
	yt=c(yt,yt2,.55,.55,.9,.9)
	polygon(xt,yt,lwd=3)
	points(hs,col="grey")
	points(D, pch=19, cex = 1.5, col = 4)
	dev.off()


	write.table(D, paste(wd.write, "/Horseshoe_Size_", design.size, "_Euclidean_Maximin_Design.txt", sep = ''))
	write.table(cands.embedded[rownames(D),], paste(wd.write, "/Horseshoe_Size_", design.size, 	"_Embedded_Euclidean_Maximin_Design.txt", sep = ''))
}

