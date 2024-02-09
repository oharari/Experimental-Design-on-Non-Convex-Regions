#***************************************************************************************************
#***************************************************************************************************
#**            This code generates IMSPE-Optimal designs for the Annulr Ring Example              **
#***************************************************************************************************
#***************************************************************************************************



#***************************************************************************************************
#**    Before proceeding, make sure you have the following packages installed on your machine:    **
#***************************************************************************************************
#install.packages("digest")
library(digest) #includes the "digest" function
#***************************************************************************************************


#***************************************************************************************************
#** Insert the working directory, where all the files and tables are stored:                      **
#***************************************************************************************************
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
#***************************************************************************************************


#***************************************************************************************************
#** Importing the functions file for the IMSPE-optimal designs:                                   **
#***************************************************************************************************
source("IMSPE_Optimal_Design_Functions.R")
#***************************************************************************************************



horseshoe<-function(candidates)
{
  n=nrow(candidates)
  idx=as.vector(rep(TRUE,n))
  
  for(i in 1:n)
  {
    if(candidates[i,1]<=.4)
    {
      di=sqrt((candidates[i,1]-.4)^2+(candidates[i,2]-.5)^2)
      if(di>.4 || di<.05-.Machine$double.eps)
        idx[i]=FALSE
    }
    else
    {
      di=abs(candidates[i,2]-.5)
      if(di<.05-.Machine$double.eps || di>.4)
        idx[i]=FALSE
    }
  }
  newcands=candidates[idx,]

  return(list(newcands=newcands,idx=idx))
}




#***************************************************************************************************
#***************************************************************************************************
#**					                 Script starts here!!! 					  **
#***************************************************************************************************
#***************************************************************************************************
wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/2D Horsehsoe"
wd.maximin <- paste(wd, "/Simulation/Maximin_Designs", sep='')
wd.imspe <- paste(wd, "/Simulation/IMSPE_Designs", sep='')
wd.euclid <- paste(wd, "/Simulation/Euclidean_Maximin_Designs", sep='')
wd.aux <- paste(wd, "/Simulation/Auxiliary_Files", sep='')


#***************************************************************************************************
#**		                          list of necessary arguments             		                  **		
#***************************************************************************************************
design.size <- 16 #your desired number of runs
grid.density <- 0.05 #your desired distance between grid points
embedding.dim <- 4 #your deired dimension of the embeddibg space
rho <- c(.5,.02,.5,.5) #your desired correlation parameters
N.Iter <- 2 #your desired number of iterations for the IMSPE optimal design search
sigmatau2 <- .01
#***************************************************************************************************



#***************************************************************************************************
#**		                          importing pre-prepared tables            		            	  **
#***************************************************************************************************
Init_Design <- read.table(paste(wd.imspe,"/Horseshoe_Size_", design.size,"_Embedded_IMSPE_Design_2.txt", sep = ''))
#importing a maximin design of identical size, to be used as an initial guess

Orig_Cands <- read.table(paste(wd,"/Horseshoe_Candidates_", grid.density, ".txt", sep = ''))
#importing the original (Cartesian) grid of candidates

Embed_Cands <- read.table(paste(wd.aux, "/Horseshoe_Embedded_Candidates_", embedding.dim, "Dim_", 
                                grid.density, ".txt", sep = ''))
#importing the embedded grid of candidates
                                
Embed_Integ_Pts <- read.table(paste(wd.aux, "/Horseshoe_Embedded_Integration_", embedding.dim, "Dim_", 
                                grid.density, ".txt", sep = ''))
#importing the embedded grid of points for numerical integration   
#***************************************************************************************************



#***************************************************************************************************
#**		                          generating and exporting the design                    		  **
#***************************************************************************************************
#Design <- Non_Convex_IMSPE_Optim_Design(Init_Design, Embed_Cands, rho,Embed_Integ_Pts,
#										  N.Iter, Orig_Cands)$design
Design <- imsedesign(Init_Design, Embed_Cands, Embed_Integ_Pts, rho, sigmatau2)$design

D <- matrix(0, nrow = design.size, ncol = ncol(Orig_Cands))
Orig_Cands <- as.matrix(Orig_Cands)
Design <- as.matrix(Design)
Embed_Cands <- as.matrix(Embed_Cands)
for(i in 1:design.size)
{
	for(j in 1:nrow(Embed_Cands))
	{
		if(sum(Design[i,]-Embed_Cands[j,])==0) D[i,] <- Orig_Cands[j,]
	}
}

#***************************************************************************************************
#**		          ploting the resultant design along with the grid of candidates          		  **
#***************************************************************************************************
	tit <- substitute(paste('Size ', nn, ' IMSPE-Optimal Design', sep = ''), list(nn=design.size))
	subtit <- substitute(paste(rho, "= (",r1,",",r2,",", r3,",",r4, ")", sep=''), 
	                     list(r1=rho[1], r2=rho[2], r3=rho[3], r4=rho[4]))
	pdf(paste(wd.imspe, "/Horseshoe_Size_", design.size, "_IMSPE_Design_", 4, 
	         ".pdf", sep = ''), width = 6.2, height = 6)
	par(mar=c(5,5,3,3))
	plot.cands=as.matrix(expand.grid(seq(0,1,.05),seq(0,1,.05)))
	hs=horseshoe(plot.cands)$newcands
	plot(plot.cands,xlab=expression(X[1]), ylab=expression(X[2]), col="white", xlim = c(0,1),
 	    ylim = c(0,1), cex.lab = 1.6)
	title(main =  tit, cex.main = 1.5, sub=subtit)
	
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


#***************************************************************************************************

write.table(D, paste(wd.imspe,"/Horseshoe_Size_", design.size, "_IMSPE_Design_4.txt", sep = ''))
write.table(Design, paste(wd.imspe,"/Horseshoe_Size_", design.size, "_Embedded_IMSPE_Design_4.txt", sep = ''))

