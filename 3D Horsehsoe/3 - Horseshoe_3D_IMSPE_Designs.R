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
wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/R Code and Datasets/3D Horsehsoe"
wd.maximin <- paste(wd, "/Simulation/Maximin_Designs", sep='')
wd.imspe <- paste(wd, "/Simulation/IMSPE_Designs", sep='')
wd.aux <- paste(wd, "/Simulation/Auxiliary_Files", sep='')
wd.fw <- paste(wd, "/Simulation/Fast_Ward_Designs", sep='')
#***************************************************************************************************


#***************************************************************************************************
#** Importing the functions file for the IMSPE-optimal designs:                                   **
#***************************************************************************************************
source(paste(wd, "/IMSPE_Optimal_Design_Functions.R", sep=''))
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



#***************************************************************************************************
#**		                          list of necessary arguments             		                  **		
#***************************************************************************************************
grid.density <- 0.1 #your desired distance between grid points
embedding.dim <- 5 #your deired dimension of the embeddibg space
N.Iter <- 1 #your desired number of iterations for the IMSPE optimal design search
sigmatau2 <- .01
#***************************************************************************************************



#***************************************************************************************************
#**		                          importing pre-prepared tables            		            	  **
#***************************************************************************************************
Orig_Cands <- read.table(paste(wd, "/Horseshoe_3D_Candidates_", grid.density, ".txt", sep = ''))
#importing the original (Cartesian) grid of candidates

Embed_Cands <- read.table(paste(wd.aux, "/Horseshoe_3D_Embedded_Candidates_", embedding.dim, "Dim_", 
                                grid.density, ".txt", sep = ''))
#importing the embedded grid of candidates
                                
Embed_Integ_Pts <- read.table(paste(wd.aux, "/Horseshoe_3D_Embedded_Integration_", embedding.dim, "Dim_", 
                                grid.density, ".txt", sep = ''))
#importing the embedded grid of points for numerical integration   

Design_List <- as.matrix(read.table(paste(wd, "/Simulation/Simulation_Design_List.txt", sep='')))
#***************************************************************************************************



#***************************************************************************************************
#**		                          generating and exporting the design                    		  **
#***************************************************************************************************

for(k in 1:nrow(Design_List))
{
	Des <- Design_List[k,]
	design.size <- Des[1] #your desired number of runs
	rho <- c(Des[2:6])
	#Init_Design <- read.table(paste(wd.maximin, "/Horseshoe_3D_Size_", design.size, "_Embedded_Maximin_Design.txt", 
                              #  sep = ''))
	Init_Design <- read.table(paste(wd.fw, "/Embedded_Fast_Ward_Size_", design.size, ".txt", sep = ''))
	#importing a maximin design of identical size, to be used as an initial guess

	Design <- imsedesign(Init_Design, Embed_Cands, Embed_Integ_Pts, rho, sigmatau2=sigmatau2)$design

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

#***************************************************************************************************
	write.table(Design, paste(wd.imspe, "/Horseshoe_3D_Size_", design.size, "_Embedded_IMSPE_Design_", 
	            ((k-1)%/%3)+1, ".txt", sep = ''))
	write.table(D, paste(wd.imspe, "/Horseshoe_3D_Size_", design.size, "_IMSPE_Design_", ((k-1)%/%3)+1,
	            ".txt", sep = ''))

}


