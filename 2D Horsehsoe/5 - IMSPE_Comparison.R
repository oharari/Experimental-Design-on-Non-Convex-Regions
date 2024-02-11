#*******************************************************************************
#* This code generates IMSPE-Optimal designs for the Annular Ring Example              
#*******************************************************************************



#*******************************************************************************
#* Before proceeding, make sure you have the following packages installed on 
#*  machine:
#*******************************************************************************
#install.packages("digest")
library(digest) #includes the "digest" function
library(proxy)
#*******************************************************************************


#*******************************************************************************
#* Insert the working directory, where all the files and tables are stored: 
#*******************************************************************************
wd = dirname(rstudioapi::getActiveDocumentContext()$path)
wd.maximin = paste(wd, "/Simulation/Maximin_Designs", sep='')
wd.imspe = paste(wd, "/Simulation/IMSPE_Designs", sep='')
wd.euclid = paste(wd, "/Simulation/Euclidean_Maximin_Designs", sep='')
wd.aux = paste(wd, "/Simulation/Auxiliary_Files", sep='')
wd.fw = paste(wd, "/Simulation/Fast_Ward_Designs", sep='')
#*******************************************************************************


#*******************************************************************************
#* Importing the functions file for the IMSPE-optimal designs:                  
#*******************************************************************************
setwd(wd)
source("IMSPE_Optimal_Design_Functions.R")
source("Pre_Optimization_Functions.R")
#*******************************************************************************




#*******************************************************************************
#* Script starts here!!! 		
#*******************************************************************************



#*******************************************************************************
#* list of necessary arguments             		                 
#*******************************************************************************
grid.density = 0.05 #your desired distance between grid points
embedding.dim = 4 #your desired dimension of the embedding space
N.Iter = 1 #your desired number of iterations for the IMSPE optimal design search
sigmatau2 = .01
epsilon = 0.05
#*******************************************************************************



#*******************************************************************************
#* importing pre-prepared tables            		            	  
#*******************************************************************************

#**importing the original (Cartesian) grid of candidates
Orig_Cands = read.table(paste(wd, "/Horseshoe_Candidates_", grid.density, 
                              ".txt", sep = ''))

#** importing the embedded grid of candidates
Embed_Cands = read.table(paste(wd.aux, "/Horseshoe_Embedded_Candidates_", 
                               embedding.dim, "Dim_", grid.density, ".txt", 
                               sep = ''))

#** importing the embedded grid of points for numerical integration                                 
Num_Integ_Pts = read.table(paste(wd.aux, "/Horseshoe_Embedded_Integration_", 
                                 embedding.dim, "Dim_", grid.density, ".txt", 
                                 sep = ''))

#** reading the simulation design  
Design_List = as.matrix(read.table(paste(wd, 
                                         "/Simulation/Simulation_Design_List.txt", 
                                         sep='')))
#*******************************************************************************



#** calculating the IMSPEs of the maximin geodesic designs and their relative 
#** efficiencies with respect to the "true" IMSPE-optimal design
Eff.Maximin = rep(0, nrow(Design_List))

for(i in 1:nrow(Design_List)){
  des = Design_List[i,]
  size = des[1]
  rho = des[2:5]
  design = read.table(paste(wd.maximin,"/Horseshoe_Size_", size,
                            "_Embedded_Maximin_Design.txt", 
                            sep=''))	
  
  opt.design = read.table(paste(wd.imspe, "/Horseshoe_Size_", size, 
                                "_Embedded_IMSPE_Design_", (i-1)%/%3+1, 
                                ".txt", sep=''))
  
  Eff.Maximin[i] = tempmse(rho, 
                           opt.design, 
                           Num_Integ_Pts, 
                           sigmatau2)/tempmse(rho, 
                                              design, 
                                              Num_Integ_Pts, 
                                              sigmatau2)
}



#** calculating the IMSPEs of the maximin Euclidean designs and their relative 
#** efficiencies with respect to the "true" IMSPE-optimal design
Eff.Euclid.Maximin = rep(0, nrow(Design_List))
for(i in 1:nrow(Design_List)){
  des = Design_List[i,]
  size = des[1]
  rho = des[2:5]
  design = read.table(paste(wd.euclid,"/Horseshoe_Size_", size ,
                            "_Embedded_Euclidean_Maximin_Design.txt", sep=''))	
  
  opt.design = read.table(paste(wd.imspe, "/Horseshoe_Size_", size, 
                                "_Embedded_IMSPE_Design_", (i-1)%/%3+1, 
                                ".txt", sep=''))
  
  Eff.Euclid.Maximin[i] = tempmse(rho,
                                  opt.design,Num_Integ_Pts,
                                  sigmatau2)/tempmse(rho,
                                                     design,
                                                     Num_Integ_Pts,
                                                     sigmatau2)
}


Eff.FW = rep(0, nrow(Design_List))
for(i in 1:nrow(Design_List)){
  des = Design_List[i,]
  size = des[1]
  rho = des[2:5]
  design = as.matrix(read.table(paste(wd.fw,"/Fast_Ward_Size_", size, 
                                      ".txt", sep='')))	

  design = read.table(paste(wd.fw, "/Embedded_Fast_Ward_Size_", size, 
                            ".txt", sep=''))
  
  opt.design = read.table(paste(wd.imspe, "/Horseshoe_Size_", size, 
                                "_Embedded_IMSPE_Design_", (i-1)%/%3+1, 
                                ".txt", sep=''))
  
  Eff.FW[i] = tempmse(rho,
                      opt.design,
                      Num_Integ_Pts,
                      sigmatau2)/tempmse(rho,
                                         design,
                                         Num_Integ_Pts,
                                         sigmatau2)
}

Maximin.Eff = data.frame(cbind(Design_List, 
                               Geo_Maximin_Eff = round(Eff.Maximin, 3), 
                               Euc_Maximin_Eff = round(Eff.Euclid.Maximin, 3), 
                               Fast_Ward_Eff = round(Eff.FW, 3)))
Maximin.Eff

Ave.Maximin.Eff = aggregate(cbind(Eff.Maximin, Eff.Euclid.Maximin, Eff.FW), 
                            by = list(Maximin.Eff$n), FUN = "mean")
Ave.Maximin.Eff


rho = rbind(c(.9, .9, .9, .9), 
            c(.6, .6, .6, .6), 
            c(.3, .3, .3, .3), 
            c(.9, .3, .9, .9))

sizes = c(8, 12, 16)
IMSPE.Iopt = c()
Eff = c()

for(n in sizes){
  for(i in 1:nrow(rho)){
    opt.design = read.table(paste(wd.imspe, "/Horseshoe_Size_", n, 
                                  "_Embedded_IMSPE_Design_", i, ".txt",
                                  sep=''))
    
    opt.IMSPE = tempmse(rho[i,], opt.design, Num_Integ_Pts,sigmatau2)
    
    for(j in 1:nrow(rho)){
      design = read.table(paste(wd.imspe,"/Horseshoe_Size_", n, 
                                "_Embedded_IMSPE_Design_", j, ".txt", 
                                sep=''))
      
      temp = tempmse(rho[i,], design, Num_Integ_Pts, sigmatau2)
      IMSPE.Iopt = c(IMSPE.Iopt, temp)
      Eff = c(Eff, opt.IMSPE/temp)
    }
  }	
}



one.col = function(n,k,l){
  u = rep(0, n)
  u[k:(l + k - 1)] = 1
  u
}

n = (nrow(rho))^2
l = nrow(rho)
temp = sapply((l*(0:(l-1))+1), one.col, n=n, l=l)
true.rho = temp%*%rho
temp = matrix(rep(diag(l),l), ncol=l, byrow=T)
assumed.rho = temp%*%rho
a = cbind(true.rho, assumed.rho)

temp = matrix(rep(diag(nrow(a)),length(sizes)), ncol=nrow(a), byrow=T)
all = temp%*%a
IMSPE.all = data.frame(cbind(rep(sizes, each = n), all, 
                             round(IMSPE.Iopt, 3), round(Eff, 3)))
names(IMSPE.all) = c("n", paste("true_rho", 1:l, sep=''), 
                     paste("assumed_rho", 1:l, sep=''), "IMSPE", "Eff")
IMSPE.all

IMSPE.Eff = aggregate(IMSPE.all$Eff, 
                      by = list(IMSPE.all$n, 
                                IMSPE.all$assumed_rho1, 
                                IMSPE.all$assumed_rho2, 
                                IMSPE.all$assumed_rho3, 
                                IMSPE.all$assumed_rho4), 
                      FUN = "mean")

IMSPE.Eff$x = round(IMSPE.Eff$x,3)

names(IMSPE.Eff) = names(IMSPE.Eff) = c("n", "rho1", "rho2", "rho3", "rho4", 
                                        "Efficiency")
IMSPE.Eff
