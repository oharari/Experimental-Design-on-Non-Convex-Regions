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




#***************************************************************************************************
#***************************************************************************************************
#** This function receives a grid of points in the [0,1]X[0,1] square and retains only the        **
#** that are inside the horseshoe                                                              **
#**                                                                                               **
#** called by: dilution                                                                           **
#** arguments: candidates - a dense grid in the [0,1]X[0,1] square                                **
#** output: a table consisting of grid points inside the annular ring                             **
#***************************************************************************************************
#***************************************************************************************************
horseshoe3d<-function(candidates)
{
  n=nrow(candidates)
  idx=as.vector(rep(TRUE,n))

  # first, rotate to the reference plane
  cands=candidates
  for(i in 1:n)
  {
    r=sqrt( (candidates[i,1])^2 + (candidates[i,3])^2 )
    # right quadrant postive, left quadrant negative
    # add .Machine$double.eps is hack for sign(0) to give us 1 not 0.
    cands[i,1]=r*sign(candidates[i,1]+.Machine$double.eps) 
    cands[i,3]=0
  }
  
  for(i in 1:n)
  {
    if(cands[i,1]>1 || cands[i,1]< -1)
      idx[i]=FALSE
    else if(cands[i,1]<=.5 && cands[i,1]>=0)
    {
      di=sqrt((cands[i,1]-.5)^2+(cands[i,2]-.5)^2)
      if(di>=.4-.Machine$double.eps || di<=.05+.Machine$double.eps)
        idx[i]=FALSE
    }
    else if(cands[i,1]>= -.5 && cands[i,1]<=0)
    {
      di=sqrt((cands[i,1]+.5)^2+(cands[i,2]-.5)^2)
      if(di>=.4-.Machine$double.eps || di<=.05+.Machine$double.eps)
        idx[i]=FALSE
    }
    else
    {
      di=abs(cands[i,2]-.5)
      if(di<=.05-.Machine$double.eps || di>=.4+.Machine$double.eps)
        idx[i]=FALSE
    }
  }
  newcands=candidates[idx,]
  newcands[,1]=newcands[,1]/2+0.5 #move it to [0,1]

  return(list(newcands=newcands,idx=idx))
}
#***************************************************************************************************
#***************************************************************************************************



#***********************************************************************************************
#***********************************************************************************************
#** Modified version of the cover.design function from the FIELDS package.                    **
#** we add the param dmat to contain the matrix of pairwise distances                         **
#** and change the code to use this matrix instead of calling a DIST function.                **
#***********************************************************************************************
#***********************************************************************************************
MaximinDesign <- function(R, nd, dmat, nruns = 1, nn = TRUE, num.nn = 100, fixed = NULL, 
    scale.type = "unscaled", R.center, R.scale, P = -20, Q = 20, 
    start = NULL, DIST = NULL, return.grid = TRUE, return.transform = TRUE, 
    max.loop = 20, verbose = FALSE) 
{
    if (!is.null(start) && is.matrix(start)) {
        if (any(duplicated.array(start))) 
            stop("Error: start must not have duplicate rows")
        start <- rowmatch(start, R)
        if (any(is.na(start))) 
            stop("Error: Starting design must be a subset of R")
    }
    R.orig <- R
    R <- as.matrix(R)
    if (nd >= nrow(R)) {
        stop(" number of design points >= the number of candidates")
    }
    if (any(duplicated.array(R))) 
        stop("Error: R must not have duplicate rows")
    if (num.nn >= (nrow(R) - nd)) {
        nn <- FALSE
        warning("Number of nearst neighbors (nn) reduced to the actual number of candidates")
    }
    if (is.null(DIST)) 
        DIST <- function(x, y) {
            rdist(x, y)
        }
    id <- 1:nrow(R)
    if (!is.null(start)) 
        nd <- length(start)
    if (is.null(fixed)) 
        n <- nd
    else {
        n <- nd + length(fixed)
    }
    R <- transformx(R, scale.type, R.center, R.scale)
    transform <- attributes(R)
    saved.crit <- rep(NA, nruns)
    saved.designs <- matrix(NA, nrow = nruns, ncol = n)
    saved.hist <- list(1:nruns)
    if (verbose) {
        cat(dim(R), fill = TRUE)
    }
    for (RUNS in 1:nruns) {
        if (is.null(start)) {
            if (!is.null(fixed)) {
                Dset <- sample((1:nrow(R))[-fixed], nd)
                Dset <- c(Dset, fixed)
            }
            else Dset <- sample(1:nrow(R), nd)
        }
        else {
            if (length(start) > nd) 
                stop("Error: the start matrix must have nd rows")
            Dset <- start
            if (!is.null(fixed)) 
                Dset <- c(Dset, fixed)
        }
        design.original <- R.orig[Dset, ]
        Dset.orginal <- Dset
        Cset <- id[-Dset]
        #dist.mat <- DIST(R[Cset, ], R[Dset, ])
        dist.mat <- dmat[Cset,Dset]
        rs <- dist.mat^P %*% rep(1, n)
        crit.i <- crit.original <- sum(rs^(Q/P))^(1/Q)
        CRIT <- rep(NA, length(Cset))
        CRIT.temp <- rep(NA, length(Cset))
        hist <- matrix(c(0, 0, crit.i), ncol = 3, nrow = 1)
        loop.counter <- 1
        repeat {
            for (i in 1:nd) {
                Dset.i <- matrix(R[Dset[i], ], nrow = 1)
                if (verbose) {
                  cat("design point", i, Dset.i, fill = T)
                }
                #partial.newrow <- sum(DIST(Dset.i, R[Dset[-i], 
                #  ])^P)
                partial.newrow <- sum(dmat[Dset[i], Dset[-i] 
                  ]^P)
                #rs.without.i <- rs - c(DIST(Dset.i, R[-Dset, 
                #  ])^P)
                rs.without.i <- rs - c(dmat[Dset[i], -Dset 
                  ]^P)
                if (nn) 
                  vec <- (1:length(Cset))[order(dist.mat[, i])[1:num.nn]]
                else vec <- 1:length(Cset)
                for (j in vec) {
                  Cset.j <- matrix(R[Cset[j], ], nrow = 1)
                  #newcol <- c(DIST(Cset.j, R[c(-Dset, -Cset[j]), 
                  #  ])^P)
                  newcol <- c(dmat[Cset[j], c(-Dset, -Cset[j]) 
                    ]^P)
                  CRIT[j] <- (sum((rs.without.i[-j] + newcol)^(Q/P)) + 
                  #  (DIST(Cset.j, Dset.i)^P + partial.newrow)^(Q/P))^(1/Q)
                    (dmat[Cset[j], Dset[i]]^P + partial.newrow)^(Q/P))^(1/Q)
                  if (verbose) {
                    cat(j, " ")
                  }
                }
                best <- min(CRIT[!is.na(CRIT)])
                best.spot <- Cset[CRIT == best][!is.na(Cset[CRIT == 
                  best])][1]
                if (verbose) {
                  cat(i, "best found ", best, " at", best.spot, 
                    fill = T)
                }
                crit.old <- crit.i
                if (best < crit.i) {
                  if (verbose) {
                    cat(i, "best swapped ", fill = T)
                  }
                  crit.i <- best
                  hist <- rbind(hist, c(Dset[i], best.spot, crit.i))
                  Dset[i] <- best.spot
                  Cset <- id[-Dset]
                  #dist.mat <- DIST(R[Cset, ], R[Dset, ])
                  dist.mat <- dmat[Cset, Dset]
                  rs <- (dist.mat^P) %*% rep(1, n)
                }
            }
            if ((crit.i == crit.old) | (loop.counter >= max.loop)) 
                break
            loop.counter <- loop.counter + 1
        }
        saved.crit[RUNS] <- crit.i
        saved.designs[RUNS, ] <- Dset
        saved.hist[[RUNS]] <- hist
    }
    ret <- (1:nruns)[saved.crit == min(saved.crit)]
    if (length(ret) > 1) {
        print("Greater than 1 optimal design; keeping first one......")
        ret <- ret[1]
    }
    crit.i <- saved.crit[ret]
    hist <- saved.hist[[ret]]
    nhist <- nrow(hist)
    nloop <- nruns
    hist <- cbind(c(0:(nrow(hist) - 1)), hist)
    dimnames(hist) <- list(NULL, c("step", "swap.out", "swap.in", 
        "new.crit"))
    out.des <- R[saved.designs[ret, ], ]
    out.des <- unscale(out.des, transform$x.center, transform$x.scale)
    out <- list(design = out.des, call = match.call(), best.id = c(saved.designs[ret, 
        ]), fixed = fixed, opt.crit = crit.i, start.design = design.original, 
        start.crit = crit.original, history = hist, other.designs = saved.designs, 
        other.crit = saved.crit, DIST = DIST, nn = nn, num.nn = num.nn, 
        P = P, Q = Q, nhist = nhist - 1, nloop = (nloop - 1)/n)
    if (return.grid) 
        out$grid <- R.orig
    if (return.transform) 
        out$transform <- transform
    class(out) <- "spatial.design"
    out
}
#***********************************************************************************************
#***********************************************************************************************





#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance parameter, creates a regular grid in the unit square with   **
#** said distance between points, retains the ones that are in the annular ring and writes them   **
#** into a text file                                                                              **         
#** calls: horseshoe                                                                              **
#** arguments: epsilon - a distance parameter                                                     **
#** output: none                                                                                  **
#***************************************************************************************************
#***************************************************************************************************
dilution <- function(epsilon)
{
	candidates <- as.matrix(expand.grid(seq(-1,1,by=epsilon),seq(0,1,by=epsilon),seq(0,1,by=epsilon)))
	candidates=horseshoe3d(candidates)$newcands
	write.table(candidates,paste("Horseshoe_3D_Candidates_", as.character(epsilon), ".txt", sep=''))
}
#***************************************************************************************************
#***************************************************************************************************




wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/3D Horsehsoe"
setwd(wd)
wd.write <- paste(wd, "/Simulation/Auxiliary_Files", sep='')


source("Pre_Optimization_Functions.R")




epsilon <- 0.1# epsilon distance between
pts1 <- -1+2*runif.sobol(16,1)
pts2 <- pts3 <- runif.sobol(16,1)
pts <- as.matrix(expand.grid(pts1,pts2,pts3))
pts=horseshoe3d(pts)
dilution(epsilon)
pts=pts$newcands
plot3d(pts)



candidates <- read.table(paste("Horseshoe_3D_Candidates_" , as.character(epsilon), ".txt", sep = ''))
plot3d(candidates)


Iso <- Isomap(candidates, epsilon)
ma <- Iso$EmbedPts
GeoDist <- Iso$GeoDist
write.table(ma, paste(wd.write, "/Horseshoe_3D_Embedded_Candidates_Full_", as.character(epsilon), ".txt", sep = ''))
write.table(GeoDist, paste(wd.write, "/Horseshoe_3D_Geodesic_Distances_", as.character(epsilon), ".txt", sep = ''))




numapproxevals <- 5

ma <- Iso$EmbedPts[,1:numapproxevals]
write.table(ma, paste(wd.write, "/Horseshoe_3D_Embedded_Candidates_", as.character(numapproxevals), "Dim_", as.character(epsilon), ".txt", sep = ''))




temp <- Embed_Integ_Pts(pts, candidates, Iso, epsilon) 
temp <- temp$Embed_Integ_Pts
D_Integ <- as.matrix(dist(temp))
Emb_Integ <- MaximinDesign(temp, 1000, D_Integ, nruns = 1, max.loop = 12000)
#
#
#
write.table(Emb_Integ$design, paste(wd.write, "/Horseshoe_Embedded_Integration_", as.character(numapproxevals), "Dim_", as.character(epsilon), ".txt", sep = ''))






