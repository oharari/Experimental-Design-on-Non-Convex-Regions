#***********************************************************************************************
#***** Before proceeding, make sure the following packages are installed on your machine: ******
#***********************************************************************************************
library(fields) #includes the function transformx()
#***********************************************************************************************
#***********************************************************************************************




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



wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
wd.read <- paste(wd, "/Simulation/Auxiliary_Files", sep='')
wd.write.geo <- paste(wd, "/Simulation/Maximin_Designs", sep = '')
wd.write.euc <- paste(wd, "/Simulation/Euclidean_Maximin_Designs", sep = '')



source("Pre_Optimization_Functions.R")


epsilon <- .1


cands <-  read.table(paste("Horseshoe_3D_Candidates_", epsilon,".txt", sep = ''))
GeoDist <- read.table(paste(wd.read, "/Horseshoe_3D_Geodesic_Distances_", epsilon,".txt", sep = ''))


design.size <- 24
D <- MaximinDesign(cands, design.size, GeoDist, nruns = 10, max.loop = 1250)$design
cands.embedded <- read.table(paste(wd.read, "/Horseshoe_3D_Embedded_Candidates_5Dim_",epsilon,".txt",sep =''))






write.table(D, paste(wd.write.geo, "/Horseshoe_3D_Size_", design.size, "_Maximin_Design.txt", sep = ''))
write.table(cands.embedded[rownames(D),], paste(wd.write.geo, "/Horseshoe_3D_Size_", design.size, "_Embedded_Maximin_Design.txt", sep = ''))


library(rgl)

plot3d(D, size = 10, col = 4, xlab='x1', ylab='x2', zlab='x3')
#rgl.snapshot("3D_Horseshoe_Maximin_Design.png")


Dist <- as.matrix(dist(cands))
D2 <- MaximinDesign(cands, design.size, Dist, nruns = 10, max.loop = 1250)$design
write.table(D2, paste(wd.write.euc, "/Horseshoe_3D_Size_", design.size, "_Euclidean_Maximin_Design.txt", sep = ''))
write.table(cands.embedded[rownames(D2),], paste(wd.write.euc, "/Horseshoe_3D_Size_", design.size, "_Embedded_Euclidean_Maximin_Design.txt", sep = ''))

