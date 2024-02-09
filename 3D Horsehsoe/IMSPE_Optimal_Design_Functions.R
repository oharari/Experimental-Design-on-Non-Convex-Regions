



#*****************************************************************************************
#*****************************************************************************************
#** This function calculates the correlation matrix at a design matrix based on the     **
#** sqaured-exponential correlation function of an isotropic Gaussian process           **
#**  																					**
#** is called by: tempmse; 																**
#** arguments: X - a design matrix; rho - a correlation parameter						**
#** values: a correlation matrix														**
#*****************************************************************************************
#*****************************************************************************************
corr.matrix <- function(X, rho)
{
	X <- as.matrix(X)
	n <- dim(X)[1]
	Theta <- diag(log(rho))	

	U <- matrix(apply(X^2%*%Theta, 1, sum), nrow = n, ncol = n ,byrow=F)
	V <- -2*X%*%Theta%*%t(X)
	Dist <- U + t(U) + V
	return(exp(4*Dist))
}
#*****************************************************************************************
#*****************************************************************************************

corr.matrix2 <- function(X, theta)
{
	X <- as.matrix(X)
	n <- dim(X)[1]
	Theta <- diag(theta)	

	U <- matrix(apply(X^2%*%Theta, 1, sum), nrow = n, ncol = n ,byrow=F)
	V <- -2*X%*%Theta%*%t(X)
	Dist <- U + t(U) + V
	return(exp(-Dist))
}



#*****************************************************************************************
#*****************************************************************************************
#** This function calculates the cross-correlation matrix between a design and a new    **
#** set of points																		**
#**  																					**
#** is called by: tempmse; 																**
#** arguments: X.new - new set of points; X.old - a design matrix; 						**
#** 				   rho - a correlation parameter									**
#** values: a cross-correlation matrix													**
#*****************************************************************************************
#*****************************************************************************************
cross.corr.matrix <- function(rho, X.new, X.old)
{
	X.new <- as.matrix(X.new)
	X.old <- as.matrix(X.old)
	n1 <- dim(X.new)[1]
	n2 <- dim(X.old)[1]
	d <- dim(X.old)[2]
	theta <- log(rho)
	Theta <- diag(theta)
	
	U <- matrix(apply(X.new^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=F)
	V <- -2*X.new%*%Theta%*%t(X.old)
	W <- matrix(apply(X.old^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=T)
	Dist <- U + V + W
	return(exp(4*Dist))
}
#*****************************************************************************************
#*****************************************************************************************


cross.corr.matrix2 <- function(theta, X.new, X.old)
{
	X.new <- as.matrix(X.new)
	X.old <- as.matrix(X.old)
	n1 <- dim(X.new)[1]
	n2 <- dim(X.old)[1]
	d <- dim(X.old)[2]

	Theta <- diag(theta)
	
	U <- matrix(apply(X.new^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=F)
	V <- -2*X.new%*%Theta%*%t(X.old)
	W <- matrix(apply(X.old^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=T)
	Dist <- U + V + W
	return(exp(-Dist))
}



#*****************************************************************************************
#*****************************************************************************************
#** This function calculates the Integrated Mean-Squared Prediction Error at a design   **
#**  																					**
#** calls: corr.matrix; cross.corr.matrix;                                              **
#** is called by: imsedesign; 															**
#** arguments: rho - a correlation parameter; design - a design matrix; 				**
#** 				 Embed_Integ_Pts - embedded points for numerical intergration;		**
#**                  sigmatau2 - a nugget parameter										**
#** values: The IMSPE																	**
#*****************************************************************************************
#*****************************************************************************************
tempmse <- function(rho,design,Embed_Integ_Pts,sigmatau2=1e-5)
{
	  N=nrow(design)
      one=rep(1,N)
      R=corr.matrix(design,rho)
			    # IMSE about 8% of the time for rho=.99, ss=10.  This would
			    # dissapear around 1e-7.  I set it to 1e-5 and it seems
			    # to work all the way to rho=.99999999, ss=10.  Hopefully
			    # will also be sufficient for larger designs (ss>10).
        Rinv=0
        ntow=length(sigmatau2)
   	    E=matrix(1,nrow=N+1,ncol=N+1)
        E[1,1]=0
        for(i in 1:ntow)
        {
          E[2:(N+1),2:(N+1)]=R+diag(N)*sigmatau2[i]
          Rinv=(Rinv*(i-1)+solve(E))/i
        }

        numpoints=nrow(Embed_Integ_Pts)

	    Rt <- cross.corr.matrix(rho, Embed_Integ_Pts, design)

        rt=rep(0,N)
        rtrt=matrix(0,nrow=N,ncol=N)
        for(i in 1:numpoints)
        {
   	      rt=rt+Rt[i,]
            rtrt=rtrt+Rt[i,]%*%t(Rt[i,])
        }
  	    rt=rt/numpoints
 	    rtrt=rtrt/numpoints
   	    RT=matrix(1,nrow=N+1,ncol=N+1)
        RT[2:(N+1),2:(N+1)]=rtrt
        RT[1,2:(N+1)]=t(rt)
        RT[2:(N+1),1]=rt
        imsehat= (1+sigmatau2-sum(diag(Rinv%*%RT)))*2/3
    
     	imsehat
}
#*****************************************************************************************
#*****************************************************************************************




tempmse2 <- function(theta,design,Embed_Integ_Pts,sigmatau2=1e-5)
{
	  N=nrow(design)
      one=rep(1,N)
      R=corr.matrix2(design,theta)
			    # IMSE about 8% of the time for rho=.99, ss=10.  This would
			    # dissapear around 1e-7.  I set it to 1e-5 and it seems
			    # to work all the way to rho=.99999999, ss=10.  Hopefully
			    # will also be sufficient for larger designs (ss>10).
        Rinv=0
        ntow=length(sigmatau2)
   	    E=matrix(1,nrow=N+1,ncol=N+1)
        E[1,1]=0
        for(i in 1:ntow)
        {
          E[2:(N+1),2:(N+1)]=R+diag(N)*sigmatau2[i]
          Rinv=(Rinv*(i-1)+solve(E))/i
        }

        numpoints=nrow(Embed_Integ_Pts)

	    Rt <- cross.corr.matrix2(theta, Embed_Integ_Pts, design)

        rt=rep(0,N)
        rtrt=matrix(0,nrow=N,ncol=N)
        for(i in 1:numpoints)
        {
   	      rt=rt+Rt[i,]
            rtrt=rtrt+Rt[i,]%*%t(Rt[i,])
        }
  	    rt=rt/numpoints
 	    rtrt=rtrt/numpoints
   	    RT=matrix(1,nrow=N+1,ncol=N+1)
        RT[2:(N+1),2:(N+1)]=rtrt
        RT[1,2:(N+1)]=t(rt)
        RT[2:(N+1),1]=rt
        imsehat= (1+sigmatau2-sum(diag(Rinv%*%RT)))*2/3
    
     	imsehat
}




#*****************************************************************************************
#*****************************************************************************************
#** This function finds an (approximate) IMSPE-optimal design - in the alternative      **
#** (embedded) coordinate system - using row exchanges with an hashing algorithm to     **
#** spare repeated comparisons of similar designs                                    	**
#** 							                            							**
#** calls: tempmse; digest (from the "digest" package);                                  **
#** is called by: IMSPE_Optim_Design; 					     							**
#** arguments: Init_Design - an initial "guess"; Embed_Cands - the embedded grid;       **
#** 				 Embed_Integ_Pts - embedded points for numerical intergration;		**
#**                  sigmatau2 - a nugget parameter;									**
#**                  hshindex; hashlist; imselist; maxhash - hashing parameters         **
#** values: imspe0 - the IMSPE of the final design; design - the chosen design;      	**
#**         hshindex; hashlist - lists of hash values; imselist - list of IMSPEs        **
#*****************************************************************************************
#*****************************************************************************************
imsedesign <- function(Init_Design, Embed_Cands, Embed_Integ_Pts, rho, sigmatau2=1e-5,
                       hshindex=NULL, hashlist=NULL, imselist=NULL, maxhash=16000)
{
	numcandidates=nrow(Embed_Cands)
	N=nrow(Init_Design)

	imse0=tempmse(rho,Init_Design,Embed_Integ_Pts,sigmatau2)
	oldimse=imse0+1 #dummy initializer
	imse=rep(1e50,N)
	cand=rep(0,N)
    imsehist=rep(0,N) #for gfx output only

    if(is.null(hshindex)) #default, do hash initialization
    {
          hashlist=rep(NA,maxhash)
          imselist=rep(NA,maxhash)# setup the hash with initial entry
          hashlist[1]=digest(Init_Design,algo="crc32")
          imselist[1]=imse0
          hshindex=2
    }

	for(i in 1:N)
	{
		for(j in 1:numcandidates)
		{
			design=Init_Design
			design[i,]=Embed_Cands[j,]
            curhsh=digest(design,algo="crc32")
            curhdx=which(hashlist==curhsh)
			curimse=imselist[curhdx]
            if(length(curhdx)==0)
            {
                 curimse=tempmse(rho,design,Embed_Integ_Pts,sigmatau2)
                 hashlist[hshindex]=curhsh
                 imselist[hshindex]=curimse
                 hshindex=hshindex%%maxhash+1
            }


			if(curimse<imse[i])
			{
				imse[i]=curimse
				cand[i]=j
			}

		}
		Init_Design[imse==min(imse),]=Embed_Cands[cand[imse==min(imse)],]
		oldimse=imse0

        imse0=imselist[which(hashlist==digest(Init_Design,algo="crc32"))]
	}

	return(list(imse=imse0,design=Init_Design,hshindex=hshindex,
	            hashlist=hashlist,imselist=imselist))

}
#*****************************************************************************************
#*****************************************************************************************




imsedesign2 <- function(Init_Design, Embed_Cands, Embed_Integ_Pts, theta, sigmatau2=1e-5,
                       hshindex=NULL, hashlist=NULL, imselist=NULL, maxhash=16000)
{
	numcandidates=nrow(Embed_Cands)
	N=nrow(Init_Design)

	imse0=tempmse2(theta,Init_Design,Embed_Integ_Pts,sigmatau2)
	oldimse=imse0+1 #dummy initializer
	imse=rep(1e50,N)
	cand=rep(0,N)
    imsehist=rep(0,N) #for gfx output only

    if(is.null(hshindex)) #default, do hash initialization
    {
          hashlist=rep(NA,maxhash)
          imselist=rep(NA,maxhash)# setup the hash with initial entry
          hashlist[1]=digest(Init_Design,algo="crc32")
          imselist[1]=imse0
          hshindex=2
    }

	for(i in 1:N)
	{
		for(j in 1:numcandidates)
		{
			design=Init_Design
			design[i,]=Embed_Cands[j,]
            curhsh=digest(design,algo="crc32")
            curhdx=which(hashlist==curhsh)
			curimse=imselist[curhdx]
            if(length(curhdx)==0)
            {
                 curimse=tempmse2(theta,design,Embed_Integ_Pts,sigmatau2)
                 hashlist[hshindex]=curhsh
                 imselist[hshindex]=curimse
                 hshindex=hshindex%%maxhash+1
            }


			if(curimse<imse[i])
			{
				imse[i]=curimse
				cand[i]=j
			}

		}
		Init_Design[imse==min(imse),]=Embed_Cands[cand[imse==min(imse)],]
		oldimse=imse0

        imse0=imselist[which(hashlist==digest(Init_Design,algo="crc32"))]
	}

	return(list(imse=imse0,design=Init_Design,hshindex=hshindex,
	            hashlist=hashlist,imselist=imselist))

}




#*****************************************************************************************
#*****************************************************************************************
#** This function calls "imsedesigns" potentially multiple times, and retains the best  **  
#** design                                   											**
#** 							                            							**
#** calls: imsedesign;                                                                  **
#** is called by: Non_Convex_IMSPE_Optim_Design;					   					**
#** arguments: Init_Design - an initial "guess"; Embed_Cands - the embedded grid;       **
#** 				 Embed_Integ_Pts - embedded points for numerical intergration;		**
#** 				 N.Iter - number of iterations;										**
#**                  sigmatau2 - a nugget parameter;									**
#** values: design - the chosen design; IMSPE - the IMSPE of the chosen design		 	**
#*****************************************************************************************
#*****************************************************************************************
IMSPE_Optim_Design <- function(Init_Design,Embed_Cands,rho,Embed_Integ_Pts,N.Iter,
                               sigmatau2=1e-5)
{
	imsevals <- rep(0, N.Iter)
	bestimspe <- 1e50
	
	for(i in 1:N.Iter) 
	{
 		d=imsedesign(Init_Design,Embed_Cands,Embed_Integ_Pts,rho,sigmatau2)
  
  		imsevals[i]=d$imse
 	    design=d$design
        if(d$imse<bestimspe) 
	    {
    	 	bestdesign=design
    		bestimspe=d$imse
  		}
	}
	return(list(design=bestdesign, IMSPE=bestimspe))
}
#*****************************************************************************************
#*****************************************************************************************



IMSPE_Optim_Design2 <- function(Init_Design,Embed_Cands,theta,Embed_Integ_Pts,N.Iter,
                               sigmatau2=1e-5)
{
	imsevals <- rep(0, N.Iter)
	bestimspe <- 1e50
	
	for(i in 1:N.Iter) 
	{
 		d=imsedesign2(Init_Design,Embed_Cands,Embed_Integ_Pts,theta,sigmatau2)
  
  		imsevals[i]=d$imse
 	    design=d$design
        if(d$imse<bestimspe) 
	    {
    	 	bestdesign=design
    		bestimspe=d$imse
  		}
	}
	return(list(design=bestdesign, IMSPE=bestimspe))
}



#*****************************************************************************************
#*****************************************************************************************
#** This function calls translates the design obtained by "IMSPE_Optim_Design" back to  **  
#** the original Cartesian coordinates													**
#** 							                            							**
#** calls: IMSPE_Optim_Design;                                                          **
#** arguments: Init_Design - an initial "guess"; Embed_Cands - the embedded grid;       **
#** 				 Embed_Integ_Pts - embedded points for numerical intergration;		**
#** 				 N.Iter - number of iterations;										**
#**                  sigmatau2 - a nugget parameter;									**
#**                  Orig_Cands - the original Cartesian grid;                          **
#** values: design - the chosen design; IMSPE - the IMSPE of the chosen design		 	**
#*****************************************************************************************
#*****************************************************************************************
Non_Convex_IMSPE_Optim_Design <- function(Init_Design, Embed_Cands, rho,Embed_Integ_Pts,
										  N.Iter, Orig_Cands, sigmatau2=1e-5)
{
	imspe.obj <- IMSPE_Optim_Design(Init_Design,Embed_Cands,rho,Embed_Integ_Pts,N.Iter)
	rows <- rownames(imspe.obj$design)
	Non.Conv.Design <- Orig_Cands[rows,]
	IMSPE <- imspe.obj$IMSPE
	
	return(list(design=Non.Conv.Design, IMSPE=IMSPE))
}
#*****************************************************************************************
#*****************************************************************************************




Non_Convex_IMSPE_Optim_Design2 <- function(Init_Design, Embed_Cands,theta,Embed_Integ_Pts,
										  N.Iter, Orig_Cands, sigmatau2=1e-5)
{
	imspe.obj <- IMSPE_Optim_Design2(Init_Design,Embed_Cands,theta,Embed_Integ_Pts,N.Iter)
	rows <- rownames(imspe.obj$design)
	Non.Conv.Design <- Orig_Cands[rows,]
	IMSPE <- imspe.obj$IMSPE
	
	return(list(design=Non.Conv.Design, IMSPE=IMSPE))
}



Non_Convex_IMSPE_Optim_Estimated <- function(X,y,Init_Design, Embed_Cands,Embed_Integ_Pts,
										  N.Iter, Orig_Cands)
{
	mod <- km(design=X, response=y, covtype="gauss", nugget.estim=TRUE, 
                control = list(trace=FALSE))
	theta <- 1/(2*(mod@covariance)@range.val^2)
	sigmatau2 <- (mod@covariance)@nugget
	imspe.obj <- IMSPE_Optim_Design2(Init_Design,Embed_Cands,theta,Embed_Integ_Pts,N.Iter,sigmatau2)
	IMSPE <- imspe.obj$IMSPE
	cat("Initial Design Criterion Value: ")
	print(tempmse2(theta,Init_Design,Embed_Integ_Pts,sigmatau2))
	cat("\n")
	cat("Final Design Criterion Value: ")
	print(IMSPE)
	cat("\n")
	
	
	return(list(design=imspe.obj$design, IMSPE=IMSPE))
}







