#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance matrix and returns an adjacency matrix of pairs of points   **
#** located "close" together                                                                      **
#**                                                                                               **
#** Called by: Isomap                                                                             **
#** arguments: dist - a distance matrix; delta - closeness parameter                              **
#** output: a logical matrix with '1' indicating an edge                                          **
#***************************************************************************************************
#***************************************************************************************************
makeGdelta<-function(distance,delta)
{
	N=nrow(distance)
	adj=matrix(0,nrow=N,ncol=N) #our adjacency matrix

	for(i in 1:N)
		  adj[i,]=as.numeric(distance[i,]<(delta+.Machine$double.eps))

	return(adj)
}
#***************************************************************************************************
#***************************************************************************************************



#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance matrix and returns an adjacency matrix of pairs of points   **
#** located "close" together                                                                      **
#**                                                                                               **
#** Called by: Isomap                                                                             **
#** arguments: dist - a distance matrix; k - number of neighbors;    							  **
#**            boundary - the adjacency matrix, filled for the boundary points                    ** 
#** output: a logical matrix with '1' indicating an edge                                          **
#***************************************************************************************************
#***************************************************************************************************
makeGknn<-function(distance,k, boundary)
{
	N=nrow(distance)
	s <- apply(boundary, 1, sum)
	bound.indices <- which(s>0)
	adj=boundary #our adjacency matrix

	for(i in (1:N)[-bound.indices])
	{
		  d <- distance[i,]
		  ord <- rank(d)
		  adj[i,] <- as.numeric(ord>1&ord<=(k+1))
	}
		  
	return(adj)
}
#***************************************************************************************************
#***************************************************************************************************




#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance matrix and an adjacency matrix, and returns a matrix of     **
#** approximated geodesic distances, using Floyd's algorithm                                      **
#**                                                                                               **
#** Called by: Isomap                                                                             **
#** arguments: dist - a distance matrix; adj - an adjecancy matrix                                **
#** output: distG - a geodesic distance matrix                                                    **
#***************************************************************************************************
#***************************************************************************************************
floyd<-function(distance,adj)
{
	print("Into Floyd\n")
	flush.console()
	N=nrow(distance)
	if(ncol(distance)!=N) stop("floyd():  Square matrix required!\n")

	distG=matrix(10000,nrow=N,ncol=N) #-1 => infinity

	for(i in 1:N)
	{
		for(j in 1:N)
		{
			if(adj[i,j])
			{
				distG[i,j]=distance[i,j]
			}
		}
	}
	print("Done Initializing\n")
	flush.console()                
	for(k in 1:N)
	{
		for(i in 1:N)
		{
			for(j in which(adj[i,]!=1))
			{
       			if(distG[i,k]+distG[k,j]<distG[i,j])
       			{
					distG[i,j]=distG[i,k]+distG[k,j]
				}
			}
		}
		if((i%%100)==0) print(paste(i, "out of ", N, " completed\n", sep = ''))
		flush.console()
	}
	diag(distG) <- 0
	
	return(distG)
}
#***************************************************************************************************
#***************************************************************************************************



#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance matrix and an an integer number indicating the number of    **
#** eigen-thingies to retain, and performs multidimensional scaling                               **
#**                                                                                               **
#** Called by: Isomap                                                                             **
#** arguments: dist - a distance matrix; k - the number of eigenvalues/eigenvectors to retain     **                                
#** output: evals - first k eigenvalues; evecs - first k eigenvectors;                            **
#**         embedding - a representation of the n points forming 'dist' in the new, k-dimensional **
#**                     coordinate system; allevals - the full set of n eigenvalues;              ** 
#**                     allevecs - the full set of n eigenvectors                                 **                                                 
#***************************************************************************************************
#***************************************************************************************************
mds<-function(D,k=NULL)
{
		A=-1/2*D^2
                n=nrow(A)
                one=as.vector(rep(1,n))
                H=diag(n)-1/n*one%*%t(one)
                B=H%*%A%*%H
         
		ev=eigen(B,symmetric=TRUE)
		evals=ev$values
		evecs=ev$vectors
                if(is.null(k)) k=sum(evals>0)
        evals=evals[1:k]
        evecs=evecs[,1:k]

        neweucpts=t(t(evecs)*sqrt(evals))  # the points transformed in the new 2-dimensional space.

        return(list(evals=evals,evecs=evecs,embedding=neweucpts,allevals=ev$values,allevecs=ev$vectors))
}
#***************************************************************************************************
#***************************************************************************************************




#***************************************************************************************************
#***************************************************************************************************
#** This function receives a distance matrix, an original set of points, a new set of points,     **
#** a closeness parameter and the MDS of the original distance matrix, and embeds the new points  ** 
#** in the k-dimensional embedding space of the original set, using Nystrom's approximation       **
#**                                                                                               **
#** arguments: Dorig - a distance matrix; origpts - the original set of points;                   ** 
#**            newpts - the new set of points;                                                    **
#**             delta - the closeness parameter for the nearest neighbors;                        **
#**             evals - the eigenvalues from the MDS performed on the original points;            **
#**             evecs - the eigenvectors from the MDS performed on the original points            **     
#** output: the embedding of the new points in the new coordinate system                          **                                                 
#***************************************************************************************************
#***************************************************************************************************
mdsemapping<-function(Dorig,origpts,newpts,delta,evals,evecs)
{
  k=length(evals) # dimesion of the embedding space
  n=nrow(Dorig)
  N=nrow(newpts)
  mu=Dorig^2%*%rep(1,n)/n # the mean vector

  # Calculate approx geo distances from all N new points to orig n points
  distance=dist(newpts,origpts) #N rows x n cols
  Dnew=matrix(0,nrow=N,ncol=n)
  for(i in 1:N)
  {
    di=sort(distance[i,],index.return=TRUE)
    knearest=di$ix[di$x<=delta]
    kndists=distance[i,knearest]
    if(kndists[1]==0) #our i'th point *is* a landmark point, so we just take the appropriate
                      #part out of our original geodesic distance matrix
      Dnew[i,]=Dorig[knearest[1],]
    else
    for(j in 1:n)
        Dnew[i,j]=min(Dorig[j,knearest]+kndists) # shortest path from a new point to an old point
                                                 # is approximated by the shortest path from the
                                                 # old point to the old points which are k-nearest
                                                 # to the new point, plus the distance to that
                                                 # k-nearest old point to the new point.
    if(i%%10==0) print(paste(i, " Completed", sep = ''))
    flush.console()
  }

  Dnew=Dnew^2 #squared dist matrix.
  Lp=t(evecs)*(1/sqrt(evals))
  embedding=matrix(0,nrow=k,ncol=N)
  for(j in 1:N) embedding[,j]=-Lp%*%(Dnew[j,]-mu)/2

  t(embedding)  #return the new embedding.
}
#***************************************************************************************************
#***************************************************************************************************




#***************************************************************************************************
#***************************************************************************************************
#** This function receives a set of points, a closeness parameter and the embedding dimension,    **
#** and performs MDS on the matrix of approximated geodesic distances between said points         ** 
#**                                                                                               **
#** Calls: makeGdelta; floyd; mds                                                                 **
#** arguments: OrigPts - the original set of points, in the Euclidean space;                      ** 
#**            epsilon - the distance between points in the original grid;                        ** 
#**            k - the embedding dimension                                                        **                                                                 
#** output: the embedding of the points in the new coordinate system                              **                                                 
#***************************************************************************************************
#***************************************************************************************************
Isomap <- function(OrigPts, epsilon, k=NA, n.neighbors=NA, boundary=NA)
{
	d=as.matrix(dist(OrigPts))
	if(is.na(n.neighbors))
	a=makeGdelta(d,sqrt(epsilon^2+epsilon^2+epsilon^2)) # 0.1 distance between the candidate points
	else 
	a=makeGknn(distance, n.neighbors, boundary)
	
	gd=floyd(d,a)
	if(is.na(k)) m=mds(gd)
	else m=mds(gd,k)
	EmbedPts=m$embedding
	Lambda <- m$evals
	V <- m$evecs

	return(list(EmbedPts=EmbedPts, Lambda=Lambda, V=V, GeoDist = gd))
}
#***************************************************************************************************
#***************************************************************************************************




#***************************************************************************************************
#***************************************************************************************************
#** This function receives a matrix of integration points, a matrix of grid points, their MDS     **
#** Object and closeness parameter, and returns the embedding of the integration points           ** 
#**                                                                                               **
#** Calls: mdsemapping; calcdist2                                                                 **
#** arguments: Integ_Pts - points used for numerical integration, in Euclidean representation;    ** 
#**            Orig_Pts - the original grid;                                                      **                                                                 
#**            Orig_MDS - the MDS object of the original grid;                                    ** 
#**            epsilon - the distance between points in the original grid;                        **                                                                 
#** output: the embedding of the integration points in the new coordinate system                  **                                                 
#***************************************************************************************************
#***************************************************************************************************
Embed_Integ_Pts <- function(Integ_Pts, Orig_Pts, Orig_MDS, epsilon)
{
	geoD = as.matrix(dist(Orig_MDS$EmbedPts))  #our *approximated* geodesic distance
	emb <- mdsemapping(geoD, Orig_Pts, Integ_Pts, sqrt(epsilon^2+epsilon^2),Orig_MDS$Lambda,Orig_MDS$V)
	dis <- dist(emb,Orig_MDS$EmbedPts)
	
	return(list(Embed_Integ_Pts = emb, Dist_Orig_Integ = dis))
}
#***************************************************************************************************
#***************************************************************************************************

