
# New 3D horseshoe region
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

horseshoe3dsurf<-function(candidates)
{
  n=nrow(candidates)
  eps=.01 # works with .025 grid
  idx=as.vector(rep(FALSE,n))

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
    di=abs(cands[i,2]-.5)
    if(abs(cands[i,1])<=1 && abs(cands[i,1])>=.5 && 
      ((di>=.05-eps && di<=.05+eps) || (di>=.4-eps && di<=.4+eps)))
      idx[i]=TRUE
    else if(abs(cands[i,1])<=1+eps && abs(cands[i,1])>=1-eps &&
            di>=.05-eps && di <=.4-eps)
      idx[i]=TRUE
    else if(cands[i,1]<=.5 && cands[i,1]>=0)
    {
      di=sqrt((cands[i,1]-.5)^2+(cands[i,2]-.5)^2)
      if( (di>=.4-eps && di <=.4+eps) ||
          (di<=.05+eps && di >=.05-eps))
        idx[i]=TRUE
    }
    else if(cands[i,1]>= -.5 && cands[i,1]<=0)
    {
      di=sqrt((cands[i,1]+.5)^2+(cands[i,2]-.5)^2)
      if( (di>=.4-eps && di<=.4+eps) ||
          (di<=.05+eps && di>=.05-eps))
          idx[i]=TRUE
    }
    else
    {
      di=abs(cands[i,2]-.5)
      if((di<=.05-eps && di>=.05-2*eps) ||
         (di>=.4+eps && di <=.4+2*eps))
        idx[i]=TRUE
    }
  }
  newcands=candidates[idx,]
  newcands[,1]=newcands[,1]/2+0.5 #move it to [0,1]

  return(list(newcands=newcands,idx=idx))
}


# Example:
#
cands=as.matrix(expand.grid(seq(-1,1,.075),seq(0,1,.075),seq(0,1,.075)))
hs=horseshoe3d(cands)$newcands

pdf("horseshoe3d_075.pdf",width=5,height=10)
par(mfrow=c(3,1))
plot(hs[,1:2],ylim=c(0,1),xlab=expression(X[1]),ylab=expression(X[2]))
plot(hs[,c(1,3)],xlab=expression(X[1]),ylab=expression(X[3]))
plot(hs[,c(2,3)],xlab=expression(X[2]),ylab=expression(X[3]))
dev.off()

# plot surface
surfcands=as.matrix(expand.grid(seq(-1,1,.025),seq(0,1,.025),seq(0,1,.025)))
hssurf=horseshoe3dsurf(surfcands)$newcands
plot3d(hssurf,xlab="X1",ylab="X2",zlab="X3")
#rgl.snapshot("hssurf.png")

# Generate settings for simulation study
# Settings are the same for both horseshoe simstudy and annular
# ring simstudy.
#
n=c(8,16,24)
num_n=length(n)
K=5 #embedding dimension
run=1:1
num_rho_settings=5
rho.true=matrix(NA,nrow=12,ncol=5)
rho.true[1:3,]=rep(.5,K)
rho.true[4:6,]=rep(.1,K)
rho.true[7:9,]=rep(.02,K)
rho.true[10:12,]=matrix(rep(c(.5,.02,rep(.5,2),.02),3), ncol=5, byrow=T)
rho.true <- cbind(rep(n,4), rho.true)

wd <- "/Users/oharari/Dropbox/Research/SFU/Non_Convex_Regions/New 3D Horsehsoe/Simulation"
setwd(wd)
rho.true <- data.frame(rho.true)
names(rho.true) <- c("n", paste("rho", 1:5, sep=''))
write.table(rho.true, "Simulation_Design_List.txt")


