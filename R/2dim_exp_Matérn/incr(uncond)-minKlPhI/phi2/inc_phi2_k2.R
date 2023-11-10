#Incremental unconditional

library(mvtnorm)
rm(list=ls())
par(mfrow=c(1,1))

#---------------------------------generate the obs from under each kernel--------------------------
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)

NRUN <- 100
nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)
N <- nrow(grid)

actin=c(1,25,601,625)
#grid[actin,]

n1 <- length(actin)
xstart <- grid[actin,]
xstart
xcan <- grid[-actin,]
nleft <- nexc-n1

M.COUNT <- matrix(0,nleft,NRUN)
dim(M.COUNT)

E11 <- E12 <- E22 <- E21 <- E1.dif <- E2.dif <- SUM.E1E2 <- numeric(length=nleft)

CRIT.vec <- numeric(nexc)
#--------------------------------------------------------------------------------------
generate.y.kernels <- function(theta,sigma2,des,model){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  if(model==0) y <- rmvnorm(1,sigma=kernel1)
  if(model==1) y <- rmvnorm(1,sigma=kernel2)
  
  return(list(y=y,DisM=DisM,kernel1=kernel1,kernel2=kernel2))
  
}

#--------------------------------------------------------------------------------------
max.likelihood.fun <- function(theta,sigma2,des,y){
  
  theta1=theta[1]
  theta2=theta[2]
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta1) * exp(-sqrt(3)*DisM*theta1)
  log.lik1 <- -n/2*log(2*pi) - 1/2*log(det(kernel1)) - 1/2*y%*%solve(kernel1)%*%t(y)
  
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta2 + 5*DisM^2*theta2^2/3 ) * exp(-sqrt(5)*DisM*theta2)
  log.lik2 <- -n/2*log(2*pi) - 1/2*log(det(kernel2)) - 1/2*y%*%solve(kernel2)%*%t(y)
  
  
  if(abs(log.lik1 - log.lik2) < 1e-6){
    lik.res <- sample(c(0, 1), 1)
  }else{
    lik.res <- as.numeric(log.lik2 > log.lik1)
  }
  
  # lik.res is  0 or 1 (prediction for which model is the TRUE model)
  return(lik.res)
}
#--------------------------------------------------------------------------------------
crit.fun.old <- function(theta,sigma2,des){
  
  n <- nrow(des)
  
  if(n==1){
    phi.crit <- 0 
  }else{
    
    DisM <- as.matrix(dist(des))
    
    kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
    kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
    
    cov.past.1 <- kernel1[1:(n-1),1:(n-1)]
    #cov.past.inv.1 <- solve(cov.past.1)
    #--------------------
    cov.past.2 <- kernel2[1:(n-1),1:(n-1)]
    #cov.past.inv.2 <- solve(cov.past.2)
    
    if(det(cov.past.1)==0 | det(cov.past.2)==0){
      
      phi.crit <- 0
      
      if(det(cov.past.1)==0 & det(cov.past.2)==0){
        warning("both kernels are singular")
      }else{
        if( det(cov.past.1)==0 ){
          warning("first kernel is singular")
        }else{
          warning("second kernel is singular")
        } 
      }
      
      
    }else{
      
      cov.past.inv.1 <- solve(cov.past.1)
      cov.past.new.1 <- kernel1[n,1:(n-1)]
      sigma.new.1 <- kernel1[n,n]
      pred.error.11 <- sigma.new.1 - cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.new.1
      
      cov.past.inv.2 <- solve(cov.past.2)
      cov.past.new.2 <- kernel2[n,1:(n-1)]
      sigma.new.2 <- kernel2[n,n]
      pred.error.22 <- sigma.new.2 - cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.new.2
      
      #----------
      pred.errorP1.12 <- sigma.new.1 + cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.1 %*% cov.past.inv.2 %*% cov.past.new.2 
      pred.errorP2.12 <- -2*cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.new.1
      pred.error.12 <- pred.errorP1.12 + pred.errorP2.12
      
      pred.errorP1.21 <- sigma.new.2 + cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.2 %*% cov.past.inv.1 %*% cov.past.new.1 
      pred.errorP2.21 <- -2*cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.new.2
      pred.error.21 <- pred.errorP1.21 + pred.errorP2.21
      
      #----------
      
      phi.crit <- pred.error.12/pred.error.11 + pred.error.21/pred.error.22 -2
    }
    
    
  }
  
  return(phi.crit)
  
}
#--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
kernel.inv.fun <- function(theta, sigma2, des) {
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  if(det(kernel1)==0 | det(kernel2)==0){
    
    kernel1.inv <- kernel2.inv <- NULL
    
    if(det(kernel1)==0 & det(kernel2)==0){
      warning("both kernels are singular")
    }else{
      if( det(kernel1)==0 ){
        warning("first kernel is singular")
      }else{
        warning("second kernel is singular")
      } 
    }
  } else {
    kernel1.inv <- solve(kernel1)
    kernel2.inv <- solve(kernel2)
    kernel1.n <- kernel1
    kernel2.n <- kernel2
    
  }	
  
  return(list(kernel1.n = kernel1.n, kernel2.n = kernel2.n, kernel1.inv = kernel1.inv, kernel2.inv = kernel2.inv))
}
#-------------------------------------------------------------------------------------------------------
crit.fun <- function(theta,sigma2,des,cov.past.1,cov.past.2,cov.past.inv.1,cov.past.inv.2){
  
  n <- nrow(des)
  
  DisVec <- sqrt( (des[1:(n-1),1]-des[n,1])^2 + (des[1:(n-1),2]-des[n,2])^2 )
  
  cov.past.new.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
  cov.past.new.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
  
  sigma.new.1 <- sigma2
  sigma.new.2 <- sigma2
  
  pred.error.11 <- sigma.new.1 - cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.new.1
  pred.error.22 <- sigma.new.2 - cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.new.2
  
  #----------
  pred.errorP1.12 <- sigma.new.1 + cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.1 %*% cov.past.inv.2 %*% cov.past.new.2 
  pred.errorP2.12 <- -2*cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.new.1
  pred.error.12 <- pred.errorP1.12 + pred.errorP2.12
  
  pred.errorP1.21 <- sigma.new.2 + cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.2 %*% cov.past.inv.1 %*% cov.past.new.1 
  pred.errorP2.21 <- -2*cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.new.2
  pred.error.21 <- pred.errorP1.21 + pred.errorP2.21
  
  #----------
  
  phi.crit <- pred.error.12/pred.error.11 + pred.error.21/pred.error.22 -2
  phi.crit
  
}
#-------------------------------------------------------------------------------------------------------
pred.error.fun <- function(theta,sigma2,pastin,leftin){
  
  n1 <- length(pastin)
  nleft <- length(leftin)
  
  DisM <- as.matrix(dist(grid))
  kernel1.grid <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2.grid <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  K1 <- kernel1.grid[pastin,pastin]
  K2 <- kernel2.grid[pastin,pastin]
  K1.inv <- solve(K1)
  K2.inv <- solve(K2)
  
  k1mat <- kernel1.grid[leftin,pastin]
  k2mat <- kernel2.grid[leftin,pastin]
  
  A1 <- crossprod(k1mat) # equivalent to t(k1mat)%*%k1mat, slightly more efficient
  A2 <- crossprod(k2mat)
  A12 <- crossprod(k1mat,k2mat) # equivalent to t(k1mat)%*%k2mat
  # remark: mat1%*%t(mat2) can be computed with tcrossprod(mat1,mat2)
  
  e1.1 <- sigma2 - sum(K1.inv * A1)/nleft
  e1.2 <- sigma2 + (sum(K2.inv%*%K1%*%K2.inv * A2) - 2*sum(K2.inv * A12))/nleft
  e2.2 <- sigma2 - sum(K2.inv * A2)/nleft
  e2.1 <- sigma2 + (sum(K1.inv%*%K2%*%K1.inv * A1) - 2*sum(K1.inv * A12))/nleft
  
  e1.dif <- e1.2 - e1.1
  e2.dif <- e2.1 - e2.2
  
  sum.e1e2 <- e1.2 - e1.1 + e2.1 - e2.2
  
  return(list(e1.1=e1.1,e1.2=e1.2,e2.2=e2.2,e2.1=e2.1,e1.dif=e1.dif,e2.dif=e2.dif,sum.e1e2=sum.e1e2))
  
}

#-------------------------------------------------------------------------------------------------------
kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart[1:3,])
kernel1.n <- kernel.inv.list$kernel1.n
kernel2.n <- kernel.inv.list$kernel2.n
kernel1.inv <- kernel.inv.list$kernel1.inv
kernel2.inv <- kernel.inv.list$kernel2.inv


CRIT.vec[1:n1] <- crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)
#CRIT.vec[1:n1] <- crit.fun.old(est.theta,sigma2.ini,xstart)
#----------------------cunstructing the rest (of design) up to the required size-------------------
set.seed(123456789)
for( k in 1: nleft){
  
  print(k)
  
  filter.can.index <- c()    
  for(t in 1:nrow(xcan)){
    test.cand <- rbind(xstart,xcan[t,])
    n <- nrow(test.cand)
    DisM <- as.matrix(dist(test.cand))
    
    kernel1 <- sigma2.ini*(1 + sqrt(3)*DisM*est.theta[1]) * exp(-sqrt(3)*DisM*est.theta[1])
    kernel2 <- sigma2.ini*(1 + sqrt(5)*DisM*est.theta[2] + 5*DisM^2*est.theta[2]^2/3 ) * exp(-sqrt(5)*DisM*est.theta[2])
    
    if(det(kernel1)==0 | det(kernel2)==0) {
      filter.can.index <- c(filter.can.index,t)
    }
    
  }
  
  
  if(length(filter.can.index) >= 1){
    xcan <- xcan[-filter.can.index,]
    print(length(filter.can.index))
  }else{
    xcan <- xcan
  }
  
  nc <- nrow(xstart)
  Ncan <- nrow(xcan)
  
  if(Ncan==0){
    xstart <- xstart
    print(k-1)
    stop("No more candidate set left")
  }
  
  kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart)
  kernel1.n <- kernel.inv.list$kernel1.n
  kernel2.n <- kernel.inv.list$kernel2.n
  kernel1.inv <- kernel.inv.list$kernel1.inv
  kernel2.inv <- kernel.inv.list$kernel2.inv
  
  if (is.null(kernel1.inv) || is.null(kernel2.inv)) {
    warning(paste0("current design singular (step ", k, ")"))
    break
  }
  
  ad.crit <- numeric(Ncan)
  for(m in 1: Ncan){
    test.des <- rbind(xstart,xcan[m,])
    ad.crit[m] <- crit.fun(est.theta,sigma2.ini,test.des,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)
  }
  
  ad.ind <- which.max(ad.crit)
  grid.max.ind <- which(grid[,1]==xcan[ad.ind,1]&grid[,2]==xcan[ad.ind,2])
  actin <- c(actin,grid.max.ind)
  canin <- setdiff(1:N, actin)
  
  ad.point <- xcan[ad.ind,]
  #xstart <- rbind(xstart,ad.point)
  #xcan <- xcan[-ad.ind,]
  xstart <- grid[actin,]
  xcan <- grid[canin,]
  
  CRIT.vec[k+n1] <- max(ad.crit)
  
  for(NN in 1:NRUN){
    y.obs.K <- generate.y.kernels(est.theta,sigma2.ini,xstart,1)$y #true pars so fix, only model changes
    M.COUNT[k,NN] <- max.likelihood.fun(est.theta,sigma2.ini,xstart,y.obs.K) 
  }
  
  pred.error.list <- pred.error.fun(est.theta,sigma2.ini,actin,canin)
  E11[k] <- pred.error.list$e1.1
  E12[k] <- pred.error.list$e1.2
  E22[k] <- pred.error.list$e2.2
  E21[k] <- pred.error.list$e2.1
  E1.dif[k] <- pred.error.list$e1.dif
  E2.dif[k] <- pred.error.list$e2.dif
  SUM.E1E2[k] <- pred.error.list$sum.e1e2
  
}

#-----------------------------------------------------------------------------------------------
nrow(xstart)
nrow(xstart)-n1 
nleft
M.COUNT[nleft,]
M.COUNT <- M.COUNT[1:(nrow(xstart)-n1),]

dim(M.COUNT)
M.COUNT[nrow(xstart)-n1,]


av.hit.iter <- apply(M.COUNT,1,mean)
av.hit.run <- apply(M.COUNT,2,mean)

av.hit.iter
av.hit.run

dim(M.COUNT)
nrow(M.COUNT)
length(av.hit.iter)

E11
E12
E22
E21
E1.dif
E2.dif
SUM.E1E2


int.hit.n <- c(1,2,3,4,5,6,16,36,76,116,156,length(av.hit.iter))
int.hit.n
av.hit.iter[int.hit.n]

E11[int.hit.n]
E12[int.hit.n]
E22[int.hit.n]
E21[int.hit.n]
E1.dif[int.hit.n]
E2.dif[int.hit.n]
SUM.E1E2[int.hit.n]

round(E11[int.hit.n],digits=4)
round(E12[int.hit.n],digits=4)
round(E22[int.hit.n],digits=4)
round(E21[int.hit.n],digits=4)
round(E1.dif[int.hit.n],digits=4)
round(E2.dif[int.hit.n],digits=4)
round(SUM.E1E2[int.hit.n],digits=4)


length(CRIT.vec)
CRIT.vec <- CRIT.vec[1:nrow(xstart)]
length(CRIT.vec)
CRIT.vec
xstart
#---------------#---------------
av.sd.iter <- apply(M.COUNT,1,function(x) sqrt((NRUN-1)/NRUN*var(x)))

av.sd.iter

av.sd.iter[int.hit.n]
round(av.sd.iter[int.hit.n],digits = 4)

count2 <- c(1,2,3,4,5,6,16,26,36,46)

av.sd.iter[count2]
round(av.sd.iter[count2],digits = 4)

av.hit.iter[count2]
##############################################################################################
#design <- data.frame(x1=xstart[,1],x2=xstart[,2])

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(xstart[,1],xstart[,2], type = "p",pch=16,col= "black");
points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");

##############################################################################################
#xstart
nrow(xstart)
size.seg <- c(20,40,80,120,160,nrow(xstart))
size.seg
seg <- length(size.seg)
seg


op <- par(mfrow=c(2,3))
for(kk in 1:seg){
  
  plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
       main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(size.seg[kk])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);
  
  xtick <- seq(0,10,by=2)
  ytick=seq(0,10,by=2)
  axis(1,at=xtick,labels=T,cex.axis=1.5)
  axis(2,at=ytick,labels=T,cex.axis=1.5)
  box()
  abline(h=seq(0,10,length=21),col="gray")
  abline(v=seq(0,10,length=21),col="gray")
  points(xstart[1:size.seg[kk],1],xstart[1:size.seg[kk],2], type = "p",pch=16,col= "black");
  points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");
  
}
par(op)

##############################################################################################
# since lower distances are more important (we are interested to analyse those in detail),
# we exclude the 3 starting points to reduce their effect on the distances
nrow(xstart)
size.seg <- c(20,40,80,120,160,nrow(xstart))
size.seg
seg <- length(size.seg)
seg


op <- par(mfrow=c(2,3))
for(l in 1:seg){
  
  Trimat.length <- choose(size.seg[l]-n1,2) #50 is the design size
  Trimat.vec <- numeric(Trimat.length)
  length(Trimat.vec)
  DisM <- as.matrix(dist(xstart[(n1+1):size.seg[l],]))
  Trimat.vec <- DisM[lower.tri(DisM)]
  hist(Trimat.vec ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("size=", .(size.seg[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

#####################################################################################################
max(CRIT.vec)

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vec),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vec),by=0.04)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vec,col="blue",lwd = 2)



