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
nexc <- 20 # size of sequential design
x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)


actin=c(1,25,601,625)
#grid[actin,]

n1 <- length(actin)
xstart <- grid[actin,]
xstart
xcan <- grid[-actin,]
nleft <- nexc-n1

M.COUNT <- matrix(0,nleft,NRUN)
dim(M.COUNT)

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
kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart[1:3,])
kernel1.n <- kernel.inv.list$kernel1.n
kernel2.n <- kernel.inv.list$kernel2.n
kernel1.inv <- kernel.inv.list$kernel1.inv
kernel2.inv <- kernel.inv.list$kernel2.inv


CRIT.vec[1:n1] <- crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)
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
  ad.point <- xcan[ad.ind,]
  xstart <- rbind(xstart,ad.point)
  xcan <- xcan[-ad.ind,]
  
  CRIT.vec[k+n1] <- max(ad.crit)
  
  for(NN in 1:NRUN){
    y.obs.K <- generate.y.kernels(est.theta,sigma2.ini,xstart,0)$y #true pars so fix, only model changes
    M.COUNT[k,NN] <- max.likelihood.fun(est.theta,sigma2.ini,xstart,y.obs.K) 
  }
  
}

#-----------------------------------------------------------------------------------------------
nrow(xstart)
nrow(xstart)-n1 
nleft
M.COUNT[nleft,]
M.COUNT <- M.COUNT[1:(nrow(xstart)-n1),]

dim(M.COUNT)
M.COUNT[nrow(xstart)-n1,]


M.COUNT.N <- 1-M.COUNT
M.COUNT.N[nrow(xstart)-n1,]

av.hit.iter <- apply(M.COUNT.N,1,mean)

av.hit.iter

##############################################################################################

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(xstart[,1],xstart[,2], type = "p",pch=16,col= "black")
points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue")

