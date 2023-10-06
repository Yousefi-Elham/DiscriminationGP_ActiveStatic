#Incremental

library(plgp)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(plotfunctions)
rm(list=ls())
set.seed(123456789)
par(mfrow=c(1,1))
indd <- 0
#-----------------------------------------------------------

grid <- seq(0,1,by=0.01)
N <- length(grid)


KL.onedim.fun <- function(sigma2,rho,actin){
  
  CRIT.vec <- numeric(N)
  
  
  par(mar=c(5.1, 4.1, 4.1, 4.1))
  plot(seq(1:N),seq(0,1,length=N), type = "n",xlab="iter",ylab="Incremental design",
       main = bquote( paste("Starting design:", .(paste(grid[actin], collapse=", "))
                            ,", ", rho, ":", .(rho),", ", sigma^2, ":",.(sigma2)) ),pch=2,col="blue");
  
  
  
  
  KL.crit.fun <- function(des){
    n <- length(des)
    DisM <- matrix(0,n,n)
    for(i in 1:n){
      for(j in 1:n){
        DisM[i,j] <- abs(des[i]-des[j]) 
      }
    }
    
    kernel1 <- sigma2*exp(-DisM/rho)
    kernel2 <- sigma2*(matrix(rep(1,n*n),n,n)+DisM/rho) * exp(-DisM/rho)
    kl.crit <- 1/4*( sum(diag(solve(kernel1) %*% kernel2)) + sum(diag(solve(kernel2) %*% kernel1)) ) - n/2
    kl.crit
  }
  
  
  n1 <- length(actin)
  xstart <- grid[actin]
  xcan <- grid[-actin]
  xcan <- sample(xcan)
  CRIT.vec[1:n1] <- KL.crit.fun(xstart)
  
  
  for(ii in 1:n1){
    points(ii, xstart[ii],col="darkblue", pch=16,cex=1)
  }
  
  #-----------------------------------------------------------
  
  
  for( k in 1:length(xcan)){
    
    nc <- length(xstart)
    Ncan <- length(xcan)
    
    
    ad.crit <- numeric(Ncan)
    for(m in 1: Ncan){
      test.des <- c(xstart,xcan[m])
      ad.crit[m] <- KL.crit.fun(test.des)
    }
    ad.ind <- which.max(ad.crit)
    ad.point <- xcan[ad.ind]
    xhat <- ad.point
    xstart <- c(xstart,ad.point)
    xcan <- xcan[-ad.ind]
    xcan <- sample(xcan)
    CRIT.vec[k+n1] <- KL.crit.fun(xstart)
    
    if(CRIT.vec[k+n1] > CRIT.vec[k+n1-1]){
      points(k+n1,xhat , type = "p",col="darkblue", pch=16,cex=1 )   
    }else{
      points(k+n1,xhat , type = "p",col="red", pch=16,cex=1 ) 
      indd <- indd+1
    }
    
    
  }
  
  check.incr <- sort(CRIT.vec)- CRIT.vec
  
  lines(seq(1:N),xstart,col="black", pch=16,cex=1)
  
  return(list(xstart=xstart,CRIT.vec= CRIT.vec,check.incr=check.incr,indd=indd)) 
  
}
###########################################################################################
###########################################################################################
#################################### non random start #####################################




try1 <- KL.onedim.fun(sigma2=2,rho=0.6,actin=51)

try1$CRIT.vec
diff(try1$CRIT.vec)

length(diff(try1$CRIT.vec))
plot(seq(1:100),diff(try1$CRIT.vec))


try2 <- KL.onedim.fun(sigma2=0.01,rho=0.2,actin=c(21,81))
try3 <- KL.onedim.fun(sigma2=0.01,rho=0.2,actin=c(21,51,81))



try4 <- KL.onedim.fun(sigma2=2,rho=0.2,actin=51)
try5 <- KL.onedim.fun(sigma2=2,rho=0.2,actin=c(21,81))
try6 <- KL.onedim.fun(sigma2=2,rho=0.6,actin=c(11,51,91))


try7 <- KL.onedim.fun(sigma2=2,rho=0.6,actin=51)
try7$CRIT.vec
diff(try7$CRIT.vec)
#plot(seq(1:100),diff(try7$CRIT.vec))
diff(try7$CRIT.vec)<0
sum(diff(try7$CRIT.vec)<0)
try7$indd

try7$xstart
try7$xstart[2:N][diff(try7$CRIT.vec)<0]


grid[c(11,51,91)]

try8 <- KL.onedim.fun(sigma2=0.01,rho=0.6,actin=c(21,81))
try9 <- KL.onedim.fun(sigma2=0.01,rho=0.6,actin=c(21,51,81))

#################################### random start ####################################
try10 <- KL.onedim.fun(sigma2=0.01,rho=0.2,actin=sample(1:N,1))
try11 <- KL.onedim.fun(sigma2=0.01,rho=0.2,actin=sample(1:N,2))
try12 <- KL.onedim.fun(sigma2=0.01,rho=0.2,actin=sample(1:N,3))



try13 <- KL.onedim.fun(sigma2=0.01,rho=0.6,actin=sample(1:N,1))
try14 <- KL.onedim.fun(sigma2=0.01,rho=0.6,actin=sample(1:N,2))
try15 <- KL.onedim.fun(sigma2=0.01,rho=0.9,actin=sample(1:N,3))
