library(plgp)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(plotfunctions)
rm(list=ls())
par(mfrow=c(1,1))


#--------------------------design size and par setting--------------------------
sigma2.ini=1
rho.ini1=0.6
rho.ini2=0.2
grid <- seq(0,1,by=0.01)
nexc <- N <- length(grid)

#------------------------------Random starting design---------------------------
#actin <- 51
#actin=c(11,51,91)
set.seed(123)
actin=sample(1:N,3) 

est.rhos <- c(0.6,0.2)
indd <- 0
#---------------------------------plot (and construct) starting design-----------------

plot(seq(1:nexc),seq(0,1,length=nexc), type = "n",xlab="iter",ylab="Incremental design",
     main = bquote( paste("size:", .(nexc),", ","start:", .(paste(grid[actin], collapse=", "))
                          ,", ", rho, ":", .(paste(round(est.rhos,digits = 2), collapse=", ")),", ", sigma^2, ":",.(sigma2.ini)) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


CRIT.vec <- numeric(nexc)

KL.crit.fun <- function(rho,sigma2,des){
  
  n <- length(des)
  DisM <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      DisM[i,j] <- abs(des[i]-des[j]) 
    }
  }
  
  
  
  kernel1 <- sigma2*exp(-DisM/rho[1])
  kernel2 <- sigma2*(matrix(rep(1,n*n),n,n)+DisM/rho[2]) * exp(-DisM/rho[2])
  
  
  kl.crit <- 1/4*( sum(diag(solve(kernel1) %*% kernel2)) + sum(diag(solve(kernel2) %*% kernel1)) ) - n/2
  kl.crit
}


n1 <- length(actin)
xstart <- grid[actin]
xcan <- grid[-actin]
#xcan <- sample(xcan)
CRIT.vec[1:n1] <- KL.crit.fun(est.rhos,sigma2.ini,xstart)


for(ii in 1:n1){
  points(ii, xstart[ii],col="darkblue", pch=16,cex=1)
}

#------------constructing the rest (of design) up to the required size----------

nleft <- nexc-n1
for( k in 1: nleft){
  
  nc <- length(xstart)
  Ncan <- length(xcan)
  
  
  ad.crit <- numeric(Ncan)
  for(m in 1: Ncan){
    test.des <- c(xstart,xcan[m])
    ad.crit[m] <- KL.crit.fun(est.rhos,sigma2.ini,test.des)
  }
  ad.ind <- which.max(ad.crit)
  ad.point <- xcan[ad.ind]
  xhat <- ad.point
  xstart <- c(xstart,ad.point)
  xcan <- xcan[-ad.ind]
  xcan <- sample(xcan)
  CRIT.vec[k+n1] <- KL.crit.fun(est.rhos,sigma2.ini,xstart)
  
  if(CRIT.vec[k+n1]>CRIT.vec[k+n1-1]){
    points(k+n1,xhat , type = "p",col="darkblue", pch=16,cex=1 )   
  }else{
    points(k+n1,xhat , type = "p",col="red", pch=16,cex=1 ) 
    indd <- indd+1
  }
  
}

lines(seq(1:N),xstart,col="black", pch=16,cex=1)