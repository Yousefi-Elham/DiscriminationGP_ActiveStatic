#exchange

library(plgp)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(plotfunctions)
rm(list=ls())
set.seed(123456789)
par(mfrow=c(1,1))
#--------------------------------------------------------------------------------------

grid <- seq(0,1,by=0.01)
N <- length(grid)

#---------------------------------set pars and starting point--------------------------
sigma2=2
rho=3
#actin=c(11,51,91)
#actin=c(11,31,51,71,91)
#actin=c(11,31,71,91)
actin=sample(1:N,3)
#actin=51
nexc=30
indd <- 0
XMAT <- matrix(0,nexc,ncol=nexc)
AD.PN <- DR.PN <- c()
epsilon <- 0.03
#---------------------------------plot (and construct) starting design-----------------

CRIT.vec <- numeric(nexc)


plot(seq(1:nexc),seq(0,1,length=nexc), type = "n",xlab="iter",ylab="Exact design",
     main = bquote( paste("size:", .(nexc),", ","starting design:", .(paste(grid[actin], collapse=", "))
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
#xcan <- sample(xcan)
CRIT.vec[1:n1] <- KL.crit.fun(xstart)


for(ii in 1:n1){
  points(ii, xstart[ii],col="blue", pch=16,cex=1)
}

#----------------------cunstructing the rest (of design) up to the required size-------------------

nleft <- nexc-n1
for( k in 1: nleft){
  
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
  
  if(CRIT.vec[k+n1]>CRIT.vec[k+n1-1]){
    points(k+n1,xhat , type = "p",col="darkblue", pch=16,cex=1 )   
  }else{
    points(k+n1,xhat , type = "p",col="red", pch=16,cex=1 ) 
    indd <- indd+1
  }
  
}
#lines(seq(1:nexc),xstart,col="black", pch=16,cex=1)


XMAT[1,] <- xstart

plot(c(0:5),seq(0,1,length=6), type = "n",xlab="exchanges, 0 is before exchange",ylab="design",
     main = bquote( paste("size:", .(nexc),", ","starting design:", .(paste(grid[actin], collapse=", "))
                          ,", ", rho, ":", .(rho),", ", sigma^2, ":",.(sigma2)) ),pch=2,col="blue");

abline(h=seq(0,1,length=101),col="lightgray")

points(rep(0,nexc),sort(xstart), type = "p",pch=16,col="blue");
text(0,1,"Exc: 0",col="red")

#-------------------------------------exchange algorithm----------------------------------------

CRIT.val <- KL.crit.fun(xstart)

exnum <- 0
finish.all <- FALSE
while(!finish.all){
  
  nc <- length(xstart)
  dr.crit <- numeric(nc)
  for(k in 1: nc){
    dr.crit[k] <- KL.crit.fun(xstart[-k])
  }
  dr.ind <- which.max(dr.crit)
  dr.point <- xstart[dr.ind]
  xcan <- c(xcan, dr.point)
  xstart <- xstart[-dr.ind]
  
  
  
  Ncan <- length(xcan)
  ad.crit <- numeric(Ncan)
  for(m in 1: Ncan){
    test.des <- c(xstart,xcan[m])
    ad.crit[m] <- KL.crit.fun(test.des)
  }
  ad.ind <- which.max(ad.crit)
  ad.point <- xcan[ad.ind]
  
  xstart <- c(xstart,ad.point)
  xcan <- xcan[-ad.ind]
  
  CRIT.exc <- KL.crit.fun(xstart)
  
  if(CRIT.exc > CRIT.val & dr.point!=ad.point){
    #if(round(CRIT.exc,1) > round(CRIT.val,1)){
    exnum <- exnum+1
    CRIT.val <- CRIT.exc 
    print(c(dr.point,ad.point))
    points(rep(exnum,nexc),xstart, type = "p",pch=16,col="blue");
    text(exnum,1-exnum*epsilon,bquote(paste( "Exc: ",.(exnum),", ", .(dr.point)," :droped, ", .(ad.point)," :added")), col="red")
    XMAT[exnum+1,] <- xstart
    DR.PN[exnum] <- dr.point
    AD.PN[exnum] <- ad.point
  }else{
    finish.all <- TRUE
  }
  
}

xstart
exnum
#-------------------------------------plotting final points----------------------------------------

plot(seq(1:nexc),seq(0,1,length=nexc), type = "n",xlab="iter",ylab="Exact design",
     main = bquote( paste("size:", .(nexc),", ","starting design:", .(paste(grid[actin], collapse=", "))
                          ,", ", rho, ":", .(rho),", ", sigma^2, ":",.(sigma2)) ),pch=2,col="blue");

points(seq(1:nexc),xstart,col="darkblue", pch=16,cex=1)
#lines(seq(1:nexc),xstart,col="black", pch=16,cex=1)

#####################################################################################################
#####################################################################################################

max.num <- exnum+1


plot(c(0:max.num),seq(0,1,length=max.num+1), type = "n",xlab="exchanges, 0 is before exchange",ylab="design",
     main = bquote( paste("size:", .(nexc),", ","starting design:", .(paste(grid[actin], collapse=", "))
                          ,", ", rho, ":", .(rho),", ", sigma^2, ":",.(sigma2)) ),pch=2,col="blue");
#xtick=seq(0,max.num,by=1)
ytick=seq(0.1,0.9,by=0.2)
#axis(1,at=xtick,labels=F)
axis(2,at=ytick,labels=T)
ytick1=seq(0.05,0.95,by=0.1)
axis(2,at=ytick1,labels=F)
abline(h=seq(0,1,length=21),col="gray")

points(rep(0,nexc),XMAT[1,], type = "p",pch=16,col="blue");
text(0,1,"Exc: 0",col="red")
if(exnum > 0){
  for(kk in 1:exnum){
    points(rep(kk,nexc),XMAT[kk+1,], type = "p",pch=16,col="blue");
    #points(kk,DR.PN[kk], type = "p",pch=16,col="red")
    points(kk,AD.PN[kk], type = "p",pch=16,col="green")
    text(kk,1-kk*epsilon,bquote(paste( "Exc: ",.(kk),", ", .(DR.PN[kk])," :dropped, ", .(AD.PN[kk])," :added")), col="red")
  }
}

