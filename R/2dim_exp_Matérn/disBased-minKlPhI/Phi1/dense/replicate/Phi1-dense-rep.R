#Distance based


library(mvtnorm)
library(dplyr) #anti_join
library(mgcv)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)
tau2 <- 0.01


NRUN <- 100
nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 51)
x2 <- seq(0,10,length.out = 51)
grid <- expand.grid(x1,x2)
N <- nrow(grid)

actin=c(1,51,2551,2601)
#grid[actin,]

n1 <- length(actin)
xstart <- grid[actin,]
xstart
xcan <- grid
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
  
  L1 <- kernel1 + tau2*diag(n)
  L2 <- kernel2 + tau2*diag(n)
  
  if(model==0) y <- rmvnorm(1,sigma=L1)
  if(model==1) y <- rmvnorm(1,sigma=L2)
  
  return(list(y=y,DisM=DisM,L1=L1,L2=L2))
}

#--------------------------------------------------------------------------------------
max.likelihood.fun <- function(theta,sigma2,des,y){
  
  
  theta1=theta[1]
  theta2=theta[2]
  
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta1) * exp(-sqrt(3)*DisM*theta1)
  L1 <- kernel1 + tau2*diag(n)
  log.lik1 <- -n/2*log(2*pi) - 1/2*log(det(L1)) - 1/2*y%*%solve(L1)%*%t(y)
  
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta2 + 5*DisM^2*theta2^2/3 ) * exp(-sqrt(5)*DisM*theta2)
  L2 <- kernel2 + tau2*diag(n)
  log.lik2 <- -n/2*log(2*pi) - 1/2*log(det(L2)) - 1/2*y%*%solve(L2)%*%t(y)
  
  
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
    Phi1 <- 0 
  }else{
    
    DisM <- as.matrix(dist(des))
    
    kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
    kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
    
    L1 <- kernel1 + tau2*diag(n)
    L2 <- kernel2 + tau2*diag(n)
    
    abs.differ <- abs(L1 - L2)
    #In <- rep(1,n)
    #Phi1 <- In %*% abs.differ %*% In
    Phi1 <- sum(abs.differ) #same formulation
    
  }
  
  return(Phi1)
  
}
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
kernel.inv.fun <- function(theta, sigma2, des) {
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  
  #kernel1.inv <- solve(kernel1)
  #kernel2.inv <- solve(kernel2)
  kernel1.n <- kernel1+ tau2*diag(n)
  kernel2.n <- kernel2+ tau2*diag(n)
  
  return(list(kernel1.n = kernel1.n, kernel2.n = kernel2.n))
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
crit.fun <- function(theta,sigma2,des,cov.past.1,cov.past.2){
  
  
  n <- nrow(des)
  DisVec <- sqrt( (des[1:(n-1),1]-des[n,1])^2 + (des[1:(n-1),2]-des[n,2])^2 )
  
  cov.past.new.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
  cov.past.new.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
  
  
  kernel1.n <- kernel1.inv <- matrix(0,n,n)
  kernel1.n[1:(n-1),1:(n-1)] <- cov.past.1
  kernel1.n[n,1:(n-1)] <- kernel1.n[1:(n-1),n] <- cov.past.new.1
  kernel1.n[n,n] <- sigma2 + tau2
  
  #S1inv.s12 <- cov.past.inv.1 %*% cov.past.new.1
  #B1 <- drop(sigma2 - cov.past.new.1 %*% S1inv.s12)
  #kernel1.inv[1:(n-1),1:(n-1)] <- cov.past.inv.1 + 1/B1 * S1inv.s12 %*% t(S1inv.s12)
  #kernel1.inv[n,1:(n-1)] <- kernel1.inv[1:(n-1),n] <- -1/B1*S1inv.s12
  #kernel1.inv[n,n] <- 1/B1
  #---------------
  kernel2.n <- kernel2.inv <- matrix(0,n,n)
  kernel2.n[1:(n-1),1:(n-1)] <- cov.past.2
  kernel2.n[n,1:(n-1)] <- kernel2.n[1:(n-1),n] <- cov.past.new.2
  kernel2.n[n,n] <- sigma2 + tau2
  
  #S2inv.s12 <- cov.past.inv.2 %*% cov.past.new.2
  #B2 <- drop(sigma2 - cov.past.new.2 %*% S2inv.s12)
  #kernel2.inv[1:(n-1),1:(n-1)] <- cov.past.inv.2 + 1/B2* S2inv.s12 %*% t(S2inv.s12)
  #kernel2.inv[n,1:(n-1)] <- kernel2.inv[1:(n-1),n] <- -1/B2*S2inv.s12
  #kernel2.inv[n,n] <- 1/B2
  
  
  abs.differ <- abs(kernel1.n - kernel2.n)
  Phi1 <- sum(abs.differ) 
  
  
  return(Phi1)
  
}
#-------------------------------------------------------------------------------------------------------
pred.error.fun <- function(theta,sigma2,pastin,leftin){
  
  n1 <- length(pastin)
  nleft <- length(leftin)
  
  DisM <- as.matrix(dist(grid))
  kernel1.grid <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2.grid <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  L1.grid <- kernel1.grid + tau2*diag(N)
  L2.grid <- kernel2.grid + tau2*diag(N)
  
  K1 <- L1.grid[pastin,pastin]
  K2 <- L2.grid[pastin,pastin]
  K1.inv <- solve(K1)
  K2.inv <- solve(K2)
  
  k1mat <- L1.grid[leftin,pastin]
  k2mat <- L2.grid[leftin,pastin]
  
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
#kernel1.inv <- kernel.inv.list$kernel1.inv
#kernel2.inv <- kernel.inv.list$kernel2.inv

crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n)

crit.fun.old(est.theta,sigma2.ini,xstart)

CRIT.vec[1:n1] <- crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n)
CRIT.vec[1:n1]
#----------------------cunstructing the rest (of design) up to the required size-------------------


set.seed(123456789)
for( k in 1: nleft){
  
  print(k)
  
  nc <- nrow(xstart)
  Ncan <- nrow(xcan)
  
  kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart)
  kernel1.n <- kernel.inv.list$kernel1.n
  kernel2.n <- kernel.inv.list$kernel2.n
  #kernel1.inv <- kernel.inv.list$kernel1.inv
  #kernel2.inv <- kernel.inv.list$kernel2.inv
  
  
  
  ad.crit <- numeric(Ncan)
  for(m in 1: Ncan){
    test.des <- rbind(xstart,xcan[m,])
    ad.crit[m] <- crit.fun(est.theta,sigma2.ini,test.des,kernel1.n,kernel2.n)
    #crit.fun.old(est.theta,sigma2.ini,test.des)
  }
  
  ad.ind <- which.max(ad.crit)
  ad.point <- xcan[ad.ind,]
  xstart <- rbind(xstart,ad.point)
  
  CRIT.vec[k+n1] <- max(ad.crit)
  #--------------------------------
  CRIT.val <- max(ad.crit) 
  #CRIT.val <-  crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n)
  
  
  exnum <- 0
  finish.all <- FALSE
  
  while(!finish.all){
    
    nc1 <- nrow(xstart)
    dr.crit <- numeric(nc1)
    
    for(a in 1: nc1){
      dr.crit[a] <- crit.fun.old(est.theta,sigma2.ini,xstart[-a,])
      #we can't use the crit.fun here while the dimensions (for past covariances) match but the designs not!
    }
    dr.ind <- which.max(dr.crit)
    dr.point <- xstart[dr.ind,]
    xstart.dr <- xstart[-dr.ind,]
    
    
    Ncan1 <- nrow(xcan)
    add.crit <- numeric(Ncan1)
    
    kernel.inv.list.dr <- kernel.inv.fun(est.theta,sigma2.ini,xstart.dr)
    kernel1.n.dr <- kernel.inv.list.dr$kernel1.n
    kernel2.n.dr <- kernel.inv.list.dr$kernel2.n
    #kernel1.inv.dr <- kernel.inv.list.dr$kernel1.inv
    #kernel2.inv.dr <- kernel.inv.list.dr$kernel2.inv
    
    for(b in 1: Ncan1){
      test.des <- rbind(xstart.dr,xcan[b,])
      add.crit[b] <- crit.fun(est.theta,sigma2.ini,test.des,kernel1.n.dr,kernel2.n.dr)
      ##crit.fun.old(est.theta,sigma2.ini,test.des)
      
    }
    add.ind <- which.max(add.crit)
    add.point <- xcan[add.ind,]
    xstart.ad <- rbind(xstart.dr,add.point)
    
    CRIT.exc <- max(add.crit)
    
    
    if( CRIT.exc > CRIT.val & dr.point[,1]!=add.point[,1] & dr.point[,2]!=add.point[,2]){
      xstart <- xstart.ad
      exnum <- exnum+1
      CRIT.val <- CRIT.exc 
      print(rbind(dr.point,add.point))
      #as.numeric(row.names(xstart))
    }else{
      finish.all <- TRUE
    }
    
  }
  
  CRIT.vec[k+n1] <- crit.fun.old(est.theta,sigma2.ini,xstart)
  
  
  #actin <- as.numeric(row.names(xstart))
  #canin <- 1:N
  
  #pred.error.list <- pred.error.fun(est.theta,sigma2.ini,actin,canin)
  #E11[k] <- pred.error.list$e1.1
  #E12[k] <- pred.error.list$e1.2
  #E22[k] <- pred.error.list$e2.2
  #E21[k] <- pred.error.list$e2.1
  #E1.dif[k] <- pred.error.list$e1.dif
  #E2.dif[k] <- pred.error.list$e2.dif
  #SUM.E1E2[k] <- pred.error.list$sum.e1e2
  
}


#-----------------------------------------------------------------------------------------------
as.numeric(row.names(xstart))

nrow(xstart)
nrow(xstart)-n1 
nleft


#E11 <- E11[1:(nrow(xstart)-n1)]
#E12 <- E12[1:(nrow(xstart)-n1)]
#E22 <- E22[1:(nrow(xstart)-n1)]
#E21 <- E21[1:(nrow(xstart)-n1)]
#E1.dif <- E1.dif[1:(nrow(xstart)-n1)]
#E2.dif <- E2.dif[1:(nrow(xstart)-n1)]
#SUM.E1E2 <- SUM.E1E2[1:(nrow(xstart)-n1)]


#int.hit.n <- c(1,2,3,4,5,6,16,36,76,116,nrow(xstart)-n1)
#int.hit.n

#E11[int.hit.n]
#E12[int.hit.n]
#E22[int.hit.n]
#E21[int.hit.n]
#E1.dif[int.hit.n]
#E2.dif[int.hit.n]
#SUM.E1E2[int.hit.n]

#round(E11[int.hit.n],digits=5)
#round(E12[int.hit.n],digits=5)
#round(E22[int.hit.n],digits=5)
#round(E21[int.hit.n],digits=5)
#round(E1.dif[int.hit.n],digits=5)
#round(E2.dif[int.hit.n],digits=5)
#round(SUM.E1E2[int.hit.n],digits=5)
##########################################################
#sel.iter <- c(1,2,3,4,5,6,16,26,36,46)

#E11[sel.iter]
#E12[sel.iter]
#E22[sel.iter]
#E21[sel.iter]
#E1.dif[sel.iter]
#E2.dif[sel.iter]
#SUM.E1E2[sel.iter]


#round(E11[sel.iter],digits=5)
#round(E12[sel.iter],digits=5)
#round(E22[sel.iter],digits=5)
#round(E21[sel.iter],digits=5)
#round(E1.dif[sel.iter],digits=5)
#round(E2.dif[sel.iter],digits=5)
#round(SUM.E1E2[sel.iter],digits=5)
##########################################################

length(CRIT.vec)
CRIT.vec <- CRIT.vec[1:nrow(xstart)]
length(CRIT.vec)
CRIT.vec
xstart
#####################################################################################################
nrow(xstart)
design1.uniq <- uniquecombs(xstart)
counter1=attr(design1.uniq,"index")
counter2=unique(sort(counter1))
replic1 <- replic2 <- numeric(nrow(design1.uniq))
#----------------------------------------------------------------
#counter <- counter2 <- numeric(length(design1.uniq))
for(z in 1:length(counter2)){
  x.replic <- which(counter1==counter2[z])
  replic1[z] <- length(x.replic)
  replic2[z] <- length(x.replic)/nrow(xstart)
  
}

design.opt <- cbind(design1.uniq,replic1,replic2)
design.opt
sum(design.opt[,3])
sum(design.opt[,4])

#####################################################################################################
subdes <- data.frame(x=xstart[1:6,1],y=xstart[1:6,2], label=as.character(1:6))

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/subdesign.jpeg")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(subdes$x,subdes$y, type = "p",pch=16,col= "black");
text(subdes$x,subdes$y,lables=subdes$label,pos=1,col = "red")

dev.off()    
#----------------------------------------------------------------
setEPS()     
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/subdesign.eps")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(subdes$x,subdes$y, type = "p",pch=16,col= "black");
text(subdes$x,subdes$y,lables=subdes$label,pos=1,col = "red")

dev.off()
#----------------------------------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/design.jpeg")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(design.opt[,1],design.opt[,2],cex=0.3*design.opt[,3], type = "p",pch=16,col= "black")

dev.off()
#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/design.eps")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(design.opt[,1],design.opt[,2],cex=0.3*design.opt[,3], type = "p",pch=16,col= "black")

dev.off()
##############################################################################################
max(CRIT.vec)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/criterionval.jpeg")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vec),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vec),by=200)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vec,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/criterionval.eps")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vec),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vec),by=200)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vec,col="blue",lwd = 2)

dev.off()

