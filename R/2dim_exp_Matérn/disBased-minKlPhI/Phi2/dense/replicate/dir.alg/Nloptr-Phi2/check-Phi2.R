library(mgcv)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)

#NRUN <- 100
#nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 51)
x2 <- seq(0,10,length.out = 51)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

opt.loc.3 <- matrix(round(c(10,0,9.502019,1.858492,-1.858492+10,-9.502019+10),digits=4),3,2,byrow = TRUE)
opt.loc.3

noopt.loc <- matrix(c(10,0,9.2,2,8,0.8),3,2,byrow = TRUE) #with unequal or equal weights, both: 98 points dont fulfill
noopt.loc

opt.w.uneq <- round(c(0.3018613,0.3490694,0.3490694),digits=4) #80
opt.w.uneq

opt.w.eq <- round(c(0.3333333,0.3333333,0.3333333),digits=4) #78
opt.w.eq

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
dir.deriv <- function(theta,sigma2,des,weight){
  
  n <- nrow(des)
  
  Phi.crit <- numeric(N)
  for(i in 1:N){
    
    DisVec <- sqrt( (grid[i,1]-des[,1])^2 + (grid[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)^2*weight) 
    
  }
  return(Phi.crit)
}
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
Phi.optim <- function(theta, sigma2, des, weight) {
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  #abs.differ <- matrix(NA,n,n)
  #for(j in 1:n){
  #  for(k in 1:n){
  #    abs.differ[j,k] <- abs(kernel1[j,k] - kernel2[j,k])^2*weight[j]*weight[k]
  
  
  #  }
  #}
  
  #opt.val <- sum(abs.differ)
  
  opt.val <- sum(abs(kernel1 - kernel2)^2*(weight%*%t(weight)))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------------------------

Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,opt.loc.3,opt.w.eq)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,opt.loc.3,opt.w.eq) 
Phi.optopt #0.0005147528


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#78

length(which(Phi.checkopt ==0))/N
#0.02998847

poi <- grid[!Phi.checkopt , ]
#poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill3.eq.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.3[,1],opt.loc.3[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill3.eq.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.3[,1],opt.loc.3[,2],pch = 17 ,cex = 2, col="red")

dev.off()

#------------------------------#------------------------------#------------------------------
#------------------------------#------------------------------#------------------------------
Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,opt.loc.3,opt.w.uneq)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,opt.loc.3,opt.w.uneq) 
Phi.optopt #0.0005138117


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#80

length(which(Phi.checkopt ==0))/N
#0.0307574

poi <- grid[!Phi.checkopt , ]
#poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill3.uneq.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.3[,1],opt.loc.3[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill3.uneq.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.3[,1],opt.loc.3[,2],pch = 17 ,cex = 2, col="red")

dev.off()

#------------------------------#------------------------------#------------------------------
#------------------------------#------------------------------#------------------------------
opt.loc.ch3 <- matrix(round(c(9,1,9.502019,1.858492,-1.858492+10,-9.502019+10),digits=4),3,2,byrow = TRUE)
opt.loc.ch3 

Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,opt.loc.ch3,opt.w.eq)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,opt.loc.3,opt.w.eq) 
Phi.optopt #0.0005147528


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#70

length(which(Phi.checkopt ==0))/N
#0.0307574

poi <- grid[!Phi.checkopt , ]
#poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfillch3.eq.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.ch3[,1],opt.loc.ch3[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfillch3.eq.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.ch3[,1],opt.loc.ch3[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#------------------------------#------------------------------#------------------------------
#------------------------------#------------------------------#------------------------------
opt.loc.4 <- matrix(c(7.5,2.5,9.502019,1.858492,-1.858492+10,-9.502019+10,10,0),4,2,byrow = TRUE) #89 points
opt.loc.4

opt.w.4 <- rep(1/4,4)
opt.w.4

Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,opt.loc.4,opt.w.4)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,opt.loc.4,opt.w.4) 
Phi.optopt #0.000482456


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#89

length(which(Phi.checkopt ==0))/N
#0.03421761

poi <- grid[!Phi.checkopt , ]
#poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill4.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.4[,1],opt.loc.4[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill4.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.4[,1],opt.loc.4[,2],pch = 17 ,cex = 2, col="red")

dev.off()

#------------------------------#------------------------------#------------------------------
#------------------------------#------------------------------#------------------------------
opt.loc.6 <- matrix(c(9.8,0.2,9.502019,1.858492,9.502019+0.2,1.858492+0.2,-1.858492+10,-9.502019+10,-1.858492+10+0.2,-9.502019+10+0.2,10,0),6,2,byrow = TRUE) #89 points
opt.loc.6

opt.w.6 <- rep(1/6,6)
opt.w.6

Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,opt.loc.6,opt.w.6)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,opt.loc.6,opt.w.6) 
Phi.optopt #0.0005264361


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#54

length(which(Phi.checkopt ==0))/N
#0.02076125

poi <- grid[!Phi.checkopt , ]
#poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill6.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.6[,1],opt.loc.6[,2],pch = 17 ,cex = 2, col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Nloptr-Phi2/nofulfill6.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=1,pch = 1 ,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(opt.loc.6[,1],opt.loc.6[,2],pch = 17 ,cex = 2, col="red")

dev.off()