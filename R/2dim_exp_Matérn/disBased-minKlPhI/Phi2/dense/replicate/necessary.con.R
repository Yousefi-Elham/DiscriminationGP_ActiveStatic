rm(list=ls())

sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)
tau2 <- 0.01

x1 <- seq(0,10,length.out = 51)
x2 <- seq(0,10,length.out = 51)
grid <- expand.grid(x1,x2)
N <- nrow(grid)


phi.optim <- function(theta,sigma2,des){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  #L1 <- kernel1+ tau2*diag(n)
  #L2 <- kernel2+ tau2*diag(n)
  
  abs.differ <- abs(kernel1 - kernel2)^2*1/n*1/n
  Phi.opt <- sum(abs.differ)
  #----------------------------------
  # Phi1.crit <- Phi1.check <- matrix(NA,n,N)
  Phi.crit <- Phi.check <- numeric(length = N)
  for(i in 1:N){
    
    DisVec <- sqrt( (grid[i,1]-des[,1])^2 + (grid[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)^2*1/n)
    
    if(Phi.crit[i] <= Phi.opt){
      Phi.check[i] <- 1
    }else{
      Phi.check[i] <- 0
    }
    
  }
  Phi.sum <- sum(Phi.check)
  return(list(Phi.opt=Phi.opt,Phi.sum=Phi.sum,Phi.check=Phi.check))
}


#--------------------------------#--------------------------------#--------------------------------

des.opt1 <- matrix(c(10,0,9.2,2,8,0.8),3,2,byrow = TRUE)

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt1)
Phi.opt <- check.ret$Phi.opt #0.0004713091
Phi.sum <- check.ret$Phi.sum #2503
Phi.check <- check.ret$Phi.check

length(Phi.check)
Phi.opt
Phi.sum  # is Phi.sum equal N?
N


length(which(Phi.check==0))
#98

length(which(Phi.check==0))/N
#0.03767782
#--------------------------------
# Points not fulfilling the criterion (points of interest):

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt1, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt1, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
des.opt11 <- matrix(c(10,0,9.2,2,8,0.8,7.5,2.5),4,2,byrow = TRUE)

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt11)
Phi.opt <- check.ret$Phi.opt #0.0004589277
Phi.sum <- check.ret$Phi.sum #2496
Phi.check <- check.ret$Phi.check

length(Phi.check)
Phi.opt
Phi.sum  # is Phi.sum equal N?
N


length(which(Phi.check==0))
#105

length(which(Phi.check==0))/N
#0.04036909
#--------------------------------
# Points not fulfilling the criterion (points of interest):

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill11.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt11, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill11.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt11, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
# A 9 point design

des.opt2 <- matrix(c(10.0,0.0,9.6,1.8,9.2,2.0,10.0,0.4,9.6,2.2,8.0,0.8,8.2,0.4,7.8,0.4,9.6,0.0),9,2,byrow = TRUE)
des.opt2

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt2)
Phi.opt <- check.ret$Phi.opt #0.0005776323
Phi.sum <- check.ret$Phi.sum #2590
Phi.check <- check.ret$Phi.check

length(Phi.check)
Phi.opt
Phi.sum  # is Phi.sum equal N?
N


length(which(Phi.check==0))
#11

length(which(Phi.check==0))/N
#0.004229143
#--------------------------------
# Points not fulfilling the criterion (points of interest):

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill2.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt2, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill2.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt2, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------

#3 point optim + a non optim design

des.opt3 <- matrix(c(10,0,9.2,2,8,0.8,5,5),4,2,byrow = TRUE)
des.opt3

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt3)
Phi.opt <- check.ret$Phi.opt #0.0002653743
Phi.sum <- check.ret$Phi.sum #2372
Phi.check <- check.ret$Phi.check

length(Phi.check)
Phi.opt
Phi.sum  # is Phi.sum equal N?
N


length(which(Phi.check==0))
#229

length(which(Phi.check==0))/N
#0.08804306
#--------------------------------
# Points not fulfilling the criterion (points of interest):

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill3.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt3, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill3.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt3, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------
#3 point optim + a non optim design

des.opt4 <- matrix(c(10,0,9.2,2,8,0.8,0.5,8),4,2,byrow = TRUE)
des.opt4

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt4)
Phi.opt <- check.ret$Phi.opt #0.0002651114
Phi.sum <- check.ret$Phi.sum # 2418
Phi.check <- check.ret$Phi.check

length(Phi.check)
Phi.opt
Phi.sum  # is Phi.sum equal N?
N


length(which(Phi.check==0))
#183

length(which(Phi.check==0))/N
#0.07035755
#--------------------------------
# Points not fulfilling the criterion (points of interest):

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill4.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt3, pch = 17, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/nofulfill4.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt3, pch = 17, col = "red", cex = 2)

dev.off()

