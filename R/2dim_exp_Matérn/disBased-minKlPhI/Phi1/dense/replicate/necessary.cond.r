#rm(list=ls())

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
  L1 <- kernel1
  L2 <- kernel2
  
  abs.differ <- abs(L1 - L2)*1/n*1/n
  Phi1.opt <- sum(abs.differ)
  #----------------------------------
  # Phi1.crit <- Phi1.check <- matrix(NA,n,N)
  Phi1.crit <- Phi1.check <- numeric(length = N)
  for(i in 1:N){
    
    DisVec <- sqrt( (grid[i,1]-des[,1])^2 + (grid[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi1.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)*1/n)
    
    if(Phi1.crit[i] <= Phi1.opt){
      Phi1.check[i] <- 1
    }else{
      Phi1.check[i] <- 0
    }
    
  }
  Phi1.sum <- sum(Phi1.check)
  return(list(Phi1.opt=Phi1.opt,Phi1.sum=Phi1.sum,Phi1.check=Phi1.check))
}

#---------------------------------------------------------------------
des.opt1 <- matrix(c(9.6,0.4,9.4,2,8,1.2),3,2,byrow = TRUE)

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt1)
Phi1.opt <- check.ret$Phi1.opt
Phi1.sum <- check.ret$Phi1.sum
Phi1.check <- check.ret$Phi1.check

length(Phi1.check)
Phi1.opt
Phi1.sum  # is Phi1.sum equal N?
N

# Points not fulfilling the criterion (points of interest):
poi <- grid[!Phi1.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt1, pch = 17, col = "red", cex = 2)

length(which(Phi1.check==0))
#241

length(which(Phi1.check==0))/N
#0.09265667

#---------------------------------------------------------------------
des.opt2 <- matrix(c(9.6,0.4,9.4,2,8,1.2,9.8,0.4,10,0.4,9.2,2.2,9.8,2,7.8,1.2,8.2,1),9,2,byrow = TRUE)
des.opt2

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt2)
Phi1.opt <- check.ret$Phi1.opt
Phi1.sum <- check.ret$Phi1.sum
Phi1.check <- check.ret$Phi1.check

length(Phi1.check)
Phi1.opt
Phi1.sum  # is Phi1.sum equal N?
N

# Points not fulfilling the criterion (points of interest):
poi <- grid[!Phi1.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt2, pch = 17, col = "red", cex = 2)

length(which(Phi1.check==0))
#77

length(which(Phi1.check==0))/N
#0.029604

#---------------------------------------------------------------------
des.opt3 <- matrix(c(9.6,0.4,9.4,2,8,1.2,8,3.5),4,2,byrow = TRUE)
des.opt3

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt3)
Phi1.opt <- check.ret$Phi1.opt
Phi1.sum <- check.ret$Phi1.sum
Phi1.check <- check.ret$Phi1.check

length(Phi1.check)
Phi1.opt
Phi1.sum  # is Phi1.sum equal N?
N

# Points not fulfilling the criterion (points of interest):
poi <- grid[!Phi1.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt3, pch = 17, col = "red", cex = 2)

length(which(Phi1.check==0))
#239

length(which(Phi1.check==0))/N
#0.09188774

#---------------------------------------------------------------------
des.opt4 <- matrix(c(9.6,0.4,9.4,2,8,1.2,2,8),4,2,byrow = TRUE)
des.opt4

check.ret <- phi.optim(est.theta,sigma2.ini,des.opt4)
Phi1.opt <- check.ret$Phi1.opt
Phi1.sum <- check.ret$Phi1.sum
Phi1.check <- check.ret$Phi1.check

length(Phi1.check)
Phi1.opt
Phi1.sum  # is Phi1.sum equal N?
N

# Points not fulfilling the criterion (points of interest):
poi <- grid[!Phi1.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des.opt4, pch = 17, col = "red", cex = 2)

length(which(Phi1.check==0))
#405

length(which(Phi1.check==0))/N
#0.1557093