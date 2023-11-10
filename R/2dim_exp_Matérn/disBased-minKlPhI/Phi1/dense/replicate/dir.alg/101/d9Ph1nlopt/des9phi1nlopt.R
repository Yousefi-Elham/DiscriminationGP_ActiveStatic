#Distance based

library(mgcv)
library(numDeriv)
library(nloptr)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)


# first location is included in the function
# starting values: just 2nd location + weights of first and second locations

des9phi1 <- matrix(c(9.6, 0.1, 0.1110002,
                     9.3, 2.2, 0.1110001,
                     8.0, 1.5, 0.1113331,
                     9.3, 2.4, 0.1110001,
                     9.2, 0.4, 0.1110001,
                     7.9, 1.1, 0.1113331,
                     9.7, 2.4, 0.1110001,
                     8.0, 1.3, 0.1113331,
                     9.7, 0.4, 0.1110001),9,3,byrow = TRUE)

des9phi1
sigma2 <- sigma2.ini
theta <- est.theta

#-------------------------------------------------------------------------------------
Phi.optim0 <- function(des) {
  
  x11 <- des[1]
  x12 <- des[2]
  x13 <- des[3]
  x14 <- des[4]
  x15 <- des[5]
  x16 <- des[6]
  x17 <- des[7]
  x18 <- des[8]
  x19 <- des[9]
  x1 <- c(x11,x12,x13,x14,x15,x16,x17,x18,x19)
  
  x21 <- des[10]
  x22 <- des[11]
  x23 <- des[12]
  x24 <- des[13]
  x25 <- des[14]
  x26 <- des[15]
  x27 <- des[16]
  x28 <- des[17]
  x29 <- des[18]
  x2 <- c(x21,x22,x23,x24,x25,x26,x27,x28,x29)
  
  w1 <- des[19]
  w2 <- des[20]
  w3 <- des[21]
  w4 <- des[22]
  w5 <- des[23]
  w6 <- des[24]
  w7 <- des[25]
  w8 <- des[26]
  w9 <- des[27]
  
  des.loc <- cbind(x1,x2)
  weight <- c(w1,w2,w3,w4,w5,w6,w7,w8,w9)
  
  n <- nrow(des.loc)
  DisM <- as.matrix(dist(des.loc))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  opt.val <- -(sum(abs(kernel1 - kernel2)*(weight%*%t(weight))))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------

eval_grad_f0 <- function( des1 ) {
  jacob=jacobian(func=Phi.optim0, x=des1, method="Richardson")
  jacob
}
#-------------------------------------------------------------------------------------

eval_g_eq0 <- function( des ) {
  #des[1:9] <- 0
  #des[10:18] <- 0
  constr <- c( des[19]+des[20]+des[21]+des[22]+des[23]+des[24]+des[25]+des[26]+des[27]-1 )
  grad <- c(rep(0,18),rep(1,9))
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-7,
              "maxeval" = 1000,
              "local_opts" = local_opts )



#the following algorithm results into equal weights, I have checked the results from both.

#opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
#              "xtol_rel" = 1.0e-7,
#              "maxeval" = 1000 )
#-------------------------------------------------------------------------------------
res <- nloptr( x0=des9phi1,eval_f=Phi.optim0 ,eval_grad_f=eval_grad_f0,lb = rep(0,27),
               ub=c(rep(10,18),rep(1,9)), eval_g_eq =eval_g_eq0 , opts=opts)


print( res )


des9nlop <- matrix(c(9.401833, 9.398308, 7.734748, 9.329245, 9.446734, 7.749855, 9.76342, 8.054958, 9.8209, 0.139857, 2.114654,
                     1.519873, 2.506829, 0.5355275, 1.027254, 2.273603, 1.283152, 0.3992495, 0.1150836, 0.1031662, 0.1150836,
                     0.1150836, 0.1031662, 0.1150836, 0.1150836, 0.1031662, 0.1150836),9,3)

des9nlop


des9nlop <- matrix(c(9.368572, 9.387661, 7.708963, 9.324074, 9.414124, 7.718722, 9.757445, 8.025477, 9.790891, 0.1188328, 2.108122,
                     1.528665, 2.501075, 0.5136339, 1.032288, 2.258624, 1.286894, 0.3794695, 0.1115604, 0.109568, 0.1118426,
                     0.1116497, 0.1091081, 0.1120211, 0.1122195, 0.1090617, 0.1129689),9,3)


des9nlop

#[,1]      [,2]      [,3]
#[1,] 9.401833 0.1398570 0.1150836
#[2,] 9.398308 2.1146540 0.1031662
#[3,] 7.734748 1.5198730 0.1150836
#[4,] 9.329245 2.5068290 0.1150836
#[5,] 9.446734 0.5355275 0.1031662
#[6,] 7.749855 1.0272540 0.1150836
#[7,] 9.763420 2.2736030 0.1150836
#[8,] 8.054958 1.2831520 0.1031662
#[9,] 9.820900 0.3992495 0.1150836

sum(des9nlop[,3])
#--------------------------------------------------------------------------------------
x1 <- seq(0,10,length.out = 501)
x2 <- seq(0,10,length.out = 501)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N
#--------------------------------------------------------------------------------------
dir.deriv <- function(theta,sigma2,des,weight){
  
  n <- nrow(des)
  
  Phi.crit <- numeric(N)
  for(i in 1:N){
    
    DisVec <- sqrt( (gridp[i,1]-des[,1])^2 + (gridp[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)*weight) 
    
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
  #    abs.differ[j,k] <- abs(kernel1[j,k] - kernel2[j,k])*weight[j]*weight[k]
  
  
  #  }
  #}
  
  #opt.val <- sum(abs.differ)
  
  opt.val <- sum(abs(kernel1 - kernel2)*(weight%*%t(weight)))
  
  return(opt.val)
}

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#CRIT.vec[max.iter]

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#5824

length(which(Phi.check==0))/N
#0.02320309

poi <- gridp[!Phi.check, ]
poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/d9Ph1nlopt/nofulfillG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/d9Ph1nlopt/nofulfillG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#CRIT.vec[max.iter]

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#241

length(which(Phi.check==0))/N
#0.02362513

poi <- gridp[!Phi.check, ]
poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/d9Ph1nlopt/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/d9Ph1nlopt/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()