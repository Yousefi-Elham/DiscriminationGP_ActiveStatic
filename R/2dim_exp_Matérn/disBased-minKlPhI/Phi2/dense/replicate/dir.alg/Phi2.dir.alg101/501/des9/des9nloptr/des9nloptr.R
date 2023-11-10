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


des9Phi2 <- matrix(c(10.00, 0.00, 0.07126207,
                 9.64, 0.04, 0.13186813,
                 9.98, 0.36, 0.12387612,
                 7.86, 0.92, 0.12820513,
                 7.98, 0.46, 0.13086913,
                 8.22, 0.78, 0.07792208,
                 9.48, 2.28, 0.07292707,
                 9.62, 1.90, 0.12787213,
                 9.18, 2.04, 0.13519814),9,3,byrow = TRUE)

des9Phi2
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
  
  opt.val <- -(sum(abs(kernel1 - kernel2)^2*(weight%*%t(weight))))
  
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
#local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
#                    "xtol_rel" = 1.0e-7 )
#opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
#              "xtol_rel" = 1.0e-7,
#              "maxeval" = 1000,
#              "local_opts" = local_opts )



#the following algorithm results into equal weights, I have checked the results from both.

opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
              "xtol_rel" = 1.0e-7,
              "maxeval" = 1000 )
#-------------------------------------------------------------------------------------
res <- nloptr( x0=des9Phi2,eval_f=Phi.optim0 ,eval_grad_f=eval_grad_f0,lb = rep(0,27),
               ub=c(rep(10,18),rep(1,9)), eval_g_eq =eval_g_eq0 , opts=opts)


print( res )


des9nlop <- matrix(res$solution,9,3)


des9nlop

#[,1]         [,2]       [,3]
#[1,] 9.997721 4.686826e-13 0.06350604
#[2,] 9.604483 6.770008e-02 0.13502750
#[3,] 9.933476 3.891591e-01 0.13373990
#[4,] 7.828677 9.693908e-01 0.12750910
#[5,] 7.936548 5.056054e-01 0.12750910
#[6,] 8.177987 8.061984e-01 0.08043498
#[7,] 9.455165 2.332699e+00 0.06350604
#[8,] 9.569251 1.955132e+00 0.13373990
#[9,] 9.132166 2.098411e+00 0.13502750

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

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des9nlop[,c(1,2)],des9nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#0.0005857937

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#0

length(which(Phi.check==0))/N
#0

poi <- gridp[!Phi.check, ]
poi
#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/nofulfillG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 8*sqrt(des9nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/nofulfillG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 8*sqrt(des9nlop[,3]))

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#plot the translated plot in the middle

findesPhi2.trl.loc <- des9nlop[,c(1,2)]-matrix(rep(c(4,-4)),9,2,byrow = TRUE)
findesPhi2.trl.loc

findesPhi2.trl <- cbind(findesPhi2.trl.loc,des9nlop[,3])
findesPhi2.trl


Phi.opt2 <-  Phi.optim(est.theta,sigma2.ini,findesPhi2.trl[,c(1,2)],findesPhi2.trl[,3])
Phi.opt2
#------------------------------------------------------------
dist(findesPhi2.trl[,c(1,2)])

#------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/Phi2_EX1_transl.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(findesPhi2.trl[,c(1,2)], pch = 19, col = "red", cex = 8*sqrt(findesPhi2.trl[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/Phi2_EX1_transl.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(findesPhi2.trl[,c(1,2)], pch = 19, col = "red", cex = 8*sqrt(findesPhi2.trl[,3]))

dev.off()


#-------------------------------------------------------------------------------------------------
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
#0

length(which(Phi.check==0))/N
#0

poi <- gridp[!Phi.check, ]
poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/des9nloptr/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9nlop[,c(1,2)], pch = 19, col = "red", cex = 20*des9nlop[,3])

dev.off()