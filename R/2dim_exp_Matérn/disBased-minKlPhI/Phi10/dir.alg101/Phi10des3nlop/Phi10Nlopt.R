#Distance based

library(mgcv)
library(numDeriv)
library(nloptr)
#rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)




des8Phi10 <- matrix(c(10.0, 0.0, 0.0003330003,
                     9.2, 2.0, 0.0003330003,
                     8.0, 0.8, 0.0003330003,
                     8.1, 0.4, 0.0009990010,
                     10.0, 0.3, 0.0009990010,
                     9.3, 2.0, 0.3326673327,
                     8.1, 0.5, 0.3326673327,
                     10.0, 0.2, 0.3316683317),8,3,byrow = TRUE)

nrow(des8Phi10)

des3Phi10 <- matrix(c(9.3, 2.0, sum(des8Phi10[c(2,6),3]),
                     8.1, 0.5, sum(des8Phi10[c(3,4,7),3]),
                     10.0, 0.2, sum(des8Phi10[c(1,5,8),3])),3,3,byrow = TRUE)
des3Phi10
sigma2 <- sigma2.ini
theta <- est.theta

#-------------------------------------------------------------------------------------
Phi.optim0 <- function(des) {
  
  x11 <- des[1]
  x12 <- des[2]
  x13 <- des[3]
  
  
  x1 <- c(x11,x12,x13)
  
  x21 <- des[4]
  x22 <- des[5]
  x23 <- des[6]
  
  
  x2 <- c(x21,x22,x23)
  
  w1 <- des[7]
  w2 <- des[8]
  w3 <- des[9]
  
  
  
  des.loc <- cbind(x1,x2)
  weight <- c(w1,w2,w3)
  
  n <- nrow(des.loc)
  DisM <- as.matrix(dist(des.loc))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  opt.val <- -(sum(abs(kernel1 - kernel2)^10*(weight%*%t(weight))))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------

eval_grad_f0 <- function( des1 ) {
  jacob=jacobian(func=Phi.optim0, x=des1, method="Richardson")
  jacob
}
#-------------------------------------------------------------------------------------

eval_g_eq0 <- function( des ) {
  #des[1:3] <- 0
  #des[4:6] <- 0
  constr <- c( des[7]+des[8]+des[9]-1 )
  grad <- c(rep(0,6),rep(1,3))
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
res <- nloptr( x0=des3Phi10,eval_f=Phi.optim0 ,eval_grad_f=eval_grad_f0,lb = rep(0,9),
               ub=c(rep(10,6),rep(1,3)), eval_g_eq =eval_g_eq0 , opts=opts)


print( res )


des3nlop <- matrix(res$solution,3,3)


des3nlop

#[,1] [,2]      [,3]
#[1,]  9.3  2.0 0.3330003
#[2,]  8.1  0.5 0.3339993
#[3,] 10.0  0.2 0.3330003

sum(des3nlop[,3])
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
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)^10*weight) 
    
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
  
  opt.val <- sum(abs(kernel1 - kernel2)^10*(weight%*%t(weight)))
  
  return(opt.val)
}

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des3nlop[,c(1,2)],des3nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des3nlop[,c(1,2)],des3nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.831108e-16

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#6

length(which(Phi.check==0))/N
#2.390429e-05

poi <- gridp[!Phi.check, ]
#poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des3nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des3nlop[,3]))

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------#-------------------------------------
# sacing on smaller scales

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501s.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 15*sqrt(des3nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501s.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 15*sqrt(des3nlop[,3]))

dev.off()
#-------------------------------------------------------------------------------------------------
ww3 <- rep(1/3,3)
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des3nlop[,c(1,2)],ww3)
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des3nlop[,c(1,2)],ww3) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.831109e-16

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#1

length(which(Phi.check==0))/N
#3.984048e-06

poi <- gridp[!Phi.check, ]
#poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501w.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(ww3))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG501w.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(ww3))

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des3nlop[,c(1,2)],ww3)
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des3nlop[,c(1,2)],ww3) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.831109e-16

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#1

length(which(Phi.check==0))/N
#9.80296e-05

poi <- gridp[!Phi.check, ]
poi

#-------------------------------------#-------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des3nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg101/Phi10des3nlop/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des3nlop[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des3nlop[,3]))

dev.off()