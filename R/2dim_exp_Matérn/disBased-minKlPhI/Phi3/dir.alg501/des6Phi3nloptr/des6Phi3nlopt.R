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



des6Phi3 <- matrix(c(9.76, 0.00, 0.1671662,
                     10.00, 0.34, 0.1708292,
                     8.02, 0.82, 0.1991342,
                     8.08, 0.38, 0.1428571,
                     9.24, 2.08, 0.1841492,
                     9.66, 1.96, 0.1358641),6,3,byrow = TRUE)

des6Phi3
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

  x1 <- c(x11,x12,x13,x14,x15,x16)
  
  x21 <- des[7]
  x22 <- des[8]
  x23 <- des[9]
  x24 <- des[10]
  x25 <- des[11]
  x26 <- des[12]

  x2 <- c(x21,x22,x23,x24,x25,x26)
  
  w1 <- des[13]
  w2 <- des[14]
  w3 <- des[15]
  w4 <- des[16]
  w5 <- des[17]
  w6 <- des[18]

  
  des.loc <- cbind(x1,x2)
  weight <- c(w1,w2,w3,w4,w5,w6)
  
  n <- nrow(des.loc)
  DisM <- as.matrix(dist(des.loc))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  opt.val <- -(sum(abs(kernel1 - kernel2)^3*(weight%*%t(weight))))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------

eval_grad_f0 <- function( des1 ) {
  jacob=jacobian(func=Phi.optim0, x=des1, method="Richardson")
  jacob
}
#-------------------------------------------------------------------------------------

eval_g_eq0 <- function( des ) {
  #des[1:6] <- 0
  #des[7:12] <- 0
  constr <- c( des[13]+des[14]+des[15]+des[16]+des[17]+des[18]-1 )
  grad <- c(rep(0,12),rep(1,6))
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
res <- nloptr( x0=des6Phi3,eval_f=Phi.optim0 ,eval_grad_f=eval_grad_f0,lb = rep(0,18),
               ub=c(rep(10,12),rep(1,6)), eval_g_eq =eval_g_eq0 , opts=opts)


print( res )


des6nlop <- matrix(res$solution,6,3)


des6nlop

#[,1]         [,2]       [,3]
#[1,] 9.686995 0.01971544 0.1666667
#[2,] 9.999666 0.32671712 0.1666667
#[3,] 7.932820 0.90084236 0.1666667
#[4,] 8.042356 0.47656076 0.1666667
#[5,] 9.150779 2.09671861 0.1666667
#[6,] 9.572986 1.97943868 0.1666667

sum(des6nlop[,3])
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
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)^3*weight) 
    
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
  
  opt.val <- sum(abs(kernel1 - kernel2)^3*(weight%*%t(weight)))
  
  return(opt.val)
}

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des6nlop[,c(1,2)],des6nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des6nlop[,c(1,2)],des6nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.543191e-05

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

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 20*sqrt(des6nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 20*sqrt(des6nlop[,3]))

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------#-------------------------------------
# sacing on smaller scales

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG501s.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 15*sqrt(des6nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG501s.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 15*sqrt(des6nlop[,3]))

dev.off()
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des6nlop[,c(1,2)],des6nlop[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des6nlop[,c(1,2)],des6nlop[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.543191e-05

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

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 20*sqrt(des6nlop[,3]))

dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/des6Phi3nloptr/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6nlop[,c(1,2)], pch = 19, col = "red", cex = 20*sqrt(des6nlop[,3]))

dev.off()