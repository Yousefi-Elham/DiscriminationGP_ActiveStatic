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


# first location is included in the function
# starting values: just 2nd location + weights of first and second locations
wstart <- des9[,3]
wstart

sigma2 <- sigma2.ini
theta <- est.theta

#-------------------------------------------------------------------------------------
Phi.optim0 <- function(w) {
  
  w1 <- w[1]
  w2 <- w[2]
  w3 <- w[3]
  w4 <- w[4]
  w5 <- w[5]
  w6 <- w[6]
  w7 <- w[7]
  w8 <- w[8]
  w9 <- w[9]
  
  des <- des9[,c(1,2)]
  weight <- w
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  opt.val <- -(sum(abs(kernel1 - kernel2)^2*(weight%*%t(weight))))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------

eval_grad_f0 <- function( w1 ) {
  jacob=jacobian(func=Phi.optim0, x=w1, method="Richardson")
  jacob
}
#-------------------------------------------------------------------------------------

eval_g_eq0 <- function( w ) {
  constr <- c( w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7]+w[8]+w[9]-1 )
  grad <- c( 1,1,1,1,1,1,1,1,1 )
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
res <- nloptr( x0=wstart,eval_f=Phi.optim0 ,eval_grad_f=eval_grad_f0,lb = c(0,0,0,0,0,0,0,0,0),
               ub=c(1,1,1,1,1,1,1,1,1), eval_g_eq =eval_g_eq0 , opts=opts)


print( res )




