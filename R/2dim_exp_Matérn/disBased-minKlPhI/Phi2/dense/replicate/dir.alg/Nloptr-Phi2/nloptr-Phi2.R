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
dstart <- c(9.2,2,0.3,0.35)
dstart

sigma2 <- sigma2.ini
theta <- est.theta

#-------------------------------------------------------------------------------------
Phi.optim <- function(dfull) {
  
  x0 <- 10
  y0 <- 0
  x1 <- dfull[1]
  y1 <- dfull[2]
  x2 <- -y1 + 10 
  y2 <- -x1 + 10
  
  w0 <- dfull[3]
  w1 <- dfull[4]
  w2 = w1
  
  des <- matrix(c(x0,y0,x1,y1,x2,y2),nrow=3,byrow = TRUE)
  weight <- c(w0,w1,w2)
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  opt.val <- -(sum(abs(kernel1 - kernel2)^2*(weight%*%t(weight))))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------

eval_grad_f0 <- function( dfull1 ) {
  jacob=jacobian(func=Phi.optim, x=dfull1, method="Richardson")
  jacob
}
#-------------------------------------------------------------------------------------

eval_g_eq0 <- function( dfull ) {
  constr <- c( dfull[3]+2*dfull[4]-1 )
  grad <- c( 0,0,1,2 )
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
#-------------------------------------------------------------------------------------

eval_g_ineq0 <- function( dfull ) {
  constr <- rbind(10-dfull[1]-dfull[2] ,dfull[3]- dfull[4] )
  grad <- rbind( c(-1,-1,0,0) ,c(0,0,1,-1) )
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
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
res <- nloptr( x0=dstart,eval_f=Phi.optim ,eval_grad_f=eval_grad_f0,lb = c(0,0,0,0),
               ub=c(10,10,1,1), eval_g_ineq =eval_g_ineq0 , eval_g_eq =eval_g_eq0 , opts=opts)


print( res )




