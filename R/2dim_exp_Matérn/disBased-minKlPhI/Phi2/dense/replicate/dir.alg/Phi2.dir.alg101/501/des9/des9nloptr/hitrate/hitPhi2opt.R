library(mvtnorm)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)

NRUN <- 100
tau2 <- 0.01

xstartw <- matrix(c(9.997721, 4.686826e-13, 0.06350604,
                   9.604483, 6.770008e-02, 0.13502750,
                   9.933476, 3.891591e-01, 0.13373990,
                   7.828677, 9.693908e-01, 0.12750910,
                   7.936548, 5.056054e-01, 0.12750910,
                   8.177987, 8.061984e-01, 0.08043498,
                   9.455165, 2.332699e+00, 0.06350604,
                   9.569251, 1.955132e+00, 0.13373990,
                   9.132166, 2.098411e+00, 0.13502750),9,3,byrow=TRUE)

round(xstartw[,3]*50,digits=2)

xstart <- matrix(c(rep(xstartw[1,c(1,2)],3),rep(xstartw[2,c(1,2)],7),rep(xstartw[3,c(1,2)],7),
         rep(xstartw[4,c(1,2)],6),rep(xstartw[5,c(1,2)],6),rep(xstartw[6,c(1,2)],4),
         rep(xstartw[7,c(1,2)],3),rep(xstartw[8,c(1,2)],7),rep(xstartw[9,c(1,2)],7)),50,2,byrow=TRUE)

xstart
M.COUNT <- matrix(0,1,NRUN)
dim(M.COUNT)

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




#-----------------------------------------------------------
set.seed(123456788)
for(NN in 1:NRUN){
  y.obs.K <- generate.y.kernels(est.theta,sigma2.ini,xstart,0)$y #true pars so fix, only model changes
  M.COUNT[1,NN] <- max.likelihood.fun(est.theta,sigma2.ini,xstart,y.obs.K)
}

M.COUNT.N <- 1-M.COUNT
res0 <- mean(M.COUNT.N)
res0 #0.51
#-----------------------------------------------------------
for(NN in 1:NRUN){
  y.obs.K <- generate.y.kernels(est.theta,sigma2.ini,xstart,1)$y #true pars so fix, only model changes
  M.COUNT[1,NN] <- max.likelihood.fun(est.theta,sigma2.ini,xstart,y.obs.K)
}

res1 <- mean(M.COUNT)
res1 #0.65

ress <- (res0+res1)/2
ress #0.58
