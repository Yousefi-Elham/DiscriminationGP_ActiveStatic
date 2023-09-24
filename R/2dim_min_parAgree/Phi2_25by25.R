# For agreement on parameter setting to use: 
# using the grid of points with 25*25=625 design locations, we want to find the parameter
# in the second covariance kernel (Matern 5/2) which minimize the \Phi_2 criterion
# when the parameter in the first kernel is fixed at \theta_1=1


sigma2.ini=1
rho.ini1=1
rho.ini2=1
x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)
#---------------------------------------------------------------------------------------------

Phi.func <- function(param,fixrho,sigma2,des){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*fixrho) * exp(-sqrt(3)*DisM*fixrho)
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*param + 5*DisM^2*param^2/3 ) * exp(-sqrt(5)*DisM*param)
  
  
  abs.differ <- abs(kernel1 - kernel2)^2
  Phi2 <- sum(abs.differ)
  return(Phi2)
}

max.log.lik1 <- optim(rho.ini2, Phi.func, fixrho=rho.ini1 ,sigma2=sigma2.ini, des=grid, method = "Brent",
                      lower=0,upper=4,control = list(maxit = 5000))




est.rho <- max.log.lik1$par
est.rho


