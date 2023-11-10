library(numDeriv)

sigma2.ini=1
theta1=1
theta2=1.07

x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)
N <- nrow(grid)

set.seed(123)
actin=sample(1:N,4)
#grid[actin,]

n1 <- length(actin)
xstart <- grid[actin,]
xstart

des <- xstart
sigma2 <- sigma2.ini

n <- nrow(des)
DisM <- as.matrix(dist(des))

DisM


nu <- 5/2
theta <- theta2


kernel <- sigma2*(1 + sqrt(5)*DisM*theta + 5*DisM^2*theta^2/3 ) * exp(-sqrt(5)*DisM*theta)
kernel

h <-  sqrt(2*nu)*DisM*theta #h in your paper, for the   2nd par
partial.h <- sqrt(2*nu)*DisM

h
partial.h


kernel.nu <- sigma2*2^(1-nu)/gamma(nu)*h^(nu)*besselK(h,nu) #compare with kernel2
diag(kernel.nu) <- 1
kernel.nu

all.equal(kernel,kernel.nu)

partial.theta <- -sigma2*2^(1-nu)/gamma(nu)*h^nu*besselK(h,nu-1)*partial.h
partial.theta
diag(partial.theta) <- 0 #??  not sure 


bessel.fun <- function(ka){
  B <- besselK(sqrt(2*ka)*DisM*theta,ka)
  return(B)
}
bessel.par <- jacobian(func=bessel.fun, x=nu, method="Richardson")
bessel.par.mat <- matrix(bessel.par,n,n)
diag(bessel.par.mat) <- 0
isSymmetric(bessel.par.mat)

logh.half <- log(h/2)
diag(logh.half) <- 0 #??  not sure 


partial.nu <- sigma2*2^(1-nu)/gamma(nu)*h^nu *( ((logh.half + 1) - digamma(nu))*besselK(h,nu) + bessel.par.mat )
partial.nu
diag(partial.nu) <- 0


kernelnu.inv <- solve(kernel.nu)

kernel.inv_th <- kernelnu.inv%*%partial.theta
kernel.inv_nu <- kernelnu.inv%*%partial.nu


FIM <- matrix(nrow = 2, ncol = 2)
FIM[1,1] <- 0.5 * sum(t(kernel.inv_th)*kernel.inv_th) # double derivative for theta
FIM[2,2] <- 0.5 * sum(t(kernel.inv_nu)*kernel.inv_nu) # double derivative for nu
FIM[1,2] <- FIM[2,1] <- 0.5 * sum(t(kernel.inv_th)*kernel.inv_nu) # cross derivative between theta and nu


FIM
FIMinv <- solve(FIM) # is the FIM invertible?

#FIM[1,1] refers to the information of the nuisance parameter

ds.crit <- det(FIM)/FIM[1,1] # this should be maximized w.r.t the design
ds.crit
