rm(list=ls())
#-----------------------------------------------------------
lenth.g <- 6
lenth.g

grid <- seq(0,lenth.g,length.out = 501)
N <- length(grid)
sigma2=1

D <- seq(0,lenth.g,length.out = 501)

rho1 <- 1


kernel0 <- sigma2*exp(-D*rho1)
kernel1 <- sigma2*(1+sqrt(3)*D*rho1) * exp(-sqrt(3)*D*rho1)
kernel2 <- sigma2*(1 + sqrt(5)*D*rho1 + 5*D^2*rho1^2/3) * exp(-sqrt(5)*D*rho1)
kernel3 <- sigma2*exp(-D^2*rho1^2/2)

plot(seq(0,lenth.g,length=N),kernel0,lty=1,type = "l",xlab="distance",ylab=expression(paste("covariance"," ,", theta,"=1")),ylim=c(0,sigma2),
     cex.main=2, cex.lab=1.3, cex.axis=1.3,lwd=2,col="blue");

lines(seq(0,lenth.g,length=N),kernel1,lty=2,lwd=2,col="red")
lines(seq(0,lenth.g,length=N),kernel2,lty=3,lwd=2,col="black")
lines(seq(0,lenth.g,length=N),kernel3,lty=4,lwd=2,col="magenta")

legend(4.2,0.8,legend=c(expression(paste(nu, "=", "1/2")),expression(paste(nu, "=", "3/2")),
                      expression(paste(nu, "=", "5/2")),expression(nu %->% infinity)),lty=c(1,2,3,4),lwd=c(2,2,2,2),col=c("blue","red","black","magenta"),bty="n",cex=1.5) 


##----------------------------------------------------------------------------------------------
gridp <- seq(0,10,length.out = 501)
THETA.FIX <- c(1, 1.07)
sigma2.ini=1



#generate.y.kernels <- function(theta,sigma2,des,model){
theta <- THETA.FIX[1]
sigma2 <- sigma2.ini
des <- gridp

n <- length(des)
DisM <- matrix(0,n,n)
for(i in 1:n){
  for(j in 1:n){
    DisM[i,j] <- abs(des[i]-des[j]) 
  }
}

kernel0 <- sigma2*exp(-DisM*theta)
kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta) * exp(-sqrt(3)*DisM*theta)
kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta + 5*DisM^2*theta^2/3 ) * exp(-sqrt(5)*DisM*theta)
kernel3 <- sigma2*exp(-DisM^2*theta^2/2)


plot(gridp,rmvnorm(1,sigma=kernel0),type = "l",col="blue",lty=1,lwd=2,xlab = "x",ylab=expression(paste("Y(x)"," ,", nu,"=1/2")),cex.main=2, cex.lab=1.2, cex.axis=1.2)
plot(gridp,rmvnorm(1,sigma=kernel1),type = "l",col="red",lty=2,lwd=2,xlab = "x",ylab=expression(paste("Y(x)"," ,", nu,"=3/2")),cex.main=2, cex.lab=1.2, cex.axis=1.2)
plot(gridp,rmvnorm(1,sigma=kernel2),type = "l",col="black",lty=3,lwd=2,xlab = "x",ylab=expression(paste("Y(x)"," ,", nu,"=5/2")),cex.main=2, cex.lab=1.2, cex.axis=1.2)
plot(gridp,rmvnorm(1,sigma=kernel3),type = "l",col="magenta",lty=4,lwd=2,xlab = "x",ylab=expression(paste("Y(x)"," ,", nu %->% infinity)),cex.main=2, cex.lab=1.2, cex.axis=1.2)

