rm(list=ls())
#-----------------------------------------------------------

grid <- seq(0,1,by=0.01)
N <- length(grid)
sigma2=2

D <- seq(0,1,by=0.01)

Matern.plot <- function(rho){
  
  kernel1 <- sigma2*exp(-D/rho)
  kernel2 <- sigma2*(1+D/rho) * exp(-D/rho)
  
  plot(seq(0,1,length=N),kernel1, type = "l",xlab="Distance",ylab="Kernel",
       main = bquote( paste( rho, ":", .(rho),", ", sigma^2, ":",.(sigma2)) ),pch=2,lwd=2,col="blue");
  
  lines(seq(0,1,length=N),kernel2,pch=2,lwd=2,col="red")
  
  legend("topright",legend=c("Matern:1/2","Matern: 3/2"),lwd=c(2,2),col=c("blue","red"),bty="n") 
  
  
}

Matern.plot(0.2)
Matern.plot(0.6)
Matern.plot(3)