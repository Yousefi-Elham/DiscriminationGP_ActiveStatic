sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)
#-------------------------------------------------------------------------------------------------
Phi.optim <- function(theta, sigma2, des, weight) {
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  
  opt.val <- sum(abs(kernel1 - kernel2)^10*(weight%*%t(weight)))
  
  return(opt.val)
}

#-------------------------------------------------------------------------------------------------
findesPhi10 <- matrix(c(9.999607, 0.00297305, 0.3333333,
                        8.159711, 0.56013476, 0.3333333,
                        9.560681, 1.87689219, 0.3333333),3,3,byrow = TRUE)


Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,findesPhi10[,c(1,2)],findesPhi10[,3])
Phi.opt1
#------------------------------------------------------------
findesPhi10.trl.loc <- findesPhi10[,c(1,2)]-matrix(rep(c(4,-4)),3,2,byrow = TRUE)
findesPhi10.trl.loc

findesPhi10.trl <- cbind(findesPhi10.trl.loc,rep(1/3,3))
findesPhi10.trl


Phi.opt2 <-  Phi.optim(est.theta,sigma2.ini,findesPhi10.trl.loc,rep(1/3,3))
Phi.opt2
#------------------------------------------------------------
dist(findesPhi10[,c(1,2)])
#1        2
#2 1.922406         
#3 1.924637 1.922646

poi <- c()
#------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg501/Phi10des3nlop/Phi10_EX1_transl.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(findesPhi10.trl[,c(1,2)], pch = 19, col = "red", cex = 12*sqrt(findesPhi10.trl[,3]))
points(findesPhi10.trl[,c(1,2)], pch = 16, cex = 1)
polygon(findesPhi10.trl[,1], findesPhi10.trl[,2], border = "blue")


dev.off()

#-------------------------------------#-------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi10/dir.alg501/Phi10des3nlop/Phi10_EX1_transl.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(findesPhi10.trl[,c(1,2)], pch = 19, col = "red", cex = 12*sqrt(findesPhi10.trl[,3]))
points(findesPhi10.trl[,c(1,2)], pch = 16, cex = 1)
polygon(findesPhi10.trl[,1], findesPhi10.trl[,2], border = "blue")

dev.off()

