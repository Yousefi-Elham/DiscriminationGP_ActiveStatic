#Distance based

library(mgcv)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)
tau2 <- 0.01


#NRUN <- 100
#nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 51)
x2 <- seq(0,10,length.out = 51)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

xstart <- matrix(c(10,0,9.2,2,8,0.8),3,2,byrow = TRUE)
xstart

n1 <- nrow(xstart)
wstart <- rep(1/n1,n1)
wstart
max.iter <- 1000

num.max <- CRIT.vec <- numeric(max.iter)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
dir.deriv <- function(theta,sigma2,des,weight){
  
  n <- nrow(des)
  
  Phi.crit <- numeric(N)
  for(i in 1:N){
    
    DisVec <- sqrt( (grid[i,1]-des[,1])^2 + (grid[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)^2*weight) 
    
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
  
  opt.val <- sum(abs(kernel1 - kernel2)^2*(weight%*%t(weight)))
  
  return(opt.val)
}
#-------------------------------------------------------------------------------------------------------
CRIT.vec0 <- Phi.optim(est.theta,sigma2.ini,xstart,wstart)
CRIT.vec0
#-------------------------------------------------------------------------------------------------------

for( k in 1: max.iter){
  
  print(k)
  
  dir.vec <- dir.deriv(est.theta,sigma2.ini,xstart,wstart)
  num.max[k] <- which.max(dir.vec)
  x.hat <- c(grid[num.max[k],1],grid[num.max[k],2])
  
  xstart <- rbind(xstart,x.hat)
  alpha <- 1/(k+1)
  wstart <- c((1-alpha)*wstart,alpha)
  design <- cbind(xstart,wstart)
  
  CRIT.vec[k] <- Phi.optim(est.theta,sigma2.ini,xstart,wstart)
  
}


#-----------------------------------------------------------------------------------------------

nrow(xstart)
nrow(xstart)-n1 


##########################################################

length(CRIT.vec)
CRIT.vec
xstart
design
#####################################################################################################
nrow(xstart)
design1.uniq <- uniquecombs(xstart)
counter1 <- attr(design1.uniq,"index")
counter2 <- unique(sort(counter1))
opt.weight <- numeric(nrow(design1.uniq))

for(z in 1:length(counter2)){
  
  x.replic <- which(counter1==counter2[z])
  opt.weight[z] <- sum(design[x.replic,3])
  
}

design.opt <- cbind(design1.uniq,opt.weight)
design.opt

nrow(design.opt)
sum(design.opt[,3])
#----------------------------------------------------------------
block <- numeric(nrow(design.opt))
for(i in 1:nrow(design.opt)){
  if(design.opt[i,2]>=0.0&design.opt[i,2]<=0.41&design.opt[i,1]>=9.6&design.opt[i,1]<=10.01) block[i] <- 2
  if(design.opt[i,2]>=0.4&design.opt[i,2]<=0.81&design.opt[i,1]>=7.8&design.opt[i,1]<=8.21) block[i] <- 1
  if(design.opt[i,2]>=1.8&design.opt[i,2]<=2.21) block[i] <- 3
}
block

disc.meas <- numeric(3)
for(j in 1:3){
  disc.meas[j] <- sum(design.opt[which(block==j),3])
}
disc.meas

disc.des <- matrix(c(8.2,0.4,10,0.4,9.0,2),3,2,byrow = TRUE)
disc.des
#--------------------------------

Phi.critopt <-  dir.deriv(est.theta,sigma2.ini,disc.des,disc.meas)
Phi.optopt  <-  Phi.optim(est.theta,sigma2.ini,disc.des,disc.meas) 
Phi.optopt 


Phi.checkopt  <- numeric(N)
for(i in 1:N){
  if(Phi.critopt [i] <= Phi.optopt ){
    Phi.checkopt [i] <- 1
  }else{
    Phi.checkopt [i] <- 0
  }
  
}

length(which(Phi.checkopt ==0))
#84

length(which(Phi.checkopt ==0))/N
#0.03229527

poi <- grid[!Phi.checkopt , ]
poi
#----------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/desmeasure.jpeg")

plot(disc.des, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=20*disc.meas,pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
text(disc.des[,1],disc.des[,2]+0.7,round(disc.meas,digits=3),cex = 1.2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/desmeasure.eps")

plot(disc.des, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=20*disc.meas,pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
text(disc.des[,1],disc.des[,2]+0.7,round(disc.meas,digits=3),cex = 1.2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#CRIT.vec[max.iter]

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#10

length(which(Phi.check==0))/N
#0.003844675

poi <- grid[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/nofulfill.jpeg")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/nofulfill.eps")

poi <- grid[!Phi.check, ]
plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------
#restricting optimal points

mm=which(design.opt[,3]>0.01)
mm
design.opt11 <- design.opt[mm,]
design.opt11
nrow(design.opt11)

Phi.crit11 <- dir.deriv(est.theta,sigma2.ini,design.opt11[,c(1,2)],design.opt11[,3])
Phi.opt11 <-  Phi.optim(est.theta,sigma2.ini,design.opt11[,c(1,2)],design.opt11[,3])  
Phi.opt11

Phi.check11 <- numeric(N)
for(i in 1:N){
  if(Phi.crit11[i] <= Phi.opt11){
    Phi.check11[i] <- 1
  }else{
    Phi.check11[i] <- 0
  }
  
}

length(which(Phi.check11==0))
# 22

length(which(Phi.check11==0))/N
#

poi11 <- grid[!Phi.check11, ]
poi11

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/nofulfillrest.jpeg")

plot(poi11, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt11[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#----------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/nofulfillrest.eps")

plot(poi11, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt11[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#####################################################################################################
CRIT.vec1 <- c(CRIT.vec0,CRIT.vec)

max(CRIT.vec1)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/criterionval.jpeg")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.0001)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/criterionval.eps")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.0001)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()
