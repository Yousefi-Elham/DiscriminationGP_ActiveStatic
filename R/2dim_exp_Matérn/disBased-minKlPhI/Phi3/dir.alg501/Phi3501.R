#Distance based

library(mgcv)
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)

#NRUN <- 100
#nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 501)
x2 <- seq(0,10,length.out = 501)
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
#-------------------------------------------------------------------------------------------------------
CRIT.vec0 <- Phi.optim(est.theta,sigma2.ini,xstart,wstart)
CRIT.vec0
#-------------------------------------------------------------------------------------------------------
start.time <- Sys.time()


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


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

nrow(xstart)
nrow(xstart)-n1 


##########################################################

length(CRIT.vec)
CRIT.vec
#xstart
#design
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

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/desloc.jpeg")

plot(0:10,0:10,type="n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(design.opt[,3]))

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/desloc.eps")

plot(0:10,0:10,type="n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(design.opt[,3]))

dev.off()

#----------------------------------------------------------------
#block <- numeric(nrow(design.opt))
#for(i in 1:nrow(design.opt)){
#  if(design.opt[i,2]>=0.0&design.opt[i,2]<=0.51&design.opt[i,1]>=9.61&design.opt[i,1]<=10.01) block[i] <- 2
#  if(design.opt[i,2]>=0.3&design.opt[i,2]<=0.81&design.opt[i,1]>=7.7&design.opt[i,1]<=8.11) block[i] <- 1
#  if(design.opt[i,2]>=1.8&design.opt[i,2]<=2.41) block[i] <- 3
#}
#block

#--------------------------------#--------------------------------#--------------------------------
design.opt[,c(1,2)]
design.opt[,3]
nrow(design.opt)
sum(design.opt[,3])
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.537331e-05

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
# 14

length(which(Phi.check==0))/N
#5.577667e-05

poi <- grid[!Phi.check, ]
poi


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/nofulfill.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(design.opt[,3]))

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/nofulfill.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(design.opt[,3]))

dev.off()
#--------------------------------#--------------------------------#--------------------------------
a <- which(design.opt[,1]>=9.39&design.opt[,1]<=10.09&design.opt[,2]>=0.0&design.opt[,2]<=0.69)
b <- which(design.opt[,1]>=7.7&design.opt[,1]<=8.39&design.opt[,2]>=0.3&design.opt[,2]<=0.99)
cc <- which(design.opt[,1]>=9.09&design.opt[,1]<=9.99&design.opt[,2]>=1.5&design.opt[,2]<=2.5)

sort(c(a,b,cc))
length(sort(c(a,b,cc)))
all.equal(1:40,sort(c(a,b,cc)))
#------------------------------------
design.opt[a,]
nrow(design.opt[a,]) #11
sum(design.opt[a,3]) ##0.3379953

al<- design.opt[a[c(9,11)],c(1,2)]

aw1 <- sum(design.opt[a[c(1,4,5,6,7,9,10)],3])
aw2 <- sum(design.opt[a[c(2,3,8,11)],3])

ad <- cbind(al,c(aw1,aw2))
ad
#------------------------------------
design.opt[b,]
nrow(design.opt[b,]) # 14
sum(design.opt[b,3]) ##0.3419913

bl<- design.opt[b[c(12,14)],c(1,2)]

bw1 <- sum(design.opt[b[c(1,3,5,6,9,12)],3])
bw2 <- sum(design.opt[b[-c(1,3,5,6,9,12)],3])

bd <- cbind(bl,c(bw1,bw2))
bd
#------------------------------------
design.opt[cc,]
nrow(design.opt[cc,])
sum(design.opt[cc,3]) ##0.3200133

cl<- design.opt[cc[c(13,14)],c(1,2)]

cw1 <- sum(design.opt[cc[c(1,3,4,6,7,10,13)],3])
cw2 <- sum(design.opt[cc[-c(1,3,4,6,7,10,13)],3])

cd <- cbind(cl,c(cw1,cw2))
cd


des6 <- rbind(ad,bd,cd)
des6

sum(des6[,3])

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/desloc6.jpeg")

plot(0:10,0:10,type="n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des6[,3]))

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/desloc6.eps")

plot(0:10,0:10,type="n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des6[,3]))

dev.off()
#--------------------------------#--------------------------------#--------------------------------
des6[,c(1,2)]
des6[,3]
nrow(des6)
sum(des6[,3])
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des6[,c(1,2)],des6[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des6[,c(1,2)],des6[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#1.538758e-05 (increased a little compared to the whole des40)

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
# 15

length(which(Phi.check==0))/N
#5.976072e-05

poi <- grid[!Phi.check, ]
poi
#Var1 Var2
#9922  8.04 0.38
#9923  8.06 0.38
#9924  8.08 0.38
#9925  8.10 0.38
#10424 8.06 0.40
#10425 8.08 0.40
#10426 8.10 0.40
#20942 8.00 0.82
#20943 8.02 0.82
#49081 9.66 1.94
#49581 9.64 1.96
#49582 9.66 1.96
#50083 9.66 1.98
#50084 9.68 1.98
#50584 9.66 2.00

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/nofulfill6.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des6[,3]))

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/nofulfill6.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des6[,c(1,2)], pch = 1, col = "red", cex = 20*sqrt(des6[,3]))

dev.off()
#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------

#####################################################################################################
CRIT.vec1 <- c(CRIT.vec0,CRIT.vec)

max(CRIT.vec1)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/criterionval.jpeg")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.00001)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi3/dir.alg501/criterionval.eps")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.00001)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()

