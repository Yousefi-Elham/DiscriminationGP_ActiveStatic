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
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

xstart <- matrix(c(9.6,0.4,9.4,2,8,1.2),3,2,byrow = TRUE)
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
    
    DisVec <- sqrt( (gridp[i,1]-des[,1])^2 + (gridp[i,2]-des[,2])^2 )
    cov.test.opt.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
    cov.test.opt.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
    Phi.crit[i] <- sum(abs(cov.test.opt.1 - cov.test.opt.2)*weight) 
    
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
  #    abs.differ[j,k] <- abs(kernel1[j,k] - kernel2[j,k])*weight[j]*weight[k]
  
  
  #  }
  #}
  
  #opt.val <- sum(abs.differ)
  
  opt.val <- sum(abs(kernel1 - kernel2)*(weight%*%t(weight)))
  
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
  x.hat <- c(gridp[num.max[k],1],gridp[num.max[k],2])
  
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
a <- which(design.opt[,1]>=8.79&design.opt[,1]<=10.01&design.opt[,2]>=0.0&design.opt[,2]<=0.69)
b <- which(design.opt[,1]>=7.3&design.opt[,1]<=8.39&design.opt[,2]>=0.5&design.opt[,2]<=2)
cc <- which(design.opt[,1]>=8.5&design.opt[,1]<=10.01&design.opt[,2]>=1.6&design.opt[,2]<=2.9)

sort(c(a,b,cc))
length(sort(c(a,b,cc)))

length(a)
length(b)
length(cc)
all.equal(1:nrow(design.opt),sort(c(a,b,cc)))

sum(design.opt[a,3])
#0.3339993
sum(design.opt[b,3])
#0.3330003
sum(design.opt[cc,3])
#0.3330003

loc3 <- matrix(c(9.6,0.08,7.76, 0.70,9.42, 2.04),3,2,byrow = TRUE)
ww3 <- c(0.3339993,0.3330003,0.3330003)
#--------------------------------#--------------------------------#--------------------------------
#----------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/desmeasure.jpeg")

plot(loc3, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=20*ww3,pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/desmeasure.eps")

plot(loc3, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=20*ww3,pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
design.opt[,c(1,2)]
design.opt[,3]
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
#1275

length(which(Phi.check==0))/N
#0.005079661

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill.eps")


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/poi.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/poi.eps")


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
sub9 <- sort(design.opt[,3])[(nrow(design.opt)-17):nrow(design.opt)]
sub9
length(sub9)
which(design.opt[,3]==sub9[1])
which(design.opt[,3]==sub9[2])
which(design.opt[,3]==sub9[3])
which(design.opt[,3]==sub9[4])
which(design.opt[,3]==sub9[5])
which(design.opt[,3]==sub9[6])
which(design.opt[,3]==sub9[7])
which(design.opt[,3]==sub9[8])
which(design.opt[,3]==sub9[9])
which(design.opt[,3]==sub9[10])
which(design.opt[,3]==sub9[11])
which(design.opt[,3]==sub9[12])
which(design.opt[,3]==sub9[13])
which(design.opt[,3]==sub9[14])
which(design.opt[,3]==sub9[15])
which(design.opt[,3]==sub9[16])
which(design.opt[,3]==sub9[17])
which(design.opt[,3]==sub9[18])

indx <- c(93, 105, 107, 133, 232, 241, 284, 303, 421, 466, 159, 179, 196, 286,287,393, 184,298, 320, 160, 244)
design.opt[indx,]

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx,c(1,2)])

indx1 <- indx[-c(1,3,19)]

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx1,c(1,2)])

des18 <- design.opt[indx1,]

des18[,c(1,2)]
des18[,3]
eqw <- rep(1/18,18)
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des18[,c(1,2)],eqw)
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des18[,c(1,2)],eqw) #equal to CRIT.vec[max.iter] 
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
#4268

length(which(Phi.check==0))/N
#0.01700392

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill18.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des18[,c(1,2)], pch = 16, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill18.eps")


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des18[,c(1,2)], pch = 16, col = "red", cex = 2)

dev.off()
#--------------------------------
Phi.optim(est.theta,sigma2.ini,des18[,c(1,2)],eqw)
#0.02244555

Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
#0.02321241

effdes18 <- Phi.optim(est.theta,sigma2.ini,des18[,c(1,2)],eqw)/Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
effdes18
#0.9669634
#--------------------------------#--------------------------------#--------------------------------
plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des18[c(11,12,15),c(1,2)])
points(des18[c(5,17,18),c(1,2)])
points(des18[c(3,14,16),c(1,2)])

loc9 <- des18[c(3,5,11,12,14,15,16,17,18),c(1,2)]
w9 <- rep(1/9,9)
des9 <- cbind(loc9,w9)
des9[,c(1,2)]
des9[,3]

Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3]) #equal to CRIT.vec[max.iter] 
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
# 5038

length(which(Phi.check==0))/N
#0.02007163

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill9.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 16, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfill9.eps")


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 16, col = "red", cex = 2)

dev.off()
#--------------------------------
Phi.optim(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3])
#0.02224422

Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
#0.02321241

effdes9 <- Phi.optim(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3])/Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
effdes9
#0.95829
#####################################################################################################
CRIT.vec1 <- c(CRIT.vec0,CRIT.vec)

max(CRIT.vec1)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/criterionval.jpeg")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.01)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/criterionval.eps")

plot(seq(1:length(CRIT.vec1)),seq(0,max(CRIT.vec1),length=length(CRIT.vec1)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,length(CRIT.vec1),by = 100),length(CRIT.vec1))
ytick=seq(0,max(CRIT.vec1),by=0.01)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:length(CRIT.vec1)),CRIT.vec1,col="blue",lwd = 2)

dev.off()

#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------

# just changing the grid (check the condition for the whole design on a sparser grid)

x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

design.opt[,c(1,2)]
design.opt[,3]
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
#51

length(which(Phi.check==0))/N
#0.00499951

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
# just changing the grid (check the condition for the whole design on a sparser grid)

x1 <- seq(0,10,length.out = 51)
x2 <- seq(0,10,length.out = 51)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

design.opt[,c(1,2)]
design.opt[,3]
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
#11

length(which(Phi.check==0))/N
#0.004229143

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfillG51.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/nofulfillG51.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
#start the same alg. with the new grid: supports of design.opt


gridp <- design.opt[,c(1,2)]
N <- nrow(gridp)
N

xstart <- des9[,c(1,2)]
xstart

#x.hat 9.02 2.14
#x.hat 9.60 0.08
#x.hat 7.76 0.70
#x.hat 7.90 1.20
#x.hat 9.06 2.52
#x.hat 7.52 1.30
#x.hat 9.42 2.04
#x.hat 9.84 0.38
#x.hat 9.60 0.60

n1 <- nrow(xstart)
wstart <- rep(1/n1,n1)
wstart
max.iter <- 1000

num.max <- CRIT.vec <- numeric(max.iter)

CRIT.vec0 <- Phi.optim(est.theta,sigma2.ini,xstart,wstart)
CRIT.vec0
#-------------------------------------------------------------------------------------------------------
start.time <- Sys.time()


for( k in 1: max.iter){
  
  print(k)
  
  dir.vec <- dir.deriv(est.theta,sigma2.ini,xstart,wstart)
  num.max[k] <- which.max(dir.vec)
  x.hat <- c(gridp[num.max[k],1],gridp[num.max[k],2])
  
  xstart <- rbind(xstart,x.hat)
  alpha <- 1/(k+1)
  wstart <- c((1-alpha)*wstart,alpha)
  design <- cbind(xstart,wstart)
  
  CRIT.vec[k] <- Phi.optim(est.theta,sigma2.ini,xstart,wstart)
  
}


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

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

#----------------------------------

x1 <- seq(0,10,length.out = 501)
x2 <- seq(0,10,length.out = 501)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

design.opt[,c(1,2)]
design.opt[,3]
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#0.02321264

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
# 1141

length(which(Phi.check==0))/N
#0.004545799

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/prognewG/nofulfillNG.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/prognewG/nofulfillNG.eps")


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#----------------------------------

x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
gridp <- expand.grid(x1,x2)
N <- nrow(gridp)
N

nrow(design.opt[,c(1,2)])
sum(design.opt[,3])
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#0.02321264

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#45

length(which(Phi.check==0))/N
#0.004411332

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/prognewG/nofulfillNG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/501/prognewG/nofulfillNG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()