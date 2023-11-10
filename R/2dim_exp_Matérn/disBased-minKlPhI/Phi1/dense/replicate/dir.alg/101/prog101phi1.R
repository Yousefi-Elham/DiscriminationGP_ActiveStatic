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
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
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
#[1] 108
sum(design.opt[,3])

plot(1:10,1:10, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",type = "n",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

#----------------------------------------------------------------
a <- which(design.opt[,1]>=8.59&design.opt[,1]<=10.01&design.opt[,2]>=0.0&design.opt[,2]<=0.69)
b <- which(design.opt[,1]>=7.3&design.opt[,1]<=8.39&design.opt[,2]>=0.5&design.opt[,2]<=2.5)
cc <- which(design.opt[,1]>=8.5&design.opt[,1]<=10.01&design.opt[,2]>=1.6&design.opt[,2]<=2.9)

sort(c(a,b,cc))
length(sort(c(a,b,cc)))

length(a)
#[1] 33

length(b)
#[1] 37

length(cc)
#[1] 38

#all.equal(1:nrow(design.opt),sort(c(a,b,cc)))

sum(design.opt[a,3])
#0.3330003
sum(design.opt[b,3])
#0.3339993
sum(design.opt[cc,3])
#0.3330003

loc3 <- matrix(c(9.6,0.08,7.76, 0.70,9.42, 2.04),3,2,byrow = TRUE)
ww3 <- c(0.3330003,0.3339994,0.3330003)
sum(ww3)

plot(1:10,1:10, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",type = "n",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

#--------------------------------#--------------------------------#--------------------------------
#----------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/desmeasure.jpeg")

plot(loc3, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=15*sqrt(ww3),pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/desmeasure.eps")

plot(loc3, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex=15*sqrt(ww3),pch = 16 ,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#--------------------------------#--------------------------------
#--------------------------------#--------------------------------#--------------------------------
design.opt[,c(1,2)]
design.opt[,3]
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,design.opt[,c(1,2)],design.opt[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#0.02320711

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#50

length(which(Phi.check==0))/N
#0.00490148

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2,pch=16, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2,pch=16, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

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

design.opt[c(39,52,54,60,65,7,23,50,64,51,10,29,14,72,9,61,76,21,18),]

#opt.weight
#x.hat 8.0 1.4 0.01698302
#x.hat 9.4 0.0 0.01698302
#x.hat 9.1 0.3 0.01698302
#x.hat 9.7 0.2 0.01698302
#x.hat 9.1 0.1 0.01698302
#x.hat 9.8 2.0 0.01798202
#x.hat 7.7 1.7 0.01798202
#x.hat 7.9 1.7 0.01798202
#x.hat 9.5 2.0 0.01898102
#x.hat 9.8 2.3 0.01998002
#x.hat 9.3 2.2 0.01998002
#x.hat 7.9 1.1 0.01998002
#x.hat 8.0 1.5 0.02097902
#x.hat 8.0 1.3 0.02097902
#x.hat 9.6 0.1 0.02097902
#x.hat 9.7 2.4 0.02097902
#x.hat 9.7 0.4 0.02197802
#x.hat 9.2 0.4 0.02297702
#x.hat 9.3 2.4 0.02497502

des18phi1 <- matrix(c(8.0, 1.4, 0.0555001,
                      9.1, 0.3, 0.0555001,
                      9.7, 0.2, 0.0555001,
                      9.1, 0.1, 0.0555001,
                      9.8, 2.0, 0.05566655,
                      7.7, 1.7, 0.0555001,
                      7.9, 1.7, 0.0555001,
                      9.5, 2.0, 0.05566655,
                      9.8, 2.3, 0.05566655,
                      9.3, 2.2, 0.055666552,
                      7.9, 1.1, 0.0555001,
                      8.0, 1.5, 0.0555001,
                      8.0, 1.3, 0.0555001,
                      9.6, 0.1, 0.0555001,
                      9.7, 2.4, 0.05566655,
                      9.7, 0.4, 0.0555001,
                      9.2, 0.4, 0.0555001,
                      9.3, 2.4, 0.05566654),18,3,byrow = TRUE)
sum(des18phi1[,3])

des18phi1[,c(1,2)]
des18phi1[,3]
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des18phi1[,c(1,2)],des18phi1[,3])
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des18phi1[,c(1,2)],des18phi1[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt1
#0.02257146

#efficiency of this 18 point design w.r.t the rich design (with 108 locations)
#0.02257146/0.02320711
#[1] 0.9726097

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#169

length(which(Phi.check==0))/N
#0.016567

poi <- gridp[!Phi.check, ]
#poi

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des18phi1[,c(1,2)], pch = 16, col = "red", cex = 2)

#-------------------------------
design.opt[c(9,10,14,18,21,29,61,72,76),]

des9phi1 <- matrix(c(9.6, 0.1, 0.1110002,
                     9.3, 2.2, 0.1110001,
                     8.0, 1.5, 0.1113331,
                     9.3, 2.4, 0.1110001,
                     9.2, 0.4, 0.1110001,
                     7.9, 1.1, 0.1113331,
                     9.7, 2.4, 0.1110001,
                     8.0, 1.3, 0.1113331,
                     9.7, 0.4, 0.1110001),9,3,byrow = TRUE)
sum(des9phi1[,3])
#--------------------------------#--------------------------------#--------------------------------

#------------------------------------
design.opt[a,]
sum(design.opt[a,3])

#------------------------------------
design.opt[b,]
sum(design.opt[b,3])

#bl<- design.opt[b[c(11,12,14)],c(1,2)]

#bw1 <- sum(design.opt[b[c(1,5,7,11,13)],3])
#bw2 <- sum(design.opt[b[c(2,4,6,9,12,15)],3])
#bw3 <- sum(design.opt[b[c(3,8,10,14)],3])

#bd <- cbind(bl,c(bw1,bw2,bw3))
#bd
#------------------------------------
design.opt[cc,]
sum(design.opt[cc,3])

#cl<- design.opt[cc[c(11,12,13)],c(1,2)]

#cw1 <- sum(design.opt[cc[c(3,7,10,11)],3])
#cw2 <- sum(design.opt[cc[c(2,5,8,9,12)],3])
#cw3 <- sum(design.opt[cc[c(1,4,6,13)],3])

#cd <- cbind(cl,c(cw1,cw2,cw3))
#cd


#des9 <- rbind(ad,bd,cd)
#des9

#sum(des9[,3])
#--------------------------------#--------------------------------#--------------------------------
sub9 <- sort(design.opt[,3])[(nrow(design.opt)-20+1):nrow(design.opt)]
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
which(design.opt[,3]==sub9[19])
which(design.opt[,3]==sub9[20])

indx <- c(39,52,54,60,65,7,23,50,64,51,10,29,14,72,9,61,76,21,18)
length(indx)
design.opt[indx,]

#x.hat 8.0 1.4 0.01698302
#x.hat 9.4 0.0 0.01698302
#x.hat 9.1 0.3 0.01698302
#x.hat 9.7 0.2 0.01698302
#x.hat 9.1 0.1 0.01698302
#x.hat 9.8 2.0 0.01798202
#x.hat 7.7 1.7 0.01798202
#x.hat 7.9 1.7 0.01798202
#x.hat 9.5 2.0 0.01898102
#x.hat 9.8 2.3 0.01998002
#x.hat 9.3 2.2 0.01998002
#x.hat 7.9 1.1 0.01998002
#x.hat 8.0 1.5 0.02097902
#x.hat 8.0 1.3 0.02097902
#x.hat 9.6 0.1 0.02097902
#x.hat 9.7 2.4 0.02097902
#x.hat 9.7 0.4 0.02197802
#x.hat 9.2 0.4 0.02297702
#x.hat 9.3 2.4 0.02497502

nrow(design.opt[indx,])
#[1] 19

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx,c(1,2)],col="red")

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/19optloc.jpeg")

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx,c(1,2)],col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/19optloc.eps")

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx,c(1,2)],col="red")

dev.off()
#--------------------------------

des19 <- design.opt[indx,c(1,2)]
sum(design.opt[indx,3])
eqw19 <- rep(1/length(indx),length(indx))
eqw19
sum(eqw19)
length(eqw19)
nrow(des19)
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des19,eqw19)
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des19,eqw19) #equal to CRIT.vec[max.iter] 
Phi.opt1
#[1] 0.02266815

#efficiency of the 19 point design with respect to the best found (with 108 points)
#0.02266815/0.02320711
#[1] 0.9767761

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#165

length(which(Phi.check==0))/N
#0.01617488

poi <- gridp[!Phi.check, ]
#poi


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill19.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des19, pch = 16, col = "red", cex = 1)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill19.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des19, pch = 16, col = "red", cex = 1)

dev.off()
#--------------------------------#--------------------------------
sub9 <- sort(design.opt[,3])[(nrow(design.opt)-50+1):nrow(design.opt)]
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
which(design.opt[,3]==sub9[19])
which(design.opt[,3]==sub9[20])
which(design.opt[,3]==sub9[21])
which(design.opt[,3]==sub9[22])
which(design.opt[,3]==sub9[23])
which(design.opt[,3]==sub9[24])
which(design.opt[,3]==sub9[25])
which(design.opt[,3]==sub9[26])
which(design.opt[,3]==sub9[27])
which(design.opt[,3]==sub9[28])
which(design.opt[,3]==sub9[29])
which(design.opt[,3]==sub9[30])
which(design.opt[,3]==sub9[31])
which(design.opt[,3]==sub9[32])
which(design.opt[,3]==sub9[33])
which(design.opt[,3]==sub9[34])
which(design.opt[,3]==sub9[35])
which(design.opt[,3]==sub9[36])
which(design.opt[,3]==sub9[37])
which(design.opt[,3]==sub9[38])
which(design.opt[,3]==sub9[39])
which(design.opt[,3]==sub9[40])
which(design.opt[,3]==sub9[41])
which(design.opt[,3]==sub9[42])
which(design.opt[,3]==sub9[43])
which(design.opt[,3]==sub9[44])
which(design.opt[,3]==sub9[45])
which(design.opt[,3]==sub9[46])
which(design.opt[,3]==sub9[47])
which(design.opt[,3]==sub9[48])
which(design.opt[,3]==sub9[49])
which(design.opt[,3]==sub9[50])

indx21 <- c(16,42,56,73,32,45,53,75,4,15,98,46,47,36,68,28,48,3,67,20,31,55,11,85,17,71,2,27,77,79)
indx22 <- c(43,39,52,54,60,65,7,23,50,64,51,10,29,14,72,9,61,76,21,18)
indx2 <- c(indx21,indx22)
length(indx2)
design.opt[indx2,]

design.opt[indx2,]
opt.weight
#x.hat 7.7 1.0 0.00999001
#x.hat 7.7 1.1 0.00999001
#x.hat 7.8 1.1 0.00999001
#x.hat 7.8 1.0 0.00999001
#x.hat 8.0 1.6 0.00999001
#x.hat 9.1 0.0 0.00999001
#x.hat 7.6 1.3 0.00999001
#x.hat 9.9 2.2 0.01098901
#x.hat 9.6 2.4 0.01098901
#x.hat 9.6 1.9 0.01098901
#x.hat 9.4 2.1 0.01198801
#x.hat 9.3 0.4 0.01198801
#x.hat 9.9 2.0 0.01198801
#x.hat 9.7 1.9 0.01198801
#x.hat 7.6 1.1 0.01298701
#x.hat 9.9 2.1 0.01298701
#x.hat 9.5 2.5 0.01298701
#      8.0 1.2 0.01332001
#x.hat 9.3 2.3 0.01398601
#x.hat 9.7 0.3 0.01398601
#x.hat 9.5 0.0 0.01398601
#x.hat 9.5 0.5 0.01398601
#x.hat 9.4 0.5 0.01498501
#x.hat 7.6 1.2 0.01498501
#x.hat 9.3 0.0 0.01498501
#x.hat 7.6 1.5 0.01498501
#      9.4 2.0 0.01531802
#x.hat 9.2 0.0 0.01598402
#x.hat 7.8 1.7 0.01598402
#x.hat 9.1 0.2 0.01598402
#x.hat 7.6 1.6 0.01698302
#x.hat 8.0 1.4 0.01698302
#x.hat 9.4 0.0 0.01698302
#x.hat 9.1 0.3 0.01698302
#x.hat 9.7 0.2 0.01698302
#x.hat 9.1 0.1 0.01698302
#x.hat 9.8 2.0 0.01798202
#x.hat 7.7 1.7 0.01798202
#x.hat 7.9 1.7 0.01798202
#x.hat 9.5 2.0 0.01898102
#x.hat 9.8 2.3 0.01998002
#x.hat 9.3 2.2 0.01998002
#x.hat 7.9 1.1 0.01998002
#x.hat 8.0 1.5 0.02097902
#x.hat 8.0 1.3 0.02097902
#x.hat 9.6 0.1 0.02097902
#x.hat 9.7 2.4 0.02097902
#x.hat 9.7 0.4 0.02197802
#x.hat 9.2 0.4 0.02297702
#x.hat 9.3 2.4 0.02497502

nrow(design.opt[indx2,])
#[1] 50

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx2,c(1,2)],col="red")

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/50optloc.jpeg")

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx2,c(1,2)],col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/50optloc.eps")

plot(0:10,0:10,type = "n", xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[indx2,c(1,2)],col="red")

dev.off()
#--------------------------------

des50 <- design.opt[indx2,c(1,2)]
sum(design.opt[indx2,3])
eqw50 <- rep(1/length(indx2),length(indx2))
eqw50
sum(eqw50)
length(eqw50)
nrow(des50)
Phi.crit1 <- dir.deriv(est.theta,sigma2.ini,des50,eqw50)
Phi.opt1 <-  Phi.optim(est.theta,sigma2.ini,des50,eqw50) #equal to CRIT.vec[max.iter] 
Phi.opt1
#[1] 0.02311141

#efficiency of the 50 point design with respect to the best found (with 108 points)
#0.02311141/0.02320711
#[1] 0.9958763

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#108

length(which(Phi.check==0))/N
#0.0105872

poi <- gridp[!Phi.check, ]
nrow(poi)
#poi


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill50.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des50, pch = 16, col = "red", cex = 1)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfill50.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des50, pch = 16, col = "red", cex = 1)

dev.off()

#####################################################################################################
CRIT.vec1 <- c(CRIT.vec0,CRIT.vec)

max(CRIT.vec1)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/criterionval.jpeg")

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
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/criterionval.eps")

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
#0.02320711

Phi.check <- numeric(N)
for(i in 1:N){
  if(Phi.crit1[i] <= Phi.opt1){
    Phi.check[i] <- 1
  }else{
    Phi.check[i] <- 0
  }
  
}

length(which(Phi.check==0))
#12

length(which(Phi.check==0))/N
#0.00461361

poi <- gridp[!Phi.check, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfillG51.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",pch = 16,cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi1/dense/replicate/dir.alg/101/nofulfillG51.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",pch = 16,cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#--------------------------------#--------------------------------#--------------------------------