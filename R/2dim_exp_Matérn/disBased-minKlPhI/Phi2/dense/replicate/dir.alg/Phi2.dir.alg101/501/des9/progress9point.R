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

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
design.opt1 <- matrix(c(10.00, 0.00,
                    9.20, 2.00,
                    8.00, 0.80,
                    9.66, 1.86,
                    8.28, 0.42,
                    8.22, 0.86,
                    9.46, 2.24,
                    10.00, 0.40,
                    7.92, 0.52,
                    9.18, 2.00,
                    9.66, 0.04,
                    7.86, 0.90,
                    9.64, 1.88,
                    10.00, 0.36,
                    7.96, 0.46,
                    9.18, 2.02,
                    9.66, 0.02,
                    7.88, 0.90,
                    9.50, 2.24,
                    8.22, 0.80,
                    9.66, 1.92,
                    8.00, 0.44,
                    9.64, 0.04,
                    9.64, 1.90,
                    9.48, 2.26,
                    8.20, 0.80,
                    7.86, 0.92,
                    7.98, 0.46,
                    7.88, 0.92,
                    9.64, 0.06,
                    8.22, 0.78,
                    8.00, 0.46,
                    9.98, 0.36,
                    9.48, 2.28,
                    9.62, 1.90,
                    9.18, 2.04),nrow=36,byrow=TRUE)
design.opt1
nrow(design.opt1)

weight <- c(0.0712620713, 0.0013320013, 0.0003330003, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010,
            0.0019980020, 0.0009990010, 0.0099900100, 0.0019980020, 0.0409590410, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 0.0009990010, 
            0.1288711289, 0.0159840160, 0.0079920080, 0.0019980020, 0.1238761239, 0.1218781219, 0.0009990010, 0.0009990010, 0.0739260739, 0.0039960040, 0.1128871129,
            0.0629370629, 0.1088911089, 0.0919080919) 



weight
length(weight)
design.opt <- cbind(design.opt1,weight)
nrow(design.opt)

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
#15

length(which(Phi.check==0))/N
#5.976072e-05

poi <- grid[!Phi.check, ]
poi


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)
#---------------------------------------------------------------
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

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
# 1

length(which(Phi.check==0))/N
#9.80296e-05

poi <- grid[!Phi.check, ]
poi


plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(design.opt[,c(1,2)], pch = 1, col = "red", cex = 2)
#------------------------------------
a <- which(design.opt[,1]>=9.49&design.opt[,1]<=10.09&design.opt[,2]>=0.0&design.opt[,2]<=0.59)
b <- which(design.opt[,1]>=7.7&design.opt[,1]<=8.39&design.opt[,2]>=0.3&design.opt[,2]<=0.99)
cc <- which(design.opt[,2]>=1.5&design.opt[,2]<=2.5)

sort(c(a,b,cc))
length(sort(c(a,b,cc)))
all.equal(1:36,sort(c(a,b,cc)))
#------------------------------------
design.opt[a,]
sum(design.opt[a,3])

al<- design.opt[a[c(1,6,8)],c(1,2)]

aw1 <- sum(design.opt[a[1],3])
aw2 <- sum(design.opt[a[c(3,5,6,7)],3])
aw3 <- sum(design.opt[a[c(2,4,8)],3])

ad <- cbind(al,c(aw1,aw2,aw3))
ad
#------------------------------------
design.opt[b,]
sum(design.opt[b,3])

bl<- design.opt[b[c(11,12,14)],c(1,2)]

bw1 <- sum(design.opt[b[c(1,5,7,11,13)],3])
bw2 <- sum(design.opt[b[c(2,4,6,9,12,15)],3])
bw3 <- sum(design.opt[b[c(3,8,10,14)],3])

bd <- cbind(bl,c(bw1,bw2,bw3))
bd
#------------------------------------
design.opt[cc,]
sum(design.opt[cc,3])

cl<- design.opt[cc[c(11,12,13)],c(1,2)]

cw1 <- sum(design.opt[cc[c(3,7,10,11)],3])
cw2 <- sum(design.opt[cc[c(2,5,8,9,12)],3])
cw3 <- sum(design.opt[cc[c(1,4,6,13)],3])

cd <- cbind(cl,c(cw1,cw2,cw3))
cd


des9 <- rbind(ad,bd,cd)
des9

sum(des9[,3])
#-------------------------------------

x1 <- seq(0,10,length.out = 501)
x2 <- seq(0,10,length.out = 501)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

Phi.crit9 <- dir.deriv(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3])
Phi.opt9 <-  Phi.optim(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt9
#CRIT.vec[max.iter]

Phi.check9 <- numeric(N)
for(i in 1:N){
  if(Phi.crit9[i] <= Phi.opt9){
    Phi.check9[i] <- 1
  }else{
    Phi.check9[i] <- 0
  }
  
}

length(which(Phi.check9==0))
#24

length(which(Phi.check9==0))/N
#9.561715e-05

poi <- grid[!Phi.check9, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfillG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()


setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfillG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()
#-------------------------------------#-------------------------------------
x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

Phi.crit9 <- dir.deriv(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3])
Phi.opt9 <-  Phi.optim(est.theta,sigma2.ini,des9[,c(1,2)],des9[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt9
#CRIT.vec[max.iter]

Phi.check9 <- numeric(N)
for(i in 1:N){
  if(Phi.crit9[i] <= Phi.opt9){
    Phi.check9[i] <- 1
  }else{
    Phi.check9[i] <- 0
  }
  
}

length(which(Phi.check9==0))
#2

length(which(Phi.check9==0))/N
#0.0001960592

poi <- grid[!Phi.check9, ]
poi

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfillG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()


setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfillG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#-------------------------------------#-------------------------------------
#the weights are resulted from nloptr in R
des9j <- cbind(des9[,c(1,2)],round(c(0.07006983, 0.1315291, 0.1241078, 0.129419, 0.1299574, 0.07911616, 0.07187044, 0.128963, 0.1349673),digits = 4))
des9j
sum(des9j[,3])

x1 <- seq(0,10,length.out = 101)
x2 <- seq(0,10,length.out = 101)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

Phi.crit9 <- dir.deriv(est.theta,sigma2.ini,des9j[,c(1,2)],des9j[,3])
Phi.opt9 <-  Phi.optim(est.theta,sigma2.ini,des9j[,c(1,2)],des9j[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt9
#CRIT.vec[max.iter]

Phi.check9 <- numeric(N)
for(i in 1:N){
  if(Phi.crit9[i] <= Phi.opt9){
    Phi.check9[i] <- 1
  }else{
    Phi.check9[i] <- 0
  }
  
}

length(which(Phi.check9==0))
#0

length(which(Phi.check9==0))/N
#0

poi <- grid[!Phi.check9, ]
poi



jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfilljG101.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()


setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfilljG101.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#-------------------------------------#-------------------------------------
#the weights are resulted from nloptr in R
des9j <- cbind(des9[,c(1,2)],round(c(0.07006983, 0.1315291, 0.1241078, 0.129419, 0.1299574, 0.07911616, 0.07187044, 0.128963, 0.1349673),digits=4))
des9j
sum(des9j[,3])


x1 <- seq(0,10,length.out = 501)
x2 <- seq(0,10,length.out = 501)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
N

Phi.crit9 <- dir.deriv(est.theta,sigma2.ini,des9j[,c(1,2)],des9j[,3])
Phi.opt9 <-  Phi.optim(est.theta,sigma2.ini,des9j[,c(1,2)],des9j[,3]) #equal to CRIT.vec[max.iter] 
Phi.opt9
#CRIT.vec[max.iter]

Phi.check9 <- numeric(N)
for(i in 1:N){
  if(Phi.crit9[i] <= Phi.opt9){
    Phi.check9[i] <- 1
  }else{
    Phi.check9[i] <- 0
  }
  
}

length(which(Phi.check9==0))
#1 or 7 if weights rounded or not

length(which(Phi.check9==0))/N
#3.984048e-06 or 2.788834e-05

poi <- grid[!Phi.check9, ]
poi



jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfilljG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()


setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfilljG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#---------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfill1jG501.jpeg")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()


setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Phi2/dense/replicate/dir.alg/Phi2.dir.alg101/501/des9/nofulfill1jG501.eps")

plot(poi, xlim = c(0,10), ylim = c(0,10),xlab="X1",ylab="X2",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(des9j[,c(1,2)], pch = 1, col = "red", cex = 2)

dev.off()

#-------------------------------------#-------------------------------------
