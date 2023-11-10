# 2dim sequential

library(parallel)
library(mvtnorm)

rm(list=ls())
par(mfrow=c(1,1))
options(nwarnings = 10000)

#--------------------------------------set pars and starting point-------------------------------------
NRUN <- 50 # number of random runs
nexc <- 50 # size of sequential design

x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)
N <- nrow(grid)
THETA.FIX <- c(1, 1.07)

actin.ini=c(1,25,601,625)
grid[actin.ini,]


sigma2.ini=1
ACTIN.MAT <- CRIT.MAT <- matrix(0, nexc, NRUN)
X.list <- list()

n1 <- length(actin.ini)
nleft <- nexc-n1

MIN.EIG.K1H <- MIN.EIG.K2H <-  M.COUNT <- THETA1.EST <- THETA2.EST <- matrix(0,nleft,NRUN)
E11 <- E12 <- E22 <- E21 <- E1.dif <- E2.dif <- SUM.E1E2 <- matrix(0,nleft,NRUN)

#-------------------------------------------------------------------------------------------------------
generate.y.kernels <- function(theta,sigma2,des,model){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta) * exp(-sqrt(3)*DisM*theta)
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta + 5*DisM^2*theta^2/3 ) * exp(-sqrt(5)*DisM*theta)
  
  if(model==0) y <- rmvnorm(1,sigma=kernel1)
  if(model==1) y <- rmvnorm(1,sigma=kernel2)
  
  return(list(y=y,DisM=DisM,kernel1=kernel1,kernel2=kernel2))
}

#-------------------------------------------------------------------------------------------------------
generate.y.kernels.cond <- function(theta,sigma2,des,ypast,model){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  if (model == 0) {
    kernelmat <- sigma2*(1 + sqrt(3)*DisM*theta) * exp(-sqrt(3)*DisM*theta)
  } else {
    kernelmat <- sigma2*(1 + sqrt(5)*DisM*theta + 5*DisM^2*theta^2/3 ) * exp(-sqrt(5)*DisM*theta)
    # NB: thanks to R's recycling capabilities, 1 is added to each element of DisM
  }
  
  cov.past <- kernelmat[1:(n-1),1:(n-1)]
  cov.past.inv <- solve(cov.past)
  #cov.past.new: specific column (row) of a matrix: a vector so this drops the dimensionality (adding drop=F, drop dim)
  cov.past.new <- kernelmat[n, 1:(n-1)]
  sigma2.new <- kernelmat[n, n]
  
  
  mu.cond <- cov.past.new %*% cov.past.inv %*% t(ypast)
  sigma2.cond <- sigma2.new - cov.past.new %*% cov.past.inv %*% cov.past.new
  
  
  if(sigma2.cond <= 0){
    ynew <- mu.cond
  }else{ 
    ynew <- rnorm(1, mu.cond, sqrt(sigma2.cond))
  }
  
  y <- cbind(ypast, ynew)
  
  return(list(y=y,DisM=DisM))
}
#-------------------------------------------------------------------------------------------------------
max.likelihood.fun <- function(theta,sigma2,des,y){
  
  
  theta1=theta[1]
  theta2=theta[2]
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  
  # log likelihood function under kernel1
  log.lik1 <- function(param1,sigma2,y){
    
    kernel1 <- sigma2*(1 + sqrt(3)*DisM*exp(param1)) * exp(-sqrt(3)*DisM*exp(param1))
    L <- -n/2*log(2*pi) - 1/2*log(det(kernel1)) - 1/2*y%*%solve(kernel1)%*%t(y)
    return(L)
  }
  
  
  # log likelihood function under kernel2
  log.lik2 <- function(param2,sigma2,y){
    
    kernel2 <- sigma2*(1 + sqrt(5)*DisM*exp(param2) + 5*DisM^2*exp(param2)^2/3 ) * exp(-sqrt(5)*DisM*exp(param2))
    L <- -n/2*log(2*pi) - 1/2*log(det(kernel2)) - 1/2*y%*%solve(kernel2)%*%t(y)
    return(L)
  }
  
  
  # Calculate the ML under kernel1
  max.log.lik1 <- optim(log(theta1), log.lik1, sigma2=sigma2, y=y, method = "Brent",
                        lower=-5.5,upper=7,control = list(maxit = 5000, fnscale=-1))
  
  # Calculate the ML under kernel2
  max.log.lik2 <- optim(log(theta2), log.lik2, sigma2=sigma2, y=y, method = "Brent",
                        lower=-1.5,upper=7,control = list(maxit = 5000, fnscale=-1))
  
  lik.eval.1 <- max.log.lik1
  lik.eval.2 <- max.log.lik2
  
  
  par.res <- c(exp(lik.eval.1$par), exp(lik.eval.2$par))
  
  kernelhat1 <- sigma2*(1 + sqrt(3)*DisM*par.res[1]) * exp(-sqrt(3)*DisM*par.res[1])
  kernelhat2 <- sigma2*(1 + sqrt(5)*DisM*par.res[2] + 5*DisM^2*par.res[2]^2/3 ) * exp(-sqrt(5)*DisM*par.res[2])
  
  min.eig.k1hat <- (min(eigen(kernelhat1)$values))^(-1)
  min.eig.k2hat <- (min(eigen(kernelhat2)$values))^(-1)
  
  
  if(abs(lik.eval.1$value - lik.eval.2$value) < 1e-6){
    lik.res <- sample(c(0, 1), 1)
  }else{
    lik.res <- as.numeric(lik.eval.2$value > lik.eval.1$value)
  }
  
  # lik.res is  0 or 1 (prediction for which model is the TRUE model)
  
  return(list(par.res=par.res,lik.res=lik.res, min.eig.k1hat=min.eig.k1hat, min.eig.k2hat=min.eig.k2hat))
}

#-------------------------------------------------------------------------------------------------------
kernel.inv.fun <- function(theta, sigma2, des) {
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  if(det(kernel1)==0 | det(kernel2)==0){
    
    kernel1.inv <- kernel2.inv <- NULL
    
    if(det(kernel1)==0 & det(kernel2)==0){
      warning("both kernels are singular")
    }else{
      if( det(kernel1)==0 ){
        warning("first kernel is singular")
      }else{
        warning("second kernel is singular")
      } 
    }
  } else {
    kernel1.inv <- solve(kernel1)
    kernel2.inv <- solve(kernel2)
  }	
  
  return(list(kernel1.inv = kernel1.inv, kernel2.inv = kernel2.inv))
}
#-------------------------------------------------------------------------------------------------------
crit.fun <- function(theta,sigma2,des,ypast,cov.past.inv.1,cov.past.inv.2){
  
  n <- nrow(des)
  
  DisVec <- sqrt( (des[1:(n-1),1]-des[n,1])^2 + (des[1:(n-1),2]-des[n,2])^2 )
  
  cov.past.new.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
  cov.past.new.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
  
  sigma.new.1 <- sigma2
  sigma.new.2 <- sigma2
  
  pred.error.11 <- sigma.new.1 - cov.past.new.1 %*% cov.past.inv.1 %*% cov.past.new.1
  yhat.1 <- cov.past.new.1 %*% cov.past.inv.1 %*% t(ypast)
  
  pred.error.22 <- sigma.new.2 - cov.past.new.2 %*% cov.past.inv.2 %*% cov.past.new.2
  yhat.2 <- cov.past.new.2 %*% cov.past.inv.2 %*% t(ypast)
  
  Prediction.dis <- pred.error.11/pred.error.22 + pred.error.22/pred.error.11 + (yhat.1 - yhat.2)^2 *(1/pred.error.11 + 1/pred.error.22)
  Prediction.dis
}
#-------------------------------------------------------------------------------------------------------
pred.error.fun <- function(theta,sigma2,pastin,leftin){
  
  n1 <- length(pastin)
  nleft <- length(leftin)
  
  DisM <- as.matrix(dist(grid))
  kernel1.grid <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2.grid <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  K1 <- kernel1.grid[pastin,pastin]
  K2 <- kernel2.grid[pastin,pastin]
  K1.inv <- solve(K1)
  K2.inv <- solve(K2)
  
  k1mat <- kernel1.grid[leftin,pastin]
  k2mat <- kernel2.grid[leftin,pastin]
  
  A1 <- crossprod(k1mat) # equivalent to t(k1mat)%*%k1mat, slightly more efficient
  A2 <- crossprod(k2mat)
  A12 <- crossprod(k1mat,k2mat) # equivalent to t(k1mat)%*%k2mat
  # remark: mat1%*%t(mat2) can be computed with tcrossprod(mat1,mat2)
  
  e1.1 <- sigma2 - sum(K1.inv * A1)/nleft
  e1.2 <- sigma2 + (sum(K2.inv%*%K1%*%K2.inv * A2) - 2*sum(K2.inv * A12))/nleft
  e2.2 <- sigma2 - sum(K2.inv * A2)/nleft
  e2.1 <- sigma2 + (sum(K1.inv%*%K2%*%K1.inv * A1) - 2*sum(K1.inv * A12))/nleft
  
  e1.dif <- e1.2 - e1.1
  e2.dif <- e2.1 - e2.2
  
  sum.e1e2 <- e1.2 - e1.1 + e2.1 - e2.2
  
  return(list(e1.1=e1.1,e1.2=e1.2,e2.2=e2.2,e2.1=e2.1,e1.dif=e1.dif,e2.dif=e2.dif,sum.e1e2=sum.e1e2))
  
}

#-------------------------------------------------------------------------------------------------------
sequential.run <- function(NN) {
  
  theta.ini <- c(1, 1.07) # starting pars (also true ones)
  THETA.FIX <- c(1, 1.07) #TRUE pars
  
  MIN.EIG.K1H <- MIN.EIG.K2H <-  M.COUNT <- THETA1.EST <- THETA2.EST <- numeric(length=nleft)
  E11 <- E12 <- E22 <- E21 <- E1.dif <- E2.dif <- SUM.E1E2 <- numeric(length=nleft)
  
  #starting designs
  actin <- actin.ini
  canin <- setdiff(1:N, actin)
  CRIT.vec <- numeric(nexc)
  
  xstart <- grid[actin,]
  xcan <- grid[canin,]
  #xcan <- grid[-actin,] CHECKED
  
  y.obs.K <- generate.y.kernels(THETA.FIX[1],sigma2.ini,xstart,0)$y
  lik.evaluation <- max.likelihood.fun(theta.ini,sigma2.ini,xstart,y.obs.K)
  theta.ini <- lik.evaluation$par.res
  
  for( k in 1: nleft){
    
    nc <- nrow(xstart)
    Ncan <- nrow(xcan)
    
    kernel.inv.list <- kernel.inv.fun(theta.ini,sigma2.ini,xstart)
    kernel1.inv <- kernel.inv.list$kernel1.inv
    kernel2.inv <- kernel.inv.list$kernel2.inv
    
    if (is.null(kernel1.inv) || is.null(kernel2.inv)) {
      warning(paste0("current design singular (run ", NN, ", step ", k, ")"))
      break
    }
    
    ad.crit <- numeric(Ncan)
    for(m in 1: Ncan){
      test.des <- rbind(xstart,xcan[m,])
      ad.crit[m] <- as.numeric(crit.fun(theta.ini,sigma2.ini,test.des,y.obs.K,kernel1.inv,kernel2.inv))
    }
    
    if(max(ad.crit)==0){
      warning("best design chosen has a zero crit. value (run ", NN, ", step ", k, ")")
      break
    }
    
    ad.ind <- which.max(ad.crit)
    grid.max.ind <- which(grid[,1]==xcan[ad.ind,1]&grid[,2]==xcan[ad.ind,2])
    actin <- c(actin,grid.max.ind)
    canin <- setdiff(1:N, actin)
    
    ad.point <- xcan[ad.ind,]
    xhat <- ad.point
    #xstart <- rbind(xstart,ad.point) CHECKED
    #xcan <- xcan[-ad.ind,] CHECKED
    xstart <- grid[actin,]
    xcan <- grid[canin,]
    
    CRIT.vec[k+n1] <- max(ad.crit)
    
    y.obs.K <- generate.y.kernels.cond(THETA.FIX[1],sigma2.ini,xstart,y.obs.K,0)$y #true pars so fix, only model changes
    
    lik.evaluation <- max.likelihood.fun(theta.ini,sigma2.ini,xstart,y.obs.K)
    
    theta.ini <- lik.evaluation$par.res #updated parameters
    THETA1.EST[k] <- theta.ini[1]
    THETA2.EST[k] <- theta.ini[2]
    
    M.COUNT[k] <- dc <- lik.evaluation$lik.res # 0 or 1
    
    MIN.EIG.K1H[k] <- lik.evaluation$min.eig.k1hat
    MIN.EIG.K2H[k] <- lik.evaluation$min.eig.k2hat 
    
    pred.error.list <- pred.error.fun(THETA.FIX,sigma2.ini,actin,canin)
    E11[k] <- pred.error.list$e1.1
    E12[k] <- pred.error.list$e1.2
    E22[k] <- pred.error.list$e2.2
    E21[k] <- pred.error.list$e2.1
    E1.dif[k] <- pred.error.list$e1.dif
    E2.dif[k] <- pred.error.list$e2.dif
    SUM.E1E2[k] <- pred.error.list$sum.e1e2
    
    
  }
  
  return(list(THETA1.EST = THETA1.EST, THETA2.EST = THETA2.EST, M.COUNT = M.COUNT, 
              MIN.EIG.K1H = MIN.EIG.K1H, MIN.EIG.K2H = MIN.EIG.K2H, X = xstart, ACTIN = actin, CRIT = CRIT.vec,
              E11=E11, E12=E12, E22=E22, E21=E21, E1.dif=E1.dif, E2.dif=E2.dif,SUM.E1E2=SUM.E1E2))
  
}
#-------------------------------------------------------------------------------------------------------


detectCores()  # how many cores are available at the local machine?
cl <- makeCluster(6) # create a socket cluster with 6 workers; should be changed to number of cores available
# On each worker, a new empty R session is started (independent from the master process)
clusterEvalQ(cl, library(mvtnorm))  # load package mvtnorm in all workers
clusterExport(cl, list("grid","actin.ini","N","sigma2.ini","THETA.FIX","nexc","n1","nleft", "generate.y.kernels",
                       "generate.y.kernels.cond","max.likelihood.fun","kernel.inv.fun","crit.fun","pred.error.fun"))
# objects that need to be exported to the workers so they are defined there
clusterSetRNGStream(cl, 123456789) # set random number seeds in all workers according to some RNG stream (for reproducibility)

system.time({
  ret <- parLapply(cl, 1:NRUN, sequential.run)
})

stopCluster(cl) # stop worker processes

# example: sapply(ret, function(x) x$THETA1.EST) --> same as do.call below
THETA1.EST <- do.call("cbind", lapply(ret, function(x) x$THETA1.EST)) 
THETA2.EST <- do.call("cbind", lapply(ret, function(x) x$THETA2.EST))
M.COUNT <- do.call("cbind", lapply(ret, function(x) x$M.COUNT))
MIN.EIG.K1H <- do.call("cbind", lapply(ret, function(x) x$MIN.EIG.K1H))
MIN.EIG.K2H <- do.call("cbind", lapply(ret, function(x) x$MIN.EIG.K2H))
X.list <- lapply(ret, function(x) x$X)
ACTIN.MAT <- do.call("cbind", lapply(ret, function(x) x$ACTIN))
CRIT.MAT <- do.call("cbind", lapply(ret, function(x) x$CRIT))
E11 <- do.call("cbind", lapply(ret, function(x) x$E11))
E12 <- do.call("cbind", lapply(ret, function(x) x$E12))
E22 <- do.call("cbind", lapply(ret, function(x) x$E22))
E21 <- do.call("cbind", lapply(ret, function(x) x$E21))
E1.dif <- do.call("cbind", lapply(ret, function(x) x$E1.dif))
E2.dif <- do.call("cbind", lapply(ret, function(x) x$E2.dif))
SUM.E1E2 <- do.call("cbind", lapply(ret, function(x) x$SUM.E1E2))
#-------------------------------------------------------------------------------------------------------
dim(M.COUNT)

###now we should sum up zeros (indication for model 0)
###instead we compute the mean for 1-M.COUNT so that we have 1 to get the mean

M.COUNT.N <- 1-M.COUNT
#-----------------------------
M.COUNT.N[nleft,]

E11[nleft,]
E12[nleft,]
E22[nleft,]
E21[nleft,]
E1.dif[nleft,]
E2.dif[nleft,]
SUM.E1E2[nleft,]

dim(M.COUNT.N)
dim(E11)

hit.iter <- apply(M.COUNT.N,1,mean)
hit.iter

E11.iter <- apply(E11,1,mean)
E12.iter <- apply(E12,1,mean)
E22.iter <- apply(E22,1,mean)
E21.iter <- apply(E21,1,mean)
E1.dif.iter <- apply(E1.dif,1,mean)
E2.dif.iter <- apply(E2.dif,1,mean)
SUM.E1E2.iter <- apply(SUM.E1E2,1,mean)

E11.iter
E12.iter
E22.iter
E21.iter
E1.dif.iter
E2.dif.iter
SUM.E1E2.iter
###########################################################
sel.iter <- c(1,2,3,4,5,6,16,26,36,nleft)

hit.iter[sel.iter]

E11.iter[sel.iter]
E12.iter[sel.iter]
E22.iter[sel.iter]
E21.iter[sel.iter]
E1.dif.iter[sel.iter]
E2.dif.iter[sel.iter]
SUM.E1E2.iter[sel.iter]

round(E11.iter[sel.iter], digits = 4)
round(E12.iter[sel.iter], digits = 4)
round(E22.iter[sel.iter], digits = 4)
round(E21.iter[sel.iter], digits = 4)
round(E1.dif.iter[sel.iter], digits = 4)
round(E2.dif.iter[sel.iter], digits = 4)
round(SUM.E1E2.iter[sel.iter], digits = 4)
###########################################################

dim(THETA1.EST)
dim(THETA2.EST)

THETA1.EST[,1:20]
THETA1.EST[,21:40]
THETA1.EST[,41:50]

THETA2.EST[,1:20]
THETA2.EST[,21:40]
THETA2.EST[,41:50]

mean(THETA1.EST[nleft,])
mean(THETA2.EST[nleft,])

min(THETA1.EST)
min(THETA2.EST)
max(THETA1.EST)
max(THETA2.EST)
max(THETA1.EST,THETA2.EST)


dim(MIN.EIG.K1H)
dim(MIN.EIG.K2H)

min(MIN.EIG.K1H)
min(MIN.EIG.K2H)

max(MIN.EIG.K1H)
max(MIN.EIG.K2H)
#---------------#---------------
av.sd.iter <- apply(M.COUNT.N,1,function(x) sqrt((NRUN-1)/NRUN*var(x)))

av.sd.iter


count2 <- c(1,2,3,4,5,6,16,26,36,46)

av.sd.iter[count2]
round(av.sd.iter[count2],digits = 4)

#---------------#---------------
########################################################################################################
#-----------------------------#-----------------------------#-----------------------------#-----------------------------
### FIRST CHOOSE THE RIGHT ADDRESS!!!!
### FIRST CHOOSE THE RIGHT ADDRESS!!!!
### FIRST CHOOSE THE RIGHT ADDRESS!!!!
### FIRST CHOOSE THE RIGHT ADDRESS!!!!

#-----------------------------#-----------------------------#-----------------------------#--------------
library(ggplot2)
library(mgcv)

des0 <- do.call(rbind.data.frame,X.list[1:3])
des1 <- data.frame(x1=des0[,1],x2=des0[,2],run=rep(1:3,each=50))
des1.unq <- uniquecombs(des1[,-3])
dim(des1.unq)
des1.ind <- attr(des1.unq,"index")
unq.ind <- unique(sort(des1.ind))
des1.runFac <- des1.rep <- numeric(length(unq.ind))
for(s in 1:length(des1.rep)){
  count <- which(des1.ind==unq.ind[s])
  des1.rep[s] <- length(count)
  if(length(count)==1){
    des1.runFac[s] <- des1$run[count]
  }else{
    des1.runFac[s] <- 4
  }
  
}

des2 <- data.frame(x1=des1.unq[,1],x2=des1.unq[,2],replicate=factor(des1.rep),runPlus1=factor(des1.runFac))

#--------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/ranRunRep.jpeg")

gfirst <- ggplot(des2,aes(x=x1, y=x2, color=runPlus1,size=replicate))
gfirst + geom_point()

dev.off()
#------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/ranRunRep.eps")

gfirst <- ggplot(des2,aes(x=x1, y=x2, color=runPlus1,size=replicate))
gfirst + geom_point()

dev.off()
#--------------------------------------------------
des11 <- data.frame(x1=des0[,1],x2=des0[,2],run=factor(rep(1:3,each=50)))

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/ranRun.jpeg")

ggfirst <- ggplot(des11, aes(x=x1, y=x2, color=run))
ggfirst  + geom_point()

dev.off()
#--------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/ranRun.eps")

ggfirst <- ggplot(des11, aes(x=x1, y=x2, color=run))
ggfirst  + geom_point()

dev.off()

#---------------------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/histogram.jpeg")

all.des.b <- do.call(rbind.data.frame,X.list)
all.des <- data.frame(x1=all.des.b[,1],x2=all.des.b[,2])
dim(all.des)

ff <- ggplot(all.des, aes(x=x1,y=x2) )
ff + geom_bin2d(bins=25) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/histogram.eps")

all.des.b <- do.call(rbind.data.frame,X.list)
all.des <- data.frame(x1=all.des.b[,1],x2=all.des.b[,2])
dim(all.des)

ff <- ggplot(all.des, aes(x=x1,y=x2) )
ff + geom_bin2d(bins=25) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


dev.off()

#-----------------------------#-----------------------------#-----------------------------#-----------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/plotparest1.jpeg")

plot(seq(1:(nleft)),seq(0,1,length=nleft),axes = F,ylim=c(0,max(THETA1.EST,THETA2.EST)), type = "n",xlab="Iter",ylab="Parameter estimate",
     main= bquote( paste(theta, ":",.(THETA.FIX[1])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for(NN in 1:NRUN){
  
  points(1:(nleft), THETA1.EST[,NN] ,ylim=c(0,max(THETA1.EST,THETA2.EST)),col=NN, pch=16,cex=1)
  
}
abline(h=THETA.FIX[1],lty=2,col="red") #true par of true model
xtick=seq(1,nleft,by=2)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()


dev.off()
#-----------------------------

#-----------------------------#-----------------------------#-----------------------------#-----------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/plotparest2.jpeg")

plot(seq(1:(nleft)),seq(0,1,length=nleft),axes = F,ylim=c(0,max(THETA1.EST,THETA2.EST)), type = "n",xlab="Iter",ylab="Parameter estimate",
     main= bquote( paste(theta, ":",.(round(THETA.FIX[2],digits=4))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for(NN in 1:NRUN){
  
  points(1:(nleft), THETA2.EST[,NN] ,ylim=c(0,max(THETA1.EST,THETA2.EST)),col=NN, pch=16,cex=1)
  
}
abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=seq(1,nleft,by=2)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
#-----------------------------#-----------------------------#-----------------------------#-----------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest1.jpeg")

boxplot(t(THETA1.EST),axes=F ,ylim=c(0,max(THETA1.EST,THETA2.EST)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[1],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[1],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest1.eps")

boxplot(t(THETA1.EST),axes=F ,ylim=c(0,max(THETA1.EST,THETA2.EST)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[1],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[1],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()

#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest1cut.eps")

boxplot(t(THETA1.EST),axes=F ,ylim=c(0,6), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[1],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[1],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=1)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest2.jpeg")

boxplot(t(THETA2.EST),axes=F ,ylim=c(0,max(THETA1.EST,THETA2.EST)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[2],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest2.eps")

boxplot(t(THETA2.EST),axes=F ,ylim=c(0,max(THETA1.EST,THETA2.EST)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[2],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=5)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxparest2cut.eps")

boxplot(t(THETA2.EST),axes=F ,ylim=c(0,6), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Parameter estimate")
title(bquote( paste(theta, ":",.(round(THETA.FIX[2],digits=4))) ),cex.main=2 )
abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(THETA1.EST,THETA2.EST),by=1)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------#-----------------------------#-----------------------------#-----------------------------
min(MIN.EIG.K1H)
max(MIN.EIG.K1H)

min(MIN.EIG.K2H)
max(MIN.EIG.K2H)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/plotmineig1.jpeg")

plot(seq(1:(nleft)),seq(0,max(MIN.EIG.K1H),length=nleft),axes = F,ylim=c(0,max(MIN.EIG.K1H)), type = "n",xlab="Iter",ylab="Inverse min eigen values of K0",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for(NN in 1:nleft){
  
  points(rep(NN,NRUN), MIN.EIG.K1H[NN,] ,ylim=c(0,max(MIN.EIG.K1H)),col=NN, pch=16,cex=1)
  
}
#abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=seq(1,nleft,by=2)
ytick=seq(0,max(MIN.EIG.K1H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/plotmineig2.jpeg")

plot(seq(1:(nleft)),seq(0,max(MIN.EIG.K2H),length=nleft),axes = F,ylim=c(0,max(MIN.EIG.K2H)), type = "n",xlab="Iter",ylab="Inverse min eigen values of K1",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for(NN in 1:nleft){
  
  points(rep(NN,NRUN), MIN.EIG.K2H[NN,] ,ylim=c(0,max(MIN.EIG.K2H)),col=NN, pch=16,cex=1)
  
}
#abline(h=THETA.FIX[2],lty=2,col="red") #true par of true model
xtick=seq(1,nleft,by=2)
ytick=seq(0,max(MIN.EIG.K2H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------------------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig1.jpeg")

boxplot(t(MIN.EIG.K1H),axes=F ,ylim=c(0,max(MIN.EIG.K1H)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K0")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K1H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig1.eps")

boxplot(t(MIN.EIG.K1H),axes=F ,ylim=c(0,max(MIN.EIG.K1H)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K0")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K1H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig1(cut).eps")

boxplot(t(MIN.EIG.K1H),axes=F ,ylim=c(0,600), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K0")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K1H),by=200)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig2.jpeg")

boxplot(t(MIN.EIG.K2H),axes=F ,ylim=c(0,max(MIN.EIG.K2H)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K1")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K2H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig2.eps")

boxplot(t(MIN.EIG.K2H),axes=F ,ylim=c(0,max(MIN.EIG.K2H)), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K1")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K2H),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxmineig2(cut).eps")

boxplot(t(MIN.EIG.K2H),axes=F ,ylim=c(0,600), cex.main=2, cex.lab=1.5, cex.axis=1.5,
        xlab="Iteration", ylab="Inverse minimum eigen values of K1")
xtick=c(seq(1,nleft,by=2),nleft)
ytick=seq(0,max(MIN.EIG.K2H),by=200)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#--------------------------------- Plot of criterion values (for all random runs)------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/criterionval.jpeg")

plot(seq(1:nexc),seq(0,max(CRIT.MAT),length=nexc), type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for (NN in 1:NRUN){
  lines(seq(1:nexc),CRIT.MAT[,NN],col=NN,cex=1)
}

dev.off()

#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/criterionval.eps")

plot(seq(1:nexc),seq(0,max(CRIT.MAT),length=nexc), type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for (NN in 1:NRUN){
  lines(seq(1:nexc),CRIT.MAT[,NN],col=NN,cex=1)
}

dev.off()

#-----------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/criterionval(cut).jpeg")

plot(seq(1:nexc),seq(0,max(CRIT.MAT),length=nexc),ylim=c(0,600), type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for (NN in 1:NRUN){
  lines(seq(1:nexc),CRIT.MAT[,NN],col=NN,cex=1)
}

dev.off()
#-----------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/criterionval(cut).eps")

plot(seq(1:nexc),seq(0,max(CRIT.MAT),length=nexc),ylim=c(0,600), type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


for (NN in 1:NRUN){
  lines(seq(1:nexc),CRIT.MAT[,NN],col=NN,cex=1)
}

dev.off()
#--------------------------------------------------------------------------------------
#------------------------
max(CRIT.MAT)
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot.jpeg")

boxplot(t(CRIT.MAT),axes=F , cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=1)
ytick=seq(0,max(CRIT.MAT),by=2000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot.eps")

boxplot(t(CRIT.MAT),axes=F , cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=1)
ytick=seq(0,max(CRIT.MAT),by=2000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot(cut).jpeg")

boxplot(t(CRIT.MAT),ylim=c(0,400),axes = F,cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=1)
ytick=seq(0,400,by=100)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot(cut).eps")

boxplot(t(CRIT.MAT),ylim=c(0,400),axes = F,cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=1)
ytick=seq(0,400,by=100)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#----------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot(cut2).eps")

boxplot(t(CRIT.MAT),ylim=c(0,10),axes = F,cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=1)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#---------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/boxplot(cut3).eps")

CRIT.MAT1 <- CRIT.MAT[-5,]
boxplot(t(CRIT.MAT1),ylim=c(0,6),axes = F,cex.main=2, cex.lab=1.5, cex.axis=1.5,xlab="Iteration", ylab="Criterion value")
xtick=seq(1,nexc,by=3)
ytick=seq(0,6,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

dev.off()
#-------------------------------------------------------------------------------------------------------
#--------------------------------------JUST FOR ONE RANDOM RUN------------------------------------------
########################################################################################################
# since lower distances are more important (we are interested to analyse those in detail),
# we exclude the 3 starting points to reduce their effect on the distances
NN <- 1 #one random run
xseq <- c(6,7,8,9,10,20,30,40,nrow(X.list[[NN]]))
xseq
length(xseq)
X.list[[NN]]
n1

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis1.jpeg")

op <- par(mfrow=c(3,3))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
##---------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis1.eps")

op <- par(mfrow=c(3,3))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
########################################################################################################
NN <- 2 #one random run
xseq <- c(6,7,8,9,10,20,30,40,nrow(X.list[[NN]]))
xseq
length(xseq)
X.list[[NN]]
n1

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis2.jpeg")

op <- par(mfrow=c(3,3))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
##---------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis2.eps")

op <- par(mfrow=c(3,3))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
######------------------------------------------
xseq <- c(10,20,30,40)
xseq
length(xseq)
NN <- 2 #one random run
X.list[[NN]]
n1


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis0.jpeg")

op <- par(mfrow=c(1,4))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
##---------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/DesDis0.eps")

op <- par(mfrow=c(1,4))
for(l in 1:length(xseq)){
  
  Trimat.lenth <- choose(xseq[l]-n1,2) #50 is the design size
  Trimat.all <- matrix(0,Trimat.lenth,NN)
  dim(Trimat.all)
  
  DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq[l],]))
  
  Trimat.all[,NN] <- DisM[lower.tri(DisM)]
  
  Vec.alldis <- as.vector(Trimat.all)
  hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
##############################################################################################
xseq <- nexc
xseq
NN <- 1
op <- par(mfrow=c(1,1))

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/onefindesDis1.jpeg")

Trimat.lenth <- choose(xseq-n1,2) #50 is the design size
Trimat.all <- matrix(0,Trimat.lenth,NN)
dim(Trimat.all)
DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq,]))
Trimat.all[,NN] <- DisM[lower.tri(DisM)]

Vec.alldis <- as.vector(Trimat.all)
hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)


dev.off()
#-----------------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/onefindesDis1.eps")

Trimat.lenth <- choose(xseq-n1,2) #50 is the design size
Trimat.all <- matrix(0,Trimat.lenth,NN)
dim(Trimat.all)
DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq,]))
Trimat.all[,NN] <- DisM[lower.tri(DisM)]

Vec.alldis <- as.vector(Trimat.all)
hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)


dev.off()
#########################################################################################################
##############################################################################################
xseq <- nexc
xseq
NN <- 2
op <- par(mfrow=c(1,1))

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/onefindesDis2.jpeg")

Trimat.lenth <- choose(xseq-n1,2) #50 is the design size
Trimat.all <- matrix(0,Trimat.lenth,NN)
dim(Trimat.all)
DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq,]))
Trimat.all[,NN] <- DisM[lower.tri(DisM)]

Vec.alldis <- as.vector(Trimat.all)
hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)


dev.off()
#-----------------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/seq(cond)-minKlPhI-Parll/K1/onefindesDis2.eps")

Trimat.lenth <- choose(xseq-n1,2) #50 is the design size
Trimat.all <- matrix(0,Trimat.lenth,NN)
dim(Trimat.all)
DisM <- as.matrix(dist(X.list[[NN]][(n1+1):xseq,]))
Trimat.all[,NN] <- DisM[lower.tri(DisM)]

Vec.alldis <- as.vector(Trimat.all)
hist(Vec.alldis ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("iter=", .(xseq))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)


dev.off()
#########################################################################################################
