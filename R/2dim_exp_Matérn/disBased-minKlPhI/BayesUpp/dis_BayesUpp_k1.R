#Distance based


library(mvtnorm)
library(dplyr) #anti_join
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)

NRUN <- 100
nexc <- 200 # size of sequential design
x1 <- seq(0,10,length.out = 25)
x2 <- seq(0,10,length.out = 25)
grid <- expand.grid(x1,x2)
N <- nrow(grid)

actin=c(1,25,601,625)
#grid[actin,]

n1 <- length(actin)
xstart <- grid[actin,]
xstart
xcan <- grid[-actin,]
nleft <- nexc-n1

all.equal(as.numeric(row.names(xcan)),setdiff(1:N,as.numeric(row.names(xstart))))

M.COUNT <- matrix(0,nleft,NRUN)
dim(M.COUNT)

E11 <- E12 <- E22 <- E21 <- E1.dif <- E2.dif <- SUM.E1E2 <- numeric(length=nleft)

CRIT.vecU <- CRIT.vecL <- numeric(nexc)
#--------------------------------------------------------------------------------------
generate.y.kernels <- function(theta,sigma2,des,model){
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
  
  if(model==0) y <- rmvnorm(1,sigma=kernel1)
  if(model==1) y <- rmvnorm(1,sigma=kernel2)
  
  return(list(y=y,DisM=DisM,kernel1=kernel1,kernel2=kernel2))
}

#--------------------------------------------------------------------------------------
max.likelihood.fun <- function(theta,sigma2,des,y){
  
  
  theta1=theta[1]
  theta2=theta[2]
  
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta1) * exp(-sqrt(3)*DisM*theta1)
  log.lik1 <- -n/2*log(2*pi) - 1/2*log(det(kernel1)) - 1/2*y%*%solve(kernel1)%*%t(y)
  
  kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta2 + 5*DisM^2*theta2^2/3 ) * exp(-sqrt(5)*DisM*theta2)
  log.lik2 <- -n/2*log(2*pi) - 1/2*log(det(kernel2)) - 1/2*y%*%solve(kernel2)%*%t(y)
  
  
  if(abs(log.lik1 - log.lik2) < 1e-6){
    lik.res <- sample(c(0, 1), 1)
  }else{
    lik.res <- as.numeric(log.lik2 > log.lik1)
  }
  
  # lik.res is  0 or 1 (prediction for which model is the TRUE model)
  return(lik.res)
}
#--------------------------------------------------------------------------------------
crit.fun.old <- function(theta,sigma2,des){
  
  
  n <- nrow(des)
  
  if(n==1){
    BayesL <- 0 
  }else{
    
    DisM <- as.matrix(dist(des))
    
    kernel1 <- sigma2*(1 + sqrt(3)*DisM*theta[1]) * exp(-sqrt(3)*DisM*theta[1])
    kernel2 <- sigma2*(1 + sqrt(5)*DisM*theta[2] + 5*DisM^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisM*theta[2])
    
    
    if(det(kernel1)<=0 | det(kernel2)<=0){
      BayesL <- 0
    }else{
      
      BayesL.p1 <- -1/2*log( 1/2 + 1/2*exp( -1/2*(sum(diag(kernel1 %*% solve(kernel2)))-log(det(kernel1 %*% solve(kernel2))) -n) ) )       
      BayesL.p2 <- -1/2*log( 1/2 + 1/2*exp( -1/2*(sum(diag(kernel2 %*% solve(kernel1)))-log(det(kernel2 %*% solve(kernel1))) -n) ) )
      BayesL <- BayesL.p1 + BayesL.p2
      
      BayesU <- 1/8*(sum(diag(solve(kernel1) %*% kernel2)) + sum(diag(solve(kernel2) %*% kernel1)) -2*n)
      
      
    }
    
  }
  
  return(list(BayesL=BayesL, BayesU=BayesU))
  
}

#-------------------------------------------------------------------------------------------------
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
    kernel1.n <- kernel1
    kernel2.n <- kernel2
    
  }	
  
  return(list(kernel1.n = kernel1.n, kernel2.n = kernel2.n, kernel1.inv = kernel1.inv, kernel2.inv = kernel2.inv))
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
crit.fun <- function(theta,sigma2,des,cov.past.1,cov.past.2,cov.past.inv.1,cov.past.inv.2){
  
  
  n <- nrow(des)
  DisVec <- sqrt( (des[1:(n-1),1]-des[n,1])^2 + (des[1:(n-1),2]-des[n,2])^2 )
  
  cov.past.new.1 <- sigma2*(1 + sqrt(3)*DisVec*theta[1]) * exp(-sqrt(3)*DisVec*theta[1])
  cov.past.new.2 <- sigma2*(1 + sqrt(5)*DisVec*theta[2] + 5*DisVec^2*theta[2]^2/3 ) * exp(-sqrt(5)*DisVec*theta[2])
  
  
  kernel1.n <- kernel1.inv <- matrix(0,n,n)
  kernel1.n[1:(n-1),1:(n-1)] <- cov.past.1
  kernel1.n[n,1:(n-1)] <- kernel1.n[1:(n-1),n] <- cov.past.new.1
  kernel1.n[n,n] <- sigma2
  
  S1inv.s12 <- cov.past.inv.1 %*% cov.past.new.1
  B1 <- drop(sigma2 - cov.past.new.1 %*% S1inv.s12)
  kernel1.inv[1:(n-1),1:(n-1)] <- cov.past.inv.1 + 1/B1 * S1inv.s12 %*% t(S1inv.s12)
  kernel1.inv[n,1:(n-1)] <- kernel1.inv[1:(n-1),n] <- -1/B1*S1inv.s12
  kernel1.inv[n,n] <- 1/B1
  #---------------
  kernel2.n <- kernel2.inv <- matrix(0,n,n)
  kernel2.n[1:(n-1),1:(n-1)] <- cov.past.2
  kernel2.n[n,1:(n-1)] <- kernel2.n[1:(n-1),n] <- cov.past.new.2
  kernel2.n[n,n] <- sigma2
  
  S2inv.s12 <- cov.past.inv.2 %*% cov.past.new.2
  B2 <- drop(sigma2 - cov.past.new.2 %*% S2inv.s12)
  kernel2.inv[1:(n-1),1:(n-1)] <- cov.past.inv.2 + 1/B2* S2inv.s12 %*% t(S2inv.s12)
  kernel2.inv[n,1:(n-1)] <- kernel2.inv[1:(n-1),n] <- -1/B2*S2inv.s12
  kernel2.inv[n,n] <- 1/B2
  
  
  if(det(kernel1.n)<=0 | det(kernel2.n)<=0){
    BayesL <- 0
  }else{
    
    BayesL.p1 <- -1/2*log( 1/2 + 1/2*exp( -1/2*(sum(diag(kernel1.n %*% kernel2.inv))-log(det(kernel1.n %*% kernel2.inv)) -n) ) )       
    BayesL.p2 <- -1/2*log( 1/2 + 1/2*exp( -1/2*(sum(diag(kernel2.n %*% kernel1.inv))-log(det(kernel2.n %*% kernel1.inv)) -n) ) )
    BayesL <- BayesL.p1 + BayesL.p2
    
    BayesU <- 1/8*( sum(diag(kernel1.inv %*% kernel2.n)) + sum(diag(kernel2.inv %*% kernel1.n))-2*n )
    
  }
  
  return(list(BayesL=BayesL, BayesU=BayesU))
  
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
kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart[1:3,])
kernel1.n <- kernel.inv.list$kernel1.n
kernel2.n <- kernel.inv.list$kernel2.n
kernel1.inv <- kernel.inv.list$kernel1.inv
kernel2.inv <- kernel.inv.list$kernel2.inv

crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)$BayesU

crit.fun.old(est.theta,sigma2.ini,xstart)$BayesU

CRIT.vecU[1:n1] <- crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)$BayesU
CRIT.vecU[1:n1]

CRIT.vecL[1:n1] <- crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)$BayesL
CRIT.vecL[1:n1]

#----------------------cunstructing the rest (of design) up to the required size-------------------
system.time({
  
  set.seed(123456789)
  for( k in 1: nleft){
    
    print(k)
    
    nc <- nrow(xstart)
    Ncan <- nrow(xcan)
    
    kernel.inv.list <- kernel.inv.fun(est.theta,sigma2.ini,xstart)
    kernel1.n <- kernel.inv.list$kernel1.n
    kernel2.n <- kernel.inv.list$kernel2.n
    kernel1.inv <- kernel.inv.list$kernel1.inv
    kernel2.inv <- kernel.inv.list$kernel2.inv
    
    
    if (is.null(kernel1.inv) || is.null(kernel2.inv)) {
      warning(paste0("current design singular (step ", k, ")"))
      break
    }
    
    
    ad.crit <- numeric(Ncan)
    for(m in 1: Ncan){
      test.des <- rbind(xstart,xcan[m,])
      ad.crit[m] <- crit.fun(est.theta,sigma2.ini,test.des,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)$BayesU
      #crit.fun.old(est.theta,sigma2.ini,test.des)$BayesU
    }
    
    if(max(ad.crit)==0){
      print(k-1)
      stop("best design chosen has a zero crit. value")
    }
    
    ad.ind <- which.max(ad.crit)
    ad.point <- xcan[ad.ind,]
    xstart <- rbind(xstart,ad.point)
    xcan <- xcan[-ad.ind,]
    
    CRIT.vecU[k+n1] <- max(ad.crit)
    #--------------------------------
    CRIT.val <- max(ad.crit) 
    #CRIT.val <-  crit.fun(est.theta,sigma2.ini,xstart,kernel1.n,kernel2.n,kernel1.inv,kernel2.inv)$BayesU
    
    
    exnum <- 0
    finish.all <- FALSE
    
    while(!finish.all){
      
      nc1 <- nrow(xstart)
      dr.crit <- numeric(nc1)
      
      for(a in 1: nc1){
        dr.crit[a] <- crit.fun.old(est.theta,sigma2.ini,xstart[-a,])$BayesU
        #we can't use the crit.fun here while the dimensions (for past covariances) match but the designs not!
      }
      dr.ind <- which.max(dr.crit)
      dr.point <- xstart[dr.ind,]
      xcan <- rbind(xcan, dr.point)
      xstart.dr <- xstart[-dr.ind,]
      
      
      Ncan1 <- nrow(xcan)
      add.crit <- numeric(Ncan1)
      
      kernel.inv.list.dr <- kernel.inv.fun(est.theta,sigma2.ini,xstart.dr)
      kernel1.n.dr <- kernel.inv.list.dr$kernel1.n
      kernel2.n.dr <- kernel.inv.list.dr$kernel2.n
      kernel1.inv.dr <- kernel.inv.list.dr$kernel1.inv
      kernel2.inv.dr <- kernel.inv.list.dr$kernel2.inv
      
      for(b in 1: Ncan1){
        test.des <- rbind(xstart.dr,xcan[b,])
        add.crit[b] <- crit.fun(est.theta,sigma2.ini,test.des,kernel1.n.dr,kernel2.n.dr,kernel1.inv.dr,kernel2.inv.dr)$BayesU
        ##crit.fun.old(est.theta,sigma2.ini,test.des)$BayesU
        
      }
      add.ind <- which.max(add.crit)
      add.point <- xcan[add.ind,]
      xstart.ad <- rbind(xstart.dr,add.point)
      xcan <- xcan[-add.ind,]
      
      CRIT.exc <- max(add.crit)
      
      
      if( CRIT.exc > CRIT.val & dr.point[,1]!=add.point[,1] & dr.point[,2]!=add.point[,2]){
        xstart <- xstart.ad
        exnum <- exnum+1
        CRIT.val <- CRIT.exc 
        print(rbind(dr.point,add.point))
        #as.numeric(row.names(xstart))
      }else{
        finish.all <- TRUE
        xcan <- grid[setdiff(1:N, as.numeric(row.names(xstart))),]
        #xcan <- anti_join(grid, xstart, by = c("Var1", "Var2"))
      }
      
    }
    
    CRIT.vecU[k+n1] <- crit.fun.old(est.theta,sigma2.ini,xstart)$BayesU
    CRIT.vecL[k+n1] <- crit.fun.old(est.theta,sigma2.ini,xstart)$BayesL
    
    for(NN in 1:NRUN){
      y.obs.K <- generate.y.kernels(est.theta,sigma2.ini,xstart,0)$y #true pars so fix, only model changes
      M.COUNT[k,NN] <- max.likelihood.fun(est.theta,sigma2.ini,xstart,y.obs.K)
    }
    
    
    actin <- as.numeric(row.names(xstart))
    canin <- setdiff(1:N, actin)
    
    pred.error.list <- pred.error.fun(est.theta,sigma2.ini,actin,canin)
    E11[k] <- pred.error.list$e1.1
    E12[k] <- pred.error.list$e1.2
    E22[k] <- pred.error.list$e2.2
    E21[k] <- pred.error.list$e2.1
    E1.dif[k] <- pred.error.list$e1.dif
    E2.dif[k] <- pred.error.list$e2.dif
    SUM.E1E2[k] <- pred.error.list$sum.e1e2
    
  }
  
})

#-----------------------------------------------------------------------------------------------
as.numeric(row.names(xstart))
setdiff(1:N,as.numeric(row.names(xstart)))
s.all <- sort(c(as.numeric(row.names(xstart)),setdiff(1:N,as.numeric(row.names(xstart)))))
s.all2 <- sort(c(as.numeric(row.names(xstart)),as.numeric(row.names(xcan))))
all.equal(1:N,s.all)
all.equal(1:N,s.all2)

nrow(xstart)
nrow(xstart)-n1 
nleft
M.COUNT[nleft,]
M.COUNT <- M.COUNT[1:(nrow(xstart)-n1),]

dim(M.COUNT)
M.COUNT[nrow(xstart)-n1,]


M.COUNT.N <- 1-M.COUNT
M.COUNT.N[nrow(xstart)-n1,]

av.hit.iter <- apply(M.COUNT.N,1,mean)
av.hit.run <- apply(M.COUNT.N,2,mean)

av.hit.iter
av.hit.run

dim(M.COUNT.N)
nrow(M.COUNT.N)
length(av.hit.iter)

E11
E12
E22
E21
E1.dif
E2.dif
SUM.E1E2


int.hit.n <- c(1,2,3,4,5,6,16,36,76,116,156,length(av.hit.iter))
int.hit.n
av.hit.iter[int.hit.n]

E11[int.hit.n]
E12[int.hit.n]
E22[int.hit.n]
E21[int.hit.n]
E1.dif[int.hit.n]
E2.dif[int.hit.n]
SUM.E1E2[int.hit.n]

round(E11[int.hit.n],digits=5)
round(E12[int.hit.n],digits=5)
round(E22[int.hit.n],digits=5)
round(E21[int.hit.n],digits=5)
round(E1.dif[int.hit.n],digits=5)
round(E2.dif[int.hit.n],digits=5)
round(SUM.E1E2[int.hit.n],digits=5)
##########################################################
sel.iter <- c(1,2,3,4,5,6,16,26,36,46)

av.hit.iter[sel.iter]
E11[sel.iter]
E12[sel.iter]
E22[sel.iter]
E21[sel.iter]
E1.dif[sel.iter]
E2.dif[sel.iter]
SUM.E1E2[sel.iter]

round(av.hit.iter[sel.iter],digits=5)
round(E11[sel.iter],digits=5)
round(E12[sel.iter],digits=5)
round(E22[sel.iter],digits=5)
round(E21[sel.iter],digits=5)
round(E1.dif[sel.iter],digits=5)
round(E2.dif[sel.iter],digits=5)
round(SUM.E1E2[sel.iter],digits=5)
##########################################################

length(CRIT.vecU)
CRIT.vecU <- CRIT.vecU[1:nrow(xstart)]
length(CRIT.vecU)
CRIT.vecU

CRIT.vecL <- CRIT.vecL[1:nrow(xstart)]
length(CRIT.vecL)
CRIT.vecL

xstart
#####################################################################################################
#####################################################################################################


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/design.jpeg")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(xstart[,1],xstart[,2], type = "p",pch=16,col= "black");
#points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");


dev.off()
#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/design.eps")

plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
     main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(nrow(xstart))) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);

xtick <- seq(0,10,by=2)
ytick=seq(0,10,by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()
abline(h=seq(0,10,length=21),col="gray")
abline(v=seq(0,10,length=21),col="gray")
points(xstart[,1],xstart[,2], type = "p",pch=16,col= "black");
#points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");


dev.off()
##############################################################################################
#xstart
nrow(xstart)
size.seg <- c(20,40,80,120,160,nrow(xstart))
size.seg
seg <- length(size.seg)
seg


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/des.seg.jpeg")

op <- par(mfrow=c(2,3))
for(kk in 1:seg){
  
  plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
       main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(size.seg[kk])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);
  
  xtick <- seq(0,10,by=2)
  ytick=seq(0,10,by=2)
  axis(1,at=xtick,labels=T,cex.axis=1.5)
  axis(2,at=ytick,labels=T,cex.axis=1.5)
  box()
  abline(h=seq(0,10,length=21),col="gray")
  abline(v=seq(0,10,length=21),col="gray")
  points(xstart[1:size.seg[kk],1],xstart[1:size.seg[kk],2], type = "p",pch=16,col= "black");
  #points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");
  
}
par(op)

dev.off()
#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/des.seg.eps")

op <- par(mfrow=c(2,3))
for(kk in 1:seg){
  
  plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
       main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(size.seg[kk])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);
  
  xtick <- seq(0,10,by=2)
  ytick=seq(0,10,by=2)
  axis(1,at=xtick,labels=T,cex.axis=1.5)
  axis(2,at=ytick,labels=T,cex.axis=1.5)
  box()
  abline(h=seq(0,10,length=21),col="gray")
  abline(v=seq(0,10,length=21),col="gray")
  points(xstart[1:size.seg[kk],1],xstart[1:size.seg[kk],2], type = "p",pch=16,col= "black");
  #points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");
  
}
par(op)

dev.off()

#---------------

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/des.seg1.jpeg", width=650, height=500)

op <- par(mfrow=c(2,3))
for(kk in 1:seg){
  
  plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
       main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(size.seg[kk])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);
  
  xtick <- seq(0,10,by=2)
  ytick=seq(0,10,by=2)
  axis(1,at=xtick,labels=T,cex.axis=1.5)
  axis(2,at=ytick,labels=T,cex.axis=1.5)
  box()
  abline(h=seq(0,10,length=21),col="gray")
  abline(v=seq(0,10,length=21),col="gray")
  points(xstart[1:size.seg[kk],1],xstart[1:size.seg[kk],2], type = "p",pch=16,col= "black");
  #points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");
  
}
par(op)

dev.off()

#----------------------
setEPS()
pdf("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/des.seg1.eps", width = 11, height = 7)

op <- par(mfrow=c(2,3))
for(kk in 1:seg){
  
  plot(seq(0,10,length=20),seq(0,10,length=20),axes = F, type = "n",xlab="X1",ylab="X2",
       main = bquote( paste(sigma^2, ":",.(sigma2.ini), ", ", size, ":",.(size.seg[kk])) ),cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2);
  
  xtick <- seq(0,10,by=2)
  ytick=seq(0,10,by=2)
  axis(1,at=xtick,labels=T,cex.axis=1.5)
  axis(2,at=ytick,labels=T,cex.axis=1.5)
  box()
  abline(h=seq(0,10,length=21),col="gray")
  abline(v=seq(0,10,length=21),col="gray")
  points(xstart[1:size.seg[kk],1],xstart[1:size.seg[kk],2], type = "p",pch=16,col= "black");
  #points(xstart[1:n1,1],xstart[1:n1,2], type = "p",pch=16,col="blue");
  
}
par(op)

dev.off()
##############################################################################################
# since lower distances are more important (we are interested to analyse those in detail),
# we exclude the 3 starting points to reduce their effect on the distances
nrow(xstart)
size.seg <- c(20,40,80,120,160,nrow(xstart))
size.seg
seg <- length(size.seg)
seg


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/DesDis.jpeg")

op <- par(mfrow=c(2,3))
for(l in 1:seg){
  Trimat.length <- choose(size.seg[l],2) #50 is the design size
  Trimat.vec <- numeric(Trimat.length)
  length(Trimat.vec)
  DisM <- as.matrix(dist(xstart[1:size.seg[l],]))
  Trimat.vec <- DisM[lower.tri(DisM)]
  hist(Trimat.vec ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("size=", .(size.seg[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/DesDis.eps")

op <- par(mfrow=c(2,3))
for(l in 1:seg){
  Trimat.length <- choose(size.seg[l],2) #50 is the design size
  Trimat.vec <- numeric(Trimat.length)
  length(Trimat.vec)
  DisM <- as.matrix(dist(xstart[1:size.seg[l],]))
  Trimat.vec <- DisM[lower.tri(DisM)]
  hist(Trimat.vec ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("size=", .(size.seg[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
#--------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/DesDis1.eps", width = 11, height = 7)

op <- par(mfrow=c(2,3))
for(l in 1:seg){
  Trimat.length <- choose(size.seg[l],2) #50 is the design size
  Trimat.vec <- numeric(Trimat.length)
  length(Trimat.vec)
  DisM <- as.matrix(dist(xstart[1:size.seg[l],]))
  Trimat.vec <- DisM[lower.tri(DisM)]
  hist(Trimat.vec ,xlim = c(0,15), xlab = "Distances",main=bquote(paste("size=", .(size.seg[l]))) ,breaks = 10, prob=TRUE,cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
}
par(op)

dev.off()
#####################################################################################################
min(CRIT.vecU)
max(CRIT.vecU)

min(CRIT.vecL)
max(CRIT.vecL)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionval.jpeg")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecU,col="blue",lwd = 2)

dev.off()
#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionval.eps")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecU,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionval2.jpeg")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecL,col="red",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU,col="blue",lwd = 2)

legend("topleft",legend=c("Bayes. Lowerbound","Bayes. Upperbound"),lwd=c(2,2),col=c("red","blue"),bty="n",cex=1.5) 


dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionval2.eps")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecL,col="red",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU,col="blue",lwd = 2)

legend("topleft",legend=c("Bayes. Lowerbound","Bayes. Upperbound"),lwd=c(2,2),col=c("red","blue"),bty="n",cex=1.5) 

dev.off()
#----------------------------------------------------------------------------
CRIT.vecL1 <- c(  2.934319e-13 ,2.934319e-13 ,2.934319e-13 ,2.934319e-13 ,7.202116e-02 ,1.190139e-01 ,1.642150e-01 ,2.065471e-01 ,2.458811e-01 ,2.822949e-01
                  ,3.159167e-01 ,3.468926e-01 ,3.753781e-01 ,4.015328e-01 ,4.255173e-01 ,4.474905e-01 ,4.676073e-01 ,4.860173e-01 ,5.028642e-01 ,5.182885e-01
                  ,5.324429e-01 ,5.455474e-01 ,5.580181e-01 ,5.707850e-01 ,5.867255e-01 ,5.873295e-01 ,5.925554e-01 ,6.289599e-01 ,6.334299e-01 ,6.398571e-01
                  ,6.435380e-01 ,6.488219e-01 ,6.518625e-01 ,6.562209e-01 ,6.587394e-01 ,6.623449e-01 ,6.644356e-01 ,6.674254e-01 ,6.691644e-01 ,6.716489e-01
                  ,6.730976e-01 ,6.751656e-01 ,6.763741e-01 ,6.780980e-01 ,6.791072e-01 ,6.805459e-01 ,6.813894e-01 ,6.825913e-01 ,6.832968e-01 ,6.843017e-01
                  ,6.846413e-01 ,6.853105e-01 ,6.886756e-01 ,6.892145e-01 ,6.894775e-01 ,6.898210e-01 ,6.902215e-01 ,6.904170e-01 ,6.906724e-01 ,6.909701e-01
                  ,6.911156e-01 ,6.913056e-01 ,6.915270e-01 ,6.916352e-01 ,6.917766e-01 ,6.919413e-01 ,6.920219e-01 ,6.921271e-01 ,6.922496e-01 ,6.923096e-01
                  ,6.923879e-01 ,6.924791e-01 ,6.925237e-01 ,6.925820e-01 ,6.926499e-01 ,6.926795e-01 ,6.927212e-01 ,6.927689e-01 ,6.928142e-01 ,6.929662e-01
                  ,6.929783e-01 ,6.929940e-01 ,6.930117e-01 ,6.930289e-01 ,6.930368e-01 ,6.930471e-01 ,6.930586e-01 ,6.930698e-01 ,6.930750e-01 ,6.930817e-01
                  ,6.930893e-01 ,6.930966e-01 ,6.931000e-01 ,6.931044e-01 ,6.931093e-01 ,6.931141e-01 ,6.931163e-01 ,6.931192e-01 ,6.931225e-01 ,6.931256e-01
                  ,6.931270e-01 ,6.931289e-01 ,6.931310e-01 ,6.931329e-01 ,6.931347e-01 ,6.931363e-01 ,6.931376e-01 ,6.931388e-01 ,6.931399e-01 ,6.931408e-01
                  ,6.931416e-01 ,6.931423e-01 ,6.931429e-01 ,6.931434e-01 ,6.931439e-01 ,6.931443e-01 ,6.931447e-01 ,6.931450e-01 ,6.931452e-01 ,6.931455e-01
                  ,6.931457e-01 ,6.931459e-01 ,6.931460e-01 ,6.931462e-01 ,6.931463e-01 ,6.931464e-01 ,6.931464e-01 ,6.931465e-01 ,6.931466e-01 ,6.931467e-01
                  ,6.931467e-01 ,6.931468e-01 ,6.931468e-01 ,6.931469e-01 ,6.931469e-01 ,6.931470e-01 ,6.931470e-01 ,6.931470e-01 ,6.931470e-01 ,6.931470e-01
                  ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01 ,6.931471e-01
                  ,6.931471e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01
                  ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01
                  ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01
                  ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01
                  ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01 ,6.931472e-01)


CRIT.vecU1 <- c(  2.935430e-13 ,2.935430e-13 ,2.935430e-13 ,2.935430e-13 ,7.497161e-02 ,1.276895e-01 ,1.819157e-01 ,2.364586e-01 ,2.910638e-01 ,3.456807e-01
                  ,4.002999e-01 ,4.549194e-01 ,5.095391e-01 ,5.641590e-01 ,6.187796e-01 ,6.734019e-01 ,7.280286e-01 ,7.826656e-01 ,8.373256e-01 ,8.920449e-01
                  ,9.469613e-01 ,1.002600e+00 ,1.060649e+00 ,1.126214e+00 ,1.219594e+00 ,1.222816e+00 ,1.253880e+00 ,1.587188e+00 ,1.637232e+00 ,1.720961e+00
                  ,1.771005e+00 ,1.854734e+00 ,1.904778e+00 ,1.988507e+00 ,2.038551e+00 ,2.122280e+00 ,2.172324e+00 ,2.256053e+00 ,2.306097e+00 ,2.389826e+00
                  ,2.439870e+00 ,2.523599e+00 ,2.573643e+00 ,2.657371e+00 ,2.707416e+00 ,2.791144e+00 ,2.841189e+00 ,2.924917e+00 ,2.974961e+00 ,3.058690e+00
                  ,3.082829e+00 ,3.140910e+00 ,3.603739e+00 ,3.707176e+00 ,3.757741e+00 ,3.831957e+00 ,3.935394e+00 ,3.985959e+00 ,4.060175e+00 ,4.163612e+00
                  ,4.214177e+00 ,4.288393e+00 ,4.391830e+00 ,4.442395e+00 ,4.516611e+00 ,4.620048e+00 ,4.670613e+00 ,4.744829e+00 ,4.848266e+00 ,4.898831e+00
                  ,4.973048e+00 ,5.076484e+00 ,5.127049e+00 ,5.201266e+00 ,5.304702e+00 ,5.346533e+00 ,5.415687e+00 ,5.506103e+00 ,5.605815e+00 ,6.118790e+00
                  ,6.169404e+00 ,6.243646e+00 ,6.338763e+00 ,6.449917e+00 ,6.500531e+00 ,6.574773e+00 ,6.669890e+00 ,6.781044e+00 ,6.831659e+00 ,6.905900e+00
                  ,7.001018e+00 ,7.112171e+00 ,7.162786e+00 ,7.237027e+00 ,7.332145e+00 ,7.443299e+00 ,7.493913e+00 ,7.568155e+00 ,7.663272e+00 ,7.774426e+00
                  ,7.822755e+00 ,7.895647e+00 ,7.989345e+00 ,8.092252e+00 ,8.197639e+00 ,8.303576e+00 ,8.409629e+00 ,8.515706e+00 ,8.621787e+00 ,8.727870e+00
                  ,8.833952e+00 ,8.940035e+00 ,9.046117e+00 ,9.152200e+00 ,9.258282e+00 ,9.364365e+00 ,9.470447e+00 ,9.576530e+00 ,9.682612e+00 ,9.788695e+00
                  ,9.894775e+00 ,1.000085e+01 ,1.010690e+01 ,1.021324e+01 ,1.032635e+01 ,1.037643e+01 ,1.045034e+01 ,1.054496e+01 ,1.064877e+01 ,1.075506e+01
                  ,1.086189e+01 ,1.096884e+01 ,1.107581e+01 ,1.118279e+01 ,1.128977e+01 ,1.139675e+01 ,1.150373e+01 ,1.161071e+01 ,1.171769e+01 ,1.182467e+01
                  ,1.193165e+01 ,1.203863e+01 ,1.214561e+01 ,1.225258e+01 ,1.235956e+01 ,1.246654e+01 ,1.257352e+01 ,1.268045e+01 ,1.278766e+01 ,1.290118e+01
                  ,1.295169e+01 ,1.302586e+01 ,1.312071e+01 ,1.322475e+01 ,1.333127e+01 ,1.343833e+01 ,1.354551e+01 ,1.365271e+01 ,1.375992e+01 ,1.386712e+01
                  ,1.397433e+01 ,1.408154e+01 ,1.418875e+01 ,1.429595e+01 ,1.440316e+01 ,1.451037e+01 ,1.461758e+01 ,1.472479e+01 ,1.483199e+01 ,1.493920e+01
                  ,1.504641e+01 ,1.515361e+01 ,1.526077e+01 ,1.536819e+01 ,1.548179e+01 ,1.553240e+01 ,1.560663e+01 ,1.570154e+01 ,1.580564e+01 ,1.591221e+01
                  ,1.601932e+01 ,1.612656e+01 ,1.623381e+01 ,1.634108e+01 ,1.644834e+01 ,1.655560e+01 ,1.666286e+01 ,1.677013e+01 ,1.687739e+01 ,1.698465e+01
                  ,1.709192e+01 ,1.719918e+01 ,1.730644e+01 ,1.741370e+01 ,1.752097e+01 ,1.762823e+01 ,1.773549e+01 ,1.784270e+01 ,1.795017e+01 ,1.806378e+01)


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionvalall4.jpeg")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecL1,col="blue",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU1,col="red",lwd = 2)
#lines(seq(1:nrow(xstart)),CRIT.vecL,col="green",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU,col="darkgreen",lwd = 2)

legend("topleft",legend=c("Bayes. Lowerbounds (both)","Bayes. Upperbound","Bayes. Upperbound1"),lwd=c(2,2,2),col=c("blue","red","darkgreen"),bty="n",cex=1.5) 


dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/BayesUpp/criterionvalall4.eps")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vecU),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");


xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vecU),by=2)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vecL1,col="blue",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU1,col="red",lwd = 2)
#lines(seq(1:nrow(xstart)),CRIT.vecL,col="green",lwd = 2)
lines(seq(1:nrow(xstart)),CRIT.vecU,col="darkgreen",lwd = 2)

legend("topleft",legend=c("Bayes. Lowerbounds (both)","Bayes. Upperbound","Bayes. Upperbound1"),lwd=c(2,2,2),col=c("blue","red","darkgreen"),bty="n",cex=1.5) 

dev.off()