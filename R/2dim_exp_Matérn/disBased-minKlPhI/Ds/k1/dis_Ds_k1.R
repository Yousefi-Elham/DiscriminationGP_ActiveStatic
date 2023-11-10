#Distance based


library(mvtnorm)
library(numDeriv)
library(dplyr) #anti_join
rm(list=ls())
par(mfrow=c(1,1))

####################################################################################
sigma2.ini=1
theta1=1
theta2=1.07 ##NEW par
est.theta <- c(theta1,theta2)
nu.set <- c(3/2,5/2)

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

CRIT.vec <- numeric(nexc)
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
crit.fun.old <- function(nu,theta,sigma2,des){
  
  
  n <- nrow(des)
  DisM <- as.matrix(dist(des))
  
  
  h <-  sqrt(2*nu)*DisM*theta 
  partial.h <- sqrt(2*nu)*DisM
  
  
  kernel.nu <- sigma2*2^(1-nu)/gamma(nu)*h^(nu)*besselK(h,nu) #compare with kernel2
  diag(kernel.nu) <- 1
  
  
  partial.theta <- -sigma2*2^(1-nu)/gamma(nu)*h^nu*besselK(h,nu-1)*partial.h
  diag(partial.theta) <- 0 
  
  
  bessel.fun <- function(ka){
    B <- besselK(sqrt(2*ka)*DisM*theta,ka)
    return(B)
  }
  bessel.par <- jacobian(func=bessel.fun, x=nu, method="Richardson")
  bessel.par.mat <- matrix(bessel.par,n,n)
  diag(bessel.par.mat) <- 0
  
  
  logh.half <- log(h/2)
  diag(logh.half) <- 0 
  
  
  partial.nu <- sigma2*2^(1-nu)/gamma(nu)*h^nu *( ((logh.half + 1) - digamma(nu))*besselK(h,nu) + bessel.par.mat )
  diag(partial.nu) <- 0
  
  
  if(det(kernel.nu)<=0){
    Ds.crit <- 0
  }else{
    
    kernelnu.inv <- solve(kernel.nu)
    
    kernel.inv_th <- kernelnu.inv%*%partial.theta
    kernel.inv_nu <- kernelnu.inv%*%partial.nu
    
    
    FIM <- matrix(nrow = 2, ncol = 2)
    FIM[1,1] <- 0.5 * sum(t(kernel.inv_th)*kernel.inv_th) # double derivative for theta
    FIM[2,2] <- 0.5 * sum(t(kernel.inv_nu)*kernel.inv_nu) # double derivative for nu
    FIM[1,2] <- FIM[2,1] <- 0.5 * sum(t(kernel.inv_th)*kernel.inv_nu) # cross derivative between theta and nu
    
    Ds.crit <- det(FIM)/FIM[1,1]
    
  }
  
  return(Ds.crit)
  
}
#-------------------------------------------------------------------------------------------------
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

crit.fun.old(nu.set[1],est.theta[1],sigma2.ini,xstart)

#----------------------cunstructing the rest (of design) up to the required size-------------------
system.time({
  
  set.seed(123456789)
  for( k in 1: nleft){
    
    
    print(k)
    
    nc <- nrow(xstart)
    Ncan <- nrow(xcan)
    
    
    ad.crit <- numeric(Ncan)
    for(m in 1: Ncan){
      
      test.des <- rbind(xstart,xcan[m,])
      ad.crit[m] <- crit.fun.old(nu.set[1],est.theta[1],sigma2.ini,test.des)
      
    }
    
    if(max(ad.crit)==0){
      print(k-1)
      stop("best design chosen has a zero crit. value")
    }
    
    ad.ind <- which.max(ad.crit)
    ad.point <- xcan[ad.ind,]
    xstart <- rbind(xstart,ad.point)
    xcan <- xcan[-ad.ind,]
    
    CRIT.vec[k+n1] <- max(ad.crit)
    #--------------------------------
    CRIT.val <- max(ad.crit)
    
    
    exnum <- 0
    finish.all <- FALSE
    
    while(!finish.all){
      
      nc1 <- nrow(xstart)
      dr.crit <- numeric(nc1)
      
      for(a in 1: nc1){
        dr.crit[a] <- crit.fun.old(nu.set[1],est.theta[1],sigma2.ini,xstart[-a,])
        
      }
      dr.ind <- which.max(dr.crit)
      dr.point <- xstart[dr.ind,]
      xcan <- rbind(xcan, dr.point)
      xstart.dr <- xstart[-dr.ind,]
      
      
      Ncan1 <- nrow(xcan)
      add.crit <- numeric(Ncan1)
      
      
      
      for(b in 1: Ncan1){
        test.des <- rbind(xstart.dr,xcan[b,])
        add.crit[b] <- crit.fun.old(nu.set[1],est.theta[1],sigma2.ini,test.des)
        
        
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
    
    CRIT.vec[k+n1] <- crit.fun.old(nu.set[1],est.theta[1],sigma2.ini,xstart)
    
    
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

length(CRIT.vec)
CRIT.vec <- CRIT.vec[1:nrow(xstart)]
length(CRIT.vec)
CRIT.vec
xstart
#---------------#---------------
av.sd.iter <- apply(M.COUNT.N,1,function(x) sqrt((NRUN-1)/NRUN*var(x)))

av.sd.iter

av.sd.iter[int.hit.n]
round(av.sd.iter[int.hit.n],digits = 4)

count2 <- c(1,2,3,4,5,6,16,26,36,46)

av.sd.iter[count2]
round(av.sd.iter[count2],digits = 4)

av.hit.iter[count2]
#---------------#---------------
#####################################################################################################
#####################################################################################################


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/design.jpeg")

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
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/design.eps")

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


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/des.seg.jpeg")

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
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/des.seg.eps")

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

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/des.seg1.jpeg", width=650, height=500)

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
pdf("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/des.seg1.eps", width = 11, height = 7)

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


jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/DesDis.jpeg")

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
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/DesDis.eps")

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
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/DesDis1.eps", width = 11, height = 7)

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

max(CRIT.vec)

jpeg("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/criterionval.jpeg")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vec),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vec),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vec,col="blue",lwd = 2)

dev.off()

#----------------------------------------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2022/2-Feb/GP-project-R/Matern1.5-2.5-newpar/2dim_exp/disBased-minKlPhI/Ds/k1/criterionval.eps")

plot(seq(1:nrow(xstart)),seq(0,max(CRIT.vec),length=nrow(xstart)),axes = F, type = "n",xlab="Iter",ylab="Criterion value",
     main = bquote( paste("Criterion values for optimum designs",", ", sigma^2, ":",.(sigma2.ini)) ),
     cex.main=2, cex.lab=1.5, cex.axis=1.5,pch=2,col="blue");

xtick=c(seq(1,nrow(xstart),by = 20),nrow(xstart))
ytick=seq(0,max(CRIT.vec),by=1000)
axis(1,at=xtick,labels=T,cex.axis=1.5)
axis(2,at=ytick,labels=T,cex.axis=1.5)
box()

lines(seq(1:nrow(xstart)),CRIT.vec,col="blue",lwd = 2)

dev.off()
