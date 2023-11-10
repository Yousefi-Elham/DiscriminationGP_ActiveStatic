DKL.1giv2 <- 1/2*(sum(diag(kernel1 %*% solve(kernel2)))-log(det(kernel1 %*% solve(kernel2)))-n)
DKL.2giv1 <- 1/2*(sum(diag(kernel2 %*% solve(kernel1)))-log(det(kernel2 %*% solve(kernel1)))-n)


pp1 <- -1/2*log( 1/2 + 1/2*exp( -DKL.1giv2 ) )       
pp2 <- -1/2*log( 1/2 + 1/2*exp( -DKL.2giv1 ) ) 
pp1+pp2