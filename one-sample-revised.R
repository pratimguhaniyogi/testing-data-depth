rm(list = ls())
library(mvnTest) #for existing testing methods
library(ddalpha) #depth calculation
library(mvtnorm) #normal and t distribution generation

multivzscore=function(x) {  #x is the data set
  x <- as.matrix(x)
  n <- nrow(x); p=ncol(x) #dimension of the data
  s <- cov(x) #covariance of the data
  m <- matrix(rep(colMeans(x),n),byrow=TRUE,ncol=p) #mean of the data
  eigen.cov <- eigen(s)
  lam <- eigen.cov$values #eigenvalues of covariance matrix
  vec <- eigen.cov$vector #eigenvectors of covariance matrix
  B <- vec%*%diag(1/sqrt(lam))%*%t(vec)
  z <- (x-m)%*%B #standardized data
  return(z)
}

B <- 500 #bootstrap replicate
n.approx <- 1000 #number of sample size that are used to build the true distribution
M <- 1000 #number of observations for to generate values of points at which statistics is calculated

n <- 100 #sample size
d <- 2 #dimension

# declare the path and the file name
path <- paste0("/Users/pratimguhaniyogi/Desktop/data-depth/github-code/one-sample-n", n, "-d", d, ".rds")

# the points at which we calculate depth. These are generated uniformly
radius <- 1 #radius
mult.dat <- matrix(rnorm(M*d), M, d) #multivariate data from N(0, I)
mult.dat.s <- t(sapply(1:M, function(aa){
  mult.dat[aa, ]/sqrt(sum(mult.dat[aa, ]^2))
}))
u <- runif(M) #uniform(0, 1)
r <- radius*u^(1/d) #r is proportional to d-th sqrt of U
t0 <- r*mult.dat.s

dat <- mvrnorm(n, rep(0, d), diag(d)) ## generate the data (change for different cases)
dat <- multivzscore(dat)

x <- rmvnorm(n.approx, mean =  rep(0, d), sigma = diag(d)) #null model N(0, I)
x <-  multivzscore(x)

s0 <- rmvnorm(M, mean =  rep(0, d), sigma = diag(d)) # from the null model
s0 <- multivzscore(s0)

Dn.KS <- depth.halfspace(t0, dat, exact = FALSE)
D0.KS <- depth.halfspace(t0, x, exact = FALSE)

Dn.CM <- depth.halfspace(s0, dat, exact = FALSE)
D0.CM <- depth.halfspace(s0, x, exact = FALSE)

depth.diff.KS <- Dn.KS - D0.KS
depth.diff.CM <- Dn.CM - D0.CM

depth.KS <- sqrt(n)*max(abs(depth.diff.KS))
depth.CM <- n*mean(depth.diff.CM^2)

depth.diff.boot.KS <- sapply(1:B, function(bb){
  set.seed(bb)
  ## parametric bootstrap
  dat.boot <- rmvnorm(n, mean =  colMeans(dat), sigma = cov(dat))
  dat.boot <- multivzscore(dat.boot)
  Dn.boot.KS <- depth.halfspace(t0, dat.boot, exact = FALSE)
  return(Dn.boot.KS-D0.KS)
})
depth.diff.boot.CM <- sapply(1:B, function(bb){
  set.seed(bb)
  ## parametric bootstrap
  dat.boot <- rmvnorm(n, mean =  colMeans(dat), sigma = cov(dat))
  dat.boot <- multivzscore(dat.boot)
  Dn.boot.CM <- depth.halfspace(s0, dat.boot, exact = FALSE)
  return(Dn.boot.CM-D0.CM)
})

depth.KS.boot <- sqrt(n)*apply(abs(depth.diff.boot.KS), 2, max)
depth.CM.boot <- n*apply(depth.diff.boot.CM^2, 2, mean)

pval.depth.KS <- mean(depth.KS.boot > depth.KS)
pval.depth.CM <- mean(depth.CM.boot > depth.CM)

multi.ad <- AD.test(dat)
multi.cm <- CM.test(dat) 
multi.dh <- DH.test(dat)
multi.hz <- HZ.test(dat)
multi.r <- R.test(dat)
multi.chisq <- S2.test(dat)

output <- list(n = n, d = d, 
               dat = dat, 
               x = x, 
               
               depth.KS = depth.KS, 
               depth.CM = depth.CM,
               
               depth.KS.boot = depth.KS.boot, 
               depth.CM.boot = depth.CM.boot,
               
               pval.KS.depth = pval.depth.KS, 
               pval.CM.depth = pval.depth.CM, 
               pval.AD = multi.ad@p.value, 
               pval.CM = multi.cm@p.value,
               pval.DH = multi.dh@p.value, 
               pval.HZ = multi.hz@p.value, 
               pval.R = multi.r@p.value, 
               pval.McCulloch = multi.chisq@p.value.s2, 
               pval.NRR = multi.chisq@p.value.y2, #Nikulin-Rao-Robson
               pval.DN = multi.chisq@p.value.u2#Dzhaparidze-Nikulin
)
saveRDS(output, path)

