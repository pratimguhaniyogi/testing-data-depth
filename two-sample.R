rm(list = ls())
library(ddalpha) #depth calculation
library(mvtnorm) #normal and t distribution generation
#This code produce the p-value of the proposed tests
#X ~ F and Y ~ G where F = N(0,I) and G = N(mu x 1, I)
#H_0: F = G against H_1: F != G (Model B)
#The distributions of the data sets can be changed easily

B <- 500 #bootstrap replicate
M <- 1000 #number of observations for to generate values of points at which statistics is calculated

n <- 100; m <- 50 #sample size
d <- 2 #dimension
mu <- 0.1 #parameter under alternative

path <- paste0("/Users/pratimguhaniyogi/Desktop/data-depth/github-code/two-sample-n", n, "-m", m,"-d", d, "-mu", mu, ".rds")
if(d==2){
  r <- sqrt(runif(M))
  theta <- runif(M, 0, 2*pi)
  t0 <- cbind(r*cos(theta), r*sin(theta))  
}

if(d>3){
  r <- sqrt(runif(M))
  theta <- matrix(runif(M*(d-2), 0, pi), M, d-2)
  theta <- cbind(theta, runif(M, 0, 2*pi))
  sin.theta <- sin(theta)
  cos.theta <- cos(theta)  
  
  u <- list()
  u[[1]] <- r*cos.theta[, 1]
  for(jj in 2:(d-1)){ 
    u[[jj]] <- r*apply(as.matrix(sin.theta[, 1:(jj-1)]), 1, prod)*cos.theta[, jj]
  }
  u[[d]] <- r*apply(as.matrix(sin.theta), 1, prod)
  t0 <- do.call(cbind, u)
}  
  
x <- rmvnorm(n, mean = rep(0, d), sigma = diag(d)) # sample from F
y <- rmvnorm(m, mean = rep(mu, d), sigma = diag(d)) # sample from G
dat <- rbind(x, y) # combined data

Dx <- depth.halfspace(t0, x, exact = FALSE) # data-depth value based on data-1 evaluated at t0
Dy <- depth.halfspace(t0, y, exact = FALSE) # data-depth value based on data-2 evaluated at t0
depth.diff <- Dx - Dy # data-depth difference 
  
depth.KS <- sqrt(n+m)*max(abs(depth.diff)) # KS-type test-statistics
depth.CM <- (n+m)*mean(depth.diff^2) # CvM-type test-statistics
  
depth.diff.boot <- sapply(1:B, function(bb){ # boostrap 
  set.seed(bb)
  index <- sample(1:(n+m), replace = TRUE, size = n+m)
  x.boot <- dat[index[1:n], ]
  y.boot <- dat[index[(n+1):(n+m)], ]
  Dx.boot <- depth.halfspace(t0, x.boot, exact = FALSE)
  Dy.boot <- depth.halfspace(t0, y.boot, exact = FALSE)
  return(Dx.boot - Dy.boot)
})  
  
depth.KS.boot <- sqrt(n+m)*apply(abs(depth.diff.boot), 2, max) # KS-type test-statistics (bootstrap)
depth.CM.boot <- (n+m)*apply(depth.diff.boot^2, 2, mean) # CvM-type test-statistics (bootstrap)
  
pval.depth.KS <- mean(depth.KS.boot > depth.KS)
pval.depth.CM <- mean(depth.CM.boot > depth.CM)
  
output <- list(seed = aa, n = n, m = m, d = d,
               x = x, y = y,
               
               depth.KS.boot = depth.KS.boot, 
               depth.CM.boot = depth.CM.boot,
               
               pval.KS.depth = pval.depth.KS, 
               pval.CM.depth = pval.depth.CM
)  
  
saveRDS(output, path)
