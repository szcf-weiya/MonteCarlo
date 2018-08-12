## get data
## you can found it on Peter Hoff's website
## https://www.stat.washington.edu/~pdhoff/book.php

#yX.sparrow = dget(yX.sparrow)

y = yX.sparrow[, 1]
X = yX.sparrow[, -1]
n = length(y)
p = dim(X)[2]

# prior expectation
pmn.beta = rep(0, p)
psd.beta = rep(10, p)
# proposal var
var.prop = var(log(y+1/2))*solve(t(X) %*% X)
S = 10000
beta = rep(0, p)
acs = 0
BETA = matrix(0, nrow = S, ncol = p)
set.seed(1)

library(MASS)
for (s in 1:S)
{
  beta.p = mvrnorm(1, beta, var.prop)
  lhr = sum(dpois(y, exp(X %*% beta.p), log = T)) -
    sum(dpois(y, exp(X %*% beta), log = T)) + 
    sum(dnorm(beta.p, pmn.beta, psd.beta, log = T)) -
    sum(dnorm(beta, pmn.beta, psd.beta, log = T))
  if (log(runif(1)) < lhr){
    beta = beta.p
    acs = acs + 1
  }
  BETA[s, ] = beta
}
# plot for beta3
plot(BETA[,3])