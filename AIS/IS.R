## target function f(x, y)
f <- function(x, y)
{
  return(0.5*exp(-90*(x-0.5)^2-45*(y+0.1)^4) + exp(-45*(x+0.4)^2-60*(y-0.5)^2))
}

x = seq(-1, 1, by = 0.1)
y = x
fxy = outer(x, y, f)
persp(x, y, z = fxy, theta = 30)

## naive
estimate.simple <- function()
{
  m = 2500
  xi = runif(m, -1, 1)
  yi = runif(m ,-1, 1)
  fi = f(xi, yi)
  mu.hat = 4/m*sum(fi)
  mu.hat.std = 4*sqrt(sum((fi-mu.hat)^2)/(m*(m-1)))
  return(c(mu.hat, mu.hat.std))
}

## repeat n times
n = 100
res = sapply(1:n, function(x) estimate.simple())
summary(res[1,])
summary(res[2,])

## Importance Sampling
mu1 = c(0.5, -0.1)
sigma1 = diag(c(1/180, 1/20))
mu2 = c(-0.4, 0.5)
sigma2 = diag(c(1/90, 1/120))

library(mvtnorm)
w1 = 0.5*sqrt(det(sigma1))
w2 = sqrt(det(sigma2))
w = c(w1, w2)
w = w/sum(w)
w # 0.464 N1 + 0.536 N2

## choose g(x, y)
g <- function(x, y)
{
  return(0.5*exp(-90*(x-0.5)^2-45*(y+0.1)^2) + exp(-45*(x+0.4)^2-60*(y-0.5)^2))
}

estimate.is <- function()
{
  m = 2500
  ## construct dataset
  flag = (runif(m) < w[1])
  data1 = rmvnorm(m, mu1, sigma1)
  data2 = rmvnorm(m, mu2, sigma2)
  data = rbind(data1[flag,], data2[!flag, ])
  w = numeric(m)
  for (i in 1:m)
  {
    if (abs(data[i, 1]) > 1 || abs(data[i, 2]) > 1)
      w[i] = 0
    else
      w[i] = f(data[i, 1], data[i, 2])/g(data[i, 1], data[i, 2])
  }
  fxy = sapply(1:m, function(i) f(data[i, 1], data[i, 2]))
  mu.hat = sum(w*fxy)/sum(w)
  ## how to calculate stderr 
  return(mu.hat)
}


## repeat n times
n = 100
res = sapply(1:n, function(x) estimate.is())
summary(res)
