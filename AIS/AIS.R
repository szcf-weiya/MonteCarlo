## multivariate
library(mvtnorm)

## parameters
e = c(1, 1, 1, 1)
I4 = diag(4)
D4 = diag(c(2, 1, 1, .5))

## generate data; d individuals
d = 100
data = rmvnorm(d, e, 2*I4) + 2*rmvnorm(d, 3*e, I4) + 1.5*rmvnorm(d, -3*e, D4)

## AIS
m = 200
nu = 2
pix <- function(x)
{
  return(1/(sqrt(2*2))*exp(-0.5*0.5*t(x)%*%x) + 
         2* 1/sqrt(2*1)*exp(-0.5*t(x-3*e)%*%(x-3*e)) +
           1.5*1/sqrt(2*det(D4))*exp(-0.5*t(x+3*e)%*%solve(D4)%*%(x+3*e)))
}

calculateCriteria <- function(mu, sigma)
{
  samples.g0 = rmvt(m, delta = mu, sigma = sigma, nu)
  weights.g0 = apply(samples.g0, 1, function(x) pix(x))
  weights.g0 = weights.g0/sum(weights.g0)*m
  mu1 = apply(samples.g0, 2, function(x) sum(x*weights.g0))/m
  sigma1 = cov(apply(samples.g0, 2, function(x) x*weights.g0))
  measure = var(weights.g0)
  return(list(mu = mu1, 
              sigma = sigma1,
              criteria = measure))
}

## start a trial density
sigma = diag(4)
mu = c(0, 0, 0, 0)
AIS.simple <- function(mu, sigma)
{
  criteria = 1e10
  while(TRUE)
  {
    res = calculateCriteria(mu, sigma)
    if (res$criteria < criteria)
    {
      criteria = res$criteria
      cat(paste0("criteria = ", criteria, "\n"))
      mu = res$mu
      sigma = res$sigma
    }
    else
      break
  }
  return(list(mu = mu,
              sigma = sigma))
  
}