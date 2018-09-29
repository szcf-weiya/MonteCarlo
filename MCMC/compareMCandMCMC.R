# R code for MCMC diagnostics
# author: szcfweiya
library(mvtnorm)
# discrete variable
delta = c(.45, .10, .45)

# continuous variable
mu = c(-3, 0, 3)
sigma2 = c(1/3, 1/3, 1/3)

# exact marginal density of theta
ext_margin_den <- function(x)
{
  dnorm(x, mu[1], sqrt(sigma2[1])) * delta[1] +
    dnorm(x, mu[2], sqrt(sigma2[2])) * delta[2] +
    dnorm(x, mu[3], sqrt(sigma2[3])) * delta[3]
}

theta = seq(-6, 6, length.out = 1000)
ptheta = ext_margin_den(theta)
plot(theta, ptheta, type = "l")

# independent Monte Carlo samplers
# step 1: sample delta
# step 2: sample theta
n = 1000
sample.delta = sample(1:3, n, prob = delta, replace = TRUE)
sample_theta <- function(cl)
{
  if (cl == 1)
    rnorm(1, mu[1], sqrt(sigma2[1]))
  else if (cl == 2)
    rnorm(1, mu[2], sqrt(sigma2[2]))
  else
    rnorm(1, mu[3], sqrt(sigma2[3]))
}
sample.theta = sapply(sample.delta, sample_theta)
hist(sample.theta, probability = TRUE, add=T, breaks = 20)

sample_delta <- function(x)
{
  denx = c(delta[1] * dnorm(x, mu[1], sqrt(sigma2[1])),
           delta[2] * dnorm(x, mu[2], sqrt(sigma2[2])),
           delta[3] * dnorm(x, mu[3], sqrt(sigma2[3])))
  condenx = denx / sum(denx)
  sample(1:3, 1, prob = condenx)
}
# gibbs sampler
gibbs <- function(m, delta = 3)
{
  DELTA = delta
  THETA = NULL
  for (i in 1:m)
  {
    # sample theta from full conditional distribution
    theta = sample_theta(delta)
    THETA = c(THETA, theta)
    delta = sample_delta(theta)
    DELTA = c(DELTA, delta)
  }
  return(THETA)
}
res = gibbs(10000, delta = 1)
par(mfrow=c(1,2))
plot(res)
hist(res, breaks = 20)
#res2 = ext_margin_den(res)
#plot(log(res2))