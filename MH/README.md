# Metropolis Algorithm

1. Monte Carlo method is sufficient for conjugate prior distribution and Gibbs sampler also can handle semiconjugate prior distribution, but when the conjugate or semiconjugate distribution is unavailable, we need the Metropolis-Hastings algorithm which is a generic method of approximating the posterior distribution.

## Procedure

1. sample $$\theta^*\sim J(\theta\mid \theta^{(s)})$$;
2. Compute the acceptance ratio
$$
r=\frac{p(\theta^*\mid y)}{p(\theta^{(s)}\mid y)}=\frac{p(y\mid\theta^*)p(\theta^*)}{p(y\mid\theta^{(s)})p(\theta^{(s)})}
$$
3. $$theta^{(s+1)=\theta^*\; \mathrm{w.p.}\; \mathrm{min}(r,1)$$, and $$theta^{(s+1)=\theta^{(s)}\; \mathrm{w.p.}\; 1-\mathrm{min}(r,1)$$

In practice, the first $$B$$ iterations are called "burn-in" period, which should be discarded.

## Example: Normal distribution with known variance

Let $$\theta\sim N(\mu,\tau^2),\;\{y_1,\ldots,y_n\mid \theta\}\sim i.i.d. N(\theta,\sigma^2)$$, then the posterior distribution of $$\theta$$ is $N(\mu_n,\tau_n^2)$.

The simulation is as follows:

```r
s2 = 1; t2 = 10; mu = 5;
y = c(9.37, 10.18, 9.16, 11.60, 10.33)
theta = 0; delta2 = 2; S = 10000; THETA = NULL;
set.seed(1)

for (s in 1:S)
{
  theta.star = rnorm(1, theta, sqrt(delta2))
  log.r = (sum(dnorm(y, theta.star, sqrt(s2), log = TRUE)) +
             dnorm(theta.star, mu, sqrt(t2), log = TRUE)) -
          (sum(dnorm(y, theta, sqrt(s2), log = TRUE)) + 
             dnorm(theta, mu, sqrt(t2), log = TRUE))
  if (log(runif(1)) < log.r){
    theta = theta.star
  }
  THETA = c(THETA, theta)
}

# figure 10.3 
## left
plot(THETA)
## right
hist(THETA, probability = TRUE, breaks = 70, xlim = c(8, 12))
## true posterior distribution
lines(sort(THETA), dnorm(sort(THETA), 10.03, 0.44))
```

## For Poisson regression

The data can be found on [Peter Hoff's personal page](https://www.stat.washington.edu/~pdhoff/book.php).

```r
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
```