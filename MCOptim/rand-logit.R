##
## adapt from mcsm's randogit.R (https://github.com/cran/mcsm/blob/master/R/randogit.R)
##               and randogibs.R (https://github.com/cran/mcsm/blob/master/R/randogibs.R)
##

rm(list = ls())

## Random effect logit model
beta0 = -3
sigma0 = 1
n = 20 # num of individuals
m = 35 # num of replicas

## simulated data
x = matrix(sample(c(-1, 0, 1), n*m, replace = TRUE), nrow = n)
y = matrix(0, ncol = m, nrow = n)
tru = rnorm(n)
for (i in 1:n)
  y[i,] = (runif(m) < 1/(1 + exp(-beta0 * x[i,] - tru[i]))) # trick: e^z / (1 + e^z) = 1 / (1 + e^(-z))

## reference value without random effects
## i.e., the MLE of beta when there is no random effect
mlan = as.numeric(glm(as.vector(y) ~ as.vector(x) - 1, family = binomial)$coef) # force no intercept via y ~ x - 1 or y ~ 0 + x (see ?glm -> ?formula)

## target function
targ = function(beta, omegas, Tmc) {
  xs = exp(beta * x)
  xxs = x * xs
  losprod = 0
  for (j in 1:m) for (t in 1:Tmc)
    losprod = losprod + sum(xxs[,j] * omegas[, t] / (1 + xs[,j] * omegas[,t])) # sum all replicas
  
  # find the root
  losprod - Tmc*sum(x*y) # Tmc is due to summation over all replicas
}

## function for MCMC (eq. 5.15)
gu = function(mu, i, beta, sigma) {
  sum((y[i,] * (beta * x[i,] + mu)) - 
      log(1 + exp(beta * x[i,] + mu))) - 0.5*mu^2 / sigma^2
}

## MCEM iterations
MCEM <- function(){
  ## MC likelihood function
  likecomp = function(beta, sigma) {
    # treat u as obs.
    exp(sum(as.vector(y) * (beta * as.vector(x) + rep(tru, m))) - 
          sum(log(1 + exp(beta * as.vector(x) + rep(tru, m)))))
  }
  
  Tmc = 10^2
  beta = mlan
  sigma = sigma0
  
  diff = iter = factor = 1
  lval = likecomp(beta, sigma)
  while (diff > 10^(-2)) {
    # simulate u's by MCMC
    samplu = matrix(tru, ncol = Tmc, nrow = n)
    acpt = 0
    for (t in 2:Tmc){
      u = rnorm(n)
      for (i in 1:n) {
        # proposal
        u[i] = factor * sigma[iter] * u[i] + samplu[i, t-1]
        # NB: different from the original code
        # factor seems only related to the proposal step
        if (log(runif(1)) > gu(u[i], i, beta[iter], sigma[iter]) - gu(samplu[i,t-1], i, beta[iter], sigma[iter]))
          u[i] = samplu[i, t-1]
        else
          acpt = acpt + 1
      }
      samplu[, t] = u
    }
    
    if (acpt < .1*Tmc) factor = factor/3
    if (acpt > .9*Tmc) factor = factor*3
    
    ## EM equation
    sigma = c(sigma, sd(as.vector(samplu)))
    
    omegas = exp(samplu)
    beta = c(beta, uniroot(targ, omegas = omegas, Tmc = Tmc, interval = mlan + 10*sigma[iter+1]*c(-1, 1))$root)
    lval = c(lval, likecomp(beta[iter + 1], sigma[iter + 1]))
    
    diff = max(abs(diff(beta[iter:(iter+1)])), abs(diff(sigma[iter:(iter+1)]))) # interestingly! diff can be numeric and function simultaneously!
    iter = iter + 1
    Tmc = Tmc*2
    cat("iter =", iter, ", diff =", diff, "\n")
  }
  
  par(mfrow = c(2, 1))
  plot(beta, sigma, pch = 19, xlab = expression(beta), ylab = expression(sigma))
  plot(lval, type = "l", xlab = "iter", ylab = expression(L^c)) 
}

## Simple MCMC
SimpleMCMC <- function() {

  # TODO: why?  
  pro = function(beta, u) {
    exp(as.vector(y) * (beta * as.vector(x) + rep(u, m))) / (1 + exp(as.vector(y) * (beta * as.vector(x) + rep(u, m))))^2
  }
  
  # u is random
  likecomp = function(beta, sigma, u) {
    sum(as.vector(y) * (beta * as.vector(x) + rep(u, m))) - 
      sum(log(1 + exp(beta * as.vector(x) + rep(u, m)))) - 
      sum(u^2) / (2*sigma^2) - log(sigma)
  }
  beta = mlan
  sigma = factor = 1
  acpt = bcpt = 0
  u = rnorm(n)
  Tmc = 10^3
  samplu = matrix(u, nrow = n, ncol = Tmc)
  for (iter in 2:Tmc) {
    # sample u
    u = rnorm(n)
    for (i in 1:n) {
      # proposal
      u[i] = factor * sigma[iter-1] * u[i] + samplu[i, iter-1]
      if (log(runif(1)) > gu(u[i], i, beta[iter-1], sigma[iter-1]) - gu(samplu[i,iter-1], i, beta[iter-1], sigma[iter-1]))
        u[i] = samplu[i, iter-1]
      else
        acpt = acpt + 1
    }
    samplu[, iter] = u
    
    # sample beta and sigma
    ## NB: 1/sigma^2 ~ Gamma(n/2+1, u^2/2), but the original code use Gamma(n/2, u^2/2)
    sigma = c(sigma, 1 / sqrt(2*rgamma(1, 0.5*n+1)/sum(u^2)))
    
    tau = sigma[iter] / sqrt(sum(as.vector(x^2) * pro(beta[iter-1], u)))
    betaprop = beta[iter-1] + rnorm(1) * factor * tau
    if (log(runif(1)) > likecomp(betaprop, sigma[iter], u) - likecomp(beta[iter-1], sigma[iter], u))
      betaprop = beta[iter-1]
    else
      bcpt = bcpt + 1
    
    beta = c(beta, betaprop)
    
    if (iter > 100) {
      if (bcpt < .1*iter) factor = factor / 3
      if (bcpt > .9*iter) factor = factor * 3
    }
  }
  library(coda)
  plot(mcmc(cbind(beta, sigma)))
  #browser()
  cumuplot(mcmc(cbind(beta, sigma)))
  list(beta = beta, sigma = sigma)
}