#!/usr/bin/env Rscript

## R program for Comparing two groups

## get data
## you can found it on Peter Hoff's website
## https://www.stat.washington.edu/~pdhoff/book.php

Y.school.mathscore<-dget("Y.school.mathscore")
y.school1<-dget("y.school1")
y.school2<-dget("y.school2")

y1 = y.school1; n1 = length(y1)
y2 = y.school2; n2 = length(y2)

## prior parameters
mu0 = 50; g02 = 625
del0 = 0; t02 = 625
s20 = 100; nu0 = 1

## starting value
mu = (mean(y1) + mean(y2)) / 2
del = (mean(y1) - mean(y2)) / 2

## gibbs sampler
MU = DEL = S2 = NULL
set.seed(1)
for (s in 1:5000)
{
  # update s2
  s2 = 1/rgamma(1, (nu0+n1+n2)/2, (nu0*s20+sum((y1-mu-del)^2)+sum((y2-mu+del)^2))/2)
  
  # update mu
  var.mu = 1/(1/g02 + (n1+n2)/s2)
  mean.mu = var.mu*(mu0/g02+sum(y1-del)/s2+sum(y2+del)/s2)
  mu = rnorm(1, mean.mu, sqrt(var.mu))
  
  # update del
  var.del = 1/(1/t02 + (n1+n2)/s2)
  mean.del = var.del*(del0/t02 + sum(y1-mu)/s2 - sum(y2 - mu)/s2)
  del = rnorm(1, mean.del, sqrt(var.del))
  
  # save parameter values
  MU = c(MU, mu)
  DEL = c(DEL, del)
  S2 = c(S2, s2)
}

# plot
png("reproduce-fig-8-2l.png")
plot(density(MU), main="", xlab=expression(mu), lwd=2, col="black")
lines(density(rnorm(5000, 50, 25)), col="red", lwd=2)
legend("topleft", c("prior", "posterior"), col=c("black","red"), lwd=2)
dev.off()

png("reproduce-fig-8-2r.png")
plot(density(DEL), main="", xlab=expression(delta), lwd=2, col="black")
lines(density(rnorm(5000, 0, 25)), col="red", lwd=2)
legend("topleft", c("prior", "posterior"), col=c("black","red"), lwd=2)
dev.off()


## comparing multiple groups
Y = Y.school.mathscore
## weakly informative priors
nu0 = 1; s20 = 100
eta0 = 1; t20 = 100
mu0 = 50; g20 = 25

## starting values
m = length(unique(Y[, 1]))
n = sv = ybar = rep(NA, m)

for (j in 1:m)
{
  ybar[j] = mean(Y[Y[, 1]==j, 2])
  sv[j] = var(Y[Y[, 1]==j, 2])
  n[j] = sum(Y[, 1]==j)
}
theta = ybar
sigma2 = mean(sv)
mu = mean(theta)
tau2 = var(theta)

## setup MCMC
set.seed(1)
S = 5000
THETA = matrix(nrow = S, ncol = m)
SMT = matrix(nrow = S, ncol = 3)

## MCMC algorithm
for (s in 1:S)
{
  # sample new values of the thetas
  for (j in 1:m)
  {
    vtheta = 1/(n[j]/sigma2+1/tau2)
    etheta = vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j] = rnorm(1, etheta, sqrt(vtheta))
  }
  # sample new value of sigma2
  nun = nu0 + sum(n)
  ss = nu0*s20
  for (j in 1:m)
  {
    ss = ss + sum(Y[Y[,1]==j, 2]-theta[j])^2
  }
  sigma2 = 1/rgamma(1, nun/2, ss/2)
  
  # sample new value of mu
  vmu = 1/(m/tau2+1/g20)
  emu = vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu = rnorm(1, emu, sqrt(vmu))
  
  # sample a new value of tau2
  etam = eta0 + m
  ss = eta0*t20 + sum((theta-mu)^2)
  tau2 = 1/rgamma(1, etam/2, ss/2)
  
  # store results
  THETA[s, ] = theta
  SMT[s, ] = c(sigma2, mu, tau2)
}
