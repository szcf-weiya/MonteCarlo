#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 20:40:53 2017

@author: root
"""
## ###############################
## ref: http://twiecki.github.io/blog/2015/11/10/mcmc-sampling/
## ###############################

## ###############################
## import modules
## ###############################
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import norm

sns.set_style('white')
sns.set_context('talk')

np.random.seed(123)

## ###############################
## generate some data
## ###############################
data = np.random.randn(20)
ax = plt.subplot()
sns.distplot(data, kde = False, ax = ax)
_ = ax.set(title='Histogram of observed data', xlabel='x', ylabel='# observations')

## ###############################
## if our prior was not conjugate
## compute the parameters of the posterior
## 
## 
## mu_0 prior
## sigma_0 prior
## sigma post    
## ###############################
def calc_posterior_analytical(data, x, mu_0, sigma_0):
    sigma = 1.
    n = len(data)
    mu_post = (mu_0/sigma_0**2 + data.sum()/sigma**2) / (1./sigma_0**2 + n/sigma**2)
    sigma_post = (1./sigma_0**2 + n/sigma**2)**-1
    return norm(mu_post, np.sqrt(sigma_post)).pdf(x)
    
ax = plt.subplot()
x = np.linspace(-1, 1, 500)
posterior_analytical = calc_posterior_analytical(data, x, 0., 1.)
ax.plot(x, posterior_analytical)
ax.set(xlabel = 'mu', ylabel = 'belief', title='Analytical posterior')
sns.despine()


## ################################
## MCMC sampling
##
## Goal:
##   get posterior of mu
## 
## approximate posterior analytical
##
## ################################

def sampler(data, samples = 4, mu_init=.5, proposal_width = .5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
   mu_current = mu_init
   posterior = [mu_current]
   for i in range(samples):
       # sugget new position
       #
       # i-th state(xi) of markov: mu_current
       # proposal: sample from q(x|xi)
       # 
       # sample u from uniform(0, 1)
       # if u < alpha(xi, x) = min(p_proposal/p_current, 1)
       #     accept!!
       #
       mu_proposal = norm(mu_current, proposal_width).rvs() # random choose
       
       # Compute likelihood by multiplying probabilities of each data point
       likelihood_current = norm(mu_current, 1).pdf(data).prod()
       likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
       
       # Compute prior probability of current and proposed mu 
       prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
       prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
       
       p_current = likelihood_current * prior_current
       p_proposal = likelihood_proposal * prior_proposal
       
       # Accept proposal?
       p_accept = p_proposal / p_current
       
       # Usually would include prior probability, which we neglect here for simplicity
       accept = np.random.rand() < p_accept

   #    if plot:
    #       plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i)
           
       if accept:
           mu_current = mu_proposal
          
       posterior.append(mu_current)
   return posterior
## ############################################  
## function to display
## http://twiecki.github.io/blog/2015/11/10/mcmc-sampling/
## ############################################    

posterior = sampler(data, samples=5000, mu_init=.1, plot=False)
fig, ax = plt.subplots()
ax.plot(posterior)
_ = ax.set(xlabel='sample', ylabel='mu')

ax = plt.subplot()

sns.distplot(posterior[500:], ax=ax, label='estimated posterior')

x = np.linspace(-.5, .5, 500)
post = calc_posterior_analytical(data, x, 0, 1)
ax.plot(x, post, 'g', label='analytic posterior')
_ = ax.set(xlabel='mu', ylabel='belief');
ax.legend();

## ##########################################
## by following the above procedure, we get samples from the same distribution 
## as what we derived analytically.
## ##########################################

## ##########################################
## Proposal width
## ##########################################

posterior_small = sampler(data, samples=5000, mu_init=1., proposal_width=.01)
fig, ax = plt.subplots()
ax.plot(posterior_small);
_ = ax.set(xlabel='sample', ylabel='mu');

posterior_large = sampler(data, samples=5000, mu_init=1., proposal_width=3.)
fig, ax = plt.subplots()
ax.plot(posterior_large); plt.xlabel('sample'); plt.ylabel('mu');
_ = ax.set(xlabel='sample', ylabel='mu');

sns.distplot(posterior_small[1000:], label='Small step size')
sns.distplot(posterior_large[1000:], label='Large step size');
_ = plt.legend();

## ##########################################
## one common metric to evaluate the efficiency of our sampler is the autocorrelation
## ##########################################

from pymc3.stats import autocorr
lags = np.arange(1, 100)
fig, ax = plt.subplots()
ax.plot(lags, [autocorr(posterior_large, l) for l in lags], label='large step size')
ax.plot(lags, [autocorr(posterior_small, l) for l in lags], label='small step size')
ax.plot(lags, [autocorr(posterior, l) for l in lags], label='medium step size')
ax.legend(loc=0)
_ = ax.set(xlabel='lag', ylabel='autocorrelation', ylim=(-.1, 1))

## ###########################################
## pymc3
## ###########################################

import pymc3 as pm

with pm.Model():
    mu = pm.Normal('mu', 0, 1)
    sigma = 1.
    returns = pm.Normal('returns', mu=mu, sd=sigma, observed=data)
    
    step = pm.Metropolis()
    trace = pm.sample(15000, step)
    
sns.distplot(trace[2000:]['mu'], label='PyMC3 sampler');
sns.distplot(posterior[500:], label='Hand-written sampler');
plt.legend();
