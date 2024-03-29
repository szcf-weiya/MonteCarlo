# MCMC

[Robert and Casella (2013)](https://www.springer.com/gp/book/9781475730715) provides the following definition:

![](mcmc-def.png)

For illustration, you can visit the jupyter notebook on [MCMC_example](http://nbviewer.jupyter.org/github/szcf-weiya/MonteCarlo/blob/master/MCMC/MCMC_example.ipynb)

Working principle: For an arbitrary starting value $$x^{(0)}$$, a chain $$X^{(t)}$$ is generated using a transition kernel with stationary distribution $$f$$, which ensures the convergence in distribution of $$X^{(t)}$$ to a random variable from $$f$$.

Some good materials about MCMC.

1. https://cosx.org/2013/01/lda-math-mcmc-and-gibbs-sampling
