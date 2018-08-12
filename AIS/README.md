# Adaptive Importance Sampling

## One Simple Way

1. start with a trial density, say $$g_0(x) = t_\alpha (x; \mu_0,\Sigma_0)$$
2. with weighted Monte Carlo samples, estimate the parameters, mean and covariate matrix, and construct new trial density, say $g_1(x) = t_\alpha (x; \mu_1,\Sigma_1)$
3. construct a certain measure of discrepancy between the trial distribution and the target distribution, such as the coefficient of variation of importance weights, does not improve any more.

## One Example
Implement an Adaptive Importance Sampling algorithm to evaluate mean and variance of a density

$$
\pi(\mathbf{x})\propto N(\mathbf{x;0}, 2I_4) + 2N(\mathbf{x; 3e}, I_4) + 1.5 N(\mathbf{x; -3e}, D_4)
$$

where $\mathbf{e} = (1,1,1,1), I_4 = diag(1,1,1,1), D_4 = diag(2,1,1,.5)$.

A possible procedure is as follows:

1. start with a trial density $g_0 =t_{\nu}(0, \Sigma)$
2. Recursively, build
$$
g_k(\mathbf{x})=(1-\epsilon)g_{k-1}(\mathbf {x}) + \epsilon t_{\nu}(\mu, \Sigma)
$$
in which one chooses $(\epsilon, \mu, \Sigma)$ to minimize the variation of coefficient of the importance weights.
