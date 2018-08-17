# Metropolis-Hastings Algorithm

Gibbs and Metropolis algoritgms are special cases of Metropolis-Hastings.

Consider a bivariate distribution for two random variables $$U$$ and $$V$$, $$p_0(u,v)$$.

## Gibbs

In Gibbs, given $$x^{(s)}=(u^{(s)},v^{(s)})$$, sample $$x^{s+1}$$ as follows.

1. Update $$U$$: sample $$u^{(s+1)}\sim p_0(u\mid v^{(s)})$$
2. Update $$V$$: sample $$v^{(s+1)}\sim p_0(v\mid u^{(s+1)})$$

Of course, we can change the sampling order.

## Metropolis-Hastings

No need to require the acceptance ratio to be symmetric.

1. Update $$U$$:
    1. sample $$u^*\sim J_u(u\mid u^{(s)}, v^{(s)})$$
    2. compute the acceptance ratio
    $$
    r=\frac{p_0(u^*, v^{(s)})}{p_0(u^{(s)},v^{(s)})}\times \frac{J_u(u^{s}\mid u^*, v^{(s)})}{J_u(u^*\mid u^{(s)},v^{(s)})}
    $$
    3. set $$u^{(s+1)}$$ to $$u^*$$ w.p. $$\mathrm{min}(1,r)$$.
2. Update $$V$$: similarly.

