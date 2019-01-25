## Essential properties of Markov Chain

In the setup of MCMC algorithms, we construct Markov chains from a transition kernel $$K$$, a conditional probability density such that $$X_{n+1}\sim K(X_n,X_{n+1})$$.

The chain encountered in MCMC settings enjoy a very strong stability property, namely a **stationary probability distribution**; that is, a distribution $$\pi$$ such that if $$X_n\sim\pi$$, then $$X_{n+1}\sim \pi$$, if the kernel $$K$$ allows for free moves all over the state space. This freedom is called **irreducibility** in the theory of Markov chains and is formalized as the existence of $$n\in\mathbb{N}$$ such that $$P(X_n\in A\mid X_0)>0$$ for every $$A$$ such that $$\pi(A)>0$$. This property also ensures that most of the chains involved in MCMC algorithms are **recurrent** (that is, that the average number of visits to an arbitrary set $A$ is infinite), or even **Harris recurrent** (that is, such that the probability of an infinite number of returns to $A$ is 1).

**Harris recurrence** ensures that the chain has the same limiting behavior for every starting value instead of almost every starting value.

The stationary distribution is also a limiting distribution in the sense that the limiting distribution of $$X_{n+1}$$ is $$\pi$$ under the total variation norm, notwithstanding the initial value of $$X_0$$.

Strong forms of convergence are also encountered in MCMC settings, like geometric and uniform convergences. 

If the marginals are proper, for convergence we only need our chain to be aperiodic. A sufficient condition is  that $$K(x_n,\cdot)>0$$ (or, equivalently, $$f(\cdot\mid x_n)>0$$) in a neighborhood of $$x_n$$.

If the marginal are not proper, or if they do not exist, then the chain is not positive recurrent. It is either null recurrent, and both cases are bad.

The [detailed balance condition](detailed_balance_condition.png) is **not necessary** for $f$ to be a stationary measure associated with the transition kernel $$K$$, but it provides a sufficient condition that is often easy to check and that can be used for most MCMC algorithms.

Ergodicity: independence of initial conditions

- [geometrically h-ergodic](def_geometrical_h_ergodic.png): decreasing at least at a geometric speed.
- [uniform ergodicity](def_uniform_ergodic.png): stronger than geometric ergodicity in the sense that the rate of geometric convergence must be uniform over the whole space.


## Basic Notions

![](def6.2.png)


![](def6.4.png)

## Bernoulli-Laplace Model

![](ex6.3.png)

The finite chain is indeed irreducible since it is possible to connect the status $$x$$ and $$y$$ in $$\vert x-y\vert$$ steps with probability

$$
\prod_{i=x\land y}^{x\vee y-1}\Big(\frac{M-i}{M}\Big)^2\,.
$$

The Bernoulli-Laplace chain is [aperiodic](def_period_discrete.png) and even [strongly aperiodic](def_period.png) since the diagonal terms satisfy $$P_{xx}>0$$ for every $$x\in \{0,\ldots,K\}$$.

Given the quasi-diagonal shape of the transition matrix, it is possible to directly determine the invariant distribution, $$\pi=(\pi_0,\ldots,\pi_K)$$. From the equation $$\pi P = \pi$$,

$$
\begin{aligned}
\pi_0 &= P_{00}\pi_0 + P_{10}\pi_1\\
\pi_1 &= P_{01}\pi_1 + P_{11}\pi_1 + P_{21}\pi_2\\
\cdots &=\cdots\\
\pi_K &= P_{(K-1)K}\pi_{K-1} + P_{KK}\pi_K\,.
\end{aligned}
$$

Thus,

$$
\pi_k=\binom{K}{k}^2\pi_0\,,\qquad k=0,\ldots,K\,,
$$

and through normalization,

$$
\pi_k=\frac{\binom{K}{k}^2}{\binom{2K}{K}}\,,
$$

by using [Chu-Vandermonde identity](https://en.wikipedia.org/wiki/Vandermonde%27s_identity)

$$
\binom{m+n}{r}=\sum_{k=0}^r\binom{m}{k}\binom{n}{r-k}
$$

with $$m=n=r=K$$. It turns out that the hypergeometric distribution $H(2K,K,1/2)$ is the [invariant distribution](def_invariant.png) for the Bernoulli-Laplace model.



## AR(1) Models

A simple illustration of Markov chains on continuous state-space. 

$$
X_n = \theta X_{n-1}+\varepsilon_n\;\theta\in \mathrm{I\!R}\,,
$$

with $$\varepsilon_n\in N(0,\sigma^2)$$, and if the $$\varepsilon_n$$'s are independent, $$X_n$$ is indeed independent from $$X_{n-2},X_{n-3},\ldots$$ conditionally on $$X_{n-1}$$.

- The Markovian properties of an AR(q) can be derived from $$(X_n,\ldots,X_{n-q+1})$$.
- ARMA(p, q) doesn't fit in the Markovian framework.

Since $$X_n\mid x_{n-1}\sim N(\theta x_{n-1},\sigma^2)$$, consider the lower bound of the transition kernel ($$\theta > 0$$):

$$
\begin{aligned}
K(x_{n-1},x_n) &= \frac{1}{\sqrt{2\pi}}\exp\Big\{-\frac{1}{2\sigma^2}(x_n-\theta x_{n-1})^2\Big\}\\
&\ge \frac{1}{\sqrt{2\pi}\sigma}\exp\Big\{-\frac{1}{2\sigma^2}\max\{(x_n-\theta \underline w)^2, (x_n-\theta \bar w)^2\}\Big\}\\
&\ge \frac{1}{\sqrt{2\pi}\sigma}\exp\Big\{-\frac{1}{2\sigma^2}\Big[ \max\{-2\theta x_n\underline w,-2\theta x_n\bar w\}+x_n^2 + \theta^2\underline w^2\land \bar w^2 \Big]\Big\}\\
&=\begin{cases}
\frac{1}{\sqrt{2\pi}\sigma}\exp\Big\{-\frac{1}{2\sigma^2}\Big[ -2\theta x_n\underline w+x_n^2 + \theta^2\underline w^2\land \bar w^2 \Big]\Big\}& \text{if }x_n>0\\
\frac{1}{\sqrt{2\pi}\sigma}\exp\Big\{-\frac{1}{2\sigma^2}\Big[ -2\theta x_n\bar w+x_n^2 + \theta^2\underline w^2\land \bar w^2 \Big]\Big\}& \text{if }x_n<0
\end{cases}\,,
\end{aligned}
$$

when $$x_{n-1}\in[\underline w, \bar w]$$. The set $$C = [\underline w, \bar w]$$ is a [small set](def_small.png), as the measure $$\nu_1$$ with density 
$$
\frac{\exp\{(-x^2+2\theta x\underline w)/2\sigma^2\}I_{x>0} + \exp\{(-x^2+2\theta x\bar w)/2\sigma^2\}I_{x<0} }{\sqrt{2\pi}\sigma\{[1-\Phi(-\theta\underline w/\sigma^2)]\exp(\theta^2\underline w^2/2\sigma^2)+\Phi(-\theta\bar w/\sigma)\exp(\theta^2\bar w^2/2\sigma^2)\}}\,,
$$
satisfy 
$$
K^1(x,A)\ge \nu_1(A),\qquad \forall x\in C, \forall A\in {\cal B(X)}\,.
$$

Given that the transition kernel corresponds to the $N(\theta x_{n-1},\sigma^2)$ distribution, a normal distribution $N(\mu,\tau^2)$ is stationary for the AR(1) chain only if 

$$
\mu=\theta\mu\qquad\text{and}\qquad \tau^2=\tau^2\theta^2+\sigma^2\,.
$$

These conditions imply that $$\mu=0$$ and that $$\tau^2=\sigma^2/(1-\theta^2)$$, which can only occur for $$\vert \theta\vert < 1$$. In this case, $$N(0,\sigma^2/(1-\theta^2))$$ is indeed the unique stationary distribution of the AR(1) chain.



## Branching process

![](ex_branching_process.png)

- If $$\phi$$ doesn't have a constant term, i.e., $$P(X_1=0)=0$$, then chain $$S_t$$ is necessarily transient since it is increasing.
- If $$P(X_1=0)>0$$, the probability of a return to 0 at time $$t$$ is $$\rho(t)=P(S_t=0)=g_t(0)$$, which thus satisfies the recurrence equation $$\rho_t=\phi(\rho_{t-1})$$. There exists a limit $$\rho$$ different from 1, such that $$\rho=\phi(\rho)$$, iff $$\phi'(1)>1$$; namely if $$E[X]>1$$. The chain is thus transient when the average number of siblings per individual is larger than 1. If there exists a restarting mechanism in 0, $$S_{t+1}\mid S_t=0\sim\phi$$, it is easily shown that when $$\phi'(1)>1$$, the number of returns to 0 follows a geometric distribution with parameter $$\rho$$.

![](i1_branch_process.png)

- If $$\phi'(1)\le 1$$, the chain is recurrent.

