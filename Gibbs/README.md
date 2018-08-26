# Gibbs Sampler

[Robert and Casella (2013)](https://www.springer.com/gp/book/9781475730715) gives the following algorithm,

![](gibbs.png)

while [Liu, Jun S. (2008)](https://www.springer.com/gp/book/9780387763699) introduces two types of Gibbs sampling strategy.

![](two-types-gibbs.png)

## Bivariate Gibbs sampler

![](bi-gibbs.png)

It is easy to implement this sampler:

```julia
## Julia program for Bivariate Gibbs sampler
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-22

function bigibbs(T, rho)
    x = ones(T+1)
    y = ones(T+1)
    for t = 1:T
        x[t+1] = randn() * sqrt(1-rho^2) + rho*y[t]
        y[t+1] = randn() * sqrt(1-rho^2) + rho*x[t+1]
    end
    return x, y
end

## example
bigibbs(100, 0.5)
```

## Completion Gibbs Sampler

![](c-gibbs.png)

Example:

![](ex-7-1-5.png)

We can use the following Julia program to implement this algorithm.

```julia
## Julia program for Truncated normal distribution
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-22

# Truncated normal distribution
function rtrunormal(T, mu, sigma, mu_down)
    x = ones(T)
    z = ones(T+1)
    # set initial value of z
    z[1] = rand()
    if mu < mu_down
        z[1] = z[1] * exp(-0.5 * (mu - mu_down)^2 / sigma^2)
    end
    for t = 1:T
        x[t] = rand() * (mu - mu_down + sqrt(-2*sigma^2*log(z[t]))) + mu_down
        z[t+1] = rand() * exp(-(x[t] - mu)^2/(2*sigma^2))
    end
    return(x)
end

## example
rtrunormal(1000, 1.0, 1.0, 1.2)
```

## Slice Sampler

[Robert and Casella (2013)](https://www.springer.com/gp/book/9781475730715) introduces the following slice sampler algorithm,

![](slice.png)

and [Liu, Jun S. (2008)](https://www.springer.com/gp/book/9780387763699) also presents the slice sampler with slightly different expression:

![](slice-lj.png)

In my opinion, we can illustrate this algorithm with one dimensioanl case. Suppose we want to sample from normal distribution (or uniform distribution), we can sample uniformly from the region encolsed by the coordinate axis and the density function, that is a bell shape (or a square).

Consider the normal distribution as an instance.

![](ex-7-1-7.png)

It is also easy to write the following Julia program.

```julia
## Julia program for Slice sampler
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-22

function rnorm_slice(T)
    x = ones(T+1)
    w = ones(T+1)
    for t = 1:T
        w[t+1] = rand() * exp(-1.0 * x[t]^2/2)
        x[t+1] = rand() * 2 * sqrt(-2*log(w[t+1])) - sqrt(-2*log(w[t+1]))
    end
    return x[2:end]
end

## example
rnorm_slice(100)
```

## Data Augmentation

A special case of Completion Gibbs Sampler.

![](dataaug.png)

Let's illustrate the scheme with grouped counting data.

![](ex-7-2-2.png)

And we can obtain the following algorithm,

![](pois_gamma_gibbs.png)

But it seems to be not obvious to derive the above algorithm, so I wrote some more details

![](details_alg35_fix.jpg)

Another example:

![](ex-7-2-3.png)

It seems that we do not need to derive the explicit form of $$g(x, z)$$, if we can directly obtain the conditional distribution. We can use the following Julia program to sample.

```julia
## Julia program for Grouped Multinomial Data (Ex. 7.2.3)
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-26

# call gamma function
#using SpecialFunctions


# sample from Dirichlet distributions
using Distributions

function gmulti(T, x, a, b, alpha1 = 0.5, alpha2 = 0.5, alpha3 = 0.5)
    z = ones(T+1, size(x, 1)-1) # initial z satisfy `z <= x`
    mu = ones(T+1)
    eta = ones(T+1)
    for t = 1:T
        # sample from g_1(theta | y)
        dir = Dirichlet([z[t, 1] + z[t, 2] + alpha1, z[t, 3] + z[t, 4] + alpha2, x[5] + alpha3])
        sample = rand(dir, 1)
        mu[t+1] = sample[1]
        eta[t+1] = sample[2]
        # sample from g_2(z | x, theta)
        for i = 1:2
            bi = Binomial(x[i], a[i]*mu[t+1]/(a[i]*mu[t+1]+b[i]))
            z[t+1, i] = rand(bi, 1)[1]
        end 
        for i = 3:4
            bi = Binomial(x[i], a[i]*eta[t+1]/(a[i]*eta[t+1]+b[i]))
            z[t+1, i] = rand(bi, 1)[1]
        end
    end
    return mu, eta
end

# example

## data
a = [0.06, 0.14, 0.11, 0.09];
b = [0.17, 0.24, 0.19, 0.20];
x = [9, 15, 12, 7, 8]; 

gmulti(100, x, a, b)
```

## Reversible Data Augmentation

![](rdataaug.png)

## Reversible Gibbs Sampler

![](rgibbs.png)

## Random Sweep Gibbs Sampler

![](rsgibbs.png)

## Random Gibbs Sampler

![](randgibbs.png)

## Hybrid Gibbs Samplers

![](hybridMCMC.png)

## Metropolization of the Gibbs Sampler

![](metrogibbs.png)

Let us illuatrate this algorithm with the following example.

![](ex-7-2-3.png)

![](ex-7-3-8.png)


