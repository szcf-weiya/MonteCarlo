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