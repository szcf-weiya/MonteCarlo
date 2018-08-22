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