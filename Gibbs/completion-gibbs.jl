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