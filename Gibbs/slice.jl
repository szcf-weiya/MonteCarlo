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