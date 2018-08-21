## Julia program for Gamma Metropolis-Hastings
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-21

## import function gamma_int
include("../GenRV/gamma.jl")

function mh_gamma(T = 100, alpha = 1.5)
    a = Int(floor(alpha))
    b = a/alpha
    x = ones(T+1) # initial value: 0 
    for t = 1:T
        yt = rgamma_int(a, 1)
        rt = (yt / x[t] * exp((x[t] - yt) / alpha))^(alpha-a)
        if rt >= 1
            x[t+1] = yt
        else
            u = rand()
            if u < rt
                x[t+1] = yt
            else
                x[t+1] = x[t]
            end
        end   
    end
    return(x)
end

# example
mh_gamma()