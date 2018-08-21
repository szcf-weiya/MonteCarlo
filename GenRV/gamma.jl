## Julia program for Generation of Gamma Random Variable
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-21

## sample from Ga(a, beta)
function rgamma_int(a::Int, beta)
    u = rand(a)
    return(-1.0 * beta * sum(log.(u)))
end

## density of Ga(alpha, beta)
include("func_gamma.jl")
function dgamma(x, alpha, beta)
    return(beta^alpha / lanczos_gamma(alpha) * x^(alpha-1) * exp(-1*beta*x))
end

## accept-reject algorithm
function rgamma(alpha = 1.5, beta = 2.1)
    a = Int(floor(alpha))
    b = a * beta / alpha
    M = exp(a * (log(a) - 1) - alpha * (log(alpha) - 1))
    while true
        x = rgamma_int(a, b)
        u = rand()
        cutpoint = dgamma(x, alpha, beta)/(M*dgamma(x, a, b))
        if u <= cutpoint
            return x
        end
    end
end

# example
println(rgamma())

